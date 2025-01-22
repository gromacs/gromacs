/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements the NNPot Force Provider class
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "nnpotforceprovider.h"

#include <filesystem>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/logger.h"

#include "nnpotmodel.h"
#include "nnpotoptions.h"
#include "torchmodel.h"

namespace gmx
{

NNPotForceProvider::NNPotForceProvider(const NNPotParameters& nnpotParameters, const MDLogger* logger) :
    params_(nnpotParameters),
    positions_(params_.numAtoms_, RVec({ 0.0, 0.0, 0.0 })),
    atomNumbers_(params_.numAtoms_, -1),
    idxLookup_(params_.numAtoms_, -1),
    box_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
    logger_(logger),
    cr_(params_.cr_)
{
    // initialize the neural network model
    std::filesystem::path modelPath(params_.modelFileName_);
    if (!std::filesystem::exists(modelPath))
    {
        GMX_THROW(FileIOError("Model file does not exist: " + params_.modelFileName_));
    }
    else if (modelPath.extension() == ".pt")
    {
        model_ = std::make_shared<TorchModel>(params_.modelFileName_, logger_);
    }
    else
    {
        GMX_THROW(FileIOError("Unrecognized extension for model file: " + params_.modelFileName_));
    }
    model_->initModel();
}

NNPotForceProvider::~NNPotForceProvider() {}

void NNPotForceProvider::calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput)
{
    // store a pointer to the communication record
    cr_ = &(fInput.cr_);
    model_->setCommRec(cr_);

    // prepare inputs for NN model
    // order in input vector is the same as in mdp file
    for (const std::string& input : params_.modelInput_)
    {
        if (input.empty())
        {
            continue;
        }
        else if (input == "atom-positions")
        {
            gatherAtomPositions(fInput.x_);
            model_->prepareAtomPositions(positions_);
        }
        else if (input == "atom-numbers")
        {
            model_->prepareAtomNumbers(atomNumbers_);
        }
        else if (input == "box")
        {
            copy_mat(fInput.box_, box_);
            model_->prepareBox(box_);
        }
        else if (input == "pbc")
        {
            copy_mat(fInput.box_, box_);
            t_pbc pbc;
            set_pbc(&pbc, *(params_.pbcType_), box_); // might not be necessary
            model_->preparePbcType(*(params_.pbcType_));
        }
        else
        {
            GMX_THROW(InconsistentInputError("Unknown input to NN model: " + input));
        }
    }

    model_->evaluateModel();

    model_->getOutputs(idxLookup_, fOutput->enerd_, fOutput->forceWithVirial_.force_);
}

void NNPotForceProvider::gatherAtomNumbersIndices()
{
    // this might not be the most efficient solution, since we are throwing away most of the
    // vectors here in case of NNP/MM

    // create lookup table for local atom indices needed for hybrid ML/MM
    // -1 is used as a flag for atoms that are not local / not in the input
    // used to distribute forces to correct local indices as the NN input tensor does not contain all atoms
    idxLookup_.assign(params_.numAtoms_, -1);
    atomNumbers_.assign(params_.numAtoms_, 0);

    int lIdx, gIdx;
    for (size_t i = 0; i < params_.inpAtoms_->numAtomsLocal(); i++)
    {
        lIdx = params_.inpAtoms_->localIndex()[i];
        gIdx = params_.inpAtoms_->globalIndex()[params_.inpAtoms_->collectiveIndex()[i]];
        // TODO: make sure that atom number indexing is correct
        atomNumbers_[gIdx] = params_.atoms_.atom[gIdx].atomnumber;
        idxLookup_[gIdx]   = lIdx;
    }

    // distribute atom numbers to all ranks
    if (havePPDomainDecomposition(cr_))
    {
        gmx_sumi(params_.numAtoms_, atomNumbers_.data(), cr_);
    }

    // remove unused elements in atomNumbers_, and idxLookup
    auto atIt  = atomNumbers_.begin();
    auto idxIt = idxLookup_.begin();
    while (atIt != atomNumbers_.end() && idxIt != idxLookup_.end())
    {
        if (*atIt == 0)
        {
            atIt  = atomNumbers_.erase(atIt);
            idxIt = idxLookup_.erase(idxIt);
        }
        else
        {
            ++atIt;
            ++idxIt;
        }
    }
}

void NNPotForceProvider::gatherAtomPositions(ArrayRef<const RVec> pos)
{
    // collect atom positions
    // at this point, we already have the atom numbers and indices, so we can fill the positions
    size_t numInput = idxLookup_.size();

    // reset positions to zero, because we might not have all atoms in the input
    positions_.assign(numInput, RVec({ 0.0, 0.0, 0.0 }));

    for (size_t i = 0; i < numInput; i++)
    {
        // if value in lookup table is -1, the atom is not local to this rank
        if (idxLookup_[i] != -1)
        {
            positions_[i] = pos[idxLookup_[i]];
        }
    }

    // in case of dom dec, distribute positions to all ranks
    if (havePPDomainDecomposition(cr_))
    {
        gmx_sum(3 * numInput, positions_.data()->as_vec(), cr_);
    }
}

} // namespace gmx
