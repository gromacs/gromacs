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
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/embedded_system_preprocessing.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/vectypes.h"

#include "nnpotmodel.h"
#include "nnpotoptions.h"
#include "torchmodel.h"

/*! \internal Helper function to find the index of a value in a vector
 *
 * Returns the index of the first occurrence of \p val in \p vec, or -1 if not found.
 * Can be used as an inversion for the index lookup tables.
 * \param[in] vec vector to search in
 * \param[in] val value to search for
 */
static std::optional<ptrdiff_t> indexOf(gmx::ArrayRef<const int> vec, const int val)
{
    auto it = std::find(vec.begin(), vec.end(), val);
    if (it == vec.end())
    {
        return std::nullopt;
    }
    return std::distance(vec.begin(), it);
}

namespace gmx
{

NNPotForceProvider::NNPotForceProvider(const NNPotParameters& nnpotParameters,
                                       const MDLogger&        logger,
                                       const MpiComm&         mpiComm) :
    params_(nnpotParameters),
    positions_(params_.numAtoms_, RVec({ 0.0, 0.0, 0.0 })),
    atomNumbers_(params_.numAtoms_, -1),
    inputToLocalIndex_(params_.numAtoms_, -1),
    box_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
    logger_(logger),
    mpiComm_(mpiComm)
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

    // set communication record
    model_->setComm(mpiComm_);
}

NNPotForceProvider::~NNPotForceProvider() {}

void NNPotForceProvider::calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput)
{
    // make sure inputs are available
    if (params_.modelNeedsInput("atom-positions"))
    {
        gatherAtomPositions(fInput.x_);
    }
    // copy box
    copy_mat(fInput.box_, box_);

    // get link atom info
    std::vector<LinkFrontierAtom> linkFrontier(constructLinkFrontier(params_.linkFrontier_));

    // prepare inputs for NN model
    model_->evaluateModel(&(fOutput->enerd_),
                          fOutput->forceWithVirial_.force_,
                          inputToLocalIndex_,
                          params_.modelInput_,
                          positions_,
                          atomNumbers_,
                          linkFrontier,
                          &box_,
                          params_.pbcType_.get());
}

void NNPotForceProvider::gatherAtomNumbersIndices(const MDModulesAtomsRedistributedSignal& signal)
{
    // create lookup table for local atom indices needed for hybrid ML/MM
    // -1 is used as a flag for atoms that are not local / part of the NN input
    // used to distribute forces to correct local indices as the NN input tensor does not contain all atoms
    const int numInput = params_.nnpAtoms_->numAtomsGlobal();

    inputToLocalIndex_.assign(numInput, -1);
    atomNumbers_.assign(numInput, 0);

    if (mpiComm_.isParallel())
    {
        GMX_RELEASE_ASSERT(signal.globalAtomIndices_.has_value(),
                           "Global atom indices must be provided when using domain decomposition.");
        auto      globalAtomIndices = signal.globalAtomIndices_.value();
        const int numLocalPlusHalo  = globalAtomIndices.size(); // includes halo atoms
        const int numLocal          = signal.x_.size();         // only local atoms on this rank
        localToInputIndex_.assign(numLocalPlusHalo, -1);

        for (int i = 0; i < numLocalPlusHalo; i++)
        {
            int globalIdx = globalAtomIndices[i];
            for (int j = 0; j < numInput; j++)
            {
                if (params_.nnpAtoms_->globalIndex()[j] == globalIdx)
                {
                    // only map input indices for home atoms
                    if (i < numLocal)
                    {
                        inputToLocalIndex_[j] = i;
                        atomNumbers_[j]       = params_.atoms_.atom[globalIdx].atomnumber;
                    }
                    localToInputIndex_[i] = j;
                    break;
                }
            }
        }
        mpiComm_.sumReduce(numInput, atomNumbers_.data());
    }
    else
    {
        localToInputIndex_.assign(params_.numAtoms_, -1);
        for (int i = 0; i < numInput; i++)
        {
            int localIndex = params_.nnpAtoms_->localIndex()[i];
            int globalIdx = params_.nnpAtoms_->globalIndex()[params_.nnpAtoms_->collectiveIndex()[i]];
            inputToLocalIndex_[i]          = localIndex;
            localToInputIndex_[localIndex] = i;
            atomNumbers_[i]                = params_.atoms_.atom[globalIdx].atomnumber;
        }
    }
    // sanity check: make sure all atom numbers have been set
    GMX_RELEASE_ASSERT(std::count(atomNumbers_.begin(), atomNumbers_.end(), 0) == 0,
                       "Some atom numbers have not been set correctly.");
}

void NNPotForceProvider::gatherAtomPositions(ArrayRef<const RVec> pos)
{
    // collect atom positions
    // at this point, we already have the atom numbers and indices, so we can fill the positions
    size_t numInput = inputToLocalIndex_.size();

    // reset positions to zero, because we might not have all atoms in the input
    positions_.assign(numInput, RVec({ 0.0, 0.0, 0.0 }));

    for (size_t i = 0; i < numInput; i++)
    {
        // if value in lookup table is -1, the atom is not local to this rank
        if (inputToLocalIndex_[i] != -1)
        {
            positions_[i] = pos[inputToLocalIndex_[i]];
        }
    }

    // in case of dom dec, distribute positions to all ranks
    if (mpiComm_.isParallel())
    {
        mpiComm_.sumReduce(3 * numInput, positions_.data()->as_vec());
    }
}

std::vector<LinkFrontierAtom>
NNPotForceProvider::constructLinkFrontier(const std::vector<LinkFrontierAtom>& inputLinkFrontier)
{
    std::vector<LinkFrontierAtom> linkFrontier;
    linkFrontier.reserve(inputLinkFrontier.size());
    // the MM atom of the link is part of the input
    // we replace it with the link atom
    for (auto& ilink : inputLinkFrontier)
    {
        const int        globalIdxMM  = ilink.getMMIndex();
        const int        globalIdxNNP = ilink.getEmbeddedIndex();
        LinkFrontierAtom link(globalIdxNNP, globalIdxMM);

        // save input indices for force redistribution later
        const std::optional<ptrdiff_t> inputIdxMM  = indexOf(inputToGlobalIndex_, globalIdxMM);
        const std::optional<ptrdiff_t> inputIdxNNP = indexOf(inputToGlobalIndex_, globalIdxNNP);
        GMX_RELEASE_ASSERT(inputIdxNNP.has_value() && inputIdxMM.has_value(),
                           "Link frontier atoms must be part of the NNPot atom set");
        link.setInputIndices(inputIdxNNP.value(), inputIdxMM.value());

        // calculate and set link atom position
        const RVec posMM  = positions_[inputIdxMM.value()];
        const RVec posNNP = positions_[inputIdxNNP.value()];
        link.setPositions(posNNP, posMM);
        // update position of MM atom to link atom position
        positions_[inputIdxMM.value()] = link.getLinkPosition();

        // overwrite atomic number of MM atom
        atomNumbers_[inputIdxMM.value()] = ilink.linkAtomNumber();
        linkFrontier.push_back(link);
    }
    return linkFrontier;
}

} // namespace gmx
