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
#include "gromacs/pbcutil/ishift.h"
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

/*! \brief Helper function to gather the positions of the NNP atoms from the input.
 *
 * \param[in] mmAtoms set of MM atoms
 * \param[in] positions positions of all atoms
 * \param[in] charges charges of all atoms
 * \returns positions and charges of the MM atoms
 */
static std::tuple<std::vector<RVec>, std::vector<real>>
getMMPositionsAndCharges(const std::unique_ptr<LocalAtomSet>& mmAtoms,
                         ArrayRef<const RVec>                 positions,
                         ArrayRef<const real>                 charges)
{
    // collect MM atom positions and charges
    // can't use lookup table here, because we need all MM atoms
    size_t numMM = mmAtoms->numAtomsLocal();

    // resize positions and charges vectors
    std::vector<RVec> mmPositions(numMM);
    std::vector<real> mmCharges(numMM);

    for (size_t i = 0; i < numMM; i++)
    {
        mmPositions[i] = positions[mmAtoms->localIndex()[i]];
        mmCharges[i]   = charges[mmAtoms->localIndex()[i]];
    }
    return std::make_tuple(mmPositions, mmCharges);
}

/*! Center positions of NN and MM atoms in the box.
 *
 * For treatment of the NNP-MM interactions, we assume that the embedding model expects
 * all atom positions to be centered around the NNP region.
 */
static void centerAtomPositions(ArrayRef<RVec> nnPositions,
                                ArrayRef<RVec> mmPositions,
                                const matrix&  box,
                                const PbcType& pbcType)
{
    t_pbc pbc;
    set_pbc(&pbc, pbcType, box);

    // center atom positions in the box around the NNP region
    // compute center of the NNP region
    RVec nnpCenter{ 0.0, 0.0, 0.0 };
    RVec dx;
    for (const auto& pos : nnPositions)
    {
        pbc_dx(&pbc, pos, nnPositions[0], dx);
        nnpCenter += dx;
    }
    nnpCenter        = nnpCenter / nnPositions.size() + nnPositions[0];
    RVec translation = RVec(box[0]) + RVec(box[1]) + RVec(box[2]);
    translation /= 2.0;
    translation -= nnpCenter;

    // apply translation to NNP and MM positions
    for (auto& pos : nnPositions)
    {
        pos += translation;
    }
    for (auto& pos : mmPositions)
    {
        pos += translation;
    }
    // put all atoms into the central box (they might be shifted out of it because of the translation)
    put_atoms_in_box(pbcType, box, nnPositions);
    put_atoms_in_box(pbcType, box, mmPositions);
}

NNPotForceProvider::NNPotForceProvider(const NNPotParameters& nnpotParameters,
                                       const MDLogger&        logger,
                                       const MpiComm&         mpiComm) :
    params_(nnpotParameters), inputToLocalIndex_(params_.numAtoms_, -1), logger_(logger), mpiComm_(mpiComm)
{
    // for now, dom dec is disabled with pair list input
    if (mpiComm_.isParallel()
        && (params_.modelNeedsInput("atom-pairs") || params_.modelNeedsInput("pair-shifts")))
    {
        GMX_THROW(NotImplementedError(
                "NNPot with pair list input does not support domain decomposition"));
    }
    if (params_.embeddingScheme_ == NNPotEmbedding::ElectrostaticModel && mpiComm_.isParallel())
    {
        GMX_THROW(NotImplementedError(
                "Electrostatic embedding scheme not yet implemented for domain decomposition."));
    }

    // initialize the neural network model
    std::filesystem::path modelPath(params_.modelFileName_);
    if (!std::filesystem::exists(modelPath))
    {
        GMX_THROW(FileIOError("Model file does not exist: " + params_.modelFileName_));
    }
    else if (modelPath.extension() == ".pt")
    {
        model_ = std::make_shared<TorchModel>(
                params_.modelFileName_, params_.embeddingScheme_, logger_, mpiComm_);
    }
    else
    {
        GMX_THROW(FileIOError("Unrecognized extension for model file: " + params_.modelFileName_));
    }
}

NNPotForceProvider::~NNPotForceProvider() {}

void NNPotForceProvider::calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput)
{
    // make sure inputs are available
    std::vector<RVec> positions;
    std::vector<RVec> mmPositions;
    std::vector<real> mmCharges;
    std::vector<int>  mmIndices;
    matrix            box;
    copy_mat(fInput.box_, box);
    if (params_.modelNeedsInput("atom-positions"))
    {
        positions = gatherAtomPositions(fInput.x_);
    }
    if (params_.modelNeedsInput("atom-positions-mm") || params_.modelNeedsInput("atom-charges-mm"))
    {
        mmIndices.assign(params_.mmIndices_.begin(), params_.mmIndices_.end());
        std::tie(mmPositions, mmCharges) =
                getMMPositionsAndCharges(params_.mmAtoms_, fInput.x_, fInput.chargeA_);
    }
    // check that pairlist is available if needed
    if (params_.modelNeedsInput("atom-pairs") || params_.modelNeedsInput("pair-shifts"))
    {
        preparePairlistInput(box);
    }

    // get link atom info
    std::vector<LinkFrontierAtom> linkFrontier(constructLinkFrontier(params_.linkFrontier_, positions));

    if (params_.modelNeedsInput("atom-positions-mm"))
    {
        // center MM atom positions in the box
        centerAtomPositions(positions, mmPositions, box, *(params_.pbcType_));
    }

    model_->evaluateModel(&(fOutput->enerd_),
                          fOutput->forceWithVirial_.force_,
                          inputToLocalIndex_,
                          mmIndices,
                          params_.modelInput_,
                          positions,
                          atomNumbers_,
                          pairlistForModel_,
                          shiftVectors_,
                          mmPositions,
                          mmCharges,
                          params_.nnpCharge_,
                          linkFrontier,
                          box,
                          *(params_.pbcType_));
}

void NNPotForceProvider::gatherAtomNumbersIndices(const MDModulesAtomsRedistributedSignal& signal)
{
    // create lookup table for local atom indices needed for hybrid ML/MM
    // -1 is used as a flag for atoms that are not local / part of the NN input
    // used to distribute forces to correct local indices as the NN input tensor does not contain all atoms
    const int numInput = params_.nnpAtoms_->numAtomsGlobal();

    inputToLocalIndex_.assign(numInput, -1);
    inputToGlobalIndex_.assign(numInput, -1);
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
                        inputToLocalIndex_[j]  = i;
                        inputToGlobalIndex_[j] = globalIdx;
                        atomNumbers_[j]        = params_.atoms_.atom[globalIdx].atomnumber;
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
            inputToGlobalIndex_[i]         = globalIdx;
            localToInputIndex_[localIndex] = i;
            atomNumbers_[i]                = params_.atoms_.atom[globalIdx].atomnumber;
        }
    }
    // sanity check: make sure all atom numbers have been set
    GMX_RELEASE_ASSERT(std::count(atomNumbers_.begin(), atomNumbers_.end(), 0) == 0,
                       "Some atom numbers have not been set correctly.");
}

std::vector<RVec> NNPotForceProvider::gatherAtomPositions(ArrayRef<const RVec> positions) const
{
    // collect atom positions
    // at this point, we already have the atom numbers and indices, so we can fill the positions
    size_t numInput = inputToLocalIndex_.size();

    // initialize NN atom positions vector and fill according to lookup table
    std::vector<RVec> nnPositions(numInput);
    for (size_t i = 0; i < numInput; i++)
    {
        // if value in lookup table is -1, the atom is not local to this rank
        if (inputToLocalIndex_[i] != -1)
        {
            nnPositions[i] = positions[inputToLocalIndex_[i]];
        }
    }
    // in case of dom dec, distribute positions to all ranks
    if (mpiComm_.isParallel())
    {
        // TODO: this should use DIM or its successor, but using DIM here causes conflicts with torch
        mpiComm_.sumReduce(3 * numInput, nnPositions.data()->as_vec());
    }

    return nnPositions;
}

std::vector<LinkFrontierAtom>
NNPotForceProvider::constructLinkFrontier(const std::vector<LinkFrontierAtom>& inputLinkFrontier,
                                          ArrayRef<RVec>                       positions)
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
        const RVec posMM  = positions[inputIdxMM.value()];
        const RVec posNNP = positions[inputIdxNNP.value()];
        link.setPositions(posNNP, posMM);
        // update position of MM atom to link atom position
        positions[inputIdxMM.value()] = link.getLinkPosition();
        // overwrite atomic number of MM atom
        atomNumbers_[inputIdxMM.value()] = ilink.linkAtomNumber();
        linkFrontier.push_back(link);
    }
    return linkFrontier;
}

void NNPotForceProvider::setPairlist(const MDModulesPairlistConstructedSignal& signal)
{
    // We don't have the updated box vectors yet here, so we can't compute the shift vectors.
    // So, we just store the pairlist and prepare the input later in calculateForces().
    fullPairlist_.assign(signal.excludedPairlist_.begin(), signal.excludedPairlist_.end());
    doPairlist_ = true;
}

void NNPotForceProvider::preparePairlistInput(const matrix& box)
{
    // New pair list constructed: Indices in the pairlist correspond to local atom indices.
    // For now, find all pairs of NNP atoms within the cutoff (which were excluded from the
    // short-range calculation in GROMACS) so the NN potential does not have to re-compute them.
    // Thus, we're interested in the excluded pairlist, because all NNP-atom pairs are excluded
    // pairs, but not all excluded pairs are necessarily NNP-atom pairs.
    // TODO: figure out dom.dec. case.
    if (!doPairlist_)
    {
        return;
    }
    // this assert should only trigger if something breaks in the mdmodules notifications
    GMX_ASSERT(!fullPairlist_.empty(), "Pairlist for NNP model is empty!");

    const int numPairs = gmx::ssize(fullPairlist_);
    pairlistForModel_.clear();
    pairlistForModel_.reserve(2 * numPairs);
    shiftVectors_.clear();
    shiftVectors_.reserve(numPairs);

    for (int i = 0; i < numPairs; i++)
    {
        const auto [atomPair, shiftIndex] = fullPairlist_[i];

        // we only want to add pairs where both atoms are part of the NNP input
        // if any of the two (or more commonly both) is not, we skip this pair
        if (const std::optional<ptrdiff_t> inputIdxA = indexOf(inputToGlobalIndex_, atomPair.first);
            inputIdxA.has_value())
        {
            if (const std::optional<ptrdiff_t> inputIdxB = indexOf(inputToGlobalIndex_, atomPair.second);
                inputIdxB.has_value())
            {
                // no need to check cutoff: pairlist already comes filtered by cutoff
                RVec       shift;
                const IVec unitShift = shiftIndexToXYZ(shiftIndex);
                mvmul_ur0(box, unitShift.toRVec(), shift);

                pairlistForModel_.push_back(inputIdxA.value());
                pairlistForModel_.push_back(inputIdxB.value());
                shiftVectors_.push_back(shift);
            }
        }
    }
    // sanity check
    GMX_RELEASE_ASSERT(pairlistForModel_.size() == shiftVectors_.size() * 2,
                       "Inconsistent pairlist and shift vector sizes");
    doPairlist_ = false;
}

} // namespace gmx
