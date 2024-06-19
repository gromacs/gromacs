/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/topology.h"

#include <cassert>

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
#include <type_traits>

#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/smalloc.h"

#include "nblib/exception.h"
#include "nblib/particletype.h"
#include "nblib/util/util.hpp"

#include "sequencing.hpp"
#include "topologyhelpers.h"

namespace nblib
{

TopologyBuilder::TopologyBuilder() : numParticles_(0) {}

ExclusionLists<int> TopologyBuilder::createExclusionsLists() const
{
    const auto& moleculesList = molecules_;

    std::vector<gmx::ExclusionBlock> exclusionBlockGlobal;
    exclusionBlockGlobal.reserve(numParticles_);

    size_t particleNumberOffset = 0;
    for (const auto& molNumberTuple : moleculesList)
    {
        const Molecule& molecule   = std::get<0>(molNumberTuple);
        size_t          numMols    = std::get<1>(molNumberTuple);
        const auto&     exclusions = molecule.getExclusions();

        // Note this is a programming error as all particles should exclude at least themselves and empty topologies are not allowed.
        const std::string message =
                "No exclusions found in the " + molecule.name().value() + " molecule.";
        assert((!exclusions.empty() && message.c_str()));

        std::vector<gmx::ExclusionBlock> exclusionBlockPerMolecule = toGmxExclusionBlock(exclusions);

        // duplicate the exclusionBlockPerMolecule for the number of Molecules of (numMols)
        for (size_t i = 0; i < numMols; ++i)
        {
            auto offsetExclusions = offsetGmxBlock(exclusionBlockPerMolecule, particleNumberOffset);

            std::copy(std::begin(offsetExclusions),
                      std::end(offsetExclusions),
                      std::back_inserter(exclusionBlockGlobal));

            particleNumberOffset += molecule.numParticlesInMolecule();
        }
    }

    gmx::ListOfLists<int> exclusionsListOfListsGlobal;
    for (const auto& block : exclusionBlockGlobal)
    {
        exclusionsListOfListsGlobal.pushBack(block.atomNumber);
    }

    std::vector<int>    listRanges(exclusionsListOfListsGlobal.listRangesView().begin(),
                                exclusionsListOfListsGlobal.listRangesView().end());
    std::vector<int>    listElements(exclusionsListOfListsGlobal.elementsView().begin(),
                                  exclusionsListOfListsGlobal.elementsView().end());
    ExclusionLists<int> exclusionListsGlobal;
    exclusionListsGlobal.ListRanges   = std::move(listRanges);
    exclusionListsGlobal.ListElements = std::move(listElements);

    return exclusionListsGlobal;
}

ListedInteractionData TopologyBuilder::createInteractionData(const ParticleSequencer& particleSequencer)
{
    ListedInteractionData interactionData;

    // this code is doing the compile time equivalent of
    // for (int i = 0; i < interactionData.size(); ++i)
    //     create(get<i>(interactionData));

    auto create = [this, &particleSequencer](auto& interactionDataElement) {
        using InteractionType = typename std::decay_t<decltype(interactionDataElement)>::type;

        // first compression stage: each bond per molecule listed once,
        // eliminates duplicates from multiple identical molecules
        auto  compressedDataStage1 = detail::collectInteractions<InteractionType>(this->molecules_);
        auto& expansionArrayStage1 = std::get<0>(compressedDataStage1);
        auto& moleculeInteractions = std::get<1>(compressedDataStage1);

        // second compression stage: recognize bond duplicates among bonds from all molecules put together
        auto  compressedDataStage2 = detail::eliminateDuplicateInteractions(moleculeInteractions);
        auto& expansionArrayStage2 = std::get<0>(compressedDataStage2);
        auto& uniqueInteractionInstances = std::get<1>(compressedDataStage2);

        // combine stage 1 + 2 expansion arrays
        std::vector<size_t> expansionArray(expansionArrayStage1.size());
        std::transform(begin(expansionArrayStage1),
                       end(expansionArrayStage1),
                       begin(expansionArray),
                       [&S2 = expansionArrayStage2](size_t S1Element) { return S2[S1Element]; });

        // add data about InteractionType instances
        interactionDataElement.parameters = std::move(uniqueInteractionInstances);

        interactionDataElement.indices.resize(expansionArray.size());
        // coordinateIndices contains the particle sequence IDs of all interaction coordinates of type <BondType>
        auto coordinateIndices = detail::sequenceIDs<InteractionType>(this->molecules_, particleSequencer);
        // zip coordinateIndices(i,j,...) + expansionArray(k) -> interactionDataElement.indices(i,j,...,k)
        std::transform(begin(coordinateIndices),
                       end(coordinateIndices),
                       begin(expansionArray),
                       begin(interactionDataElement.indices),
                       [](auto coordinateIndex, auto interactionIndex) {
                           std::array<int, coordinateIndex.size() + 1> ret{ 0 };
                           for (int i = 0; i < int(coordinateIndex.size()); ++i)
                           {
                               ret[i] = coordinateIndex[i];
                           }
                           ret[coordinateIndex.size()] = interactionIndex;
                           return ret;
                       });
    };

    for_each_tuple(create, interactionData);

    return interactionData;
}

template<typename T, class Extractor>
std::vector<T> TopologyBuilder::extractParticleTypeQuantity(Extractor&& extractor)
{
    auto& moleculesList = molecules_;

    // returned object
    std::vector<T> ret;
    ret.reserve(numParticles_);

    for (auto& molNumberTuple : moleculesList)
    {
        Molecule& molecule = std::get<0>(molNumberTuple);
        size_t    numMols  = std::get<1>(molNumberTuple);

        for (size_t i = 0; i < numMols; ++i)
        {
            for (auto& particleData : molecule.particleData())
            {
                auto particleTypesMap = molecule.particleTypesMap();
                ret.push_back(extractor(particleData, particleTypesMap));
            }
        }
    }

    return ret;
}

Topology TopologyBuilder::buildTopology()
{
    assert((!(numParticles_ < 0) && "It should not be possible to have negative particles"));
    if (numParticles_ == 0)
    {
        throw InputException("You cannot build a topology with no particles");
    }
    topology_.numParticles_ = numParticles_;

    topology_.exclusionLists_ = createExclusionsLists();
    topology_.charges_        = extractParticleTypeQuantity<real>(
            [](const auto& data, [[maybe_unused]] auto& map) { return data.charge_; });

    // map unique ParticleTypes to IDs
    std::unordered_map<std::string, int> nameToId;
    for (auto& name_particleType_tuple : particleTypes_)
    {
        topology_.particleTypes_.push_back(name_particleType_tuple.second);
        nameToId[name_particleType_tuple.first] = nameToId.size();
    }

    topology_.particleTypeIdOfAllParticles_ =
            extractParticleTypeQuantity<int>([&nameToId](const auto& data, [[maybe_unused]] auto& map) {
                return nameToId[data.particleTypeName_];
            });

    ParticleSequencer particleSequencer;
    particleSequencer.build(molecules_);
    topology_.particleSequencer_ = std::move(particleSequencer);

    topology_.combinationRule_         = particleTypesInteractions_.getCombinationRule();
    topology_.nonBondedInteractionMap_ = particleTypesInteractions_.generateTable();

    topology_.interactionData_ = createInteractionData(topology_.particleSequencer_);

    // Check whether there is any missing term in the particleTypesInteractions compared to the
    // list of particletypes
    for (const auto& particleType1 : particleTypes_)
    {
        for (const auto& particleType2 : particleTypes_)
        {
            auto interactionKey = std::make_tuple(ParticleTypeName(particleType1.first),
                                                  ParticleTypeName(particleType2.first));
            if (topology_.nonBondedInteractionMap_.count(interactionKey) == 0)
            {
                std::string message =
                        formatString("Missing nonbonded interaction parameters for pair {} {}",
                                     particleType1.first,
                                     particleType2.first);
                throw InputException(message);
            }
        }
    }

    return topology_;
}

TopologyBuilder& TopologyBuilder::addMolecule(const Molecule& molecule, const int nMolecules)
{
    /*
     * 1. Push-back a tuple of molecule type and nMolecules
     * 2. Append exclusion list into the data structure
     */

    molecules_.emplace_back(molecule, nMolecules);
    numParticles_ += nMolecules * molecule.numParticlesInMolecule();

    auto particleTypesInMolecule = molecule.particleTypesMap();

    for (const auto& name_type_tuple : particleTypesInMolecule)
    {
        // If we already have the particleType, we need to make
        // sure that the type's parameters are actually the same
        // otherwise we would overwrite them
        if (particleTypes_.count(name_type_tuple.first) > 0)
        {
            if (!(particleTypes_.at(name_type_tuple.first) == name_type_tuple.second))
            {
                throw InputException("Differing ParticleTypes with identical names encountered");
            }
        }
    }

    // Note: insert does nothing if the key already exists
    particleTypes_.insert(particleTypesInMolecule.begin(), particleTypesInMolecule.end());

    return *this;
}

void TopologyBuilder::addParticleTypesInteractions(const ParticleTypesInteractions& particleTypesInteractions)
{
    particleTypesInteractions_.merge(particleTypesInteractions);
}

int Topology::numParticles() const
{
    return numParticles_;
}

std::vector<real> Topology::getCharges() const
{
    return charges_;
}

std::vector<ParticleType> Topology::getParticleTypes() const
{
    return particleTypes_;
}

std::vector<int> Topology::getParticleTypeIdOfAllParticles() const
{
    return particleTypeIdOfAllParticles_;
}

int Topology::sequenceID(MoleculeName moleculeName, int moleculeNr, ResidueName residueName, ParticleName particleName) const
{
    return particleSequencer_(moleculeName, moleculeNr, residueName, particleName);
}

NonBondedInteractionMap Topology::getNonBondedInteractionMap() const
{
    return nonBondedInteractionMap_;
}

ListedInteractionData Topology::getInteractionData() const
{
    return interactionData_;
}

CombinationRule Topology::getCombinationRule() const
{
    return combinationRule_;
}

ExclusionLists<int> Topology::exclusionLists() const
{
    return exclusionLists_;
}

} // namespace nblib
