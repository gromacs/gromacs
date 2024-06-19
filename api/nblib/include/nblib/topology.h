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
/*! \inpublicapi \file
 * \brief
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_TOPOLOGY_H
#define NBLIB_TOPOLOGY_H

#include <cstddef>

#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "nblib/basicdefinitions.h"
#include "nblib/interactions.h"
#include "nblib/kerneloptions.h"
#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/molecules.h"
#include "nblib/particlesequencer.h"
#include "nblib/particletype.h"

namespace nblib
{

/*! \brief A list of exclusion lists, optimized for performance
 *
 * This struct is to hold data that comes from, and will be copied into,
 * gmx::ListOfLists. Since the ListElements and ListRanges are sufficient
 * to create a gmx::ListOfLists, they are stored in the topoology instead
 * of a gmx::ListOfLists to avoid dependencies.
 */
template<typename T>
struct ExclusionLists
{
    //! Ranges of the concatenated lists
    std::vector<int> ListRanges;
    //! Elements of the concatenated lists
    std::vector<T> ListElements;
};

/*! \inpublicapi
 * \ingroup nblib
 * \brief System Topology
 *
 * Contains all topology information meant to be used by the simulation
 * engine internally. Private constructor ensures that a Topology object
 * exists in a scope in a valid state after it has been built using a
 * Topology Builder.
 */
class Topology final
{

public:
    //! Returns the total number of particles in the system
    int numParticles() const;

    //! Returns a vector of particle types
    std::vector<ParticleType> getParticleTypes() const;

    //! Return the ParticleType ID of all particles
    std::vector<int> getParticleTypeIdOfAllParticles() const;

    //! Returns a vector of particles partial charges
    std::vector<real> getCharges() const;

    //! Returns exclusions in a format that can be used to create gmx::ListOfLists
    ExclusionLists<int> exclusionLists() const;

    //! Returns the unique ID of a specific particle belonging to a molecule in the global space
    int sequenceID(MoleculeName moleculeName, int moleculeNr, ResidueName residueName, ParticleName particleName) const;

    //! Returns a map of non-bonded force parameters indexed by ParticleType names
    NonBondedInteractionMap getNonBondedInteractionMap() const;

    //! Returns the interaction data
    ListedInteractionData getInteractionData() const;

    //! Returns the combination rule used to generate the NonBondedInteractionMap
    CombinationRule getCombinationRule() const;

private:
    Topology() = default;

    friend class TopologyBuilder;

    //! Total number of particles in the system
    int numParticles_;
    //! unique collection of ParticleTypes
    std::vector<ParticleType> particleTypes_;
    //! store an ID of each particles's type
    std::vector<int> particleTypeIdOfAllParticles_;
    //! Storage for particles partial charges
    std::vector<real> charges_;
    //! Information about exclusions.
    ExclusionLists<int> exclusionLists_;
    //! Associate molecule, residue and particle names with sequence numbers
    ParticleSequencer particleSequencer_;
    //! Map that should hold all nonbonded interactions for all particle types
    NonBondedInteractionMap nonBondedInteractionMap_;
    //! data about bonds for all supported types
    ListedInteractionData interactionData_;
    //! Combination Rule used to generate the nonbonded interactions
    CombinationRule combinationRule_;
};

/*! \brief Topology Builder
 *
 * \libinternal
 * \ingroup nblib
 *
 * A helper class to assist building of topologies. They also ensure that
 * topologies only exist in a valid state within the scope of the
 * simulation program.
 *
 */
class TopologyBuilder final
{
public:
    //! Constructor
    TopologyBuilder();

    /*! \brief
     * Builds and Returns a valid Topology
     *
     * This function accounts for all the molecules added along with their
     * exclusions and returns a topology with a valid state that is usable
     * by the GROMACS back-end.
     */
    Topology buildTopology();

    //! Adds a molecules of a certain type into the topology
    TopologyBuilder& addMolecule(const Molecule& moleculeType, int nMolecules);

    //! Add non-bonded interaction map to the topology
    void addParticleTypesInteractions(const ParticleTypesInteractions& particleTypesInteractions);

private:
    //! Internally stored topology
    Topology topology_;

    //! Total number of particles in the system
    int numParticles_;

    //! List of molecule types and number of molecules
    std::vector<std::tuple<Molecule, int>> molecules_;

    //! Builds an exclusions list aggregating exclusions from all molecules
    ExclusionLists<int> createExclusionsLists() const;

    //! Gather interaction data from molecules
    ListedInteractionData createInteractionData(const ParticleSequencer&);

    //! Helper function to extract quantities like mass, charge, etc from the system
    template<typename T, class Extractor>
    std::vector<T> extractParticleTypeQuantity(Extractor&& extractor);

    //! Distinct collection of ParticleTypes
    std::unordered_map<std::string, ParticleType> particleTypes_;

    //! ParticleType nonbonded parameters
    ParticleTypesInteractions particleTypesInteractions_;
};

//! utility function to extract Particle quantities and expand them to the full
//! array of length numParticles()
template<class F>
inline auto expandQuantity(const Topology& topology, F&& particleTypeExtractor)
{
    using ValueType = decltype((std::declval<ParticleType>().*std::declval<F>())());

    std::vector<ValueType> ret;
    ret.reserve(topology.numParticles());

    const std::vector<ParticleType>& particleTypes = topology.getParticleTypes();

    for (size_t id : topology.getParticleTypeIdOfAllParticles())
    {
        ret.push_back((particleTypes[id].*particleTypeExtractor)());
    }

    return ret;
}

} // namespace nblib

#endif // NBLIB_TOPOLOGY_H
