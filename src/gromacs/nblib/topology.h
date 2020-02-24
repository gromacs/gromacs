/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \inpublicapi \file
 * \brief
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/nblib/molecules.h"
#include "gromacs/utility/listoflists.h"

namespace gmx
{
    struct ExclusionBlock;
}

namespace nblib
{

    namespace detail
    {

// Converts tuples of atom indices to exclude to the gmx::ExclusionBlock format
        std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList);

// Add offset to all indices in inBlock
        std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset);

    } // namespace detail

/*! \inpublicapi
 * \ingroup nblib
 * \brief System Topology
 *
 * Contains all topology information meant to be used by the simulation
 * engine internally. Private constructor ensures that a Topology object
 * exists in a scope in a valid state after it has been built using a
 * Topology Builder.
 */
    class Topology
    {
    public:
        //! Returns the total number of atoms in the system
        const int& numAtoms() const;

        //! Returns a vector of atom names
        const std::vector<AtomType>& getAtomTypes() const;

        //! Return the AtomType ID of all atoms
        const std::vector<int>& getAtomTypeIdOfAllAtoms() const;

        //! Returns a vector of atom partial charges
        const std::vector<real>& getCharges() const;

        //! Returns exclusions in proper, performant, gromacs layout
        const gmx::ListOfLists<int>& getGmxExclusions() const { return exclusions_; }

    private:
        Topology() = default;

        friend class TopologyBuilder;

        //! Total number of atoms in the system
        int numAtoms_;
        //! unique collection of AtomTypes
        std::vector<AtomType> atomTypes_;
        //! store an ID of each atom's type
        std::vector<int> atomTypeIdOfAllAtoms_;
        //! Storage for atom partial charges
        std::vector<real> charges_;
        //! Information about exclusions.
        gmx::ListOfLists<int> exclusions_;
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
    class TopologyBuilder
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

        // Adds a molecules of a certain type into the topology
        TopologyBuilder& addMolecule(const Molecule& moleculeType, int nMolecules);

    private:
        //! Internally stored topology
        Topology topology_;

        //! Total number of atoms in the system
        int numAtoms_;

        //! List of molecule types and number of molecules
        std::vector<std::tuple<Molecule, int>> molecules_;

        // Builds a GROMACS-compliant performant exclusions list aggregating exclusions from all molecules
        gmx::ListOfLists<int> createExclusionsListOfLists() const;

        // Helper function to extract quantities like mass, charge, etc from the system
        template<typename T, class Extractor>
        std::vector<T> extractAtomTypeQuantity(Extractor extractor);

        //! distinct collection of AtomTypes
        std::unordered_map<std::string, AtomType> atomTypes_;
    };

//! utility function to extract AtomType quantities and expand them to the full
//! array of length numAtoms()
    template<class F>
    inline auto expandQuantity(const Topology& topology, F atomTypeExtractor)
    {
        using ValueType = decltype((std::declval<AtomType>().*std::declval<F>())());

        std::vector<ValueType> ret;
        ret.reserve(topology.numAtoms());

        const std::vector<AtomType>& atomTypes = topology.getAtomTypes();

        for (size_t id : topology.getAtomTypeIdOfAllAtoms())
        {
            ret.push_back((atomTypes[id].*atomTypeExtractor)());
        }

        return ret;
    }

} // namespace nblib

#endif // GROMACS_TOPOLOGY_H