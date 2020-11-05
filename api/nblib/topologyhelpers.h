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
 * Implements helper functions needed for the nblib topology
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_TOPOLOGY_HELPERS_H
#define NBLIB_TOPOLOGY_HELPERS_H

#include <tuple>
#include <unordered_map>
#include <vector>

#include "gromacs/utility/listoflists.h"
#include "nblib/listed_forces/traits.h"
#include "nblib/molecules.h"

namespace gmx
{
struct ExclusionBlock;
}

namespace nblib
{

namespace detail
{

//! Converts tuples of particle indices to exclude to the gmx::ExclusionBlock format
std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList);

//! Add offset to all indices in inBlock
std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset);

/*!
 * \brief
 * Extract all interactions of type I from a vector of molecules. The second argument tuple element
 * specifies multiples of the molecule given as first tuple element. Let (S, I) denote the return
 * value tuple. Then J[i] = I[S[i]] for all i in 0...S.size() is the full sequence of BondType
 * instances as they occur in the input tuple
 *
 */
template<class I>
std::tuple<std::vector<size_t>, std::vector<I>>
collectInteractions(const std::vector<std::tuple<Molecule, int>>&);

#define COLLECT_BONDS_EXTERN_TEMPLATE(x)                                                 \
    extern template std::tuple<std::vector<size_t>, std::vector<x>> collectInteractions( \
            const std::vector<std::tuple<Molecule, int>>&);
MAP(COLLECT_BONDS_EXTERN_TEMPLATE, SUPPORTED_TWO_CENTER_TYPES)
#undef COLLECT_BONDS_EXTERN_TEMPLATE

/*!
 * \brief
 * Return a list of unique BondType instances U and an index list S of size aggregatedBonds.size()
 * such that the BondType instance at aggregatedBonds[i] is equal to U[S[i]]
 * returns std::tuple(S, U)
 *
 */
template<class I>
std::tuple<std::vector<size_t>, std::vector<I>> eliminateDuplicateInteractions(const std::vector<I>& collectedBonds);

/// \cond DO_NOT_DOCUMENT
#define ELIMINATE_DUPLICATE_EXTERN_TEMPLATE(x)                                                      \
    extern template std::tuple<std::vector<size_t>, std::vector<x>> eliminateDuplicateInteractions( \
            const std::vector<x>& collectedBonds);
MAP(ELIMINATE_DUPLICATE_EXTERN_TEMPLATE, SUPPORTED_LISTED_TYPES)
#undef ELIMINATE_DUPLICATE_EXTERN_TEMPLATE
/// \endcond

//! Helper class for Topology to keep track of particle IDs
class ParticleSequencer
{
    //! Alias for storing by (molecule name, molecule nr, residue name, particle name)
    using DataType = std::unordered_map<
            std::string,
            std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, int>>>>;

public:
    //! Build sequence from a list of molecules
    void build(const std::vector<std::tuple<Molecule, int>>& moleculesList);

    //! Access ID by (molecule name, molecule nr, residue name, particle name)
    int operator()(const MoleculeName&, int, const ResidueName&, const ParticleName&) const;

private:
    DataType data_;
};

//!
template<class B>
std::vector<CoordinateIndex<B>> sequenceIDs(const std::vector<std::tuple<Molecule, int>>&,
                                            const detail::ParticleSequencer&);

/// \cond DO_NOT_DOCUMENT
#define SEQUENCE_PAIR_ID_EXTERN_TEMPLATE(x)                         \
    extern template std::vector<CoordinateIndex<x>> sequenceIDs<x>( \
            const std::vector<std::tuple<Molecule, int>>&, const detail::ParticleSequencer&);
MAP(SEQUENCE_PAIR_ID_EXTERN_TEMPLATE, SUPPORTED_LISTED_TYPES)
#undef SEQUENCE_PAIR_ID_EXTERN_TEMPLATE
/// \endcond

} // namespace detail

} // namespace nblib

#endif // NBLIB_TOPOLOGY_HELPERS_H
