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
 * Implements helper functions needed for the nblib topology that
 * are part of the translation of string-based particle identifiers
 * used in Molecule to the sequence of integer IDs that the topology
 * stores. This file is included only in the topology.cpp translation unit.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_SEQUENCING_H
#define NBLIB_SEQUENCING_H

#include <numeric>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "listed_forces/transformations.h"

#include "gromacs/utility/listoflists.h"

#include "nblib/molecules.h"
#include "nblib/particlesequencer.h"

namespace nblib
{

namespace detail
{

/*! \brief Extract all interactions of type I from a vector of molecules.
 *
 * The second argument tuple element specifies multiples of the molecule given as first tuple
 * element. Let (S, I) denote the return value tuple. Then J[i] = I[S[i]] for all i in 0...S.size()
 * is the full sequence of BondType instances as they occur in the input tuple
 */
template<class I>
std::tuple<std::vector<size_t>, std::vector<I>>
collectInteractions(const std::vector<std::tuple<Molecule, int>>& molecules)
{
    std::vector<I>      collectedBonds;
    std::vector<size_t> expansionArray;
    for (const auto& molNumberTuple : molecules)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        auto& interactions = pickType<I>(molecule.interactionData()).interactionTypes_;

        std::vector<size_t> moleculeExpansion(interactions.size());
        // assign indices to the bonds in the current molecule, continue counting from
        // the number of bonds seen so far (=collectedBonds.size())
        std::iota(begin(moleculeExpansion), end(moleculeExpansion), collectedBonds.size());

        std::copy(begin(interactions), end(interactions), std::back_inserter(collectedBonds));

        for (size_t i = 0; i < numMols; ++i)
        {
            std::copy(begin(moleculeExpansion), end(moleculeExpansion), std::back_inserter(expansionArray));
        }
    }
    return std::make_tuple(expansionArray, collectedBonds);
}

/*! \brief Return a list of unique BondType instances
 *
 * U is the unique list of bonds and S is an index list of size aggregatedBonds.size(),
 * such that the BondType instance at aggregatedBonds[i] is equal to U[S[i]]
 * returns std::tuple(S, U)
 */
template<class I>
std::tuple<std::vector<size_t>, std::vector<I>> eliminateDuplicateInteractions(const std::vector<I>& aggregatedInteractions)
{
    std::vector<size_t> uniqueIndices(aggregatedInteractions.size());
    std::vector<I>      uniquInteractionsInstances;
    // if there are no interactions of type B we're done now
    if (aggregatedInteractions.empty())
    {
        return std::make_tuple(uniqueIndices, uniquInteractionsInstances);
    }

    // create 0,1,2,... sequence
    std::iota(begin(uniqueIndices), end(uniqueIndices), 0);

    std::vector<std::tuple<I, size_t>> enumeratedBonds(aggregatedInteractions.size());
    // append each interaction with its index
    std::transform(begin(aggregatedInteractions),
                   end(aggregatedInteractions),
                   begin(uniqueIndices),
                   begin(enumeratedBonds),
                   [](I b, size_t i) { return std::make_tuple(b, i); });

    auto sortKey = [](const auto& t1, const auto& t2) { return std::get<0>(t1) < std::get<0>(t2); };
    // sort w.r.t bonds. the result will contain contiguous segments of identical bond instances
    // the associated int indicates the original index of each BondType instance in the input vector
    std::sort(begin(enumeratedBonds), end(enumeratedBonds), sortKey);

    // initialize it1 and it2 to delimit first range of equal BondType instances
    auto range = std::equal_range(begin(enumeratedBonds), end(enumeratedBonds), enumeratedBonds[0], sortKey);
    auto it1 = range.first;
    auto it2 = range.second;

    // number of unique instances of BondType B = number of contiguous segments in enumeratedBonds =
    //         number of iterations in the outer while loop below
    while (it1 != end(enumeratedBonds))
    {
        uniquInteractionsInstances.push_back(std::get<0>(*it1));

        // loop over all identical BondType instances;
        for (; it1 != it2; ++it1)
        {
            // we note down that the BondType instance at index <interactionIndex>
            // can be found in the uniqueBondInstances container at index <uniqueBondInstances.size()>
            int interactionIndex            = std::get<1>(*it1);
            uniqueIndices[interactionIndex] = uniquInteractionsInstances.size() - 1;
        }

        // Note it1 has been incremented and is now equal to it2
        if (it1 != end(enumeratedBonds))
        {
            it2 = std::upper_bound(it1, end(enumeratedBonds), *it1, sortKey);
        }
    }

    return make_tuple(uniqueIndices, uniquInteractionsInstances);
}


namespace sequence_detail
{

//! \brief Helper function to convert a tuple of strings into a particle index sequence
template<class Tuple, class F, class... Args, size_t... Is>
auto stringsToIndices_impl(const Tuple& tuple, [[maybe_unused]] std::index_sequence<Is...> is, F&& f, Args... args)
{
    return std::array<int, sizeof...(Is)>{ f(
            args..., std::get<2 * Is + 1>(tuple), std::get<2 * Is>(tuple))... };
}

/*! \brief
 *  This takes a tuple<(string, string) * nCenter> from molecule
 *  where nCenter = 2 for bonds, 3 for angles and 4 for dihedrals
 *  each (ResidueName, ParticleName)-pair is converted to a particle sequence index
 *  by calling the supplied function object f, containing the particleSequencer at the call site
 *  Therefore, the return type is tuple<int * nCenter>
 *
 */
template<class Tuple, class F, class... Args>
auto stringsToIndices(const Tuple& tuple, F&& f, Args... args)
{
    auto is = std::make_index_sequence<std::tuple_size<Tuple>::value / 2>{};
    return stringsToIndices_impl(tuple, is, std::forward<F>(f), args...);
}

} // namespace sequence_detail

/*! \brief Convert string-based interactions to integer-index based ones
 *
 * For each interaction, translate particle identifiers (moleculeName, nr, residueName,
 * particleName) to particle coordinate indices
 */
template<class B>
std::vector<CoordinateIndex<B>> sequenceIDs(const std::vector<std::tuple<Molecule, int>>& molecules,
                                            const ParticleSequencer& particleSequencer)
{
    std::vector<CoordinateIndex<B>> coordinateIndices;

    auto callSequencer = [&particleSequencer](const MoleculeName& moleculeName,
                                              int                 i,
                                              const ResidueName&  residueName,
                                              const ParticleName& particleName) {
        return particleSequencer(moleculeName, i, residueName, particleName);
    };

    // loop over all molecules
    for (const auto& molNumberTuple : molecules)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        for (size_t i = 0; i < numMols; ++i)
        {
            auto& interactions = pickType<B>(molecule.interactionData()).interactions_;
            for (const auto& interactionString : interactions)
            {
                CoordinateIndex<B> index = sequence_detail::stringsToIndices(
                        interactionString, callSequencer, molecule.name(), i);
                coordinateIndices.push_back(nblibOrdering(index));
            }
        }
    }
    return coordinateIndices;
}

} // namespace detail

} // namespace nblib

#endif // NBLIB_TOPOLOGY_HELPERS_H
