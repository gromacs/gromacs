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
/*! \internal \file
 * \brief
 * Implements nblib Topology helpers
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include <numeric>

#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/utility/smalloc.h"
#include "nblib/exception.h"
#include "nblib/listed_forces/transformations.h"
#include "nblib/topologyhelpers.h"
#include "nblib/util/internal.h"

namespace nblib
{

namespace detail
{

std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList)
{
    std::vector<gmx::ExclusionBlock> ret;

    auto firstLowerThan = [](auto const& tup1, auto const& tup2) {
        return std::get<0>(tup1) < std::get<0>(tup2);
    };

    // Note this is a programming error as all particles should exclude at least themselves and empty topologies are not allowed.
    // Note also that this is also checked in the parent function with a more informative error message.
    assert((!tupleList.empty() && "No exclusions found.\n"));

    // initialize pair of iterators delimiting the range of exclusions for
    // the first particle in the list
    auto range = std::equal_range(std::begin(tupleList), std::end(tupleList), tupleList[0], firstLowerThan);
    auto it1 = range.first;
    auto it2 = range.second;

    // loop over all exclusions in molecule, linear in tupleList.size()
    while (it1 != std::end(tupleList))
    {
        gmx::ExclusionBlock localBlock;
        // loop over all exclusions for current particle
        for (; it1 != it2; ++it1)
        {
            localBlock.atomNumber.push_back(std::get<1>(*it1));
        }

        ret.push_back(localBlock);

        // update the upper bound of the range for the next particle
        if (it1 != end(tupleList))
        {
            it2 = std::upper_bound(it1, std::end(tupleList), *it1, firstLowerThan);
        }
    }

    return ret;
}

std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset)
{
    // shift particle numbers by offset
    for (auto& localBlock : inBlock)
    {
        std::transform(std::begin(localBlock.atomNumber),
                       std::end(localBlock.atomNumber),
                       std::begin(localBlock.atomNumber),
                       [offset](auto i) { return i + offset; });
    }

    return inBlock;
}

int ParticleSequencer::operator()(const MoleculeName& moleculeName,
                                  int                 moleculeNr,
                                  const ResidueName&  residueName,
                                  const ParticleName& particleName) const
{
    try
    {
        return data_.at(moleculeName).at(moleculeNr).at(residueName).at(particleName);
    }
    catch (const std::out_of_range& outOfRange)
    {
        // TODO: use string format function once we have it
        if (moleculeName.value() == residueName.value())
        {
            printf("No particle %s in residue %s in molecule %s found\n",
                   particleName.value().c_str(),
                   residueName.value().c_str(),
                   moleculeName.value().c_str());
        }
        else
        {
            printf("No particle %s in molecule %s found\n",
                   particleName.value().c_str(),
                   moleculeName.value().c_str());
        }

        throw InputException(outOfRange.what());
    };
}

void ParticleSequencer::build(const std::vector<std::tuple<Molecule, int>>& moleculesList)
{
    int currentID = 0;
    for (auto& molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        const size_t    numMols  = std::get<1>(molNumberTuple);

        auto& moleculeMap = data_[molecule.name()];

        for (size_t i = 0; i < numMols; ++i)
        {
            auto& moleculeNrMap = moleculeMap[i];
            for (int j = 0; j < molecule.numParticlesInMolecule(); ++j)
            {
                moleculeNrMap[molecule.residueName(j)][molecule.particleName(j)] = currentID++;
            }
        }
    }
}

template<class I>
std::tuple<std::vector<size_t>, std::vector<I>>
collectInteractions(const std::vector<std::tuple<Molecule, int>>& molecules)
{
    std::vector<I>      collectedBonds;
    std::vector<size_t> expansionArray;
    for (auto& molNumberTuple : molecules)
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

/// \cond DO_NOT_DOCUMENT
#define COLLECT_BONDS_INSTANTIATE_TEMPLATE(x)                                     \
    template std::tuple<std::vector<size_t>, std::vector<x>> collectInteractions( \
            const std::vector<std::tuple<Molecule, int>>&);
MAP(COLLECT_BONDS_INSTANTIATE_TEMPLATE, SUPPORTED_LISTED_TYPES)
/// \endcond

namespace sequence_detail
{

//! Helper function to convert a tuple of strings into a particle index sequence
template<class Tuple, class F, class... Args, size_t... Is>
auto stringsToIndices_impl(const Tuple& tuple, std::index_sequence<Is...> is, F&& f, Args... args)
{
    // return std::make_tuple(f(args..., std::get<2 * Is + 1>(tuple), std::get<2 * Is>(tuple))...);
    ignore_unused(is);
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

//! Tuple ordering for two center interactions
[[maybe_unused]] static std::array<int, 2> nblibOrdering(const std::array<int, 2>& t)
{
    // for bonds (two coordinate indices),
    // we choose to store the lower sequence ID first. this allows for better unit tests
    // that are agnostic to how the input was set up
    int id1 = std::min(std::get<0>(t), std::get<1>(t));
    int id2 = std::max(std::get<0>(t), std::get<1>(t));

    return std::array<int, 2>{ id1, id2 };
}

//! Tuple ordering for three center interactions
[[maybe_unused]] static std::array<int, 3> nblibOrdering(const std::array<int, 3>& t)
{
    // for angles (three coordinate indices),
    // we choose to store the two non-center coordinate indices sorted.
    // such that ret[0] < ret[2] always (ret = returned tuple)
    int id1 = std::min(std::get<0>(t), std::get<2>(t));
    int id3 = std::max(std::get<0>(t), std::get<2>(t));

    return std::array<int, 3>{ id1, std::get<1>(t), id3 };
}

//! Tuple ordering for four center interactions
[[maybe_unused]] static std::array<int, 4> nblibOrdering(const std::array<int, 4>& t)
{
    return t;
}

//! Tuple ordering for five center interactions
[[maybe_unused]] static std::array<int, 5> nblibOrdering(const std::array<int, 5>& t)
{
    return t;
}

} // namespace sequence_detail

//! For each interaction, translate particle identifiers (moleculeName, nr, residueName,
//! particleName) to particle coordinate indices
template<class B>
std::vector<CoordinateIndex<B>> sequenceIDs(const std::vector<std::tuple<Molecule, int>>& molecules,
                                            const detail::ParticleSequencer& particleSequencer)
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

/// \cond DO_NOT_DOCUMENT
#define SEQUENCE_PAIR_ID_INSTANTIATE_TEMPLATE(x)             \
    template std::vector<CoordinateIndex<x>> sequenceIDs<x>( \
            const std::vector<std::tuple<Molecule, int>>&, const detail::ParticleSequencer&);
MAP(SEQUENCE_PAIR_ID_INSTANTIATE_TEMPLATE, SUPPORTED_LISTED_TYPES)
#undef SEQUENCE_PAIR_ID_INSTANTIATE_TEMPLATE
/// \endcond

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

/// \cond DO_NOT_DOCUMENT
#define ELIMINATE_DUPLICATE_BONDS_INSTANTIATE_TEMPLATE(x)                                    \
    template std::tuple<std::vector<size_t>, std::vector<x>> eliminateDuplicateInteractions( \
            const std::vector<x>& aggregatedBonds);
MAP(ELIMINATE_DUPLICATE_BONDS_INSTANTIATE_TEMPLATE, SUPPORTED_LISTED_TYPES)
#undef ELIMINATE_DUPLICATE_BONDS_INSTANTIATE_TEMPLATE
/// \endcond

} // namespace detail

} // namespace nblib
