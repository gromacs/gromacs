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
 * Transformations to combine elementary interaction types into aggregate types
 * are implemented here
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GMX_NBLIB_LISTED_FORCES_AGGREGATE_TRANSFORMATIONS_HPP
#define GMX_NBLIB_LISTED_FORCES_AGGREGATE_TRANSFORMATIONS_HPP

#include <unordered_map>

#include "nblib/listed_forces/traits.h"
#include "nblib/listed_forces/transformations.h"

namespace nblib
{

template<class InteractionData>
void deleteInteractions(InteractionData& interactions, const std::vector<int>& deleteList)
{
    std::vector<int> keep(interactions.indices.size(), 1);

    for (int deletion : deleteList)
    {
        keep[deletion] = 0;
    }

    auto trimmedInteractions = interactions.indices;
    trimmedInteractions.clear();

    for (size_t i = 0; i < keep.size(); ++i)
    {
        if (keep[i])
        {
            trimmedInteractions.push_back(interactions.indices[i]);
        }
    }

    swap(trimmedInteractions, interactions.indices);

    // Todo: possible optimization: also remove parameters not needed anymore
}

template<class CarrierType, class AggregateType>
void migrateParameters(ListedTypeData<CarrierType>& source, ListedTypeData<AggregateType>& destination)
{
    destination.indices.reserve(source.indices.size());
    destination.parametersA.reserve(source.indices.size());
    for (size_t i = 0; i < source.indices.size(); ++i)
    {
        auto sourceIndex = source.indices[i];
        auto sourceParam = source.parametersA[sourceIndex.back()];

        AggregateType destinationParam;
        destinationParam.carrier() = sourceParam;

        auto destinationIndex   = sourceIndex;
        destinationIndex.back() = i;

        destination.indices.push_back(destinationIndex);
        destination.parametersA.push_back(destinationParam);
    }

    // delete source
    source = ListedTypeData<CarrierType>();
}

template<class TwoCenterType, class AggregateType>
bool integrateTwoInThree(IndexArray<3>                        index,
                         int                                  searchIndex,
                         const ListedTypeData<TwoCenterType>& twoCSource,
                         ListedTypeData<AggregateType>&       aggregate)
{
    auto sortKeyObj = [](const auto& lhs, const auto& rhs) { return interactionSortKey(lhs, rhs); };
    IndexArray<4> carrierSearch{ 0, searchIndex, 0, 0 };

    auto range = std::equal_range(
            begin(aggregate.indices), end(aggregate.indices), carrierSearch, sortKeyObj);

    int  otherIndex = (searchIndex == index[0]) ? index[1] : index[0];
    bool integrated = false;
    for (auto it = range.first; it != range.second; ++it)
    {
        int            aggregateIndex = std::distance(begin(aggregate.indices), it);
        const auto&    carrierIndex   = *it;
        AggregateType& aggregateParam = aggregate.parametersA[aggregateIndex];

        if (aggregateParam.manifest & AggregateType::has_bond)
        {
            continue;
        }

        if (otherIndex == carrierIndex[0])
        {
            aggregateParam.manifest |= AggregateType::bond_ij;
            aggregateParam.twoCenter() = twoCSource.parametersA[index[2]];
            integrated                 = true;
            break;
        }
        if (otherIndex == carrierIndex[2])
        {
            aggregateParam.manifest |= AggregateType::bond_jk;
            aggregateParam.twoCenter() = twoCSource.parametersA[index[2]];
            integrated                 = true;
            break;
        }
    }
    return integrated;
}

template<class TwoCenterType, class ThreeCenterType>
void migrateTwoToThree(ListedTypeData<TwoCenterType>&   twoCSource,
                       ListedTypeData<ThreeCenterType>& threeCSource,
                       ListedTypeData<ThreeCenterAggregate<TwoCenterType, ThreeCenterType>>& aggregateDestination)
{
    migrateParameters(threeCSource, aggregateDestination);

    std::vector<int> eliminatedBonds;

    // go through 2-center source and try to integrate into aggregates
    for (size_t bondIndex = 0; bondIndex < twoCSource.indices.size(); ++bondIndex)
    {
        auto index      = twoCSource.indices[bondIndex];
        bool integrated = false;
        integrated      = integrateTwoInThree(index, index[0], twoCSource, aggregateDestination);
        if (!integrated)
        {
            integrated = integrateTwoInThree(index, index[1], twoCSource, aggregateDestination);
        }
        if (integrated)
        {
            eliminatedBonds.push_back(bondIndex);
        }
    }

    deleteInteractions(twoCSource, eliminatedBonds);
}

template<class TwoCenterType, class AggregateType>
bool integrateTwoInFour(IndexArray<3>                        index,
                        int                                  searchIndex,
                        const ListedTypeData<TwoCenterType>& twoCSource,
                        ListedTypeData<AggregateType>&       aggregate)
{
    auto sortKeyObj = [](const auto& lhs, const auto& rhs) { return interactionSortKey(lhs, rhs); };
    IndexArray<5> carrierSearch{ 0, searchIndex, 0, 0 };

    auto range = std::equal_range(
            begin(aggregate.indices), end(aggregate.indices), carrierSearch, sortKeyObj);

    int  otherIndex = (searchIndex == index[0]) ? index[1] : index[0];
    bool integrated = false;
    for (auto it = range.first; it != range.second; ++it)
    {
        int            aggregateIndex = std::distance(begin(aggregate.indices), it);
        const auto&    carrierIndex   = *it;
        AggregateType& aggregateParam = aggregate.parametersA[aggregateIndex];

        if (aggregateParam.manifest & AggregateType::has_bond)
        {
            continue;
        }
        if (otherIndex == carrierIndex[0])
        {
            aggregateParam.manifest |= AggregateType::bond_ij;
            aggregateParam.twoCenter() = twoCSource.parametersA[index[2]];
            integrated                 = true;
            break;
        }
        if (otherIndex == carrierIndex[2])
        {
            aggregateParam.manifest |= AggregateType::bond_jk;
            aggregateParam.twoCenter() = twoCSource.parametersA[index[2]];
            integrated                 = true;
            break;
        }
        if (otherIndex == carrierIndex[3])
        {
            aggregateParam.manifest |= AggregateType::bond_jl;
            aggregateParam.twoCenter() = twoCSource.parametersA[index[2]];
            integrated                 = true;
            break;
        }
    }
    return integrated;
}

template<class ThreeCenterType, class AggregateType>
bool integrateThreeInFour(IndexArray<4>                          index,
                          int                                    searchIndex,
                          const ListedTypeData<ThreeCenterType>& threeCSource,
                          ListedTypeData<AggregateType>&         aggregate)
{
    auto sortKeyObj = [](const auto& lhs, const auto& rhs) { return interactionSortKey(lhs, rhs); };
    IndexArray<5> carrierSearch{ 0, searchIndex, 0, 0, 0 };

    auto range = std::equal_range(
            begin(aggregate.indices), end(aggregate.indices), carrierSearch, sortKeyObj);

    bool integrated = false;
    for (auto it = range.first; it != range.second; ++it)
    {
        int            aggregateIndex = std::distance(begin(aggregate.indices), it);
        const auto&    carrierIndex   = *it;
        AggregateType& aggregateParam = aggregate.parametersA[aggregateIndex];

        if (aggregateParam.manifest & AggregateType::has_angle)
        {
            continue;
        }
        if (searchIndex == index[1]
            && ((index[0] == carrierIndex[0] && index[2] == carrierIndex[2])
                || (index[2] == carrierIndex[0] && index[0] == carrierIndex[2])))
        {
            aggregateParam.manifest |= AggregateType::angle_j;
            aggregateParam.threeCenter() = threeCSource.parametersA[index[3]];
            integrated                   = true;
            break;
        }
        if ((searchIndex == index[0] && (index[1] == carrierIndex[2] && index[2] == carrierIndex[3]))
            || (searchIndex == index[2] && (index[1] == carrierIndex[2] && index[0] == carrierIndex[3])))
        {
            aggregateParam.manifest |= AggregateType::angle_k;
            aggregateParam.threeCenter() = threeCSource.parametersA[index[3]];
            integrated                   = true;
            break;
        }
    }
    return integrated;
}

template<class PairType, class AggregateType>
bool integratePair(IndexArray<3>                     index,
                   int                               searchIndex,
                   const ListedTypeData<PairType>&   pairSource,
                   const std::vector<IndexArray<5>>& aggregateIndices,
                   ListedTypeData<AggregateType>&    aggregates)
{
    auto sortKeyObj = [](const auto& lhs, const auto& rhs) { return lhs[0] < rhs[0]; };

    bool          integrated = false;
    IndexArray<5> carrierSearch{ searchIndex };

    int otherIndex = (searchIndex == index[0]) ? index[1] : index[0];

    auto range = std::equal_range(begin(aggregateIndices), end(aggregateIndices), carrierSearch, sortKeyObj);

    for (auto it = range.first; it != range.second; ++it)
    {
        const auto&    carrierIndex   = *it;
        AggregateType& aggregateParam = aggregates.parametersA[carrierIndex[4]];
        if (aggregateParam.manifest & AggregateType::pair_14)
        {
            continue;
        }
        if (otherIndex == carrierIndex[3])
        {
            aggregateParam.manifest |= AggregateType::pair_14;
            aggregateParam.pair() = pairSource.parametersA[index[2]];
            integrated            = true;
            break;
        }
    }
    return integrated;
}

template<class TwoCenterType, class AggregateType>
void migrateTwoToFour(ListedTypeData<TwoCenterType>& twoCSource, ListedTypeData<AggregateType>& destination)
{
    std::vector<int> eliminatedBonds;
    for (size_t bondIndex = 0; bondIndex < twoCSource.indices.size(); ++bondIndex)
    {
        IndexArray<3> index = twoCSource.indices[bondIndex];

        bool integrated = false;
        integrated      = integrateTwoInFour(index, index[0], twoCSource, destination);
        if (!integrated)
        {
            integrated = integrateTwoInFour(index, index[1], twoCSource, destination);
        }
        if (integrated)
        {
            eliminatedBonds.push_back(bondIndex);
        }
    }
    deleteInteractions(twoCSource, eliminatedBonds);
}

template<class ThreeCenterType, class AggregateType>
void migrateThreeToFour(ListedTypeData<ThreeCenterType>& threeCSource, ListedTypeData<AggregateType>& destination)
{
    std::vector<int> eliminatedAngles;
    for (size_t angleIndex = 0; angleIndex < threeCSource.indices.size(); ++angleIndex)
    {
        auto index = threeCSource.indices[angleIndex];

        bool integrated = false;
        integrated      = integrateThreeInFour(index, index[1], threeCSource, destination);
        if (!integrated)
        {
            integrated = integrateThreeInFour(index, index[0], threeCSource, destination);
        }
        if (!integrated)
        {
            integrated = integrateThreeInFour(index, index[2], threeCSource, destination);
        }
        if (integrated)
        {
            eliminatedAngles.push_back(angleIndex);
        }
    }
    deleteInteractions(threeCSource, eliminatedAngles);
}

template<class PairType, class AggregateType>
void migratePairs(ListedTypeData<PairType>& pairSource, ListedTypeData<AggregateType>& aggregates)
{
    auto aggregateIndices = aggregates.indices;
    auto sortKeyObj       = [](const auto& lhs, const auto& rhs) { return lhs[0] < rhs[0]; };
    std::sort(begin(aggregateIndices), end(aggregateIndices), sortKeyObj);

    std::vector<int> eliminatedPairs;

    for (size_t pairIndex = 0; pairIndex < pairSource.indices.size(); ++pairIndex)
    {
        IndexArray<3> index = pairSource.indices[pairIndex];

        bool integrated = false;
        integrated      = integratePair(index, index[0], pairSource, aggregateIndices, aggregates);
        if (!integrated)
        {
            integrated = integratePair(index, index[1], pairSource, aggregateIndices, aggregates);
        }
        if (integrated)
        {
            eliminatedPairs.push_back(pairIndex);
        }
    }

    deleteInteractions(pairSource, eliminatedPairs);
}

inline void createAggregates(ListedInteractionData& interactions)
{
    sortInteractions(interactions);

    auto migrate = [&interactions](auto& interactionElement) {
        using AggregateType = typename std::decay_t<decltype(interactionElement)>::type;

        auto& bonds      = pickType<typename AggregateType::TwoCenterAggregateType>(interactions);
        auto& angles     = pickType<typename AggregateType::ThreeCenterAggregateType>(interactions);
        auto& dihedrals  = pickType<typename AggregateType::CarrierType>(interactions);
        auto& pairs      = pickType<typename AggregateType::PairAggregateType>(interactions);
        auto& aggregates = pickType<AggregateType>(interactions);

        migrateParameters(dihedrals, aggregates);
        migrateTwoToFour(bonds, aggregates);
        migrateThreeToFour(angles, aggregates);
        migratePairs(pairs, aggregates);
    };
    auto computeIndices = subsetIndices(AggregateTypes{}, AllListedTypes{});
    for_each_tuple(migrate, tieElements(interactions, computeIndices));
}

/*! \brief
 * Combine multiple array elements into a single hash
 *
 * From https://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html
 */
template<int N>
struct ArrayHasher
{
    std::size_t operator()(const util::array<int, N>& a) const
    {
        std::size_t hash = 0;

        for (auto e : a)
        {
            hash ^= std::hash<int>{}(e) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

inline bool hasBond(int config)
{
    return config & 0b1111;
}

inline bool hasAngle(int config)
{
    return config & 0b1111110000;
}

inline bool hasPair(int config)
{
    return config & (0b1111 << 10);
}

inline int setBondIndices(util::array<int, 2> grouping, int config)
{
    assert(!hasBond(config));
    config |= grouping[1];
    config |= (grouping[0] << 2);

    return config;
}

inline util::array<int, 2> getBondIndices(int config)
{
    return { (config & 0b1100) >> 2, config & 0b11 };
}

inline int setAngleIndices(util::array<int, 3> grouping, int config)
{
    assert(!hasAngle(config));
    config |= (grouping[2] << 4);
    config |= (grouping[1] << 6);
    config |= (grouping[0] << 8);

    return config;
}

inline util::array<int, 3> getAngleIndices(int config)
{
    return { (config & 0b1100000000) >> 8, (config & 0b11000000) >> 6, (config & 0b110000) >> 4 };
}

inline int setPairIndices(util::array<int, 2> grouping, int config)
{
    assert(!hasPair(config));
    config |= (grouping[1] << 10);
    config |= (grouping[0] << 12);

    return config;
}

inline util::array<int, 2> getPairIndices(int config)
{
    return { (config & (0b11 << 12)) >> 12, (config & (0b11 << 10)) >> 10 };
}

template<size_t N>
void removeIntegrated(std::vector<IndexArray<N>>& indices)
{
    auto paramPositive = [](auto i) { return i.back() >= 0; };
    int  numRemaining  = std::count_if(begin(indices), end(indices), paramPositive);

    std::vector<IndexArray<N>> remaining;
    remaining.reserve(numRemaining);
    std::copy_if(begin(indices), end(indices), std::back_inserter(remaining), paramPositive);

    swap(indices, remaining);
}

inline void integrateBonds(std::vector<IndexArray<3>>&       twoCSource,
                           const std::vector<IndexArray<5>>& carriers,
                           std::vector<IndexArray<4>>&       aggregates)
{
    std::unordered_map<util::array<int, 2>, int, ArrayHasher<2>> bonds;
    for (size_t i = 0; i < twoCSource.size(); ++i)
    {
        auto          e = twoCSource[i];
        IndexArray<2> key{ e[0], e[1] };
        bonds[key] = i;
    }

    util::array<IndexArray<2>, 12> queries{
        { { 0, 1 }, { 1, 0 }, { 0, 2 }, { 2, 0 }, { 1, 2 }, { 2, 1 }, { 1, 3 }, { 3, 1 }, { 2, 3 }, { 3, 2 }, { 0, 3 }, { 3, 0 } }
    };

    for (size_t i = 0; i < carriers.size(); ++i)
    {
        auto host = carriers[i];
        for (auto grouping : queries)
        {
            auto it = bonds.find({ host[grouping[0]], host[grouping[1]] });
            // if a bond matching these indices exists
            if (it != bonds.end())
            {
                // the index in twoCSource
                int bondIndex = it->second;
                int parameter = twoCSource[bondIndex][2];

                // if bond not already integrated elsewhere
                if (parameter >= 0)
                {
                    // store the bond orientation within the host indices
                    aggregates[i][0] = setBondIndices(grouping, aggregates[i][0]);
                    // store the bond parameter index in the aggregate key
                    aggregates[i][1] = parameter;
                    // set bond as integrated
                    twoCSource[bondIndex][2] = -1;
                    // move to next carrier 4-center index
                    break;
                }
            }
        }
    }

    removeIntegrated(twoCSource);
}

inline void integratePairs(std::vector<IndexArray<3>>&       twoCSource,
                           const std::vector<IndexArray<5>>& carriers,
                           std::vector<IndexArray<4>>&       aggregates)
{
    std::unordered_map<util::array<int, 2>, int, ArrayHasher<2>> pairs;
    for (size_t i = 0; i < twoCSource.size(); ++i)
    {
        auto          e = twoCSource[i];
        IndexArray<2> key{ e[0], e[1] };
        pairs[key] = i;
    }

    util::array<IndexArray<2>, 12> queries{
        { { 0, 3 }, { 3, 0 }, { 0, 1 }, { 1, 0 }, { 0, 2 }, { 2, 0 }, { 1, 2 }, { 2, 1 }, { 1, 3 }, { 3, 1 }, { 2, 3 }, { 3, 2 } }
    };

    for (size_t i = 0; i < carriers.size(); ++i)
    {
        auto host = carriers[i];
        for (auto grouping : queries)
        {
            auto it = pairs.find({ host[grouping[0]], host[grouping[1]] });
            // if a bond matching these indices exists
            if (it != pairs.end())
            {
                // the index in twoCSource
                int pairIndex = it->second;
                int parameter = twoCSource[pairIndex][2];

                // if bond not already integrated elsewhere
                if (parameter >= 0)
                {
                    // store the bond orientation within the host indices
                    aggregates[i][0] = setPairIndices(grouping, aggregates[i][0]);
                    // store the bond parameter index in the aggregate key
                    aggregates[i][3] = parameter;
                    // set bond as integrated
                    twoCSource[pairIndex][2] = -1;
                    // move to next carrier 4-center index
                    break;
                }
            }
        }
    }

    removeIntegrated(twoCSource);
}

inline void integrateAngles(std::vector<IndexArray<4>>&       threeCSource,
                            const std::vector<IndexArray<5>>& carriers,
                            std::vector<IndexArray<4>>&       aggregates)
{
    std::unordered_map<util::array<int, 3>, int, ArrayHasher<3>> angles;
    for (size_t i = 0; i < threeCSource.size(); ++i)
    {
        auto          e = threeCSource[i];
        IndexArray<3> key{ e[0], e[1], e[2] };
        angles[key] = i;
    }

    util::array<IndexArray<3>, 12> queries{ {
            { 1, 0, 2 },
            { 2, 0, 1 },
            { 0, 1, 2 },
            { 2, 1, 0 },
            { 3, 1, 2 },
            { 2, 1, 3 },
            { 0, 2, 1 },
            { 1, 2, 0 },
            { 1, 2, 3 },
            { 3, 2, 1 },
            { 2, 3, 1 },
            { 1, 3, 2 },
    } };

    for (size_t i = 0; i < carriers.size(); ++i)
    {
        auto host = carriers[i];
        for (auto grouping : queries)
        {
            auto it = angles.find({ host[grouping[0]], host[grouping[1]], host[grouping[2]] });
            // if a bond matching these indices exists
            if (it != angles.end())
            {
                // the index in twoCSource
                int bondIndex = it->second;
                int parameter = threeCSource[bondIndex][3];

                // if bond not already integrated elsewhere
                if (parameter >= 0)
                {
                    // store the bond orientation within the host indices
                    aggregates[i][0] = setAngleIndices(grouping, aggregates[i][0]);
                    // store the bond parameter index in the aggregate key
                    aggregates[i][2] = parameter;
                    // set bond as integrated
                    threeCSource[bondIndex][3] = -1;
                    // move to next carrier 4-center index
                    break;
                }
            }
        }
    }

    removeIntegrated(threeCSource);
}

} // namespace nblib
#endif // GMX_NBLIB_LISTED_FORCES_AGGREGATE_TRANSFORMATIONS_HPP
