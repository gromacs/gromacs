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
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "topology.h"

#include <numeric>

#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/util.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

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

    //! initialize pair of iterators delimiting the range of exclusions for
    //! the first atom in the list
    GMX_ASSERT(!tupleList.empty(), "tupleList must not be empty\n");
    auto range = std::equal_range(std::begin(tupleList), std::end(tupleList), tupleList[0], firstLowerThan);
    auto it1 = range.first;
    auto it2 = range.second;

    //! loop over all exclusions in molecule, linear in tupleList.size()
    while (it1 != std::end(tupleList))
    {
        gmx::ExclusionBlock localBlock;
        //! loop over all exclusions for current atom
        for (; it1 != it2; ++it1)
        {
            localBlock.atomNumber.push_back(std::get<1>(*it1));
        }

        ret.push_back(localBlock);

        //! update the upper bound of the range for the next atom
        if (it1 != end(tupleList))
        {
            it2 = std::upper_bound(it1, std::end(tupleList), *it1, firstLowerThan);
        }
    }

    return ret;
}

std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset)
{
    //! shift atom numbers by offset
    for (auto& localBlock : inBlock)
    {
        std::transform(std::begin(localBlock.atomNumber), std::end(localBlock.atomNumber),
                       std::begin(localBlock.atomNumber), [offset](auto i) { return i + offset; });
    }

    return inBlock;
}

} // namespace detail

TopologyBuilder::TopologyBuilder() : numAtoms_(0) {}

gmx::ListOfLists<int> TopologyBuilder::createExclusionsListOfLists() const
{
    const auto& moleculesList = molecules_;

    std::vector<gmx::ExclusionBlock> exclusionBlockGlobal;
    exclusionBlockGlobal.reserve(numAtoms_);

    size_t atomNumberOffset = 0;
    for (const auto& molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t          numMols  = std::get<1>(molNumberTuple);

        std::vector<gmx::ExclusionBlock> exclusionBlockPerMolecule =
                detail::toGmxExclusionBlock(molecule.getExclusions());

        //! duplicate the exclusionBlockPerMolecule for the number of Molecules of (numMols)
        for (size_t i = 0; i < numMols; ++i)
        {
            auto offsetExclusions = detail::offsetGmxBlock(exclusionBlockPerMolecule, atomNumberOffset);

            std::copy(std::begin(offsetExclusions), std::end(offsetExclusions),
                      std::back_inserter(exclusionBlockGlobal));

            atomNumberOffset += molecule.numAtomsInMolecule();
        }
    }

    gmx::ListOfLists<int> exclusionsListOfListsGlobal;
    for (const auto& block : exclusionBlockGlobal)
    {
        exclusionsListOfListsGlobal.pushBack(block.atomNumber);
    }

    return exclusionsListOfListsGlobal;
}

template<typename T, class Extractor>
std::vector<T> TopologyBuilder::extractAtomTypeQuantity(Extractor extractor)
{
    auto& moleculesList = molecules_;

    //! returned object
    std::vector<T> ret;
    ret.reserve(numAtoms_);

    for (auto& molNumberTuple : moleculesList)
    {
        Molecule& molecule = std::get<0>(molNumberTuple);
        size_t    numMols  = std::get<1>(molNumberTuple);

        for (size_t i = 0; i < numMols; ++i)
        {
            for (auto& atomData : molecule.atoms_)
            {
                ret.push_back(extractor(atomData, molecule.particleTypes_));
            }
        }
    }

    return ret;
}

Topology TopologyBuilder::buildTopology()
{
    topology_.numAtoms_ = numAtoms_;

    topology_.exclusions_ = createExclusionsListOfLists();
    topology_.charges_    = extractAtomTypeQuantity<real>([](const auto& data, auto& map) {
        ignore_unused(map);
        return data.charge_;
    });

    std::unordered_map<std::string, int> nameToId;
    for (auto& name_atomType_tuple : atomTypes_)
    {
        topology_.atomTypes_.push_back(name_atomType_tuple.second);
        nameToId[name_atomType_tuple.first] = nameToId.size();
    }

    topology_.atomTypeIdOfAllAtoms_ =
            extractAtomTypeQuantity<int>([&nameToId](const auto& data, auto& map) {
                ignore_unused(map);
                return nameToId[data.particleTypeName_];
            });

    return topology_;
}

TopologyBuilder& TopologyBuilder::addMolecule(const Molecule& molecule, const int nMolecules)
{
    /*!
     * 1. Push-back a tuple of molecule type and nMolecules
     * 2. Append exclusion list into the data structure
     */

    molecules_.emplace_back(std::make_tuple(molecule, nMolecules));
    numAtoms_ += nMolecules * molecule.numAtomsInMolecule();

    for (const auto& name_type_tuple : molecule.particleTypes_)
    {
        //! If we already have the atomType, we need to make
        //! sure that the type's parameters are actually the same
        //! otherwise we would overwrite them
        if (atomTypes_.count(name_type_tuple.first) > 0)
        {
            if (!(atomTypes_[name_type_tuple.first] == name_type_tuple.second))
            {
                GMX_THROW(gmx::InvalidInputError(
                        "Differing AtomTypes with identical names encountered"));
            }
        }
    }

    // Note: insert does nothing if the key already exists
    atomTypes_.insert(molecule.particleTypes_.begin(), molecule.particleTypes_.end());

    return *this;
}

const int& Topology::numAtoms() const
{
    return numAtoms_;
}

const std::vector<real>& Topology::getCharges() const
{
    return charges_;
}

const std::vector<ParticleType>& Topology::getAtomTypes() const
{
    return atomTypes_;
}

const std::vector<int>& Topology::getAtomTypeIdOfAllAtoms() const
{
    return atomTypeIdOfAllAtoms_;
}

} // namespace nblib