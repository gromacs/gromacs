//
// Created by sebkelle on 19.11.19.
//

#include <vector>
#include <array>

#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"

#include "topology.h"

#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/block.h"

namespace nblib {

t_blocka TopologyBuilder::createExclusionsList() const
{
    const auto &moleculesList = molecules_;

    std::vector<gmx::ExclusionBlock> exclusionBlockGlobal;
    exclusionBlockGlobal.reserve(numAtoms_);

    //! compare tuples by comparing the first element
    auto firstLowerThan = [](auto tup1, auto tup2) { return std::get<0>(tup1) < std::get<0>(tup2); };

    for (const auto &molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t numMols = std::get<1>(molNumberTuple);

        std::vector<std::tuple<int, int>> exclusionList = molecule.exclusions_;
        //! sorted (by first tuple element as key) list of exclusion pairs
        std::sort(std::begin(exclusionList), std::end(exclusionList), firstLowerThan);
        std::vector<gmx::ExclusionBlock> exclusionBlockPerMolecule;

        //! initialize pair of iterators delimiting the range of exclusions for
        //! the first atom in the list
        GMX_ASSERT(!exclusionList.empty(), "exclusionList must not be empty\n");
        auto range = std::equal_range(std::begin(exclusionList), std::end(exclusionList),
                                  exclusionList[0], firstLowerThan);
        auto it1 = range.first;
        auto it2 = range.second;

        //! loop over all exclusions in molecule, linear in exclusionList.size()
        while (it1 != std::end(exclusionList))
        {
            gmx::ExclusionBlock localBlock;
            //! loop over all exclusions for current atom
            for ( ; it1 != it2; ++it1)
            {
                localBlock.atomNumber.push_back(std::get<1>(*it1));
            }

            exclusionBlockPerMolecule.push_back(localBlock);

            //! update the upper bound of the range for the next atom
            if (it1 != end(exclusionList))
                it2 = std::upper_bound(it1, std::end(exclusionList), *it1, firstLowerThan);
        }

        //! duplicate the exclusionBlockPerMolecule for the number of Molecules of (numMols)
        for (size_t i = 0; i < numMols; ++i)
        {
            std::copy(std::begin(exclusionBlockPerMolecule), std::end(exclusionBlockPerMolecule),
                      std::back_inserter(exclusionBlockGlobal));
        }
    }

    //! At the very end, convert the exclusionBlockGlobal into
    //! a massive t_blocka and return
    t_blocka tBlockGlobal;
    gmx::exclusionBlocksToBlocka(exclusionBlockGlobal, &tBlockGlobal);

    return tBlockGlobal;
}

template <class Extractor>
std::vector<real> TopologyBuilder::extractQuantity(Extractor extractor)
{
    auto &moleculesList = molecules_;

    //! returned object
    std::vector<real> ret;
    ret.reserve(numAtoms_);

    for (auto &molNumberTuple : moleculesList)
    {
        Molecule &molecule = std::get<0>(molNumberTuple);
        size_t numMols = std::get<1>(molNumberTuple);

        for (size_t i = 0; i < numMols; ++i)
        {
            for (auto &atomTuple : molecule.atoms_)
            {
                std::string atomTypeName = std::get<1>(atomTuple);

                Atom &atomType = molecule.atomTypes_[atomTypeName];
                ret.push_back(extractor(atomType));
            }
        }
    }

    return ret;
}

Topology TopologyBuilder::buildTopology()
{
    topology_.excls = createExclusionsList();
    topology_.masses = extractQuantity([](const Atom &atomType){ return atomType.mass(); });
    topology_.charges = extractQuantity([](const Atom &atomType){ return atomType.charge(); });

    return topology_;
}

TopologyBuilder& TopologyBuilder::addMolecule(const Molecule molecule, const int nMolecules)
{
    /*!
     * 1. Push-back a tuple of molecule type and nMolecules
     * 2. Append exclusion list into the data structure
     */

    molecules_.emplace_back(std::make_tuple(molecule, nMolecules));
    numAtoms_ += nMolecules*molecule.numAtomsInMolecule();

    return *this;
}

const std::vector<real>& Topology::getMasses() const
{
    return masses;
}

const std::vector<real>& Topology::getCharges() const
{
    return charges;
}

const std::vector<int>& Topology::getAtoms() const
{
    return atomTypes;
}

const std::vector<real>& Topology::getNonbondedParameters() const
{
    return nonbondedParameters;
}

const std::vector<int>& Topology::getAtomInfoAllVdw() const
{
   return atomInfoAllVdw;
}



}
