//
// Created by sebkelle on 19.11.19.
//

#include <numeric>

#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"

#include "atomtype.h"
#include "topology.h"

#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/block.h"

namespace nblib {

TopologyBuilder::TopologyBuilder() : numAtoms_(0)
{}

t_blocka TopologyBuilder::createExclusionsList() const
{
    const auto &moleculesList = molecules_;

    std::vector<gmx::ExclusionBlock> exclusionBlockGlobal;
    exclusionBlockGlobal.reserve(numAtoms_);

    //! compare tuples by comparing the first element
    auto firstLowerThan = [](auto tup1, auto tup2) { return std::get<0>(tup1) < std::get<0>(tup2); };

    size_t atomNumberOffset = 0;
    for (const auto &molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        size_t numMols = std::get<1>(molNumberTuple);

        std::vector<std::tuple<int, int>> exclusionList = molecule.exclusions_;

        //! sorted (by first tuple element as key) list of exclusion pairs
        std::sort(std::begin(exclusionList), std::end(exclusionList), firstLowerThan);

        //! remove duplicates
        auto unique_end = std::unique(std::begin(exclusionList), std::end(exclusionList));
        if (unique_end != std::end(exclusionList))
            printf("[nblib] Warning: exclusionList of molecule \"%s\" contained duplicates",
                    molecule.name_.c_str());

        exclusionList.erase(unique_end, std::end(exclusionList));

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
            auto shiftedExclusions = exclusionBlockPerMolecule;
            //! shift atom numbers by atomNumberOffset
            for (auto& localBlock : shiftedExclusions)
            {
                std::transform(std::begin(localBlock.atomNumber), std::end(localBlock.atomNumber),
                               std::begin(localBlock.atomNumber),
                               [atomNumberOffset](auto i){ return i + atomNumberOffset; });
            }

            std::copy(std::begin(shiftedExclusions), std::end(shiftedExclusions),
                           std::back_inserter(exclusionBlockGlobal));

            atomNumberOffset += molecule.numAtomsInMolecule();
        }
    }

    size_t numberOfExclusions = std::accumulate(std::begin(exclusionBlockGlobal),
                                                std::end(exclusionBlockGlobal), size_t(0),
                                                [](size_t acc, const auto &block)
                                                { return acc + block.atomNumber.size();});

    //! At the very end, convert the exclusionBlockGlobal into
    //! a massive t_blocka and return
    t_blocka tBlockGlobal;
    snew(tBlockGlobal.index, numberOfExclusions + 1);
    snew(tBlockGlobal.a, numberOfExclusions + 1);

    gmx::exclusionBlocksToBlocka(exclusionBlockGlobal, &tBlockGlobal);

    return tBlockGlobal;
}

template <class Extractor>
std::vector<real> TopologyBuilder::extractAtomTypeQuantity(Extractor extractor)
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

                AtomType &atomType = std::get<0>(molecule.atomTypes_[atomTypeName]);
                ret.push_back(extractor(atomType));
            }
        }
    }

    return ret;
}

std::vector<real> TopologyBuilder::extractCharge()
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

                real charge = std::get<1>(molecule.atomTypes_[atomTypeName]);
                ret.push_back(charge);
            }
        }
    }

    return ret;
}

Topology TopologyBuilder::buildTopology()
{
    topology_.excls = createExclusionsList();
    topology_.masses = extractAtomTypeQuantity([](const AtomType &atomType){ return atomType.mass(); });
    topology_.charges = extractCharge();

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
