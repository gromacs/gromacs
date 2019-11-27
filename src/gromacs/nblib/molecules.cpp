/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Implements nblib Molecule
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "atomtype.h"
#include "molecules.h"

#include <tuple>

namespace nblib {


Molecule::Molecule(std::string moleculeName) : name_(std::move(moleculeName)) {}

void Molecule::addAtomSelfExclusion(std::string atomName, std::string resName)
{
    bool found = false;
    int atomIndex = atomNameAndResidueToIndex(std::make_tuple(atomName, resName));

    for(auto &tuple : exclusions_)
    {
        bool found = false;
        if(std::get<0>(tuple) == atomIndex) {
            if(std::get<1>(tuple) == atomIndex) {
                found = true;
            }
        }
    }
    if(!found) {
        exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndex));
    }
}

Molecule& Molecule::addAtom(const AtomName& atomName, const ResidueName& residueName, const Charge& charge, AtomType const &atomType)
{
    // this check can only ensure that we don't use the same atomName twice, not that an (AtomType, charge) pair is not repeated
    if (!atomTypes_.count(atomName))
    {
        atomTypes_[atomName] = std::make_tuple(atomType, charge);
    }

    atoms_.emplace_back(std::make_tuple(atomName, residueName));
    addAtomSelfExclusion(atomName, residueName);

    return *this;
}

Molecule& Molecule::addAtom(const AtomName& atomName, const ResidueName& residueName, AtomType const &atomType)
{
    real charge = 0;
    addAtom(atomName, residueName, charge, atomType);

    return *this;
}

Molecule& Molecule::addAtom(const AtomName& atomName, const Charge& charge, AtomType const &atomType)
{
    addAtom(atomName, name_, charge, atomType);

    return *this;
}

Molecule& Molecule::addAtom(const AtomName& atomName, const AtomType& atomType)
{
    real charge = 0;
    addAtom(atomName, name_, charge, atomType);

    return *this;
}

int Molecule::numAtomsInMolecule() const
{
    return atoms_.size();
}

void Molecule::addHarmonicBond(HarmonicType harmonicBond)
{
    harmonicInteractions_.push_back(harmonicBond);
}


int Molecule::atomNameAndResidueToIndex(std::tuple<std::string, std::string> atomResNameTuple)
{
    auto equal = [](auto tup1, auto tup2) { return (std::get<0>(tup1) == std::get<0>(tup2) and std::get<1>(tup1) == std::get<1>(tup2)); };

    auto posIter = std::find_if(begin(atoms_), end(atoms_),
        [&](std::tuple<std::string, std::string> & atomAndResName)
        {
            return equal(atomAndResName, atomResNameTuple);
        });

    if(posIter != end(atoms_)) {
        return posIter - begin(atoms_);
    } else {
        // TODO throw exception
        return -1;
    }
}

void Molecule::addExclusion(const int atomIndex, const int atomIndexToExclude)
{
    // We do not need to add exclusion in case the atom indexes are the same
    // because self exclusion are added by addAtom
    if(atomIndex != atomIndexToExclude){
        exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndexToExclude));
        exclusions_.emplace_back(std::make_tuple(atomIndexToExclude, atomIndex));
    }
}

void Molecule::addExclusion(std::tuple<std::string, std::string> atom, std::tuple<std::string, std::string> atomToExclude)
{
    auto atomNameIndex = atomNameAndResidueToIndex(atom);
    auto atomToExcludeIndex = atomNameAndResidueToIndex(atomToExclude);

    addExclusion(atomNameIndex, atomToExcludeIndex);
}

void Molecule::addExclusion(std::string atomName, std::string atomNameToExclude)
{
    auto atomNameIndex = atomNameAndResidueToIndex(std::make_tuple(atomName, name_));
    auto atomToExcludeIndex = atomNameAndResidueToIndex(std::make_tuple(atomNameToExclude, name_));

    addExclusion(atomNameIndex, atomToExcludeIndex);
}

} // namespace nblib
