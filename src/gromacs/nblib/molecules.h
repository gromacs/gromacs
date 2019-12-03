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
#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "interactions.h"

#include "gromacs/math/vectypes.h"

class AtomType;
class TopologyBuilder;

namespace nblib
{

class Molecule {
public:
    Molecule(std::string moleculeName);

    Molecule& addAtom(const std::string &atomName, const std::string &residueName, AtomType const &atomType);

    Molecule& addAtom(const std::string &atomName, AtomType const &atomType);

    void addHarmonicBond(HarmonicType harmonicBond);

    // TODO: add exclusions based on the unique ID given to the atom of the molecule
    void addExclusion(const int atomIndex, const int atomIndexToExclude);

    void addExclusion(std::tuple<std::string, std::string> atom, std::tuple<std::string, std::string> atomToExclude);

    void addExclusion(std::string atomName, std::string atomNameToExclude);

    int numAtomsInMolecule() const;

    friend class TopologyBuilder;

private:
    std::string name_;

    //! one entry per atom in molecule
    std::vector<std::tuple<std::string, std::string>> atoms_;
    //! collection of distinct Atoms in molecule
    std::unordered_map<std::string, AtomType> atomTypes_;

    std::vector<std::tuple<int, int>> exclusions_;

    std::vector<HarmonicType> harmonicInteractions_;

    int atomNameAndResidueToIndex(std::tuple<std::string, std::string> atomResNameTuple);

    void addAtomSelfExclusion(std::string atomName, std::string resName);

};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
