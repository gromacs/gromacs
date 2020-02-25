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
/*! \file
 * \brief
 * Implements nblib Molecule
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \inpublicapi
 * \ingroup nblib
 */
#ifndef GMX_NBLIB_MOLECULES_H
#define GMX_NBLIB_MOLECULES_H

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/nblib/atomtype.h"

#include "bondtypes.h"

namespace nblib
{
class TopologyBuilder;

using ParticleName = std::string;
using Charge       = real;
using ResidueName  = std::string;

class Molecule
{
public:
    Molecule(std::string moleculeName);

    // Add an atom to the molecule with full specification of parameters.
    Molecule& addAtom(const ParticleName&    particleName,
                      const ResidueName& residueName,
                      const Charge&      charge,
                      AtomType const&    particleType);

    // Force explicit use of correct types
    template<typename T, typename U, typename V>
    Molecule& addAtom(const T& particleName, const U& residueName, const V& charge, AtomType const& particleType) = delete;

    // Add an atom to the molecule with implicit charge of 0
    Molecule& addAtom(const ParticleName& particleName, const ResidueName& residueName, AtomType const& particleType);

    // Add an atom to the molecule with residueName set using particleName
    Molecule& addAtom(const ParticleName& particleName, const Charge& charge, AtomType const& particleType);

    // Force explicit use of correct types, covers both implicit charge and residueName
    template<typename T, typename U>
    Molecule& addAtom(const T& particleName, const U& charge, AtomType const& particleType) = delete;

    // Add an atom to the molecule with residueName set using particleName with implicit charge of 0
    Molecule& addAtom(const ParticleName& particleName, AtomType const& particleType);

    // Force explicit use of correct types
    template<typename T>
    Molecule& addAtom(const T& particleName, AtomType const& particleType) = delete;

    // TODO: add exclusions based on the unique ID given to the atom of the molecule
    void addExclusion(int atomIndex, int atomIndexToExclude);

    // Specify an exclusion with atom and residue names that have been added to molecule
    void addExclusion(std::tuple<std::string, std::string> atom,
                      std::tuple<std::string, std::string> atomToExclude);

    // Specify an exclusion with atoms names that have been added to molecule
    void addExclusion(const std::string& particleName, const std::string& particleNameToExclude);

    void addInteraction(ParticleName particleNameI, ParticleName particleNameJ, HarmonicBondType bondType);

    void addInteraction(ParticleName particleNameI, ParticleName particleNameJ, G96BondType bondType);

    void addInteraction(ParticleName particleNameI, ParticleName particleNameJ, HalfAttractiveQuarticBondType bondType);

    // The number of molecules
    int numAtomsInMolecule() const;

    // Return the AtomType data for a specific atom name that has been added to the molecule
    const AtomType& at(const std::string& particleTypeName) const;

    // convert exclusions given by name to indices and unify with exclusions given by indices
    // returns a sorted vector containing no duplicates of atoms to exclude by indices
    std::vector<std::tuple<int, int>> getExclusions() const;

    friend class TopologyBuilder;

private:
    //! Name of the molecule
    std::string name_;

    struct AtomData
    {
        std::string particleName_;
        std::string residueName_;
        std::string particleTypeName_;
        real        charge_;
    };

    template <class Bond>
    struct BondData {
        std::unordered_map<BondName, Bond> bondTypes_;
        std::vector<std::tuple<ParticleName, ParticleName, BondName>> bonds_;
    };

    // template <class BondTypes...>
    // struct BondTypesData {

    // }

    //! one entry per atom in molecule
    std::vector<AtomData> atoms_;

    //! collection of distinct Atoms in molecule
    std::unordered_map<std::string, AtomType> particleTypes_;

    //! Used for calculated exclusions based on atom indices in molecule
    std::vector<std::tuple<int, int>> exclusions_;

    //! we cannot efficiently compute indices during the build-phase
    //! so we delay the conversion until TopologyBuilder requests it
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> exclusionsByName_;

    BondData<HarmonicBondType> harmonicBonds_;
    BondData<G96BondType> g96Bonds_;
    BondData<HalfAttractiveQuarticBondType> halfAttractiveBonds_;

};

} // namespace nblib
#endif // GMX_NBLIB_MOLECULES_H
