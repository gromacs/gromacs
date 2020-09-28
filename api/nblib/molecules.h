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
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inpublicapi
 * \ingroup nblib
 */
#ifndef NBLIB_MOLECULES_H
#define NBLIB_MOLECULES_H

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "nblib/particletype.h"

namespace nblib
{
//! Named type for unique identifier for a particle in a molecule
using ParticleName = StrongType<std::string, struct ParticleNameParameter>;

//! Named type for charges on a particle within a molecule
using Charge = StrongType<real, struct ChargeParameter>;

//! Named type for residue name used to diffentiate between sections of a molecule
using ResidueName = StrongType<std::string, struct ResidueNameParameter>;

//! Named type for the name of a molecule
using MoleculeName = StrongType<std::string, struct MoleculeNameParameter>;

struct ParticleData
{
    std::string particleName_;
    std::string residueName_;
    std::string particleTypeName_;
    real        charge_;
};

class Molecule final
{
public:
    explicit Molecule(MoleculeName moleculeName);

    //! Add a particle to the molecule with full specification of parameters.
    Molecule& addParticle(const ParticleName& particleName,
                          const ResidueName&  residueName,
                          const Charge&       charge,
                          ParticleType const& particleType);

    //! Add a particle to the molecule with implicit charge of 0
    Molecule& addParticle(const ParticleName& particleName,
                          const ResidueName&  residueName,
                          ParticleType const& particleType);

    //! Add a particle to the molecule with residueName set using particleName
    Molecule& addParticle(const ParticleName& particleName, const Charge& charge, ParticleType const& particleType);

    //! Add a particle to the molecule with residueName set using particleName with implicit charge of 0
    Molecule& addParticle(const ParticleName& particleName, ParticleType const& particleType);

    // TODO: add exclusions based on the unique ID given to the particle of the molecule
    void addExclusion(int particleIndex, int particleIndexToExclude);

    //! Specify an exclusion with particle and residue names that have been added to molecule
    void addExclusion(std::tuple<std::string, std::string> particle,
                      std::tuple<std::string, std::string> particleToExclude);

    //! Specify an exclusion with particle names that have been added to molecule
    void addExclusion(const std::string& particleName, const std::string& particleNameToExclude);

    //! The number of molecules
    int numParticlesInMolecule() const;

    //! Return the ParticleType data for a specific particle name that has been added to the molecule
    const ParticleType& at(const std::string& particlesTypeName) const;

    //! Convert exclusions given by name to indices and unify with exclusions given by indices
    //! returns a sorted vector containing no duplicates of particles to exclude by indices
    std::vector<std::tuple<int, int>> getExclusions() const;

    //! Return name of ith particle
    ParticleName particleName(int i) const;

    //! Return name of ith residue
    ResidueName residueName(int i) const;

    //! Return array of data structs on particle types
    std::vector<ParticleData> particleData() const;

    //! Return map of particle types and their names
    std::unordered_map<std::string, ParticleType> particleTypesMap() const;

    //! The molecule name
    MoleculeName name() const;

private:
    //! Name of the molecule
    MoleculeName name_;

    //! one entry per particle in molecule
    std::vector<ParticleData> particles_;

    //! collection of distinct particle types in molecule
    std::unordered_map<std::string, ParticleType> particleTypes_;

    //! Used for calculated exclusions based on particle indices in molecule
    std::vector<std::tuple<int, int>> exclusions_;

    //! we cannot efficiently compute indices during the build-phase
    //! so we delay the conversion until TopologyBuilder requests it
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> exclusionsByName_;
};

} // namespace nblib
#endif // NBLIB_MOLECULES_H
