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

#include "nblib/listed_forces/definitions.h"
#include "nblib/particletype.h"

namespace nblib
{
class TopologyBuilder;

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

//! \brief uniquely identifies a particle within a Molecule
class ParticleIdentifier final
{
public:
    //! \brief construct form a ParticleName, allow implicit conversion
    ParticleIdentifier(ParticleName particleName) :
        particleName_(std::move(particleName)), residueName_()
    {
    }

    //! \brief construct with a non-default ResidueName
    ParticleIdentifier(ParticleName particleName, ResidueName residueName) :
        particleName_(std::move(particleName)), residueName_(std::move(residueName))
    {
    }

    [[nodiscard]] const ParticleName& particleName() const { return particleName_; }
    [[nodiscard]] const ResidueName&  residueName() const { return residueName_; }

private:
    ParticleName particleName_;
    ResidueName  residueName_;

    friend inline bool operator==(const ParticleIdentifier& lhs, const ParticleIdentifier& rhs)
    {
        return lhs.particleName_ == rhs.particleName_ && lhs.residueName_ == rhs.residueName_;
    }
};

//! \brief Molecule class that holds particles and their bonded interactions
class Molecule final
{
    //! \brief per-InteractionType storage container for listed interaction with string-based particle IDs
    template<class InteractionType>
    struct InteractionTypeData
    {
        using type           = InteractionType;
        using IdentifierType = Repeat<TypeList<ParticleName, ResidueName>, NCenter<InteractionType>{}>;

        std::vector<InteractionType>                    interactionTypes_;
        std::vector<Reduce<std::tuple, IdentifierType>> interactions_;
    };

    //! this creates a tuple containing an instance of InteractionType data for each supported listed type
    using InteractionTuple = Reduce<std::tuple, Map<InteractionTypeData, SupportedListedTypes>>;

    //! \brief returns the default residue name if necessary
    ResidueName residueName(const ParticleIdentifier& particleIdentifier);

    //! \brief adds an interaction to the InteractionTuple
    template<class ListedVariant, class... ParticleIdentifiers>
    void addInteractionImpl(const ListedVariant& interaction, const ParticleIdentifiers&... particles);

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

    //! \brief Add a particle to the molecule with residueName set using particleName
    Molecule& addParticle(const ParticleName& particleName, const Charge& charge, ParticleType const& particleType);

    //! \brief Add a particle to the molecule with residueName set using particleName with implicit charge of 0
    Molecule& addParticle(const ParticleName& particleName, const ParticleType& particleType);

    /*! \brief Specify an exclusion between two particles that have been added to the molecule
     *
     * Exclusion of a particle with itself is detected and handled correctly.
     * Note that adding an exclusion between particles not present in the Molecule will \throw an
     * exception.
     */
    void addExclusion(const ParticleIdentifier& particle, const ParticleIdentifier& particleToExclude);

    /*! \brief Add 2-particle interactions such as harmonic bonds
     *
     * Note that adding an interaction between particles not present in the Molecule will \throw an
     * exception.
     */
    void addInteraction(const ParticleIdentifier&   particleI,
                        const ParticleIdentifier&   particleJ,
                        const TwoCenterInteraction& interaction);

    //! \brief Add 3-particle interactions such as angles
    void addInteraction(const ParticleIdentifier&     particleI,
                        const ParticleIdentifier&     particleJ,
                        const ParticleIdentifier&     particleK,
                        const ThreeCenterInteraction& interaction);

    //! \brief Add 4-particle interactions such as (im)proper-dihedrals
    void addInteraction(const ParticleIdentifier&    particleI,
                        const ParticleIdentifier&    particleJ,
                        const ParticleIdentifier&    particleK,
                        const ParticleIdentifier&    particleL,
                        const FourCenterInteraction& interaction);

    //! \brief Add 5-particle interactions such as CMAP
    void addInteraction(const ParticleIdentifier&    particleI,
                        const ParticleIdentifier&    particleJ,
                        const ParticleIdentifier&    particleK,
                        const ParticleIdentifier&    particleL,
                        const ParticleIdentifier&    particleM,
                        const FiveCenterInteraction& interaction);

    //! \brief The number of molecules
    int numParticlesInMolecule() const;

    //! \brief Return the ParticleType data for a specific particle name that has been added to the molecule
    const ParticleType& at(const std::string& particlesTypeName) const;

    /*! \brief access integer-based exclusions
     *
     * Convert exclusions given by name to indices and unify with exclusions given by indices
     * returns a sorted vector containing no duplicates of particles to exclude by indices
     */
    std::vector<std::tuple<int, int>> getExclusions() const;

    //! \brief Return all interactions stored in Molecule
    const InteractionTuple& interactionData() const;

    //! \brief Return name of i-th particle
    ParticleName particleName(int i) const;

    //! \brief Return name of i-th residue
    ResidueName residueName(int i) const;

    //! \brief Return array of data structs on particle types
    std::vector<ParticleData> particleData() const;

    //! \brief Return map of particle types and their names
    std::unordered_map<std::string, ParticleType> particleTypesMap() const;

    //! \brief The molecule name
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

    //! collection of data for all types of interactions
    InteractionTuple interactionData_;
};

/*! \brief Helper function to add harmonic angle and harmonic bonds for Urey-Bradley term.
 *
 * Urey-Bradley consist of two harmonic terms:
 *   1. Harmonic angle, connecting all three particles.
 *   2. Harmonic correction to the distance between two non-central particles (particles 1 and 3)
 * This function creates theese terms and adds them to the \c molecule as independent harmonic
 * angle and harmonic bond.
 *
 * \todo This should be moved to another location (e.g. to TPR reader).
 *
 * \param[in,out] molecule   The molecule to add Urey-bradley to.
 * \param[in]     particleI  First interacting particle.
 * \param[in]     particleJ  Second (central) interacting particle.
 * \param[in]     particleK  Third interacting particle.
 * \param[in]     theta0     Equilibrium angle (in radians).
 * \param[in]     kTheta     Force-constant for angle.
 * \param[in]     r130       Equilibrium distance between particles 1 and 3.
 * \param[in]     kUB        Force constant for bond correction term.
 */
static inline void addUreyBradleyInteraction(Molecule&                 molecule,
                                             const ParticleIdentifier& particleI,
                                             const ParticleIdentifier& particleJ,
                                             const ParticleIdentifier& particleK,
                                             const Radians             theta0,
                                             const ForceConstant       kTheta,
                                             const EquilConstant       r130,
                                             const ForceConstant       kUB)
{
    HarmonicAngle    ubAngle(kTheta, theta0);
    HarmonicBondType ubBond(kUB, r130);
    molecule.addInteraction(particleI, particleJ, particleK, ubAngle);
    molecule.addInteraction(particleI, particleK, ubBond);
}

} // namespace nblib
#endif // NBLIB_MOLECULES_H
