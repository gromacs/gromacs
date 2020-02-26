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
 * Implements nblib Molecule
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "molecules.h"

#include <tuple>

#include "gromacs/nblib/particletype.h"

namespace nblib
{


Molecule::Molecule(std::string moleculeName) : name_(std::move(moleculeName)) {}

Molecule& Molecule::addParticle(const ParticleName& particleName,
                                const ResidueName&  residueName,
                                const Charge&       charge,
                                ParticleType const& particleType)
{
    if (particleTypes_.count(particleType.name()) == 0)
    {
        particleTypes_[particleType.name()] = particleType;
    }

    particles_.emplace_back(ParticleData{ particleName, residueName, particleType.name(), charge });

    //! Add self exclusion. We just added the particle, so we know its index and that the exclusion doesn't exist yet
    std::size_t id = particles_.size() - 1;
    exclusions_.emplace_back(std::make_tuple(id, id));

    return *this;
}

Molecule& Molecule::addParticle(const ParticleName& particleName,
                                const ResidueName&  residueName,
                                ParticleType const& particleType)
{
    real charge = 0;
    addParticle(particleName, residueName, charge, particleType);

    return *this;
}

Molecule& Molecule::addParticle(const ParticleName& particleName,
                                const Charge&       charge,
                                ParticleType const& particleType)
{
    addParticle(particleName, name_, charge, particleType);

    return *this;
}

Molecule& Molecule::addParticle(const ParticleName& particleName, const ParticleType& particleType)
{
    real charge = 0;
    addParticle(particleName, name_, charge, particleType);

    return *this;
}

int Molecule::numParticlesInMolecule() const
{
    return particles_.size();
}

void Molecule::addExclusion(const int particleIndex, const int particleIndexToExclude)
{
    // We do not need to add exclusion in case the particle indexes are the same
    // because self exclusion are added by addParticle
    if (particleIndex != particleIndexToExclude)
    {
        exclusions_.emplace_back(particleIndex, particleIndexToExclude);
        exclusions_.emplace_back(particleIndexToExclude, particleIndex);
    }
}

void Molecule::addExclusion(std::tuple<std::string, std::string> particle,
                            std::tuple<std::string, std::string> particleToExclude)
{
    //! duplication for the swapped pair happens in getExclusions()
    exclusionsByName_.emplace_back(std::make_tuple(std::get<0>(particle), std::get<1>(particle),
                                                   std::get<0>(particleToExclude),
                                                   std::get<1>(particleToExclude)));
}

void Molecule::addExclusion(const std::string& particleName, const std::string& particleNameToExclude)
{
    addExclusion(std::make_tuple(particleName, name_), std::make_tuple(particleNameToExclude, name_));
}

const Molecule::InteractionTuple& Molecule::interactionData() const { return interactionData_; }

const ParticleType& Molecule::at(const std::string& particleTypeName) const
{
    return particleTypes_.at(particleTypeName);
}

std::vector<std::tuple<int, int>> Molecule::getExclusions() const
{
    //! tuples of (particleName, residueName, index)
    std::vector<std::tuple<std::string, std::string, int>> indexKey;
    indexKey.reserve(numParticlesInMolecule());

    for (int i = 0; i < numParticlesInMolecule(); ++i)
    {
        indexKey.emplace_back(particles_[i].particleName_, particles_[i].residueName_, i);
    }

    std::sort(std::begin(indexKey), std::end(indexKey));

    std::vector<std::tuple<int, int>> ret = exclusions_;
    ret.reserve(exclusions_.size() + exclusionsByName_.size());

    //! normal operator<, except ignore third element
    auto sortKey = [](const auto& tup1, const auto& tup2) {
        if (std::get<0>(tup1) < std::get<0>(tup2))
        {
            return true;
        }
        else
        {
            return std::get<1>(tup1) < std::get<1>(tup2);
        }
    };

    //! convert exclusions given by names to indices and append
    for (auto& tup : exclusionsByName_)
    {
        const std::string& particleName1 = std::get<0>(tup);
        const std::string& residueName1  = std::get<1>(tup);
        const std::string& particleName2 = std::get<2>(tup);
        const std::string& residueName2  = std::get<3>(tup);

        //! look up first index (binary search)
        auto it1 = std::lower_bound(std::begin(indexKey), std::end(indexKey),
                                    std::make_tuple(particleName1, residueName2, 0), sortKey);

        //! make sure we have the (particleName,residueName) combo
        if (it1 == std::end(indexKey) or std::get<0>(*it1) != particleName1 or std::get<1>(*it1) != residueName1)
        {
            throw std::runtime_error((std::string("Particle ") += particleName1 + std::string(" in residue ") +=
                                      residueName1 + std::string(" not found in list of particles\n"))
                                             .c_str());
        }

        int firstIndex = std::get<2>(*it1);

        //! look up second index (binary search)
        auto it2 = std::lower_bound(std::begin(indexKey), std::end(indexKey),
                                    std::make_tuple(particleName2, residueName2, 0), sortKey);

        //! make sure we have the (particleName,residueName) combo
        if (it2 == std::end(indexKey) or std::get<0>(*it2) != particleName2 or std::get<1>(*it2) != residueName2)
        {
            throw std::runtime_error((std::string("Particle ") += particleName2 + std::string(" in residue ") +=
                                      residueName2 + std::string(" not found in list of particles\n"))
                                             .c_str());
        }

        int secondIndex = std::get<2>(*it2);

        ret.emplace_back(firstIndex, secondIndex);
        ret.emplace_back(secondIndex, firstIndex);
    }

    std::sort(std::begin(ret), std::end(ret));

    auto uniqueEnd = std::unique(std::begin(ret), std::end(ret));
    if (uniqueEnd != std::end(ret))
    {
        printf("[nblib] Warning: exclusionList for molecule %s contained duplicates", name_.c_str());
    }

    ret.erase(uniqueEnd, std::end(ret));
    return ret;
}

} // namespace nblib
