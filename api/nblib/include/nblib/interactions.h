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
/*! \inpublicapi \file
 * \brief
 * Implements nblib particle-types interactions
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_INTERACTIONS_H
#define NBLIB_INTERACTIONS_H

#include <cstddef>

#include <map>
#include <string>
#include <tuple>
#include <unordered_map>

#include "nblib/kerneloptions.h"
#include "nblib/particletype.h"
#include "nblib/util/util.hpp"

namespace nblib
{

using NonBondedInteractionMapImpl =
        std::map<std::tuple<ParticleTypeName, ParticleTypeName>, std::tuple<C6, C12>>;

//! Map used for looking up non-bonded parameters using particle types
class NonBondedInteractionMap final
{
private:
    using NamePairTuple   = std::tuple<ParticleTypeName, ParticleTypeName>;
    using ComboParameters = std::tuple<C6, C12>;
    using InteractionMap  = std::map<NamePairTuple, ComboParameters>;
    InteractionMap interactionMap_;

public:
    void   setInteractions(const ParticleTypeName&, const ParticleTypeName&, const C6, const C12);
    size_t count(const NamePairTuple&);

    [[nodiscard]] C6  getC6(const ParticleTypeName&, const ParticleTypeName&) const;
    [[nodiscard]] C12 getC12(const ParticleTypeName&, const ParticleTypeName&) const;

    InteractionMap::iterator begin() { return interactionMap_.begin(); }
    InteractionMap::iterator end() { return interactionMap_.end(); }
};

/*! \brief Non-Bonded Interactions between Particle Types
 *
 * \inpublicapi
 * \ingroup nblib
 *
 * A class to hold a mapping between pairs of particle types and the non-bonded
 * interactions between them. One may add the non-bonded parameters, namely the
 * C6/C12 params for each particle type individually and construct a pair-wise
 * mapping using combination rules or manually specify the parameters between
 * a specific pair.
 *
 */
class ParticleTypesInteractions final
{
public:
    //! Initialized with the default geometric combination rule
    explicit ParticleTypesInteractions(CombinationRule = CombinationRule::Geometric);

    //! Specify non-bonded params of a particle type
    ParticleTypesInteractions& add(const ParticleTypeName& particleTypeName, C6 c6, C12 c12);

    //! Specify the non-bonded params of a specific pair of particle types
    ParticleTypesInteractions& add(const ParticleTypeName& particleTypeName1,
                                   const ParticleTypeName& particleTypeName2,
                                   C6                      c6,
                                   C12                     c12);

    //! Generate table based on the parameters stored
    [[nodiscard]] NonBondedInteractionMap generateTable() const;

    //! Get combination rule enabled in this object
    [[nodiscard]] CombinationRule getCombinationRule() const;

    //! Merge with the information stored in another ParticleTypesInteractions object
    void merge(const ParticleTypesInteractions&);

private:
    CombinationRule combinationRule_;

    std::map<ParticleTypeName, std::tuple<C6, C12>> singleParticleInteractionsMap_;
    std::map<std::tuple<ParticleTypeName, ParticleTypeName>, std::tuple<C6, C12>> twoParticlesInteractionsMap_;
};

} // namespace nblib
#endif // NBLIB_INTERACTIONS_H
