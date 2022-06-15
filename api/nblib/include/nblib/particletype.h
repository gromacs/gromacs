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
 * Declares nblib ParticleTypes
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
#ifndef NBLIB_PARTICLETYPE_H
#define NBLIB_PARTICLETYPE_H

#include <string>

#include "nblib/basicdefinitions.h"
#include "nblib/util/util.hpp"

namespace nblib
{
class TopologyBuilder;

//! Named type for particle type name
using ParticleTypeName = StrongType<std::string, struct ParticleTypeNameParameter>;
//! Named type for particle mass
using Mass = StrongType<real, struct MassParameter>;

//! Shorthand for a map used for looking up non-bonded parameters using particle types
//! Named type for the C6 parameter in the Lennard-Jones potential
using C6 = StrongType<real, struct C6Parameter>;
//! Named type for the C12 parameter in the Lennard-Jones potential
using C12 = StrongType<real, struct C12Parameter>;

/*! \brief Class that represents the particle type.
 *
 * The particle type is used in lookup tables for masses, non-bonded parameters, etc.
 * Every particle has to assigned an atom type.
 */
class ParticleType final
{
public:
    /*! \brief Constructor with explicit name and mass specification.
     *
     * \param[in] name The unique name to reference the particle type.
     * \param[in] mass The mass of the particle of this type.
     */
    ParticleType(ParticleTypeName name, Mass mass);

    //! Get the type name
    [[nodiscard]] ParticleTypeName name() const;

    //! Get the mass
    [[nodiscard]] Mass mass() const;

private:
    //! The name
    ParticleTypeName name_;
    //! The mass
    Mass mass_;
};

/*! \brief Comparison operator
 *
 * \param[in] a First type.
 * \param[in] b Second type.
 * \returns If the types are identical.
 */
bool operator==(const ParticleType& a, const ParticleType& b);

} // namespace nblib
#endif // NBLIB_PARTICLETYPE_H
