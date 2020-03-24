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
#ifndef GMX_NBLIB_PARTICLETYPE_H
#define GMX_NBLIB_PARTICLETYPE_H

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "gromacs/math/vectypes.h"

namespace nblib
{
class TopologyBuilder;

using ParticleTypeName = std::string;
using Mass             = real;
using C6               = real;
using C12              = real;

/*! \brief Class that represents the particle type.
 *
 * The particle type is used in lookup tables for masses, non-bonded parameters, etc.
 * Every particle has to assigned an atom type.
 */
class ParticleType
{
public:
    /*! \brief Constructor with explicit name and mass specification.
     *
     * \param[in] name The unique name to reference the particle type.
     * \param[in] mass The mass of the particle of this type.
     */
    ParticleType(ParticleTypeName name, Mass mass);

    /*! \brief Constructor with explicit type specification
     *
     * \param[in] name The unique name to reference the particle type.
     * \param[in] mass The mass of the particle of this type.
     * \param[in] c6   The C6 Lennard-Jones parameter.
     * \param[in] c12  The C12 Lennard-Jones parameter.
     */
    ParticleType(ParticleTypeName name, Mass mass, C6 c6, C12 c12);

    //! Force explicit use of correct types
    template<typename T, typename U>
    ParticleType(T name, U mass) = delete;

    //! Force explicit use of correct types
    template<typename T, typename U, typename V, typename W>
    ParticleType(T name, U mass, V c6, W c12) = delete;

    //! Get the type name
    ParticleTypeName name() const;

    //! Get the mass
    Mass mass() const;

    //! Get the c6 param
    C6 c6() const;

    //! Get the c12 param
    C12 c12() const;

private:
    //! The name
    ParticleTypeName name_;
    //! The mass
    Mass mass_;
    //! The c12 param
    C6 c6_;
    //! The c12 param
    C12 c12_;
};

/*! \brief Comparison operator
 *
 * \param[in] a First type.
 * \param[in] b Second type.
 * \returns If the types are identical.
 */
bool operator==(const ParticleType& a, const ParticleType& b);

} // namespace nblib
#endif // GMX_NBLIB_PARTICLETYPE_H
