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
 * Implements nblib ParticleType
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include "gmxpre.h"

#include "particletype.h"

namespace nblib
{

ParticleType::ParticleType() noexcept :
    name_(ParticleTypeName("")),
    mass_(Mass(0)),
    c6_(C6(0)),
    c12_(C12(0))
{
}

ParticleType::ParticleType(ParticleTypeName name, Mass mass, C6 c6, C12 c12) :
    name_(std::move(name)),
    mass_(mass),
    c6_(c6),
    c12_(c12)
{
}

ParticleTypeName ParticleType::name() const
{
    return name_;
}

Mass ParticleType::mass() const
{
    return mass_;
}


C6 ParticleType::c6() const
{
    return c6_;
}

C12 ParticleType::c12() const
{
    return c12_;
}

bool operator==(const ParticleType& a, const ParticleType& b)
{
    return a.name() == b.name() && a.mass() == b.mass() && a.c6() == b.c6() && a.c12() == b.c12();
}

} // namespace nblib
