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
 * Implements nblib integrator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_INTEGRATOR_H
#define NBLIB_INTEGRATOR_H

#include <vector>

#include "nblib/box.h"
#include "nblib/vector.h"

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx

namespace nblib
{

class Topology;

/*! \brief Simple integrator
 *
 */
class LeapFrog final
{
public:
    /*! \brief Constructor.
     *
     * \param[in] topology  Topology object to build list of inverse masses.
     * \param[in] box       Box object for ensuring that coordinates remain within bounds
     */
    LeapFrog(const Topology& topology, const Box& box);

    /*! \brief Constructor.
     *
     * \param[in] masses  List of inverse masses.
     * \param[in] box     Box object for ensuring that coordinates remain within bounds
     */
    LeapFrog(gmx::ArrayRef<const real> inverseMasses, const Box& box);

    /*! \brief Integrate
     *
     * Integrates the equation of motion using Leap-Frog algorithm.
     * Updates coordinates and velocities.
     *
     * \param[in]  dt          Timestep.
     * \param[out] coordinates Coordinate array that would be modified in-place.
     * \param[out] velocities  Velocity array that would be modified in-place.
     * \param[in]  forces      Force array to be read.
     *
     */
    void integrate(real                      dt,
                   gmx::ArrayRef<Vec3>       coordinates,
                   gmx::ArrayRef<Vec3>       velocities,
                   gmx::ArrayRef<const Vec3> forces);

private:
    //! 1/mass for all atoms
    std::vector<real> inverseMasses_;
    //! Box for PBC conformity
    Box box_;
};

} // namespace nblib

#endif // NBLIB_INTEGRATOR_H
