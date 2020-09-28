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
/*! \inpublicapi \file
 * \brief
 * Implements nblib ForceCalculator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_FORCECALCULATOR_H
#define NBLIB_FORCECALCULATOR_H

#include "nblib/interactions.h"
#include "nblib/kerneloptions.h"
#include "nblib/simulationstate.h"

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx

namespace nblib
{
class NbvSetupUtil;
class GmxForceCalculator;

/*! \brief Setups up and computes forces using gromacs backend.
 *
 * The ForceCalculator uses the data in the SimulationState and NBKernelOptions to opaquely
 * construct all gromacs data structures needed to perform nonbonded force calculations. It is
 * costly to create this object since much of the SimulationState and NBKernelOptions has to be
 * passed to the gromacs backend. However, once constructed, compute can be called repeatedly only
 * paying the cost of the actual nonbonded force calculation. Repeated calls to compute on the same
 * coordinated will always return the same forces (within precision), so the user must update the
 * positions using the forces generated here to advance a simulation. If the coordinates move
 * sufficiently far from their positions at construction time, the efficiency of the calculation
 * will suffer. To alleviate this, the user can call updatePairList.
 *
 */
class ForceCalculator final
{
public:
    ForceCalculator(const SimulationState& system, const NBKernelOptions& options);

    ~ForceCalculator();

    /*! \brief Dispatch the nonbonded force kernels and reduce the forces
     *
     * This function zeroes out all values in the passed in forces buffer, so it can be regarded as
     * an output only param.
     *
     * \param[in] coordinates to be used for the force calculation
     * \param[out] forces buffer to store the output forces
     */
    void compute(gmx::ArrayRef<const Vec3> coordinates, gmx::ArrayRef<Vec3> forces);

    /*! \brief Puts particles on a grid based on bounds specified by the box
     *
     * As compute is called repeatedly, the particles drift apart and the force computation becomes
     * progressively less efficient. Calling this function recomputes the particle-particle pair
     * lists so that computation can proceed efficiently. Should be called around every 100 steps.
     *
     * \param particleInfoAllVdW The types of the particles to be placed on grids
     * \param coordinates The coordinates to be placed on grids
     * \param[in] box The system simulation box
     */
    void updatePairList(gmx::ArrayRef<const int> particleInfoAllVdW,
                        gmx::ArrayRef<Vec3>      coordinates,
                        const Box&               box);

private:
    //! GROMACS force calculator to compute forces
    std::unique_ptr<GmxForceCalculator> gmxForceCalculator_;
};

} // namespace nblib

#endif // NBLIB_FORCECALCULATOR_H
