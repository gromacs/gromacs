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
/*! \libinternal \file
 * \brief
 * Implements a force calculator based on GROMACS data structures.
 *
 * Intended for internal use inside the ForceCalculator.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_GMXCALCULATORCPU_H
#define NBLIB_GMXCALCULATORCPU_H

#include <memory>
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
struct NBKernelOptions;
class Topology;
struct TprReader;

class GmxNBForceCalculatorCpu final
{
public:
    GmxNBForceCalculatorCpu(gmx::ArrayRef<int>     particleTypeIdOfAllParticles,
                            gmx::ArrayRef<real>    nonBondedParams,
                            gmx::ArrayRef<real>    charges,
                            gmx::ArrayRef<int64_t> particleInteractionFlags,
                            gmx::ArrayRef<int>     exclusionRanges,
                            gmx::ArrayRef<int>     exclusionElements,
                            const NBKernelOptions& options);

    ~GmxNBForceCalculatorCpu();

    //! calculates a new pair list based on new coordinates (for every NS step)
    void updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box);

    //! Compute forces and return
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput);

    //! Compute forces and virial tensor
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput,
                 gmx::ArrayRef<real>            virialOutput);

    //! Compute forces, virial tensor and potential energies
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput,
                 gmx::ArrayRef<real>            virialOutput,
                 gmx::ArrayRef<real>            energyOutput);

private:
    //! Private implementation
    class CpuImpl;
    std::unique_ptr<CpuImpl> impl_;
};

std::unique_ptr<GmxNBForceCalculatorCpu> setupGmxForceCalculatorCpu(const Topology&        topology,
                                                                    const NBKernelOptions& options);

//! Sets up and returns a GmxForceCalculatorCpu based on a TPR file as input
std::unique_ptr<GmxNBForceCalculatorCpu> setupGmxForceCalculatorCpu(TprReader& tprReader,
                                                                    const NBKernelOptions& options);

} // namespace nblib

#endif // NBLIB_GMXCALCULATORCPU_H
