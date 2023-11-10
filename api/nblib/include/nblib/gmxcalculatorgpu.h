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

#ifndef NBLIB_GMXCALCULATORGPU_H
#define NBLIB_GMXCALCULATORGPU_H

#include <memory>
#include <vector>

#include "gromacs/gpu_utils/devicebuffer.h"

#include "nblib/box.h"
#include "nblib/tpr.h"
#include "nblib/vector.h"

struct DeviceInformation;
class DeviceContext;

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx

namespace nblib
{
struct NBKernelOptions;
class Topology;

class GmxNBForceCalculatorGpu final
{
public:
    GmxNBForceCalculatorGpu(gmx::ArrayRef<int>       particleTypeIdOfAllParticles,
                            gmx::ArrayRef<real>      nonBondedParams,
                            gmx::ArrayRef<real>      charges,
                            gmx::ArrayRef<int64_t>   particleInteractionFlags,
                            gmx::ArrayRef<int>       exclusionRanges,
                            gmx::ArrayRef<int>       exclusionElements,
                            const NBKernelOptions&   options,
                            const DeviceInformation& deviceInfo);

    ~GmxNBForceCalculatorGpu();
    /*
        //! Compute forces & virials and return forces in the forceOutput buffer
        void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                     gmx::ArrayRef<gmx::RVec>       forceOutput,
                     VirialTensor&                  virialOutput) const;
    */

    //! calculates a new pair list based on new coordinates (for every NS step)
    void updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box);

    /*! \brief reorder the input coordinates into nbnxm ordering and return them in coordinateOutput
     *
     * \param coordinateInput    in XYZ format
     * \param coordinateOutput   in XYZQ format (Note: Q)
     */
    void reorder(gmx::ArrayRef<const gmx::RVec> coordinateInput, gmx::ArrayRef<real> coordinateOutput);

    /*! \brief convert from nbnxm ordering back to the original ordering
     *
     * \param[in]  input   coordinates or forces in nbnxm ordering in XYZ format
     * \param[out] output  output in original ordering in XYZ format
     */
    void undoReorder(gmx::ArrayRef<const gmx::RVec> input, gmx::ArrayRef<gmx::RVec> output);

    /*! \brief Compute forces on the device, copying from/to CPU buffers
     *
     * \param[in]    coordinateInput  in XYZ format
     * \param[in]    box              the coordinate bounding box
     * \param[inout] forceOutput      in XYZ format, forces are added to \p forceOutput
     */
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput) const;

    /*! \brief Compute forces based on GPU buffers directly
     *
     * \param[in]    coordinateInput in XYZQ format and in nbnxm ordering, length = nbnxmBufferSize()
     * \param[in]    box             the coordinate bounding box
     * \param[inout] forceOutput     in XYZ format and nbxnxm ordering, length = nbnxmBufferSize()
     *                               forces are added to \p forceOutput
     */
    void compute(DeviceBuffer<Float4> coordinateInput, const Box& box, DeviceBuffer<Float3> forceOutput) const;

    //! \brief return buffer size for the DeviceBuffer compute interface
    [[nodiscard]] std::size_t nbnxmBufferSize() const;

    [[nodiscard]] const DeviceContext& deviceContext() const;
    [[nodiscard]] const DeviceStream&  deviceStream() const;

private:
    //! Private implementation
    class GpuImpl;
    std::unique_ptr<GpuImpl> impl_;
};

//! \brief reorders a scalar real-valued array into the nbnxm ordering defined in \p calculator
void reorderScalarArray(GmxNBForceCalculatorGpu&  calculator,
                        gmx::ArrayRef<const real> input,
                        gmx::ArrayRef<real>       output);

/*! \brief Sets up and returns a GmxForceCalculatorGpu
 *
 * This is a convenience wrapper to drive the underlying builder object
 * with a NBLIB Topology. Note however that the builder is agnostic to the topology
 * type and that it is possible to implement the same functionality based on
 * GROMACS internal data structures, such as mtop.
 *
 */
std::unique_ptr<GmxNBForceCalculatorGpu> setupGmxForceCalculatorGpu(const Topology&        topology,
                                                                    const NBKernelOptions& options,
                                                                    const DeviceInformation& deviceInfo);

} // namespace nblib

#endif // NBLIB_GMXCALCULATORGPU_H
