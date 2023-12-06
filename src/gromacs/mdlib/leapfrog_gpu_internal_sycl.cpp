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
/*! \internal \file
 *
 * \brief Implements Leap-Frog using SYCL
 *
 * This file contains SYCL implementation of back-end specific code for Leap-Frog.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/leapfrog_gpu.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/template_mp.h"

#include "leapfrog_gpu_internal.h"

//! \brief Class name for leap-frog kernel
template<gmx::NumTempScaleValues numTempScaleValues, gmx::ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling>
class LeapFrogKernel;

namespace gmx
{

using mode = sycl::access_mode;

/*! \brief Main kernel for the Leap-Frog integrator.
 *
 *  The coordinates and velocities are updated on the GPU. Also saves the intermediate values of the coordinates for
 *   further use in constraints.
 *
 *  Each GPU thread works with a single particle.
 *
 * \tparam        numTempScaleValues                The number of different T-couple values.
 * \tparam        parrinelloRahmanVelocityScaling   The properties of the Parrinello-Rahman velocity scaling matrix.
 * \param[in,out] gm_x                              Coordinates to update upon integration.
 * \param[out]    gm_x0                             A copy of the coordinates before the integration (for constraints).
 * \param[in,out] gm_v                              Velocities to update.
 * \param[in]     gm_f                              Atomic forces.
 * \param[in]     gm_inverseMasses                  Reciprocal masses.
 * \param[in]     dt                                Timestep.
 * \param[in]     gm_lambdas                        Temperature scaling factors (one per group).
 * \param[in]     gm_tempScaleGroups                Mapping of atoms into groups.
 * \param[in]     prVelocityScalingMatrixDiagonal   Diagonal elements of Parrinello-Rahman velocity scaling matrix.
 */
template<NumTempScaleValues numTempScaleValues, ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling>
auto leapFrogKernel(Float3* __restrict__ gm_x,
                    Float3* __restrict__ gm_x0,
                    Float3* __restrict__ gm_v,
                    const Float3* __restrict__ gm_f,
                    const float* __restrict__ gm_inverseMasses,
                    float dt,
                    const float* __restrict__ gm_lambdas /* used iff numTempScaleValues != None */,
                    const unsigned short* __restrict__ gm_tempScaleGroups /* used iff numTempScaleValues == Multiple */,
                    Float3 prVelocityScalingMatrixDiagonal)
{
    return [=](sycl::id<1> itemIdx) {
        const Float3 x    = gm_x[itemIdx];
        const Float3 v    = gm_v[itemIdx];
        const Float3 f    = gm_f[itemIdx];
        const float  im   = gm_inverseMasses[itemIdx];
        const float  imdt = im * dt;

        gm_x0[itemIdx] = x;

        const float lambda = [=]() {
            if constexpr (numTempScaleValues == NumTempScaleValues::None)
            {
                return 1.0F;
            }
            else if constexpr (numTempScaleValues == NumTempScaleValues::Single)
            {
                return gm_lambdas[0];
            }
            else if constexpr (numTempScaleValues == NumTempScaleValues::Multiple)
            {
                const int tempScaleGroup = gm_tempScaleGroups[itemIdx];
                return gm_lambdas[tempScaleGroup];
            }
        }();

        const Float3 prVelocityDelta = [=]() {
            if constexpr (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal)
            {
                return Float3{ prVelocityScalingMatrixDiagonal[0] * v[0],
                               prVelocityScalingMatrixDiagonal[1] * v[1],
                               prVelocityScalingMatrixDiagonal[2] * v[2] };
            }
            else if constexpr (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::No)
            {
                return Float3{ 0, 0, 0 };
            }
            else
            {
                // Other kinds of scaling not yet
                // implemented. Assertions higher in the call stack
                // prevent reaching this code.
                return Float3{ 1e9, 1e9, 1e9 };
            }
        }();

        const Float3 v_new = v * lambda - prVelocityDelta + f * imdt;
        gm_v[itemIdx]      = v_new;
        gm_x[itemIdx]      = x + v_new * dt;
    };
}

//! \brief Leap Frog SYCL kernel launch code.
template<NumTempScaleValues numTempScaleValues, ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling, class... Args>
static void launchLeapFrogKernel(const DeviceStream& deviceStream, int numAtoms, Args&&... args)
{
    // Should not be needed for SYCL2020.
    using kernelNameType = LeapFrogKernel<numTempScaleValues, parrinelloRahmanVelocityScaling>;

    const sycl::range<1> rangeAllAtoms(numAtoms);
    sycl::queue          q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = leapFrogKernel<numTempScaleValues, parrinelloRahmanVelocityScaling>(
                std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(rangeAllAtoms, kernel);
    });
}

//! Convert \p doTemperatureScaling and \p numTempScaleValues to \ref NumTempScaleValues.
static NumTempScaleValues getTempScalingType(bool doTemperatureScaling, int numTempScaleValues)
{
    if (!doTemperatureScaling)
    {
        return NumTempScaleValues::None;
    }
    else if (numTempScaleValues == 1)
    {
        return NumTempScaleValues::Single;
    }
    else if (numTempScaleValues > 1)
    {
        return NumTempScaleValues::Multiple;
    }
    else
    {
        gmx_incons("Temperature coupling was requested with no temperature coupling groups.");
    }
}

/*! \brief Select templated kernel and launch it. */
template<class... Args>
static inline void launchLeapFrogKernel(NumTempScaleValues              tempScalingType,
                                        ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
                                        Args&&... args)
{
    GMX_ASSERT(parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::No
                       || parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal,
               "Only isotropic Parrinello-Rahman pressure coupling is supported.");

    dispatchTemplatedFunction(
            [&](auto tempScalingType_, auto prScalingType_) {
                return launchLeapFrogKernel<tempScalingType_, prScalingType_>(std::forward<Args>(args)...);
            },
            tempScalingType,
            parrinelloRahmanVelocityScaling);
}

void launchLeapFrogKernel(int                                   numAtoms,
                          DeviceBuffer<Float3>                  d_x,
                          DeviceBuffer<Float3>                  d_x0,
                          DeviceBuffer<Float3>                  d_v,
                          const DeviceBuffer<Float3>            d_f,
                          const DeviceBuffer<float>             d_inverseMasses,
                          const float                           dt,
                          const bool                            doTemperatureScaling,
                          const int                             numTempScaleValues,
                          const DeviceBuffer<unsigned short>    d_tempScaleGroups,
                          const DeviceBuffer<float>             d_lambdas,
                          const ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
                          const Float3                          prVelocityScalingMatrixDiagonal,
                          const DeviceStream&                   deviceStream)
{
    NumTempScaleValues tempParrinelloRahmanVelocityScaling =
            getTempScalingType(doTemperatureScaling, numTempScaleValues);


    launchLeapFrogKernel(tempParrinelloRahmanVelocityScaling,
                         parrinelloRahmanVelocityScaling,
                         deviceStream,
                         numAtoms,
                         d_x.get_pointer(),
                         d_x0.get_pointer(),
                         d_v.get_pointer(),
                         d_f.get_pointer(),
                         d_inverseMasses.get_pointer(),
                         dt,
                         d_lambdas.get_pointer(),
                         d_tempScaleGroups.get_pointer(),
                         prVelocityScalingMatrixDiagonal);
}

} // namespace gmx
