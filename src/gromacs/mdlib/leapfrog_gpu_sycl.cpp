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
 *
 * \brief Implements Leap-Frog using SYCL
 *
 * This file contains implementation of basic Leap-Frog integrator
 * using SYCL, including class initialization, data-structures management
 * and GPU kernel.
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

namespace gmx
{

using cl::sycl::access::mode;

/*! \brief Main kernel for the Leap-Frog integrator.
 *
 *  The coordinates and velocities are updated on the GPU. Also saves the intermediate values of the coordinates for
 *   further use in constraints.
 *
 *  Each GPU thread works with a single particle.
 *
 * \tparam        numTempScaleValues               The number of different T-couple values.
 * \tparam        velocityScaling                  Type of the Parrinello-Rahman velocity rescaling.
 * \param         cgh                              SYCL's command group handler.
 * \param[in,out] a_x                              Coordinates to update upon integration.
 * \param[out]    a_xp                             A copy of the coordinates before the integration (for constraints).
 * \param[in,out] a_v                              Velocities to update.
 * \param[in]     a_f                              Atomic forces.
 * \param[in]     a_inverseMasses                  Reciprocal masses.
 * \param[in]     dt                               Timestep.
 * \param[in]     a_lambdas                        Temperature scaling factors (one per group).
 * \param[in]     a_tempScaleGroups                Mapping of atoms into groups.
 * \param[in]     prVelocityScalingMatrixDiagonal  Diagonal elements of Parrinello-Rahman velocity scaling matrix
 */
template<NumTempScaleValues numTempScaleValues, VelocityScalingType velocityScaling>
auto leapFrogKernel(
        cl::sycl::handler&                          cgh,
        DeviceAccessor<float3, mode::read_write>    a_x,
        DeviceAccessor<float3, mode::discard_write> a_xp,
        DeviceAccessor<float3, mode::read_write>    a_v,
        DeviceAccessor<float3, mode::read>          a_f,
        DeviceAccessor<float, mode::read>           a_inverseMasses,
        float                                       dt,
        OptionalAccessor<float, mode::read, numTempScaleValues != NumTempScaleValues::None> a_lambdas,
        OptionalAccessor<unsigned short, mode::read, numTempScaleValues == NumTempScaleValues::Multiple> a_tempScaleGroups,
        float3 prVelocityScalingMatrixDiagonal)
{
    cgh.require(a_x);
    cgh.require(a_xp);
    cgh.require(a_v);
    cgh.require(a_f);
    cgh.require(a_inverseMasses);
    if constexpr (numTempScaleValues != NumTempScaleValues::None)
    {
        cgh.require(a_lambdas);
    }
    if constexpr (numTempScaleValues == NumTempScaleValues::Multiple)
    {
        cgh.require(a_tempScaleGroups);
    }

    return [=](cl::sycl::id<1> itemIdx) {
        const float3 x    = a_x[itemIdx];
        const float3 v    = a_v[itemIdx];
        const float3 f    = a_f[itemIdx];
        const float  im   = a_inverseMasses[itemIdx];
        const float  imdt = im * dt;

        // Swapping places for xp and x so that the x will contain the updated coordinates and xp -
        // the coordinates before update. This should be taken into account when (if) constraints
        // are applied after the update: x and xp have to be passed to constraints in the 'wrong'
        // order. See Issue #3727
        a_xp[itemIdx] = x;

        const float lambda = [=]() {
            if constexpr (numTempScaleValues == NumTempScaleValues::None)
            {
                return 1.0F;
            }
            else if constexpr (numTempScaleValues == NumTempScaleValues::Single)
            {
                return a_lambdas[0];
            }
            else if constexpr (numTempScaleValues == NumTempScaleValues::Multiple)
            {
                const int tempScaleGroup = a_tempScaleGroups[itemIdx];
                return a_lambdas[tempScaleGroup];
            }
        }();

        const float3 prVelocityDelta = [=]() {
            if constexpr (velocityScaling == VelocityScalingType::Diagonal)
            {
                return float3{ prVelocityScalingMatrixDiagonal[0] * v[0],
                               prVelocityScalingMatrixDiagonal[1] * v[1],
                               prVelocityScalingMatrixDiagonal[2] * v[2] };
            }
            else if constexpr (velocityScaling == VelocityScalingType::None)
            {
                return float3{ 0, 0, 0 };
            }
        }();

        const float3 v_new = v * lambda - prVelocityDelta + f * imdt;
        a_v[itemIdx]       = v_new;
        a_x[itemIdx]       = x + v_new * dt;
    };
}

// SYCL 1.2.1 requires providing a unique type for a kernel. Should not be needed for SYCL2020.
template<NumTempScaleValues numTempScaleValues, VelocityScalingType velocityScaling>
class LeapFrogKernelName;

template<NumTempScaleValues numTempScaleValues, VelocityScalingType velocityScaling, class... Args>
static cl::sycl::event launchLeapFrogKernel(const DeviceStream& deviceStream, int numAtoms, Args&&... args)
{
    // Should not be needed for SYCL2020.
    using kernelNameType = LeapFrogKernelName<numTempScaleValues, velocityScaling>;

    const cl::sycl::range<1> rangeAllAtoms(numAtoms);
    cl::sycl::queue          q = deviceStream.stream();

    cl::sycl::event e = q.submit([&](cl::sycl::handler& cgh) {
        auto kernel =
                leapFrogKernel<numTempScaleValues, velocityScaling>(cgh, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(rangeAllAtoms, kernel);
    });

    return e;
}

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
static inline cl::sycl::event launchLeapFrogKernel(NumTempScaleValues  tempScalingType,
                                                   VelocityScalingType prVelocityScalingType,
                                                   Args&&... args)
{
    GMX_ASSERT(prVelocityScalingType == VelocityScalingType::None
                       || prVelocityScalingType == VelocityScalingType::Diagonal,
               "Only isotropic Parrinello-Rahman pressure coupling is supported.");

    return dispatchTemplatedFunction(
            [&](auto tempScalingType_, auto prScalingType_) {
                return launchLeapFrogKernel<tempScalingType_, prScalingType_>(std::forward<Args>(args)...);
            },
            tempScalingType, prVelocityScalingType);
}

void LeapFrogGpu::integrate(DeviceBuffer<float3>              d_x,
                            DeviceBuffer<float3>              d_xp,
                            DeviceBuffer<float3>              d_v,
                            DeviceBuffer<float3>              d_f,
                            const real                        dt,
                            const bool                        doTemperatureScaling,
                            gmx::ArrayRef<const t_grp_tcstat> tcstat,
                            const bool                        doParrinelloRahman,
                            const float                       dtPressureCouple,
                            const matrix                      prVelocityScalingMatrix)
{
    if (doTemperatureScaling)
    {
        GMX_ASSERT(checkDeviceBuffer(d_lambdas_, numTempScaleValues_),
                   "Number of temperature scaling factors changed since it was set for the "
                   "last time.");
        { // Explicitly limiting the scope of host accessor. Not strictly necessary here.
            auto ha_lambdas_ = d_lambdas_.buffer_->get_access<mode::discard_write>();
            for (int i = 0; i < numTempScaleValues_; i++)
            {
                ha_lambdas_[i] = tcstat[i].lambda;
            }
        }
    }
    NumTempScaleValues tempVelocityScalingType =
            getTempScalingType(doTemperatureScaling, numTempScaleValues_);

    VelocityScalingType prVelocityScalingType = VelocityScalingType::None;
    if (doParrinelloRahman)
    {
        prVelocityScalingType = VelocityScalingType::Diagonal;
        GMX_ASSERT(prVelocityScalingMatrix[YY][XX] == 0 && prVelocityScalingMatrix[ZZ][XX] == 0
                           && prVelocityScalingMatrix[ZZ][YY] == 0 && prVelocityScalingMatrix[XX][YY] == 0
                           && prVelocityScalingMatrix[XX][ZZ] == 0 && prVelocityScalingMatrix[YY][ZZ] == 0,
                   "Fully anisotropic Parrinello-Rahman pressure coupling is not yet supported "
                   "in GPU version of Leap-Frog integrator.");
        prVelocityScalingMatrixDiagonal_ =
                dtPressureCouple
                * float3{ prVelocityScalingMatrix[XX][XX], prVelocityScalingMatrix[YY][YY],
                          prVelocityScalingMatrix[ZZ][ZZ] };
    }

    launchLeapFrogKernel(tempVelocityScalingType, prVelocityScalingType, deviceStream_, numAtoms_,
                         d_x, d_xp, d_v, d_f, d_inverseMasses_, dt, d_lambdas_, d_tempScaleGroups_,
                         prVelocityScalingMatrixDiagonal_);
}

LeapFrogGpu::LeapFrogGpu(const DeviceContext& deviceContext, const DeviceStream& deviceStream) :
    deviceContext_(deviceContext),
    deviceStream_(deviceStream),
    numAtoms_(0)
{
}

LeapFrogGpu::~LeapFrogGpu()
{
    freeDeviceBuffer(&d_inverseMasses_);
}

void LeapFrogGpu::set(const int             numAtoms,
                      const real*           inverseMasses,
                      const int             numTempScaleValues,
                      const unsigned short* tempScaleGroups)
{
    numAtoms_           = numAtoms;
    numTempScaleValues_ = numTempScaleValues;

    reallocateDeviceBuffer(&d_inverseMasses_, numAtoms_, &numInverseMasses_,
                           &numInverseMassesAlloc_, deviceContext_);
    copyToDeviceBuffer(&d_inverseMasses_, inverseMasses, 0, numAtoms_, deviceStream_,
                       GpuApiCallBehavior::Sync, nullptr);

    // Temperature scale group map only used if there are more then one group
    if (numTempScaleValues_ > 1)
    {
        reallocateDeviceBuffer(&d_tempScaleGroups_, numAtoms_, &numTempScaleGroups_,
                               &numTempScaleGroupsAlloc_, deviceContext_);
        copyToDeviceBuffer(&d_tempScaleGroups_, tempScaleGroups, 0, numAtoms_, deviceStream_,
                           GpuApiCallBehavior::Sync, nullptr);
    }

    // If the temperature coupling is enabled, we need to make space for scaling factors
    if (numTempScaleValues_ > 0)
    {
        reallocateDeviceBuffer(&d_lambdas_, numTempScaleValues_, &numLambdas_, &numLambdasAlloc_,
                               deviceContext_);
    }
}

} // namespace gmx
