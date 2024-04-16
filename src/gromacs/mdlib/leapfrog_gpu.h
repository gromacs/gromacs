/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 *
 * \brief Declarations for GPU implementation of Leap-Frog.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_LEAPFROG_GPU_H
#define GMX_MDLIB_LEAPFROG_GPU_H

#include "config.h"

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gputraits.cuh"
#endif
#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gputraits_sycl.h"
#endif

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/matrix.h"
#include "gromacs/utility/arrayref.h"

class DeviceContext;
class DeviceStream;
struct t_grp_tcstat;

namespace gmx
{


/*! \brief Sets the number of different temperature coupling values
 *
 *  This is needed to template the kernel
 *  \todo Unify with similar enum in CPU update module
 */
enum class NumTempScaleValues
{
    None     = 0, //!< No temperature coupling
    Single   = 1, //!< Single T-scaling value (one group)
    Multiple = 2, //!< Multiple T-scaling values, need to use T-group indices
    Count    = 3  //!< Number of valid values
};

// Avoid check-source warnings about the duplicate class enum
#ifndef DOXYGEN
/*! \brief Describes the properties of the Parrinello-Rahman pressure
 * scaling matrix
 *
 * This is needed to template the kernel
 */
enum class ParrinelloRahmanVelocityScaling
{
    No,          //!< Do not apply velocity scaling (not a PR-coupling run or step)
    Diagonal,    //!< Apply velocity scaling using a diagonal matrix
    Anisotropic, //!< Apply velocity scaling using a matrix with off-diagonal elements
    Count        //!< Number of valid values
};
#endif

class LeapFrogGpu
{

public:
    /*! \brief Constructor.
     *
     * \param[in] deviceContext       Device context.
     * \param[in] deviceStream        Device stream to use.
     * \param[in] numTempScaleValues  Number of temperature scale groups.
     */
    LeapFrogGpu(const DeviceContext& deviceContext, const DeviceStream& deviceStream, int numTempScaleValues);
    ~LeapFrogGpu();

    /*! \brief Integrate
     *
     * Integrates the equation of motion using Leap-Frog algorithm.
     * Updates coordinates and velocities on the GPU. The current coordinates are saved for constraints.
     *
     * \param[in,out] d_x                      Coordinates to update
     * \param[out]    d_x0                     Place to save the values of initial coordinates coordinates to.
     * \param[in,out] d_v                      Velocities (will be updated).
     * \param[in]     d_f                      Forces.
     * \param[in]     dt                       Timestep.
     * \param[in]     doTemperatureScaling     If velocities should be scaled for temperature coupling.
     * \param[in]     tcstat                   Temperature coupling data.
     * \param[in]     doParrinelloRahman       If current step is a Parrinello-Rahman pressure coupling step.
     * \param[in]     dtPressureCouple         Period between pressure coupling steps
     * \param[in]     prVelocityScalingMatrix  Parrinello-Rahman velocity scaling matrix
     */
    void integrate(DeviceBuffer<Float3>              d_x,
                   DeviceBuffer<Float3>              d_x0,
                   DeviceBuffer<Float3>              d_v,
                   DeviceBuffer<Float3>              d_f,
                   float                             dt,
                   bool                              doTemperatureScaling,
                   gmx::ArrayRef<const t_grp_tcstat> tcstat,
                   bool                              doParrinelloRahman,
                   float                             dtPressureCouple,
                   const Matrix3x3&                  prVelocityScalingMatrix);

    /*! \brief Set the integrator
     *
     * Allocates memory for inverse masses, and, if needed for temperature scaling factor(s)
     * and temperature coupling groups. Copies inverse masses and temperature coupling groups
     * to the GPU.
     *
     * \param[in] numAtoms        Number of atoms in the system.
     * \param[in] inverseMasses   Inverse masses of atoms.
     * \param[in] tempScaleGroups Maps the atom index to temperature scale value.
     */
    void set(int numAtoms, ArrayRef<const real> inverseMasses, ArrayRef<const unsigned short> tempScaleGroups);

    /*! \brief Class with hardware-specific interfaces and implementations.*/
    class Impl;

private:
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! GPU stream
    const DeviceStream& deviceStream_;

    //! Number of atoms
    int numAtoms_;

    //! 1/mass for all atoms (GPU)
    DeviceBuffer<float> d_inverseMasses_;
    //! Current size of the reciprocal masses array
    int numInverseMasses_ = -1;
    //! Maximum size of the reciprocal masses array
    int numInverseMassesAlloc_ = -1;

    //! Number of temperature coupling groups (zero = no coupling)
    int numTempScaleValues_ = 0;
    /*! \brief Array with temperature scaling factors.
     * This is temporary solution to remap data from t_grp_tcstat into plain array.
     * \todo Replace with better solution.
     */
    gmx::HostVector<float> h_lambdas_;
    //! Device-side temperature scaling factors
    DeviceBuffer<float> d_lambdas_;
    //! Current size of the array with temperature scaling factors (lambdas)
    int numLambdas_ = -1;
    //! Maximum size of the array with temperature scaling factors (lambdas)
    int numLambdasAlloc_ = -1;


    //! Array that maps atom index onto the temperature scaling group to get scaling parameter
    DeviceBuffer<unsigned short> d_tempScaleGroups_;
    //! Current size of the temperature coupling groups array
    int numTempScaleGroups_ = -1;
    //! Maximum size of the temperature coupling groups array
    int numTempScaleGroupsAlloc_ = -1;

    //! Vector with diagonal elements of the Parrinello-Rahman pressure coupling velocity rescale factors
    Float3 prVelocityScalingMatrixDiagonal_;
};

} // namespace gmx

#endif
