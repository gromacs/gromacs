/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief Declarations for CUDA implementation of Leap-Frog.
 *
 * \todo Reconsider naming towards using "gpu" suffix instead of "cuda".
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_LEAPFROG_CUDA_CUH
#define GMX_MDLIB_LEAPFROG_CUDA_CUH

#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class LeapFrogCuda
{

public:
    /*! \brief Constructor.
     *
     * \param[in] commandStream  Device command stream to use.
     */
    LeapFrogCuda(CommandStream commandStream);
    ~LeapFrogCuda();

    /*! \brief
     * Update PBC data.
     *
     * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \param[in] pbc The PBC data in t_pbc format.
     */
    void setPbc(const t_pbc* pbc);

    /*! \brief Integrate
     *
     * Integrates the equation of motion using Leap-Frog algorithm.
     * Updates coordinates and velocities on the GPU. The current coordinates are saved for constraints.
     *
     * \param[in,out] d_x                      Coordinates to update
     * \param[out]    d_xp                     Place to save the values of initial coordinates coordinates to.
     * \param[in,out] d_v                      Velocities (will be updated).
     * \param[in]     d_f                      Forces.
     * \param[in]     dt                       Timestep.
     * \param[in]     doTemperatureScaling     If velocities should be scaled for temperature coupling.
     * \param[in]     tcstat                   Temperature coupling data.
     * \param[in]     doParrinelloRahman       If current step is a Parrinello-Rahman pressure coupling step.
     * \param[in]     dtPressureCouple         Period between pressure coupling steps
     * \param[in]     prVelocityScalingMatrix  Parrinello-Rahman velocity scaling matrix
     */
    void integrate(const float3*                     d_x,
                   float3*                           d_xp,
                   float3*                           d_v,
                   const float3*                     d_f,
                   const real                        dt,
                   const bool                        doTemperatureScaling,
                   gmx::ArrayRef<const t_grp_tcstat> tcstat,
                   const bool                        doParrinelloRahman,
                   const float                       dtPressureCouple,
                   const matrix                      prVelocityScalingMatrix);

    /*! \brief Set the integrator
     *
     * Allocates memory for inverse masses, and, if needed for temperature scaling factor(s)
     * and temperature coupling groups. Copies inverse masses and temperature coupling groups
     * to the GPU.
     *
     * \param[in] md                  MD atoms, from which inverse masses are taken.
     * \param[in] numTempScaleValues  Number of temperature scale groups.
     * \param[in] tempScaleGroups     Maps the atom index to temperature scale value.
     */
    void set(const t_mdatoms& md, int numTempScaleValues, const unsigned short* tempScaleGroups);

    /*! \brief Class with hardware-specific interfaces and implementations.*/
    class Impl;

private:
    //! CUDA stream
    CommandStream commandStream_;
    //! CUDA kernel launch config
    KernelLaunchConfig kernelLaunchConfig_;
    //! Periodic boundary data
    PbcAiuc pbcAiuc_;
    //! Number of atoms
    int numAtoms_;

    //! 1/mass for all atoms (GPU)
    real* d_inverseMasses_;
    //! Current size of the reciprocal masses array
    int numInverseMasses_ = -1;
    //! Maximum size of the reciprocal masses array
    int numInverseMassesAlloc_ = -1;

    //! Number of temperature coupling groups (zero = no coupling)
    int numTempScaleValues_ = 0;
    /*! \brief Array with temperature scaling factors.
     * This is temporary solution to remap data from t_grp_tcstat into plain array
     * \todo Replace with better solution.
     */
    gmx::HostVector<float> h_lambdas_;
    //! Device-side temperature scaling factors
    float* d_lambdas_;
    //! Current size of the array with temperature scaling factors (lambdas)
    int numLambdas_ = -1;
    //! Maximum size of the array with temperature scaling factors (lambdas)
    int numLambdasAlloc_ = -1;


    //! Array that maps atom index onto the temperature scaling group to get scaling parameter
    unsigned short* d_tempScaleGroups_;
    //! Current size of the temperature coupling groups array
    int numTempScaleGroups_ = -1;
    //! Maximum size of the temperature coupling groups array
    int numTempScaleGroupsAlloc_ = -1;

    //! Vector with diagonal elements of the Parrinello-Rahman pressure coupling velocity rescale factors
    float3 prVelocityScalingMatrixDiagonal_;
};

} // namespace gmx

#endif
