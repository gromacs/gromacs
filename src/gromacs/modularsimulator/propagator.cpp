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
/*! \libinternal
 * \brief Defines the propagator element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "propagator.h"

#include "gromacs/utility.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/timing/wallcycle.h"

#include "statepropagatordata.h"

namespace gmx
{
//! Update velocities
static void inline
updateVelocities(int                       a,
                 real                      dt,
                 const rvec * gmx_restrict invMassPerDim,
                 rvec       * gmx_restrict v,
                 const rvec * gmx_restrict f)
{
    for (int d = 0; d < DIM; d++)
    {
        v[a][d] += f[a][d]*invMassPerDim[a][d]*dt;
    }
}

//! Update velocities - 2nd version. The parentheses allows to exactly replicate do_md results under VV
static void inline
updateVelocities2(int                       a,
                  real                      dt,
                  const rvec * gmx_restrict invMassPerDim,
                  rvec       * gmx_restrict v,
                  const rvec * gmx_restrict f)
{
    for (int d = 0; d < DIM; d++)
    {
        v[a][d] += f[a][d]*(invMassPerDim[a][d]*dt);
    }
}

//! Update positions
static void inline
updatePositions(int                       a,
                real                      dt,
                const rvec * gmx_restrict x,
                rvec       * gmx_restrict xprime,
                const rvec * gmx_restrict v)
{
    for (int d = 0; d < DIM; d++)
    {
        xprime[a][d]  = x[a][d] + v[a][d]*dt;
    }
}

//! Propagation (position only)
template <>
void Propagator<IntegrationStep::PositionsOnly>::run()
{
    wallcycle_start(wcycle_, ewcUPDATE);

    auto xp = as_rvec_array(statePropagatorData_->positionsView().paddedArrayRef().data());
    auto x  = as_rvec_array(statePropagatorData_->constPreviousPositionsView().paddedArrayRef().data());
    auto v  = as_rvec_array(statePropagatorData_->constVelocitiesView().paddedArrayRef().data());

    int  nth    = gmx_omp_nthreads_get(emntUpdate);
    int  homenr = mdAtoms_->mdatoms()->homenr;

    #pragma omp parallel for num_threads(nth) schedule(static) default(none) \
    shared(nth, homenr, x, xp, v)
    for (int th = 0; th < nth; th++)
    {
        try
        {

            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                updatePositions(a, timestep_, x, xp, v);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, ewcUPDATE);
}

//! Propagation (velocity only)
template <>
void Propagator<IntegrationStep::VelocitiesOnly>::run()
{
    wallcycle_start(wcycle_, ewcUPDATE);

    auto v             = as_rvec_array(statePropagatorData_->velocitiesView().paddedArrayRef().data());
    auto f             = as_rvec_array(statePropagatorData_->constForcesView().paddedArrayRef().data());
    auto invMassPerDim = mdAtoms_->mdatoms()->invMassPerDim;

    int  nth           = gmx_omp_nthreads_get(emntUpdate);
    int  homenr        = mdAtoms_->mdatoms()->homenr;

    #pragma omp parallel for num_threads(nth) schedule(static) default(none) \
    shared(nth, homenr, v, f, invMassPerDim)
    for (int th = 0; th < nth; th++)
    {
        try
        {

            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                updateVelocities(a, timestep_, invMassPerDim, v, f);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, ewcUPDATE);
}

//! Propagation (leapfrog case - position and velocity)
template <>
void Propagator<IntegrationStep::LeapFrog>::run()
{
    wallcycle_start(wcycle_, ewcUPDATE);

    auto xp            = as_rvec_array(statePropagatorData_->positionsView().paddedArrayRef().data());
    auto x             = as_rvec_array(statePropagatorData_->constPreviousPositionsView().paddedArrayRef().data());
    auto v             = as_rvec_array(statePropagatorData_->velocitiesView().paddedArrayRef().data());
    auto f             = as_rvec_array(statePropagatorData_->constForcesView().paddedArrayRef().data());
    auto invMassPerDim = mdAtoms_->mdatoms()->invMassPerDim;

    int  nth           = gmx_omp_nthreads_get(emntUpdate);
    int  homenr        = mdAtoms_->mdatoms()->homenr;

    #pragma omp parallel for num_threads(nth) schedule(static) default(none) \
    shared(nth, homenr, x, xp, v, f, invMassPerDim)
    for (int th = 0; th < nth; th++)
    {
        try
        {

            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                updateVelocities(a, timestep_, invMassPerDim, v, f);
                updatePositions(a, timestep_, x, xp, v);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, ewcUPDATE);
}

//! Propagation (velocity verlet stage 2 - velocity and position)
template <>
void Propagator<IntegrationStep::VelocityVerletPositionsAndVelocities>::run()
{
    wallcycle_start(wcycle_, ewcUPDATE);

    auto xp            = as_rvec_array(statePropagatorData_->positionsView().paddedArrayRef().data());
    auto x             = as_rvec_array(statePropagatorData_->constPreviousPositionsView().paddedArrayRef().data());
    auto v             = as_rvec_array(statePropagatorData_->velocitiesView().paddedArrayRef().data());
    auto f             = as_rvec_array(statePropagatorData_->constForcesView().paddedArrayRef().data());
    auto invMassPerDim = mdAtoms_->mdatoms()->invMassPerDim;

    int  nth           = gmx_omp_nthreads_get(emntUpdate);
    int  homenr        = mdAtoms_->mdatoms()->homenr;

    #pragma omp parallel for num_threads(nth) schedule(static) default(none) \
    shared(nth, homenr, x, xp, v, f, invMassPerDim)
    for (int th = 0; th < nth; th++)
    {
        try
        {

            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                updateVelocities2(a, timestep_*0.5, invMassPerDim, v, f);
                updatePositions(a, timestep_, x, xp, v);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, ewcUPDATE);
}

template <IntegrationStep algorithm>
Propagator<algorithm>::Propagator(
        double               timestep,
        StatePropagatorData *statePropagatorData,
        const MDAtoms       *mdAtoms,
        gmx_wallcycle       *wcycle) :
    timestep_(timestep),
    statePropagatorData_(statePropagatorData),
    mdAtoms_(mdAtoms),
    wcycle_(wcycle)
{}

template <IntegrationStep algorithm>
void Propagator<algorithm>::scheduleTask(
        Step gmx_unused step, Time gmx_unused time,
        const RegisterRunFunctionPtr &registerRunFunction)
{
    (*registerRunFunction)(
            std::make_unique<SimulatorRunFunction>(
                    [this](){run(); }));
}

//! Explicit template initialization
//! @{
template class Propagator<IntegrationStep::PositionsOnly>;
template class Propagator<IntegrationStep::VelocitiesOnly>;
template class Propagator<IntegrationStep::LeapFrog>;
template class Propagator<IntegrationStep::VelocityVerletPositionsAndVelocities>;
//! @}

}  // namespace gmx
