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
/*! \internal \file
 * \brief Defines the propagator element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "propagator.h"

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/timing/wallcycle.h"

#include "modularsimulator.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

namespace gmx
{
namespace
{
// Names of integration steps, only used locally for error messages
constexpr EnumerationArray<IntegrationStage, const char*> integrationStepNames = {
    "IntegrationStage::PositionsOnly",   "IntegrationStage::VelocitiesOnly",
    "IntegrationStage::LeapFrog",        "IntegrationStage::VelocityVerletPositionsAndVelocities",
    "IntegrationStage::ScaleVelocities", "IntegrationStage::ScalePositions"
};
} // namespace

/*! \brief Update velocities
 *
 * To maximize the ability of the compiler to optimize, all the arrays
 * of RVec should be annotated with gmx_restrict, so the compiler knows
 * there is no aliasing, and for the same reason we do not use
 * ArrayRef<RVec> for them. */
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues>
static void inline updateVelocities(int                      a,
                                    real                     dt,
                                    real                     lambdaStart,
                                    real                     lambdaEnd,
                                    const RVec* gmx_restrict invMassPerDim,
                                    RVec* gmx_restrict       v,
                                    const RVec* gmx_restrict f,
                                    const RVec&              diagPR,
                                    const Matrix3x3&         matrixPR)
{
    RVec parrinelloRahmanScaledVelocity;
    if (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Anisotropic)
    {
        parrinelloRahmanScaledVelocity = matrixPR * v[a];
    }
    for (int d = 0; d < DIM; d++)
    {
        // TODO: Extract this into policy classes
        if (numStartVelocityScalingValues != NumVelocityScalingValues::None
            && parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::No)
        {
            v[a][d] *= lambdaStart;
        }
        if (numStartVelocityScalingValues != NumVelocityScalingValues::None
            && parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal)
        {
            v[a][d] *= (lambdaStart - diagPR[d]);
        }
        if (numStartVelocityScalingValues != NumVelocityScalingValues::None
            && parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Anisotropic)
        {
            v[a][d] = lambdaStart * v[a][d] - parrinelloRahmanScaledVelocity[d];
        }
        if (numStartVelocityScalingValues == NumVelocityScalingValues::None
            && parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal)
        {
            v[a][d] *= (1 - diagPR[d]);
        }
        if (numStartVelocityScalingValues == NumVelocityScalingValues::None
            && parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Anisotropic)
        {
            v[a][d] -= parrinelloRahmanScaledVelocity[d];
        }
        v[a][d] += f[a][d] * invMassPerDim[a][d] * dt;
        if (numEndVelocityScalingValues != NumVelocityScalingValues::None)
        {
            v[a][d] *= lambdaEnd;
        }
    }
}

/*! \brief Update positions
 *
 * To maximize the ability of the compiler to optimize, all the arrays
 * of RVec should be annotated with gmx_restrict, so the compiler knows
 * there is no aliasing, and for the same reason we do not use
 * ArrayRef<RVec> for them. */
static void inline updatePositions(int                      a,
                                   real                     dt,
                                   const RVec* gmx_restrict x,
                                   RVec* gmx_restrict       xprime,
                                   const RVec* gmx_restrict v)
{
    for (int d = 0; d < DIM; d++)
    {
        xprime[a][d] = x[a][d] + v[a][d] * dt;
    }
}

/*! \brief Scale velocities
 *
 * To maximize the ability of the compiler to optimize, all the arrays
 * of RVec should be annotated with gmx_restrict, so the compiler knows
 * there is no aliasing, and for the same reason we do not use
 * ArrayRef<RVec> for them. */
template<NumVelocityScalingValues numStartVelocityScalingValues>
static void inline scaleVelocities(int a, real lambda, RVec* gmx_restrict v)
{
    if (numStartVelocityScalingValues != NumVelocityScalingValues::None)
    {
        for (int d = 0; d < DIM; d++)
        {
            v[a][d] *= lambda;
        }
    }
}

/*! \brief Scale positions
 *
 * To maximize the ability of the compiler to optimize, all the arrays
 * of RVec should be annotated with gmx_restrict, so the compiler knows
 * there is no aliasing, and for the same reason we do not use
 * ArrayRef<RVec> for them. */
template<NumPositionScalingValues numPositionScalingValues>
static void inline scalePositions(int a, real lambda, RVec* gmx_restrict x)
{
    if (numPositionScalingValues != NumPositionScalingValues::None)
    {
        for (int d = 0; d < DIM; d++)
        {
            x[a][d] *= lambda;
        }
    }
}

//! Is the PR matrix diagonal?
template<ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling>
static inline bool canTreatPRScalingMatrixAsDiagonal(const Matrix3x3& matrixPR)
{
    if (parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::Anisotropic)
    {
        return false;
    }
    else
    {
        return (matrixPR(YY, XX) == 0 && matrixPR(ZZ, XX) == 0 && matrixPR(ZZ, YY) == 0);
    }
}

//! Propagation (position only)
template<>
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues,
         NumPositionScalingValues        numPositionScalingValues>
void Propagator<IntegrationStage::PositionsOnly>::run()
{
    wallcycle_start(wcycle_, WallCycleCounter::Update);

    RVec*       xp = statePropagatorData_->positionsView().paddedArrayRef().data();
    const RVec* x  = statePropagatorData_->constPositionsView().paddedArrayRef().data();
    const RVec* v  = statePropagatorData_->constVelocitiesView().paddedArrayRef().data();

    int nth    = gmx_omp_nthreads_get(ModuleMultiThread::Update);
    int homenr = mdAtoms_->mdatoms()->homenr;

#pragma omp parallel for num_threads(nth) schedule(static) default(none) shared(nth, homenr, x, xp, v)
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
    wallcycle_stop(wcycle_, WallCycleCounter::Update);
}

//! Propagation (scale position only)
template<>
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues,
         NumPositionScalingValues        numPositionScalingValues>
void Propagator<IntegrationStage::ScalePositions>::run()
{
    wallcycle_start(wcycle_, WallCycleCounter::Update);

    RVec* x = statePropagatorData_->positionsView().paddedArrayRef().data();

    const real lambda =
            (numPositionScalingValues == NumPositionScalingValues::Single) ? positionScaling_[0] : 1.0;

    int nth    = gmx_omp_nthreads_get(ModuleMultiThread::Update);
    int homenr = mdAtoms_->mdatoms()->homenr;

#pragma omp parallel for num_threads(nth) schedule(static) default(none) shared(nth, homenr, x) \
        firstprivate(lambda)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                scalePositions<numPositionScalingValues>(
                        a,
                        (numPositionScalingValues == NumPositionScalingValues::Multiple)
                                ? positionScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                : lambda,
                        x);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, WallCycleCounter::Update);
}

//! Propagation (velocity only)
template<>
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues,
         NumPositionScalingValues        numPositionScalingValues>
void Propagator<IntegrationStage::VelocitiesOnly>::run()
{
    wallcycle_start(wcycle_, WallCycleCounter::Update);

    RVec*                      v = statePropagatorData_->velocitiesView().paddedArrayRef().data();
    const RVec*                f = statePropagatorData_->constForcesView().force().data();
    const ArrayRef<const RVec> invMassPerDim = mdAtoms_->mdatoms()->invMassPerDim;

    const real lambdaStart = (numStartVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? startVelocityScaling_[0]
                                     : 1.0;
    const real lambdaEnd   = (numEndVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? endVelocityScaling_[0]
                                     : 1.0;

    const bool treatPRScalingMatrixAsDiagonal =
            canTreatPRScalingMatrixAsDiagonal<parrinelloRahmanVelocityScaling>(matrixPR_);
    const RVec diagonalOfPRScalingMatrix = treatPRScalingMatrixAsDiagonal ? diagonal(matrixPR_) : RVec{};

    const int nth    = gmx_omp_nthreads_get(ModuleMultiThread::Update);
    const int homenr = mdAtoms_->mdatoms()->homenr;

#pragma omp parallel for num_threads(nth) schedule(static) default(none) shared(v, f, invMassPerDim) \
        shared(nth, homenr, lambdaStart, lambdaEnd, treatPRScalingMatrixAsDiagonal, diagonalOfPRScalingMatrix)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                if (treatPRScalingMatrixAsDiagonal)
                {
                    updateVelocities<numStartVelocityScalingValues, ParrinelloRahmanVelocityScaling::Diagonal, numEndVelocityScalingValues>(
                            a,
                            timestep_,
                            numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaStart,
                            numEndVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? endVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaEnd,
                            invMassPerDim.data(),
                            v,
                            f,
                            diagonalOfPRScalingMatrix,
                            matrixPR_);
                }
                else
                {
                    updateVelocities<numStartVelocityScalingValues, parrinelloRahmanVelocityScaling, numEndVelocityScalingValues>(
                            a,
                            timestep_,
                            numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaStart,
                            numEndVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? endVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaEnd,
                            invMassPerDim.data(),
                            v,
                            f,
                            diagonalOfPRScalingMatrix,
                            matrixPR_);
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, WallCycleCounter::Update);
}

//! Propagation (leapfrog case - position and velocity)
template<>
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues,
         NumPositionScalingValues        numPositionScalingValues>
void Propagator<IntegrationStage::LeapFrog>::run()
{
    wallcycle_start(wcycle_, WallCycleCounter::Update);

    RVec*       xp = statePropagatorData_->positionsView().paddedArrayRef().data();
    const RVec* x  = statePropagatorData_->constPositionsView().paddedArrayRef().data();
    RVec*       v  = statePropagatorData_->velocitiesView().paddedArrayRef().data();
    const RVec* f  = statePropagatorData_->constForcesView().force().data();
    const ArrayRef<const RVec> invMassPerDim = mdAtoms_->mdatoms()->invMassPerDim;

    const real lambdaStart = (numStartVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? startVelocityScaling_[0]
                                     : 1.0;
    const real lambdaEnd   = (numEndVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? endVelocityScaling_[0]
                                     : 1.0;

    const bool treatPRScalingMatrixAsDiagonal =
            canTreatPRScalingMatrixAsDiagonal<parrinelloRahmanVelocityScaling>(matrixPR_);
    const RVec diagonalOfPRScalingMatrix = treatPRScalingMatrixAsDiagonal ? diagonal(matrixPR_) : RVec{};

    const int nth    = gmx_omp_nthreads_get(ModuleMultiThread::Update);
    const int homenr = mdAtoms_->mdatoms()->homenr;

#pragma omp parallel for num_threads(nth) schedule(static) default(none) \
        shared(x, xp, v, f, invMassPerDim)                               \
        firstprivate(nth, homenr, lambdaStart, lambdaEnd, treatPRScalingMatrixAsDiagonal, diagonalOfPRScalingMatrix)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                if (treatPRScalingMatrixAsDiagonal)
                {
                    updateVelocities<numStartVelocityScalingValues, ParrinelloRahmanVelocityScaling::Diagonal, numEndVelocityScalingValues>(
                            a,
                            timestep_,
                            numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaStart,
                            numEndVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? endVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaEnd,
                            invMassPerDim.data(),
                            v,
                            f,
                            diagonalOfPRScalingMatrix,
                            matrixPR_);
                }
                else
                {
                    updateVelocities<numStartVelocityScalingValues, parrinelloRahmanVelocityScaling, numEndVelocityScalingValues>(
                            a,
                            timestep_,
                            numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaStart,
                            numEndVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? endVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaEnd,
                            invMassPerDim.data(),
                            v,
                            f,
                            diagonalOfPRScalingMatrix,
                            matrixPR_);
                }
                updatePositions(a, timestep_, x, xp, v);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, WallCycleCounter::Update);
}

//! Propagation (velocity verlet stage 2 - velocity and position)
template<>
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues,
         NumPositionScalingValues        numPositionScalingValues>
void Propagator<IntegrationStage::VelocityVerletPositionsAndVelocities>::run()
{
    wallcycle_start(wcycle_, WallCycleCounter::Update);

    RVec*       xp = statePropagatorData_->positionsView().paddedArrayRef().data();
    const RVec* x  = statePropagatorData_->constPositionsView().paddedArrayRef().data();
    RVec*       v  = statePropagatorData_->velocitiesView().paddedArrayRef().data();
    const RVec* f  = statePropagatorData_->constForcesView().force().data();
    const ArrayRef<const RVec> invMassPerDim = mdAtoms_->mdatoms()->invMassPerDim;

    const real lambdaStart = (numStartVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? startVelocityScaling_[0]
                                     : 1.0;
    const real lambdaEnd   = (numEndVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? endVelocityScaling_[0]
                                     : 1.0;

    const bool treatPRScalingMatrixAsDiagonal =
            canTreatPRScalingMatrixAsDiagonal<parrinelloRahmanVelocityScaling>(matrixPR_);
    const RVec diagonalOfPRScalingMatrix = treatPRScalingMatrixAsDiagonal ? diagonal(matrixPR_) : RVec{};

    const int nth    = gmx_omp_nthreads_get(ModuleMultiThread::Update);
    const int homenr = mdAtoms_->mdatoms()->homenr;

#pragma omp parallel for num_threads(nth) schedule(static) default(none) \
        shared(x, xp, v, f, invMassPerDim)                               \
        firstprivate(nth, homenr, lambdaStart, lambdaEnd, treatPRScalingMatrixAsDiagonal, diagonalOfPRScalingMatrix)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                if (treatPRScalingMatrixAsDiagonal)
                {
                    updateVelocities<numStartVelocityScalingValues, ParrinelloRahmanVelocityScaling::Diagonal, numEndVelocityScalingValues>(
                            a,
                            0.5 * timestep_,
                            numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaStart,
                            numEndVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? endVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaEnd,
                            invMassPerDim.data(),
                            v,
                            f,
                            diagonalOfPRScalingMatrix,
                            matrixPR_);
                }
                else
                {
                    updateVelocities<numStartVelocityScalingValues, parrinelloRahmanVelocityScaling, numEndVelocityScalingValues>(
                            a,
                            0.5 * timestep_,
                            numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaStart,
                            numEndVelocityScalingValues == NumVelocityScalingValues::Multiple
                                    ? endVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                    : lambdaEnd,
                            invMassPerDim.data(),
                            v,
                            f,
                            diagonalOfPRScalingMatrix,
                            matrixPR_);
                }
                updatePositions(a, timestep_, x, xp, v);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, WallCycleCounter::Update);
}

//! Scaling (velocity scaling only)
template<>
template<NumVelocityScalingValues        numStartVelocityScalingValues,
         ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
         NumVelocityScalingValues        numEndVelocityScalingValues,
         NumPositionScalingValues        numPositionScalingValues>
void Propagator<IntegrationStage::ScaleVelocities>::run()
{
    if (numStartVelocityScalingValues == NumVelocityScalingValues::None)
    {
        return;
    }
    wallcycle_start(wcycle_, WallCycleCounter::Update);

    RVec* v = statePropagatorData_->velocitiesView().paddedArrayRef().data();

    const real lambdaStart = (numStartVelocityScalingValues == NumVelocityScalingValues::Single)
                                     ? startVelocityScaling_[0]
                                     : 1.0;

    const int nth    = gmx_omp_nthreads_get(ModuleMultiThread::Update);
    const int homenr = mdAtoms_->mdatoms()->homenr;

#pragma omp parallel for num_threads(nth) schedule(static) default(none) \
        shared(v, lambdaStart, nth, homenr)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th = 0;
            int end_th   = 0;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            for (int a = start_th; a < end_th; a++)
            {
                scaleVelocities<numStartVelocityScalingValues>(
                        a,
                        numStartVelocityScalingValues == NumVelocityScalingValues::Multiple
                                ? startVelocityScaling_[mdAtoms_->mdatoms()->cTC[a]]
                                : lambdaStart,
                        v);
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    wallcycle_stop(wcycle_, WallCycleCounter::Update);
}

template<IntegrationStage integrationStage>
Propagator<integrationStage>::Propagator(double               timestep,
                                         StatePropagatorData* statePropagatorData,
                                         const MDAtoms*       mdAtoms,
                                         gmx_wallcycle*       wcycle) :
    timestep_(timestep),
    statePropagatorData_(statePropagatorData),
    doSingleStartVelocityScaling_(false),
    doGroupStartVelocityScaling_(false),
    doSingleEndVelocityScaling_(false),
    doGroupEndVelocityScaling_(false),
    scalingStepVelocity_(-1),
    matrixPR_{ 0 },
    scalingStepPR_(-1),
    mdAtoms_(mdAtoms),
    wcycle_(wcycle)
{
}

template<IntegrationStage integrationStage>
void Propagator<integrationStage>::scheduleTask(Step                       step,
                                                Time gmx_unused            time,
                                                const RegisterRunFunction& registerRunFunction)
{
    const bool doSingleVScalingThisStep =
            (doSingleStartVelocityScaling_ && (step == scalingStepVelocity_));
    const bool doGroupVScalingThisStep = (doGroupStartVelocityScaling_ && (step == scalingStepVelocity_));

    if (integrationStage == IntegrationStage::ScaleVelocities)
    {
        // IntegrationStage::ScaleVelocities only needs to run if some kind of
        // velocity scaling is needed on the current step.
        if (!doSingleVScalingThisStep && !doGroupVScalingThisStep)
        {
            return;
        }
    }

    if (integrationStage == IntegrationStage::ScalePositions)
    {
        // IntegrationStage::ScalePositions only needs to run if
        // position scaling is needed on the current step.
        if (step != scalingStepPosition_)
        {
            return;
        }
        // Since IntegrationStage::ScalePositions is the only stage for which position scaling
        // is implemented we handle it here to avoid enlarging the decision tree below.
        if (doSinglePositionScaling_)
        {
            registerRunFunction(
                    [this]()
                    {
                        run<NumVelocityScalingValues::None,
                            ParrinelloRahmanVelocityScaling::No,
                            NumVelocityScalingValues::None,
                            NumPositionScalingValues::Single>();
                    });
        }
        else if (doGroupPositionScaling_)
        {
            registerRunFunction(
                    [this]()
                    {
                        run<NumVelocityScalingValues::None,
                            ParrinelloRahmanVelocityScaling::No,
                            NumVelocityScalingValues::None,
                            NumPositionScalingValues::Multiple>();
                    });
        }
    }

    const bool doParrinelloRahmanThisStep = (step == scalingStepPR_);

    if (doSingleVScalingThisStep)
    {
        if (doParrinelloRahmanThisStep)
        {
            if (doSingleEndVelocityScaling_)
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Single,
                                ParrinelloRahmanVelocityScaling::Anisotropic,
                                NumVelocityScalingValues::Single,
                                NumPositionScalingValues::None>();
                        });
            }
            else
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Single,
                                ParrinelloRahmanVelocityScaling::Anisotropic,
                                NumVelocityScalingValues::None,
                                NumPositionScalingValues::None>();
                        });
            }
        }
        else
        {
            if (doSingleEndVelocityScaling_)
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Single,
                                ParrinelloRahmanVelocityScaling::No,
                                NumVelocityScalingValues::Single,
                                NumPositionScalingValues::None>();
                        });
            }
            else
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Single,
                                ParrinelloRahmanVelocityScaling::No,
                                NumVelocityScalingValues::None,
                                NumPositionScalingValues::None>();
                        });
            }
        }
    }
    else if (doGroupVScalingThisStep)
    {
        if (doParrinelloRahmanThisStep)
        {
            if (doGroupEndVelocityScaling_)
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Multiple,
                                ParrinelloRahmanVelocityScaling::Anisotropic,
                                NumVelocityScalingValues::Multiple,
                                NumPositionScalingValues::None>();
                        });
            }
            else
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Multiple,
                                ParrinelloRahmanVelocityScaling::Anisotropic,
                                NumVelocityScalingValues::None,
                                NumPositionScalingValues::None>();
                        });
            }
        }
        else
        {
            if (doGroupEndVelocityScaling_)
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Multiple,
                                ParrinelloRahmanVelocityScaling::No,
                                NumVelocityScalingValues::Multiple,
                                NumPositionScalingValues::None>();
                        });
            }
            else
            {
                registerRunFunction(
                        [this]()
                        {
                            run<NumVelocityScalingValues::Multiple,
                                ParrinelloRahmanVelocityScaling::No,
                                NumVelocityScalingValues::None,
                                NumPositionScalingValues::None>();
                        });
            }
        }
    }
    else
    {
        if (doParrinelloRahmanThisStep)
        {
            registerRunFunction(
                    [this]()
                    {
                        run<NumVelocityScalingValues::None,
                            ParrinelloRahmanVelocityScaling::Anisotropic,
                            NumVelocityScalingValues::None,
                            NumPositionScalingValues::None>();
                    });
        }
        else
        {
            registerRunFunction(
                    [this]()
                    {
                        run<NumVelocityScalingValues::None,
                            ParrinelloRahmanVelocityScaling::No,
                            NumVelocityScalingValues::None,
                            NumPositionScalingValues::None>();
                    });
        }
    }
}

template<IntegrationStage integrationStage>
constexpr bool hasStartVelocityScaling()
{
    return (integrationStage == IntegrationStage::VelocitiesOnly
            || integrationStage == IntegrationStage::LeapFrog
            || integrationStage == IntegrationStage::VelocityVerletPositionsAndVelocities
            || integrationStage == IntegrationStage::ScaleVelocities);
}

template<IntegrationStage integrationStage>
constexpr bool hasEndVelocityScaling()
{
    return (hasStartVelocityScaling<integrationStage>()
            && integrationStage != IntegrationStage::ScaleVelocities);
}

template<IntegrationStage integrationStage>
constexpr bool hasPositionScaling()
{
    return (integrationStage == IntegrationStage::ScalePositions);
}

template<IntegrationStage integrationStage>
constexpr bool hasParrinelloRahmanScaling()
{
    return (integrationStage == IntegrationStage::VelocitiesOnly
            || integrationStage == IntegrationStage::LeapFrog
            || integrationStage == IntegrationStage::VelocityVerletPositionsAndVelocities);
}

template<IntegrationStage integrationStage>
void Propagator<integrationStage>::setNumVelocityScalingVariables(int numVelocityScalingVariables,
                                                                  ScaleVelocities scaleVelocities)
{
    GMX_RELEASE_ASSERT(
            hasStartVelocityScaling<integrationStage>() || hasEndVelocityScaling<integrationStage>(),
            formatString("Velocity scaling not implemented for %s", integrationStepNames[integrationStage])
                    .c_str());
    GMX_RELEASE_ASSERT(startVelocityScaling_.empty(),
                       "Number of velocity scaling variables cannot be changed once set.");

    const bool scaleEndVelocities = (scaleVelocities == ScaleVelocities::PreStepAndPostStep);
    startVelocityScaling_.resize(numVelocityScalingVariables, 1.);
    if (scaleEndVelocities)
    {
        endVelocityScaling_.resize(numVelocityScalingVariables, 1.);
    }
    doSingleStartVelocityScaling_ = numVelocityScalingVariables == 1;
    doGroupStartVelocityScaling_  = numVelocityScalingVariables > 1;
    doSingleEndVelocityScaling_   = doSingleStartVelocityScaling_ && scaleEndVelocities;
    doGroupEndVelocityScaling_    = doGroupStartVelocityScaling_ && scaleEndVelocities;
}

template<IntegrationStage integrationStage>
void Propagator<integrationStage>::setNumPositionScalingVariables(int numPositionScalingVariables)
{
    GMX_RELEASE_ASSERT(hasPositionScaling<integrationStage>(),
                       formatString("Position scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());
    GMX_RELEASE_ASSERT(positionScaling_.empty(),
                       "Number of position scaling variables cannot be changed once set.");
    positionScaling_.resize(numPositionScalingVariables, 1.);
    doSinglePositionScaling_ = (numPositionScalingVariables == 1);
    doGroupPositionScaling_  = (numPositionScalingVariables > 1);
}

template<IntegrationStage integrationStage>
ArrayRef<real> Propagator<integrationStage>::viewOnStartVelocityScaling()
{
    GMX_RELEASE_ASSERT(hasStartVelocityScaling<integrationStage>(),
                       formatString("Start velocity scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());
    GMX_RELEASE_ASSERT(!startVelocityScaling_.empty(),
                       "Number of velocity scaling variables not set.");

    return startVelocityScaling_;
}

template<IntegrationStage integrationStage>
ArrayRef<real> Propagator<integrationStage>::viewOnEndVelocityScaling()
{
    GMX_RELEASE_ASSERT(hasEndVelocityScaling<integrationStage>(),
                       formatString("End velocity scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());
    GMX_RELEASE_ASSERT(!endVelocityScaling_.empty(),
                       "Number of velocity scaling variables not set.");

    return endVelocityScaling_;
}

template<IntegrationStage integrationStage>
ArrayRef<real> Propagator<integrationStage>::viewOnPositionScaling()
{
    GMX_RELEASE_ASSERT(hasPositionScaling<integrationStage>(),
                       formatString("Position scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());
    GMX_RELEASE_ASSERT(!positionScaling_.empty(), "Number of position scaling variables not set.");

    return positionScaling_;
}

template<IntegrationStage integrationStage>
PropagatorCallback Propagator<integrationStage>::velocityScalingCallback()
{
    GMX_RELEASE_ASSERT(
            hasStartVelocityScaling<integrationStage>() || hasEndVelocityScaling<integrationStage>(),
            formatString("Velocity scaling not implemented for %s", integrationStepNames[integrationStage])
                    .c_str());

    return [this](Step step) { scalingStepVelocity_ = step; };
}

template<IntegrationStage integrationStage>
PropagatorCallback Propagator<integrationStage>::positionScalingCallback()
{
    GMX_RELEASE_ASSERT(hasPositionScaling<integrationStage>(),
                       formatString("Position scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());

    return [this](Step step) { scalingStepPosition_ = step; };
}

template<IntegrationStage integrationStage>
Matrix3x3* Propagator<integrationStage>::viewOnPRScalingMatrix()
{
    GMX_RELEASE_ASSERT(hasParrinelloRahmanScaling<integrationStage>(),
                       formatString("Parrinello-Rahman scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());
    return &matrixPR_;
}

template<IntegrationStage integrationStage>
PropagatorCallback Propagator<integrationStage>::prScalingCallback()
{
    GMX_RELEASE_ASSERT(hasParrinelloRahmanScaling<integrationStage>(),
                       formatString("Parrinello-Rahman scaling not implemented for %s",
                                    integrationStepNames[integrationStage])
                               .c_str());

    return [this](Step step) { scalingStepPR_ = step; };
}

template<IntegrationStage integrationStage>
static PropagatorConnection getConnection(Propagator<integrationStage>* propagator,
                                          const PropagatorTag&          propagatorTag)
{
    PropagatorConnection propagatorConnection{ propagatorTag };

    if constexpr (hasStartVelocityScaling<integrationStage>() || hasEndVelocityScaling<integrationStage>())
    {
        propagatorConnection.setNumVelocityScalingVariables =
                [propagator](int num, ScaleVelocities scaleVelocities)
        { propagator->setNumVelocityScalingVariables(num, scaleVelocities); };
        propagatorConnection.getVelocityScalingCallback = [propagator]()
        { return propagator->velocityScalingCallback(); };
    }
    if constexpr (hasStartVelocityScaling<integrationStage>())
    {
        propagatorConnection.getViewOnStartVelocityScaling = [propagator]()
        { return propagator->viewOnStartVelocityScaling(); };
    }
    if constexpr (hasEndVelocityScaling<integrationStage>())
    {
        propagatorConnection.getViewOnEndVelocityScaling = [propagator]()
        { return propagator->viewOnEndVelocityScaling(); };
    }
    if constexpr (hasPositionScaling<integrationStage>())
    {
        propagatorConnection.setNumPositionScalingVariables = [propagator](int num)
        { propagator->setNumPositionScalingVariables(num); };
        propagatorConnection.getViewOnPositionScaling = [propagator]()
        { return propagator->viewOnPositionScaling(); };
        propagatorConnection.getPositionScalingCallback = [propagator]()
        { return propagator->positionScalingCallback(); };
    }
    if constexpr (hasParrinelloRahmanScaling<integrationStage>())
    {
        propagatorConnection.getViewOnPRScalingMatrix = [propagator]()
        { return propagator->viewOnPRScalingMatrix(); };
        propagatorConnection.getPRScalingCallback = [propagator]()
        { return propagator->prScalingCallback(); };
    }

    return propagatorConnection;
}

// doxygen is confused by the two definitions
//! \cond
template<IntegrationStage integrationStage>
ISimulatorElement* Propagator<integrationStage>::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData gmx_unused*                  energyData,
        FreeEnergyPerturbationData gmx_unused*  freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused*   globalCommunicationHelper,
        ObservablesReducer* /* observablesReducer */,
        const PropagatorTag& propagatorTag,
        TimeStep             timestep)
{
    GMX_RELEASE_ASSERT(!(integrationStage == IntegrationStage::ScaleVelocities
                         || integrationStage == IntegrationStage::ScalePositions)
                               || (timestep == 0.0),
                       "Scaling elements don't propagate the system.");
    auto* element    = builderHelper->storeElement(std::make_unique<Propagator<integrationStage>>(
            timestep, statePropagatorData, legacySimulatorData->mdAtoms_, legacySimulatorData->wallCycleCounters_));
    auto* propagator = static_cast<Propagator<integrationStage>*>(element);
    builderHelper->registerPropagator(getConnection<integrationStage>(propagator, propagatorTag));
    return element;
}

template<IntegrationStage integrationStage>
ISimulatorElement* Propagator<integrationStage>::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData*             freeEnergyPerturbationData,
        GlobalCommunicationHelper*              globalCommunicationHelper,
        ObservablesReducer*                     observablesReducer,
        const PropagatorTag&                    propagatorTag)
{
    GMX_RELEASE_ASSERT(
            integrationStage == IntegrationStage::ScaleVelocities
                    || integrationStage == IntegrationStage::ScalePositions,
            "Adding a propagator without time step is only allowed for scaling elements");
    return getElementPointerImpl(legacySimulatorData,
                                 builderHelper,
                                 statePropagatorData,
                                 energyData,
                                 freeEnergyPerturbationData,
                                 globalCommunicationHelper,
                                 observablesReducer,
                                 propagatorTag,
                                 TimeStep(0.0));
}
//! \endcond

// Explicit template initializations
template class Propagator<IntegrationStage::PositionsOnly>;
template class Propagator<IntegrationStage::VelocitiesOnly>;
template class Propagator<IntegrationStage::LeapFrog>;
template class Propagator<IntegrationStage::VelocityVerletPositionsAndVelocities>;
template class Propagator<IntegrationStage::ScaleVelocities>;
template class Propagator<IntegrationStage::ScalePositions>;

} // namespace gmx
