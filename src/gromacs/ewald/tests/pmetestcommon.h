/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * \brief
 * Describes common routines and types for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_PME_TEST_COMMON_H
#define GMX_EWALD_PME_TEST_COMMON_H

#include <array>
#include <map>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/scoped_cptr.h"

struct gmx_gpu_opt_t;
struct gmx_hw_info_t;
struct gmx_pme_t;
struct t_inputrec;

namespace gmx
{

//! \internal \brief This class performs one-time test initialization (enumerating the hardware)
class PmeTestEnvironment : public ::testing::Environment
{
    private:
        //! General hardware info
        scoped_cptr<gmx_hw_info_t, gmx_hardware_info_free> hardwareInfo_;
        //TODO: Should gpu_options live here as well? Storing gmx_gpu_opt_t is wrong (runtime per-test-case info)!
        //! GPU assignment information
        std::unique_ptr<gmx_gpu_opt_t>       gpuOptions_;
        //! Dummy communication structure which the test does not really care about
        scoped_cptr<t_commrec, done_commrec> commrec_;

        //! Simple GPU initialization, allowing for PME to work on GPU
        //! \todo: Do we want to run tests on several GPUs, maybe?
        void hardwareInit();

    public:
        //! Default
        ~PmeTestEnvironment() = default;
        //! Is called once to query the hardware
        void SetUp();
        //! Get hardware information
        const gmx_hw_info_t *getHardwareInfo(){return hardwareInfo_.get(); }
        //! Get GPU information
        const gmx_gpu_opt_t *getGpuOptions(){return gpuOptions_.get(); }
};

//! The test environment
//extern PmeTestEnvironment *const pmeEnv;

// Convenience typedefs
//! A safe pointer type for PME.
typedef std::unique_ptr<gmx_pme_t, void(*)(gmx_pme_t *)> PmeSafePointer;
//! Charges
typedef std::vector<real> ChargesVector;
//! Coordinates
typedef std::vector<RVec> CoordinatesVector;
//! Gridline indices
typedef std::vector<IVec> GridLineIndicesVector;
//! Spline parameters (theta or dtheta)
typedef std::vector<real> SplineParamsVector;
//! Non-zero grid values
typedef std::map<IVec, real> SparseGridValues;
//! TODO: make proper C++ matrix for the whole Gromacs, get rid of this
typedef std::array<real, DIM * DIM> Matrix3x3;
//! PME code path being tested
enum class PmeCodePath
{
    CPU, // serial CPU code
    CUDA
};

// Misc.

//! Tells if this generally valid PME input is supported for this mode
bool PmeSupportsInputForMode(const t_inputrec *inputRec, PmeCodePath mode);

// PME stages

//! Simple PME initialization based on input, no atom data; only good for testing the initialization stage
PmeSafePointer PmeInitEmpty(const t_inputrec *inputRec);
//! PME initialization with atom data and system box
PmeSafePointer PmeInitWithAtoms(const t_inputrec        *inputRec, PmeCodePath mode,
                                const CoordinatesVector &coordinates,
                                const ChargesVector     &charges,
                                const gmx::Matrix3x3    &box
                                );
//! PME spline computation and charge spreading
void PmePerformSplineAndSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                               bool computeSplines, bool spreadCharges);

// PME stage outputs
//! Fetching the spline computation outputs of PmePerformSplineAndSpread()
void PmeFetchOutputsSpline(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SplineParamsVector &splineValues,
                           SplineParamsVector &splineDerivatives,
                           GridLineIndicesVector &gridLineIndices);

//! Fetching the spreading output of PmePerformSplineAndSpread()
void PmeFetchOutputsSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SparseGridValues &gridValues);
}
#endif
