/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#include <list>
#include <map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/unique_cptr.h"

//FIXME included for contexts
#include "gromacs/utility/stringutil.h"
#include "gromacs/ewald/pme-gpu-internal.h"

namespace gmx
{

//! An interface to describe a hardware context
class ITestHardwareContext
{
    public:
        //! Activates the context
        virtual void activate() = 0;
        //! Gets a human-readable context description line
        virtual std::string getDescription() = 0;
};

//! A default empty context - should be used when we don't actually have any differing contexts
class EmptyTestHardwareContext : public ITestHardwareContext
{
    public:
        void activate(){}
        std::string getDescription() {return std::string(""); }
};

//! A GPU context
class GpuTestHardwareContext : public ITestHardwareContext
{
    // a GPU id - such as CUDA device id
    int id_;
    public:
        GpuTestHardwareContext(int id) {id_ = id; }
        void activate(){myInitGpu(id_); }
        std::string getDescription() {return formatString("GPU %d", id_); } //TODO maybe somethign from gpuutils
};

// Convenience typedefs
//! A safe pointer type for PME.
typedef gmx::unique_cptr<gmx_pme_t, gmx_pme_destroy> PmeSafePointer;
//! Charges
typedef std::vector<real> ChargesVector;
//! Coordinates
typedef std::vector<RVec> CoordinatesVector;
//! Forces
typedef std::vector<RVec> ForcesVector;
//! Gridline indices
typedef std::vector<IVec> GridLineIndicesVector;
//! Spline parameters (theta or dtheta)
typedef std::vector<real> SplineParamsVector;
//! Non-zero grid values accessed by their gridline indices
template<typename ValueType>using SparseGridValues = std::map<IVec, ValueType>;
//! Non-zero real grid values accessed by their gridline indices
typedef SparseGridValues<real> SparseRealGridValues;
//! Non-zero complex grid values accessed by their gridline indices
typedef SparseGridValues<t_complex> SparseComplexGridValues;
//! TODO: make proper C++ matrix for the whole Gromacs, get rid of this
typedef std::array<real, DIM * DIM> Matrix3x3;
//! PME code path being tested
enum class PmeCodePath
{
    CPU, // serial CPU code
    CUDA
};
//! A list of hardware contexts
typedef std::list<std::shared_ptr<ITestHardwareContext> > TestHardwareContexts;
//! Type of spline data
enum class PmeSplineDataType
{
    Values,      // theta
    Derivatives, // dtheta
};
//! PME solver type
enum class PmeSolveAlgorithm
{
    Normal,
    LennardJones,
};
//! PME gathering input forces treatment
enum class PmeGatherInputHandling
{
    Overwrite,
    ReduceWith,
};

// Misc.

//! \internal \brief This class performs one-time test initialization (enumerating the hardware)
class PmeTestEnvironment : public ::testing::Environment
{
    private:
        //! General hardware info
        unique_cptr<gmx_hw_info_t, gmx_hardware_info_free> hardwareInfo_;
        //TODO: Should gpu_options live here as well? Storing gmx_gpu_opt_t is wrong (runtime per-test-case info)!
        //! GPU assignment information
        std::unique_ptr<gmx_gpu_opt_t>              gpuOptions_;
        //! Dummy communication structure which the tests do not really care about currently
        unique_cptr<t_commrec, done_commrec>        commrec_;
        //! Storage of hardware contexts
        std::map<PmeCodePath, TestHardwareContexts> hardwareContextsByMode_;

        //! Simple GPU initialization, allowing for PME to work on GPU
        void hardwareInit();

    public:
        //! Default
        ~PmeTestEnvironment() = default; //destryoing contexts?
        //! Is called once to query the hardware
        void SetUp();
        //! Get hardware information
        const gmx_hw_info_t *getHardwareInfo(){return hardwareInfo_.get(); }
        //! Get GPU information
        const gmx_gpu_opt_t *getGpuOptions(){return gpuOptions_.get(); }
        //! Get hardware contexts for given code path
        const TestHardwareContexts &getHardwareContexts(PmeCodePath mode){return hardwareContextsByMode_.at(mode); }
};

/*! \brief
 * Returns all available contexts for given code path.
 * Currently returns single empty context for CPU, and CUDA contexts for all visible CUDA-capable GPU's.
 */
const TestHardwareContexts &GetContextsForMode(PmeCodePath mode);
//! Tells if this generally valid PME input is supported for this mode
bool PmeSupportsInputForMode(const t_inputrec *inputRec, PmeCodePath mode);

// PME stages

//! Simple PME initialization based on input, no atom data
PmeSafePointer PmeInitEmpty(const t_inputrec *inputRec,
                            const Matrix3x3 &box = {{1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f}},
                            real ewaldCoeff_q = 0.0f, real ewaldCoeff_lj = 0.0f,
                            PmeCodePath mode = PmeCodePath::CPU);
//! PME initialization with atom data
PmeSafePointer PmeInitWithAtoms(const t_inputrec        *inputRec,
                                PmeCodePath              mode,
                                const CoordinatesVector &coordinates,
                                const ChargesVector     &charges,
                                const Matrix3x3         &box);
//! PME spline computation and charge spreading
void PmePerformSplineAndSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                               bool computeSplines, bool spreadCharges);
//! PME solving
void PmePerformSolve(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                     PmeSolveAlgorithm method, real cellVolume);
//! PME force gathering
void PmePerformGather(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                      PmeGatherInputHandling inputTreatment, ForcesVector &forces);

// PME stage inputs - setters for skipping stages

//! Setting atom spline values or derivatives to be used in spread/gather
void PmeSetSplineData(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                      const SplineParamsVector &splineValues, PmeSplineDataType type);
//! Setting gridline indices be used in spread/gather
void PmeSetGridLineIndices(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           const GridLineIndicesVector &gridLineIndices);
//! Setting real grid to be used in gather
void PmeSetRealGrid(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                    const SparseRealGridValues &gridValues);
//! Setting complex grid to be used in solve
void PmeSetComplexGrid(const PmeSafePointer          &pmeSafe,
                       PmeCodePath                    mode,
                       const SparseComplexGridValues &gridValues);

// PME stage outputs
//! Fetching the spline computation outputs of PmePerformSplineAndSpread()
void PmeFetchOutputsSpline(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SplineParamsVector &splineValues,
                           SplineParamsVector &splineDerivatives,
                           GridLineIndicesVector &gridLineIndices);

//! Fetching the outputs of PmePerformSolve()
void PmeFetchOutputsSolve(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                          PmeSolveAlgorithm method,
                          SparseComplexGridValues &gridValues,
                          real &energy,
                          Matrix3x3 &virial);

//! Fetching the spreading output of PmePerformSplineAndSpread()
void PmeFetchOutputsSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SparseRealGridValues &gridValues);

// Fetching the output of PmePerformGather() is not needed since it gets
// passed the forces buffer (the only output)
}
#endif
