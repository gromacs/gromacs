/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief
 * Describes common routines and types for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_PME_TEST_COMMON_H
#define GMX_EWALD_PME_TEST_COMMON_H

#include <cstdint>

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_gpu_internal.h"
#include "gromacs/ewald/pme_gpu_program.h"
#include "gromacs/ewald/pme_output.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/test_device.h"

class DeviceContext;
class DeviceStream;
struct gmx_pme_t;
struct t_inputrec;

namespace gmx
{

template<typename>
class ArrayRef;

namespace test
{

//! Hardware code path being tested
enum class CodePath : int
{
    //! CPU code path
    CPU,
    //! GPU code path
    GPU,
    //! Total number of code paths
    Count
};

// Convenience typedefs
//! A safe pointer type for PME.
typedef gmx::unique_cptr<gmx_pme_t, gmx_pme_destroy> PmeSafePointer;
//! Charges
typedef std::vector<real> ChargesVector;
//! Coordinates
typedef std::vector<RVec> CoordinatesVector;
//! Forces
typedef ArrayRef<RVec> ForcesVector;
//! Gridline indices
typedef std::vector<IVec> GridLineIndicesVector;
/*! \brief Spline parameters (theta or dtheta).
 * A reference to a single dimension's spline data; this means (atomCount * pmeOrder) values or derivatives.
 */
typedef ArrayRef<const real> SplineParamsDimVector;
/*! \brief Spline parameters (theta or dtheta) in all 3 dimensions
 */
typedef std::array<SplineParamsDimVector, DIM> SplineParamsVector;

//! Non-zero grid values for test input; keys are 3d indices (IVec)
template<typename ValueType>
using SparseGridValuesInput = std::map<IVec, ValueType>;
//! Non-zero real grid values
typedef SparseGridValuesInput<real> SparseRealGridValuesInput;
//! Non-zero complex grid values
typedef SparseGridValuesInput<t_complex> SparseComplexGridValuesInput;
//! Non-zero grid values for test output; keys are string representations of the cells' 3d indices (IVec); this allows for better sorting.
template<typename ValueType>
using SparseGridValuesOutput = std::map<std::string, ValueType>;
//! Non-zero real grid values
typedef SparseGridValuesOutput<real> SparseRealGridValuesOutput;
//! Non-zero complex grid values
typedef SparseGridValuesOutput<t_complex> SparseComplexGridValuesOutput;
//! TODO: make proper C++ matrix for the whole Gromacs, get rid of this
typedef std::array<real, DIM * DIM> Matrix3x3;
//! PME solver type
enum class PmeSolveAlgorithm : int
{
    //! Coulomb electrostatics
    Coulomb,
    //! Lennard-Jones
    LennardJones,
    //! Total number of solvers
    Count
};

// Misc.

/*! \brief Returns message describing why PME in this \c mode is not
 * supported for this \c inputRec, or empty of messages when PME is supported. */
MessageStringCollector getSkipMessagesIfNecessary(const t_inputrec& inputRec, CodePath mode);

//! Spline moduli are computed in double precision, so they're very good in single precision
constexpr int64_t c_splineModuliSinglePrecisionUlps = 1;
/*! \brief For double precision checks, the recursive interpolation
 * and use of trig functions in make_dft_mod require a lot more flops,
 * and thus opportunity for deviation between implementations. */
uint64_t getSplineModuliDoublePrecisionUlps(int splineOrder);

// PME stages

//! PME initialization
PmeSafePointer pmeInitWrapper(const t_inputrec*    inputRec,
                              CodePath             mode,
                              const DeviceContext* deviceContext,
                              const DeviceStream*  deviceStream,
                              const PmeGpuProgram* pmeGpuProgram,
                              const Matrix3x3&     box,
                              real                 ewaldCoeff_q  = 1.0F,
                              real                 ewaldCoeff_lj = 1.0F);

//! Simple PME initialization based on inputrec only
PmeSafePointer pmeInitEmpty(const t_inputrec* inputRec);

//! Make a GPU state-propagator manager
std::unique_ptr<StatePropagatorDataGpu> makeStatePropagatorDataGpu(const gmx_pme_t&     pme,
                                                                   const DeviceContext* deviceContext,
                                                                   const DeviceStream* deviceStream);
//! PME initialization with atom data and system box
void pmeInitAtoms(gmx_pme_t*               pme,
                  StatePropagatorDataGpu*  stateGpu,
                  CodePath                 mode,
                  const CoordinatesVector& coordinates,
                  const ChargesVector&     charges);
//! PME spline computation and charge spreading
void pmePerformSplineAndSpread(gmx_pme_t* pme, CodePath mode, bool computeSplines, bool spreadCharges);
//! PME solving
void pmePerformSolve(gmx_pme_t*        pme,
                     CodePath          mode,
                     PmeSolveAlgorithm method,
                     real              cellVolume,
                     GridOrdering      gridOrdering,
                     bool              computeEnergyAndVirial);
//! PME force gathering
void pmePerformGather(gmx_pme_t*    pme,
                      CodePath      mode,
                      ForcesVector& forces); //NOLINT(google-runtime-references)
//! PME test finalization before fetching the outputs
void pmeFinalizeTest(const gmx_pme_t* pme, CodePath mode);

// PME state setters

//! Setting atom spline values or derivatives to be used in spread/gather
void pmeSetSplineData(const gmx_pme_t*             pme,
                      CodePath                     mode,
                      const SplineParamsDimVector& splineValues,
                      PmeSplineDataType            type,
                      int                          dimIndex);

//! Setting gridline indices be used in spread/gather
void pmeSetGridLineIndices(gmx_pme_t* pme, CodePath mode, const GridLineIndicesVector& gridLineIndices);
//! Setting real grid to be used in gather
void pmeSetRealGrid(gmx_pme_t* pme, CodePath mode, const SparseRealGridValuesInput& gridValues);
void pmeSetComplexGrid(gmx_pme_t*                          pme,
                       CodePath                            mode,
                       GridOrdering                        gridOrdering,
                       const SparseComplexGridValuesInput& gridValues);

// PME state getters

//! Getting the single dimension's spline values or derivatives
SplineParamsDimVector pmeGetSplineData(const gmx_pme_t* pme, CodePath mode, PmeSplineDataType type, int dimIndex);
//! Getting the gridline indices
GridLineIndicesVector pmeGetGridlineIndices(const gmx_pme_t* pme, CodePath mode);
//! Getting the real grid (spreading output of pmePerformSplineAndSpread())
SparseRealGridValuesOutput pmeGetRealGrid(gmx_pme_t* pme, CodePath mode);
//! Getting the complex grid output of pmePerformSolve()
SparseComplexGridValuesOutput pmeGetComplexGrid(gmx_pme_t* pme, CodePath mode, GridOrdering gridOrdering);
//! Getting the reciprocal energy and virial
PmeOutput pmeGetReciprocalEnergyAndVirial(const gmx_pme_t* pme, CodePath mode, PmeSolveAlgorithm method);

struct PmeTestHardwareContext
{
    //! Hardware path for the code being tested.
    CodePath codePath_;
    //! Returns a human-readable context description line
    std::string description() const;
    //! Returns an optional GPU ID, with a valid value when the context is for a GPU
    std::optional<int> gpuId() const;
    //! Pointer to the global test hardware device (if on GPU)
    TestDevice* testDevice_ = nullptr;
    //! PME GPU program if needed
    PmeGpuProgramStorage pmeGpuProgram_ = nullptr;
    //! Constructor for CPU context
    PmeTestHardwareContext();
    //! Constructor for GPU context
    explicit PmeTestHardwareContext(TestDevice* testDevice);

    //! Get the code path
    CodePath codePath() const { return codePath_; }
    //! Get the PME GPU program
    const PmeGpuProgram* pmeGpuProgram() const
    {
        return codePath() == CodePath::GPU ? pmeGpuProgram_.get() : nullptr;
    }

    const DeviceContext* deviceContext() const
    {
        return codePath() == CodePath::GPU ? &testDevice_->deviceContext() : nullptr;
    }

    const DeviceStream* deviceStream() const
    {
        return codePath() == CodePath::GPU ? &testDevice_->deviceStream() : nullptr;
    }

    //! Activate the context (set the device)
    void activate() const;
};

//! Return a view of the current PME test hardware contexts
ArrayRef<const PmeTestHardwareContext> getPmeTestHardwareContexts();

/*! \brief Construct a refdata filename for a test
 *
 * We want the same reference data to apply to every hardware context
 * for which we test PME. That means we need to store it in a file
 * whose name relates to the name of the test, but excluding the part
 * related to the context. */
std::string makeRefDataFileName();

/*! \brief Make a terse description of the hardware context suitable
 * for use in naming the test case.
 *
 * The full hardware device description is used in a SCOPED_TRACE
 * message, but that is too long for the test name and has too many
 * possible characters that might break GoogleTest. */
std::string makeHardwareContextName(int hardwareContextIndex);

/*! \brief Functions that dynamically register test cases
 *
 * This are called by registerTestsDynamically and customize the
 * range of test cases to suit the available hardware. */
//!\{
void registerDynamicalPmeSplineSpreadTests(Range<int> contextIndexRange);
void registerDynamicalPmeSolveTests(Range<int> contextIndexRange);
void registerDynamicalPmeGatherTests(Range<int> contextIndexRange);
//!\}

//! A couple of valid inputs for boxes.
extern const std::map<std::string, Matrix3x3> c_inputBoxes;

//! Valid PME orders for testing
extern std::vector<int> c_inputPmeOrders;

} // namespace test
} // namespace gmx

#endif // GMX_EWALD_PME_TEST_COMMON_H
