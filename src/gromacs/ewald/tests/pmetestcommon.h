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
#include <map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-gpu-internal.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/unique_cptr.h"

#include "testhardwarecontexts.h"

namespace gmx
{
namespace test
{

// Convenience typedefs
//! A safe pointer type for PME.
typedef gmx::unique_cptr<gmx_pme_t, gmx_pme_destroy> PmeSafePointer;
//! Charges
typedef ConstArrayRef<real> ChargesVector;
//! Coordinates
typedef std::vector<RVec> CoordinatesVector;
//! Forces
typedef ArrayRef<RVec> ForcesVector;
//! Gridline indices
typedef ConstArrayRef<IVec> GridLineIndicesVector;
/*! \brief Spline parameters (theta or dtheta).
 * A reference to a single dimension's spline data; this means (atomCount * pmeOrder) values or derivatives.
 */
typedef ConstArrayRef<real> SplineParamsDimVector;
/*! \brief Spline parameters (theta or dtheta) in all 3 dimensions
 */
typedef std::array<SplineParamsDimVector, DIM> SplineParamsVector;

//! Non-zero grid values for test input; keys are 3d indices (IVec)
template<typename ValueType>using SparseGridValuesInput = std::map<IVec, ValueType>;
//! Non-zero real grid values
typedef SparseGridValuesInput<real> SparseRealGridValuesInput;
//! Non-zero complex grid values
typedef SparseGridValuesInput<t_complex> SparseComplexGridValuesInput;
//! Non-zero grid values for test output; keys are string representations of the cells' 3d indices (IVec); this allows for better sorting.
template<typename ValueType>using SparseGridValuesOutput = std::map<std::string, ValueType>;
//! Non-zero real grid values
typedef SparseGridValuesOutput<real> SparseRealGridValuesOutput;
//! Non-zero complex grid values
typedef SparseGridValuesOutput<t_complex> SparseComplexGridValuesOutput;
//! TODO: make proper C++ matrix for the whole Gromacs, get rid of this
typedef std::array<real, DIM * DIM> Matrix3x3;
//! PME gathering input forces treatment
enum class PmeGatherInputHandling
{
    Overwrite,
    ReduceWith,
};
//! PME solver type
enum class PmeSolveAlgorithm
{
    Normal,
    LennardJones,
};
//! PME solve grid ordering //TODO rename
enum class PmeSolveOptions
{
    YZXGridOrder,
    XYZGridOrder
};
//! PME solver results - reciprocal energy and virial
typedef std::tuple<real, Matrix3x3> PmeSolveOutput;

// Misc.

//! Tells if this generally valid PME input is supported for this mode
bool pmeSupportsInputForMode(const t_inputrec *inputRec, CodePath mode);

// PME stages

//! Simple PME initialization based on input, no atom data
PmeSafePointer pmeInitEmpty(const t_inputrec *inputRec,
                            CodePath mode = CodePath::CPU,
                            const Matrix3x3 &box = {{1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f}},
                            real ewaldCoeff_q = 0.0f, real ewaldCoeff_lj = 0.0f);
//! PME initialization with atom data and system box
PmeSafePointer pmeInitWithAtoms(const t_inputrec         *inputRec,
                                CodePath                  mode,
                                const CoordinatesVector  &coordinates,
                                const ChargesVector      &charges,
                                const Matrix3x3          &box
                                );
//! PME spline computation and charge spreading
void pmePerformSplineAndSpread(gmx_pme_t *pme, CodePath mode,
                               bool computeSplines, bool spreadCharges);
//! PME solving
void pmePerformSolve(const gmx_pme_t *pme, CodePath mode,
                     PmeSolveAlgorithm method, real cellVolume,
                     PmeSolveOptions option);
//! PME force gathering
void pmePerformGather(gmx_pme_t *pme, CodePath mode,
                      PmeGatherInputHandling inputTreatment, ForcesVector &forces);
//! PME test finalization before fetching the outputs
void pmeFinalizeTest(const gmx_pme_t *pme, CodePath mode);

// PME state setters

//! Setting atom spline values or derivatives to be used in spread/gather
void pmeSetSplineData(const gmx_pme_t *pme, CodePath mode,
                      const SplineParamsDimVector &splineValues, PmeSplineDataType type, int dimIndex);
//! Setting gridline indices be used in spread/gather
void pmeSetGridLineIndices(const gmx_pme_t *pme, CodePath mode,
                           const GridLineIndicesVector &gridLineIndices);
//! Setting real grid to be used in gather
void pmeSetRealGrid(const gmx_pme_t *pme, CodePath mode,
                    const SparseRealGridValuesInput &gridValues);
void pmeSetComplexGrid(const gmx_pme_t *pme, CodePath mode, PmeSolveOptions option,
                       const SparseComplexGridValuesInput &gridValues);

// PME state getters

//! Getting the single dimension's spline values or derivatives
SplineParamsDimVector pmeGetSplineData(const gmx_pme_t *pme, CodePath mode,
                                       PmeSplineDataType type, int dimIndex);
//! Getting the gridline indices
GridLineIndicesVector pmeGetGridlineIndices(const gmx_pme_t *pme, CodePath mode);
//! Getting the real grid (spreading output of pmePerformSplineAndSpread())
SparseRealGridValuesOutput pmeGetRealGrid(const gmx_pme_t *pme, CodePath mode);
//! Getting the complex grid output of pmePerformSolve()
SparseComplexGridValuesOutput pmeGetComplexGrid(const gmx_pme_t *pme, CodePath mode,
                                                PmeSolveOptions option);
//! Getting the reciprocal energy and virial
PmeSolveOutput pmeGetReciprocalEnergyAndVirial(const gmx_pme_t *pme, CodePath mode,
                                               PmeSolveAlgorithm method);
}
}

#endif
