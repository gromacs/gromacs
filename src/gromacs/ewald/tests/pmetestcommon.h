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

#include "gromacs/ewald/pme.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/unique_cptr.h"

struct t_inputrec;

namespace gmx
{

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
//! Non-zero grid values by their gridline indices
typedef SparseGridValues<real> SparseRealGridValues;
//! Non-zero complex grid values by their gridline indices
typedef SparseGridValues<t_complex> SparseComplexGridValues;
//! TODO: make proper C++ matrix for the whole Gromacs, get rid of this
typedef std::array<real, DIM * DIM> Matrix3x3;
//! PME code path being tested
enum class PmeCodePath
{
    CPU, // serial CPU code
};
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

// PME stages

//! Simple PME initialization based on input, no atom data
PmeSafePointer PmeInitEmpty(const t_inputrec *inputRec,
                            const Matrix3x3 &box = {1, 0, 0, 0, 1, 0, 0, 0, 1},
                            real ewaldCoeff_q = 0.0f, real ewaldCoeff_lj = 0.0f);
//! PME initialization with atom data
PmeSafePointer PmeInitWithAtoms(const t_inputrec        *inputRec,
                                const CoordinatesVector &coordinates,
                                const ChargesVector     &charges,
                                const Matrix3x3         &box
                                );
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
void PmeSetRealGrid(const PmeSafePointer       &pmeSafe,
                    PmeCodePath                 mode,
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
