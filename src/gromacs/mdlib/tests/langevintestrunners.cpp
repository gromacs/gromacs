/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Runner for CPU-based implementation of the Langevin integrator.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "langevintestrunners.h"

#include <array>
#include <memory>
#include <utility>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/matrix.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/arrayref.h"

#include "langevintestdata.h"


namespace gmx
{
namespace test
{

void LangevinHostTestRunner::integrate(LangevinTestData* testData, int numSteps)
{
    testData->state_.x.resizeWithPadding(testData->numAtoms_);
    testData->state_.v.resizeWithPadding(testData->numAtoms_);
    for (int i = 0; i < testData->numAtoms_; i++)
    {
        testData->state_.x[i] = testData->x_[i];
        testData->state_.v[i] = testData->v_[i];
    }

    gmx_omp_nthreads_set(ModuleMultiThread::Update, 1);

    Matrix3x3 parrinelloRahmanM{ { 0. } };

    for (int step = 0; step < numSteps; step++)
    {
        testData->update_->update_coords(testData->inputRecord_,
                                         step,
                                         testData->mdAtoms_.homenr,
                                         testData->mdAtoms_.havePartiallyFrozenAtoms,
                                         testData->mdAtoms_.ptype,
                                         testData->mdAtoms_.invmass,
                                         testData->mdAtoms_.invMassPerDim,
                                         &testData->state_,
                                         testData->f_,
                                         &testData->forceCalculationData_,
                                         &testData->kineticEnergyData_,
                                         parrinelloRahmanM,
                                         etrtNONE,
                                         nullptr,
                                         false);
        testData->update_->finish_update(testData->inputRecord_,
                                         testData->mdAtoms_.havePartiallyFrozenAtoms,
                                         testData->mdAtoms_.homenr,
                                         &testData->state_,
                                         nullptr,
                                         false);
    }
    const auto xp = makeArrayRef(*testData->update_->xp()).subArray(0, testData->numAtoms_);
    for (int i = 0; i < testData->numAtoms_; i++)
    {
        for (int d = 0; d < DIM; d++)
        {
            testData->x_[i][d]      = testData->state_.x[i][d];
            testData->v_[i][d]      = testData->state_.v[i][d];
            testData->xPrime_[i][d] = xp[i][d];
        }
    }
}

} // namespace test
} // namespace gmx
