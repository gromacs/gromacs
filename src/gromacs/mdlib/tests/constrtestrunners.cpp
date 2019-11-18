/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Function to run SHAKE and LINCS on CPU.
 *
 * Functions used in the test to apply constraints on the test data:
 * CPU-based implementation and a stub for GPU-based implementation.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "constrtestrunners.h"

#include "config.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/shake.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/*! \brief
 * Initialize and apply SHAKE constraints.
 *
 * \param[in] testData        Test data structure.
 * \param[in] pbc             Periodic boundary data.
 */
void applyShake(ConstraintsTestData* testData, t_pbc gmx_unused pbc)
{
    shakedata* shaked = shake_init();
    make_shake_sblock_serial(shaked, &testData->idef_, testData->md_);
    bool success = constrain_shake(
            nullptr, shaked, testData->invmass_.data(), testData->idef_, testData->ir_,
            as_rvec_array(testData->x_.data()), as_rvec_array(testData->xPrime_.data()),
            as_rvec_array(testData->xPrime2_.data()), &testData->nrnb_, testData->md_.lambda,
            &testData->dHdLambda_, testData->invdt_, as_rvec_array(testData->v_.data()),
            testData->computeVirial_, testData->virialScaled_, false, gmx::ConstraintVariable::Positions);
    EXPECT_TRUE(success) << "Test failed with a false return value in SHAKE.";
    done_shake(shaked);
}

/*! \brief
 * Initialize and apply LINCS constraints.
 *
 * \param[in] testData        Test data structure.
 * \param[in] pbc             Periodic boundary data.
 */
void applyLincs(ConstraintsTestData* testData, t_pbc pbc)
{

    Lincs* lincsd;
    int    maxwarn         = 100;
    int    warncount_lincs = 0;
    gmx_omp_nthreads_set(emntLINCS, 1);

    // Make blocka structure for faster LINCS setup
    std::vector<t_blocka> at2con_mt;
    at2con_mt.reserve(testData->mtop_.moltype.size());
    for (const gmx_moltype_t& moltype : testData->mtop_.moltype)
    {
        // This function is in constr.cpp
        at2con_mt.push_back(make_at2con(moltype, testData->mtop_.ffparams.iparams,
                                        flexibleConstraintTreatment(EI_DYNAMICS(testData->ir_.eI))));
    }
    // Initialize LINCS
    lincsd = init_lincs(nullptr, testData->mtop_, testData->nflexcon_, at2con_mt, false,
                        testData->ir_.nLincsIter, testData->ir_.nProjOrder);
    set_lincs(testData->idef_, testData->md_, EI_DYNAMICS(testData->ir_.eI), &testData->cr_, lincsd);

    // Evaluate constraints
    bool success = constrain_lincs(
            false, testData->ir_, 0, lincsd, testData->md_, &testData->cr_, &testData->ms_,
            as_rvec_array(testData->x_.data()), as_rvec_array(testData->xPrime_.data()),
            as_rvec_array(testData->xPrime2_.data()), pbc.box, &pbc, testData->md_.lambda,
            &testData->dHdLambda_, testData->invdt_, as_rvec_array(testData->v_.data()),
            testData->computeVirial_, testData->virialScaled_, gmx::ConstraintVariable::Positions,
            &testData->nrnb_, maxwarn, &warncount_lincs);
    EXPECT_TRUE(success) << "Test failed with a false return value in LINCS.";
    EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS.";
    for (auto& moltype : at2con_mt)
    {
        sfree(moltype.index);
        sfree(moltype.a);
    }
    done_lincs(lincsd);
}

#if GMX_GPU != GMX_GPU_CUDA
/*! \brief
 * Stub for LINCS constraints on CUDA-enabled GPU to satisfy compiler.
 *
 * \param[in] testData        Test data structure.
 * \param[in] pbc             Periodic boundary data.
 */
void applyLincsCuda(ConstraintsTestData gmx_unused* testData, t_pbc gmx_unused pbc)
{
    FAIL() << "Dummy LINCS CUDA function was called instead of the real one.";
}
#endif

} // namespace test
} // namespace gmx
