/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include <cassert>
#include <cmath>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/shake.h"
#include "gromacs/mdlib/tests/constrtestdata.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/testasserts.h"

namespace gmx
{
class Lincs;

namespace test
{

void ShakeConstraintsRunner::applyConstraints(ConstraintsTestData* testData, t_pbc /* pbc */)
{
    shakedata shaked;
    make_shake_sblock_serial(&shaked, testData->idef_.get(), testData->numAtoms_);
    bool success = constrain_shake(nullptr,
                                   &shaked,
                                   testData->invmass_,
                                   *testData->idef_,
                                   testData->ir_,
                                   testData->x_,
                                   testData->xPrime_,
                                   testData->xPrime2_,
                                   nullptr,
                                   &testData->nrnb_,
                                   testData->lambda_,
                                   &testData->dHdLambda_,
                                   testData->invdt_,
                                   testData->v_,
                                   testData->computeVirial_,
                                   testData->virialScaled_,
                                   false,
                                   gmx::ConstraintVariable::Positions);
    EXPECT_TRUE(success) << "Test failed with a false return value in SHAKE.";
}

void LincsConstraintsRunner::applyConstraints(ConstraintsTestData* testData, t_pbc pbc)
{

    Lincs* lincsd;
    int    maxwarn         = 100;
    int    warncount_lincs = 0;
    gmx_omp_nthreads_set(ModuleMultiThread::Lincs, 1);

    // Communication record
    t_commrec cr;
    cr.nnodes = 1;
    cr.dd     = nullptr;

    // Multi-sim record
    gmx_multisim_t ms{ 1, 0, MPI_COMM_NULL, MPI_COMM_NULL };

    // Make blocka structure for faster LINCS setup
    std::vector<ListOfLists<int>> at2con_mt;
    at2con_mt.reserve(testData->mtop_.moltype.size());
    for (const gmx_moltype_t& moltype : testData->mtop_.moltype)
    {
        // This function is in constr.cpp
        at2con_mt.push_back(make_at2con(moltype,
                                        testData->mtop_.ffparams.iparams,
                                        flexibleConstraintTreatment(EI_DYNAMICS(testData->ir_.eI))));
    }
    // Initialize LINCS
    lincsd = init_lincs(nullptr,
                        testData->mtop_,
                        testData->nflexcon_,
                        at2con_mt,
                        false,
                        testData->ir_.nLincsIter,
                        testData->ir_.nProjOrder,
                        nullptr);
    set_lincs(*testData->idef_,
              testData->numAtoms_,
              testData->invmass_,
              testData->lambda_,
              EI_DYNAMICS(testData->ir_.eI),
              &cr,
              lincsd);

    // Evaluate constraints
    bool success = constrain_lincs(false,
                                   testData->ir_,
                                   0,
                                   lincsd,
                                   testData->invmass_,
                                   &cr,
                                   &ms,
                                   testData->x_.arrayRefWithPadding(),
                                   testData->xPrime_.arrayRefWithPadding(),
                                   testData->xPrime2_.arrayRefWithPadding().unpaddedArrayRef(),
                                   pbc.box,
                                   &pbc,
                                   testData->hasMassPerturbed_,
                                   testData->lambda_,
                                   &testData->dHdLambda_,
                                   testData->invdt_,
                                   testData->v_.arrayRefWithPadding().unpaddedArrayRef(),
                                   testData->computeVirial_,
                                   testData->virialScaled_,
                                   gmx::ConstraintVariable::Positions,
                                   &testData->nrnb_,
                                   maxwarn,
                                   &warncount_lincs,
                                   nullptr);
    EXPECT_TRUE(success) << "Test failed with a false return value in LINCS.";
    EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS.";
    done_lincs(lincsd);
}

} // namespace test
} // namespace gmx
