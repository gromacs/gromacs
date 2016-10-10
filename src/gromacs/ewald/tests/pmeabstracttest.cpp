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
 * Implements abstract base class and tester functions for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pmeabstracttest.h"

#include <gmock/gmock.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/gmxassert.h"

// Testers

//! \brief A tester function which tests for graceful termination with gmx_fatal()
void TestForFatalError(ITest *test)
{
    EXPECT_DEATH(test->RunTest(), "Fatal error");
}
//! \brief A global pointer to the death tester
TesterFuncPtr g_TestForFatalErrorPtr = &TestForFatalError;

//! \brief A tester function which tests the output with CheckOutputs()
void TestForCorrectness(ITest *test)
{
    test->RunTest();
    test->CheckOutputs();
}
//! \brief A global pointer to the correctness tester
TesterFuncPtr g_TestForCorrectnessPtr = &TestForCorrectness;

// PME-specific stuff

//! Environment for getting the t_inputrec structure easily
static gmx::MDModules s_mdModules;
t_inputrec           *AbstractPmeTest::inputRec_ = s_mdModules.inputrec();

/*! \brief PME destructor function wrapper for the safe pointer. */
void pme_free_wrapper(gmx_pme_t *pme)
{
    if (pme)
    {
        gmx_pme_destroy(&pme);
    }
}

//! PME initialization based on input
pmeSafePointer PmeDummyInit(const t_inputrec *inputRec)
{
    gmx_pme_t *pmeDataRaw = NULL;
    gmx_pme_init(&pmeDataRaw, NULL, 1, 1, inputRec,
                 0, FALSE, FALSE, TRUE, 0.0, 0.0, 1);
    pmeSafePointer pme(pmeDataRaw, pme_free_wrapper); // taking ownership
    return pme;
}

void AbstractPmeTest::PmeInit()
{
    pme_ = PmeDummyInit(inputRec_);
}

void AbstractPmeTest::Check()
{
    GMX_ASSERT((void *)testerFunction_ != NULL, "No tester function specified");
    testerFunction_(static_cast<ITest *>(this));
}
