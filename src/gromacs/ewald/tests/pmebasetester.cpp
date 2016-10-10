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
 * Implements tester functions and common routines for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pmebasetester.h"

#include <gmock/gmock.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/utility/gmxassert.h"

// Testers

//! A tester function which tests RunTest() for graceful termination with gmx_fatal()
void TestForFatalError(const TestFunc &RunTest)
{
    EXPECT_DEATH(RunTest(), "Fatal error");
    /* Note that this tests only for a general gmx_fatal error message
     * (looking for "Fatal error" substring in the program output).
     * The function being tested can fail not where we expect it to,
     * and still pass the death test.
     * We would ideally want to pass the specific error message here,
     * but first we would like to have a single program-wide error message database,
     * to not introduce horrible duplication.
     */
}
//! A global pointer to the death tester
TesterFuncPtr g_TestForFatalErrorPtr = &TestForFatalError;

//! A tester function which tests RunTest() for successful run with no failures
void TestForSuccess(const TestFunc &RunTest)
{
    EXPECT_NO_FATAL_FAILURE(RunTest());
    /* It seems that at least gmx_fatal bypasses EXPECT_NO_FATAL_FAILURE anyway,
     * halting the whole test execution.
     * (Tested by mixing obviously invalid input among valid ones).
     */
}
//! A global pointer to the success tester
TesterFuncPtr g_TestForSuccessPtr = &TestForSuccess;

// Consider having a segfault tester as well


// PME-specific stuff

//! PME destructor function wrapper for the safe pointer.
void pme_free_wrapper(gmx_pme_t *pme)
{
    if (pme)
    {
        gmx_pme_destroy(&pme);
    }
}

//! Simple PME initialization based on input
pmeSafePointer PmeDummyInit(const t_inputrec *inputRec)
{
    gmx_pme_t *pmeDataRaw = NULL;
    gmx_pme_init(&pmeDataRaw, NULL, 1, 1, inputRec,
                 0, FALSE, FALSE, TRUE, 0.0, 0.0, 1);
    pmeSafePointer pme(pmeDataRaw, pme_free_wrapper); // taking ownership
    return pme;
}
