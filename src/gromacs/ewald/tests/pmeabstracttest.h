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
 * Describes abstract base classes for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#ifndef PMEABSTRACTTEST_H
#define PMEABSTRACTTEST_H

#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/scoped_cptr.h"

#include "testutils/refdata.h"

struct gmx_pme_t;

/*! \brief PME destructor function wrapper. */
void pme_free_wrapper(gmx_pme_t *pme);
/*! \brief A safe scoped pointer type for PME. */
typedef gmx::scoped_cptr<gmx_pme_t, pme_free_wrapper> pmeSafePointer;

//! Environment for getting the t_inputrec structure easily
static gmx::MDModules s_mdModules;
struct t_inputrec;

/*! \internal \brief
 * Abstract test fixture for initializing the PME test.
 */
class AbstractPmeTestBase
{
    //! Allocates the PME data
    void Init();

    public:
        //! Pointer to the static input structure
        t_inputrec           *inputRec_;
        //! PME data safe pointer
        pmeSafePointer        pme_;

        AbstractPmeTestBase() : inputRec_(s_mdModules.inputrec()){}

        //! Should contain the whole testing procedure and checks, called after SetUp()
        virtual void Check() = 0;

        //! This should run the functionality to be tested and is called by Check() in some shape or form.
        virtual void RunTest()
        {
            Init();
        }
};

/*! \internal \brief
 * Abstract test fixture for a PME death test (test that is supposed to terminate with gmx_fatal()).
 */
class AbstractPmeDeathTestBase : virtual public AbstractPmeTestBase
{
    public:
        //! Runs the test and checks for the fatal error message on the output
        void Check();
};

/*! \internal \brief
 * Abstract test fixture for a PME test that checks the outputs.
 */
class AbstractPmeCorrectnessTestBase : virtual public AbstractPmeTestBase
{
    //! Dummy reference data
    gmx::test::TestReferenceData        data_;

    public:
        //! Reference data checker
        gmx::test::TestReferenceChecker checker_;

        AbstractPmeCorrectnessTestBase() : checker_(data_.rootChecker()){}
        //! This function should check the correctness of the output, using the checker_
        virtual void CheckOutputs() = 0;
        //! Runs the test and checks the correctness of the output by calling CheckOutputs()
        void Check();
};

#endif // PMEABSTRACTTEST_H
