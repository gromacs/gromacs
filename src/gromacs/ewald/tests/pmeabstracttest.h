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
 * Describes abstract base class and tester functions for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_PME_TEST_BASE_H
#define GMX_EWALD_PME_TEST_BASE_H

#include <memory>
#include <string>

struct gmx_pme_t;
struct t_inputrec;

/*! \internal \brief
 * A basic test interface.
 */
class ITest
{
    public:
        //! This should run the functionality to be tested.
        virtual void RunTest() = 0;
        //! This should check the results of RunTest() if it succeeded
        virtual void CheckOutputs() = 0;
};

//! \brief A tester function pointer; a tester function should accept the ITest interface and run with it
typedef void (*TesterFuncPtr)(ITest *);
//! \brief A global pointer to the death tester
extern TesterFuncPtr g_TestForFatalErrorPtr;
//! \brief A global pointer to the correctness tester
extern TesterFuncPtr g_TestForCorrectnessPtr;

/*! \brief A safe pointer type for PME. */
typedef std::unique_ptr<gmx_pme_t, void(*)(gmx_pme_t *)> pmeSafePointer;

/*! \internal \brief
 * An abstract PME test.
 */
class AbstractPmeTest : public ITest
{
    protected:
        //! Unique test ID for the reference data
        std::string   testId_;
        //! Tester function pointer
        TesterFuncPtr testerFunction_;

        //! Input data structure pointer
        static t_inputrec *inputRec_;
        //! PME data
        pmeSafePointer     pme_;
        //! PME initialization based on input
        void PmeInit();
    public:
        AbstractPmeTest() : testerFunction_(NULL), pme_(NULL, NULL) {}
        //! Should be called by unit tests - calls testerFunction_() inside
        void Check();
};

#endif // GMX_EWALD_PME_TEST_BASE_H
