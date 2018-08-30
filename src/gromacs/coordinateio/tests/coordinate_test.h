/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!\file
 * \libinternal
 * \brief
 * Helper classes for outputmanager and coordinateio tests
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */

#ifndef GMX_COORDINATEIO_TESTS_COORDINATEIO_H
#define GMX_COORDINATEIO_TESTS_COORDINATEIO_H

#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"
#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/coordinateio/outputmanager.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

class DummyOutputModule : public ICoordinateOutput
{
    public:
        explicit DummyOutputModule(unsigned long flag) :
            ICoordinateOutput(flag)
        {
        }
        DummyOutputModule(DummyOutputModule &&old) noexcept :
            ICoordinateOutput(old.moduleFlags_)
        {}

        ~DummyOutputModule() {}

        virtual void processFrame(int /*framenumber*/, t_trxframe * /*input*/)
        {}
};

//! Convenience typedef for dummy module
typedef std::shared_ptr<DummyOutputModule>
    DummyOutputModulePointer;

/*! \libinternal \brief
 * Tests for internal functions of the outputmanager.
 *
 * This test fixture is used for testing the internal functions of
 * the OutputManager class for file opening and module adding.
 *
 * \inlibraryapi
 * \ingroup module_coordinateio
 */
class OutputManagerTest : public ::testing::Test
{
    protected:
        //! Empty outptutmanager used for testing.
        OutputManagerPointer output_;
        //! Empty topology used for initializing outputmanager.
        gmx_mtop_t           dummyTopology_;
        //! Empty selection used for initializing outputmanager.
        Selection            dummySelection_;

        OutputManagerTest()
        {
            init_mtop(&dummyTopology_);
            output_ = compat::make_unique<OutputManager>("test.pdb", dummySelection_, &dummyTopology_);
        }
        ~OutputManagerTest()
        {
        }

    public:
        /*! \brief
         * Method to create a dummy module to add to the outputmanager during tests.
         *
         * Allocates a new minimal outputadapter module with a specified \p flag
         * and tries to add it to the outputmanager. Used for testing that adding
         * works and that only matching modules are accepted.
         *
         * \param[in] flag Flag setting propagated to the dummy module to check
         *                 that mismatched flags are rejected and matching flags accepted.
         */
        void addEmptyModule(unsigned long flag)
        {
            DummyOutputModulePointer basicModule_ = std::make_shared<DummyOutputModule>(flag);
            output_->addFlagModule(basicModule_);
        }
};

} // namespace test

} // namespace gmx

#endif
