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
#include "gromacs/coordinateio/builder.h"
#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/coordinateio/outputmanager.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

class DummyOutputModule : public ICoordinateOutput
{
    public:
        explicit DummyOutputModule(unsigned long requirementsFlag) :
            moduleRequirements(requirementsFlag) {}

        DummyOutputModule(DummyOutputModule &&old) noexcept = default;

        ~DummyOutputModule() override {}

        void processFrame(int /*framenumber*/, t_trxframe * /*input*/) override
        {}

        unsigned long getModuleFlag() override { return moduleRequirements; }

        //! Local requirements
        unsigned long moduleRequirements;
};

//! Convenience typedef for dummy module
typedef std::unique_ptr<DummyOutputModule>
    DummyOutputModulePointer;

/*!\brief
 * Create minimal OutputManager using the provided builder.
 *
 * \param[in] filename Name of file to create OutputManager for.
 * \param[in] dummyTopology Pointer to input top or null.
 * \param[in] adapters Container for the outputadapters.
 * \throws InconsistentInputError When builder can not create the OutputManager.
 * \returns unique_ptr to new OutputManager object.
 */
OutputManagerPointer createMinimalOutputManager(const std::string       &filename,
                                                const gmx_mtop_t        *dummyTopology,
                                                CoordinateOutputAdapters adapters);

/*! \libinternal \brief
 * Test fixture to test matching file types for modules.
 */
class ModuleTest : public gmx::test::CommandLineTestBase,
                   public ::testing::WithParamInterface<const char *>
{
    public:
        /*! \brief
         * Run the builder to create an OutputManager during tests.
         *
         * \param[in] filename Name for output file, to determine filetype during construction.
         * \param[in] adapters Container of outputadapters to add to the OutputManager.
         */
        void runTest(const char *filename, CoordinateOutputAdapters adapters)
        {
            output_ = createMinimalOutputManager(filename,
                                                 &dummyTopology_,
                                                 std::move(adapters));
        }
        //! Pointer to new OutputManager object.
        OutputManagerPointer output_;
        //! Dummy topology to use to create OutputManager.
        gmx_mtop_t           dummyTopology_;
};

} // namespace test

} // namespace gmx

#endif
