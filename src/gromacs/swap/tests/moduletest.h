/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
/*! \defgroup module_mdrun_integration_tests Integration test utilities
 * \ingroup group_mdrun
 *
 * \brief Functionality for testing computational elecrophysiology setups
 * in mdrun
 */

#include "testutils/integrationtests.h"

#include <gtest/gtest.h>

namespace gmx
{

namespace test
{

class SwapTestFixture : public IntegrationTestFixture
{
    protected:
        SwapTestFixture();
        ~SwapTestFixture();

        //! Calls grompp to prepare for the mdrun test
        int callGrompp();
        //! Calls mdrun for testing
        int callMdrun();
        //! Calls mdrun and expects to read from checkpoint
        int callMdrunAppend();


        //@{
        /*! \name Names for grompp and mdrun output files
         *
         * These strings can be set to point to files present in the
         * source tree, or to temporary files created for the test
         * fixture. In the latter case,
         * IntegrationTestFixture::fileManager_ should be used to fill
         * these strings with paths to files, so that they are created
         * in a temporary directory and (by default behaviour of
         * TestFileManager) deleted when the test is complete.
         */
        std::string topFileName;
        std::string groInputFileName;
        std::string xtcFileName;
        std::string mdpInputFileName;
        std::string mdpOutputFileName;
        std::string tprFileName;
        std::string logFileName;
        std::string ndxFileName;
        std::string edrFileName;
        std::string cptFileName;
        std::string swapFileName;
        std::string groOutputFileName;
        //@}
};

} // namespace test
} // namespace gmx
