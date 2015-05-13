/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_MDRUN_TESTS_MULTISIMTEST_H
#define GMX_MDRUN_TESTS_MULTISIMTEST_H

/*! \internal \file
 * \brief
 * Declares test fixture for the mdrun multi-simulation functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/uniqueptr.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

//! Convenience typedef
typedef gmx_unique_ptr<CommandLine>::type CommandLinePointer;

class MultiSimTest : public gmx::test::ParameterizedMdrunTestFixture
{
    public:
        //! Constructor
        MultiSimTest();

        /*! \brief Organize the .mdp file for this rank
         *
         * For testing multi-simulation, this .mdp file is more
         * complicated than it needs to be, but it does little harm,
         * and doing it this way allows this function to be re-used
         * for testing replica-exchange.
         *
         * \param controlVariable Allows parameterization to work with
         * T, P or (later) lambda as the control variable, by passing a
         * string with "mdp-param = value" such that different paths
         * in init_replica_exchange() are followed.
         */
        void organizeMdpFile(const char *controlVariable);

        //! Number of MPI ranks
        int                size_;
        //! MPI rank of this process
        int                rank_;
        //! Object for building the mdrun command line
        CommandLinePointer mdrunCaller_;
        //! Name of .tpr file to be used by mdrun
        std::string        mdrunTprFileName_;
};

} // namespace
} // namespace

#endif
