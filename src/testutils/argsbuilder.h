/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Declares test fixture for mdrun tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_INTEGRATION_TESTS_ARGSBUILDER_H
#define GMX_INTEGRATION_TESTS_ARGSBUILDER_H

#include <string>
#include <vector>
#include "gromacs/options/optionsvisitor.h"
#include "testutils/cmdlinetest.h"

namespace gmx
{

namespace test
{

/*! \internal \brief This class constructs a C-style command line from a gmx::Options object.
 *
 * It does so by implementing the Visitor pattern, by inheriting the
 * gmx::OptionsVisitor interface. This means it can traverse a
 * gmx::Options data structure and build a CommandLine.
 *
 * Does not throw.
 *
 * \ingroup module_testutils
 */

class ArgsBuilder : public gmx::OptionsVisitor
{
    public:
        //! Construct an ArgsBuilder for the named program.
        ArgsBuilder(std::string const &programName);
        //! Clean up after the GROMACS tool's main().
        ~ArgsBuilder();
        //! Visit any (sub)sections of the gmx::Options
        void visitSubSection(const Options &section);
        //! Visit an option and pass it to commandLine_
        void visitOption(const OptionInfo &option);
        //! Returns argv for passing into C-style command-line handling.
        int getArgc();
        //! Returns argv for passing into C-style command-line handling.
        char **getArgv();

    private:
        /*! /brief Helps build a command line
         *
         * This gets filled during the visiting process.
         */
        gmx::test::CommandLine commandLine_;
};

} // namespace test
} // namespace gmx

#endif
