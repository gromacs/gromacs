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

namespace gmx
{

namespace test
{

/*! \internal \brief Class for building a command-line-style argv
 * array to pass to a GROMACS tool.
 *
 * This class implements the Visitor pattern by inheriting the
 * gmx::OptionsVisitor interface. This means it can traverse a
 * gmx::Options data structure and build an argv array from the
 * options that have been set to real values, while ignoring the
 * options that are either missing, or defined and not set.
 *
 * It also handles the fact that parse_file_args() used in all the
 * GROMACS legacy tools feels the need to modify argv, and has done so
 * since the original GROMACS CVS commit. (Presumably nobody has any
 * idea whether there's a sound reason for that.) That fact prevents
 * us from calling a GROMACS main() and then freeing memory in the
 * argv we passed. So, we have to keep copies of the pointers that
 * were in argv before they were passed into the main(), so we can
 * free them once main() returns.
 *
 * // TODO think about this
 * //  Any method in this class may throw std::bad_alloc if out of
 * // memory.
 *
 * \ingroup module_integration_tests
 */
class ArgsBuilder : public gmx::OptionsVisitor
{
public:
    /*! /brief Construct an ArgsBuilder for the named program. */
    ArgsBuilder(const char* programName);
    /*! /brief Clean up after the GROMACS tool's main(). */
    ~ArgsBuilder();
    /*! /brief Visit any (sub)sections of the gmx::Options */
    void visitSubSection(const Options &section);
    /*! /brief Visit an option and convert it to argv arry entries. */
    void visitOption(const OptionInfo &option);
    /*! /brief Return the length of the argv array, such as can be
     *  used for the argc argument to a main() function. */
    int getArgc();
    /*! /brief Return a C-style argv array, such as can be used for
     *  the argv argument to a main() function. */
    char **getArgv();

private:
    /*! /brief Container for the arguments that will be copied to
     *  argv.
     *
     * This gets filled during the visiting process.
     */
    std::vector<std::string> arguments;
    /*! /brief The vector that will actually be passed to a GROMACS
     *  tool's main(). */
    std::vector<char *> theArgv;
    /*! /brief Hold handles to dynamically allocated C strings that
     *  will be passed to a GROMACS legacy tool, so that we can free
     *  them after the tool's main() returns. */
    std::vector<char *> theArgvHolder;
};
 
} // namespace test
} // namespace gmx

#endif
