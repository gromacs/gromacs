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
#ifndef GMX_TESTUTILS_PROGRAM_CALLER_H
#define GMX_TESTUTILS_PROGRAM_CALLER_H

#include "gromacs/options/options.h"
#include "testutils/cmdlinetest.h"

namespace gmx
{

namespace test
{

class CommandLine;

typedef gmx_unique_ptr<CommandLine>::type CommandLinePointer;

/*! \internal \brief Class for constructing a command line to call
 * a GROMACS tool.
 *
 * An object of this type can be instantiated directly. Any options
 * will then need to be added, with a value for that option already
 * set, because the command-line contribution of that option is
 * constructed immediately.
 *
 * \todo Example usage:
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_testutils
 */
class ProgramCaller
{
    public:
        //! Function pointer type for a C main function.
        // TODO decide whether this is better or worse than #including
        // cmdlinemodulemanager.h
        typedef int (*CMainFunction)(int argc, char *argv[]);
        //! Constructor
        // TODO can registerLegacyModules be used here instead?
        explicit ProgramCaller(const char *programName,
                               CMainFunction cmain);
        //! Copy constructor
        void addOption(const AbstractOption &settings);
        //! Run the program and return its result.
        int run();

    private:
        //! Holds the options that have been added for this program
        gmx::Options  options_;
        //! Holds the command line that will run the program
        CommandLinePointer commandLine_;
        /*! Permit lazy visiting of options to create command line for
         *  the object to use, and block adding options after that. */
        bool bDoneOptionsVisiting_;
        //! Pointer to the GROMACS-tool function to actually run
        CMainFunction cmainToRun_;
};

} // namespace test
} // namespace gmx

#endif
