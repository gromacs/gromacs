/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \defgroup module_commandline Command Line Program Management (commandline)
 * \ingroup group_utilitymodules
 * \brief
 * Provides functionality for managing command line programs.
 *
 * This module provides utility classes and functions for implementing command
 * line programs.  They are mainly used within \Gromacs, but can also be used
 * from external programs if they want to get a similar user experience to
 * \Gromacs tools.
 *
 * The classes exposed from this module can be roughly divided into two groups:
 *
 *  - Helper classes/functions for implementing the %main() function.
 *    See \ref page_usinglibrary for an overview of those available for user
 *    programs.  These are declared in cmdlineinit.h
 *    (gmx::ICommandLineModule is declared in cmdlinemodule.h and
 *    gmx::ICommandLineOptions in cmdlineoptionsmodule.h).
 *    \if libapi
 *
 *    Additionally, for internal \Gromacs use, gmx::CommandLineModuleManager
 *    provides the functionality to implement the `gmx` wrapper binary, as well
 *    as command line options common to all \Gromacs programs (such as
 *    `-version`).
 *    This is described in more detail at \ref page_wrapperbinary.
 *    \endif
 *
 *  - Helper classes for particular command line tasks:
 *     - gmx::CommandLineParser implements command line parsing to assign
 *       values to gmx::Options (see \ref module_options).
 *     - gmx::CommandLineHelpWriter writes help text for a program that uses
 *       the parser.
 *     - parse_common_args() is an old interface to \Gromacs command line
 *       parsing.  This is still used by many parts of \Gromacs.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Declares gmx::ICommandLineModule and supporting classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULE_H
#define GMX_COMMANDLINE_CMDLINEMODULE_H

namespace gmx
{

class CommandLineHelpContext;
class CommandLineModuleSettings;

/*! \brief
 * Module that can be run from command line using CommandLineModuleManager.
 *
 * \see CommandLineModuleManager
 * \see CommandLineOptionsModule
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class ICommandLineModule
{
public:
    virtual ~ICommandLineModule() {}

    //! Returns the name of the module.
    virtual const char* name() const = 0;
    //! Returns a one-line description of the module.
    virtual const char* shortDescription() const = 0;

    /*! \brief
     * Initializes the module and provides settings for the runner.
     *
     * This will be called before run(), and can be used to adjust
     * initialization that the runner does.
     *
     * This method is currently not called when writing the help.
     */
    virtual void init(CommandLineModuleSettings* settings) = 0;
    /*! \brief
     * Runs the module with the given arguments.
     *
     * \param[in] argc  Number of elements in \p argv.
     * \param[in] argv  Command-line arguments.
     * \throws   unspecified  May throw exceptions to indicate errors.
     * \returns  Exit code for the program.
     * \retval   0 on successful termination.
     *
     * \p argv[0] is the name of the module, i.e., the arguments are as if
     * the module was run as a standalone executable.
     */
    virtual int run(int argc, char* argv[]) = 0;
    /*! \brief
     * Prints help for the module.
     *
     * \param[in] context  Context object for writing the help.
     * \throws    std::bad_alloc if out of memory.
     * \throws    FileIOError on any I/O error.
     *
     * Note that for MPI-enabled builds, this is called only on the main
     * rank.
     */
    virtual void writeHelp(const CommandLineHelpContext& context) const = 0;
};

//! \cond libapi
/*! \libinternal \brief
 * Helper to implement ICommandLineModule::writeHelp() with a C-like
 * main() function that calls parse_common_args().
 *
 * \param[in] context      Context object for writing the help.
 * \param[in] name         Name of the module.
 * \param[in] mainFunction C-like main() function that calls parse_common_args().
 *
 * \ingroup module_commandline
 */
void writeCommandLineHelpCMain(const CommandLineHelpContext& context,
                               const char*                   name,
                               int (*mainFunction)(int argc, char* argv[]));
//! \endcond

} // namespace gmx

#endif
