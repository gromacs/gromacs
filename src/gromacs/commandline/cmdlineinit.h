/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares functions for initializing the \Gromacs library for command line use.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEINIT_H
#define GMX_COMMANDLINE_CMDLINEINIT_H

#ifdef __cplusplus

// Forward declaration of class CommandLineProgramContext is not sufficient for
// MSVC if the return value of initForCommandLine() is ignored(!)
#include "gromacs/commandline/cmdlineprogramcontext.h"

namespace gmx
{

class CommandLineModuleInterface;
class CommandLineOptionsModuleInterface;

/*! \brief
 * Initializes the \Gromacs library for command-line use.
 *
 * \param[in] argc  argc value passed to main().
 * \param[in] argv  argv array passed to main().
 * \returns   Reference to initialized program context object.
 *
 * This function is tailored for use in command line applications.
 * For other usage, combination of gmx::init() and gmx::setProgramContext()
 * provides more flexible initialization alternatives.
 * Unlike gmx::init(), calls to this method cannot be nested.
 *
 * The command line arguments are communicated so that they can be
 * parsed on each processor.
 * \p argc and \p argv are passed to gmx::init(); see there for additional
 * discussion.  This method does not place any additional limitations, but
 * generally there should be no need to pass NULL values.
 *
 * Does not throw. Terminates the program on out-of-memory error.
 *
 * This method is not thread-safe, since it is intended to be the first method
 * called.  See setProgramContext() for additional discussion.
 *
 * \see gmx::init()
 * \see setProgramContext()
 * \ingroup module_commandline
 */
CommandLineProgramContext &initForCommandLine(int *argc, char ***argv);
/*! \brief
 * Deinitializes the \Gromacs library after initForCommandLine().
 *
 * Calls gmx::finalize() and additionally undoes the work done by
 * initForCommandLine().
 *
 * \see gmx::finalize()
 * \ingroup module_commandline
 */
void finalizeForCommandLine();
/*! \brief
 * Handles an exception and deinitializes after initForCommandLine.
 *
 * \param[in] ex  Exception that is the cause for terminating the program.
 * \returns   Return code to return from main().
 *
 * This method should be called as the last thing before terminating the
 * program because of an exception. See processExceptionAtExit() for details.
 * Additionally this method undoes the work done by initForCommandLine.
 *
 * Does not throw.
 */
int processExceptionAtExitForCommandLine(const std::exception &ex);
/*! \brief
 * Implements a main() method that runs a single module.
 *
 * \param argc   \c argc passed to main().
 * \param argv   \c argv passed to main().
 * \param module Module to run.
 *
 * This method allows for uniform behavior for binaries that only
 * contain a single module without duplicating any of the
 * implementation from CommandLineModuleManager (startup headers,
 * common options etc.).
 *
 * The signature assumes that \p module construction does not throw
 * (because otherwise the caller would need to duplicate all the
 * exception handling code).  It is possible to move the construction
 * inside the try/catch in this method using an indirection similar to
 * TrajectoryAnalysisCommandLineRunner::runAsMain(), but until that is
 * necessary, the current approach leads to simpler code.
 *
 * Usage:
 * \code
   int main(int argc, char *argv[])
   {
       CustomCommandLineModule module;
       return gmx::runCommandLineModule(argc, argv, &module);
   }
   \endcode
 *
 * Does not throw.  All exceptions are caught and handled internally.
 */
int runCommandLineModule(int argc, char *argv[],
                         CommandLineModuleInterface *module);
/*! \brief
 * Implements a main() method that runs a single module.
 *
 * \param     argc        \c argc passed to main().
 * \param     argv        \c argv passed to main().
 * \param[in] name        Name for the module.
 * \param[in] description Short description for the module.
 * \param     factory Factory method that creates the module to run.
 *
 * This method allows for uniform behavior for binaries that only
 * contain a single module without duplicating any of the
 * implementation from CommandLineModuleManager (startup headers,
 * common options etc.).
 *
 * Usage:
 * \code
   class CustomCommandLineOptionsModule : public CommandLineOptionsModuleInterface
   {
       // <...>
   };

   static CommandLineOptionsModuleInterface *create()
   {
       return new CustomCommandLineOptionsModule();
   }

   int main(int argc, char *argv[])
   {
       return gmx::runCommandLineModule(
               argc, argv, "mymodule", "short description", &create);
   }
   \endcode
 *
 * Does not throw.  All exceptions are caught and handled internally.
 */
int runCommandLineModule(int argc, char *argv[],
                         const char *name, const char *description,
                         CommandLineOptionsModuleInterface *(*factory)());

} // namespace gmx

extern "C"
{
#endif

/*! \brief
 * Implements a main() method that runs a given C main function.
 *
 * \param argc         \c argc passed to main().
 * \param argv         \c argv passed to main().
 * \param mainFunction The main()-like method to wrap.
 *
 * This method creates a dummy command line module that does its
 * processing by calling \p mainFunction.  It then runs this module as with
 * gmx::runCommandLineModule().
 * This allows the resulting executable to handle common options and do
 * other common actions (e.g., startup headers) without duplicate code
 * in the main methods.
 *
 * \p mainFunction should call parse_common_args() to process its command-line
 * arguments.
 *
 * Usage:
 * \code
   int my_main(int argc, char *argv[])
   {
       // <...>
   }

   int main(int argc, char *argv[])
   {
       return gmx_run_cmain(argc, argv, &my_main);
   }
   \endcode
 *
 * Does not throw.  All exceptions are caught and handled internally.
 */
int gmx_run_cmain(int argc, char *argv[], int (*mainFunction)(int, char *[]));

#ifdef __cplusplus
}
#endif

#endif
