/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULEMANAGER_H
#define GMX_COMMANDLINE_CMDLINEMODULEMANAGER_H

#include "../onlinehelp/helptopicinterface.h"
#include "../utility/common.h"
#include "../utility/uniqueptr.h"

namespace gmx
{

class CommandLineModuleInterface;
class ProgramInfo;

//! Smart pointer type for managing a CommandLineModuleInterface.
typedef gmx_unique_ptr<CommandLineModuleInterface>::type
    CommandLineModulePointer;

/*! \brief
 * Implements a wrapper command-line interface for multiple modules.
 *
 * Typical usage:
 * \code
   int main(int argc, char *argv[])
   {
       const gmx::ProgramInfo &programInfo =
           gmx::ProgramInfo::init("gmx", argc, argv);
       try
       {
           gmx::CommandLineModuleManager manager(programInfo);
           // <register all necessary modules>
           return manager.run(argc, argv);
       }
       catch (const std::exception &ex)
       {
           gmx::printFatalErrorMessage(stderr, ex);
           return 1;
       }
   }
 * \endcode
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class CommandLineModuleManager
{
    public:
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
               return gmx::CommandLineModuleManager::runAsMainSingleModule(argc, argv, &module);
           }
         * \endcode
         *
         * Does not throw.  All exceptions are caught and handled internally.
         */
        static int runAsMainSingleModule(int argc, char *argv[],
                                         CommandLineModuleInterface *module);

        /*! \brief
         * Initializes a command-line module manager.
         *
         * \param     programInfo  Program information for the running binary.
         * \throws    std::bad_alloc if out of memory.
         *
         * The binary name is used to detect when the binary is run through a
         * symlink, and automatically invoke a matching module in such a case.
         *
         * \p programInfo is non-const to allow the manager to amend it based
         * on the actual module that is getting executed.
         */
        explicit CommandLineModuleManager(ProgramInfo *programInfo);
        ~CommandLineModuleManager();

        /*! \brief
         * Sets the module manager to quiet mode: don't print anything.
         *
         * \param[in] bQuiet  Whether the module manager should remain silent.
         *
         * Normally, the module manager prints out some information to stderr
         * before it starts the module and after it finishes.  This removes
         * that output, which is useful in particular for unit tests so that
         * they don't spam stderr.
         */
        void setQuiet(bool bQuiet);

        /*! \brief
         * Adds a given module to this manager.
         *
         * \param   module  Module to add.
         * \throws  std::bad_alloc if out of memory.
         *
         * The manager takes ownership of the object.
         *
         * This method is public mostly for testing purposes; for typical uses,
         * registerModule() is a more convenient way of adding modules.
         *
         * \see registerModule()
         */
        void addModule(CommandLineModulePointer module);
        /*! \brief
         * Registers a module of a certain type to this manager.
         *
         * \tparam  Module  Type of module to register.
         * \throws  std::bad_alloc if out of memory.
         *
         * \p Module must be default-constructible and implement
         * CommandLineModuleInterface.
         *
         * This method is provided as a convenient alternative to addModule()
         * for cases where each module is implemented by a different type
         * (which should be the case for typical situations outside unit
         * tests).
         */
        template <class Module>
        void registerModule()
        {
            addModule(CommandLineModulePointer(new Module));
        }

        /*! \brief
         * Make given help topic available through the manager's help module.
         *
         * \param[in]  topic  Help topic to add.
         * \throws     std::bad_alloc if out of memory.
         *
         * The manager takes ownership of the help topic.
         */
        void addHelpTopic(HelpTopicPointer topic);

        /*! \brief
         * Runs a module based on given command line.
         *
         * \param[in] argc  Number of elements in \p argv.
         * \param[in] argv  Command-line arguments.
         * \throws   unspecified  Throws any exception that the selected module
         *      throws.
         * \returns  Exit code for the program.
         * \retval   0 on successful termination.
         * \retval   2 if no module is specified, or if the module is not found.
         *
         * Runs the module whose name matches \p argv[1].
         */
        int run(int argc, char *argv[]);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
