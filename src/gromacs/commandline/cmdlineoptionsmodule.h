/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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
 * Declares gmx::ICommandLineOptionsModule and supporting routines.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEOPTIONSMODULE_H
#define GMX_COMMANDLINE_CMDLINEOPTIONSMODULE_H

#include <functional>
#include <memory>

#include "gromacs/commandline/cmdlinemodule.h"

namespace gmx
{

template <typename T> class ArrayRef;

class CommandLineModuleManager;
class ICommandLineModule;
class ICommandLineOptionsModule;
class IOptionsBehavior;
class IOptionsContainer;

//! Smart pointer to manage an ICommandLineOptionsModule.
typedef std::unique_ptr<ICommandLineOptionsModule>
    ICommandLineOptionsModulePointer;

/*! \brief
 * Settings to pass information between a CommandLineOptionsModule and generic
 * code that runs it.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class ICommandLineOptionsModuleSettings
{
    public:
        /*! \brief
         * Sets the help text for the module from string array.
         *
         * \param[in] help  String array to set as the description.
         * \throws    std::bad_alloc if out of memory.
         *
         * Formatting for the help text is described on \ref page_onlinehelp.
         *
         * Example usage:
         * \code
           const char *const desc[] = {
               "This is the description",
               "for the options"
           };

           settings->setHelpText(desc);
           \endcode
         */
        virtual void setHelpText(const ArrayRef<const char *const> &help) = 0;
        /*! \brief
         * Adds an option behavior that performs actions before
         * ICommandLineOptionsModule::run() is called.
         *
         * For now, this takes a shared_ptr to make it easier for the caller to
         * keep a reference to the behavior, but the behavior should be treated
         * as owned by the options module after this call.
         */
        virtual void addOptionsBehavior(
            const std::shared_ptr<IOptionsBehavior> &behavior) = 0;

    protected:
        // Disallow deletion through the interface.
        // (no need for the virtual, but some compilers warn otherwise)
        virtual ~ICommandLineOptionsModuleSettings();
};

/*! \brief
 * Module that can be run from a command line and uses gmx::Options for
 * argument processing.
 *
 * This class provides a higher-level interface on top of
 * gmx::ICommandLineModule for cases where gmx::Options will be used
 * for declaring the command-line arguments.  The module only needs to declare
 * the options it uses, and the framework takes care of command-line parsing
 * and help output.  The module typically consists of the following parts:
 *  - init() allows for some interaction between the module and the framework
 *    when running the module; see ICommandLineModule::init().  If no
 *    such customization is necessary, an empty implementation is sufficient.
 *  - initOptions() is called both for running the module and for printing help
 *    for the module, and it should add the options that the module
 *    understands.  Values provided for the options are typically stored in
 *    member variables.
 *  - optionsFinished() can be implemented in case additional processing is
 *    needed (e.g., checking whether an option was set by the user).
 *  - run() is called when running the module, after command-line options have
 *    been parsed and their values stored in the corresponding member
 *    variables.
 *
 * registerModule(), runAsMain(), or createModule() can be used to use modules
 * of this type in all contexts where a gmx::ICommandLineModule is
 * expected.  These methods create a gmx::ICommandLineModule
 * implementation that contains the common code needed to parse command-line
 * options and write help, based on the information provided from the methods
 * in this class.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class ICommandLineOptionsModule
{
    public:
        /*! \brief
         * Function pointer to a factory method that returns an interface of
         * this type.
         *
         * \returns Module to run.
         * \throws  std::bad_alloc if out of memory.
         */
        typedef std::function<ICommandLineOptionsModulePointer()> FactoryMethod;

        /*! \brief
         * Creates a ICommandLineModule to run the specified module.
         *
         * \param[in] name        Name for the module.
         * \param[in] description Short description for the module.
         * \param[in] module      Module to run.
         * \returns ICommandLineModule object that runs \p module module.
         * \throws  std::bad_alloc if out of memory.
         */
        static std::unique_ptr<ICommandLineModule>
        createModule(const char *name, const char *description,
                     ICommandLineOptionsModulePointer module);
        /*! \brief
         * Implements a main() method that runs a single module.
         *
         * \param     argc    \c argc passed to main().
         * \param     argv    \c argv passed to main().
         * \param[in] name        Name for the module.
         * \param[in] description Short description for the module.
         * \param[in] factory     Factory that returns the module to run.
         *
         * This method allows for uniform behavior for binaries that only
         * contain a single module without duplicating any of the
         * implementation from CommandLineModuleManager (startup headers,
         * common options etc.).
         *
         * \see runCommandLineModule()
         */
        static int
        runAsMain(int argc, char *argv[], const char *name,
                  const char *description, FactoryMethod factory);
        /*! \brief
         * Registers a module of a certain type to this manager.
         *
         * \param     manager     Manager to register to.
         * \param[in] name        Name for the module.
         * \param[in] description Short description for the module.
         * \param[in] factory     Factory that returns the module to register.
         * \throws  std::bad_alloc if out of memory.
         *
         * This method internally creates a ICommandLineModule module
         * with the given \p name and \p description, and adds that to
         * \p manager.  When run or asked to write the help, the module calls
         * \p factory to get the actual module, and forwards the necessary
         * calls.
         */
        static void
        registerModuleFactory(CommandLineModuleManager *manager,
                              const char *name, const char *description,
                              FactoryMethod factory);
        /*! \brief
         * Registers a module to this manager.
         *
         * \param     manager     Manager to register to.
         * \param[in] name        Name for the module.
         * \param[in] description Short description for the module.
         * \param[in] module      Module to register.
         * \throws  std::bad_alloc if out of memory.
         *
         * This method internally creates a ICommandLineModule module
         * with the given \p name and \p description, and adds that to
         * \p manager.
         *
         * This method is mainly used by tests that need to have a reference to
         * the ICommandLineOptionsModule instance (e.g., for mocking).
         */
        static void
        registerModuleDirect(CommandLineModuleManager *manager,
                             const char *name, const char *description,
                             ICommandLineOptionsModulePointer module);

        virtual ~ICommandLineOptionsModule();

        //! \copydoc gmx::ICommandLineModule::init()
        virtual void init(CommandLineModuleSettings *settings) = 0;
        /*! \brief
         * Initializes command-line arguments understood by the module.
         *
         * \param[in,out] options  Options object to add the options to.
         * \param[in,out] settings Settings to communicate information
         *     to/from generic code running the module.
         *
         * When running the module, this method is called after init().
         * When printing help, there is no call to init(), and this is the only
         * method called.
         * In both cases, the implementation should add options understood by
         * the module to \p options.  Output values from options should be
         * stored in member variables.
         */
        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings) = 0;
        /*! \brief
         * Called after all option values have been set.
         *
         * When running the module, this method is called after all
         * command-line arguments have been parsed.
         */
        virtual void optionsFinished() = 0;

        /*! \brief
         * Runs the module.
         *
         * \throws   unspecified  May throw exceptions to indicate errors.
         * \returns  Exit code for the program.
         * \retval   0 on successful termination.
         *
         * This method is called after optionsFinished() when running the
         * module, and should do all the processing for the module.
         */
        virtual int run() = 0;
};

} // namespace gmx

#endif
