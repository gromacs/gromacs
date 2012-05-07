/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \file
 * \brief
 * Declares gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULEMANAGER_H
#define GMX_COMMANDLINE_CMDLINEMODULEMANAGER_H

#include "../utility/common.h"
#include "../utility/uniqueptr.h"

namespace gmx
{

class CommandLineModuleInterface;

//! Smart pointer type for managing a CommandLineModuleInterface.
typedef gmx_unique_ptr<CommandLineModuleInterface>::type
        CommandLineModulePointer;

/*! \brief
 * Implements a wrapper command-line interface for multiple modules.
 *
 * Typical usage:
 * \code
int
main(int argc, char *argv[])
{
    CopyRight(stderr, argv[0]);
    try
    {
        gmx::CommandLineModuleManager manager;
        // <register all necessary modules>
        return manager.run(argc, argv);
    }
    catch (const std::exception &ex)
    {
        fprintf(stderr, "%s", gmx::formatErrorMessage(ex).c_str());
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
        CommandLineModuleManager();
        ~CommandLineModuleManager();

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
