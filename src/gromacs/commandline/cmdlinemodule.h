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
 * Declares gmx::CommandLineModuleInterface.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULE_H
#define GMX_COMMANDLINE_CMDLINEMODULE_H

namespace gmx
{

class File;

/*! \brief
 * Module that can be run from command line using CommandLineModuleManager.
 *
 * \see CommandLineModuleManager
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class CommandLineModuleInterface
{
    public:
        virtual ~CommandLineModuleInterface() {}

        //! Returns the name of the module.
        virtual const char *name() const = 0;
        //! Returns a one-line description of the module.
        virtual const char *shortDescription() const = 0;

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
        virtual int run(int argc, char *argv[]) = 0;
        /*! \brief
         * Prints help for the module.
         *
         * \param[in] file  File to write the help to.
         * \throws    std::bad_alloc if out of memory.
         */
        virtual void writeHelp(File *file) const = 0;
};

} // namespace gmx

#endif
