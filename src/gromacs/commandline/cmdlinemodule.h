/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Declares gmx::CommandLineModuleInterface and supporting classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULE_H
#define GMX_COMMANDLINE_CMDLINEMODULE_H

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class CommandLineHelpContext;

/*! \brief
 * Settings to pass information between a module and the general runner.
 *
 * Methods in this class do not throw, except that construction may throw
 * std::bad_alloc.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class CommandLineModuleSettings
{
    public:
        CommandLineModuleSettings();
        ~CommandLineModuleSettings();

        //! Returns the default nice level for this module.
        int defaultNiceLevel() const;

        /*! \brief
         * Sets the default nice level for this module.
         *
         * If not called, the module will be niced.
         */
        void setDefaultNiceLevel(int niceLevel);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \brief
 * Module that can be run from command line using CommandLineModuleManager.
 *
 * \see CommandLineModuleManager
 * \see CommandLineOptionsModule
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
         * Initializes the module and provides settings for the runner.
         *
         * This will be called before run(), and can be used to adjust
         * initialization that the runner does.
         *
         * This method is currently not called when writing the help.
         */
        virtual void init(CommandLineModuleSettings *settings) = 0;
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
         * \param[in] context  Context object for writing the help.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Note that for MPI-enabled builds, this is called only on the master
         * rank.
         */
        virtual void writeHelp(const CommandLineHelpContext &context) const = 0;
};

//! \cond libapi
/*! \libinternal \brief
 * Helper to implement CommandLineModuleInterface::writeHelp() with a C-like
 * main() function that calls parse_common_args().
 *
 * \param[in] context      Context object for writing the help.
 * \param[in] name         Name of the module.
 * \param[in] mainFunction C-like main() function that calls parse_common_args().
 *
 * \ingroup module_commandline
 */
void writeCommandLineHelpCMain(
        const CommandLineHelpContext &context, const char *name,
        int (*mainFunction)(int argc, char *argv[]));
//! \endcond

} // namespace gmx

#endif
