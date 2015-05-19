/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares gmx::CommandLineHelpModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEHELPMODULE_H
#define GMX_COMMANDLINE_CMDLINEHELPMODULE_H

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/onlinehelp/helptopicinterface.h"
#include "gromacs/utility/classhelpers.h"

#include "cmdlinemodulemanager-impl.h"

namespace gmx
{

class CommandLineHelpContext;
class FileOutputRedirectorInterface;
class ProgramContextInterface;

class CommandLineHelpModuleImpl;

/*! \internal
 * \brief
 * Command-line module for producing help.
 *
 * This module implements the 'help' subcommand that is automatically added by
 * CommandLineModuleManager.
 *
 * \ingroup module_commandline
 */
class CommandLineHelpModule : public CommandLineModuleInterface
{
    public:
        /*! \brief
         * Creates a command-line help module.
         *
         * \param[in] programContext Information about the running binary.
         * \param[in] binaryName     Name of the running binary
         *     (without Gromacs binary suffix or .exe on Windows).
         * \param[in] modules  List of modules for to use for module listings.
         * \param[in] groups   List of module groups.
         * \throws    std::bad_alloc if out of memory.
         */
        CommandLineHelpModule(const ProgramContextInterface    &programContext,
                              const std::string                &binaryName,
                              const CommandLineModuleMap       &modules,
                              const CommandLineModuleGroupList &groups);
        ~CommandLineHelpModule();

        /*! \brief
         * Creates a help topic for a command-line module.
         *
         * \param[in] module  Module the create the help topic for.
         * \throws    std::bad_alloc if out of memory.
         *
         * The caller should add the topic using addTopic() if that is desired.
         * This method is provided separately to allow for strong exception
         * safety in CommandLineModuleManager::addModule().
         */
        HelpTopicPointer
        createModuleHelpTopic(const CommandLineModuleInterface &module) const;
        /*! \brief
         * Adds a top-level help topic.
         *
         * \param[in] topic     Help topic to add.
         * \param[in] bExported Whether this topic will be directly exported to
         *     the user guide.
         * \throws    std::bad_alloc if out of memory.
         */
        void addTopic(HelpTopicPointer topic, bool bExported);
        //! Sets whether hidden options will be shown in help.
        void setShowHidden(bool bHidden);
        /*! \brief
         * Sets an override to show the help for the given module.
         *
         * If called, the help module directly prints the help for the given
         * module when called, skipping any other processing.
         */
        void setModuleOverride(const CommandLineModuleInterface &module);

        /*! \brief
         * Sets a file redirector for writing help output.
         *
         * Used for unit testing; see
         * CommandLineModuleManager::setOutputRedirector() for more details.
         */
        void setOutputRedirector(FileOutputRedirectorInterface *output);

        virtual const char *name() const { return "help"; }
        virtual const char *shortDescription() const
        {
            return "Print help information";
        }

        virtual void init(CommandLineModuleSettings *settings)
        {
            settings->setDefaultNiceLevel(0);
        }
        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const CommandLineHelpContext &context) const;

    private:
        typedef CommandLineHelpModuleImpl Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
