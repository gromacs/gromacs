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
/*! \internal \file
 * \brief
 * Declares implementation types for gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULEMANAGER_IMPL_H
#define GMX_COMMANDLINE_CMDLINEMODULEMANAGER_IMPL_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

//! \addtogroup module_commandline
//! \{

//! Container type for mapping module names to module objects.
typedef std::map<std::string, CommandLineModulePointer> CommandLineModuleMap;

/*! \internal
 * \brief
 * Internal data for a CommandLineModuleManager module group.
 *
 * This class contains the state of a module group.  CommandLineModuleGroup
 * provides the public interface to construct/alter the state, and
 * CommandLineModuleManager and its associated classes use it for help output.
 */
class CommandLineModuleGroupData
{
    public:
        /*! \brief
         * Shorthand for a list of modules contained in the group.
         *
         * The first element in the contained pair contains the tag
         * (gmx-something) of the module, and the second element contains the
         * description.  The second element is never NULL.
         */
        typedef std::vector<std::pair<std::string, const char *> > ModuleList;

        /*! \brief
         * Constructs an empty module group.
         *
         * \param[in] modules     List of all modules
         *     (used for checking and default descriptions).
         * \param[in] binaryName  Name of the binary containing the modules.
         * \param[in] title       Title of the group.
         *
         * Does not throw.
         */
        CommandLineModuleGroupData(const CommandLineModuleMap &modules,
                                   const char                 *binaryName,
                                   const char                 *title)
            : allModules_(modules), binaryName_(binaryName), title_(title)
        {
        }

        //! Returns the title for the group.
        const char *title() const { return title_; }
        //! Returns the list of modules in the group.
        const ModuleList &modules() const { return modules_; }

        /*! \brief
         * Adds a module to the group.
         *
         * \param[in] name        Name of the module.
         * \param[in] description Description of the module in this group.
         * \throws    std::bad_alloc if out of memory.
         *
         * If \p description is NULL, the description returned by the module is
         * used.
         */
        void addModule(const char *name, const char *description);

    private:
        const CommandLineModuleMap &allModules_;
        const char                 *binaryName_;
        const char                 *title_;
        ModuleList                  modules_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineModuleGroupData);
};

//! Smart pointer type for managing a CommandLineModuleGroup.
typedef gmx_unique_ptr<CommandLineModuleGroupData>::type
    CommandLineModuleGroupDataPointer;
//! Container type for keeping a list of module groups.
typedef std::vector<CommandLineModuleGroupDataPointer>
    CommandLineModuleGroupList;

/*! \internal
 * \brief
 * Encapsulates some handling of common options to the wrapper binary.
 */
class CommandLineCommonOptionsHolder
{
    public:
        CommandLineCommonOptionsHolder();
        ~CommandLineCommonOptionsHolder();

        //! Initializes the common options.
        void initOptions();
        /*! \brief
         * Finishes option parsing.
         *
         * \returns `false` if the wrapper binary should quit without executing
         *     any module.
         */
        bool finishOptions();

        //! Adjust defaults based on module settings.
        void adjustFromSettings(const CommandLineModuleSettings &settings);

        //! Returns the internal Options object.
        Options *options() { return &options_; }
        //! Returns the settings for printing startup information.
        const BinaryInformationSettings &binaryInfoSettings() const
        {
            return binaryInfoSettings_;
        }

        /*! \brief
         * Returns `true` if common options are set such that the wrapper
         * binary should quit, without running the actual module.
         */
        bool shouldIgnoreActualModule() const
        {
            return bHelp_ || bVersion_;
        }
        //! Returns whether common options specify showing help.
        bool shouldShowHelp() const { return bHelp_; }
        //! Returns whether common options specify showing hidden options in help.
        bool shouldShowHidden() const { return bHidden_; }
        //! Returns whether common options specify quiet execution.
        bool shouldBeQuiet() const
        {
            return bQuiet_ && !bVersion_;
        }
        //! Returns whether backups should be made.
        bool shouldBackup() const { return bBackup_; }

        //! Returns the nice level.
        int niceLevel() const { return niceLevel_; }
        //! Returns whether floating-point exception should be enabled
        bool enableFPExceptions() const { return bFpexcept_; }
        //! Returns the debug level.
        int debugLevel() const { return debugLevel_; }

        //! Returns the file to which startup information should be printed.
        FILE *startupInfoFile() const { return (bVersion_ ? stdout : stderr); }

    private:
        Options                      options_;
        //! Settings for what to write in the startup header.
        BinaryInformationSettings    binaryInfoSettings_;
        bool                         bHelp_;
        bool                         bHidden_;
        bool                         bQuiet_;
        bool                         bVersion_;
        bool                         bCopyright_;
        int                          niceLevel_;
        bool                         bBackup_;
        bool                         bFpexcept_;
        int                          debugLevel_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineCommonOptionsHolder);
};

//! \}

} // namespace gmx

#endif
