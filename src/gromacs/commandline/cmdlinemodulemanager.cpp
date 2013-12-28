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
 * Implements gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "cmdlinemodulemanager.h"

#include <cstdio>

#include <string>
#include <utility>

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/network.h"

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpmodule.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager-impl.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/programinfo.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/********************************************************************
 * CMainCommandLineModule
 */

/*! \internal \brief
 * Implements a CommandLineModuleInterface, given a function with C/C++ main()
 * signature.
 *
 * \ingroup module_commandline
 */
class CMainCommandLineModule : public CommandLineModuleInterface
{
    public:
        //! \copydoc gmx::CommandLineModuleManager::CMainFunction
        typedef CommandLineModuleManager::CMainFunction CMainFunction;

        /*! \brief
         * Creates a wrapper module for the given main function.
         *
         * \param[in] name             Name for the module.
         * \param[in] shortDescription One-line description for the module.
         * \param[in] mainFunction     Main function to wrap.
         *
         * Does not throw.  This is essential for correct implementation of
         * CommandLineModuleManager::runAsMainCMain().
         */
        CMainCommandLineModule(const char *name, const char *shortDescription,
                               CMainFunction mainFunction)
            : name_(name), shortDescription_(shortDescription),
              mainFunction_(mainFunction)
        {
        }

        virtual const char *name() const
        {
            return name_;
        }
        virtual const char *shortDescription() const
        {
            return shortDescription_;
        }

        virtual int run(int argc, char *argv[])
        {
            return mainFunction_(argc, argv);
        }
        virtual void writeHelp(const CommandLineHelpContext &context) const
        {
            char *argv[2];
            int   argc = 1;
            // TODO: The constness should not be cast away.
            argv[0] = const_cast<char *>(name_);
            argv[1] = NULL;
            GlobalCommandLineHelpContext global(context);
            mainFunction_(argc, argv);
        }

    private:
        const char             *name_;
        const char             *shortDescription_;
        CMainFunction           mainFunction_;

};

}   // namespace

/********************************************************************
 * CommandLineModuleManager::Impl
 */

/*! \internal \brief
 * Private implementation class for CommandLineModuleManager.
 *
 * \ingroup module_commandline
 */
class CommandLineModuleManager::Impl
{
    public:
        /*! \brief
         * Initializes the implementation class.
         *
         * \param     programInfo  Program information for the running binary.
         */
        explicit Impl(ProgramInfo *programInfo);

        /*! \brief
         * Helper method that adds a given module to the module manager.
         *
         * \throws    std::bad_alloc if out of memory.
         */
        void addModule(CommandLineModulePointer module);
        /*! \brief
         * Creates the help module if it does not yet exist.
         *
         * \throws    std::bad_alloc if out of memory.
         *
         * This method should be called before accessing \a helpModule_.
         */
        void ensureHelpModuleExists();

        /*! \brief
         * Finds a module that matches a name.
         *
         * \param[in] name  Module name to find.
         * \returns   Iterator to the found module, or
         *      \c modules_.end() if not found.
         *
         * Does not throw.
         */
        CommandLineModuleMap::const_iterator
        findModuleByName(const std::string &name) const;
        /*! \brief
         * Finds a module that the name of the binary.
         *
         * \param[in] programInfo  Program information object to use.
         * \throws    std::bad_alloc if out of memory.
         * \returns   Iterator to the found module, or
         *      \c modules_.end() if not found.
         *
         * Checks whether the program is invoked through a symlink whose name
         * is different from ProgramInfo::realBinaryName(), and if so, checks
         * if a module name matches the name of the symlink.
         *
         * Note that the \p programInfo parameter is currently not necessary
         * (as the program info object is also contained as a member), but it
         * clarifies the control flow.
         */
        CommandLineModuleMap::const_iterator
        findModuleFromBinaryName(const ProgramInfo &programInfo) const;

        /*! \brief
         * Processes command-line options for the wrapper binary.
         *
         * \param[in,out] argc On input, argc passed to run().
         *     On output, argc to be passed to the module.
         * \param[in,out] argv On input, argv passed to run().
         *     On output, argv to be passed to the module.
         * \throws    InvalidInputError if there are invalid options.
         * \returns   The module that should be run.
         *
         * Handles command-line options that affect the wrapper binary
         * (potentially changing the members of \c this in response to the
         * options).  Also finds the module that should be run and the
         * arguments that should be passed to it.
         */
        CommandLineModuleInterface *
        processCommonOptions(int *argc, char ***argv);

        /*! \brief
         * Maps module names to module objects.
         *
         * Owns the contained modules.
         */
        CommandLineModuleMap         modules_;
        /*! \brief
         * List of groupings for modules for help output.
         *
         * Owns the contained module group data objects.
         * CommandLineModuleGroup objects point to the data objects contained
         * here.
         */
        CommandLineModuleGroupList   moduleGroups_;
        //! Information about the currently running program.
        ProgramInfo                 &programInfo_;
        /*! \brief
         * Module that implements help for the binary.
         *
         * The pointed module is owned by the \a modules_ container.
         */
        CommandLineHelpModule       *helpModule_;
        //! Settings for what to write in the startup header.
        BinaryInformationSettings    binaryInfoSettings_;
        //! If non-NULL, run this module in single-module mode.
        CommandLineModuleInterface  *singleModule_;
        //! Whether all stderr output should be suppressed.
        bool                         bQuiet_;
        //! Whether to write the startup information to stdout iso stderr.
        bool                         bStdOutInfo_;

    private:
        GMX_DISALLOW_COPY_AND_ASSIGN(Impl);
};

CommandLineModuleManager::Impl::Impl(ProgramInfo *programInfo)
    : programInfo_(*programInfo), helpModule_(NULL), singleModule_(NULL),
      bQuiet_(false), bStdOutInfo_(false)
{
    binaryInfoSettings_.copyright(true);
}

void CommandLineModuleManager::Impl::addModule(CommandLineModulePointer module)
{
    GMX_ASSERT(modules_.find(module->name()) == modules_.end(),
               "Attempted to register a duplicate module name");
    ensureHelpModuleExists();
    HelpTopicPointer helpTopic(helpModule_->createModuleHelpTopic(*module));
    modules_.insert(std::make_pair(std::string(module->name()),
                                   move(module)));
    helpModule_->addTopic(move(helpTopic));
}

void CommandLineModuleManager::Impl::ensureHelpModuleExists()
{
    if (helpModule_ == NULL)
    {
        helpModule_ = new CommandLineHelpModule(programInfo_, modules_,
                                                moduleGroups_);
        addModule(CommandLineModulePointer(helpModule_));
    }
}

CommandLineModuleMap::const_iterator
CommandLineModuleManager::Impl::findModuleByName(const std::string &name) const
{
    // TODO: Accept unambiguous prefixes?
    return modules_.find(name);
}

CommandLineModuleMap::const_iterator
CommandLineModuleManager::Impl::findModuleFromBinaryName(
        const ProgramInfo &programInfo) const
{
    std::string binaryName = programInfo.invariantProgramName();
    if (binaryName == programInfo.realBinaryName())
    {
        return modules_.end();
    }
    if (binaryName.compare(0, 2, "g_") == 0)
    {
        binaryName.erase(0, 2);
    }
    if (binaryName.compare(0, 3, "gmx") == 0)
    {
        binaryName.erase(0, 3);
    }
    return findModuleByName(binaryName);
}

CommandLineModuleInterface *
CommandLineModuleManager::Impl::processCommonOptions(int *argc, char ***argv)
{
    // Check if we are directly invoking a certain module.
    CommandLineModuleInterface *module = singleModule_;
    if (module == NULL)
    {
        // Also check for invokation through named symlinks.
        CommandLineModuleMap::const_iterator moduleIter
            = findModuleFromBinaryName(programInfo_);
        if (moduleIter != modules_.end())
        {
            module = moduleIter->second.get();
        }
    }

    bool bHelp      = false;
    bool bHidden    = false;
    bool bVersion   = false;
    bool bCopyright = true;
    // TODO: Print the common options into the help.
    // TODO: It would be nice to propagate at least the -quiet option to
    // the modules so that they can also be quiet in response to this.
    Options options(NULL, NULL);
    options.addOption(BooleanOption("h").store(&bHelp));
    options.addOption(BooleanOption("hidden").store(&bHidden));
    options.addOption(BooleanOption("quiet").store(&bQuiet_));
    options.addOption(BooleanOption("version").store(&bVersion));
    options.addOption(BooleanOption("copyright").store(&bCopyright));

    if (module == NULL)
    {
        // If not in single-module mode, process options to the wrapper binary.
        // TODO: Ideally, this could be done by CommandLineParser.
        int argcForWrapper = 1;
        while (argcForWrapper < *argc && (*argv)[argcForWrapper][0] == '-')
        {
            ++argcForWrapper;
        }
        if (argcForWrapper > 1)
        {
            CommandLineParser(&options).parse(&argcForWrapper, *argv);
        }
        // If no action requested and there is a module specified, process it.
        if (argcForWrapper < *argc && !bHelp && !bVersion)
        {
            const char *moduleName = (*argv)[argcForWrapper];
            CommandLineModuleMap::const_iterator moduleIter
                = findModuleByName(moduleName);
            if (moduleIter == modules_.end())
            {
                std::string message =
                    formatString("'%s' is not a GROMACS command.", moduleName);
                GMX_THROW(InvalidInputError(message));
            }
            module = moduleIter->second.get();
            *argc -= argcForWrapper;
            *argv += argcForWrapper;
            // After this point, argc and argv are the same independent of
            // which path is taken: (*argv)[0] is the module name.
        }
    }
    if (module != NULL)
    {
        if (singleModule_ == NULL)
        {
            programInfo_.setDisplayName(
                    programInfo_.realBinaryName() + " " + module->name());
        }
        // Recognize the common options also after the module name.
        // TODO: It could be nicer to only recognize -h/-hidden if module is not
        // null.
        CommandLineParser(&options).skipUnknown(true).parse(argc, *argv);
    }
    options.finish();
    binaryInfoSettings_.extendedInfo(bVersion);
    binaryInfoSettings_.copyright(bCopyright);
    if (bVersion)
    {
        bQuiet_      = false;
        bStdOutInfo_ = true;
        return NULL;
    }
    // If no module specified and no other action, show the help.
    // Also explicitly specifying -h for the wrapper binary goes here.
    if (module == NULL || bHelp)
    {
        ensureHelpModuleExists();
        if (module != NULL)
        {
            helpModule_->setModuleOverride(*module);
        }
        *argc  = 1;
        module = helpModule_;
    }
    if (module == helpModule_)
    {
        helpModule_->setShowHidden(bHidden);
    }
    return module;
}

/********************************************************************
 * CommandLineModuleManager
 */

CommandLineModuleManager::CommandLineModuleManager(ProgramInfo *programInfo)
    : impl_(new Impl(programInfo))
{
}

CommandLineModuleManager::~CommandLineModuleManager()
{
}

void CommandLineModuleManager::setQuiet(bool bQuiet)
{
    impl_->bQuiet_ = bQuiet;
}

void CommandLineModuleManager::setSingleModule(CommandLineModuleInterface *module)
{
    impl_->singleModule_ = module;
}

void CommandLineModuleManager::addModule(CommandLineModulePointer module)
{
    impl_->addModule(move(module));
}

void CommandLineModuleManager::addModuleCMain(
        const char *name, const char *shortDescription,
        CMainFunction mainFunction)
{
    CommandLineModulePointer module(
            new CMainCommandLineModule(name, shortDescription, mainFunction));
    addModule(move(module));
}

CommandLineModuleGroup CommandLineModuleManager::addModuleGroup(
        const char *title)
{
    const char *const                 binaryName =
        impl_->programInfo_.invariantProgramName().c_str();
    CommandLineModuleGroupDataPointer group(
            new CommandLineModuleGroupData(impl_->modules_, binaryName, title));
    impl_->moduleGroups_.push_back(move(group));
    return CommandLineModuleGroup(impl_->moduleGroups_.back().get());
}

void CommandLineModuleManager::addHelpTopic(HelpTopicPointer topic)
{
    impl_->ensureHelpModuleExists();
    impl_->helpModule_->addTopic(move(topic));
}

int CommandLineModuleManager::run(int argc, char *argv[])
{
    CommandLineModuleInterface *module;
    const bool                  bMaster = (!gmx_mpi_initialized() || gmx_node_rank() == 0);
    try
    {
        module = impl_->processCommonOptions(&argc, &argv);
    }
    catch (const std::exception &)
    {
        if (bMaster && !impl_->bQuiet_)
        {
            printBinaryInformation(stderr, impl_->programInfo_,
                                   impl_->binaryInfoSettings_);
        }
        throw;
    }
    if (!bMaster)
    {
        impl_->bQuiet_ = true;
    }
    if (!impl_->bQuiet_)
    {
        FILE *out = (impl_->bStdOutInfo_ ? stdout : stderr);
        printBinaryInformation(out, impl_->programInfo_,
                               impl_->binaryInfoSettings_);
        fprintf(out, "\n");
    }
    if (module == NULL)
    {
        return 0;
    }
    int rc = module->run(argc, argv);
    if (!impl_->bQuiet_)
    {
        gmx_thanx(stderr);
    }
    return rc;
}

// static
int CommandLineModuleManager::runAsMainSingleModule(
        int argc, char *argv[], CommandLineModuleInterface *module)
{
    ProgramInfo &programInfo = gmx::init(&argc, &argv);
    try
    {
        CommandLineModuleManager manager(&programInfo);
        manager.setSingleModule(module);
        int rc = manager.run(argc, argv);
        gmx::finalize();
        return rc;
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        return processExceptionAtExit(ex);
    }
}

// static
int CommandLineModuleManager::runAsMainCMain(
        int argc, char *argv[], CMainFunction mainFunction)
{
    CMainCommandLineModule module(argv[0], NULL, mainFunction);
    return runAsMainSingleModule(argc, argv, &module);
}

/********************************************************************
 * CommandLineModuleGroupData
 */

void CommandLineModuleGroupData::addModule(const char *name,
                                           const char *description)
{
    CommandLineModuleMap::const_iterator moduleIter = allModules_.find(name);
    GMX_RELEASE_ASSERT(moduleIter != allModules_.end(),
                       "Non-existent module added to a group");
    if (description == NULL)
    {
        description = moduleIter->second->shortDescription();
        GMX_RELEASE_ASSERT(description != NULL,
                           "Module without a description added to a group");
    }
    std::string       tag(formatString("%s-%s", binaryName_, name));
    modules_.push_back(std::make_pair(tag, description));
}

/********************************************************************
 * CommandLineModuleGroup
 */

void CommandLineModuleGroup::addModule(const char *name)
{
    impl_->addModule(name, NULL);
}

void CommandLineModuleGroup::addModuleWithDescription(const char *name,
                                                      const char *description)
{
    impl_->addModule(name, description);
}

} // namespace gmx
