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
 * Implements gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlinemodulemanager.h"

#include <cstdio>

#include <string>
#include <utility>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#include "cmdlinehelpmodule.h"
#include "cmdlinemodulemanager-impl.h"

namespace gmx
{

namespace
{

//! \addtogroup module_commandline
//! \{

/********************************************************************
 * CMainCommandLineModule
 */

/*! \brief
 * Implements a CommandLineModuleInterface, given a function with C/C++ main()
 * signature.
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

        virtual void init(CommandLineModuleSettings * /*settings*/)
        {
        }
        virtual int run(int argc, char *argv[])
        {
            return mainFunction_(argc, argv);
        }
        virtual void writeHelp(const CommandLineHelpContext &context) const
        {
            writeCommandLineHelpCMain(context, name_, mainFunction_);
        }

    private:
        const char             *name_;
        const char             *shortDescription_;
        CMainFunction           mainFunction_;
};

//! \}

}   // namespace

/********************************************************************
 * CommandLineCommonOptionsHolder
 */

CommandLineCommonOptionsHolder::CommandLineCommonOptionsHolder()
    : options_(NULL, NULL), bHelp_(false), bHidden_(false),
      bQuiet_(false), bVersion_(false), bCopyright_(true),
      niceLevel_(19), bBackup_(true), bFpexcept_(false), debugLevel_(0)
{
    binaryInfoSettings_.copyright(true);
}

CommandLineCommonOptionsHolder::~CommandLineCommonOptionsHolder()
{
}

void CommandLineCommonOptionsHolder::initOptions()
{
    options_.addOption(BooleanOption("h").store(&bHelp_)
                           .description("Print help and quit"));
    options_.addOption(BooleanOption("hidden").store(&bHidden_)
                           .hidden()
                           .description("Show hidden options in help"));
    options_.addOption(BooleanOption("quiet").store(&bQuiet_)
                           .description("Do not print common startup info or quotes"));
    options_.addOption(BooleanOption("version").store(&bVersion_)
                           .description("Print extended version information and quit"));
    options_.addOption(BooleanOption("copyright").store(&bCopyright_)
                           .description("Print copyright information on startup"));
    options_.addOption(IntegerOption("nice").store(&niceLevel_)
                           .description("Set the nicelevel (default depends on command)"));
    options_.addOption(BooleanOption("backup").store(&bBackup_)
                           .description("Write backups if output files exist"));
    options_.addOption(BooleanOption("fpexcept").store(&bFpexcept_)
                           .hidden().description("Enable floating-point exceptions"));
    options_.addOption(IntegerOption("debug").store(&debugLevel_)
                           .hidden().defaultValueIfSet(1)
                           .description("Write file with debug information, "
                                        "1: short (default), 2: also x and f"));
}

bool CommandLineCommonOptionsHolder::finishOptions()
{
    options_.finish();
    binaryInfoSettings_.extendedInfo(bVersion_);
    // The latter condition suppresses the copyright with
    // -quiet -version.
    binaryInfoSettings_.copyright(bCopyright_ && !bQuiet_);
    return !bVersion_;
}

void CommandLineCommonOptionsHolder::adjustFromSettings(
        const CommandLineModuleSettings &settings)
{
    if (!options_.isSet("nice"))
    {
        niceLevel_ = settings.defaultNiceLevel();
    }
}

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
         * \param[in] binaryName     Name of the running binary
         *     (without Gromacs binary suffix or .exe on Windows).
         * \param     programContext Program information for the running binary.
         */
        Impl(const char *binaryName, CommandLineProgramContext *programContext);

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
         * Processes command-line options for the wrapper binary.
         *
         * \param[in,out] optionsHolder Common options.
         * \param[in,out] argc          On input, argc passed to run().
         *     On output, argc to be passed to the module.
         * \param[in,out] argv          On input, argv passed to run().
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
        processCommonOptions(CommandLineCommonOptionsHolder *optionsHolder,
                             int *argc, char ***argv);

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
        CommandLineProgramContext   &programContext_;
        //! Name of the binary.
        std::string                  binaryName_;
        /*! \brief
         * Module that implements help for the binary.
         *
         * The pointed module is owned by the \a modules_ container.
         */
        CommandLineHelpModule       *helpModule_;
        //! If non-NULL, run this module in single-module mode.
        CommandLineModuleInterface  *singleModule_;
        //! Stores the value set with setQuiet().
        bool                         bQuiet_;

    private:
        GMX_DISALLOW_COPY_AND_ASSIGN(Impl);
};

CommandLineModuleManager::Impl::Impl(const char                *binaryName,
                                     CommandLineProgramContext *programContext)
    : programContext_(*programContext),
      binaryName_(binaryName != NULL ? binaryName : ""),
      helpModule_(NULL), singleModule_(NULL),
      bQuiet_(false)
{
    GMX_RELEASE_ASSERT(binaryName_.find('-') == std::string::npos,
                       "Help export does not currently work with binary names with dashes");
}

void CommandLineModuleManager::Impl::addModule(CommandLineModulePointer module)
{
    GMX_ASSERT(modules_.find(module->name()) == modules_.end(),
               "Attempted to register a duplicate module name");
    ensureHelpModuleExists();
    HelpTopicPointer helpTopic(helpModule_->createModuleHelpTopic(*module));
    modules_.insert(std::make_pair(std::string(module->name()),
                                   move(module)));
    helpModule_->addTopic(move(helpTopic), false);
}

void CommandLineModuleManager::Impl::ensureHelpModuleExists()
{
    if (helpModule_ == NULL)
    {
        helpModule_ = new CommandLineHelpModule(programContext_, binaryName_,
                                                modules_, moduleGroups_);
        addModule(CommandLineModulePointer(helpModule_));
    }
}

CommandLineModuleMap::const_iterator
CommandLineModuleManager::Impl::findModuleByName(const std::string &name) const
{
    // TODO: Accept unambiguous prefixes?
    return modules_.find(name);
}

CommandLineModuleInterface *
CommandLineModuleManager::Impl::processCommonOptions(
        CommandLineCommonOptionsHolder *optionsHolder, int *argc, char ***argv)
{
    // Check if we are directly invoking a certain module.
    CommandLineModuleInterface *module = singleModule_;

    // TODO: It would be nice to propagate at least the -quiet option to
    // the modules so that they can also be quiet in response to this.

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
            CommandLineParser(optionsHolder->options())
                .parse(&argcForWrapper, *argv);
        }
        // If no action requested and there is a module specified, process it.
        if (argcForWrapper < *argc && !optionsHolder->shouldIgnoreActualModule())
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
            programContext_.setDisplayName(binaryName_ + " " + module->name());
        }
        // Recognize the common options also after the module name.
        // TODO: It could be nicer to only recognize -h/-hidden if module is not
        // null.
        CommandLineParser(optionsHolder->options())
            .skipUnknown(true).parse(argc, *argv);
    }
    if (!optionsHolder->finishOptions())
    {
        return NULL;
    }
    // If no module specified and no other action, show the help.
    // Also explicitly specifying -h for the wrapper binary goes here.
    if (module == NULL || optionsHolder->shouldShowHelp())
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
        helpModule_->setShowHidden(optionsHolder->shouldShowHidden());
    }
    return module;
}

/********************************************************************
 * CommandLineModuleManager
 */

CommandLineModuleManager::CommandLineModuleManager(
        const char *binaryName, CommandLineProgramContext *programContext)
    : impl_(new Impl(binaryName, programContext))
{
}

CommandLineModuleManager::~CommandLineModuleManager()
{
}

void CommandLineModuleManager::setQuiet(bool bQuiet)
{
    impl_->bQuiet_ = bQuiet;
}

void CommandLineModuleManager::setOutputRedirector(
        FileOutputRedirectorInterface *output)
{
    impl_->ensureHelpModuleExists();
    impl_->helpModule_->setOutputRedirector(output);
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
    const char *const                 binaryName = impl_->binaryName_.c_str();
    CommandLineModuleGroupDataPointer group(
            new CommandLineModuleGroupData(impl_->modules_, binaryName, title));
    impl_->moduleGroups_.push_back(move(group));
    return CommandLineModuleGroup(impl_->moduleGroups_.back().get());
}

void CommandLineModuleManager::addHelpTopic(HelpTopicPointer topic)
{
    impl_->ensureHelpModuleExists();
    impl_->helpModule_->addTopic(move(topic), true);
}

int CommandLineModuleManager::run(int argc, char *argv[])
{
    CommandLineModuleInterface    *module;
    const bool                     bMaster = (gmx_node_rank() == 0);
    bool                           bQuiet  = impl_->bQuiet_ || !bMaster;
    CommandLineCommonOptionsHolder optionsHolder;
    try
    {
        optionsHolder.initOptions();
        module = impl_->processCommonOptions(&optionsHolder, &argc, &argv);
    }
    catch (const std::exception &)
    {
        bQuiet |= optionsHolder.shouldBeQuiet();
        if (!bQuiet)
        {
            printBinaryInformation(stderr, impl_->programContext_,
                                   optionsHolder.binaryInfoSettings());
        }
        throw;
    }
    bQuiet |= optionsHolder.shouldBeQuiet();
    if (!bQuiet)
    {
        FILE *out = optionsHolder.startupInfoFile();
        printBinaryInformation(out, impl_->programContext_,
                               optionsHolder.binaryInfoSettings());
        fprintf(out, "\n");
    }
    if (module == NULL)
    {
        return 0;
    }

    CommandLineModuleSettings settings;
    module->init(&settings);
    optionsHolder.adjustFromSettings(settings);

    gmx_set_max_backup_count(optionsHolder.shouldBackup() ? -1 : 0);

    // Open the debug file.
    if (optionsHolder.debugLevel() > 0)
    {
        std::string filename(impl_->programContext_.programName());
        if (gmx_node_num() > 1)
        {
            filename.append(formatString("%d", gmx_node_rank()));
        }
        filename.append(".debug");

        fprintf(stderr, "Will write debug log file: %s\n", filename.c_str());
        gmx_init_debug(optionsHolder.debugLevel(), filename.c_str());
    }
    // Set the nice level unless disabled in the configuration.
    if (optionsHolder.niceLevel() != 0)
    {
        static bool bNiceSet = false; // Only set it once.
        if (!bNiceSet)
        {
            // TODO: Diagnostic if this fails and the user explicitly requested it.
            gmx_set_nice(optionsHolder.niceLevel());
            bNiceSet = true;
        }
    }
    if (optionsHolder.enableFPExceptions())
    {
        //TODO: currently it is always enabled for mdrun (verlet) and tests.
        gmx_feenableexcept();
    }

    int rc = 0;
    if (!(module == impl_->helpModule_ && !bMaster))
    {
        rc = module->run(argc, argv);
    }
    if (!bQuiet)
    {
        gmx_thanx(stderr);
    }
    return rc;
}

// static
int CommandLineModuleManager::runAsMainSingleModule(
        int argc, char *argv[], CommandLineModuleInterface *module)
{
    CommandLineProgramContext &programContext = gmx::initForCommandLine(&argc, &argv);
    try
    {
        CommandLineModuleManager manager(NULL, &programContext);
        manager.setSingleModule(module);
        int rc = manager.run(argc, argv);
        gmx::finalizeForCommandLine();
        return rc;
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        return processExceptionAtExitForCommandLine(ex);
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
