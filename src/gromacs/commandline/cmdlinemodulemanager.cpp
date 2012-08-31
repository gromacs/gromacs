/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "cmdlinemodulemanager.h"

#include <cstdio>
#include <cstring>

#include <map>
#include <string>
#include <utility>

#include <boost/scoped_ptr.hpp>

#include "gromacs/legacyheaders/copyrite.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/helpmanager.h"
#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programinfo.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Container type for mapping module names to module objects.
typedef std::map<std::string, CommandLineModulePointer> CommandLineModuleMap;

namespace
{

/********************************************************************
 * RootHelpTopic
 */

struct RootHelpText
{
    static const char        name[];
    static const char        title[];
    static const char *const text[];
};

// The first two are not used.
const char        RootHelpText::name[]  = "";
const char        RootHelpText::title[] = "";
const char *const RootHelpText::text[]  = {
    "Usage: [PROGRAM] <command> [<args>]",
};

/*! \internal \brief
 * Help topic that forms the root of the help tree for the help subcommand.
 *
 * \ingroup module_commandline
 */
class RootHelpTopic : public CompositeHelpTopic<RootHelpText>
{
    public:
        /*! \brief
         * Creates a root help topic.
         *
         * \param[in] modules  List of modules for to use for module listings.
         *
         * Does not throw.
         */
        explicit RootHelpTopic(const CommandLineModuleMap &modules)
            : modules_(modules)
        {
        }

        virtual void writeHelp(const HelpWriterContext &context) const;

    private:
        void printModuleList(const HelpWriterContext &context) const;

        const CommandLineModuleMap &modules_;

        GMX_DISALLOW_COPY_AND_ASSIGN(RootHelpTopic);
};

void RootHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    if (context.outputFormat() == eHelpOutputFormat_Console)
    {
        context.writeTextBlock(helpText());
        // TODO: If/when this list becomes long, it may be better to only print
        // "common" commands here, and have a separate topic (e.g.,
        // "help commands") that prints the full list.
        printModuleList(context);
        context.writeTextBlock(
                "For additional help on a command, use '[PROGRAM] help <command>'");
    }
    writeSubTopicList(context,
                      "\nAdditional help is available on the following topics:");
    if (context.outputFormat() == eHelpOutputFormat_Console)
    {
        context.writeTextBlock(
                "To access the help, use '[PROGRAM] help <topic>'.");
    }
}

void RootHelpTopic::printModuleList(const HelpWriterContext &context) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        GMX_THROW(NotImplementedError(
                          "Module list is not implemented for this output format"));
    }
    int maxNameLength = 0;
    CommandLineModuleMap::const_iterator module;
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        int nameLength = static_cast<int>(module->first.length());
        if (module->second->shortDescription() != NULL
            && nameLength > maxNameLength)
        {
            maxNameLength = nameLength;
        }
    }
    File              &file = context.outputFile();
    TextTableFormatter formatter;
    formatter.addColumn(NULL, maxNameLength + 1, false);
    formatter.addColumn(NULL, 72 - maxNameLength, true);
    formatter.setFirstColumnIndent(4);
    file.writeLine();
    file.writeLine("Available commands:");
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        const char *name        = module->first.c_str();
        const char *description = module->second->shortDescription();
        if (description != NULL)
        {
            formatter.clear();
            formatter.addColumnLine(0, name);
            formatter.addColumnLine(1, description);
            file.writeString(formatter.formatRow());
        }
    }
}

/********************************************************************
 * ModuleHelpTopic
 */

/*! \internal \brief
 * Help topic wrapper for a command-line module.
 *
 * This class implements HelpTopicInterface such that it wraps a
 * CommandLineModuleInterface, allowing subcommand "help <command>"
 * to produce the help for "<command>".
 *
 * \ingroup module_commandline
 */
class ModuleHelpTopic : public HelpTopicInterface
{
    public:
        //! Constructs a help topic for a specific module.
        explicit ModuleHelpTopic(const CommandLineModuleInterface &module)
            : module_(module)
        {
        }

        virtual const char *name() const { return module_.name(); }
        virtual const char *title() const { return NULL; }
        virtual bool hasSubTopics() const { return false; }
        virtual const HelpTopicInterface *findSubTopic(const char * /*name*/) const
        {
            return NULL;
        }
        virtual void writeHelp(const HelpWriterContext &context) const;

    private:
        const CommandLineModuleInterface &module_;

        GMX_DISALLOW_COPY_AND_ASSIGN(ModuleHelpTopic);
};

void ModuleHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    module_.writeHelp(context);
}

/********************************************************************
 * HelpExportInterface
 */

class HelpExportInterface
{
    public:
        virtual ~HelpExportInterface() {};

        virtual void startModuleExport() = 0;
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module) = 0;
        virtual void finishModuleExport() = 0;
        virtual void exportHelpTopics(const RootHelpTopic &root) = 0;
};

/********************************************************************
 * HelpExportReStructuredText
 */

class HelpExportReStructuredText : public HelpExportInterface
{
    public:
        virtual void startModuleExport() {}
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module);
        virtual void finishModuleExport() {}
        virtual void exportHelpTopics(const RootHelpTopic &root);
};

void HelpExportReStructuredText::exportModuleHelp(
        const std::string &tag, const CommandLineModuleInterface &module)
{
    File              moduleHelpFile(tag + ".txt", "w");
    HelpWriterContext context(&moduleHelpFile, eHelpOutputFormat_Export);
    HelpWriterContext titleSec(context.createSubsection(tag));
    std::string       desc = module.shortDescription();
    HelpWriterContext subtitleSec(titleSec.createSubsection(desc));
    //moduleHelpFile.writeLine(formatString(":Author:  Gromacs team"));
    // TODO: Implement date generation
    moduleHelpFile.writeLine(formatString(":Date:    2012-08-24"));
    // TODO: The current version string is too long
    moduleHelpFile.writeLine(formatString(":Version: %s", GromacsVersion()));
    moduleHelpFile.writeLine(formatString(":Manual section: 1"));
    moduleHelpFile.writeLine(formatString(":Manual group:   Gromacs Manual"));
    moduleHelpFile.writeLine();
    module.writeHelp(context);
}

void HelpExportReStructuredText::exportHelpTopics(const RootHelpTopic &root)
{
    File              topicsHelpFile("help-topics.txt", "w");
    HelpWriterContext context(&topicsHelpFile, eHelpOutputFormat_Export);
    root.writeHelp(context);
}

}   // namespace

/********************************************************************
 * CommandLineHelpModule
 */

/*! \internal \brief
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
         * \param[in] modules  List of modules for to use for module listings.
         * \throws    std::bad_alloc if out of memory.
         */
        explicit CommandLineHelpModule(const CommandLineModuleMap &modules);

        /*! \brief
         * Adds a top-level help topic.
         *
         * \param[in] topic  Help topic to add.
         * \throws    std::bad_alloc if out of memory.
         */
        void addTopic(HelpTopicPointer topic);

        virtual const char *name() const { return "help"; }
        virtual const char *shortDescription() const
        {
            return "Print help information";
        }

        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const HelpWriterContext &context) const;

        //! Prints usage message to stderr.
        void printUsage() const;

    private:
        void exportHelp(HelpExportInterface *exporter) const;

        boost::scoped_ptr<RootHelpTopic>  rootTopic_;
        const CommandLineModuleMap       &modules_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineHelpModule);
};

CommandLineHelpModule::CommandLineHelpModule(const CommandLineModuleMap &modules)
    : rootTopic_(new RootHelpTopic(modules)), modules_(modules)
{
}

void CommandLineHelpModule::addTopic(HelpTopicPointer topic)
{
    rootTopic_->addSubTopic(move(topic));
}

int CommandLineHelpModule::run(int argc, char *argv[])
{
    // TODO: It would be nicer to use a CommandLineParser here.
    if (argc == 3 && std::strcmp(argv[1], "-export") == 0)
    {
        boost::scoped_ptr<HelpExportInterface> exporter;
        if (std::strcmp(argv[2], "rst") == 0)
        {
            exporter.reset(new HelpExportReStructuredText);
        }
        else
        {
            GMX_THROW(InvalidInputError(
                              formatString("Unknown help export format '%s'",
                                           argv[2])));
        }
        exportHelp(exporter.get());
        return 0;
    }
    HelpWriterContext context(&File::standardOutput(),
                              eHelpOutputFormat_Console);
    HelpManager       helpManager(*rootTopic_, context);
    try
    {
        for (int i = 1; i < argc; ++i)
        {
            helpManager.enterTopic(argv[i]);
        }
    }
    catch (const InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 2;
    }
    helpManager.writeCurrentTopic();
    fprintf(stderr, "\n");
    return 0;
}

void CommandLineHelpModule::writeHelp(const HelpWriterContext &context) const
{
    context.writeTextBlock(
            "Usage: [PROGRAM] help [<command>|<topic> [<subtopic> [...]]]");
    // TODO: More information.
}

void CommandLineHelpModule::printUsage() const
{
    HelpWriterContext context(&File::standardError(),
                              eHelpOutputFormat_Console);
    rootTopic_->writeHelp(context);
}

void CommandLineHelpModule::exportHelp(HelpExportInterface *exporter) const
{
    // TODO: Would be nicer to have the file names supplied by the build system
    // and/or export a list of files from here.
    const char *program = ProgramInfo::getInstance().programName().c_str();
    CommandLineModuleMap::const_iterator module;
    exporter->startModuleExport();
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        const char *moduleName = module->first.c_str();
        // For testing, only export the select module.
        if (std::strcmp(moduleName, "select") != 0)
        {
            continue;
        }
        std::string tag(formatString("%s-%s", program, moduleName));
        exporter->exportModuleHelp(tag, *module->second);
    }
    exporter->finishModuleExport();
    exporter->exportHelpTopics(*rootTopic_);
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
         * \param[in] programInfo  Program information for the running binary.
         */
        explicit Impl(const ProgramInfo &programInfo);

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
         * Maps module names to module objects.
         *
         * Owns the contained modules.
         */
        CommandLineModuleMap    modules_;
        //! Information about the currently running program.
        const ProgramInfo      &programInfo_;
        /*! \brief
         * Module that implements help for the binary.
         *
         * The pointed module is owned by the \a modules_ container.
         */
        CommandLineHelpModule  *helpModule_;
};

CommandLineModuleManager::Impl::Impl(const ProgramInfo &programInfo)
    : programInfo_(programInfo)
{
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
    return findModuleByName(binaryName);
}

/********************************************************************
 * CommandLineModuleManager
 */

CommandLineModuleManager::CommandLineModuleManager(const ProgramInfo &programInfo)
    : impl_(new Impl(programInfo))
{
    impl_->helpModule_ = new CommandLineHelpModule(impl_->modules_);
    addModule(CommandLineModulePointer(impl_->helpModule_));
}

CommandLineModuleManager::~CommandLineModuleManager()
{
}

void CommandLineModuleManager::addModule(CommandLineModulePointer module)
{
    GMX_ASSERT(impl_->modules_.find(module->name()) == impl_->modules_.end(),
               "Attempted to register a duplicate module name");
    HelpTopicPointer helpTopic(new ModuleHelpTopic(*module));
    impl_->modules_.insert(std::make_pair(std::string(module->name()),
                                          move(module)));
    addHelpTopic(move(helpTopic));
}

void CommandLineModuleManager::addHelpTopic(HelpTopicPointer topic)
{
    impl_->helpModule_->addTopic(move(topic));
}

int CommandLineModuleManager::run(int argc, char *argv[])
{
    int argOffset = 0;
    CommandLineModuleMap::const_iterator module
        = impl_->findModuleFromBinaryName(impl_->programInfo_);
    if (module == impl_->modules_.end())
    {
        if (argc < 2)
        {
            impl_->helpModule_->printUsage();
            return 2;
        }
        module    = impl_->findModuleByName(argv[1]);
        argOffset = 1;
    }
    if (module == impl_->modules_.end())
    {
        fprintf(stderr, "Unknown command: '%s'\n\n", argv[1]);
        impl_->helpModule_->printUsage();
        return 2;
    }
    return module->second->run(argc - argOffset, argv + argOffset);
}

} // namespace gmx
