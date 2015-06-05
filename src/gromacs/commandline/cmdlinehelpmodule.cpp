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
 * Implements gmx::CommandLineHelpModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlinehelpmodule.h"

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/helpmanager.h"
#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/fileredirector.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

#include "shellcompletions.h"

namespace gmx
{

namespace
{
class HelpExportInterface;
class RootHelpTopic;
}   // namespace

/********************************************************************
 * CommandLineHelpModuleImpl declaration
 */

class CommandLineHelpModuleImpl
{
    public:
        CommandLineHelpModuleImpl(const ProgramContextInterface    &programContext,
                                  const std::string                &binaryName,
                                  const CommandLineModuleMap       &modules,
                                  const CommandLineModuleGroupList &groups);

        void exportHelp(HelpExportInterface *exporter);

        boost::scoped_ptr<RootHelpTopic>  rootTopic_;
        const ProgramContextInterface    &programContext_;
        std::string                       binaryName_;
        const CommandLineModuleMap       &modules_;
        const CommandLineModuleGroupList &groups_;

        CommandLineHelpContext           *context_;
        const CommandLineModuleInterface *moduleOverride_;
        bool                              bHidden_;

        FileOutputRedirectorInterface    *outputRedirector_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineHelpModuleImpl);
};

namespace
{

/********************************************************************
 * HelpExportInterface
 */

/*! \brief
 * Callbacks for exporting help information for command-line modules.
 *
 * \ingroup module_commandline
 */
class HelpExportInterface
{
    public:
        //! Shorthand for a list of modules contained in a group.
        typedef CommandLineModuleGroupData::ModuleList ModuleGroupContents;

        virtual ~HelpExportInterface() {};

        /*! \brief
         * Called once before exporting individual modules.
         *
         * Can, e.g., open shared output files (e.g., if the output is written
         * into a single file, or if a separate index is required) and write
         * headers into them.
         */
        virtual void startModuleExport() = 0;
        /*! \brief
         * Called to export the help for each module.
         *
         * \param[in] module      Module for which the help should be exported.
         * \param[in] tag         Unique tag for the module (gmx-something).
         * \param[in] displayName Display name for the module (gmx something).
         */
        virtual void exportModuleHelp(
            const CommandLineModuleInterface &module,
            const std::string                &tag,
            const std::string                &displayName) = 0;
        /*! \brief
         * Called after all modules have been exported.
         *
         * Can close files opened in startModuleExport(), write footers to them
         * etc.
         */
        virtual void finishModuleExport() = 0;

        /*! \brief
         * Called once before exporting module groups.
         *
         * Can, e.g., open a single output file for listing all the groups.
         */
        virtual void startModuleGroupExport() = 0;
        /*! \brief
         * Called to export the help for each module group.
         *
         * \param[in] title    Title for the group.
         * \param[in] modules  List of modules in the group.
         */
        virtual void exportModuleGroup(const char                *title,
                                       const ModuleGroupContents &modules) = 0;
        /*! \brief
         * Called after all module groups have been exported.
         *
         * Can close files opened in startModuleGroupExport(), write footers to them
         * etc.
         */
        virtual void finishModuleGroupExport() = 0;

        /*! \brief
         * Called to export the help for a top-level topic.
         *
         * \param[in] topic   Topic to export.
         */
        virtual void exportTopic(const HelpTopicInterface &topic) = 0;
};

/********************************************************************
 * RootHelpTopic
 */

struct RootHelpText
{
    static const char        title[];
    static const char *const text[];
};

// These are used for the gmx.1 man page.
// TODO: Do not hardcode them here, but pass them from the outside to make this
// code more generic.
const char        RootHelpText::title[] = "molecular dynamics simulation suite";
const char *const RootHelpText::text[]  = {
    "|Gromacs| is a full-featured suite of programs to perform molecular",
    "dynamics simulations, i.e., to simulate the behavior of systems with",
    "hundreds to millions of particles using Newtonian equations of motion.",
    "It is primarily used for research on proteins, lipids, and polymers, but",
    "can be applied to a wide variety of chemical and biological research",
    "questions.",
};

/*! \brief
 * Help topic that forms the root of the help tree for the help subcommand.
 *
 * \ingroup module_commandline
 */
class RootHelpTopic : public AbstractCompositeHelpTopic
{
    public:
        /*! \brief
         * Creates a root help topic.
         *
         * Does not throw.
         */
        explicit RootHelpTopic(const CommandLineHelpModuleImpl &helpModule)
            : helpModule_(helpModule)
        {
        }

        virtual const char *name() const { return helpModule_.binaryName_.c_str(); }
        virtual const char *title() const { return title_.c_str(); }

        //! Adds a top-level topic and optionally marks it as exported.
        void addTopic(HelpTopicPointer topic, bool bExported)
        {
            if (bExported)
            {
                exportedTopics_.push_back(topic->name());
            }
            addSubTopic(move(topic));
        }
        //! Exports all the top-level topics with the given exporter.
        void exportHelp(HelpExportInterface *exporter);

        virtual void writeHelp(const HelpWriterContext &context) const;

    private:
        // unused because of the writeHelp() override
        virtual std::string helpText() const { return ""; }

        const CommandLineHelpModuleImpl  &helpModule_;
        std::string                       title_;
        std::vector<std::string>          exportedTopics_;

        GMX_DISALLOW_COPY_AND_ASSIGN(RootHelpTopic);
};

void RootHelpTopic::exportHelp(HelpExportInterface *exporter)
{
    std::vector<std::string>::const_iterator topicName;
    for (topicName = exportedTopics_.begin();
         topicName != exportedTopics_.end();
         ++topicName)
    {
        const HelpTopicInterface *topic = findSubTopic(topicName->c_str());
        GMX_RELEASE_ASSERT(topic != NULL, "Exported help topic no longer found");
        exporter->exportTopic(*topic);
    }
    // For now, the title is only set for the export to make it not appear in
    // console output, which makes things consistent for 'gmx help' and
    // 'gmx help <command>'.
    title_ = RootHelpText::title;
    exporter->exportTopic(*this);
}

void RootHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    {
        CommandLineCommonOptionsHolder            optionsHolder;
        boost::scoped_ptr<CommandLineHelpContext> cmdlineContext;
        if (helpModule_.context_ != NULL)
        {
            cmdlineContext.reset(new CommandLineHelpContext(*helpModule_.context_));
        }
        else
        {
            cmdlineContext.reset(new CommandLineHelpContext(context));
        }
        cmdlineContext->setModuleDisplayName(helpModule_.binaryName_);
        optionsHolder.initOptions();
        Options &options = *optionsHolder.options();
        options.setDescription(RootHelpText::text);
        // TODO: Add <command> [<args>] into the synopsis.
        CommandLineHelpWriter(options)
            .setShowDescriptions(context.outputFormat() != eHelpOutputFormat_Console)
            .writeHelp(*cmdlineContext);
    }
    if (context.outputFormat() == eHelpOutputFormat_Console)
    {
        // TODO: Consider printing a list of "core" commands. Would require someone
        // to determine such a set...
        writeSubTopicList(context,
                          "Additional help is available on the following topics:");
        context.writeTextBlock("To access the help, use '[PROGRAM] help <topic>'.");
        context.writeTextBlock("For help on a command, use '[PROGRAM] help <command>'.");
    }
    else
    {
        // TODO: This should not really end up on the HTML page.
        context.writeTitle(formatString("%s commands", helpModule_.binaryName_.c_str()));
        context.writeTextBlock(
                "The following commands are available. Please refer to their "
                "individual man pages or [TT][PROGRAM] help <command>[tt] "
                "for further details.");
        context.writeTextBlock("");
        context.writeTextBlock(".. include:: /fragments/bytopic-man.rst");
    }
}

/********************************************************************
 * CommandsHelpTopic
 */

/*! \brief
 * Help topic for listing the commands.
 *
 * \ingroup module_commandline
 */
class CommandsHelpTopic : public HelpTopicInterface
{
    public:
        /*! \brief
         * Creates a command list help topic.
         *
         * \param[in]     helpModule Help module to get module information from.
         *
         * Does not throw.
         */
        explicit CommandsHelpTopic(const CommandLineHelpModuleImpl &helpModule)
            : helpModule_(helpModule)
        {
        }

        virtual const char *name() const { return "commands"; }
        virtual const char *title() const { return "List of available commands"; }
        virtual bool hasSubTopics() const { return false; }
        virtual const HelpTopicInterface *findSubTopic(const char * /*name*/) const
        {
            return NULL;
        }

        virtual void writeHelp(const HelpWriterContext &context) const;

    private:
        const CommandLineHelpModuleImpl &helpModule_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandsHelpTopic);
};

void CommandsHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        GMX_THROW(NotImplementedError(
                          "Module list is not implemented for this output format"));
    }
    int maxNameLength = 0;
    const CommandLineModuleMap           &modules = helpModule_.modules_;
    CommandLineModuleMap::const_iterator  module;
    for (module = modules.begin(); module != modules.end(); ++module)
    {
        int nameLength = static_cast<int>(module->first.length());
        if (module->second->shortDescription() != NULL
            && nameLength > maxNameLength)
        {
            maxNameLength = nameLength;
        }
    }
    context.writeTextBlock(
            "Usage: [PROGRAM] [<options>] <command> [<args>][PAR]"
            "Available commands:");
    File              &file = context.outputFile();
    TextTableFormatter formatter;
    formatter.addColumn(NULL, maxNameLength + 1, false);
    formatter.addColumn(NULL, 72 - maxNameLength, true);
    formatter.setFirstColumnIndent(4);
    for (module = modules.begin(); module != modules.end(); ++module)
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
    context.writeTextBlock(
            "For help on a command, use '[PROGRAM] help <command>'.");
}

/********************************************************************
 * ModuleHelpTopic
 */

/*! \brief
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
        ModuleHelpTopic(const CommandLineModuleInterface &module,
                        const CommandLineHelpModuleImpl  &helpModule)
            : module_(module), helpModule_(helpModule)
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
        const CommandLineHelpModuleImpl  &helpModule_;

        GMX_DISALLOW_COPY_AND_ASSIGN(ModuleHelpTopic);
};

void ModuleHelpTopic::writeHelp(const HelpWriterContext & /*context*/) const
{
    CommandLineHelpContext context(*helpModule_.context_);
    const char *const      program = helpModule_.binaryName_.c_str();
    context.setModuleDisplayName(formatString("%s %s", program, module_.name()));
    module_.writeHelp(context);
}

/********************************************************************
 * HelpExportReStructuredText
 */

/*! \internal \brief
 * Adds hyperlinks to modules within this binary.
 *
 * \param[in,out] links      Links are added here.
 * \param[in]     helpModule Help module to get module information from.
 * \throws        std::bad_alloc if out of memory.
 *
 * Initializes a HelpLinks object with links to modules defined in
 * \p helpModule.
 *
 * \ingroup module_commandline
 */
void initProgramLinks(HelpLinks *links, const CommandLineHelpModuleImpl &helpModule)
{
    const char *const                    program = helpModule.binaryName_.c_str();
    CommandLineModuleMap::const_iterator module;
    for (module = helpModule.modules_.begin();
         module != helpModule.modules_.end();
         ++module)
    {
        if (module->second->shortDescription() != NULL)
        {
            std::string linkName("[gmx-" + module->first + "]");
            const char *name = module->first.c_str();
            std::string reference(
                    formatString(":doc:`%s %s <%s-%s>`", program, name, program, name));
            std::string displayName(
                    formatString("[TT]%s %s[tt]", program, name));
            links->addLink(linkName, reference, displayName);
        }
    }
}

/*! \internal \brief
 * Implements export for web pages as reStructuredText.
 *
 * \ingroup module_commandline
 */
class HelpExportReStructuredText : public HelpExportInterface
{
    public:
        //! Initializes reST exporter.
        explicit HelpExportReStructuredText(
            const CommandLineHelpModuleImpl &helpModule);

        virtual void startModuleExport();
        virtual void exportModuleHelp(
            const CommandLineModuleInterface &module,
            const std::string                &tag,
            const std::string                &displayName);
        virtual void finishModuleExport();

        virtual void startModuleGroupExport();
        virtual void exportModuleGroup(const char                *title,
                                       const ModuleGroupContents &modules);
        virtual void finishModuleGroupExport();

        virtual void exportTopic(const HelpTopicInterface &topic);

    private:
        FileOutputRedirectorInterface  *outputRedirector_;
        const std::string              &binaryName_;
        HelpLinks                       links_;
        boost::scoped_ptr<File>         indexFile_;
        boost::scoped_ptr<File>         manPagesFile_;
};

HelpExportReStructuredText::HelpExportReStructuredText(
        const CommandLineHelpModuleImpl &helpModule)
    : outputRedirector_(helpModule.outputRedirector_),
      binaryName_(helpModule.binaryName_),
      links_(eHelpOutputFormat_Rst)
{
    File             linksFile("links.dat", "r");
    std::string      line;
    while (linksFile.readLine(&line))
    {
        links_.addLink("[REF]." + line + "[ref]",
                       formatString(":ref:`.%s <%s>`", line.c_str(), line.c_str()),
                       line);
        links_.addLink("[REF]" + line + "[ref]", formatString(":ref:`%s`", line.c_str()), line);
    }
    linksFile.close();
    initProgramLinks(&links_, helpModule);
}

void HelpExportReStructuredText::startModuleExport()
{
    indexFile_.reset(
            new File(outputRedirector_->openFileForWriting("fragments/byname.rst")));
    indexFile_->writeLine(formatString("* :doc:`%s </onlinehelp/%s>` - %s",
                                       binaryName_.c_str(), binaryName_.c_str(),
                                       RootHelpText::title));
    manPagesFile_.reset(
            new File(outputRedirector_->openFileForWriting("conf-man.py")));
    manPagesFile_->writeLine("man_pages = [");
}

void HelpExportReStructuredText::exportModuleHelp(
        const CommandLineModuleInterface &module,
        const std::string                &tag,
        const std::string                &displayName)
{
    // TODO: Ideally, the file would only be touched if it really changes.
    // This would make Sphinx reruns much faster.
    File file(outputRedirector_->openFileForWriting("onlinehelp/" + tag + ".rst"));
    file.writeLine(formatString(".. _%s:", displayName.c_str()));
    if (0 == displayName.compare(binaryName_ + " mdrun"))
    {
        // Make an extra link target for the convenience of
        // MPI-specific documentation
        file.writeLine(".. _mdrun_mpi:");
    }
    file.writeLine();

    CommandLineHelpContext context(&file, eHelpOutputFormat_Rst, &links_, binaryName_);
    context.enterSubSection(displayName);
    context.setModuleDisplayName(displayName);
    module.writeHelp(context);

    file.writeLine();
    file.writeLine(".. only:: man");
    file.writeLine();
    file.writeLine("   See also");
    file.writeLine("   --------");
    file.writeLine();
    file.writeLine(formatString("   :manpage:`%s(1)`", binaryName_.c_str()));
    file.writeLine();
    file.writeLine("   More information about |Gromacs| is available at <http://www.gromacs.org/>.");
    file.close();

    indexFile_->writeLine(formatString("* :doc:`%s </onlinehelp/%s>` - %s",
                                       displayName.c_str(), tag.c_str(),
                                       module.shortDescription()));
    manPagesFile_->writeLine(
            formatString("    ('onlinehelp/%s', '%s', \"%s\", '', 1),",
                         tag.c_str(), tag.c_str(), module.shortDescription()));
}

void HelpExportReStructuredText::finishModuleExport()
{
    indexFile_->close();
    indexFile_.reset();
    // TODO: Generalize.
    manPagesFile_->writeLine(
            formatString("    ('onlinehelp/%s', '%s', '%s', '', 1)",
                         binaryName_.c_str(), binaryName_.c_str(),
                         RootHelpText::title));
    manPagesFile_->writeLine("]");
    manPagesFile_->close();
    manPagesFile_.reset();
}

void HelpExportReStructuredText::startModuleGroupExport()
{
    indexFile_.reset(
            new File(outputRedirector_->openFileForWriting("fragments/bytopic.rst")));
    manPagesFile_.reset(
            new File(outputRedirector_->openFileForWriting("fragments/bytopic-man.rst")));
}

void HelpExportReStructuredText::exportModuleGroup(
        const char                *title,
        const ModuleGroupContents &modules)
{
    indexFile_->writeLine(title);
    indexFile_->writeLine(std::string(std::strlen(title), '^'));
    manPagesFile_->writeLine(title);
    manPagesFile_->writeLine(std::string(std::strlen(title), '^'));

    ModuleGroupContents::const_iterator module;
    for (module = modules.begin(); module != modules.end(); ++module)
    {
        const std::string     &tag(module->first);
        std::string            displayName(tag);
        // TODO: This does not work if the binary name would contain a dash,
        // but that is not currently the case.
        const size_t           dashPos = displayName.find('-');
        GMX_RELEASE_ASSERT(dashPos != std::string::npos,
                           "There should always be at least one dash in the tag");
        displayName[dashPos] = ' ';
        indexFile_->writeLine(formatString(":doc:`%s </onlinehelp/%s>`\n  %s",
                                           displayName.c_str(), tag.c_str(),
                                           module->second));
        manPagesFile_->writeLine(formatString(":manpage:`%s(1)`\n  %s",
                                              tag.c_str(),
                                              module->second));
    }
    indexFile_->writeLine();
    manPagesFile_->writeLine();
}

void HelpExportReStructuredText::finishModuleGroupExport()
{
    indexFile_->close();
    indexFile_.reset();
    manPagesFile_->close();
    manPagesFile_.reset();
}

void HelpExportReStructuredText::exportTopic(const HelpTopicInterface &topic)
{
    const std::string      path("onlinehelp/" + std::string(topic.name()) + ".rst");
    File                   file(outputRedirector_->openFileForWriting(path));
    CommandLineHelpContext context(&file, eHelpOutputFormat_Rst, &links_,
                                   binaryName_);
    HelpManager            manager(topic, context.writerContext());
    manager.writeCurrentTopic();
}

/********************************************************************
 * HelpExportCompletion
 */

/*! \internal \brief
 * Implements export for command-line completion.
 *
 * \ingroup module_commandline
 */
class HelpExportCompletion : public HelpExportInterface
{
    public:
        //! Initializes completion exporter.
        explicit HelpExportCompletion(const CommandLineHelpModuleImpl &helpModule);

        virtual void startModuleExport();
        virtual void exportModuleHelp(
            const CommandLineModuleInterface &module,
            const std::string                &tag,
            const std::string                &displayName);
        virtual void finishModuleExport();

        virtual void startModuleGroupExport() {}
        virtual void exportModuleGroup(const char                * /*title*/,
                                       const ModuleGroupContents & /*modules*/) {}
        virtual void finishModuleGroupExport() {}

        virtual void exportTopic(const HelpTopicInterface & /*topic*/) {}

    private:
        ShellCompletionWriter    bashWriter_;
        std::vector<std::string> modules_;
};

HelpExportCompletion::HelpExportCompletion(
        const CommandLineHelpModuleImpl &helpModule)
    : bashWriter_(helpModule.binaryName_, eShellCompletionFormat_Bash)
{
}

void HelpExportCompletion::startModuleExport()
{
    bashWriter_.startCompletions();
}

void HelpExportCompletion::exportModuleHelp(
        const CommandLineModuleInterface &module,
        const std::string                 & /*tag*/,
        const std::string                 & /*displayName*/)
{
    modules_.push_back(module.name());
    {
        CommandLineHelpContext context(&bashWriter_);
        // We use the display name to pass the name of the module to the
        // completion writer.
        context.setModuleDisplayName(module.name());
        module.writeHelp(context);
    }
}

void HelpExportCompletion::finishModuleExport()
{
    CommandLineCommonOptionsHolder optionsHolder;
    optionsHolder.initOptions();
    bashWriter_.writeWrapperCompletions(modules_, *optionsHolder.options());
    bashWriter_.finishCompletions();
}

}   // namespace

/********************************************************************
 * CommandLineHelpModuleImpl implementation
 */

CommandLineHelpModuleImpl::CommandLineHelpModuleImpl(
        const ProgramContextInterface    &programContext,
        const std::string                &binaryName,
        const CommandLineModuleMap       &modules,
        const CommandLineModuleGroupList &groups)
    : rootTopic_(new RootHelpTopic(*this)), programContext_(programContext),
      binaryName_(binaryName), modules_(modules), groups_(groups),
      context_(NULL), moduleOverride_(NULL), bHidden_(false),
      outputRedirector_(&defaultFileOutputRedirector())
{
}

void CommandLineHelpModuleImpl::exportHelp(HelpExportInterface *exporter)
{
    // TODO: Would be nicer to have the file names supplied by the build system
    // and/or export a list of files from here.
    const char *const program = binaryName_.c_str();

    exporter->startModuleExport();
    CommandLineModuleMap::const_iterator module;
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        if (module->second->shortDescription() != NULL)
        {
            const char *const moduleName = module->first.c_str();
            std::string       tag(formatString("%s-%s", program, moduleName));
            std::string       displayName(formatString("%s %s", program, moduleName));
            exporter->exportModuleHelp(*module->second, tag, displayName);
        }
    }
    exporter->finishModuleExport();

    exporter->startModuleGroupExport();
    CommandLineModuleGroupList::const_iterator group;
    for (group = groups_.begin(); group != groups_.end(); ++group)
    {
        exporter->exportModuleGroup((*group)->title(), (*group)->modules());
    }
    exporter->finishModuleGroupExport();

    rootTopic_->exportHelp(exporter);
}

/********************************************************************
 * CommandLineHelpModule
 */

CommandLineHelpModule::CommandLineHelpModule(
        const ProgramContextInterface    &programContext,
        const std::string                &binaryName,
        const CommandLineModuleMap       &modules,
        const CommandLineModuleGroupList &groups)
    : impl_(new Impl(programContext, binaryName, modules, groups))
{
}

CommandLineHelpModule::~CommandLineHelpModule()
{
}

HelpTopicPointer CommandLineHelpModule::createModuleHelpTopic(
        const CommandLineModuleInterface &module) const
{
    return HelpTopicPointer(new ModuleHelpTopic(module, *impl_));
}

void CommandLineHelpModule::addTopic(HelpTopicPointer topic, bool bExported)
{
    impl_->rootTopic_->addTopic(move(topic), bExported);
}

void CommandLineHelpModule::setShowHidden(bool bHidden)
{
    impl_->bHidden_ = bHidden;
}

void CommandLineHelpModule::setModuleOverride(
        const CommandLineModuleInterface &module)
{
    impl_->moduleOverride_ = &module;
}

void CommandLineHelpModule::setOutputRedirector(
        FileOutputRedirectorInterface *output)
{
    impl_->outputRedirector_ = output;
}

int CommandLineHelpModule::run(int argc, char *argv[])
{
    // Add internal topics lazily here.
    addTopic(HelpTopicPointer(new CommandsHelpTopic(*impl_)), false);

    const char *const exportFormats[] = { "rst", "completion" };
    std::string       exportFormat;
    Options           options(NULL, NULL);
    options.addOption(StringOption("export").store(&exportFormat)
                          .enumValue(exportFormats));
    CommandLineParser(&options).parse(&argc, argv);
    if (!exportFormat.empty())
    {
        boost::scoped_ptr<HelpExportInterface> exporter;
        if (exportFormat == "rst")
        {
            exporter.reset(new HelpExportReStructuredText(*impl_));
        }
        else if (exportFormat == "completion")
        {
            exporter.reset(new HelpExportCompletion(*impl_));
        }
        else
        {
            GMX_THROW(NotImplementedError("This help format is not implemented"));
        }
        impl_->exportHelp(exporter.get());
        return 0;
    }

    File                  &outputFile = impl_->outputRedirector_->standardOutput();
    HelpLinks              links(eHelpOutputFormat_Console);
    initProgramLinks(&links, *impl_);
    CommandLineHelpContext context(&outputFile, eHelpOutputFormat_Console, &links,
                                   impl_->binaryName_);
    context.setShowHidden(impl_->bHidden_);
    if (impl_->moduleOverride_ != NULL)
    {
        context.setModuleDisplayName(impl_->programContext_.displayName());
        impl_->moduleOverride_->writeHelp(context);
        return 0;
    }
    impl_->context_ = &context;

    HelpManager helpManager(*impl_->rootTopic_, context.writerContext());
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
    return 0;
}

void CommandLineHelpModule::writeHelp(const CommandLineHelpContext &context) const
{
    const HelpWriterContext &writerContext = context.writerContext();
    // TODO: Implement.
    if (writerContext.outputFormat() != eHelpOutputFormat_Console)
    {
        return;
    }
    writerContext.writeTextBlock(
            "Usage: [PROGRAM] help [<command>|<topic> [<subtopic> [...]]]");
    // TODO: More information.
}

} // namespace gmx
