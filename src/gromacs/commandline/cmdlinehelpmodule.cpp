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
 * Implements gmx::CommandLineHelpModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gromacs/commandline/cmdlinehelpmodule.h"

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/shellcompletions.h"
#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/helpmanager.h"
#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

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

        void exportHelp(HelpExportInterface *exporter) const;

        boost::scoped_ptr<RootHelpTopic>  rootTopic_;
        const ProgramContextInterface    &programContext_;
        std::string                       binaryName_;
        const CommandLineModuleMap       &modules_;
        const CommandLineModuleGroupList &groups_;

        CommandLineHelpContext           *context_;
        const CommandLineModuleInterface *moduleOverride_;
        bool                              bHidden_;

        File                             *outputOverride_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineHelpModuleImpl);
};

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
const char *const RootHelpText::text[]  = { "" };

/*! \brief
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
         * Does not throw.
         */
        explicit RootHelpTopic(const std::string &binaryName)
            : binaryName_(binaryName)
        {
        }

        virtual void writeHelp(const HelpWriterContext &context) const;

    private:
        std::string                 binaryName_;

        GMX_DISALLOW_COPY_AND_ASSIGN(RootHelpTopic);
};

void RootHelpTopic::writeHelp(const HelpWriterContext &context) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                          "Root help is not implemented for this output format"));
    }
    {
        CommandLineCommonOptionsHolder optionsHolder;
        CommandLineHelpContext         cmdlineContext(context);
        optionsHolder.initOptions();
        cmdlineContext.setModuleDisplayName(binaryName_);
        // TODO: Add <command> [<args>] into the synopsis.
        // TODO: Propagate the -hidden option here.
        CommandLineHelpWriter(*optionsHolder.options())
            .writeHelp(cmdlineContext);
    }
    // TODO: Consider printing a list of "core" commands. Would require someone
    // to determine such a set...
    writeSubTopicList(context,
                      "Additional help is available on the following topics:");
    context.writeTextBlock(
            "To access the help, use '[PROGRAM] help <topic>'.[BR]"
            "For help on a command, use '[PROGRAM] help <command>'.");
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
};

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
            std::string targetName(
                    formatString("%s-%s", program, module->first.c_str()));
            std::string displayName(
                    formatString("[TT]%s %s[tt]", program, module->first.c_str()));
            links->addLink(linkName, targetName, displayName);
        }
    }
}

/********************************************************************
 * HelpExportMan
 */

/*! \internal \brief
 * Implements export for man pages.
 *
 * \ingroup module_commandline
 */
class HelpExportMan : public HelpExportInterface
{
    public:
        //! Initializes man page exporter.
        explicit HelpExportMan(const CommandLineHelpModuleImpl &helpModule)
            : links_(eHelpOutputFormat_Man)
        {
            initProgramLinks(&links_, helpModule);
        }

        virtual void startModuleExport() {}
        virtual void exportModuleHelp(
            const CommandLineModuleInterface &module,
            const std::string                &tag,
            const std::string                &displayName);
        virtual void finishModuleExport() {}

        virtual void startModuleGroupExport();
        virtual void exportModuleGroup(const char                *title,
                                       const ModuleGroupContents &modules);
        virtual void finishModuleGroupExport();

    private:
        HelpLinks                links_;
        boost::scoped_ptr<File>  man7File_;
        std::string              man7Footer_;
};

void HelpExportMan::exportModuleHelp(
        const CommandLineModuleInterface &module,
        const std::string                &tag,
        const std::string                &displayName)
{
    File file("man1/" + tag + ".1", "w");

    // TODO: It would be nice to remove the VERSION prefix from the version
    // string to make it shorter.
    file.writeLine(formatString(".TH %s 1 \"\" \"%s\" \"GROMACS Manual\"\n",
                                tag.c_str(), gmx_version()));
    file.writeLine(".SH NAME");
    file.writeLine(formatString("%s - %s", tag.c_str(),
                                module.shortDescription()));
    file.writeLine();

    CommandLineHelpContext context(&file, eHelpOutputFormat_Man, &links_);
    context.setModuleDisplayName(displayName);
    module.writeHelp(context);

    file.writeLine(".SH SEE ALSO");
    file.writeLine(".BR gromacs(7)");
    file.writeLine();
    file.writeLine("More information about \\fBGROMACS\\fR is available at <\\fIhttp://www.gromacs.org/\\fR>.");
}

void HelpExportMan::startModuleGroupExport()
{
    const char *const programListPlaceholder = "@PROGMANPAGES@";

    const std::string man7Template = gmx::File::readToString("man7/gromacs.7.in");
    const size_t      index        = man7Template.find(programListPlaceholder);
    GMX_RELEASE_ASSERT(index != std::string::npos,
                       "gromacs.7.in must contain a @PROGMANPAGES@ line");
    std::string header = man7Template.substr(0, index);
    man7Footer_ = man7Template.substr(index + std::strlen(programListPlaceholder));
    header      = replaceAll(header, "@VERSION@", gmx_version());
    man7File_.reset(new File("man7/gromacs.7", "w"));
    man7File_->writeLine(header);
}

void HelpExportMan::exportModuleGroup(const char                *title,
                                      const ModuleGroupContents &modules)
{
    man7File_->writeLine(formatString(".Sh \"%s\"", title));
    man7File_->writeLine(formatString(".IX Subsection \"%s\"", title));
    man7File_->writeLine(".Vb");
    man7File_->writeLine(".ta 16n");

    ModuleGroupContents::const_iterator module;
    for (module = modules.begin(); module != modules.end(); ++module)
    {
        const std::string &tag(module->first);
        man7File_->writeLine(formatString("\\&  %s\t%s",
                                          tag.c_str(), module->second));
    }

    man7File_->writeLine(".Ve");
}

void HelpExportMan::finishModuleGroupExport()
{
    man7File_->writeLine(man7Footer_);
    man7File_->close();
}

/********************************************************************
 * HelpExportHtml
 */

/*! \internal \brief
 * Implements export for HTML help.
 *
 * \ingroup module_commandline
 */
class HelpExportHtml : public HelpExportInterface
{
    public:
        //! Initializes HTML exporter.
        explicit HelpExportHtml(const CommandLineHelpModuleImpl &helpModule);

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

    private:
        void setupHeaderAndFooter();

        void writeHtmlHeader(File *file, const std::string &title) const;
        void writeHtmlFooter(File *file) const;

        HelpLinks                links_;
        boost::scoped_ptr<File>  indexFile_;
        std::string              header_;
        std::string              footer_;
};

HelpExportHtml::HelpExportHtml(const CommandLineHelpModuleImpl &helpModule)
    : links_(eHelpOutputFormat_Html)
{
    File             linksFile("links.dat", "r");
    std::string      line;
    while (linksFile.readLine(&line))
    {
        links_.addLink(line, "../online/" + line, line);
    }
    linksFile.close();
    initProgramLinks(&links_, helpModule);
    setupHeaderAndFooter();
}

void HelpExportHtml::setupHeaderAndFooter()
{
    header_ = gmx::File::readToString("header.html.in");
    header_ = replaceAll(header_, "@VERSION@", gmx_version());
    gmx::File::writeFileFromString("header.html", header_);
    header_ = replaceAll(header_, "@ROOTPATH@", "../");
    footer_ = gmx::File::readToString("footer.html");
}

void HelpExportHtml::startModuleExport()
{
    indexFile_.reset(new File("final/programs/byname.html", "w"));
    writeHtmlHeader(indexFile_.get(), "GROMACS Programs by Name");
    indexFile_->writeLine("<H3>GROMACS Programs Alphabetically</H3>");
}

void HelpExportHtml::exportModuleHelp(
        const CommandLineModuleInterface &module,
        const std::string                &tag,
        const std::string                &displayName)
{
    File file("final/programs/" + tag + ".html", "w");
    writeHtmlHeader(&file, displayName);

    CommandLineHelpContext context(&file, eHelpOutputFormat_Html, &links_);
    context.setModuleDisplayName(displayName);
    module.writeHelp(context);

    writeHtmlFooter(&file);
    file.close();

    indexFile_->writeLine(formatString("<a href=\"%s.html\">%s</a> - %s<br>",
                                       tag.c_str(), displayName.c_str(),
                                       module.shortDescription()));
}

void HelpExportHtml::finishModuleExport()
{
    writeHtmlFooter(indexFile_.get());
    indexFile_->close();
}

void HelpExportHtml::startModuleGroupExport()
{
    indexFile_.reset(new File("final/programs/bytopic.html", "w"));
    writeHtmlHeader(indexFile_.get(), "GROMACS Programs by Topic");
    indexFile_->writeLine("<H3>GROMACS Programs by Topic</H3>");
}

void HelpExportHtml::exportModuleGroup(const char                *title,
                                       const ModuleGroupContents &modules)
{
    indexFile_->writeLine(formatString("<H4>%s</H4>", title));

    ModuleGroupContents::const_iterator module;
    for (module = modules.begin(); module != modules.end(); ++module)
    {
        const std::string     &tag(module->first);
        std::string            displayName(tag);
        // TODO: This does not work if the binary name would contain a dash,
        // but that is not currently the case.
        size_t                 dashPos = displayName.find('-');
        GMX_RELEASE_ASSERT(dashPos != std::string::npos,
                           "There should always be at least one dash in the tag");
        displayName[dashPos] = ' ';
        indexFile_->writeLine(formatString("<a href=\"%s.html\">%s</a> - %s<br>",
                                           tag.c_str(), displayName.c_str(),
                                           module->second));
    }
}

void HelpExportHtml::finishModuleGroupExport()
{
    writeHtmlFooter(indexFile_.get());
    indexFile_->close();
}

void HelpExportHtml::writeHtmlHeader(File *file, const std::string &title) const
{
    file->writeLine(replaceAll(header_, "@TITLE@", title));
}

void HelpExportHtml::writeHtmlFooter(File *file) const
{
    file->writeLine(footer_);
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
    : rootTopic_(new RootHelpTopic(binaryName)), programContext_(programContext),
      binaryName_(binaryName), modules_(modules), groups_(groups),
      context_(NULL), moduleOverride_(NULL), bHidden_(false),
      outputOverride_(NULL)
{
}

void CommandLineHelpModuleImpl::exportHelp(HelpExportInterface *exporter) const
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

void CommandLineHelpModule::addTopic(HelpTopicPointer topic)
{
    impl_->rootTopic_->addSubTopic(move(topic));
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

void CommandLineHelpModule::setOutputRedirect(File *output)
{
    impl_->outputOverride_ = output;
}

int CommandLineHelpModule::run(int argc, char *argv[])
{
    // Add internal topics lazily here.
    addTopic(HelpTopicPointer(new CommandsHelpTopic(*impl_)));

    const char *const exportFormats[] = { "man", "html", "completion" };
    std::string       exportFormat;
    Options           options(NULL, NULL);
    options.addOption(StringOption("export").store(&exportFormat)
                          .enumValue(exportFormats));
    CommandLineParser(&options).parse(&argc, argv);
    if (!exportFormat.empty())
    {
        boost::scoped_ptr<HelpExportInterface> exporter;
        if (exportFormat == "man")
        {
            exporter.reset(new HelpExportMan(*impl_));
        }
        else if (exportFormat == "html")
        {
            exporter.reset(new HelpExportHtml(*impl_));
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

    File *outputFile = &File::standardOutput();
    if (impl_->outputOverride_ != NULL)
    {
        outputFile = impl_->outputOverride_;
    }
    HelpLinks                                 links(eHelpOutputFormat_Console);
    initProgramLinks(&links, *impl_);
    boost::scoped_ptr<CommandLineHelpContext> context(
            new CommandLineHelpContext(outputFile,
                                       eHelpOutputFormat_Console, &links));
    context->setShowHidden(impl_->bHidden_);
    if (impl_->moduleOverride_ != NULL)
    {
        context->setModuleDisplayName(impl_->programContext_.displayName());
        impl_->moduleOverride_->writeHelp(*context);
        return 0;
    }
    impl_->context_ = context.get();

    HelpManager helpManager(*impl_->rootTopic_, context->writerContext());
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
