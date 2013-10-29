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

#include <algorithm>
#include <map>
#include <string>
#include <utility>

#include <boost/scoped_ptr.hpp>

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/smalloc.h"

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/helpmanager.h"
#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/programinfo.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Container type for mapping module names to module objects.
typedef std::map<std::string, CommandLineModulePointer> CommandLineModuleMap;

class CommandLineHelpModule;

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
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                          "Root help is not implemented for this output format"));
    }
    writeBasicHelpTopic(context, *this, helpText());
    // TODO: If/when this list becomes long, it may be better to only print
    // "common" commands here, and have a separate topic (e.g.,
    // "help commands") that prints the full list.
    printModuleList(context);
    context.writeTextBlock(
            "For additional help on a command, use '[PROGRAM] help <command>'");
    writeSubTopicList(context,
                      "\nAdditional help is available on the following topics:");
    context.writeTextBlock(
            "To access the help, use '[PROGRAM] help <topic>'.");
}

void RootHelpTopic::printModuleList(const HelpWriterContext &context) const
{
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
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
 * ModuleHelpTopic declaration
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
                        const CommandLineHelpModule      &helpModule)
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
        const CommandLineHelpModule      &helpModule_;

        GMX_DISALLOW_COPY_AND_ASSIGN(ModuleHelpTopic);
};

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
         * \param[in] tag     Unique tag for the module (gmx-something).
         * \param[in] module  Module for which the help should be exported.
         */
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module) = 0;
        /*! \brief
         * Called after all modules have been exported.
         *
         * Can close files opened in startModuleExport(), write footers to them
         * etc.
         */
        virtual void finishModuleExport() = 0;
};

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
        virtual void startModuleExport() {}
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module);
        virtual void finishModuleExport() {}
};

void HelpExportMan::exportModuleHelp(const std::string                &tag,
                                     const CommandLineModuleInterface &module)
{
    File file("man1/" + tag + ".1", "w");

    // TODO: It would be nice to remove the VERSION prefix from the version
    // string to make it shorter.
    file.writeLine(formatString(".TH %s 1 \"\" \"%s\" \"GROMACS Manual\"\n",
                                tag.c_str(),
                                GromacsVersion()));
    file.writeLine(".SH NAME");
    file.writeLine(formatString("%s - %s", tag.c_str(),
                                module.shortDescription()));
    file.writeLine();

    CommandLineHelpContext context(&file, eHelpOutputFormat_Man);
    std::string            displayName(tag);
    std::replace(displayName.begin(), displayName.end(), '-', ' ');
    context.setModuleDisplayName(displayName);
    module.writeHelp(context);

    file.writeLine(".SH SEE ALSO");
    file.writeLine(".BR gromacs(7)");
    file.writeLine();
    file.writeLine("More information about \\fBGROMACS\\fR is available at <\\fIhttp://www.gromacs.org/\\fR>.");
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
        HelpExportHtml();

        virtual void startModuleExport();
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module);
        virtual void finishModuleExport();

    private:
        void writeHtmlHeader(File *file, const std::string &title) const;
        void writeHtmlFooter(File *file) const;

        boost::scoped_ptr<File>  byNameFile_;
        HelpLinks                links_;
};

HelpExportHtml::HelpExportHtml()
{
    char *linksFilename = low_gmxlibfn("links.dat", FALSE, FALSE);
    if (linksFilename != NULL)
    {
        scoped_ptr_sfree guard(linksFilename);
        File             linksFile(linksFilename, "r");
        std::string      line;
        while (linksFile.readLine(&line))
        {
            links_.addLink(line, "../online/" + line);
        }
    }
}

void HelpExportHtml::startModuleExport()
{
    byNameFile_.reset(new File("byname.html", "w"));
    writeHtmlHeader(byNameFile_.get(), "GROMACS Programs by Name");
    byNameFile_->writeLine("<H3>GROMACS Programs Alphabetically</H3>");
}

void HelpExportHtml::exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module)
{
    File file(tag + ".html", "w");
    writeHtmlHeader(&file, tag);

    CommandLineHelpContext context(&file, eHelpOutputFormat_Html);
    std::string            displayName(tag);
    std::replace(displayName.begin(), displayName.end(), '-', ' ');
    context.setModuleDisplayName(displayName);
    context.setLinks(links_);
    module.writeHelp(context);

    writeHtmlFooter(&file);
    file.close();

    byNameFile_->writeLine(formatString("<a href=\"%s.html\">%s</a> - %s<br>",
                                        tag.c_str(), displayName.c_str(),
                                        module.shortDescription()));
}

void HelpExportHtml::finishModuleExport()
{
    writeHtmlFooter(byNameFile_.get());
    byNameFile_->close();
}

void HelpExportHtml::writeHtmlHeader(File *file, const std::string &title) const
{
    file->writeLine("<HTML>");
    file->writeLine("<HEAD>");
    file->writeLine(formatString("<TITLE>%s</TITLE>", title.c_str()));
    file->writeLine("<LINK rel=stylesheet href=\"../online/style.css\" type=\"text/css\">");
    file->writeLine("<BODY text=\"#000000\" bgcolor=\"#FFFFFF\" link=\"#0000FF\" vlink=\"#990000\" alink=\"#FF0000\">");
    file->writeLine("<TABLE WIDTH=\"98%%\" NOBORDER><TR>");
    file->writeLine("<TD WIDTH=400><TABLE WIDTH=400 NOBORDER>");
    file->writeLine("<TD WIDTH=116>");
    file->writeLine("<A HREF=\"http://www.gromacs.org/\">"
                    "<IMG SRC=\"../images/gmxlogo_small.jpg\" BORDER=0>"
                    "</A>");
    file->writeLine("</TD>");
    file->writeLine(formatString("<TD ALIGN=LEFT VALIGN=TOP WIDTH=280>"
                                 "<BR><H2>%s</H2>", title.c_str()));
    file->writeLine("<FONT SIZE=-1><A HREF=\"../online.html\">Main Table of Contents</A></FONT>");
    file->writeLine("</TD>");
    file->writeLine("</TABLE></TD>");
    file->writeLine("<TD WIDTH=\"*\" ALIGN=RIGHT VALIGN=BOTTOM>");
    file->writeLine(formatString("<P><B>%s</B>", GromacsVersion()));
    file->writeLine("</TD>");
    file->writeLine("</TR></TABLE>");
    file->writeLine("<HR>");
}

void HelpExportHtml::writeHtmlFooter(File *file) const
{
    file->writeLine("<P>");
    file->writeLine("<HR>");
    file->writeLine("<DIV ALIGN=RIGHT><FONT SIZE=\"-1\">");
    file->writeLine("<A HREF=\"http://www.gromacs.org\">http://www.gromacs.org</A><BR>");
    file->writeLine("</FONT></DIV>");
    file->writeLine("</BODY>");
    file->writeLine("</HTML>");
}

}   // namespace

/********************************************************************
 * CommandLineHelpModule
 */

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
        //! Sets whether hidden options will be shown in help.
        void setShowHidden(bool bHidden) { bHidden_ = bHidden; }
        /*! \brief
         * Sets an override to show the help for the given module.
         *
         * If called, the help module directly prints the help for the given
         * module when called, skipping any other processing.
         */
        void setModuleOverride(const CommandLineModuleInterface &module)
        {
            moduleOverride_ = &module;
        }

        //! Returns the context object for help output.
        const CommandLineHelpContext &context() const
        {
            return *context_;
        }

        virtual const char *name() const { return "help"; }
        virtual const char *shortDescription() const
        {
            return "Print help information";
        }

        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const CommandLineHelpContext &context) const;

    private:
        void exportHelp(HelpExportInterface *exporter) const;

        boost::scoped_ptr<RootHelpTopic>  rootTopic_;
        const CommandLineModuleMap       &modules_;

        CommandLineHelpContext           *context_;
        const CommandLineModuleInterface *moduleOverride_;
        bool                              bHidden_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineHelpModule);
};

CommandLineHelpModule::CommandLineHelpModule(const CommandLineModuleMap &modules)
    : rootTopic_(new RootHelpTopic(modules)), modules_(modules),
      context_(NULL), moduleOverride_(NULL), bHidden_(false)
{
}

void CommandLineHelpModule::addTopic(HelpTopicPointer topic)
{
    rootTopic_->addSubTopic(move(topic));
}

int CommandLineHelpModule::run(int argc, char *argv[])
{
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
            exporter.reset(new HelpExportMan);
        }
        else if (exportFormat == "html")
        {
            exporter.reset(new HelpExportHtml);
        }
        else
        {
            GMX_THROW(NotImplementedError("This help format is not implemented"));
        }
        exportHelp(exporter.get());
        return 0;
    }

    boost::scoped_ptr<CommandLineHelpContext> context(
            new CommandLineHelpContext(&File::standardOutput(),
                                       eHelpOutputFormat_Console));
    context->setShowHidden(bHidden_);
    context_ = context.get();
    if (moduleOverride_ != NULL)
    {
        ModuleHelpTopic(*moduleOverride_, *this).writeHelp(context->writerContext());
        return 0;
    }

    HelpManager       helpManager(*rootTopic_, context->writerContext());
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

void CommandLineHelpModule::exportHelp(HelpExportInterface *exporter) const
{
    // TODO: Would be nicer to have the file names supplied by the build system
    // and/or export a list of files from here.
    const char *const program =
        ProgramInfo::getInstance().invariantProgramName().c_str();

    exporter->startModuleExport();
    CommandLineModuleMap::const_iterator module;
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        if (module->second->shortDescription() != NULL)
        {
            const char *const moduleName = module->first.c_str();
            std::string       tag(formatString("%s-%s", program, moduleName));
            exporter->exportModuleHelp(tag, *module->second);
        }
    }
    exporter->finishModuleExport();
}

namespace
{

/********************************************************************
 * ModuleHelpTopic implementation
 */

void ModuleHelpTopic::writeHelp(const HelpWriterContext & /*context*/) const
{
    module_.writeHelp(helpModule_.context());
}

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
            const HelpOutputFormat format = context.writerContext().outputFormat();
            const char            *type;
            switch (format)
            {
                case eHelpOutputFormat_Console:
                    type = "help";
                    break;
                case eHelpOutputFormat_Man:
                    type = "nroff";
                    break;
                case eHelpOutputFormat_Html:
                    type = "html";
                    break;
                default:
                    GMX_THROW(NotImplementedError(
                                      "Command-line help is not implemented for this output format"));
            }
            char *argv[4];
            int   argc = 3;
            // TODO: The constness should not be cast away.
            argv[0] = const_cast<char *>(name_);
            argv[1] = const_cast<char *>("-man");
            argv[2] = const_cast<char *>(type);
            argv[3] = NULL;
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
    HelpTopicPointer helpTopic(new ModuleHelpTopic(*module, *helpModule_));
    modules_.insert(std::make_pair(std::string(module->name()),
                                   move(module)));
    helpModule_->addTopic(move(helpTopic));
}

void CommandLineModuleManager::Impl::ensureHelpModuleExists()
{
    if (helpModule_ == NULL)
    {
        helpModule_ = new CommandLineHelpModule(modules_);
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
            programInfo_.setDisplayName(
                    programInfo_.realBinaryName() + "-" + moduleIter->first);
            *argc -= argcForWrapper;
            *argv += argcForWrapper;
            // After this point, argc and argv are the same independent of
            // which path is taken: (*argv)[0] is the module name.
        }
    }
    if (module != NULL)
    {
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
        manager.impl_->singleModule_ = module;
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

} // namespace gmx
