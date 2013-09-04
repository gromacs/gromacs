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

#include <map>
#include <string>
#include <utility>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/network.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineparser.h"
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

/*! \internal \brief
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

    // TODO: Add date? The current version string contains the date, except for
    // releases; with the date, it is too long to fit in a 80-column footer
    // together with the date.
    file.writeLine(formatString(".TH %s 1 \"\" \"%s\" \"Gromacs Manual\"\n",
                                tag.c_str(),
                                GromacsVersion()));
    file.writeLine(".SH NAME");
    file.writeLine(formatString("%s - %s", tag.c_str(),
                                module.shortDescription()));
    file.writeLine();

    HelpWriterContext context(&file, eHelpOutputFormat_Man);
    module.writeHelp(context);

    file.writeLine(".SH SEE ALSO");
    file.writeLine(".BR gromacs(7)");
    file.writeLine();
    file.writeLine("More information about \\fBGROMACS\\fR is available at <\\fIhttp://www.gromacs.org/\\fR>.");

    file.close();
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
        virtual void startModuleExport() {}
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module);
        virtual void finishModuleExport() {}
};

void HelpExportHtml::exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module)
{
    File file(tag + ".html", "w");

    file.writeLine("<HTML>");
    file.writeLine("<HEAD>");
    file.writeLine(formatString("<TITLE>%s</TITLE>", tag.c_str()));
    file.writeLine("<LINK rel=stylesheet href=\"style.css\" type=\"text/css\">");
    file.writeLine("<BODY text=\"#000000\" bgcolor=\"#FFFFFF\" link=\"#0000FF\" vlink=\"#990000\" alink=\"#FF0000\">");
    file.writeLine("<TABLE WIDTH=\"98%%\" NOBORDER><TR>");
    file.writeLine("<TD WIDTH=400><TABLE WIDTH=400 NOBORDER>");
    file.writeLine("<TD WIDTH=116>");
    file.writeLine("<A HREF=\"http://www.gromacs.org/\">"
                   "<IMG SRC=\"../images/gmxlogo_small.png\" BORDER=0>"
                   "</A>");
    file.writeLine("</TD>");
    file.writeLine(formatString("<TD ALIGN=LEFT VALIGN=TOP WIDTH=280>"
                                "<BR><H2>%s</H2>", tag.c_str()));
    file.writeLine("<FONT SIZE=-1><A HREF=\"../online.html\">Main Table of Contents</A></FONT>");
    file.writeLine("</TD>");
    file.writeLine("</TABLE></TD>");
    file.writeLine("<TD WIDTH=\"*\" ALIGN=RIGHT VALIGN=BOTTOM>");
    file.writeLine(formatString("<P><B>%s</B>", GromacsVersion()));
    file.writeLine("</TD>");
    file.writeLine("</TR></TABLE>");
    file.writeLine("<HR>");

    HelpWriterContext context(&file, eHelpOutputFormat_Html);
    module.writeHelp(context);

    file.writeLine("<P>");
    file.writeLine("<HR>");
    file.writeLine("<DIV ALIGN=RIGHT><FONT SIZE=\"-1\">");
    file.writeLine("<A HREF=\"http://www.gromacs.org\">http://www.gromacs.org</A><BR>");
    file.writeLine("<A HREF=\"mailto:gromacs@gromacs.org\">gromacs@gromacs.org</A>");
    file.writeLine("</FONT></DIV>");
    file.writeLine("</BODY>");
    file.writeLine("</HTML>");

    file.close();
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
        virtual void startModuleExport();
        virtual void exportModuleHelp(const std::string                &tag,
                                      const CommandLineModuleInterface &module);
        virtual void finishModuleExport();

    private:
        boost::scoped_ptr<File>  bashFile_;
        boost::scoped_ptr<File>  cshFile_;
        boost::scoped_ptr<File>  zshFile_;
        std::vector<std::string> modules_;
};

void HelpExportCompletion::startModuleExport()
{
    bashFile_.reset(new File("completion.bash", "w"));
    bashFile_->writeLine("shopt -s extglob");
    cshFile_.reset(new File("completion.csh", "w"));
    zshFile_.reset(new File("completion.zsh", "w"));
}

void HelpExportCompletion::exportModuleHelp(const std::string                 & /*tag*/,
                                            const CommandLineModuleInterface &module)
{
    modules_.push_back(module.name());
    {
        HelpWriterContext context(bashFile_.get(), eHelpOutputFormat_CompletionBash);
        context.setModuleDisplayName(std::string("gmx_") + module.name());
        module.writeHelp(context);
    }
    {
        HelpWriterContext context(cshFile_.get(), eHelpOutputFormat_CompletionCsh);
        module.writeHelp(context);
    }
    {
        HelpWriterContext context(zshFile_.get(), eHelpOutputFormat_CompletionZsh);
        module.writeHelp(context);
    }
}

void HelpExportCompletion::finishModuleExport()
{
    const char *const programName = ProgramInfo::getInstance().programName().c_str();

    bashFile_->writeLine("_gmx_compl() {");
    bashFile_->writeLine("local i c m");
    bashFile_->writeLine("COMPREPLY=()");
    bashFile_->writeLine("unset COMP_WORDS[0]");
    bashFile_->writeLine("for ((i=1;i<COMP_CWORD;++i)) ; do");
    bashFile_->writeLine("if [[ \"${COMP_WORDS[i]}\" != -* ]]; then break ; fi");
    bashFile_->writeLine("unset COMP_WORDS[i]");
    bashFile_->writeLine("done");
    bashFile_->writeLine("if (( i == COMP_CWORD )); then");
    bashFile_->writeLine("c=${COMP_WORDS[COMP_CWORD]}");
    std::string gmxCompletions("-h -quiet -version -nocopyright");
    for (std::vector<std::string>::const_iterator i = modules_.begin();
         i != modules_.end(); ++i)
    {
        gmxCompletions.append(" ");
        gmxCompletions.append(*i);
    }
    bashFile_->writeLine("COMPREPLY=( $(compgen -W '" + gmxCompletions + "' -- $c) )");
    bashFile_->writeLine("return 0");
    bashFile_->writeLine("fi");
    bashFile_->writeLine("m=${COMP_WORDS[i]}");
    bashFile_->writeLine("COMP_WORDS=( \"${COMP_WORDS[@]}\" )");
    bashFile_->writeLine("COMP_CWORD=$((COMP_CWORD-i))");
    bashFile_->writeLine("case \"$m\" in");
    for (std::vector<std::string>::const_iterator i = modules_.begin();
         i != modules_.end(); ++i)
    {
        const char *const name = i->c_str();
        bashFile_->writeLine(formatString("%s) _gmx_%s_compl ;;", name, name));
    }
    bashFile_->writeLine("esac }");
    bashFile_->writeLine(formatString("complete -F _gmx_compl %s", programName));
    bashFile_->close();

    cshFile_->close();

    zshFile_->close();
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
        else if (exportFormat == "completion")
        {
            exporter.reset(new HelpExportCompletion);
        }
        else
        {
            GMX_RELEASE_ASSERT(false,
                               "Export format list inconsistent with the implementation");
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
    return 0;
}

void CommandLineHelpModule::writeHelp(const HelpWriterContext &context) const
{
    // TODO: Implement.
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        return;
    }
    context.writeTextBlock(
            "Usage: [PROGRAM] help [<command>|<topic> [<subtopic> [...]]]");
    // TODO: More information.
}

void CommandLineHelpModule::exportHelp(HelpExportInterface *exporter) const
{
    // TODO: Would be nicer to have the file names supplied by the build system
    // and/or export a list of files from here.
    const char *const program =
        ProgramInfo::getInstance().programName().c_str();

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
        virtual void writeHelp(const HelpWriterContext &context) const
        {
            char *argv[3];
            int   argc = 1;
            // TODO: The constness should not be cast away.
            argv[0] = const_cast<char *>(name_);
            if (context.outputFormat() == eHelpOutputFormat_Console)
            {
                argv[1] = const_cast<char *>("-h");
                argc    = 2;
            }
            else
            {
                const char *type;
                switch (context.outputFormat())
                {
                    case eHelpOutputFormat_Man:
                        type = "nroff";
                        break;
                    case eHelpOutputFormat_Html:
                        type = "html";
                        break;
                    case eHelpOutputFormat_CompletionBash:
                        type = "completion-bash";
                        break;
                    case eHelpOutputFormat_CompletionCsh:
                        type = "completion-csh";
                        break;
                    case eHelpOutputFormat_CompletionZsh:
                        type = "completion-zsh";
                        break;
                    default:
                        GMX_THROW(NotImplementedError(
                                          "Command-line help is not implemented for this output format"));
                }
                argv[1] = const_cast<char *>("-man");
                argv[2] = const_cast<char *>(type);
                argc    = 3;
            }
            HelpWriterGlobalContext global(context);
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
    bool bVersion   = false;
    bool bCopyright = true;
    // TODO: Print the common options into the help.
    // TODO: It would be nice to propagate at least the -quiet option to
    // the modules so that they can also be quiet in response to this.
    // TODO: Consider handling -h and related options here instead of in the
    // modules (also -hidden needs to be transfered here to make that work).
    // That would mean that with -h, all module-specific options would get
    // ignored.  This means that the help output would not depend on the
    // command line, but would always show the default values (making it
    // possible to simplify it further), but also that mdrun -h could not be
    // used for option validation in g_tune_pme.
    Options options(NULL, NULL);
    if (module == NULL)
    {
        options.addOption(BooleanOption("h").store(&bHelp));
    }
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
    else
    {
        // In single-module mode, recognize the common options also after the
        // module name.
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
    if (module == NULL)
    {
        *argc = 1;
        return helpModule_;
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
    GMX_ASSERT(impl_->modules_.find(module->name()) == impl_->modules_.end(),
               "Attempted to register a duplicate module name");
    HelpTopicPointer helpTopic(new ModuleHelpTopic(*module));
    impl_->modules_.insert(std::make_pair(std::string(module->name()),
                                          move(module)));
    addHelpTopic(move(helpTopic));
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
    if (impl_->helpModule_ == NULL)
    {
        impl_->helpModule_ = new CommandLineHelpModule(impl_->modules_);
        addModule(CommandLineModulePointer(impl_->helpModule_));
    }
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
