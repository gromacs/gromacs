/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements gmx::CommandLineHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlinehelpwriter.h"

#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "shellcompletions.h"

namespace gmx
{

namespace
{

//! \addtogroup module_commandline
//! \{

/********************************************************************
 * IOptionsFormatter
 */

/*! \brief
 * Interface for output format specific formatting of options.
 *
 * \see OptionsFilter
 */
class IOptionsFormatter
{
    public:
        virtual ~IOptionsFormatter() {}

        //! Formats a single option option.
        virtual void formatOption(const OptionInfo &option) = 0;
};

/********************************************************************
 * OptionsFilter
 */

/*! \brief
 * Output format independent processing of options.
 *
 * Together with code in CommandLineHelpWriter::writeHelp(), this class
 * implements the common logic for writing out the help.
 * An object implementing the IOptionsFormatter must be provided to the
 * constructor, and does the actual formatting that is specific to the output
 * format.
 */
class OptionsFilter : public OptionsVisitor
{
    public:
        //! Specifies the type of output that formatSelected() produces.
        enum FilterType
        {
            eSelectInputFileOptions,
            eSelectInputOutputFileOptions,
            eSelectOutputFileOptions,
            eSelectOtherOptions
        };

        /*! \brief
         * Creates a new filtering object.
         *
         * Does not throw.
         */
        OptionsFilter()
            : formatter_(nullptr), filterType_(eSelectOtherOptions),
              bShowHidden_(false)
        {
        }

        //! Sets whether hidden options will be shown.
        void setShowHidden(bool bShowHidden)
        {
            bShowHidden_ = bShowHidden;
        }

        //! Formats selected options using the formatter.
        void formatSelected(FilterType                 type,
                            IOptionsFormatter         *formatter,
                            const Options             &options);

        virtual void visitSection(const OptionSectionInfo &section);
        virtual void visitOption(const OptionInfo &option);

    private:
        IOptionsFormatter              *formatter_;
        FilterType                      filterType_;
        bool                            bShowHidden_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsFilter);
};

void OptionsFilter::formatSelected(FilterType                 type,
                                   IOptionsFormatter         *formatter,
                                   const Options             &options)
{
    formatter_  = formatter;
    filterType_ = type;
    visitSection(options.rootSection());
}

void OptionsFilter::visitSection(const OptionSectionInfo &section)
{
    OptionsIterator iterator(section);
    iterator.acceptSections(this);
    iterator.acceptOptions(this);
}

void OptionsFilter::visitOption(const OptionInfo &option)
{
    if (!bShowHidden_ && option.isHidden())
    {
        return;
    }
    const FileNameOptionInfo *const fileOption = option.toType<FileNameOptionInfo>();
    if (fileOption != nullptr && fileOption->isInputFile())
    {
        if (filterType_ == eSelectInputFileOptions)
        {
            formatter_->formatOption(option);
        }
        return;
    }
    if (fileOption != nullptr && fileOption->isInputOutputFile())
    {
        if (filterType_ == eSelectInputOutputFileOptions)
        {
            formatter_->formatOption(option);
        }
        return;
    }
    if (fileOption != nullptr && fileOption->isOutputFile())
    {
        if (filterType_ == eSelectOutputFileOptions)
        {
            formatter_->formatOption(option);
        }
        return;
    }
    if (filterType_ == eSelectOtherOptions)
    {
        formatter_->formatOption(option);
        return;
    }
}

/********************************************************************
 * CommonFormatterData
 */

class CommonFormatterData
{
    public:
        explicit CommonFormatterData(const char *timeUnit)
            : timeUnit(timeUnit)
        {
        }

        const char             *timeUnit;
};

/********************************************************************
 * Helper functions
 */

//! Formats option name and value.
void formatOptionNameAndValue(const OptionInfo &option, std::string *name,
                              std::string *value)
{
    *name  = option.name();
    *value = "<" + option.type() + ">";
    if (option.isType<FileNameOptionInfo>())
    {
        // TODO: Make these work also for other option types.
        if (option.maxValueCount() != 1)
        {
            *value += " [...]";
        }
        if (option.minValueCount() == 0)
        {
            *value = "[" + *value + "]";
        }
    }
    if (option.isType<BooleanOptionInfo>())
    {
        *name = "[no]" + *name;
        // Old command-line parser doesn't accept any values for these.
        // value = "[" + value + "]";
        value->clear();
    }
}

//! Formats the default option value as a string.
std::string defaultOptionValue(const OptionInfo &option)
{
    return joinStrings(option.defaultValuesAsStrings(), " ");
}

//! Formats the flags for a file option as a string.
std::string
fileOptionFlagsAsString(const FileNameOptionInfo &option, bool bAbbrev)
{
    std::string type;
    if (!option.isRequired())
    {
        type = bAbbrev ? "Opt." : "Optional";
    }
    if (option.isLibraryFile())
    {
        if (!type.empty())
        {
            type.append(", ");
        }
        type.append(bAbbrev ? "Lib." : "Library");
    }
    return type;
}

//! Formats the description for an option as a string.
std::string
descriptionWithOptionDetails(const CommonFormatterData &common,
                             const OptionInfo          &option)
{
    std::string             description(option.formatDescription());

    const FloatOptionInfo  *floatOption  = option.toType<FloatOptionInfo>();
    const DoubleOptionInfo *doubleOption = option.toType<DoubleOptionInfo>();
    if ((floatOption != nullptr && floatOption->isTime())
        || (doubleOption != nullptr && doubleOption->isTime()))
    {
        // TODO: It could be nicer to have this in basicoptions.cpp.
        description = replaceAll(description, "%t", common.timeUnit);
    }

    return description;
}

/********************************************************************
 * OptionsSynopsisFormatter
 */

/*! \brief
 * Formatter implementation for synopsis.
 */
class SynopsisFormatter : public IOptionsFormatter
{
    public:
        //! Creates a helper object for formatting the synopsis.
        explicit SynopsisFormatter(const HelpWriterContext &context)
            : context_(context), bFormatted_(false), lineLength_(0), indent_(0),
              currentLength_(0)
        {
        }

        //! Starts formatting the synopsis.
        void start(const char *name);
        //! Finishes formatting the synopsis.
        void finish();

        virtual void formatOption(const OptionInfo &option);

    private:
        const HelpWriterContext &context_;
        bool                     bFormatted_;
        int                      lineLength_;
        int                      indent_;
        int                      currentLength_;

        GMX_DISALLOW_COPY_AND_ASSIGN(SynopsisFormatter);
};

void SynopsisFormatter::start(const char *name)
{
    currentLength_ = std::strlen(name) + 1;
    indent_        = std::min(currentLength_, 13);
    TextWriter &file = context_.outputFile();
    switch (context_.outputFormat())
    {
        case eHelpOutputFormat_Console:
            lineLength_ = 78;
            file.writeString(name);
            break;
        case eHelpOutputFormat_Rst:
            bFormatted_ = true;
            lineLength_ = 74;
            indent_    += 4;
            file.writeLine(".. parsed-literal::");
            file.writeLine();
            file.writeString("    ");
            file.writeString(name);
            break;
        default:
            GMX_THROW(NotImplementedError("Synopsis formatting not implemented for this output format"));
    }
}

void SynopsisFormatter::finish()
{
    context_.outputFile().ensureLineBreak();
}

void SynopsisFormatter::formatOption(const OptionInfo &option)
{
    std::string name, value;
    formatOptionNameAndValue(option, &name, &value);
    int         totalLength = name.length() + 4;
    std::string fullOptionText
        = formatString(" [%s-%s", bFormatted_ ? ":strong:`" : "", name.c_str());
    if (!value.empty())
    {
        fullOptionText.append(bFormatted_ ? "` :emphasis:`" : " ");
        fullOptionText.append(value);
        totalLength += value.length() + 1;
    }
    fullOptionText.append(bFormatted_ ? "`]" : "]");

    TextWriter &file = context_.outputFile();
    currentLength_ += totalLength;
    if (currentLength_ >= lineLength_)
    {
        file.writeString(formatString("\n%*c", indent_ - 1, ' '));
        currentLength_ = indent_ - 1 + totalLength;
    }
    file.writeString(fullOptionText);
}

/********************************************************************
 * OptionsListFormatter
 */

/*! \brief
 * Formatter implementation for help export.
 */
class OptionsListFormatter : public IOptionsFormatter
{
    public:
        //! Creates a helper object for formatting options.
        OptionsListFormatter(const HelpWriterContext   &context,
                             const CommonFormatterData &common,
                             const char                *title);

        //! Initiates a new section in the options list.
        void startSection(const char *header)
        {
            header_     = header;
            bDidOutput_ = false;
        }
        //! Finishes a section in the options list.
        void finishSection()
        {
            if (bDidOutput_)
            {
                context_.writeOptionListEnd();
            }
        }

        virtual void formatOption(const OptionInfo &option);

    private:
        void writeSectionStartIfNecessary()
        {
            if (title_ != nullptr)
            {
                context_.writeTitle(title_);
                title_ = nullptr;
            }
            if (!bDidOutput_)
            {
                if (header_ != nullptr)
                {
                    context_.paragraphBreak();
                    context_.writeTextBlock(header_);
                    context_.paragraphBreak();
                }
                context_.writeOptionListStart();
            }
            bDidOutput_ = true;
        }

        const HelpWriterContext               &context_;
        const CommonFormatterData             &common_;
        const char                            *title_;
        const char                            *header_;
        bool                                   bDidOutput_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsListFormatter);
};

OptionsListFormatter::OptionsListFormatter(
        const HelpWriterContext   &context,
        const CommonFormatterData &common,
        const char                *title)
    : context_(context), common_(common),
      title_(title), header_(nullptr), bDidOutput_(false)
{
}

void OptionsListFormatter::formatOption(const OptionInfo &option)
{
    writeSectionStartIfNecessary();
    std::string               name, value;
    formatOptionNameAndValue(option, &name, &value);
    std::string               defaultValue(defaultOptionValue(option));
    std::string               info;
    const FileNameOptionInfo *fileOption = option.toType<FileNameOptionInfo>();
    if (fileOption != nullptr)
    {
        const bool bAbbrev = (context_.outputFormat() == eHelpOutputFormat_Console);
        info = fileOptionFlagsAsString(*fileOption, bAbbrev);
    }
    std::string description(descriptionWithOptionDetails(common_, option));
    context_.writeOptionItem("-" + name, value, defaultValue, info, description);
}

//! \}

}   // namespace

/********************************************************************
 * CommandLineHelpWriter::Impl
 */

/*! \internal \brief
 * Private implementation class for CommandLineHelpWriter.
 *
 * \ingroup module_commandline
 */
class CommandLineHelpWriter::Impl
{
    public:
        //! Sets the Options object to use for generating help.
        explicit Impl(const Options &options)
            : options_(options)
        {
        }

        //! Format the list of known issues.
        void formatBugs(const HelpWriterContext &context);

        //! Options object to use for generating help.
        const Options               &options_;
        //! Help text.
        std::string                  helpText_;
        //! List of bugs/knows issues.
        ArrayRef<const char *const>  bugs_;
};

void CommandLineHelpWriter::Impl::formatBugs(const HelpWriterContext &context)
{
    if (bugs_.empty())
    {
        return;
    }
    context.writeTitle("Known Issues");
    ArrayRef<const char *const>::const_iterator i;
    for (i = bugs_.begin(); i != bugs_.end(); ++i)
    {
        const char *const       bug = *i;
        context.writeTextBlock(formatString("* %s", bug));
    }
}


/********************************************************************
 * CommandLineHelpWriter
 */

CommandLineHelpWriter::CommandLineHelpWriter(const Options &options)
    : impl_(new Impl(options))
{
}

CommandLineHelpWriter::~CommandLineHelpWriter()
{
}

CommandLineHelpWriter &
CommandLineHelpWriter::setHelpText(const std::string &help)
{
    impl_->helpText_ = help;
    return *this;
}

CommandLineHelpWriter &
CommandLineHelpWriter::setHelpText(const ArrayRef<const char *const> &help)
{
    impl_->helpText_ = joinStrings(help, "\n");
    return *this;
}

CommandLineHelpWriter &
CommandLineHelpWriter::setKnownIssues(const ArrayRef<const char *const> &bugs)
{
    impl_->bugs_ = bugs;
    return *this;
}

void CommandLineHelpWriter::writeHelp(const CommandLineHelpContext &context)
{
    if (context.isCompletionExport())
    {
        context.shellCompletionWriter().writeModuleCompletions(
                context.moduleDisplayName(), impl_->options_);
        return;
    }
    const HelpWriterContext &writerContext = context.writerContext();
    OptionsFilter            filter;
    filter.setShowHidden(context.showHidden());

    {
        writerContext.writeTitle("Synopsis");
        SynopsisFormatter synopsisFormatter(writerContext);
        synopsisFormatter.start(context.moduleDisplayName());
        filter.formatSelected(OptionsFilter::eSelectInputFileOptions,
                              &synopsisFormatter, impl_->options_);
        filter.formatSelected(OptionsFilter::eSelectInputOutputFileOptions,
                              &synopsisFormatter, impl_->options_);
        filter.formatSelected(OptionsFilter::eSelectOutputFileOptions,
                              &synopsisFormatter, impl_->options_);
        filter.formatSelected(OptionsFilter::eSelectOtherOptions,
                              &synopsisFormatter, impl_->options_);
        synopsisFormatter.finish();
    }

    if (!impl_->helpText_.empty())
    {
        writerContext.writeTitle("Description");
        writerContext.writeTextBlock(impl_->helpText_);
    }
    CommonFormatterData    common(TimeUnitManager().timeUnitAsString());
    OptionsListFormatter   formatter(writerContext, common, "Options");
    formatter.startSection("Options to specify input files:");
    filter.formatSelected(OptionsFilter::eSelectInputFileOptions,
                          &formatter, impl_->options_);
    formatter.finishSection();
    formatter.startSection("Options to specify input/output files:");
    filter.formatSelected(OptionsFilter::eSelectInputOutputFileOptions,
                          &formatter, impl_->options_);
    formatter.finishSection();
    formatter.startSection("Options to specify output files:");
    filter.formatSelected(OptionsFilter::eSelectOutputFileOptions,
                          &formatter, impl_->options_);
    formatter.finishSection();
    formatter.startSection("Other options:");
    filter.formatSelected(OptionsFilter::eSelectOtherOptions,
                          &formatter, impl_->options_);
    formatter.finishSection();

    impl_->formatBugs(writerContext);
}

} // namespace gmx
