/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
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
 * Implements gmx::CommandLineHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "cmdlinehelpwriter.h"

#include <string>

#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/********************************************************************
 * OptionsFormatterInterface
 */

/*! \internal \brief
 * Interface for output format specific formatting of options.
 *
 * \see OptionsFilter
 *
 * \ingroup module_commandline
 */
class OptionsFormatterInterface
{
    public:
        virtual ~OptionsFormatterInterface() {}

        //! Formats the description text block for a section.
        virtual void formatDescription(const HelpWriterContext &context,
                                       const Options           &section) = 0;
        //! Formats a single file option.
        virtual void formatFileOption(const HelpWriterContext  &context,
                                      const FileNameOptionInfo &option) = 0;
        //! Formats a single non-file, non-selection option.
        virtual void formatOption(const HelpWriterContext &context,
                                  const OptionInfo        &option) = 0;
        //! Formats a single selection option.
        virtual void formatSelectionOption(const HelpWriterContext &context,
                                           const OptionInfo        &option) = 0;
};

/********************************************************************
 * OptionsFilter
 */

/*! \internal \brief
 * Output format independent processing of options.
 *
 * Together with code in CommandLineHelpWriter::writeHelp(), this class
 * implements the common logic for writing out the help.
 * An object implementing the OptionsFormatterInterface must be provided to the
 * constructor, and does the actual formatting that is specific to the output
 * format.
 *
 * \ingroup module_commandline
 */
class OptionsFilter : public OptionsVisitor
{
    public:
        //! Specifies the type of output that formatSelected() produces.
        enum FilterType
        {
            eSelectDescriptions,
            eSelectFileOptions,
            eSelectSelectionOptions,
            eSelectOtherOptions
        };

        /*! \brief
         * Creates a new filtering object.
         *
         * \param[in] context   Help context to use to write the help.
         * \param[in] formatter Output format specific formatter.
         *
         * Does not throw.
         */
        OptionsFilter(const HelpWriterContext   &context,
                      OptionsFormatterInterface *formatter)
            : context_(context), formatter_(*formatter),
              filterType_(eSelectOtherOptions), bShowHidden_(false),
              bDidOutput_(false)
        {
        }

        //! Sets whether hidden options will be shown.
        void setShowHidden(bool bShowHidden)
        {
            bShowHidden_ = bShowHidden;
        }

        //! Formats selected options using the formatter.
        void formatSelected(FilterType type, const Options &options);

        virtual void visitSubSection(const Options &section);
        virtual void visitOption(const OptionInfo &option);

    private:
        const HelpWriterContext        &context_;
        OptionsFormatterInterface      &formatter_;
        FilterType                      filterType_;
        bool                            bShowHidden_;
        bool                            bDidOutput_;
};

void OptionsFilter::formatSelected(FilterType type, const Options &options)
{
    filterType_ = type;
    bDidOutput_ = false;
    visitSubSection(options);
    if (bDidOutput_)
    {
        context_.outputFile().writeLine();
    }
}

void OptionsFilter::visitSubSection(const Options &section)
{
    if (filterType_ == eSelectDescriptions)
    {
        if (bDidOutput_)
        {
            context_.outputFile().writeLine();
        }
        formatter_.formatDescription(context_, section);
        bDidOutput_ = true;
    }

    OptionsIterator iterator(section);
    iterator.acceptSubSections(this);
    iterator.acceptOptions(this);
}

void OptionsFilter::visitOption(const OptionInfo &option)
{
    if (!bShowHidden_ && option.isHidden())
    {
        return;
    }
    if (option.isType<FileNameOptionInfo>())
    {
        if (filterType_ == eSelectFileOptions)
        {
            formatter_.formatFileOption(context_,
                                        *option.toType<FileNameOptionInfo>());
            bDidOutput_ = true;
        }
        return;
    }
    if (option.isType<SelectionFileOptionInfo>()
        || option.isType<SelectionOptionInfo>())
    {
        if (filterType_ == eSelectSelectionOptions)
        {
            formatter_.formatSelectionOption(context_, option);
            bDidOutput_ = true;
        }
        return;
    }
    if (filterType_ == eSelectOtherOptions)
    {
        formatter_.formatOption(context_, option);
        bDidOutput_ = true;
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

std::string
descriptionWithOptionDetails(const CommonFormatterData &common,
                             const OptionInfo          &option)
{
    std::string             description(option.description());

    const DoubleOptionInfo *doubleOption = option.toType<DoubleOptionInfo>();
    if (doubleOption != NULL && doubleOption->isTime())
    {
        description = replaceAll(description, "%t", common.timeUnit);
    }

    const StringOptionInfo *stringOption = option.toType<StringOptionInfo>();
    if (stringOption != NULL && stringOption->isEnumerated())
    {
        const std::vector<std::string> &allowedValues
            = stringOption->allowedValues();
        description.append(": ");
        for (size_t i = 0; i < allowedValues.size(); ++i)
        {
            if (i > 0)
            {
                description.append(i + 1 < allowedValues.size()
                                   ? ", " : ", or ");
            }
            description.append(allowedValues[i]);
        }
    }
    return description;
}

std::string
fileOptionFlagsAsString(const FileNameOptionInfo &option, bool bAbbrev)
{
    std::string type;
    if (option.isInputOutputFile())
    {
        type = bAbbrev ? "In/Out" : "Input/Output";
    }
    else if (option.isInputFile())
    {
        type = "Input";
    }
    else if (option.isOutputFile())
    {
        type = "Output";
    }
    if (!option.isRequired())
    {
        type += bAbbrev ? ", Opt." : ", Optional";
    }
    if (option.isLibraryFile())
    {
        type += bAbbrev ? ", Lib." : ", Library";
    }
    return type;
}

/********************************************************************
 * OptionsConsoleFormatter
 */

/*! \internal \brief
 * Formatter implementation for console help.
 *
 * \ingroup module_commandline
 */
class OptionsConsoleFormatter : public OptionsFormatterInterface
{
    public:
        //! Creates a helper object for formatting options help for console.
        explicit OptionsConsoleFormatter(const CommonFormatterData &common);

        virtual void formatDescription(const HelpWriterContext &context,
                                       const Options           &section);
        virtual void formatFileOption(const HelpWriterContext  &context,
                                      const FileNameOptionInfo &option);
        virtual void formatOption(const HelpWriterContext &context,
                                  const OptionInfo        &option);
        virtual void formatSelectionOption(const HelpWriterContext &context,
                                           const OptionInfo        &option);

    private:
        const CommonFormatterData &common_;
        TextTableFormatter         fileOptionFormatter_;
        TextTableFormatter         genericOptionFormatter_;
        TextTableFormatter         selectionOptionFormatter_;
};

OptionsConsoleFormatter::OptionsConsoleFormatter(const CommonFormatterData &common)
    : common_(common)
{
    fileOptionFormatter_.addColumn("Option",      6, false);
    fileOptionFormatter_.addColumn("Filename",    12, false);
    fileOptionFormatter_.addColumn("Type",        12, false);
    fileOptionFormatter_.addColumn("Description", 45, true);

    genericOptionFormatter_.addColumn("Option",      12, false);
    genericOptionFormatter_.addColumn("Type",         6, false);
    genericOptionFormatter_.addColumn("Value",        6, false);
    genericOptionFormatter_.addColumn("Description", 51, true);

    selectionOptionFormatter_.addColumn("Selection",   10, false);
    selectionOptionFormatter_.addColumn("Description", 67, true);
}

void OptionsConsoleFormatter::formatDescription(
        const HelpWriterContext &context, const Options &section)
{
    if (!section.description().empty())
    {
        File              &file  = context.outputFile();
        const std::string &title = section.title();
        if (!title.empty())
        {
            file.writeLine(title);
            file.writeLine();
        }
        context.writeTextBlock(section.description());
    }
}

void OptionsConsoleFormatter::formatFileOption(
        const HelpWriterContext &context, const FileNameOptionInfo &option)
{
    int firstShortValue = 0;  // The first value after which the type fits.
    int firstLongValue  = -1; // First value that overlaps description column.
    int lastLongValue   = -1; // Last value like the above.

    // Get the values to write and check where text overflows the columns.
    fileOptionFormatter_.clear();
    std::string name(formatString("-%s", option.name().c_str()));
    fileOptionFormatter_.addColumnLine(0, name);
    for (int i = 0; i < option.valueCount() || i == 0; ++i)
    {
        std::string value;
        if (option.valueCount() == 0
            || (option.valueCount() == 1 && option.formatValue(0).empty()))
        {
            value = option.formatDefaultValueIfSet();
        }
        else
        {
            value = option.formatValue(i);
        }
        fileOptionFormatter_.addColumnLine(1, value);
        if (value.length() > 12U && i == firstShortValue)
        {
            firstShortValue = i + 1;
        }
        if (value.length() > 25U)
        {
            if (firstLongValue == -1)
            {
                firstLongValue = i;
            }
            lastLongValue = i;
        }
    }
    std::string type(fileOptionFlagsAsString(option, true));
    bool        bLongType = (type.length() > 12U);
    fileOptionFormatter_.addColumnLine(2, type);
    fileOptionFormatter_.addColumnHelpTextBlock(3, context, option.description());

    // Compute layout.
    if (name.length() > 6U || firstShortValue > 0)
    {
        fileOptionFormatter_.setColumnFirstLineOffset(1, 1);
        // Assume that the name is <20 chars, so that the type fits
        if (firstLongValue >= 0)
        {
            ++firstLongValue;
            ++lastLongValue;
        }
    }
    int firstDescriptionLine = 0;
    if (bLongType)
    {
        firstDescriptionLine = 1;
    }
    fileOptionFormatter_.setColumnFirstLineOffset(3, firstDescriptionLine);
    if (firstLongValue >= 0
        && fileOptionFormatter_.lastColumnLine(3) >= firstLongValue)
    {
        firstDescriptionLine = lastLongValue + 1;
        fileOptionFormatter_.setColumnFirstLineOffset(3, firstDescriptionLine);
    }

    // Do the formatting.
    context.outputFile().writeString(fileOptionFormatter_.formatRow());
}

void OptionsConsoleFormatter::formatOption(
        const HelpWriterContext &context, const OptionInfo &option)
{
    genericOptionFormatter_.clear();
    bool        bIsBool = option.isType<BooleanOptionInfo>();
    std::string name(formatString("-%s%s", bIsBool ? "[no]" : "",
                                  option.name().c_str()));
    genericOptionFormatter_.addColumnLine(0, name);
    genericOptionFormatter_.addColumnLine(1, option.type());
    if (name.length() > 12U)
    {
        genericOptionFormatter_.setColumnFirstLineOffset(1, 1);
    }

    // TODO: Better handling of multiple long values
    std::string values;
    for (int i = 0; i < option.valueCount(); ++i)
    {
        if (i != 0)
        {
            values.append(" ");
        }
        values.append(option.formatValue(i));
    }
    genericOptionFormatter_.addColumnLine(2, values);

    std::string description(descriptionWithOptionDetails(common_, option));
    genericOptionFormatter_.addColumnHelpTextBlock(3, context, description);
    if (values.length() > 6U)
    {
        genericOptionFormatter_.setColumnFirstLineOffset(3, 1);
    }

    context.outputFile().writeString(genericOptionFormatter_.formatRow());
}

void OptionsConsoleFormatter::formatSelectionOption(
        const HelpWriterContext &context, const OptionInfo &option)
{
    File &file = context.outputFile();

    selectionOptionFormatter_.clear();
    std::string name(formatString("-%s", option.name().c_str()));
    selectionOptionFormatter_.addColumnLine(0, name);
    selectionOptionFormatter_.addColumnHelpTextBlock(1, context,
                                                     option.description());
    file.writeString(selectionOptionFormatter_.formatRow());

    TextLineWrapper wrapper;
    wrapper.settings().setLineLength(77);
    wrapper.settings().setFirstLineIndent(4);
    wrapper.settings().setIndent(8);
    wrapper.settings().setContinuationChar('\\');
    // TODO: What to do with selection variables?
    // They are not printed as values for any option.
    for (int i = 0; i < option.valueCount(); ++i)
    {
        std::string value(option.formatValue(i));
        file.writeLine(wrapper.wrapToString(value));
    }
}

/********************************************************************
 * OptionsExportFormatter
 */

/*! \internal \brief
 * Formatter implementation for help export.
 *
 * \ingroup module_commandline
 */
class OptionsExportFormatter : public OptionsFormatterInterface
{
    public:
        //! Creates a helper object for formatting options help for man pages.
        explicit OptionsExportFormatter(const CommonFormatterData &common)
            : common_(common)
        {
        }

        virtual void formatDescription(const HelpWriterContext &context,
                                       const Options           &section);
        virtual void formatFileOption(const HelpWriterContext  &context,
                                      const FileNameOptionInfo &option);
        virtual void formatOption(const HelpWriterContext &context,
                                  const OptionInfo        &option);
        virtual void formatSelectionOption(const HelpWriterContext &context,
                                           const OptionInfo        &option);

    private:
        const CommonFormatterData &common_;
};

void OptionsExportFormatter::formatDescription(
        const HelpWriterContext &context, const Options &section)
{
    if (!section.description().empty())
    {
        context.writeTextBlock(section.description());
    }
}

void OptionsExportFormatter::formatFileOption(
        const HelpWriterContext &context, const FileNameOptionInfo &option)
{
    std::string defaultValue;
    if (option.valueCount() == 0
        || (option.valueCount() == 1 && option.formatValue(0).empty()))
    {
        defaultValue = option.formatDefaultValueIfSet();
    }
    else
    {
        defaultValue = option.formatValue(0);
    }
    std::string description(option.description());
    std::string type(fileOptionFlagsAsString(option, false));
    description.append(formatString(" (%s)", type.c_str()));
    context.writeOptionItem(option.name(),
                            formatString("<%s>", defaultValue.c_str()),
                            description);
}

void OptionsExportFormatter::formatOption(
        const HelpWriterContext &context, const OptionInfo &option)
{
    if (option.isType<BooleanOptionInfo>())
    {
        // Implementation depends on approach.
    }
    else
    {
        std::string description(descriptionWithOptionDetails(common_, option));
        context.writeOptionItem(option.name(),
                                formatString("<%s>", option.type()),
                                description);
    }
}

void OptionsExportFormatter::formatSelectionOption(
        const HelpWriterContext &context, const OptionInfo &option)
{
    const char *args = (option.isType<SelectionFileOptionInfo>()
                        ? "<selection.dat>"
                        : "SELECTION");
    context.writeOptionItem(option.name(), args, option.description());
}

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
        explicit Impl(const Options &options);

        //! Options object to use for generating help.
        const Options          &options_;
        //! Time unit to show in descriptions.
        std::string             timeUnit_;
        //! Whether to write descriptions to output.
        bool                    bShowDescriptions_;
        //! Whether to write hidden options to output.
        bool                    bShowHidden_;
};

CommandLineHelpWriter::Impl::Impl(const Options &options)
    : options_(options), timeUnit_(TimeUnitManager().timeUnitAsString()),
      bShowDescriptions_(false), bShowHidden_(false)
{
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

CommandLineHelpWriter &CommandLineHelpWriter::setShowHidden(bool bSet)
{
    impl_->bShowHidden_ = bSet;
    return *this;
}

CommandLineHelpWriter &CommandLineHelpWriter::setShowDescriptions(bool bSet)
{
    impl_->bShowDescriptions_ = bSet;
    return *this;
}

CommandLineHelpWriter &CommandLineHelpWriter::setTimeUnitString(const char *timeUnit)
{
    impl_->timeUnit_ = timeUnit;
    return *this;
}

void CommandLineHelpWriter::writeHelp(const HelpWriterContext &context)
{
    boost::scoped_ptr<OptionsFormatterInterface> formatter;
    CommonFormatterData common(impl_->timeUnit_.c_str());
    switch (context.outputFormat())
    {
        case eHelpOutputFormat_Console:
            formatter.reset(new OptionsConsoleFormatter(common));
            break;
        default:
            GMX_THROW(NotImplementedError(
                              "Command-line help is not implemented for this output format"));
    }
    OptionsFilter filter(context, formatter.get());
    filter.setShowHidden(impl_->bShowHidden_);

    File &file = context.outputFile();
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        context.writeTitle("Synopsis");
        context.writeTextBlock(formatString("[PROGRAM] %s", impl_->options_.name().c_str()));
        file.writeLine("\n\n");
    }

    if (impl_->bShowDescriptions_)
    {
        context.writeTitle("Description");
        filter.formatSelected(OptionsFilter::eSelectDescriptions,
                              impl_->options_);
    }
    if (context.outputFormat() != eHelpOutputFormat_Console)
    {
        context.writeTitle("Options");
    }
    filter.formatSelected(OptionsFilter::eSelectFileOptions,
                          impl_->options_);
    filter.formatSelected(OptionsFilter::eSelectOtherOptions,
                          impl_->options_);
    filter.formatSelected(OptionsFilter::eSelectSelectionOptions,
                          impl_->options_);
}

} // namespace gmx
