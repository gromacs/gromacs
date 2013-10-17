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

#include "gromacs/commandline/cmdlinehelpcontext.h"
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

//! \addtogroup module_commandline
//! \{

/********************************************************************
 * OptionsFormatterInterface
 */

/*! \brief
 * Interface for output format specific formatting of options.
 *
 * \see OptionsFilter
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

/*! \brief
 * Output format independent processing of options.
 *
 * Together with code in CommandLineHelpWriter::writeHelp(), this class
 * implements the common logic for writing out the help.
 * An object implementing the OptionsFormatterInterface must be provided to the
 * constructor, and does the actual formatting that is specific to the output
 * format.
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
              bShowHidden_(false), filterType_(eSelectOtherOptions),
              title_(NULL), bDidOutput_(false)
        {
        }

        //! Sets whether hidden options will be shown.
        void setShowHidden(bool bShowHidden)
        {
            bShowHidden_ = bShowHidden;
        }

        //! Formats selected options using the formatter.
        bool formatSelected(FilterType type, const Options &options,
                            const char *title);

        virtual void visitSubSection(const Options &section);
        virtual void visitOption(const OptionInfo &option);

    private:
        const HelpWriterContext        &context_;
        OptionsFormatterInterface      &formatter_;
        bool                            bShowHidden_;
        FilterType                      filterType_;
        const char                     *title_;
        bool                            bDidOutput_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsFilter);
};

bool OptionsFilter::formatSelected(FilterType type, const Options &options,
                                   const char *title)
{
    filterType_ = type;
    title_      = title;
    bDidOutput_ = false;
    if (type != eSelectDescriptions)
    {
        context_.writeOptionListStart();
    }
    visitSubSection(options);
    if (type != eSelectDescriptions)
    {
        context_.writeOptionListEnd();
    }
    if (bDidOutput_)
    {
        context_.outputFile().writeLine();
    }
    return bDidOutput_;
}

void OptionsFilter::visitSubSection(const Options &section)
{
    if (filterType_ == eSelectDescriptions)
    {
        if (bDidOutput_)
        {
            context_.outputFile().writeLine();
        }
        else if (title_ != NULL)
        {
            context_.writeTitle(title_);
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
            if (!bDidOutput_ && title_ != NULL)
            {
                context_.writeTitle(title_);
            }
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
            if (!bDidOutput_ && title_ != NULL)
            {
                context_.writeTitle(title_);
            }
            formatter_.formatSelectionOption(context_, option);
            bDidOutput_ = true;
        }
        return;
    }
    if (filterType_ == eSelectOtherOptions)
    {
        if (!bDidOutput_ && title_ != NULL)
        {
            context_.writeTitle(title_);
        }
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

//! Formats the default option value as a string.
std::string
defaultOptionValue(const OptionInfo &option)
{
    if (option.valueCount() == 0
        || (option.valueCount() == 1 && option.formatValue(0).empty()))
    {
        return option.formatDefaultValueIfSet();
    }
    else
    {
        std::string result;
        for (int i = 0; i < option.valueCount(); ++i)
        {
            if (i != 0)
            {
                result.append(" ");
            }
            result.append(option.formatValue(i));
        }
        return result;
    }
}

//! Formats the flags for a file option as a string.
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

//! Formats the description for an option as a string.
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

/********************************************************************
 * OptionsExportFormatter
 */

/*! \brief
 * Formatter implementation for help export.
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

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsExportFormatter);
};

void OptionsExportFormatter::formatDescription(
        const HelpWriterContext &context, const Options &section)
{
    if (!section.description().empty())
    {
        // TODO: Print title for the section?
        context.writeTextBlock(section.description());
    }
}

void OptionsExportFormatter::formatFileOption(
        const HelpWriterContext &context, const FileNameOptionInfo &option)
{
    std::string defaultValue(defaultOptionValue(option));
    std::string type(fileOptionFlagsAsString(option, false));
    std::string description(option.description());
    context.writeOptionItem("-" + option.name(),
                            formatString("<%s> (%s)", defaultValue.c_str(), type.c_str()),
                            description);
}

void OptionsExportFormatter::formatOption(
        const HelpWriterContext &context, const OptionInfo &option)
{
    std::string name(option.name());
    std::string value(formatString("<%s>", option.type()));
    std::string defaultValue(defaultOptionValue(option));
    std::string description(descriptionWithOptionDetails(common_, option));
    if (option.isType<BooleanOptionInfo>())
    {
        name  = "[no]" + name;
        value = "[" + value + "]";
    }
    if (!defaultValue.empty())
    {
        value += " (" + defaultValue + ")";
    }
    context.writeOptionItem("-" + name, value, description);
}

void OptionsExportFormatter::formatSelectionOption(
        const HelpWriterContext &context, const OptionInfo &option)
{
    const char *args = (option.isType<SelectionFileOptionInfo>()
                        ? "<selection.dat>"
                        : "<selection(s)>");
    context.writeOptionItem("-" + option.name(), args, option.description());
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
        explicit Impl(const Options &options);

        //! Options object to use for generating help.
        const Options          &options_;
        //! Time unit to show in descriptions.
        std::string             timeUnit_;
        //! Whether to write descriptions to output.
        bool                    bShowDescriptions_;
};

CommandLineHelpWriter::Impl::Impl(const Options &options)
    : options_(options), timeUnit_(TimeUnitManager().timeUnitAsString()),
      bShowDescriptions_(false)
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

void CommandLineHelpWriter::writeHelp(const CommandLineHelpContext &context)
{
    const HelpWriterContext &writerContext = context.writerContext();
    CommonFormatterData      common(impl_->timeUnit_.c_str());
    OptionsExportFormatter   formatter(common);
    OptionsFilter            filter(writerContext, &formatter);
    filter.setShowHidden(context.showHidden());

    File &file = writerContext.outputFile();
    if (writerContext.outputFormat() != eHelpOutputFormat_Console)
    {
        writerContext.writeTitle("Synopsis");
        writerContext.writeTextBlock(context.moduleDisplayName());
        file.writeLine("\n\n");
    }

    if (impl_->bShowDescriptions_)
    {
        filter.formatSelected(OptionsFilter::eSelectDescriptions,
                              impl_->options_, "Description");
    }
    filter.formatSelected(OptionsFilter::eSelectFileOptions,
                          impl_->options_, "File Options");
    const bool bAnyOptionsPrinted =
        filter.formatSelected(OptionsFilter::eSelectOtherOptions,
                              impl_->options_, "Options");
    filter.formatSelected(OptionsFilter::eSelectSelectionOptions,
                          impl_->options_,
                          !bAnyOptionsPrinted ? "Options" : NULL);
}

} // namespace gmx
