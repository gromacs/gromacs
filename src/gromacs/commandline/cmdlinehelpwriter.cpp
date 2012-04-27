/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::CommandLineHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_commandline
 */
#include "cmdlinehelpwriter.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include "gromacs/options/basicoptioninfo.h"
#include "gromacs/options/filenameoptioninfo.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"

#include "cmdlinehelpwriter-impl.h"

namespace gmx
{

namespace
{

/********************************************************************
 * DescriptionWriter
 */

/*! \internal \brief
 * Helper object for writing section descriptions to help.
 *
 * \ingroup module_commandline
 */
class DescriptionWriter : public OptionsVisitor
{
    public:
        //! Creates a helper object for writing section descriptions.
        explicit DescriptionWriter(FILE *fp) : fp_(fp) {}

        virtual void visitSubSection(const Options &section);
        virtual void visitOption(const OptionInfo & /*option*/) { }

    private:
        FILE                   *fp_;
};

void DescriptionWriter::visitSubSection(const Options &section)
{
    if (!section.description().empty())
    {
        fprintf(fp_, "\n");
        const std::string &title = section.title();
        if (!title.empty())
        {
            fprintf(fp_, "%s\n\n", title.c_str());
        }
        // TODO: Wrap lines and do markup substitutions.
        fprintf(fp_, "%s\n\n", section.description().c_str());
    }
    OptionsIterator(section).acceptSubSections(this);
}


/********************************************************************
 * FileParameterWriter
 */

/*! \internal \brief
 * Helper object for writing help for file parameters.
 *
 * \ingroup module_commandline
 */
class FileParameterWriter : public OptionsTypeVisitor<FileNameOptionInfo>
{
    public:
        //! Creates a helper object for writing file parameters.
        explicit FileParameterWriter(FILE *fp)
            : fp_(fp), bFirst_(true)
        {
        }

        //! Returns true if anything was written out.
        bool didOutput() const { return !bFirst_; }

        virtual void visitSubSection(const Options &section);
        virtual void visitOptionType(const FileNameOptionInfo &option);

    private:
        FILE                   *fp_;
        bool                    bFirst_;
};

void FileParameterWriter::visitSubSection(const Options &section)
{
    OptionsIterator iterator(section);
    iterator.acceptSubSections(this);
    iterator.acceptOptions(this);
}

void FileParameterWriter::visitOptionType(const FileNameOptionInfo &option)
{
    if (bFirst_)
    {
        fprintf(fp_, "%6s %12s  %-12s %s\n",
                "Option", "Filename", "Type", "Description");
        fprintf(fp_, "------------------------------------------------------------\n");
        bFirst_ = false;
    }

    std::string optionLine("-");
    optionLine.reserve(30 + option.description().size());
    optionLine.append(option.name()).append(" ");
    if (optionLine.size() < 11)
    {
        optionLine.resize(11, ' ');
    }
    bool bTypePrinted = false;
    size_t lineStart = 0;
    for (int i = 0; i < option.valueCount(); ++i)
    {
        if (i > 0)
        {
            optionLine.append("\n");
            lineStart = optionLine.size();
            optionLine.append(11, ' ');
        }
        optionLine.append(option.formatValue(i)).append(" ");
        // TODO: Do eliding
        if (optionLine.size() <= lineStart + 21)
        {
            optionLine.resize(lineStart + 21, ' ');
            if (!bTypePrinted)
            {
                optionLine.append(option.type()).append(" ");
                bTypePrinted = true;
            }
        }
    }
    if (!bTypePrinted)
    {
        optionLine.append("\n");
        lineStart = optionLine.size();
        optionLine.append(21, ' ');
        optionLine.append(option.type()).append(" ");
    }
    if (optionLine.size() > lineStart + 34)
    {
        if (!option.description().empty())
        {
            optionLine.append("\n");
            optionLine.append(34, ' ');
        }
    }
    else
    {
        optionLine.resize(lineStart + 34, ' ');
    }
    // TODO: Markup substitution.
    optionLine.append(option.description());
    // TODO: Wrap lines.
    fprintf(fp_, "%s\n", optionLine.c_str());
}


/********************************************************************
 * ParameterWriter
 */

/*! \internal \brief
 * Helper object for writing help for non-file parameters.
 *
 * \ingroup module_commandline
 */
class ParameterWriter : public OptionsVisitor
{
    public:
        //! Creates a helper object for writing non-file parameters.
        explicit ParameterWriter(FILE *fp)
            : fp_(fp), bFirst_(true), bShowHidden_(false)
        {
        }

        //! Sets the writer to show hidden options.
        void setShowHidden(bool bSet) { bShowHidden_ = bSet; }
        //! Returns true if anything was written out.
        bool didOutput() const { return !bFirst_; }

        virtual void visitSubSection(const Options &section);
        virtual void visitOption(const OptionInfo &option);

    private:
        FILE                   *fp_;
        bool                    bFirst_;
        bool                    bShowHidden_;
};

void ParameterWriter::visitSubSection(const Options &section)
{
    OptionsIterator iterator(section);
    iterator.acceptSubSections(this);
    iterator.acceptOptions(this);
}

void ParameterWriter::visitOption(const OptionInfo &option)
{
    if (option.isType<FileNameOptionInfo>()
        || (!bShowHidden_ && option.isHidden()))
    {
        return;
    }

    if (bFirst_)
    {
        fprintf(fp_, "%-12s %-6s %-6s  %s\n",
                "Option", "Type", "Value", "Description");
        fprintf(fp_, "----------------------------------------------------\n");
        bFirst_ = false;
    }

    std::string optionLine("-");
    optionLine.reserve(30 + option.description().size());
    if (option.isType<BooleanOptionInfo>())
    {
        optionLine.append("[no]");
    }
    optionLine.append(option.name()).append(" ");
    if (optionLine.size() < 13)
    {
        optionLine.resize(13, ' ');
    }
    optionLine.append(option.type()).append(" ");
    if (optionLine.size() < 20)
    {
        optionLine.resize(20, ' ');
    }
    optionLine.append(option.formatValues()).append(" ");
    if (optionLine.size() > 28)
    {
        // TODO: Wrap lines / do eliding
        if (!option.description().empty())
        {
            optionLine.append("\n");
            optionLine.append(28, ' ');
        }
    }
    else
    {
        optionLine.resize(28, ' ');
    }
    // TODO: Markup substitution.
    optionLine.append(option.description());
    // TODO: Wrap lines.
    fprintf(fp_, "%s\n", optionLine.c_str());
}

} // namespace

/********************************************************************
 * CommandLineHelpWriter::Impl
 */

CommandLineHelpWriter::Impl::Impl(const Options &options)
    : options_(options), bShowDescriptions_(false), bShowHidden_(false)
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

void CommandLineHelpWriter::writeHelp(FILE *fp)
{
    if (impl_->bShowDescriptions_)
    {
        fprintf(fp, "DESCRIPTION\n"
                    "-----------\n");
        DescriptionWriter(fp).visitSubSection(impl_->options_);
    }
    {
        FileParameterWriter writer(fp);
        writer.visitSubSection(impl_->options_);
        if (writer.didOutput())
        {
            fprintf(fp, "\n");
        }
    }
    {
        ParameterWriter writer(fp);
        writer.setShowHidden(impl_->bShowHidden_);
        writer.visitSubSection(impl_->options_);
        if (writer.didOutput())
        {
            fprintf(fp, "\n");
        }
    }
}

} // namespace gmx
