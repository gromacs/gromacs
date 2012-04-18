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
 * \ingroup module_options
 */
#include "gromacs/options/cmdlinehelpwriter.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"

#include "cmdlinehelpwriter-impl.h"

namespace gmx
{

/********************************************************************
 * CommandLineHelpWriter::Impl
 */

CommandLineHelpWriter::Impl::Impl(const Options &options)
    : _options(options)
{
}


/********************************************************************
 * AsciiDescriptionWriter
 */

void AsciiDescriptionWriter::visitSubSection(const Options &section)
{
    if (!section.description().empty())
    {
        fprintf(_fp, "\n");
        const std::string &title = section.title();
        if (!title.empty())
        {
            fprintf(_fp, "%s\n\n", title.c_str());
        }
        // TODO: Wrap lines and do markup substitutions.
        fprintf(_fp, "%s\n\n", section.description().c_str());
    }
    OptionsIterator(section).acceptSubSections(this);
}


/********************************************************************
 * AsciiFileParameterWriter
 */

void AsciiFileParameterWriter::visitSubSection(const Options &section)
{
    OptionsIterator iterator(section);
    iterator.acceptSubSections(this);
    iterator.acceptOptions(this);
}

void AsciiFileParameterWriter::visitOptionType(const FileNameOptionInfo &option)
{
    if (_bFirst)
    {
        fprintf(_fp, "%6s %12s  %-12s %s\n",
                "Option", "Filename", "Type", "Description");
        fprintf(_fp, "------------------------------------------------------------\n");
        _bFirst = false;
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
    fprintf(_fp, "%s\n", optionLine.c_str());
}


/********************************************************************
 * AsciiParameterWriter
 */

void AsciiParameterWriter::visitSubSection(const Options &section)
{
    OptionsIterator iterator(section);
    iterator.acceptSubSections(this);
    iterator.acceptOptions(this);
}

void AsciiParameterWriter::visitOption(const OptionInfo &option)
{
    if (option.isType<FileNameOptionInfo>()
        || (!_bShowHidden && option.isHidden()))
    {
        return;
    }

    if (_bFirst)
    {
        fprintf(_fp, "%-12s %-6s %-6s  %s\n",
                "Option", "Type", "Value", "Description");
        fprintf(_fp, "----------------------------------------------------\n");
        _bFirst = false;
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
    fprintf(_fp, "%s\n", optionLine.c_str());
}

/********************************************************************
 * CommandLineHelpWriter
 */

CommandLineHelpWriter::CommandLineHelpWriter(const Options &options)
    : _impl(new Impl(options))
{
}

CommandLineHelpWriter::~CommandLineHelpWriter()
{
}

CommandLineHelpWriter &CommandLineHelpWriter::setShowHidden(bool bSet)
{
    _impl->_flags.set(Impl::efShowHidden, bSet);
    return *this;
}

CommandLineHelpWriter &CommandLineHelpWriter::setShowDescriptions(bool bSet)
{
    _impl->_flags.set(Impl::efShowDescriptions, bSet);
    return *this;
}

void CommandLineHelpWriter::writeHelp(FILE *fp)
{
    if (_impl->_flags.test(Impl::efShowDescriptions))
    {
        fprintf(fp, "DESCRIPTION\n"
                    "-----------\n");
        AsciiDescriptionWriter(fp).visitSubSection(_impl->_options);
    }
    {
        AsciiFileParameterWriter writer(fp);
        writer.visitSubSection(_impl->_options);
        if (writer.didOutput())
        {
            fprintf(fp, "\n");
        }
    }
    {
        AsciiParameterWriter writer(fp);
        writer.setShowHidden(_impl->_flags.test(Impl::efShowHidden));
        writer.visitSubSection(_impl->_options);
        if (writer.didOutput())
        {
            fprintf(fp, "\n");
        }
    }
}

} // namespace gmx
