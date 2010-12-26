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
 * Implements gmx::AsciiHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/asciihelpwriter.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"

#include "asciihelpwriter-impl.h"

namespace gmx
{

/********************************************************************
 * AsciiHelpWriter::Impl
 */

AsciiHelpWriter::Impl::Impl(const Options &options)
    : _options(options)
{
}

class AsciiDescriptionWriter : public OptionsVisitor
{
    public:
        AsciiDescriptionWriter(FILE *fp) : _fp(fp) { }

        void visitSubSection(const Options &section)
        {
            const std::string &title = section.title();
            if (!title.empty())
            {
                fprintf(_fp, "\n%s\n\n", title.c_str());
            }
            // TODO: Wrap lines and do markup substitutions.
            fprintf(_fp, "%s\n\n", section.description().c_str());
            OptionsIterator(section).acceptSubSections(this);
        }
        void visitOption(const OptionInfo & /*option*/) { }

    private:
        FILE                   *_fp;
};

class AsciiFileParameterWriter : public OptionsVisitor
{
    public:
        AsciiFileParameterWriter(FILE *fp) : _fp(fp) { }

        void visitSubSection(const Options &section)
        {
            OptionsIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }

        void visitOption(const OptionInfo &option)
        {
            if (!option.isFile())
            {
                return;
            }
            // TODO: Write the information out.
        }

    private:
        FILE                   *_fp;
};

class AsciiParameterWriter : public OptionsVisitor
{
    public:
        AsciiParameterWriter(FILE *fp) : _fp(fp) { }

        void visitSubSection(const Options &section)
        {
            OptionsIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }

        void visitOption(const OptionInfo &option)
        {
            if (option.isFile())
            {
                return;
            }

            std::string optionLine("-");
            optionLine.reserve(30 + option.description().size());
            if (option.isBoolean())
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
            if (optionLine.size() < 28)
            {
                optionLine.resize(28, ' ');
            }
            // TODO: Markup substitution.
            optionLine.append(option.description());
            // TODO: Wrap lines.
            fprintf(_fp, "%s\n", optionLine.c_str());
        }

    private:
        FILE                   *_fp;
};

/********************************************************************
 * AsciiHelpWriter
 */

AsciiHelpWriter::AsciiHelpWriter(const Options &options)
    : _impl(new Impl(options))
{
}

AsciiHelpWriter::~AsciiHelpWriter()
{
    delete _impl;
}

int AsciiHelpWriter::writeHelp(FILE *fp)
{
    fprintf(fp, "DESCRIPTION\n"
                "-----------\n\n");
    AsciiDescriptionWriter(fp).visitSubSection(_impl->_options);
    if (_impl->_options.hasFileOptions())
    {
        fprintf(fp, "%6s %12s  %-12s %s\n",
                "Option", "Filename", "Type", "Description");
        fprintf(fp, "------------------------------------------------------------\n");
        AsciiFileParameterWriter(fp).visitSubSection(_impl->_options);
        fprintf(fp, "\n");
    }
    if (_impl->_options.hasNonFileOptions())
    {
        fprintf(fp, "%-12s %-6s %-6s  %s\n",
                "Option", "Type", "Value", "Description");
        fprintf(fp, "----------------------------------------------------\n");
        AsciiParameterWriter(fp).visitSubSection(_impl->_options);
    }
    return 0;
}

} // namespace gmx
