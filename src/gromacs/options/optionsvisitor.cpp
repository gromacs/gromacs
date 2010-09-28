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
 * Implements classes in optionsvisitor.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/optionsvisitor.h"

#include "gromacs/options/options.h"

#include "option.h"
#include "options-impl.h"

namespace gmx
{

/********************************************************************
 * OptionInfo
 */

OptionInfo::OptionInfo(const Option &option)
    : _option(option)
{
}

bool OptionInfo::isBoolean() const
{
    return _option.isBoolean();
}

bool OptionInfo::isFile() const
{
    return _option.isFile();
}

bool OptionInfo::isHidden() const
{
    return _option.isHidden();
}

const std::string &OptionInfo::name() const
{
    return _option.name();
}

const std::string &OptionInfo::description() const
{
    return _option.description();
}

const char *OptionInfo::type() const
{
    return _option.type();
}

int OptionInfo::valueCount() const
{
    return _option.valueCount();
}

std::string OptionInfo::formatValue(int i) const
{
    return _option.formatValue(i);
}

std::string OptionInfo::formatValues() const
{
    return _option.formatValues();
}

/********************************************************************
 * OptionsIterator
 */

OptionsIterator::OptionsIterator(const Options &options)
    : _options(options)
{
}

void OptionsIterator::acceptSubSections(OptionsVisitor *visitor) const
{
    const Options::Impl::SubSectionList &subSectionList =
        _options._impl->_subSections;
    Options::Impl::SubSectionList::const_iterator i;
    for (i = subSectionList.begin(); i != subSectionList.end(); ++i)
    {
        visitor->visitSubSection(*(*i));
    }
}

void OptionsIterator::acceptOptions(OptionsVisitor *visitor) const
{
    const Options::Impl::OptionList &optionList =
        _options._impl->_options;
    Options::Impl::OptionList::const_iterator i;
    for (i = optionList.begin(); i != optionList.end(); ++i)
    {
        visitor->visitOption(OptionInfo(*(*i)));
    }
}

} // namespace gmx
