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
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gromacs/options/optionsvisitor.h"

#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/options.h"

#include "options-impl.h"

namespace gmx
{

/********************************************************************
 * OptionsIterator
 */

OptionsIterator::OptionsIterator(const Options &options)
    : options_(options)
{
}

void OptionsIterator::acceptSubSections(OptionsVisitor *visitor) const
{
    const Options::Impl::SubSectionList          &subSectionList =
        options_.impl_->subSections_;
    Options::Impl::SubSectionList::const_iterator i;
    for (i = subSectionList.begin(); i != subSectionList.end(); ++i)
    {
        visitor->visitSubSection(*(*i));
    }
}

void OptionsIterator::acceptOptions(OptionsVisitor *visitor) const
{
    const Options::Impl::OptionList          &optionList =
        options_.impl_->options_;
    Options::Impl::OptionList::const_iterator i;
    for (i = optionList.begin(); i != optionList.end(); ++i)
    {
        // This is not strictly const-correct, since optionInfo() is
        // not const (while the input options is), but this makes things much
        // simpler.
        visitor->visitOption((*i)->optionInfo());
    }
}

/********************************************************************
 * OptionsModifyingIterator
 */

OptionsModifyingIterator::OptionsModifyingIterator(Options *options)
    : options_(*options)
{
}

void OptionsModifyingIterator::acceptSubSections(OptionsModifyingVisitor *visitor) const
{
    const Options::Impl::SubSectionList          &subSectionList =
        options_.impl_->subSections_;
    Options::Impl::SubSectionList::const_iterator i;
    for (i = subSectionList.begin(); i != subSectionList.end(); ++i)
    {
        visitor->visitSubSection(*i);
    }
}

void OptionsModifyingIterator::acceptOptions(OptionsModifyingVisitor *visitor) const
{
    const Options::Impl::OptionList          &optionList =
        options_.impl_->options_;
    Options::Impl::OptionList::const_iterator i;
    for (i = optionList.begin(); i != optionList.end(); ++i)
    {
        visitor->visitOption(&(*i)->optionInfo());
    }
}

} // namespace gmx
