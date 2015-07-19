/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2014,2015, by the GROMACS development team, led by
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
 * Implements classes in optionsvisitor.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "optionsvisitor.h"

#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/options.h"

#include "options-impl.h"

namespace gmx
{

using internal::OptionsImpl;

namespace
{

//! Helper function to call visitOptions() and handle correct indirection.
void visitOption(OptionsVisitor *visitor, OptionInfo &optionInfo)
{
    visitor->visitOption(optionInfo);
}
//! Helper function to call visitOptions() and handle correct indirection.
void visitOption(OptionsModifyingVisitor *visitor, OptionInfo &optionInfo)
{
    visitor->visitOption(&optionInfo);
}

//! Helper function to recursively visit all options in a group.
template <class VisitorType>
void acceptOptionsGroup(const OptionsImpl::Group &group, VisitorType *visitor)
{
    OptionsImpl::Group::OptionList::const_iterator option;
    for (option = group.options_.begin(); option != group.options_.end(); ++option)
    {
        visitOption(visitor, (*option)->optionInfo());
    }
    OptionsImpl::Group::SubgroupList::const_iterator subgroup;
    for (subgroup = group.subgroups_.begin(); subgroup != group.subgroups_.end(); ++subgroup)
    {
        acceptOptionsGroup(*subgroup, visitor);
    }
}

}   // namespace

/********************************************************************
 * OptionsIterator
 */

OptionsIterator::OptionsIterator(const Options &options)
    : options_(options)
{
}

void OptionsIterator::acceptSubSections(OptionsVisitor *visitor) const
{
    const OptionsImpl::SubSectionList          &subSectionList =
        options_.impl_->subSections_;
    OptionsImpl::SubSectionList::const_iterator i;
    for (i = subSectionList.begin(); i != subSectionList.end(); ++i)
    {
        visitor->visitSubSection(*(*i));
    }
}

void OptionsIterator::acceptOptions(OptionsVisitor *visitor) const
{
    acceptOptionsGroup(options_.impl_->rootGroup_, visitor);
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
    const OptionsImpl::SubSectionList          &subSectionList =
        options_.impl_->subSections_;
    OptionsImpl::SubSectionList::const_iterator i;
    for (i = subSectionList.begin(); i != subSectionList.end(); ++i)
    {
        visitor->visitSubSection(*i);
    }
}

void OptionsModifyingIterator::acceptOptions(OptionsModifyingVisitor *visitor) const
{
    acceptOptionsGroup(options_.impl_->rootGroup_, visitor);
}

} // namespace gmx
