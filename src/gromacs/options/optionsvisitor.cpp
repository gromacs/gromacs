/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <list>

#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"

#include "options_impl.h"

namespace gmx
{

namespace
{

//! Helper function to call visitOptions() and handle correct indirection.
void visitOption(OptionsVisitor* visitor, OptionInfo& optionInfo) //NOLINT(google-runtime-references)
{
    visitor->visitOption(optionInfo);
}
//! Helper function to call visitOptions() and handle correct indirection.
void visitOption(OptionsModifyingVisitor* visitor, OptionInfo& optionInfo) //NOLINT(google-runtime-references)
{
    visitor->visitOption(&optionInfo);
}

//! Helper function to recursively visit all options in a group.
template<class VisitorType>
void acceptOptionsGroup(const internal::OptionSectionImpl::Group& group, VisitorType* visitor)
{
    for (const auto& option : group.options_)
    {
        visitOption(visitor, option->optionInfo());
    }
    for (const auto& subgroup : group.subgroups_)
    {
        acceptOptionsGroup(subgroup, visitor);
    }
}

} // namespace

/********************************************************************
 * OptionsIterator
 */

OptionsIterator::OptionsIterator(const Options& options) : section_(options.rootSection().section())
{
}

OptionsIterator::OptionsIterator(const OptionSectionInfo& section) : section_(section.section()) {}

void OptionsIterator::acceptSections(OptionsVisitor* visitor) const
{
    for (const auto& section : section_.subsections_)
    {
        visitor->visitSection(section->info());
    }
}

void OptionsIterator::acceptOptions(OptionsVisitor* visitor) const
{
    acceptOptionsGroup(section_.rootGroup_, visitor);
}

/********************************************************************
 * OptionsModifyingIterator
 */

OptionsModifyingIterator::OptionsModifyingIterator(Options* options) :
    section_(options->rootSection().section())
{
}

OptionsModifyingIterator::OptionsModifyingIterator(OptionSectionInfo* section) :
    section_(section->section())
{
}

void OptionsModifyingIterator::acceptSections(OptionsModifyingVisitor* visitor) const
{
    for (auto& section : section_.subsections_)
    {
        visitor->visitSection(&section->info());
    }
}

void OptionsModifyingIterator::acceptOptions(OptionsModifyingVisitor* visitor) const
{
    acceptOptionsGroup(section_.rootGroup_, visitor);
}

} // namespace gmx
