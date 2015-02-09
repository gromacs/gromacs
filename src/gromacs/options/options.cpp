/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::Options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "options.h"

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "options-impl.h"

namespace gmx
{

/********************************************************************
 * OptionManagerInterface
 */

OptionManagerInterface::~OptionManagerInterface()
{
}

/********************************************************************
 * Options::Impl
 */

Options::Impl::Impl(const char *name, const char *title)
    : name_(name != NULL ? name : ""), title_(title != NULL ? title : ""),
      parent_(NULL)
{
}

Options::Impl::~Impl()
{
}

Options *Options::Impl::findSubSection(const char *name) const
{
    SubSectionList::const_iterator i;
    for (i = subSections_.begin(); i != subSections_.end(); ++i)
    {
        if ((*i)->name() == name)
        {
            return *i;
        }
    }
    return NULL;
}

AbstractOptionStorage *Options::Impl::findOption(const char *name) const
{
    OptionList::const_iterator i;
    for (i = options_.begin(); i != options_.end(); ++i)
    {
        if ((*i)->name() == name)
        {
            return i->get();
        }
    }
    return NULL;
}

void Options::Impl::startSource()
{
    OptionList::const_iterator i;
    for (i = options_.begin(); i != options_.end(); ++i)
    {
        AbstractOptionStorage &option = **i;
        option.startSource();
    }
    SubSectionList::const_iterator j;
    for (j = subSections_.begin(); j != subSections_.end(); ++j)
    {
        Options &section = **j;
        section.impl_->startSource();
    }
}

/********************************************************************
 * Options
 */

Options::Options(const char *name, const char *title)
    : impl_(new Impl(name, title))
{
}

Options::~Options()
{
}

const std::string &Options::name() const
{
    return impl_->name_;
}

const std::string &Options::title() const
{
    return impl_->title_;
}

const std::string &Options::description() const
{
    return impl_->description_;
}

void Options::setDescription(const std::string &desc)
{
    impl_->description_ = desc;
}

void Options::setDescription(const ConstArrayRef<const char *> &descArray)
{
    impl_->description_ = joinStrings(descArray, "\n");
}

void Options::addManager(OptionManagerInterface *manager)
{
    GMX_RELEASE_ASSERT(impl_->parent_ == NULL,
                       "Can only add a manager in a top-level Options object");
    // This ensures that all options see the same set of managers.
    GMX_RELEASE_ASSERT(impl_->options_.empty(),
                       "Can only add a manager before options");
    // This check could be relaxed if we instead checked that the subsections
    // do not have options.
    GMX_RELEASE_ASSERT(impl_->subSections_.empty(),
                       "Can only add a manager before subsections");
    impl_->managers_.add(manager);
}

void Options::addSubSection(Options *section)
{
    // This is required, because managers are used from the root Options
    // object, so they are only seen after the subsection has been added.
    GMX_RELEASE_ASSERT(section->impl_->options_.empty(),
                       "Can only add a subsection before it has any options");
    GMX_RELEASE_ASSERT(section->impl_->managers_.empty(),
                       "Can only have managers in a top-level Options object");
    // Make sure that section is not already inserted somewhere.
    GMX_RELEASE_ASSERT(section->impl_->parent_ == NULL,
                       "Cannot add as subsection twice");
    // Make sure that there are no duplicate sections.
    GMX_RELEASE_ASSERT(impl_->findSubSection(section->name().c_str()) == NULL,
                       "Duplicate subsection name");
    impl_->subSections_.push_back(section);
    section->impl_->parent_ = this;
}

OptionInfo *Options::addOption(const AbstractOption &settings)
{
    Options::Impl *root = impl_.get();
    while (root->parent_ != NULL)
    {
        root = root->parent_->impl_.get();
    }
    Impl::AbstractOptionStoragePointer option(settings.createStorage(root->managers_));
    if (impl_->findOption(option->name().c_str()) != NULL)
    {
        GMX_THROW(APIError("Duplicate option: " + option->name()));
    }
    impl_->options_.push_back(move(option));
    return &impl_->options_.back()->optionInfo();
}

bool Options::isSet(const char *name) const
{
    AbstractOptionStorage *option = impl_->findOption(name);
    return (option != NULL ? option->isSet() : false);
}

void Options::finish()
{
    // TODO: Consider how to customize these error messages based on context.
    ExceptionInitializer             errors("Invalid input values");
    Impl::OptionList::const_iterator i;
    for (i = impl_->options_.begin(); i != impl_->options_.end(); ++i)
    {
        AbstractOptionStorage &option = **i;
        try
        {
            option.finish();
        }
        catch (UserInputError &ex)
        {
            ex.prependContext("In option " + option.name());
            errors.addCurrentExceptionAsNested();
        }
    }
    Impl::SubSectionList::const_iterator j;
    for (j = impl_->subSections_.begin(); j != impl_->subSections_.end(); ++j)
    {
        Options &section = **j;
        try
        {
            section.finish();
        }
        catch (const UserInputError &)
        {
            errors.addCurrentExceptionAsNested();
        }
    }
    if (errors.hasNestedExceptions())
    {
        // TODO: This exception type may not always be appropriate.
        GMX_THROW(InvalidInputError(errors));
    }
}

} // namespace gmx
