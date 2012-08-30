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
 * Implements gmx::Options.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/options.h"

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/messagestringcollector.h"
#include "gromacs/utility/stringutil.h"

#include "options-impl.h"

namespace gmx
{

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

void Options::addSubSection(Options *section)
{
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
    AbstractOptionStoragePointer option(settings.createStorage());
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
    MessageStringCollector errors;
    Impl::OptionList::const_iterator i;
    for (i = impl_->options_.begin(); i != impl_->options_.end(); ++i)
    {
        AbstractOptionStorage &option = **i;
        try
        {
            option.finish();
        }
        catch (const UserInputError &ex)
        {
            MessageStringContext context(&errors, "In option " + option.name());
            errors.append(ex.what());
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
        catch (const UserInputError &ex)
        {
            errors.append(ex.what());
        }
    }
    if (!errors.isEmpty())
    {
        // TODO: This exception type may not always be appropriate.
        GMX_THROW(InvalidInputError(errors.toString()));
    }
}

} // namespace gmx
