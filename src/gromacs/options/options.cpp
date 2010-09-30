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

#include <cassert>
#include <cctype>
#include <cstring>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/globalproperties.h"

#include "options-impl.h"

namespace gmx
{

static std::string composeString(const char *const *sarray)
{
    std::string result;

    for (int i = 0; sarray[i] != NULL; ++i)
    {
        if (sarray[i][0] != '\0')
        {
            result.append(sarray[i]);
            char lastchar = sarray[i][std::strlen(sarray[i])-1];
            if (!std::isspace(lastchar))
            {
                result.append(" ");
            }
        }
    }
    result.resize(result.find_last_not_of(" \n\r\t") + 1);
    return result;
}

/********************************************************************
 * Options::Impl
 */

Options::Impl::Impl(const char *name, const char *title)
    : _name(name != NULL ? name : ""), _title(title != NULL ? title : ""),
      _parent(NULL), _globalProperties(new OptionsGlobalProperties)
{
}

Options::Impl::~Impl()
{
    OptionList::const_iterator i;
    for (i = _options.begin(); i != _options.end(); ++i)
    {
        delete *i;
    }
    delete _globalProperties;
}

Options *Options::Impl::findSubSection(const char *name) const
{
    SubSectionList::const_iterator i;
    for (i = _subSections.begin(); i != _subSections.end(); ++i)
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
    for (i = _options.begin(); i != _options.end(); ++i)
    {
        if ((*i)->name() == name)
        {
            return *i;
        }
    }
    return NULL;
}

int Options::Impl::startSource()
{
    int rc = 0;
    OptionList::const_iterator i;
    for (i = _options.begin(); i != _options.end(); ++i)
    {
        AbstractOptionStorage *option = *i;
        int rc1 = option->startSource();
        rc = (rc != 0 ? rc : rc1);
    }
    SubSectionList::const_iterator j;
    for (j = _subSections.begin(); j != _subSections.end(); ++j)
    {
        Options *section = *j;
        int rc1 = section->_impl->startSource();
        rc = (rc != 0 ? rc : rc1);
    }
    return rc;
}

/********************************************************************
 * Options
 */

Options::Options(const char *name, const char *title)
    : _impl(new Impl(name, title))
{
}

Options::~Options()
{
    delete _impl;
}

const std::string &Options::name() const
{
    return _impl->_name;
}

const std::string &Options::title() const
{
    return _impl->_title;
}

const std::string &Options::description() const
{
    return _impl->_description;
}

void Options::setDescription(const char *const *desc)
{
    _impl->_description = composeString(desc);
}

void Options::addSubSection(Options *section)
{
    // Make sure that section is not already inserted somewhere.
    assert(section->_impl->_parent == NULL);
    // Make sure that there are no duplicate sections.
    assert(_impl->findSubSection(section->name().c_str()) == NULL);
    _impl->_subSections.push_back(section);
    section->_impl->_parent = this;

    globalProperties()._usedProperties |=
        section->_impl->_globalProperties->_usedProperties;
    delete section->_impl->_globalProperties;
    section->_impl->_globalProperties = NULL;
}

void Options::addOption(const AbstractOption &settings)
{
    AbstractOptionStorage *option = NULL;
    int rc = settings.createDefaultStorage(this, &option);
    // Caller code should be fixed if option initialization fails.
    assert(rc == 0);
    // Make sure that there are no duplicate options.
    assert(_impl->findOption(option->name().c_str()) == NULL);
    _impl->_options.push_back(option);
}

void Options::addDefaultOptions()
{
    globalProperties().addDefaultOptions(this);
}

bool Options::isSet(const char *name) const
{
    AbstractOptionStorage *option = _impl->findOption(name);
    return (option != NULL ? option->isSet() : false);
}

int Options::finish(AbstractErrorReporter *errors)
{
    int rc = 0;
    Impl::OptionList::const_iterator i;
    for (i = _impl->_options.begin(); i != _impl->_options.end(); ++i)
    {
        AbstractOptionStorage *option = *i;
        int rc1 = option->finish(errors);
        rc = (rc != 0 ? rc : rc1);
    }
    Impl::SubSectionList::const_iterator j;
    for (j = _impl->_subSections.begin(); j != _impl->_subSections.end(); ++j)
    {
        Options *section = *j;
        int rc1 = section->finish(errors);
        rc = (rc != 0 ? rc : rc1);
    }
    return rc;
}

OptionsGlobalProperties &Options::globalProperties()
{
    Options *section = this;
    while (section->_impl->_parent != NULL)
    {
        section = section->_impl->_parent;
    }
    return *section->_impl->_globalProperties;
}

} // namespace gmx
