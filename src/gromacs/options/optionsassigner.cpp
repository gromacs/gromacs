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
 * Implements gmx::OptionsAssigner.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/optionsassigner.h"

#include <deque>

#include <cassert>

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/options.h"

#include "optionsassigner-impl.h"
#include "options-impl.h"

namespace gmx
{

/********************************************************************
 * OptionsAssigner::Impl
 */

OptionsAssigner::Impl::Impl(Options *options, AbstractErrorReporter *errors)
    : _options(*options), _errors(errors), _flags(0), _currentOption(NULL),
      _errorCode(0), _currentValueCount(0), _reverseBoolean(false)
{
    _sectionStack.push_back(&_options);
}

void OptionsAssigner::Impl::setFlag(OptionsAssigner::Impl::Flag flag, bool bSet)
{
    if (bSet)
    {
        _flags |= flag;
    }
    else
    {
        _flags &= ~flag;
    }
}

AbstractOptionStorage *
OptionsAssigner::Impl::findOption(const char *name)
{
    assert(_currentOption == NULL);
    AbstractOptionStorage *option = NULL;
    Options *section = NULL;
    Options *root = &currentSection();
    Options *oldRoot = NULL;
    int      upcount = 0;
    std::deque<Options *> searchList;
    searchList.push_back(root);
    while (option == NULL && !searchList.empty())
    {
        section = searchList.front();
        option = section->_impl->findOption(name);
        if (option == NULL && hasFlag(efAcceptBooleanNoPrefix))
        {
            if (name[0] == 'n' && name[1] == 'o')
            {
                option = section->_impl->findOption(name + 2);
                if (option != NULL && option->isBoolean())
                {
                    _reverseBoolean = true;
                }
                else
                {
                    option = NULL;
                }
            }
        }
        searchList.pop_front();
        if (hasFlag(efNoStrictSectioning))
        {
            Options::Impl::SubSectionList::const_iterator i;
            for (i = section->_impl->_subSections.begin();
                 i != section->_impl->_subSections.end(); ++i)
            {
                if (*i != oldRoot)
                {
                    searchList.push_back(*i);
                }
            }
            if (searchList.empty() && root != &_options)
            {
                Options *oldRoot = root;
                root = root->_impl->_parent;
                ++upcount;
                searchList.push_back(root);
            }
        }
    }
    if (hasFlag(efNoStrictSectioning) && option != NULL)
    {
        while (upcount > 0)
        {
            _sectionStack.pop_back();
            --upcount;
        }
        std::vector<Options *> sections;
        while (section != &currentSection())
        {
            sections.push_back(section);
            section = section->_impl->_parent;
        }
        while (!sections.empty())
        {
            _sectionStack.push_back(sections.back());
            sections.pop_back();
        }
    }
    return option;
}

/********************************************************************
 * OptionsAssigner
 */

OptionsAssigner::OptionsAssigner(Options *options, AbstractErrorReporter *errors)
    : _impl(new Impl(options, errors))
{
}

OptionsAssigner::~OptionsAssigner()
{
    delete _impl;
}

AbstractErrorReporter *OptionsAssigner::errorReporter() const
{
    return _impl->_errors;
}

void OptionsAssigner::setAcceptBooleanNoPrefix(bool enabled)
{
    _impl->setFlag(Impl::efAcceptBooleanNoPrefix, enabled);
}

void OptionsAssigner::setNoStrictSectioning(bool enabled)
{
    _impl->setFlag(Impl::efNoStrictSectioning, enabled);
}

int OptionsAssigner::start()
{
    return _impl->keepError(_impl->_options._impl->startSource());
}

int OptionsAssigner::startSubSection(const char *name)
{
    if (_impl->_currentOption != NULL)
    {
        // The return code is ignored to keep on assigning, but any error is
        // stored to be returned in finish().
        finishOption();
    }

    Options *section = _impl->currentSection()._impl->findSubSection(name);
    if (section == NULL)
    {
        // TODO: Print an error
        return _impl->keepError(eeInvalidInput);
    }
    _impl->_sectionStack.push_back(section);
    return 0;
}

int OptionsAssigner::startOption(const char *name)
{
    if (_impl->_currentOption != NULL)
    {
        // The return code is ignored to keep on assigning, but any error is
        // stored to be returned in finish().
        finishOption();
    }

    AbstractOptionStorage *option = _impl->findOption(name);
    if (option == NULL)
    {
        _impl->_errors->error("Unknown option");
        return _impl->keepError(eeInvalidInput);
    }
    int rc = option->startSet(_impl->_errors);
    if (rc != 0)
    {
        return _impl->keepError(rc);
    }
    _impl->_currentOption = option;
    _impl->_currentValueCount = 0;
    return 0;
}

int OptionsAssigner::appendValue(const std::string &value)
{
    AbstractOptionStorage *option = _impl->_currentOption;
    // The option should have been successfully started.
    assert(option != NULL);
    ++_impl->_currentValueCount;
    return _impl->keepError(option->appendValue(value, _impl->_errors));
}

int OptionsAssigner::finishOption()
{
    AbstractOptionStorage *option = _impl->_currentOption;
    // The option should have been successfully started.
    assert(option != NULL);
    int rc = 0;
    if (option->isBoolean())
    {
        if (_impl->_currentValueCount == 0)
        {
            // TODO: Get rid of the hard-coded strings.
            rc = option->appendValue(_impl->_reverseBoolean ? "0" : "1",
                                     _impl->_errors);
            // If the above fails, there is something wrong.
            assert(rc == 0);
        }
        else if (_impl->_reverseBoolean)
        {
            _impl->_errors->error("Cannot specify a value together with 'no' prefix");
            rc = eeInvalidInput;
        }
    }
    int rc1 = _impl->_currentOption->finishSet(_impl->_errors);
    rc = (rc != 0 ? rc : rc1);
    _impl->_currentOption = NULL;
    _impl->_reverseBoolean = false;
    return _impl->keepError(rc);
}

int OptionsAssigner::finishSubSection()
{
    // Should only be called if we are in a subsection.
    assert(_impl->inSubSection());
    if (_impl->_currentOption != NULL)
    {
        // Possible error codes are stored and returned in the end.
        finishOption();
    }
    _impl->_sectionStack.pop_back();
    return 0;
}

int OptionsAssigner::finish()
{
    if (_impl->_currentOption != NULL)
    {
        // Possible error codes are stored and returned in the end.
        finishOption();
    }
    while (_impl->inSubSection())
    {
        // Possible error codes are stored and returned in the end.
        finishSubSection();
    }
    return _impl->_errorCode;
}

} // namespace gmx
