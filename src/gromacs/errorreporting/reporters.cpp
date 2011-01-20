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
 * Implements error reporter classes.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_errorreporting
 */
#include <cstdio>

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/errorreporting/emptyerrorreporter.h"
#include "gromacs/errorreporting/standarderrorreporter.h"

#include "standarderrorreporter-impl.h"

static const char *const error_type_names[] = { "note", "warning", "error" };

namespace gmx
{

/********************************************************************
 * EmptyErrorReporter
 */

void EmptyErrorReporter::startContext(const char * /*name*/)
{
}

void EmptyErrorReporter::add(ErrorType type, const char * /*reason*/)
{
    incrementCount(type);
}

void EmptyErrorReporter::finishContext()
{
}

/********************************************************************
 * StandardErrorReporter::Impl
 */

StandardErrorReporter::Impl::Impl()
    : _prevContext(-1)
{
}

/********************************************************************
 * StandardErrorReporter
 */

StandardErrorReporter::StandardErrorReporter()
    : _impl(new Impl)
{
}

StandardErrorReporter::~StandardErrorReporter()
{
    delete _impl;
}

void StandardErrorReporter::startContext(const char *name)
{
    _impl->_contexts.push_back(name);
}

void StandardErrorReporter::add(ErrorType type, const char *reason)
{
    incrementCount(type);
    int indent = (_impl->_prevContext + 1) * 2;
    if (!_impl->_contexts.empty())
    {
        std::vector<std::string>::const_iterator ci;
        for (ci = _impl->_contexts.begin() + _impl->_prevContext + 1;
             ci != _impl->_contexts.end(); ++ci)
        {
            std::fprintf(stderr, "%*s%s\n", indent, "", ci->c_str());
            indent += 2;
        }
    }
    _impl->_prevContext = _impl->_contexts.size() - 1;
    std::fprintf(stderr, "%*s%s: %s\n", indent, "",
                 error_type_names[type], reason);
}

void StandardErrorReporter::finishContext()
{
    _impl->_contexts.pop_back();
    if (_impl->_prevContext >= static_cast<int>(_impl->_contexts.size()))
    {
        _impl->_prevContext = _impl->_contexts.size() - 1;
    }
}

} // namespace gmx
