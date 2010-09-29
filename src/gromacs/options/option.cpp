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
 * Implements gmx::Option.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "option.h"

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/optionflags.h"

namespace gmx
{

/********************************************************************
 * Option
 */

Option::Option()
    : _flags(0), _minValueCount(1), _maxValueCount(1), _storage(NULL),
      _currentValueCount(-1)
{
}

Option::~Option()
{
    delete _storage;
}

int Option::init(const AbstractOption &settings, Options *options)
{
    // Copy cheap values first to make them easy to access in checks.
    _flags = settings._flags;
    _minValueCount = settings._minValueCount;
    _maxValueCount = settings._maxValueCount;

    // If the maximum number of values is not known, storage to
    // caller-allocated memory is unsafe.
    if (_maxValueCount < 0 && hasFlag(efExternalStore))
    {
        GMX_ERROR(eeInvalidValue,
                  "Cannot set user-allocated storage for arbitrary number of values");
    }
    // Check that user has not provided incorrect values for vectors.
    if (hasFlag(efVector) && (_minValueCount > 1 || _maxValueCount < 1))
    {
        GMX_ERROR(eeInvalidValue,
                  "Inconsistent value counts for vector values");
    }

    int rc = settings.createDefaultStorage(options, &_storage);
    if (rc != 0)
    {
        return rc;
    }

    if (settings._name != NULL)
    {
        _name  = settings._name;
    }
    _descr = settings.createDescription();
    _flags |= efHasDefaultValue;
    return 0;
}

const char *Option::type() const
{
    return _storage->typeString();
}

int Option::valueCount() const
{
    return _storage->valueCount();
}

std::string Option::formatValue(int i) const
{
    return _storage->formatValue(i);
}

std::string Option::formatValues() const
{
    return _storage->formatValues();
}

int Option::startSource()
{
    _flags |= efHasDefaultValue;
    return 0;
}

int Option::startSet(AbstractErrorReporter *errors)
{
    if (hasFlag(efHasDefaultValue))
    {
        _flags &= ~efHasDefaultValue;
        _storage->clear();
    }
    else if (isSet() && !hasFlag(efMulti))
    {
        errors->error("Option specified multiple times");
        return eeInvalidInput;
    }
    _currentValueCount = 0;
    return 0;
}

int Option::appendValue(const std::string &value,
                        AbstractErrorReporter *errors)
{
    if (_currentValueCount == -1)
    {
        int rc = startSet(errors);
        if (rc != 0)
        {
            return rc;
        }
    }
    if (_maxValueCount >= 0 && _currentValueCount >= _maxValueCount)
    {
        errors->warning("Ignoring extra value: " + value);
        return eeInvalidInput;
    }
    int oldCount = _storage->valueCount();
    int rc = _storage->appendValue(value, errors);
    if (rc == 0)
    {
        _currentValueCount += _storage->valueCount() - oldCount;
    }
    return rc;
}

int Option::finishSet(AbstractErrorReporter *errors)
{
    int nvalues = _currentValueCount;

    _flags |= efSet;
    _currentValueCount = -1;

    // TODO: We probably should always call finishSet() to keep storage state
    // more consistent.
    if (nvalues < _minValueCount)
    {
        errors->error("Too few (valid) values");
        return eeInvalidInput;
    }
    if (_maxValueCount >= 0 && nvalues > _maxValueCount)
    {
        errors->error("Too many values");
        return eeInvalidInput;
    }

    return _storage->finishSet(nvalues, errors);
}

int Option::finish(AbstractErrorReporter *errors)
{
    if (_currentValueCount >= 0)
    {
        int rc = finishSet(errors);
        if (rc != 0)
        {
            return rc;
        }
    }
    if (hasFlag(efRequired) && !isSet())
    {
        errors->error("Required option '" + _name + "' not set");
        return eeInconsistentInput;
    }
    return _storage->finish(errors);
}

} // namespace gmx
