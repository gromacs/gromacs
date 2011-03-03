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
 * Implements classes in abstractoption.h and abstractoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/abstractoption.h"

#include <cassert>

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/optionflags.h"

namespace gmx
{

/********************************************************************
 * AbstractOptionStorage
 */

AbstractOptionStorage::AbstractOptionStorage()
    : _minValueCount(0), _maxValueCount(0), _currentValueCount(-1),
      _options(NULL)
{
}

AbstractOptionStorage::~AbstractOptionStorage()
{
}

int AbstractOptionStorage::init(const AbstractOption &settings,
                                Options *options)
{
    // We add user-provided flags to the ones possibly set by the subclass.
    _flags |= settings._flags;
    _minValueCount = settings._minValueCount;
    _maxValueCount = settings._maxValueCount;
    _options = options;

    // If the maximum number of values is not known, storage to
    // caller-allocated memory is unsafe.
    if ((_maxValueCount < 0 || hasFlag(efMulti)) && hasFlag(efExternalStore))
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

    if (settings._name != NULL)
    {
        _name  = settings._name;
    }
    _descr = settings.createDescription();

    return 0;
}

int AbstractOptionStorage::startSource()
{
    setFlag(efHasDefaultValue);
    return 0;
}

int AbstractOptionStorage::startSet(AbstractErrorReporter *errors)
{
    if (hasFlag(efHasDefaultValue))
    {
        clearFlag(efHasDefaultValue);
        clear();
    }
    else if (isSet() && !hasFlag(efMulti))
    {
        errors->error("Option specified multiple times");
        return eeInvalidInput;
    }
    _currentValueCount = 0;
    return 0;
}

int AbstractOptionStorage::appendValue(const std::string &value,
                                       AbstractErrorReporter *errors)
{
    assert(_currentValueCount >= 0);
    if (!hasFlag(efConversionMayNotAddValues) && _maxValueCount >= 0
        && _currentValueCount >= _maxValueCount)
    {
        errors->warning("Ignoring extra value: " + value);
        return eeInvalidInput;
    }
    return convertValue(value, errors);
}

int AbstractOptionStorage::finishSet(AbstractErrorReporter *errors)
{
    assert(_currentValueCount >= 0);

    setFlag(efSet);
    // TODO: Remove invalid values if there are too few
    int rc = processSet(_currentValueCount, errors);
    if (!hasFlag(efDontCheckMinimumCount)
        && _currentValueCount < _minValueCount)
    {
        errors->error("Too few (valid) values");
        rc = eeInvalidInput;
    }
    _currentValueCount = -1;
    return rc;
}

int AbstractOptionStorage::finish(AbstractErrorReporter *errors)
{
    assert(_currentValueCount == -1);
    int rc = processAll(errors);
    if (hasFlag(efRequired) && !isSet())
    {
        errors->error("Required option '" + _name + "' not set");
        rc = eeInconsistentInput;
    }
    return rc;
}

int AbstractOptionStorage::setMinValueCount(int count,
                                            AbstractErrorReporter *errors)
{
    assert(!hasFlag(efMulti));
    assert(count >= 0);
    _minValueCount = count;
    if (isSet()
        && !hasFlag(efDontCheckMinimumCount) && valueCount() < _minValueCount)
    {
        if (_maxValueCount == -1 || valueCount() <= _maxValueCount)
        {
            errors->error("Too few values");
        }
        return eeInconsistentInput;
    }
    return 0;
}

int AbstractOptionStorage::setMaxValueCount(int count,
                                            AbstractErrorReporter *errors)
{
    assert(!hasFlag(efMulti));
    assert(count >= -1);
    _maxValueCount = count;
    if (isSet() && _maxValueCount >= 0 && valueCount() > _maxValueCount)
    {
        if (hasFlag(efDontCheckMinimumCount) || valueCount() >= _minValueCount)
        {
            errors->error("Too many values");
        }
        return eeInconsistentInput;
    }
    return 0;
}

int AbstractOptionStorage::incrementValueCount()
{
    if (_currentValueCount == -1)
    {
        return 0;
    }
    if (_maxValueCount >= 0 && _currentValueCount >= _maxValueCount)
    {
        return eeInvalidInput;
    }
    ++_currentValueCount;
    return 0;
}

} // namespace gmx
