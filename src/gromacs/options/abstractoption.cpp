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

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/optionflags.h"

namespace gmx
{

/********************************************************************
 * AbstractOptionStorage
 */

AbstractOptionStorage::AbstractOptionStorage(const AbstractOption &settings,
                                             Options *options,
                                             OptionFlags staticFlags)
    : _flags(settings._flags | staticFlags),
      _minValueCount(settings._minValueCount),
      _maxValueCount(settings._maxValueCount),
      _currentValueCount(-1),
      _options(options)
{
    // If the maximum number of values is not known, storage to
    // caller-allocated memory is unsafe.
    if ((_maxValueCount < 0 || hasFlag(efMulti)) && hasFlag(efExternalStore))
    {
        GMX_THROW(APIError("Cannot set user-allocated storage for arbitrary number of values"));
    }
    // Check that user has not provided incorrect values for vectors.
    if (hasFlag(efVector) && (_minValueCount > 1 || _maxValueCount < 1))
    {
        GMX_THROW(APIError("Inconsistent value counts for vector values"));
    }

    if (settings._name != NULL)
    {
        _name  = settings._name;
    }
    _descr = settings.createDescription();
}

AbstractOptionStorage::~AbstractOptionStorage()
{
}

void AbstractOptionStorage::startSource()
{
    setFlag(efHasDefaultValue);
}

void AbstractOptionStorage::startSet()
{
    if (hasFlag(efHasDefaultValue))
    {
        clearFlag(efHasDefaultValue);
        clear();
    }
    else if (isSet() && !hasFlag(efMulti))
    {
        GMX_THROW(InvalidInputError("Option specified multiple times"));
    }
    _currentValueCount = 0;
}

void AbstractOptionStorage::appendValue(const std::string &value)
{
    GMX_RELEASE_ASSERT(_currentValueCount >= 0, "startSet() not called");
    if (!hasFlag(efConversionMayNotAddValues) && _maxValueCount >= 0
        && _currentValueCount >= _maxValueCount)
    {
        GMX_THROW(InvalidInputError("Ignoring extra value: " + value));
    }
    convertValue(value);
}

class CurrentCountClearer
{
    public:
        explicit CurrentCountClearer(int *currentValueCount)
            : value_(currentValueCount)
        {
        }
        ~CurrentCountClearer()
        {
            *value_ = -1;
        }

    private:
        int                    *value_;
};

void AbstractOptionStorage::finishSet()
{
    GMX_RELEASE_ASSERT(_currentValueCount >= 0, "startSet() not called");

    setFlag(efSet);

    CurrentCountClearer clearOnExit(&_currentValueCount);
    // TODO: Remove invalid values
    processSet(_currentValueCount);
    // TODO: Should this also be checked if processSet throws?
    if (!hasFlag(efDontCheckMinimumCount)
        && _currentValueCount < _minValueCount)
    {
        GMX_THROW(InvalidInputError("Too few (valid) values"));
    }
}

void AbstractOptionStorage::finish()
{
    GMX_RELEASE_ASSERT(_currentValueCount == -1, "finishSet() not called");
    processAll();
    // TODO: Should this also be checked if processAll throws?
    if (hasFlag(efRequired) && !isSet())
    {
        GMX_THROW(InvalidInputError("Option is required, but not set"));
    }
}

void AbstractOptionStorage::setMinValueCount(int count)
{
    GMX_RELEASE_ASSERT(!hasFlag(efMulti),
                       "setMinValueCount() not supported with efMulti");
    GMX_RELEASE_ASSERT(count >= 0, "Invalid value count");
    _minValueCount = count;
    if (isSet()
        && !hasFlag(efDontCheckMinimumCount) && valueCount() < _minValueCount)
    {
        GMX_THROW(InvalidInputError("Too few values"));
    }
}

void AbstractOptionStorage::setMaxValueCount(int count)
{
    GMX_RELEASE_ASSERT(!hasFlag(efMulti),
                       "setMaxValueCount() not supported with efMulti");
    GMX_RELEASE_ASSERT(count >= -1, "Invalid value count");
    _maxValueCount = count;
    if (isSet() && _maxValueCount >= 0 && valueCount() > _maxValueCount)
    {
        GMX_THROW(InvalidInputError("Too many values"));
    }
}

void AbstractOptionStorage::incrementValueCount()
{
    if (_currentValueCount == -1)
    {
        return;
    }
    if (_maxValueCount >= 0 && _currentValueCount >= _maxValueCount)
    {
        GMX_THROW(InvalidInputError("Too many values"));
    }
    ++_currentValueCount;
}

} // namespace gmx
