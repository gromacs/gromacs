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
 * Implements classes in basicoptions.h, basicoptioninfo.h and
 * basicoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/basicoptions.h"

#include <cstdio>
#include <cstdlib>

#include <string>
#include <vector>

#include "gromacs/options/basicoptioninfo.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/format.h"

#include "basicoptionstorage.h"

template <typename T> static
void expandVector(size_t length, std::vector<T> *values)
{
    if (length > 0 && values->size() > 0 && values->size() != length)
    {
        if (values->size() != 1)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString(
                      "Expected 1 or %d values, got %d", length, values->size())));
        }
        const T &value = (*values)[0];
        values->resize(length, value);
    }
}

namespace gmx
{

/********************************************************************
 * BooleanOptionStorage
 */

std::string BooleanOptionStorage::formatValue(int i) const
{
    bool value = values()[i];
    return value ? "yes" : "no";
}

void BooleanOptionStorage::convertValue(const std::string &value)
{
    // TODO: Case-independence
    if (value == "1" || value == "yes" || value == "true")
    {
        addValue(true);
        return;
    }
    else if (value == "0" || value == "no" || value == "false")
    {
        addValue(false);
        return;
    }
    GMX_THROW(InvalidInputError("Invalid value: '" + value + "'; supported values are: 1, 0, yes, no, true, false"));
}

/********************************************************************
 * BooleanOptionInfo
 */

BooleanOptionInfo::BooleanOptionInfo(BooleanOptionStorage *option)
    : OptionInfo(option)
{
}

/********************************************************************
 * BooleanOption
 */

AbstractOptionStoragePointer BooleanOption::createStorage() const
{
    return AbstractOptionStoragePointer(new BooleanOptionStorage(*this));
}


/********************************************************************
 * IntegerOptionStorage
 */

std::string IntegerOptionStorage::formatValue(int i) const
{
    int value = values()[i];
    return formatString("%d", value);
}

void IntegerOptionStorage::convertValue(const std::string &value)
{
    const char *ptr = value.c_str();
    char *endptr;
    long int ival = std::strtol(ptr, &endptr, 10);
    if (*endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: " + value));
    }
    addValue(ival);
}

void IntegerOptionStorage::processSetValues(ValueList *values)
{
    if (hasFlag(efVector))
    {
        expandVector(maxValueCount(), values);
    }
}

/********************************************************************
 * IntegerOptionInfo
 */

IntegerOptionInfo::IntegerOptionInfo(IntegerOptionStorage *option)
    : OptionInfo(option)
{
}

/********************************************************************
 * IntegerOption
 */

AbstractOptionStoragePointer IntegerOption::createStorage() const
{
    return AbstractOptionStoragePointer(new IntegerOptionStorage(*this));
}


/********************************************************************
 * DoubleOptionStorage
 */

DoubleOptionStorage::DoubleOptionStorage(const DoubleOption &settings)
    : MyBase(settings), info_(this), bTime_(settings._bTime), factor_(1.0)
{
}

const char *DoubleOptionStorage::typeString() const
{
    return hasFlag(efVector) ? "vector" : (isTime() ? "time" : "double");
}

std::string DoubleOptionStorage::formatValue(int i) const
{
    return formatString("%g", values()[i] / factor_);
}

void DoubleOptionStorage::convertValue(const std::string &value)
{
    const char *ptr = value.c_str();
    char *endptr;
    double dval = std::strtod(ptr, &endptr);
    if (*endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: " + value));
    }
    addValue(dval * factor_);
}

void DoubleOptionStorage::processSetValues(ValueList *values)
{
    if (hasFlag(efVector))
    {
        expandVector(maxValueCount(), values);
    }
}

void DoubleOptionStorage::processAll()
{
}

void DoubleOptionStorage::setScaleFactor(double factor)
{
    GMX_RELEASE_ASSERT(factor > 0.0, "Invalid scaling factor");
    if (!hasFlag(efHasDefaultValue))
    {
        double scale = factor / factor_;
        ValueList::iterator i;
        for (i = values().begin(); i != values().end(); ++i)
        {
            (*i) *= scale;
        }
        refreshValues();
    }
    factor_ = factor;
}

/********************************************************************
 * DoubleOptionInfo
 */

DoubleOptionInfo::DoubleOptionInfo(DoubleOptionStorage *option)
    : OptionInfo(option)
{
}

DoubleOptionStorage &DoubleOptionInfo::option()
{
    return static_cast<DoubleOptionStorage &>(OptionInfo::option());
}

const DoubleOptionStorage &DoubleOptionInfo::option() const
{
    return static_cast<const DoubleOptionStorage &>(OptionInfo::option());
}

bool DoubleOptionInfo::isTime() const
{
    return option().isTime();
}

void DoubleOptionInfo::setScaleFactor(double factor)
{
    option().setScaleFactor(factor);
}

/********************************************************************
 * DoubleOption
 */

AbstractOptionStoragePointer DoubleOption::createStorage() const
{
    return AbstractOptionStoragePointer(new DoubleOptionStorage(*this));
}


/********************************************************************
 * StringOptionStorage
 */

StringOptionStorage::StringOptionStorage(const StringOption &settings)
    : MyBase(settings), _info(this), _enumIndexStore(NULL)
{
    if (settings._defaultEnumIndex >= 0 && settings._enumValues == NULL)
    {
        GMX_THROW(APIError("Cannot set default enum index without enum values"));
    }
    if (settings._enumIndexStore != NULL && settings._enumValues == NULL)
    {
        GMX_THROW(APIError("Cannot set enum index store without enum values"));
    }
    if (settings._enumIndexStore != NULL && settings._maxValueCount < 0)
    {
        GMX_THROW(APIError("Cannot set enum index store with arbitrary number of values"));
    }
    if (settings._enumValues != NULL)
    {
        _enumIndexStore = settings._enumIndexStore;
        const std::string *defaultValue = settings.defaultValue();
        int match = -1;
        for (int i = 0; settings._enumValues[i] != NULL; ++i)
        {
            if (defaultValue != NULL && settings._enumValues[i] == *defaultValue)
            {
                match = i;
            }
            _allowed.push_back(settings._enumValues[i]);
        }
        if (defaultValue != NULL)
        {
            if (match < 0)
            {
                GMX_THROW(APIError("Default value is not one of allowed values"));
            }
        }
        if (settings._defaultEnumIndex >= 0)
        {
            if (settings._defaultEnumIndex >= static_cast<int>(_allowed.size()))
            {
                GMX_THROW(APIError("Default enumeration index is out of range"));
            }
            if (defaultValue != NULL && *defaultValue != _allowed[settings._defaultEnumIndex])
            {
                GMX_THROW(APIError("Conflicting default values"));
            }
        }
        // If there is no default value, match is still -1.
        if (_enumIndexStore != NULL)
        {
            *_enumIndexStore = match;
        }
    }
    if (settings._defaultEnumIndex >= 0)
    {
        clear();
        addValue(_allowed[settings._defaultEnumIndex]);
        commitValues();
    }
}

std::string StringOptionStorage::formatValue(int i) const
{
    return values()[i];
}

void StringOptionStorage::convertValue(const std::string &value)
{
    if (_allowed.size() == 0)
    {
        addValue(value);
    }
    else
    {
        ValueList::const_iterator  i;
        ValueList::const_iterator  match = _allowed.end();
        for (i = _allowed.begin(); i != _allowed.end(); ++i)
        {
            // TODO: Case independence.
            if (i->find(value) == 0)
            {
                if (match == _allowed.end() || i->size() < match->size())
                {
                    match = i;
                }
            }
        }
        if (match == _allowed.end())
        {
            GMX_THROW(InvalidInputError("Invalid value: " + value));
        }
        addValue(*match);
    }
}

void StringOptionStorage::refreshValues()
{
    MyBase::refreshValues();
    if (_enumIndexStore != NULL)
    {
        for (size_t i = 0; i < values().size(); ++i)
        {
            ValueList::const_iterator match =
                std::find(_allowed.begin(), _allowed.end(), values()[i]);
            GMX_ASSERT(match != _allowed.end(),
                       "Enum value not found (internal error)");
            _enumIndexStore[i] = static_cast<int>(match - _allowed.begin());
        }
    }
}

/********************************************************************
 * StringOptionInfo
 */

StringOptionInfo::StringOptionInfo(StringOptionStorage *option)
    : OptionInfo(option)
{
}

/********************************************************************
 * StringOption
 */

AbstractOptionStoragePointer StringOption::createStorage() const
{
    return AbstractOptionStoragePointer(new StringOptionStorage(*this));
}

std::string StringOption::createDescription() const
{
    std::string value(MyBase::createDescription());

    if (_enumValues != NULL)
    {
        value.append(": ");
        for (int i = 0; _enumValues[i] != NULL; ++i)
        {
            value.append(_enumValues[i]);
            if (_enumValues[i + 1] != NULL)
            {
                value.append(_enumValues[i + 2] != NULL ? ", " : ", or ");
            }
        }
    }
    return value;
}

} // namespace gmx
