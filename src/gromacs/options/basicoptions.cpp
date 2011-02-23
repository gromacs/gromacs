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
 * Implements classes in basicoptions.h and basicoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/basicoptions.h"

#include <cstdio>
#include <cstdlib>

#include <string>
#include <vector>

#include "gromacs/errorreporting/abstracterrorreporter.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"

#include "basicoptionstorage.h"

template <typename T> static
int expandVector(int length, int *nvalues, gmx::AbstractErrorReporter *errors,
                 std::vector<T> *values)
{
    if (length > 0 && *nvalues > 0 && *nvalues != length)
    {
        if (*nvalues != 1)
        {
            char err_buf[256];
            std::sprintf(err_buf, "Expected 1 or %d values, got %d",
                         length, *nvalues);
            errors->error(err_buf);
            values->resize(values->size() - *nvalues);
            return gmx::eeInvalidInput;
        }
        const T &value = (*values)[values->size() - 1];
        values->resize(values->size() + length - 1, value);
        *nvalues = length;
    }
    return 0;
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

int BooleanOptionStorage::convertValue(const std::string &value,
                                       AbstractErrorReporter *errors)
{
    // TODO: Case-independence
    if (value == "1" || value == "yes" || value == "true")
    {
        addValue(true);
        return 0;
    }
    else if (value == "0" || value == "no" || value == "false")
    {
        addValue(false);
        return 0;
    }
    errors->error("Invalid value: '" + value + "'; supported values are: 1, 0, yes, no, true, false");
    return eeInvalidInput;
}

/********************************************************************
 * BooleanOption
 */

int BooleanOption::createDefaultStorage(Options *options,
                                        AbstractOptionStorage **storage) const
{
    return createOptionStorage<BooleanOption, BooleanOptionStorage>(this, options, storage);
}


/********************************************************************
 * IntegerOptionStorage
 */

IntegerOptionStorage::IntegerOptionStorage()
{
}

std::string IntegerOptionStorage::formatValue(int i) const
{
    char buf[64];
    int value = values()[i];
    std::sprintf(buf, "%d", value);
    return std::string(buf);
}

int IntegerOptionStorage::convertValue(const std::string &value,
                                       AbstractErrorReporter *errors)
{
    const char *ptr = value.c_str();
    char *endptr = NULL;
    long int ival = std::strtol(ptr, &endptr, 10);
    if (*endptr == '\0')
    {
        addValue(ival);
        return 0;
    }
    errors->error("Invalid value: " + value);
    return eeInvalidInput;
}

int IntegerOptionStorage::processSet(int nvalues, AbstractErrorReporter *errors)
{
    if (hasFlag(efVector))
    {
        int rc = expandVector(maxValueCount(), &nvalues, errors, &values());
        if (rc != 0)
        {
            return rc;
        }
    }
    return MyBase::processSet(nvalues, errors);
}

/********************************************************************
 * IntegerOption
 */

int IntegerOption::createDefaultStorage(Options *options,
                                        AbstractOptionStorage **storage) const
{
    return createOptionStorage<IntegerOption, IntegerOptionStorage>(this, options, storage);
}


/********************************************************************
 * DoubleOptionStorage
 */

DoubleOptionStorage::DoubleOptionStorage()
    : _bTime(false)
{
}

int DoubleOptionStorage::init(const DoubleOption &settings, Options *options)
{
    _bTime = settings._bTime;
    if (_bTime)
    {
        options->globalProperties().request(eogpTimeScaleFactor);
    }
    return MyBase::init(settings, options);
}

const char *DoubleOptionStorage::typeString() const
{
    return hasFlag(efVector) > 0 ? "vector" : (_bTime ? "time" : "double");
}

std::string DoubleOptionStorage::formatValue(int i) const
{
    char buf[64];
    double value = values()[i];
    if (_bTime)
    {
        double factor = hostOptions().globalProperties().timeScaleFactor();
        value /= factor;
    }
    std::sprintf(buf, "%g", value);
    return std::string(buf);
}

int DoubleOptionStorage::convertValue(const std::string &value,
                                      AbstractErrorReporter *errors)
{
    const char *ptr = value.c_str();
    char *endptr = NULL;
    double dval = std::strtod(ptr, &endptr);
    if (*endptr == '\0')
    {
        addValue(dval);
        return 0;
    }
    errors->error("Invalid value: " + value);
    return eeInvalidInput;
}

int DoubleOptionStorage::processSet(int nvalues, AbstractErrorReporter *errors)
{
    if (hasFlag(efVector))
    {
        int rc = expandVector(maxValueCount(), &nvalues, errors, &values());
        if (rc != 0)
        {
            return rc;
        }
    }
    return MyBase::processSet(nvalues, errors);
}

int DoubleOptionStorage::processAll(AbstractErrorReporter *errors)
{
    if (_bTime)
    {
        double factor = hostOptions().globalProperties().timeScaleFactor();
        ValueList::iterator i;
        for (i = values().begin(); i != values().end(); ++i)
        {
            (*i) *= factor;
        }
    }
    return MyBase::processAll(errors);
}

/********************************************************************
 * DoubleOption
 */

int DoubleOption::createDefaultStorage(Options *options,
                                       AbstractOptionStorage **storage) const
{
    return createOptionStorage<DoubleOption, DoubleOptionStorage>(this, options, storage);
}


/********************************************************************
 * StringOptionStorage
 */

StringOptionStorage::StringOptionStorage()
    : _enumIndexStore(NULL)
{
}

int StringOptionStorage::init(const StringOption &settings, Options *options)
{
    if (settings._defaultEnumIndex >= 0 && settings._enumValues == NULL)
    {
        GMX_ERROR(eeInvalidValue,
                  "Cannot set default enum index without enum values");
    }
    if (settings._enumIndexStore != NULL && settings._enumValues == NULL)
    {
        GMX_ERROR(eeInvalidValue,
                  "Cannot set enum index store without enum values");
    }
    if (settings._enumIndexStore != NULL && settings._maxValueCount < 0)
    {
        GMX_ERROR(eeInvalidValue,
                  "Cannot set enum index store with arbitrary number of values");
    }
    if (settings._enumValues != NULL)
    {
        _enumIndexStore = settings._enumIndexStore;
        const std::string *defaultValue = settings.defaultValue();
        int match = -1;
        for (int i = 0; settings._enumValues[i] != NULL; ++i)
        {
            if (defaultValue && settings._enumValues[i] == *defaultValue)
            {
                match = i;
            }
            _allowed.push_back(settings._enumValues[i]);
        }
        if (defaultValue)
        {
            if (match < 0)
            {
                GMX_ERROR(eeInvalidValue,
                          "Default value is not one of allowed values");
            }
        }
        if (settings._defaultEnumIndex >= 0)
        {
            if (settings._defaultEnumIndex >= static_cast<int>(_allowed.size()))
            {
                GMX_ERROR(eeInvalidValue,
                          "Default enumeration index is out of range");
            }
            if (defaultValue && *defaultValue != _allowed[settings._defaultEnumIndex])
            {
                GMX_ERROR(eeInvalidValue,
                          "Conflicting default values");
            }
        }
        // If there is no default value, match is still -1.
        if (_enumIndexStore != NULL)
        {
            *_enumIndexStore = match;
        }
    }
    int rc = MyBase::init(settings, options);
    if (rc == 0)
    {
        if (settings._defaultEnumIndex >= 0)
        {
            clear();
            addValue(_allowed[settings._defaultEnumIndex]);
            if (_enumIndexStore != NULL)
            {
                *_enumIndexStore = settings._defaultEnumIndex;
            }
            processValues(1, false);
        }
    }
    return rc;
}

std::string StringOptionStorage::formatValue(int i) const
{
    return values()[i];
}

int StringOptionStorage::convertValue(const std::string &value,
                                      AbstractErrorReporter *errors)
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
            errors->error("Invalid value: " + value);
            return eeInvalidInput;
        }
        addValue(*match);
        if (_enumIndexStore)
        {
            _enumIndexStore[valueCount() - 1] = (match - _allowed.begin());
        }
    }
    return 0;
}

/********************************************************************
 * StringOption
 */

int StringOption::createDefaultStorage(Options *options,
                                       AbstractOptionStorage **storage) const
{
    return createOptionStorage<StringOption, StringOptionStorage>(this, options, storage);
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


/********************************************************************
 * FileNameOptionStorage
 */

FileNameOptionStorage::FileNameOptionStorage()
    : _filetype(eftUnknown)
{
}

int FileNameOptionStorage::init(const FileNameOption &settings, Options *options)
{
    _filetype = settings._filetype;
    if (_filetype == eftPlot)
    {
        options->globalProperties().request(eogpPlotFormat);
    }
    return MyBase::init(settings, options);
}

std::string FileNameOptionStorage::formatValue(int i) const
{
    return values()[i];
}

int FileNameOptionStorage::convertValue(const std::string &value,
                                        AbstractErrorReporter * /*errors*/)
{
    // TODO: Proper implementation.
    addValue(value);
    return 0;
}

/********************************************************************
 * FileNameOption
 */

int FileNameOption::createDefaultStorage(Options *options,
                                         AbstractOptionStorage **storage) const
{
    return createOptionStorage<FileNameOption, FileNameOptionStorage>(this, options, storage);
}

} // namespace gmx
