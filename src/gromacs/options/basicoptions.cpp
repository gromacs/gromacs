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
#include "basicoptions.h"
#include "basicoptionstorage.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>

#include <limits>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace
{

/*! \brief
 * Expands a single value to a vector by copying the value.
 *
 * \tparam        ValueType  Type of values to process.
 * \param[in]     length     Length of the resulting vector.
 * \param[in,out] values     Values to process.
 * \throws   std::bad_alloc    if out of memory.
 * \throws   InvalidInputError if \p values has an invalid number of values.
 *
 * \p values should have 0, 1, or \p length values.
 * If \p values has 1 value, it is expanded such that it has \p length
 * identical values.  In other valid cases, nothing is done.
 */
template <typename ValueType>
void expandVector(size_t length, std::vector<ValueType> *values)
{
    if (length > 0 && values->size() > 0 && values->size() != length)
    {
        if (values->size() != 1)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString(
                                                     "Expected 1 or %d values, got %d", length, values->size())));
        }
        const ValueType &value = (*values)[0];
        values->resize(length, value);
    }
}

} // namespace

namespace gmx
{

/********************************************************************
 * BooleanOptionStorage
 */

std::string BooleanOptionStorage::formatSingleValue(const bool &value) const
{
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

std::string IntegerOptionStorage::formatSingleValue(const int &value) const
{
    return formatString("%d", value);
}

void IntegerOptionStorage::convertValue(const std::string &value)
{
    const char *ptr = value.c_str();
    char       *endptr;
    errno = 0;
    long int    ival = std::strtol(ptr, &endptr, 10);
    if (errno == ERANGE
        || ival < std::numeric_limits<int>::min()
        || ival > std::numeric_limits<int>::max())
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + value
                                    + "'; it causes an integer overflow"));
    }
    if (*ptr == '\0' || *endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + value
                                    + "'; expected an integer"));
    }
    addValue(ival);
}

void IntegerOptionStorage::processSetValues(ValueList *values)
{
    if (isVector())
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
    : MyBase(settings), info_(this), bTime_(settings.bTime_), factor_(1.0)
{
}

const char *DoubleOptionStorage::typeString() const
{
    return isVector() ? "vector" : (isTime() ? "time" : "double");
}

std::string DoubleOptionStorage::formatSingleValue(const double &value) const
{
    return formatString("%g", value / factor_);
}

void DoubleOptionStorage::convertValue(const std::string &value)
{
    const char *ptr = value.c_str();
    char       *endptr;
    errno = 0;
    double      dval = std::strtod(ptr, &endptr);
    if (errno == ERANGE)
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + value
                                    + "'; it causes an overflow/underflow"));
    }
    if (*ptr == '\0' || *endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + value
                                    + "'; expected a number"));
    }
    addValue(dval * factor_);
}

void DoubleOptionStorage::processSetValues(ValueList *values)
{
    if (isVector())
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
    if (!hasFlag(efOption_HasDefaultValue))
    {
        double              scale = factor / factor_;
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
    : MyBase(settings), info_(this), enumIndexStore_(NULL)
{
    if (settings.defaultEnumIndex_ >= 0 && settings.enumValues_ == NULL)
    {
        GMX_THROW(APIError("Cannot set default enum index without enum values"));
    }
    if (settings.enumIndexStore_ != NULL && settings.enumValues_ == NULL)
    {
        GMX_THROW(APIError("Cannot set enum index store without enum values"));
    }
    if (settings.enumIndexStore_ != NULL && settings.maxValueCount_ < 0)
    {
        GMX_THROW(APIError("Cannot set enum index store with arbitrary number of values"));
    }
    if (settings.enumValues_ != NULL)
    {
        enumIndexStore_ = settings.enumIndexStore_;
        const std::string *defaultValue = settings.defaultValue();
        int                match        = -1;
        for (int i = 0; settings.enumValues_[i] != NULL; ++i)
        {
            if (defaultValue != NULL && settings.enumValues_[i] == *defaultValue)
            {
                match = i;
            }
            allowed_.push_back(settings.enumValues_[i]);
        }
        if (defaultValue != NULL)
        {
            if (match < 0)
            {
                GMX_THROW(APIError("Default value is not one of allowed values"));
            }
        }
        if (settings.defaultEnumIndex_ >= 0)
        {
            if (settings.defaultEnumIndex_ >= static_cast<int>(allowed_.size()))
            {
                GMX_THROW(APIError("Default enumeration index is out of range"));
            }
            if (defaultValue != NULL && *defaultValue != allowed_[settings.defaultEnumIndex_])
            {
                GMX_THROW(APIError("Conflicting default values"));
            }
        }
        // If there is no default value, match is still -1.
        if (enumIndexStore_ != NULL)
        {
            *enumIndexStore_ = match;
        }
    }
    if (settings.defaultEnumIndex_ >= 0)
    {
        clear();
        addValue(allowed_[settings.defaultEnumIndex_]);
        commitValues();
    }
}

std::string StringOptionStorage::formatSingleValue(const std::string &value) const
{
    return value;
}

void StringOptionStorage::convertValue(const std::string &value)
{
    if (allowed_.size() == 0)
    {
        addValue(value);
    }
    else
    {
        ValueList::const_iterator  i;
        ValueList::const_iterator  match = allowed_.end();
        for (i = allowed_.begin(); i != allowed_.end(); ++i)
        {
            // TODO: Case independence.
            if (i->find(value) == 0)
            {
                if (match == allowed_.end() || i->size() < match->size())
                {
                    match = i;
                }
            }
        }
        if (match == allowed_.end())
        {
            GMX_THROW(InvalidInputError("Invalid value: " + value));
        }
        addValue(*match);
    }
}

void StringOptionStorage::refreshValues()
{
    MyBase::refreshValues();
    if (enumIndexStore_ != NULL)
    {
        for (size_t i = 0; i < values().size(); ++i)
        {
            ValueList::const_iterator match =
                std::find(allowed_.begin(), allowed_.end(), values()[i]);
            GMX_ASSERT(match != allowed_.end(),
                       "Enum value not found (internal error)");
            enumIndexStore_[i] = static_cast<int>(match - allowed_.begin());
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

    if (enumValues_ != NULL)
    {
        value.append(": ");
        for (int i = 0; enumValues_[i] != NULL; ++i)
        {
            value.append(enumValues_[i]);
            if (enumValues_[i + 1] != NULL)
            {
                value.append(enumValues_[i + 2] != NULL ? ", " : ", or ");
            }
        }
    }
    return value;
}

} // namespace gmx
