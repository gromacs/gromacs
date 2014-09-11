/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements classes in basicoptions.h and basicoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "basicoptions.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>

#include <limits>
#include <string>
#include <vector>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "basicoptionstorage.h"

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
    if (length > 0 && !values->empty() && values->size() != length)
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

const BooleanOptionStorage &BooleanOptionInfo::option() const
{
    return static_cast<const BooleanOptionStorage &>(OptionInfo::option());
}

bool BooleanOptionInfo::defaultValue() const
{
    return option().defaultValue();
}

/********************************************************************
 * BooleanOption
 */

AbstractOptionStorage *
BooleanOption::createStorage(const OptionManagerContainer & /*managers*/) const
{
    return new BooleanOptionStorage(*this);
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

AbstractOptionStorage *
IntegerOption::createStorage(const OptionManagerContainer & /*managers*/) const
{
    return new IntegerOptionStorage(*this);
}


/********************************************************************
 * Int64OptionStorage
 */

std::string Int64OptionStorage::formatSingleValue(const gmx_int64_t &value) const
{
    return formatString("%" GMX_PRId64, value);
}

void Int64OptionStorage::convertValue(const std::string &value)
{
    const char       *ptr = value.c_str();
    char             *endptr;
    errno = 0;
    const gmx_int64_t ival = str_to_int64_t(ptr, &endptr);
    if (errno == ERANGE)
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

/********************************************************************
 * Int64OptionInfo
 */

Int64OptionInfo::Int64OptionInfo(Int64OptionStorage *option)
    : OptionInfo(option)
{
}

/********************************************************************
 * Int64Option
 */

AbstractOptionStorage *
Int64Option::createStorage(const OptionManagerContainer & /*managers*/) const
{
    return new Int64OptionStorage(*this);
}


/********************************************************************
 * DoubleOptionStorage
 */

DoubleOptionStorage::DoubleOptionStorage(const DoubleOption &settings)
    : MyBase(settings), info_(this), bTime_(settings.bTime_), factor_(1.0)
{
}

std::string DoubleOptionStorage::typeString() const
{
    return isVector() ? "vector" : (isTime() ? "time" : "real");
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

AbstractOptionStorage *
DoubleOption::createStorage(const OptionManagerContainer & /*managers*/) const
{
    return new DoubleOptionStorage(*this);
}


/********************************************************************
 * FloatOptionStorage
 */

FloatOptionStorage::FloatOptionStorage(const FloatOption &settings)
    : MyBase(settings), info_(this), bTime_(settings.bTime_), factor_(1.0)
{
}

std::string FloatOptionStorage::typeString() const
{
    return isVector() ? "vector" : (isTime() ? "time" : "real");
}

std::string FloatOptionStorage::formatSingleValue(const float &value) const
{
    return formatString("%g", value / factor_);
}

void FloatOptionStorage::convertValue(const std::string &value)
{
    const char *ptr = value.c_str();
    char       *endptr;
    errno = 0;
    double      dval = std::strtod(ptr, &endptr);
    if (errno == ERANGE
        || dval * factor_ < -std::numeric_limits<float>::max()
        || dval * factor_ >  std::numeric_limits<float>::max())
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

void FloatOptionStorage::processSetValues(ValueList *values)
{
    if (isVector())
    {
        expandVector(maxValueCount(), values);
    }
}

void FloatOptionStorage::setScaleFactor(double factor)
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
 * FloatOptionInfo
 */

FloatOptionInfo::FloatOptionInfo(FloatOptionStorage *option)
    : OptionInfo(option)
{
}

FloatOptionStorage &FloatOptionInfo::option()
{
    return static_cast<FloatOptionStorage &>(OptionInfo::option());
}

const FloatOptionStorage &FloatOptionInfo::option() const
{
    return static_cast<const FloatOptionStorage &>(OptionInfo::option());
}

bool FloatOptionInfo::isTime() const
{
    return option().isTime();
}

void FloatOptionInfo::setScaleFactor(double factor)
{
    option().setScaleFactor(factor);
}

/********************************************************************
 * FloatOption
 */

AbstractOptionStorage *
FloatOption::createStorage(const OptionManagerContainer & /*managers*/) const
{
    return new FloatOptionStorage(*this);
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
        int                count        = settings.enumValuesCount_;
        if (count < 0)
        {
            count = 0;
            while (settings.enumValues_[count] != NULL)
            {
                ++count;
            }
        }
        for (int i = 0; i < count; ++i)
        {
            if (settings.enumValues_[i] == NULL)
            {
                GMX_THROW(APIError("Enumeration value cannot be NULL"));
            }
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
    }
    if (settings.defaultEnumIndex_ >= 0)
    {
        clear();
        addValue(allowed_[settings.defaultEnumIndex_]);
        commitValues();
    }
    // Somewhat subtly, this does not update the stored enum index if the
    // caller has not provided store() or storeVector(), because values()
    // will be empty in such a case.  This leads to (desired) behavior of
    // preserving the existing value in the enum index store variable in such
    // cases.
    refreshEnumIndexStore();
}

std::string StringOptionStorage::formatExtraDescription() const
{
    std::string result;
    if (!allowed_.empty())
    {
        result.append(": ");
        ValueList::const_iterator i;
        for (i = allowed_.begin(); i != allowed_.end(); ++i)
        {
            if (i != allowed_.begin())
            {
                result.append(", ");
            }
            result.append(*i);
        }
    }
    return result;
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
    refreshEnumIndexStore();
}

void StringOptionStorage::refreshEnumIndexStore()
{
    if (enumIndexStore_ != NULL)
    {
        for (size_t i = 0; i < values().size(); ++i)
        {
            if (values()[i].empty())
            {
                enumIndexStore_[i] = -1;
            }
            else
            {
                ValueList::const_iterator match =
                    std::find(allowed_.begin(), allowed_.end(), values()[i]);
                GMX_ASSERT(match != allowed_.end(),
                           "Enum value not found (internal error)");
                enumIndexStore_[i] = static_cast<int>(match - allowed_.begin());
            }
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

StringOptionStorage &StringOptionInfo::option()
{
    return static_cast<StringOptionStorage &>(OptionInfo::option());
}

const StringOptionStorage &StringOptionInfo::option() const
{
    return static_cast<const StringOptionStorage &>(OptionInfo::option());
}

bool StringOptionInfo::isEnumerated() const
{
    return !allowedValues().empty();
}

const std::vector<std::string> &StringOptionInfo::allowedValues() const
{
    return option().allowedValues();
}

/********************************************************************
 * StringOption
 */

AbstractOptionStorage *
StringOption::createStorage(const OptionManagerContainer & /*managers*/) const
{
    return new StringOptionStorage(*this);
}

} // namespace gmx
