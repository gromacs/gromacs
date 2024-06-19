/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements classes in basicoptions.h and basicoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/basicoptions.h"

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/ivaluestore.h"
#include "gromacs/options/optionflags.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "basicoptionstorage.h"

namespace gmx
{
class AbstractOptionStorage;
class OptionManagerContainer;
} // namespace gmx

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
 *
 * \ingroup module_options
 */
template<typename ValueType>
void expandVector(size_t length, std::vector<ValueType>* values)
{
    if (length > 0 && !values->empty() && values->size() != length)
    {
        if (values->size() != 1)
        {
            GMX_THROW(gmx::InvalidInputError(
                    gmx::formatString("Expected 1 or %zu values, got %zu", length, values->size())));
        }
        const ValueType& value = (*values)[0];
        values->resize(length, value);
    }
}

/*! \brief
 * Finds an enumerated value from the list of allowed values.
 *
 * \param[in] allowedValues  List of allowed values.
 * \param[in] value          Value to search for.
 * \throws    gmx::InvalidInputError if \p value does not match anything in
 *     \p allowedValues.
 * \returns   Iterator to the found value.
 *
 * \ingroup module_options
 */
std::vector<std::string>::const_iterator findEnumValue(const std::vector<std::string>& allowedValues,
                                                       const std::string&              value)
{
    std::vector<std::string>::const_iterator i;
    std::vector<std::string>::const_iterator match = allowedValues.end();
    for (i = allowedValues.begin(); i != allowedValues.end(); ++i)
    {
        // TODO: Case independence.
        if (gmx::startsWith(*i, value))
        {
            if (match == allowedValues.end() || i->size() < match->size())
            {
                match = i;
            }
        }
    }
    if (match == allowedValues.end())
    {
        GMX_THROW(gmx::InvalidInputError("Invalid value: " + value));
    }
    return match;
}

} // namespace

namespace gmx
{

/********************************************************************
 * BooleanOptionStorage
 */

std::string BooleanOptionStorage::formatSingleValue(const bool& value) const
{
    return value ? "yes" : "no";
}

void BooleanOptionStorage::initConverter(ConverterType* converter)
{
    converter->addConverter<std::string>(&fromStdString<bool>);
}

/********************************************************************
 * BooleanOptionInfo
 */

BooleanOptionInfo::BooleanOptionInfo(BooleanOptionStorage* option) : OptionInfo(option) {}

const BooleanOptionStorage& BooleanOptionInfo::option() const
{
    return static_cast<const BooleanOptionStorage&>(OptionInfo::option());
}

bool BooleanOptionInfo::defaultValue() const
{
    return option().defaultValue();
}

/********************************************************************
 * BooleanOption
 */

AbstractOptionStorage* BooleanOption::createStorage(const OptionManagerContainer& /*managers*/) const
{
    return new BooleanOptionStorage(*this);
}


/********************************************************************
 * IntegerOptionStorage
 */

std::string IntegerOptionStorage::formatSingleValue(const int& value) const
{
    return toString(value);
}

void IntegerOptionStorage::initConverter(ConverterType* converter)
{
    converter->addConverter<std::string>(&fromStdString<int>);
}

void IntegerOptionStorage::processSetValues(ValueList* values)
{
    if (isVector())
    {
        expandVector(maxValueCount(), values);
    }
}

/********************************************************************
 * IntegerOptionInfo
 */

IntegerOptionInfo::IntegerOptionInfo(IntegerOptionStorage* option) : OptionInfo(option) {}

/********************************************************************
 * IntegerOption
 */

AbstractOptionStorage* IntegerOption::createStorage(const OptionManagerContainer& /*managers*/) const
{
    return new IntegerOptionStorage(*this);
}


/********************************************************************
 * Int64OptionStorage
 */

std::string Int64OptionStorage::formatSingleValue(const int64_t& value) const
{
    return toString(value);
}

void Int64OptionStorage::initConverter(ConverterType* converter)
{
    converter->addConverter<std::string>(&fromStdString<int64_t>);
}

/********************************************************************
 * Int64OptionInfo
 */

Int64OptionInfo::Int64OptionInfo(Int64OptionStorage* option) : OptionInfo(option) {}

/********************************************************************
 * Int64Option
 */

AbstractOptionStorage* Int64Option::createStorage(const OptionManagerContainer& /*managers*/) const
{
    return new Int64OptionStorage(*this);
}


/********************************************************************
 * DoubleOptionStorage
 */

DoubleOptionStorage::DoubleOptionStorage(const DoubleOption& settings) :
    MyBase(settings), info_(this), bTime_(settings.bTime_), factor_(1.0)
{
}

std::string DoubleOptionStorage::typeString() const
{
    return isVector() ? "vector" : (isTime() ? "time" : "real");
}

std::string DoubleOptionStorage::formatSingleValue(const double& value) const
{
    return toString(value / factor_);
}

void DoubleOptionStorage::initConverter(ConverterType* converter)
{
    converter->addConverter<std::string>(&fromStdString<double>);
    converter->addCastConversion<float>();
}

double DoubleOptionStorage::processValue(const double& value) const
{
    // TODO: Consider testing for overflow when scaling with factor_.
    return value * factor_;
}

void DoubleOptionStorage::processSetValues(ValueList* values)
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
        double scale = factor / factor_;
        for (double& value : values())
        {
            value *= scale;
        }
    }
    factor_ = factor;
}

/********************************************************************
 * DoubleOptionInfo
 */

DoubleOptionInfo::DoubleOptionInfo(DoubleOptionStorage* option) : OptionInfo(option) {}

DoubleOptionStorage& DoubleOptionInfo::option()
{
    return static_cast<DoubleOptionStorage&>(OptionInfo::option());
}

const DoubleOptionStorage& DoubleOptionInfo::option() const
{
    return static_cast<const DoubleOptionStorage&>(OptionInfo::option());
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

AbstractOptionStorage* DoubleOption::createStorage(const OptionManagerContainer& /*managers*/) const
{
    return new DoubleOptionStorage(*this);
}


/********************************************************************
 * FloatOptionStorage
 */

FloatOptionStorage::FloatOptionStorage(const FloatOption& settings) :
    MyBase(settings), info_(this), bTime_(settings.bTime_), factor_(1.0)
{
}

std::string FloatOptionStorage::typeString() const
{
    return isVector() ? "vector" : (isTime() ? "time" : "real");
}

std::string FloatOptionStorage::formatSingleValue(const float& value) const
{
    return toString(value / factor_);
}

void FloatOptionStorage::initConverter(ConverterType* converter)
{
    converter->addConverter<std::string>(&fromStdString<float>);
    converter->addCastConversion<double>();
}

float FloatOptionStorage::processValue(const float& value) const
{
    // TODO: Consider testing for overflow when scaling with factor_.
    return value * factor_;
}

void FloatOptionStorage::processSetValues(ValueList* values)
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
        float scale = factor / factor_;
        for (float& value : values())
        {
            value *= scale;
        }
    }
    factor_ = factor;
}

/********************************************************************
 * FloatOptionInfo
 */

FloatOptionInfo::FloatOptionInfo(FloatOptionStorage* option) : OptionInfo(option) {}

FloatOptionStorage& FloatOptionInfo::option()
{
    return static_cast<FloatOptionStorage&>(OptionInfo::option());
}

const FloatOptionStorage& FloatOptionInfo::option() const
{
    return static_cast<const FloatOptionStorage&>(OptionInfo::option());
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

AbstractOptionStorage* FloatOption::createStorage(const OptionManagerContainer& /*managers*/) const
{
    return new FloatOptionStorage(*this);
}


/********************************************************************
 * StringOptionStorage
 */

StringOptionStorage::StringOptionStorage(const StringOption& settings) :
    MyBase(settings), info_(this)
{
    if (settings.defaultEnumIndex_ >= 0 && settings.enumValues_ == nullptr)
    {
        GMX_THROW(APIError("Cannot set default enum index without enum values"));
    }
    if (settings.enumValues_ != nullptr)
    {
        int count = settings.enumValuesCount_;
        if (count < 0)
        {
            count = 0;
            while (settings.enumValues_[count] != nullptr)
            {
                ++count;
            }
        }
        for (int i = 0; i < count; ++i)
        {
            if (settings.enumValues_[i] == nullptr)
            {
                GMX_THROW(APIError("Enumeration value cannot be NULL"));
            }
            allowed_.emplace_back(settings.enumValues_[i]);
        }
        if (settings.defaultEnumIndex_ >= 0)
        {
            if (settings.defaultEnumIndex_ >= count)
            {
                GMX_THROW(APIError("Default enumeration index is out of range"));
            }
            const std::string* defaultValue = settings.defaultValue();
            if (defaultValue != nullptr && *defaultValue != allowed_[settings.defaultEnumIndex_])
            {
                GMX_THROW(APIError("Conflicting default values"));
            }
            setDefaultValue(allowed_[settings.defaultEnumIndex_]);
        }
    }
}

std::string StringOptionStorage::formatExtraDescription() const
{
    std::string result;
    if (!allowed_.empty())
    {
        result.append(": ");
        result.append(joinStrings(allowed_, ", "));
    }
    return result;
}

std::string StringOptionStorage::formatSingleValue(const std::string& value) const
{
    return value;
}

void StringOptionStorage::initConverter(ConverterType* /*converter*/) {}

std::string StringOptionStorage::processValue(const std::string& value) const
{
    if (!allowed_.empty())
    {
        return *findEnumValue(this->allowed_, value);
    }
    return value;
}

/********************************************************************
 * StringOptionInfo
 */

StringOptionInfo::StringOptionInfo(StringOptionStorage* option) : OptionInfo(option) {}

const StringOptionStorage& StringOptionInfo::option() const
{
    return static_cast<const StringOptionStorage&>(OptionInfo::option());
}

bool StringOptionInfo::isEnumerated() const
{
    return !allowedValues().empty();
}

const std::vector<std::string>& StringOptionInfo::allowedValues() const
{
    return option().allowedValues();
}

/********************************************************************
 * StringOption
 */

AbstractOptionStorage* StringOption::createStorage(const OptionManagerContainer& /*managers*/) const
{
    return new StringOptionStorage(*this);
}


/********************************************************************
 * EnumOptionStorage
 */

EnumOptionStorage::EnumOptionStorage(const AbstractOption& settings,
                                     const char* const*    enumValues,
                                     int                   count,
                                     int                   defaultValue,
                                     int                   defaultValueIfSet,
                                     StorePointer          store) :
    MyBase(settings, std::move(store)), info_(this)
{
    if (enumValues == nullptr)
    {
        GMX_THROW(APIError("Allowed values must be provided to EnumOption"));
    }

    if (count < 0)
    {
        count = 0;
        while (enumValues[count] != nullptr)
        {
            ++count;
        }
    }
    for (int i = 0; i < count; ++i)
    {
        if (enumValues[i] == nullptr)
        {
            GMX_THROW(APIError("Enumeration value cannot be NULL"));
        }
        allowed_.emplace_back(enumValues[i]);
    }

    GMX_ASSERT(defaultValue < count, "Default enumeration value is out of range");
    GMX_ASSERT(defaultValueIfSet < count, "Default enumeration value is out of range");
    setFlag(efOption_HasDefaultValue);
    if (defaultValue >= 0)
    {
        setDefaultValue(defaultValue);
    }
    if (defaultValueIfSet >= 0)
    {
        setDefaultValueIfSet(defaultValueIfSet);
    }
}

std::string EnumOptionStorage::formatExtraDescription() const
{
    std::string result;
    result.append(": ");
    result.append(joinStrings(allowed_, ", "));
    return result;
}

std::string EnumOptionStorage::formatSingleValue(const int& value) const
{
    if (value < 0 || value >= gmx::ssize(allowed_))
    {
        return std::string();
    }
    return allowed_[value];
}

Any EnumOptionStorage::normalizeValue(const int& value) const
{
    return Any::create<std::string>(formatSingleValue(value));
}

void EnumOptionStorage::initConverter(ConverterType* converter)
{
    converter->addConverter<std::string>([this](const std::string& value) {
        return static_cast<int>(findEnumValue(this->allowed_, value) - this->allowed_.begin());
    });
}

/********************************************************************
 * EnumOptionInfo
 */

EnumOptionInfo::EnumOptionInfo(EnumOptionStorage* option) : OptionInfo(option) {}

const EnumOptionStorage& EnumOptionInfo::option() const
{
    return static_cast<const EnumOptionStorage&>(OptionInfo::option());
}

const std::vector<std::string>& EnumOptionInfo::allowedValues() const
{
    return option().allowedValues();
}

/********************************************************************
 * EnumOption helpers
 */

namespace internal
{

//! \cond internal
AbstractOptionStorage* createEnumOptionStorage(const AbstractOption& option,
                                               const char* const*    enumValues,
                                               int                   count,
                                               int                   defaultValue,
                                               int                   defaultValueIfSet,
                                               std::unique_ptr<IOptionValueStore<int>> store)
{
    return new EnumOptionStorage(
            option, enumValues, count, defaultValue, defaultValueIfSet, std::move(store));
}
//! \endcond

} // namespace internal

} // namespace gmx
