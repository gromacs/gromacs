/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Declares internal helper classes for the reference data framework to check
 * reference data values of different types.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_REFDATA_CHECKERS_H
#define GMX_TESTUTILS_REFDATA_CHECKERS_H

#include <cstdlib>

#include <limits>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata-impl.h"
#include "testutils/testasserts.h"
#include "testutils/testexceptions.h"

namespace gmx
{
namespace test
{

class IReferenceDataEntryChecker
{
    public:
        virtual void fillEntry(ReferenceDataEntry *entry) const = 0;
        virtual ::testing::AssertionResult
        checkEntry(const ReferenceDataEntry &entry, const std::string &fullId) const = 0;

    protected:
        virtual ~IReferenceDataEntryChecker() {}
};

class NullChecker : public IReferenceDataEntryChecker
{
    public:
        virtual void fillEntry(ReferenceDataEntry *) const {}
        virtual ::testing::AssertionResult
        checkEntry(const ReferenceDataEntry &, const std::string &) const
        {
            return ::testing::AssertionSuccess();
        }
};

class ExactStringChecker : public IReferenceDataEntryChecker
{
    public:
        explicit ExactStringChecker(const std::string &value)
            : value_(value)
        {
        }

        virtual void fillEntry(ReferenceDataEntry *entry) const
        {
            entry->setValue(value_);
        }
        virtual ::testing::AssertionResult
        checkEntry(const ReferenceDataEntry &entry, const std::string &fullId) const
        {
            if (entry.value() == value_)
            {
                return ::testing::AssertionSuccess();
            }
            return ::testing::AssertionFailure()
                   << "  In item: " << fullId << std::endl
                   << "   Actual: '" << value_ << "'" << std::endl
                   << "Reference: '" << entry.value() << "'";
        }

    private:
        std::string  value_;
};

class ExactStringBlockChecker : public IReferenceDataEntryChecker
{
    public:
        explicit ExactStringBlockChecker(const std::string &value)
            : value_(value)
        {
        }

        virtual void fillEntry(ReferenceDataEntry *entry) const
        {
            entry->setTextBlockValue(value_);
        }
        virtual ::testing::AssertionResult
        checkEntry(const ReferenceDataEntry &entry, const std::string &fullId) const
        {
            if (entry.value() == value_)
            {
                return ::testing::AssertionSuccess();
            }
            return ::testing::AssertionFailure()
                   << "  In item: " << fullId << std::endl
                   << "   Actual: '" << value_ << "'" << std::endl
                   << "Reference: '" << entry.value() << "'";
        }

    private:
        std::string  value_;
};

//! Helper function to parse a floating-point value.
// TODO: Move this into src/gromacs/utility/, and consolidate with similar code
// elsewhere.
double convertDouble(const std::string &value)
{
    char   *endptr;
    double  convertedValue = std::strtod(value.c_str(), &endptr);
    // TODO: Check for overflow
    if (*endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid floating-point value: " + value));
    }
    return convertedValue;
}

//! Helper function to parse a floating-point reference data value.
double convertDoubleReferenceValue(const std::string &value)
{
    try
    {
        return convertDouble(value);
    }
    catch (const InvalidInputError &ex)
    {
        GMX_THROW_WRAPPER_TESTEXCEPTION(ex);
    }
}

template <typename FloatType>
class FloatingPointChecker : public IReferenceDataEntryChecker
{
    public:
        FloatingPointChecker(FloatType value, const FloatingPointTolerance &tolerance)
            : value_(value), tolerance_(tolerance)
        {
        }

        virtual void fillEntry(ReferenceDataEntry *entry) const
        {
            const int prec = std::numeric_limits<FloatType>::digits10 + 2;
            entry->setValue(formatString("%.*g", prec, value_));
        }
        virtual ::testing::AssertionResult
        checkEntry(const ReferenceDataEntry &entry, const std::string &fullId) const
        {
            FloatType               refValue = static_cast<FloatType>(convertDoubleReferenceValue(entry.value()));
            FloatingPointDifference diff(refValue, value_);
            if (tolerance_.isWithin(diff))
            {
                return ::testing::AssertionSuccess();
            }
            return ::testing::AssertionFailure()
                   << "   In item: " << fullId << std::endl
                   << "    Actual: " << value_ << std::endl
                   << " Reference: " << refValue << std::endl
                   << "Difference: " << diff.toString() << std::endl
                   << " Tolerance: " << tolerance_.toString(diff);
        }

    private:
        FloatType               value_;
        FloatingPointTolerance  tolerance_;
};

template <typename FloatType>
class FloatingPointFromStringChecker : public IReferenceDataEntryChecker
{
    public:
        FloatingPointFromStringChecker(
            const std::string &value, const FloatingPointTolerance &tolerance)
            : value_(value), tolerance_(tolerance)
        {
        }

        virtual void fillEntry(ReferenceDataEntry *entry) const
        {
            entry->setValue(value_);
        }
        virtual ::testing::AssertionResult
        checkEntry(const ReferenceDataEntry &entry, const std::string &fullId) const
        {
            FloatType               value    = static_cast<FloatType>(convertDouble(value_));
            FloatType               refValue = static_cast<FloatType>(convertDoubleReferenceValue(entry.value()));
            FloatingPointDifference diff(refValue, value);
            if (tolerance_.isWithin(diff))
            {
                return ::testing::AssertionSuccess();
            }
            return ::testing::AssertionFailure()
                   << "   In item: " << fullId << std::endl
                   << "    Actual: " << value << std::endl
                   << " Reference: " << entry.value() << std::endl
                   << "Difference: " << diff.toString() << std::endl
                   << " Tolerance: " << tolerance_.toString(diff);
        }

    private:
        std::string             value_;
        FloatingPointTolerance  tolerance_;
};

} // namespace test
} // namespace gmx

#endif
