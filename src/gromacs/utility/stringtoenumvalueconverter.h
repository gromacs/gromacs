/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
/*! \file
 * \brief Defines helper function object for class enumerations
 *
 * This helper type facilitates efficient lookup of values from
 * an enumeration identified by a string, given a pre-declared mapping of
 * values to such strings.
 *
 * It is separated from enumerationhelpers.h because it is not as
 * widely used and brings in several extra dependencies not needed by
 * that header.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_STRINGTOENUMVALUECONVERTER_H
#define GMX_UTILITY_STRINGTOENUMVALUECONVERTER_H

#include <map>
#include <optional>
#include <string>

#include "enumerationhelpers.h"
#include "stringcompare.h"
#include "stringutil.h"

namespace gmx
{

/*! \brief Enum class for whether StringToEnumValueConverter will strip strings
 * of leading and trailing whitespace before comparison. */
enum class StripStrings : int
{
    //! Do not strip strings
    No,
    //! Strip strings
    Yes
};

/*! \brief A class to convert a string to an enum value of type \c EnumType.
 *
 * It can be configured:
 *  - to match different enumerations,
 *  - to convert enumerations to strings in a custom way,
 *  - either strip strings of leading and trailing whitespace before
 *    comparison or not,
 *  - with different handling of how the string characters are compared.
 *
 * Usage example:
 *
 *    enum class Foo : int {
 *       Fizz, Buzz, Count, Default = Fizz
 *    };
 *    StringToEnumValueConverter<Foo, enumValueToString> converter;
 *    Foo type = converter.valueFrom(theString);
 *
 * \tparam EnumType                    A class enum for which enumValueToString
 *                                     is defined and maps all values
 *                                     (except EnumType::Count) to a string.
 * \tparam enumValueToStringFunction   Function to convert EnumValue to string, which
 *                                     is typically enumValueToString, per convention
 * \tparam stringCompareType           Indicates how the string should be compared
 *                                     with respect to case, hyphens, underscores, etc.
 * \tparam stripStrings                Indicates whether strings should have leading
 *                                     and trailing whitespace removed before comparison
 */
template<typename EnumType,
         const char*(enumValueToStringFunction)(EnumType),
         StringCompareType stringCompareType = StringCompareType::Exact,
         StripStrings      stripStrings      = StripStrings::No>
class StringToEnumValueConverter
{
public:
    StringToEnumValueConverter() : stringToEnumValue_(stringCompareType)
    {
        for (const auto type : EnumerationWrapper<EnumType>{})
        {
            GMX_RELEASE_ASSERT(type != EnumType::Count,
                               "EnumerationWrapper<EnumType> should never return EnumType::Count");
            std::string stringFromType = enumValueToStringFunction(type);
            if (stripStrings == StripStrings::Yes)
            {
                // Ensure leading and trailing spaces are removed
                stringFromType = stripString(stringFromType);
            }
            stringToEnumValue_[stringFromType] = type;
        }
    }

    //! Return an optional enum value identified from the \c s (which is never EnumType::Count)
    std::optional<EnumType> valueFrom(const std::string& s) const
    {
        typename decltype(stringToEnumValue_)::const_iterator typeIt;
        if (stripStrings == StripStrings::Yes)
        {
            // Ensure leading and trailing spaces are removed
            typeIt = stringToEnumValue_.find(stripString(s));
        }
        else
        {
            typeIt = stringToEnumValue_.find(s);
        }
        return (typeIt != stringToEnumValue_.end()) ? std::make_optional(typeIt->second) : std::nullopt;
    }

private:
    /*! \brief Map of strings to enumeration values.
     *
     * By construction, those values never include
     * EnumType::Count. */
    std::map<std::string, EnumType, StringCompare> stringToEnumValue_;
};

} // namespace gmx

#endif
