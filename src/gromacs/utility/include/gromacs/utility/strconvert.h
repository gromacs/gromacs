/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares common utility functions for conversions to and from strings.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_STRCONVERT_H
#define GMX_UTILITY_STRCONVERT_H

#include <algorithm>
#include <array>
#include <optional>
#include <string>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! \cond libapi
//! \addtogroup module_utility
//! \{

/*! \brief
 * Parses a boolean from a string.
 *
 * \throws  InvalidInputError if `str` is not recognized as a boolean value.
 */
bool boolFromString(const char* str);
/*! \brief
 * Parses an integer from a string.
 *
 * \throws  InvalidInputError if `str` is not a valid integer.
 *
 * Also checks for overflow.
 */
int intFromString(const char* str);
/*! \brief
 * Parses a 64-bit integer from a string.
 *
 * \throws  InvalidInputError if `str` is not a valid integer.
 *
 * Also checks for overflow.
 */
int64_t int64FromString(const char* str);
/*! \brief
 * Parses a float value from a string.
 *
 * \throws  InvalidInputError if `str` is not a valid number.
 *
 * Also checks for overflow.
 */
float floatFromString(const char* str);
/*! \brief
 * Parses a double value from a string.
 *
 * \throws  InvalidInputError if `str` is not a valid number.
 *
 * Also checks for overflow.
 */
double doubleFromString(const char* str);

/*! \brief
 * Parses a value from a string to a given type.
 *
 * \tparam T Type of value to parse.
 *
 * `T` can only be one of the types that is explicity supported.
 * The main use for this function is to write `fromString<real>(value)`,
 * but it can also be used for other types for consistency.
 */
template<typename T>
static inline T fromString(const char* str);
//! \copydoc fromString(const char *)
template<typename T>
static inline T fromString(const std::string& str)
{
    return fromString<T>(str.c_str());
}
/*! \copydoc fromString(const char *)
 *
 * Provided for situations where overload resolution cannot easily resolve the
 * desired std::string parameter.
 */
template<typename T>
static inline T fromStdString(const std::string& str)
{
    return fromString<T>(str.c_str());
}

//! Implementation for boolean values.
template<>
inline bool fromString<bool>(const char* str)
{
    return boolFromString(str);
}
//! Implementation for integer values.
template<>
inline int fromString<int>(const char* str)
{
    return intFromString(str);
}
//! Implementation for 64-bit integer values.
template<>
inline int64_t fromString<int64_t>(const char* str)
{
    return int64FromString(str);
}
//! Implementation for float values.
template<>
inline float fromString<float>(const char* str)
{
    return floatFromString(str);
}
//! Implementation for double values.
template<>
inline double fromString<double>(const char* str)
{
    return doubleFromString(str);
}

/*! \brief
 * Converts a boolean to a "true"/"false" string.
 *
 * Does not throw.
 */
static inline const char* boolToString(bool value)
{
    return value ? "true" : "false";
}
/*! \brief
 * Returns a string containing the value of \c t.
 *
 * \throws std::bad_alloc if out of memory.
 */
static inline std::string intToString(int t)
{
    return formatString("%d", t);
}
//! \copydoc intToString(int)
static inline std::string int64ToString(int64_t t)
{
    return formatString("%" PRId64, t);
}
//! \copydoc intToString(int)
static inline std::string doubleToString(double t)
{
    return formatString("%g", t);
}

/*! \name
 * Overloads for converting a value of a given type to a string.
 *
 * \throws std::bad_alloc if out of memory.
 * \{
 */
static inline std::string toString(bool t)
{
    return boolToString(t);
}
static inline std::string toString(int t)
{
    return intToString(t);
}
static inline std::string toString(int64_t t)
{
    return int64ToString(t);
}
static inline std::string toString(float t)
{
    return doubleToString(t);
}
static inline std::string toString(double t)
{
    return doubleToString(t);
}
static inline std::string toString(std::string t)
{
    return t;
}
//! \}

//! \}
//! \endcond

/*! \brief Convert a string into an array of values.
 *
 * \tparam ValueType array element type to convert into
 * \tparam NumExpectedValues number of values of the array
 *
 * \returns an array containing the converted string, optionally null if
 *          the white-space stripped string is empty
 *
 * \throws InvalidInputError if splitting the string at whitespaces does not
 *                           result in NumExpectedValues or zero substrings
 *
 * \throws InvalidInputError if conversion of any of the NumExpectedValues
 *                           substrings of the splitted input string fails
 *
 * Converts a string into an array of type ValueType with exactly NumExpectedValues.
 *
 * No result is returned if the string is empty or contains only whitespace .
 *
 */
template<typename ValueType, int NumExpectedValues>
static inline std::optional<std::array<ValueType, NumExpectedValues>>
parsedArrayFromInputString(const std::string& str)
{
    // return nullopt right away if the string is just whitespace or empty
    {
        const std::string& strippedString = stripString(str);
        if (strippedString.empty())
        {
            return std::nullopt;
        }
    }

    const std::vector<std::string>& valuesAsStrings = splitString(str);

    // throw right away if we don't have the expected number of string entries
    if (valuesAsStrings.size() != NumExpectedValues)
    {
        const std::string errorMessage =
                "Expected empty string or string with " + intToString(NumExpectedValues)
                + " elements to convert, but received " + intToString(valuesAsStrings.size())
                + " elements instead.";
        GMX_THROW(InvalidInputError(errorMessage));
    }

    // will throw if any conversion from string to value fails
    std::array<ValueType, NumExpectedValues> valuesAsArray;
    for (int i = 0; i < NumExpectedValues; i++)
    {
        valuesAsArray[i] = fromString<ValueType>(valuesAsStrings[i]);
    }

    return { valuesAsArray };
}

/*! \brief Returns the input string, throwing an excpetion if the demanded
 *         conversion to an array will not succeed.
 *
 * \tparam ValueType array element type to convert into
 * \tparam NumExpectedValues number of values of the array
 *
 * \param[in] toConvert the string to convert
 * \param[in] errorContextMessage the message to add to the thrown exceptions if
 *                                conversion of the string is bound to fail at
 *                                some point
 * \returns the input string
 *
 * \throws InvalidInputError if splitting the string at whitespaces does not
 *                           result in NumExpectedValues or zero substrings
 *
 * \throws InvalidInputError if conversion of any of the NumExpectedValues
 *                           substrings of the splitted input string fails
 *
 * A typical use of this function would be in .mdp string option parsing
 * where information in the .mdp file is transformed into the data that is
 * stored in the .tpr file.
 */
template<typename ValueType, int NumExpectedValues>
static inline std::string stringIdentityTransformWithArrayCheck(const std::string& toConvert,
                                                                const std::string& errorContextMessage)
{
    // Attempt the conversion to an array so that the string parsing routine
    // will throw an InvalidInputError if the string is not fit for conversion.
    try
    {
        // The converted array is discarded.
        gmx_unused const auto& val = parsedArrayFromInputString<ValueType, NumExpectedValues>(toConvert);
    }
    catch (const InvalidInputError& e)
    {
        InvalidInputError toThrow(errorContextMessage + std::string(e.what()));
        GMX_THROW(toThrow);
    }

    return toConvert;
}

} // namespace gmx

#endif
