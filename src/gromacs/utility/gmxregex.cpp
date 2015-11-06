/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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
 * Implements regular expression wrappers.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gmxregex.h"

#include <regex>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Helper function to transform the std error code into a string that gmx::Regex can report
std::string getRegexErrorString(std::regex_constants::error_type code)
{
    using namespace std::regex_constants;

    std::string errorString;
    switch (code)
    {
        case error_collate:
            errorString = "The expression contained an invalid collating element name.";
            break;
        case error_ctype:
            errorString = "The expression contained an invalid character class name.";
            break;
        case error_escape:
            errorString = "The expression contained an invalid escaped character, or a trailing escape.";
            break;
        case error_backref:
            errorString = "The expression contained an invalid back reference.";
            break;
        case error_brack:
            errorString = "The expression contained mismatched brackets ([ and ]).";
            break;
        case error_paren:
            errorString = "The expression contained mismatched parentheses (( and )).";
            break;
        case error_brace:
            errorString = "The expression contained mismatched braces ({ and }).";
            break;
        case error_badbrace:
            errorString = "The expression contained an invalid range between braces ({ and }).";
            break;
        case error_range:
            errorString = "The expression contained an invalid character range.";
            break;
        case error_space:
            errorString = "There was insufficient memory to convert the expression into a finite state machine.";
            break;
        case error_badrepeat:
            errorString = "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression.";
            break;
        case error_complexity:
            errorString = "The complexity of an attempted match against a regular expression exceeded a pre-set level.";
            break;
        case error_stack:
            errorString = "There was insufficient memory to determine whether the regular expression could match the specified character sequence.";
            break;
        default:
            break;
    }
    return errorString;
}

class Regex::Impl
{
    public:
        explicit Impl(const char *value)
        try : regex_(value, std::regex::nosubs | std::regex::extended)
        {
        }
        catch (const std::regex_error &ex)
        {
            GMX_THROW(InvalidInputError(formatString("Error in regular expression \"%s\": %s",
                                                     value, getRegexErrorString(ex.code()).c_str())));
        }

        explicit Impl(const std::string &value)
        try : regex_(value, std::regex::nosubs | std::regex::extended)
        {
        }
        catch (const std::regex_error &ex)
        {
            GMX_THROW(InvalidInputError(formatString("Error in regular expression \"%s\": %s",
                                                     value.c_str(), getRegexErrorString(ex.code()).c_str())));
        }

        bool match(const char *value) const
        {
            try
            {
                return std::regex_match(value, regex_);
            }
            catch (const std::regex_error &)
            {
                // TODO: Handle errors.
                return false;
            }
        }

    private:
        std::regex              regex_;
};

Regex::Regex()
    : impl_(nullptr)
{
}

Regex::Regex(const char *value)
    : impl_(new Impl(value))
{
}

Regex::Regex(const std::string &value)
    : impl_(new Impl(value))
{
}

Regex::~Regex()
{
}

/*! \cond libapi */
bool regexMatch(const char *str, const Regex &regex)
{
    if (regex.impl_ == nullptr)
    {
        return false;
    }
    return regex.impl_->match(str);
}

bool regexMatch(const std::string &str, const Regex &regex)
{
    return regexMatch(str.c_str(), regex);
}
//! \endcond

} // namespace gmx
