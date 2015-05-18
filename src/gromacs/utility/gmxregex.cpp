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

#include "config.h"

#if defined(HAVE_POSIX_REGEX)
#include <sys/types.h>
// old Mac needs sys/types.h before regex.h
#include <regex.h>
#define USE_POSIX_REGEX
#elif defined(HAVE_CXX11_REGEX)
#include <regex>
#define USE_CXX11_REGEX
#endif

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

// static
bool Regex::isSupported()
{
#if defined(USE_POSIX_REGEX) || defined(USE_CXX11_REGEX)
    return true;
#else
    return false;
#endif
}

#if defined(USE_POSIX_REGEX)
class Regex::Impl
{
    public:
        explicit Impl(const char *value)
        {
            compile(value);
        }
        explicit Impl(const std::string &value)
        {
            compile(value.c_str());
        }
        ~Impl()
        {
            regfree(&regex_);
        }

        bool match(const char *value) const
        {
            int rc = regexec(&regex_, value, 0, NULL, 0);
            if (rc != 0 && rc != REG_NOMATCH)
            {
                // TODO: Handle errors.
            }
            return (rc == 0);
        }

    private:
        void compile(const char *value)
        {
            std::string buf(formatString("^%s$", value));
            int         rc = regcomp(&regex_, buf.c_str(), REG_EXTENDED | REG_NOSUB);
            if (rc != 0)
            {
                // TODO: Better error messages.
                GMX_THROW(InvalidInputError(formatString(
                                                    "Error in regular expression \"%s\"", value)));
            }
        }

        regex_t                 regex_;
};
#elif defined(USE_CXX11_REGEX)
class Regex::Impl
{
    public:
        explicit Impl(const char *value)
        try : regex_(value, std::regex::nosubs | std::regex::extended)
        {
        }
        catch (const std::regex_error &)
        {
            // TODO: Better error messages.
            GMX_THROW(InvalidInputError(formatString(
                                                "Error in regular expression \"%s\"", value)));
        }
        explicit Impl(const std::string &value)
        try : regex_(value, std::regex::nosubs | std::regex::extended)
        {
        }
        catch (const std::regex_error &)
        {
            // TODO: Better error messages.
            GMX_THROW(InvalidInputError(formatString(
                                                "Error in regular expression \"%s\"", value)));
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
#else
class Regex::Impl
{
    public:
        explicit Impl(const char * /*value*/)
        {
            GMX_THROW(NotImplementedError(
                              "GROMACS is compiled without regular expression support"));
        }
        explicit Impl(const std::string & /*value*/)
        {
            GMX_THROW(NotImplementedError(
                              "GROMACS is compiled without regular expression support"));
        }

        bool match(const char * /*value*/) const
        {
            // Should never be reached.
            GMX_THROW(NotImplementedError(
                              "GROMACS is compiled without regular expression support"));
        }
};
#endif

Regex::Regex()
    : impl_(NULL)
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
    if (regex.impl_.get() == NULL)
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
