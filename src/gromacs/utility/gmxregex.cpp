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
 *
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
 * Implements regular expression wrappers.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "gmxregex.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(HAVE_POSIX_REGEX)
// old Mac needs sys/types.h before regex.h
#include <sys/types.h>
#include <regex.h>
#define USE_POSIX_REGEX
#elif defined(HAVE_CXX11_REGEX)
#include <regex>
#define USE_CXX11_REGEX
#endif

#include "exceptions.h"
#include "stringutil.h"

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

        regex_t regex_;
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
        std::regex regex_;
};
#else
class Regex::Impl
{
    public:
        explicit Impl(const char * /*value*/)
        {
            GMX_THROW(NotImplementedError(
                          "Gromacs is compiled without regular expression support"));
        }
        explicit Impl(const std::string & /*value*/)
        {
            GMX_THROW(NotImplementedError(
                          "Gromacs is compiled without regular expression support"));
        }

        bool match(const char * /*value*/) const
        {
            // Should never be reached.
            GMX_THROW(NotImplementedError(
                          "Gromacs is compiled without regular expression support"));
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
