/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares simple wrapper for regular expression functionality.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GMXREGEX_H
#define GMX_UTILITY_GMXREGEX_H

#include <string>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class Regex;

/*! \cond libapi */
/*! \libinternal \brief
 * Matches a string with a regular expression.
 *
 * \param[in] str   String to match.
 * \param[in] regex Regular expression to match.
 * \returns   true if \p regex matches the whole \p str.
 *
 * Does not throw currently, but this is subject to change if/when better error
 * handling is implemented (currently, it returns false if the matching fails,
 * e.g., because of out-of-memory).
 *
 * \ingroup module_utility
 */
bool regexMatch(const char *str, const Regex &regex);
//! \copydoc regexMatch(const char *, const Regex &)
bool regexMatch(const std::string &str, const Regex &regex);
//! \endcond

/*! \libinternal \brief
 * Represents a regular expression.
 *
 * This class provides a simple interface for regular expression construction.
 * regexMatch() is used to match the regular expression against a string.
 * POSIX extended regular expression syntax is used.
 *
 * Currently, isSupported() will return true if either
 *
 *  -# POSIX regular expression header <regex.h> is available, or
 *  -# C++11 header \<regex> is available (e.g., new enough MSVC has this).
 *
 * In other cases, isSupported() returns false and calling other
 * constructors than the default constructor throws an exception.
 *
 * \see regexMatch()
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class Regex
{
    public:
        /*! \brief
         * Returns true if regular expression support has been compiled in.
         *
         * Does not throw.
         */
        static bool isSupported();

        /*! \brief
         * Constructs a regular expression that matches nothing.
         *
         * Does not throw.
         */
        Regex();
        /*! \brief
         * Constructs a regular expression from a string.
         *
         * \param[in] value  String to compile into a regular expression.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if \p value is not a valid regular
         *      expression.
         *
         * \todo
         * Consider whether some other exception type would be better.
         */
        explicit Regex(const char *value);
        //! \copydoc Regex(const char *)
        explicit Regex(const std::string &value);
        //! Frees memory allocated for the regular expression.
        ~Regex();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        friend bool regexMatch(const char *str, const Regex &regex);
        friend bool regexMatch(const std::string &str, const Regex &regex);
};

} // namespace gmx

#endif
