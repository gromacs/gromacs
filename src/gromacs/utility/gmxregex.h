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
/*! \libinternal \file
 * \brief
 * Declares simple wrapper for regular expression functionality.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GMXREGEX_H
#define GMX_UTILITY_GMXREGEX_H

#include <string>

#include "common.h"

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
 * Currently, isSupported() will return false if POSIX regular expression
 * header is not available (i.e., on Windows).  In this case, calling other
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

