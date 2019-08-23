/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \file
 * \brief C++ exceptions for the gmxapi compatibility tools.
 *
 * Used internally for the gmxapi compatibility helpers that manage type
 * mappings of older GROMACS structures. The long-term disposition of this
 * code is uncertain, but the headers are not likely to be public. If they
 * do persist in some form, we can integrate the exception hierarchy into
 * whatever module takes ownership of this code.
 *
 * Exceptions defined here should only be caught by code that understands the
 * implementation details of these compatibility tools. Exposure of these
 * exceptions outside of the installed object files should be treated as a bug.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_compat
 */

#ifndef GMXAPICOMPAT_EXCEPTIONS_H
#define GMXAPICOMPAT_EXCEPTIONS_H

#include <exception>
#include <string>

namespace gmxapicompat
{

/*!
 * \brief Generic exception class for gmxapicompat.
 */
class Exception : public std::exception
{
    public:
        using std::exception::exception;

        explicit Exception(const std::string &message) :
            message_ {message}
        {}
        explicit Exception(const char* message) : Exception(std::string(message)) {}

        const char *what() const noexcept override
        {
            return message_.c_str();
        }

    private:
        std::string message_;
};

/*!
 * \brief The key name provided for a key-value mapping is invalid.
 */
class KeyError : public Exception
{
    using Exception::Exception;
};

/*!
 * \brief The value provided for a key-value mapping is invalid.
 */
class ValueError : public Exception
{
    using Exception::Exception;
};

/*!
 * \brief Value provided for a key-value mapping is of an incompatible type.
 */
class TypeError : public Exception
{
    using Exception::Exception;
};

}      // end namespace gmxapicompat
#endif //GMXAPICOMPAT_EXCEPTIONS_H
