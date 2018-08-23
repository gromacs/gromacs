/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMXAPI_EXCEPTIONS_H
#define GMXAPI_EXCEPTIONS_H
/*! \defgroup gmxapi_exceptions Exceptions
 *
 * \brief Exceptions thrown by gmxapi components.
 *
 * \ingroup gmxapi
 */
/*! \file
 * \brief Declare exception classes for operations in gmxapi.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup gmxapi_exceptions
 */

#include <exception>

namespace gmxapi
{

/*! \brief Base exception for gmxapi library.
 *
 * Exceptions thrown in the gmxapi namespace are descended from gmxapi::Exception
 * or there is a bug.
 *
 * \ingroup gmxapi_exceptions
 */
class Exception : public std::exception
{
    public:
        //! \cond
        Exception();
        ~Exception() override;
        Exception(const Exception &);
        Exception &operator=(const Exception &);

        Exception(Exception &&) noexcept;
        Exception &operator=(Exception &&) noexcept;

        const char* what() const noexcept override;
        //! \endcond
};

/*! \brief Behavioral protocol violated.
 *
 * Indicates that a behavioral protocol specified in the API is not being followed. The class
 * throwing this exception expects certain methods to be called in a certain order.
 *
 * If this exception is encountered in client code, the API is being misused or there is a bug.
 * (Generally, required behaviors should be implemented in templates or base classes rather than
 * exposing and requiring complete implementation of the protocol in client code.)
 *
 * \ingroup gmxapi_exceptions
 */
class ProtocolError;

/*!
 * \brief Intended feature is not implemented.
 *
 * Indicates a bug in the API implementation. Either a version mismatch between the client
 * and library has gone undetected, or the API has purported to offer functionality that does
 * not exist.
 *
 * \ingroup gmxapi_exceptions
 */
class NotImplementedError;

}      // end namespace gmxapi

#endif // header guard
