/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#include <string>
#include <utility>


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
    Exception(const Exception& /*unused*/);
    Exception& operator=(const Exception& /*unused*/);

    Exception(Exception&& /*unused*/) noexcept;
    Exception& operator=(Exception&& /*unused*/) noexcept;

    [[nodiscard]] const char* what() const noexcept override;
    //! \endcond
};

/*!
 * \brief Basic implementation mix-in for exceptions.
 *
 * Allow exceptions to be defined with minimal syntax when their primary function is
 * to exist as distinct named types.
 *
 * \tparam E the class using this template as a base class.
 *
 * Use in the "curiously recurring template pattern".
 *
 * Example:
 *
 *     class DerivedException : public BasicException<DerivedException>
 *     {
 *          public:
 *              using BasicException::BasicException;
 *     };
 *
 * \note Current implementation only provides constructors and no specialized or dispatched behavior.
 *
 * \ingroup gmxapi_exceptions
 */
template<class E>
class BasicException : public Exception
{
private:
    //! Store the usual exception message as a std::string instead of C string.
    std::string what_;

public:
    //! Initialize with empty message.
    BasicException() : BasicException{ std::string() } {}

    /*!
     * \brief Copy a string to use for the exception message.
     *
     * \param message
     * \{
     */
    explicit BasicException(std::string message) noexcept : what_{ std::move(message) } {}

    explicit BasicException(const char* message) { what_ = message; }
    //! \}

    /*!
     * \brief Get message.
     *
     * \return pointer to C string.
     *
     * It is the responsibility of the caller to keep the Exception object alive while the char
     * pointer is in use.
     */
    [[nodiscard]] const char* what() const noexcept override { return what_.c_str(); }
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
class ProtocolError : public BasicException<ProtocolError>
{
public:
    using BasicException<ProtocolError>::BasicException;
};

/*!
 * \deprecated \brief Intended feature is not implemented (Deprecated).
 *
 * This exception is no longer used.
 * See, instead, MissingImplementationError and \issue{4337}.
 *
 * \ingroup gmxapi_exceptions
 */
class NotImplementedError : public BasicException<NotImplementedError>
{
public:
    using BasicException<NotImplementedError>::BasicException;
};

/*!
 * \brief Intended feature is not implemented.
 *
 * Indicates a bug in the API implementation. Either a version mismatch between the client
 * and library has gone undetected, or the API has purported to offer functionality that does
 * not exist.
 *
 * \version 0.3.1
 * This exception replaces `gmxapi::NotImplementedError`. See \issue{4337}.
 *
 * \ingroup gmxapi_exceptions
 */
class MissingImplementationError : public BasicException<MissingImplementationError>
{
public:
    using BasicException<MissingImplementationError>::BasicException;
};

/*!
 * \brief Unacceptable API usage.
 *
 * Usage error. Examples include calling a function with invalid arguments.
 *
 * \ingroup gmxapi_exceptions
 */
class UsageError : public BasicException<UsageError>
{
public:
    using BasicException<UsageError>::BasicException;
};

} // end namespace gmxapi

#endif // header guard
