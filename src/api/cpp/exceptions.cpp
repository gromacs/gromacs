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

#include "gmxapi/exceptions.h"

#include <string>

namespace gmxapi
{

//! \cond
Exception::Exception() = default;

Exception::~Exception() = default;

Exception::Exception(const Exception &) = default;

Exception &Exception::operator=(const Exception &) = default;

Exception::Exception(Exception &&) noexcept = default;

Exception &Exception::operator=(Exception &&) noexcept = default;

const char* Exception::what() const noexcept
{
    return "Gromacs API error";
}
//! \endcond

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
        //! Store the usual exception method as a std::string instead of C string.
        std::string what_;
    public:
        //! Initialize with empty message.
        BasicException() : BasicException{std::string()}
        {}

        /*!
         * \brief Copy a string to use for the exception message.
         *
         * \param message
         * \{
         */
        explicit BasicException(std::string message) noexcept :
            what_ {std::move(message)}
        {}

        explicit BasicException(const char* message)
        {
            what_ = message;
        }
        //! \}

        /*!
         * \brief Get message.
         *
         * \return pointer to C string.
         *
         * It is the responsibility of the caller to keep the Exception object alive while the char
         * pointer is in use.
         */
        const char* what() const noexcept override
        {
            return what_.c_str();
        }
};

class ProtocolError : public BasicException<ProtocolError>
{
    public:
        using BasicException<ProtocolError>::BasicException;
};

class NotImplementedError : public BasicException<NotImplementedError>
{
    public:
        using BasicException<NotImplementedError>::BasicException;
};

} // end namespace gmxapi
