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
/*! \file
 * \brief
 * Declares common exception classes for fatal error handling.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_fatalerror
 */
#ifndef GMX_FATALERROR_EXCEPTIONS_H
#define GMX_FATALERROR_EXCEPTIONS_H

#include <exception>
#include <string>

#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/throw_exception.hpp>

namespace gmx
{

/*! \brief
 * Stores a user-friendly explanation for the reason of an exception.
 *
 * Typically, should not be used directly, but through the GromacsException
 * class: it is initialized by the constructor, and can be accessed with
 * GromacsException::what().
 *
 * \inlibraryapi
 */
typedef boost::error_info<struct errinfo_message_, std::string> errinfo_message;

/*! \addtopublicapi
 * \{
 */

/*! \brief
 * Base class for all exception objects in Gromacs.
 *
 * Although boost recommends using virtual inheritance in exception hiearchies,
 * it is not used here for two reasons:
 * -# It is only useful when there is diamond inheritance, and that should
 *    never occur in this exception hierarchy because this class is the only
 *    instance of multiple inheritance (Gromacs programming guidelines prohibit
 *    multiple inheritance from concrete classes, but it is unavoidable here
 *    because of the design of boost::exception).
 * -# Because the constructor takes an argument, virtual inheritance would
 *    complicate any classes that inherit indirectly from this class.
 *
 * \ingroup module_fatalerror
 */
class GromacsException : public std::exception, public boost::exception
{
    public:
        /*! \brief
         * Returns the reason string for the exception.
         *
         * The return value is the string that was passed to the constructor.
         */
        virtual const char *what() const throw();
        /*! \brief
         * Returns the error code corresponding to the exception type.
         */
        virtual int errorCode() const = 0;

    protected:
        /*! \brief
         * Creates an exception object with the provided detailed reason.
         *
         * \param[in] reason Detailed reason for the exception.
         */
        explicit GromacsException(const std::string &reason);
};

/*! \brief
 * Exception class for file I/O errors.
 *
 * \ingroup module_fatalerror
 */
class FileIOError : public GromacsException
{
    public:
        /*! \brief
         * Creates an exception object with the provided detailed reason.
         *
         * \param[in] reason Detailed reason for the exception.
         */
        explicit FileIOError(const std::string &reason)
            : GromacsException(reason) {}

        virtual int errorCode() const;
};

/*! \brief
 * Exception class for user input errors.
 *
 * Derived classes should be used to indicate the nature of the error instead
 * of throwing this class directly.
 *
 * \ingroup module_fatalerror
 */
class UserInputError : public GromacsException
{
    protected:
        //! \copydoc FileIOError::FileIOError()
        explicit UserInputError(const std::string &reason)
            : GromacsException(reason) {}
};

/*! \brief
 * Exception class for situations where user input cannot be parsed/understood.
 *
 * \ingroup module_fatalerror
 */
class InvalidInputError : public UserInputError
{
    public:
        //! \copydoc FileIOError::FileIOError()
        explicit InvalidInputError(const std::string &reason)
            : UserInputError(reason) {}

        virtual int errorCode() const;
};

/*! \brief
 * Exception class for situations where user input is inconsistent.
 *
 * \ingroup module_fatalerror
 */
class InconsistentInputError : public UserInputError
{
    public:
        //! \copydoc FileIOError::FileIOError()
        explicit InconsistentInputError(const std::string &reason)
            : UserInputError(reason) {}

        virtual int errorCode() const;
};

/*! \brief
 * Exception class for simulation instabilities.
 *
 * \ingroup module_fatalerror
 */
class SimulationInstabilityError : public GromacsException
{
    public:
        //! \copydoc FileIOError::FileIOError()
        explicit SimulationInstabilityError(const std::string &reason)
            : GromacsException(reason) {}

        virtual int errorCode() const;
};

/*! \brief
 * Exception class for internal errors.
 *
 * \ingroup module_fatalerror
 */
class InternalError : public GromacsException
{
    public:
        //! \copydoc FileIOError::FileIOError()
        explicit InternalError(const std::string &reason)
            : GromacsException(reason) {}

        virtual int errorCode() const;
};

/*! \brief
 * Exception class for incorrect use of an API.
 *
 * \ingroup module_fatalerror
 */
class APIError : public GromacsException
{
    public:
        //! \copydoc FileIOError::FileIOError()
        explicit APIError(const std::string &reason)
            : GromacsException(reason) {}

        virtual int errorCode() const;
};

/*! \brief
 * Exception class for use of an unimplemented feature.
 *
 * \ingroup module_fatalerror
 */
class NotImplementedError : public APIError
{
    public:
        //! \copydoc FileIOError::FileIOError()
        explicit NotImplementedError(const std::string &reason)
            : APIError(reason) {}

        virtual int errorCode() const;
};


/*! \brief
 * Macro for throwing an exception.
 *
 * \param[in] e    Exception object to throw.
 *
 * Using this macro instead of \c throw directly makes it possible to uniformly
 * attach information into the exception objects.
 * \p e should evaluate to an instance of an object derived from
 * GromacsException.
 *
 * Basic usage:
 * \code
if (value < 0)
{
    GMX_THROW(InconsistentUserInput("Negative values not allowed for value"));
}
 * \endcode
 */
#define GMX_THROW(e) \
    BOOST_THROW_EXCEPTION((e))

/*! \brief
 * Formats a standard error message for reporting an error.
 *
 * Normal usage in Gromacs command-line programs is like this:
 * \code
int main(int argc, char *argv[])
{
    try
    {
        // The actual code for the program
        return 0;
    }
    catch (const std::exception &ex)
    {
        fprintf(stderr, "%s", gmx::formatErrorMessage(ex).c_str());
        return 1;
    }
}
 * \endcode
 */
std::string formatErrorMessage(const std::exception &ex);

/*! \brief
 * Converts an exception into a return code.
 */
int translateException(const std::exception &ex);

/*!\}*/

} // namespace gmx

#endif
