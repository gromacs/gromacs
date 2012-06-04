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
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_EXCEPTIONS_H
#define GMX_UTILITY_EXCEPTIONS_H

#include <cstdio>
#include <cstdlib>

#include <exception>
#include <string>

#include <boost/exception/errinfo_api_function.hpp>
#include <boost/exception/errinfo_errno.hpp>
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * \ingroup module_utility
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
 * Macro for throwing an exception based on errno.
 *
 * \param[in] e       Exception object to throw.
 * \param[in] syscall Name of the syscall that returned the error.
 * \param[in] err     errno value returned by the syscall.
 *
 * This macro provides a convenience interface for throwing an exception to
 * report an error based on a errno value.  In addition to adding the necessary
 * information to the exception object, the macro also ensures that \p errno is
 * evaluated before, e.g., the constructor of \p e may call other functions
 * that could overwrite the errno value.
 * \p e should evaluate to an instance of an object derived from
 * GromacsException.
 *
 * Typical usage (note that gmx::File wraps this particular case):
 * \code
FILE *fp = fopen("filename.txt", "r");
if (fp == NULL)
{
    GMX_THROW(FileIOError("Could not open file"), "fopen", errno);
}
 * \endcode
 */
#define GMX_THROW_WITH_ERRNO(e, syscall, err) \
    do { \
        int stored_errno_ = (err); \
        GMX_THROW((e) << boost::errinfo_errno(stored_errno_) \
                      << boost::errinfo_api_function(syscall)); \
    } while(0)

/*! \brief
 * Formats a standard fatal error message for reporting an exception.
 *
 * Does not throw.  If memory allocation fails or some other error occurs
 * while formatting the error, tries to print a reasonable alternative message.
 *
 * Normal usage in Gromacs command-line programs is like this:
 * \code
int main(int argc, char *argv[])
{
    gmx::ProgramInfo::init(argc, argv);
    try
    {
        // The actual code for the program
        return 0;
    }
    catch (const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return 1;
    }
}
 * \endcode
 */
void printFatalErrorMessage(FILE *fp, const std::exception &ex);

/*! \brief
 * Converts an exception into a return code.
 */
int translateException(const std::exception &ex);

/*!\}*/

/*! \cond libapi */
/*! \libinternal \brief
 * Macro for catching exceptions at C++ -> C boundary.
 *
 * This macro is intended for uniform handling of exceptions when C++ code is
 * called from C code within Gromacs.  Since most existing code is written
 * using the assumption that fatal errors terminate the program, this macro
 * implements this behavior for exceptions.  It should only be used in cases
 * where the error cannot be propagated upwards using return values or such.
 *
 * Having this as a macro instead of having the same code in each place makes
 * it easy to 1) find all such locations in the code, and 2) change the exact
 * behavior if needed.
 *
 * Usage:
 * \code
try
{
    // C++ code
}
GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
 * \endcode
 *
 * \inlibraryapi
 */
#define GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR \
    catch (const std::exception &ex) { \
        ::gmx::printFatalErrorMessage(stderr, ex); \
        std::exit(1); \
    }
//! \endcond

} // namespace gmx

#endif
