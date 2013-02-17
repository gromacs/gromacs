/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Implements classes and functions in exceptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "exceptions.h"

#include <cstring>

#include <new>
#include <stdexcept>
#include <typeinfo>

#include <boost/exception/get_error_info.hpp>
#include <boost/shared_ptr.hpp>

#include "gromacs/legacyheaders/thread_mpi/system_error.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/gmxassert.h"

#include "errorformat.h"

namespace gmx
{

namespace
{

/********************************************************************
 * ErrorMessage
 */

class ErrorMessage
{
    public:
        /*! \brief
         * Creates an error message object with the specified text.
         *
         * \param[in] text  Text for the message.
         */
        explicit ErrorMessage(const std::string &text);

        //! Whether this object is a context string.
        bool isContext() const { return child_; }
        //! Returns the text for this object.
        const std::string &text() const { return text_; }
        /*! \brief
         * Returns the child object for a context object.
         *
         * Must not be called if isContext() returns false.
         */
        const ErrorMessage &child() const
        {
            GMX_ASSERT(isContext(),
                       "Attempting to access nonexistent message object");
            return *child_;
        }

        /*! \brief
         * Creates a new message object with context prepended.
         *
         * \param[in] context  Context string to add.
         * \returns   New error message object that has \p context as its text
         *      and \c this as its child.
         * \throws    std::bad_alloc if out of memory.
         */
        ErrorMessage prependContext(const std::string &context) const;

    private:
        std::string                     text_;
        boost::shared_ptr<ErrorMessage> child_;
};

/*! \internal \brief
 * Stores a reason or the top-most context string of an exception.
 *
 * \ingroup module_utility
 */
typedef boost::error_info<struct errinfo_message_, ErrorMessage>
    errinfo_message;

ErrorMessage::ErrorMessage(const std::string &text)
    : text_(text)
{
    size_t length = text_.find_last_not_of(" \n");
    if (length == std::string::npos)
    {
        length = text_.length() - 1;
    }
    text_.resize(length + 1);
}

ErrorMessage
ErrorMessage::prependContext(const std::string &context) const
{
    ErrorMessage newMessage(context);
    newMessage.child_.reset(new ErrorMessage(*this));
    return newMessage;
}

/*! \brief
 * Stores list of nested exceptions for Gromacs exceptions.
 *
 * \ingroup module_utility
 */
typedef boost::error_info<struct errinfo_message_, internal::NestedExceptionList>
    errinfo_nested_exceptions;

}   // namespace

/********************************************************************
 * GromacsException
 */

GromacsException::GromacsException(const ExceptionInitializer &details)
{
    *this << errinfo_message(ErrorMessage(details.reason_));
    if (details.hasNestedExceptions())
    {
        *this << errinfo_nested_exceptions(details.nested_);
    }
}

const char *GromacsException::what() const throw()
{
    const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*this);
    while (msg != NULL && msg->isContext())
    {
        msg = &msg->child();
    }
    return msg != NULL ? msg->text().c_str() : "No reason provided";
}

void GromacsException::prependContext(const std::string &context)
{
    const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*this);
    GMX_RELEASE_ASSERT(msg != NULL, "Message should always be set");
    *this << errinfo_message(msg->prependContext(context));
}

/********************************************************************
 * Derived exception classes
 */

int FileIOError::errorCode() const
{
    return eeFileIO;
}

int InvalidInputError::errorCode() const
{
    return eeInvalidInput;
}

int InconsistentInputError::errorCode() const
{
    return eeInconsistentInput;
}

int SimulationInstabilityError::errorCode() const
{
    return eeInstability;
}

int InternalError::errorCode() const
{
    return eeInternalError;
}

int APIError::errorCode() const
{
    return eeAPIError;
}

int NotImplementedError::errorCode() const
{
    return eeNotImplemented;
}


/********************************************************************
 * Global functions
 */

namespace
{

/*! \brief
 * Prints error information for an exception object.
 *
 * \param[in] fp      File to write the information out to (typically stderr).
 * \param[in] ex      Exception object to print.
 * \param[in] indent  Indentation for the information.
 *
 * If the exception contains nested exceptions, information from them is
 * recursively printed.
 */
void printExceptionMessage(FILE *fp, const std::exception &ex, int indent)
{
    const boost::exception *boostEx = dynamic_cast<const boost::exception *>(&ex);
    if (boostEx != NULL)
    {
        // TODO: Remove duplicate context if present in multiple nested exceptions.
        const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*boostEx);
        if (msg != NULL)
        {
            while (msg != NULL && msg->isContext())
            {
                internal::printFatalErrorMessageLine(fp, msg->text().c_str(), indent*2);
                ++indent;
                msg = &msg->child();
            }
            if (msg != NULL && !msg->text().empty())
            {
                internal::printFatalErrorMessageLine(fp, msg->text().c_str(), indent*2);
            }
        }
        else
        {
            internal::printFatalErrorMessageLine(fp, ex.what(), indent);
        }

        const int *errorNumber
            = boost::get_error_info<boost::errinfo_errno>(*boostEx);
        if (errorNumber != NULL)
        {
            std::fprintf(fp, "%*sReason: %s\n", (indent+1)*2, "",
                         std::strerror(*errorNumber));
            const char * const *funcName
                = boost::get_error_info<boost::errinfo_api_function>(*boostEx);
            if (funcName != NULL)
            {
                std::fprintf(fp, "%*s(call to %s() returned error code %d)\n",
                             (indent+1)*2, "", *funcName, *errorNumber);
            }
        }

        // TODO: Treat also boost::nested_exception (not currently used, though)

        const internal::NestedExceptionList *nested
            = boost::get_error_info<errinfo_nested_exceptions>(*boostEx);
        if (nested != NULL)
        {
            internal::NestedExceptionList::const_iterator ni;
            for (ni = nested->begin(); ni != nested->end(); ++ni)
            {
                try
                {
                    rethrow_exception(*ni);
                }
                catch (const std::exception &nestedEx)
                {
                    printExceptionMessage(fp, nestedEx, indent + 1);
                }
            }
        }
    }
    else
    {
        internal::printFatalErrorMessageLine(fp, ex.what(), indent);
    }
}

}   // namespace

void printFatalErrorMessage(FILE *fp, const std::exception &ex)
{
    const char             *title      = "Unknown exception";
    bool                    bPrintType = false;
    const GromacsException *gmxEx      = dynamic_cast<const GromacsException *>(&ex);
    // TODO: Treat more of the standard exceptions
    if (gmxEx != NULL)
    {
        title = getErrorCodeString(gmxEx->errorCode());
    }
    else if (dynamic_cast<const tMPI::system_error *>(&ex) != NULL)
    {
        title = "System error in thread synchronization";
    }
    else if (dynamic_cast<const std::bad_alloc *>(&ex) != NULL)
    {
        title = "Memory allocation failed";
    }
    else if (dynamic_cast<const std::logic_error *>(&ex) != NULL)
    {
        title      = "Standard library logic error (bug)";
        bPrintType = true;
    }
    else if (dynamic_cast<const std::runtime_error *>(&ex) != NULL)
    {
        title      = "Standard library runtime error (possible bug)";
        bPrintType = true;
    }
    else
    {
        bPrintType = true;
    }
    // We can't call get_error_info directly on ex since our internal boost
    // needs to be compiled with BOOST_NO_RTTI. So we do the dynamic_cast
    // here instead.
    const char *const      *funcPtr = NULL;
    const char *const      *filePtr = NULL;
    const int              *linePtr = NULL;
    const boost::exception *boostEx = dynamic_cast<const boost::exception *>(&ex);
    if (boostEx != NULL)
    {
        funcPtr = boost::get_error_info<boost::throw_function>(*boostEx);
        filePtr = boost::get_error_info<boost::throw_file>(*boostEx);
        linePtr = boost::get_error_info<boost::throw_line>(*boostEx);
    }
    internal::printFatalErrorHeader(fp, title,
                                    funcPtr != NULL ? *funcPtr : NULL,
                                    filePtr != NULL ? *filePtr : NULL,
                                    linePtr != NULL ? *linePtr : 0);
    if (bPrintType)
    {
        std::fprintf(fp, "(exception type: %s)\n", typeid(ex).name());
    }
    printExceptionMessage(fp, ex, 0);
    internal::printFatalErrorFooter(fp);
}

} // namespace gmx
