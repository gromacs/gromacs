/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements classes and functions in exceptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "exceptions.h"

#include "config.h"

#include <cstring>

#include <new>
#include <stdexcept>
#include <typeinfo>

#include <boost/shared_ptr.hpp>
#include <boost/exception/get_error_info.hpp>

#include "thread_mpi/system_error.h"

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "errorformat.h"

namespace gmx
{

namespace
{

/********************************************************************
 * ErrorMessage
 */

/*! \brief
 * Error message or error context text item.
 *
 * Error messages for an exception are represented as a chain of ErrorMessage
 * objects: the elements at the bottom of the chain (with no children) is the
 * error message, and other elements are the context strings added.
 *
 * \ingroup module_utility
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
        bool isContext() const { return static_cast<bool>(child_); }
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
    if (msg == NULL)
    {
        return "No reason provided";
    }
    while (msg->isContext())
    {
        msg = &msg->child();
    }
    return msg->text().c_str();
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

//! \addtogroup module_utility
//! \{

/*! \brief
 * Abstracts actual output from the other logic in exception formatting.
 *
 * Object that implements this interface is passed to
 * formatExceptionMessageInternal(), and is responsible for composing the
 * output.  This allows using the same implementation of interpreting the
 * exceptions while still supporting output to different formats (e.g., to a
 * string or to \c stderr).
 */
class MessageWriterInterface
{
    public:
        virtual ~MessageWriterInterface() {}

        /*! \brief
         * Writes a single line of text into the output.
         *
         * \param[in] text    Text to write on the line.
         * \param[in] indent  Suggested number of spaces to indent the line.
         */
        virtual void writeLine(const char *text, int indent) = 0;
        /*! \brief
         * Writes information about a system error (errno-based).
         *
         * \param[in] errorNumber  errno value
         * \param[in] funcName     Name of the system call (can be NULL).
         * \param[in] indent       Suggested number of spaces to indent the output.
         */
        virtual void writeErrNoInfo(int errorNumber, const char *funcName,
                                    int indent) = 0;
};

/*! \brief
 * Exception information writer for cases where exceptions should be avoided.
 *
 * Formats the messages into the provided FILE handle without checking for
 * errors in std::fprintf() calls.
 */
class MessageWriterFileNoThrow : public MessageWriterInterface
{
    public:
        //! Initializes a writer that writes to the given file handle.
        explicit MessageWriterFileNoThrow(FILE *fp) : fp_(fp) {}

        virtual void writeLine(const char *text, int indent)
        {
            internal::printFatalErrorMessageLine(fp_, text, indent);
        }
        virtual void writeErrNoInfo(int errorNumber, const char *funcName,
                                    int indent)
        {
            std::fprintf(fp_, "%*sReason: %s\n", indent, "",
                         std::strerror(errorNumber));
            if (funcName != NULL)
            {
                std::fprintf(fp_, "%*s(call to %s() returned error code %d)\n",
                             indent, "", funcName, errorNumber);
            }
        }

    private:
        FILE                   *fp_;
};

/*! \brief
 * Exception information writer to format into an std::string.
 */
class MessageWriterString : public MessageWriterInterface
{
    public:
        //! Post-processes the output string to not end in a line feed.
        void removeTerminatingLineFeed()
        {
            if (result_.size() > 0U)
            {
                result_.erase(result_.size() - 1);
            }
        }
        //! Returns the constructed string.
        const std::string &result() const { return result_; }

        virtual void writeLine(const char *text, int indent)
        {
            result_.append(indent, ' ');
            result_.append(text);
            result_.append("\n");
        }
        virtual void writeErrNoInfo(int errorNumber, const char *funcName,
                                    int indent)
        {
            writeLine(formatString("Reason: %s", std::strerror(errorNumber)).c_str(),
                      indent);
            if (funcName != NULL)
            {
                writeLine(formatString("(call to %s() returned error code %d)",
                                       funcName, errorNumber).c_str(),
                          indent);
            }
        }

    private:
        std::string             result_;
};

/*! \brief
 * Prints error information for an exception object.
 *
 * \param[in] writer  Writer to write out the information.
 * \param[in] ex      Exception object to print.
 * \param[in] indent  Indentation for the information.
 *
 * If the exception contains nested exceptions, information from them is
 * recursively printed.
 *
 * Does not throw unless the writer throws.
 */
void formatExceptionMessageInternal(MessageWriterInterface *writer,
                                    const std::exception &ex, int indent)
{
    const boost::exception *boostEx = dynamic_cast<const boost::exception *>(&ex);
    if (boostEx != NULL)
    {
        // TODO: Add an option to print this information for the tests
        // const char *const *funcPtr =
        //     boost::get_error_info<boost::throw_function>(*boostEx);
        // const char *const *filePtr =
        //     boost::get_error_info<boost::throw_file>(*boostEx);
        // const int         *linePtr =
        //     boost::get_error_info<boost::throw_line>(*boostEx);

        // std::string        result;
        // if (filePtr != NULL && linePtr != NULL)
        // {
        //     result = formatString("%s:%d: %s\n", *filePtr, *linePtr,
        //                           funcPtr != NULL ? *funcPtr : "");
        // }

        bool                bAnythingWritten = false;
        // TODO: Remove duplicate context if present in multiple nested exceptions.
        const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*boostEx);
        if (msg != NULL)
        {
            while (msg != NULL && msg->isContext())
            {
                writer->writeLine(msg->text().c_str(), indent*2);
                ++indent;
                msg = &msg->child();
            }
            if (msg != NULL && !msg->text().empty())
            {
                writer->writeLine(msg->text().c_str(), indent*2);
                bAnythingWritten = true;
            }
        }
        else
        {
            writer->writeLine(ex.what(), indent*2);
            bAnythingWritten = true;
        }

        const int *errorNumber
            = boost::get_error_info<boost::errinfo_errno>(*boostEx);
        if (errorNumber != NULL)
        {
            const char * const *funcName
                = boost::get_error_info<boost::errinfo_api_function>(*boostEx);
            writer->writeErrNoInfo(*errorNumber,
                                   funcName != NULL ? *funcName : NULL,
                                   (indent+1)*2);
            bAnythingWritten = true;
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
                    const int newIndent = indent + (bAnythingWritten ? 1 : 0);
                    formatExceptionMessageInternal(writer, nestedEx, newIndent);
                }
            }
        }
    }
    else
    {
        writer->writeLine(ex.what(), indent*2);
    }
}

//! \}

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
    MessageWriterFileNoThrow writer(fp);
    formatExceptionMessageInternal(&writer, ex, 0);
    internal::printFatalErrorFooter(fp);
}

std::string formatExceptionMessageToString(const std::exception &ex)
{
    MessageWriterString writer;
    formatExceptionMessageInternal(&writer, ex, 0);
    writer.removeTerminatingLineFeed();
    return writer.result();
}

void formatExceptionMessageToFile(FILE *fp, const std::exception &ex)
{
    MessageWriterFileNoThrow writer(fp);
    formatExceptionMessageInternal(&writer, ex, 0);
}

int processExceptionAtExit(const std::exception & /*ex*/)
{
    int returnCode = 1;
#ifdef GMX_LIB_MPI
    // TODO: Consider moving the output done in gmx_abort() into the message
    // printing routine above, so that this could become a simple MPI_Abort().
    gmx_abort(returnCode);
#endif
    return returnCode;
}

} // namespace gmx
