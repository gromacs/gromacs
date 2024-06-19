/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements classes and functions in exceptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/exceptions.h"

#include <cstdio>
#include <cstring>

#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <utility>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "errorcodes.h"
#include "errorformat.h"

namespace gmx
{

namespace internal
{

IExceptionInfo::~IExceptionInfo() {}

class ExceptionData
{
public:
    std::map<std::type_index, ExceptionInfoPointer> infos_;
};

} // namespace internal

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
    explicit ErrorMessage(const std::string& text);

    //! Whether this object is a context string.
    bool isContext() const { return static_cast<bool>(child_); }
    //! Returns the text for this object.
    const std::string& text() const { return text_; }
    /*! \brief
     * Returns the child object for a context object.
     *
     * Must not be called if isContext() returns false.
     */
    const ErrorMessage& child() const
    {
        GMX_ASSERT(isContext(), "Attempting to access nonexistent message object");
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
    ErrorMessage prependContext(const std::string& context) const;

private:
    std::string                   text_;
    std::shared_ptr<ErrorMessage> child_;
};

/*! \internal \brief
 * Stores a reason or the top-most context string of an exception.
 *
 * \ingroup module_utility
 */
typedef ExceptionInfo<struct ExceptionInfoMessage_, ErrorMessage> ExceptionInfoMessage;

ErrorMessage::ErrorMessage(const std::string& text) : text_(text)
{
    size_t length = text_.find_last_not_of(" \n");
    if (length == std::string::npos)
    {
        length = text_.length() - 1;
    }
    text_.resize(length + 1);
}

ErrorMessage ErrorMessage::prependContext(const std::string& context) const
{
    ErrorMessage newMessage(context);
    newMessage.child_ = std::make_shared<ErrorMessage>(*this);
    return newMessage;
}

/*! \brief
 * Stores list of nested exceptions for Gromacs exceptions.
 *
 * \ingroup module_utility
 */
typedef ExceptionInfo<struct ExceptionInfoNestedExceptions_, internal::NestedExceptionList> ExceptionInfoNestedExceptions;

} // namespace

/********************************************************************
 * GromacsException
 */

GromacsException::GromacsException(const ExceptionInitializer& details) :
    data_(new internal::ExceptionData)
{
    setInfo(ExceptionInfoMessage(ErrorMessage(details.reason_)));
    if (details.hasNestedExceptions())
    {
        setInfo(ExceptionInfoNestedExceptions(details.nested_));
    }
}

const char* GromacsException::what() const noexcept
{
    const ErrorMessage* msg = getInfo<ExceptionInfoMessage>();
    if (msg == nullptr)
    {
        return "No reason provided";
    }
    while (msg->isContext())
    {
        msg = &msg->child();
    }
    return msg->text().c_str();
}

void GromacsException::prependContext(const std::string& context)
{
    const ErrorMessage* msg = getInfo<ExceptionInfoMessage>();
    GMX_RELEASE_ASSERT(msg != nullptr, "Message should always be set");
    setInfo(ExceptionInfoMessage(msg->prependContext(context)));
}

const internal::IExceptionInfo* GromacsException::getInfo(const std::type_index& index) const
{
    auto iter = data_->infos_.find(index);
    if (iter != data_->infos_.end())
    {
        return iter->second.get();
    }
    return nullptr;
}

void GromacsException::setInfo(const std::type_index& index, internal::ExceptionInfoPointer&& item)
{
    data_->infos_[index] = std::move(item);
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

int ToleranceError::errorCode() const
{
    return eeTolerance;
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

int RangeError::errorCode() const
{
    return eeRange;
}

int NotImplementedError::errorCode() const
{
    return eeNotImplemented;
}

int ParallelConsistencyError::errorCode() const
{
    return eeParallelConsistency;
}

int ModularSimulatorError::errorCode() const
{
    return eeModularSimulator;
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
class IMessageWriter
{
public:
    virtual ~IMessageWriter() {}

    /*! \brief
     * Writes a single line of text into the output.
     *
     * \param[in] text    Text to write on the line.
     * \param[in] indent  Suggested number of spaces to indent the line.
     */
    virtual void writeLine(const char* text, int indent) = 0;
    /*! \brief
     * Writes information about a system error (errno-based).
     *
     * \param[in] errorNumber  errno value
     * \param[in] funcName     Name of the system call (can be NULL).
     * \param[in] indent       Suggested number of spaces to indent the output.
     */
    virtual void writeErrNoInfo(int errorNumber, const char* funcName, int indent) = 0;
};

/*! \brief
 * Exception information writer for cases where exceptions should be avoided.
 *
 * Formats the messages into the provided FILE handle without checking for
 * errors in std::fprintf() calls.
 */
class MessageWriterFileNoThrow : public IMessageWriter
{
public:
    //! Initializes a writer that writes to the given file handle.
    explicit MessageWriterFileNoThrow(FILE* fp) : fp_(fp) {}

    void writeLine(const char* text, int indent) override
    {
        internal::printFatalErrorMessageLine(fp_, text, indent);
    }
    void writeErrNoInfo(int errorNumber, const char* funcName, int indent) override
    {
        std::fprintf(fp_, "%*sReason: %s\n", indent, "", std::strerror(errorNumber));
        if (funcName != nullptr)
        {
            std::fprintf(fp_, "%*s(call to %s() returned error code %d)\n", indent, "", funcName, errorNumber);
        }
    }

private:
    FILE* fp_;
};

/*! \brief
 * Exception information writer to format into a TextOutputStream.
 */
class MessageWriterTextWriter : public IMessageWriter
{
public:
    //! Initializes a writer that writes to the given stream.
    explicit MessageWriterTextWriter(TextWriter* writer) : writer_(writer) {}

    void writeLine(const char* text, int indent) override
    {
        writer_->wrapperSettings().setIndent(indent);
        writer_->writeLine(text);
    }
    void writeErrNoInfo(int errorNumber, const char* funcName, int indent) override
    {
        writer_->wrapperSettings().setIndent(indent);
        writer_->writeLine(formatString("Reason: %s", std::strerror(errorNumber)));
        if (funcName != nullptr)
        {
            writer_->writeLine(formatString("(call to %s() returned error code %d)", funcName, errorNumber));
        }
    }

private:
    TextWriter* writer_;
};

/*! \brief
 * Exception information writer to format into an std::string.
 */
class MessageWriterString : public IMessageWriter
{
public:
    //! Post-processes the output string to not end in a line feed.
    void removeTerminatingLineFeed()
    {
        if (!result_.empty())
        {
            result_.erase(result_.size() - 1);
        }
    }
    //! Returns the constructed string.
    const std::string& result() const { return result_; }

    void writeLine(const char* text, int indent) override
    {
        result_.append(indent, ' ');
        result_.append(text);
        result_.append("\n");
    }
    void writeErrNoInfo(int errorNumber, const char* funcName, int indent) override
    {
        writeLine(formatString("Reason: %s", std::strerror(errorNumber)).c_str(), indent);
        if (funcName != nullptr)
        {
            writeLine(formatString("(call to %s() returned error code %d)", funcName, errorNumber).c_str(),
                      indent);
        }
    }

private:
    std::string result_;
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
void formatExceptionMessageInternal(IMessageWriter* writer, const std::exception& ex, int indent)
{
    const GromacsException* gmxEx = dynamic_cast<const GromacsException*>(&ex);
    if (gmxEx != nullptr)
    {
        // TODO: Add an option to print location information for the tests

        // std::string        result;
        // if (filePtr != NULL && linePtr != NULL)
        // {
        //     result = formatString("%s:%d: %s\n", *filePtr, *linePtr,
        //                           funcPtr != NULL ? *funcPtr : "");
        // }

        bool bAnythingWritten = false;
        // TODO: Remove duplicate context if present in multiple nested exceptions.
        const ErrorMessage* msg = gmxEx->getInfo<ExceptionInfoMessage>();
        if (msg != nullptr)
        {
            while (msg != nullptr && msg->isContext())
            {
                writer->writeLine(msg->text().c_str(), indent * 2);
                ++indent;
                msg = &msg->child();
            }
            if (msg != nullptr && !msg->text().empty())
            {
                writer->writeLine(msg->text().c_str(), indent * 2);
                bAnythingWritten = true;
            }
        }
        else
        {
            writer->writeLine(ex.what(), indent * 2);
            bAnythingWritten = true;
        }

        const int* errorNumber = gmxEx->getInfo<ExceptionInfoErrno>();
        if (errorNumber != nullptr && *errorNumber != 0)
        {
            const char* const* funcName = gmxEx->getInfo<ExceptionInfoApiFunction>();
            writer->writeErrNoInfo(
                    *errorNumber, funcName != nullptr ? *funcName : nullptr, (indent + 1) * 2);
            bAnythingWritten = true;
        }

        const internal::NestedExceptionList* nested = gmxEx->getInfo<ExceptionInfoNestedExceptions>();
        if (nested != nullptr)
        {
            internal::NestedExceptionList::const_iterator ni;
            for (ni = nested->begin(); ni != nested->end(); ++ni)
            {
                try
                {
                    std::rethrow_exception(*ni);
                }
                catch (const std::exception& nestedEx)
                {
                    const int newIndent = indent + (bAnythingWritten ? 1 : 0);
                    formatExceptionMessageInternal(writer, nestedEx, newIndent);
                }
            }
        }
    }
    else
    {
        writer->writeLine(ex.what(), indent * 2);
    }
}

//! \}

} // namespace

void printFatalErrorMessage(FILE* fp, const std::exception& ex)
{
    const char*             title      = "Unknown exception";
    bool                    bPrintType = false;
    const GromacsException* gmxEx      = dynamic_cast<const GromacsException*>(&ex);
    // TODO: Treat more of the standard exceptions
    if (gmxEx != nullptr)
    {
        title = getErrorCodeString(gmxEx->errorCode());
    }
    else if (dynamic_cast<const std::bad_alloc*>(&ex) != nullptr)
    {
        title = "Memory allocation failed";
    }
    else if (dynamic_cast<const std::logic_error*>(&ex) != nullptr)
    {
        title      = "Standard library logic error (bug)";
        bPrintType = true;
    }
    else if (dynamic_cast<const std::runtime_error*>(&ex) != nullptr)
    {
        title      = "Standard library runtime error (possible bug)";
        bPrintType = true;
    }
    else
    {
        bPrintType = true;
    }
    const char* func = nullptr;
    const char* file = nullptr;
    int         line = 0;
    if (gmxEx != nullptr)
    {
        const ThrowLocation* loc = gmxEx->getInfo<ExceptionInfoLocation>();
        if (loc != nullptr)
        {
            func = loc->func;
            file = loc->file;
            line = loc->line;
        }
    }
    internal::printFatalErrorHeader(fp, title, func, file, line);
    if (bPrintType)
    {
        std::fprintf(fp, "(exception type: %s)\n", typeid(ex).name());
    }
    MessageWriterFileNoThrow writer(fp);
    formatExceptionMessageInternal(&writer, ex, 0);
    internal::printFatalErrorFooter(fp);
}

std::string formatExceptionMessageToString(const std::exception& ex)
{
    MessageWriterString writer;
    formatExceptionMessageInternal(&writer, ex, 0);
    writer.removeTerminatingLineFeed();
    return writer.result();
}

void formatExceptionMessageToFile(FILE* fp, const std::exception& ex)
{
    MessageWriterFileNoThrow writer(fp);
    formatExceptionMessageInternal(&writer, ex, 0);
}

void formatExceptionMessageToWriter(TextWriter* writer, const std::exception& ex)
{
    MessageWriterTextWriter messageWriter(writer);
    formatExceptionMessageInternal(&messageWriter, ex, 0);
}

int processExceptionAtExit(const std::exception& /*ex*/)
{
    int returnCode = 1;
    // If we have more than one rank (whether real MPI or thread-MPI),
    // we cannot currently know whether just one rank or all ranks
    // actually threw the error, so we need to exit here.
    // Returning would mean graceful cleanup, which is not possible if
    // some code is still executing on other ranks/threads.
    if (gmx_node_num() > 1)
    {
        gmx_exit_on_fatal_error(ExitType_Abort, returnCode);
    }
    return returnCode;
}

void processExceptionAsFatalError(const std::exception& ex)
{
    printFatalErrorMessage(stderr, ex);
    gmx_exit_on_fatal_error(ExitType_Abort, 1);
}

} // namespace gmx
