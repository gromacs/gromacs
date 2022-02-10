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
/*! \file
 * \brief
 * Declares common exception classes and macros for fatal error handling.
 *
 * The basic approach is the same as in boost::exception for storing additional
 * context information to exceptions, but since that functionality is a very
 * small and simple part of boost::exception, the code is duplicated here.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_EXCEPTIONS_H
#define GMX_UTILITY_EXCEPTIONS_H

#include <cstdio>
#include <cstdlib>

#include <exception>
#include <memory>
#include <string>
#include <type_traits>
#include <typeindex>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

class TextWriter;

namespace internal
{
//! Internal container type for storing a list of nested exceptions.
typedef std::vector<std::exception_ptr> NestedExceptionList;

/*! \internal
 * \brief
 * Base class for ExceptionInfo.
 *
 * This class only provides a way to store different ExceptionInfo objects in
 * the same container.  Actual access to the ExceptionInfo items is handled by
 * downcasting, after looking up the correct item based on its type.
 *
 * \ingroup module_utility
 */
class IExceptionInfo
{
public:
    virtual ~IExceptionInfo();
    IExceptionInfo()                          = default;
    IExceptionInfo(const IExceptionInfo&)     = default;
    IExceptionInfo(IExceptionInfo&&) noexcept = default;
    IExceptionInfo& operator=(const IExceptionInfo&) = default;
    IExceptionInfo& operator=(IExceptionInfo&&) noexcept = default;
};

//! Smart pointer to manage IExceptionInfo ownership.
typedef std::unique_ptr<IExceptionInfo> ExceptionInfoPointer;

class ExceptionData;

} // namespace internal

//! \addtogroup module_utility
//! \{

/*! \brief
 * Stores additional context information for exceptions.
 *
 * \tparam  Tag  Tag type (typically, a forward-declared struct that is not
 *     defined anywhere) that makes all ExceptionInfo types unique, even if
 *     they have the same value type.
 * \tparam  T    Type of value this object stores.
 *     Needs to be copy-constructible.
 *
 * Example of declaring a new info type that stores an integer:
 * \code
   typedef ExceptionInfo<struct ExceptionInfoMyInfo_, int> ExceptionInfoMyInfo;
   \endcode
 *
 * \inpublicapi
 */
template<class Tag, typename T>
class ExceptionInfo : public internal::IExceptionInfo
{
public:
    //! The type of value stored in this object.
    typedef T value_type;

    //! Creates an info object from given value.
    explicit ExceptionInfo(const T& value) : value_(value) {}

    //! Returns the stored value.
    const T& value() const { return value_; }

private:
    T value_;
};

/*! \internal
 * \brief
 * Stores the location from which an exception was thrown.
 */
struct ThrowLocation
{
    //! Creates an object for storing the throw location.
    ThrowLocation(const char* func, const char* file, int line) : func(func), file(file), line(line)
    {
    }

    //! Function where the throw occurred.
    const char* func;
    //! File where the throw occurred.
    const char* file;
    //! Line number where the throw occurred.
    int line;
};

//! Stores `errno` value that triggered the exception.
typedef ExceptionInfo<struct ExceptionInfoErrno_, int> ExceptionInfoErrno;
//! Stores the function name that returned the `errno` in ExceptionInfoErrno.
typedef ExceptionInfo<struct ExceptionInfoApiFunc_, const char*> ExceptionInfoApiFunction;
//! Stores the location where the exception was thrown.
typedef ExceptionInfo<struct ExceptionInfoLocation_, ThrowLocation> ExceptionInfoLocation;

/*! \brief
 * Provides information for Gromacs exception constructors.
 *
 * This class exists to implement common functionality for initializing all
 * Gromacs exceptions without having extra code in each exception class.
 * In simple cases, it can be implicitly constructed by passing a simple string
 * to an exception constructor.
 * If more complex initialization is necessary, it is possible to explicitly
 * construct an object of this type and then call other methods to add
 * information before actually creating the exception object.
 *
 * \todo
 * With the exception of the reason string, information added with this class
 * is not currently accessible through any public API, except for calling
 * printFatalErrorMessage(), formatExceptionMessageToString() or
 * formatExceptionMessageToFile().  This is not implemented as there is not yet
 * need for it, and it is not clear what would be the best alternative for the
 * access.  It should be possible to refactor the internal implementation to
 * suit the needs of such external access without requiring changes in code
 * that throws these exceptions.
 *
 * \ingroup module_utility
 */
class ExceptionInitializer
{
public:
    /*! \brief
     * Creates an initialized with the given string as the reason.
     *
     * \param[in] reason  Detailed reason for the exception.
     * \throw     std::bad_alloc if out of memory.
     *
     * This constructor is not explicit to allow constructing exceptions
     * with a plain string argument given to the constructor without adding
     * extra code to each exception class.
     */
    ExceptionInitializer(const char* reason) : reason_(reason) {}
    //! \copydoc ExceptionInitializer(const char *)
    ExceptionInitializer(const std::string& reason) : reason_(reason) {}

    /*! \brief
     * Returns true if addCurrentExceptionAsNested() has been called.
     *
     * Provided for convenience for cases where exceptions will be added
     * conditionally, and the caller wants to check whether any excetions
     * were actually added.
     */
    bool hasNestedExceptions() const { return !nested_.empty(); }
    /*! \brief
     * Adds the currently caught exception as a nested exception.
     *
     * May be called multiple times; all provided exceptions will be added
     * in a list of nested exceptions.
     *
     * Must not be called outside a catch block.
     */
    void addCurrentExceptionAsNested() { nested_.push_back(std::current_exception()); }
    /*! \brief
     * Adds the specified exception as a nested exception.
     *
     * May be called multiple times; all provided exceptions will be added
     * in a list of nested exceptions.
     *
     * This is equivalent to throwing \p ex and calling
     * addCurrentExceptionAsNested() in the catch block, but potentially
     * more efficient.
     */
    template<class Exception>
    void addNested(const Exception& ex)
    {
        nested_.push_back(std::make_exception_ptr(ex));
    }

private:
    std::string                   reason_;
    internal::NestedExceptionList nested_;

    friend class GromacsException;
};

/*! \brief
 * Base class for all exception objects in Gromacs.
 *
 * \inpublicapi
 */
class GromacsException : public std::exception
{
public:
    // Explicitly declared because some compiler/library combinations warn
    // about missing noexcept otherwise.
    ~GromacsException() noexcept override {}

    GromacsException()                            = default;
    GromacsException(const GromacsException&)     = default;
    GromacsException(GromacsException&&) noexcept = default;
    GromacsException& operator=(const GromacsException&) = default;
    GromacsException& operator=(GromacsException&&) noexcept = default;

    /*! \brief
     * Returns the reason string for the exception.
     *
     * The return value is the string that was passed to the constructor.
     */
    const char* what() const noexcept override;
    /*! \brief
     * Returns the error code corresponding to the exception type.
     */
    virtual int errorCode() const = 0;

    /*! \brief
     * Returns the value associated with given ExceptionInfo.
     *
     * \tparam  InfoType  ExceptionInfo type to get the value for.
     * \returns Value set for `InfoType`, or `nullptr` if such info has not
     *     been set.
     *
     * Does not throw.
     */
    template<class InfoType>
    const typename InfoType::value_type* getInfo() const
    {
        const internal::IExceptionInfo* item = getInfo(typeid(InfoType));
        if (item != nullptr)
        {
            GMX_ASSERT(dynamic_cast<const InfoType*>(item) != nullptr,
                       "Invalid exception info item found");
            return &static_cast<const InfoType*>(item)->value();
        }
        return nullptr;
    }

    /*! \brief
     * Associates extra information with the exception.
     *
     * \tparam  Tag  ExceptionInfo tag type.
     * \tparam  T          ExceptionInfo value type.
     * \param[in] item  ExceptionInfo to associate.
     * \throws std::bad_alloc if out of memory.
     * \throws unspecified    any exception thrown by `T` copy construction.
     *
     * If an item of this type is already associated, it is overwritten.
     */
    template<class Tag, typename T>
    void setInfo(const ExceptionInfo<Tag, T>& item)
    {
        typedef ExceptionInfo<Tag, T>  ItemType;
        internal::ExceptionInfoPointer itemPtr(new ItemType(item));
        setInfo(typeid(ItemType), std::move(itemPtr));
    }

    /*! \brief
     * Adds context information to this exception.
     *
     * \param[in] context  Context string to add.
     * \throws    std::bad_alloc if out of memory.
     *
     * Typical use is to add additional information higher up in the call
     * stack using this function in a catch block and the rethrow the
     * exception.
     *
     * \todo
     * The added information is currently not accessible through what(),
     * nor through any other means except for calling
     * printFatalErrorMessage(), formatExceptionMessageToString() or
     * formatExceptionMessageToFile(). See ExceptionInitializer for more
     * discussion.
     */
    void prependContext(const std::string& context);

protected:
    /*! \brief
     * Creates an exception object with the provided initializer/reason.
     *
     * \param[in] details  Initializer for the exception.
     * \throws    std::bad_alloc if out of memory.
     */
    explicit GromacsException(const ExceptionInitializer& details);

private:
    const internal::IExceptionInfo* getInfo(const std::type_index& index) const;
    void setInfo(const std::type_index& index, internal::ExceptionInfoPointer&& item);

    std::shared_ptr<internal::ExceptionData> data_;
};

/*! \brief
 * Associates extra information with an exception.
 *
 * \tparam  Exception  Exception type (must be derived from GromacsException).
 * \tparam  Tag        ExceptionInfo tag.
 * \tparam  T          ExceptionInfo value type.
 * \param[in,out] ex   Exception to associate the information to.
 * \param[in]     item Information to associate.
 *
 * \internal
 * The association is done with a templated non-member operator of exactly this
 * form to make the simple syntax of GMX_THROW() possible.  To support this,
 * this operation needs to:
 *  - Allow setting information in a temporary to support
 *    `GMX_THROW(InvalidInputError(ex))`.
 *  - Return a copy of the same class it takes in.  The compiler needs
 *    this information to throw the correct type of exception.  This
 *    would be tedious to achieve with a member function (without a
 *    lot of code duplication).  Generally, \c ex will be a temporary,
 *    copied twice and returned by value, which the compiler will
 *    typically elide away (and anyway performance is not important
 *    when throwing).  We are not using the typical
 *    return-by-const-reference idiom for this operator so that
 *    tooling can reliably see that we are throwing by value.
 *  - Provide convenient syntax for adding multiple items.  A non-member
 *    function that would require nested calls would look ugly for such cases.
 *
 * The reason for the enable_if is that this way, it does not conflict with
 * other overloads of `operator<<` for ExceptionInfo objects, in case someone
 * would like to declare those.  But currently we do not have such overloads, so
 * if the enable_if causes problems with some compilers, it can be removed.
 *
 * \todo Use std::is_base_of_v when CUDA 11 is a requirement.
 */
template<class Exception, class Tag, class T>
inline std::enable_if_t<std::is_base_of<GromacsException, Exception>::value, Exception>
operator<<(Exception ex, const ExceptionInfo<Tag, T>& item)
{
    ex.setInfo(item);
    return ex;
}

/*! \brief
 * Exception class for file I/O errors.
 *
 * \inpublicapi
 */
class FileIOError : public GromacsException
{
public:
    /*! \brief
     * Creates an exception object with the provided initializer/reason.
     *
     * \param[in] details  Initializer for the exception.
     * \throws    std::bad_alloc if out of memory.
     *
     * It is possible to call this constructor either with an explicit
     * ExceptionInitializer object (useful for more complex cases), or
     * a simple string if only a reason string needs to be provided.
     */
    explicit FileIOError(const ExceptionInitializer& details) : GromacsException(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for user input errors.
 *
 * Derived classes should be used to indicate the nature of the error instead
 * of throwing this class directly.
 *
 * \inpublicapi
 */
class UserInputError : public GromacsException
{
protected:
    //! \copydoc FileIOError::FileIOError()
    explicit UserInputError(const ExceptionInitializer& details) : GromacsException(details) {}
};

/*! \brief
 * Exception class for situations where user input cannot be parsed/understood.
 *
 * \inpublicapi
 */
class InvalidInputError : public UserInputError
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit InvalidInputError(const ExceptionInitializer& details) : UserInputError(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for situations where user input is inconsistent.
 *
 * \inpublicapi
 */
class InconsistentInputError : public UserInputError
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit InconsistentInputError(const ExceptionInitializer& details) : UserInputError(details)
    {
    }

    int errorCode() const override;
};

/*! \brief
 * Exception class when a specified tolerance cannot be achieved.
 *
 * \inpublicapi
 */
class ToleranceError : public GromacsException
{
public:
    /*! \brief
     * Creates an exception object with the provided initializer/reason.
     *
     * \param[in] details  Initializer for the exception.
     * \throws    std::bad_alloc if out of memory.
     *
     * It is possible to call this constructor either with an explicit
     * ExceptionInitializer object (useful for more complex cases), or
     * a simple string if only a reason string needs to be provided.
     */
    explicit ToleranceError(const ExceptionInitializer& details) : GromacsException(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for simulation instabilities.
 *
 * \inpublicapi
 */
class SimulationInstabilityError : public GromacsException
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit SimulationInstabilityError(const ExceptionInitializer& details) :
        GromacsException(details)
    {
    }

    int errorCode() const override;
};

/*! \brief
 * Exception class for internal errors.
 *
 * \inpublicapi
 */
class InternalError : public GromacsException
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit InternalError(const ExceptionInitializer& details) : GromacsException(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for incorrect use of an API.
 *
 * \inpublicapi
 */
class APIError : public GromacsException
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit APIError(const ExceptionInitializer& details) : GromacsException(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for out-of-range values or indices
 *
 * \inpublicapi
 */
class RangeError : public GromacsException
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit RangeError(const ExceptionInitializer& details) : GromacsException(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for use of an unimplemented feature.
 *
 * \inpublicapi
 */
class NotImplementedError : public APIError
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit NotImplementedError(const ExceptionInitializer& details) : APIError(details) {}

    int errorCode() const override;
};

/*! \brief Exception class for use when ensuring that MPI ranks to throw
 * in a coordinated fashion.
 *
 * Generally all ranks that can throw would need to check for whether
 * an exception has been caught, communicate whether any rank caught,
 * then all throw one of these, with either a string that describes
 * any exception caught on that rank, or a generic string.
 *
 * \inpublicapi
 */
class ParallelConsistencyError : public APIError
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit ParallelConsistencyError(const ExceptionInitializer& details) : APIError(details) {}

    int errorCode() const override;
};

/*! \brief
 * Exception class for modular simulator.
 *
 * \inpublicapi
 */
class ModularSimulatorError : public GromacsException
{
public:
    //! \copydoc FileIOError::FileIOError()
    explicit ModularSimulatorError(const ExceptionInitializer& details) : GromacsException(details)
    {
    }

    [[nodiscard]] int errorCode() const override;
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
   \endcode
 */
#define GMX_THROW(e) \
    throw(e) << gmx::ExceptionInfoLocation(gmx::ThrowLocation(GMX_CURRENT_FUNCTION, __FILE__, __LINE__))

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
   \endcode
 */
#define GMX_THROW_WITH_ERRNO(e, syscall, err)                     \
    do                                                            \
    {                                                             \
        int stored_errno_ = (err);                                \
        GMX_THROW((e) << gmx::ExceptionInfoErrno(stored_errno_)   \
                      << gmx::ExceptionInfoApiFunction(syscall)); \
    } while (0)
// TODO: Add an equivalent macro for Windows GetLastError

/*! \brief
 * Formats a standard fatal error message for reporting an exception.
 *
 * \param[in] fp  %File to format the message to.
 * \param[in] ex  Exception to format.
 *
 * Does not throw.  If memory allocation fails or some other error occurs
 * while formatting the error, tries to print a reasonable alternative message.
 *
 * Normal usage in Gromacs command-line programs is like this:
 * \code
   int main(int argc, char *argv[])
   {
       gmx::init(&argc, &argv);
       try
       {
           // The actual code for the program
           return 0;
       }
       catch (const std::exception &ex)
       {
           gmx::printFatalErrorMessage(stderr, ex);
           return gmx::processExceptionAtExit(ex);
       }
   }
   \endcode
 */
void printFatalErrorMessage(FILE* fp, const std::exception& ex);
/*! \brief
 * Formats an error message for reporting an exception.
 *
 * \param[in] ex  Exception to format.
 * \returns   Formatted string containing details of \p ex.
 * \throws    std::bad_alloc if out of memory.
 */
std::string formatExceptionMessageToString(const std::exception& ex);
/*! \brief
 * Formats an error message for reporting an exception.
 *
 * \param     fp  %File to write the message to.
 * \param[in] ex  Exception to format.
 * \throws    std::bad_alloc if out of memory.
 */
void formatExceptionMessageToFile(FILE* fp, const std::exception& ex);
/*! \brief
 * Formats an error message for reporting an exception.
 *
 * \param     writer  Writer to use for writing the message.
 * \param[in] ex      Exception to format.
 * \throws    std::bad_alloc if out of memory.
 */
void formatExceptionMessageToWriter(TextWriter* writer, const std::exception& ex);
/*! \brief
 * Handles an exception that is causing the program to terminate.
 *
 * \param[in] ex  Exception that is the cause for terminating the program.
 * \returns   Return code to return from main().
 *
 * This method should be called as the last thing before terminating the
 * program because of an exception.  It exists to terminate the program as
 * gracefully as possible in the case of MPI processing (but the current
 * implementation always calls MPI_Abort()).
 *
 * See printFatalErrorMessage() for example usage.
 *
 * Does not throw.
 */
int processExceptionAtExit(const std::exception& ex);

/*! \brief
 * Helper function for terminating the program on an exception.
 *
 * \param[in] ex  Exception that is the cause for terminating the program.
 *
 * Does not throw, and does not return.
 */
[[noreturn]] void processExceptionAsFatalError(const std::exception& ex);

/*! \brief
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
   \code
   try
   {
       // C++ code
   }
   GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
   \endcode
 */
#define GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR \
    catch (const std::exception& ex) { ::gmx::processExceptionAsFatalError(ex); }

//! \}

} // namespace gmx

#endif
