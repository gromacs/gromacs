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
/*! \internal \file
 * \brief
 * Implements classes and functions in exceptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "exceptions.h"

#include <boost/exception/get_error_info.hpp>
#include <boost/shared_ptr.hpp>

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
        bool isContext() const { return child().impl_; }
        const std::string &text() const;
        const ErrorMessage &child() const;

        ErrorMessage prependContext(const std::string &context) const;

    private:
        ErrorMessage() {}
        explicit ErrorMessage(const std::string &reason);

        struct Impl;

        boost::shared_ptr<Impl> impl_;

        friend class gmx::GromacsException;
};

/*! \internal \brief
 * Stores a user-friendly explanation for the reason of an exception.
 *
 * Typically, should not be used directly, but through the GromacsException
 * class: it is initialized by the constructor, and can be accessed with
 * GromacsException::what().
 *
 * \ingroup module_utility
 */
typedef boost::error_info<struct errinfo_message_, ErrorMessage>
        errinfo_message;

struct ErrorMessage::Impl
{
    std::string         text_;
    ErrorMessage        child_;
};

ErrorMessage::ErrorMessage(const std::string &reason)
    : impl_(new Impl)
{
    size_t length = reason.find_last_not_of(" \n");
    if (length == std::string::npos)
    {
        length = reason.length() - 1;
    }
    impl_->text_ = reason.substr(0, length + 1);
}

const std::string &ErrorMessage::text() const
{
    return impl_->text_;
}

const ErrorMessage &ErrorMessage::child() const
{
    return impl_->child_;
}

ErrorMessage
ErrorMessage::prependContext(const std::string &context) const
{
    ErrorMessage newMessage(context);
    newMessage.impl_->child_ = *this;
    return newMessage;
}

typedef boost::error_info<struct errinfo_message_, internal::NestedExceptionList>
        errinfo_nested_exceptions;

} // namespace

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
            internal::printFatalErrorMessageLine(fp, ex.what(), 0);
        }

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
                    printExceptionMessage(fp, nestedEx, indent);
                }
            }
        }
    }
    else
    {
        internal::printFatalErrorMessageLine(fp, ex.what(), 0);
    }
    // TODO: Treat errno information in boost exceptions
}

} // namespace

void printFatalErrorMessage(FILE *fp, const std::exception &ex)
{
    const char *title = "Unknown exception";
    const GromacsException *gmxEx = dynamic_cast<const GromacsException *>(&ex);
    // TODO: Also treat common standard exceptions
    if (gmxEx != NULL)
    {
        title = getErrorCodeString(gmxEx->errorCode());
    }
    else if (dynamic_cast<const std::bad_alloc *>(&ex) != NULL)
    {
        title = "Memory allocation failed";
    }
    // We can't call get_error_info directly on ex since our internal boost
    // needs to be compiled with BOOST_NO_RTTI. So we do the dynamic_cast
    // here instead.
    const char *const *funcPtr = NULL;
    const char *const *filePtr = NULL;
    const int         *linePtr = NULL;
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
    printExceptionMessage(fp, ex, 0);
    internal::printFatalErrorFooter(fp);
}

} // namespace gmx
