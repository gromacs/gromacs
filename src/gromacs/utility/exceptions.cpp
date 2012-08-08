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

#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/messagestringcollector.h"

#include "errorformat.h"

namespace gmx
{

namespace
{

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

} // namespace

/********************************************************************
 * ErrorMessage
 */

ErrorMessage::ErrorMessage(const char *reason)
    : impl_(new Impl)
{
    impl_->reason_ = reason;
}

ErrorMessage::ErrorMessage(const std::string &reason)
    : impl_(new Impl)
{
    impl_->reason_ = reason;
}

ErrorMessage::ErrorMessage(const std::string &context, const std::string &reason)
    : impl_(new Impl)
{
    impl_->context_ = context;
    impl_->reason_ = reason;
}

ErrorMessage::~ErrorMessage()
{
}

ErrorMessage
ErrorMessage::prependContext(const std::string &context) const
{
    return prependContext(context, std::string());
}

ErrorMessage
ErrorMessage::prependContext(const std::string &context,
                             const std::string &reason) const
{
    ErrorMessage newMessage(context, reason);
    newMessage.addDetails(*this);
    return newMessage;
}

/********************************************************************
 * GromacsException
 */

GromacsException::GromacsException(const ErrorMessage &reason)
{
    *this << errinfo_message(reason);
}

const char *GromacsException::what() const throw()
{
    const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*this);
    while (msg != NULL && msg->reason().empty())
    {
        msg = msg->hasDetails() ? &msg->details().front() : NULL;
    }
    return msg != NULL ? msg->reason().c_str() : "No reason provided";
}

const ErrorMessage &GromacsException::message() const
{
    const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*this);
    GMX_RELEASE_ASSERT(msg != NULL, "Message should always be set");
    return *msg;
}

void GromacsException::prependContext(const std::string &context)
{
    const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*this);
    GMX_RELEASE_ASSERT(msg != NULL, "Message should always be set");
    *this << errinfo_message(msg->prependContext(context));
}

void GromacsException::prependContext(const std::string &context,
                                      const std::string &reason)
{
    const ErrorMessage *msg = boost::get_error_info<errinfo_message>(*this);
    GMX_RELEASE_ASSERT(msg != NULL, "Message should always be set");
    *this << errinfo_message(msg->prependContext(context, reason));
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

void printGromacsExceptionMessage(FILE *fp, const ErrorMessage &msg, int indent)
{
    // TODO: Line wrapping in the whole function.
    // TODO: Remove duplicate context if present in multiple details objects.
    if (!msg.context().empty())
    {
        std::fprintf(fp, "%*s%s\n", indent*2, "", msg.context().c_str());
        ++indent;
    }
    if (!msg.reason().empty())
    {
        std::fprintf(fp, "%*s%s\n", indent*2, "", msg.reason().c_str());
    }

    const std::vector<ErrorMessage> &details = msg.details();
    std::vector<ErrorMessage>::const_iterator di;
    for (di = details.begin(); di != details.end(); ++di)
    {
        printGromacsExceptionMessage(fp, *di, indent);
    }
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
    if (gmxEx != NULL)
    {
        printGromacsExceptionMessage(fp, gmxEx->message(), 0);
    }
    else
    {
        // TODO: Line wrapping
        std::fprintf(fp, "%s\n", ex.what());
    }
    // TODO: Treat errno information in boost exceptions
    internal::printFatalErrorFooter(fp);
}

} // namespace gmx
