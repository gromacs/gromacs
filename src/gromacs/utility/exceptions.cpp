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
#include "gromacs/utility/exceptions.h"

#include <boost/exception/get_error_info.hpp>

#include "gromacs/utility/errorcodes.h"

#include "errorformat.h"

namespace gmx
{

/********************************************************************
 * GromacsException
 */

GromacsException::GromacsException(const std::string &reason)
{
    *this << errinfo_message(reason);
}

const char *GromacsException::what() const throw()
{
    const std::string *msg = boost::get_error_info<errinfo_message>(*this);
    return msg != NULL ? msg->c_str() : "No reason provided";
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
    internal::printFatalErrorMessageLine(fp, ex.what(), 0);
    // TODO: Treat errno information in boost exceptions
    internal::printFatalErrorFooter(fp);
}

} // namespace gmx
