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

std::string formatErrorMessage(const std::exception &ex)
{
    const char *title = "Unknown exception";
    const char *func = NULL;
    const char *file = NULL;
    int line = 0;
    const GromacsException *gmxEx = dynamic_cast<const GromacsException *>(&ex);
    // TODO: Also treat common standard exceptions
    if (gmxEx != NULL)
    {
        title = getErrorCodeString(gmxEx->errorCode());
        func = *boost::get_error_info<boost::throw_function>(*gmxEx);
        file = *boost::get_error_info<boost::throw_file>(*gmxEx);
        line = *boost::get_error_info<boost::throw_line>(*gmxEx);
    }
    // TODO: Treat errno information in boost exceptions
    return internal::formatFatalError(title, ex.what(), func, file, line);
}

} // namespace gmx
