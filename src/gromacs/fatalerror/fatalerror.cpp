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
 * Implements functions in fatalerror.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_fatalerror
 */
#include "gromacs/fatalerror/fatalerror.h"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

static const char *const error_names[] =
{
    "No error",
    "Out of memory",
    "File not found",
    "System I/O error",
    "Error in user input",
    "Inconsistency in user input",
    "Simulation instability detected",

    "Feature not implemented",
    "Invalid value (bug)",
    "Invalid call (bug)",
    "Internal error (bug)",
    "Range checking error (possible bug)",
    "Communication error (possible bug)",

    "Unknown error",
};

namespace gmx
{

static void standardErrorHandler(int retcode, const char *msg,
                                 const char *file, int line)
{
    using std::fprintf;

    fprintf(stderr, "\n-------------------------------------------------------\n");
    fprintf(stderr, "Program %s, %s\n", "TEST", "VERSION 0.1");
    fprintf(stderr, "In source file %s, line %d\n\n", file, line);
    if (retcode < 0 || retcode >= eeUnknownError)
    {
        retcode = eeUnknownError;
    }
    fprintf(stderr, "%s:\n", error_names[retcode]);
    fprintf(stderr, "%s\n", msg);
    fprintf(stderr, "For more information and tips for troubleshooting, please check the GROMACS\n"
                    "website at http://www.gromacs.org/Documentation/Errors");
    fprintf(stderr, "\n-------------------------------------------------------\n");
    std::exit(1);
}

static ErrorHandlerFunc error_handler = standardErrorHandler;

ErrorHandlerFunc setFatalErrorHandler(ErrorHandlerFunc handler)
{
    // TODO: Acquire a mutex here
    ErrorHandlerFunc old_handler = error_handler;
    error_handler = handler;
    // TODO: Release the mutex here
    return old_handler;
}

void fatalError(int retcode, const char *msg, const char *file, int line)
{
    // TODO: Acquire a mutex here
    ErrorHandlerFunc handler = error_handler;
    // TODO: Release the mutex here
    if (handler != NULL)
    {
        handler(retcode, msg, file, line);
    }
}

void fatalErrorFormatted(int retcode, const char *file, int line,
                         const char *fmt, ...)
{
    // TODO: Acquire a mutex here
    ErrorHandlerFunc handler = error_handler;
    // TODO: Release the mutex here
    if (handler != NULL)
    {
        char    errorbuf[10000];
        va_list va;

        va_start(va, fmt);
        // FIXME: vsnprintf should really be used here
        std::vsprintf(errorbuf, fmt, va);
        va_end(va);
        handler(retcode, errorbuf, file, line);
    }
}

} // namespace gmx
