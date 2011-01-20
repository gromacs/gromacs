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
/*! \defgroup module_fatalerror Handling of Fatal Errors
 * \ingroup group_utilitymodules
 * \brief
 * Provides functions for handling fatal errors.
 *
 * Facilities for handling fatal errors are provided by the fatalerror.h header
 * file.  It provides a set of error codes (the enum ::ErrorCode) that should
 * be used for return codes in functions.  It also provides function fatalError()
 * for reporting the cause of the error to the user, and convenience macros
 * ::GMX_ERROR and ::GMX_ERROR_NORET for calling fatalError().  If the reason
 * string needs formatting, fatalErrorFormatted() is also provided.
 *
 * For users of the library, setFatalErrorHandler() is provided to alter the
 * behavior of fatalError() and fatalErrorFormatted().  fatalError() simply
 * calls the provided handler, while fatalErrorFormatted() does the formatting
 * internally and then calls the same handler.  The default handler prints the
 * reason of the error to standard error and aborts the execution.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
/*! \file
 * \brief
 * Public API convenience header for fatal error handling.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_fatalerror
 */
#ifndef GMX_FATALERROR_H
#define GMX_FATALERROR_H

#include "fatalerror/fatalerror.h"

#endif
