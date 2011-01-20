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
/*! \defgroup module_errorreporting Reporting of Non-Fatal Errors
 * \ingroup group_utilitymodules
 * \brief
 * Provides functions and classes for reporting non-fatal errors.
 *
 * Facilities for customizable reporting of non-fatal errors are provided by an
 * abstract class AbstractErrorReporter.  Objects of this class can be passed
 * to functions for reporting non-fatal errors.  The caller can then create a
 * reporter object derived from AbstractErrorReporter to implement desired
 * behavior for reporting the errors to the user.  Two basic reporter classes
 * are provided by this module: EmptyErrorReporter that only counts the number
 * of errors occurred, and StandardErrorReporter that writes the errors to
 * standard error.
 *
 * An ErrorContext class is also provided to ease implementation of functions
 * that need to report errors that originate deeper in the call stack.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
/*! \file
 * \brief
 * Public API convenience header for non-fatal error reporting.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_errorreporting
 */
#ifndef GMX_ERRORREPORTING_H
#define GMX_ERRORREPORTING_H

#include "errorreporting/abstracterrorreporter.h"
#include "errorreporting/emptyerrorreporter.h"
#include "errorreporting/errorcontext.h"
#include "errorreporting/standarderrorreporter.h"

#endif
