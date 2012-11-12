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
 * Declares an internal helper function for formatting standard error messages.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ERRORFORMAT_H
#define GMX_UTILITY_ERRORFORMAT_H

#include <cstdio>

namespace gmx
{

/*! \cond internal */
namespace internal
{

/*! \internal \brief
 * Formats a common header for fatal error messages.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
void printFatalErrorHeader(FILE *fp, const char *title,
                           const char *func, const char *file, int line);
/*! \internal \brief
 * Formats a line of fatal error message text.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
void printFatalErrorMessageLine(FILE *fp, const char *text, int indent);
/*! \internal \brief
 * Formats a common footer for fatal error messages.
 *
 * Does not throw.
 *
 * \ingroup module_utility
 */
void printFatalErrorFooter(FILE *fp);

}   // namespace internal
//! \endcond

} // namespace gmx

#endif
