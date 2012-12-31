/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _md_logging_h
#define _md_logging_h

#include "types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

    void md_print_info(const t_commrec *cr, FILE *fplog,
                       const char *fmt, ...);
    /* Print an general information message to stderr on the master node
     * and to fplog if fplog!=NULL.
     * fmt is a standard printf formatting string which should end in \n,
     * the arguments after that contain the values to be printed, as in printf.
     */

    void md_print_warn(const t_commrec *cr, FILE *fplog,
                       const char *fmt, ...);
    /* As md_print_info above, but for important notices or warnings.
     * The only difference with md_print_info is that a newline is printed
     * before and after the message such that it stands out.
     */

#ifdef __cplusplus
}
#endif

#endif  /* _md_logging_h */
