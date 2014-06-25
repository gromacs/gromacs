/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _md_logging_h
#define _md_logging_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

struct t_commrec;

void md_print_info(const struct t_commrec *cr, FILE *fplog,
                   const char *fmt, ...);
/* Print an general information message to stderr on the master node
 * and to fplog if fplog!=NULL.
 * fmt is a standard printf formatting string which should end in \n,
 * the arguments after that contain the values to be printed, as in printf.
 */

void md_print_warn(const struct t_commrec *cr, FILE *fplog,
                   const char *fmt, ...);
/* As md_print_info above, but for important notices or warnings.
 * The only difference with md_print_info is that a newline is printed
 * before and after the message such that it stands out.
 */

#ifdef __cplusplus
}
#endif

#endif  /* _md_logging_h */
