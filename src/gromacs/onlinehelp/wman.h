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
#ifndef GMX_ONLINEHELP_WMAN_H
#define GMX_ONLINEHELP_WMAN_H

#include "readinp.h"
#include "gromacs/fileio/filenm.h"

#ifdef __cplusplus
extern "C" {
#endif

void write_man(const char *mantp, const char *program,
               int nldesc, const char **desc,
               int nfile, t_filenm *fnm,
               int npargs, t_pargs *pa,
               int nbug, const char **bugs);

char *fileopt(unsigned long flag, char buf[], int maxsize);
/* Return a string describing the file type in flag.
 * flag should the flag field of a filenm struct.
 * You have to provide a buffer and buffer length in which
 * the result will be written. The returned pointer is just
 * a pointer to this buffer.
 */

char *check_nroff(const char *s);
char *check_html(const char *s);
char *check_tty(const char *s);

#ifdef __cplusplus
}
#endif

#endif
