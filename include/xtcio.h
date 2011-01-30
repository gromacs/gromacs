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

#ifndef _xtcio_h
#define _xtcio_h

#include "typedefs.h"
#include "gmxfio.h"
#include "xdrf.h"

#ifdef __cplusplus
extern "C" {
#endif

/* All functions return 1 if successful, 0 otherwise 
 * bOK tells if a frame is not corrupted 
 */  

t_fileio *open_xtc(const char *filename,const char *mode);
/* Open a file for xdr I/O */
  
void close_xtc(t_fileio *fio);
/* Close the file for xdr I/O */
  
int read_first_xtc(t_fileio *fio,
			  int *natoms,int *step,real *time,
			  matrix box,rvec **x,real *prec,gmx_bool *bOK);
/* Open xtc file, read xtc file first time, allocate memory for x */

int read_next_xtc(t_fileio *fio,
			 int natoms,int *step,real *time,
			 matrix box,rvec *x,real *prec,gmx_bool *bOK);
/* Read subsequent frames */

int write_xtc(t_fileio *fio,
		     int natoms,int step,real time,
		     matrix box,rvec *x,real prec);
/* Write a frame to xtc file */

int xtc_check(const char *str,gmx_bool bResult,const char *file,int line);
#define XTC_CHECK(s,b) xtc_check(s,b,__FILE__,__LINE__)

void xtc_check_fat_err(const char *str,gmx_bool bResult,const char *file,int line);
#define XTC_CHECK_FAT_ERR(s,b) xtc_check_fat_err(s,b,__FILE__,__LINE__)

#ifdef __cplusplus
}
#endif

#endif
