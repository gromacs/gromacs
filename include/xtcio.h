/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _xtcio_h
#define _xtcio_h

static char *SRCID_xtcio_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "typedefs.h"
#include "xdrf.h"

/* All functions return 1 if succesfull, 0 otherwise 
 * bOK tells if a frame is not corrupted 
 */  

extern int open_xtc(char *filename,char *mode);
/* Open a file for xdr I/O */
  
extern void close_xtc(int fp);
/* Close the file for xdr I/O */
  
extern int read_first_xtc(int fp,
			  int *natoms,int *step,real *time,
			  matrix box,rvec **x,real *prec,bool *bOK);
/* Open xtc file, read xtc file first time, allocate memory for x */

extern int read_next_xtc(int fp,
			 int *natoms,int *step,real *time,
			 matrix box,rvec *x,real *prec,bool *bOK);
/* Read subsequent frames */

extern int write_xtc(int fp,
		     int natoms,int step,real time,
		     matrix box,rvec *x,real prec);
/* Write a frame to xtc file */

extern int xtc_check(char *str,bool bResult,char *file,int line);
#define XTC_CHECK(s,b) xtc_check(s,b,__FILE__,__LINE__)

extern void xtc_check_fat_err(char *str,bool bResult,char *file,int line);
#define XTC_CHECK_FAT_ERR(s,b) xtc_check_fat_err(s,b,__FILE__,__LINE__)

#ifdef CPLUSPLUS
}
#endif

#endif
