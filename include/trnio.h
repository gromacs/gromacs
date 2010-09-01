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

#ifndef _trnio_h
#define _trnio_h

/**************************************************************
 *
 * These routines handle trj (trajectory) I/O, they read and
 * write trj/trr files. The routines should be able to read single
 * and double precision files without the user noting it.
 * The files are backward compatible, therefore the header holds
 * some unused variables.
 *
 * The routines in the corresponding c-file trnio.c
 * are based on the lower level routines in gmxfio.c
 * The integer file pointer returned from open_trn
 * can also be used with the routines in gmxfio.h
 *
 **************************************************************/
	
#include "typedefs.h"
#include "gmxfio.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct		/* This struct describes the order and the	*/
  /* sizes of the structs in a trjfile, sizes are given in bytes.	*/
{
  gmx_bool  bDouble;        /* Double precision?                            */
  int	ir_size;	/* Backward compatibility		        */
  int	e_size;		/* Backward compatibility		        */
  int	box_size;	/* Non zero if a box is present			*/
  int   vir_size;       /* Backward compatibility		        */
  int   pres_size;      /* Backward compatibility		        */
  int	top_size;	/* Backward compatibility		        */
  int	sym_size;	/* Backward compatibility		        */
  int	x_size;		/* Non zero if coordinates are present		*/
  int	v_size;		/* Non zero if velocities are present		*/
  int	f_size;		/* Non zero if forces are present		*/

  int	natoms;		/* The total number of atoms			*/
  int	step;		/* Current step number				*/
  int	nre;		/* Backward compatibility		        */
  real	t;		/* Current time					*/
  real	lambda;		/* Current value of lambda			*/
} t_trnheader;

t_fileio *open_trn(const char *fn,const char *mode);
/* Open a trj / trr file */

void close_trn(t_fileio *fio);
/* Close it */

gmx_bool fread_trnheader(t_fileio *fio,t_trnheader *trn,gmx_bool *bOK);
/* Read the header of a trn file. Return FALSE if there is no frame.
 * bOK will be FALSE when the header is incomplete.
 */

void read_trnheader(const char *fn,t_trnheader *header);
/* Read the header of a trn file from fn, and close the file afterwards. 
 */

void pr_trnheader(FILE *fp,int indent,char *title,t_trnheader *sh);
/* Print the header of a trn file to fp */

gmx_bool is_trn(FILE *fp);
/* Return true when the file is a trn file. File will be rewound
 * afterwards.
 */

void fwrite_trn(t_fileio *fio,int step,real t,real lambda,
		       rvec *box,int natoms,rvec *x,rvec *v,rvec *f);
/* Write a trn frame to file fp, box, x, v, f may be NULL */

gmx_bool fread_htrn(t_fileio *fio,t_trnheader *sh,
		       rvec *box,rvec *x,rvec *v,rvec *f);
/* Extern read a frame except the header (that should be pre-read,
 * using routine read_trnheader, see above) from a trn file.
 * Return FALSE on error
 */
 
gmx_bool fread_trn(t_fileio *fio,int *step,real *t,real *lambda,
		      rvec *box,int *natoms,rvec *x,rvec *v,rvec *f);
/* Read a trn frame, including the header from fp. box, x, v, f may
 * be NULL, in which case the data will be skipped over.
 * return FALSE on error
 */
 
void write_trn(const char *fn,int step,real t,real lambda,
		      rvec *box,int natoms,rvec *x,rvec *v,rvec *f);
/* Write a single trn frame to file fn, which is closed afterwards */

void read_trn(const char *fn,int *step,real *t,real *lambda,
		     rvec *box,int *natoms,rvec *x,rvec *v,rvec *f);
/* Read a single trn frame from file fn, which is closed afterwards 
 */

#ifdef __cplusplus
}
#endif


#endif
