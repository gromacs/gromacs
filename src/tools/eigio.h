/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * S  C  A  M  O  R  G
 */

#ifndef _eigio_h
#define _eigio_h

static char *SRCID_eigio_h = "$Id$";
#include "typedefs.h"

enum { eWXR_NO, eWXR_YES, eWXR_NOFIT };

extern void read_eigenvectors(char *file,int *natoms,bool *bFit,
			      rvec **xref,bool *bDMR,
			      rvec **xav,bool *bDMA,
			      int *nvec, int **eignr, rvec ***eigvec);
/* Read eigenvectors from file into eigvec, the eigenvector numbers   */
/* are stored in eignr.                                               */
/* When bFit=FALSE no fitting was used in the covariance analysis.    */
/* xref is the reference structure, can be NULL if not present.       */
/* bDMR indicates mass weighted fit.                                  */
/* xav is the average/minimum structure is written (t=0).             */
/* bDMA indicates mass weighted analysis/eigenvectors.                */ 

extern void write_eigenvectors(char *trnname,int natoms,real mat[],
			       bool bReverse,int begin,int end,
			       int WriteXref,rvec *xref,bool bDMR,
			       rvec xav[],bool bDMA);
/* Write eigenvectors in mat to a TRN file.                           */
/* The reference structure is written (t=-1) when WriteXref=eWXR_YES. */
/* When WriteXref==eWXR_NOFIT a zero frame is written (t=-1),         */
/* with lambda=-1.                                                    */ 
/* bDMR indicates mass weighted fit.                                  */ 
/* The average/minimum structure is written (t=0).                    */
/* bDMA indicates mass weighted analysis/eigenvectors.                */
/* eigenvectors with begin <= num <= end are written (num is base-1), */
/* the timestamp of eigenvector num is num.                           */
/* If bReverse==TRUE, num=1 is the last vector in mat.                */

#endif
