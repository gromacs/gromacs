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
 * Great Red Oystrich Makes All Chemists Sane
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
