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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _nrjac_h
#define _nrjac_h

static char *SRCID_nrjac_h = "$Id$";

extern void jacobi(double **a,int n,double d[],double **v,int *nrot);
/* 
 * real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
 * int     natoms = number of rows and columns
 * real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
 * real       **v = v[0..n-1][0..n-1] contains the vectors in columns
 * int      *irot = number of jacobi rotations
 */
#endif
