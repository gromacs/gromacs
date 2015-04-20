/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _eigio_h
#define _eigio_h

#include "gromacs/legacyheaders/typedefs.h"

enum {
    eWXR_NO, eWXR_YES, eWXR_NOFIT
};

extern void read_eigenvectors(const char *file, int *natoms, gmx_bool *bFit,
                              rvec **xref, gmx_bool *bDMR,
                              rvec **xav, gmx_bool *bDMA,
                              int *nvec, int **eignr, rvec ***eigvec, real **eigval);
/* Read eigenvectors from file into eigvec, the eigenvector numbers   */
/* are stored in eignr.                                               */
/* When bFit=FALSE no fitting was used in the covariance analysis.    */
/* xref is the reference structure, can be NULL if not present.       */
/* bDMR indicates mass weighted fit.                                  */
/* xav is the average/minimum structure is written (t=0).             */
/* bDMA indicates mass weighted analysis/eigenvectors.                */

extern void write_eigenvectors(const char *trnname, int natoms, real mat[],
                               gmx_bool bReverse, int begin, int end,
                               int WriteXref, rvec *xref, gmx_bool bDMR,
                               rvec xav[], gmx_bool bDMA, real *eigval);
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


/* Read up to nmax eigenvalues from file fn, store the values in eigval[],
 * and the corresponding indices (start counting on 0) in eigvalnr[].
 * Returns the number of values read.
 */
int read_eigval  (const char *          fn,
                  int                   nmax,
                  int                   eigvalnr[],
                  real                  eigval[]);

#endif
