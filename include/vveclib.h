/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */

#ifndef _vveclib_h
#define _vveclib_h

static char *SRCID_vveclib_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) vveclib.h 1.11 2/2/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"

/* 
 * Interface module for vector library.
 * This is a set of routines that uses 
 * - libfv  on the i860
 * - veclib on the convex (future)
 * - some general C-routines
 *
 * The routines provided here just provide a standard header for
 * all systems.
 * All routine-names start with "v_"
 * All routines take type real which can be single or double
 * depending on options set and machine.
 * All routines take an integer N as first parameter which is
 * the vector length.
 */

extern double Flops;		/* The number of flops to date 	*/

#define FV_SAFE 12   /* Make your arrays FV_SAFE longer allways */

extern void v_scopy(int N,real s,real Z[]);
/* copy scalar to vector (Z[i] := s) */

extern void v_vcopy(int N,real X[],real Z[]);
/* copy vector X to vector Z (Z[i] := X[i]) */

extern real v_inner(int N,real X[],real Y[]);
/* calculate the inner product of X and Y */

extern void v_scopy_s(int N,real s,real Z[],int stride);
/* copy scalar to each stride's element of Z (Z[i x stride] := s */

extern void v_gather(int N,real B[],int IX[],real Z[]);
/* gather vector: Z[i] := B[IX[i]] */

extern void v_scatter(int N,real B[],int IX[],real Z[]);
/* scatter vector: Z[IX[i]] := B[i] */

extern void v_gather_us(int N,real B[],unsigned short IX[],real Z[]);
/* gather vector: Z[i] := B[IX[i]] */

extern void v_scatter_us(int N,real B[],unsigned short IX[],real Z[]);
/* scatter vector: Z[IX[i]] := B[i] */

extern void v_gather_us_s(int N,real B[],unsigned short IX[],real Z[],
			  int stride);
/* gather vector: Z[stride*i] := B[IX[stride*i]] */

extern void v_scatter_us_s(int N,real B[],unsigned short IX[],real Z[],
			   int stride);
/* scatter vector: Z[IX[stride*i]] := B[stride*i] */

extern void v_sub(int N,real X[],real Y[],real Z[]);
/* vector subtract: Z[i] := X[i]-Y[i] */

extern void v_mul(int N,real X[],real Y[],real Z[]);
/* vector multiply: Z[i] := X[i]*Y[i] */

extern void v_add(int N,real X[],real Y[],real Z[]);
/* vector add: Z[i] := X[i]+Y[i] */

/********************* Other routines ****************************/

void dprod(tensor T, rvec V, rvec W);
/* calculate the direct product: T = V(x)W */

void fac_dprod(tensor T, rvec V, rvec W, real fac);
/* calculate the direct product and every component fac
 * times to T: Tij += fac * [V(x)W]ij 
 */

#endif	/* _vveclib_h */
