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
static char *SRCID_fftw_threads_h = "$Id$";

#ifndef FFTW_THREADS_H
#define FFTW_THREADS_H

#include <fftw.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/******************** User Interface *********************/

extern void fftw_threads(int nthreads,
		  fftw_plan plan, int howmany, fftw_complex *in, int istride,
		  int idist, fftw_complex *out, int ostride, int odist);
extern void fftwnd_threads(int nthreads,
			   fftwnd_plan plan, int howmany,
			   fftw_complex *in, int istride, int idist,
			   fftw_complex *out, int ostride, int odist);

extern void fftw_threads_one(int nthreads,
			     fftw_plan plan,
			     fftw_complex *in, fftw_complex *out);
extern void fftwnd_threads_one(int nthreads,
			       fftwnd_plan plan,
			       fftw_complex *in, fftw_complex *out);

extern int fftw_threads_init(void);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTW_THREADS_H */
