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
#ifndef RFFTW_THREADS_H
#define RFFTW_THREADS_H

static char *SRCID_rfftw_threads_h = "$Id$";

#include "rfftw.h"
#include "fftw_threads.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/******************** User Interface *********************/

extern void rfftw_threads(int nthreads,
                   fftw_plan plan, int howmany, fftw_real *in, int istride,
                   int idist, fftw_real *out, int ostride, int odist);
extern void rfftw_threads_one(int nthread, fftw_plan plan,
			      fftw_real *in, fftw_real *out);

extern void rfftwnd_threads_real_to_complex(int nthreads, fftwnd_plan p,
					    int howmany,
					    fftw_real *in,
					    int istride, int idist,
					    fftw_complex *out,
					    int ostride, int odist);
extern void rfftwnd_threads_complex_to_real(int nthreads, fftwnd_plan p,
					    int howmany,
					    fftw_complex *in,
					    int istride, int idist,
					    fftw_real *out,
					    int ostride, int odist);
extern void rfftwnd_threads_one_real_to_complex(int nthreads, fftwnd_plan p,
						fftw_real *in,
						fftw_complex *out);
extern void rfftwnd_threads_one_complex_to_real(int nthreads, fftwnd_plan p,
						fftw_complex *in,
						fftw_real *out);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_THREADS_H */
