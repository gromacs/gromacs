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
#ifndef RFFTW_MPI_THREADS_H
#define RFFTW_MPI_THREADS_H

static char *SRCID_rfftw_mpi_threads_h = "$Id$";

#include "rfftw_mpi.h"
#include "fftw_threads.h"
#include "rfftw.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***********************************************************************/

extern void rfftwnd_mpi_threads(int nthreads,
				    rfftwnd_mpi_plan p,
				    int n_fields,
				    fftw_real *local_data,
				    fftw_real *work,
				    fftwnd_mpi_output_order
				    output_order);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_MPI_THREADS_H */
