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
#ifndef RFFTW_MPI_H
#define RFFTW_MPI_H

static char *SRCID_rfftw_mpi_h = "$Id$";

#include <fftw_mpi.h>
#include <rfftw.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***********************************************************************/

typedef struct {
     fftw_plan p_fft_x;  /* plan for first dimension */
     rfftwnd_plan p_fft;  /* plan for subsequent dimensions */
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_complex *work; /* extra workspace, if needed */
} rfftwnd_mpi_plan_data;

typedef rfftwnd_mpi_plan_data *rfftwnd_mpi_plan;

extern rfftwnd_mpi_plan rfftwnd_mpi_create_plan(MPI_Comm comm,
					      int rank, const int *n,
					      fftw_direction dir,
					      int flags);
extern rfftwnd_mpi_plan rfftw2d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny,
					  fftw_direction dir, int flags);
extern rfftwnd_mpi_plan rfftw3d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny, int nz,
					  fftw_direction dir, int flags);

extern void rfftwnd_mpi_destroy_plan(rfftwnd_mpi_plan p);

extern void rfftwnd_mpi_local_sizes(rfftwnd_mpi_plan p,
				   int *local_nx,
				   int *local_x_start,
				   int *local_ny_after_transpose,
				   int *local_y_start_after_transpose,
				   int *total_local_size);

extern void rfftwnd_mpi(rfftwnd_mpi_plan p,
		       int n_fields,
		       fftw_real *local_data, fftw_real *work,
		       fftwnd_mpi_output_order output_order);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_MPI_H */
