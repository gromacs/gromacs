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
static char *SRCID_fftwnd_mpi_h = "$Id$";

#ifndef FFTWND_MPI_H
#define FFTWND_MPI_H

#include "fftw.h"
#include "transpose_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct {
    fftwnd_plan p_fft;
    transpose_mpi_plan p_transpose, p_transpose_inv;
} fftwnd_mpi_aux_data;

typedef fftwnd_mpi_aux_data *fftwnd_mpi_plan;

typedef enum {
    FFTW_NORMAL_ORDER,
    FFTW_TRANSPOSED_ORDER
} fftwnd_mpi_output_order;

extern fftwnd_mpi_plan fftwnd_mpi_create_plan(MPI_Comm comm,
					      int rank, const int *n,
					      fftw_direction dir,
					      int flags);
extern fftwnd_mpi_plan fftw2d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny,
					      fftw_direction dir, int flags);
extern fftwnd_mpi_plan fftw3d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny, int nz,
					      fftw_direction dir, int flags);

extern void fftwnd_mpi_destroy_plan(fftwnd_mpi_plan p);

extern void fftwnd_mpi_local_sizes(fftwnd_mpi_plan p,
				   int *local_nx,
				   int *local_x_start,
				   int *local_ny_after_transpose,
				   int *local_y_start_after_transpose,
				   int *total_local_size);

extern void fftwnd_mpi(fftwnd_mpi_plan p,
		       int n_fields, fftw_complex * local_data,
		       fftwnd_mpi_output_order output_order);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTWND_MPI_H */
