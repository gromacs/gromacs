/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef FFTW_MPI_H
#define FFTW_MPI_H

#include <fftw.h>
#include <mpi.h> /* need access to the MPI type definitions */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***********************************************************************/

typedef fftw_real TRANSPOSE_EL_TYPE;

typedef struct {
     int block_num, dest_pe, send_size, recv_size;
} transpose_mpi_exchange;

typedef struct {
     MPI_Comm comm;
     int n_pes, my_pe;
     
     int nx,ny,local_nx,local_ny;

     transpose_mpi_exchange *exchange;
     int num_steps, send_block_size, recv_block_size;

     MPI_Datatype el_type;

     MPI_Request request[2];

     int *perm_block_dest;
     int num_perm_blocks, perm_block_size;

     int all_blocks_equal;
     int *send_block_sizes, *send_block_offsets;
     int *recv_block_sizes, *recv_block_offsets;

     char *move;
     int move_size;
} transpose_mpi_plan_struct;

typedef transpose_mpi_plan_struct *transpose_mpi_plan;

extern void transpose_mpi_get_local_size(int n, int my_pe, int n_pes,
					 int *local_n, int *local_start);
extern int transpose_mpi_get_local_storage_size(int nx, int ny,
						int my_pe, int n_pes);

extern transpose_mpi_plan transpose_mpi_create_plan(int nx, int ny,
						    MPI_Comm comm);
extern void transpose_mpi_destroy_plan(transpose_mpi_plan p);

extern void transpose_mpi(transpose_mpi_plan p, int el_size,
			  TRANSPOSE_EL_TYPE *local_data,
			  TRANSPOSE_EL_TYPE *work);

typedef enum { BEFORE_TRANSPOSE, AFTER_TRANSPOSE } transpose_in_place_which;

typedef enum { TRANSPOSE_SYNC, TRANSPOSE_ASYNC } transpose_sync_type;

extern void transpose_in_place_local(transpose_mpi_plan p,
                              int el_size, TRANSPOSE_EL_TYPE *local_data,
                              transpose_in_place_which which);

extern TRANSPOSE_EL_TYPE *transpose_allocate_send_buf(transpose_mpi_plan p,
						      int el_size);
extern void transpose_get_send_block(transpose_mpi_plan p, int step,
				     int *block_y_start, int *block_ny);
extern void transpose_start_exchange_step(transpose_mpi_plan p,
					  int el_size,
					  TRANSPOSE_EL_TYPE *local_data,
					  TRANSPOSE_EL_TYPE *send_buf,
					  int step,
					  transpose_sync_type sync_type);
extern void transpose_finish_exchange_step(transpose_mpi_plan p, int step);

/***********************************************************************/

typedef struct {
     fftw_plan p_fft_x;  /* plan for first dimension */
     fftwnd_plan p_fft;  /* plan for subsequent dimensions */
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_complex *work; /* extra workspace, if needed */
} fftwnd_mpi_plan_data;

typedef fftwnd_mpi_plan_data *fftwnd_mpi_plan;

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
		       int n_fields,
		       fftw_complex *local_data, fftw_complex *work,
		       fftwnd_mpi_output_order output_order);

extern void fftw_mpi_die(const char *error_string);

/***********************************************************************/

typedef struct fftw_mpi_twiddle_struct {
     int rows, rowstart, cols, n;
     fftw_complex *W;
     int refcount;
     struct fftw_mpi_twiddle_struct *next;
} fftw_mpi_twiddle;

typedef struct fftw_mpi_plan_struct {
     int n, m, r, local_m, local_m_start, local_r, local_r_start;
     fftw_complex *fft_work;
     fftw_mpi_twiddle *tw;
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_plan pm, pr;
     int flags;
     fftw_direction dir;
} *fftw_mpi_plan;

/* new flags for the MPI planner: */
#define FFTW_SCRAMBLED_INPUT (8192)
#define FFTW_SCRAMBLED_OUTPUT (16384)

extern void fftw_mpi_local_sizes(fftw_mpi_plan p,
				 int *local_n,
				 int *local_start,
				 int *local_n_after_transform,
				 int *local_start_after_transform,
				 int *total_local_size);

extern fftw_mpi_plan fftw_mpi_create_plan(MPI_Comm comm,
					  int n,
					  fftw_direction dir, int flags);

extern void fftw_mpi_destroy_plan(fftw_mpi_plan p);

extern void fftw_mpi(fftw_mpi_plan p, int n_fields,
		     fftw_complex *local_data, fftw_complex *work);

extern void fftw_mpi_print_plan(fftw_mpi_plan p);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTW_MPI_H */
