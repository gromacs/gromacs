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

#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <fftw_mpi.h>

#include "fftw_sched.h"
#include "TOMS_transpose.h"

/**************************************************************************/

static int transpose_mpi_get_block_size(int n, int n_pes)
{
     return((n + n_pes - 1) / n_pes);
}

void transpose_mpi_get_local_size(int n, int my_pe, int n_pes,
				  int *local_n, int *local_start)
{
     int block_size;

     block_size = transpose_mpi_get_block_size(n,n_pes);
     n_pes = (n + block_size - 1) / block_size;

     if (my_pe >= n_pes) {
	  *local_n = 0;
	  *local_start = 0;
     }
     else {
	  *local_start = my_pe * block_size;
	  if (my_pe == n_pes - 1)
	  *local_n = n - *local_start;
     else
	  *local_n = block_size;
     }
}

#define MAX2(a,b) ((a) > (b) ? (a) : (b))

int transpose_mpi_get_local_storage_size(int nx, int ny,
					 int my_pe, int n_pes)
{
     int local_nx, local_ny, local_x_start, local_y_start;

     transpose_mpi_get_local_size(nx,my_pe,n_pes,&local_nx,&local_x_start);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,&local_ny,&local_y_start);

     return MAX2(1, MAX2(local_nx*ny, nx*local_ny));
}

static int gcd(int a, int b)
{
        int r;
        do {
                r = a % b;
                a = b;
                b = r;
        } while (r != 0);
        return a;
}

/**************************************************************************/

transpose_mpi_plan transpose_mpi_create_plan(int nx, int ny,
					     MPI_Comm transpose_comm)
{
     transpose_mpi_plan p;
     int my_pe, n_pes, pe;
     int x_block_size, local_nx, local_x_start;
     int y_block_size, local_ny, local_y_start;
     transpose_mpi_exchange *exchange = 0;
     int step, send_block_size = 0, recv_block_size = 0, num_steps = 0;
     int **sched, sched_npes, sched_sortpe, sched_sort_ascending = 0;
     int *perm_block_dest = NULL;
     int num_perm_blocks = 0, perm_block_size = 0, perm_block;
     char *move = NULL;
     int move_size = 0;
     int *send_block_sizes = 0, *send_block_offsets = 0;
     int *recv_block_sizes = 0, *recv_block_offsets = 0;
     MPI_Comm comm;

     /* create a new "clone" communicator so that transpose
        communications do not interfere with caller communications. */
     MPI_Comm_dup(transpose_comm, &comm);

     MPI_Comm_rank(comm,&my_pe);
     MPI_Comm_size(comm,&n_pes);

     /* work space for in-place transpose routine: */
     move_size = (nx + ny) / 2;
     move = (char *) fftw_malloc(sizeof(char) * move_size);

     x_block_size = transpose_mpi_get_block_size(nx,n_pes);
     transpose_mpi_get_local_size(nx,my_pe,n_pes,&local_nx,&local_x_start);
     y_block_size = transpose_mpi_get_block_size(ny,n_pes);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,&local_ny,&local_y_start);

     /* allocate pre-computed post-transpose permutation: */
     perm_block_size = gcd(nx,x_block_size);
     num_perm_blocks = (nx / perm_block_size) * local_ny;
     perm_block_dest = (int *) fftw_malloc(sizeof(int) * num_perm_blocks);
     for (perm_block = 0; perm_block < num_perm_blocks; ++perm_block)
	  perm_block_dest[perm_block] = num_perm_blocks; 

     /* allocate block sizes/offsets arrays for out-of-place transpose: */
     send_block_sizes = (int *) fftw_malloc(n_pes * sizeof(int));
     send_block_offsets = (int *) fftw_malloc(n_pes * sizeof(int));
     recv_block_sizes = (int *) fftw_malloc(n_pes * sizeof(int));
     recv_block_offsets = (int *) fftw_malloc(n_pes * sizeof(int));
     for (step = 0; step < n_pes; ++step)
	  send_block_sizes[step] = send_block_offsets[step] =
	       recv_block_sizes[step] = recv_block_offsets[step] = 0;

     if (local_nx > 0 || local_ny > 0) {

     sched_npes = n_pes;
     sched_sortpe = -1;
     for (pe = 0; pe < n_pes; ++pe) {
	  int pe_nx, pe_x_start, pe_ny, pe_y_start;

	  transpose_mpi_get_local_size(nx,pe,n_pes,
				       &pe_nx,&pe_x_start);
	  transpose_mpi_get_local_size(ny,pe,n_pes,
				       &pe_ny,&pe_y_start);

	  if (pe_nx == 0 && pe_ny == 0) {
	       sched_npes = pe;
	       break;
	  }
	  else if (pe_nx * y_block_size != pe_ny * x_block_size
		   && pe_nx != 0 && pe_ny != 0) {
	       if (sched_sortpe != -1)
		    fftw_mpi_die("BUG: More than one PE needs sorting!\n");
	       sched_sortpe = pe;
	       sched_sort_ascending = 
		    pe_nx * y_block_size > pe_ny * x_block_size;
	  }
     }

     sched = make_comm_schedule(sched_npes);
     if (!sched) {
	  MPI_Comm_free(&comm);
	  return 0;
     }

     if (sched_sortpe != -1) {
	  sort_comm_schedule(sched,sched_npes,sched_sortpe);
	  if (!sched_sort_ascending)
	       invert_comm_schedule(sched,sched_npes);
     }

     send_block_size = local_nx * y_block_size;
     recv_block_size = local_ny * x_block_size;

     num_steps = sched_npes;

     exchange = (transpose_mpi_exchange *)
	  fftw_malloc(num_steps * sizeof(transpose_mpi_exchange));
     if (!exchange) {
	  free_comm_schedule(sched,sched_npes);
	  MPI_Comm_free(&comm);
	  return 0;
     }

     for (step = 0; step < sched_npes; ++step) {
	  int dest_pe;
	  int dest_nx, dest_x_start;
	  int dest_ny, dest_y_start;
	  int num_perm_blocks_received, num_perm_rows_received;
	  
	  exchange[step].dest_pe = dest_pe =
	       exchange[step].block_num = sched[my_pe][step];

	  if (exchange[step].block_num == -1)
	       fftw_mpi_die("BUG: schedule ended too early.\n");

	       transpose_mpi_get_local_size(nx,dest_pe,n_pes,
					    &dest_nx,&dest_x_start);
	  transpose_mpi_get_local_size(ny,dest_pe,n_pes,
				       &dest_ny,&dest_y_start);
	       
	  exchange[step].send_size = local_nx * dest_ny;
	  exchange[step].recv_size = dest_nx * local_ny;
	  
	  send_block_sizes[dest_pe] = exchange[step].send_size;
	  send_block_offsets[dest_pe] = dest_pe * send_block_size;
	  recv_block_sizes[dest_pe] = exchange[step].recv_size;
	  recv_block_offsets[dest_pe] = dest_pe * recv_block_size;

	  /* Precompute the post-transpose permutation (ugh): */
          if (exchange[step].recv_size > 0) {
               num_perm_blocks_received =
		    exchange[step].recv_size / perm_block_size;
	       num_perm_rows_received = num_perm_blocks_received / local_ny;

	       for (perm_block = 0; perm_block < num_perm_blocks_received;
		    ++perm_block) {
		    int old_block, new_block;
		    
                    old_block = perm_block + exchange[step].block_num *
                         (recv_block_size / perm_block_size);
		    new_block = perm_block % num_perm_rows_received +
                         dest_x_start / perm_block_size +
                         (perm_block / num_perm_rows_received)
			 * (nx / perm_block_size);
		    
		    if (old_block >= num_perm_blocks ||
			new_block >= num_perm_blocks)
			 fftw_mpi_die("bad block index in permutation!");

		    perm_block_dest[old_block] = new_block;
	       }
	  }
     }

     free_comm_schedule(sched,sched_npes);

     } /* if (local_nx > 0 || local_ny > 0) */

     p = (transpose_mpi_plan) 
	  fftw_malloc(sizeof(transpose_mpi_plan_struct));
     if (!p) {
          fftw_free(exchange);
	  MPI_Comm_free(&comm);
	  return 0;
     }

     p->comm = comm;
     p->nx = nx; p->ny = ny;
     p->local_nx = local_nx; p->local_ny = local_ny;

     p->my_pe = my_pe; p->n_pes = n_pes;

     p->exchange = exchange;
     p->send_block_size = send_block_size;
     p->recv_block_size = recv_block_size;
     p->num_steps = num_steps;

     p->perm_block_dest = perm_block_dest;
     p->num_perm_blocks = num_perm_blocks;
     p->perm_block_size = perm_block_size;

     p->move = move;
     p->move_size = move_size;

     p->send_block_sizes = send_block_sizes;
     p->send_block_offsets = send_block_offsets;
     p->recv_block_sizes = recv_block_sizes;
     p->recv_block_offsets = recv_block_offsets;

     p->all_blocks_equal = send_block_size * n_pes * n_pes == nx * ny &&
	                   recv_block_size * n_pes * n_pes == nx * ny;
     if (p->all_blocks_equal)
	  for (step = 0; step < n_pes; ++step)
	       if (send_block_sizes[step] != send_block_size ||
		   recv_block_sizes[step] != recv_block_size) {
		    p->all_blocks_equal = 0;
		    break;
	       }
     if (nx % n_pes == 0 && ny % n_pes == 0 && !p->all_blocks_equal)
	  fftw_mpi_die("n_pes divided dimensions but blocks are unequal!");

     /* Set the type constant for passing to the MPI routines;
	here, we assume that TRANSPOSE_EL_TYPE is one of the
	floating-point types. */
     if (sizeof(TRANSPOSE_EL_TYPE) == sizeof(double))
	  p->el_type = MPI_DOUBLE;
     else if (sizeof(TRANSPOSE_EL_TYPE) == sizeof(float))
	  p->el_type = MPI_FLOAT;
     else
	  fftw_mpi_die("Unknown TRANSPOSE_EL_TYPE!\n");

     return p;
}

/**************************************************************************/

void transpose_mpi_destroy_plan(transpose_mpi_plan p)
{
     if (p) {
	  if (p->exchange)
	       fftw_free(p->exchange);
	  if (p->perm_block_dest)
               fftw_free(p->perm_block_dest);
	  if (p->move)
	       fftw_free(p->move);
	  if (p->send_block_sizes)
	       fftw_free(p->send_block_sizes);
	  if (p->send_block_offsets)
	       fftw_free(p->send_block_offsets);
	  if (p->recv_block_sizes)
	       fftw_free(p->recv_block_sizes);
	  if (p->recv_block_offsets)
	       fftw_free(p->recv_block_offsets);
          MPI_Comm_free(&p->comm);
	  fftw_free(p);
     }
}

/**************************************************************************/

static void exchange_elements(TRANSPOSE_EL_TYPE *buf1,
			      TRANSPOSE_EL_TYPE *buf2, int n)
{
     int i;
     TRANSPOSE_EL_TYPE t0,t1,t2,t3;

     for (i = 0; i < (n & 3); ++i) {
	  t0 = buf1[i];
	  buf1[i] = buf2[i];
	  buf2[i] = t0;
     }
     for (; i < n; i += 4) {
	  t0 = buf1[i];
	  t1 = buf1[i+1];
	  t2 = buf1[i+2];
	  t3 = buf1[i+3];
	  buf1[i] = buf2[i];
	  buf1[i+1] = buf2[i+1];
	  buf1[i+2] = buf2[i+2];
	  buf1[i+3] = buf2[i+3];
	  buf2[i] = t0;
	  buf2[i+1] = t1;
	  buf2[i+2] = t2;
	  buf2[i+3] = t3;
     }
}

static void do_permutation(TRANSPOSE_EL_TYPE *data,
			   int *perm_block_dest,
			   int num_perm_blocks,
			   int perm_block_size)
{
     int start_block;

     /* Perform the permutation in the perm_block_dest array, following
	the cycles around and *changing* the perm_block_dest array
	to reflect the permutations that have already been performed.
	At the end of this routine, we change the perm_block_dest
        array back to its original state. (ugh) */

     for (start_block = 0; start_block < num_perm_blocks; ++start_block) {
	  int cur_block = start_block;
	  int new_block = perm_block_dest[start_block];
	  
	  while (new_block > 0 && new_block < num_perm_blocks &&
		 new_block != start_block) {
	       exchange_elements(data + perm_block_size*start_block,
				data + perm_block_size*new_block,
				perm_block_size);
	       perm_block_dest[cur_block] = -1 - new_block;
	       cur_block = new_block;
	       new_block = perm_block_dest[cur_block];
	  }

          if (new_block == start_block)
	       perm_block_dest[cur_block] = -1 - new_block;
     }

     /* reset the permutation array (ugh): */
     for (start_block = 0; start_block < num_perm_blocks; ++start_block)
	  perm_block_dest[start_block] = -1 - perm_block_dest[start_block];
}

TRANSPOSE_EL_TYPE *transpose_allocate_send_buf(transpose_mpi_plan p,
					       int el_size)
{
     TRANSPOSE_EL_TYPE *send_buf = 0;

     /* allocate the send buffer: */
     if (p->send_block_size > 0) {
          send_buf = (TRANSPOSE_EL_TYPE *)
               fftw_malloc(p->send_block_size * el_size
                                * sizeof(TRANSPOSE_EL_TYPE));
          if (!send_buf)
               fftw_mpi_die("Out of memory!\n");
     }

     return send_buf;
}

void transpose_in_place_local(transpose_mpi_plan p,
			      int el_size, TRANSPOSE_EL_TYPE *local_data,
			      transpose_in_place_which which)
{
     switch (which) {
	 case BEFORE_TRANSPOSE:
	      if (el_size == 1)
		   TOMS_transpose_2d(local_data,
				     p->local_nx, p->ny,
				     p->move, p->move_size);
	      else
		   TOMS_transpose_2d_arbitrary(local_data,
					       p->local_nx, p->ny,
					       el_size,
					       p->move, p->move_size);
	      break;
	 case AFTER_TRANSPOSE:
	      do_permutation(local_data, p->perm_block_dest,
			     p->num_perm_blocks, p->perm_block_size * el_size);
	      break;
     }
}			      

/**************************************************************************/

static void local_transpose_copy(TRANSPOSE_EL_TYPE *src, 
				 TRANSPOSE_EL_TYPE *dest,
				 int el_size, int nx, int ny)
{
     int x, y;

	  if (el_size == 1)
	  for (x = 0; x < nx; ++x)
	       for (y = 0; y < ny; ++y)
		    dest[y * nx + x] = src[x * ny + y];
     else if (el_size == 2)
	  for (x = 0; x < nx; ++x)
	       for (y = 0; y < ny; ++y) {
		    dest[y * (2 * nx) + 2*x]     = src[x * (2 * ny) + 2*y];
		    dest[y * (2 * nx) + 2*x + 1] = src[x * (2 * ny) + 2*y + 1];
	       }
	  else
	  for (x = 0; x < nx; ++x)
	       for (y = 0; y < ny; ++y)
		    memcpy(&dest[y * (el_size*nx) + (el_size*x)],
			   &src[x * (el_size*ny) + (el_size*y)],
			   el_size * sizeof(TRANSPOSE_EL_TYPE));

}

/* Out-of-place version of transpose_mpi (or rather, in place using
   a scratch array): */
static void transpose_mpi_out_of_place(transpose_mpi_plan p, int el_size,
				       TRANSPOSE_EL_TYPE *local_data,
				       TRANSPOSE_EL_TYPE *work)
{
     local_transpose_copy(local_data, work, el_size, p->local_nx, p->ny);

     if (p->all_blocks_equal)
	  MPI_Alltoall(work, p->send_block_size * el_size, p->el_type,
		       local_data, p->recv_block_size * el_size, p->el_type,
		       p->comm);
     else {
	  int i, n_pes = p->n_pes;

	  for (i = 0; i < n_pes; ++i) {
	       p->send_block_sizes[i] *= el_size;
	       p->recv_block_sizes[i] *= el_size;
	       p->send_block_offsets[i] *= el_size;
	       p->recv_block_offsets[i] *= el_size;
     }
	  MPI_Alltoallv(work, p->send_block_sizes, p->send_block_offsets,
			p->el_type,
			local_data, p->recv_block_sizes, p->recv_block_offsets,
			p->el_type,
			p->comm);
	  for (i = 0; i < n_pes; ++i) {
	       p->send_block_sizes[i] /= el_size;
	       p->recv_block_sizes[i] /= el_size;
	       p->send_block_offsets[i] /= el_size;
	       p->recv_block_offsets[i] /= el_size;
	  }
     }

     do_permutation(local_data, p->perm_block_dest, p->num_perm_blocks,
		    p->perm_block_size * el_size);
}

/**************************************************************************/

void transpose_mpi(transpose_mpi_plan p, int el_size,
		   TRANSPOSE_EL_TYPE *local_data,
		   TRANSPOSE_EL_TYPE *work)
{
     /* if local_data and work are both NULL, we have no way of knowing
	whether we should use in-place or out-of-place transpose routine;
	if we guess wrong, MPI_Alltoall will block.  We prevent this
	by making sure that transpose_mpi_get_local_storage_size returns
	at least 1. */
     if (!local_data && !work)
	  fftw_mpi_die("local_data and work are both NULL!");

     if (work)
	  transpose_mpi_out_of_place(p, el_size, local_data, work);
     else if (p->local_nx > 0 || p->local_ny > 0) {
	  int step;
	  TRANSPOSE_EL_TYPE *send_buf = transpose_allocate_send_buf(p,el_size);

	  transpose_in_place_local(p, el_size, local_data, BEFORE_TRANSPOSE);

	  for (step = 0; step < p->num_steps; ++step) {
               transpose_finish_exchange_step(p, step - 1);
	       
               transpose_start_exchange_step(p, el_size, local_data,
                                             send_buf, step, TRANSPOSE_SYNC);
     }

	  transpose_finish_exchange_step(p, step - 1);

	  transpose_in_place_local(p, el_size, local_data, AFTER_TRANSPOSE);
     
	  if (send_buf)
	       fftw_free(send_buf);
     } /* if (local_nx > 0 || local_ny > 0) */
}
	       
/**************************************************************************/

/* non-blocking routines for overlapping communication and computation: */

#define USE_SYNCHRONOUS_ISEND 1
#if USE_SYNCHRONOUS_ISEND
#define ISEND MPI_Issend
#else
#define ISEND MPI_Isend
#endif

void transpose_get_send_block(transpose_mpi_plan p, int step,
			      int *block_y_start, int *block_ny)
{
     if (p->local_nx > 0) {
	  *block_y_start = 
	       p->send_block_size / p->local_nx * p->exchange[step].block_num;
	  *block_ny = p->exchange[step].send_size / p->local_nx;
     }
     else {
	  *block_y_start = 0;
	  *block_ny = 0;
     }
}

void transpose_start_exchange_step(transpose_mpi_plan p,
				   int el_size,
				   TRANSPOSE_EL_TYPE *local_data,
				   TRANSPOSE_EL_TYPE *send_buf,
				   int step,
				   transpose_sync_type sync_type)
{
     if (p->local_nx > 0 || p->local_ny > 0) {
	  transpose_mpi_exchange *exchange = p->exchange;
	  int block = exchange[step].block_num;
	  int send_block_size = p->send_block_size;
	  int recv_block_size = p->recv_block_size;
	  
          if (exchange[step].dest_pe != p->my_pe) {

	       /* first, copy to send buffer: */
	       if (exchange[step].send_size > 0)
		    memcpy(send_buf,
			   local_data + el_size*send_block_size*block,
			   el_size * exchange[step].send_size *
			   sizeof(TRANSPOSE_EL_TYPE));
				 
#define DO_ISEND  \
               if (exchange[step].send_size > 0) {  \
			 ISEND(send_buf, \
			       exchange[step].send_size * el_size, \
			       p->el_type, \
			       exchange[step].dest_pe, 0, \
			       p->comm, \
			       &p->request[0]); \
	       }
 
	       p->request[0] = MPI_REQUEST_NULL;
	       p->request[1] = MPI_REQUEST_NULL;

	       if (sync_type == TRANSPOSE_ASYNC) {
		    /* Note that we impose an ordering on the sends and
		       receives (lower pe sends first) so that we won't
		       have deadlock if Isend & Irecv are blocking in some
		       MPI implementation: */
	  
		    if (p->my_pe < exchange[step].dest_pe)
			 DO_ISEND;
		    
		    if (exchange[step].recv_size > 0) {
			 MPI_Irecv(local_data + el_size*recv_block_size*block,
				   exchange[step].recv_size * el_size,
				   p->el_type,
				   exchange[step].dest_pe, MPI_ANY_TAG,
				   p->comm,
				   &p->request[1]);
	       }
	       
		    if (p->my_pe > exchange[step].dest_pe)
			 DO_ISEND;
	  }
	       else /* (sync_type == TRANSPOSE_SYNC) */ {
		    MPI_Status status;

		    MPI_Sendrecv(send_buf,
				 exchange[step].send_size * el_size,
				 p->el_type,
				 exchange[step].dest_pe, 0,

				 local_data + el_size*recv_block_size*block,
				 exchange[step].recv_size * el_size,
				 p->el_type,
				 exchange[step].dest_pe, MPI_ANY_TAG,

				 p->comm, &status);
	       }
	  }
	  else if (exchange[step].recv_size > 0 &&
		   recv_block_size != send_block_size)
	       memmove(local_data + el_size*recv_block_size*block,
		       local_data + el_size*send_block_size*block,
		       exchange[step].recv_size * el_size *
		       sizeof(TRANSPOSE_EL_TYPE));
     }
}

void transpose_finish_exchange_step(transpose_mpi_plan p, int step)
{
     if ((p->local_nx > 0 || p->local_ny > 0) && step >= 0
	 && p->exchange[step].dest_pe != p->my_pe) {
	  MPI_Status status[2];

	  MPI_Waitall(2,p->request,status);
     }
}

