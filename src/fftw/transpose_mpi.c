/*
 * Copyright (c) 1997,1998 Massachusetts Institute of Technology
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

#include "transpose_mpi.h"

#include "TOMS_transpose.h"

static int transpose_mpi_get_block_size(int n, int n_pes)
{
     if (n_pes > n)
	  n_pes = n;

     return((n + n_pes - 1) / n_pes);
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

void transpose_mpi_get_local_size(int n, int my_pe, int n_pes,
				  int *local_n, int *local_start)
{
     int block_size;

     if (n_pes > n)
	  n_pes = n;

     block_size = transpose_mpi_get_block_size(n,n_pes);
     *local_start = my_pe * block_size;

     if (*local_start >= n)
	  *local_n = 0;
     else if (*local_start + block_size >= n)
	  *local_n = n - *local_start;
     else
	  *local_n = block_size;
}

int transpose_mpi_get_local_storage_size(int nx, int ny,
					 int my_pe, int n_pes)
{
     int local_nx, local_ny, local_x_start, local_y_start;

     /* bug fix 3/26/98: force all processors to use the maximum
	local storage (the storage used by process 0).  Otherwise,
	the pad_blocks routine would go past the end of the array
	in some cases (e.g. 5x24 on 4 processors) 
	-- this is a quick hack, but oh well */
     my_pe = 0;

     transpose_mpi_get_local_size(nx,my_pe,n_pes,&local_nx,&local_x_start);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,&local_ny,&local_y_start);

     if (local_nx*ny > nx*local_ny)
	  return (local_nx*ny);
     else
	  return (nx*local_ny);
}

static int sort_block_compare(int src_pe, int dest_pe1, int dest_pe2)
{
     int xor1,xor2;

     xor1 = src_pe ^ dest_pe1;
     xor2 = src_pe ^ dest_pe2;
     if (xor1 > xor2)
	  return 1;
     else if (xor1 < xor2)
	  return -1;
     else
	  return 0;
}

static void sort_block_dests(int my_pe,
			     transpose_mpi_block_dest *block_dest,
			     int num_blocks)
{
     transpose_mpi_block_dest swap;
     int i,j;

     /* Do a simple N^2 sort of the blocks: */
     for (i = 0; i < num_blocks-1; ++i)
	  for (j = i+1; j < num_blocks; ++j)
	       if (sort_block_compare(my_pe,
				      block_dest[i].dest_pe,
				      block_dest[j].dest_pe) > 0) {
		    swap = block_dest[i];
		    block_dest[i] = block_dest[j];
		    block_dest[j] = swap;
	       }
}

#define WHICH_PE(dest_i,i_block_size) ((dest_i)/(i_block_size))

transpose_mpi_plan transpose_mpi_create_plan(int nx, int ny,
					     MPI_Comm transpose_comm)
{
     transpose_mpi_plan p;
     int my_pe, n_pes;
     int x_block_size, local_nx, local_x_start;
     int y_block_size, local_ny, local_y_start;
     transpose_mpi_block_dest *block_dest = 0;
     int block, block_size = 0, num_blocks = 0;
     int block_size_padded = 0;
     int *perm_block_dest = 0;
     int num_perm_blocks = 0, perm_block_size = 0, perm_block;
     int move_size;
     char *move = 0;
     MPI_Comm comm;

     /* create a new "clone" communicator so that transpose
        communications do not interfere with caller communications. */
     MPI_Comm_dup(transpose_comm, &comm);

     MPI_Comm_rank(comm,&my_pe);
     MPI_Comm_size(comm,&n_pes);

     x_block_size = transpose_mpi_get_block_size(nx,n_pes);
     transpose_mpi_get_local_size(nx,my_pe,n_pes,&local_nx,&local_x_start);
     y_block_size = transpose_mpi_get_block_size(ny,n_pes);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,&local_ny,&local_y_start);

     if (local_nx > 0 || local_ny > 0) {

     perm_block_size = gcd(nx,x_block_size);
     num_perm_blocks = (nx / perm_block_size) * local_ny;
     {
	  int tmp = (local_nx / perm_block_size) * ny;
	  if (tmp > num_perm_blocks)
	       num_perm_blocks = tmp;
     }

     perm_block_dest = (int *) malloc(sizeof(int) * num_perm_blocks);
     if (!perm_block_dest)
	  return 0;

     for (perm_block = 0; perm_block < num_perm_blocks; ++perm_block)
	  perm_block_dest[perm_block] = num_perm_blocks; 

     block_size = local_nx * y_block_size;

     num_blocks = n_pes;

     block_dest = (transpose_mpi_block_dest *)
	  malloc(num_blocks * sizeof(transpose_mpi_block_dest));
     if (!block_dest) {
	  free(perm_block_dest);
	  return 0;
     }

     if (local_nx < x_block_size && block_size < local_ny * x_block_size)
	  block_size_padded = local_ny * x_block_size;
     else
	  block_size_padded = block_size;

     for (block = 0; block < num_blocks; ++block) {
	  int block_x, block_y, block_new_x, block_new_y;
	  int num_perm_blocks_received, num_perm_rows_received;
	  int dest_pe;
	  
	  block_x = local_x_start;
	  block_y = block * y_block_size;

	  block_dest[block].block_num = block;

	  if (block_y >= ny)
	       block_dest[block].send_size = 0;
	  else if (block_y + y_block_size >= ny) /* the last block */
	       block_dest[block].send_size = local_nx*ny - block*block_size;
	  else
	       block_dest[block].send_size = block_size;

	  dest_pe = block;

	  if (dest_pe == my_pe) {
	       block_dest[block].dest_pe = dest_pe;
	       block_dest[block].recv_size = 
		    block_dest[block].send_size;
	       block_new_x = block_x;
	       block_new_y = block_y;
	  }
	  else {
	       int dest_nx, dest_x_start;

	       transpose_mpi_get_local_size(nx,dest_pe,n_pes,
					    &dest_nx,&dest_x_start);
	       
	       block_dest[block].dest_pe = dest_pe;
	       block_new_x = dest_x_start;
	       block_new_y = local_y_start;
	       block_dest[block].recv_size = dest_nx * local_ny;
	  }
	  
	  if (block_dest[block].recv_size > 0) {
	       num_perm_blocks_received = block_dest[block].recv_size /
		    perm_block_size;
	       num_perm_rows_received = num_perm_blocks_received / local_ny;

	       for (perm_block = 0; perm_block < num_perm_blocks_received;
		    ++perm_block) {
		    int old_block, new_block;
		    
		    old_block = perm_block + 
			 block * (block_size_padded / perm_block_size);
		    new_block = perm_block % num_perm_rows_received +
			 block_new_x / perm_block_size +
			 (block_new_y - local_y_start + 
			  perm_block / num_perm_rows_received) 
			 * (nx / perm_block_size);
		    
		    perm_block_dest[old_block] = new_block;
	       }
	  }
     }

     sort_block_dests(my_pe,block_dest,num_blocks);

     move_size = (local_nx + ny)/ 2;
     move = (char *) malloc(move_size*sizeof(char));
     if (!move) {
          free(block_dest);
          free(perm_block_dest);
          return 0;
     }

     }

     p = (transpose_mpi_plan) malloc(sizeof(transpose_mpi_plan_struct));
     if (!p) {
          free(block_dest);
	  free(perm_block_dest);
	  free(move);
	  return 0;
     }

     p->comm = comm;
     p->nx = nx; p->ny = ny;
     p->local_nx = local_nx; p->local_ny = local_ny;

     p->my_pe = my_pe; p->n_pes = n_pes;

     p->block_dest = block_dest;
     p->block_size = block_size;
     p->block_size_padded = block_size_padded;
     p->num_blocks = num_blocks;

     p->perm_block_dest = perm_block_dest;
     p->num_perm_blocks = num_perm_blocks;
     p->perm_block_size = perm_block_size;

     p->move_size = move_size;
     p->move = move;

     return p;
}

void transpose_mpi_destroy_plan(transpose_mpi_plan p)
{
     if (p) {
	  if (p->block_dest)
	       free(p->block_dest);
	  if (p->perm_block_dest)
	       free(p->perm_block_dest);
	  if (p->move)
	    free(p->move);
          MPI_Comm_free(&p->comm);
	  free(p);
     }
}

static void exchange_elements(TRANSPOSE_EL_TYPE *buf1, TRANSPOSE_EL_TYPE *buf2, int n)
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
	array back to its original state. */

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

     /* reset the permutation array: */
     for (start_block = 0; start_block < num_perm_blocks; ++start_block)
	  perm_block_dest[start_block] = -1 - perm_block_dest[start_block];
}

static void move_elements_forward(TRANSPOSE_EL_TYPE *buf1, TRANSPOSE_EL_TYPE *buf2, int n)
{
     int i;

     for (i = 0; i < (n & 3); ++i)
	  buf2[n-1-i] = buf1[n-1-i];
     for (i=n-1-i; i >= 0; i -= 4) {
	  buf2[i] = buf1[i];
	  buf2[i-1] = buf1[i-1];
	  buf2[i-2] = buf1[i-2];
	  buf2[i-3] = buf1[i-3];
     }
}

static void pad_blocks(TRANSPOSE_EL_TYPE *data, int block_size, int total_size, int shift)
{
     int num_blocks = total_size / block_size;
     int last_block_size = total_size % block_size;
     int block;

     if (last_block_size != 0)
	  move_elements_forward(data+block_size*num_blocks,
			       data+block_size*num_blocks+shift*num_blocks,
			       last_block_size);
     for (block = num_blocks - 1; block > 0; --block)
	  move_elements_forward(data+block_size*block,
			       data+block_size*block+shift*block,
			       block_size);
}

void transpose_mpi(transpose_mpi_plan p, TRANSPOSE_EL_TYPE *local_data, int el_size)
{
     int my_pe, n_pes;
     int nx, ny, local_nx, local_ny;
     transpose_mpi_block_dest *block_dest;
     int block_size, num_blocks, block, blk_num;
     int block_size_padded;
     int *perm_block_dest, num_perm_blocks, perm_block_size;
     MPI_Comm comm;
     TRANSPOSE_EL_TYPE *send_buf;

     my_pe = p->my_pe;
     n_pes = p->n_pes;
     nx = p->nx;
     ny = p->ny;
     local_nx = p->local_nx;
     local_ny = p->local_ny;

     if (local_nx > 0 || local_ny > 0) {

     block_dest = p->block_dest;
     block_size = p->block_size;
     block_size_padded = p->block_size_padded;
     num_blocks = p->num_blocks;
     perm_block_dest = p->perm_block_dest;
     num_perm_blocks = p->num_perm_blocks;
     perm_block_size = p->perm_block_size;
     comm = p->comm;

     /* first, transpose the local data: */

     if (local_nx > 0) {
	  if (el_size == 1)
	       TOMS_transpose_2d(local_data,local_nx,ny,p->move,p->move_size);
	  else
	       TOMS_transpose_2d_arbitrary(local_data,local_nx,ny,el_size,
					   p->move,p->move_size);
     }

     /* on last PE, may need to add padding in between blocks */
     if (block_size != block_size_padded && block_size > 0) {
	  int shift_size = block_size_padded - block_size;

	  pad_blocks(local_data,
		     block_size*el_size,
		     local_nx*ny*el_size,
		     shift_size*el_size);
     }

     /* now, exchange the blocks between processors*/

     send_buf = 0;
     
     for (block = 0; block < num_blocks; ++block)
	  if (block_dest[block].dest_pe != my_pe) {
	       MPI_Status status;
	       
	       blk_num = block_dest[block].block_num;

	       if (block_dest[block].send_size > 0 &&
		   block_dest[block].send_size == block_dest[block].recv_size)
		    MPI_Sendrecv_replace(local_data + 
					 el_size*block_size_padded*blk_num,
				    block_dest[block].send_size * el_size,
				    MPI_TRANSPOSE_EL_TYPE_CONSTANT,
				    block_dest[block].dest_pe, blk_num,
				    block_dest[block].dest_pe, MPI_ANY_TAG,
				    comm,
				    &status);
	       else if (block_dest[block].send_size > 0 &&
			block_dest[block].recv_size == 0)
		    MPI_Send(local_data +
			     el_size*block_size_padded*blk_num,
			     block_dest[block].send_size * el_size,
			     MPI_TRANSPOSE_EL_TYPE_CONSTANT,
			     block_dest[block].dest_pe, blk_num,
			     comm);
	       else if (block_dest[block].send_size == 0 &&
			block_dest[block].recv_size > 0)
		    MPI_Recv(local_data +
			     el_size*block_size_padded*blk_num,
			     block_dest[block].recv_size * el_size,
			     MPI_TRANSPOSE_EL_TYPE_CONSTANT,
			     block_dest[block].dest_pe, MPI_ANY_TAG,
			     comm,
			     &status);
	       else if (block_dest[block].send_size > 0 &&
			block_dest[block].recv_size > 0) {
		    if (send_buf == 0)
			 send_buf = (TRANSPOSE_EL_TYPE *)
			      malloc(block_dest[block].send_size * 
				     el_size * sizeof(TRANSPOSE_EL_TYPE));

		    memcpy(send_buf,local_data + 
			   el_size*block_size_padded*blk_num,
			   block_dest[block].send_size * el_size *
			   sizeof(TRANSPOSE_EL_TYPE));

		    MPI_Sendrecv(send_buf,
				 block_dest[block].send_size * el_size,
				 MPI_TRANSPOSE_EL_TYPE_CONSTANT,
				 block_dest[block].dest_pe, blk_num,

				 local_data + 
				 el_size*block_size_padded*blk_num,
				 block_dest[block].recv_size * el_size,
                                 MPI_TRANSPOSE_EL_TYPE_CONSTANT,
                                 block_dest[block].dest_pe, MPI_ANY_TAG,

				 comm,
				 &status);
				 
	       }
	  }

     if (send_buf)
	  free(send_buf);

     if (local_ny > 0)
	  do_permutation(local_data,perm_block_dest,num_perm_blocks,
			 perm_block_size*el_size);

     }
}

