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
#include <stdio.h>
#include <mpi.h>

#include "transpose_mpi.h"

#define NX 14
#define NY 10

#define OUTPUT_PERMUTATION 0

int main(int argc, char **argv)
{
     int my_pe,n_pes;
     int i,j;
     int nx,ny;
     int local_x_start,local_nx,local_y_start,local_ny, xbsize,ybsize;
     int total_local_size;
     TRANSPOSE_EL_TYPE *data;
     transpose_mpi_plan tp;
     int el_size = 1, el;

     MPI_Init(&argc,&argv);
     MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
     MPI_Comm_size(MPI_COMM_WORLD,&n_pes);

     if (argc < 3) {
	  if (my_pe == 1)
	       printf("Syntax: test_transpose_mpi nx ny [el_size]\n");
	  MPI_Finalize();
	  return 1;
     }

     nx = atoi(argv[1]);
     ny = atoi(argv[2]);

     if (argc > 3)
	  el_size = atoi(argv[3]);

     transpose_mpi_get_local_size(nx,my_pe,n_pes,
				  &local_nx,&local_x_start);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,
				  &local_ny,&local_y_start);
     total_local_size = transpose_mpi_get_local_storage_size(nx,ny,
							     my_pe,n_pes);


     if (total_local_size == 0)
	  data = 0;
     else
	  data = (TRANSPOSE_EL_TYPE *)
	       malloc(total_local_size*el_size*sizeof(TRANSPOSE_EL_TYPE));

     tp = transpose_mpi_create_plan(nx,ny,MPI_COMM_WORLD);
     
     /* output plan data, one process at a time: */
     for (j = 0; j < n_pes; ++j) {
	  if (j == my_pe) {
	       printf("Plan for process %d:\n",j);
	       printf("   nx = %d, ny = %d\n",tp->nx,tp->ny);
	       printf("   local_nx = %d, local_ny = %d\n",
		      tp->local_nx,tp->local_ny);

	       if (local_nx > 0 || local_ny > 0) {
	       printf("   block_size = %d (padded = %d), num_blocks = %d\n",
		      tp->block_size, tp->block_size_padded, tp->num_blocks);
	       for (i = 0; i < tp->num_blocks; ++i)
		    printf("      block_dest[%d]: block = %d, dest_pe = %d, "
			   "send_size = %d, recv_size = %d\n",
			   i,tp->block_dest[i].block_num,
			   tp->block_dest[i].dest_pe,
			   tp->block_dest[i].send_size,
			   tp->block_dest[i].recv_size);
	       printf("   perm_block_size = %d, num_perm_blocks = %d\n",
		      tp->perm_block_size, tp->num_perm_blocks);
	       #if OUTPUT_PERMUTATION
	       for (i = 0; i < tp->num_perm_blocks; ++i)
		    printf("      perm_block_dest[%d] = %d\n",
			   i,tp->perm_block_dest[i]);
	       #endif
	       }

	       printf("\n");
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
     }
	      
     for (i = 0; i < local_nx; ++i)
	  for (j = 0; j < ny; ++j)
	       for (el = 0; el < el_size; ++el)
		    data[el+el_size*(i*ny+j)] = ny * (i + local_x_start) + j;

     if (my_pe == 0)
	  printf("\nComputing transpose of %dx%d matrix (%d elements)...\n",
		 nx,ny,el_size);

     transpose_mpi(tp,data,el_size);

     MPI_Barrier(MPI_COMM_WORLD);
     
     printf("Checking result on process %d...\n",my_pe);

     for (i = 0; i < local_ny; ++i)
	  for (j = 0; j < nx; ++j) 
	       for (el = 0; el < el_size; ++el) {
		    if (data[el+el_size*(i*nx + j)] != 
			(TRANSPOSE_EL_TYPE) (j*ny + i+local_y_start)) {
			 printf("Error with x=%d, y=%d on process %d!\n"
				"  -- is %g rather than %g\n",
				j,i+local_y_start,my_pe,
				data[el+el_size*(i*nx + j)], 
				(TRANSPOSE_EL_TYPE) (j*ny + i+local_y_start));
			 MPI_Finalize();
			 exit(1);
	       }
	  }

     printf("Process %d okay!\n",my_pe);
     
     /* Initialize data and check a second time, to make sure that
	we can reuse plans.*/

     MPI_Barrier(MPI_COMM_WORLD);
     
     for (i = 0; i < local_nx; ++i)
	  for (j = 0; j < ny; ++j)
	       for (el = 0; el < el_size; ++el)
		    data[el+el_size*(i*ny+j)] = ny * (i + local_x_start) + j;

     if (my_pe == 0)
	  printf("\nTrying again, just to be sure...\n"
		 "Computing transpose of %dx%d matrix (%d elements)...\n",
		 nx,ny,el_size);

     transpose_mpi(tp,data,el_size);

     MPI_Barrier(MPI_COMM_WORLD);
     
     printf("Checking result on process %d...\n",my_pe);

     for (i = 0; i < local_ny; ++i)
	  for (j = 0; j < nx; ++j) 
	       for (el = 0; el < el_size; ++el) {
		    if (data[el+el_size*(i*nx + j)] != 
			(TRANSPOSE_EL_TYPE) (j*ny + i+local_y_start)) {
			 printf("Error with x=%d, y=%d on process %d!\n"
				"  -- is %g rather than %g\n",
				j,i+local_y_start,my_pe,
				data[el+el_size*(i*nx + j)], 
				(TRANSPOSE_EL_TYPE) (j*ny + i+local_y_start));
			 MPI_Finalize();
			 exit(1);
	       }
	  }

     transpose_mpi_destroy_plan(tp);

     if (data)
	  free(data);

     printf("Process %d okay!\n",my_pe);
     
     MPI_Finalize();

     return 0;
}

