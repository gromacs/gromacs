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
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include <fftw_mpi.h>

#define NX 14
#define NY 10

void perform_transpose_async(transpose_mpi_plan tp, int el_size,
			     TRANSPOSE_EL_TYPE *data,
			     int local_nx, int local_x_start, int ny)
{
     int i,j,el;
     TRANSPOSE_EL_TYPE *send_buf;
     int step;

     for (i = 0; i < local_nx * ny * el_size; ++i)
	  data[i] = 0.0;

     send_buf = transpose_allocate_send_buf(tp,el_size);

     transpose_in_place_local(tp, el_size, data, BEFORE_TRANSPOSE);

     for (step = 0; step < tp->num_steps; ++step) {
	  int block_y_start, block_ny;

	  /* initialize block for step here (overlap initialization
	     and communication): */
	  transpose_get_send_block(tp, step,
				   &block_y_start, &block_ny);
	  for (i = 0; i < local_nx; ++i)
	       for (j = block_y_start; j < block_y_start + block_ny; ++j)
		    for (el = 0; el < el_size; ++el)
			 data[el+el_size*(i+j*local_nx)]
			      = ny * (i + local_x_start) + j;

	  transpose_finish_exchange_step(tp, step - 1);
	  transpose_start_exchange_step(tp,el_size,data,send_buf,step,
					TRANSPOSE_ASYNC);
     }

     transpose_finish_exchange_step(tp, step - 1);

     transpose_in_place_local(tp, el_size, data, AFTER_TRANSPOSE);

     if (send_buf)
	  fftw_free(send_buf);
}

void test_transpose_plan(int nx, int ny, int el_size, int verbose,
			 transpose_mpi_plan tp,
			 TRANSPOSE_EL_TYPE *data,
			 TRANSPOSE_EL_TYPE *work,
			 int async)
{
     int i,j,el;
     int my_pe,n_pes;
     int local_x_start,local_nx,local_y_start,local_ny, xbsize,ybsize;

     MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
     MPI_Comm_size(MPI_COMM_WORLD,&n_pes);

     transpose_mpi_get_local_size(nx,my_pe,n_pes,
                                  &local_nx,&local_x_start);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,
                                  &local_ny,&local_y_start);

     if (async) {
	  perform_transpose_async(tp,el_size,data,local_nx,local_x_start,ny);
     }
     else {
	  for (i = 0; i < local_nx; ++i)
	       for (j = 0; j < ny; ++j)
		    for (el = 0; el < el_size; ++el)
			 data[el+el_size*(i*ny+j)] 
			      = ny * (i + local_x_start) + j;
	  
	  transpose_mpi(tp,el_size,data,work);
     }

     /* print out result matrix if it is small: */
     if (verbose && local_ny < 16 && nx < 16) {
	  int pe;

	  for (pe = 0; pe < n_pes; ++pe) {
	       MPI_Barrier(MPI_COMM_WORLD);
	       if (pe == my_pe) {
		    printf("\nprocess %d result: \n", my_pe);
		    for (j = 0; j < nx; ++j) {
			 for (i = 0; i < local_ny; ++i)
			      printf("%4.0f",data[el_size*(i*nx + j)]);
			 printf("\n");
		    }
	       }
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
     }

     if (verbose)
       printf("Checking result on process %d...\n",my_pe);

     for (i = 0; i < local_ny; ++i)
	  for (j = 0; j < nx; ++j) 
	       for (el = 0; el < el_size; ++el) {
		    if (data[el+el_size*(i*nx + j)] != 
			(TRANSPOSE_EL_TYPE) (j*ny + i+local_y_start)) {
			 fprintf(stderr,
				 "Error with x=%d, y=%d on process %d!\n"
				 "  -- is %g rather than %g\n",
				 j,i+local_y_start,my_pe,
				 data[el+el_size*(i*nx + j)], 
				 (TRANSPOSE_EL_TYPE) (j*ny + i+local_y_start));
			 fftw_die("incorrect result.\n");
	       }
	  }

     if (verbose)
       printf("Process %d okay!\n",my_pe);
}

void test_transpose(int nx, int ny, int el_size, int verbose)
{
     int my_pe,n_pes;
     int i,j;
     int local_x_start,local_nx,local_y_start,local_ny, xbsize,ybsize;
     int total_local_size;
     TRANSPOSE_EL_TYPE *data, *work;
     transpose_mpi_plan tp;
     int el;

     MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
     MPI_Comm_size(MPI_COMM_WORLD,&n_pes);

     if (my_pe == 0)
	  printf("\nTesting transpose of %dx%d matrix (%d elements)...\n",
		 nx,ny,el_size);

     transpose_mpi_get_local_size(nx,my_pe,n_pes,
				  &local_nx,&local_x_start);
     transpose_mpi_get_local_size(ny,my_pe,n_pes,
				  &local_ny,&local_y_start);
     total_local_size = transpose_mpi_get_local_storage_size(nx,ny,
							     my_pe,n_pes);

     if (total_local_size == 0)
	  data = work = 0;
     else {
	  data = (TRANSPOSE_EL_TYPE *)
	       fftw_malloc(total_local_size*el_size*sizeof(TRANSPOSE_EL_TYPE));
	  work = (TRANSPOSE_EL_TYPE *)
	       fftw_malloc(total_local_size*el_size*sizeof(TRANSPOSE_EL_TYPE));
     }

     tp = transpose_mpi_create_plan(nx,ny,MPI_COMM_WORLD);
     
     /* output plan data, one process at a time: */
     if (verbose)
     for (j = 0; j < n_pes; ++j) {
	  if (j == my_pe) {
	       printf("Plan for process %d:\n",j);
	       printf("   nx = %d, ny = %d\n",tp->nx,tp->ny);
	       printf("   local_nx = %d, local_ny = %d\n",
		      tp->local_nx,tp->local_ny);

	       if (local_nx > 0 || local_ny > 0) {
		      printf("  send/recv block_size = %d/%d, "
			     "num_steps = %d\n",
			     tp->send_block_size, tp->recv_block_size, 
			     tp->num_steps);
		      for (i = 0; i < tp->num_steps; ++i)
			printf("      exchange[%d]: block = %d, dest_pe = %d, "
			   "send_size = %d, recv_size = %d\n",
			       i,tp->exchange[i].block_num,
			       tp->exchange[i].dest_pe,
			       tp->exchange[i].send_size,
			       tp->exchange[i].recv_size);
	       printf("   perm_block_size = %d, num_perm_blocks = %d\n",
		      tp->perm_block_size, tp->num_perm_blocks);
#if 0
	       for (i = 0; i < tp->num_perm_blocks; ++i)
		    printf("      perm_block_dest[%d] = %d\n",
			   i,tp->perm_block_dest[i]);
#endif
	       }

	       printf("\n");
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
     }
	      
     /* testing blocking, in-place transpose */
     if (verbose && my_pe == 0)
	  printf("\nTesting blocking in-place transpose:\n");
     test_transpose_plan(nx,ny,el_size,verbose,tp,data,NULL,0);
     MPI_Barrier(MPI_COMM_WORLD);

     /* Initialize data and check a second time, to make sure that
	we can reuse plans:*/     
     if (verbose && my_pe == 0) printf("\nTesting again:\n");
     test_transpose_plan(nx,ny,el_size,verbose,tp,data,NULL,0);     
     
     /* testing blocking, out-of-place transpose */
     if (verbose && my_pe == 0)
	  printf("\nTesting blocking out-of-place transpose:\n");
     test_transpose_plan(nx,ny,el_size,verbose,tp,data,work,0);
     MPI_Barrier(MPI_COMM_WORLD);
     
     /* Initialize data and check a second time, to make sure that
	we can reuse plans:*/     
     if (verbose && my_pe == 0) printf("\nTesting again:\n");
     test_transpose_plan(nx,ny,el_size,verbose,tp,data,work,0);     

     /* Test non-blocking (asynchronous) transpose calls: */
     if (verbose && my_pe == 0) printf("\nTesting non-blocking transpose:\n");

     test_transpose_plan(nx,ny,el_size,verbose,tp,data,NULL,1);
     MPI_Barrier(MPI_COMM_WORLD);
     
     /* Initialize data and check a second time, to make sure that
	we can reuse plans:*/     
     if (verbose && my_pe == 0) printf("\nTesting again:\n");
     test_transpose_plan(nx,ny,el_size,verbose,tp,data,NULL,1);

     transpose_mpi_destroy_plan(tp);
     
     if (data)
	  fftw_free(data);
     if (work)
	  fftw_free(work);

     if (verbose)
       printf("Process %d okay!\n",my_pe);

     fftw_check_memory_leaks();

     if (my_pe == 0)
       printf("okay!\n");
}

int main(int argc, char **argv)
{
     int my_pe;
     int nx = -1,ny = -1;
     int el_size = 1;
     int seed;

     MPI_Init(&argc,&argv);
     
     fftw_die_hook = fftw_mpi_die;

     seed = time(NULL);
     MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
     srand(seed);

     MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);

     /* only process 0 is guaranteed to have the correct args */
     if (my_pe == 0) {
	  if (argc >= 3) {
	       nx = atoi(argv[1]);
	       ny = atoi(argv[2]);
	  }

	  if (argc >= 4)
	       el_size = atoi(argv[3]);
     }
     /* broadcast args to other processes: */
     MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&el_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

     if (nx > 0 && ny > 0 && el_size > 0)
       test_transpose(nx,ny,el_size,1);
     else {
	  if (my_pe == 0)
	    printf("Doing random tests (does not terminate)...\n");
	  while (1) {
	       test_transpose(rand() % 10 + 1, rand() % 10 + 1,
			      rand() % 10 + 1, 0);
	       test_transpose(rand() % 100 + 1, rand() % 100 + 1,
			      rand() % 10 + 1, 0);
	       test_transpose(rand() % 1000 + 1, rand() % 1000 + 1,
			      rand() % 10 + 1, 0);
	  }
     }
     
     MPI_Finalize();

     return 0;
}

