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

#include <stdio.h>
#include <stdlib.h>

#include "sched.h"

int main(int argc, char **argv)
{
     int **sched;
     int npes = -1, mype, sortpe = -1, steps;

     if (argc >= 2) {
	  npes = atoi(argv[1]);
	  if (npes <= 0) {
	       fprintf(stderr,"npes must be positive!");
	       return 1;
	  }
     }
     if (argc >= 3) {
	  sortpe = atoi(argv[2]);
	  if (sortpe < 0 || sortpe >= npes) {
	       fprintf(stderr,"sortpe must be between 0 and npes-1.\n");
	       return 1;
	  }
     }

     if (npes != -1) {
	  printf("Computing schedule for npes = %d:\n",npes);
	  sched = make_comm_schedule(npes);
	  if (!sched) {
	       fprintf(stderr,"Out of memory!");
	       return 6;
	  }
	  
	  if (steps = check_comm_schedule(sched,npes))
	       printf("schedule OK (takes %d steps to complete).\n", steps);
	  else
	       printf("schedule not OK.\n");

	  print_comm_schedule(sched, npes);
	  
	  if (sortpe != -1) {
	       printf("\nSorting schedule for sortpe = %d...\n", sortpe);
	       sort_comm_schedule(sched,npes,sortpe);
	       
	       if (steps = check_comm_schedule(sched,npes))
		    printf("schedule OK (takes %d steps to complete).\n", 
			   steps);
	       else
		    printf("schedule not OK.\n");

	       print_comm_schedule(sched, npes);
	       
	       printf("\nInverting schedule...\n");
	       invert_comm_schedule(sched,npes);
	       
	       if (steps = check_comm_schedule(sched,npes))
		    printf("schedule OK (takes %d steps to complete).\n", 
			   steps);
	       else
		    printf("schedule not OK.\n");

	       print_comm_schedule(sched, npes);
	       
	       free_comm_schedule(sched,npes);
	  }
     }
     else {
	  printf("Doing infinite tests...\n");
	  for (npes = 1; ; ++npes) {
	       printf("npes = %d...",npes);
	       sched = make_comm_schedule(npes);
	       if (!sched) {
		    fprintf(stderr,"Out of memory!\n");
		    return 5;
	       }
	       for (sortpe = 0; sortpe < npes; ++sortpe) {
		    empty_comm_schedule(sched,npes);
		    fill_comm_schedule(sched,npes);
		    if (!check_comm_schedule(sched,npes)) {
			 fprintf(stderr,
				 "\n -- fill error for sortpe = %d!\n",sortpe);
			 return 2;
		    }
		    sort_comm_schedule(sched,npes,sortpe);
		    if (!check_comm_schedule(sched,npes)) {
			 fprintf(stderr,
				 "\n -- sort error for sortpe = %d!\n",sortpe);
			 return 3;
		    }
		    invert_comm_schedule(sched,npes);
		    if (!check_comm_schedule(sched,npes)) {
			 fprintf(stderr,
				 "\n -- invert error for sortpe = %d!\n",
				 sortpe);
			 return 4;
		    }
	       }
	       free_comm_schedule(sched,npes);
	       printf("OK\n");
	       if (npes % 50 == 0)
		    printf("(...Hit Ctrl-C to stop...)\n");
	  }
     }

     return 0;
}
