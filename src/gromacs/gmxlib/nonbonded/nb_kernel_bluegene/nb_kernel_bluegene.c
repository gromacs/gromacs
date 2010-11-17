/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Must come directly after config.h */
#ifdef GMX_THREAD_SHM_FDECOMP
#include <thread_mpi.h>
#endif

#include <types/simple.h>
#include <types/nrnb.h>

#include "nb_kernel_bluegene.h"

/* Include Bluegene kernel headers in local directory */
#include "nb_kernel010_bluegene.h"
#include "nb_kernel020_bluegene.h"
#include "nb_kernel030_bluegene.h"
#include "nb_kernel100_bluegene.h"
#include "nb_kernel101_bluegene.h"
#include "nb_kernel102_bluegene.h"
#include "nb_kernel103_bluegene.h"
#include "nb_kernel104_bluegene.h"
#include "nb_kernel110_bluegene.h"
#include "nb_kernel111_bluegene.h"
#include "nb_kernel112_bluegene.h"
#include "nb_kernel113_bluegene.h"
#include "nb_kernel114_bluegene.h"
#include "nb_kernel120_bluegene.h"
#include "nb_kernel121_bluegene.h"
#include "nb_kernel122_bluegene.h"
#include "nb_kernel123_bluegene.h"
#include "nb_kernel124_bluegene.h"
#include "nb_kernel130_bluegene.h"
#include "nb_kernel131_bluegene.h"
#include "nb_kernel132_bluegene.h"
#include "nb_kernel133_bluegene.h"
#include "nb_kernel134_bluegene.h"
#include "nb_kernel200_bluegene.h"
#include "nb_kernel201_bluegene.h"
#include "nb_kernel202_bluegene.h"
#include "nb_kernel203_bluegene.h"
#include "nb_kernel204_bluegene.h"
#include "nb_kernel210_bluegene.h"
#include "nb_kernel211_bluegene.h"
#include "nb_kernel212_bluegene.h"
#include "nb_kernel213_bluegene.h"
#include "nb_kernel214_bluegene.h"
#include "nb_kernel220_bluegene.h"
#include "nb_kernel221_bluegene.h"
#include "nb_kernel222_bluegene.h"
#include "nb_kernel223_bluegene.h"
#include "nb_kernel224_bluegene.h"
#include "nb_kernel230_bluegene.h"
#include "nb_kernel231_bluegene.h"
#include "nb_kernel232_bluegene.h"
#include "nb_kernel233_bluegene.h"
#include "nb_kernel234_bluegene.h"
#include "nb_kernel300_bluegene.h"
#include "nb_kernel301_bluegene.h"
#include "nb_kernel302_bluegene.h"
#include "nb_kernel303_bluegene.h"
#include "nb_kernel304_bluegene.h"
#include "nb_kernel310_bluegene.h"
#include "nb_kernel311_bluegene.h"
#include "nb_kernel312_bluegene.h"
#include "nb_kernel313_bluegene.h"
#include "nb_kernel314_bluegene.h"
#include "nb_kernel320_bluegene.h"
#include "nb_kernel321_bluegene.h"
#include "nb_kernel322_bluegene.h"
#include "nb_kernel323_bluegene.h"
#include "nb_kernel324_bluegene.h"
#include "nb_kernel330_bluegene.h"
#include "nb_kernel331_bluegene.h"
#include "nb_kernel332_bluegene.h"
#include "nb_kernel333_bluegene.h"
#include "nb_kernel334_bluegene.h"
#include "nb_kernel400_bluegene.h"
#include "nb_kernel410_bluegene.h"
#include "nb_kernel420_bluegene.h"
#include "nb_kernel430_bluegene.h"



#include <stdlib.h>
#include <stdio.h>
/* Necessary headers for POSIX-style long jumps. */
#include <signal.h>
#include <setjmp.h>



#include "../nb_kerneltype.h"
#include "nb_kernel_bluegene.h"


static nb_kernel_t *
kernellist_bluegene[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_bluegene,
    nb_kernel020_bluegene,
    nb_kernel030_bluegene,
    nb_kernel100_bluegene,
    nb_kernel101_bluegene,
    nb_kernel102_bluegene,
    nb_kernel103_bluegene,
    nb_kernel104_bluegene,
    nb_kernel110_bluegene,
    nb_kernel111_bluegene,
    nb_kernel112_bluegene,
    nb_kernel113_bluegene,
    nb_kernel114_bluegene,
    nb_kernel120_bluegene,
    nb_kernel121_bluegene,
    nb_kernel122_bluegene,
    nb_kernel123_bluegene,
    nb_kernel124_bluegene,
    nb_kernel130_bluegene,
    nb_kernel131_bluegene,
    nb_kernel132_bluegene,
    nb_kernel133_bluegene,
    nb_kernel134_bluegene,
    nb_kernel200_bluegene,
    nb_kernel201_bluegene,
    nb_kernel202_bluegene,
    nb_kernel203_bluegene,
    nb_kernel204_bluegene,
    nb_kernel210_bluegene,
    nb_kernel211_bluegene,
    nb_kernel212_bluegene,
    nb_kernel213_bluegene,
    nb_kernel214_bluegene,
    nb_kernel220_bluegene,
    nb_kernel221_bluegene,
    nb_kernel222_bluegene,
    nb_kernel223_bluegene,
    nb_kernel224_bluegene,
    nb_kernel230_bluegene,
    nb_kernel231_bluegene,
    nb_kernel232_bluegene,
    nb_kernel233_bluegene,
    nb_kernel234_bluegene,
    nb_kernel300_bluegene,
    nb_kernel301_bluegene,
    nb_kernel302_bluegene,
    nb_kernel303_bluegene,
    nb_kernel304_bluegene,
    nb_kernel310_bluegene,
    nb_kernel311_bluegene,
    nb_kernel312_bluegene,
    nb_kernel313_bluegene,
    nb_kernel314_bluegene,
    nb_kernel320_bluegene,
    nb_kernel321_bluegene,
    nb_kernel322_bluegene,
    nb_kernel323_bluegene,
    nb_kernel324_bluegene,
    nb_kernel330_bluegene,
    nb_kernel331_bluegene,
    nb_kernel332_bluegene,
    nb_kernel333_bluegene,
    nb_kernel334_bluegene,
    nb_kernel400_bluegene,
    nb_kernel410_bluegene,
    nb_kernel430_bluegene
};


#ifdef GMX_THREAD_SHM_FDECOMP
static tMPI_Thread_mutex_t 
nb_kernel_bluegene_test_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif


/*! Posix long jump label */
static jmp_buf         
nb_kernel_bluegene_testprog;

/*! Result of bluegene test */
static gmx_bool      
nb_kernel_bluegene_present;


static void 
nb_kernel_bluegene_sigill_handler(int n)
{
  nb_kernel_bluegene_present=FALSE;
  longjmp(nb_kernel_bluegene_testprog,n);
}



/* Return 0 if Bluegene support is present, or
 * non-zero on failure.
 */
int 
nb_kernel_bluegene_test(FILE *                log)
{
	/* 
	 * This should NOT be called from threads, 
	 * but just in case you still try to do it...
	 */
#ifdef GMX_THREAD_SHM_FDECOMP
	tMPI_Thread_mutex_lock(&nb_kernel_bluegene_test_mutex);
#endif
    
    if(log)
        fprintf(log,"Testing BlueGene support...");

	nb_kernel_bluegene_present = TRUE;
	signal(SIGILL,nb_kernel_bluegene_sigill_handler);

	/* return to this point after executing the signal handler
	 * if we catch a SIGILL
	 */
	setjmp(nb_kernel_bluegene_testprog); 

	if(nb_kernel_bluegene_present)
    {
        double _Complex a,b;
        double c = 2.0;
        
        a = __cmplx(c,c);
        b = __fpmadd(a,a,a);
    }
	
	/* If Bluegene worked, then success is still 1.
     * If we got SIGILL, it was set to 0 in sigill_handler().
     */

	if(log)
		fprintf(log," %spresent.\n", 
				nb_kernel_bluegene_present ? "":"not ");
	
#ifdef GMX_THREAD_SHM_FDECOMP
	tMPI_Thread_mutex_unlock(&nb_kernel_bluegene_test_mutex);
#endif
    
	return ((nb_kernel_bluegene_present) ? 0 : -1);
}

				
			

void
nb_kernel_setup_bluegene(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_bluegene_test(log) != 0)
        return;
    
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_bluegene[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    






	
	

	
