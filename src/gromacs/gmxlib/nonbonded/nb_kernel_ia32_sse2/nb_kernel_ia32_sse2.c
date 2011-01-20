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

#include "nb_kernel_ia32_sse2.h"

/* Include ia32 SSE kernel headers in local directory */
#include "nb_kernel010_ia32_sse2.h"
#include "nb_kernel030_ia32_sse2.h"
#include "nb_kernel100_ia32_sse2.h"
#include "nb_kernel101_ia32_sse2.h"
#include "nb_kernel102_ia32_sse2.h"
#include "nb_kernel103_ia32_sse2.h"
#include "nb_kernel104_ia32_sse2.h"
#include "nb_kernel110_ia32_sse2.h"
#include "nb_kernel111_ia32_sse2.h"
#include "nb_kernel112_ia32_sse2.h"
#include "nb_kernel113_ia32_sse2.h"
#include "nb_kernel114_ia32_sse2.h"
#include "nb_kernel130_ia32_sse2.h"
#include "nb_kernel131_ia32_sse2.h"
#include "nb_kernel132_ia32_sse2.h"
#include "nb_kernel133_ia32_sse2.h"
#include "nb_kernel134_ia32_sse2.h"
#include "nb_kernel200_ia32_sse2.h"
#include "nb_kernel201_ia32_sse2.h"
#include "nb_kernel202_ia32_sse2.h"
#include "nb_kernel203_ia32_sse2.h"
#include "nb_kernel204_ia32_sse2.h"
#include "nb_kernel210_ia32_sse2.h"
#include "nb_kernel211_ia32_sse2.h"
#include "nb_kernel212_ia32_sse2.h"
#include "nb_kernel213_ia32_sse2.h"
#include "nb_kernel214_ia32_sse2.h"
#include "nb_kernel230_ia32_sse2.h"
#include "nb_kernel231_ia32_sse2.h"
#include "nb_kernel232_ia32_sse2.h"
#include "nb_kernel233_ia32_sse2.h"
#include "nb_kernel234_ia32_sse2.h"
#include "nb_kernel300_ia32_sse2.h"
#include "nb_kernel301_ia32_sse2.h"
#include "nb_kernel302_ia32_sse2.h"
#include "nb_kernel303_ia32_sse2.h"
#include "nb_kernel304_ia32_sse2.h"
#include "nb_kernel310_ia32_sse2.h"
#include "nb_kernel311_ia32_sse2.h"
#include "nb_kernel312_ia32_sse2.h"
#include "nb_kernel313_ia32_sse2.h"
#include "nb_kernel314_ia32_sse2.h"
#include "nb_kernel330_ia32_sse2.h"
#include "nb_kernel331_ia32_sse2.h"
#include "nb_kernel332_ia32_sse2.h"
#include "nb_kernel333_ia32_sse2.h"
#include "nb_kernel334_ia32_sse2.h"
#include "nb_kernel400_ia32_sse2.h"
#include "nb_kernel410_ia32_sse2.h"
#include "nb_kernel430_ia32_sse2.h"


#include <stdlib.h>
#include <stdio.h>
/* Necessary headers for POSIX-style long jumps. */
#include <signal.h>
#include <setjmp.h>




#include "../nb_kerneltype.h"
#include "nb_kernel_ia32_sse2.h"
#include "nb_kernel_ia32_sse2_test_asm.h"

static nb_kernel_t *
kernellist_ia32_sse2[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_ia32_sse2,
    NULL,
    nb_kernel030_ia32_sse2,
    nb_kernel100_ia32_sse2,
    nb_kernel101_ia32_sse2,
    nb_kernel102_ia32_sse2,
    nb_kernel103_ia32_sse2,
    nb_kernel104_ia32_sse2,
    nb_kernel110_ia32_sse2,
    nb_kernel111_ia32_sse2,
    nb_kernel112_ia32_sse2,
    nb_kernel113_ia32_sse2,
    nb_kernel114_ia32_sse2,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel130_ia32_sse2,
    nb_kernel131_ia32_sse2,
    nb_kernel132_ia32_sse2,
    nb_kernel133_ia32_sse2,
    nb_kernel134_ia32_sse2,
    nb_kernel200_ia32_sse2,
    nb_kernel201_ia32_sse2,
    nb_kernel202_ia32_sse2,
    nb_kernel203_ia32_sse2,
    nb_kernel204_ia32_sse2,
    nb_kernel210_ia32_sse2,
    nb_kernel211_ia32_sse2,
    nb_kernel212_ia32_sse2,
    nb_kernel213_ia32_sse2,
    nb_kernel214_ia32_sse2,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel230_ia32_sse2,
    nb_kernel231_ia32_sse2,
    nb_kernel232_ia32_sse2,
    nb_kernel233_ia32_sse2,
    nb_kernel234_ia32_sse2,
    nb_kernel300_ia32_sse2,
    nb_kernel301_ia32_sse2,
    nb_kernel302_ia32_sse2,
    nb_kernel303_ia32_sse2,
    nb_kernel304_ia32_sse2,
    nb_kernel310_ia32_sse2,
    nb_kernel311_ia32_sse2,
    nb_kernel312_ia32_sse2,
    nb_kernel313_ia32_sse2,
    nb_kernel314_ia32_sse2,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel330_ia32_sse2,
    nb_kernel331_ia32_sse2,
    nb_kernel332_ia32_sse2,
    nb_kernel333_ia32_sse2,
    nb_kernel334_ia32_sse2,
    nb_kernel400_ia32_sse2,
    nb_kernel410_ia32_sse2,
    nb_kernel430_ia32_sse2
};

#ifdef GMX_THREAD_SHM_FDECOMP
static tMPI_Thread_mutex_t 
nb_kernel_ia32_sse2_test_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif

/*! Posix long jump label */
static jmp_buf         
nb_kernel_ia32_sse2_testprog;

/*! Result of ia32 SSE2 test */
static gmx_bool
nb_kernel_ia32_sse2_present;


static void 
nb_kernel_ia32_sse2_sigill_handler(int n)
{
  nb_kernel_ia32_sse2_present=FALSE;
  longjmp(nb_kernel_ia32_sse2_testprog,n);
}





/* Return 0 if SSE2 support is present, or
 * or non-zero on failure.
 */
int 
nb_kernel_ia32_sse2_test(FILE *                log)
{
	/* 
	 * This should NOT be called from threads, 
	 * but just in case you still try to do it...
	 */
#ifdef GMX_THREAD_SHM_FDECOMP
	tMPI_Thread_mutex_lock(&nb_kernel_ia32_sse2_test_mutex);
#endif

    if(log)
        fprintf(log,"Testing ia32 SSE2 support...");
    
	nb_kernel_ia32_sse2_present = TRUE;
	signal(SIGILL,nb_kernel_ia32_sse2_sigill_handler);

	/* return to this point after executing the signal handler
	 * if we catch a SIGILL
	 */
	setjmp(nb_kernel_ia32_sse2_testprog); 

	if(nb_kernel_ia32_sse2_present)
		nb_kernel_ia32_sse2_test_asm();
	
	/* If SSE2 worked, then success is still 1.
     * If we got SIGILL, it was set to 0 in sigill_handler().
     */

	if(log)
		fprintf(log," %spresent.\n", 
				nb_kernel_ia32_sse2_present ? "":"not ");
	
#ifdef GMX_THREAD_SHM_FDECOMP
	tMPI_Thread_mutex_unlock(&nb_kernel_ia32_sse2_test_mutex);
#endif
    
	return ((nb_kernel_ia32_sse2_present) ? 0 : -1);
}

				
			

void
nb_kernel_setup_ia32_sse2(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_ia32_sse2_test(log) != 0)
        return;
    
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_ia32_sse2[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    





