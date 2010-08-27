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

#include "nb_kernel_x86_64_sse2.h"

/* Include x86_64 SSE kernel headers in local directory */
#include "nb_kernel010_x86_64_sse2.h"
#include "nb_kernel030_x86_64_sse2.h"
#include "nb_kernel100_x86_64_sse2.h"
#include "nb_kernel101_x86_64_sse2.h"
#include "nb_kernel102_x86_64_sse2.h"
#include "nb_kernel103_x86_64_sse2.h"
#include "nb_kernel104_x86_64_sse2.h"
#include "nb_kernel110_x86_64_sse2.h"
#include "nb_kernel111_x86_64_sse2.h"
#include "nb_kernel112_x86_64_sse2.h"
#include "nb_kernel113_x86_64_sse2.h"
#include "nb_kernel114_x86_64_sse2.h"
#include "nb_kernel130_x86_64_sse2.h"
#include "nb_kernel131_x86_64_sse2.h"
#include "nb_kernel132_x86_64_sse2.h"
#include "nb_kernel133_x86_64_sse2.h"
#include "nb_kernel134_x86_64_sse2.h"
#include "nb_kernel200_x86_64_sse2.h"
#include "nb_kernel201_x86_64_sse2.h"
#include "nb_kernel202_x86_64_sse2.h"
#include "nb_kernel203_x86_64_sse2.h"
#include "nb_kernel204_x86_64_sse2.h"
#include "nb_kernel210_x86_64_sse2.h"
#include "nb_kernel211_x86_64_sse2.h"
#include "nb_kernel212_x86_64_sse2.h"
#include "nb_kernel213_x86_64_sse2.h"
#include "nb_kernel214_x86_64_sse2.h"
#include "nb_kernel230_x86_64_sse2.h"
#include "nb_kernel231_x86_64_sse2.h"
#include "nb_kernel232_x86_64_sse2.h"
#include "nb_kernel233_x86_64_sse2.h"
#include "nb_kernel234_x86_64_sse2.h"
#include "nb_kernel300_x86_64_sse2.h"
#include "nb_kernel301_x86_64_sse2.h"
#include "nb_kernel302_x86_64_sse2.h"
#include "nb_kernel303_x86_64_sse2.h"
#include "nb_kernel304_x86_64_sse2.h"
#include "nb_kernel310_x86_64_sse2.h"
#include "nb_kernel311_x86_64_sse2.h"
#include "nb_kernel312_x86_64_sse2.h"
#include "nb_kernel313_x86_64_sse2.h"
#include "nb_kernel314_x86_64_sse2.h"
#include "nb_kernel330_x86_64_sse2.h"
#include "nb_kernel331_x86_64_sse2.h"
#include "nb_kernel332_x86_64_sse2.h"
#include "nb_kernel333_x86_64_sse2.h"
#include "nb_kernel334_x86_64_sse2.h"
#include "nb_kernel400_x86_64_sse2.h"
#include "nb_kernel410_x86_64_sse2.h"
#include "nb_kernel430_x86_64_sse2.h"



#include <stdlib.h>
#include <stdio.h>
/* Necessary headers for POSIX-style long jumps. */
#include <signal.h>
#include <setjmp.h>


#include "../nb_kerneltype.h"
#include "nb_kernel_x86_64_sse2.h"
#include "nb_kernel_x86_64_sse2_test_asm.h"


static nb_kernel_t *
kernellist_x86_64_sse2[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_x86_64_sse2,
    NULL,
    nb_kernel030_x86_64_sse2,
    nb_kernel100_x86_64_sse2,
    nb_kernel101_x86_64_sse2,
    nb_kernel102_x86_64_sse2,
    nb_kernel103_x86_64_sse2,
    nb_kernel104_x86_64_sse2,
    nb_kernel110_x86_64_sse2,
    nb_kernel111_x86_64_sse2,
    nb_kernel112_x86_64_sse2,
    nb_kernel113_x86_64_sse2,
    nb_kernel114_x86_64_sse2,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel130_x86_64_sse2,
    nb_kernel131_x86_64_sse2,
    nb_kernel132_x86_64_sse2,
    nb_kernel133_x86_64_sse2,
    nb_kernel134_x86_64_sse2,
    nb_kernel200_x86_64_sse2,
    nb_kernel201_x86_64_sse2,
    nb_kernel202_x86_64_sse2,
    nb_kernel203_x86_64_sse2,
    nb_kernel204_x86_64_sse2,
    nb_kernel210_x86_64_sse2,
    nb_kernel211_x86_64_sse2,
    nb_kernel212_x86_64_sse2,
    nb_kernel213_x86_64_sse2,
    nb_kernel214_x86_64_sse2,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel230_x86_64_sse2,
    nb_kernel231_x86_64_sse2,
    nb_kernel232_x86_64_sse2,
    nb_kernel233_x86_64_sse2,
    nb_kernel234_x86_64_sse2,
    nb_kernel300_x86_64_sse2,
    nb_kernel301_x86_64_sse2,
    nb_kernel302_x86_64_sse2,
    nb_kernel303_x86_64_sse2,
    nb_kernel304_x86_64_sse2,
    nb_kernel310_x86_64_sse2,
    nb_kernel311_x86_64_sse2,
    nb_kernel312_x86_64_sse2,
    nb_kernel313_x86_64_sse2,
    nb_kernel314_x86_64_sse2,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel330_x86_64_sse2,
    nb_kernel331_x86_64_sse2,
    nb_kernel332_x86_64_sse2,
    nb_kernel333_x86_64_sse2,
    nb_kernel334_x86_64_sse2,
    nb_kernel400_x86_64_sse2,
    nb_kernel410_x86_64_sse2,
    nb_kernel430_x86_64_sse2
};

#ifdef GMX_THREAD_SHM_FDECOMP
static tMPI_Thread_mutex_t 
nb_kernel_x86_64_sse2_test_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif

/*! Posix long jump label */
static jmp_buf         
nb_kernel_x86_64_sse2_testprog;

/*! Result of x86_64 SSE2 test */
static gmx_bool
nb_kernel_x86_64_sse2_present;


static void 
nb_kernel_x86_64_sse2_sigill_handler(int n)
{
  nb_kernel_x86_64_sse2_present=FALSE;
  longjmp(nb_kernel_x86_64_sse2_testprog,n);
}




/* Return GMX_SUCCESS (0) if SSE2 support is present, or
 * general error GMX_EFAILURE.
 */
int 
nb_kernel_x86_64_sse2_test(FILE *                log)
{
	/* 
	 * This should NOT be called from threads, 
	 * but just in case you still try to do it...
	 */
#ifdef GMX_THREAD_SHM_FDECOMP
	tMPI_Thread_mutex_lock(&nb_kernel_x86_64_sse2_test_mutex);
#endif
    
    if(log)
        fprintf(log,"Testing x86_64 SSE2 support...");

	nb_kernel_x86_64_sse2_present = TRUE;
	signal(SIGILL,nb_kernel_x86_64_sse2_sigill_handler);

	/* return to this point after executing the signal handler
	 * if we catch a SIGILL
	 */
	setjmp(nb_kernel_x86_64_sse2_testprog); 

	if(nb_kernel_x86_64_sse2_present)
		nb_kernel_x86_64_sse2_test_asm();
	
	/* If SSE2 worked, then success is still 1.
     * If we got SIGILL, it was set to 0 in sigill_handler().
     */

	if(log)
		fprintf(log," %spresent.\n", 
				nb_kernel_x86_64_sse2_present ? "":"not ");
	
#ifdef GMX_THREAD_SHM_FDECOMP
	tMPI_Thread_mutex_unlock(&nb_kernel_x86_64_sse2_test_mutex);
#endif
    
	return ((nb_kernel_x86_64_sse2_present) ? 0 : -1);
}

				

void
nb_kernel_setup_x86_64_sse2(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_x86_64_sse2_test(log) != 0)
        return;
    
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_x86_64_sse2[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    
