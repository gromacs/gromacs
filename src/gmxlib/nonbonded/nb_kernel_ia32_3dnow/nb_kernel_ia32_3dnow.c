/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- 
 *
 * $Id$
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
#include <gmx_thread.h>

#include <types/simple.h>
#include <types/nrnb.h>

/* Include ia32_3dnow kernel headers in local directory */
#include "nb_kernel010_ia32_3dnow.h"
#include "nb_kernel030_ia32_3dnow.h"
#include "nb_kernel100_ia32_3dnow.h"
#include "nb_kernel101_ia32_3dnow.h"
#include "nb_kernel102_ia32_3dnow.h"
#include "nb_kernel110_ia32_3dnow.h"
#include "nb_kernel111_ia32_3dnow.h"
#include "nb_kernel112_ia32_3dnow.h"
#include "nb_kernel300_ia32_3dnow.h"
#include "nb_kernel301_ia32_3dnow.h"
#include "nb_kernel302_ia32_3dnow.h"
#include "nb_kernel310_ia32_3dnow.h"
#include "nb_kernel311_ia32_3dnow.h"
#include "nb_kernel312_ia32_3dnow.h"
#include "nb_kernel330_ia32_3dnow.h"
#include "nb_kernel331_ia32_3dnow.h"
#include "nb_kernel332_ia32_3dnow.h"



#include <stdlib.h>
#include <stdio.h>
/* Necessary headers for POSIX-style long jumps. */
#include <signal.h>
#include <setjmp.h>

#include "../nb_kerneltype.h"
#include "nb_kernel_ia32_3dnow.h"
#include "nb_kernel_ia32_3dnow_test_asm.h"



static nb_kernel_t *
kernellist_ia32_3dnow[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_ia32_3dnow,
    NULL,
    nb_kernel030_ia32_3dnow,
    nb_kernel100_ia32_3dnow,
    nb_kernel101_ia32_3dnow,
    nb_kernel102_ia32_3dnow,
    NULL,
    NULL,
    nb_kernel110_ia32_3dnow,
    nb_kernel111_ia32_3dnow,
    nb_kernel112_ia32_3dnow,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel300_ia32_3dnow,
    nb_kernel301_ia32_3dnow,
    nb_kernel302_ia32_3dnow,
    NULL,
    NULL,
    nb_kernel310_ia32_3dnow,
    nb_kernel311_ia32_3dnow,
    nb_kernel312_ia32_3dnow,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel330_ia32_3dnow,
    nb_kernel331_ia32_3dnow,
    nb_kernel332_ia32_3dnow,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
};



#ifdef GMX_THREADS
static gmx_thread_mutex_t 
nb_kernel_ia32_3dnow_test_mutex = GMX_THREAD_MUTEX_INITIALIZER;
#endif

/*! Posix long jump label */
static jmp_buf         
nb_kernel_ia32_3dnow_testprog;

/*! Result of ia32 3DNow test */
static bool   
nb_kernel_ia32_3dnow_present;


static void 
nb_kernel_ia32_3dnow_sigill_handler(int n)
{
  nb_kernel_ia32_3dnow_present=FALSE;
  longjmp(nb_kernel_ia32_3dnow_testprog,n);
}





/* Return 0 if 3DNow support is present, or
 * non-zero upon failure.
 */
int 
nb_kernel_ia32_3dnow_test(FILE *                log)
{
	/* 
	 * This should NOT be called from threads, 
	 * but just in case you still try to do it...
	 */
#ifdef GMX_THREADS
	gmx_thread_mutex_lock(&nb_kernel_ia32_3dnow_test_mutex);
#endif
    
    if(log)
        fprintf(log,"Testing AMD 3DNow support...");

	nb_kernel_ia32_3dnow_present = TRUE;
	signal(SIGILL,nb_kernel_ia32_3dnow_sigill_handler);

	/* return to this point after executing the signal handler
	 * if we catch a SIGILL
	 */
	setjmp(nb_kernel_ia32_3dnow_testprog); 

	if(nb_kernel_ia32_3dnow_present)
		nb_kernel_ia32_3dnow_test_asm();
	
	/* If 3DNow worked, then success is still 1.
     * If we got SIGILL, it was set to 0 in sigill_handler().
     */

	if(log) 
		fprintf(log," %spresent.\n", 
				nb_kernel_ia32_3dnow_present ? "":"not ");
	
#ifdef GMX_THREADS
	pthread_mutex_unlock(&nb_kernel_ia32_3dnow_test_mutex);
#endif
	
	return ((nb_kernel_ia32_3dnow_present) ? 0 : -1);
}

				


void
nb_kernel_setup_ia32_3dnow(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_ia32_3dnow_test(log) != 0)
        return;
    
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_ia32_3dnow[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    

