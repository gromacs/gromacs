/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
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

#include "nb_kernel_sse2_double.h"

/* Include double precision SSE intrinsics kernel headers in local directory */
#include "nb_kernel400_sse2_double.h"
#include "nb_kernel410_sse2_double.h"
#include "nb_kernel430_sse2_double.h"

#include <stdlib.h>
#include <stdio.h>

#include "../nb_kerneltype.h"
#include "nb_kernel_sse2_double.h"

static nb_kernel_t *
kernellist_sse2_double[eNR_NBKERNEL_NR] = 
{
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
    nb_kernel400_sse2_double,
    nb_kernel410_sse2_double,
    nb_kernel430_sse2_double
};


/* Return 0 if SSE support is present, or
 * non-zero on failure.
 */
int 
nb_kernel_sse2_double_test(FILE *                log)
{
	unsigned int level;
	unsigned int _eax,_ebx,_ecx,_edx;
	int status;
	int CPUInfo[4];
	
	if(NULL != log)
    {
		fprintf(log,"Checking CPU SSE2 support... ");
    }
	
	level = 1;
#ifdef _MSC_VER
	__cpuid(CPUInfo,1);
	
	_eax=CPUInfo[0];
	_ebx=CPUInfo[1];
	_ecx=CPUInfo[2];
	_edx=CPUInfo[3];
	
#elif defined(__x86_64__)
	/* GCC 64-bit inline asm */
	__asm__ ("push %%rbx\n\tcpuid\n\tpop %%rbx\n"                 \
			 : "=a" (_eax), "=S" (_ebx), "=c" (_ecx), "=d" (_edx) \
			 : "0" (level));
#elif defined(__i386__)
	__asm__ ("push %%ebx\n\tcpuid\n\tpop %%ebx\n"                 \
			 : "=a" (_eax), "=S" (_ebx), "=c" (_ecx), "=d" (_edx) \
			 : "0" (level));
#else
	if(NULL != log)
	{
		fprintf(log,"Don't know how to call cpuid() on this system!\n");
	}
	_eax=_ebx=_ecx=_edx=0;
#endif
	
	/* Features:                                                                                                       
	 *                                                                                                                 
	 * SSE      Bit 25 of edx should be set                                                                            
	 * SSE2     Bit 26 of edx should be set                                                                            
	 * SSE3     Bit  0 of ecx should be set                                                                            
	 * SSE4.1   Bit 19 of ecx should be set                                                                            
	 */
	status =  (_edx & (1 << 26)) != 0;
	
	if(NULL != log)
	{
		fprintf(log,"%s present.", (status==0) ? "not" : "");
	}
	
	/* Return SSE2 status */
	return status;
}




void
nb_kernel_setup_sse2_double(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_sse2_double_test(log) == 0)
    {
		return;
    }
	
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_sse2_double[i];
        if(p!=NULL)
		{
			list[i] = p; 
		}
    }
}    


	

	
