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


#include <types/nrnb.h>
#include "types/simple.h"

#include "../nb_kerneltype.h"
#include "nb_kernel_ppc_altivec.h"
#include "nb_kernel_ppc_altivec_test.h"

/* Include altivec kernel headers in local directory */
#include "nb_kernel010_ppc_altivec.h"
#include "nb_kernel030_ppc_altivec.h"
#include "nb_kernel100_ppc_altivec.h"
#include "nb_kernel101_ppc_altivec.h"
#include "nb_kernel102_ppc_altivec.h"
#include "nb_kernel103_ppc_altivec.h"
#include "nb_kernel104_ppc_altivec.h"
#include "nb_kernel110_ppc_altivec.h"
#include "nb_kernel111_ppc_altivec.h"
#include "nb_kernel112_ppc_altivec.h"
#include "nb_kernel113_ppc_altivec.h"
#include "nb_kernel114_ppc_altivec.h"
#include "nb_kernel130_ppc_altivec.h"
#include "nb_kernel131_ppc_altivec.h"
#include "nb_kernel132_ppc_altivec.h"
#include "nb_kernel133_ppc_altivec.h"
#include "nb_kernel134_ppc_altivec.h"
#include "nb_kernel200_ppc_altivec.h"
#include "nb_kernel201_ppc_altivec.h"
#include "nb_kernel202_ppc_altivec.h"
#include "nb_kernel203_ppc_altivec.h"
#include "nb_kernel204_ppc_altivec.h"
#include "nb_kernel210_ppc_altivec.h"
#include "nb_kernel211_ppc_altivec.h"
#include "nb_kernel212_ppc_altivec.h"
#include "nb_kernel213_ppc_altivec.h"
#include "nb_kernel214_ppc_altivec.h"
#include "nb_kernel230_ppc_altivec.h"
#include "nb_kernel231_ppc_altivec.h"
#include "nb_kernel232_ppc_altivec.h"
#include "nb_kernel233_ppc_altivec.h"
#include "nb_kernel234_ppc_altivec.h"
#include "nb_kernel300_ppc_altivec.h"
#include "nb_kernel301_ppc_altivec.h"
#include "nb_kernel302_ppc_altivec.h"
#include "nb_kernel303_ppc_altivec.h"
#include "nb_kernel304_ppc_altivec.h"
#include "nb_kernel310_ppc_altivec.h"
#include "nb_kernel311_ppc_altivec.h"
#include "nb_kernel312_ppc_altivec.h"
#include "nb_kernel313_ppc_altivec.h"
#include "nb_kernel314_ppc_altivec.h"
#include "nb_kernel330_ppc_altivec.h"
#include "nb_kernel331_ppc_altivec.h"
#include "nb_kernel332_ppc_altivec.h"
#include "nb_kernel333_ppc_altivec.h"
#include "nb_kernel334_ppc_altivec.h"
#include "nb_kernel400_ppc_altivec.h"
#include "nb_kernel410_ppc_altivec.h"
#include "nb_kernel430_ppc_altivec.h"

/* Necessary headers for POSIX-style long jumps. */
#include <signal.h>
#include <setjmp.h>


static nb_kernel_t *
kernellist_ppc_altivec[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_ppc_altivec,
    NULL,
    nb_kernel030_ppc_altivec,
    nb_kernel100_ppc_altivec,
    nb_kernel101_ppc_altivec,
    nb_kernel102_ppc_altivec,
    nb_kernel103_ppc_altivec,
    nb_kernel104_ppc_altivec,
    nb_kernel110_ppc_altivec,
    nb_kernel111_ppc_altivec,
    nb_kernel112_ppc_altivec,
    nb_kernel113_ppc_altivec,
    nb_kernel114_ppc_altivec,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel130_ppc_altivec,
    nb_kernel131_ppc_altivec,
    nb_kernel132_ppc_altivec,
    nb_kernel133_ppc_altivec,
    nb_kernel134_ppc_altivec,
    nb_kernel200_ppc_altivec,
    nb_kernel201_ppc_altivec,
    nb_kernel202_ppc_altivec,
    nb_kernel203_ppc_altivec,
    nb_kernel204_ppc_altivec,
    nb_kernel210_ppc_altivec,
    nb_kernel211_ppc_altivec,
    nb_kernel212_ppc_altivec,
    nb_kernel213_ppc_altivec,
    nb_kernel214_ppc_altivec,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel230_ppc_altivec,
    nb_kernel231_ppc_altivec,
    nb_kernel232_ppc_altivec,
    nb_kernel233_ppc_altivec,
    nb_kernel234_ppc_altivec,
    nb_kernel300_ppc_altivec,
    nb_kernel301_ppc_altivec,
    nb_kernel302_ppc_altivec,
    nb_kernel303_ppc_altivec,
    nb_kernel304_ppc_altivec,
    nb_kernel310_ppc_altivec,
    nb_kernel311_ppc_altivec,
    nb_kernel312_ppc_altivec,
    nb_kernel313_ppc_altivec,
    nb_kernel314_ppc_altivec,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel330_ppc_altivec,
    nb_kernel331_ppc_altivec,
    nb_kernel332_ppc_altivec,
    nb_kernel333_ppc_altivec,
    nb_kernel334_ppc_altivec,
    nb_kernel400_ppc_altivec,
    nb_kernel410_ppc_altivec,
    nb_kernel430_ppc_altivec
};


static jmp_buf         
nb_kernel_ppc_altivec_testprog;


/*! Result of Altivec test */
static gmx_bool
nb_kernel_ppc_altivec_present;


static void 
nb_kernel_ppc_altivec_sigill_handler(int n)
{
    nb_kernel_ppc_altivec_present=FALSE;
    longjmp(nb_kernel_ppc_altivec_testprog,n);
}



/* Return 0 if Altivec support is present, or
 *  non-zero upon error.
 */
int 
nb_kernel_ppc_altivec_test(FILE * log)
{
    
    if(log)
        fprintf(log,"Testing Altivec/VMX support...");

    nb_kernel_ppc_altivec_present = TRUE;
    signal(SIGILL,nb_kernel_ppc_altivec_sigill_handler);
    
    /* return to this point after executing the signal handler
     * if we catch a SIGILL
     */
    setjmp(nb_kernel_ppc_altivec_testprog); 
    
    if(nb_kernel_ppc_altivec_present)
        nb_kernel_ppc_altivec_issue_instructions();
    
    /* If altivec worked, then success is still 1.
        * If we got SIGILL, it was set to 0 in sigill_handler().
        */
    
    if(log)
        fprintf(log," %spresent.\n", 
                nb_kernel_ppc_altivec_present ? "":"not ");
    
    
    return ((nb_kernel_ppc_altivec_present) ? 0 : -1);
}





void
nb_kernel_setup_ppc_altivec(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_ppc_altivec_test(log) != 0)
        return;
    
	if(log)
		fprintf(log,"Configuring PPC/Altivec nonbonded kernels...\n");

    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_ppc_altivec[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    





