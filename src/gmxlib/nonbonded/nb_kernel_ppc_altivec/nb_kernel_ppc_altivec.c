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

#include <types/nrnb.h>

#include "../nb_kerneltype.h"

#include "nb_kernel_ppc_altivec.h"

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

#ifdef GMX_THREADS
static gmx_thread_mutex_t 
nb_kernel_ppc_altivec_test_mutex = GMX_THREAD_MUTEX_INITIALIZER;
#endif

/*! Posix long jump label */
static jmp_buf         
nb_kernel_ppc_altivec_testprog;

/*! Result of Altivec test */
static bool      
nb_kernel_ppc_altivec_present;


static void 
nb_kernel_ppc_altivec_sigill_handler(int n)
{
    nb_kernel_ppc_altivec_present=FALSE;
    longjmp(nb_kernel_ppc_altivec_testprog,n);
}




/*! \brief Issue Altivec instructions - can cause SIGILL!
*
*  Utility function for nb_kernel_ppc_altivec_test.
*  This code issues a couple of altivec instructions; if Altivec
*  support is not present it will generate an illegal instruction
*  signal which we should capture in nb_kernel_ppc_altivec_sigill_handler().
*
*  Since we need to do something Altivec-related in it anyway, we try to 
*  set the vector unit to non-java mode which is slightly faster.
* 
*  \warning     This HAS to be a separate function, since any function
*               with Altivec instructions in it will cause stack-related
*               Altivec instructions to be issued, even if we never
*               access the "real" Altivec instructions.
*
*  \threadsafe  No. The code itself is clean, but on many operating 
*               systems the signal handling does not work for threads.
*              
*/ 
static void
nb_kernel_ppc_altivec_issue_instructions(void)
{
    vector unsigned short vsr1,vsr2;
    vector unsigned int tmp1,tmp2;
    
    vsr1=vec_mfvscr();
    tmp1=vec_splat_u32(1);
    tmp2=vec_splat_u32(8);
    tmp1=vec_sl(tmp1,tmp2);
    vsr2=(vector unsigned short)vec_sl(tmp1,tmp2);
    vsr1=vec_or(vsr1,vsr2);
    vec_mtvscr(vsr1);
}



/* Return 0 if Altivec support is present, or
 *  non-zero upon error.
 */
int 
nb_kernel_ppc_altivec_test(FILE * log)
{
    /* 
     * This should NOT be called from threads, 
     * but just in case you still try to do it...
     */
#ifdef GMX_THREADS
    gmx_thread_mutex_lock(&nb_kernel_ppc_altivec_test_mutex);
#endif
    
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
    
#ifdef GMX_THREADS
    gmx_thread_mutex_unlock(&nb_kernel_ppc_altivec_test_mutex);
#endif
    
    return ((nb_kernel_ppc_altivec_present) ? 0 : -1);
}





void
nb_kernel_setup_ppc_altivec(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
    
    if(nb_kernel_ppc_altivec_test(log) != 0)
        return;
    
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_ppc_altivec[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    





