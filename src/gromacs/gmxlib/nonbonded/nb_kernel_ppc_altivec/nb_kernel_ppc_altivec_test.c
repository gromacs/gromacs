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

#ifdef HAVE_ALTIVEC_H
#include <altivec.h>
#endif


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
void
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


