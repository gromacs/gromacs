/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_detectcpu_c = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <detectcpu.h>


#if (defined USE_X86_ASM || defined USE_PPC_ALTIVEC)
/* We use some special posix-style long jumps in
 * this routine to catch exceptions, but for portability
 * we only include them when actually needed.
 */
#include <signal.h>
#include <setjmp.h>

/* These are actually only local variables to x86cpu,
 * but then gcc complains hard and says they might
 * be clobbered by the longjmp. We get rid of this
 * warning by making them global...
 */
static int cpuSSE,cpu3DNow,doSSE,do3DNow,cpuflags;
static int cpuAltivec;
  
static jmp_buf mainloop;
static int success=1;

static void sigill_handler(int n)
{
  success=0;
  longjmp(mainloop,n);
}


/* The general detection routine at the end of this
 * file calls one of the cpu-specific routines below.
 */

#ifdef USE_PPC_ALTIVEC
/* separate small routine that performs an altivec instruction -
 * but make sure we use the result so it is not optimized away.
 * We use the oppurtunity to set the control register to non-java mode.
 */
static int check_altivec(void)
{
#ifdef __VEC__
  vector unsigned short vsr1,vsr2;
  vsr1=vec_mfvscr();
  vsr2=(vector unsigned short)vec_sl(vec_splat_u32(1),vec_splat_u32(16));
  vsr1=vec_or(vsr1,vsr2);
  vec_mtvscr(vsr1);

  return 1;
#else
  return 0;
#endif
}
  
static int detect_altivec(FILE *log)
{
  /* dont know how to detect ppc, but this
   * will only be called on ppc, so never mind for now...
   * Even linux 2.2 should support altivec, so we dont need
   * to perform the extensive tests like for x86 to recommend
   * users to upgrade. Just check if an instruction works!
   */
  cpuflags=cpuAltivec=0;
  
  /* if we got here, its ppc! */
  cpuflags |= PPC_CPU;

  if(log)
    fprintf(log,"\nTesting PPC Altivec support...\n");
  success=1;
  signal(SIGILL,sigill_handler);
  setjmp(mainloop); /* return to this point if we get SIGILL */
  if(success)
    success=check_altivec();
  
  if(success) {
    cpuflags |= PPC_ALTIVEC_SUPPORT;
    fprintf(log,"CPU supports Altivec.\n"
	    "Using Gromacs Altivec innerloops.\n\n");
  } else {
    fprintf(log,"No Altivec support found.\n");
  }
  return cpuflags;
}
#endif


#ifdef USE_X86_ASM
static int detect_sse3dnow(FILE *log)
{
  unsigned long eax,ebx,ecx,edx;
  
  cpuSSE=cpu3DNow=doSSE=do3DNow=cpuflags=0;

  if(log)    
    fprintf(log,"\nTesting x86 processor CPUID...\n");
  /* start by trying to issue the cpuid instruction */
  success=1;
  signal(SIGILL,sigill_handler);
  setjmp(mainloop); /* return to this point if we get SIGILL */

  if(success) 
    x86_cpuid(0,&eax,&ebx,&ecx,&edx);
  else if(log)
    fprintf(log,"This CPU doesn't support CPUID.\n");

  if(eax>0) {
    cpuflags |= X86_CPU;
    fprintf(log,"\nTesting x86 SSE capabilities...\n");
    if(ebx==VENDOR_INTEL) {
      /* intel - we need SSE support, bit 25 of edx should be set */
      x86_cpuid(1,&eax,&ebx,&ecx,&edx);
      if(edx & FLAGS_SUPPORT_SSE) {
	cpuSSE=1;
	/* try it */
	success=1;
	setjmp(mainloop); /* return to this point if we get SIGILL */
	if(success) {
	  checksse();
	  doSSE=1;
	  cpuflags |= X86_SSE_SUPPORT;
	}
      }
    }

    fprintf(log,"\nTesting x86 3DNow capabilities...\n");
    if(ebx==VENDOR_AMD) {
      /* Newer athlons support SSE which might be faster (or at least IEEE=better),
       * but I havent been able to test one yet, so it is not turned on yet.
       * In theory it should be easy - just insert a test for bit 25 of edx in
       * the same way we do for intel above. 
       */
      /* amd - start by checking for extended functions */
      x86_cpuid(0x80000000,&eax,&ebx,&ecx,&edx);
      if(eax>=0x80000001) {
	x86_cpuid(0x80000001,&eax,&ebx,&ecx,&edx);
	if(edx & FLAGS_SUPPORT_EXT_3DNOW) {
	  cpu3DNow=1;
	  /* try it */
	  success=1;
	  setjmp(mainloop); /* return to this point if we get SIGILL */
	  if(success) {
	    check3dnow();
	    do3DNow=1;
	    cpuflags |= X86_3DNOW_SUPPORT;
	  }
	}
      }
    }
  }
  if(doSSE) {
    if(log)
      fprintf(log,"CPU and OS support SSE.\n"
	      "Using Gromacs SSE assembly innerloops.\n\n");
  } else if(do3DNow) {
    if(log)
      fprintf(log,"CPU and OS support extended 3DNow.\n"
	      "Using Gromacs 3DNow assembly innerloops.\n\n");
  } else if(log) {
    if(cpuSSE) {
      fprintf(log,"CPU supports SSE, but your OS doesn't.\n");
      fprintf(stderr,"NOTE: This version of gromacs is compiled with assembly innerloops\n" 
	      "      using Intel SSE instructions. Your processor supports this,\n"
	      "      but not your OS. Fixing this (e.g., upgrade to linux kernel 2.4)\n"
	      "      will boost your gromacs performance SIGNIFICANTLY.\n");
    } else if(cpu3DNow) {      
      fprintf(log,"CPU supports extended 3DNow, but your OS doesn't.\n");
    } else if(!cpuSSE && !cpu3DNow) {
      fprintf(log,"No SSE/3DNow support found.\n");
    }
    fprintf(log,"Using normal Gromacs innerloops.\n\n");
  }
  return cpuflags; 
} 
#endif


int detect_cpu(FILE *log)
{
#ifdef USE_PPC_ALTIVEC
  return detect_altivec(log);
#elif defined USE_X86_ASM
  return detect_sse3dnow(log);
#endif
  fprintf(log,"Not checking cpu support for SSE/3DNow/Altivec\n");
  return UNKNOWN_CPU;
}



#endif /* detectcpu.h */

