/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/* This file is completely threadsafe - keep it that way! */

#ifdef USE_THREADS
#include <pthread.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <detectcpu.h>

/* This file uses some global variables to catch illegal instruction
 * signals from the operating system. To make it thread-safe, we use
 * a global mutex which is locked immediately when the detect_cpu()
 * routine is called (if we use threads). Do NOT call any of the 
 * lower-level CPU detection routines directly, since these do not
 * lock the mutex.
 */


#if (defined USE_X86_SSE_AND_3DNOW || defined USE_X86_SSE2 || defined USE_PPC_ALTIVEC)
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

#ifdef USE_THREADS
static pthread_mutex_t detectcpu_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif
  
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
void check_altivec(void);


static int detect_altivec(FILE *log)
{
  /* dont know how to detect ppc, but this
   * will only be called on ppc, so never mind for now...
   * Even linux 2.2 should support altivec, so we dont need
   * to perform the extensive tests like for x86 to recommend
   * users to upgrade. Just check if an instruction works!
   */
  int cpuflags;

  /* BUG warning: Due to problems in the Apple gcc compiler,
   * the cpuflags variable MUST be assigned inside the
   * if-statement below, and returned inside the same
   * statement.
   */
  if(log)
    fprintf(log,"\nTesting PPC Altivec support...\n");
  signal(SIGILL,sigill_handler);
  if(sigsetjmp(mainloop,0)) { /* return here with val=1 if we get SIGILL */
    if(log)
      fprintf(log,"No Altivec support found.\n");
    cpuflags = PPC_CPU;
    return cpuflags;
  } else {
    check_altivec();

    cpuflags = PPC_CPU | PPC_ALTIVEC_SUPPORT;
    if(log)
      fprintf(log,"CPU supports Altivec.\n"
              "Using Gromacs Altivec innerloops.\n\n");
    return cpuflags;
  }
  /* We used to return here, outside the if-statement, but
   * that exposed a bug in Apple gcc-3.x, so be VERY
   * careful and test it if you ever change this code...
   */
}
#endif /* PPC_ALTIVEC */



#ifdef USE_X86_SSE2
static int detect_sse2(FILE *log)
{
  unsigned long eax,ebx,ecx,edx;
  int cpuflags=0;
  
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
    if(log)
      fprintf(log,"\nTesting x86 SSE2 capabilities...\n");
    /* for SSE support, bit 25 of edx should be set */
    x86_cpuid(1,&eax,&ebx,&ecx,&edx);
    if(edx & FLAGS_SUPPORT_SSE2) {
      cpuflags |= CPU_SSE2_SUPPORT;/* CPU ok, but OS might not be */
      /* try it */
      success=1;
      setjmp(mainloop); /* return to this point if we get SIGILL */
      if(success) {
	checksse2();
	cpuflags |= X86_SSE2_SUPPORT;
      }
    }
  }
  return cpuflags;
}
#endif /* sse2 */

#ifdef USE_X86_SSE_AND_3DNOW
static int detect_sse(FILE *log)
{
  unsigned long eax,ebx,ecx,edx;
  
  int cpuflags=0;
  
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
    if(log)
      fprintf(log,"\nTesting x86 SSE capabilities...\n");
    /* for SSE support, bit 25 of edx should be set */
    x86_cpuid(1,&eax,&ebx,&ecx,&edx);
    if(edx & FLAGS_SUPPORT_SSE) {
      cpuflags |= CPU_SSE_SUPPORT;/* CPU ok, but OS might not be */
      /* try it */
      success=1;
      setjmp(mainloop); /* return to this point if we get SIGILL */
      if(success) {
	checksse();
	cpuflags |= X86_SSE_SUPPORT;
      }
    }
  }
  return cpuflags;
}

static int detect_3dnow(FILE *log)
{
  unsigned long eax,ebx,ecx,edx;
  
  int cpuflags=0;
  
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
    if(log)
      fprintf(log,"\nTesting x86 3DNow capabilities...\n");
    if(ebx==VENDOR_AMD) {
      /* amd - start by checking for extended functions */
      x86_cpuid(0x80000000,&eax,&ebx,&ecx,&edx);
      if(eax>=0x80000001) {
	x86_cpuid(0x80000001,&eax,&ebx,&ecx,&edx);
	if(edx & FLAGS_SUPPORT_EXT_3DNOW) {
	  /* try it */
	  success=1;
	  setjmp(mainloop); /* return to this point if we get SIGILL */
	  if(success) {
	    check3dnow();
	    cpuflags |= X86_3DNOW_SUPPORT;
	  }
	}
      }
    }
  }
  return cpuflags;
}
#endif /* sse and 3dnow */

#if (defined USE_X86_SSE_AND_3DNOW || defined USE_X86_SSE2)
static int detect_x86(FILE *log)
{
  unsigned long eax,ebx,ecx,edx;

  int cpuflags=0;

  /* First we check if this is an Intel or AMD cpu */
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
  } else {
    if(log)
      fprintf(log,"Didn't recognize this as an x86-type CPU.\n");
    return cpuflags;
  }
  
  if(ebx==VENDOR_INTEL) {
    if(log) 
      fprintf(log,"CPU manufactured by Intel.\n");
    cpuflags |= INTEL_CPU;
  }
  else if(ebx==VENDOR_AMD) {
     if(log) 
      fprintf(log,"CPU manufactured by AMD.\n");
    cpuflags |= AMD_CPU;
  }
    
#ifdef USE_X86_SSE2
  /* For double precision we need SSE2 */
  cpuflags=detect_sse2(log);
  if(cpuflags & CPU_SSE2_SUPPORT) {
    if (cpuflags & X86_SSE2_SUPPORT)
      fprintf(log,"CPU and OS support SSE2.\n"
	      "Using Gromacs SSE2 double precision assembly innerloops.\n\n");
    else {
      fprintf(log,"CPU supports SSE2, but your OS doesn't.\n");
      fprintf(stderr,"NOTE: This version of gromacs is compiled with assembly innerloops\n" 
	      "      using Intel SSE2 instructions. Your processor supports this,\n"
	      "      but not your OS. Fixing this (e.g., upgrade your linux kernel)\n"
	      "      will boost your gromacs performance SIGNIFICANTLY.\n");
    }
  } else 
    fprintf(log,"No SSE2 support found for this CPU.\n");
  
#endif   
#ifdef USE_X86_SSE_AND_3DNOW
  /* Either SSE or 3DNow will do for single precision.
   * We first check for SSE since it is IEEE-compliant
   */
  cpuflags=detect_sse(log);
  if(cpuflags & X86_SSE_SUPPORT) {
    fprintf(log,"CPU and OS support SSE.\n"
	    "Using Gromacs SSE single precision assembly innerloops.\n\n");
    return cpuflags;
  } else {
    if(cpuflags & CPU_SSE_SUPPORT)
      fprintf(log,"Your processor supports SSE, but not your OS.\n"
	      "Checking for Extended 3DNow support instead...\n");
    else /* No sse support at all */
      fprintf(log,"No SSE support in CPU. Trying Extended 3DNow...\n");
    cpuflags=detect_3dnow(log);
    if(cpuflags & X86_3DNOW_SUPPORT) {
      fprintf(log,"CPU and OS support extended 3DNow.\n"
	      "Using Gromacs 3DNow single precision assembly innerloops.\n\n");
      return cpuflags;
    }
  }
  fprintf(log,"SSE or Extended 3DNow not supported.\n");
  fprintf(log,"Using normal Gromacs innerloops.\n\n");
  /* Tell the user if he has SSE support but cannot use it ... */
  if(cpuflags & CPU_SSE_SUPPORT)
    fprintf(stderr,"NOTE: This version of gromacs is compiled with assembly innerloops\n" 
	    "      using Intel SSE instructions. Your processor supports this,\n"
	    "      but not your OS. Fixing this (e.g., upgrade your linux kernel)\n"
	    "      will boost your gromacs performance SIGNIFICANTLY.\n");
#endif
  return cpuflags;
}
#endif


int 
detect_cpu(FILE *log)
{
  int cpuflags=0;

  /* Use a global lock for all CPU detection */
#ifdef USE_THREADS
  pthread_mutex_lock(&detectcpu_mutex);
#endif

#ifdef USE_PPC_ALTIVEC
  cpuflags=detect_altivec(log);
#elif (defined USE_X86_SSE_AND_3DNOW || defined USE_X86_SSE2)
  cpuflags=detect_x86(log);
#else  
  if(log)
    fprintf(log,"Not checking cpu support for SSE/SSE2/3DNow/Altivec\n");
  cpuflags=UNKNOWN_CPU;
#endif

  /* release lock */
#ifdef USE_THREADS
  pthread_mutex_unlock(&detectcpu_mutex);
#endif

  if (getenv("NOASSEMBLYLOOPS") != NULL) {
    cpuflags &= (~X86_SSE_SUPPORT) & (~X86_SSE2_SUPPORT) & (~X86_3DNOW_SUPPORT) & (~PPC_ALTIVEC_SUPPORT);
    
    if(log)
      fprintf(log,
	      "Found environment variable NOASSEMBLYLOOPS.\n"
	      "Disabling all SSE/SSE2/3DNow/Altivec support.\n");
  }

  return cpuflags;
}
#endif /* detectcpu */




