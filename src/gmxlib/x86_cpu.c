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
static char *SRCID_x86_cpu_c = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <x86_cpu.h>

#ifdef USE_X86_ASM
/* there is a fallback routine at the end
 * of this file which returns UNKNOWN_CPU
 * if we do not use x86 support
 */
#include <signal.h>
#include <setjmp.h>

/* These are actually only local variables to x86cpu,
 * but then gcc complains hard and says they might
 * be clobbered by the longjmp. We get rid of this
 * warning by making them global...
 */
static int cpuSSE,cpu3DNow,doSSE,do3DNow,cpuflags;
  
static jmp_buf mainloop;
static int success=1;

static void sigill_handler(int n)
{
  success=0;
  longjmp(mainloop,n);
}

int check_x86cpu(FILE *log)
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
#ifdef USE_X86_ASM
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
#endif    
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
#ifdef USE_X86_ASM
      fprintf(log,"No SSE/3DNow support found.\n");
#else
      fprintf(log,"Not checking SEE/3DNow support.\n");
#endif      
    }
    fprintf(log,"Using normal Gromacs innerloops.\n\n");
  }
  return cpuflags; 
}

#else /* no x86 support */
int check_x86cpu(FILE *log)
{
  return UNKNOWN_CPU;
}
#endif

