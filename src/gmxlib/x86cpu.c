#include <stdio.h>
#include <signal.h>
#include <setjmp.h>
#include "x86cpu.h"



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
  int cpuSSE;
  int cpu3DNow;
  int doSSE;
  int do3DNow;

  cpuSSE=cpu3DNow=doSSE=do3DNow=0;

  if(log)
    fprintf(log,"\nTesting x86 processor and OS capabilities:\n");
  /* start by trying to issue the cpuid instruction */
  success=1;
  signal(SIGILL,sigill_handler);
  setjmp(mainloop); /* return to this point if we get SIGILL */

  if(success) 
    gmxcpuid(0,&eax,&ebx,&ecx,&edx);
  else if(log)
    fprintf(log,"This CPU doesn't support CPUID.\n");

  if(eax>0) {
    if(ebx==VENDOR_INTEL) {
      /* intel - we need SSE support, bit 25 of edx should be set */
      gmxcpuid(1,&eax,&ebx,&ecx,&edx);
      if(edx & FLAGS_SUPPORT_SSE) {
	cpuSSE=1;
	/* try it */
	success=1;
	setjmp(mainloop); /* return to this point if we get SIGILL */
	if(success) {
	  checksse();
	  doSSE=1;
	}
      }
    } else if(ebx==VENDOR_AMD) {
      /* amd - start by checking for extended functions */
      gmxcpuid(0x80000000,&eax,&ebx,&ecx,&edx);
      if(eax>=0x80000001) {
	gmxcpuid(0x80000001,&eax,&ebx,&ecx,&edx);
	if(edx & FLAGS_SUPPORT_EXT_3DNOW) {
	  cpu3DNow=1;
	  /* try it */
	  success=1;
	  setjmp(mainloop); /* return to this point if we get SIGILL */
	  if(success) {
	    check3dnow();
	    do3DNow=1;
	  }
	}
      }
    }
  }
  if(doSSE) {
    if(log)
      fprintf(log,"CPU and OS support SSE.\n"
	      "Using Gromacs SSE assembly innerloops.\n\n");
    return X86_SSE;
  } else if(do3DNow) {
    if(log)
      fprintf(log,"CPU and OS support extended 3DNow.\n"
	      "Using Gromacs 3DNow assembly innerloops.\n\n");
    return X86_3DNOW;
  } else if(log) {
    if(!cpuSSE && !cpu3DNow)
      fprintf(log,"No SSE or 3DNow support found.\n");
    else if(cpuSSE)
      fprintf(log,"CPU supports SSE, but your OS doesn't.\n");
    else if(cpu3DNow)
      fprintf(log,"CPU supports extended 3DNow, but your OS doesn't.\n");
    fprintf(log,"Using normal Gromacs innerloops.\n\n");
  }
  return X86_NOOPT; 
}

