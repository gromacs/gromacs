#ifndef _x86_cpu_h
#define _x86_cpu_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern int cpu_capabilities;
#define UNKNOWN_CPU           0

#if (defined USE_SSE || defined USE_3DNOW)

#define VENDOR_AMD   0x68747541
#define VENDOR_INTEL 0x756e6547
#define FLAGS_SUPPORT_SSE 0x02000000
#define FLAGS_SUPPORT_EXT_3DNOW 0xc0000000
/* Flags for x86 and future processor capabilities */
#define X86_CPU               1
#define X86_SSE_SUPPORT       (1 << 1)
#define X86_3DNOW_SUPPORT     (1 << 2)

#include <x86_sse.h>
#include <x86_3dnow.h>

int check_x86cpu(FILE *log);

/* Assembly routines in gmxcpuid.s */
void x86_cpuid(int,unsigned long *,unsigned long *,unsigned long *,unsigned long *);

#endif

#endif
