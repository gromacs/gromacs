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
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _x86_cpu_h
#define _x86_cpu_h

static char *SRCID_x86_cpu_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

extern int cpu_capabilities;
#define UNKNOWN_CPU           0

#ifdef USE_X86_ASM

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
