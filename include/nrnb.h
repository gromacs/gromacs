/*
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _nrnb_h
#define _nrnb_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_nrnb(t_nrnb *nrnb);

void cp_nrnb(t_nrnb *dest, t_nrnb *src);

void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2);

void print_nrnb(FILE *out, t_nrnb *nrnb);

void _inc_nrnb(t_nrnb *nrnb,int enr,int inc,char *file,int line);

#if DEBUG_NRNB
#define inc_nrnb(nrnb,enr,inc) _inc_nrnb(nrnb,enr,inc,__FILE__,__LINE__)
#else
#define inc_nrnb(nrnb,enr,inc) (nrnb)->n[enr] += inc
#endif

 
void print_flop(FILE *out,t_nrnb *nrnb,double *nbfs,double *mflop);
/* Calculates the non-bonded forces and flop count.
 * When out!=NULL also prints the full count table.
 */

void print_perf(FILE *out,double nodetime,double realtime,int nprocs,
		       gmx_large_int_t nsteps,real delta_t,
		       double nbfs,double mflop);
/* Prints the performance, nbfs and mflop come from print_flop */

void pr_load(FILE *log,t_commrec *cr,t_nrnb nrnb[]);
/* Print detailed load balancing info */

int cost_nrnb(int enr);
/* Cost in i860 cycles of this component of MD */

const char *nrnb_str(int enr);
/* Name of this component */

#ifdef __cplusplus
}
#endif

#endif	/* _nrnb_h */
