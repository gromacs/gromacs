/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _nrnb_h
#define _nrnb_h

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_nrnb(t_nrnb *nrnb);

void cp_nrnb(t_nrnb *dest, t_nrnb *src);

void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2);

void print_nrnb(FILE *out, t_nrnb *nrnb);

void _inc_nrnb(t_nrnb *nrnb, int enr, int inc, char *file, int line);

#ifdef DEBUG_NRNB
#define inc_nrnb(nrnb, enr, inc) _inc_nrnb(nrnb, enr, inc, __FILE__, __LINE__)
#else
#define inc_nrnb(nrnb, enr, inc) (nrnb)->n[enr] += inc
#endif


void print_flop(FILE *out, t_nrnb *nrnb, double *nbfs, double *mflop);
/* Calculates the non-bonded forces and flop count.
 * When out!=NULL also prints the full count table.
 */

void print_perf(FILE *out, double nodetime, double realtime,
                gmx_int64_t nsteps, double delta_t,
                double nbfs, double mflop);
/* Prints the performance, nbfs and mflop come from print_flop */

void pr_load(FILE *log, t_commrec *cr, t_nrnb nrnb[]);
/* Print detailed load balancing info */

int cost_nrnb(int enr);
/* Cost in i860 cycles of this component of MD */

const char *nrnb_str(int enr);
/* Name of this component */

#ifdef __cplusplus
}
#endif

#endif  /* _nrnb_h */
