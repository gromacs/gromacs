/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Grunge ROck MAChoS
 */

#ifndef _nrnb_h
#define _nrnb_h

static char *SRCID_nrnb_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) nrnb.h 1.9 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"

extern void init_nrnb(t_nrnb *nrnb);

extern void cp_nrnb(t_nrnb *dest, t_nrnb *src);

extern void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2);

extern void print_nrnb(FILE *out, t_nrnb *nrnb);

extern void _inc_nrnb(t_nrnb *nrnb,int enr,int inc,char *file,int line);

#ifdef DEBUG_NRNB
#define inc_nrnb(nrnb,enr,inc) _inc_nrnb(nrnb,enr,inc,__FILE__,__LINE__)
#else
#define inc_nrnb(nrnb,enr,inc) (nrnb)->n[enr] += inc
#endif

extern void print_perf(FILE *out,double cputime,double realtime,real runtime,
		       t_nrnb *nrnb,int nprocs);

extern void pr_load(FILE *log,int nprocs,t_nrnb nrnb[]);
/* Print detailed load balancing info */

extern int cost_nrnb(int enr);
/* Cost in i860 cycles of this component of MD */

extern char *nrnb_str(int enr);
/* Name of this component */

#endif	/* _nrnb_h */
