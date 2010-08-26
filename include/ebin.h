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

#ifndef _ebin_h
#define _ebin_h

#include "sysstuff.h"
#include "typedefs.h"
#include "enxio.h"	

#ifdef __cplusplus
extern "C" {
#endif


/* This is a running averaging structure ('energy bin') for use during mdrun. */
typedef struct {
  int        nener;
  gmx_enxnm_t *enm;
  gmx_large_int_t nsteps;
  gmx_large_int_t nsum;
  t_energy   *e;
  gmx_large_int_t nsteps_sim;
  gmx_large_int_t nsum_sim;
  t_energy   *e_sim;
} t_ebin;

enum { eprNORMAL, eprAVER, eprRMS, eprNR };

t_ebin *mk_ebin(void);
/* Create an energy bin */

int get_ebin_space(t_ebin *eb,int nener,const char *enm[],const char *unit);

/* Create space in the energy bin and register names.
 * The enm array must be static, because the contents are not copied,
 * but only the pointers.
 * Function returns an index number that must be used in subsequent
 * calls to add_ebin.
 */

void add_ebin(t_ebin *eb,int index,int nener,real ener[],gmx_bool bSum);
/* Add nener reals (eg. energies, box-lengths, pressures) to the
 * energy bin at position index. 
 * If bSum is TRUE then the reals are also added to the sum
 * and sum of squares.
 */

void ebin_increase_count(t_ebin *eb,gmx_bool bSum);
/* Increase the counters for the sums.
 * This routine should be called AFTER all add_ebin calls for this step.
 */

void reset_ebin_sums(t_ebin *eb);
/* Reset the average and fluctuation sums */

void pr_ebin(FILE *fp,t_ebin *eb,int index,int nener,int nperline,
		    int prmode,gmx_bool bPrHead);
/* Print the contents of the energy bin. If nener = -1 ALL energies from
 * index to the end will be printed. We will print nperline entries on a text
 * line (advisory <= 5). prmode may be any of the above listed enum values.
 * tsteps is used only when eprAVER or eprRMS is set.
 * If bPrHead than the header is printed.
 */

#ifdef __cplusplus
}
#endif

#endif	/* _ebin_h */
