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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _nblist_h
#define _nblist_h

#ifdef __cplusplus
extern "C" {
#endif

/* Neighborlist type */
enum {
  enlistATOM_ATOM,
  enlistSPC_ATOM,   enlistSPC_SPC,
  enlistTIP4P_ATOM, enlistTIP4P_TIP4P,
  enlistCG_CG,
  enlistNR
};

typedef unsigned long t_excl;

/* The maximum charge group size because of minimum size of t_excl
 * could be 32 bits.
 */
#define MAX_CHARGEGROUP_SIZE 32

/* The maximum charge group size for CG-CG nblists.
 * The excl entry in t_nblist uses blocks of this size.
 */
#define MAX_CGCGSIZE 32

typedef struct 
{
  int             enlist;      /* The type of nblist, enum, see above    */
  int             il_code;      /* Innerloop index from nrnb.h, used     */
                                /* for flop accounting.                  */
  int             icoul;        /* Coulomb loop type index for kernels   */
  int             ivdw;         /* VdW loop type index for kernels       */
  int             free_energy;  /* Free energy setting for this list     */

  int             nri,maxnri;   /* Current/max number of i particles	 */
  int             nrj,maxnrj;   /* Current/max number of j particles	 */
  int             maxlen;       /* maxnr of j atoms for a single i atom  */
  int *           iinr;	        /* The i-elements	   	         */
  int *           iinr_end;     /* The end atom, only with enlistCG      */
  int *           gid;          /* Index in energy arrays                */
  int *           shift;        /* Shift vector index                    */
  int *           jindex;       /* Index in jjnr                         */
  int *           jjnr;	        /* The j-atom list                       */
  int *           jjnr_end;     /* The end atom, only with enltypeCG     */
  t_excl *        excl;         /* Exclusions, only with enltypeCG       */
  int             count;        /* counter to multithread the innerloops */
  void *          mtx;          /* mutex to lock the counter             */
} t_nblist;


/* For atom I =  nblist->iinr[N] (0 <= N < nblist->nri) there can be
 * several neighborlists (N's), for different energy groups (gid) and
 * different shifts (shift).
 * For corresponding J atoms for each list start at:
 * nblist->jjnr[JI]
 * with nblist->jindex[N] <= JI < nblist->jindex[N+1]
 *
 * enlist is of the form enlistUNIT1_UNIT2:
 * UNIT ATOM:  there is one atom: iinr[N] or jjnr[JI]
 * UNIT SPC:   there are 3 atoms: iinr[N],iinr[N]+1,iinr[N]+2, jjnr analog.
 * UNIT TIP4P: there are 4 atoms: iinr[N],...,iinr[N]+3, jjnr analog.
 * UNIT CG:    there are N atoms: iinr[N],...,iinr_end[N]-1, jjnr analog.
 *
 * Clear?
 */

#ifdef __cplusplus
}
#endif

#endif
