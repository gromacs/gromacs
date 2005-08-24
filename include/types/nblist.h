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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct 
{
  int             il_code;      /* Innerloop index from nrnb.h, used     */
                                /* for flop accounting.                  */
  int             icoul;        /* Coulomb loop type index for kernels   */
  int             ivdw;         /* VdW loop type index for kernels       */
  int             free_energy;  /* Free energy setting for this list     */
  int             solvent_opt;  /* Atom, water, or water-water list      */

  int             nri,maxnri;   /* Current/max number of i particles	 */
  int             nrj,maxnrj;   /* Current/max number of j particles	 */
  int             maxlen;       /* maxnr of j atoms for a single i atom  */
  int *           iinr;	        /* The i-elements	   	         */
  int *           gid;          /* Index in energy arrays                */
  int *           shift;        /* Shift vector index                    */
  int *           jindex;       /* Index in jjnr                         */
  int *           jjnr;	        /* The j-atom list                       */
  int             count;        /* counter to multithread the innerloops */
  void *          mtx;          /* mutex to lock the counter             */
} t_nblist;


/* For atom I =  nblist->iinr[N] (0 <= N < nblist->nri) there can be
 * several neighborlists (N's), for different energy groups (gid) and
 * different shifts (shift).
 * For corresponding J atoms for each list are are:
 * nblist->jjnr[JI]
 * with nblist->jindex[N] <= JI < nblist->jindex[N+1]
 *
 * Clear?
 */










