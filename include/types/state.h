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
#ifndef _state_h_
#define _state_h_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 * The t_state struct should contain all the (possibly) non-static
 * information required to define the state of the system.
 * Currently the random seeds for SD and BD are missing.
 */

/* These enums are used in flags as (1<<est...).
 * The order of these enums should not be changed,
 * since that affects the checkpoint (.cpt) file format.
 */
enum { estLAMBDA,
       estBOX, estBOX_REL, estBOXV, estPRES_PREV, estNH_XI,  estTC_INT,
       estX,   estV,       estSDX,  estCGP,       estLD_RNG, estLD_RNGI,
       estDISRE_INITF, estDISRE_RM3TAV,
       estORIRE_INITF, estORIRE_DTAV,
       estNR };

/* The names of the state entries, defined in src/gmxib/checkpoint.c */
extern const char *est_names[estNR];

typedef struct
{
  real disre_initf;    /* The scaling factor for initializing the time av. */
  int  ndisrepairs;    /* The number of distance restraints                */
  real *disre_rm3tav;  /* The r^-3 time averaged pair distances            */
  real orire_initf;    /* The scaling factor for initializing the time av. */
  int  norire_Dtav;    /* The number of matrix element in dtav (npair*5)   */
  real *orire_Dtav;    /* The time averaged orientation tensors            */
} history_t;

/* Struct used for checkpointing only.
 * This struct would not be required with unlimited precision.
 * But because of limited precision, the COM motion removal implementation
 * can cause the kinetic energy in the MD loop to differ by a few bits from
 * the kinetic energy one would determine from state.v.
 */
typedef struct
{
  int      ekinh_n;
  tensor * ekinh;
  real     dekindl;
  real     mvcos;
} ekinstate_t;

typedef struct
{
  double *   ener_ave;  /* Energy term history sum to get fluctuations      */
  double *   ener_sum;  /* Energy term history sum to get fluctuations      */
  int        nener;     /* Number of energy terms in two previous arrays    */
}
energyhistory_t;

typedef struct
{
  int           natoms;
  int           ngtc;
  int           nrng;
  int           nrngi;
  int           flags;  /* Flags telling which entries are present      */
  real          lambda; /* the free energy switching parameter          */
  matrix 	    box;    /* box vector coordinates                      	*/
  matrix     	box_rel; /* Relitaive box vectors to preserve shape    	*/
  matrix 	    boxv;   /* box velocitites for Parrinello-Rahman pcoupl */
  matrix        pres_prev; /* Pressure of the previous step for pcoupl  */
  real          *nosehoover_xi;  /* for Nose-Hoover tcoupl (ngtc)       */
  double        *therm_integral; /* for N-H/V-rescale tcoupl (ngtc)     */
  int           nalloc; /* Allocation size for x, v and sd_x when !=NULL*/
  rvec          *x;     /* the coordinates (natoms)                     */
  rvec          *v;     /* the velocities (natoms)                      */
  rvec          *sd_X;  /* random part of the x update for stoch. dyn.  */
  rvec          *cg_p;  /* p vector for conjugate gradient minimization */

  unsigned int  *ld_rng;  /* RNG random state                           */
  int           *ld_rngi; /* RNG index                                  */

  history_t     hist;   /* Time history for restraints                  */

  ekinstate_t   ekinstate; /* The state of the kinetic energy data      */

  energyhistory_t  enerhist; /* Energy history for statistics           */
	
  int           ddp_count; /* The DD partitioning count for this state  */
  int           ddp_count_cg_gl; /* The DD part. count for index_gl     */
  int           ncg_gl; /* The number of local charge groups            */
  int           *cg_gl; /* The global cg number of the local cgs        */
  int           cg_gl_nalloc; /* Allocation size of cg_gl;              */
} t_state;

#endif /* _state_h_ */
