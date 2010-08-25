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

#ifndef _vcm_h
#define _vcm_h

#include "sysstuff.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  int    nr;                   /* Number of groups                    */
  int    mode;                 /* One of the enums above              */
  gmx_bool   ndim;                 /* The number of dimensions for corr.  */
  real   *group_ndf;           /* Number of degrees of freedom        */
  rvec   *group_p;             /* Linear momentum per group           */
  rvec   *group_v;             /* Linear velocity per group           */
  rvec   *group_x;             /* Center of mass per group            */
  rvec   *group_j;             /* Angular momentum per group          */
  rvec   *group_w;             /* Angular velocity (omega)            */
  tensor *group_i;             /* Moment of inertia per group         */
  real   *group_mass;          /* Mass per group                      */
  char   **group_name;         /* These two are copies to pointers in */
} t_vcm;

t_vcm *init_vcm(FILE *fp,gmx_groups_t *groups,t_inputrec *ir);

/* Do a per group center of mass things */
void calc_vcm_grp(FILE *fp,int start,int homenr,t_mdatoms *md,
			 rvec x[],rvec v[],t_vcm *vcm);

void do_stopcm_grp(FILE *fp,int start,int homenr,
			  unsigned short *group_id,
			  rvec x[],rvec v[],t_vcm *vcm);

void check_cm_grp(FILE *fp,t_vcm *vcm,t_inputrec *ir,real Temp_Max);


#ifdef __cplusplus
}
#endif

			 
#endif /* _vcm_h */
