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

enum { eNL_VDWQQ, eNL_VDW, eNL_QQ, 
       eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE, 
       eNL_VDWQQ_SOLMNO, eNL_VDW_SOLMNO, eNL_QQ_SOLMNO, 
       eNL_VDWQQ_WATER, eNL_QQ_WATER, 
       eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER, 
       eNL_NR };
 
typedef struct {
  /* Cut-Off stuff */
  int  ePBC;
  real rlist,rlistlong;
  
  /* Dielectric constant resp. multiplication factor for charges */
  real zsquare,temp;
  real epsilon_r,epsfac;  
  
  /* Constants for reaction fields */
  bool bRF;
  real kappa,k_rf,c_rf;
  
  /* Integrated differces for energy and virial with cut-off functions */
  real enerdiffsix;
  real enerdifftwelve;
  real virdiffsix;
  real virdifftwelve;
  /* Constant for long range dispersion correction (average dispersion) */
  real avcsix;
  /* Constant for long range repulsion term. Relative difference of about 
   * 0.1 percent with 0.8 nm cutoffs. But hey, it's cheap anyway...
   */
  real avctwelve;
  
  /* Fudge factors */
  real fudgeQQ;

  /* Table stuff */
  bool bcoultab;
  bool bvdwtab;
  real rtab;
  int  ntab;
  real tabscale;
  /* We duplicate tables for cache optimization purposes */
  real *coultab;      /* Coul only */
  real *vdwtab;       /* Vdw only   */
  real *coulvdwtab;   /* Both      */
  real *coulvdw14tab; /* 1,4 table with both */

  /* PPPM & Shifting stuff */
  real rcoulomb_switch,rcoulomb;
  real *phi;

  /* VdW stuff */
  real rvdw_switch,rvdw;
  real bham_b_max;
  real tabscale_exp;

  /* Free energy ? */
  int  efep;
  real sc_alpha;
  real sc_sigma6;
  bool bSepDVDL;

  /* NS Stuff */
  int  eeltype;
  int  vdwtype;
  int  cg0,hcg;
  int  ndelta;
  bool bSolvOpt;
  int  nMNOMol;
  real nMNOav[3];
  int  nWatMol;
  int  Dimension;
  bool bGrid,bDomDecomp;
  int  *solvent_type;
  int  *mno_index;
  rvec *cg_cm;
  rvec *shift_vec;
  
  /* The actual neighbor lists, short and long range, see enum above
   * for definition of neighborlist indices.
   */
  t_nblist nlist_sr[eNL_NR];
  t_nblist nlist_lr[eNL_NR];
  
  /* This mask array of length nn determines whether or not this bit of the
   * neighbourlists should be computed. Usually all these are true of course,
   * but not when shells are used. During minimisation all the forces that 
   * include shells are done, then after minimsation is converged the remaining
   * forces are computed.
   */
  /* bool *bMask; */
    
  /* Twin Range stuff. */
  bool bTwinRange;
  int  nlr;
  rvec *f_twin;
  rvec *fshift_twin;
  
  /* PME/Ewald stuff */
  rvec *f_pme;
  bool bEwald;
  real ewaldcoeff;

  /* Virial Stuff */
  rvec *fshift;
  
  /* Free energy stuff */
  int     nmol;
  atom_id *mol_nr;
  real    *mol_epot;
  int     nstcalc;
  
  /* Non bonded Parameter lists */
  int  ntype; /* Number of atom types */
  bool bBHAM;
  real *nbfp;

  /* Energy group exclusions */
  bool *eg_excl;

  /* xmdrun flexible constraints */
  real fc_stepsize;

  /* Generalized born stuff */
  /* VdW radius for each atomtype (dim is thus ntype) */
  real *atype_radius;
  /* Effective radius (derived from effective volume) for each type */
  real *atype_vol;
  /* Implicit solvent - surface tension for each atomtype */
  real *atype_surftens;

} t_forcerec;

#define C6(nbfp,ntp,ai,aj)     (nbfp)[2*((ntp)*(ai)+(aj))]
#define C12(nbfp,ntp,ai,aj)    (nbfp)[2*((ntp)*(ai)+(aj))+1]
#define BHAMC(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))]
#define BHAMA(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))+1]
#define BHAMB(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))+2]
