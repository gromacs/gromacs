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
  real r;         /* range of the table */
  int  n;         /* n+1 is the number of points */
  real scale;     /* distance between two points */
  real scale_exp; /* distance for exponential Buckingham table */
  real *tab;      /* the actual tables, per point there are  4 numbers for
		   * Coulomb, dispersion and repulsion (in total 12 numbers)
		   */
} t_forcetable;
 
typedef struct {
  /* Cut-Off stuff */
  int  ePBC;
  real rlist,rlistlong;
  
  /* Dielectric constant resp. multiplication factor for charges */
  real zsquare,temp;
  real epsilon_r,epsfac;  
  
  /* Constants for reaction fields */
  real kappa,k_rf,c_rf;

  /* Charge sum for topology A/B ([0]/[1]) for Ewald corrections */
  double qsum[2];

  /* The shift of the shift or user potentials */
  real enershiftsix;
  real enershifttwelve;
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
  t_forcetable tab;
  /* We duplicate tables for cache optimization purposes */
  real *coultab;      /* Coul only */
  real *vdwtab;       /* Vdw only   */
  t_forcetable tab14; /* for 1-4 interactions only */

  /* PPPM & Shifting stuff */
  real rcoulomb_switch,rcoulomb;
  real *phi;

  /* VdW stuff */
  real rvdw_switch,rvdw;
  real bham_b_max;

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

  /* Long-range forces and virial for PPPM/PME/Ewald */
  rvec *f_el_recip;
  tensor vir_el_recip;

  /* PME/Ewald stuff */
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

  /* User determined parameters, copied from the inputrec */
  int  userint1;
  int  userint2;
  int  userint3;
  int  userint4;
  real userreal1;
  real userreal2;
  real userreal3;
  real userreal4;
} t_forcerec;

#define C6(nbfp,ntp,ai,aj)     (nbfp)[2*((ntp)*(ai)+(aj))]
#define C12(nbfp,ntp,ai,aj)    (nbfp)[2*((ntp)*(ai)+(aj))+1]
#define BHAMC(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))]
#define BHAMA(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))+1]
#define BHAMB(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))+2]
