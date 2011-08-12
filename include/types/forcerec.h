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

#include "ns.h"
#include "genborn.h"
#include "qmmmrec.h"
#include "idef.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Abstract type for PME that is defined only in the routine that use them. */
typedef struct gmx_pme *gmx_pme_t;

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
  t_forcetable tab;
  /* We duplicate tables for cache optimization purposes */
  real *coultab;      /* Coul only */
  real *vdwtab;       /* Vdw only   */
  /* The actual neighbor lists, short and long range, see enum above
   * for definition of neighborlist indices.
   */
  t_nblist nlist_sr[eNL_NR];
  t_nblist nlist_lr[eNL_NR];
} t_nblists;

/* macros for the cginfo data in forcerec */
/* The maximum cg size in cginfo is 255,
 * because we only have space for 8 bits in cginfo,
 * this cg size entry is actually only read with domain decomposition.
 * But there is a smaller limit due to the t_excl data structure
 * which is defined in nblist.h.
 */
#define SET_CGINFO_GID(cgi,gid)      (cgi) = (((cgi)  &  ~65535)  |  (gid)   )
#define GET_CGINFO_GID(cgi)        ( (cgi)            &   65535)
#define SET_CGINFO_EXCL_INTRA(cgi)   (cgi) =  ((cgi)  |  (1<<16))
#define GET_CGINFO_EXCL_INTRA(cgi) ( (cgi)            &  (1<<16))
#define SET_CGINFO_EXCL_INTER(cgi)   (cgi) =  ((cgi)  |  (1<<17))
#define GET_CGINFO_EXCL_INTER(cgi) ( (cgi)            &  (1<<17))
#define SET_CGINFO_SOLOPT(cgi,opt)   (cgi) = (((cgi)  & ~(15<<18)) | ((opt)<<18))
#define GET_CGINFO_SOLOPT(cgi)     (((cgi)>>18)       &   15)
/* This bit is only used with bBondComm in the domain decomposition */
#define SET_CGINFO_BOND_INTER(cgi)   (cgi) =  ((cgi)  |  (1<<22))
#define GET_CGINFO_BOND_INTER(cgi) ( (cgi)            &  (1<<22))
#define SET_CGINFO_NATOMS(cgi,opt)   (cgi) = (((cgi)  & ~(255<<23)) | ((opt)<<23))
#define GET_CGINFO_NATOMS(cgi)     (((cgi)>>23)       &   255)


/* Value to be used in mdrun for an infinite cut-off.
 * Since we need to compare with the cut-off squared,
 * this value should be slighlty smaller than sqrt(GMX_FLOAT_MAX).
 */
#define GMX_CUTOFF_INF 1E+18


enum { egCOULSR, egLJSR, egBHAMSR, egCOULLR, egLJLR, egBHAMLR,
       egCOUL14, egLJ14, egGB, egNR };

typedef struct {
  int  nener;        /* The number of energy group pairs     */
  real *ener[egNR];  /* Energy terms for each pair of groups */
} gmx_grppairener_t;

typedef struct {
  real term[F_NRE];    /* The energies for all different interaction types */
  gmx_grppairener_t grpp;
  double dvdl_lin;     /* Contributions to dvdl with linear lam-dependence */
  double dvdl_nonlin;  /* Idem, but non-linear dependence                  */
  int    n_lambda;
  double *enerpart_lambda; /* Partial energy for lambda and flambda[] */
} gmx_enerdata_t;
/* The idea is that dvdl terms with linear lambda dependence will be added
 * automatically to enerpart_lambda. Terms with non-linear lambda dependence
 * should explicitly determine the energies at foreign lambda points
 * when n_lambda > 0.
 */

typedef struct {
  int cg_start;
  int cg_end;
  int cg_mod;
  int *cginfo;
} cginfo_mb_t;


/* ewald table type */
typedef struct ewald_tab *ewald_tab_t; 

typedef struct {
  /* Domain Decomposition */
  gmx_bool bDomDec;

  /* PBC stuff */
  int  ePBC;
  gmx_bool bMolPBC;
  int  rc_scaling;
  rvec posres_com;
  rvec posres_comB;

  gmx_bool UseOptimizedKernels;

  /* Use special N*N kernels? */
  gmx_bool bAllvsAll;
  /* Private work data */
  void *AllvsAll_work;
  void *AllvsAll_workgb;

  /* Cut-Off stuff.
   * Infinite cut-off's will be GMX_CUTOFF_INF (unlike in t_inputrec: 0).
   */
  real rlist,rlistlong;
  
  /* Dielectric constant resp. multiplication factor for charges */
  real zsquare,temp;
  real epsilon_r,epsilon_rf,epsfac;  
  
  /* Constants for reaction fields */
  real kappa,k_rf,c_rf;

  /* Charge sum and dipole for topology A/B ([0]/[1]) for Ewald corrections */
  double qsum[2];
  rvec   mu_tot[2];

  /* Dispersion correction stuff */
  int  eDispCorr;
  /* The shift of the shift or user potentials */
  real enershiftsix;
  real enershifttwelve;
  /* Integrated differces for energy and virial with cut-off functions */
  real enerdiffsix;
  real enerdifftwelve;
  real virdiffsix;
  real virdifftwelve;
  /* Constant for long range dispersion correction (average dispersion)
   * for topology A/B ([0]/[1]) */
  real avcsix[2];
  /* Constant for long range repulsion term. Relative difference of about 
   * 0.1 percent with 0.8 nm cutoffs. But hey, it's cheap anyway...
   */
  real avctwelve[2];
  
  /* Fudge factors */
  real fudgeQQ;

  /* Table stuff */
  gmx_bool bcoultab;
  gmx_bool bvdwtab;
  /* The normal tables are in the nblists struct(s) below */
  t_forcetable tab14; /* for 1-4 interactions only */

  /* PPPM & Shifting stuff */
  real rcoulomb_switch,rcoulomb;
  real *phi;

  /* VdW stuff */
  double reppow;
  real rvdw_switch,rvdw;
  real bham_b_max;

  /* Free energy ? */
  int  efep;
  real sc_alpha;
  int  sc_power;
  real sc_sigma6_def;
  real sc_sigma6_min;
  gmx_bool bSepDVDL;

  /* NS Stuff */
  int  eeltype;
  int  vdwtype;
  int  cg0,hcg;
  /* solvent_opt contains the enum for the most common solvent
   * in the system, which will be optimized.
   * It can be set to esolNO to disable all water optimization */
  int  solvent_opt;
  int  nWatMol;
  gmx_bool bGrid;
  cginfo_mb_t *cginfo_mb;
  int  *cginfo;
  rvec *cg_cm;
  int  cg_nalloc;
  rvec *shift_vec;

  /* The neighborlists including tables */
  int  nnblists;
  int  *gid2nblists;
  t_nblists *nblists;

  /* The wall tables (if used) */
  int  nwall;
  t_forcetable **wall_tab;

  /* This mask array of length nn determines whether or not this bit of the
   * neighbourlists should be computed. Usually all these are true of course,
   * but not when shells are used. During minimisation all the forces that 
   * include shells are done, then after minimsation is converged the remaining
   * forces are computed.
   */
  /* gmx_bool *bMask; */

  /* The number of charge groups participating in do_force_lowlevel */
  int ncg_force;
  /* The number of atoms participating in do_force_lowlevel */
  int natoms_force;
  /* The number of atoms participating in force and constraints */
  int natoms_force_constr;
  /* The allocation size of vectors of size natoms_force */
  int nalloc_force;

  /* Twin Range stuff, f_twin has size natoms_force */
  gmx_bool bTwinRange;
  int  nlr;
  rvec *f_twin;

  /* Forces that should not enter into the virial summation:
   * PPPM/PME/Ewald/posres
   */
  gmx_bool bF_NoVirSum;
  int  f_novirsum_n;
  int  f_novirsum_nalloc;
  rvec *f_novirsum_alloc;
  /* Pointer that points to f_novirsum_alloc when pressure is calcaluted,
   * points to the normal force vectors wen pressure is not requested.
   */
  rvec *f_novirsum;

  /* Long-range forces and virial for PPPM/PME/Ewald */
  gmx_pme_t pmedata;
  tensor    vir_el_recip;

  /* PME/Ewald stuff */
  gmx_bool bEwald;
  real ewaldcoeff;
  ewald_tab_t ewald_table;

  /* Virial Stuff */
  rvec *fshift;
  rvec vir_diag_posres;
  dvec vir_wall_z;

  /* Non bonded Parameter lists */
  int  ntype; /* Number of atom types */
  gmx_bool bBHAM;
  real *nbfp;

  /* Energy group pair flags */
  int *egp_flags;

  /* xmdrun flexible constraints */
  real fc_stepsize;

  /* Generalized born implicit solvent */
  gmx_bool bGB;
  /* Generalized born stuff */
  real gb_epsilon_solvent;
  /* Table data for GB */
  t_forcetable gbtab;
  /* VdW radius for each atomtype (dim is thus ntype) */
  real *atype_radius;
  /* Effective radius (derived from effective volume) for each type */
  real *atype_vol;
  /* Implicit solvent - surface tension for each atomtype */
  real *atype_surftens;
  /* Implicit solvent - radius for GB calculation */
  real *atype_gb_radius;
  /* Implicit solvent - overlap for HCT model */
  real *atype_S_hct;
  /* Generalized born interaction data */
  gmx_genborn_t *born;

  /* Table scale for GB */
  real gbtabscale;
  /* Table range for GB */
  real gbtabr;
  /* GB neighborlists (the sr list will contain for each atom all other atoms                                            
   * (for use in the SA calculation) and the lr list will contain                                                  
   * for each atom all atoms 1-4 or greater (for use in the GB calculation)                                        
   */
  t_nblist gblist_sr;
  t_nblist gblist_lr;
  t_nblist gblist;

  /* Inverse square root of the Born radii for implicit solvent */
  real *invsqrta;
  /* Derivatives of the potential with respect to the Born radii */
  real *dvda;
  /* Derivatives of the Born radii with respect to coordinates */
  real *dadx;
  real *dadx_rawptr;
  int   nalloc_dadx; /* Allocated size of dadx */
	
  /* If > 0 signals Test Particle Insertion,
   * the value is the number of atoms of the molecule to insert
   * Only the energy difference due to the addition of the last molecule
   * should be calculated.
   */
  gmx_bool n_tpi;

  /* Neighbor searching stuff */
  gmx_ns_t ns;

  /* QMMM stuff */
  gmx_bool         bQMMM;
  t_QMMMrec    *qr;

  /* QM-MM neighborlists */
  t_nblist QMMMlist;

  /* Limit for printing large forces, negative is don't print */
  real print_force;

  /* coarse load balancing time measurement */
  double t_fnbf;
  double t_wait;
  int timesteps;

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

#ifdef __cplusplus
}
#endif

