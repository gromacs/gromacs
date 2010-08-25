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
#ifndef _inputrec_h_
#define _inputrec_h_


#include "simple.h"
#include "../sysstuff.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  int  n;		/* Number of terms				*/
  real *a;		/* Coeffients (V / nm )                  	*/
  real *phi;		/* Phase angles					*/
} t_cosines;

typedef struct {
  real E0;              /* Field strength (V/nm)                        */
  real omega;           /* Frequency (1/ps)                             */
  real t0;              /* Centre of the Gaussian pulse (ps)            */
  real sigma;           /* Width of the Gaussian pulse (FWHM) (ps)      */
} t_efield;

#define EGP_EXCL  (1<<0)
#define EGP_TABLE (1<<1)

typedef struct {
  int     ngtc;                  /* # T-Coupl groups                        */
  int     nhchainlength;         /* # of nose-hoover chains per group       */
  int     ngacc;                 /* # Accelerate groups                     */
  int     ngfrz;                 /* # Freeze groups                         */
  int     ngener;	         /* # Ener groups			    */
  real    *nrdf;	         /* Nr of degrees of freedom in a group	    */
  real    *ref_t;	         /* Coupling temperature	per group   */
  int     *annealing;            /* No/simple/periodic SA for each group    */
  int     *anneal_npoints;       /* Number of annealing time points per grp */    
  real    **anneal_time;         /* For ea. group: Time points              */
  real    **anneal_temp;         /* For ea. grp: Temperature at these times */
                                 /* Final temp after all intervals is ref_t */ 
  real    *tau_t;	         /* Tau coupling time 			    */
  rvec    *acc;		         /* Acceleration per group		    */
  ivec    *nFreeze;	         /* Freeze the group in each direction ?    */
  int     *egp_flags;            /* Exclusions/tables of energy group pairs */

  /* QMMM stuff */
  int     ngQM;         /* nr of QM groups                              */
  int     *QMmethod;    /* Level of theory in the QM calculation        */
  int     *QMbasis;     /* Basisset in the QM calculation               */
  int     *QMcharge;    /* Total charge in the QM region                */
  int     *QMmult;      /* Spin multiplicicty in the QM region          */
  gmx_bool    *bSH;         /* surface hopping (diabatic hop only)          */
  int     *CASorbitals; /* number of orbiatls in the active space       */
  int     *CASelectrons;/* number of electrons in the active space      */
  real    *SAon;        /* at which gap (A.U.) the SA is switched on    */
  real    *SAoff;
  int     *SAsteps;     /* in how many steps SA goes from 1-1 to 0.5-0.5*/
  gmx_bool    *bOPT;
  gmx_bool    *bTS;
} t_grpopts;

enum { epgrppbcNONE, epgrppbcREFAT, epgrppbcCOS };

typedef struct {
  int        nat;      /* Number of atoms in the pull group */
  atom_id    *ind;     /* The global atoms numbers */
  int        nat_loc;  /* Number of local pull atoms */
  int        nalloc_loc; /* Allocation size for ind_loc and weight_loc */ 
  atom_id    *ind_loc; /* Local pull indices */
  int        nweight;  /* The number of weights (0 or nat) */
  real       *weight;  /* Weights (use all 1 when weight==NULL) */
  real       *weight_loc; /* Weights for the local indices */
  int        epgrppbc; /* The type of pbc for this pull group, see enum above */
  atom_id    pbcatom;  /* The reference atom for pbc (global number) */
  rvec       vec;      /* The pull vector, direction or position */
  rvec       init;     /* Initial reference displacement */
  real       rate;     /* Rate of motion (nm/ps) */
  real       k;        /* force constant */
  real       kB;       /* force constant for state B */
  real       wscale;   /* scaling factor for the weights: sum w m/sum w w m */
  real       invtm;    /* inverse total mass of the group: 1/wscale sum w m */
  dvec       x;        /* center of mass before update */
  dvec       xp;       /* center of mass after update before constraining */
  dvec       dr;       /* The distance from the reference group */
  double     f_scal;   /* Scalar force for directional pulling */
  dvec       f;        /* force due to the pulling/constraining */
} t_pullgrp; 

typedef struct {
  int        ngrp;        /* number of groups */
  int        eGeom;       /* pull geometry */
  ivec       dim;         /* used to select components for constraint */
  real       cyl_r1;      /* radius of cylinder for dynamic COM */
  real       cyl_r0;      /* radius of cylinder including switch length */
  real       constr_tol;  /* absolute tolerance for constraints in (nm) */
  int        nstxout;     /* Output frequency for pull x */
  int        nstfout;     /* Output frequency for pull f */
  int        ePBC;        /* the boundary conditions */
  int        npbcdim;     /* do pbc in dims 0 <= dim < npbcdim */
  gmx_bool       bRefAt;      /* do we need reference atoms for a group COM ? */
  int        cosdim;      /* dimension for cosine weighting, -1 if none */
  gmx_bool       bVirial;     /* do we need to add the pull virial? */
  t_pullgrp  *grp;        /* groups to pull/restrain/etc/ */
  t_pullgrp  *dyna;       /* dynamic groups for use with local constraints */
  rvec       *rbuf;       /* COM calculation buffer */
  dvec       *dbuf;       /* COM calculation buffer */
  double     *dbuf_cyl;   /* cylinder ref. groups COM calculation buffer */

  FILE       *out_x;      /* output file for pull data */
  FILE       *out_f;      /* output file for pull data */
} t_pull;

typedef struct {
  int  eI;              /* Integration method 				*/
  gmx_large_int_t nsteps;	/* number of steps to be taken			*/
  int  simulation_part; /* Used in checkpointing to separate chunks */
  gmx_large_int_t init_step;	/* start at a stepcount >0 (used w. tpbconv)    */
  int  nstcalcenergy;	/* fequency of energy calc. and T/P coupl. upd.	*/
  int  ns_type;		/* which ns method should we use?               */
  int  nstlist;		/* number of steps before pairlist is generated	*/
  int  ndelta;		/* number of cells per rlong			*/
  int  nstcomm;		/* number of steps after which center of mass	*/
                        /* motion is removed				*/
  int  comm_mode;       /* Center of mass motion removal algorithm      */
  int nstcheckpoint;    /* checkpointing frequency                      */
  int nstlog;		/* number of steps after which print to logfile	*/
  int nstxout;		/* number of steps after which X is output	*/
  int nstvout;		/* id. for V					*/
  int nstfout;		/* id. for F					*/
  int nstenergy;	/* number of steps after which energies printed */
  int nstxtcout;	/* id. for compressed trj (.xtc)		*/
  double init_t;	/* initial time (ps) 				*/
  double delta_t;	/* time step (ps)				*/
  real xtcprec;         /* precision of xtc file                        */
  int  nkx,nky,nkz;     /* number of k vectors in each spatial dimension*/
                        /* for fourier methods for long range electrost.*/
  int  pme_order;       /* interpolation order for PME                  */
  real ewald_rtol;      /* Real space tolerance for Ewald, determines   */
                        /* the real/reciprocal space relative weight    */
  int  ewald_geometry;  /* normal/3d ewald, or pseudo-2d LR corrections */
  real epsilon_surface; /* Epsilon for PME dipole correction            */
  gmx_bool bOptFFT;         /* optimize the fft plan at start               */
  int  ePBC;		/* Type of periodic boundary conditions		*/
  int  bPeriodicMols;   /* Periodic molecules                           */
  gmx_bool bContinuation;   /* Continuation run: starting state is correct	*/
  int  etc;		/* temperature coupling         		*/
  int  nsttcouple;      /* interval in steps for temperature coupling   */
  int  epc;		/* pressure coupling                            */
  int  epct;		/* pressure coupling type			*/
  int  nstpcouple;      /* interval in steps for pressure coupling      */
  real tau_p;		/* pressure coupling time (ps)			*/
  tensor ref_p;		/* reference pressure (kJ/(mol nm^3))		*/
  tensor compress;	/* compressability ((mol nm^3)/kJ) 		*/
  int  refcoord_scaling;/* How to scale absolute reference coordinates  */
  rvec posres_com;      /* The COM of the posres atoms                  */
  rvec posres_comB;     /* The B-state COM of the posres atoms          */
  int  andersen_seed;   /* Random seed for Andersen thermostat.         */
  real rlist;		/* short range pairlist cut-off (nm)		*/
  real rlistlong;	/* long range pairlist cut-off (nm)		*/
  real rtpi;            /* Radius for test particle insertion           */
  int  coulombtype;	/* Type of electrostatics treatment             */
  real rcoulomb_switch; /* Coulomb switch range start (nm)		*/
  real rcoulomb;        /* Coulomb cutoff (nm)		                */
  real epsilon_r;       /* relative dielectric constant                 */ 
  real epsilon_rf;      /* relative dielectric constant of the RF       */ 
  int  implicit_solvent;/* No (=explicit water), or GBSA solvent models */
  int  gb_algorithm;    /* Algorithm to use for calculation Born radii  */
  int  nstgbradii;      /* Frequency of updating Generalized Born radii */
  real rgbradii;        /* Cutoff for GB radii calculation              */
  real gb_saltconc;     /* Salt concentration (M) for GBSA models       */
  real gb_epsilon_solvent; /* dielectric coeff. of implicit solvent     */
  real gb_obc_alpha;    /* 1st scaling factor for Bashford-Case GB      */
  real gb_obc_beta;     /* 2nd scaling factor for Bashford-Case GB      */
  real gb_obc_gamma;    /* 3rd scaling factor for Bashford-Case GB      */
  real gb_dielectric_offset; /* Dielectric offset for Still/HCT/OBC     */
  int  sa_algorithm;    /* Algorithm for SA part of GBSA                */
  real sa_surface_tension; /* Energy factor for SA part of GBSA */
  int  vdwtype;         /* Type of Van der Waals treatment              */
  real rvdw_switch;     /* Van der Waals switch range start (nm)        */
  real rvdw;		    /* Van der Waals cutoff (nm)	        */
  int  eDispCorr;       /* Perform Long range dispersion corrections    */
  real tabext;          /* Extension of the table beyond the cut-off,   *
		 	 * as well as the table length for 1-4 interac. */
  real shake_tol;	/* tolerance for shake				*/
  int  efep;   		/* free energy interpolation no/yes		*/
  double init_lambda;	/* initial value for perturbation variable	*/
  double delta_lambda;	/* change of lambda per time step (1/dt)	*/
  int  n_flambda;       /* The number of foreign lambda points          */
  double *flambda;      /* The foreign lambda values                    */
  real sc_alpha;        /* free energy soft-core parameter              */
  int  sc_power;        /* lambda power for soft-core interactions      */
  real sc_sigma;        /* free energy soft-core sigma when c6 or c12=0 */
  real sc_sigma_min;    /* minimum FE sc sigma (default: =sg_sigma)     */
  int  nstdhdl;         /* The frequency for writing to dhdl.xvg        */
  int  separate_dhdl_file; /* whether to write a separate dhdl.xvg file 
                              note: NOT a gmx_bool, but an enum */
  int  dhdl_derivatives;/* whether to calculate+write dhdl derivatives 
                              note: NOT a gmx_bool, but an enum */
  int  dh_hist_size;    /* The maximum size for the dH histogram        */
  double dh_hist_spacing; /* The spacing for the dH histogram           */
  int  eDisre;          /* Type of distance restraining                 */
  real dr_fc;		    /* force constant for ta_disre			*/
  int  eDisreWeighting; /* type of weighting of pairs in one restraints	*/
  gmx_bool bDisreMixed;     /* Use comb of time averaged and instan. viol's	*/
  int  nstdisreout;     /* frequency of writing pair distances to enx   */ 
  real dr_tau;		    /* time constant for memory function in disres 	*/
  real orires_fc;	    /* force constant for orientational restraints  */
  real orires_tau;	    /* time constant for memory function in orires 	*/
  int  nstorireout;     /* frequency of writing tr(SD) to enx           */ 
  real dihre_fc;        /* force constant for dihedral restraints	*/
  real em_stepsize;	    /* The stepsize for updating			*/
  real em_tol;		    /* The tolerance				*/
  int  niter;           /* Number of iterations for convergence of      */
                        /* steepest descent in relax_shells             */
  real fc_stepsize;     /* Stepsize for directional minimization        */
                        /* in relax_shells                              */
  int  nstcgsteep;      /* number of steps after which a steepest       */
                        /* descents step is done while doing cg         */
  int  nbfgscorr;       /* Number of corrections to the hessian to keep */
  int  eConstrAlg;      /* Type of constraint algorithm                 */
  int  nProjOrder;      /* Order of the LINCS Projection Algorithm      */
  real LincsWarnAngle;  /* If bond rotates more than %g degrees, warn   */
  int  nLincsIter;      /* Number of iterations in the final Lincs step */
  gmx_bool bShakeSOR;       /* Use successive overrelaxation for shake      */
  real bd_fric;         /* Friction coefficient for BD (amu/ps)         */
  int  ld_seed;         /* Random seed for SD and BD                    */
  int  nwall;           /* The number of walls                          */
  int  wall_type;       /* The type of walls                            */
  real wall_r_linpot;   /* The potentail is linear for r<=wall_r_linpot */
  int  wall_atomtype[2];/* The atom type for walls                      */
  real wall_density[2]; /* Number density for walls                     */
  real wall_ewald_zfac; /* Scaling factor for the box for Ewald         */
  int  ePull;           /* Type of pulling: no, umbrella or constraint  */
  t_pull *pull;         /* The data for center of mass pulling          */
  real cos_accel;       /* Acceleration for viscosity calculation       */
  tensor deform;        /* Triclinic deformation velocities (nm/ps)     */
  int  userint1;        /* User determined parameters                   */
  int  userint2;
  int  userint3;
  int  userint4;
  real userreal1;
  real userreal2;
  real userreal3;
  real userreal4;
  t_grpopts opts;	/* Group options				*/
  t_cosines ex[DIM];	/* Electric field stuff	(spatial part)		*/
  t_cosines et[DIM];	/* Electric field stuff	(time part)		*/
  gmx_bool bQMMM;           /* QM/MM calculation                            */ 
  int  QMconstraints;   /* constraints on QM bonds                      */
  int  QMMMscheme;      /* Scheme: ONIOM or normal                      */
  real scalefactor;     /* factor for scaling the MM charges in QM calc.*/
} t_inputrec;

#define DEFORM(ir) ((ir).deform[XX][XX]!=0 || (ir).deform[YY][YY]!=0 || (ir).deform[ZZ][ZZ]!=0 || (ir).deform[YY][XX]!=0 || (ir).deform[ZZ][XX]!=0 || (ir).deform[ZZ][YY]!=0)

#define DYNAMIC_BOX(ir) ((ir).epc!=epcNO || (ir).eI==eiTPI || DEFORM(ir))

#define PRESERVE_SHAPE(ir) ((ir).epc != epcNO && (ir).deform[XX][XX] == 0 && ((ir).epct == epctISOTROPIC || (ir).epct == epctSEMIISOTROPIC))

#define NEED_MUTOT(ir) (((ir).coulombtype==eelEWALD || EEL_PME((ir).coulombtype)) && ((ir).ewald_geometry==eewg3DC || (ir).epsilon_surface!=0))

#define IR_TWINRANGE(ir) ((ir).rlist > 0 && ((ir).rlistlong == 0 || (ir).rlistlong > (ir).rlist))

#define IR_ELEC_FIELD(ir) ((ir).ex[XX].n > 0 || (ir).ex[YY].n > 0 || (ir).ex[ZZ].n > 0)

#define IR_EXCL_FORCES(ir) (EEL_FULL((ir).coulombtype) || (EEL_RF((ir).coulombtype) && (ir).coulombtype != eelRF_NEC) || (ir).implicit_solvent != eisNO)
/* use pointer definitions of ir here, since that's what's usually used in the code */
#define IR_NVT_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && ((ir)->etc == etcNOSEHOOVER))

#define IR_NPT_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && ((ir)->epc == epcMTTK))

#ifdef __cplusplus
}
#endif


#endif
