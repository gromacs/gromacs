/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int  n;		/* Number of terms				*/
  real *a;		/* Coeffients (V / nm )                  	*/
  real *phi;		/* Phase angles					*/
} t_cosines;

typedef struct {
  int     ngtc;         /* # T-Coupl groups                             */
  int     ngacc;        /* # Accelerate groups                          */
  int     ngfrz;        /* # Freeze groups                              */
  int     ngener;	/* # Ener groups				*/
  real    *nrdf;	/* Nr of degrees of freedom in a group		*/
  real    *ref_t;	/* Coupling temperature	per group		*/
  real    *tau_t;	/* Tau coupling time 				*/
  rvec    *acc;		/* Acceleration per group			*/
  ivec    *nFreeze;	/* Freeze the group in each direction ?	        */
  bool    *eg_excl;     /* Exclusions of energy group pairs             */
} t_grpopts;

typedef struct {
  int  eI;              /* Integration method 				*/
  int  nsteps;		/* number of steps to be taken			*/
  int  ns_type;		/* which ns method should we use?               */
  int  nstlist;		/* number of steps before pairlist is generated	*/
  int  ndelta;		/* number of cells per rlong			*/
  bool bDomDecomp;      /* Should we do domain decomposition?           */
  int  decomp_dir;      /* Direction of decomposition (may not be opt.) */
  int  nstcomm;		/* number of steps after which center of mass	*/
                        /* motion is removed				*/
  int nstlog;		/* number of steps after which print to logfile	*/
  int nstxout;		/* number of steps after which X is output	*/
  int nstvout;		/* id. for V					*/
  int nstfout;		/* id. for F					*/
  int nstenergy;	/* number of steps after which energies printed */
  int nstxtcout;	/* id. for compressed trj (.xtc)		*/
  real init_t;		/* initial time (ps) 				*/
  real delta_t;		/* time step (ps)				*/
  real xtcprec;         /* precision of xtc file                        */
  int  nkx,nky,nkz;     /* number of k vectors in each spatial dimension*/
                        /* for fourier methods for long range electrost.*/
  int  pme_order;       /* interpolation order for PME                  */
  real ewald_rtol;      /* Real space tolerance for Ewald, determines   */
                        /* the real/reciprocal space relative weight    */
  int  ewald_geometry;  /* normal/3d ewald, or pseudo-2d LR corrections */
  bool epsilon_surface; /* Epsilon for PME dipole correction            */
  bool bOptFFT;         /* optimize the fft plan at start               */
  int  ePBC;		/* Type of periodic boundary conditions		*/
  bool bUncStart;       /* Do not constrain the start configuration	*/
  int  etc;		/* temperature coupling         		*/
  int  epc;		/* pressure coupling                            */
  int  epct;		/* pressure coupling type			*/
  real tau_p;		/* pressure coupling time (ps)			*/
  tensor ref_p;		/* reference pressure (kJ/(mol nm^3))		*/
  tensor compress;	/* compressability ((mol nm^3)/kJ) 		*/
  bool bSimAnn;         /* simulated annealing (SA)                     */
  real zero_temp_time;  /* time at which temp becomes zero in sim. ann. */
  real rlist;		/* short range pairlist cut-off (nm)		*/
  int  coulombtype;	/* Type of electrostatics treatment             */
  real rcoulomb_switch; /* Coulomb switch range start (nm)		*/
  real rcoulomb;        /* Coulomb cutoff (nm)		                */
  int  vdwtype;         /* Type of Van der Waals treatment              */
  real rvdw_switch;     /* Van der Waals switch range start (nm)        */
  real rvdw;		/* Van der Waals cutoff (nm)		        */
  real epsilon_r;       /* relative dielectric constant                 */
  int  eDispCorr;       /* Perform Long range dispersion corrections    */
  real shake_tol;	/* tolerance for shake				*/
  real fudgeQQ;		/* Id. for 1-4 coulomb interactions		*/
  int  efep;   		/* free energy interpolation no/yes		*/
  real init_lambda;	/* initial value for perturbation variable	*/
  real delta_lambda;	/* change of lambda per time step (1/dt)	*/
  real sc_alpha;        /* free energy soft-core parameter              */
  real sc_sigma;        /* free energy soft-core sigma when c6 or c12=0 */
  real dr_fc;		/* force constant for ta_disre			*/
  int  eDisreWeighting; /* type of weighting of pairs in one restraints	*/
  bool bDisreMixed;     /* Use comb of time averaged and instan. viol's	*/
  int  nstdisreout;     /* frequency of writing pair distances to enx   */ 
  real dr_tau;		/* time constant for memory function in disres 	*/
  real orires_fc;	/* force constant for orientational restraints  */
  real orires_tau;	/* time constant for memory function in orires 	*/
  int  nstorireout;     /* frequency of writing tr(SD) to enx           */ 
  real em_stepsize;	/* The stepsize for updating			*/
  real em_tol;		/* The tolerance				*/
  int  niter;           /* Number of iterations for convergence of      */
                        /* steepest descent in relax_shells             */
  int  fc_stepsize;     /* Stepsize for directional minimization        */
                        /* in relax_shells                              */
  int  nstcgsteep;      /* number of steps after which a steepest       */
                        /* descents step is done while doing cg         */
  int  eConstrAlg;      /* Type of constraint algorithm                 */
  int  nProjOrder;      /* Order of the LINCS Projection Algorithm      */
  real LincsWarnAngle;  /* If bond rotates more than %g degrees, warn   */
  bool bShakeSOR;       /* Use successive overrelaxation for shake      */
  real bd_temp;         /* Temperature for Brownian Dynamics (BD)       */
  real bd_fric;         /* Friction coefficient for BD (amu / ps)       */
  int  ld_seed;         /* Random seed for SD and BD                    */
  real cos_accel;       /* Acceleration for viscosity calculation       */
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
} t_inputrec;

