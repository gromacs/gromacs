/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

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
  int     *nrdf;	/* Nr of degrees of freedom in a group		*/
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
  int  solvent_opt;     /* atomtype of water oxygen                     */
  int  nsatoms;         /* Number of atoms in solvent molecule          */
  int  nkx,nky,nkz;     /* number of k vectors in each spatial dimension*/
                        /* for fourier methods for long range electrost.*/
  int  pme_order;       /* interpolation order for PME                  */
  real ewald_rtol;      /* Real space tolerance for Ewald, determines   */
                        /* the real/reciprocal space relative weight    */
  bool bOptFFT;         /* optimize the fft plan at start               */
  int  eBox;		/* The box type					*/
  bool bUncStart;       /* Do not constrain the start configuration	*/
  bool btc;		/* temperature coupling         		*/
  int  ntcmemory;       /* Memory (steps) for T coupling                */
  int  epc;		/* pressure coupling type			*/
  int  npcmemory;       /* Memory (steps) for P coupling                */
  real tau_p;		/* pressure coupling time (ps)			*/
  rvec ref_p;		/* reference pressure (kJ/(mol nm^3))		*/
  rvec compress;	/* compressability ((mol nm^3)/kJ) 		*/
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
  bool bDispCorr;       /* Perform Long range dispersion corrections    */
  real shake_tol;	/* tolerance for shake				*/
  real fudgeQQ;		/* Id. for 1-4 coulomb interactions		*/
  bool bPert;		/* is perturbation turned on			*/
  real init_lambda;	/* initial value for perturbation variable	*/
  real delta_lambda;	/* change of lambda per time step (1/dt)	*/
  real dr_fc;		/* force constant for ta_disre			*/
  int  eDisreWeighting; /* type of weighting of pairs in one restraints	*/
  bool bDisreMixed;     /* Use comb of time averaged and instan. viol's	*/
  int  nstdisreout;     /* frequency of writing pair distances to enx   */ 
  real dr_tau;		/* time constant for memory function in disres 	*/
  real em_stepsize;	/* The stepsize for updating			*/
  real em_tol;		/* The tolerance				*/
  int  niter;           /* Number of iterations for convergence of      */
                        /* steepest descent in relax_shells             */
  int  nstcgsteep;      /* number of steps after which a steepest       */
                        /* descents step is done while doing cg         */
  int  eConstrAlg;      /* Type of constraint algorithm                 */
  int  nProjOrder;      /* Order of the LINCS Projection Algorithm      */
  real LincsWarnAngle;  /* If bond rotates more than %g degrees, warn   */
  int  nstLincsout;     /* Frequency for output from LINCS              */
  real ld_temp;         /* Temperature for Langevin Dynamics (LD)       */
  real ld_fric;         /* Friction coefficient for LD (amu / ps)       */
  int  ld_seed;         /* Random seed for LD                           */
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

