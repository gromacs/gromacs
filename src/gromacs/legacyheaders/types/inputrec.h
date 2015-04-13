/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef _inputrec_h_
#define _inputrec_h_

#include <stdio.h>

#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/swap/enums.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
    int   n;    /* Number of terms				*/
    real *a;    /* Coeffients (V / nm )                     */
    real *phi;  /* Phase angles					*/
} t_cosines;

typedef struct {
    real E0;            /* Field strength (V/nm)                        */
    real omega;         /* Frequency (1/ps)                             */
    real t0;            /* Centre of the Gaussian pulse (ps)            */
    real sigma;         /* Width of the Gaussian pulse (FWHM) (ps)      */
} t_efield;

#define EGP_EXCL  (1<<0)
#define EGP_TABLE (1<<1)

typedef struct {
    int       ngtc;           /* # T-Coupl groups                        */
    int       nhchainlength;  /* # of nose-hoover chains per group       */
    int       ngacc;          /* # Accelerate groups                     */
    int       ngfrz;          /* # Freeze groups                         */
    int       ngener;         /* # Ener groups			    */
    real     *nrdf;           /* Nr of degrees of freedom in a group	    */
    real     *ref_t;          /* Coupling temperature	per group   */
    int      *annealing;      /* No/simple/periodic SA for each group    */
    int      *anneal_npoints; /* Number of annealing time points per grp */
    real    **anneal_time;    /* For ea. group: Time points              */
    real    **anneal_temp;    /* For ea. grp: Temperature at these times */
                              /* Final temp after all intervals is ref_t */
    real     *tau_t;          /* Tau coupling time              */
    rvec     *acc;            /* Acceleration per group		    */
    ivec     *nFreeze;        /* Freeze the group in each direction ?    */
    int      *egp_flags;      /* Exclusions/tables of energy group pairs */

    /* QMMM stuff */
    int          ngQM;         /* nr of QM groups                              */
    int         *QMmethod;     /* Level of theory in the QM calculation        */
    int         *QMbasis;      /* Basisset in the QM calculation               */
    int         *QMcharge;     /* Total charge in the QM region                */
    int         *QMmult;       /* Spin multiplicicty in the QM region          */
    gmx_bool    *bSH;          /* surface hopping (diabatic hop only)          */
    int         *CASorbitals;  /* number of orbiatls in the active space       */
    int         *CASelectrons; /* number of electrons in the active space      */
    real        *SAon;         /* at which gap (A.U.) the SA is switched on    */
    real        *SAoff;
    int         *SAsteps;      /* in how many steps SA goes from 1-1 to 0.5-0.5*/
    gmx_bool    *bOPT;
    gmx_bool    *bTS;
} t_grpopts;

typedef struct {
    int         nat;        /* Number of atoms in the pull group */
    atom_id    *ind;        /* The global atoms numbers */
    int         nweight;    /* The number of weights (0 or nat) */
    real       *weight;     /* Weights (use all 1 when weight==NULL) */
    atom_id     pbcatom;    /* The reference atom for pbc (global number) */
} t_pull_group;

typedef struct {
    int         group[4];   /* The pull groups, index in group in t_pull */
    int         eType;      /* The pull type: umbrella, constraint, ... */
    int         eGeom;      /* The pull geometry */
    ivec        dim;        /* Used to select components for constraint */
    rvec        origin;     /* The origin for the absolute reference */
    rvec        vec;        /* The pull vector, direction or position */
    gmx_bool    bStart;     /* Set init based on the initial structure */
    real        init;       /* Initial reference displacement */
    real        rate;       /* Rate of motion (nm/ps) */
    real        k;          /* force constant */
    real        kB;         /* force constant for state B */
} t_pull_coord;

typedef struct {
    int   eSimTempScale; /* simulated temperature scaling; linear or exponential */
    real  simtemp_low;   /* the low temperature for simulated tempering  */
    real  simtemp_high;  /* the high temperature for simulated tempering */
    real *temperatures;  /* the range of temperatures used for simulated tempering */
} t_simtemp;

typedef struct {
    int    nstdhdl;                 /* The frequency for calculating dhdl           */
    double init_lambda;             /* fractional value of lambda (usually will use
                                       init_fep_state, this will only be for slow growth,
                                       and for legacy free energy code. Only has a
                                       valid value if positive)   */
    int      init_fep_state;        /* the initial number of the state                   */
    double   delta_lambda;          /* change of lambda per time step (fraction of (0.1) */
    int      edHdLPrintEnergy;      /* print no, total or potential energies in dhdl    */
    int      n_lambda;              /* The number of foreign lambda points               */
    double **all_lambda;            /* The array of all lambda values                    */
    int      lambda_neighbors;      /* The number of neighboring lambda states to
                                       calculate the energy for in up and down directions
                                       (-1 for all) */
    int      lambda_start_n;        /* The first lambda to calculate energies for */
    int      lambda_stop_n;         /* The last lambda +1 to calculate energies for */
    real     sc_alpha;              /* free energy soft-core parameter                   */
    int      sc_power;              /* lambda power for soft-core interactions           */
    real     sc_r_power;            /* r power for soft-core interactions                */
    real     sc_sigma;              /* free energy soft-core sigma when c6 or c12=0      */
    real     sc_sigma_min;          /* free energy soft-core sigma for ?????             */
    gmx_bool bScCoul;               /* use softcore for the coulomb portion as well (default FALSE) */
    gmx_bool separate_dvdl[efptNR]; /* whether to print the dvdl term associated with
                                       this term; if it is not specified as separate,
                                       it is lumped with the FEP term */
    int    separate_dhdl_file;      /* whether to write a separate dhdl.xvg file
                                       note: NOT a gmx_bool, but an enum */
    int    dhdl_derivatives;        /* whether to calculate+write dhdl derivatives
                                       note: NOT a gmx_bool, but an enum */
    int    dh_hist_size;            /* The maximum table size for the dH histogram */
    double dh_hist_spacing;         /* The spacing for the dH histogram */
} t_lambda;

typedef struct {
    int      nstexpanded;         /* The frequency of expanded ensemble state changes */
    int      elamstats;           /* which type of move updating do we use for lambda monte carlo (or no for none) */
    int      elmcmove;            /* what move set will be we using for state space moves */
    int      elmceq;              /* the method we use to decide of we have equilibrated the weights */
    int      equil_n_at_lam;      /* the minumum number of samples at each lambda for deciding whether we have reached a minimum */
    real     equil_wl_delta;      /* WL delta at which we stop equilibrating weights */
    real     equil_ratio;         /* use the ratio of weights (ratio of minimum to maximum) to decide when to stop equilibrating */
    int      equil_steps;         /* after equil_steps steps we stop equilibrating the weights */
    int      equil_samples;       /* after equil_samples total samples (steps/nstfep), we stop equilibrating the weights */
    int      lmc_seed;            /* random number seed for lambda mc switches */
    gmx_bool minvar;              /* whether to use minumum variance weighting */
    int      minvarmin;           /* the number of samples needed before kicking into minvar routine */
    real     minvar_const;        /* the offset for the variance in MinVar */
    int      c_range;             /* range of cvalues used for BAR */
    gmx_bool bSymmetrizedTMatrix; /* whether to print symmetrized matrices */
    int      nstTij;              /* How frequently to print the transition matrices */
    int      lmc_repeats;         /* number of repetitions in the MC lambda jumps */  /*MRS -- VERIFY THIS */
    int      lmc_forced_nstart;   /* minimum number of samples for each state before free sampling */ /* MRS -- VERIFY THIS! */
    int      gibbsdeltalam;       /* distance in lambda space for the gibbs interval */
    real     wl_scale;            /* scaling factor for wang-landau */
    real     wl_ratio;            /* ratio between largest and smallest number for freezing the weights */
    real     init_wl_delta;       /* starting delta for wang-landau */
    gmx_bool bWLoneovert;         /* use one over t convergece for wang-landau when the delta get sufficiently small */
    gmx_bool bInit_weights;       /* did we initialize the weights? TODO: REMOVE FOR 5.0, no longer needed with new logic */
    real     mc_temp;             /* To override the main temperature, or define it if it's not defined */
    real    *init_lambda_weights; /* user-specified initial weights to start with  */
} t_expanded;

typedef struct {
    int            ngroup;         /* number of pull groups */
    int            ncoord;         /* number of pull coordinates */
    real           cylinder_r;     /* radius of cylinder for dynamic COM */
    real           constr_tol;     /* absolute tolerance for constraints in (nm) */
    gmx_bool       bPrintCOM1;     /* Print coordinates of COM 1 for each coord */
    gmx_bool       bPrintCOM2;     /* Print coordinates of COM 2 for each coord */
    gmx_bool       bPrintRefValue; /* Print the reference value for each coord */
    gmx_bool       bPrintComp;     /* Print cartesian components for each coord with geometry=distance */
    int            nstxout;        /* Output frequency for pull x */
    int            nstfout;        /* Output frequency for pull f */

    t_pull_group  *group;          /* groups to pull/restrain/etc/ */
    t_pull_coord  *coord;          /* the pull coordinates */
} pull_params_t;

/* Abstract type for COM pull caclulations only defined in the pull module */
struct pull_t;


/* Abstract types for enforced rotation only defined in pull_rotation.c       */
typedef struct gmx_enfrot *gmx_enfrot_t;
typedef struct gmx_enfrotgrp *gmx_enfrotgrp_t;

typedef struct {
    int         eType;             /* Rotation type for this group                  */
    int         bMassW;            /* Use mass-weighed positions?                   */
    int         nat;               /* Number of atoms in the group                  */
    atom_id    *ind;               /* The global atoms numbers                      */
    rvec       *x_ref;             /* The reference positions                       */
    rvec        vec;               /* The normalized rotation vector                */
    real        rate;              /* Rate of rotation (degree/ps)                  */
    real        k;                 /* Force constant (kJ/(mol nm^2)                 */
    rvec        pivot;             /* Pivot point of rotation axis (nm)             */
    int         eFittype;          /* Type of fit to determine actual group angle   */
    int         PotAngle_nstep;    /* Number of angles around the reference angle
                                      for which the rotation potential is also
                                      evaluated (for fit type 'potential' only)     */
    real            PotAngle_step; /* Distance between two angles in degrees (for
                                      fit type 'potential' only)                    */
    real            slab_dist;     /* Slab distance (nm)                            */
    real            min_gaussian;  /* Minimum value the gaussian must have so that
                                      the force is actually evaluated               */
    real            eps;           /* Additive constant for radial motion2 and
                                      flexible2 potentials (nm^2)                   */
    gmx_enfrotgrp_t enfrotgrp;     /* Stores non-inputrec rotation data per group   */
} t_rotgrp;

typedef struct {
    int          ngrp;       /* Number of rotation groups                     */
    int          nstrout;    /* Output frequency for main rotation outfile    */
    int          nstsout;    /* Output frequency for per-slab data            */
    t_rotgrp    *grp;        /* Groups to rotate                              */
    gmx_enfrot_t enfrot;     /* Stores non-inputrec enforced rotation data    */
} t_rot;

/* Abstract type for IMD only defined in IMD.c */
typedef struct gmx_IMD *t_gmx_IMD;

typedef struct {
    int         nat;         /* Number of interactive atoms                   */
    atom_id    *ind;         /* The global indices of the interactive atoms   */
    t_gmx_IMD   setup;       /* Stores non-inputrec IMD data                  */
} t_IMD;

/* Abstract types for position swapping only defined in swapcoords.c */
typedef struct t_swap *gmx_swapcoords_t;

typedef struct {
    int              nstswap;           /* Every how many steps a swap is attempted?    */
    int              nat;               /* Number of atoms in the ion group             */
    int              nat_split[2];      /* Number of atoms in the split group           */
    int              nat_sol;           /* Number of atoms in the solvent group         */
    atom_id         *ind;               /* The global ion group atoms numbers           */
    atom_id         *ind_split[2];      /* Split groups for compartment partitioning    */
    atom_id         *ind_sol;           /* The global solvent group atom numbers        */
    gmx_bool         massw_split[2];    /* Use mass-weighted positions in split group?  */
    real             cyl0r, cyl1r;      /* Split cylinders defined by radius, upper and */
    real             cyl0u, cyl1u;      /* ... lower extension. The split cylinders de- */
    real             cyl0l, cyl1l;      /* ... fine the channels and are each anchored  */
                                        /* ... in the center of the split group         */
    int              nanions[eCompNR];  /* Requested number of anions and               */
    int              nAverage;          /* Coupling constant (nr of swap attempt steps) */
    real             threshold;         /* Ion counts may deviate from the requested
                                           values by +-threshold before a swap is done  */
    int              ncations[eCompNR]; /* ... cations for both compartments            */
    gmx_swapcoords_t si_priv;           /* swap private data accessible in
                                         * swapcoords.c                                 */
} t_swapcoords;


typedef struct {
    int      type;           /* type of AdResS simulation                    */
    gmx_bool bnew_wf;        /* enable new AdResS weighting function         */
    gmx_bool bchempot_dx;    /*true:interaction table format input is F=-dmu/dx   false: dmu_dwp  */
    gmx_bool btf_full_box;   /* true: appy therm force everywhere in the box according to table false: only in hybrid region */
    real     const_wf;       /* value of weighting function for eAdressConst */
    real     ex_width;       /* center of the explicit zone                  */
    real     hy_width;       /* width of the hybrid zone                     */
    int      icor;           /* type of interface correction                 */
    int      site;           /* AdResS CG site location                      */
    rvec     refs;           /* Coordinates for AdResS reference             */
    real     ex_forcecap;    /* in the hybrid zone, cap forces large then this to adress_ex_forcecap */
    gmx_bool do_hybridpairs; /* If true pair interaction forces are also scaled in an adress way*/

    int    * tf_table_index; /* contains mapping of energy group index -> i-th adress tf table*/
    int      n_tf_grps;
    int     *group_explicit;
    int      n_energy_grps;
} t_adress;

typedef struct {
    int             eI;                      /* Integration method                 */
    gmx_int64_t     nsteps;                  /* number of steps to be taken			*/
    int             simulation_part;         /* Used in checkpointing to separate chunks */
    gmx_int64_t     init_step;               /* start at a stepcount >0 (used w. convert-tpr)    */
    int             nstcalcenergy;           /* frequency of energy calc. and T/P coupl. upd.	*/
    int             cutoff_scheme;           /* group or verlet cutoffs     */
    int             ns_type;                 /* which ns method should we use?               */
    int             nstlist;                 /* number of steps before pairlist is generated	*/
    int             ndelta;                  /* number of cells per rlong			*/
    int             nstcomm;                 /* number of steps after which center of mass	*/
                                             /* motion is removed				*/
    int             comm_mode;               /* Center of mass motion removal algorithm      */
    int             nstlog;                  /* number of steps after which print to logfile	*/
    int             nstxout;                 /* number of steps after which X is output	*/
    int             nstvout;                 /* id. for V					*/
    int             nstfout;                 /* id. for F					*/
    int             nstenergy;               /* number of steps after which energies printed */
    int             nstxout_compressed;      /* id. for compressed trj (.xtc,.tng)           */
    double          init_t;                  /* initial time (ps)              */
    double          delta_t;                 /* time step (ps)				*/
    real            x_compression_precision; /* precision of x in compressed trajectory file */
    real            fourier_spacing;         /* requested fourier_spacing, when nk? not set  */
    int             nkx, nky, nkz;           /* number of k vectors in each spatial dimension*/
                                             /* for fourier methods for long range electrost.*/
    int             pme_order;               /* interpolation order for PME                  */
    real            ewald_rtol;              /* Real space tolerance for Ewald, determines   */
                                             /* the real/reciprocal space relative weight    */
    real            ewald_rtol_lj;           /* Real space tolerance for LJ-Ewald            */
    int             ewald_geometry;          /* normal/3d ewald, or pseudo-2d LR corrections */
    real            epsilon_surface;         /* Epsilon for PME dipole correction            */
    int             ljpme_combination_rule;  /* Type of combination rule in LJ-PME          */
    int             ePBC;                    /* Type of periodic boundary conditions		*/
    int             bPeriodicMols;           /* Periodic molecules                           */
    gmx_bool        bContinuation;           /* Continuation run: starting state is correct	*/
    int             etc;                     /* temperature coupling               */
    int             nsttcouple;              /* interval in steps for temperature coupling   */
    gmx_bool        bPrintNHChains;          /* whether to print nose-hoover chains        */
    int             epc;                     /* pressure coupling                            */
    int             epct;                    /* pressure coupling type			*/
    int             nstpcouple;              /* interval in steps for pressure coupling      */
    real            tau_p;                   /* pressure coupling time (ps)			*/
    tensor          ref_p;                   /* reference pressure (kJ/(mol nm^3))		*/
    tensor          compress;                /* compressability ((mol nm^3)/kJ)        */
    int             refcoord_scaling;        /* How to scale absolute reference coordinates  */
    rvec            posres_com;              /* The COM of the posres atoms                  */
    rvec            posres_comB;             /* The B-state COM of the posres atoms          */
    int             andersen_seed;           /* Random seed for Andersen thermostat (obsolete) */
    real            verletbuf_tol;           /* Per atom pair energy drift tolerance (kJ/mol/ps/atom) for list buffer  */
    real            rlist;                   /* short range pairlist cut-off (nm)		*/
    real            rlistlong;               /* long range pairlist cut-off (nm)		*/
    int             nstcalclr;               /* Frequency of evaluating direct space long-range interactions */
    real            rtpi;                    /* Radius for test particle insertion           */
    int             coulombtype;             /* Type of electrostatics treatment             */
    int             coulomb_modifier;        /* Modify the Coulomb interaction              */
    real            rcoulomb_switch;         /* Coulomb switch range start (nm)		*/
    real            rcoulomb;                /* Coulomb cutoff (nm)		                */
    real            epsilon_r;               /* relative dielectric constant                 */
    real            epsilon_rf;              /* relative dielectric constant of the RF       */
    int             implicit_solvent;        /* No (=explicit water), or GBSA solvent models */
    int             gb_algorithm;            /* Algorithm to use for calculation Born radii  */
    int             nstgbradii;              /* Frequency of updating Generalized Born radii */
    real            rgbradii;                /* Cutoff for GB radii calculation              */
    real            gb_saltconc;             /* Salt concentration (M) for GBSA models       */
    real            gb_epsilon_solvent;      /* dielectric coeff. of implicit solvent     */
    real            gb_obc_alpha;            /* 1st scaling factor for Bashford-Case GB      */
    real            gb_obc_beta;             /* 2nd scaling factor for Bashford-Case GB      */
    real            gb_obc_gamma;            /* 3rd scaling factor for Bashford-Case GB      */
    real            gb_dielectric_offset;    /* Dielectric offset for Still/HCT/OBC     */
    int             sa_algorithm;            /* Algorithm for SA part of GBSA                */
    real            sa_surface_tension;      /* Energy factor for SA part of GBSA */
    int             vdwtype;                 /* Type of Van der Waals treatment              */
    int             vdw_modifier;            /* Modify the VdW interaction                   */
    real            rvdw_switch;             /* Van der Waals switch range start (nm)        */
    real            rvdw;                    /* Van der Waals cutoff (nm)	        */
    int             eDispCorr;               /* Perform Long range dispersion corrections    */
    real            tabext;                  /* Extension of the table beyond the cut-off,   *
                                              * as well as the table length for 1-4 interac. */
    real            shake_tol;               /* tolerance for shake				*/
    int             efep;                    /* free energy calculations                     */
    t_lambda       *fepvals;                 /* Data for the FEP state                       */
    gmx_bool        bSimTemp;                /* Whether to do simulated tempering            */
    t_simtemp      *simtempvals;             /* Variables for simulated tempering            */
    gmx_bool        bExpanded;               /* Whether expanded ensembles are used          */
    t_expanded     *expandedvals;            /* Expanded ensemble parameters              */
    int             eDisre;                  /* Type of distance restraining                 */
    real            dr_fc;                   /* force constant for ta_disre			*/
    int             eDisreWeighting;         /* type of weighting of pairs in one restraints	*/
    gmx_bool        bDisreMixed;             /* Use comb of time averaged and instan. viol's	*/
    int             nstdisreout;             /* frequency of writing pair distances to enx   */
    real            dr_tau;                  /* time constant for memory function in disres    */
    real            orires_fc;               /* force constant for orientational restraints  */
    real            orires_tau;              /* time constant for memory function in orires    */
    int             nstorireout;             /* frequency of writing tr(SD) to enx           */
    real            em_stepsize;             /* The stepsize for updating			*/
    real            em_tol;                  /* The tolerance				*/
    int             niter;                   /* Number of iterations for convergence of      */
                                             /* steepest descent in relax_shells             */
    real            fc_stepsize;             /* Stepsize for directional minimization        */
                                             /* in relax_shells                              */
    int             nstcgsteep;              /* number of steps after which a steepest       */
                                             /* descents step is done while doing cg         */
    int             nbfgscorr;               /* Number of corrections to the hessian to keep */
    int             eConstrAlg;              /* Type of constraint algorithm                 */
    int             nProjOrder;              /* Order of the LINCS Projection Algorithm      */
    real            LincsWarnAngle;          /* If bond rotates more than %g degrees, warn   */
    int             nLincsIter;              /* Number of iterations in the final Lincs step */
    gmx_bool        bShakeSOR;               /* Use successive overrelaxation for shake      */
    real            bd_fric;                 /* Friction coefficient for BD (amu/ps)         */
    gmx_int64_t     ld_seed;                 /* Random seed for SD and BD                    */
    int             nwall;                   /* The number of walls                          */
    int             wall_type;               /* The type of walls                            */
    real            wall_r_linpot;           /* The potentail is linear for r<=wall_r_linpot */
    int             wall_atomtype[2];        /* The atom type for walls                      */
    real            wall_density[2];         /* Number density for walls                     */
    real            wall_ewald_zfac;         /* Scaling factor for the box for Ewald         */
    gmx_bool        bPull;                   /* Do we do COM pulling?                        */
    pull_params_t  *pull;                    /* The data for center of mass pulling          */
    struct pull_t  *pull_work;               /* The COM pull force calculation data structure; TODO this pointer should live somewhere else */

    gmx_bool        bRot;                    /* Calculate enforced rotation potential(s)?    */
    t_rot          *rot;                     /* The data for enforced rotation potentials    */
    int             eSwapCoords;             /* Do ion/water position exchanges (CompEL)?    */
    t_swapcoords   *swap;
    gmx_bool        bIMD;                    /* Allow interactive MD sessions for this .tpr? */
    t_IMD          *imd;                     /* Interactive molecular dynamics               */
    real            cos_accel;               /* Acceleration for viscosity calculation       */
    tensor          deform;                  /* Triclinic deformation velocities (nm/ps)     */
    int             userint1;                /* User determined parameters                   */
    int             userint2;
    int             userint3;
    int             userint4;
    real            userreal1;
    real            userreal2;
    real            userreal3;
    real            userreal4;
    t_grpopts       opts;          /* Group options				*/
    t_cosines       ex[DIM];       /* Electric field stuff	(spatial part)		*/
    t_cosines       et[DIM];       /* Electric field stuff	(time part)		*/
    gmx_bool        bQMMM;         /* QM/MM calculation                            */
    int             QMconstraints; /* constraints on QM bonds                      */
    int             QMMMscheme;    /* Scheme: ONIOM or normal                      */
    real            scalefactor;   /* factor for scaling the MM charges in QM calc.*/
                                   /* parameter needed for AdResS simulation       */
    gmx_bool        bAdress;       /* Is AdResS enabled ? */
    t_adress       *adress;        /* The data for adress simulations */
} t_inputrec;

#define DEFORM(ir) ((ir).deform[XX][XX] != 0 || (ir).deform[YY][YY] != 0 || (ir).deform[ZZ][ZZ] != 0 || (ir).deform[YY][XX] != 0 || (ir).deform[ZZ][XX] != 0 || (ir).deform[ZZ][YY] != 0)

#define DYNAMIC_BOX(ir) ((ir).epc != epcNO || (ir).eI == eiTPI || DEFORM(ir))

#define PRESERVE_SHAPE(ir) ((ir).epc != epcNO && (ir).deform[XX][XX] == 0 && ((ir).epct == epctISOTROPIC || (ir).epct == epctSEMIISOTROPIC))

#define NEED_MUTOT(ir) (((ir).coulombtype == eelEWALD || EEL_PME((ir).coulombtype)) && ((ir).ewald_geometry == eewg3DC || (ir).epsilon_surface != 0))

#define IR_TWINRANGE(ir) ((ir).rlist > 0 && ((ir).rlistlong == 0 || (ir).rlistlong > (ir).rlist))

#define IR_ELEC_FIELD(ir) ((ir).ex[XX].n > 0 || (ir).ex[YY].n > 0 || (ir).ex[ZZ].n > 0)

#define IR_EXCL_FORCES(ir) (EEL_FULL((ir).coulombtype) || (EEL_RF((ir).coulombtype) && (ir).coulombtype != eelRF_NEC) || (ir).implicit_solvent != eisNO)
/* use pointer definitions of ir here, since that's what's usually used in the code */
#define IR_NPT_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && (((ir)->epc == epcMTTK) && ((ir)->etc == etcNOSEHOOVER)))

#define IR_NVT_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && ((!((ir)->epc == epcMTTK)) && ((ir)->etc == etcNOSEHOOVER)))

#define IR_NPH_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && (((ir)->epc == epcMTTK) && (!(((ir)->etc == etcNOSEHOOVER)))))

#ifdef __cplusplus
}
#endif


#endif
