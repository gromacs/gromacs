/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_MDTYPES_INPUTREC_H
#define GMX_MDTYPES_INPUTREC_H

#include <cstdio>

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#define EGP_EXCL (1 << 0)
#define EGP_TABLE (1 << 1)

struct gmx_enfrot;
struct gmx_enfrotgrp;
struct pull_params_t;

namespace gmx
{
class Awh;
struct AwhParams;
class KeyValueTreeObject;
} // namespace gmx

struct t_grpopts
{
    //! Number of T-Coupl groups
    int ngtc;
    //! Number of of Nose-Hoover chains per group
    int nhchainlength;
    //! Number of Accelerate groups
    int ngacc;
    //! Number of Freeze groups
    int ngfrz;
    //! Number of Energy groups
    int ngener;
    //! Number of degrees of freedom in a group
    real* nrdf;
    //! Coupling temperature	per group
    real* ref_t;
    //! No/simple/periodic simulated annealing for each group
    int* annealing;
    //! Number of annealing time points per group
    int* anneal_npoints;
    //! For each group: Time points
    real** anneal_time;
    //! For each group: Temperature at these times. Final temp after all intervals is ref_t
    real** anneal_temp;
    //! Tau coupling time
    real* tau_t;
    //! Acceleration per group
    rvec* acc;
    //! Whether the group will be frozen in each direction
    ivec* nFreeze;
    //! Exclusions/tables of energy group pairs
    int* egp_flags;

    /* QMMM stuff */
    //! Number of QM groups
    int ngQM;
    //! Level of theory in the QM calculation
    int* QMmethod;
    //! Basisset in the QM calculation
    int* QMbasis;
    //! Total charge in the QM region
    int* QMcharge;
    //! Spin multiplicicty in the QM region
    int* QMmult;
    //! Surface hopping (diabatic hop only)
    gmx_bool* bSH;
    //! Number of orbiatls in the active space
    int* CASorbitals;
    //! Number of electrons in the active space
    int* CASelectrons;
    //! At which gap (A.U.) the SA is switched on
    real* SAon;
    //! At which gap (A.U.) the SA is switched off
    real* SAoff;
    //! In how many steps SA goes from 1-1 to 0.5-0.5
    int* SAsteps;
};

struct t_simtemp
{
    //! Simulated temperature scaling; linear or exponential
    int eSimTempScale;
    //! The low temperature for simulated tempering
    real simtemp_low;
    //! The high temperature for simulated tempering
    real simtemp_high;
    //! The range of temperatures used for simulated tempering
    real* temperatures;
};

struct t_lambda
{
    //! The frequency for calculating dhdl
    int nstdhdl;
    //! Fractional value of lambda (usually will use init_fep_state, this will only be for slow growth, and for legacy free energy code. Only has a valid value if positive)
    double init_lambda;
    //! The initial number of the state
    int init_fep_state;
    //! Change of lambda per time step (fraction of (0.1)
    double delta_lambda;
    //! Print no, total or potential energies in dhdl
    int edHdLPrintEnergy;
    //! The number of foreign lambda points
    int n_lambda;
    //! The array of all lambda values
    double** all_lambda;
    //! The number of neighboring lambda states to calculate the energy for in up and down directions (-1 for all)
    int lambda_neighbors;
    //! The first lambda to calculate energies for
    int lambda_start_n;
    //! The last lambda +1 to calculate energies for
    int lambda_stop_n;
    //! Free energy soft-core parameter
    real sc_alpha;
    //! Lambda power for soft-core interactions
    int sc_power;
    //! R power for soft-core interactions
    real sc_r_power;
    //! Free energy soft-core sigma when c6 or c12=0
    real sc_sigma;
    //! Free energy soft-core sigma for ?????
    real sc_sigma_min;
    //! Use softcore for the coulomb portion as well (default FALSE)
    gmx_bool bScCoul;
    //! Whether to print the dvdl term associated with this term; if it is not specified as separate, it is lumped with the FEP term
    gmx_bool separate_dvdl[efptNR];
    //! Whether to write a separate dhdl.xvg file note: NOT a gmx_bool, but an enum
    int separate_dhdl_file;
    //! Whether to calculate+write dhdl derivatives note: NOT a gmx_bool, but an enum
    int dhdl_derivatives;
    //! The maximum table size for the dH histogram
    int dh_hist_size;
    //! The spacing for the dH histogram
    double dh_hist_spacing;
};

struct t_expanded
{
    //! The frequency of expanded ensemble state changes
    int nstexpanded;
    //! Which type of move updating do we use for lambda monte carlo (or no for none)
    int elamstats;
    //! What move set will be we using for state space moves
    int elmcmove;
    //! The method we use to decide of we have equilibrated the weights
    int elmceq;
    //! The minumum number of samples at each lambda for deciding whether we have reached a minimum
    int equil_n_at_lam;
    //! Wang-Landau delta at which we stop equilibrating weights
    real equil_wl_delta;
    //! Use the ratio of weights (ratio of minimum to maximum) to decide when to stop equilibrating
    real equil_ratio;
    //! After equil_steps steps we stop equilibrating the weights
    int equil_steps;
    //! After equil_samples total samples (steps/nstfep), we stop equilibrating the weights
    int equil_samples;
    //! Random number seed for lambda mc switches
    int lmc_seed;
    //! Whether to use minumum variance weighting
    gmx_bool minvar;
    //! The number of samples needed before kicking into minvar routine
    int minvarmin;
    //! The offset for the variance in MinVar
    real minvar_const;
    //! Range of cvalues used for BAR
    int c_range;
    //! Whether to print symmetrized matrices
    gmx_bool bSymmetrizedTMatrix;
    //! How frequently to print the transition matrices
    int nstTij;
    //! Number of repetitions in the MC lambda jumps MRS -- VERIFY THIS
    int lmc_repeats;
    //! Minimum number of samples for each state before free sampling MRS -- VERIFY THIS!
    int lmc_forced_nstart;
    //! Distance in lambda space for the gibbs interval
    int gibbsdeltalam;
    //! Scaling factor for Wang-Landau
    real wl_scale;
    //! Ratio between largest and smallest number for freezing the weights
    real wl_ratio;
    //! Starting delta for Wang-Landau
    real init_wl_delta;
    //! Use one over t convergence for Wang-Landau when the delta get sufficiently small
    gmx_bool bWLoneovert;
    //! Did we initialize the weights? TODO: REMOVE FOR 5.0, no longer needed with new logic
    gmx_bool bInit_weights;
    //! To override the main temperature, or define it if it's not defined
    real mc_temp;
    //! User-specified initial weights to start with
    real* init_lambda_weights;
};

struct t_rotgrp
{
    //! Rotation type for this group
    int eType;
    //! Use mass-weighed positions?
    int bMassW;
    //! Number of atoms in the group
    int nat;
    //! The global atoms numbers
    int* ind;
    //! The reference positions
    rvec* x_ref;
    //! The normalized rotation vector
    rvec inputVec;
    //! Rate of rotation (degree/ps)
    real rate;
    //! Force constant (kJ/(mol nm^2)
    real k;
    //! Pivot point of rotation axis (nm)
    rvec pivot;
    //! Type of fit to determine actual group angle
    int eFittype;
    //! Number of angles around the reference angle for which the rotation potential is also evaluated (for fit type 'potential' only)
    int PotAngle_nstep;
    //! Distance between two angles in degrees (for fit type 'potential' only)
    real PotAngle_step;
    //! Slab distance (nm)
    real slab_dist;
    //! Minimum value the gaussian must have so that the force is actually evaluated
    real min_gaussian;
    //! Additive constant for radial motion2 and flexible2 potentials (nm^2)
    real eps;
};

struct t_rot
{
    //! Number of rotation groups
    int ngrp;
    //! Output frequency for main rotation outfile
    int nstrout;
    //! Output frequency for per-slab data
    int nstsout;
    //! Groups to rotate
    t_rotgrp* grp;
};

struct t_IMD
{
    //! Number of interactive atoms
    int nat;
    //! The global indices of the interactive atoms
    int* ind;
};

struct t_swapGroup
{
    //! Name of the swap group, e.g. NA, CL, SOL
    char* molname;
    //! Number of atoms in this group
    int nat;
    //! The global ion group atoms numbers
    int* ind;
    //! Requested number of molecules of this type per compartment
    int nmolReq[eCompNR];
};

struct t_swapcoords
{
    //! Period between when a swap is attempted
    int nstswap;
    //! Use mass-weighted positions in split group
    gmx_bool massw_split[2];
    /*! \brief Split cylinders defined by radius, upper and lower
     * extension. The split cylinders define the channels and are
     * each anchored in the center of the split group */
    /**@{*/
    real cyl0r, cyl1r;
    real cyl0u, cyl1u;
    real cyl0l, cyl1l;
    /**@}*/
    //! Coupling constant (number of swap attempt steps)
    int nAverage;
    //! Ion counts may deviate from the requested values by +-threshold before a swap is done
    real threshold;
    //! Offset of the swap layer (='bulk') with respect to the compartment-defining layers
    real bulkOffset[eCompNR];
    //! Number of groups to be controlled
    int ngrp;
    //! All swap groups, including split and solvent
    t_swapGroup* grp;
};

struct t_inputrec // NOLINT (clang-analyzer-optin.performance.Padding)
{
    t_inputrec();
    explicit t_inputrec(const t_inputrec&) = delete;
    t_inputrec& operator=(const t_inputrec&) = delete;
    ~t_inputrec();

    //! Integration method
    int eI;
    //! Number of steps to be taken
    int64_t nsteps;
    //! Used in checkpointing to separate chunks
    int simulation_part;
    //! Start at a stepcount >0 (used w. convert-tpr)
    int64_t init_step;
    //! Frequency of energy calc. and T/P coupl. upd.
    int nstcalcenergy;
    //! Group or verlet cutoffs
    int cutoff_scheme;
    //! Number of steps before pairlist is generated
    int nstlist;
    //! Number of steps after which center of mass motion is removed
    int nstcomm;
    //! Center of mass motion removal algorithm
    int comm_mode;
    //! Number of steps after which print to logfile
    int nstlog;
    //! Number of steps after which X is output
    int nstxout;
    //! Number of steps after which V is output
    int nstvout;
    //! Number of steps after which F is output
    int nstfout;
    //! Number of steps after which energies printed
    int nstenergy;
    //! Number of steps after which compressed trj (.xtc,.tng) is output
    int nstxout_compressed;
    //! Initial time (ps)
    double init_t;
    //! Time step (ps)
    double delta_t;
    //! Precision of x in compressed trajectory file
    real x_compression_precision;
    //! Requested fourier_spacing, when nk? not set
    real fourier_spacing;
    //! Number of k vectors in x dimension for fourier methods for long range electrost.
    int nkx;
    //! Number of k vectors in y dimension for fourier methods for long range electrost.
    int nky;
    //! Number of k vectors in z dimension for fourier methods for long range electrost.
    int nkz;
    //! Interpolation order for PME
    int pme_order;
    //! Real space tolerance for Ewald, determines the real/reciprocal space relative weight
    real ewald_rtol;
    //! Real space tolerance for LJ-Ewald
    real ewald_rtol_lj;
    //! Normal/3D ewald, or pseudo-2D LR corrections
    int ewald_geometry;
    //! Epsilon for PME dipole correction
    real epsilon_surface;
    //! Type of combination rule in LJ-PME
    int ljpme_combination_rule;
    //! Type of periodic boundary conditions
    int ePBC;
    //! Periodic molecules
    bool bPeriodicMols;
    //! Continuation run: starting state is correct (ie. constrained)
    gmx_bool bContinuation;
    //! Temperature coupling
    int etc;
    //! Interval in steps for temperature coupling
    int nsttcouple;
    //! Whether to print nose-hoover chains
    gmx_bool bPrintNHChains;
    //! Pressure coupling
    int epc;
    //! Pressure coupling type
    int epct;
    //! Interval in steps for pressure coupling
    int nstpcouple;
    //! Pressure coupling time (ps)
    real tau_p;
    //! Reference pressure (kJ/(mol nm^3))
    tensor ref_p;
    //! Compressibility ((mol nm^3)/kJ)
    tensor compress;
    //! How to scale absolute reference coordinates
    int refcoord_scaling;
    //! The COM of the posres atoms
    rvec posres_com;
    //! The B-state COM of the posres atoms
    rvec posres_comB;
    //! Random seed for Andersen thermostat (obsolete)
    int andersen_seed;
    //! Per atom pair energy drift tolerance (kJ/mol/ps/atom) for list buffer
    real verletbuf_tol;
    //! Short range pairlist cut-off (nm)
    real rlist;
    //! Radius for test particle insertion
    real rtpi;
    //! Type of electrostatics treatment
    int coulombtype;
    //! Modify the Coulomb interaction
    int coulomb_modifier;
    //! Coulomb switch range start (nm)
    real rcoulomb_switch;
    //! Coulomb cutoff (nm)
    real rcoulomb;
    //! Relative dielectric constant
    real epsilon_r;
    //! Relative dielectric constant of the RF
    real epsilon_rf;
    //! Always false (no longer supported)
    bool implicit_solvent;
    //! Type of Van der Waals treatment
    int vdwtype;
    //! Modify the Van der Waals interaction
    int vdw_modifier;
    //! Van der Waals switch range start (nm)
    real rvdw_switch;
    //! Van der Waals cutoff (nm)
    real rvdw;
    //! Perform Long range dispersion corrections
    int eDispCorr;
    //! Extension of the table beyond the cut-off, as well as the table length for 1-4 interac.
    real tabext;
    //! Tolerance for shake
    real shake_tol;
    //! Free energy calculations
    int efep;
    //! Data for the FEP state
    t_lambda* fepvals;
    //! Whether to do simulated tempering
    gmx_bool bSimTemp;
    //! Variables for simulated tempering
    t_simtemp* simtempvals;
    //! Whether expanded ensembles are used
    gmx_bool bExpanded;
    //! Expanded ensemble parameters
    t_expanded* expandedvals;
    //! Type of distance restraining
    int eDisre;
    //! Force constant for time averaged distance restraints
    real dr_fc;
    //! Type of weighting of pairs in one restraints
    int eDisreWeighting;
    //! Use combination of time averaged and instantaneous violations
    gmx_bool bDisreMixed;
    //! Frequency of writing pair distances to enx
    int nstdisreout;
    //! Time constant for memory function in disres
    real dr_tau;
    //! Force constant for orientational restraints
    real orires_fc;
    //! Time constant for memory function in orires
    real orires_tau;
    //! Frequency of writing tr(SD) to energy output
    int nstorireout;
    //! The stepsize for updating
    real em_stepsize;
    //! The tolerance
    real em_tol;
    //! Number of iterations for convergence of steepest descent in relax_shells
    int niter;
    //! Stepsize for directional minimization in relax_shells
    real fc_stepsize;
    //! Number of steps after which a steepest descents step is done while doing cg
    int nstcgsteep;
    //! Number of corrections to the Hessian to keep
    int nbfgscorr;
    //! Type of constraint algorithm
    int eConstrAlg;
    //! Order of the LINCS Projection Algorithm
    int nProjOrder;
    //! Warn if any bond rotates more than this many degrees
    real LincsWarnAngle;
    //! Number of iterations in the final LINCS step
    int nLincsIter;
    //! Use successive overrelaxation for shake
    gmx_bool bShakeSOR;
    //! Friction coefficient for BD (amu/ps)
    real bd_fric;
    //! Random seed for SD and BD
    int64_t ld_seed;
    //! The number of walls
    int nwall;
    //! The type of walls
    int wall_type;
    //! The potentail is linear for r<=wall_r_linpot
    real wall_r_linpot;
    //! The atom type for walls
    int wall_atomtype[2];
    //! Number density for walls
    real wall_density[2];
    //! Scaling factor for the box for Ewald
    real wall_ewald_zfac;

    /* COM pulling data */
    //! Do we do COM pulling?
    gmx_bool bPull;
    //! The data for center of mass pulling
    pull_params_t* pull;

    /* AWH bias data */
    //! Whether to use AWH biasing for PMF calculations
    gmx_bool bDoAwh;
    //! AWH biasing parameters
    gmx::AwhParams* awhParams;

    /* Enforced rotation data */
    //! Whether to calculate enforced rotation potential(s)
    gmx_bool bRot;
    //! The data for enforced rotation potentials
    t_rot* rot;

    //! Whether to do ion/water position exchanges (CompEL)
    int eSwapCoords;
    //! Swap data structure.
    t_swapcoords* swap;

    //! Whether the tpr makes an interactive MD session possible.
    gmx_bool bIMD;
    //! Interactive molecular dynamics
    t_IMD* imd;

    //! Acceleration for viscosity calculation
    real cos_accel;
    //! Triclinic deformation velocities (nm/ps)
    tensor deform;
    /*! \brief User determined parameters */
    /**@{*/
    int  userint1;
    int  userint2;
    int  userint3;
    int  userint4;
    real userreal1;
    real userreal2;
    real userreal3;
    real userreal4;
    /**@}*/
    //! Group options
    t_grpopts opts;
    //! QM/MM calculation
    gmx_bool bQMMM;
    //! Constraints on QM bonds
    int QMconstraints;
    //! Scheme: ONIOM or normal
    int QMMMscheme;
    //! Factor for scaling the MM charges in QM calc.
    real scalefactor;

    /* Fields for removed features go here (better caching) */
    //! Whether AdResS is enabled - always false if a valid .tpr was read
    gmx_bool bAdress;
    //! Whether twin-range scheme is active - always false if a valid .tpr was read
    gmx_bool useTwinRange;

    //! KVT object that contains input parameters converted to the new style.
    gmx::KeyValueTreeObject* params;

    //! KVT for storing simulation parameters that are not part of the mdp file.
    std::unique_ptr<gmx::KeyValueTreeObject> internalParameters;
};

int ir_optimal_nstcalcenergy(const t_inputrec* ir);

int tcouple_min_integration_steps(int etc);

int ir_optimal_nsttcouple(const t_inputrec* ir);

int pcouple_min_integration_steps(int epc);

int ir_optimal_nstpcouple(const t_inputrec* ir);

/* Returns if the Coulomb force or potential is switched to zero */
gmx_bool ir_coulomb_switched(const t_inputrec* ir);

/* Returns if the Coulomb interactions are zero beyond the rcoulomb.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec* ir);

/* As ir_coulomb_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rcoulomb.
 */
gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec* ir);

/* Returns if the Van der Waals force or potential is switched to zero */
gmx_bool ir_vdw_switched(const t_inputrec* ir);

/* Returns if the Van der Waals interactions are zero beyond the rvdw.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec* ir);

/* As ir_vdw_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rvdw.
 */
gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec* ir);

/*! \brief Free memory from input record.
 *
 * All arrays and pointers will be freed.
 *
 * \param[in] ir The data structure
 */
void done_inputrec(t_inputrec* ir);

void pr_inputrec(FILE* fp, int indent, const char* title, const t_inputrec* ir, gmx_bool bMDPformat);

void cmp_inputrec(FILE* fp, const t_inputrec* ir1, const t_inputrec* ir2, real ftol, real abstol);

void comp_pull_AB(FILE* fp, pull_params_t* pull, real ftol, real abstol);


gmx_bool inputrecDeform(const t_inputrec* ir);

gmx_bool inputrecDynamicBox(const t_inputrec* ir);

gmx_bool inputrecPreserveShape(const t_inputrec* ir);

gmx_bool inputrecNeedMutot(const t_inputrec* ir);

gmx_bool inputrecTwinRange(const t_inputrec* ir);

gmx_bool inputrecExclForces(const t_inputrec* ir);

gmx_bool inputrecNptTrotter(const t_inputrec* ir);

gmx_bool inputrecNvtTrotter(const t_inputrec* ir);

gmx_bool inputrecNphTrotter(const t_inputrec* ir);

/*! \brief Return true if the simulation is 2D periodic with two walls. */
bool inputrecPbcXY2Walls(const t_inputrec* ir);

/*! \brief Returns true for MD integator with T and/or P-coupling that supports
 * calculating the conserved energy quantity.
 */
bool integratorHasConservedEnergyQuantity(const t_inputrec* ir);

/*! \brief Returns true when temperature is coupled or constant. */
bool integratorHasReferenceTemperature(const t_inputrec* ir);

/*! \brief Return the number of bounded dimensions
 *
 * \param[in] ir The input record with MD parameters
 * \return the number of dimensions in which
 * the coordinates of the particles are bounded, starting at X.
 */
int inputrec2nboundeddim(const t_inputrec* ir);

/*! \brief Returns the number of degrees of freedom in center of mass motion
 *
 * \param[in] ir  The inputrec structure
 * \return the number of degrees of freedom of the center of mass
 */
int ndof_com(const t_inputrec* ir);

/*! \brief Returns the maximum reference temperature over all coupled groups
 *
 * Returns 0 for energy minimization and normal mode computation.
 * Returns -1 for MD without temperature coupling.
 *
 * \param[in] ir  The inputrec structure
 */
real maxReferenceTemperature(const t_inputrec& ir);

/*! \brief Returns whether there is an Ewald surface contribution
 */
bool haveEwaldSurfaceContribution(const t_inputrec& ir);

#endif /* GMX_MDTYPES_INPUTREC_H */
