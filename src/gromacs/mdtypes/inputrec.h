/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#define EGP_EXCL  (1<<0)
#define EGP_TABLE (1<<1)

struct pull_params_t;
struct pull_t;

namespace gmx
{
class Awh;
struct AwhParams;
class KeyValueTreeObject;
}

struct t_grpopts
{
    //! # T-Coupl groups
    int    ngtc = 0;
    //! # of nose-hoover chains per group
    int    nhchainlength = 0;
    //! # Accelerate groups
    int    ngacc = 0;
    //! # Freeze groups
    int    ngfrz = 0;
    //! # Ener groups
    int    ngener = 0;
    //! Nr of degrees of freedom in a group
    real  *nrdf = nullptr;
    //! Coupling temperature	per group
    real  *ref_t = nullptr;
    //! No/simple/periodic SA for each group
    int   *annealing = nullptr;
    //! Number of annealing time points per grp
    int   *anneal_npoints = nullptr;
    //! For ea. group: Time points
    real **anneal_time = nullptr;
    //! For ea. grp: Temperature at these times. Final temp after all intervals is ref_t
    real **anneal_temp = nullptr;
    //! Tau coupling time
    real  *tau_t = nullptr;
    //! Acceleration per group
    rvec  *acc = nullptr;
    //! Freeze the group in each direction ?
    ivec  *nFreeze = nullptr;
    //! Exclusions/tables of energy group pairs
    int   *egp_flags = nullptr;

    /* QMMM stuff */
    //! nr of QM groups
    int       ngQM = 0;
    //! Level of theory in the QM calculation
    int      *QMmethod = nullptr;
    //! Basisset in the QM calculation
    int      *QMbasis = nullptr;
    //! Total charge in the QM region
    int      *QMcharge = nullptr;
    //! Spin multiplicicty in the QM region
    int      *QMmult = nullptr;
    //! surface hopping (diabatic hop only)
    gmx_bool *bSH = nullptr;
    //! number of orbiatls in the active space
    int      *CASorbitals = nullptr;
    //! number of electrons in the active space
    int      *CASelectrons = nullptr;
    //! at which gap (A.U.) the SA is switched on
    real     *SAon = nullptr;
    //! at which gap (A.U.) the SA is switched off
    real     *SAoff = nullptr;
    //! in how many steps SA goes from 1-1 to 0.5-0.5
    int      *SAsteps = nullptr;
};

struct t_simtemp
{
    //! simulated temperature scaling; linear or exponential
    int   eSimTempScale = 0;
    //! the low temperature for simulated tempering
    real  simtemp_low = 0;
    //! the high temperature for simulated tempering
    real  simtemp_high = 0;
    //! the range of temperatures used for simulated tempering
    real *temperatures = nullptr;
};

struct t_lambda
{
    //! The frequency for calculating dhdl
    int      nstdhdl = 0;
    //! fractional value of lambda (usually will use init_fep_state, this will only be for slow growth, and for legacy free energy code. Only has a valid value if positive)
    double   init_lambda = 0;
    //! the initial number of the state
    int      init_fep_state = 0;
    //! change of lambda per time step (fraction of (0.1)
    double   delta_lambda = 0;
    //! print no, total or potential energies in dhdl
    int      edHdLPrintEnergy = 0;
    //! The number of foreign lambda points
    int      n_lambda = 0;
    //! The array of all lambda values
    double **all_lambda = nullptr;
    //! The number of neighboring lambda states to calculate the energy for in up and down directions (-1 for all)
    int      lambda_neighbors = 0;
    //! The first lambda to calculate energies for
    int      lambda_start_n = 0;
    //! The last lambda +1 to calculate energies for
    int      lambda_stop_n = 0;
    //! free energy soft-core parameter
    real     sc_alpha = 0;
    //! lambda power for soft-core interactions
    int      sc_power = 0;
    //! r power for soft-core interactions
    real     sc_r_power = 0;
    //! free energy soft-core sigma when c6 or c12=0
    real     sc_sigma = 0;
    //! free energy soft-core sigma for ?????
    real     sc_sigma_min = 0;
    //! use softcore for the coulomb portion as well (default FALSE)
    gmx_bool bScCoul = 0;
    //! whether to print the dvdl term associated with this term; if it is not specified as separate, it is lumped with the FEP term
    gmx_bool separate_dvdl[efptNR] = {0};
    //! whether to write a separate dhdl.xvg file note: NOT a gmx_bool, but an enum
    int      separate_dhdl_file = 0;
    //! whether to calculate+write dhdl derivatives note: NOT a gmx_bool, but an enum
    int      dhdl_derivatives = 0;
    //! The maximum table size for the dH histogram
    int      dh_hist_size = 0;
    //! The spacing for the dH histogram
    double   dh_hist_spacing = 0;
};

struct t_expanded {
    //! The frequency of expanded ensemble state changes
    int      nstexpanded = 0;
    //! which type of move updating do we use for lambda monte carlo (or no for none)
    int      elamstats = 0;
    //! what move set will be we using for state space moves
    int      elmcmove = 0;
    //! the method we use to decide of we have equilibrated the weights
    int      elmceq = 0;
    //! the minumum number of samples at each lambda for deciding whether we have reached a minimum
    int      equil_n_at_lam = 0;
    //! WL delta at which we stop equilibrating weights
    real     equil_wl_delta = 0;
    //! use the ratio of weights (ratio of minimum to maximum) to decide when to stop equilibrating
    real     equil_ratio = 0;
    //! after equil_steps steps we stop equilibrating the weights
    int      equil_steps = 0;
    //! after equil_samples total samples (steps/nstfep), we stop equilibrating the weights
    int      equil_samples = 0;
    //! random number seed for lambda mc switches
    int      lmc_seed = 0;
    //! whether to use minumum variance weighting
    gmx_bool minvar = 0;
    //! the number of samples needed before kicking into minvar routine
    int      minvarmin = 0;
    //! the offset for the variance in MinVar
    real     minvar_const = 0;
    //! range of cvalues used for BAR
    int      c_range = 0;
    //! whether to print symmetrized matrices
    gmx_bool bSymmetrizedTMatrix = 0;
    //! How frequently to print the transition matrices
    int      nstTij = 0;
    //! number of repetitions in the MC lambda jumps MRS -- VERIFY THIS
    int      lmc_repeats = 0;
    //! minimum number of samples for each state before free sampling MRS -- VERIFY THIS!
    int      lmc_forced_nstart = 0;
    //! distance in lambda space for the gibbs interval
    int      gibbsdeltalam = 0;
    //! scaling factor for wang-landau
    real     wl_scale = 0;
    //! ratio between largest and smallest number for freezing the weights
    real     wl_ratio = 0;
    //! starting delta for wang-landau
    real     init_wl_delta = 0;
    //! use one over t convergece for wang-landau when the delta get sufficiently small
    gmx_bool bWLoneovert = 0;
    //! did we initialize the weights? TODO: REMOVE FOR 5.0, no longer needed with new logic
    gmx_bool bInit_weights = 0;
    //! To override the main temperature, or define it if it's not defined
    real     mc_temp = 0;
    //! user-specified initial weights to start with
    real    *init_lambda_weights = nullptr;
};


/* Abstract types for enforced rotation only defined in pull_rotation.c       */
typedef struct gmx_enfrot *gmx_enfrot_t;
typedef struct gmx_enfrotgrp *gmx_enfrotgrp_t;

struct t_rotgrp
{
    //! Rotation type for this group
    int              eType = 0;
    //! Use mass-weighed positions?
    int              bMassW = 0;
    //! Number of atoms in the group
    int              nat = 0;
    //! The global atoms numbers
    std::vector<int> ind;
    //! The reference positions
    rvec            *x_ref = nullptr;
    //! The normalized rotation vector
    rvec             vec = {0};
    //! Rate of rotation (degree/ps)
    real             rate = 0;
    //! Force constant (kJ/(mol nm^2)
    real             k = 0;
    //! Pivot point of rotation axis (nm)
    rvec             pivot = {0};
    //! Type of fit to determine actual group angle
    int              eFittype = 0;
    //! Number of angles around the reference angle for which the rotation potential is also evaluated (for fit type 'potential' only)
    int              PotAngle_nstep = 0;
    //! Distance between two angles in degrees (for fit type 'potential' only)
    real             PotAngle_step = 0;
    //! Slab distance (nm)
    real             slab_dist = 0;
    //! Minimum value the gaussian must have so that the force is actually evaluated
    real             min_gaussian = 0;
    //! Additive constant for radial motion2 and flexible2 potentials (nm^2)
    real             eps = 0;
    //! Stores non-inputrec rotation data per group
    gmx_enfrotgrp_t  enfrotgrp = nullptr;
};

struct t_rot
{
    //! Number of rotation groups
    int                   ngrp = 0;
    //! Output frequency for main rotation outfile
    int                   nstrout = 0;
    //! Output frequency for per-slab data
    int                   nstsout = 0;
    //! Groups to rotate
    std::vector<t_rotgrp> grp;
    //! Stores non-inputrec enforced rotation data
    gmx_enfrot_t          enfrot = nullptr;
};

/* Abstract type for IMD only defined in IMD.c */
struct t_gmx_IMD;

struct t_IMD
{
    //! Number of interactive atoms
    int              nat = 0;
    //! The global indices of the interactive atoms
    std::vector<int> ind;
    struct //! Stores non-inputrec IMD data
    t_gmx_IMD       *setup = nullptr;
};

/* Abstract types for position swapping only defined in swapcoords.cpp */
typedef struct t_swap *gmx_swapcoords_t;

struct t_swapGroup
{
    //! Name of the swap group, e.g. NA, CL, SOL
    char            *molname = nullptr;
    //! Number of atoms in this group
    int              nat = 0;
    //! The global ion group atoms numbers
    std::vector<int> ind;
    //! Requested number of molecules of this type per compartment
    int              nmolReq[eCompNR] = {0};
};

struct t_swapcoords
{
    //! Every how many steps a swap is attempted?
    int      nstswap = 0;
    //! Use mass-weighted positions in split group?
    gmx_bool massw_split[2] = {0, 0};
    /*! \brief Split cylinders defined by radius, upper and lower
     * extension. The split cylinders define the channels and are
     * each anchored in the center of the split group */
    /**@{*/
    real cyl0r, cyl1r = 0;
    real cyl0u, cyl1u = 0;
    real cyl0l, cyl1l = 0;
    /**@}*/
    //! Coupling constant (nr of swap attempt steps)
    int                      nAverage = 0;
    //! Ion counts may deviate from the requested values by +-threshold before a swap is done
    real                     threshold = 0;
    //! Offset of the swap layer (='bulk') with respect to the compartment-defining layers
    real                     bulkOffset[eCompNR] = {0};
    //! Number of groups to be controlled
    int                      ngrp = 0;
    //! All swap groups, including split and solvent
    std::vector<t_swapGroup> grp;
    //! swap private data accessible in swapcoords.cpp
    gmx_swapcoords_t         si_priv = nullptr;
};

struct t_inputrec
{
    t_inputrec();
    explicit t_inputrec(const t_inputrec &) = delete;
    t_inputrec &operator=(const t_inputrec &) = delete;
    ~t_inputrec();

    //! Integration method
    int                         eI = 0;
    //! number of steps to be taken
    gmx_int64_t                 nsteps = 0;
    //! Used in checkpointing to separate chunks
    int                         simulation_part = 0;
    //! start at a stepcount >0 (used w. convert-tpr)
    gmx_int64_t                 init_step = 0;
    //! frequency of energy calc. and T/P coupl. upd.
    int                         nstcalcenergy = 0;
    //! group or verlet cutoffs
    int                         cutoff_scheme = 0;
    //! which ns method should we use?
    int                         ns_type = 0;
    //! number of steps before pairlist is generated
    int                         nstlist = 0;
    //! number of cells per rlong
    int                         ndelta = 0;
    //! number of steps after which center of mass motion is removed
    int                         nstcomm = 0;
    //! Center of mass motion removal algorithm
    int                         comm_mode = 0;
    //! number of steps after which print to logfile
    int                         nstlog = 0;
    //! number of steps after which X is output
    int                         nstxout = 0;
    //! number of steps after which V is output
    int                         nstvout = 0;
    //! number of steps after which F is output
    int                         nstfout = 0;
    //! number of steps after which energies printed
    int                         nstenergy = 0;
    //! number of steps after which compressed trj (.xtc,.tng) is output
    int                         nstxout_compressed = 0;
    //! initial time (ps)
    double                      init_t = 0;
    //! time step (ps)
    double                      delta_t = 0;
    //! precision of x in compressed trajectory file
    real                        x_compression_precision = 0;
    //! requested fourier_spacing, when nk? not set
    real                        fourier_spacing = 0;
    //! number of k vectors in x dimension for fourier methods for long range electrost.
    int                         nkx = 0;
    //! number of k vectors in y dimension for fourier methods for long range electrost.
    int                         nky = 0;
    //! number of k vectors in z dimension for fourier methods for long range electrost.
    int                         nkz = 0;
    //! interpolation order for PME
    int                         pme_order = 0;
    //! Real space tolerance for Ewald, determines the real/reciprocal space relative weight
    real                        ewald_rtol = 0;
    //! Real space tolerance for LJ-Ewald
    real                        ewald_rtol_lj = 0;
    //! normal/3d ewald, or pseudo-2d LR corrections
    int                         ewald_geometry = 0;
    //! Epsilon for PME dipole correction
    real                        epsilon_surface = 0;
    //! Type of combination rule in LJ-PME
    int                         ljpme_combination_rule = 0;
    //! Type of periodic boundary conditions
    int                         ePBC = 0;
    //! Periodic molecules
    int                         bPeriodicMols = 0;
    //! Continuation run: starting state is correct
    gmx_bool                    bContinuation = 0;
    //! temperature coupling
    int                         etc = 0;
    //! interval in steps for temperature coupling
    int                         nsttcouple = 0;
    //! whether to print nose-hoover chains
    gmx_bool                    bPrintNHChains = 0;
    //! pressure coupling
    int                         epc = 0;
    //! pressure coupling type
    int                         epct = 0;
    //! interval in steps for pressure coupling
    int                         nstpcouple = 0;
    //! pressure coupling time (ps)
    real                        tau_p = 0;
    //! reference pressure (kJ/(mol nm^3))
    tensor                      ref_p = {{0}};
    //! compressability ((mol nm^3)/kJ)
    tensor                      compress = {{0}};
    //! How to scale absolute reference coordinates
    int                         refcoord_scaling = 0;
    //! The COM of the posres atoms
    rvec                        posres_com = {0};
    //! The B-state COM of the posres atoms
    rvec                        posres_comB = {0};
    //! Random seed for Andersen thermostat (obsolete)
    int                         andersen_seed = 0;
    //! Per atom pair energy drift tolerance (kJ/mol/ps/atom) for list buffer
    real                        verletbuf_tol = 0;
    //! short range pairlist cut-off (nm)
    real                        rlist = 0;
    //! Radius for test particle insertion
    real                        rtpi = 0;
    //! Type of electrostatics treatment
    int                         coulombtype = 0;
    //! Modify the Coulomb interaction
    int                         coulomb_modifier = 0;
    //! Coulomb switch range start (nm)
    real                        rcoulomb_switch = 0;
    //! Coulomb cutoff (nm)
    real                        rcoulomb = 0;
    //! relative dielectric constant
    real                        epsilon_r = 0;
    //! relative dielectric constant of the RF
    real                        epsilon_rf = 0;
    //! Always false (no longer supported
    bool                        implicit_solvent = 0;
    //! Type of Van der Waals treatment
    int                         vdwtype = 0;
    //! Modify the VdW interaction
    int                         vdw_modifier = 0;
    //! Van der Waals switch range start (nm)
    real                        rvdw_switch = 0;
    //! Van der Waals cutoff (nm)
    real                        rvdw = 0;
    //! Perform Long range dispersion corrections
    int                         eDispCorr = 0;
    //! Extension of the table beyond the cut-off, as well as the table length for 1-4 interac.
    real                        tabext = 0;
    //! tolerance for shake
    real                        shake_tol = 0;
    //! free energy calculations
    int                         efep = 0;
    //! Data for the FEP state
    std::unique_ptr<t_lambda>   fepvals;
    //! Whether to do simulated tempering
    gmx_bool                    bSimTemp = 0;
    //! Variables for simulated tempering
    std::unique_ptr<t_simtemp>  simtempvals;
    //! Whether expanded ensembles are used
    gmx_bool                    bExpanded = 0;
    //! Expanded ensemble parameters
    std::unique_ptr<t_expanded> expandedvals;
    //! Type of distance restraining
    int                         eDisre = 0;
    //! force constant for ta_disre
    real                        dr_fc = 0;
    //! type of weighting of pairs in one restraints
    int                         eDisreWeighting = 0;
    //! Use combination of time averaged and instantaneous violations
    gmx_bool                    bDisreMixed = 0;
    //! frequency of writing pair distances to enx
    int                         nstdisreout = 0;
    //! time constant for memory function in disres
    real                        dr_tau = 0;
    //! force constant for orientational restraints
    real                        orires_fc = 0;
    //! time constant for memory function in orires
    real                        orires_tau = 0;
    //! frequency of writing tr(SD) to energy output
    int                         nstorireout = 0;
    //! The stepsize for updating
    real                        em_stepsize = 0;
    //! The tolerance
    real                        em_tol = 0;
    //! Number of iterations for convergence of steepest descent in relax_shells
    int                         niter = 0;
    //! Stepsize for directional minimization in relax_shells
    real                        fc_stepsize = 0;
    //! number of steps after which a steepest descents step is done while doing cg
    int                         nstcgsteep = 0;
    //! Number of corrections to the hessian to keep
    int                         nbfgscorr = 0;
    //! Type of constraint algorithm
    int                         eConstrAlg = 0;
    //! Order of the LINCS Projection Algorithm
    int                         nProjOrder = 0;
    //! If bond rotates more than %g degrees, warn
    real                        LincsWarnAngle = 0;
    //! Number of iterations in the final Lincs step
    int                         nLincsIter = 0;
    //! Use successive overrelaxation for shake
    gmx_bool                    bShakeSOR = 0;
    //! Friction coefficient for BD (amu/ps)
    real                        bd_fric = 0;
    //! Random seed for SD and BD
    gmx_int64_t                 ld_seed = 0;
    //! The number of walls
    int                         nwall = 0;
    //! The type of walls
    int                         wall_type = 0;
    //! The potentail is linear for r<=wall_r_linpot
    real                        wall_r_linpot = 0;
    //! The atom type for walls
    int                         wall_atomtype[2] = {0, 0};
    //! Number density for walls
    real                        wall_density[2] = {0, 0};
    //! Scaling factor for the box for Ewald
    real                        wall_ewald_zfac = 0;

    /* COM pulling data */
    //! Do we do COM pulling?
    gmx_bool       bPull = 0;
    //! The data for center of mass pulling
    pull_params_t *pull = nullptr;
    // TODO: Remove this by converting pull into a ForceProvider
    //! The COM pull force calculation data structure
    pull_t *pull_work = nullptr;

    /* AWH bias data */
    //! Use awh biasing for PMF calculations?
    gmx_bool        bDoAwh = 0;
    //! AWH biasing parameters
    gmx::AwhParams *awhParams = nullptr;

    /* Enforced rotation data */
    //! Calculate enforced rotation potential(s)?
    gmx_bool               bRot = 0;
    //! The data for enforced rotation potentials
    std::unique_ptr<t_rot> rot;

    //! Do ion/water position exchanges (CompEL)?
    int                           eSwapCoords = 0;
    //! Swap data structure.
    std::unique_ptr<t_swapcoords> swap;

    //! Allow interactive MD sessions for this .tpr?
    gmx_bool               bIMD = 0;
    //! Interactive molecular dynamics
    std::unique_ptr<t_IMD> imd;

    //! Acceleration for viscosity calculation
    real   cos_accel = 0;
    //! Triclinic deformation velocities (nm/ps)
    tensor deform = {{0}};
    /*! \brief User determined parameters */
    /**@{*/
    int  userint1  = 0;
    int  userint2  = 0;
    int  userint3  = 0;
    int  userint4  = 0;
    real userreal1 = 0;
    real userreal2 = 0;
    real userreal3 = 0;
    real userreal4 = 0;
    /**@}*/
    //! Group options
    t_grpopts opts;
    //! QM/MM calculation
    gmx_bool  bQMMM = 0;
    //! constraints on QM bonds
    int       QMconstraints = 0;
    //! Scheme: ONIOM or normal
    int       QMMMscheme = 0;
    //! factor for scaling the MM charges in QM calc.
    real      scalefactor = 0;

    /* Fields for removed features go here (better caching) */
    //! Whether AdResS is enabled - always false if a valid .tpr was read
    gmx_bool                 bAdress = 0;
    //! Whether twin-range scheme is active - always false if a valid .tpr was read
    gmx_bool                 useTwinRange = 0;

    //! KVT object that contains input parameters converted to the new style.
    gmx::KeyValueTreeObject *params = nullptr;
};

int ir_optimal_nstcalcenergy(const t_inputrec *ir);

int tcouple_min_integration_steps(int etc);

int ir_optimal_nsttcouple(const t_inputrec *ir);

int pcouple_min_integration_steps(int epc);

int ir_optimal_nstpcouple(const t_inputrec *ir);

/* Returns if the Coulomb force or potential is switched to zero */
gmx_bool ir_coulomb_switched(const t_inputrec *ir);

/* Returns if the Coulomb interactions are zero beyond the rcoulomb.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec *ir);

/* As ir_coulomb_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rcoulomb.
 */
gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec *ir);

/* Returns if the Van der Waals force or potential is switched to zero */
gmx_bool ir_vdw_switched(const t_inputrec *ir);

/* Returns if the Van der Waals interactions are zero beyond the rvdw.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec *ir);

/* As ir_vdw_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rvdw.
 */
gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec *ir);

/*! \brief Free memory from input record.
 *
 * All arrays and pointers will be freed.
 *
 * \param[in] ir The data structure
 */
void done_inputrec(t_inputrec *ir);

void pr_inputrec(FILE *fp, int indent, const char *title, const t_inputrec *ir,
                 gmx_bool bMDPformat);

void cmp_inputrec(FILE *fp, const t_inputrec *ir1, const t_inputrec *ir2, real ftol, real abstol);

void comp_pull_AB(FILE *fp, pull_params_t *pull, real ftol, real abstol);


gmx_bool inputrecDeform(const t_inputrec *ir);

gmx_bool inputrecDynamicBox(const t_inputrec *ir);

gmx_bool inputrecPreserveShape(const t_inputrec *ir);

gmx_bool inputrecNeedMutot(const t_inputrec *ir);

gmx_bool inputrecTwinRange(const t_inputrec *ir);

gmx_bool inputrecExclForces(const t_inputrec *ir);

gmx_bool inputrecNptTrotter(const t_inputrec *ir);

gmx_bool inputrecNvtTrotter(const t_inputrec *ir);

gmx_bool inputrecNphTrotter(const t_inputrec *ir);

/*! \brief Return true if the simulation is 2D periodic with two walls. */
bool     inputrecPbcXY2Walls(const t_inputrec *ir);

/*! \brief Returns true for MD integator with T and/or P-coupling that supports
 * calculating the conserved energy quantity.
 */
bool integratorHasConservedEnergyQuantity(const t_inputrec *ir);

/*! \brief Returns true when temperature is coupled or constant. */
bool integratorHasReferenceTemperature(const t_inputrec *ir);

/*! \brief Return the number of bounded dimensions
 *
 * \param[in] ir The input record with MD parameters
 * \return the number of dimensions in which
 * the coordinates of the particles are bounded, starting at X.
 */
int inputrec2nboundeddim(const t_inputrec *ir);

/*! \brief Returns the number of degrees of freedom in center of mass motion
 *
 * \param[in] ir the inputrec structure
 * \return the number of degrees of freedom of the center of mass
 */
int ndof_com(const t_inputrec *ir);

#endif /* GMX_MDTYPES_INPUTREC_H */
