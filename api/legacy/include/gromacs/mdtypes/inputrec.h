/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_MDTYPES_INPUTREC_H
#define GMX_MDTYPES_INPUTREC_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#define EGP_EXCL (1 << 0)
#define EGP_TABLE (1 << 1)

struct gmx_enfrot;
struct gmx_enfrotgrp;
struct pull_params_t;

namespace gmx
{
class Awh;
class AwhParams;
template<typename T1, typename T2, T1 U>
struct EnumerationArray;
class KeyValueTreeObject;
struct MtsLevel;
} // namespace gmx

struct t_grpopts
{
    //! Number of T-Coupl groups
    int ngtc = 0;
    //! Number of of Nose-Hoover chains per group
    int nhchainlength = 0;
    //! Number of Accelerate groups
    int ngacc;
    //! Number of Freeze groups
    int ngfrz = 0;
    //! Number of Energy groups
    int ngener = 0;
    //! Number of degrees of freedom in a temperature-coupling group
    real* nrdf = nullptr;
    //! Coupling temperature	per group
    real* ref_t = nullptr;
    //! No/simple/periodic simulated annealing for each group
    SimulatedAnnealing* annealing = nullptr;
    //! Number of annealing time points per group
    int* anneal_npoints = nullptr;
    //! For each group: Time points
    real** anneal_time = nullptr;
    //! For each group: Temperature at these times. Final temp after all intervals is ref_t
    real** anneal_temp = nullptr;
    //! Tau coupling time
    real* tau_t = nullptr;
    //! Acceleration per group
    rvec* acceleration = nullptr;
    //! Whether the group will be frozen in each direction
    ivec* nFreeze = nullptr;
    //! Exclusions/tables of energy group pairs
    int* egp_flags = nullptr;

    /* QMMM stuff */
    //! Number of QM groups
    int ngQM = 0;
};

struct t_simtemp
{
    //! Simulated temperature scaling; linear or exponential
    SimulatedTempering eSimTempScale = SimulatedTempering::Default;
    //! The low temperature for simulated tempering
    real simtemp_low = 0;
    //! The high temperature for simulated tempering
    real simtemp_high = 0;
    //! The range of temperatures used for simulated tempering
    std::vector<real> temperatures;
};

struct t_lambda
{
    //! Return the initial lambda
    double initialLambda(const FreeEnergyPerturbationCouplingType couplingType) const
    {
        if (init_lambda_without_states >= 0)
        {
            return init_lambda_without_states;
        }
        else
        {
            return all_lambda[couplingType][init_fep_state];
        }
    }

    //! The frequency for calculating dhdl
    int nstdhdl = 0;
    //! Fractional value of lambda (usually will use init_fep_state, this will only be for slow growth, and for legacy free energy code. Only has a valid value if positive)
    double init_lambda_without_states = -1;
    //! The initial number of the state
    int init_fep_state = -1;
    //! Change of lambda per time step (fraction of (0.1)
    double delta_lambda = 0;
    //! Print no, total or potential energies in dhdl
    FreeEnergyPrintEnergy edHdLPrintEnergy = FreeEnergyPrintEnergy::Default;
    //! The number of foreign lambda points
    int n_lambda = 0;
    //! The array of all lambda values
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>> all_lambda;
    //! The number of neighboring lambda states to calculate the energy for in up and down directions (-1 for all)
    int lambda_neighbors = 0;
    //! The first lambda to calculate energies for
    int lambda_start_n = 0;
    //! The last lambda +1 to calculate energies for
    int lambda_stop_n = 0;
    //! Free energy soft-core parameter
    real sc_alpha = 0;
    //! Lambda power for soft-core interactions
    int sc_power = 0;
    //! R power for soft-core interactions
    real sc_r_power = 0;
    //! Free energy soft-core sigma when c6 or c12=0
    real sc_sigma = 0;
    //! Free energy soft-core sigma for ?????
    real sc_sigma_min = 0;
    //! Use softcore for the coulomb portion as well (default FALSE)
    bool bScCoul = false;
    //! The specific soft-core function to use
    SoftcoreType softcoreFunction = SoftcoreType::Beutler;
    //! scale for the linearization point for the vdw interaction with gapsys soft-core
    real scGapsysScaleLinpointLJ = 0.85;
    //! scale for the linearization point for the coulomb interaction with gapsys soft-core
    real scGapsysScaleLinpointQ = 0.3;
    //! lower bound for c12/c6 in gapsys soft-core
    real scGapsysSigmaLJ = 0.3;
    //! Whether to print the dvdl term associated with this term; if it is not specified as separate, it is lumped with the FEP term
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, bool> separate_dvdl;
    //! Whether to write a separate dhdl.xvg file note: NOT a gmx_bool, but an enum
    SeparateDhdlFile separate_dhdl_file = SeparateDhdlFile::Default;
    //! Whether to calculate+write dhdl derivatives note: NOT a gmx_bool, but an enum
    DhDlDerivativeCalculation dhdl_derivatives = DhDlDerivativeCalculation::Default;
    //! The maximum table size for the dH histogram
    int dh_hist_size = 0;
    //! The spacing for the dH histogram
    double dh_hist_spacing = 0;
};

struct t_expanded
{
    //! The frequency of expanded ensemble state changes
    int nstexpanded = 0;
    //! Which type of move updating do we use for lambda monte carlo (or no for none)
    LambdaWeightCalculation elamstats = LambdaWeightCalculation::Default;
    //! What move set will be we using for state space moves
    LambdaMoveCalculation elmcmove = LambdaMoveCalculation::Default;
    //! The method we use to decide of we have equilibrated the weights
    LambdaWeightWillReachEquilibrium elmceq = LambdaWeightWillReachEquilibrium::Default;
    //! The minimum number of samples at each lambda for deciding whether we have reached a minimum
    int equil_n_at_lam = 0;
    //! Wang-Landau delta at which we stop equilibrating weights
    real equil_wl_delta = 0;
    //! Use the ratio of weights (ratio of minimum to maximum) to decide when to stop equilibrating
    real equil_ratio = 0;
    //! After equil_steps steps we stop equilibrating the weights
    int equil_steps = 0;
    //! After equil_samples total samples (steps/nstfep), we stop equilibrating the weights
    int equil_samples = 0;
    //! Random number seed for lambda mc switches
    int lmc_seed = 0;
    //! Whether to use minimum variance weighting
    bool minvar = false;
    //! The number of samples needed before kicking into minvar routine
    int minvarmin = 0;
    //! The offset for the variance in MinVar
    real minvar_const = 0;
    //! Range of cvalues used for BAR
    int c_range = 0;
    //! Whether to print symmetrized matrices
    bool bSymmetrizedTMatrix = false;
    //! How frequently to print the transition matrices
    int nstTij = 0;
    //! Number of repetitions in the MC lambda jumps MRS -- VERIFY THIS
    int lmc_repeats = 0;
    //! Minimum number of samples for each state before free sampling MRS -- VERIFY THIS!
    int lmc_forced_nstart = 0;
    //! Distance in lambda space for the gibbs interval
    int gibbsdeltalam = 0;
    //! Scaling factor for Wang-Landau
    real wl_scale = 0;
    //! Ratio between largest and smallest number for freezing the weights
    real wl_ratio = 0;
    //! Starting delta for Wang-Landau
    real init_wl_delta = 0;
    //! Use one over t convergence for Wang-Landau when the delta get sufficiently small
    bool bWLoneovert = false;
    //! Did we initialize the weights? TODO: REMOVE FOR 5.0, no longer needed with new logic
    bool bInit_weights = false;
    //! To override the main temperature, or define it if it's not defined
    real mc_temp = 0;
    //! User-specified initial weights to start with
    std::vector<real> init_lambda_weights;
};

struct t_rotgrp
{
    //! Rotation type for this group
    EnforcedRotationGroupType eType = EnforcedRotationGroupType::Default;
    //! Use mass-weighed positions?
    bool bMassW = false;
    //! Number of atoms in the group
    int nat = 0;
    //! The global atoms numbers
    int* ind = nullptr;
    //! The reference positions (which have not been centered)
    std::vector<gmx::RVec> x_ref_original;
    //! The normalized rotation vector
    rvec inputVec = { 0, 0, 0 };
    //! Rate of rotation (degree/ps)
    real rate = 0;
    //! Force constant (kJ/(mol nm^2)
    real k = 0;
    //! Pivot point of rotation axis (nm)
    rvec pivot = { 0, 0, 0 };
    //! Type of fit to determine actual group angle
    RotationGroupFitting eFittype = RotationGroupFitting::Default;
    //! Number of angles around the reference angle for which the rotation potential is also evaluated (for fit type 'potential' only)
    int PotAngle_nstep = 0;
    //! Distance between two angles in degrees (for fit type 'potential' only)
    real PotAngle_step = 0;
    //! Slab distance (nm)
    real slab_dist = 0;
    //! Minimum value the gaussian must have so that the force is actually evaluated
    real min_gaussian = 0;
    //! Additive constant for radial motion2 and flexible2 potentials (nm^2)
    real eps = 0;
};

struct t_rot
{
    //! Output frequency for main rotation outfile
    int nstrout;
    //! Output frequency for per-slab data
    int nstsout;
    //! Groups to rotate
    std::vector<t_rotgrp> grp;
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
    gmx::EnumerationArray<Compartment, int> nmolReq;
};

struct t_swapcoords
{
    //! Period between when a swap is attempted
    int nstswap;
    //! Use mass-weighted positions in split group
    bool massw_split[2];
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
    gmx::EnumerationArray<Compartment, real> bulkOffset;
    //! Number of groups to be controlled
    int ngrp;
    //! All swap groups, including split and solvent
    t_swapGroup* grp;
};

struct PressureCouplingOptions
{
    //! Pressure coupling algorithm
    PressureCoupling epc = PressureCoupling::Default;
    //! Pressure coupling isotropy
    PressureCouplingType epct = PressureCouplingType::Default;
    //! Interval in steps for pressure coupling
    int nstpcouple = 0;
    //! Pressure coupling time (ps)
    real tau_p = 0;
    //! Reference pressure (kJ/(mol nm^3))
    tensor ref_p = { { 0 } };
    //! Compressibility ((mol nm^3)/kJ)
    tensor compress = { { 0 } };
    //! How to scale absolute reference coordinates
    RefCoordScaling refcoord_scaling = RefCoordScaling::Default;
};

struct t_inputrec // NOLINT (clang-analyzer-optin.performance.Padding)
{
    t_inputrec();
    explicit t_inputrec(const t_inputrec&) = delete;
    t_inputrec& operator=(const t_inputrec&) = delete;
    ~t_inputrec();

    //! The tpx version number this inputrec was read from, -1 when not read from tpr
    int tpxFileVersion = -1;

    //! Integration method
    IntegrationAlgorithm eI = IntegrationAlgorithm::Default;
    //! Number of steps to be taken
    int64_t nsteps = 0;
    //! Used in checkpointing to separate chunks
    int simulation_part = 0;
    //! Start at a stepcount >0 (used w. convert-tpr)
    int64_t init_step = 0;
    //! Frequency of energy calc. and T/P coupl. upd.
    int nstcalcenergy = 0;
    //! Group or verlet cutoffs
    CutoffScheme cutoff_scheme = CutoffScheme::Default;
    //! Number of steps before pairlist is generated
    int nstlist = 0;
    //! Number of steps after which center of mass motion is removed
    int nstcomm = 0;
    //! Center of mass motion removal algorithm
    ComRemovalAlgorithm comm_mode = ComRemovalAlgorithm::Default;
    //! Number of steps after which print to logfile
    int nstlog = 0;
    //! Number of steps after which X is output
    int nstxout = 0;
    //! Number of steps after which V is output
    int nstvout = 0;
    //! Number of steps after which F is output
    int nstfout = 0;
    //! Number of steps after which energies printed
    int nstenergy = 0;
    //! Number of steps after which compressed trj (.xtc,.tng) is output
    int nstxout_compressed = 0;
    //! Initial time (ps)
    double init_t = 0;
    //! Time step (ps)
    double delta_t = 0;
    //! Whether we use multiple time stepping
    bool useMts = false;
    //! The multiple time stepping levels
    std::vector<gmx::MtsLevel> mtsLevels;

    //! The factor for repartitioning atom masses
    real massRepartitionFactor = 1;

    //! Precision of x in compressed trajectory file
    real x_compression_precision = 0;

    //! Requested fourier_spacing, when nk? not set
    real fourier_spacing = 0;
    //! Number of k vectors in x dimension for fourier methods for long range electrost.
    int nkx = 0;
    //! Number of k vectors in y dimension for fourier methods for long range electrost.
    int nky = 0;
    //! Number of k vectors in z dimension for fourier methods for long range electrost.
    int nkz = 0;
    //! Interpolation order for PME
    int pme_order = 0;
    //! Real space tolerance for Ewald, determines the real/reciprocal space relative weight
    real ewald_rtol = 0;
    //! Real space tolerance for LJ-Ewald
    real ewald_rtol_lj = 0;
    //! Normal/3D ewald, or pseudo-2D LR corrections
    EwaldGeometry ewald_geometry = EwaldGeometry::Default;
    //! Epsilon for PME dipole correction
    real epsilon_surface = 0;
    //! Type of combination rule in LJ-PME
    LongRangeVdW ljpme_combination_rule = LongRangeVdW::Default;
    //! Type of periodic boundary conditions
    PbcType pbcType = PbcType::Default;
    //! Periodic molecules
    bool bPeriodicMols = false;
    //! Continuation run: starting state is correct (ie. constrained)
    bool bContinuation = false;
    //! Whether and what kind of ensemble temperature we have for the system
    EnsembleTemperatureSetting ensembleTemperatureSetting;
    //! The ensemble temperature of the system, see ensembleTemperatureSetting for validity
    real ensembleTemperature;
    //! Temperature coupling
    TemperatureCoupling etc = TemperatureCoupling::Default;
    //! Interval in steps for temperature coupling
    int nsttcouple = 0;
    //! Whether to print nose-hoover chains
    bool bPrintNHChains = false;
    //! Pressure coupling
    PressureCouplingOptions pressureCouplingOptions;
    //! The COM of the posres atoms
    rvec posres_com = { 0, 0, 0 };
    //! The B-state COM of the posres atoms
    rvec posres_comB = { 0, 0, 0 };
    //! Random seed for Andersen thermostat (obsolete)
    int andersen_seed = 0;
    //! Per atom pair energy drift tolerance (kJ/mol/ps/atom) for the pairlist buffer
    real verletbuf_tol = 0;
    //! The tolerance for the average LJ pressure deviation for the pairlist buffer
    real verletBufferPressureTolerance = 0;
    //! Short range pairlist cut-off (nm)
    real rlist = 0;
    //! Radius for test particle insertion
    real rtpi = 0;
    //! Type of electrostatics treatment
    CoulombInteractionType coulombtype = CoulombInteractionType::Default;
    //! Modify the Coulomb interaction
    InteractionModifiers coulomb_modifier = InteractionModifiers::Default;
    //! Coulomb switch range start (nm)
    real rcoulomb_switch = 0;
    //! Coulomb cutoff (nm)
    real rcoulomb = 0;
    //! Relative dielectric constant
    real epsilon_r = 0;
    //! Relative dielectric constant of the RF
    real epsilon_rf = 0;
    //! Always false (no longer supported)
    bool implicit_solvent = false;
    //! Type of Van der Waals treatment
    VanDerWaalsType vdwtype = VanDerWaalsType::Default;
    //! Modify the Van der Waals interaction
    InteractionModifiers vdw_modifier = InteractionModifiers::Default;
    //! Van der Waals switch range start (nm)
    real rvdw_switch = 0;
    //! Van der Waals cutoff (nm)
    real rvdw = 0;
    //! Perform Long range dispersion corrections
    DispersionCorrectionType eDispCorr = DispersionCorrectionType::Default;
    //! Extension of the table beyond the cut-off, as well as the table length for 1-4 interac.
    real tabext = 0;
    //! Tolerance for shake
    real shake_tol = 0;
    //! Free energy calculations
    FreeEnergyPerturbationType efep = FreeEnergyPerturbationType::Default;
    //! Data for the FEP state
    std::unique_ptr<t_lambda> fepvals;
    //! Whether to do simulated tempering
    bool bSimTemp = false;
    //! Variables for simulated tempering
    std::unique_ptr<t_simtemp> simtempvals;
    //! Whether expanded ensembles are used
    bool bExpanded = false;
    //! Expanded ensemble parameters
    std::unique_ptr<t_expanded> expandedvals;
    //! Type of distance restraining
    DistanceRestraintRefinement eDisre = DistanceRestraintRefinement::Default;
    //! Force constant for time averaged distance restraints
    real dr_fc = 0;
    //! Type of weighting of pairs in one restraints
    DistanceRestraintWeighting eDisreWeighting = DistanceRestraintWeighting::Default;
    //! Use combination of time averaged and instantaneous violations
    bool bDisreMixed = false;
    //! Frequency of writing pair distances to enx
    int nstdisreout = 0;
    //! Time constant for memory function in disres
    real dr_tau = 0;
    //! Force constant for orientational restraints
    real orires_fc = 0;
    //! Time constant for memory function in orires
    real orires_tau = 0;
    //! Frequency of writing tr(SD) to energy output
    int nstorireout = 0;
    //! The stepsize for updating
    real em_stepsize = 0;
    //! The tolerance
    real em_tol = 0;
    //! Number of iterations for convergence of steepest descent in relax_shells
    int niter = 0;
    //! Stepsize for directional minimization in relax_shells
    real fc_stepsize = 0;
    //! Number of steps after which a steepest descents step is done while doing cg
    int nstcgsteep = 0;
    //! Number of corrections to the Hessian to keep
    int nbfgscorr = 0;
    //! Type of constraint algorithm
    ConstraintAlgorithm eConstrAlg = ConstraintAlgorithm::Default;
    //! Order of the LINCS Projection Algorithm
    int nProjOrder = 0;
    //! Warn if any bond rotates more than this many degrees
    real LincsWarnAngle = 0;
    //! Number of iterations in the final LINCS step
    int nLincsIter = 0;
    //! Use successive overrelaxation for shake
    bool bShakeSOR = false;
    //! Friction coefficient for BD (amu/ps)
    real bd_fric = 0;
    //! Random seed for SD and BD
    int64_t ld_seed = 0;
    //! The number of walls
    int nwall = 0;
    //! The type of walls
    WallType wall_type = WallType::Default;
    //! The potentail is linear for r<=wall_r_linpot
    real wall_r_linpot = 0;
    //! The atom type for walls
    int wall_atomtype[2] = { 0, 0 };
    //! Number density for walls
    real wall_density[2] = { 0, 0 };
    //! Scaling factor for the box for Ewald
    real wall_ewald_zfac = 0;

    /* COM pulling data */
    //! Do we do COM pulling?
    bool bPull = false;
    //! The data for center of mass pulling
    std::unique_ptr<pull_params_t> pull;

    /* AWH bias data */
    //! Whether to use AWH biasing for PMF calculations
    bool bDoAwh = false;
    //! AWH biasing parameters
    std::unique_ptr<gmx::AwhParams> awhParams;

    /* Enforced rotation data */
    //! Whether to calculate enforced rotation potential(s)
    bool bRot = false;
    //! The data for enforced rotation potentials
    std::unique_ptr<t_rot> rot;

    //! Whether to do ion/water position exchanges (CompEL)
    SwapType eSwapCoords = SwapType::Default;
    //! Swap data structure.
    t_swapcoords* swap = nullptr;

    //! Whether the tpr makes an interactive MD session possible.
    bool bIMD = false;
    //! Interactive molecular dynamics
    t_IMD* imd = nullptr;

    //! Acceleration for viscosity calculation
    real cos_accel = 0;
    //! Triclinic deformation velocities (nm/ps)
    tensor deform = { { 0 } };
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
    bool bQMMM = false;

    /* Fields for removed features go here (better caching) */
    //! Whether AdResS is enabled - always false if a valid .tpr was read
    bool bAdress = false;
    //! Whether twin-range scheme is active - always false if a valid .tpr was read
    bool useTwinRange = false;
    //! Whether we have constant acceleration
    bool useConstantAcceleration = false;

    //! KVT object that contains input parameters converted to the new style.
    gmx::KeyValueTreeObject* params = nullptr;

    //! KVT for storing simulation parameters that are not part of the mdp file.
    std::unique_ptr<gmx::KeyValueTreeObject> internalParameters;
};

int tcouple_min_integration_steps(TemperatureCoupling etc);

int ir_optimal_nsttcouple(const t_inputrec* ir);

int pcouple_min_integration_steps(PressureCoupling epc);

int ir_optimal_nstpcouple(const t_inputrec* ir);

/* Returns if the Coulomb force or potential is switched to zero */
bool ir_coulomb_switched(const t_inputrec* ir);

/* Returns if the Coulomb interactions are zero beyond the rcoulomb.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
bool ir_coulomb_is_zero_at_cutoff(const t_inputrec* ir);

/* As ir_coulomb_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rcoulomb.
 */
bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec* ir);

/* Returns if the Van der Waals force or potential is switched to zero */
bool ir_vdw_switched(const t_inputrec* ir);

/* Returns if the Van der Waals interactions are zero beyond the rvdw.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
bool ir_vdw_is_zero_at_cutoff(const t_inputrec* ir);

/* As ir_vdw_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rvdw.
 */
bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec* ir);

/*! \brief Free memory from input record.
 *
 * All arrays and pointers will be freed.
 *
 * \param[in] ir The data structure
 */
void done_inputrec(t_inputrec* ir);

void pr_inputrec(FILE* fp, int indent, const char* title, const t_inputrec* ir, bool bMDPformat);

void cmp_inputrec(FILE* fp, const t_inputrec* ir1, const t_inputrec* ir2, real ftol, real abstol);

void comp_pull_AB(FILE* fp, const pull_params_t& pull, real ftol, real abstol);


bool inputrecDeform(const t_inputrec* ir);

bool inputrecDynamicBox(const t_inputrec* ir);

bool shouldPreserveBoxShape(const PressureCouplingOptions& pressureCoupling, const tensor deform);

bool inputrecNeedMutot(const t_inputrec* ir);

bool inputrecExclForces(const t_inputrec* ir);

bool inputrecNptTrotter(const t_inputrec* ir);

bool inputrecNvtTrotter(const t_inputrec* ir);

bool inputrecNphTrotter(const t_inputrec* ir);

/*! \brief Return true if the simulation is 2D periodic with two walls. */
bool inputrecPbcXY2Walls(const t_inputrec* ir);

//! \brief Return true if the simulation has frozen atoms (non-trivial freeze groups).
bool inputrecFrozenAtoms(const t_inputrec* ir);

/*! \brief Returns true when a constant ensemble temperature is available for the system. */
bool haveConstantEnsembleTemperature(const t_inputrec& ir);

/*! \brief Returns the constant ensemble temperature for the system. */
real constantEnsembleTemperature(const t_inputrec& ir);

/*! \brief Returns true when a constant or variable ensemble temperature is available for the system.
 *
 * \note The current ensemble temperature can be obtained from \p gmx_ekindata_t.
 */
bool haveEnsembleTemperature(const t_inputrec& ir);

/*! \brief Returns true for MD integator with T and/or P-coupling that supports
 * calculating a conserved energy quantity.
 *
 * Note that dynamical integrators without T and P coupling (ie NVE)
 * return false, i.e. the return value refers to whether there
 * is a conserved quantity different than the total energy.
 */
bool integratorHasConservedEnergyQuantity(const t_inputrec* ir);

/*! \brief Returns true when the integrator, and possibly T-coupling, has a reference temperature. */
bool integratorHasReferenceTemperature(const t_inputrec& ir);

/*! \brief Returns whether we are doing simulated annealing */
bool doSimulatedAnnealing(const t_inputrec& ir);

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

/*! \brief Check if calculation of the specific FEP type was requested.
 *
 * \param[in] ir       Input record.
 * \param[in] fepType  Free-energy perturbation type to check for.
 *
 * \returns If the \p fepType is perturbed in this run.
 */
bool haveFreeEnergyType(const t_inputrec& ir, int fepType);

/*! \brief Checks whether all lambda components change at the same or zero rate
 *
 * Returns true when over each lambda interval the change in lambda is the same
 * for all components. The change can differ by a factor -1 and it can be zero
 * for some or all components. When this function returns true, the composed,
 * total dH/lambda value will be correct and meaningful.
 *
 * \param[in] allLambdas  Sets of lambda values for all components, as used in \p t_lambda
 *
 * \returns Whether all lambda components change at the same or zero rate
 */
bool fepLambdasChangeAtSameRate(
        const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>>& allLambdas);

//! Return whether the box is continuously deformed
bool ir_haveBoxDeformation(const t_inputrec& ir);

#endif /* GMX_MDTYPES_INPUTREC_H */
