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
/*! \file
 * \brief
 * Declares enumerated types used throughout the code.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_MD_ENUMS_H
#define GMX_MDTYPES_MD_ENUMS_H

/*! \brief Return a string from a list of strings
 *
 * If index if within 0 .. max_index-1 returns the corresponding string
 * or "no name defined" otherwise, in other words this is a range-check that does
 * not crash.
 * \param[in] index     The index in the array
 * \param[in] max_index The length of the array
 * \param[in] names     The array
 * \return the correct string or "no name defined"
 */
const char* enum_name(int index, int max_index, const char* const names[]);

/*! \brief Enum for setting answer to yes or no
 */
enum class Boolean : int
{
    No,
    Yes,
    Count,
    Default = No
};

//! Return name of boolean selection.
const char* enumValueToString(Boolean enumValue);
//! Return name of boolean selection for actual bool.
const char* booleanValueToString(bool value);

//! \brief The two compartments for CompEL setups.
enum class Compartment : int
{
    A,
    B,
    Count
};

/*! \brief The channels that define with their COM the compartment boundaries in CompEL setups.
 *
 * In principle one could also use modified setups with more than two channels.
 */
enum class Channel : int
{
    Zero,
    One,
    Count
};

/*! \brief Temperature coupling type
 *
 * yes is an alias for berendsen
 *
 * Note: Keep `Count` as the second-to-last entry, and `Default` as the last entry -
 *       this is needed to keep EnumerationWrapper, EnumerationArray and (de)serialization
 *       working.
 */
enum class TemperatureCoupling : int
{
    No,
    Berendsen,
    NoseHoover,
    Yes,
    Andersen,
    AndersenMassive,
    VRescale,
    Count,
    Default = No
};
//! Return names of temperature coupling schemes
const char* enumValueToString(TemperatureCoupling enumValue);
//! Return whether this is andersen coupling
#define ETC_ANDERSEN(e) \
    (((e) == TemperatureCoupling::AndersenMassive) || ((e) == TemperatureCoupling::Andersen))

/*! \brief Pressure coupling types
 *
 * isotropic is an alias for berendsen
 *
 * Note: Keep `Count` as the second-to-last entry, and `Default` as the last entry -
 *       this is needed to keep EnumerationWrapper, EnumerationArray and (de)serialization
 *       working.
 */
enum class PressureCoupling : int
{
    No,
    Berendsen,
    ParrinelloRahman,
    Isotropic,
    Mttk,
    CRescale,
    Count,
    Default = No
};
//! Return names of pressure coupling schemes
const char* enumValueToString(PressureCoupling enumValue);

//! Flat-bottom posres geometries
enum
{
    efbposresZERO,
    efbposresSPHERE,
    efbposresCYLINDER,
    efbposresX,
    efbposresY,
    efbposresZ,
    efbposresCYLINDERX,
    efbposresCYLINDERY,
    efbposresCYLINDERZ,
    efbposresNR
};

//! Relative coordinate scaling type for position restraints.
enum class RefCoordScaling : int
{
    No,
    All,
    Com,
    Count,
    Default = No
};

//! String corresponding to relative coordinate scaling.
const char* enumValueToString(RefCoordScaling enumValue);

//! Trotter decomposition extended variable parts.
enum
{
    etrtNONE,
    etrtNHC,
    etrtBAROV,
    etrtBARONHC,
    etrtNHC2,
    etrtBAROV2,
    etrtBARONHC2,
    etrtVELOCITY1,
    etrtVELOCITY2,
    etrtPOSITION,
    etrtSKIPALL,
    etrtNR
};

//! Sequenced parts of the trotter decomposition.
enum class TrotterSequence : int
{
    Zero,
    One,
    Two,
    Three,
    Four,
    Count
};

//! Pressure coupling type
enum class PressureCouplingType : int
{
    Isotropic,
    SemiIsotropic,
    Anisotropic,
    SurfaceTension,
    Count,
    Default = Isotropic
};
//! String corresponding to pressure coupling type
const char* enumValueToString(PressureCouplingType enumValue);

//! \\brief Cutoff scheme
enum class CutoffScheme : int
{
    Verlet,
    Group,
    Count,
    Default = Verlet
};
//! String corresponding to cutoff scheme
const char* enumValueToString(CutoffScheme enumValue);

/*! \brief Coulomb / VdW interaction modifiers.
 *
 * grompp replaces eintmodPOTSHIFT_VERLET_UNSUPPORTED by eintmodPOTSHIFT.
 * Exactcutoff is only used by Reaction-field-zero, and is not user-selectable.
 */
enum class InteractionModifiers : int
{
    PotShiftVerletUnsupported,
    PotShift,
    None,
    PotSwitch,
    ExactCutoff,
    ForceSwitch,
    Count,
    Default = PotShiftVerletUnsupported
};
//! String corresponding to interaction modifiers
const char* enumValueToString(InteractionModifiers enumValue);

/*! \brief Cut-off treatment for Coulomb */
enum class CoulombInteractionType : int
{
    Cut,
    RF,
    GRFNotused,
    Pme,
    Ewald,
    P3mAD,
    Poisson,
    Switch,
    Shift,
    User,
    GBNotused,
    RFNecUnsupported,
    EncadShiftNotused,
    PmeUser,
    PmeSwitch,
    PmeUserSwitch,
    RFZero,
    Count,
    Default = Cut
};
//! String corresponding to Coulomb treatment
const char* enumValueToString(CoulombInteractionType enumValue);

//! Ewald geometry.
enum class EwaldGeometry : int
{
    ThreeD,
    ThreeDC,
    Count,
    Default = ThreeD
};
//! String corresponding to Ewald geometry
const char* enumValueToString(EwaldGeometry enumValue);

//! Returns whether we use reaction field
static inline bool usingRF(const CoulombInteractionType& cit)
{
    return (cit == CoulombInteractionType::RF || cit == CoulombInteractionType::GRFNotused
            || cit == CoulombInteractionType::RFNecUnsupported || cit == CoulombInteractionType::RFZero);
};

//! Returns whether we use PME
static inline bool usingPme(const CoulombInteractionType& cit)
{
    return (cit == CoulombInteractionType::Pme || cit == CoulombInteractionType::PmeSwitch
            || cit == CoulombInteractionType::PmeUser
            || cit == CoulombInteractionType::PmeUserSwitch || cit == CoulombInteractionType::P3mAD);
}

//! Returns whether we use PME or full Ewald
static inline bool usingPmeOrEwald(const CoulombInteractionType& cit)
{
    return (usingPme(cit) || cit == CoulombInteractionType::Ewald);
};

//! Returns whether we use full electrostatics of any sort
static inline bool usingFullElectrostatics(const CoulombInteractionType& cit)
{
    return (usingPmeOrEwald(cit) || cit == CoulombInteractionType::Poisson);
}

//! Returns whether we use user defined electrostatics
static inline bool usingUserTableElectrostatics(const CoulombInteractionType& cit)
{
    return (cit == CoulombInteractionType::User || cit == CoulombInteractionType::PmeUser
            || cit == CoulombInteractionType::PmeUserSwitch);
}

//! Van der Waals interaction treatment
enum class VanDerWaalsType : int
{
    Cut,
    Switch,
    Shift,
    User,
    EncadShiftUnused,
    Pme,
    Count,
    Default = Cut
};
//! String corresponding to Van der Waals treatment
const char* enumValueToString(VanDerWaalsType enumValue);

//! Type of long-range VdW treatment of combination rules
enum class LongRangeVdW : int
{
    Geom,
    LB,
    Count,
    Default = Geom
};
//! String for LJPME combination rule treatment
const char* enumValueToString(LongRangeVdW enumValue);

//! Returns whether we use LJPME
static inline bool usingLJPme(const VanDerWaalsType& vanDerWaalsType)
{
    return vanDerWaalsType == VanDerWaalsType::Pme;
};

/*! \brief Integrator algorithm
 *
 * eiSD2 has been removed, but we keep a renamed enum entry,
 * so we can refuse to do MD with such .tpr files.
 * eiVV is normal velocity verlet
 * eiVVAK uses 1/2*(KE(t-dt/2)+KE(t+dt/2)) as the kinetic energy,
 * and the half step kinetic energy for temperature control
 */
enum class IntegrationAlgorithm : int
{
    MD,
    Steep,
    CG,
    BD,
    SD2Removed,
    NM,
    LBFGS,
    TPI,
    TPIC,
    SD1,
    VV,
    VVAK,
    Mimic,
    Count,
    Default = MD
};
//! Name of the integrator algorithm
const char* enumValueToString(IntegrationAlgorithm enumValue);
//! Do we use MiMiC QM/MM?
#define EI_MIMIC(e) ((e) == IntegrationAlgorithm::Mimic)
//! Do we use velocity Verlet
#define EI_VV(e) ((e) == IntegrationAlgorithm::VV || (e) == IntegrationAlgorithm::VVAK)
//! Do we use molecular dynamics
#define EI_MD(e) ((e) == IntegrationAlgorithm::MD || EI_VV(e) || EI_MIMIC(e))
//! Do we use stochastic dynamics
#define EI_SD(e) ((e) == IntegrationAlgorithm::SD1)
//! Do we use any stochastic integrator
#define EI_RANDOM(e) (EI_SD(e) || (e) == IntegrationAlgorithm::BD)
/*above integrators may not conserve momenta*/
//! Do we use any type of dynamics
#define EI_DYNAMICS(e) (EI_MD(e) || EI_RANDOM(e))
//! Or do we use minimization
#define EI_ENERGY_MINIMIZATION(e)                                          \
    ((e) == IntegrationAlgorithm::Steep || (e) == IntegrationAlgorithm::CG \
     || (e) == IntegrationAlgorithm::LBFGS)
//! Do we apply test particle insertion
#define EI_TPI(e) ((e) == IntegrationAlgorithm::TPI || (e) == IntegrationAlgorithm::TPIC)
//! Do we deal with particle velocities
#define EI_STATE_VELOCITY(e) (EI_MD(e) || EI_SD(e))

//! Constraint algorithm
enum class ConstraintAlgorithm : int
{
    Lincs,
    Shake,
    Count,
    Default = Lincs
};
//! String corresponding to constraint algorithm
const char* enumValueToString(ConstraintAlgorithm enumValue);

//! Distance restraint refinement algorithm
enum class DistanceRestraintRefinement : int
{
    None,
    Simple,
    Ensemble,
    Count,
    Default = None
};
//! String corresponding to distance restraint algorithm
const char* enumValueToString(DistanceRestraintRefinement enumValue);

//! Distance restraints weighting type
enum class DistanceRestraintWeighting : int
{
    Conservative,
    Equal,
    Count,
    Default = Conservative
};
//! String corresponding to distance restraint weighting
const char* enumValueToString(DistanceRestraintWeighting enumValue);

//! Combination rule algorithm.
enum class CombinationRule : int
{
    None,
    Geometric,
    Arithmetic,
    GeomSigEps,
    Count,
    Default = Geometric
};
//! String for combination rule algorithm
const char* enumValueToString(CombinationRule enumValue);

//! Van der Waals potential.
enum class VanDerWaalsPotential : int
{
    None,
    LJ,
    Buckingham,
    Count,
    Default = LJ
};
//! String corresponding to Van der Waals potential
const char* enumValueToString(VanDerWaalsPotential enumValue);

//! Simulated tempering methods.
enum class SimulatedTempering : int
{
    Geometric,
    Exponential,
    Linear,
    Count,
    Default = Geometric
};
//! String corresponding to simulated tempering
const char* enumValueToString(SimulatedTempering enumValue);

/*! \brief Free energy perturbation type
 */
enum class FreeEnergyPerturbationType : int
{
    //! there are no evaluations at other states
    No,
    //! treated equivalently to Static
    Yes,
    //! then lambdas do not change during the simulation
    Static,
    //! then the states change monotonically throughout the simulation
    SlowGrowth,
    //! then expanded ensemble simulations are occurring
    Expanded,
    Count,
    Default = No
};
//! String corresponding to FEP type.
const char* enumValueToString(FreeEnergyPerturbationType enumValue);

//! Free energy perturbation coupling types.
enum class FreeEnergyPerturbationCouplingType : int
{
    Fep,
    Mass,
    Coul,
    Vdw,
    Bonded,
    Restraint,
    Temperature,
    Count,
    Default = Fep
};
//! String for FEP coupling type
const char* enumValueToString(FreeEnergyPerturbationCouplingType enumValue);
//! String for FEP coupling type, singular mention.
const char* enumValueToStringSingular(FreeEnergyPerturbationCouplingType enumValue);

/*! \brief What to print for free energy calculations
 *
 * Printing the energy to the free energy dhdl file.
 * Yes is an alias to Total, and
 * will be converted in readir, so we never have to account for it in code.
 */
enum class FreeEnergyPrintEnergy : int
{
    No,
    Total,
    Potential,
    Yes,
    Count,
    Default = No
};
//! String corresponding to printing of free energy
const char* enumValueToString(FreeEnergyPrintEnergy enumValue);

/*! \brief How the lambda weights are calculated
 */
enum class LambdaWeightCalculation : int
{
    //! don't calculate
    No,
    //! using the metropolis criteria
    Metropolis,
    //! using the Barker critera for transition weights, also called unoptimized Bennett
    Barker,
    //! using Barker + minimum variance for weights
    Minvar,
    //! Wang-Landu (using visitation counts)
    WL,
    //! Weighted Wang-Landau (using optimized Gibbs weighted visitation counts)
    WWL,
    Count,
    Default = No
};
//! String corresponding to lambda weights
const char* enumValueToString(LambdaWeightCalculation enumValue);
//! Macro telling us whether we use expanded ensemble
#define ELAMSTATS_EXPANDED(e) ((e) > LambdaWeightCalculation::No)
//! Macro telling us whether we use some kind of Wang-Landau
#define EWL(e) ((e) == LambdaWeightCalculation::WL || (e) == LambdaWeightCalculation::WWL)

/*! \brief How moves in lambda are calculated
 */
enum class LambdaMoveCalculation : int
{
    //! don't calculate move
    No,
    //! using the Metropolis criteria, and 50% up and down
    Metropolis,
    //! using the Barker criteria, and 50% up and down
    Barker,
    //! computing the transition using the marginalized probabilities of the lambdas
    Gibbs,
    /*! \brief
     * using the metropolized version of Gibbs
     *
     * Monte Carlo Strategies in Scientific computing, Liu, p. 134
     */
    MetropolisGibbs,
    Count,
    Default = No
};
//! String corresponding to lambda moves
const char* enumValueToString(LambdaMoveCalculation enumValue);

/*! \brief How we decide whether weights have reached equilibrium
 */
enum class LambdaWeightWillReachEquilibrium : int
{
    //! never stop, weights keep going
    No,
    //! fix the weights from the beginning; no movement
    Yes,
    //! stop when the WL-delta falls below a certain level
    WLDelta,
    //! stop when we have a certain number of samples at every step
    NumAtLambda,
    //! stop when we've run a certain total number of steps
    Steps,
    //! stop when we've run a certain total number of samples
    Samples,
    //! stop when the ratio of samples (lowest to highest) is sufficiently large
    Ratio,
    Count,
    Default = No
};
//! String corresponding to equilibrium algorithm
const char* enumValueToString(LambdaWeightWillReachEquilibrium enumValue);

/*! \brief separate_dhdl_file selection
 *
 * NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool
 * Why was this done this way, just .........
 */
enum class SeparateDhdlFile : int
{
    Yes,
    No,
    Count,
    Default = Yes
};
//! String corresponding to separate DHDL file selection
const char* enumValueToString(SeparateDhdlFile enumValue);

/*! \brief dhdl_derivatives selection \
 *
 * NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool
 * Why was this done this way, just .........
 */
enum class DhDlDerivativeCalculation : int
{
    Yes,
    No,
    Count,
    Default = Yes
};
//! String for DHDL derivatives
const char* enumValueToString(DhDlDerivativeCalculation enumValue);

/*! \brief soft-core function \
 *
 * Distinguishes between soft-core functions in the input.
 */
enum class SoftcoreType : int
{
    Beutler,
    Gapsys,
    Count,
    Default = Beutler
};
//! Strings for softcore function names
const char* enumValueToString(SoftcoreType enumValue);

/*! \brief soft-core function as parameter to the nb-fep kernel/14-interaction.\
 *
 * Distinguishes between soft-core functions internally. This is different
 * from SoftcoreType in that it offers 'None' which is not exposed to the user.
 */
enum class KernelSoftcoreType : int
{
    Beutler,
    Gapsys,
    None,
    Count,
    Default = Beutler
};
//! Strings for softcore function names
const char* enumValueToString(KernelSoftcoreType enumValue);

/*! \brief Solvent model
 *
 * Distinguishes classical water types with 3 or 4 particles
 */
enum class SolventModel : int
{
    No,
    Spc,
    Tip4p,
    Count,
    Default = Spc
};
//! String corresponding to solvent type
const char* enumValueToString(SolventModel enumValue);

//! Dispersion correction.
enum class DispersionCorrectionType : int
{
    No,
    EnerPres,
    Ener,
    AllEnerPres,
    AllEner,
    Count,
    Default = No
};
//! String corresponding to dispersion corrections
const char* enumValueToString(DispersionCorrectionType enumValue);

//! Algorithm for simulated annealing.
enum class SimulatedAnnealing : int
{
    No,
    Single,
    Periodic,
    Count,
    Default = No
};
//! String for simulated annealing
const char* enumValueToString(SimulatedAnnealing enumValue);

//! Wall types.
enum class WallType : int
{
    NineThree,
    TenFour,
    Table,
    TwelveSix,
    Count,
    Default = NineThree
};
//! String corresponding to wall type
const char* enumValueToString(WallType enumValue);

//! Pulling algorithm.
enum class PullingAlgorithm : int
{
    Umbrella,
    Constraint,
    ConstantForce,
    FlatBottom,
    FlatBottomHigh,
    External,
    Count,
    Default = Umbrella
};
//! String for pulling algorithm
const char* enumValueToString(PullingAlgorithm enumValue);

//! Control of pull groups
enum class PullGroupGeometry : int
{
    Distance,
    Direction,
    Cylinder,
    DirectionPBC,
    DirectionRelative,
    Angle,
    Dihedral,
    AngleAxis,
    Transformation,
    Count,
    Default = Distance
};
//! String for pull groups
const char* enumValueToString(PullGroupGeometry enumValue);

//! Enforced rotation groups.
enum class EnforcedRotationGroupType : int
{
    Iso,
    Isopf,
    Pm,
    Pmpf,
    Rm,
    Rmpf,
    Rm2,
    Rm2pf,
    Flex,
    Flext,
    Flex2,
    Flex2t,
    Count,
    Default = Iso
};
//! Rotation group names
const char* enumValueToString(EnforcedRotationGroupType enumValue);
//! String for rotation group origin names
const char* enumValueToLongString(EnforcedRotationGroupType enumValue);

//! Rotation group fitting type
enum class RotationGroupFitting : int
{
    Rmsd,
    Norm,
    Pot,
    Count,
    Default = Rmsd
};
//! String corresponding to rotation group fitting
const char* enumValueToString(RotationGroupFitting enumValue);

/*! \brief Direction along which ion/water swaps happen
 *
 * Part of "Computational Electrophysiology" (CompEL) setups
 */
enum class SwapType : int
{
    No,
    X,
    Y,
    Z,
    Count,
    Default = No
};
//! Names for swapping
const char* enumValueToString(SwapType enumValue);

/*! \brief Swap group splitting type
 *
 * These are just the fixed groups we need for any setup. In t_swap's grp
 * entry after that follows the variable number of swap groups.
 */
enum class SwapGroupSplittingType : int
{
    Split0,
    Split1,
    Solvent,
    Count,
    Default = Solvent
};
//! String for swap group splitting
const char* enumValueToString(SwapGroupSplittingType enumValue);

/*! \brief Types of electrostatics calculations
 *
 * Types of electrostatics calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections
 * in the MDP file; many interactions for instance map to tabulated kernels.
 */
enum class NbkernelElecType : int
{
    None,
    Coulomb,
    ReactionField,
    CubicSplineTable,
    Ewald,
    Count,
    Default = None
};
//! String corresponding to electrostatics kernels
const char* enumValueToString(NbkernelElecType enumValue);

/*! \brief Types of vdw calculations available
 *
 * Types of vdw calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections
 * in the MDP file; many interactions for instance map to tabulated kernels.
 */
enum class NbkernelVdwType : int
{
    None,
    LennardJones,
    Buckingham,
    CubicSplineTable,
    LJEwald,
    Count,
    Default = None
};
//! String corresponding to VdW kernels
const char* enumValueToString(NbkernelVdwType enumValue);

//! Center of mass motion removal algorithm.
enum class ComRemovalAlgorithm : int
{
    Linear,
    Angular,
    No,
    LinearAccelerationCorrection,
    Count,
    Default = Linear
};
//! String corresponding to COM removal
const char* enumValueToString(ComRemovalAlgorithm enumValue);

//! Enumeration that contains all supported periodic boundary setups.
enum class PbcType : int
{
    Xyz     = 0, //!< Periodic boundaries in all dimensions.
    No      = 1, //!< No periodic boundaries.
    XY      = 2, //!< Only two dimensions are periodic.
    Screw   = 3, //!< Screw.
    Unset   = 4, //!< The type of PBC is not set or invalid.
    Count   = 5,
    Default = Xyz
};

#endif /* GMX_MDTYPES_MD_ENUMS_H */
