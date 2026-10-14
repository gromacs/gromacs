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
 * Declares interaction function types and related things.
 *
 * \ingroup module_topology
 * \inlibraryapi
 */
#ifndef GMX_TOPOLOGY_IFUNC_H
#define GMX_TOPOLOGY_IFUNC_H

#include "gromacs/libgromacs_export.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/vectypes.h"

/* TODO: Remove this typedef when t_ilist is removed */
typedef int t_iatom;

/* Real vector type with an additional, unused 4th element.
 * This type is used to allow aligned 4-wide SIMD loads and stores.
 */
typedef real rvec4[4];

/*! \brief These flags tell some of the routines what can be done with
 * particular kinds of interaction functions.
 *
 * These are used to define the values used for
 * \c InteractionFunctionCategoryFlags.
 */
enum class InteractionFunctionCategory : int
{
    //!< Indicates the intent to set no category
    Null,
    //!< A bonded force interaction will be calculated, e.g. bond, angle, dihedral, etc.
    BondedForce,
    //!< A virtual site
    VirtualSite,
    //!< Either a constraint or a SETTLE
    Constraint,
    //!< A bond, constraint, or SETTLE that expresses connectivity
    ChemicalBond,
    /*!< \brief A two-body bond that grompp can manipulate in topology
     * preparation.
     *
     * Originally IF_BYTPE was introduced to describe all kinds of bonds that
     * were then supported in GROMACS, and used in preparing virtual sites
     * and in conversions to either constraints or Morse bonds.  Later
     * other kinds of two-body interactions were introduced, and sometimes
     * they were added to sc_ifcBType. Then the descriptor sc_ifcChemicalBond
     * was introduced to describe all kinds of connective two-body
     * associations, including constraints.
     *
     * Bond would be a better name, but that word was originally used
     * in IF_BOND to refer to all kinds of bonded interactions, so
     * should be avoided here.
     *
     * This category is used in preprocessing for handling the above three
     * roles.  There are other two-body bonded interactions that could
     * perhaps have this category added. Or those three roles could be
     * differentiated when we have a use case for that.
     *
     * This category should be avoided in mdrun, to avoid making the
     * mulptiple use cases even more complex. Prefer \c ChemicalBond
     * or checking the number of particles in the interaction. */
    BType,
    //!< Angle interactions
    Angle,
    //!< Dihedral interactions
    Dihedral,
    //!< Pair interactions
    Pair,
    //!< Tabulated interactions
    Tabulated,
    /*!< \brief Indicates that for a bonded interaction the potential
     * does goes to zero for large distances, thus if such an interaction
     * is not assigned to any rank by the domain decomposition, the simulation
     * may continue, when mdrun has been told to permit that.
     */
    ZeroAtLargeDistances,
    //!< Conventional final delimiter for an enum class
    Count
};

/*! \brief Type for specifying a set of interaction categories
 *
 * A std::bitset might be nice in C++23, but until then it is not
 * constexpr enough.
 */
using InteractionFunctionCategoryFlags = unsigned int;

/*! \brief Flags describing categories of interaction functions
 *
 * These values can be used (individually or in groups) in bitwise
 * operations to see e.g. if an interaction type is either a
 * constraint or a bond.
 */
//! \{
// clang-format off
static constexpr InteractionFunctionCategoryFlags sc_ifcNull =                 0;
static constexpr InteractionFunctionCategoryFlags sc_ifcBondedForce =          1 << static_cast<int>(InteractionFunctionCategory::BondedForce);
static constexpr InteractionFunctionCategoryFlags sc_ifcVirtualSite =          1 << static_cast<int>(InteractionFunctionCategory::VirtualSite);
static constexpr InteractionFunctionCategoryFlags sc_ifcConstraint =           1 << static_cast<int>(InteractionFunctionCategory::Constraint);
static constexpr InteractionFunctionCategoryFlags sc_ifcChemicalBond =         1 << static_cast<int>(InteractionFunctionCategory::ChemicalBond);
static constexpr InteractionFunctionCategoryFlags sc_ifcBType =                1 << static_cast<int>(InteractionFunctionCategory::BType);
static constexpr InteractionFunctionCategoryFlags sc_ifcAngle =                1 << static_cast<int>(InteractionFunctionCategory::Angle);
static constexpr InteractionFunctionCategoryFlags sc_ifcDihedral =             1 << static_cast<int>(InteractionFunctionCategory::Dihedral);
static constexpr InteractionFunctionCategoryFlags sc_ifcPair =                 1 << static_cast<int>(InteractionFunctionCategory::Pair);
static constexpr InteractionFunctionCategoryFlags sc_ifcTabulated =            1 << static_cast<int>(InteractionFunctionCategory::Tabulated);
static constexpr InteractionFunctionCategoryFlags sc_ifcZeroAtLargeDistances = 1 << static_cast<int>(InteractionFunctionCategory::ZeroAtLargeDistances);
// clang-format on
//! \}

//! Interaction function.
struct t_interaction_function // NOLINT(clang-analyzer-optin.performance.Padding)
{
    const char* name;                            //!< The name of this function.
    const char* longname;                        //!< The name for printing etc.
    int         nratoms;                         //!< Number of atoms needed for this function.
    int         nrfpA, nrfpB;                    /*!< \brief Number of parameters for this function.
                                                  *
                                                  * This corresponds to the number of parameters in
                                                  * iparams struct.
                                                  *
                                                  * A and B are for normal and free energy components
                                                  * respectively.
                                                  */
    InteractionFunctionCategoryFlags categories; //!< Categories for this interaction (bond, constraint, etc.)
};

/*! \brief Interaction function enums.
 *
 * This MUST correspond to the t_interaction_function array in src/gromacs/topology/ifunc.cpp.
 */
enum class InteractionFunction : int
{
    Bonds,
    GROMOS96Bonds,
    MorsePotential,
    CubicBonds,
    ConnectBonds,
    HarmonicPotential,
    FENEBonds,
    TabulatedBonds,
    TabulatedBondsNoCoupling,
    RestraintBonds,
    Angles,
    GROMOS96Angles,
    RestrictedBendingPotential,
    LinearAngles,
    CrossBondBonds,
    CrossBondAngles,
    UreyBradleyPotential,
    QuarticAngles,
    TabulatedAngles,
    ProperDihedrals,
    RyckaertBellemansDihedrals,
    RestrictedTorsionPotential,
    CombinedBendingTorsionPotential,
    FourierDihedrals,
    ImproperDihedrals,
    PeriodicImproperDihedrals,
    TabulatedDihedrals,
    DihedralEnergyCorrectionMap,
    GeneralizedBorn12PolarizationUnused,
    GeneralizedBorn13PolarizationUnused,
    GeneralizedBorn14PolarizationUnused,
    GeneralizedBornPolarizationUnused,
    NonpolarSolvationUnused,
    LennardJones14,
    Coulomb14,
    LennardJonesCoulomb14Q,
    LennardJonesCoulombNonBondedPairs,
    LennardJonesShortRange,
    BuckinghamShortRange,
    LennardJonesLongRangeUnused,
    BuckinghamLongRangeUnused,
    DispersionCorrection,
    CoulombShortRange,
    CoulombLongRangeUnused,
    ReactionFieldExclusion,
    CoulombReciprocalSpace,
    LennardJonesReciprocalSpace,
    DissipativeParticleDynamics,
    Polarization,
    WaterPolarization,
    TholePolarization,
    AnharmonicPolarization,
    PositionRestraints,
    FlatBottomedPositionRestraints,
    DistanceRestraints,
    DistanceRestraintViolations,
    OrientationRestraints,
    OrientationRestraintDeviations,
    AngleRestraints,
    AngleZAxisRestraints,
    DihedralRestraints,
    DihedralRestraintViolations,
    Constraints,
    ConstraintsNoCoupling,
    SETTLE,
    VirtualSite1,
    VirtualSite2,
    VirtualSite2FlexibleDistance,
    VirtualSite3,
    VirtualSite3FlexibleDistance,
    VirtualSite3FlexibleAngleDistance,
    VirtualSite3Outside,
    VirtualSite4FlexibleDistance,
    VirtualSite4FlexibleDistanceNormalization,
    VirtualSiteN,
    CenterOfMassPullingEnergy,
    DensityFitting,
    QuantumMechanicalRegionEnergy,
    NeuralNetworkPotentialEnergy,
    PotentialEnergy,
    KineticEnergy,
    TotalEnergy,
    ConservedEnergy,
    Temperature,
    VirialTemperatureUnused,
    PressureDispersionCorrection,
    Pressure,
    dHdLambdaConstraint,
    dVremainingdLambda,
    dEkineticdLambda,
    dVCoulombdLambda,
    dVvanderWaalsdLambda,
    dVbondeddLambda,
    dVrestraintdLambda,
    dVtemperaturedLambda, //!< Not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE).
    Count //!< This number is for the total number of energies.
};

static inline bool IS_RESTRAINT_TYPE(InteractionFunction ifunc)
{
    return ifunc == InteractionFunction::PositionRestraints
           || ifunc == InteractionFunction::FlatBottomedPositionRestraints
           || ifunc == InteractionFunction::DistanceRestraints || ifunc == InteractionFunction::RestraintBonds
           || ifunc == InteractionFunction::DistanceRestraintViolations
           || ifunc == InteractionFunction::OrientationRestraints
           || ifunc == InteractionFunction::OrientationRestraintDeviations
           || ifunc == InteractionFunction::AngleRestraints
           || ifunc == InteractionFunction::AngleZAxisRestraints
           || ifunc == InteractionFunction::DihedralRestraints;
}

//! Maximum allowed number of atoms
constexpr int MAXATOMLIST = 6;
//! Maximum allowed number of parameters
constexpr int MAXFORCEPARAM = 12;
//! Maximum terms in interaction_function. Check src/gromacs/gmxpreprocess/toppush.cpp when you change these numbers.
constexpr int NR_RBDIHS   = 6;
constexpr int NR_CBTDIHS  = 6;
constexpr int NR_FOURDIHS = 4;

//! Initialised interaction functions descriptor.
LIBGROMACS_EXPORT extern const gmx::EnumerationArray<InteractionFunction, t_interaction_function> interaction_function;

static inline int NRFPA(InteractionFunction iftype)
{
    return interaction_function[iftype].nrfpA;
}

static inline int NRFPB(InteractionFunction iftype)
{
    return interaction_function[iftype].nrfpB;
}

static inline int NRFP(InteractionFunction iftype)
{
    return NRFPA(iftype) + NRFPB(iftype);
}

static inline int NRAL(InteractionFunction iftype)
{
    return interaction_function[iftype].nratoms;
}

//! \brief Tells if function type ftype represents a two-bond Lewis-style chemical bond.
static inline bool IS_CHEMBOND(InteractionFunction iftype)
{
    // Note that SETTLE has sc_ifcChemicalBond, so we need to check for the
    // atom count.
    return interaction_function[iftype].nratoms == 2
           && (interaction_function[iftype].categories & sc_ifcChemicalBond) != 0;
}

/*! \brief Tells if a function type ftype represents an angle.
 *
 * Per Larsson, 2007-11-06.
 */
static inline bool IS_ANGLE(InteractionFunction iftype)
{
    return interaction_function[iftype].nratoms == 3
           && (interaction_function[iftype].categories & sc_ifcAngle) != 0;
}

static inline bool IS_VSITE(InteractionFunction iftype)
{
    return (interaction_function[iftype].categories & sc_ifcVirtualSite) != 0;
}

static inline bool IS_TABULATED(InteractionFunction iftype)
{
    return (interaction_function[iftype].categories & sc_ifcTabulated) != 0;
}

#endif
