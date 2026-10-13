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
 */
//! \{
constexpr unsigned int IF_NULL = 0;
//! With IF_BOND a bonded interaction will be calculated.
// TODO rename me to IF_BONDEDFORCE for clarity wrt a 2-body bond
constexpr unsigned int IF_BOND = 1 << 0;
//! A virtual site
constexpr unsigned int IF_VSITE = 1 << 1;
//! Either a constraint or a SETTLE
constexpr unsigned int IF_CONSTRAINT = 1 << 2;
//! A bond, constraint, or SETTLE that expresses connectivity
constexpr unsigned int IF_CHEMBOND = 1 << 3;
/*! \brief A two-body bond that grompp can convert to a Morse potential.
 *
 * \todo rename me to IF_CONVERTIBLE_TO_MORSE for clarity
 *
 * With IF_BTYPE or IF_ATYPE the bond/angle can be converted to
 * a constraint or used for vsite parameter determination by grompp.
 *
 * \todo express this as IF_CONVERTIBLE_TO_CONSTRAINT and IF_VSITE_BASIS
 */
constexpr unsigned int IF_BTYPE = 1 << 4;
/*! \brief Angle interactions
 *
 * With IF_BTYPE or IF_ATYPE the bond/angle can be converted to
 * a constraint or used for vsite parameter determination by grompp.
 *
 * \todo express this as IF_CONVERTIBLE_TO_CONSTRAINT and IF_VSITE_BASIS
 */
constexpr unsigned int IF_ATYPE = 1 << 5;
//! Dihedral interactions
constexpr unsigned int IF_DIHEDRAL = 1 << 6;
//! Pair interactions
constexpr unsigned int IF_PAIR = 1 << 7;
//! Tabulated interactions
constexpr unsigned int IF_TABULATED = 1 << 8;
/*! \brief IF_LIMZERO indicates that for a bonded interaction the potential
 * does goes to zero for large distances.
 *
 * Thus if such an interaction it not assigned to any domain by the
 * domain decomposition, the simulation still continues, if mdrun has
 * been told to do so.
 */
constexpr unsigned int IF_LIMZERO = 1 << 9;
//! \}

//! Interaction function.
struct t_interaction_function // NOLINT(clang-analyzer-optin.performance.Padding)
{
    const char* name;         //!< The name of this function.
    const char* longname;     //!< The name for printing etc.
    int         nratoms;      //!< Number of atoms needed for this function.
    int         nrfpA, nrfpB; /*!< \brief Number of parameters for this function.
                               *
                               * This corresponds to the number of parameters in
                               * iparams struct.
                               *
                               * A and B are for normal and free energy components
                               * respectively.
                               */
    unsigned int flags;       //!< Flags (IF_NULL, IF_BOND, etc.).
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
    // Note that SETTLE has IF_CHEMBOND, so we need to check for the
    // atom count.
    return interaction_function[iftype].nratoms == 2
           && (interaction_function[iftype].flags & IF_CHEMBOND) != 0;
}

/*! \brief Tells if a function type ftype represents an angle.
 *
 * Per Larsson, 2007-11-06.
 */
static inline bool IS_ANGLE(InteractionFunction iftype)
{
    return interaction_function[iftype].nratoms == 3
           && (interaction_function[iftype].flags & IF_ATYPE) != 0;
}

static inline bool IS_VSITE(InteractionFunction iftype)
{
    return (interaction_function[iftype].flags & IF_VSITE) != 0;
}

static inline bool IS_TABULATED(InteractionFunction iftype)
{
    return (interaction_function[iftype].flags & IF_TABULATED) != 0;
}

#endif
