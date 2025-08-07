/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * This implements common definitions used in both
 * GMX->NB-LIB and NB-LIB->GMX listed conversion utilities
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_CONVERSIONSCOMMON_H
#define NBLIB_LISTEDFORCES_CONVERSIONSCOMMON_H

#include <memory>

#include "gromacs/topology/idef.h"

#include "nblib/listed_forces/definitions.h"

#include "traits.h"

struct gmx_ffparams_t;

namespace nblib
{
/*! \brief Interactions not supported by NB-LIB but there in GROMACS
 *
 * Dummy no op interaction as placeholder
 * V(r, forceConstant, equilConstant) = 0;
 */
class NotInNblibButInGMX
{
};

/*! \brief Interactions not supported by both NB-LIB and GROMACS
 *
 * Dummy no op interaction as placeholder
 * V(r, forceConstant, equilConstant) = 0;
 */
class Unimplemented
{
};

using GmxToNblibMapping =
        std::tuple<HarmonicBondType,   //    InteractionFunction::Bonds,
                   G96BondType,        //    InteractionFunction::GROMOS96Bonds,
                   MorseBondType,      //    InteractionFunction::MorsePotential,
                   CubicBondType,      //    InteractionFunction::CubicBonds,
                   Unimplemented,      //    InteractionFunction::ConnectBonds,
                   NotInNblibButInGMX, //    InteractionFunction::HarmonicPotential,
                   FENEBondType,       //    InteractionFunction::FENEBonds,
                   NotInNblibButInGMX, //    InteractionFunction::TabulatedBonds,
                   NotInNblibButInGMX, //    InteractionFunction::TabulatedBondsNoCoupling,
                   NotInNblibButInGMX, //    InteractionFunction::RestraintBonds,
                   HarmonicAngle,      //    InteractionFunction::Angles
                   G96Angle,           //    InteractionFunction::GROMOS96Angles
                   RestrictedAngle,    //    InteractionFunction::RestrictedBendingPotential,
                   LinearAngle,        //    InteractionFunction::LinearAngles,
                   CrossBondBond,      //    InteractionFunction::CrossBondBonds,
                   CrossBondAngle,     //    InteractionFunction::CrossBondAngles,
                   // TODO: Implement custom param transfer function for this
                   NotInNblibButInGMX,       //    InteractionFunction::UreyBradleyPotential,
                   QuarticAngle,             //    InteractionFunction::QuarticAngles,
                   NotInNblibButInGMX,       //    InteractionFunction::TabulatedAngles,
                   ProperDihedral,           //    InteractionFunction::ProperDihedrals,
                   RyckaertBellemanDihedral, //    InteractionFunction::RyckaertBellemansDihedrals,
                   NotInNblibButInGMX,       //    InteractionFunction::RestrictedTorsionPotential,
                   NotInNblibButInGMX, //    InteractionFunction::CombinedBendingTorsionPotential,
                   NotInNblibButInGMX, //    InteractionFunction::FourierDihedrals,
                   ImproperDihedral,   //    InteractionFunction::ImproperDihedrals,
                   NotInNblibButInGMX, //    InteractionFunction::PeriodicImproperDihedrals,
                   NotInNblibButInGMX, //    InteractionFunction::TabulatedDihedrals,
                   Unimplemented,      //    InteractionFunction::DihedralEnergyCorrectionMap,
                   Unimplemented, //    InteractionFunction::GeneralizedBorn12PolarizationUnused,
                   Unimplemented, //    InteractionFunction::GeneralizedBorn13PolarizationUnused,
                   Unimplemented, //    InteractionFunction::GeneralizedBorn14PolarizationUnused,
                   Unimplemented, //    InteractionFunction::GeneralizedBornPolarizationUnused,
                   Unimplemented, //    InteractionFunction::NonpolarSolvationUnused,
                   PairLJType,    //    InteractionFunction::LennardJones14,
                   Unimplemented, //    InteractionFunction::Coulomb14,
                   Unimplemented, //    InteractionFunction::LennardJonesCoulomb14Q,
                   Unimplemented, //    InteractionFunction::LennardJonesCoulombNonBondedPairs,
                   Unimplemented, //    InteractionFunction::LennardJonesShortRange,
                   Unimplemented, //    InteractionFunction::BuckinghamShortRange,
                   Unimplemented, //    InteractionFunction::LennardJonesLongRangeUnused,
                   Unimplemented, //    InteractionFunction::BuckinghamShortRange_LR_NOLONGERUSED,
                   Unimplemented, //    InteractionFunction::DispersionCorrection,
                   Unimplemented, //    InteractionFunction::CoulombShortRange,
                   Unimplemented, //    InteractionFunction::CoulombLongRangeUnused,
                   Unimplemented, //    InteractionFunction::ReactionFieldExclusion,
                   Unimplemented, //    InteractionFunction::CoulombReciprocalSpace,
                   Unimplemented, //    InteractionFunction::LennardJonesReciprocalSpace,
                   Unimplemented, //    InteractionFunction::DissipativeParticleDynamics,
                   NotInNblibButInGMX, //    InteractionFunction::Polarization,
                   NotInNblibButInGMX, //    InteractionFunction::WaterPolarization,
                   NotInNblibButInGMX, //    InteractionFunction::TholePolarization,
                   NotInNblibButInGMX, //    InteractionFunction::AnharmonicPolarization,
                   Unimplemented,      //    InteractionFunction::PositionRestraints,
                   Unimplemented,      //    InteractionFunction::FlatBottomedPositionRestraints,
                   NotInNblibButInGMX, //    InteractionFunction::DistanceRestraints,
                   Unimplemented,      //    InteractionFunction::DistanceRestraintViolations,
                   NotInNblibButInGMX, //    InteractionFunction::OrientationRestraints,
                   Unimplemented,      //    InteractionFunction::OrientationRestraintDeviations,
                   NotInNblibButInGMX, //    InteractionFunction::AngleRestraints,
                   NotInNblibButInGMX, //    InteractionFunction::AngleZAxisRestraints,
                   NotInNblibButInGMX, //    InteractionFunction::DihedralRestraints,
                   Unimplemented,      //    InteractionFunction::DihedralRestraintViolations,
                   Unimplemented,      //    InteractionFunction::Constraints,
                   Unimplemented,      //    InteractionFunction::ConstraintsNoCoupling,
                   Unimplemented,      //    InteractionFunction::SETTLE,
                   Unimplemented,      //    InteractionFunction::VirtualSite1,
                   Unimplemented,      //    InteractionFunction::VirtualSite2,
                   Unimplemented,      //    InteractionFunction::VirtualSite2FlexibleDistance,
                   Unimplemented,      //    InteractionFunction::VirtualSite3,
                   Unimplemented,      //    InteractionFunction::VirtualSite3FlexibleDistance,
                   Unimplemented,      //    InteractionFunction::VirtualSite3FlexibleAngleDistance,
                   Unimplemented,      //    InteractionFunction::VirtualSite3Outside,
                   Unimplemented,      //    InteractionFunction::VirtualSite4FlexibleDistance,
                   Unimplemented, //    InteractionFunction::VirtualSite4FlexibleDistanceNormalization,
                   Unimplemented, //    InteractionFunction::VirtualSiteN,
                   Unimplemented, //    InteractionFunction::CenterOfMassPullingEnergy,
                   Unimplemented, //    InteractionFunction::DensityFitting,
                   Unimplemented, //    InteractionFunction::QuantumMechanicalRegionEnergy,
                   Unimplemented, //    InteractionFunction::NeuralNetworkPotentialEnergy,
                   Unimplemented, //    InteractionFunction::PotentialEnergy,
                   Unimplemented, //    InteractionFunction::KineticEnergy,
                   Unimplemented, //    InteractionFunction::TotalEnergy,
                   Unimplemented, //    InteractionFunction::ConservedEnergy,
                   Unimplemented, //    InteractionFunction::Temperature,
                   Unimplemented, //    InteractionFunction::VirialTemperatureUnused,
                   Unimplemented, //    InteractionFunction::PressureDispersionCorrection,
                   Unimplemented, //    InteractionFunction::Pressure,
                   Unimplemented, //    InteractionFunction::dHdLambdaConstraint,
                   Unimplemented, //    InteractionFunction::dVremainingdLambda,
                   Unimplemented, //    InteractionFunction::dEkineticdLambda,
                   Unimplemented, //    InteractionFunction::dVCoulombdLambda,
                   Unimplemented, //    InteractionFunction::dVvanderWaalsdLambda,
                   Unimplemented, //    InteractionFunction::dVbondeddLambda,
                   Unimplemented, //    InteractionFunction::dVrestraintdLambda,
                   Unimplemented, //    InteractionFunction::dVtemperaturedLambda,
                   /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
                   Unimplemented //    InteractionFunction::Count
                                 /* This number is for the total number of energies */
                   >;

template<class InteractionType>
struct ListedTypeIsImplemented : std::bool_constant<Contains<InteractionType, GmxToNblibMapping>{}>
{
};

/*! \brief Function to convert from nblib-ListedInteractionData to gmx-InteractionDefinitions
 *
 * Supports all types implemented in NB-LIB that are supported by GROMACS
 *
 * \param  interactions NB-LIB ListedInteractionData
 * \return interactionDefinitions and gmx_ffparams_t objects as used in GROMACS
 */
std::tuple<std::unique_ptr<InteractionDefinitions>, std::unique_ptr<gmx_ffparams_t>>
convertToGmxInteractions(const ListedInteractionData& interactions);

/*! \brief Function to convert from gmx-InteractionDefinitions to nblib-ListedInteractionData
 *
 * Supports only those types that are supported by both NB-LIB and GROMACS. Primary usecase
 * to read data from TPR files using the GROMACS machinery
 *
 * \param   interactionDefinitions object as used in GROMACS
 * \return  NB-LIB ListedInteractionData
 */
ListedInteractionData convertToNblibInteractions(const InteractionDefinitions& interactionDefinitions);

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERSIONSCOMMON_H
