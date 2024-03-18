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
        std::tuple<HarmonicBondType,   //    F_BONDS,
                   G96BondType,        //    F_G96BONDS,
                   MorseBondType,      //    F_MORSE,
                   CubicBondType,      //    F_CUBICBONDS,
                   Unimplemented,      //    F_CONNBONDS,
                   NotInNblibButInGMX, //    F_HARMONIC,
                   FENEBondType,       //    F_FENEBONDS,
                   NotInNblibButInGMX, //    F_TABBONDS,
                   NotInNblibButInGMX, //    F_TABBONDSNC,
                   NotInNblibButInGMX, //    F_RESTRBONDS,
                   HarmonicAngle,      //    F_ANGLES
                   G96Angle,           //    F_G96ANGLES
                   RestrictedAngle,    //    F_RESTRANGLES,
                   LinearAngle,        //    F_LINEAR_ANGLES,
                   CrossBondBond,      //    F_CROSS_BOND_BONDS,
                   CrossBondAngle,     //    F_CROSS_BOND_ANGLES,
                   // TODO: Implement custom param transfer function for this
                   NotInNblibButInGMX,       //    F_UREY_BRADLEY,
                   QuarticAngle,             //    F_QUARTIC_ANGLES,
                   NotInNblibButInGMX,       //    F_TABANGLES,
                   ProperDihedral,           //    F_PDIHS,
                   RyckaertBellemanDihedral, //    F_RBDIHS,
                   NotInNblibButInGMX,       //    F_RESTRDIHS,
                   NotInNblibButInGMX,       //    F_CBTDIHS,
                   NotInNblibButInGMX,       //    F_FOURDIHS,
                   ImproperDihedral,         //    F_IDIHS,
                   NotInNblibButInGMX,       //    F_PIDIHS,
                   NotInNblibButInGMX,       //    F_TABDIHS,
                   Unimplemented,            //    F_CMAP,
                   Unimplemented,            //    F_GB12_NOLONGERUSED,
                   Unimplemented,            //    F_GB13_NOLONGERUSED,
                   Unimplemented,            //    F_GB14_NOLONGERUSED,
                   Unimplemented,            //    F_GBPOL_NOLONGERUSED,
                   Unimplemented,            //    F_NPSOLVATION_NOLONGERUSED,
                   PairLJType,               //    F_LJ14,
                   Unimplemented,            //    F_COUL14,
                   Unimplemented,            //    F_LJC14_Q,
                   Unimplemented,            //    F_LJC_PAIRS_NB,
                   Unimplemented,            //    F_LJ,
                   Unimplemented,            //    F_BHAM,
                   Unimplemented,            //    F_LJ_LR_NOLONGERUSED,
                   Unimplemented,            //    F_BHAM_LR_NOLONGERUSED,
                   Unimplemented,            //    F_DISPCORR,
                   Unimplemented,            //    F_COUL_SR,
                   Unimplemented,            //    F_COUL_LR_NOLONGERUSED,
                   Unimplemented,            //    F_RF_EXCL,
                   Unimplemented,            //    F_COUL_RECIP,
                   Unimplemented,            //    F_LJ_RECIP,
                   Unimplemented,            //    F_DPD,
                   NotInNblibButInGMX,       //    F_POLARIZATION,
                   NotInNblibButInGMX,       //    F_WATER_POL,
                   NotInNblibButInGMX,       //    F_THOLE_POL,
                   NotInNblibButInGMX,       //    F_ANHARM_POL,
                   Unimplemented,            //    F_POSRES,
                   Unimplemented,            //    F_FBPOSRES,
                   NotInNblibButInGMX,       //    F_DISRES,
                   Unimplemented,            //    F_DISRESVIOL,
                   NotInNblibButInGMX,       //    F_ORIRES,
                   Unimplemented,            //    F_ORIRESDEV,
                   NotInNblibButInGMX,       //    F_ANGRES,
                   NotInNblibButInGMX,       //    F_ANGRESZ,
                   NotInNblibButInGMX,       //    F_DIHRES,
                   Unimplemented,            //    F_DIHRESVIOL,
                   Unimplemented,            //    F_CONSTR,
                   Unimplemented,            //    F_CONSTRNC,
                   Unimplemented,            //    F_SETTLE,
                   Unimplemented,            //    F_VSITE1,
                   Unimplemented,            //    F_VSITE2,
                   Unimplemented,            //    F_VSITE2FD,
                   Unimplemented,            //    F_VSITE3,
                   Unimplemented,            //    F_VSITE3FD,
                   Unimplemented,            //    F_VSITE3FAD,
                   Unimplemented,            //    F_VSITE3OUT,
                   Unimplemented,            //    F_VSITE4FD,
                   Unimplemented,            //    F_VSITE4FDN,
                   Unimplemented,            //    F_VSITEN,
                   Unimplemented,            //    F_COM_PULL,
                   Unimplemented,            //    F_DENSITYFITTING,
                   Unimplemented,            //    F_EQM,
                   Unimplemented,            //    F_EPOT,
                   Unimplemented,            //    F_EKIN,
                   Unimplemented,            //    F_ETOT,
                   Unimplemented,            //    F_ECONSERVED,
                   Unimplemented,            //    F_TEMP,
                   Unimplemented,            //    F_VTEMP_NOLONGERUSED,
                   Unimplemented,            //    F_PDISPCORR,
                   Unimplemented,            //    F_PRES,
                   Unimplemented,            //    F_DVDL_CONSTR,
                   Unimplemented,            //    F_DVDL,
                   Unimplemented,            //    F_DKDL,
                   Unimplemented,            //    F_DVDL_COUL,
                   Unimplemented,            //    F_DVDL_VDW,
                   Unimplemented,            //    F_DVDL_BONDED,
                   Unimplemented,            //    F_DVDL_RESTRAINT,
                   Unimplemented,            //    F_DVDL_TEMPERATURE,
                   /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
                   Unimplemented //    F_NRE
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
