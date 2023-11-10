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
 * GMX->NBLIB and NBLIB->GMX listed conversion utilities
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_CONVERSIONSCOMMON_H
#define NBLIB_LISTEDFORCES_CONVERSIONSCOMMON_H

#include <memory>
#include <tuple>
#include <type_traits>

#include "gromacs/topology/idef.h"

#include "nblib/listed_forces/definitions.h"
#include "nblib/listed_forces/traits.h"

struct gmx_ffparams_t;

namespace nblib
{
/*! \brief Interactions not supported by NBLIB but there in GROMACS
 *
 * Dummy no op interaction as placeholder
 * V(r, forceConstant, equilConstant) = 0;
 */
class NotInNblibButInGMX
{
};

/*! \brief Interactions not supported by both NBLIB and GROMACS
 *
 * Dummy no op interaction as placeholder
 * V(r, forceConstant, equilConstant) = 0;
 */
class Unimplemented
{
};

using GmxToNblibMapping =
        std::tuple<HarmonicBondType,         //    0 F_BONDS,
                   G96BondType,              //    1 F_G96BONDS,
                   MorseBondType,            //    2 F_MORSE,
                   CubicBondType,            //    3 F_CUBICBONDS,
                   Unimplemented,            //    4 F_CONNBONDS,
                   NotInNblibButInGMX,       //    5 F_HARMONIC,
                   FENEBondType,             //    6 F_FENEBONDS,
                   NotInNblibButInGMX,       //    7 F_TABBONDS,
                   NotInNblibButInGMX,       //    8 F_TABBONDSNC,
                   NotInNblibButInGMX,       //    9 F_RESTRBONDS,
                   HarmonicAngle,            //   10 F_ANGLES
                   G96Angle,                 //   11 F_G96ANGLES
                   RestrictedAngle,          //   12 F_RESTRANGLES,
                   LinearAngle,              //   13 F_LINEAR_ANGLES,
                   CrossBondBond,            //   14 F_CROSS_BOND_BONDS,
                   CrossBondAngle,           //   15 F_CROSS_BOND_ANGLES,
                   UreyBradley,              //   16 F_UREY_BRADLEY,
                   QuarticAngle,             //   17 F_QUARTIC_ANGLES,
                   NotInNblibButInGMX,       //   18 F_TABANGLES,
                   ProperDihedral,           //   19 F_PDIHS,
                   RyckaertBellemanDihedral, //   20 F_RBDIHS,
                   NotInNblibButInGMX,       //   21 F_RESTRDIHS,
                   NotInNblibButInGMX,       //   22 F_CBTDIHS,
                   NotInNblibButInGMX,       //   23 F_FOURDIHS,
                   ImproperDihedral,         //   24 F_IDIHS,
                   ImproperProperDihedral,   //   25 F_PIDIHS,
                   NotInNblibButInGMX,       //   26 F_TABDIHS,
                   Unimplemented,            //   27 F_CMAP,
                   Unimplemented,            //   28 F_GB12_NOLONGERUSED,
                   Unimplemented,            //   29 F_GB13_NOLONGERUSED,
                   Unimplemented,            //   30 F_GB14_NOLONGERUSED,
                   Unimplemented,            //   31 F_GBPOL_NOLONGERUSED,
                   Unimplemented,            //   32 F_NPSOLVATION_NOLONGERUSED,
                   PairLJType,               //   33 F_LJ14,
                   Unimplemented,            //   34 F_COUL14,
                   PairLJChargeType,         //   35 F_LJC14_Q,
                   Unimplemented,            //   36 F_LJC_PAIRS_NB,
                   Unimplemented,            //   37 F_LJ,
                   Unimplemented,            //   38 F_BHAM,
                   Unimplemented,            //   39 F_LJ_LR_NOLONGERUSED,
                   Unimplemented,            //   40 F_BHAM_LR_NOLONGERUSED,
                   Unimplemented,            //   41 F_DISPCORR,
                   Unimplemented,            //   42 F_COUL_SR,
                   Unimplemented,            //   43 F_COUL_LR_NOLONGERUSED,
                   Unimplemented,            //   44 F_RF_EXCL,
                   Unimplemented,            //   45 F_COUL_RECIP,
                   Unimplemented,            //   46 F_LJ_RECIP,
                   Unimplemented,            //   47 F_DPD,
                   SimplePolarization,       //   48 F_POLARIZATION,
                   NotInNblibButInGMX,       //   49 F_WATER_POL,
                   NotInNblibButInGMX,       //   50 F_THOLE_POL,
                   NotInNblibButInGMX,       //   51 F_ANHARM_POL,
                   PositionRestraints,       //   52 F_POSRES,
                   Unimplemented,            //   53 F_FBPOSRES,
                   NotInNblibButInGMX,       //   54 F_DISRES,
                   Unimplemented,            //   55 F_DISRESVIOL,
                   NotInNblibButInGMX,       //   56 F_ORIRES,
                   Unimplemented,            //   57 F_ORIRESDEV,
                   NotInNblibButInGMX,       //   58 F_ANGRES,
                   NotInNblibButInGMX,       //   59 F_ANGRESZ,
                   NotInNblibButInGMX,       //   60 F_DIHRES,
                   Unimplemented,            //   61 F_DIHRESVIOL,
                   Unimplemented,            //   62 F_CONSTR,
                   Unimplemented,            //   64 F_CONSTRNC,
                   Unimplemented,            //   65 F_SETTLE,
                   Unimplemented,            //   66 F_VSITE1,
                   Unimplemented,            //   67 F_VSITE2,
                   Unimplemented,            //   68 F_VSITE2FD,
                   Unimplemented,            //   69 F_VSITE3,
                   Unimplemented,            //   70 F_VSITE3FD,
                   Unimplemented,            //   71 F_VSITE3FAD,
                   Unimplemented,            //   72 F_VSITE3OUT,
                   Unimplemented,            //   73 F_VSITE4FD,
                   Unimplemented,            //   74 F_VSITE4FDN,
                   Unimplemented,            //   75 F_VSITEN,
                   Unimplemented,            //   76 F_COM_PULL,
                   Unimplemented,            //   77 F_DENSITYFITTING,
                   Unimplemented,            //   78 F_EQM,
                   Unimplemented,            //   79 F_EPOT,
                   Unimplemented,            //   80 F_EKIN,
                   Unimplemented,            //   81 F_ETOT,
                   Unimplemented,            //   82 F_ECONSERVED,
                   Unimplemented,            //   83 F_TEMP,
                   Unimplemented,            //   84 F_VTEMP_NOLONGERUSED,
                   Unimplemented,            //   85 F_PDISPCORR,
                   Unimplemented,            //   86 F_PRES,
                   Unimplemented,            //   87 F_DVDL_CONSTR,
                   Unimplemented,            //   88 F_DVDL,
                   Unimplemented,            //   89 F_DKDL,
                   Unimplemented,            //   90 F_DVDL_COUL,
                   Unimplemented,            //   91 F_DVDL_VDW,
                   Unimplemented,            //   92 F_DVDL_BONDED,
                   Unimplemented,            //   93 F_DVDL_RESTRAINT,
                   Unimplemented,            //   94 F_DVDL_TEMPERATURE,
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
 * Supports all types implemented in NBLIB that are supported by GROMACS
 *
 * \param  NBLIB ListedInteractionData
 * \return interactionDefinitions and gmx_ffparams_t objects as used in GROMACS
 */
std::tuple<std::unique_ptr<InteractionDefinitions>, std::unique_ptr<gmx_ffparams_t>>
convertToGmxInteractions(const ListedInteractionData& interactions);

/*! \brief Function to convert from gmx-InteractionDefinitions to nblib-ListedInteractionData
 *
 * Supports only those types that are supported by both NBLIB and GROMACS. Primary usecase
 * to read data from TPR files using the GROMACS machinery
 *
 * \param   interactionDefinitions object as used in GROMACS
 * \return  NBLIB ListedInteractionData
 */
ListedInteractionData convertToNblibInteractions(const InteractionDefinitions& interactionDefinitions);

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERSIONSCOMMON_H
