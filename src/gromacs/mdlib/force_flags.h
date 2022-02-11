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
#ifndef GMX_MDLIB_FORCE_FLAGS_H
#define GMX_MDLIB_FORCE_FLAGS_H

/* Flags to tell the force calculation routines what (not) to do */

/* The state has changed, always set unless TPI is used. */
#define GMX_FORCE_STATECHANGED (1u << 0u)
/* The box might have changed */
#define GMX_FORCE_DYNAMICBOX (1u << 1u)
/* Do neighbor searching */
#define GMX_FORCE_NS (1u << 2u)
/* Calculate listed energies/forces (e.g. bonds, restraints, 1-4, FEP non-bonded) */
#define GMX_FORCE_LISTED (1u << 4u)
/* Calculate non-bonded energies/forces */
#define GMX_FORCE_NONBONDED (1u << 6u)
/* Calculate forces (not only energies) */
#define GMX_FORCE_FORCES (1u << 7u)
/* Calculate the virial */
#define GMX_FORCE_VIRIAL (1u << 8u)
/* Calculate energies */
#define GMX_FORCE_ENERGY (1u << 9u)
/* Calculate dHdl */
#define GMX_FORCE_DHDL (1u << 10u)
/* Tells whether only the MTS combined force buffer is needed and not the normal force buffer */
#define GMX_FORCE_DO_NOT_NEED_NORMAL_FORCE (1u << 11u)

/* Normally one want all energy terms and forces */
#define GMX_FORCE_ALLFORCES (GMX_FORCE_LISTED | GMX_FORCE_NONBONDED | GMX_FORCE_FORCES)

#endif
