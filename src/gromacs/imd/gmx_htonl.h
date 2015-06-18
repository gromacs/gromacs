/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

/*! \libinternal
 * \defgroup module_imd Interactive molecular dynamics (IMD)
 * \ingroup group_mdrun
 *
 * \brief
 * Allows mdrun to interface with VMD via the interactive molecular dynamics
 * (IMD) protocol.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 */

/*! \libinternal \file
 *
 * \brief
 * This file contains function declarations for byte swapping.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 *
 * \inlibraryapi
 * \ingroup module_imd
 */
#ifndef GMX_IMD_HTONL_H
#define GMX_IMD_HTONL_H

#include "gromacs/utility/basedefinitions.h"

/*! \brief Swap endianness of src if needed */
gmx_int32_t gmx_htonl(gmx_int32_t src);

/*! \brief And the reverse */
gmx_int32_t gmx_ntohl(gmx_int32_t src);

#endif
