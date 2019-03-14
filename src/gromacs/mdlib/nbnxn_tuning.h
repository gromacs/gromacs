/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief Declares functions for tuning adjustable parameters for the nbnxn non-bonded search and interaction kernels
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \inlibraryapi
 * \ingroup __module_nb_verlet
 */

#ifndef NBNXN_TUNING_H
#define NBNXN_TUNING_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"

namespace gmx
{
class CpuInfo;
class MDLogger;
}

struct gmx_mtop_t;
struct interaction_const_t;
struct NbnxnListParameters;
struct t_commrec;
struct t_inputrec;

/*! \brief Try to increase nstlist when using the Verlet cut-off scheme
 *
 * \param[in,out] fplog    Log file
 * \param[in]     cr       The communication record
 * \param[in]     ir       The input parameter record
 * \param[in]     nstlistOnCmdline  The value of nstlist provided on the command line
 * \param[in]     mtop     The global topology
 * \param[in]     box      The unit cell
 * \param[in]     useOrEmulateGpuForNonbondeds  Tells if we are using a GPU for non-bondeds
 * \param[in]     cpuinfo  Information about the CPU(s)
 */
void increaseNstlist(FILE *fplog, t_commrec *cr,
                     t_inputrec *ir, int nstlistOnCmdline,
                     const gmx_mtop_t *mtop,
                     const matrix box,
                     bool useOrEmulateGpuForNonbondeds,
                     const gmx::CpuInfo &cpuinfo);

/*! \brief Set up the dynamic pairlist pruning
 *
 * \param[in,out] mdlog            MD logger
 * \param[in]     ir               The input parameter record
 * \param[in]     mtop             The global topology
 * \param[in]     box              The unit cell
 * \param[in]     nbnxnKernelType  The type of nbnxn kernel used
 * \param[in]     ic               The nonbonded interactions constants
 * \param[in,out] listParams       The list setup parameters
 */
void setupDynamicPairlistPruning(const gmx::MDLogger       &mdlog,
                                 const t_inputrec          *ir,
                                 const gmx_mtop_t          *mtop,
                                 matrix                     box,
                                 int                        nbnxnKernelType,
                                 const interaction_const_t *ic,
                                 NbnxnListParameters       *listParams);

#endif /* NBNXN_TUNING_H */
