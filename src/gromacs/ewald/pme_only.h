/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief This file contains function declarations necessary for
 * running on an MPI rank doing only PME long-ranged work.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_ONLY_H
#define GMX_EWALD_PME_ONLY_H

#include <string>

#include "gromacs/timing/walltime_accounting.h"

struct t_commrec;
struct t_inputrec;
struct t_nrnb;
struct gmx_pme_t;
struct gmx_wallcycle;

enum class PmeRunMode;
namespace gmx
{
class DeviceStreamManager;
}

/*! \brief Called on the nodes that do PME exclusively */
int gmx_pmeonly(gmx_pme_t*                      pme,
                const t_commrec*                cr,
                t_nrnb*                         mynrnb,
                gmx_wallcycle*                  wcycle,
                gmx_walltime_accounting_t       walltime_accounting,
                t_inputrec*                     ir,
                PmeRunMode                      runMode,
                const gmx::DeviceStreamManager* deviceStreamManager);

#endif
