/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \libinternal \file
 * \brief Declares functions to collect state data to the main rank.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_COLLECT_H
#define GMX_DOMDEC_COLLECT_H

#include "gromacs/math/vectypes.h"

namespace gmx
{
template<typename>
class ArrayRef;
}
struct gmx_domdec_t;
class t_state;

/*! \brief Gathers rvec arrays \p localVector to \p globalVector on the main rank */
void dd_collect_vec(gmx_domdec_t*                  dd,
                    int                            ddpCount,
                    int                            ddpCountCgGl,
                    gmx::ArrayRef<const int>       localCGNumbers,
                    gmx::ArrayRef<const gmx::RVec> localVector,
                    gmx::ArrayRef<gmx::RVec>       globalVector);

/*! \brief Gathers state \p localState to \p globalState on the main rank */
void dd_collect_state(gmx_domdec_t* dd, const t_state* localState, t_state* globalState);

#endif
