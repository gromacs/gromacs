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

#ifndef GMX_MDLIB_STAT_H
#define GMX_MDLIB_STAT_H

#include <cstdint>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

class gmx_ekindata_t;
struct gmx_enerdata_t;
struct t_vcm;
struct t_inputrec;
struct t_commrec;

namespace gmx
{
template<typename T>
class ArrayRef;
class ObservablesReducer;
} // namespace gmx

typedef struct gmx_global_stat* gmx_global_stat_t;

gmx_global_stat_t global_stat_init(const t_inputrec* ir);

void global_stat_destroy(gmx_global_stat_t gs);

/*! \brief All-reduce energy-like quantities over cr->mpi_comm_mysim  */
void global_stat(const gmx_global_stat&   gs,
                 const t_commrec*         cr,
                 gmx_enerdata_t*          enerd,
                 tensor                   fvir,
                 tensor                   svir,
                 const t_inputrec&        inputrec,
                 gmx_ekindata_t*          ekind,
                 t_vcm*                   vcm,
                 gmx::ArrayRef<real>      sig,
                 bool                     bSumEkinhOld,
                 int                      flags,
                 int64_t                  step,
                 gmx::ObservablesReducer* observablesReducer);

/*! \brief Returns TRUE if io should be done */
inline bool do_per_step(int64_t step, int64_t nstep)
{
    if (nstep != 0)
    {
        return (step % nstep) == 0;
    }
    else
    {
        return false;
    }
}

#endif // GMX_MDLIB_STAT_H
