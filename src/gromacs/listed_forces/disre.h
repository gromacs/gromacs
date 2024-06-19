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
/*! \libinternal \file
 * \brief
 * Declares functions for handling distance restraints.
 *
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_DISRE_H
#define GMX_LISTED_FORCES_DISRE_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_multisim_t;
class history_t;
struct t_commrec;
struct t_disresdata;
struct t_oriresdata;
struct t_fcdata;
struct t_inputrec;
struct t_pbc;
class t_state;
enum class DDRole;
enum class NumRanks;
union t_iparams;

namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx

//! Whether distance restraints are called from mdrun or from an analysis tool
enum class DisResRunMode
{
    MDRun,
    AnalysisTool
};

/*! \brief
 * Initiates *disresdata.
 *
 * Must be called once, nbonds is the number
 * of iatoms in the ilist of the idef struct.
 * When time averaging is used, the history is initialized in state,
 * unless it was read before from a checkpoint file.
 * The implementation of distance restraints with -multidir
 * must differ according to whether REMD is active.
 */
void init_disres(FILE*                 fplog,
                 const gmx_mtop_t&     mtop,
                 t_inputrec*           ir,
                 DisResRunMode         disResRunMode,
                 DDRole                ddRole,
                 NumRanks              numRanks,
                 MPI_Comm              communicator,
                 const gmx_multisim_t* ms,
                 t_disresdata*         disresdata,
                 t_state*              state,
                 gmx_bool              bIsREMD);

/*! \brief
 * Calculates r and r^-3 (inst. and time averaged) for all pairs
 * and the ensemble averaged r^-6 (inst. and time averaged) for all restraints
 */
void calc_disres_R_6(const t_commrec*      cr,
                     const gmx_multisim_t* ms,
                     int                   nfa,
                     const t_iatom*        fa,
                     const rvec*           x,
                     const t_pbc*          pbc,
                     t_disresdata*         disresdata,
                     const history_t*      hist);

//! Calculates the distance restraint forces, return the potential.
real ta_disres(int                       nfa,
               const t_iatom*            forceatoms,
               const t_iparams*          ip,
               const rvec*               x,
               rvec4*                    f,
               rvec*                     fshift,
               const t_pbc*              pbc,
               real                      lambda,
               real*                     dvdlambda,
               gmx::ArrayRef<const real> charge,
               t_fcdata gmx_unused* fcd,
               t_disresdata*        disresdata,
               t_oriresdata gmx_unused* oriresdata,
               int*                     global_atom_index);

//! Copies the new time averages that have been calculated in calc_disres_R_6.
void update_disres_history(const t_disresdata& disresdata, history_t* hist);

#endif
