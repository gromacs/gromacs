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
#ifndef GMX_MDLIB_TRAJECTORY_WRITING_H
#define GMX_MDLIB_TRAJECTORY_WRITING_H

#include <cstdint>
#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/utility/basedefinitions.h"

class gmx_ekindata_t;
struct gmx_mtop_t;
struct ObservablesHistory;
struct t_commrec;
struct t_filenm;
struct t_forcerec;
class t_state;
struct t_inputrec;

namespace gmx
{
class EnergyOutput;
template<typename T>
class ArrayRef;
} // namespace gmx

//! The current state of ekindata as passed to do_md_trajectory_writing()
enum class EkindataState
{
    NotUsed,               //!< ekindata is not used this step and should thus not be checkpointed
    UsedNeedToReduce,      //!< ekindata is used this step and terms need to be reduced
    UsedDoNotNeedToReduce, //!< ekindata is used this step and no reduction is needed
};

/*! \brief Wrapper routine for writing trajectories during mdrun
 *
 * This routine does communication (e.g. collecting distributed coordinates).
 *
 * The kinetic energy data \p ekind is only used at steps where energies are
 * calculated or temperature or pressure coupling is done. Thus this data only
 * needs to be written to checkpoint at such steps. It might also contain
 * local contributions that are not yet reduced over ranks. The state of
 * \p ekind is described by \p ekindataState.
 */
void do_md_trajectory_writing(FILE*                          fplog,
                              struct t_commrec*              cr,
                              int                            nfile,
                              const t_filenm                 fnm[],
                              int64_t                        step,
                              int64_t                        step_rel,
                              double                         t,
                              const t_inputrec*              ir,
                              t_state*                       state,
                              t_state*                       state_global,
                              ObservablesHistory*            observablesHistory,
                              const gmx_mtop_t&              top_global,
                              t_forcerec*                    fr,
                              gmx_mdoutf_t                   outf,
                              const gmx::EnergyOutput&       energyOutput,
                              gmx_ekindata_t*                ekind,
                              gmx::ArrayRef<const gmx::RVec> f,
                              gmx_bool                       bCPT,
                              gmx_bool                       bRerunMD,
                              gmx_bool                       bLastStep,
                              gmx_bool                       bDoConfOut,
                              EkindataState                  ekindataState);

#endif
