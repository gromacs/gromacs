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
#ifndef GMX_MDLIB_UPDATE_VV_H
#define GMX_MDLIB_UPDATE_VV_H

#include <cstdint>
#include <cstdio>

#include <array>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_global_stat;
struct gmx_wallcycle;
struct pull_t;
struct t_commrec;
struct t_extmass;
struct t_fcdata;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
class t_state;
struct t_vcm;

namespace gmx
{
class Constraints;
class ForceBuffers;
class ObservablesReducer;
class SimulationSignaller;
class Update;
enum class StartingBehavior : int;
} // namespace gmx

/*! \brief Make the first step of Velocity Verlet integration
 *
 * \param[in]  step              Current timestep.
 * \param[in]  bFirstStep        Is it a first step.
 * \param[in]  bInitStep         Is it an initialization step.
 * \param[in]  startingBehavior  Describes whether this is a restart appending to output files.
 * \param[in]  nstglobalcomm     Will globals be computed on this step.
 * \param[in]  ir                Input record.
 * \param[in]  fr                Force record.
 * \param[in]  cr                Comunication record.
 * \param[in]  state             Simulation state.
 * \param[in]  mdatoms           MD atoms data.
 * \param[in]  fcdata            Force calculation data.
 * \param[in]  MassQ             Mass/pressure data.
 * \param[in]  vcm               Center of mass motion removal.
 * \param[in]  enerd             Energy data.
 * \param[in]  observablesReducer Pointer to the \c ObservablesReducer object
 * \param[in]  ekind             Kinetic energy data.
 * \param[in]  gstat             Storage of thermodynamic parameters data.
 * \param[out] last_ekin         Kinetic energies of the last step.
 * \param[in]  bCalcVir          If the virial is computed on this step.
 * \param[in]  total_vir         Total virial tensor.
 * \param[in]  shake_vir         Constraints virial.
 * \param[in]  force_vir         Force virial.
 * \param[in]  pres              Pressure tensor.
 * \param[in]  do_log            Do logging on this step.
 * \param[in]  do_ene            Print energies on this step.
 * \param[in]  bCalcEner         Compute energies on this step.
 * \param[in]  bGStat            Collect globals this step.
 * \param[in]  bStopCM           Stop the center of mass motion on this step.
 * \param[in]  bTrotter          Do trotter routines this step.
 * \param[in]  bExchanged        If this is a replica exchange step.
 * \param[out] bSumEkinhOld      Old kinetic energies will need to be summed up.
 * \param[out] saved_conserved_quantity  Place to store the conserved energy.
 * \param[in]  f                 Force buffers.
 * \param[in]  upd               Update object.
 * \param[in]  constr            Constraints object.
 * \param[in]  nullSignaller     Simulation signaller.
 * \param[in]  trotter_seq       NPT variables.
 * \param[in]  nrnb              Cycle counters.
 * \param[in]  fplog             Another logger.
 * \param[in]  wcycle            Wall-clock cycle counter.
 */
void integrateVVFirstStep(int64_t                   step,
                          bool                      bFirstStep,
                          bool                      bInitStep,
                          gmx::StartingBehavior     startingBehavior,
                          int                       nstglobalcomm,
                          const t_inputrec*         ir,
                          t_forcerec*               fr,
                          t_commrec*                cr,
                          t_state*                  state,
                          t_mdatoms*                mdatoms,
                          t_fcdata*                 fcdata,
                          t_extmass*                MassQ,
                          t_vcm*                    vcm,
                          gmx_enerdata_t*           enerd,
                          gmx::ObservablesReducer*  observablesReducer,
                          gmx_ekindata_t*           ekind,
                          gmx_global_stat*          gstat,
                          real*                     last_ekin,
                          bool                      bCalcVir,
                          tensor                    total_vir,
                          tensor                    shake_vir,
                          tensor                    force_vir,
                          tensor                    pres,
                          bool                      do_log,
                          bool                      do_ene,
                          bool                      bCalcEner,
                          bool                      bGStat,
                          bool                      bStopCM,
                          bool                      bTrotter,
                          bool                      bExchanged,
                          bool*                     bSumEkinhOld,
                          real*                     saved_conserved_quantity,
                          gmx::ForceBuffers*        f,
                          gmx::Update*              upd,
                          gmx::Constraints*         constr,
                          gmx::SimulationSignaller* nullSignaller,
                          gmx::EnumerationArray<TrotterSequence, std::vector<int>> trotter_seq,
                          t_nrnb*                                                  nrnb,
                          FILE*                                                    fplog,
                          gmx_wallcycle*                                           wcycle);


/*! \brief Make the second step of Velocity Verlet integration
 *
 * \param[in]  step              Current timestep.
 * \param[in]  ir                Input record.
 * \param[in]  fr                Force record.
 * \param[in]  cr                Comunication record.
 * \param[in]  state             Simulation state.
 * \param[in]  mdatoms           MD atoms data.
 * \param[in]  fcdata            Force calculation data.
 * \param[in]  MassQ             Mass/pressure data.
 * \param[in]  vcm               Center of mass motion removal.
 * \param[in]  pull_work         Pulling data.
 * \param[in]  enerd             Energy data.
 * \param[in]  observablesReducer Pointer to the \c ObservablesReducer object
 * \param[in]  ekind             Kinetic energy data.
 * \param[in]  gstat             Storage of thermodynamic parameters data.
 * \param[out] dvdl_constr       FEP data for constraints.
 * \param[in]  bCalcVir          If the virial is computed on this step.
 * \param[in]  total_vir         Total virial tensor.
 * \param[in]  shake_vir         Constraints virial.
 * \param[in]  force_vir         Force virial.
 * \param[in]  pres              Pressure tensor.
 * \param[in]  lastbox           Last recorded PBC box.
 * \param[in]  do_log            Do logging on this step.
 * \param[in]  do_ene            Print energies on this step.
 * \param[in]  bGStat            Collect globals this step.
 * \param[out] bSumEkinhOld      Old kinetic energies need to be summed up.
 * \param[in]  f                 Force buffers.
 * \param[in]  cbuf              Buffer to store intermediate coordinates
 * \param[in]  upd               Update object.
 * \param[in]  constr            Constraints object.
 * \param[in]  nullSignaller     Simulation signaller.
 * \param[in]  trotter_seq       NPT variables.
 * \param[in]  nrnb              Cycle counters.
 * \param[in]  wcycle            Wall-clock cycle counter.
 */
void integrateVVSecondStep(int64_t                   step,
                           const t_inputrec*         ir,
                           t_forcerec*               fr,
                           t_commrec*                cr,
                           t_state*                  state,
                           t_mdatoms*                mdatoms,
                           t_fcdata*                 fcdata,
                           t_extmass*                MassQ,
                           t_vcm*                    vcm,
                           pull_t*                   pull_work,
                           gmx_enerdata_t*           enerd,
                           gmx::ObservablesReducer*  observablesReducer,
                           gmx_ekindata_t*           ekind,
                           gmx_global_stat*          gstat,
                           real*                     dvdl_constr,
                           bool                      bCalcVir,
                           tensor                    total_vir,
                           tensor                    shake_vir,
                           tensor                    force_vir,
                           tensor                    pres,
                           matrix                    lastbox,
                           bool                      do_log,
                           bool                      do_ene,
                           bool                      bGStat,
                           bool*                     bSumEkinhOld,
                           gmx::ForceBuffers*        f,
                           std::vector<gmx::RVec>*   cbuf,
                           gmx::Update*              upd,
                           gmx::Constraints*         constr,
                           gmx::SimulationSignaller* nullSignaller,
                           gmx::EnumerationArray<TrotterSequence, std::vector<int>> trotter_seq,
                           t_nrnb*                                                  nrnb,
                           gmx_wallcycle*                                           wcycle);


#endif // GMX_MDLIB_UPDATE_VV_H
