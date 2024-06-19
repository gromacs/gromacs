/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
#ifndef GMX_MDLIB_MDOUTF_H
#define GMX_MDLIB_MDOUTF_H

#include <cstdint>
#include <cstdio>

#include "gromacs/fileio/enxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

class energyhistory_t;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct ObservablesHistory;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
class t_state;

namespace gmx
{
enum class StartingBehavior;
class IMDOutputProvider;
struct MDModulesNotifiers;
struct MdrunOptions;
class WriteCheckpointDataHolder;
} // namespace gmx

typedef struct gmx_mdoutf* gmx_mdoutf_t;

/*! \brief Allocate and initialize object to manager trajectory writing output
 *
 * Returns a pointer to a data structure with all output file pointers
 * and names required by mdrun.
 */
gmx_mdoutf_t init_mdoutf(FILE*                          fplog,
                         int                            nfile,
                         const t_filenm                 fnm[],
                         const gmx::MdrunOptions&       mdrunOptions,
                         const t_commrec*               cr,
                         gmx::IMDOutputProvider*        outputProvider,
                         const gmx::MDModulesNotifiers& mdModulesNotifiers,
                         const t_inputrec*              ir,
                         const gmx_mtop_t&              mtop,
                         const gmx_output_env_t*        oenv,
                         gmx_wallcycle*                 wcycle,
                         gmx::StartingBehavior          startingBehavior,
                         bool                           simulationsShareState,
                         const gmx_multisim_t*          ms);

/*! \brief Getter for file pointer */
ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of);

/*! \brief Getter for file pointer */
FILE* mdoutf_get_fp_dhdl(gmx_mdoutf_t of);

/*! \brief Getter for wallcycle timer */
gmx_wallcycle* mdoutf_get_wcycle(gmx_mdoutf_t of);

/*! \brief Close TNG files if they are open.
 *
 * This also measures the time it takes to close the TNG
 * files.
 */
void mdoutf_tng_close(gmx_mdoutf_t of);

/*! \brief Close all open output files and free the of pointer */
void done_mdoutf(gmx_mdoutf_t of);

/*! \brief Routine that writes trajectory-like frames.
 *
 * Writes data to trn, xtc and/or checkpoint. What is written is
 * determined by the mdof_flags defined below. Data is collected to
 * the main node only when necessary. Without domain decomposition
 * only data from state_local is used and state_global is ignored.
 *
 * Note that the mdoutf_write_checkpoint function below allows to
 * write only the checkpoint file, and does no data gathering.
 * Except for the modular simulator checkpointing,  the
 * mdoutf_write_to_trajectory_files function is used for
 * _all trajectory / checkpoint writing_.
 *
 * \param[in] fplog                           File handler to log file.
 * \param[in] cr                              Communication record.
 * \param[in] of                              File handler to trajectory file.
 * \param[in] mdof_flags                      Flags indicating what data is written.
 * \param[in] natoms                          The total number of atoms in the system.
 * \param[in] step                            The current time step.
 * \param[in] t                               The current time.
 * \param[in] state_local                     Pointer to the local state object.
 * \param[in] state_global                    Pointer to the global state object.
 * \param[in] observablesHistory              Pointer to the ObservableHistory object.
 * \param[in] f_local                         The local forces.
 * \param[in] modularSimulatorCheckpointData  CheckpointData object used by modular simulator.
 */
void mdoutf_write_to_trajectory_files(FILE*                          fplog,
                                      const t_commrec*               cr,
                                      gmx_mdoutf_t                   of,
                                      int                            mdof_flags,
                                      int                            natoms,
                                      int64_t                        step,
                                      double                         t,
                                      t_state*                       state_local,
                                      t_state*                       state_global,
                                      ObservablesHistory*            observablesHistory,
                                      gmx::ArrayRef<const gmx::RVec> f_local,
                                      gmx::WriteCheckpointDataHolder* modularSimulatorCheckpointData);

/*! \brief Routine allowing to write checkpoint files directly.
 *
 * Unlike mdoutf_write_to_trajectory_files, this will only write a checkpoint file,
 * and it will not gather the state over the different ranks.
 *
 * This should only be called if no other trajectory file writing is needed, and
 * the global state has all required information. Currently, this is only used by
 * the modular checkpointing facility.
 *
 * \param[in] of                              File handler to trajectory file.
 * \param[in] fplog                           File handler to log file.
 * \param[in] cr                              Communication record.
 * \param[in] step                            The current time step.
 * \param[in] t                               The current time.
 * \param[in] state_global                    Pointer to the global state object.
 * \param[in] observablesHistory              Pointer to the ObservableHistory object.
 * \param[in] modularSimulatorCheckpointData  CheckpointData object used by modular simulator.
 */
void mdoutf_write_checkpoint(gmx_mdoutf_t                    of,
                             FILE*                           fplog,
                             const t_commrec*                cr,
                             int64_t                         step,
                             double                          t,
                             t_state*                        state_global,
                             ObservablesHistory*             observablesHistory,
                             gmx::WriteCheckpointDataHolder* modularSimulatorCheckpointData);

/*! \brief Get the output interval of box size of uncompressed TNG output.
 * Returns 0 if no uncompressed TNG file is open.
 */
int mdoutf_get_tng_box_output_interval(gmx_mdoutf_t of);

/*! \brief Get the output interval of lambda of uncompressed TNG output.
 * Returns 0 if no uncompressed TNG file is open.
 */
int mdoutf_get_tng_lambda_output_interval(gmx_mdoutf_t of);

/*! \brief Get the output interval of box size of compressed TNG output.
 * Returns 0 if no compressed TNG file is open.
 */
int mdoutf_get_tng_compressed_box_output_interval(gmx_mdoutf_t of);

/*! \brief Get the output interval of lambda of compressed TNG output.
 * Returns 0 if no compressed TNG file is open.
 */
int mdoutf_get_tng_compressed_lambda_output_interval(gmx_mdoutf_t of);

#define MDOF_X (1u << 0u)
#define MDOF_V (1u << 1u)
#define MDOF_F (1u << 2u)
#define MDOF_X_COMPRESSED (1u << 3u)
#define MDOF_CPT (1u << 4u)
#define MDOF_IMD (1u << 5u)
#define MDOF_BOX (1u << 6u)
#define MDOF_LAMBDA (1u << 7u)
#define MDOF_BOX_COMPRESSED (1u << 8u)
#define MDOF_LAMBDA_COMPRESSED (1u << 9u)

#endif
