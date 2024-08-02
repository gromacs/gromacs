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
#include "gmxpre.h"

#include "trajectory_writing.h"

#include <filesystem>
#include <memory>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

void do_md_trajectory_writing(FILE*                          fplog,
                              t_commrec*                     cr,
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
                              const EkindataState            ekindataState)
{
    int   mdof_flags;
    rvec* x_for_confout = nullptr;

    mdof_flags = 0;
    if (do_per_step(step, ir->nstxout))
    {
        mdof_flags |= MDOF_X;
    }
    if (do_per_step(step, ir->nstvout))
    {
        mdof_flags |= MDOF_V;
    }
    if (do_per_step(step, ir->nstfout))
    {
        mdof_flags |= MDOF_F;
    }
    if (do_per_step(step, ir->nstxout_compressed))
    {
        mdof_flags |= MDOF_X_COMPRESSED;
    }
    if (bCPT)
    {
        mdof_flags |= MDOF_CPT;
    }
    if (do_per_step(step, mdoutf_get_tng_box_output_interval(outf)))
    {
        mdof_flags |= MDOF_BOX;
    }
    if (do_per_step(step, mdoutf_get_tng_lambda_output_interval(outf)))
    {
        mdof_flags |= MDOF_LAMBDA;
    }
    if (do_per_step(step, mdoutf_get_tng_compressed_box_output_interval(outf)))
    {
        mdof_flags |= MDOF_BOX_COMPRESSED;
    }
    if (do_per_step(step, mdoutf_get_tng_compressed_lambda_output_interval(outf)))
    {
        mdof_flags |= MDOF_LAMBDA_COMPRESSED;
    }

    if (mdof_flags != 0)
    {
        wallcycle_start(mdoutf_get_wcycle(outf), WallCycleCounter::Traj);
        if (bCPT)
        {
            const bool checkpointEkindata = (ekindataState != EkindataState::NotUsed);
            if (checkpointEkindata)
            {
                update_ekinstate(MAIN(cr) ? &state_global->ekinstate : nullptr,
                                 ekind,
                                 ekindataState == EkindataState::UsedNeedToReduce,
                                 cr);
            }

            if (MAIN(cr))
            {
                state_global->ekinstate.bUpToDate = checkpointEkindata;

                energyOutput.fillEnergyHistory(observablesHistory->energyHistory.get());
            }
        }
        // The current function is only called by legacy code, while
        // mdoutf_write_to_trajectory_files is also called from modular simulator. Create a dummy
        // modular simulator checkpointing object for compatibility.
        gmx::WriteCheckpointDataHolder checkpointDataHolder;
        // Note that part of the following code is duplicated in StatePropagatorData::trajectoryWriterTeardown.
        // This duplication is needed while both legacy and modular code paths are in use.
        // TODO: Remove duplication asap, make sure to keep in sync in the meantime.
        mdoutf_write_to_trajectory_files(
                fplog, cr, outf, mdof_flags, top_global.natoms, step, t, state, state_global, observablesHistory, f, &checkpointDataHolder);
        if (bLastStep && step_rel == ir->nsteps && bDoConfOut && MAIN(cr) && !bRerunMD)
        {
            // With box deformation we would have to correct the output velocities, which is tedious
            const bool makeMoleculesWholeInConfout =
                    (fr->bMolPBC && !ir->bPeriodicMols && !ir_haveBoxDeformation(*ir));

            if (makeMoleculesWholeInConfout && state == state_global)
            {
                /* This (single-rank) run needs to allocate a
                   temporary array of size natoms so that any
                   periodicity removal for mdrun -confout does not
                   perturb the update and thus the final .edr
                   output. This makes .cpt restarts look binary
                   identical, and makes .edr restarts binary
                   identical. */
                snew(x_for_confout, state_global->numAtoms());
                copy_rvecn(state_global->x.rvec_array(), x_for_confout, 0, state_global->numAtoms());
            }
            else
            {
                /* With DD, or no bMolPBC, it doesn't matter if
                   we change state_global->x.rvec_array() */
                x_for_confout = state_global->x.rvec_array();
            }

            /* x and v have been collected in mdoutf_write_to_trajectory_files,
             * because a checkpoint file will always be written
             * at the last step.
             */
            fprintf(stderr, "\nWriting final coordinates.\n");
            if (makeMoleculesWholeInConfout)
            {
                /* Make molecules whole only for confout writing */
                do_pbc_mtop(ir->pbcType, state->box, &top_global, x_for_confout);
            }
            write_sto_conf_mtop(ftp2fn(efSTO, nfile, fnm),
                                *top_global.name,
                                top_global,
                                x_for_confout,
                                state_global->v.rvec_array(),
                                ir->pbcType,
                                state->box);
            if (makeMoleculesWholeInConfout && state == state_global)
            {
                sfree(x_for_confout);
            }
        }
        wallcycle_stop(mdoutf_get_wcycle(outf), WallCycleCounter::Traj);
    }
#if GMX_FAHCORE
    if (MAIN(cr))
    {
        fcWriteVisFrame(ir->pbcType, state_global->box, top_global, state_global->x.rvec_array());
    }
#endif
}
