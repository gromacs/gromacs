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

#include "mdoutf.h"

#include "config.h"

#include <cstdlib>
#include <cstring>

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/edsamhistory.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/swaphistory.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

struct gmx_mdoutf
{
    t_fileio*                      fp_trn;
    t_fileio*                      fp_xtc;
    gmx_tng_trajectory_t           tng;
    gmx_tng_trajectory_t           tng_low_prec;
    int                            x_compression_precision; /* only used by XTC output */
    ener_file_t                    fp_ene;
    const char*                    fn_cpt;
    gmx_bool                       bKeepAndNumCPT;
    IntegrationAlgorithm           eIntegrator;
    gmx_bool                       bExpanded;
    LambdaWeightCalculation        elamstats;
    int                            simulation_part;
    FILE*                          fp_dhdl;
    int                            natoms_global;
    int                            natoms_x_compressed;
    const SimulationGroups*        groups; /* for compressed position writing */
    gmx_wallcycle*                 wcycle;
    rvec*                          f_global;
    gmx::IMDOutputProvider*        outputProvider;
    const gmx::MDModulesNotifiers* mdModulesNotifiers;
    bool                           simulationsShareState;
    MPI_Comm                       mainRanksComm;
};


gmx_mdoutf_t init_mdoutf(FILE*                          fplog,
                         int                            nfile,
                         const t_filenm                 fnm[],
                         const gmx::MdrunOptions&       mdrunOptions,
                         const t_commrec*               cr,
                         gmx::IMDOutputProvider*        outputProvider,
                         const gmx::MDModulesNotifiers& mdModulesNotifiers,
                         const t_inputrec*              ir,
                         const gmx_mtop_t&              top_global,
                         const gmx_output_env_t*        oenv,
                         gmx_wallcycle*                 wcycle,
                         const gmx::StartingBehavior    startingBehavior,
                         bool                           simulationsShareState,
                         const gmx_multisim_t*          ms)
{
    gmx_mdoutf_t of;
    const char * appendMode = "a+", *writeMode = "w+", *filemode;
    gmx_bool     bCiteTng = FALSE;
    int          i;
    bool restartWithAppending = (startingBehavior == gmx::StartingBehavior::RestartWithAppending);

    snew(of, 1);

    of->fp_trn       = nullptr;
    of->fp_ene       = nullptr;
    of->fp_xtc       = nullptr;
    of->tng          = nullptr;
    of->tng_low_prec = nullptr;
    of->fp_dhdl      = nullptr;

    of->eIntegrator             = ir->eI;
    of->bExpanded               = ir->bExpanded;
    of->elamstats               = ir->expandedvals->elamstats;
    of->simulation_part         = ir->simulation_part;
    of->x_compression_precision = static_cast<int>(ir->x_compression_precision);
    of->wcycle                  = wcycle;
    of->f_global                = nullptr;
    of->outputProvider          = outputProvider;

    GMX_RELEASE_ASSERT(!simulationsShareState || ms != nullptr,
                       "Need valid multisim object when simulations share state");
    of->simulationsShareState = simulationsShareState;
    if (of->simulationsShareState)
    {
        of->mainRanksComm = ms->mainRanksComm_;
    }

    if (MAIN(cr))
    {
        of->bKeepAndNumCPT = mdrunOptions.checkpointOptions.keepAndNumberCheckpointFiles;

        filemode = restartWithAppending ? appendMode : writeMode;

        if (EI_DYNAMICS(ir->eI) && ir->nstxout_compressed > 0)
        {
            const char* filename;
            filename = ftp2fn(efCOMPRESSED, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efXTC: of->fp_xtc = open_xtc(filename, filemode); break;
                case efTNG:
                    gmx_tng_open(filename, filemode[0], &of->tng_low_prec);
                    if (filemode[0] == 'w')
                    {
                        gmx_tng_prepare_low_prec_writing(of->tng_low_prec, &top_global, ir);
                    }
                    bCiteTng = TRUE;
                    break;
                default: gmx_incons("Invalid reduced precision file format");
            }
        }
        if ((EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
            && (!GMX_FAHCORE
                && !(EI_DYNAMICS(ir->eI) && ir->nstxout == 0 && ir->nstvout == 0 && ir->nstfout == 0)))
        {
            const char* filename;
            filename = ftp2fn(efTRN, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efTRR:
                case efTRN:
                    /* If there is no uncompressed coordinate output and
                       there is compressed TNG output write forces
                       and/or velocities to the TNG file instead. */
                    if (ir->nstxout != 0 || ir->nstxout_compressed == 0 || !of->tng_low_prec)
                    {
                        of->fp_trn = gmx_trr_open(filename, filemode);
                    }
                    break;
                case efTNG:
                    gmx_tng_open(filename, filemode[0], &of->tng);
                    if (filemode[0] == 'w')
                    {
                        gmx_tng_prepare_md_writing(of->tng, &top_global, ir);
                    }
                    bCiteTng = TRUE;
                    break;
                default: gmx_incons("Invalid full precision file format");
            }
        }
        if (EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
        {
            of->fp_ene = open_enx(ftp2fn(efEDR, nfile, fnm), filemode);
        }
        of->fn_cpt = opt2fn("-cpo", nfile, fnm);

        if ((ir->efep != FreeEnergyPerturbationType::No || ir->bSimTemp) && ir->fepvals->nstdhdl > 0
            && (ir->fepvals->separate_dhdl_file == SeparateDhdlFile::Yes) && EI_DYNAMICS(ir->eI))
        {
            if (restartWithAppending)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-dhdl", nfile, fnm), filemode);
            }
            else
            {
                of->fp_dhdl = open_dhdl(opt2fn("-dhdl", nfile, fnm), ir, oenv);
            }
        }

        outputProvider->initOutput(fplog, nfile, fnm, restartWithAppending, oenv);
        of->mdModulesNotifiers = &mdModulesNotifiers;

        /* Set up atom counts so they can be passed to actual
           trajectory-writing routines later. Also, XTC writing needs
           to know what (and how many) atoms might be in the XTC
           groups, and how to look up later which ones they are. */
        of->natoms_global       = top_global.natoms;
        of->groups              = &top_global.groups;
        of->natoms_x_compressed = 0;
        for (i = 0; (i < top_global.natoms); i++)
        {
            if (getGroupType(*of->groups, SimulationAtomGroupType::CompressedPositionOutput, i) == 0)
            {
                of->natoms_x_compressed++;
            }
        }

        if (ir->nstfout && haveDDAtomOrdering(*cr))
        {
            snew(of->f_global, top_global.natoms);
        }
    }

    if (bCiteTng)
    {
        please_cite(fplog, "Lundborg2014");
    }

    return of;
}

ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of)
{
    return of->fp_ene;
}

FILE* mdoutf_get_fp_dhdl(gmx_mdoutf_t of)
{
    return of->fp_dhdl;
}

gmx_wallcycle* mdoutf_get_wcycle(gmx_mdoutf_t of)
{
    return of->wcycle;
}

static void mpiBarrierBeforeRename(const bool applyMpiBarrierBeforeRename, MPI_Comm mpiBarrierCommunicator)
{
    if (applyMpiBarrierBeforeRename)
    {
#if GMX_MPI
        MPI_Barrier(mpiBarrierCommunicator);
#else
        GMX_RELEASE_ASSERT(false, "Should not request a barrier without MPI");
        GMX_UNUSED_VALUE(mpiBarrierCommunicator);
#endif
    }
}
/*! \brief Write a checkpoint to the filename
 *
 * Appends the _step<step>.cpt with bNumberAndKeep, otherwise moves
 * the previous checkpoint filename with suffix _prev.cpt.
 */
static void write_checkpoint(const char*                     fn,
                             gmx_bool                        bNumberAndKeep,
                             FILE*                           fplog,
                             const t_commrec*                cr,
                             ivec                            domdecCells,
                             int                             nppnodes,
                             IntegrationAlgorithm            eIntegrator,
                             int                             simulation_part,
                             gmx_bool                        bExpanded,
                             LambdaWeightCalculation         elamstats,
                             int64_t                         step,
                             double                          t,
                             t_state*                        state,
                             ObservablesHistory*             observablesHistory,
                             const gmx::MDModulesNotifiers&  mdModulesNotifiers,
                             gmx::WriteCheckpointDataHolder* modularSimulatorCheckpointData,
                             bool                            applyMpiBarrierBeforeRename,
                             MPI_Comm                        mpiBarrierCommunicator)
{
    t_fileio* fp;
    char*     fntemp; /* the temporary checkpoint file name */
    int       npmenodes;
    char      buf[1024], suffix[5 + STEPSTRSIZE], sbuf[STEPSTRSIZE];
    t_fileio* ret;

    if (haveDDAtomOrdering(*cr))
    {
        npmenodes = cr->npmenodes;
    }
    else
    {
        npmenodes = 0;
    }

#if !GMX_NO_RENAME
    /* make the new temporary filename */
    snew(fntemp, std::strlen(fn) + 5 + STEPSTRSIZE);
    std::strcpy(fntemp, fn);
    fntemp[std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
    sprintf(suffix, "_%s%s", "step", gmx_step_str(step, sbuf));
    std::strcat(fntemp, suffix);
    std::strcat(fntemp, fn + std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1);
#else
    /* if we can't rename, we just overwrite the cpt file.
     * dangerous if interrupted.
     */
    snew(fntemp, std::strlen(fn));
    std::strcpy(fntemp, fn);
#endif
    std::string timebuf = gmx_format_current_time();

    if (fplog)
    {
        fprintf(fplog, "Writing checkpoint, step %s at %s\n\n", gmx_step_str(step, buf), timebuf.c_str());
    }

    /* Get offsets for open files */
    auto outputfiles = gmx_fio_get_output_file_positions();

    fp = gmx_fio_open(fntemp, "w");

    /* We can check many more things now (CPU, acceleration, etc), but
     * it is highly unlikely to have two separate builds with exactly
     * the same version, user, time, and build host!
     */

    int nlambda = (state->dfhist ? state->dfhist->nlambda : 0);

    edsamhistory_t* edsamhist = observablesHistory->edsamHistory.get();
    int             nED       = (edsamhist ? edsamhist->nED : 0);

    swaphistory_t* swaphist    = observablesHistory->swapHistory.get();
    SwapType       eSwapCoords = (swaphist ? swaphist->eSwapCoords : SwapType::No);

    CheckpointHeaderContents headerContents = { CheckPointVersion::UnknownVersion0,
                                                { 0 },
                                                { 0 },
                                                { 0 },
                                                { 0 },
                                                GMX_DOUBLE,
                                                { 0 },
                                                { 0 },
                                                eIntegrator,
                                                simulation_part,
                                                step,
                                                t,
                                                nppnodes,
                                                { 0 },
                                                npmenodes,
                                                state->numAtoms(),
                                                state->ngtc,
                                                state->nnhpres,
                                                state->nhchainlength,
                                                nlambda,
                                                state->flags(),
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                nED,
                                                eSwapCoords,
                                                false };
    std::strcpy(headerContents.version, gmx_version());
    std::strcpy(headerContents.fprog, gmx::getProgramContext().fullBinaryPath().string().c_str());
    std::strcpy(headerContents.ftime, timebuf.c_str());
    if (haveDDAtomOrdering(*cr))
    {
        copy_ivec(domdecCells, headerContents.dd_nc);
    }

    write_checkpoint_data(fp,
                          headerContents,
                          bExpanded,
                          elamstats,
                          state,
                          observablesHistory,
                          mdModulesNotifiers,
                          &outputfiles,
                          modularSimulatorCheckpointData);

    /* we really, REALLY, want to make sure to physically write the checkpoint,
       and all the files it depends on, out to disk. Because we've
       opened the checkpoint with gmx_fio_open(), it's in our list
       of open files.  */
    ret = gmx_fio_all_output_fsync();

    if (ret)
    {
        char buf[STRLEN];
        sprintf(buf,
                "Cannot fsync '%s'; maybe you are out of disk space?",
                gmx_fio_getname(ret).string().c_str());

        if (getenv(GMX_IGNORE_FSYNC_FAILURE_ENV) == nullptr)
        {
            gmx_file(buf);
        }
        else
        {
            gmx_warning("%s", buf);
        }
    }

    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    /* we don't move the checkpoint if the user specified they didn't want it,
       or if the fsyncs failed */
#if !GMX_NO_RENAME
    if (!bNumberAndKeep && !ret)
    {
        // Add a barrier before renaming to reduce chance to get out of sync (#2440)
        // Note: Checkpoint might only exist on some ranks, so put barrier before if clause (#3919)
        mpiBarrierBeforeRename(applyMpiBarrierBeforeRename, mpiBarrierCommunicator);
        if (gmx_fexist(fn))
        {
            /* Rename the previous checkpoint file */
            std::strcpy(buf, fn);
            buf[std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
            std::strcat(buf, "_prev");
            std::strcat(buf, fn + std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1);
            if (!GMX_FAHCORE)
            {
                /* we copy here so that if something goes wrong between now and
                 * the rename below, there's always a state.cpt.
                 * If renames are atomic (such as in POSIX systems),
                 * this copying should be unneccesary.
                 */
                if (gmx_file_copy(fn, buf, FALSE) != 0)
                {
                    GMX_THROW(gmx::FileIOError(
                            gmx::formatString("Cannot rename checkpoint file from %s to %s; maybe "
                                              "you are out of disk space?",
                                              fn,
                                              buf)));
                }
            }
            else
            {
                gmx_file_rename(fn, buf);
            }
        }

        /* Rename the checkpoint file from the temporary to the final name */
        mpiBarrierBeforeRename(applyMpiBarrierBeforeRename, mpiBarrierCommunicator);

        try
        {
            gmx_file_rename(fntemp, fn);
        }
        catch (gmx::FileIOError const&)
        {
            // In this case we can be more helpful than the generic message from gmx_file_rename
            GMX_THROW(gmx::FileIOError(
                    "Cannot rename checkpoint file; maybe you are out of disk space?"));
        }
    }
#endif /* GMX_NO_RENAME */

    sfree(fntemp);

#if GMX_FAHCORE
    /* Always FAH checkpoint immediately after a GROMACS checkpoint.
     *
     * Note that it is critical that we save a FAH checkpoint directly
     * after writing a GROMACS checkpoint. If the program dies, either
     * by the machine powering off suddenly or the process being,
     * killed, FAH can recover files that have only appended data by
     * truncating them to the last recorded length. The GROMACS
     * checkpoint does not just append data, it is fully rewritten each
     * time so a crash between moving the new Gromacs checkpoint file in
     * to place and writing a FAH checkpoint is not recoverable. Thus
     * the time between these operations must be kept as short as
     * possible.
     */
    fcCheckpoint();
#endif /* end GMX_FAHCORE block */
}

void mdoutf_write_checkpoint(gmx_mdoutf_t                    of,
                             FILE*                           fplog,
                             const t_commrec*                cr,
                             int64_t                         step,
                             double                          t,
                             t_state*                        state_global,
                             ObservablesHistory*             observablesHistory,
                             gmx::WriteCheckpointDataHolder* modularSimulatorCheckpointData)
{
    fflush_tng(of->tng);
    fflush_tng(of->tng_low_prec);
    /* Write the checkpoint file.
     * When simulations share the state, an MPI barrier is applied before
     * renaming old and new checkpoint files to minimize the risk of
     * checkpoint files getting out of sync.
     */
    gmx::IVec one_ivec = { 1, 1, 1 };
    write_checkpoint(of->fn_cpt,
                     of->bKeepAndNumCPT,
                     fplog,
                     cr,
                     haveDDAtomOrdering(*cr) ? cr->dd->numCells : one_ivec,
                     haveDDAtomOrdering(*cr) ? cr->dd->nnodes : cr->nnodes,
                     of->eIntegrator,
                     of->simulation_part,
                     of->bExpanded,
                     of->elamstats,
                     step,
                     t,
                     state_global,
                     observablesHistory,
                     *(of->mdModulesNotifiers),
                     modularSimulatorCheckpointData,
                     of->simulationsShareState,
                     of->mainRanksComm);
}

void mdoutf_write_to_trajectory_files(FILE*                           fplog,
                                      const t_commrec*                cr,
                                      gmx_mdoutf_t                    of,
                                      int                             mdof_flags,
                                      int                             natoms,
                                      int64_t                         step,
                                      double                          t,
                                      t_state*                        state_local,
                                      t_state*                        state_global,
                                      ObservablesHistory*             observablesHistory,
                                      gmx::ArrayRef<const gmx::RVec>  f_local,
                                      gmx::WriteCheckpointDataHolder* modularSimulatorCheckpointData)
{
    const rvec* f_global;

    if (haveDDAtomOrdering(*cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            dd_collect_state(cr->dd, state_local, state_global);
        }
        else
        {
            if (mdof_flags & (MDOF_X | MDOF_X_COMPRESSED))
            {
                auto globalXRef = MAIN(cr) ? state_global->x : gmx::ArrayRef<gmx::RVec>();
                dd_collect_vec(cr->dd,
                               state_local->ddp_count,
                               state_local->ddp_count_cg_gl,
                               state_local->cg_gl,
                               state_local->x,
                               globalXRef);
            }
            if (mdof_flags & MDOF_V)
            {
                auto globalVRef = MAIN(cr) ? state_global->v : gmx::ArrayRef<gmx::RVec>();
                dd_collect_vec(cr->dd,
                               state_local->ddp_count,
                               state_local->ddp_count_cg_gl,
                               state_local->cg_gl,
                               state_local->v,
                               globalVRef);
            }
        }
        f_global = of->f_global;
        if (mdof_flags & MDOF_F)
        {
            auto globalFRef = MAIN(cr) ? gmx::arrayRefFromArray(
                                      reinterpret_cast<gmx::RVec*>(of->f_global), of->natoms_global)
                                       : gmx::ArrayRef<gmx::RVec>();
            dd_collect_vec(cr->dd,
                           state_local->ddp_count,
                           state_local->ddp_count_cg_gl,
                           state_local->cg_gl,
                           f_local,
                           globalFRef);
        }
    }
    else
    {
        /* We have the whole state locally: copy the local state pointer */
        state_global = state_local;

        f_global = as_rvec_array(f_local.data());
    }

    if (MAIN(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            mdoutf_write_checkpoint(
                    of, fplog, cr, step, t, state_global, observablesHistory, modularSimulatorCheckpointData);
        }

        if (mdof_flags & (MDOF_X | MDOF_V | MDOF_F))
        {
            const rvec* x = (mdof_flags & MDOF_X) ? state_global->x.rvec_array() : nullptr;
            const rvec* v = (mdof_flags & MDOF_V) ? state_global->v.rvec_array() : nullptr;
            const rvec* f = (mdof_flags & MDOF_F) ? f_global : nullptr;

            if (of->fp_trn)
            {
                gmx_trr_write_frame(of->fp_trn,
                                    step,
                                    t,
                                    state_local->lambda[FreeEnergyPerturbationCouplingType::Fep],
                                    state_local->box,
                                    natoms,
                                    x,
                                    v,
                                    f);
                if (gmx_fio_flush(of->fp_trn) != 0)
                {
                    gmx_file("Cannot write trajectory; maybe you are out of disk space?");
                }
            }

            /* If a TNG file is open for uncompressed coordinate output also write
               velocities and forces to it. */
            else if (of->tng)
            {
                gmx_fwrite_tng(of->tng,
                               FALSE,
                               step,
                               t,
                               state_local->lambda[FreeEnergyPerturbationCouplingType::Fep],
                               state_local->box,
                               natoms,
                               x,
                               v,
                               f);
            }
            /* If only a TNG file is open for compressed coordinate output (no uncompressed
               coordinate output) also write forces and velocities to it. */
            else if (of->tng_low_prec)
            {
                gmx_fwrite_tng(of->tng_low_prec,
                               FALSE,
                               step,
                               t,
                               state_local->lambda[FreeEnergyPerturbationCouplingType::Fep],
                               state_local->box,
                               natoms,
                               x,
                               v,
                               f);
            }
        }
        if (mdof_flags & MDOF_X_COMPRESSED)
        {
            rvec* xxtc = nullptr;

            if (of->natoms_x_compressed == of->natoms_global)
            {
                /* We are writing the positions of all of the atoms to
                   the compressed output */
                xxtc = state_global->x.rvec_array();
            }
            else
            {
                /* We are writing the positions of only a subset of
                   the atoms to the compressed output, so we have to
                   make a copy of the subset of coordinates. */
                int i, j;

                snew(xxtc, of->natoms_x_compressed);
                auto x = makeArrayRef(state_global->x);
                for (i = 0, j = 0; (i < of->natoms_global); i++)
                {
                    if (getGroupType(*of->groups, SimulationAtomGroupType::CompressedPositionOutput, i) == 0)
                    {
                        copy_rvec(x[i], xxtc[j++]);
                    }
                }
            }
            if (write_xtc(of->fp_xtc, of->natoms_x_compressed, step, t, state_local->box, xxtc, of->x_compression_precision)
                == 0)
            {
                gmx_fatal(FARGS,
                          "XTC error. This indicates you are out of disk space, or a "
                          "simulation with major instabilities resulting in coordinates "
                          "that are NaN or too large to be represented in the XTC format.\n");
            }
            gmx_fwrite_tng(of->tng_low_prec,
                           TRUE,
                           step,
                           t,
                           state_local->lambda[FreeEnergyPerturbationCouplingType::Fep],
                           state_local->box,
                           of->natoms_x_compressed,
                           xxtc,
                           nullptr,
                           nullptr);
            if (of->natoms_x_compressed != of->natoms_global)
            {
                sfree(xxtc);
            }
        }
        if (mdof_flags & (MDOF_BOX | MDOF_LAMBDA) && !(mdof_flags & (MDOF_X | MDOF_V | MDOF_F)))
        {
            if (of->tng)
            {
                real  lambda = -1;
                rvec* box    = nullptr;
                if (mdof_flags & MDOF_BOX)
                {
                    box = state_local->box;
                }
                if (mdof_flags & MDOF_LAMBDA)
                {
                    lambda = state_local->lambda[FreeEnergyPerturbationCouplingType::Fep];
                }
                gmx_fwrite_tng(of->tng, FALSE, step, t, lambda, box, natoms, nullptr, nullptr, nullptr);
            }
        }
        if (mdof_flags & (MDOF_BOX_COMPRESSED | MDOF_LAMBDA_COMPRESSED)
            && !(mdof_flags & (MDOF_X_COMPRESSED)))
        {
            if (of->tng_low_prec)
            {
                real  lambda = -1;
                rvec* box    = nullptr;
                if (mdof_flags & MDOF_BOX_COMPRESSED)
                {
                    box = state_local->box;
                }
                if (mdof_flags & MDOF_LAMBDA_COMPRESSED)
                {
                    lambda = state_local->lambda[FreeEnergyPerturbationCouplingType::Fep];
                }
                gmx_fwrite_tng(of->tng_low_prec, FALSE, step, t, lambda, box, natoms, nullptr, nullptr, nullptr);
            }
        }

#if GMX_FAHCORE
        /* Write a FAH checkpoint after writing any other data.  We may end up
           checkpointing twice but it's fast so it's ok. */
        if ((mdof_flags & ~MDOF_CPT))
        {
            fcCheckpoint();
        }
#endif
    }
}

void mdoutf_tng_close(gmx_mdoutf_t of)
{
    if (of->tng || of->tng_low_prec)
    {
        wallcycle_start(of->wcycle, WallCycleCounter::Traj);
        gmx_tng_close(&of->tng);
        gmx_tng_close(&of->tng_low_prec);
        wallcycle_stop(of->wcycle, WallCycleCounter::Traj);
    }
}

void done_mdoutf(gmx_mdoutf_t of)
{
    if (of->fp_ene != nullptr)
    {
        done_ener_file(of->fp_ene);
    }
    if (of->fp_xtc)
    {
        close_xtc(of->fp_xtc);
    }
    if (of->fp_trn)
    {
        gmx_trr_close(of->fp_trn);
    }
    if (of->fp_dhdl != nullptr)
    {
        gmx_fio_fclose(of->fp_dhdl);
    }
    of->outputProvider->finishOutput();
    if (of->f_global != nullptr)
    {
        sfree(of->f_global);
    }

    gmx_tng_close(&of->tng);
    gmx_tng_close(&of->tng_low_prec);

    sfree(of);
}

int mdoutf_get_tng_box_output_interval(gmx_mdoutf_t of)
{
    if (of->tng)
    {
        return gmx_tng_get_box_output_interval(of->tng);
    }
    return 0;
}

int mdoutf_get_tng_lambda_output_interval(gmx_mdoutf_t of)
{
    if (of->tng)
    {
        return gmx_tng_get_lambda_output_interval(of->tng);
    }
    return 0;
}

int mdoutf_get_tng_compressed_box_output_interval(gmx_mdoutf_t of)
{
    if (of->tng_low_prec)
    {
        return gmx_tng_get_box_output_interval(of->tng_low_prec);
    }
    return 0;
}

int mdoutf_get_tng_compressed_lambda_output_interval(gmx_mdoutf_t of)
{
    if (of->tng_low_prec)
    {
        return gmx_tng_get_lambda_output_interval(of->tng_low_prec);
    }
    return 0;
}
