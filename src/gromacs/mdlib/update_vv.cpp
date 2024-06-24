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
#include "gmxpre.h"

#include "update_vv.h"

#include <cmath>

#include <algorithm>

#include "gromacs/domdec/localtopologychecker.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"

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
                          gmx_wallcycle*                                           wcycle)
{
    if (!bFirstStep || startingBehavior == gmx::StartingBehavior::NewSimulation)
    {
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        rvec* vbuf = nullptr;

        wallcycle_start(wcycle, WallCycleCounter::Update);
        if (ir->eI == IntegrationAlgorithm::VV && bInitStep)
        {
            /* if using velocity verlet with full time step Ekin,
             * take the first half step only to compute the
             * virial for the first step. From there,
             * revert back to the initial coordinates
             * so that the input is actually the initial step.
             */
            snew(vbuf, state->numAtoms());
            copy_rvecn(state->v.rvec_array(), vbuf, 0, state->numAtoms()); /* should make this better for parallelizing? */
        }
        else
        {
            /* this is for NHC in the Ekin(t+dt/2) version of vv */
            trotter_update(ir,
                           step,
                           ekind,
                           state,
                           total_vir,
                           mdatoms->homenr,
                           mdatoms->cTC,
                           mdatoms->invmass,
                           MassQ,
                           trotter_seq,
                           TrotterSequence::One);
        }

        // This is not used when updating under VV
        gmx::Matrix3x3 dummyParrinelloRahmanM;
        upd->update_coords(*ir,
                           step,
                           mdatoms->homenr,
                           mdatoms->havePartiallyFrozenAtoms,
                           mdatoms->ptype,
                           mdatoms->invmass,
                           mdatoms->invMassPerDim,
                           state,
                           f->view().forceWithPadding(),
                           fcdata,
                           ekind,
                           dummyParrinelloRahmanM,
                           etrtVELOCITY1,
                           cr,
                           constr != nullptr);

        wallcycle_stop(wcycle, WallCycleCounter::Update);
        constrain_velocities(constr, do_log, do_ene, step, state, nullptr, bCalcVir, shake_vir);
        wallcycle_start(wcycle, WallCycleCounter::Update);
        /* if VV, compute the pressure and constraints */
        /* For VV2, we strictly only need this if using pressure
         * control, but we really would like to have accurate pressures
         * printed out.
         * Think about ways around this in the future?
         * For now, keep this choice in comments.
         */
        /*bPres = (ir->eI==IntegrationAlgorithm::VV || inputrecNptTrotter(ir)); */
        /*bTemp = ((ir->eI==IntegrationAlgorithm::VV &&(!bInitStep)) || (ir->eI==IntegrationAlgorithm::VVAK && inputrecNptTrotter(ir)));*/
        bool bPres = TRUE;
        bool bTemp = ((ir->eI == IntegrationAlgorithm::VV && (!bInitStep))
                      || (ir->eI == IntegrationAlgorithm::VVAK));
        if (bCalcEner && ir->eI == IntegrationAlgorithm::VVAK)
        {
            *bSumEkinhOld = TRUE;
        }
        /* for vv, the first half of the integration actually corresponds to the previous step.
            So we need information from the last step in the first half of the integration */
        if (bGStat || do_per_step(step - 1, nstglobalcomm))
        {
            wallcycle_stop(wcycle, WallCycleCounter::Update);
            int cglo_flags =
                    ((bGStat ? CGLO_GSTAT : 0) | (bCalcEner ? CGLO_ENERGY : 0)
                     | (bTemp ? CGLO_TEMPERATURE : 0) | (bPres ? CGLO_PRESSURE : 0)
                     | (bPres ? CGLO_CONSTRAINT : 0) | (bStopCM ? CGLO_STOPCM : 0) | CGLO_SCALEEKIN);
            compute_globals(gstat,
                            cr,
                            ir,
                            fr,
                            ekind,
                            makeConstArrayRef(state->x),
                            makeConstArrayRef(state->v),
                            state->box,
                            mdatoms,
                            nrnb,
                            vcm,
                            wcycle,
                            enerd,
                            force_vir,
                            shake_vir,
                            total_vir,
                            pres,
                            nullSignaller,
                            state->box,
                            bSumEkinhOld,
                            cglo_flags,
                            step,
                            observablesReducer);
            /* explanation of above:
                a) We compute Ekin at the full time step
                if 1) we are using the AveVel Ekin, and it's not the
                initial step, or 2) if we are using AveEkin, but need the full
                time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                b) If we are using EkinAveEkin for the kinetic energy for the temperature control, we still feed in
                EkinAveVel because it's needed for the pressure */
            if (bStopCM)
            {
                process_and_stopcm_grp(
                        fplog, vcm, *mdatoms, makeArrayRef(state->x), makeArrayRef(state->v));
                inc_nrnb(nrnb, eNR_STOPCM, mdatoms->homenr);
            }
            wallcycle_start(wcycle, WallCycleCounter::Update);
        }
        /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
        if (!bInitStep)
        {
            if (bTrotter)
            {
                m_add(force_vir, shake_vir, total_vir); /* we need the un-dispersion corrected total vir here */
                trotter_update(ir,
                               step,
                               ekind,
                               state,
                               total_vir,
                               mdatoms->homenr,
                               mdatoms->cTC,
                               mdatoms->invmass,
                               MassQ,
                               trotter_seq,
                               TrotterSequence::Two);

                /* TODO This is only needed when we're about to write
                 * a checkpoint, because we use it after the restart
                 * (in a kludge?). But what should we be doing if
                 * the startingBehavior is NewSimulation or bInitStep are true? */
                if (inputrecNptTrotter(ir) || inputrecNphTrotter(ir))
                {
                    copy_mat(shake_vir, state->svir_prev);
                    copy_mat(force_vir, state->fvir_prev);
                }
                if ((inputrecNptTrotter(ir) || inputrecNvtTrotter(ir)) && ir->eI == IntegrationAlgorithm::VV)
                {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(
                            &(ir->opts), ekind, nullptr, (ir->eI == IntegrationAlgorithm::VV), FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            else if (bExchanged)
            {
                wallcycle_stop(wcycle, WallCycleCounter::Update);
                /* We need the kinetic energy at minus the half step for determining
                 * the full step kinetic energy and possibly for T-coupling.*/
                /* This may not be quite working correctly yet . . . . */
                compute_globals(gstat,
                                cr,
                                ir,
                                fr,
                                ekind,
                                makeConstArrayRef(state->x),
                                makeConstArrayRef(state->v),
                                state->box,
                                mdatoms,
                                nrnb,
                                vcm,
                                wcycle,
                                enerd,
                                nullptr,
                                nullptr,
                                nullptr,
                                nullptr,
                                nullSignaller,
                                state->box,
                                bSumEkinhOld,
                                CGLO_GSTAT | CGLO_TEMPERATURE,
                                step,
                                observablesReducer);
                wallcycle_start(wcycle, WallCycleCounter::Update);
            }
        }
        /* if it's the initial step, we performed this first step just to get the constraint virial */
        if (ir->eI == IntegrationAlgorithm::VV && bInitStep)
        {
            copy_rvecn(vbuf, state->v.rvec_array(), 0, state->numAtoms());
            sfree(vbuf);
        }
        wallcycle_stop(wcycle, WallCycleCounter::Update);
    }

    /* compute the conserved quantity */
    *saved_conserved_quantity = NPT_energy(ir->pressureCouplingOptions,
                                           ir->etc,
                                           gmx::constArrayRefFromArray(ir->opts.nrdf, ir->opts.ngtc),
                                           *ekind,
                                           inputrecNvtTrotter(ir) || inputrecNptTrotter(ir),
                                           state,
                                           MassQ);
    if (ir->eI == IntegrationAlgorithm::VV)
    {
        *last_ekin = enerd->term[F_EKIN];
    }
    if ((ir->eDispCorr != DispersionCorrectionType::EnerPres)
        && (ir->eDispCorr != DispersionCorrectionType::AllEnerPres))
    {
        *saved_conserved_quantity -= enerd->term[F_DISPCORR];
    }
    /* sum up the foreign kinetic energy and dK/dl terms for vv.  currently done every step so that dhdl is correct in the .edr */
    if (ir->efep != FreeEnergyPerturbationType::No)
    {
        accumulateKineticLambdaComponents(enerd, state->lambda, *ir->fepvals);
    }
}

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
                           gmx_wallcycle*                                           wcycle)
{
    // This is not used when updating under VV
    gmx::Matrix3x3 dummyParrinelloRahmanM;
    /* velocity half-step update */
    upd->update_coords(*ir,
                       step,
                       mdatoms->homenr,
                       mdatoms->havePartiallyFrozenAtoms,
                       mdatoms->ptype,
                       mdatoms->invmass,
                       mdatoms->invMassPerDim,
                       state,
                       f->view().forceWithPadding(),
                       fcdata,
                       ekind,
                       dummyParrinelloRahmanM,
                       etrtVELOCITY2,
                       cr,
                       constr != nullptr);


    /* Above, initialize just copies ekinh into ekin,
     * it doesn't copy position (for VV),
     * and entire integrator for MD.
     */

    if (ir->eI == IntegrationAlgorithm::VVAK)
    {
        cbuf->resize(state->x.size());
        std::copy(state->x.begin(), state->x.end(), cbuf->begin());
    }

    if (ir->bPull && ir->pull->bSetPbcRefToPrevStepCOM)
    {
        updatePrevStepPullCom(pull_work, state->pull_com_prev_step);
    }

    upd->update_coords(*ir,
                       step,
                       mdatoms->homenr,
                       mdatoms->havePartiallyFrozenAtoms,
                       mdatoms->ptype,
                       mdatoms->invmass,
                       mdatoms->invMassPerDim,
                       state,
                       f->view().forceWithPadding(),
                       fcdata,
                       ekind,
                       dummyParrinelloRahmanM,
                       etrtPOSITION,
                       cr,
                       constr != nullptr);

    wallcycle_stop(wcycle, WallCycleCounter::Update);

    constrain_coordinates(
            constr, do_log, do_ene, step, state, upd->xp()->arrayRefWithPadding(), dvdl_constr, bCalcVir, shake_vir);

    upd->update_sd_second_half(
            *ir, step, dvdl_constr, mdatoms->homenr, mdatoms->ptype, mdatoms->invmass, state, cr, nrnb, wcycle, constr, do_log, do_ene);
    upd->finish_update(
            *ir, mdatoms->havePartiallyFrozenAtoms, mdatoms->homenr, state, wcycle, constr != nullptr);

    if (ir->eI == IntegrationAlgorithm::VVAK)
    {
        /* erase F_EKIN and F_TEMP here? */
        /* just compute the kinetic energy at the half step to perform a trotter step */
        compute_globals(gstat,
                        cr,
                        ir,
                        fr,
                        ekind,
                        makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v),
                        state->box,
                        mdatoms,
                        nrnb,
                        vcm,
                        wcycle,
                        enerd,
                        force_vir,
                        shake_vir,
                        total_vir,
                        pres,
                        nullSignaller,
                        lastbox,
                        bSumEkinhOld,
                        (bGStat ? CGLO_GSTAT : 0) | CGLO_TEMPERATURE,
                        step,
                        observablesReducer);
        wallcycle_start(wcycle, WallCycleCounter::Update);
        trotter_update(ir,
                       step,
                       ekind,
                       state,
                       total_vir,
                       mdatoms->homenr,
                       mdatoms->cTC,
                       mdatoms->invmass,
                       MassQ,
                       trotter_seq,
                       TrotterSequence::Four);
        /* now we know the scaling, we can compute the positions again */
        std::copy(cbuf->begin(), cbuf->end(), state->x.begin());

        upd->update_coords(*ir,
                           step,
                           mdatoms->homenr,
                           mdatoms->havePartiallyFrozenAtoms,
                           mdatoms->ptype,
                           mdatoms->invmass,
                           mdatoms->invMassPerDim,
                           state,
                           f->view().forceWithPadding(),
                           fcdata,
                           ekind,
                           dummyParrinelloRahmanM,
                           etrtPOSITION,
                           cr,
                           constr != nullptr);
        wallcycle_stop(wcycle, WallCycleCounter::Update);

        /* do we need an extra constraint here? just need to copy out of as_rvec_array(state->v.data()) to upd->xp? */
        /* are the small terms in the shake_vir here due
         * to numerical errors, or are they important
         * physically? I'm thinking they are just errors, but not completely sure.
         * For now, will call without actually constraining, constr=NULL*/
        upd->finish_update(*ir, mdatoms->havePartiallyFrozenAtoms, mdatoms->homenr, state, wcycle, false);
    }
    /* this factor or 2 correction is necessary
        because half of the constraint force is removed
        in the vv step, so we have to double it.  See
        the Issue #1255.  It is not yet clear
        if the factor of 2 is exact, or just a very
        good approximation, and this will be
        investigated.  The next step is to see if this
        can be done adding a dhdl contribution from the
        rattle step, but this is somewhat more
        complicated with the current code. Will be
        investigated, hopefully for 4.6.3. However,
        this current solution is much better than
        having it completely wrong.
        */
    enerd->term[F_DVDL_CONSTR] += 2 * *dvdl_constr;
}
