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

#include "md_support.h"

#include <climits>
#include <cmath>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

static void calc_ke_part_normal(gmx::ArrayRef<const gmx::RVec> v,
                                const t_grpopts*               opts,
                                const t_mdatoms*               md,
                                gmx_ekindata_t*                ekind,
                                t_nrnb*                        nrnb,
                                gmx_bool                       bEkinAveVel)
{
    int                         g;
    gmx::ArrayRef<t_grp_tcstat> tcstat = ekind->tcstat;

    /* three main: VV with AveVel, vv with AveEkin, leap with AveEkin.  Leap with AveVel is also
       an option, but not supported now.
       bEkinAveVel: If TRUE, we sum into ekin, if FALSE, into ekinh.
     */

    // Now accumulate the partial global and groups ekin.
    for (g = 0; (g < opts->ngtc); g++)
    {
        copy_mat(tcstat[g].ekinh, tcstat[g].ekinh_old);
        if (bEkinAveVel)
        {
            clear_mat(tcstat[g].ekinf);
            tcstat[g].ekinscalef_nhc = 1.0; /* need to clear this -- logic is complicated! */
        }
        else
        {
            clear_mat(tcstat[g].ekinh);
        }
    }
    ekind->dekindl_old = ekind->dekindl;
    int nthread        = gmx_omp_nthreads_get(ModuleMultiThread::Update);

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        // This OpenMP only loops over arrays and does not call any functions
        // or memory allocation. It should not be able to throw, so for now
        // we do not need a try/catch wrapper.
        int     start_t, end_t, n;
        int     gt;
        real    hm;
        int     d, m;
        matrix* ekin_sum;
        real*   dekindl_sum;

        start_t = ((thread + 0) * md->homenr) / nthread;
        end_t   = ((thread + 1) * md->homenr) / nthread;

        ekin_sum    = ekind->ekin_work[thread];
        dekindl_sum = ekind->dekindl_work[thread];

        for (gt = 0; gt < opts->ngtc; gt++)
        {
            clear_mat(ekin_sum[gt]);
        }
        *dekindl_sum = 0.0;

        gt = 0;
        for (n = start_t; n < end_t; n++)
        {
            if (!md->cTC.empty())
            {
                gt = md->cTC[n];
            }
            hm = 0.5 * md->massT[n];

            for (d = 0; (d < DIM); d++)
            {
                for (m = 0; (m < DIM); m++)
                {
                    /* if we're computing a full step velocity, v[d] has v(t).  Otherwise, v(t+dt/2) */
                    ekin_sum[gt][m][d] += hm * v[n][m] * v[n][d];
                }
            }
            if (md->nMassPerturbed && md->bPerturbed[n])
            {
                *dekindl_sum += 0.5 * (md->massB[n] - md->massA[n]) * iprod(v[n], v[n]);
            }
        }
    }

    ekind->dekindl = 0;
    for (int thread = 0; thread < nthread; thread++)
    {
        for (g = 0; g < opts->ngtc; g++)
        {
            if (bEkinAveVel)
            {
                m_add(tcstat[g].ekinf, ekind->ekin_work[thread][g], tcstat[g].ekinf);
            }
            else
            {
                m_add(tcstat[g].ekinh, ekind->ekin_work[thread][g], tcstat[g].ekinh);
            }
        }

        ekind->dekindl += *ekind->dekindl_work[thread];
    }

    inc_nrnb(nrnb, eNR_EKIN, md->homenr);
}

static void calc_ke_part_visc(const matrix                   box,
                              gmx::ArrayRef<const gmx::RVec> x,
                              gmx::ArrayRef<const gmx::RVec> v,
                              const t_grpopts*               opts,
                              const t_mdatoms*               md,
                              gmx_ekindata_t*                ekind,
                              t_nrnb*                        nrnb,
                              gmx_bool                       bEkinAveVel)
{
    int                         start = 0, homenr = md->homenr;
    int                         g, d, n, m, gt = 0;
    rvec                        v_corrt;
    real                        hm;
    gmx::ArrayRef<t_grp_tcstat> tcstat = ekind->tcstat;
    t_cos_acc*                  cosacc = &(ekind->cosacc);
    real                        dekindl;
    real                        fac, cosz;
    double                      mvcos;

    for (g = 0; g < opts->ngtc; g++)
    {
        copy_mat(ekind->tcstat[g].ekinh, ekind->tcstat[g].ekinh_old);
        clear_mat(ekind->tcstat[g].ekinh);
    }
    ekind->dekindl_old = ekind->dekindl;

    fac     = 2 * M_PI / box[ZZ][ZZ];
    mvcos   = 0;
    dekindl = 0;
    for (n = start; n < start + homenr; n++)
    {
        if (!md->cTC.empty())
        {
            gt = md->cTC[n];
        }
        hm = 0.5 * md->massT[n];

        /* Note that the times of x and v differ by half a step */
        /* MRS -- would have to be changed for VV */
        cosz = std::cos(fac * x[n][ZZ]);
        /* Calculate the amplitude of the new velocity profile */
        mvcos += 2 * cosz * md->massT[n] * v[n][XX];

        copy_rvec(v[n], v_corrt);
        /* Subtract the profile for the kinetic energy */
        v_corrt[XX] -= cosz * cosacc->vcos;
        for (d = 0; (d < DIM); d++)
        {
            for (m = 0; (m < DIM); m++)
            {
                /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
                if (bEkinAveVel)
                {
                    tcstat[gt].ekinf[m][d] += hm * v_corrt[m] * v_corrt[d];
                }
                else
                {
                    tcstat[gt].ekinh[m][d] += hm * v_corrt[m] * v_corrt[d];
                }
            }
        }
        if (md->nPerturbed && md->bPerturbed[n])
        {
            /* The minus sign here might be confusing.
             * The kinetic contribution from dH/dl doesn't come from
             * d m(l)/2 v^2 / dl, but rather from d p^2/2m(l) / dl,
             * where p are the momenta. The difference is only a minus sign.
             */
            dekindl -= 0.5 * (md->massB[n] - md->massA[n]) * iprod(v_corrt, v_corrt);
        }
    }
    ekind->dekindl = dekindl;
    cosacc->mvcos  = mvcos;

    inc_nrnb(nrnb, eNR_EKIN, homenr);
}

static void calc_ke_part(gmx::ArrayRef<const gmx::RVec> x,
                         gmx::ArrayRef<const gmx::RVec> v,
                         const matrix                   box,
                         const t_grpopts*               opts,
                         const t_mdatoms*               md,
                         gmx_ekindata_t*                ekind,
                         t_nrnb*                        nrnb,
                         gmx_bool                       bEkinAveVel)
{
    if (ekind->cosacc.cos_accel == 0)
    {
        calc_ke_part_normal(v, opts, md, ekind, nrnb, bEkinAveVel);
    }
    else
    {
        calc_ke_part_visc(box, x, v, opts, md, ekind, nrnb, bEkinAveVel);
    }
}

/* TODO Specialize this routine into init-time and loop-time versions?
   e.g. bReadEkin is only true when restoring from checkpoint */
void compute_globals(gmx_global_stat*               gstat,
                     t_commrec*                     cr,
                     const t_inputrec*              ir,
                     t_forcerec*                    fr,
                     gmx_ekindata_t*                ekind,
                     gmx::ArrayRef<const gmx::RVec> x,
                     gmx::ArrayRef<const gmx::RVec> v,
                     const matrix                   box,
                     const t_mdatoms*               mdatoms,
                     t_nrnb*                        nrnb,
                     t_vcm*                         vcm,
                     gmx_wallcycle*                 wcycle,
                     gmx_enerdata_t*                enerd,
                     tensor                         force_vir,
                     tensor                         shake_vir,
                     tensor                         total_vir,
                     tensor                         pres,
                     gmx::SimulationSignaller*      signalCoordinator,
                     const matrix                   lastbox,
                     gmx_bool*                      bSumEkinhOld,
                     const int                      flags,
                     int64_t                        step,
                     gmx::ObservablesReducer*       observablesReducer)
{
    gmx_bool bEner, bPres, bTemp;
    gmx_bool bStopCM, bGStat, bReadEkin, bEkinAveVel, bScaleEkin, bConstrain;
    real     dvdl_ekin;

    /* translate CGLO flags to gmx_booleans */
    bStopCM    = ((flags & CGLO_STOPCM) != 0);
    bGStat     = ((flags & CGLO_GSTAT) != 0);
    bReadEkin  = ((flags & CGLO_READEKIN) != 0);
    bScaleEkin = ((flags & CGLO_SCALEEKIN) != 0);
    bEner      = ((flags & CGLO_ENERGY) != 0);
    bTemp      = ((flags & CGLO_TEMPERATURE) != 0);
    bPres      = ((flags & CGLO_PRESSURE) != 0);
    bConstrain = ((flags & CGLO_CONSTRAINT) != 0);

    /* we calculate a full state kinetic energy either with full-step velocity verlet
       or half step where we need the pressure */

    bEkinAveVel = (ir->eI == IntegrationAlgorithm::VV
                   || (ir->eI == IntegrationAlgorithm::VVAK && bPres) || bReadEkin);

    /* in initalization, it sums the shake virial in vv, and to
       sums ekinh_old in leapfrog (or if we are calculating ekinh_old) for other reasons */

    /* ########## Kinetic energy  ############## */

    if (bTemp)
    {
        if (!bReadEkin)
        {
            calc_ke_part(x, v, box, &(ir->opts), mdatoms, ekind, nrnb, bEkinAveVel);
        }
    }

    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM)
    {
        calc_vcm_grp(*mdatoms, x, v, vcm);
    }

    if (bTemp || bStopCM || bPres || bEner || bConstrain || observablesReducer->isReductionRequired())
    {
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            *bSumEkinhOld = TRUE;
        }
        else
        {
            gmx::ArrayRef<real> signalBuffer = signalCoordinator->getCommunicationBuffer();
            if (PAR(cr))
            {
                wallcycle_start(wcycle, WallCycleCounter::MoveE);
                global_stat(*gstat,
                            cr,
                            enerd,
                            force_vir,
                            shake_vir,
                            *ir,
                            ekind,
                            bStopCM ? vcm : nullptr,
                            signalBuffer,
                            *bSumEkinhOld,
                            flags,
                            step,
                            observablesReducer);
                wallcycle_stop(wcycle, WallCycleCounter::MoveE);
            }
            signalCoordinator->finalizeSignals();
            *bSumEkinhOld = FALSE;
        }
    }

    if (bEner)
    {
        /* Calculate the amplitude of the cosine velocity profile */
        ekind->cosacc.vcos = ekind->cosacc.mvcos / mdatoms->tmass;
    }

    if (bTemp)
    {
        /* Sum the kinetic energies of the groups & calc temp */
        /* compute full step kinetic energies if vv, or if vv-avek and we are computing the pressure with inputrecNptTrotter */
        /* three maincase:  VV with AveVel (md-vv), vv with AveEkin (md-vv-avek), leap with AveEkin (md).
           Leap with AveVel is not supported; it's not clear that it will actually work.
           bEkinAveVel: If TRUE, we simply multiply ekin by ekinscale to get a full step kinetic energy.
           If FALSE, we average ekinh_old and ekinh*ekinscale_nhc to get an averaged half step kinetic energy.
         */
        enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, &dvdl_ekin, bEkinAveVel, bScaleEkin);
        enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Mass] = static_cast<double>(dvdl_ekin);

        enerd->term[F_EKIN] = trace(ekind->ekin);
    }

    /* ########## Now pressure ############## */
    // TODO: For the VV integrator bConstrain is needed in the conditional. This is confusing, so get rid of this.
    if (bPres || bConstrain)
    {
        m_add(force_vir, shake_vir, total_vir);

        /* Calculate pressure and apply LR correction if PPPM is used.
         * Use the box from last timestep since we already called update().
         */

        enerd->term[F_PRES] = calc_pres(fr->pbcType, ir->nwall, lastbox, ekind->ekin, total_vir, pres);
    }
}

static void min_zero(int* n, int i)
{
    if (i > 0 && (*n == 0 || i < *n))
    {
        *n = i;
    }
}

static int lcd4(int i1, int i2, int i3, int i4)
{
    int nst;

    nst = 0;
    min_zero(&nst, i1);
    min_zero(&nst, i2);
    min_zero(&nst, i3);
    min_zero(&nst, i4);
    if (nst == 0)
    {
        gmx_incons("All 4 inputs for determining nstglobalcomm are <= 0");
    }

    while (nst > 1
           && ((i1 > 0 && i1 % nst != 0) || (i2 > 0 && i2 % nst != 0) || (i3 > 0 && i3 % nst != 0)
               || (i4 > 0 && i4 % nst != 0)))
    {
        nst--;
    }

    return nst;
}

int computeGlobalCommunicationPeriod(const t_inputrec* ir)
{
    int nstglobalcomm = 10;
    {
        // Set up the default behaviour
        if (!(ir->nstcalcenergy > 0 || ir->nstlist > 0 || ir->etc != TemperatureCoupling::No
              || ir->pressureCouplingOptions.epc != PressureCoupling::No))
        {
            /* The user didn't choose the period for anything
               important, so we just make sure we can send signals and
               write output suitably. */
            if (ir->nstenergy > 0 && ir->nstenergy < nstglobalcomm)
            {
                nstglobalcomm = ir->nstenergy;
            }
        }
        else
        {
            /* The user has made a choice (perhaps implicitly), so we
             * ensure that we do timely intra-simulation communication
             * for (possibly) each of the four parts that care.
             *
             * TODO Does the Verlet scheme (+ DD) need any
             * communication at nstlist steps? Is the use of nstlist
             * here a leftover of the twin-range scheme? Can we remove
             * nstlist when we remove the group scheme?
             */
            nstglobalcomm = lcd4(ir->nstcalcenergy,
                                 ir->nstlist,
                                 ir->etc != TemperatureCoupling::No ? ir->nsttcouple : 0,
                                 ir->pressureCouplingOptions.epc != PressureCoupling::No
                                         ? ir->pressureCouplingOptions.nstpcouple
                                         : 0);
        }
    }
    return nstglobalcomm;
}

int computeGlobalCommunicationPeriod(const gmx::MDLogger& mdlog, const t_inputrec* ir, const t_commrec* cr)
{
    const int nstglobalcomm = computeGlobalCommunicationPeriod(ir);

    if (cr->nnodes > 1)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted("Intra-simulation communication will occur every %d steps.\n",
                                     nstglobalcomm);
    }
    return nstglobalcomm;
}

void rerun_parallel_comm(t_commrec* cr, t_trxframe* fr, gmx_bool* bLastStep)
{
    rvec *xp, *vp;

    if (MASTER(cr) && *bLastStep)
    {
        fr->natoms = -1;
    }
    xp = fr->x;
    vp = fr->v;
    gmx_bcast(sizeof(*fr), fr, cr->mpi_comm_mygroup);
    fr->x = xp;
    fr->v = vp;

    *bLastStep = (fr->natoms < 0);
}

// TODO Most of this logic seems to belong in the respective modules
void set_state_entries(t_state* state, const t_inputrec* ir, bool useModularSimulator)
{
    /* The entries in the state in the tpx file might not correspond
     * with what is needed, so we correct this here.
     */
    state->flags = 0;
    if (ir->efep != FreeEnergyPerturbationType::No || ir->bExpanded)
    {
        state->flags |= enumValueToBitMask(StateEntry::Lambda);
        state->flags |= enumValueToBitMask(StateEntry::FepState);
    }
    state->flags |= enumValueToBitMask(StateEntry::X);
    GMX_RELEASE_ASSERT(state->x.size() == state->natoms,
                       "We should start a run with an initialized state->x");
    if (EI_DYNAMICS(ir->eI))
    {
        state->flags |= enumValueToBitMask(StateEntry::V);
    }

    state->nnhpres = 0;
    if (ir->pbcType != PbcType::No)
    {
        state->flags |= enumValueToBitMask(StateEntry::Box);
        if (shouldPreserveBoxShape(ir->pressureCouplingOptions, ir->deform))
        {
            state->flags |= enumValueToBitMask(StateEntry::BoxRel);
        }
        if ((ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman)
            || (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk))
        {
            state->flags |= enumValueToBitMask(StateEntry::BoxV);
            if (!useModularSimulator)
            {
                state->flags |= enumValueToBitMask(StateEntry::PressurePrevious);
            }
        }
        if (inputrecNptTrotter(ir) || (inputrecNphTrotter(ir)))
        {
            state->nnhpres = 1;
            state->flags |= enumValueToBitMask(StateEntry::Nhpresxi);
            state->flags |= enumValueToBitMask(StateEntry::Nhpresvxi);
            state->flags |= enumValueToBitMask(StateEntry::SVirPrev);
            state->flags |= enumValueToBitMask(StateEntry::FVirPrev);
            state->flags |= enumValueToBitMask(StateEntry::Veta);
            state->flags |= enumValueToBitMask(StateEntry::Vol0);
        }
        if (ir->pressureCouplingOptions.epc == PressureCoupling::Berendsen
            || ir->pressureCouplingOptions.epc == PressureCoupling::CRescale)
        {
            state->flags |= enumValueToBitMask(StateEntry::BarosInt);
        }
    }

    if (ir->etc == TemperatureCoupling::NoseHoover)
    {
        state->flags |= enumValueToBitMask(StateEntry::Nhxi);
        state->flags |= enumValueToBitMask(StateEntry::Nhvxi);
    }

    if (ir->etc == TemperatureCoupling::VRescale || ir->etc == TemperatureCoupling::Berendsen)
    {
        state->flags |= enumValueToBitMask(StateEntry::ThermInt);
    }

    init_gtc_state(state, state->ngtc, state->nnhpres, ir->opts.nhchainlength); /* allocate the space for nose-hoover chains */
    init_ekinstate(&state->ekinstate, ir);

    if (ir->bExpanded && !useModularSimulator)
    {
        snew(state->dfhist, 1);
        init_df_history(state->dfhist, ir->fepvals->n_lambda);
    }

    if (ir->pull && ir->pull->bSetPbcRefToPrevStepCOM)
    {
        state->flags |= enumValueToBitMask(StateEntry::PullComPrevStep);
    }
}
