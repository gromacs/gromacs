/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#include "gmxpre.h"

#include "md_support.h"

#include <climits>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/md_logging.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

// TODO move this to multi-sim module
bool multisim_int_all_are_equal(const gmx_multisim_t *ms,
                                gmx_int64_t           value)
{
    bool         allValuesAreEqual = true;
    gmx_int64_t *buf;

    GMX_RELEASE_ASSERT(ms, "Invalid use of multi-simulation pointer");

    snew(buf, ms->nsim);
    /* send our value to all other master ranks, receive all of theirs */
    buf[ms->sim] = value;
    gmx_sumli_sim(ms->nsim, buf, ms);

    for (int s = 0; s < ms->nsim; s++)
    {
        if (buf[s] != value)
        {
            allValuesAreEqual = false;
            break;
        }
    }

    sfree(buf);

    return allValuesAreEqual;
}

int multisim_min(const gmx_multisim_t *ms, int nmin, int n)
{
    int     *buf;
    gmx_bool bPos, bEqual;
    int      s, d;

    snew(buf, ms->nsim);
    buf[ms->sim] = n;
    gmx_sumi_sim(ms->nsim, buf, ms);
    bPos   = TRUE;
    bEqual = TRUE;
    for (s = 0; s < ms->nsim; s++)
    {
        bPos   = bPos   && (buf[s] > 0);
        bEqual = bEqual && (buf[s] == buf[0]);
    }
    if (bPos)
    {
        if (bEqual)
        {
            nmin = std::min(nmin, buf[0]);
        }
        else
        {
            /* Find the least common multiple */
            for (d = 2; d < nmin; d++)
            {
                s = 0;
                while (s < ms->nsim && d % buf[s] == 0)
                {
                    s++;
                }
                if (s == ms->nsim)
                {
                    /* We found the LCM and it is less than nmin */
                    nmin = d;
                    break;
                }
            }
        }
    }
    sfree(buf);

    return nmin;
}

void copy_coupling_state(t_state *statea, t_state *stateb,
                         gmx_ekindata_t *ekinda, gmx_ekindata_t *ekindb, t_grpopts* opts)
{

    /* MRS note -- might be able to get rid of some of the arguments.  Look over it when it's all debugged */

    int i, j, nc;

    /* Make sure we have enough space for x and v */
    if (statea->nalloc > stateb->nalloc)
    {
        stateb->nalloc = statea->nalloc;
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        srenew(stateb->x, stateb->nalloc + 1);
        srenew(stateb->v, stateb->nalloc + 1);
    }

    stateb->natoms     = statea->natoms;
    stateb->ngtc       = statea->ngtc;
    stateb->nnhpres    = statea->nnhpres;
    stateb->veta       = statea->veta;
    if (ekinda)
    {
        copy_mat(ekinda->ekin, ekindb->ekin);
        for (i = 0; i < stateb->ngtc; i++)
        {
            ekindb->tcstat[i].T  = ekinda->tcstat[i].T;
            ekindb->tcstat[i].Th = ekinda->tcstat[i].Th;
            copy_mat(ekinda->tcstat[i].ekinh, ekindb->tcstat[i].ekinh);
            copy_mat(ekinda->tcstat[i].ekinf, ekindb->tcstat[i].ekinf);
            ekindb->tcstat[i].ekinscalef_nhc =  ekinda->tcstat[i].ekinscalef_nhc;
            ekindb->tcstat[i].ekinscaleh_nhc =  ekinda->tcstat[i].ekinscaleh_nhc;
            ekindb->tcstat[i].vscale_nhc     =  ekinda->tcstat[i].vscale_nhc;
        }
    }
    copy_rvecn(statea->x, stateb->x, 0, stateb->natoms);
    copy_rvecn(statea->v, stateb->v, 0, stateb->natoms);
    copy_mat(statea->box, stateb->box);
    copy_mat(statea->box_rel, stateb->box_rel);
    copy_mat(statea->boxv, stateb->boxv);

    for (i = 0; i < stateb->ngtc; i++)
    {
        nc = i*opts->nhchainlength;
        for (j = 0; j < opts->nhchainlength; j++)
        {
            stateb->nosehoover_xi[nc+j]  = statea->nosehoover_xi[nc+j];
            stateb->nosehoover_vxi[nc+j] = statea->nosehoover_vxi[nc+j];
        }
    }
    if (stateb->nhpres_xi != NULL)
    {
        for (i = 0; i < stateb->nnhpres; i++)
        {
            nc = i*opts->nhchainlength;
            for (j = 0; j < opts->nhchainlength; j++)
            {
                stateb->nhpres_xi[nc+j]  = statea->nhpres_xi[nc+j];
                stateb->nhpres_vxi[nc+j] = statea->nhpres_vxi[nc+j];
            }
        }
    }
}

real compute_conserved_from_auxiliary(t_inputrec *ir, t_state *state, t_extmass *MassQ)
{
    real quantity = 0;
    switch (ir->etc)
    {
        case etcNO:
            break;
        case etcBERENDSEN:
            break;
        case etcNOSEHOOVER:
            quantity = NPT_energy(ir, state, MassQ);
            break;
        case etcVRESCALE:
            quantity = vrescale_energy(&(ir->opts), state->therm_integral);
            break;
        default:
            break;
    }
    return quantity;
}

/* TODO Specialize this routine into init-time and loop-time versions?
   e.g. bReadEkin is only true when restoring from checkpoint */
void compute_globals(FILE *fplog, gmx_global_stat *gstat, t_commrec *cr, t_inputrec *ir,
                     t_forcerec *fr, gmx_ekindata_t *ekind,
                     t_state *state, t_mdatoms *mdatoms,
                     t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
                     gmx_enerdata_t *enerd, tensor force_vir, tensor shake_vir, tensor total_vir,
                     tensor pres, rvec mu_tot, gmx_constr_t constr,
                     gmx::SimulationSignaller *signalCoordinator,
                     matrix box, int *totalNumberOfBondedInteractions,
                     gmx_bool *bSumEkinhOld, int flags)
{
    tensor   corr_vir, corr_pres;
    gmx_bool bEner, bPres, bTemp;
    gmx_bool bStopCM, bGStat,
             bReadEkin, bEkinAveVel, bScaleEkin, bConstrain;
    real     prescorr, enercorr, dvdlcorr, dvdl_ekin;

    /* translate CGLO flags to gmx_booleans */
    bStopCM       = flags & CGLO_STOPCM;
    bGStat        = flags & CGLO_GSTAT;
    bReadEkin     = (flags & CGLO_READEKIN);
    bScaleEkin    = (flags & CGLO_SCALEEKIN);
    bEner         = flags & CGLO_ENERGY;
    bTemp         = flags & CGLO_TEMPERATURE;
    bPres         = (flags & CGLO_PRESSURE);
    bConstrain    = (flags & CGLO_CONSTRAINT);

    /* we calculate a full state kinetic energy either with full-step velocity verlet
       or half step where we need the pressure */

    bEkinAveVel = (ir->eI == eiVV || (ir->eI == eiVVAK && bPres) || bReadEkin);

    /* in initalization, it sums the shake virial in vv, and to
       sums ekinh_old in leapfrog (or if we are calculating ekinh_old) for other reasons */

    /* ########## Kinetic energy  ############## */

    if (bTemp)
    {
        /* Non-equilibrium MD: this is parallellized, but only does communication
         * when there really is NEMD.
         */

        if (PAR(cr) && (ekind->bNEMD))
        {
            accumulate_u(cr, &(ir->opts), ekind);
        }
        if (!bReadEkin)
        {
            calc_ke_part(state, &(ir->opts), mdatoms, ekind, nrnb, bEkinAveVel);
        }
    }

    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM)
    {
        calc_vcm_grp(0, mdatoms->homenr, mdatoms,
                     state->x, state->v, vcm);
    }

    if (bTemp || bStopCM || bPres || bEner || bConstrain)
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
                wallcycle_start(wcycle, ewcMoveE);
                global_stat(gstat, cr, enerd, force_vir, shake_vir, mu_tot,
                            ir, ekind, constr, bStopCM ? vcm : NULL,
                            signalBuffer.size(), signalBuffer.data(),
                            totalNumberOfBondedInteractions,
                            *bSumEkinhOld, flags);
                wallcycle_stop(wcycle, ewcMoveE);
            }
            signalCoordinator->finalizeSignals();
            *bSumEkinhOld = FALSE;
        }
    }

    if (!ekind->bNEMD && debug && bTemp && (vcm->nr > 0))
    {
        correct_ekin(debug,
                     0, mdatoms->homenr,
                     state->v, vcm->group_p[0],
                     mdatoms->massT, mdatoms->tmass, ekind->ekin);
    }

    /* Do center of mass motion removal */
    if (bStopCM)
    {
        check_cm_grp(fplog, vcm, ir, 1);
        do_stopcm_grp(0, mdatoms->homenr, mdatoms->cVCM,
                      state->x, state->v, vcm);
        inc_nrnb(nrnb, eNR_STOPCM, mdatoms->homenr);
    }

    if (bEner)
    {
        /* Calculate the amplitude of the cosine velocity profile */
        ekind->cosacc.vcos = ekind->cosacc.mvcos/mdatoms->tmass;
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
        enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, &dvdl_ekin,
                                       bEkinAveVel, bScaleEkin);
        enerd->dvdl_lin[efptMASS] = (double) dvdl_ekin;

        enerd->term[F_EKIN] = trace(ekind->ekin);
    }

    /* ##########  Long range energy information ###### */

    if (bEner || bPres || bConstrain)
    {
        calc_dispcorr(ir, fr, box, state->lambda[efptVDW],
                      corr_pres, corr_vir, &prescorr, &enercorr, &dvdlcorr);
    }

    if (bEner)
    {
        enerd->term[F_DISPCORR]  = enercorr;
        enerd->term[F_EPOT]     += enercorr;
        enerd->term[F_DVDL_VDW] += dvdlcorr;
    }

    /* ########## Now pressure ############## */
    if (bPres || bConstrain)
    {

        m_add(force_vir, shake_vir, total_vir);

        /* Calculate pressure and apply LR correction if PPPM is used.
         * Use the box from last timestep since we already called update().
         */

        enerd->term[F_PRES] = calc_pres(fr->ePBC, ir->nwall, box, ekind->ekin, total_vir, pres);

        /* Calculate long range corrections to pressure and energy */
        /* this adds to enerd->term[F_PRES] and enerd->term[F_ETOT],
           and computes enerd->term[F_DISPCORR].  Also modifies the
           total_vir and pres tesors */

        m_add(total_vir, corr_vir, total_vir);
        m_add(pres, corr_pres, pres);
        enerd->term[F_PDISPCORR] = prescorr;
        enerd->term[F_PRES]     += prescorr;
    }
}

void check_nst_param(FILE *fplog, t_commrec *cr,
                     const char *desc_nst, int nst,
                     const char *desc_p, int *p)
{
    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p)/nst + 1)*nst;
        md_print_warn(cr, fplog,
                      "NOTE: %s changes %s to %d\n", desc_nst, desc_p, *p);
    }
}

void set_current_lambdas(gmx_int64_t step, t_lambda *fepvals, gmx_bool bRerunMD,
                         t_trxframe *rerun_fr, t_state *state_global, t_state *state, double lam0[])
/* find the current lambdas.  If rerunning, we either read in a state, or a lambda value,
   requiring different logic. */
{
    real frac;
    int  i, fep_state = 0;
    if (bRerunMD)
    {
        if (rerun_fr->bLambda)
        {
            if (fepvals->delta_lambda == 0)
            {
                state_global->lambda[efptFEP] = rerun_fr->lambda;
                for (i = 0; i < efptNR; i++)
                {
                    if (i != efptFEP)
                    {
                        state->lambda[i] = state_global->lambda[i];
                    }
                }
            }
            else
            {
                /* find out between which two value of lambda we should be */
                frac      = (step*fepvals->delta_lambda);
                fep_state = static_cast<int>(floor(frac*fepvals->n_lambda));
                /* interpolate between this state and the next */
                /* this assumes that the initial lambda corresponds to lambda==0, which is verified in grompp */
                frac = (frac*fepvals->n_lambda)-fep_state;
                for (i = 0; i < efptNR; i++)
                {
                    state_global->lambda[i] = lam0[i] + (fepvals->all_lambda[i][fep_state]) +
                        frac*(fepvals->all_lambda[i][fep_state+1]-fepvals->all_lambda[i][fep_state]);
                }
            }
        }
        else if (rerun_fr->bFepState)
        {
            state_global->fep_state = rerun_fr->fep_state;
            for (i = 0; i < efptNR; i++)
            {
                state_global->lambda[i] = fepvals->all_lambda[i][fep_state];
            }
        }
    }
    else
    {
        if (fepvals->delta_lambda != 0)
        {
            /* find out between which two value of lambda we should be */
            frac = (step*fepvals->delta_lambda);
            if (fepvals->n_lambda > 0)
            {
                fep_state = static_cast<int>(floor(frac*fepvals->n_lambda));
                /* interpolate between this state and the next */
                /* this assumes that the initial lambda corresponds to lambda==0, which is verified in grompp */
                frac = (frac*fepvals->n_lambda)-fep_state;
                for (i = 0; i < efptNR; i++)
                {
                    state_global->lambda[i] = lam0[i] + (fepvals->all_lambda[i][fep_state]) +
                        frac*(fepvals->all_lambda[i][fep_state+1]-fepvals->all_lambda[i][fep_state]);
                }
            }
            else
            {
                for (i = 0; i < efptNR; i++)
                {
                    state_global->lambda[i] = lam0[i] + frac;
                }
            }
        }
        else
        {
            /* if < 0, fep_state was never defined, and we should not set lambda from the state */
            if (state_global->fep_state > -1)
            {
                state_global->fep_state = state->fep_state; /* state->fep_state is the one updated by bExpanded */
                for (i = 0; i < efptNR; i++)
                {
                    state_global->lambda[i] = fepvals->all_lambda[i][state_global->fep_state];
                }
            }
        }
    }
    for (i = 0; i < efptNR; i++)
    {
        state->lambda[i] = state_global->lambda[i];
    }
}

static void min_zero(int *n, int i)
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

    while (nst > 1 && ((i1 > 0 && i1 % nst != 0)  ||
                       (i2 > 0 && i2 % nst != 0)  ||
                       (i3 > 0 && i3 % nst != 0)  ||
                       (i4 > 0 && i4 % nst != 0)))
    {
        nst--;
    }

    return nst;
}

int check_nstglobalcomm(FILE *fplog, t_commrec *cr,
                        int nstglobalcomm, t_inputrec *ir)
{
    if (!EI_DYNAMICS(ir->eI))
    {
        nstglobalcomm = 1;
    }

    if (nstglobalcomm == -1)
    {
        // Set up the default behaviour
        if (!(ir->nstcalcenergy > 0 ||
              ir->nstlist > 0 ||
              ir->etc != etcNO ||
              ir->epc != epcNO))
        {
            /* The user didn't choose the period for anything
               important, so we just make sure we can send signals and
               write output suitably. */
            nstglobalcomm = 10;
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
                                 ir->etc != etcNO ? ir->nsttcouple : 0,
                                 ir->epc != epcNO ? ir->nstpcouple : 0);
        }
    }
    else
    {
        // Check that the user's choice of mdrun -gcom will work
        if (ir->nstlist > 0 &&
            nstglobalcomm > ir->nstlist && nstglobalcomm % ir->nstlist != 0)
        {
            nstglobalcomm = (nstglobalcomm / ir->nstlist)*ir->nstlist;
            md_print_warn(cr, fplog, "WARNING: nstglobalcomm is larger than nstlist, but not a multiple, setting it to %d\n", nstglobalcomm);
        }
        if (ir->nstcalcenergy > 0)
        {
            check_nst_param(fplog, cr, "-gcom", nstglobalcomm,
                            "nstcalcenergy", &ir->nstcalcenergy);
        }
        if (ir->etc != etcNO && ir->nsttcouple > 0)
        {
            check_nst_param(fplog, cr, "-gcom", nstglobalcomm,
                            "nsttcouple", &ir->nsttcouple);
        }
        if (ir->epc != epcNO && ir->nstpcouple > 0)
        {
            check_nst_param(fplog, cr, "-gcom", nstglobalcomm,
                            "nstpcouple", &ir->nstpcouple);
        }

        check_nst_param(fplog, cr, "-gcom", nstglobalcomm,
                        "nstenergy", &ir->nstenergy);

        check_nst_param(fplog, cr, "-gcom", nstglobalcomm,
                        "nstlog", &ir->nstlog);
    }

    if (ir->comm_mode != ecmNO && ir->nstcomm < nstglobalcomm)
    {
        md_print_warn(cr, fplog, "WARNING: Changing nstcomm from %d to %d\n",
                      ir->nstcomm, nstglobalcomm);
        ir->nstcomm = nstglobalcomm;
    }

    if (fplog)
    {
        fprintf(fplog, "Intra-simulation communication will occur every %d steps.\n", nstglobalcomm);
    }
    return nstglobalcomm;
}

void rerun_parallel_comm(t_commrec *cr, t_trxframe *fr,
                         gmx_bool *bLastStep)
{
    rvec    *xp, *vp;

    if (MASTER(cr) && *bLastStep)
    {
        fr->natoms = -1;
    }
    xp = fr->x;
    vp = fr->v;
    gmx_bcast(sizeof(*fr), fr, cr);
    fr->x = xp;
    fr->v = vp;

    *bLastStep = (fr->natoms < 0);

}

void set_state_entries(t_state *state, const t_inputrec *ir)
{
    /* The entries in the state in the tpx file might not correspond
     * with what is needed, so we correct this here.
     */
    state->flags = 0;
    if (ir->efep != efepNO || ir->bExpanded)
    {
        state->flags |= (1<<estLAMBDA);
        state->flags |= (1<<estFEPSTATE);
    }
    state->flags |= (1<<estX);
    if (state->lambda == NULL)
    {
        snew(state->lambda, efptNR);
    }
    if (state->x == NULL)
    {
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        snew(state->x, state->nalloc + 1);
    }
    if (EI_DYNAMICS(ir->eI))
    {
        state->flags |= (1<<estV);
        if (state->v == NULL)
        {
            snew(state->v, state->nalloc + 1);
        }
    }
    if (ir->eI == eiCG)
    {
        state->flags |= (1<<estCGP);
        if (state->cg_p == NULL)
        {
            /* cg_p is not stored in the tpx file, so we need to allocate it */
            snew(state->cg_p, state->nalloc + 1);
        }
    }

    state->nnhpres = 0;
    if (ir->ePBC != epbcNONE)
    {
        state->flags |= (1<<estBOX);
        if (inputrecPreserveShape(ir))
        {
            state->flags |= (1<<estBOX_REL);
        }
        if ((ir->epc == epcPARRINELLORAHMAN) || (ir->epc == epcMTTK))
        {
            state->flags |= (1<<estBOXV);
        }
        if (ir->epc != epcNO)
        {
            if (inputrecNptTrotter(ir) || (inputrecNphTrotter(ir)))
            {
                state->nnhpres = 1;
                state->flags  |= (1<<estNHPRES_XI);
                state->flags  |= (1<<estNHPRES_VXI);
                state->flags  |= (1<<estSVIR_PREV);
                state->flags  |= (1<<estFVIR_PREV);
                state->flags  |= (1<<estVETA);
                state->flags  |= (1<<estVOL0);
            }
            else
            {
                state->flags |= (1<<estPRES_PREV);
            }
        }
    }

    if (ir->etc == etcNOSEHOOVER)
    {
        state->flags |= (1<<estNH_XI);
        state->flags |= (1<<estNH_VXI);
    }

    if (ir->etc == etcVRESCALE)
    {
        state->flags |= (1<<estTC_INT);
    }

    init_gtc_state(state, state->ngtc, state->nnhpres, ir->opts.nhchainlength); /* allocate the space for nose-hoover chains */
    init_ekinstate(&state->ekinstate, ir);
    snew(state->enerhist, 1);
    init_energyhistory(state->enerhist);
    init_df_history(&state->dfhist, ir->fepvals->n_lambda);
    state->swapstate.eSwapCoords = ir->eSwapCoords;
}
