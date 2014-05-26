/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "mdrun.h"
#include "domdec.h"
#include "mtop_util.h"
#include "vcm.h"
#include "nrnb.h"
#include "macros.h"
#include "md_logging.h"
#include "md_support.h"
#include "names.h"

#include "gromacs/timing/wallcycle.h"

/* Is the signal in one simulation independent of other simulations? */
gmx_bool gs_simlocal[eglsNR] = { TRUE, FALSE, FALSE, TRUE };

/* check which of the multisim simulations has the shortest number of
   steps and return that number of nsteps */
gmx_int64_t get_multisim_nsteps(const t_commrec *cr,
                                gmx_int64_t      nsteps)
{
    gmx_int64_t steps_out;

    if (MASTER(cr))
    {
        gmx_int64_t     *buf;
        int              s;

        snew(buf, cr->ms->nsim);

        buf[cr->ms->sim] = nsteps;
        gmx_sumli_sim(cr->ms->nsim, buf, cr->ms);

        steps_out = -1;
        for (s = 0; s < cr->ms->nsim; s++)
        {
            /* find the smallest positive number */
            if (buf[s] >= 0 && ((steps_out < 0) || (buf[s] < steps_out)) )
            {
                steps_out = buf[s];
            }
        }
        sfree(buf);

        /* if we're the limiting simulation, don't do anything */
        if (steps_out >= 0 && steps_out < nsteps)
        {
            char strbuf[255];
            snprintf(strbuf, 255, "Will stop simulation %%d after %s steps (another simulation will end then).\n", "%"GMX_PRId64);
            fprintf(stderr, strbuf, cr->ms->sim, steps_out);
        }
    }
    /* broadcast to non-masters */
    gmx_bcast(sizeof(gmx_int64_t), &steps_out, cr);
    return steps_out;
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
            nmin = min(nmin, buf[0]);
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

int multisim_nstsimsync(const t_commrec *cr,
                        const t_inputrec *ir, int repl_ex_nst)
{
    int nmin;

    if (MASTER(cr))
    {
        nmin = INT_MAX;
        nmin = multisim_min(cr->ms, nmin, ir->nstlist);
        nmin = multisim_min(cr->ms, nmin, ir->nstcalcenergy);
        nmin = multisim_min(cr->ms, nmin, repl_ex_nst);
        if (nmin == INT_MAX)
        {
            gmx_fatal(FARGS, "Can not find an appropriate interval for inter-simulation communication, since nstlist, nstcalcenergy and -replex are all <= 0");
        }
        /* Avoid inter-simulation communication at every (second) step */
        if (nmin <= 2)
        {
            nmin = 10;
        }
    }

    gmx_bcast(sizeof(int), &nmin, cr);

    return nmin;
}

void init_global_signals(globsig_t *gs, const t_commrec *cr,
                         const t_inputrec *ir, int repl_ex_nst)
{
    int i;

    if (MULTISIM(cr))
    {
        gs->nstms = multisim_nstsimsync(cr, ir, repl_ex_nst);
        if (debug)
        {
            fprintf(debug, "Syncing simulations for checkpointing and termination every %d steps\n", gs->nstms);
        }
    }
    else
    {
        gs->nstms = 1;
    }

    for (i = 0; i < eglsNR; i++)
    {
        gs->sig[i] = 0;
        gs->set[i] = 0;
    }
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
        srenew(stateb->x, stateb->nalloc);
        srenew(stateb->v, stateb->nalloc);
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

void compute_globals(FILE *fplog, gmx_global_stat_t gstat, t_commrec *cr, t_inputrec *ir,
                     t_forcerec *fr, gmx_ekindata_t *ekind,
                     t_state *state, t_state *state_global, t_mdatoms *mdatoms,
                     t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
                     gmx_enerdata_t *enerd, tensor force_vir, tensor shake_vir, tensor total_vir,
                     tensor pres, rvec mu_tot, gmx_constr_t constr,
                     globsig_t *gs, gmx_bool bInterSimGS,
                     matrix box, gmx_mtop_t *top_global,
                     gmx_bool *bSumEkinhOld, int flags)
{
    int      i, gsi;
    real     gs_buf[eglsNR];
    tensor   corr_vir, corr_pres;
    gmx_bool bEner, bPres, bTemp, bVV;
    gmx_bool bRerunMD, bStopCM, bGStat, bIterate,
             bFirstIterate, bReadEkin, bEkinAveVel, bScaleEkin, bConstrain;
    real     ekin, temp, prescorr, enercorr, dvdlcorr, dvdl_ekin;

    /* translate CGLO flags to gmx_booleans */
    bRerunMD = flags & CGLO_RERUNMD;
    bStopCM  = flags & CGLO_STOPCM;
    bGStat   = flags & CGLO_GSTAT;

    bReadEkin     = (flags & CGLO_READEKIN);
    bScaleEkin    = (flags & CGLO_SCALEEKIN);
    bEner         = flags & CGLO_ENERGY;
    bTemp         = flags & CGLO_TEMPERATURE;
    bPres         = (flags & CGLO_PRESSURE);
    bConstrain    = (flags & CGLO_CONSTRAINT);
    bIterate      = (flags & CGLO_ITERATE);
    bFirstIterate = (flags & CGLO_FIRSTITERATE);

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
        debug_gmx();
        if (bReadEkin)
        {
            restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
        }
        else
        {

            calc_ke_part(state, &(ir->opts), mdatoms, ekind, nrnb, bEkinAveVel, bIterate);
        }

        debug_gmx();
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
            if (gs != NULL)
            {
                for (i = 0; i < eglsNR; i++)
                {
                    gs_buf[i] = gs->sig[i];
                }
            }
            if (PAR(cr))
            {
                wallcycle_start(wcycle, ewcMoveE);
                global_stat(fplog, gstat, cr, enerd, force_vir, shake_vir, mu_tot,
                            ir, ekind, constr, bStopCM ? vcm : NULL,
                            gs != NULL ? eglsNR : 0, gs_buf,
                            top_global, state,
                            *bSumEkinhOld, flags);
                wallcycle_stop(wcycle, ewcMoveE);
            }
            if (gs != NULL)
            {
                if (MULTISIM(cr) && bInterSimGS)
                {
                    if (MASTER(cr))
                    {
                        /* Communicate the signals between the simulations */
                        gmx_sum_sim(eglsNR, gs_buf, cr->ms);
                    }
                    /* Communicate the signals form the master to the others */
                    gmx_bcast(eglsNR*sizeof(gs_buf[0]), gs_buf, cr);
                }
                for (i = 0; i < eglsNR; i++)
                {
                    if (bInterSimGS || gs_simlocal[i])
                    {
                        /* Set the communicated signal only when it is non-zero,
                         * since signals might not be processed at each MD step.
                         */
                        gsi = (gs_buf[i] >= 0 ?
                               (int)(gs_buf[i] + 0.5) :
                               (int)(gs_buf[i] - 0.5));
                        if (gsi != 0)
                        {
                            gs->set[i] = gsi;
                        }
                        /* Turn off the local signal */
                        gs->sig[i] = 0;
                    }
                }
            }
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
        /* compute full step kinetic energies if vv, or if vv-avek and we are computing the pressure with IR_NPT_TROTTER */
        /* three maincase:  VV with AveVel (md-vv), vv with AveEkin (md-vv-avek), leap with AveEkin (md).
           Leap with AveVel is not supported; it's not clear that it will actually work.
           bEkinAveVel: If TRUE, we simply multiply ekin by ekinscale to get a full step kinetic energy.
           If FALSE, we average ekinh_old and ekinh*ekinscale_nhc to get an averaged half step kinetic energy.
           bSaveEkinOld: If TRUE (in the case of iteration = bIterate is TRUE), we don't reset the ekinscale_nhc.
           If FALSE, we go ahead and erase over it.
         */
        enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, &dvdl_ekin,
                                       bEkinAveVel, bScaleEkin);
        enerd->dvdl_lin[efptMASS] = (double) dvdl_ekin;

        enerd->term[F_EKIN] = trace(ekind->ekin);
    }

    /* ##########  Long range energy information ###### */

    if (bEner || bPres || bConstrain)
    {
        calc_dispcorr(fplog, ir, fr, 0, top_global->natoms, box, state->lambda[efptVDW],
                      corr_pres, corr_vir, &prescorr, &enercorr, &dvdlcorr);
    }

    if (bEner && bFirstIterate)
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
                fep_state = floor(frac*fepvals->n_lambda);
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
                fep_state = floor(frac*fepvals->n_lambda);
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
            if (state->fep_state > 0)
            {
                state_global->fep_state = state->fep_state; /* state->fep is the one updated by bExpanded */
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
        gmx_incons("All 4 inputs for determininig nstglobalcomm are <= 0");
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
        if (!(ir->nstcalcenergy > 0 ||
              ir->nstlist > 0 ||
              ir->etc != etcNO ||
              ir->epc != epcNO))
        {
            nstglobalcomm = 10;
            if (ir->nstenergy > 0 && ir->nstenergy < nstglobalcomm)
            {
                nstglobalcomm = ir->nstenergy;
            }
        }
        else
        {
            /* Ensure that we do timely global communication for
             * (possibly) each of the four following options.
             */
            nstglobalcomm = lcd4(ir->nstcalcenergy,
                                 ir->nstlist,
                                 ir->etc != etcNO ? ir->nsttcouple : 0,
                                 ir->epc != epcNO ? ir->nstpcouple : 0);
        }
    }
    else
    {
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

    return nstglobalcomm;
}

void check_ir_old_tpx_versions(t_commrec *cr, FILE *fplog,
                               t_inputrec *ir, gmx_mtop_t *mtop)
{
    /* Check required for old tpx files */
    if (IR_TWINRANGE(*ir) && ir->nstlist > 1 &&
        ir->nstcalcenergy % ir->nstlist != 0)
    {
        md_print_warn(cr, fplog, "Old tpr file with twin-range settings: modifying energy calculation and/or T/P-coupling frequencies\n");

        if (gmx_mtop_ftype_count(mtop, F_CONSTR) +
            gmx_mtop_ftype_count(mtop, F_CONSTRNC) > 0 &&
            ir->eConstrAlg == econtSHAKE)
        {
            md_print_warn(cr, fplog, "With twin-range cut-off's and SHAKE the virial and pressure are incorrect\n");
            if (ir->epc != epcNO)
            {
                gmx_fatal(FARGS, "Can not do pressure coupling with twin-range cut-off's and SHAKE");
            }
        }
        check_nst_param(fplog, cr, "nstlist", ir->nstlist,
                        "nstcalcenergy", &ir->nstcalcenergy);
        if (ir->epc != epcNO)
        {
            check_nst_param(fplog, cr, "nstlist", ir->nstlist,
                            "nstpcouple", &ir->nstpcouple);
        }
        check_nst_param(fplog, cr, "nstcalcenergy", ir->nstcalcenergy,
                        "nstenergy", &ir->nstenergy);
        check_nst_param(fplog, cr, "nstcalcenergy", ir->nstcalcenergy,
                        "nstlog", &ir->nstlog);
        if (ir->efep != efepNO)
        {
            check_nst_param(fplog, cr, "nstcalcenergy", ir->nstcalcenergy,
                            "nstdhdl", &ir->fepvals->nstdhdl);
        }
    }

    if (EI_VV(ir->eI) && IR_TWINRANGE(*ir) && ir->nstlist > 1)
    {
        gmx_fatal(FARGS, "Twin-range multiple time stepping does not work with integrator %s.", ei_names[ir->eI]);
    }
}

void rerun_parallel_comm(t_commrec *cr, t_trxframe *fr,
                         gmx_bool *bNotLastFrame)
{
    gmx_bool bAlloc;
    rvec    *xp, *vp;

    bAlloc = (fr->natoms == 0);

    if (MASTER(cr) && !*bNotLastFrame)
    {
        fr->natoms = -1;
    }
    xp = fr->x;
    vp = fr->v;
    gmx_bcast(sizeof(*fr), fr, cr);
    fr->x = xp;
    fr->v = vp;

    *bNotLastFrame = (fr->natoms >= 0);

}
