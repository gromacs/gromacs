/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _isnan() */
#include <float.h>
#endif

#include "typedefs.h"
#include "string2.h"
#include "smalloc.h"
#include "mdrun.h"
#include "domdec.h"
#include "mtop_util.h"
#include "gmx_wallcycle.h"
#include "vcm.h"
#include "nrnb.h"

/* check which of the multisim simulations has the shortest number of
   steps and return that number of nsteps */
gmx_large_int_t get_multisim_nsteps(const t_commrec *cr,
                                    gmx_large_int_t nsteps)
{
    gmx_large_int_t steps_out;

    if MASTER(cr)
    {
        gmx_large_int_t *buf;
        int s;

        snew(buf,cr->ms->nsim);

        buf[cr->ms->sim] = nsteps;
        gmx_sumli_sim(cr->ms->nsim, buf, cr->ms);

        steps_out=-1;
        for(s=0; s<cr->ms->nsim; s++)
        {
            /* find the smallest positive number */
            if (buf[s]>= 0 && ((steps_out < 0) || (buf[s]<steps_out)) )
            {
                steps_out=buf[s];
            }
        }
        sfree(buf);

        /* if we're the limiting simulation, don't do anything */
        if (steps_out>=0 && steps_out<nsteps) 
        {
            char strbuf[255];
            snprintf(strbuf, 255, "Will stop simulation %%d after %s steps (another simulation will end then).\n", gmx_large_int_pfmt);
            fprintf(stderr, strbuf, cr->ms->sim, steps_out);
        }
    }
    /* broadcast to non-masters */
    gmx_bcast(sizeof(gmx_large_int_t), &steps_out, cr);
    return steps_out;
}

int multisim_min(const gmx_multisim_t *ms,int nmin,int n)
{
    int  *buf;
    gmx_bool bPos,bEqual;
    int  s,d;

    snew(buf,ms->nsim);
    buf[ms->sim] = n;
    gmx_sumi_sim(ms->nsim,buf,ms);
    bPos   = TRUE;
    bEqual = TRUE;
    for(s=0; s<ms->nsim; s++)
    {
        bPos   = bPos   && (buf[s] > 0);
        bEqual = bEqual && (buf[s] == buf[0]);
    }
    if (bPos)
    {
        if (bEqual)
        {
            nmin = min(nmin,buf[0]);
        }
        else
        {
            /* Find the least common multiple */
            for(d=2; d<nmin; d++)
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

void
init_signals(FILE *fplog,
             gmx_signal *signal,
             const t_commrec *cr,
             const t_inputrec *ir,
             int nstglobalcomm,
             int nstsignalcomm,
             real max_hours)
{
    int i;

    /* We always communicate intra-simulation signals */
    signal->nstsim = 1;

    if (MULTISIM(cr))
    {
        signal->nstms = check_nstsignalcomm(fplog,cr,nstsignalcomm,ir,max_hours);
    }
    else
    {
        signal->nstms = nstglobalcomm;
    }

    for(i = 0; i < esignalNR; i++)
    {
        signal->init[i] = 0;
        signal->set[i] = 0;
    }
}

void copy_coupling_state(t_state *statea,t_state *stateb, 
                         gmx_ekindata_t *ekinda,gmx_ekindata_t *ekindb, t_grpopts* opts) 
{
    
    /* MRS note -- might be able to get rid of some of the arguments.  Look over it when it's all debugged */
    
    int i,j,nc;

    /* Make sure we have enough space for x and v */
    if (statea->nalloc > stateb->nalloc)
    {
        stateb->nalloc = statea->nalloc;
        srenew(stateb->x,stateb->nalloc);
        srenew(stateb->v,stateb->nalloc);
    }

    stateb->natoms     = statea->natoms;
    stateb->ngtc       = statea->ngtc;
    stateb->nnhpres    = statea->nnhpres;
    stateb->veta       = statea->veta;
    if (ekinda) 
    {
        copy_mat(ekinda->ekin,ekindb->ekin);
        for (i=0; i<stateb->ngtc; i++) 
        {
            ekindb->tcstat[i].T = ekinda->tcstat[i].T;
            ekindb->tcstat[i].Th = ekinda->tcstat[i].Th;
            copy_mat(ekinda->tcstat[i].ekinh,ekindb->tcstat[i].ekinh);
            copy_mat(ekinda->tcstat[i].ekinf,ekindb->tcstat[i].ekinf);
            ekindb->tcstat[i].ekinscalef_nhc =  ekinda->tcstat[i].ekinscalef_nhc;
            ekindb->tcstat[i].ekinscaleh_nhc =  ekinda->tcstat[i].ekinscaleh_nhc;
            ekindb->tcstat[i].vscale_nhc =  ekinda->tcstat[i].vscale_nhc;
        }
    }
    copy_rvecn(statea->x,stateb->x,0,stateb->natoms);
    copy_rvecn(statea->v,stateb->v,0,stateb->natoms);
    copy_mat(statea->box,stateb->box);
    copy_mat(statea->box_rel,stateb->box_rel);
    copy_mat(statea->boxv,stateb->boxv);

    for (i = 0; i<stateb->ngtc; i++) 
    { 
        nc = i*opts->nhchainlength;
        for (j=0; j<opts->nhchainlength; j++) 
        {
            stateb->nosehoover_xi[nc+j]  = statea->nosehoover_xi[nc+j];
            stateb->nosehoover_vxi[nc+j] = statea->nosehoover_vxi[nc+j];
        }
    }
    if (stateb->nhpres_xi != NULL)
    {
        for (i = 0; i<stateb->nnhpres; i++) 
        {
            nc = i*opts->nhchainlength;
            for (j=0; j<opts->nhchainlength; j++) 
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
        quantity = NPT_energy(ir,state,MassQ);                
        break;
    case etcVRESCALE:
        quantity = vrescale_energy(&(ir->opts),state->therm_integral);
        break;
    default:
        break;
    }
    return quantity;
}

void
compute_globals(FILE *fplog,
                gmx_global_stat_t gstat,
                t_commrec *cr,
                t_inputrec *ir,
                t_forcerec *fr,
                gmx_ekindata_t *ekind,
                t_state *state,
                t_state *state_global,
                t_mdatoms *mdatoms,
                t_nrnb *nrnb,
                t_vcm *vcm,
                gmx_wallcycle_t wcycle,
                gmx_enerdata_t *enerd,
                tensor force_vir,
                tensor shake_vir,
                tensor total_vir,
                tensor pres,
                rvec mu_tot,
                gmx_constr_t constr,
                gmx_signal *signal,
                gmx_bool bIntraSimSignal,
                gmx_bool bInterSimSignal,
                matrix box,
                gmx_mtop_t *top_global,
                real *pcurr,
                int natoms,
                gmx_bool *bSumEkinhOld,
                int flags)
{
    int  i,signal_value,nsignals;
    real signal_buf[esignalNR];
    tensor corr_vir,corr_pres,shakeall_vir;
    gmx_bool bEner,bPres,bTemp, bVV;
    gmx_bool bRerunMD, bStopCM, bGStat, bIterate, 
        bFirstIterate,bReadEkin,bEkinAveVel,bScaleEkin, bConstrain;
    real ekin,temp,prescorr,enercorr,dvdlcorr;
    
    /* translate CGLO flags to gmx_booleans */
    bRerunMD = flags & CGLO_RERUNMD;
    bStopCM = flags & CGLO_STOPCM;
    bGStat = flags & CGLO_GSTAT;

    bReadEkin = (flags & CGLO_READEKIN);
    bScaleEkin = (flags & CGLO_SCALEEKIN);
    bEner = flags & CGLO_ENERGY;
    bTemp = flags & CGLO_TEMPERATURE;
    bPres  = (flags & CGLO_PRESSURE);
    bConstrain = (flags & CGLO_CONSTRAINT);
    bIterate = (flags & CGLO_ITERATE);
    bFirstIterate = (flags & CGLO_FIRSTITERATE);

    /* we calculate a full state kinetic energy either with full-step velocity verlet
       or half step where we need the pressure */
    
    bEkinAveVel = (ir->eI==eiVV || (ir->eI==eiVVAK && bPres) || bReadEkin);
    
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
            accumulate_u(cr,&(ir->opts),ekind);
        }
        debug_gmx();
        if (bReadEkin)
        {
            restore_ekinstate_from_state(cr,ekind,&state_global->ekinstate);
        }
        else 
        {

            calc_ke_part(state,&(ir->opts),mdatoms,ekind,nrnb,bEkinAveVel,bIterate);
        }
        
        debug_gmx();
        
        /* Calculate center of mass velocity if necessary, also parallellized */
        if (bStopCM && !bRerunMD && bEner) 
        {
            calc_vcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms,
                         state->x,state->v,vcm);
        }
    }

    if (bTemp || bPres || bEner || bConstrain) 
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
            /* Handle inter- and intra-simulation signals. Even when
             * there is only one simulation, the "inter-simulation"
             * signals (checkpointing, termination) still have to be
             * handled.
             */

            /* First, construct a buffer containing the signal values. */
            if (NULL != signal)
            {
                for(i = 0; i < esignalINTRA_NR; i++)
                {
                    signal_buf[i] = bIntraSimSignal ? signal->init[i] : 0.0;
                }

                /* i value carries over from above loop */
                for(; i < esignalINTER_NR; i++)
                {
                    signal_buf[i] = bInterSimSignal ? signal->init[i] : 0.0;
                }

                if (bInterSimSignal && MULTISIM(cr) && MASTER(cr))
                {
                    /* Communicate the inter-simulation signals
                     * between the simulations - but only the
                     * inter-simulation part of the buffer should be
                     * communicated. */
                    gmx_sum_sim(esignalINTER_NR - esignalINTRA_NR, signal_buf + esignalINTRA_NR, cr->ms);
                }

                if (bInterSimSignal || bIntraSimSignal)
                {
                    nsignals = esignalNR;
                }
                else
                {
                    nsignals = 0;
                }
            }
            else
            {
                nsignals = 0;
            }

            if (PAR(cr)) 
            {
                wallcycle_start(wcycle,ewcMoveE);
                GMX_MPE_LOG(ev_global_stat_start);
                global_stat(fplog,gstat,cr,enerd,force_vir,shake_vir,mu_tot,
                            ir,ekind,constr,vcm,
                            nsignals,signal_buf,
                            top_global,state,
                            *bSumEkinhOld,flags);
                GMX_MPE_LOG(ev_global_stat_finish);
                wallcycle_stop(wcycle,ewcMoveE);
            }
            if (NULL != signal)
            {
                /* We set the communicated signals only when they are
                 * non-zero, since signals might not be processed at
                 * each MD step.
                 */

                /* First, the intra-simulation signals, only when wanted: */
                if (bIntraSimSignal)
                {
                    for(i = 0; i < esignalINTRA_NR; i++)
                    {
                        signal_value = (signal_buf[i] >= 0 ?
                                        (int)(signal_buf[i] + 0.5) :
                                        (int)(signal_buf[i] - 0.5));
                        if (signal_value != 0)
                        {
                            signal->set[i] = signal_value;
                        }
                        /* Turn off the local signal */
                        signal->init[i] = 0;
                    }
                }

                /* Then, the inter-simulation signals, only when
                 * wanted. Again, the i value carries over from the
                 * above loop.*/
                if (bInterSimSignal)
                {
                    for(; i < esignalINTER_NR; i++)
                    {
                        signal_value = (signal_buf[i] >= 0 ?
                                        (int)(signal_buf[i] + 0.5) :
                                        (int)(signal_buf[i] - 0.5));
                        if (signal_value != 0)
                        {
                            signal->set[i] = signal_value;
                        }
                        /* Turn off the local signal */
                        signal->init[i] = 0;
                    }
                }
            }
            *bSumEkinhOld = FALSE;
        }
    }
    
    if (!ekind->bNEMD && debug && bTemp && (vcm->nr > 0))
    {
        correct_ekin(debug,
                     mdatoms->start,mdatoms->start+mdatoms->homenr,
                     state->v,vcm->group_p[0],
                     mdatoms->massT,mdatoms->tmass,ekind->ekin);
    }
    
    if (bEner) {
        /* Do center of mass motion removal */
        if (bStopCM && !bRerunMD) /* is this correct?  Does it get called too often with this logic? */
        {
            check_cm_grp(fplog,vcm,ir,1);
            do_stopcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms->cVCM,
                          state->x,state->v,vcm);
            inc_nrnb(nrnb,eNR_STOPCM,mdatoms->homenr);
        }

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
        enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,&(enerd->term[F_DKDL]),
                                       bEkinAveVel,bIterate,bScaleEkin);
 
        enerd->term[F_EKIN] = trace(ekind->ekin);
    }
    
    /* ##########  Long range energy information ###### */
    
    if (bEner || bPres || bConstrain) 
    {
        calc_dispcorr(fplog,ir,fr,0,top_global->natoms,box,state->lambda,
                      corr_pres,corr_vir,&prescorr,&enercorr,&dvdlcorr);
    }
    
    if (bEner && bFirstIterate) 
    {
        enerd->term[F_DISPCORR] = enercorr;
        enerd->term[F_EPOT] += enercorr;
        enerd->term[F_DVDL] += dvdlcorr;
        if (fr->efep != efepNO) {
            enerd->dvdl_lin += dvdlcorr;
        }
    }
    
    /* ########## Now pressure ############## */
    if (bPres || bConstrain) 
    {
        
        m_add(force_vir,shake_vir,total_vir);
        
        /* Calculate pressure and apply LR correction if PPPM is used.
         * Use the box from last timestep since we already called update().
         */
        
        enerd->term[F_PRES] = calc_pres(fr->ePBC,ir->nwall,box,ekind->ekin,total_vir,pres,
                                        (fr->eeltype==eelPPPM)?enerd->term[F_COUL_RECIP]:0.0);
        
        /* Calculate long range corrections to pressure and energy */
        /* this adds to enerd->term[F_PRES] and enerd->term[F_ETOT], 
           and computes enerd->term[F_DISPCORR].  Also modifies the 
           total_vir and pres tesors */
        
        m_add(total_vir,corr_vir,total_vir);
        m_add(pres,corr_pres,pres);
        enerd->term[F_PDISPCORR] = prescorr;
        enerd->term[F_PRES] += prescorr;
        *pcurr = enerd->term[F_PRES];
        /* calculate temperature using virial */
        enerd->term[F_VTEMP] = calc_temp(trace(total_vir),ir->opts.nrdf[0]);
        
    }    
}

void check_nst_param(FILE *fplog,const t_commrec *cr,
                     const char *desc_nst,int nst,
                     const char *desc_p,int *p)
{
    char buf[STRLEN];

    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p)/nst + 1)*nst;
        sprintf(buf,"NOTE: %s changes %s to %d\n",desc_nst,desc_p,*p);
        md_print_warning(cr,fplog,buf);
    }
}

void reset_all_counters(FILE *fplog,t_commrec *cr,
                        gmx_large_int_t step,
                        gmx_large_int_t *step_rel,t_inputrec *ir,
                        gmx_wallcycle_t wcycle,t_nrnb *nrnb,
                        gmx_runtime_t *runtime)
{
    char buf[STRLEN],sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    sprintf(buf,"Step %s: resetting all time and cycle counters\n",
            gmx_step_str(step,sbuf));
    md_print_warning(cr,fplog,buf);

    wallcycle_stop(wcycle,ewcRUN);
    wallcycle_reset_all(wcycle);
    if (DOMAINDECOMP(cr))
    {
        reset_dd_statistics_counters(cr->dd);
    }
    init_nrnb(nrnb);
    ir->init_step += *step_rel;
    ir->nsteps    -= *step_rel;
    *step_rel = 0;
    wallcycle_start(wcycle,ewcRUN);
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Restarted time",runtime);
}

void min_zero(int *n,int i)
{
    if (i > 0 && (*n == 0 || i < *n))
    {
        *n = i;
    }
}

int lcd4(int i1,int i2,int i3,int i4)
{
    int nst;

    nst = 0;
    min_zero(&nst,i1);
    min_zero(&nst,i2);
    min_zero(&nst,i3);
    min_zero(&nst,i4);
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

int check_nstglobalcomm(FILE *fplog,const t_commrec *cr,
                        int nstglobalcomm,t_inputrec *ir)
{
    char buf[STRLEN];

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
            sprintf(buf,"WARNING: nstglobalcomm is larger than nstlist, but not a multiple, setting it to %d\n",nstglobalcomm);
            md_print_warning(cr,fplog,buf);
        }
        if (ir->nstcalcenergy > 0)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nstcalcenergy",&ir->nstcalcenergy);
        }
        if (ir->etc != etcNO && ir->nsttcouple > 0)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nsttcouple",&ir->nsttcouple);
        }
        if (ir->epc != epcNO && ir->nstpcouple > 0)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nstpcouple",&ir->nstpcouple);
        }

        check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                        "nstenergy",&ir->nstenergy);

        check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                        "nstlog",&ir->nstlog);
    }

    if (ir->comm_mode != ecmNO && ir->nstcomm < nstglobalcomm)
    {
        sprintf(buf,"WARNING: Changing nstcomm from %d to %d\n",
                ir->nstcomm,nstglobalcomm);
        md_print_warning(cr,fplog,buf);
        ir->nstcomm = nstglobalcomm;
    }

    if (fplog)
    {
        fprintf(fplog, "Global intra-simulation communication will occur every %d steps.\n", nstglobalcomm);
    }
    return nstglobalcomm;
}

int
check_nstsignalcomm(FILE *fplog,
                    const t_commrec *cr,
                    int nstsignalcomm,
                    const t_inputrec *ir,
                    real max_hours)
{
    char buf[STRLEN];
    int new_nstsignalcomm;

    if (-1 == nstsignalcomm)
    {
        real temp_float = ir->nstcalcenergy;
        /* The inter-simulation communication is used for things like
         * synchronizing checkpointing, dealing with Unix signals, and
         * run termination. Accordingly its frequency does not need to
         * reflect nstlist, nstcalcenergy, the -gcom setting or the
         * number of processors. Ideally, we want it to take place at
         * a constant time interval, but we do not know when a
         * simulation starts how long each time step will take, and we
         * probably do not want to be polling for system time either.
         *
         * So, if we have to guess how often to do inter-simulation
         * synchronization, we may as well make the asymptotic cost of
         * synchronization constant. Assuming a tree implementation of
         * global communication, the synchronization cost has
         * components that scale
         *
         *   as O(nsim * log(nsim)) for the phase that broadcasts
         * across the simulation masters, and
         *
         *   as O(nnodes * log(nnodes)) for the subsequent
         * intra-simulation broadcast phase.
         *
         * However, we would like to be confident that a signal check
         * will occur in the interval between 0.99 * max_hours and
         * max_hours, so that that termination mechanism works. Same
         * for checkpoint intervals. However it's hard to see a
         * general way of doing that, since we can't know a priori how
         * much walltime will elapse for a given number of simulation
         * steps. The user can always set nstsignalcomm by hand if
         * they have problems with the default, and we provide a hint
         * when that might be necessary. */
        if (1 < cr->nnodes)
        {
            temp_float *= cr->nnodes * log2((double)cr->nnodes);
        }
        if (NULL != cr->ms && 1 < cr->ms->nsim)
        {
            temp_float *= cr->ms->nsim * log2((double)cr->ms->nsim);
        }
        /* So for dt = 0.002, nstcalcenergy = 10 and 512 processors, we
         * transmit inter-simulation signals every 92ps. This seems
         * like it might be a bit too infrequent.
         *
         * Hence we use a fudge factor, which we can adjust to taste
         * as hardware gets faster. */
        temp_float *= 1;
        new_nstsignalcomm = (int) temp_float;
    }
    else
    {
        new_nstsignalcomm = nstsignalcomm;
    }

    /* We want to do inter-simulation communication at times when
     * intra-simulation communication is already happening, so as not
     * to introduce barriers. */
    check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                    "-signalcom", &new_nstsignalcomm);
    if (fplog && 0 < max_hours)
    {
        fprintf(fplog, "If you have problems with the mdrun -maxh mechanism not stopping runs in time,\n inter-simulation signalling may be too infrequent");
        if (nstsignalcomm != new_nstsignalcomm)
        {
            fprintf(fplog, " (mdrun chose an interval of %d steps)", new_nstsignalcomm);
        }
        fprintf(fplog, ".\n Try a smaller value with mdrun -signalcom.\n\n");
    }

    return new_nstsignalcomm;
}

void check_ir_old_tpx_versions(t_commrec *cr,FILE *fplog,
                               t_inputrec *ir,gmx_mtop_t *mtop)
{
    /* Check required for old tpx files */
    if (IR_TWINRANGE(*ir) && ir->nstlist > 1 &&
        ir->nstcalcenergy % ir->nstlist != 0)
    {
        md_print_warning(cr,fplog,"Old tpr file with twin-range settings: modifying energy calculation and/or T/P-coupling frequencies");

        if (gmx_mtop_ftype_count(mtop,F_CONSTR) +
            gmx_mtop_ftype_count(mtop,F_CONSTRNC) > 0 &&
            ir->eConstrAlg == econtSHAKE)
        {
            md_print_warning(cr,fplog,"With twin-range cut-off's and SHAKE the virial and pressure are incorrect");
            if (ir->epc != epcNO)
            {
                gmx_fatal(FARGS,"Can not do pressure coupling with twin-range cut-off's and SHAKE");
            }
        }
        check_nst_param(fplog,cr,"nstlist",ir->nstlist,
                        "nstcalcenergy",&ir->nstcalcenergy);
        if (ir->epc != epcNO)
        {
            check_nst_param(fplog,cr,"nstlist",ir->nstlist,
                            "nstpcouple",&ir->nstpcouple);
        }
        check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                        "nstenergy",&ir->nstenergy);
        check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                        "nstlog",&ir->nstlog);
        if (ir->efep != efepNO)
        {
            check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                            "nstdhdl",&ir->nstdhdl);
        }
    }
}

void rerun_parallel_comm(t_commrec *cr,t_trxframe *fr,
                         gmx_bool *bNotLastFrame)
{
    gmx_bool bAlloc;
    rvec *xp,*vp;

    bAlloc = (fr->natoms == 0);

    if (MASTER(cr) && !*bNotLastFrame)
    {
        fr->natoms = -1;
    }
    xp = fr->x;
    vp = fr->v;
    gmx_bcast(sizeof(*fr),fr,cr);
    fr->x = xp;
    fr->v = vp;

    *bNotLastFrame = (fr->natoms >= 0);

    if (*bNotLastFrame && PARTDECOMP(cr))
    {
        /* x and v are the only variable size quantities stored in trr
         * that are required for rerun (f is not needed).
         */
        if (bAlloc)
        {
            snew(fr->x,fr->natoms);
            snew(fr->v,fr->natoms);
        }
        if (fr->bX)
        {
            gmx_bcast(fr->natoms*sizeof(fr->x[0]),fr->x[0],cr);
        }
        if (fr->bV)
        {
            gmx_bcast(fr->natoms*sizeof(fr->v[0]),fr->v[0],cr);
        }
    }
}

void md_print_warning(const t_commrec *cr,FILE *fplog,const char *buf)
{
    if (MASTER(cr))
    {
        fprintf(stderr,"\n%s\n",buf);
    }
    if (fplog)
    {
        fprintf(fplog,"\n%s\n",buf);
    }
}
