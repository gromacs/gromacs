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

#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "vsite.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "xtcio.h"
#include "mdrun.h"
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"
#include "ionize.h"
#include "disre.h"
#include "orires.h"
#include "dihre.h"
#include "pppm.h"
#include "pme.h"
#include "mdatoms.h"
#include "repl_ex.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "genborn.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

/* The following two variables and the signal_handler function
 * is used from pme.c as well 
 */
extern bool bGotTermSignal, bGotUsr1Signal;

static RETSIGTYPE signal_handler(int n)
{
    switch (n) {
        case SIGTERM:
            bGotTermSignal = TRUE;
            break;
#ifdef HAVE_SIGUSR1
        case SIGUSR1:
            bGotUsr1Signal = TRUE;
            break;
#endif
    }
}

typedef struct { 
    gmx_integrator_t *func;
} gmx_intp_t;

/* The array should match the eI array in include/types/enums.h */
const gmx_intp_t integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md}, {do_md}};

/* Static variables for temporary use with the deform option */
static int    init_step_tpx;
static matrix box_tpx;
#ifdef GMX_THREADS
static tMPI_Thread_mutex_t box_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


static void copy_coupling_state(t_state *statea,t_state *stateb, 
                                gmx_ekindata_t *ekinda,gmx_ekindata_t *ekindb) 
{
    
    /* MRS note -- might be able to get rid of some of the arguments.  Look over it when it's all debugged */
    
    int i,j,ngtc_eff;
    
    stateb->natoms     = statea->natoms;
    stateb->ngtc       = statea->ngtc;
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
    /* need extra term for the barostat */
    ngtc_eff = stateb->ngtc+1;

    for (i = 0; i < ngtc_eff; i++) 
    { 
        for (j=0; j < NNHCHAIN; j++) 
        {
            stateb->nosehoover_xi[i*NNHCHAIN + j]       = statea->nosehoover_xi[i*NNHCHAIN + j];
            stateb->nosehoover_vxi[i*NNHCHAIN + j]      = statea->nosehoover_vxi[i*NNHCHAIN + j];
        }
    }
}

static void compute_globals(FILE *fplog, gmx_global_stat_t gstat, t_commrec *cr, t_inputrec *ir, 
                            t_forcerec *fr, gmx_ekindata_t *ekind, 
                            t_state *state, t_state *state_global, t_mdatoms *mdatoms, 
                            t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
                            gmx_enerdata_t *enerd,tensor force_vir, tensor shake_vir, tensor total_vir, 
                            tensor pres, rvec mu_tot, gmx_constr_t constr, 
                            real *chkpt,real *terminate, real *terminate_now,
                            int *nabnsb, matrix box, gmx_mtop_t *top_global, real *pcurr, 
                            int natoms, bool *bSumEkinhOld, int flags)
{
    
    int i;
    tensor corr_vir,corr_pres,shakeall_vir;
    bool bEner,bPres,bTemp, bVV;
    bool bRerunMD, bEkinAveVel, bStopCM, bGStat, bNEMD, bIterate, 
        bFirstIterate,bReadEkin,bScaleEkin, bConstrain;
    real ekin,temp,prescorr,enercorr,dvdlcorr;
    
    /* translate CGLO flags to booleans */
    bRerunMD = flags & CGLO_RERUNMD;
    bEkinAveVel = flags & CGLO_EKINAVEVEL;
    bStopCM = flags & CGLO_STOPCM;
    bGStat = flags & CGLO_GSTAT;
    bNEMD = flags & CGLO_NEMD;
    bReadEkin = flags & CGLO_READEKIN;
    bScaleEkin = flags & CGLO_SCALEEKIN;
    bEner = flags & CGLO_ENERGY;
    bTemp = flags & CGLO_TEMPERATURE;
    bPres  = flags & CGLO_PRESSURE;
    bConstrain = flags & CGLO_CONSTRAINT;
    bIterate = flags & CGLO_ITERATE;
    bFirstIterate = flags & CGLO_FIRSTITERATE;

    /* in initalization, it sums the shake virial in vv, and to 
       sums ekinh_old in leapfrog (or if we are calculating ekinh_old for other reasons */
    
    bVV = (ir->eI == eiVV);
    
    /* ########## Kinetic energy  ############## */
    
    if (bTemp) 
    {
        /* Non-equilibrium MD: this is parallellized, but only does communication
         * when there really is NEMD.
         */
        
        if (PAR(cr) && (bNEMD)) 
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
            if (PAR(cr)) 
            {
                GMX_MPE_LOG(ev_global_stat_start);
                global_stat(fplog,gstat,cr,enerd,force_vir,shake_vir,mu_tot,
                            ir,ekind,constr,vcm,NULL,NULL,terminate,top_global,state,
                            *bSumEkinhOld,flags);
                GMX_MPE_LOG(ev_global_stat_finish);
            }
            
            if (terminate != 0)
            {
                terminate_now = terminate;
                terminate = 0;
            }
            wallcycle_stop(wcycle,ewcMoveE);
            *bSumEkinhOld = FALSE;
        }
    }
    
    if (!bNEMD && debug && bTemp && (vcm->nr > 0))
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
    }

    if (bTemp) 
    {
        /* Sum the kinetic energies of the groups & calc temp */
        enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,&(enerd->term[F_DKDL]),
                                       bEkinAveVel||bReadEkin,bIterate,bScaleEkin);
        /* three main: VV with AveVel, vv with AveEkin, leap with AveEkin.  Leap with AveVel is also
           an option, but not supported now.  Additionally, if we are doing iterations.  
           bCopyHalf: if TRUE, we simply copy the ekinh directly to ekin, multiplying be ekinscale_nhc.
           bEkinAveVel: If TRUE, we simply multiply ekin by ekinscale. 
           If FALSE, we average ekinh_old and ekinh*ekinscale_nhc.
           bSaveEkinOld: If TRUE (in the case of iteration = bIterate is TRUE), we don't reset the ekinscale_nhc.  
           If FALSE, we go ahead and erase.
        */
        enerd->term[F_EKIN] = trace(ekind->ekin);
    }
    
    /* ##########  Long range energy information ###### */
    
    if (bEner || bPres || bConstrain) 
    {
        calc_dispcorr(fplog,ir,fr,0,top_global,box,state->lambda,
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


#ifdef GMX_THREADS
struct mdrunner_arglist
{
    FILE *fplog;
    t_commrec *cr;
    int nfile;
    const t_filenm *fnm;
    output_env_t oenv;
    bool bVerbose;
    bool bCompact;
    int nstglobalcomm;
    ivec ddxyz;
    int dd_node_order;
    real rdd;
    real rconstr;
    const char *dddlb_opt;
    real dlb_scale;
    const char *ddcsx;
    const char *ddcsy;
    const char *ddcsz;
    int nstepout;
    int nmultisim;
    int repl_ex_nst;
    int repl_ex_seed;
    real pforce;
    real cpt_period;
    real max_hours;
    unsigned long Flags;
    int ret; /* return value */
};


static void mdrunner_start_fn(void *arg)
{
    struct mdrunner_arglist *mda=(struct mdrunner_arglist*)arg;
    struct mdrunner_arglist mc=*mda; /* copy the arg list to make sure 
                                        that it's thread-local. This doesn't
                                        copy pointed-to items, of course,
                                        but those are all const. */
    t_commrec *cr;  /* we need a local version of this */
    FILE *fplog=NULL;
    t_filenm *fnm=dup_tfn(mc.nfile, mc.fnm);

    cr=init_par_threads(mc.cr);
    if (MASTER(cr))
    {
        fplog=mc.fplog;
    }


    mda->ret=mdrunner(fplog, cr, mc.nfile, mc.fnm, mc.oenv, mc.bVerbose,
                      mc.bCompact, mc.nstglobalcomm, 
                      mc.ddxyz, mc.dd_node_order, mc.rdd,
                      mc.rconstr, mc.dddlb_opt, mc.dlb_scale, 
                      mc.ddcsx, mc.ddcsy, mc.ddcsz, mc.nstepout, mc.nmultisim,
                      mc.repl_ex_nst, mc.repl_ex_seed, mc.pforce, 
                      mc.cpt_period, mc.max_hours, mc.Flags);
}

#endif

int mdrunner_threads(int nthreads, 
                     FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
                     const output_env_t oenv, bool bVerbose,bool bCompact,
                     int nstglobalcomm,
                     ivec ddxyz,int dd_node_order,real rdd,real rconstr,
                     const char *dddlb_opt,real dlb_scale,
                     const char *ddcsx,const char *ddcsy,const char *ddcsz,
                     int nstepout,int nmultisim,int repl_ex_nst,
                     int repl_ex_seed, real pforce,real cpt_period,
                     real max_hours, unsigned long Flags)
{
    int ret;
    /* first check whether we even need to start tMPI */
    if (nthreads < 2)
    {
        ret=mdrunner(fplog, cr, nfile, fnm, oenv, bVerbose, bCompact,
                     nstglobalcomm,
                     ddxyz, dd_node_order, rdd, rconstr, dddlb_opt, dlb_scale,
                     ddcsx, ddcsy, ddcsz, nstepout, nmultisim, repl_ex_nst, 
                     repl_ex_seed, pforce, cpt_period, max_hours, Flags);
    }
    else
    {
#ifdef GMX_THREADS
        struct mdrunner_arglist mda;
        /* fill the data structure to pass as void pointer to thread start fn */
        mda.fplog=fplog;
        mda.cr=cr;
        mda.nfile=nfile;
        mda.fnm=fnm;
        mda.oenv=oenv;
        mda.bVerbose=bVerbose;
        mda.bCompact=bCompact;
        mda.nstglobalcomm=nstglobalcomm;
        mda.ddxyz[XX]=ddxyz[XX];
        mda.ddxyz[YY]=ddxyz[YY];
        mda.ddxyz[ZZ]=ddxyz[ZZ];
        mda.dd_node_order=dd_node_order;
        mda.rdd=rdd;
        mda.rconstr=rconstr;
        mda.dddlb_opt=dddlb_opt;
        mda.dlb_scale=dlb_scale;
        mda.ddcsx=ddcsx;
        mda.ddcsy=ddcsy;
        mda.ddcsz=ddcsz;
        mda.nstepout=nstepout;
        mda.nmultisim=nmultisim;
        mda.repl_ex_nst=repl_ex_nst;
        mda.repl_ex_seed=repl_ex_seed;
        mda.pforce=pforce;
        mda.cpt_period=cpt_period;
        mda.max_hours=max_hours;
        mda.Flags=Flags;

        fprintf(stderr, "Starting %d threads\n",nthreads);
        fflush(stderr);
        tMPI_Init_fn(nthreads, mdrunner_start_fn, (void*)(&mda) );
        ret=mda.ret;
#else
        ret=-1;
        gmx_comm("Multiple threads requested but not compiled with threads");
#endif
    }
    return ret;
}


int mdrunner(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, bool bVerbose,bool bCompact,
             int nstglobalcomm,
             ivec ddxyz,int dd_node_order,real rdd,real rconstr,
             const char *dddlb_opt,real dlb_scale,
             const char *ddcsx,const char *ddcsy,const char *ddcsz,
             int nstepout,int nmultisim,int repl_ex_nst,int repl_ex_seed,
             real pforce,real cpt_period,real max_hours,
             unsigned long Flags)
{
    double     nodetime=0,realtime;
    t_inputrec *inputrec;
    t_state    *state=NULL;
    matrix     box;
    gmx_ddbox_t ddbox;
    int        npme_major;
    real       tmpr1,tmpr2;
    t_nrnb     *nrnb;
    gmx_mtop_t *mtop=NULL;
    t_mdatoms  *mdatoms=NULL;
    t_forcerec *fr=NULL;
    t_fcdata   *fcd=NULL;
    real       ewaldcoeff=0;
    gmx_pme_t  *pmedata=NULL;
    gmx_vsite_t *vsite=NULL;
    gmx_constr_t constr;
    int        i,m,nChargePerturbed=-1,status,nalloc;
    char       *gro;
    gmx_wallcycle_t wcycle;
    bool       bReadRNG,bReadEkin;
    int        list;
    gmx_runtime_t runtime;
    int        rc;
    gmx_large_int_t reset_counters;
    gmx_edsam_t ed;

    /* Essential dynamics */
    if (opt2bSet("-ei",nfile,fnm)) 
    {
        /* Open input and output files, allocate space for ED data structure */
        ed = ed_open(nfile,fnm,cr);
    } 
    else
        ed=NULL;

    snew(inputrec,1);
    snew(mtop,1);

    if (bVerbose && SIMMASTER(cr))
    {
        fprintf(stderr,"Getting Loaded...\n");
    }

    if (Flags & MD_APPENDFILES) 
    {
        fplog = NULL;
    }

    if (PAR(cr))
    {
        /* The master thread on the master node reads from disk, 
         * then distributes everything to the other processors.
         */

        list = (SIMMASTER(cr) && !(Flags & MD_APPENDFILES)) ?  (LIST_SCALARS | LIST_INPUTREC) : 0;

        snew(state,1);
        init_parallel(fplog, opt2fn_master("-s",nfile,fnm,cr),cr,
                      inputrec,mtop,state,list);

    }
    else
    {
        /* Read a file for a single processor */
        snew(state,1);
        init_single(fplog,inputrec,ftp2fn(efTPX,nfile,fnm),mtop,state);
    }
    if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC))
    {
        cr->npmenodes = 0;
    }

#ifdef GMX_FAHCORE
    fcRegisterSteps(inputrec->nsteps,inputrec->init_step);
#endif

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */
    snew(fcd,1);

    /* This needs to be called before read_checkpoint to extend the state */
    init_disres(fplog,mtop,inputrec,cr,Flags & MD_PARTDEC,fcd,state);

    if (gmx_mtop_ftype_count(mtop,F_ORIRES) > 0)
    {
        if (PAR(cr) && !(Flags & MD_PARTDEC))
        {
            gmx_fatal(FARGS,"Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)");
        }
        /* Orientation restraints */
        if (MASTER(cr))
        {
            init_orires(fplog,mtop,state->x,inputrec,cr->ms,&(fcd->orires),
                        state);
        }
    }

    if (DEFORM(*inputrec))
    {
        /* Store the deform reference box before reading the checkpoint */
        if (SIMMASTER(cr))
        {
            copy_mat(state->box,box);
        }
        if (PAR(cr))
        {
            gmx_bcast(sizeof(box),box,cr);
        }
        /* Because we do not have the update struct available yet
         * in which the reference values should be stored,
         * we store them temporarily in static variables.
         * This should be thread safe, since they are only written once
         * and with identical values.
         */
#ifdef GMX_THREADS
        tMPI_Thread_mutex_lock(&box_mutex);
#endif
        init_step_tpx = inputrec->init_step;
        copy_mat(box,box_tpx);
#ifdef GMX_THREADS
        tMPI_Thread_mutex_unlock(&box_mutex);
#endif
    }

    if (opt2bSet("-cpi",nfile,fnm)) 
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        if( gmx_fexist_master(opt2fn_master("-cpi",nfile,fnm,cr),cr) )
        {
            load_checkpoint(opt2fn_master("-cpi",nfile,fnm,cr),fplog,
                            cr,Flags & MD_PARTDEC,ddxyz,
                            inputrec,state,&bReadRNG,&bReadEkin,
                            (Flags & MD_APPENDFILES));
            
            if (bReadRNG)
            {
                Flags |= MD_READ_RNG;
            }
            if (bReadEkin)
            {
                Flags |= MD_READ_EKIN;
            }
        }
    }

    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        fplog = gmx_log_open(ftp2fn(efLOG,nfile,fnm),cr,!(Flags & MD_SEPPOT),
                             Flags);
    }

    if (SIMMASTER(cr)) 
    {
        copy_mat(state->box,box);
    }

    if (PAR(cr)) 
    {
        gmx_bcast(sizeof(box),box,cr);
    }

    if (bVerbose && SIMMASTER(cr))
    {
        fprintf(stderr,"Loaded with Money\n\n");
    }

    if (PAR(cr) && !((Flags & MD_PARTDEC) || EI_TPI(inputrec->eI)))
    {
        cr->dd = init_domain_decomposition(fplog,cr,Flags,ddxyz,rdd,rconstr,
                                           dddlb_opt,dlb_scale,
                                           ddcsx,ddcsy,ddcsz,
                                           mtop,inputrec,
                                           box,state->x,
                                           &ddbox,&npme_major);

        make_dd_communicators(fplog,cr,dd_node_order);

        /* Set overallocation to avoid frequent reallocation of arrays */
        set_over_alloc_dd(TRUE);
    }
    else
    {
        cr->duty = (DUTY_PP | DUTY_PME);
        npme_major = cr->nnodes;
        
        if (inputrec->ePBC == epbcSCREW)
        {
            gmx_fatal(FARGS,
                      "pbc=%s is only implemented with domain decomposition",
                      epbc_names[inputrec->ePBC]);
        }
    }

    if (PAR(cr))
    {
        /* After possible communicator splitting in make_dd_communicators.
         * we can set up the intra/inter node communication.
         */
        gmx_setup_nodecomm(fplog,cr);
    }

    wcycle = wallcycle_init(fplog,cr);
    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes 
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast_sim(sizeof(reset_counters),&reset_counters,cr);
        wcycle_set_reset_counters(wcycle, reset_counters);
    }


    snew(nrnb,1);
    if (cr->duty & DUTY_PP)
    {
        /* For domain decomposition we allocate dynamically
         * in dd_partition_system.
         */
        if (DOMAINDECOMP(cr))
        {
            bcast_state_setup(cr,state);
        }
        else
        {
            if (PAR(cr))
            {
                if (!MASTER(cr))
                {
                    snew(state,1);
                }
                bcast_state(cr,state,TRUE);
            }
        }

        /* Dihedral Restraints */
        if (gmx_mtop_ftype_count(mtop,F_DIHRES) > 0)
        {
            init_dihres(fplog,mtop,inputrec,fcd);
        }

        /* Initiate forcerecord */
        fr = mk_forcerec();
        init_forcerec(fplog,oenv,fr,fcd,inputrec,mtop,cr,box,FALSE,
                      opt2fn("-table",nfile,fnm),
                      opt2fn("-tablep",nfile,fnm),
                      opt2fn("-tableb",nfile,fnm),FALSE,pforce);

        /* version for PCA_NOT_READ_NODE (see md.c) */
        /*init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
          "nofile","nofile","nofile",FALSE,pforce);
          */        
        fr->bSepDVDL = ((Flags & MD_SEPPOT) == MD_SEPPOT);

        /* Initialize QM-MM */
        if(fr->bQMMM)
        {
            init_QMMMrec(cr,box,mtop,inputrec,fr);
        }

        /* Initialize the mdatoms structure.
         * mdatoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdatoms = init_mdatoms(fplog,mtop,inputrec->efep!=efepNO);

        /* Initialize the virtual site communication */
        vsite = init_vsite(mtop,cr);

        calc_shifts(box,fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(cr) &&
            !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->ePBC != epbcNONE)
            {
                do_pbc_first_mtop(fplog,inputrec->ePBC,box,mtop,state->x);
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                construct_vsites_mtop(fplog,vsite,mtop,state->x);
            }
        }

        /* Initiate PPPM if necessary */
        if (fr->eeltype == eelPPPM)
        {
            if (mdatoms->nChargePerturbed)
            {
                gmx_fatal(FARGS,"Free energy with %s is not implemented",
                          eel_names[fr->eeltype]);
            }
            status = gmx_pppm_init(fplog,cr,oenv,FALSE,TRUE,box,
                                   getenv("GMXGHAT"),inputrec, (Flags & MD_REPRODUCIBLE));
            if (status != 0)
            {
                gmx_fatal(FARGS,"Error %d initializing PPPM",status);
            }
        }

        if (EEL_PME(fr->eeltype))
        {
            ewaldcoeff = fr->ewaldcoeff;
            pmedata = &fr->pmedata;
        }
        else
        {
            pmedata = NULL;
        }
    }
    else
    {
        /* This is a PME only node */

        /* We don't need the state */
        done_state(state);

        ewaldcoeff = calc_ewaldcoeff(inputrec->rcoulomb, inputrec->ewald_rtol);
        snew(pmedata,1);
    }

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype))
    {
        if (mdatoms)
        {
            nChargePerturbed = mdatoms->nChargePerturbed;
        }
        if (cr->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed */
            gmx_bcast_sim(sizeof(nChargePerturbed),&nChargePerturbed,cr);
        }
        if (cr->duty & DUTY_PME)
        {
            status = gmx_pme_init(pmedata,cr,npme_major,inputrec,
                                  mtop ? mtop->natoms : 0,nChargePerturbed,
                                  (Flags & MD_REPRODUCIBLE));
            if (status != 0) 
            {
                gmx_fatal(FARGS,"Error %d initializing PME",status);
            }
        }
    }


    if (integrator[inputrec->eI].func == do_md)
    {
        /* Turn on signal handling on all nodes */
        /*
         * (A user signal from the PME nodes (if any)
         * is communicated to the PP nodes.
         */
        if (getenv("GMX_NO_TERM") == NULL)
        {
            if (debug)
            {
                fprintf(debug,"Installing signal handler for SIGTERM\n");
            }
            signal(SIGTERM,signal_handler);
        }
#ifdef HAVE_SIGUSR1
        if (getenv("GMX_NO_USR1") == NULL)
        {
            if (debug)
            {
                fprintf(debug,"Installing signal handler for SIGUSR1\n");
            }
            signal(SIGUSR1,signal_handler);
        }
#endif
    }

    if (cr->duty & DUTY_PP)
    {
        if (inputrec->ePull != epullNO)
        {
            /* Initialize pull code */
            init_pull(fplog,inputrec,nfile,fnm,mtop,cr,oenv,
                      EI_DYNAMICS(inputrec->eI) && MASTER(cr),Flags);
        }

        constr = init_constraints(fplog,mtop,inputrec,ed,state,cr);

        if (DOMAINDECOMP(cr))
        {
            dd_init_bondeds(fplog,cr->dd,mtop,vsite,constr,inputrec,
                            Flags & MD_DDBONDCHECK,fr->cginfo_mb);

            set_dd_parameters(fplog,cr->dd,dlb_scale,inputrec,fr,&ddbox);

            setup_dd_grid(fplog,cr->dd);
        }

        /* Now do whatever the user wants us to do (how flexible...) */
        integrator[inputrec->eI].func(fplog,cr,nfile,fnm,
                                      oenv,bVerbose,bCompact,
                                      nstglobalcomm,
                                      vsite,constr,
                                      nstepout,inputrec,mtop,
                                      fcd,state,
                                      mdatoms,nrnb,wcycle,ed,fr,
                                      repl_ex_nst,repl_ex_seed,
                                      cpt_period,max_hours,
                                      Flags,
                                      &runtime);

        if (inputrec->ePull != epullNO)
        {
            finish_pull(fplog,inputrec->pull);
        }
    } 
    else 
    {
        /* do PME only */
        gmx_pmeonly(*pmedata,cr,nrnb,wcycle,ewaldcoeff,FALSE,inputrec);
    }

    if (EI_DYNAMICS(inputrec->eI) || EI_TPI(inputrec->eI))
    {
        /* Some timing stats */  
        if (MASTER(cr))
        {
            if (runtime.proc == 0)
            {
                runtime.proc = runtime.real;
            }
        }
        else
        {
            runtime.real = 0;
        }
    }

    wallcycle_stop(wcycle,ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again 
     */
    finish_run(fplog,cr,ftp2fn(efSTO,nfile,fnm),
               inputrec,nrnb,wcycle,&runtime,
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

    /* Does what it says */  
    print_date_and_time(fplog,cr->nodeid,"Finished mdrun",&runtime);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_close(fplog);
    }	

    if(bGotTermSignal)
    {
        rc = 1;
    }
    else if(bGotUsr1Signal)
    {
        rc = 2;
    }
    else
    {
        rc = 0;
    }

    return rc;
}

static void md_print_warning(const t_commrec *cr,FILE *fplog,const char *buf)
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

static bool done_iterating(const t_commrec *cr,FILE *fplog, bool *bFirstIterate, bool *bIterate, real fom, real *newf, int n) 
{


    
    /* monitor convergence, and use a secant search to propose new
       values.  
                                                                  x_{i} - x_{i-1}
       The secant method computes x_{i+1} = x_{i} - f(x_{i}) * ---------------------
                                                                f(x_{i}) - f(x_{i-1})
       
       The function we are trying to zero is fom-x, where fom is the
       "figure of merit" which is the pressure (or the veta value) we
       would get by putting in an old value of the pressure or veta into
       the incrementor function for the step or half step.  I have
       verified that this gives the same answer as self consistent
       iteration, usually in many fewer steps, especially for small tau_p.
       
       We could possibly eliminate an iteration with proper use
       of the value from the previous step, but that would take a bit
       more bookkeeping, especially for veta, since tests indicate the
       function of veta on the last step is not sufficiently close to
       guarantee convergence this step. This is
       good enough for now.  On my tests, I could use tau_p down to
       0.02, which is smaller that would ever be necessary in
       practice. Generally, 3-5 iterations will be sufficient */

/* Definitions for convergence.  Could be optimized . . .  */
#ifdef GMX_DOUBLE
#define CONVERGEITER  0.000000001
#else
#define CONVERGEITER  0.0001
#endif
#define MAXITERCONST       200

/* used to escape out of cyclic traps because of limited numberical precision  */
#define CYCLEMAX            20
#define FALLBACK          1000 

    static real f,fprev,x,xprev;  
    static int iter_i;
    static real allrelerr[MAXITERCONST+2];
    double relerr,xmin;
    char buf[256];
    int i;
    bool incycle;

    if (*bFirstIterate) 
    {
        *bFirstIterate = FALSE;
        iter_i = 0;
        x = fom;
        f = fom-x;
        *newf = fom;
    } 
    else 
    {
        f = fom-x; /* we want to zero this difference */
        if ((iter_i > 1) && (iter_i < MAXITERCONST)) 
        {
            if (f==fprev) 
            {
                *newf = f;
            } 
            else 
            {
                *newf = x - (x-xprev)*(f)/(f-fprev); 
            }
        } 
        else 
        {
            /* just use self-consistent iteration the first step to initialize, or 
               if it's not converging (which happens occasionally -- need to investigate why) */
            *newf = fom; 
        }
    }
    /* Consider a slight shortcut allowing us to exit one sooner -- we check the
       difference between the closest of x and xprev to the new
       value. To be 100% certain, we should check the difference between
       the last result, and the previous result, or
       
       relerr = (fabs((x-xprev)/fom));
       
       but this is pretty much never necessary under typical conditions.
       Checking numerically, it seems to lead to almost exactly the same
       trajectories, but there are small differences out a few decimal
       places in the pressure, and eventually in the v_eta, but it could
       save an interation.
       
       if (fabs(*newf-x) < fabs(*newf - xprev)) { xmin = x;} else { xmin = xprev;}
       relerr = (fabs((*newf-xmin) / *newf));
    */
    
    relerr = (fabs((f-fprev)/fom));
    
    allrelerr[iter_i] = relerr;

    if (iter_i > 0) 
    {
        if (debug) 
        {
            fprintf(debug,"Iterating NPT constraints #%i: %6i %20.12f%14.6g%20.12f\n",n,iter_i,fom,relerr,*newf);
        }
        
        if ((relerr < CONVERGEITER) || (fom==0))
        {
            *bIterate = FALSE;
            if (debug) 
            {
                fprintf(debug,"Iterating NPT constraints #%i: CONVERGED\n",n);
            }
            return TRUE;
        }
        if (iter_i > MAXITERCONST) 
        {
            incycle = FALSE;
            for (i=0;i<CYCLEMAX;i++) {
                if (allrelerr[(MAXITERCONST-2*CYCLEMAX)-i] == allrelerr[iter_i-1]) {
                    incycle = TRUE;
                    if (relerr > CONVERGEITER*FALLBACK) {incycle = FALSE; break;}
                }
            }
            
            if (incycle) {
                /* we are trapped in a numerical attractor, and can't converge any more, and are close to the final result.
                   Better to give up here and proceed with a warning than have the simulation die.
                */
                sprintf(buf,"numerical convergence issues with NPT, relative error only %10.5g, continuing\n",relerr);
                md_print_warning(cr,fplog,buf);
                return TRUE;
            } else {
                gmx_fatal(FARGS,"Could not converge NPT constraints\n");
            }
        }
    }
    
    xprev = x;
    x = *newf;
    fprev = f;
    iter_i++;
    
    return FALSE;
}

static void check_nst_param(FILE *fplog,t_commrec *cr,
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

static void reset_all_counters(FILE *fplog,t_commrec *cr,
                               gmx_large_int_t step,
                               gmx_large_int_t *step_rel,t_inputrec *ir,
                               gmx_wallcycle_t wcycle,t_nrnb *nrnb,
                               gmx_runtime_t *runtime)
{
    char buf[STRLEN],sbuf[22];

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

static int check_nstglobalcomm(FILE *fplog,t_commrec *cr,
                               int nstglobalcomm,t_inputrec *ir)
{
    char buf[STRLEN];

    if (!EI_DYNAMICS(ir->eI))
    {
        nstglobalcomm = 1;
    }

    if (nstglobalcomm == -1)
    {
        if (ir->nstcalcenergy == 0 && ir->nstlist == 0)
        {
            nstglobalcomm = 10;
            if (ir->nstenergy > 0 && ir->nstenergy < nstglobalcomm)
            {
                nstglobalcomm = ir->nstenergy;
            }
        }
        else
        {
            /* We assume that if nstcalcenergy > nstlist,
             * nstcalcenergy is a multiple of nstlist.
             */
            if (ir->nstcalcenergy == 0 ||
                (ir->nstlist > 0 && ir->nstlist < ir->nstcalcenergy))
            {
                nstglobalcomm = ir->nstlist;
            }
            else
            {
                nstglobalcomm = ir->nstcalcenergy;
            }
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
        if (nstglobalcomm > ir->nstcalcenergy)
        {
            check_nst_param(fplog,cr,"-gcom",nstglobalcomm,
                            "nstcalcenergy",&ir->nstcalcenergy);
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

    return nstglobalcomm;
}

static void check_ir_old_tpx_versions(t_commrec *cr,FILE *fplog,
                                      t_inputrec *ir,gmx_mtop_t *mtop)
{
    /* Check required for old tpx files */
    if (IR_TWINRANGE(*ir) && ir->nstlist > 1 &&
        ir->nstcalcenergy % ir->nstlist != 0)
    {
        md_print_warning(cr,fplog,"Old tpr file with twin-range settings: modifiying energy calculation and/or T/P-coupling frequencies");

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
    }
    check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                    "nstenergy",&ir->nstenergy);
    check_nst_param(fplog,cr,"nstcalcenergy",ir->nstcalcenergy,
                    "nstlog",&ir->nstlog);
    if (ir->efep != efepNO)
    {
        check_nst_param(fplog,cr,"nstdhdl",ir->nstdhdl,
                        "nstenergy",&ir->nstenergy);
    }
}

typedef struct {
    bool       bGStatEveryStep;
    gmx_large_int_t step_ns;
    gmx_large_int_t step_nscheck;
    gmx_large_int_t nns;
    matrix     scale_tot;
    int        nabnsb;
    double     s1;
    double     s2;
    double     ab;
    double     lt_runav;
    double     lt_runav2;
} gmx_nlheur_t;

static void reset_nlistheuristics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    nlh->lt_runav  = 0;
    nlh->lt_runav2 = 0;
    nlh->step_nscheck = step;
}

static void init_nlistheuristics(gmx_nlheur_t *nlh,
                                 bool bGStatEveryStep,gmx_large_int_t step)
{
    nlh->bGStatEveryStep = bGStatEveryStep;
    nlh->nns       = 0;
    nlh->nabnsb    = 0;
    nlh->s1        = 0;
    nlh->s2        = 0;
    nlh->ab        = 0;

    reset_nlistheuristics(nlh,step);
}

static void update_nliststatistics(gmx_nlheur_t *nlh,gmx_large_int_t step)
{
    gmx_large_int_t nl_lt;
    char sbuf[22],sbuf2[22];

    /* Determine the neighbor list life time */
    nl_lt = step - nlh->step_ns;
    if (debug)
    {
        fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %s steps\n",nlh->nabnsb,gmx_step_str(nl_lt,sbuf));
    }
    nlh->nns++;
    nlh->s1 += nl_lt;
    nlh->s2 += nl_lt*nl_lt;
    nlh->ab += nlh->nabnsb;
    if (nlh->lt_runav == 0)
    {
        nlh->lt_runav  = nl_lt;
        /* Initialize the fluctuation average
         * such that at startup we check after 0 steps.
         */
        nlh->lt_runav2 = sqr(nl_lt/2.0);
    }
    /* Running average with 0.9 gives an exp. history of 9.5 */
    nlh->lt_runav2 = 0.9*nlh->lt_runav2 + 0.1*sqr(nlh->lt_runav - nl_lt);
    nlh->lt_runav  = 0.9*nlh->lt_runav  + 0.1*nl_lt;
    if (nlh->bGStatEveryStep)
    {
        /* Always check the nlist validity */
        nlh->step_nscheck = step;
    }
    else
    {
        /* We check after:  <life time> - 2*sigma
         * The factor 2 is quite conservative,
         * but we assume that with nstlist=-1 the user
         * prefers exact integration over performance.
         */
        nlh->step_nscheck = step
                  + (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)) - 1;
    }
    if (debug)
    {
        fprintf(debug,"nlist life time %s run av. %4.1f sig %3.1f check %s check with -gcom %d\n",
                gmx_step_str(nl_lt,sbuf),nlh->lt_runav,sqrt(nlh->lt_runav2),
                gmx_step_str(nlh->step_nscheck-step+1,sbuf2),
                (int)(nlh->lt_runav - 2.0*sqrt(nlh->lt_runav2)));
    }
}

static void set_nlistheuristics(gmx_nlheur_t *nlh,bool bReset,gmx_large_int_t step)
{
    int d;

    if (bReset)
    {
        reset_nlistheuristics(nlh,step);
    }
    else
    {
        update_nliststatistics(nlh,step);
    }

    nlh->step_ns = step;
    /* Initialize the cumulative coordinate scaling matrix */
    clear_mat(nlh->scale_tot);
    for(d=0; d<DIM; d++)
    {
        nlh->scale_tot[d][d] = 1.0;
    }
}

double do_md(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, bool bVerbose,bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite,gmx_constr_t constr,
             int stepout,t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb,gmx_wallcycle_t wcycle,
             gmx_edsam_t ed,t_forcerec *fr,
             int repl_ex_nst,int repl_ex_seed,
             real cpt_period,real max_hours,
             unsigned long Flags,
             gmx_runtime_t *runtime)
{
    int        fp_trn=0,fp_xtc=0;
    ener_file_t fp_ene=NULL;
    gmx_large_int_t step,step_rel;
    const char *fn_cpt;
    FILE       *fp_dhdl=NULL,*fp_field=NULL;
    double     run_time;
    double     t,t0,lam0;
    bool       bGStatEveryStep,bGStat,bNstEner,bCalcPres,bCalcEner;
    bool       bNS,bNStList,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
               bFirstStep,bStateFromTPX,bInitStep,bLastStep,bEkinAveVel,
               bBornRadii,bStartingFromCpt,bVVAveVel,bVVAveEkin;
    bool       bDoDHDL=FALSE;
    bool       bNEMD,do_ene,do_log,do_verbose,bRerunWarnNoV=TRUE,
               bForceUpdate=FALSE,bX,bV,bF,bXTC,bCPT;
    bool       bMasterState;
    int        force_flags,cglo_flags;
    tensor     force_vir,shake_vir,total_vir,tmp_vir,pres;
    int        i,m,status;
    rvec       mu_tot;
    t_vcm      *vcm;
	t_state    *bufstate;   
    matrix     *scale_tot,pcoupl_mu,M,ebox;
    gmx_nlheur_t nlh;
    t_trxframe rerun_fr;
    gmx_repl_ex_t repl_ex=NULL;
    int        nchkpt=1;
    /* Booleans (disguised as a reals) to checkpoint and terminate mdrun */  
    real       chkpt=0,terminate=0,terminate_now=0;

    gmx_localtop_t *top;	
    t_mdebin *mdebin=NULL;
    t_state    *state=NULL;
    rvec       *f_global=NULL;
    int        n_xtc=-1;
    rvec       *x_xtc=NULL;
    gmx_enerdata_t *enerd;
    rvec       *f=NULL;
    gmx_global_stat_t gstat;
    gmx_update_t upd=NULL;
    t_graph    *graph=NULL;

    bool        bFFscan;
    gmx_groups_t *groups;
    gmx_ekindata_t *ekind, *ekind_save;
    gmx_shellfc_t shellfc;
    int         count,nconverged=0;
    real        timestep=0;
    double      tcount=0;
    bool        bIonize=FALSE;
    bool        bTCR=FALSE,bConverged=TRUE,bOK,bSumEkinhOld,bExchanged;
    bool        bAppend;
    bool        bVV,bIterations,bIterate,bFirstIterate,bTemp,bPres,bTrotter;
    real        temp0,mu_aver=0,dvdl;
    int         a0,a1,gnx=0,ii;
    atom_id     *grpindex=NULL;
    char        *grpname;
    t_coupl_rec *tcr=NULL;
    rvec        *xcopy=NULL,*vcopy=NULL,*cbuf=NULL;
    matrix      boxcopy,lastbox;
	tensor      tmpvir;
	real        fom,oldfom,veta_save,pcurr,scalevir,tracevir;
	real        vetanew = 0;
    double      cycles;
    int         reset_counters=-1;
	real        last_conserved = 0;
    real        last_ekin = 0;
	int         iter_i;
	t_extmass   MassQ;
    int         **trotter_seq; 
    char        sbuf[22],sbuf2[22];
    bool        bHandledSignal=FALSE;
#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif

    /* Check for special mdrun options */
    bRerunMD = (Flags & MD_RERUN);
    bIonize  = (Flags & MD_IONIZE);
    bFFscan  = (Flags & MD_FFSCAN);
    bAppend  = (Flags & MD_APPENDFILES);
    /* md-vv uses averaged half step velocities to determine temperature unless defined otherwise by GMX_EKIN_AVE_EKIN;
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = (ir->eI==eiVV);
    if (bVV) {
        bEkinAveVel = (getenv("GMX_EKIN_AVE_EKIN")==NULL);
    } else {
        bEkinAveVel = (getenv("GMX_EKIN_AVE_VEL")!=NULL);
    }
    bVVAveVel = bVV && bEkinAveVel;
    bVVAveEkin = bVV && !bEkinAveVel;

    if (bVV) /* to store the initial velocities while computing virial */
    {
        snew(cbuf,top_global->natoms);
    }
    /* all the iteratative cases - only if there are constraints */ 
    bIterations = ((IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || (IR_NVT_TROTTER(ir))));        
    
    if (bRerunMD)
    {
        /* Since we don't know if the frames read are related in any way,
         * rebuild the neighborlist at every step.
         */
        ir->nstlist       = 1;
        ir->nstcalcenergy = 1;
        nstglobalcomm     = 1;
    }

    check_ir_old_tpx_versions(cr,fplog,ir,top_global);

    nstglobalcomm = check_nstglobalcomm(fplog,cr,nstglobalcomm,ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
    {
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
    }

    if (bRerunMD || bFFscan)
    {
        ir->nstxtcout = 0;
    }
    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog,cr,ir,oenv,&t,&t0,&state_global->lambda,&lam0,
            nrnb,top_global,&upd,
            nfile,fnm,&fp_trn,&fp_xtc,&fp_ene,&fn_cpt,
            &fp_dhdl,&fp_field,&mdebin,
            force_vir,shake_vir,mu_tot,&bNEMD,&bSimAnn,&vcm,Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd,1);
    init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f,top_global->natoms);
    }

    /* Kinetic energy data */
    snew(ekind,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind);
    /* needed for iteration of constraints */
    snew(ekind_save,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind_save);
    /* Copy the cos acceleration to the groups struct */    
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global,n_flexible_constraints(constr),
                                 (ir->bContinuation || 
                                  (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                                 NULL : state_global->x);

    if (DEFORM(*ir))
    {
#ifdef GMX_THREADS
        tMPI_Thread_mutex_lock(&box_mutex);
#endif
        set_deform_reference_box(upd,init_step_tpx,box_tpx);
#ifdef GMX_THREADS
        tMPI_Thread_mutex_unlock(&box_mutex);
#endif
    }

    {
        double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
        if ((io > 2000) && MASTER(cr))
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
    }

    if (DOMAINDECOMP(cr)) {
        top = dd_init_local_top(top_global);

        snew(state,1);
        dd_init_local_state(cr->dd,state_global,state);

        if (DDMASTER(cr->dd) && ir->nstfout) {
            snew(f_global,state_global->natoms);
        }
    } else {
        if (PAR(cr)) {
            /* Initialize the particle decomposition and split the topology */
            top = split_system(fplog,top_global,ir,cr);

            pd_cg_range(cr,&fr->cg0,&fr->hcg);
            pd_at_range(cr,&a0,&a1);
        } else {
            top = gmx_mtop_generate_local_top(top_global,ir);

            a0 = 0;
            a1 = top_global->natoms;
        }

        state = partdec_init_local_state(cr,state_global);
        f_global = f;

        atoms2md(top_global,ir,0,NULL,a0,a1-a0,mdatoms);

        if (vsite) {
            set_vsite_top(vsite,top,mdatoms,cr);
        }

        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
            graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
        }

        if (shellfc) {
            make_local_shells(cr,mdatoms,shellfc);
        }

        if (ir->pull && PAR(cr)) {
            dd_make_local_pull_groups(NULL,ir->pull,mdatoms);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog,ir->init_step,cr,TRUE,1,
                            state_global,top_global,ir,
                            state,&f,mdatoms,top,fr,
                            vsite,shellfc,constr,
                            nrnb,wcycle,FALSE);
    }

    /* If not DD, copy gb data */
    if(ir->implicit_solvent && !DOMAINDECOMP(cr))
    {
        make_local_gb(cr,fr->born,ir->gb_algorithm);
    }

    update_mdatoms(mdatoms,state->lambda);

    if (MASTER(cr))
    {
        /* Update mdebin with energy history if appending to output files */
        if ( Flags & MD_APPENDFILES )
        {
            restore_energyhistory_from_state(mdebin,&state_global->enerhist);
        }
        /* Set the initial energy history in state to zero by updating once */
        update_energyhistory(&state_global->enerhist,mdebin);
    }	

    if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG)) {
        /* Set the random state if we read a checkpoint file */
        set_stochd_state(upd,state);
    }

    /* Initialize constraints */
    if (constr) {
        if (!DOMAINDECOMP(cr))
            set_constraints(constr,top,ir,mdatoms,cr);
    }

    /* Check whether we have to GCT stuff */
    bTCR = ftp2bSet(efGCT,nfile,fnm);
    if (bTCR) {
        if (MASTER(cr)) {
            fprintf(stderr,"Will do General Coupling Theory!\n");
        }
        gnx = top_global->mols.nr;
        snew(grpindex,gnx);
        for(i=0; (i<gnx); i++) {
            grpindex[i] = i;
        }
    }

    if (repl_ex_nst > 0 && MASTER(cr))
        repl_ex = init_replica_exchange(fplog,cr->ms,state_global,ir,
                                        repl_ex_nst,repl_ex_seed);

    if (!ir->bContinuation && !bRerunMD)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++)
            {
                for(m=0; m<DIM; m++)
                {
                    if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                    {
                        state->v[i][m] = 0;
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog,constr,ir,mdatoms,state,f,
                               graph,cr,nrnb,fr,top,shake_vir);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,NULL,
                             top->idef.iparams,top->idef.il,
                             fr->ePBC,fr->bMolPBC,graph,cr,state->box);
        }
    }

    debug_gmx();
  
    /* would need to modify if vv starts with full time step */
    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                    constr,&chkpt,&terminate,&terminate_now,&(nlh.nabnsb),state->box,
                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                    (CGLO_FIRSTCALL |
                     (bVV ? CGLO_PRESSURE:0) |
                     (CGLO_TEMPERATURE) |  
                     (bRerunMD ? CGLO_RERUNMD:0) | 
                     (bEkinAveVel ? CGLO_EKINAVEVEL:0) | 
                     (CGLO_GSTAT) | 
                     ((Flags & MD_READ_EKIN) ? CGLO_READEKIN:0)));
    /* I'm assuming we need global communication the first time */
    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT)) 
    {
        for(i=0; (i<ir->opts.ngtc); i++) 
        {
            copy_mat(ekind->tcstat[i].ekinh,ekind->tcstat[i].ekinh_old);
        } 
    }
    temp0 = enerd->term[F_TEMP];
    
    /* Initiate data for the special cases */
    if (bIterations) 
    {
        bufstate = init_bufstate(state->natoms,(ir->opts.ngtc+1)*NNHCHAIN); /* extra state for barostat */
    }
    
    if (bFFscan) 
    {
        snew(xcopy,state->natoms);
        snew(vcopy,state->natoms);
        copy_rvecn(state->x,xcopy,0,state->natoms);
        copy_rvecn(state->v,vcopy,0,state->natoms);
        copy_mat(state->box,boxcopy);
    } 
    
    /* need to make an initial call to get the Trotter variables set. */
    trotter_seq = init_trotter(ir,state,&MassQ,bTrotter);
    
    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr,FALSE));
        }
        fprintf(fplog,"Initial temperature: %g K\n",enerd->term[F_TEMP]);
        if (bRerunMD)
        {
            fprintf(stderr,"starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name),opt2fn("-rerun",nfile,fnm));
            if (bVerbose)
            {
                fprintf(stderr,"Calculated time to finish depends on nsteps from "
                        "run input file,\nwhich may not correspond to the time "
                        "needed to process input trajectory.\n\n");
            }
        }
        else
        {
            fprintf(stderr,"starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->init_step > 0)
            {
                fprintf(stderr,"%s steps, %8.1f ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps,sbuf),
                        (ir->init_step+ir->nsteps)*ir->delta_t,
                        gmx_step_str(ir->init_step,sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr,"%s steps, %8.1f ps.\n",
                        gmx_step_str(ir->nsteps,sbuf),ir->nsteps*ir->delta_t);
            }
        }
        fprintf(fplog,"\n");
    }

    /* Set and write start time */
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Started mdrun",runtime);
    wallcycle_start(wcycle,ewcRUN);
    if (fplog)
        fprintf(fplog,"\n");

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
#ifdef GMX_FAHCORE
    chkpt_ret=fcCheckPointParallel( cr->nodeid,
                                    NULL,0);
    if ( chkpt_ret == 0 ) 
        gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", 0 );
#endif

    debug_gmx();
    /***********************************************************
     *
     *             Loop over MD steps 
     *
     ************************************************************/

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD)
    {
        if (getenv("GMX_FORCE_UPDATE"))
        {
            bForceUpdate = TRUE;
        }

        bNotLastFrame = read_first_frame(oenv,&status,
                                         opt2fn("-rerun",nfile,fnm),
                                         &rerun_fr,TRX_NEED_X | TRX_READ_V);
        if (rerun_fr.natoms != top_global->natoms)
        {
            gmx_fatal(FARGS,
                      "Number of atoms in trajectory (%d) does not match the "
                      "run input file (%d)\n",
                      rerun_fr.natoms,top_global->natoms);
        }
        if (ir->ePBC != epbcNONE)
        {
            if (!rerun_fr.bBox)
            {
                gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used",rerun_fr.step,rerun_fr.time);
            }
            if (max_cutoff2(ir->ePBC,rerun_fr.box) < sqr(fr->rlistlong))
            {
                gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions",rerun_fr.step,rerun_fr.time);
            }

            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box,fr->shift_vec);
        }
    }

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
    bInitStep = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep = FALSE;
    bSumEkinhOld = FALSE;
    bExchanged = FALSE;

    step = ir->init_step;
    step_rel = 0;

    if (ir->nstlist == -1)
    {
        init_nlistheuristics(&nlh,bGStatEveryStep,step);
    }

    bLastStep = (bRerunMD || step_rel > ir->nsteps);
    while (!bLastStep || (bRerunMD && bNotLastFrame)) {

        wallcycle_start(wcycle,ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        if (bRerunMD) {
            if (rerun_fr.bStep) {
                step = rerun_fr.step;
                step_rel = step - ir->init_step;
            }
            if (rerun_fr.bTime) {
                t = rerun_fr.time;
            }
            else
            {
                t = step;
            }
        } 
        else 
        {
            bLastStep = (step_rel == ir->nsteps);
            t = t0 + step*ir->delta_t;
        }

        if (ir->efep != efepNO)
        {
            if (bRerunMD && rerun_fr.bLambda && (ir->delta_lambda!=0))
            {
                state_global->lambda = rerun_fr.lambda;
            }
            else
            {
                state_global->lambda = lam0 + step*ir->delta_lambda;
            }
            state->lambda = state_global->lambda;
            bDoDHDL = do_per_step(step,ir->nstdhdl);
        }

        if (bSimAnn) 
        {
            update_annealing_target_temp(&(ir->opts),t);
        }

        if (bRerunMD)
        {
            if (!(DOMAINDECOMP(cr) && !MASTER(cr)))
            {
                for(i=0; i<state_global->natoms; i++)
                {
                    copy_rvec(rerun_fr.x[i],state_global->x[i]);
                }
                if (rerun_fr.bV)
                {
                    for(i=0; i<state_global->natoms; i++)
                    {
                        copy_rvec(rerun_fr.v[i],state_global->v[i]);
                    }
                }
                else
                {
                    for(i=0; i<state_global->natoms; i++)
                    {
                        clear_rvec(state_global->v[i]);
                    }
                    if (bRerunWarnNoV)
                    {
                        fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
                                "         Ekin, temperature and pressure are incorrect,\n"
                                "         the virial will be incorrect when constraints are present.\n"
                                "\n");
                        bRerunWarnNoV = FALSE;
                    }
                }
            }
            copy_mat(rerun_fr.box,state_global->box);
            copy_mat(state_global->box,state->box);

            if (vsite && (Flags & MD_RERUN_VSITE))
            {
                if (DOMAINDECOMP(cr))
                {
                    gmx_fatal(FARGS,"Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
                }
                if (graph)
                {
                    /* Following is necessary because the graph may get out of sync
                     * with the coordinates if we only have every N'th coordinate set
                     */
                    mk_mshift(fplog,graph,fr->ePBC,state->box,state->x);
                    shift_self(graph,state->box,state->x);
                }
                construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                                 top->idef.iparams,top->idef.il,
                                 fr->ePBC,fr->bMolPBC,graph,cr,state->box);
                if (graph)
                {
                    unshift_self(graph,state->box,state->x);
                }
            }
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step,ir->nstcomm));

        /* Copy back starting coordinates in case we're doing a forcefield scan */
        if (bFFscan)
        {
            for(ii=0; (ii<state->natoms); ii++)
            {
                copy_rvec(xcopy[ii],state->x[ii]);
                copy_rvec(vcopy[ii],state->v[ii]);
            }
            copy_mat(boxcopy,state->box);
        }

        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            bNS = (bFirstStep || ir->nstlist != 0);
            bNStList = bNS;
        }
        else
        {
            /* Determine whether or not to do Neighbour Searching and LR */
            bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);
            
            bNS = (bFirstStep || bExchanged || bNStList ||
                   (ir->nstlist == -1 && nlh.nabnsb > 0));

            if (bNS && ir->nstlist == -1)
            {
                set_nlistheuristics(&nlh,bFirstStep || bExchanged,step);
            }
        } 

        if (terminate_now > 0 || (terminate_now < 0 && bNS))
        {
            bLastStep = TRUE;
        }

        /* Determine whether or not to update the Born radii if doing GB */
        bBornRadii=bFirstStep;
        if (ir->implicit_solvent && (step % ir->nstgbradii==0))
        {
            bBornRadii=TRUE;
        }
        
        do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
                  (step % stepout == 0 || bFirstStep || bLastStep);

        if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD))
        {
            if (bRerunMD)
            {
                bMasterState = TRUE;
            }
            else
            {
                bMasterState = FALSE;
                /* Correct the new box if it is too skewed */
                if (DYNAMIC_BOX(*ir))
                {
                    if (correct_box(fplog,step,state->box,graph))
                    {
                        bMasterState = TRUE;
                    }
                }
                if (DOMAINDECOMP(cr) && bMasterState)
                {
                    dd_collect_state(cr->dd,state,state_global);
                }
            }

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                wallcycle_start(wcycle,ewcDOMDEC);
                dd_partition_system(fplog,step,cr,
                                    bMasterState,nstglobalcomm,
                                    state_global,top_global,ir,
                                    state,&f,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,do_verbose);
                wallcycle_stop(wcycle,ewcDOMDEC);
            }
        }

        if (MASTER(cr) && do_log && !bFFscan)
        {
            print_ebin_header(fplog,step,t,state->lambda);
        }

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms,state->lambda); 
        }

        if (bRerunMD && rerun_fr.bV)
        {
            
            /* We need the kinetic energy at minus the half step for determining
             * the full step kinetic energy and possibly for T-coupling.*/
            /* This may not be quite working correctly yet . . . . */
            compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                            wcycle,enerd,NULL,NULL,NULL,NULL,mu_tot,
                            constr,&chkpt,&terminate,&terminate_now,&(nlh.nabnsb),state->box,
                            top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                            (CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE) |
                            (bEkinAveVel?CGLO_EKINAVEVEL:0));
        }
        clear_mat(force_vir);
        
        /* Ionize the atoms if necessary */
        if (bIonize)
        {
            ionize(fplog,oenv,mdatoms,top_global,t,ir,state->x,state->v,
                   mdatoms->start,mdatoms->start+mdatoms->homenr,state->box,cr);
        }
        
        /* Update force field in ffscan program */
        if (bFFscan)
        {
            if (update_forcefield(fplog,
                                  nfile,fnm,fr,
                                  mdatoms->nr,state->x,state->box)) {
                if (gmx_parallel_env())
                {
                    gmx_finalize();
                }
                exit(0);
            }
        }

        GMX_MPE_LOG(ev_timestep2);

        if ((bNS || bLastStep) && (step > ir->init_step) && !bRerunMD)
        {
            bCPT = (chkpt > 0 || (bLastStep && (Flags & MD_CONFOUT)));
            if (bCPT)
            {
                chkpt = 0;
            }
        }
        else
        {
            bCPT = FALSE;
        }

        /* Determine the pressure:
         * always when we want exact averages in the energy file,
         * at ns steps when we have pressure coupling,
         * otherwise only at energy output steps (set below).
         */
        bNstEner = (bGStatEveryStep || do_per_step(step,ir->nstcalcenergy));
        bCalcEner = bNstEner;
		bCalcPres = (bGStatEveryStep || (ir->epc != epcNO && bNS));

        /* Do we need global communication ? */
        bGStat = (bCalcEner || bStopCM ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcPres = TRUE;
            bCalcEner = TRUE;
            bGStat    = TRUE;
        }
        
        /* these CGLO_ options remain the same throughout the iteration */
        cglo_flags = ((bRerunMD ? CGLO_RERUNMD : 0) |
                      (bEkinAveVel ? CGLO_EKINAVEVEL : 0) |
                      (bStopCM ? CGLO_STOPCM : 0) |
                      (bGStat ? CGLO_GSTAT : 0) |
                      (bNEMD ? CGLO_NEMD : 0) |
                      (constr ? CGLO_CONSTRAINT : 0)
            );
        
        force_flags = (GMX_FORCE_STATECHANGED |
                       ((DYNAMIC_BOX(*ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       (bNStList ? GMX_FORCE_DOLR : 0) |
                       GMX_FORCE_SEPLRF |
                       (bCalcEner ? GMX_FORCE_VIRIAL : 0) |
                       (bDoDHDL ? GMX_FORCE_DHDL : 0));
        
        if (shellfc)
        {
            /* Now is the time to relax the shells */
            count=relax_shell_flexcon(fplog,cr,bVerbose,bFFscan ? step+1 : step,
                                      ir,bNS,force_flags,
                                      bStopCM,top,top_global,
                                      constr,enerd,fcd,
                                      state,f,force_vir,mdatoms,
                                      nrnb,wcycle,graph,groups,
                                      shellfc,fr,bBornRadii,t,mu_tot,
                                      state->natoms,&bConverged,vsite,
                                      fp_field);
            tcount+=count;

            if (bConverged)
            {
                nconverged++;
            }
        }
        else
        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too. 
             * Check comments in sim_util.c
             */

            do_force(fplog,cr,ir,step,nrnb,wcycle,top,top_global,groups,
                     state->box,state->x,&state->hist,
                     f,force_vir,mdatoms,enerd,fcd,
                     state->lambda,graph,
                     fr,vsite,mu_tot,t,fp_field,ed,bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);
        }

        GMX_BARRIER(cr->mpi_comm_mygroup);

        if (bTCR)
        {
            mu_aver = calc_mu_aver(cr,state->x,mdatoms->chargeA,
                                   mu_tot,&top_global->mols,mdatoms,gnx,grpindex);
        }

        if (bTCR && bFirstStep)
        {
            tcr=init_coupling(fplog,nfile,fnm,cr,fr,mdatoms,&(top->idef));
            fprintf(fplog,"Done init_coupling\n"); 
            fflush(fplog);
        }
        
        /*  ############### START FIRST UPDATE HALF-STEP ############### */

        if (!bStartingFromCpt && !bRerunMD) {
            if (bVVAveVel && bInitStep) {
                /* if using velocity verlet with full time step Ekin, take the first half step only to compute the 
                   virial for the first step. From there, revert back to the initial coordinates
                   so that the input is actually the initial step */
                copy_rvecn(state->v,cbuf,0,state->natoms); /* should make this better for parallelizing? */
            }

            /* this is for NHC in the Ekin(t+dt/2) version of vv */
            if (!bInitStep || !bEkinAveVel)
            {
                trotter_update(ir,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq[1],bEkinAveVel);            
            }
            
            update_coords(fplog,step,ir,mdatoms,state,
                          f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                          ekind,M,wcycle,upd,bInitStep,etrtVELOCITY,cr,nrnb,constr,&top->idef);
            
            bIterate = bIterations && (!bInitStep);
            bFirstIterate = TRUE;
            /* for iterations, we save these vectors, as we will be self-consistently iterating
               the calculations */
            /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */
            
            /* save the state */
            if (bIterate) { 
                copy_coupling_state(state,bufstate,ekind,ekind_save);
            }
            
            while (bIterate || bFirstIterate) 
            {
                if (bIterate) 
                {
                    copy_coupling_state(bufstate,state,ekind_save,ekind);
                }
                if (bIterate) 
                {
                    if (bFirstIterate && bTrotter) 
                    {
                        /* The first time, we need a decent first estimate
                           of veta(t+dt) to compute the constraints.  Do
                           this by computing the box volume part of the
                           trotter integration at this time. Nothing else
                           should be changed by this routine here.  If
                           !(first time), we start with the previous value
                           of veta.  */
                        
                        veta_save = state->veta;
                        trotter_update(ir,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq[0],FALSE);
                        vetanew = state->veta;
                        state->veta = veta_save;
                    } 
                } 
                
                bOK = TRUE;
                if ( !bRerunMD || rerun_fr.bV || bForceUpdate) {  /* Why is rerun_fr.bV here?  Unclear. */
                    wallcycle_start(wcycle,ewcUPDATE);
                    dvdl = 0;
                    
                    update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                       &top->idef,shake_vir,NULL,
                                       cr,nrnb,wcycle,upd,constr,
                                       bInitStep,TRUE,bCalcPres,vetanew);
                    
                    if (!bOK && !bFFscan)
                    {
                        gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
                    }
                    
                } 
                else if (graph)
                { /* Need to unshift here if a do_force has been
                     called in the previous step */
                    unshift_self(graph,state->box,state->x);
                }
                
                
                if (bVV) {
                    /* if EkinAveVel, compute the pressure and constraints */
                    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                                    constr,&chkpt,&terminate,&terminate_now,&(nlh.nabnsb),state->box,
                                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                                    (cglo_flags | CGLO_ENERGY) | 
                                    (bFirstIterate ? (cglo_flags | CGLO_PRESSURE):0) |
                                    (((bEkinAveVel &&(!bInitStep)) || (!bEkinAveVel && IR_NPT_TROTTER(ir)))? (cglo_flags | CGLO_TEMPERATURE):0) | 
                                    ((bEkinAveVel || (!bEkinAveVel && IR_NPT_TROTTER(ir))) ? CGLO_EKINAVEVEL : 0) |
                                    (bIterate ? CGLO_ITERATE : 0) | 
                                    (bFirstIterate ? CGLO_FIRSTITERATE : 0) | 
                                    (cglo_flags | CGLO_SCALEEKIN)
                        );
                }
                /* explanation of above: 
                   a) We compute Ekin at the full time step
                   if 1) we are using the AveVel Ekin, and it's not the
                   initial step, or 2) if we are using AveEkin, but need the full
                   time step kinetic energy for the pressure.
                   b) If we are using EkinAveEkin for the kinetic energy for the temperture control, we still feed in 
                   EkinAveVel because it's needed for the pressure */
                
                /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
                if (!bInitStep) 
                {
                    trotter_update(ir,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq[2],bEkinAveVel);
                }
                
                if (done_iterating(cr,fplog,&bFirstIterate,&bIterate,state->veta,&vetanew,0)) break;
            }
            
            if (bTrotter && !bInitStep) {
                copy_mat(shake_vir,state->vir_prev);
                if (IR_NVT_TROTTER(ir) && bEkinAveVel) {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts),ekind,NULL,bEkinAveVel,FALSE,FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (bInitStep && bVVAveVel) {
                copy_rvecn(cbuf,state->v,0,state->natoms);
            }
            
            if (fr->bSepDVDL && fplog && do_log) 
            {
                fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
            }
            enerd->term[F_DHDL_CON] += dvdl;
            
            wallcycle_stop(wcycle,ewcUPDATE);
            GMX_MPE_LOG(ev_timestep1);
            
            /* MRS -- done iterating -- compute the conserved quantity */
        }
        if (bVV) {
            last_conserved = 0;
            if (IR_NVT_TROTTER(ir) || IR_NPT_TROTTER(ir))
            {
                last_conserved = 
                    NPT_energy(ir,state->nosehoover_xi,state->nosehoover_vxi,state->veta, state->box, &MassQ); 
                if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres)) 
                {
                    last_conserved -= enerd->term[F_DISPCORR];
                }
            }
            if (bEkinAveVel) {
                last_ekin = enerd->term[F_EKIN]; /* does this get preserved through checkpointing? */
            }
        }
        
        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */
        
        /* ################## START TRAJECTORY OUTPUT ################# */
        
        /* Now we have the energies and forces corresponding to the 
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        bX   = do_per_step(step,ir->nstxout);
        bV   = do_per_step(step,ir->nstvout);
        bF   = do_per_step(step,ir->nstfout);
        bXTC = do_per_step(step,ir->nstxtcout);

#ifdef GMX_FAHCORE
        if (MASTER(cr))
            fcReportProgress( ir->nsteps, step );

        bX = bX || bLastStep; /*enforce writing positions and velocities 
                                at end of run */
        bV = bV || bLastStep;
        {
            int nthreads=(cr->nthreads==0 ? 1 : cr->nthreads);
            int nnodes=(cr->nnodes==0 ? 1 : cr->nnodes);

            bCPT = bCPT;
            /*Gromacs drives checkpointing; no ||  
              fcCheckPointPendingThreads(cr->nodeid,
              nthreads*nnodes);*/
            /* sync bCPT and fc record-keeping */
            if (bCPT && MASTER(cr))
                fcRequestCheckPoint();
        }
#endif

        
        if (bX || bV || bF || bXTC || bCPT)
        {
            wallcycle_start(wcycle,ewcTRAJ);
            if (bCPT)
            {
                if (state->flags & (1<<estLD_RNG))
                {
                    get_stochd_state(upd,state);
                }
                if (MASTER(cr))
                {
                    if (bSumEkinhOld)
                    {
                        state_global->ekinstate.bUpToDate = FALSE;
                    }
                    else
                    {
                        update_ekinstate(&state_global->ekinstate,ekind);
                        state_global->ekinstate.bUpToDate = TRUE;
                    }
                    update_energyhistory(&state_global->enerhist,mdebin);
                }
            }
            write_traj(fplog,cr,fp_trn,bX,bV,bF,fp_xtc,bXTC,ir->xtcprec,
                       fn_cpt,bCPT,top_global,ir->eI,ir->simulation_part,
                       step,t,state,state_global,f,f_global,&n_xtc,&x_xtc);
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                (Flags & MD_CONFOUT) && MASTER(cr) &&
                !bRerunMD && !bFFscan)
            {
                /* x and v have been collected in write_traj */
                fprintf(stderr,"\nWriting final coordinates.\n");
                if (ir->ePBC != epbcNONE && !ir->bPeriodicMols &&
                    DOMAINDECOMP(cr))
                {
                    /* Make molecules whole only for confout writing */
                    do_pbc_mtop(fplog,ir->ePBC,state->box,top_global,state_global->x);
                }
                write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
                                    *top_global->name,top_global,
                                    state_global->x,state_global->v,
                                    ir->ePBC,state->box);
                debug_gmx();
            }
            wallcycle_stop(wcycle,ewcTRAJ);
        }
        GMX_MPE_LOG(ev_output_finish);

        /* kludge -- virial is lost with restart for NPT control. Must restart */
        if (bStartingFromCpt && bVV) 
        {
            copy_mat(state->vir_prev,shake_vir);
        }
        /*  ################## END TRAJECTORY OUTPUT ################ */
        
        /* Determine the pressure:                                                             
         * always when we want exact averages in the energy file,                              
         * at ns steps when we have pressure coupling,                                         
         * otherwise only at energy output steps (set below).                                  
         */
        
        bNstEner = (bGStatEveryStep || do_per_step(step,ir->nstcalcenergy));
        bCalcEner = bNstEner;
        bCalcPres = (bGStatEveryStep || (ir->epc != epcNO && bNS));
        
        /* Do we need global communication ? */
        bGStat = (bGStatEveryStep || bStopCM || bNS ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));
        
        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);
        
        if (do_ene || do_log)
        {
            bCalcPres = TRUE;
            bGStat    = TRUE;
        }
        
        bIterate = bIterations;
        bFirstIterate = TRUE;
        
        /* for iterations, we save these vectors, as we will be redoing the calculations */
        if (bIterate) 
        {
            copy_coupling_state(state,bufstate,ekind,ekind_save);
        }
        
        while (bIterate || bFirstIterate) 
        {
            /* We now restore these vectors to redo the calculation with improved extended variables */    
            if (bIterate) 
            { 
                copy_coupling_state(bufstate,state,ekind_save,ekind);
            }
            
            /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
               so scroll down for that logic */
            
            /* #########   START SECOND UPDATE STEP ################# */
            
            GMX_MPE_LOG(ev_update_start);
            bOK = TRUE;
            if (!bRerunMD || rerun_fr.bV || bForceUpdate) 
            {
                wallcycle_start(wcycle,ewcUPDATE);
                dvdl = 0;
                
                
                /* Box is changed in update() when we do pressure coupling,
                 * but we should still use the old box for energy corrections and when
                 * writing it to the energy file, so it matches the trajectory files for
                 * the same timestep above. Make a copy in a separate array.
                 */
                
                copy_mat(state->box,lastbox);
                
                /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS*/
                if (bTrotter) 
                {
                    if (bFirstIterate) 
                    {
                        tracevir = trace(shake_vir);
                    }
                    if (bIterate) 
                    {
                        /* we use a new value of scalevir to converge the iterations faster */
                        scalevir = tracevir/trace(shake_vir);
                        msmul(shake_vir,scalevir,shake_vir); 
                        m_add(force_vir,shake_vir,total_vir);
                        clear_mat(shake_vir);
                    }
                    trotter_update(ir,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq[3],bEkinAveVel);
                }
                /* We can only do Berendsen coupling after we have summed the kinetic
                 * energy or virial. Since the happens in global_state after update,
                 * we should only do it at step % nstlist = 1 with bGStatEveryStep=FALSE.
                 */
                
                update_extended(fplog,step,ir,mdatoms,state,ekind,pcoupl_mu,M,wcycle,
                                upd,bInitStep,FALSE,&MassQ);
                
                /* velocity (for VV) */
                update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                              ekind,M,wcycle,upd,FALSE,etrtVELOCITY,cr,nrnb,constr,&top->idef);

                /* above, initialize just copies ekinh into ekin, it doesn't copy 
                   position (for VV), and entire integrator for MD */
                
                if (bVVAveEkin) 
                {
                    copy_rvecn(state->x,cbuf,0,state->natoms);
                }

                update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                              ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
                
                update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                   &top->idef,shake_vir,force_vir,
                                   cr,nrnb,wcycle,upd,constr,
                                   bInitStep,FALSE,bCalcPres,state->veta);  
                
                if (bVVAveEkin)
                {
                    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                                    constr,&chkpt,&terminate,&terminate_now,&(nlh.nabnsb),lastbox,
                                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                                    (cglo_flags | CGLO_TEMPERATURE) 

                        );
                    trotter_update(ir,ekind,enerd,state,total_vir,mdatoms,&MassQ,trotter_seq[4],bEkinAveVel);            
                    copy_rvecn(cbuf,state->x,0,state->natoms);
                    update_coords(fplog,step,ir,mdatoms,state,f,fr->bTwinRange && bNStList,fr->f_twin,fcd,
                                  ekind,M,wcycle,upd,bInitStep,etrtPOSITION,cr,nrnb,constr,&top->idef);
                    /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                    /* are the small terms in the shake_vir here due
                     * to numerical errors, or are they important
                     * physically? I'm thinking they are just errors, but not completely sure. 
                     * For now, will call without actually constraining, constr=NULL*/
 
                    update_constraints(fplog,step,&dvdl,ir,ekind,mdatoms,state,graph,f,
                                       &top->idef,tmp_vir,force_vir,
                                       cr,nrnb,wcycle,upd,NULL,
                                       bInitStep,FALSE,bCalcPres,state->veta);  
                }
                
                if (!bOK && !bFFscan) 
                {
                    gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
                }
                
                if (fr->bSepDVDL && fplog && do_log) 
                {
                    fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
                }
                enerd->term[F_DHDL_CON] += dvdl;
                wallcycle_stop(wcycle,ewcUPDATE);
            } 
            else if (graph) 
            {
                /* Need to unshift here */
                unshift_self(graph,state->box,state->x);
            }
            
            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_update_finish);
            
            if (vsite != NULL) 
            {
                wallcycle_start(wcycle,ewcVSITECONSTR);
                if (graph != NULL) 
                {
                    shift_self(graph,state->box,state->x);
                }
                construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                                 top->idef.iparams,top->idef.il,
                                 fr->ePBC,fr->bMolPBC,graph,cr,state->box);
                
                if (graph != NULL) 
                {
                    unshift_self(graph,state->box,state->x);
                }
                wallcycle_stop(wcycle,ewcVSITECONSTR);
            }
            
            /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints ############ */
            compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
                            wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
                            constr,&chkpt,&terminate,&terminate_now,&(nlh.nabnsb),lastbox,
                            top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
                            (!bVV ? (cglo_flags | CGLO_ENERGY):0) | 
                            (!bVV ? (cglo_flags | CGLO_PRESSURE):0) |
                            (!bVV ? (cglo_flags | CGLO_TEMPERATURE):0) | 
                            (bIterate ? (cglo_flags | CGLO_ITERATE) : 0) |
                            (bFirstIterate ? (cglo_flags | CGLO_FIRSTITERATE) : 0)
                );            
            /* bIterate is set to keep it from eliminating the old ekinh */
            /* #############  END CALC EKIN AND PRESSURE ################# */
            
            if (done_iterating(cr,fplog,&bFirstIterate,&bIterate,trace(shake_vir),&tracevir,1)) break;
        }
        
        update_box(fplog,step,ir,mdatoms,state,graph,f,
                   ir->nstlist==-1 ? &nlh.scale_tot : NULL,pcoupl_mu,nrnb,wcycle,upd,bInitStep,FALSE);

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */
        
        /* Determine the wallclock run time up till now */
        run_time = (double)time(NULL) - (double)runtime->real;

        /* Check whether everything is still allright */    
        if ((bGotTermSignal || bGotUsr1Signal) && !bHandledSignal)
        {
            if (bGotTermSignal || ir->nstlist == 0)
            {
                terminate = 1;
            }
            else
            {
                terminate = -1;
            }
            if (!PAR(cr))
            {
                terminate_now = terminate;
            }
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        bGotTermSignal ? "TERM" : "USR1",
                        terminate==-1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    bGotTermSignal ? "TERM" : "USR1",
                    terminate==-1 ? "NS " : "");
            fflush(stderr);
            bHandledSignal=TRUE;
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 terminate == 0)
        {
            /* Signal to terminate the run */
            terminate = (ir->nstlist == 0 ? 1 : -1);
            if (!PAR(cr))
            {
                terminate_now = terminate;
            }
            if (fplog)
            {
                fprintf(fplog,"\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
        }

        if (ir->nstlist == -1 && !bRerunMD)
        {
            /* When bGStatEveryStep=FALSE, global_stat is only called
             * when we check the atom displacements, not at NS steps.
             * This means that also the bonded interaction count check is not
             * performed immediately after NS. Therefore a few MD steps could
             * be performed with missing interactions.
             * But wrong energies are never written to file,
             * since energies are only written after global_stat
             * has been called.
             */
            if (step >= nlh.step_nscheck)
            {
                nlh.nabnsb = natoms_beyond_ns_buffer(ir,fr,&top->cgs,
                                                     nlh.scale_tot,state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nlh.nabnsb = 0;
            }
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 || 
                            run_time >= nchkpt*cpt_period*60.0)))
        {
            if (chkpt == 0)
            {
                nchkpt++;
            }
            chkpt = 1;
        }
        
        /* The coordinates (x) were unshifted in update */
        if (bFFscan && (shellfc==NULL || bConverged))
        {
            if (print_forcefield(fplog,enerd->term,mdatoms->homenr,
                                 f,NULL,xcopy,
                                 &(top_global->mols),mdatoms->massT,pres))
            {
                if (gmx_parallel_env()) {
                    gmx_finalize();
                }
                fprintf(stderr,"\n");
                exit(0);
            }
        }
        if (!bGStat)
        {
            /* We will not sum ekinh_old,                                                            
             * so signal that we still have to do it.                                                
             */
            bSumEkinhOld = TRUE;
        }

        if (bTCR)
        {
            /* Only do GCT when the relaxation of shells (minimization) has converged,
             * otherwise we might be coupling to bogus energies. 
             * In parallel we must always do this, because the other sims might
             * update the FF.
             */

            /* Since this is called with the new coordinates state->x, I assume
             * we want the new box state->box too. / EL 20040121
             */
            do_coupling(fplog,oenv,nfile,fnm,tcr,t,step,enerd->term,fr,
                        ir,MASTER(cr),
                        mdatoms,&(top->idef),mu_aver,
                        top_global->mols.nr,cr,
                        state->box,total_vir,pres,
                        mu_tot,state->x,f,bConverged);
            debug_gmx();
        }

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */
        
        sum_dhdl(enerd,state->lambda,ir);
        /* use the directly determined last velocity, not actually the averaged half steps */
        if (bTrotter && bEkinAveVel) {
            enerd->term[F_EKIN] = last_ekin;
        }
        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];
        
        switch (ir->etc) 
        {
        case etcNO:
            break;
        case etcBERENDSEN:
            break;
        case etcNOSEHOOVER:
            if (IR_NVT_TROTTER(ir)) {
                enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + last_conserved;
            } else {
                enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + 
                    NPT_energy(ir,state->nosehoover_xi,
                               state->nosehoover_vxi,state->veta,state->box,&MassQ);	
            }
            break;
        case etcVRESCALE:
            enerd->term[F_ECONSERVED] =
                enerd->term[F_ETOT] + vrescale_energy(&(ir->opts),
                                                      state->therm_integral);
            break;
        default:
            break;
        }
        
        /* Check for excessively large energies */
        if (bIonize) 
        {
#ifdef GMX_DOUBLE
            real etot_max = 1e200;
#else
            real etot_max = 1e30;
#endif
            if (fabs(enerd->term[F_ETOT]) > etot_max) 
            {
                fprintf(stderr,"Energy too large (%g), giving up\n",
                        enerd->term[F_ETOT]);
            }
        }
        /* #########  END PREPARING EDR OUTPUT  ###########  */
        
        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep) 
        {
            runtime_upd_proc(runtime);
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            bool do_dr,do_or;

            if (bNstEner)
            {
                upd_mdebin(mdebin,bDoDHDL ? fp_dhdl : NULL,TRUE,
                           t,mdatoms->tmass,enerd,state,lastbox,
                           shake_vir,force_vir,total_vir,pres,
                           ekind,mu_tot,constr);
            }
            else
            {
                upd_mdebin_step(mdebin);
            }

            do_dr  = do_per_step(step,ir->nstdisreout);
            do_or  = do_per_step(step,ir->nstorireout);

            print_ebin(fp_ene,do_ene,do_dr,do_or,do_log?fplog:NULL,step,t,
                       eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));

            if (ir->ePull != epullNO)
            {
                pull_print_output(ir->pull,step,t);
            }

            if (do_per_step(step,ir->nstlog))
            {
                if(fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of quota?");
                }
            }
        }


        /* Remaining runtime */
        if (MULTIMASTER(cr) && do_verbose) 
        {
            if (shellfc) 
            {
                fprintf(stderr,"\n");
            }
            print_time(stderr,runtime,step,ir);
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step,repl_ex_nst)) 
        {
            bExchanged = replica_exchange(fplog,cr,repl_ex,
                                          state_global,enerd->term,
                                          state,step,t);
        }
        if (bExchanged && PAR(cr)) 
        {
            if (DOMAINDECOMP(cr)) 
            {
                dd_partition_system(fplog,step,cr,TRUE,1,
                                    state_global,top_global,ir,
                                    state,&f,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,FALSE);
            } 
            else 
            {
                bcast_state(cr,state,FALSE);
            }
        }

        bFirstStep = FALSE;
        bInitStep = FALSE;
        bStartingFromCpt = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* Complicated conditional when bGStatEveryStep=FALSE.
         * We can not just use bGStat, since then the simulation results
         * would depend on nstenergy and nstlog or step_nscheck.
         */
        if (((state->flags & (1<<estPRES_PREV)) || (state->flags & (1<<estVIR_PREV))) &&
            (bGStatEveryStep ||
             (ir->nstlist > 0 && step % ir->nstlist == 0) ||
             (ir->nstlist < 0 && nlh.nabnsb > 0) ||
             (ir->nstlist == 0 && bGStat))) 
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            if (state->flags & (1<<estPRES_PREV))
            {
                copy_mat(pres,state->pres_prev);
            }
        }

    /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        if (bRerunMD) 
        {
            /* read next frame from input trajectory */
            bNotLastFrame = read_next_frame(oenv,status,&rerun_fr);
        }

        if (!bRerunMD || !rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }

        cycles = wallcycle_stop(wcycle,ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd,cycles,ddCyclStep);
        }

        if (step_rel == wcycle_get_reset_counters(wcycle))
        {
            /* Reset all the counters related to performance over the run */
            reset_all_counters(fplog,cr,step,&step_rel,ir,wcycle,nrnb,runtime);
            wcycle_set_reset_counters(wcycle, 0);
        }
    }
    /* End of main MD loop */
    debug_gmx();

    /* Stop the time */
    runtime_end(runtime);

    if (bRerunMD)
    {
        close_trj(status);
    }

    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_finish(cr);
    }

    if (MASTER(cr))
    {
        if (bGStatEveryStep && !bRerunMD) 
        {
            print_ebin(fp_ene,FALSE,FALSE,FALSE,fplog,step,t,
                       eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
        }
        close_enx(fp_ene);
        if (ir->nstxtcout)
        {
            close_xtc(fp_xtc);
        }
        close_trn(fp_trn);
        if (fp_dhdl)
        {
            gmx_fio_fclose(fp_dhdl);
        }
        if (fp_field)
        {
            gmx_fio_fclose(fp_field);
        }
    }
    debug_gmx();

    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
    {
        fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",nlh.s1/nlh.nns,sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
        fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",nlh.ab/nlh.nns);
    }

    if (shellfc && fplog)
    {
        fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
                (nconverged*100.0)/step_rel);
        fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
                tcount/step_rel);
    }

    if (repl_ex_nst > 0 && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog,repl_ex);
    }

    runtime->nsteps_done = step_rel;

    return 0;
}
        
