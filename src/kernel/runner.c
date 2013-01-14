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

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _isnan() */
#include <float.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "statutil.h"
#include "mdrun.h"
#include "network.h"
#include "pull.h"
#include "names.h"
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
#include "coulomb.h"
#include "constr.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sighandler.h"
#include "tpxio.h"
#include "txtdump.h"

#include "md_openmm.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

#ifdef GMX_OPENMM
#include "md_openmm.h"
#endif


typedef struct { 
    gmx_integrator_t *func;
} gmx_intp_t;

/* The array should match the eI array in include/types/enums.h */
#ifdef GMX_OPENMM  /* FIXME do_md_openmm needs fixing */
const gmx_intp_t integrator[eiNR] = { {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm}, {do_md_openmm},{do_md_openmm}};
#else
const gmx_intp_t integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md}, {do_md},{do_md}};
#endif

gmx_large_int_t     deform_init_init_step_tpx;
matrix              deform_init_box_tpx;
#ifdef GMX_THREADS
tMPI_Thread_mutex_t deform_init_box_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


#ifdef GMX_THREADS
struct mdrunner_arglist
{
    FILE *fplog;
    t_commrec *cr;
    int nfile;
    const t_filenm *fnm;
    output_env_t oenv;
    gmx_bool bVerbose;
    gmx_bool bCompact;
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
    int resetstep;
    int nmultisim;
    int repl_ex_nst;
    int repl_ex_seed;
    real pforce;
    real cpt_period;
    real max_hours;
    const char *deviceOptions;
    unsigned long Flags;
    int ret; /* return value */
};


/* The function used for spawning threads. Extracts the mdrunner() 
   arguments from its one argument and calls mdrunner(), after making
   a commrec. */
static void mdrunner_start_fn(void *arg)
{
    struct mdrunner_arglist *mda=(struct mdrunner_arglist*)arg;
    struct mdrunner_arglist mc=*mda; /* copy the arg list to make sure 
                                        that it's thread-local. This doesn't
                                        copy pointed-to items, of course,
                                        but those are all const. */
    t_commrec *cr;  /* we need a local version of this */
    FILE *fplog=NULL;
    t_filenm *fnm;

    fnm = dup_tfn(mc.nfile, mc.fnm);

    cr = init_par_threads(mc.cr);

    if (MASTER(cr))
    {
        fplog=mc.fplog;
    }

    mda->ret=mdrunner(cr->nnodes, fplog, cr, mc.nfile, fnm, mc.oenv, 
                      mc.bVerbose, mc.bCompact, mc.nstglobalcomm, 
                      mc.ddxyz, mc.dd_node_order, mc.rdd,
                      mc.rconstr, mc.dddlb_opt, mc.dlb_scale, 
                      mc.ddcsx, mc.ddcsy, mc.ddcsz, mc.nstepout, mc.resetstep, 
                      mc.nmultisim, mc.repl_ex_nst, mc.repl_ex_seed, mc.pforce, 
                      mc.cpt_period, mc.max_hours, mc.deviceOptions, mc.Flags);
}

/* called by mdrunner() to start a specific number of threads (including 
   the main thread) for thread-parallel runs. This in turn calls mdrunner()
   for each thread. 
   All options besides nthreads are the same as for mdrunner(). */
static t_commrec *mdrunner_start_threads(int nthreads, 
              FILE *fplog,t_commrec *cr,int nfile, 
              const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
              gmx_bool bCompact, int nstglobalcomm,
              ivec ddxyz,int dd_node_order,real rdd,real rconstr,
              const char *dddlb_opt,real dlb_scale,
              const char *ddcsx,const char *ddcsy,const char *ddcsz,
              int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
              int repl_ex_seed, real pforce,real cpt_period, real max_hours, 
              const char *deviceOptions, unsigned long Flags)
{
    int ret;
    struct mdrunner_arglist *mda;
    t_commrec *crn; /* the new commrec */
    t_filenm *fnmn;

    /* first check whether we even need to start tMPI */
    if (nthreads<2)
        return cr;

    /* a few small, one-time, almost unavoidable memory leaks: */
    snew(mda,1);
    fnmn=dup_tfn(nfile, fnm);

    /* fill the data structure to pass as void pointer to thread start fn */
    mda->fplog=fplog;
    mda->cr=cr;
    mda->nfile=nfile;
    mda->fnm=fnmn;
    mda->oenv=oenv;
    mda->bVerbose=bVerbose;
    mda->bCompact=bCompact;
    mda->nstglobalcomm=nstglobalcomm;
    mda->ddxyz[XX]=ddxyz[XX];
    mda->ddxyz[YY]=ddxyz[YY];
    mda->ddxyz[ZZ]=ddxyz[ZZ];
    mda->dd_node_order=dd_node_order;
    mda->rdd=rdd;
    mda->rconstr=rconstr;
    mda->dddlb_opt=dddlb_opt;
    mda->dlb_scale=dlb_scale;
    mda->ddcsx=ddcsx;
    mda->ddcsy=ddcsy;
    mda->ddcsz=ddcsz;
    mda->nstepout=nstepout;
    mda->resetstep=resetstep;
    mda->nmultisim=nmultisim;
    mda->repl_ex_nst=repl_ex_nst;
    mda->repl_ex_seed=repl_ex_seed;
    mda->pforce=pforce;
    mda->cpt_period=cpt_period;
    mda->max_hours=max_hours;
    mda->deviceOptions=deviceOptions;
    mda->Flags=Flags;

    fprintf(stderr, "Starting %d threads\n",nthreads);
    fflush(stderr);
    /* now spawn new threads that start mdrunner_start_fn(), while 
       the main thread returns */
    ret=tMPI_Init_fn(TRUE, nthreads, mdrunner_start_fn, (void*)(mda) );
    if (ret!=TMPI_SUCCESS)
        return NULL;

    /* make a new comm_rec to reflect the new situation */
    crn=init_par_threads(cr);
    return crn;
}


/* get the number of threads based on how many there were requested, 
   which algorithms we're using, and how many particles there are. */
static int get_nthreads(int nthreads_requested, t_inputrec *inputrec,
                        gmx_mtop_t *mtop)
{
    int nthreads,nthreads_new;
    int min_atoms_per_thread;
    char *env;

    nthreads = nthreads_requested;

    /* determine # of hardware threads. */
    if (nthreads_requested < 1)
    {
        if ((env = getenv("GMX_MAX_THREADS")) != NULL)
        {
            nthreads = 0;
            sscanf(env,"%d",&nthreads);
            if (nthreads < 1)
            {
                gmx_fatal(FARGS,"GMX_MAX_THREADS (%d) should be larger than 0",
                          nthreads);
            }
        }
        else
        {
            nthreads = tMPI_Thread_get_hw_number();
        }
    }

    if (inputrec->eI == eiNM || EI_TPI(inputrec->eI))
    {
        /* Steps are divided over the nodes iso splitting the atoms */
        min_atoms_per_thread = 0;
    }
    else
    {
        min_atoms_per_thread = MIN_ATOMS_PER_THREAD;
    }

    /* Check if an algorithm does not support parallel simulation.  */
    if (nthreads != 1 && 
        ( inputrec->eI == eiLBFGS ||
          inputrec->coulombtype == eelEWALD ) )
    {
        fprintf(stderr,"\nThe integration or electrostatics algorithm doesn't support parallel runs. Not starting any threads.\n");
        nthreads = 1;
    }
    else if (nthreads_requested < 1 &&
             mtop->natoms/nthreads < min_atoms_per_thread)
    {
        /* the thread number was chosen automatically, but there are too many
           threads (too few atoms per thread) */
        nthreads_new = max(1,mtop->natoms/min_atoms_per_thread);

        if (nthreads_new > 8 || (nthreads == 8 && nthreads_new > 4))
        {
            /* Use only multiples of 4 above 8 threads
             * or with an 8-core processor
             * (to avoid 6 threads on 8 core processors with 4 real cores).
             */
            nthreads_new = (nthreads_new/4)*4;
        }
        else if (nthreads_new > 4)
        {
            /* Avoid 5 or 7 threads */
            nthreads_new = (nthreads_new/2)*2;
        }

        nthreads = nthreads_new;

        fprintf(stderr,"\n");
        fprintf(stderr,"NOTE: Parallelization is limited by the small number of atoms,\n");
        fprintf(stderr,"      only starting %d threads.\n",nthreads);
        fprintf(stderr,"      You can use the -nt option to optimize the number of threads.\n\n");
    }
    return nthreads;
}
#endif


int mdrunner(int nthreads_requested, FILE *fplog,t_commrec *cr,int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm,
             ivec ddxyz,int dd_node_order,real rdd,real rconstr,
             const char *dddlb_opt,real dlb_scale,
             const char *ddcsx,const char *ddcsy,const char *ddcsz,
             int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
             int repl_ex_seed, real pforce,real cpt_period,real max_hours,
             const char *deviceOptions, unsigned long Flags)
{
    double     nodetime=0,realtime;
    t_inputrec *inputrec;
    t_state    *state=NULL;
    matrix     box;
    gmx_ddbox_t ddbox={0};
    int        npme_major,npme_minor;
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
    gmx_bool       bReadRNG,bReadEkin;
    int        list;
    gmx_runtime_t runtime;
    int        rc;
    gmx_large_int_t reset_counters;
    gmx_edsam_t ed=NULL;
    t_commrec   *cr_old=cr; 
    int         nthreads=1;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
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

    snew(state,1);
    if (MASTER(cr)) 
    {
        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPX,nfile,fnm),inputrec,state,NULL,mtop);

        /* NOW the threads will be started: */
#ifdef GMX_THREADS
        nthreads = get_nthreads(nthreads_requested, inputrec, mtop);

        if (nthreads > 1)
        {
            /* now start the threads. */
            cr=mdrunner_start_threads(nthreads, fplog, cr_old, nfile, fnm, 
                                      oenv, bVerbose, bCompact, nstglobalcomm, 
                                      ddxyz, dd_node_order, rdd, rconstr, 
                                      dddlb_opt, dlb_scale, ddcsx, ddcsy, ddcsz,
                                      nstepout, resetstep, nmultisim, 
                                      repl_ex_nst, repl_ex_seed, pforce, 
                                      cpt_period, max_hours, deviceOptions, 
                                      Flags);
            /* the main thread continues here with a new cr. We don't deallocate
               the old cr because other threads may still be reading it. */
            if (cr == NULL)
            {
                gmx_comm("Failed to spawn threads");
            }
        }
#endif
    }
    /* END OF CAUTION: cr is now reliable */

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(fplog, cr, inputrec, mtop);
    }
    if (fplog != NULL)
    {
        pr_inputrec(fplog,0,"Input Parameters",inputrec,FALSE);
    }

    /* now make sure the state is initialized and propagated */
    set_state_entries(state,inputrec,cr->nnodes);

    /* remove when vv and rerun works correctly! */
    if (PAR(cr) && EI_VV(inputrec->eI) && ((Flags & MD_RERUN) || (Flags & MD_RERUN_VSITE)))
    {
        gmx_fatal(FARGS, "Currently can't do velocity verlet with rerun in parallel.");
    }
    if (EI_VV(inputrec->eI) && etcVRESCALE == inputrec->etc)
    {
        gmx_fatal(FARGS, "In GROMACS 4.5.x, velocity-Verlet integrators do not work with velocity-rescaling temperature coupling. They do work in 4.6.x. Please ugprade your GROMACS version.");
    }

    /* A parallel command line option consistency check that we can
       only do after any threads have started. */
    if (!PAR(cr) &&
        (ddxyz[XX] > 1 || ddxyz[YY] > 1 || ddxyz[ZZ] > 1 || cr->npmenodes > 0))
    {
        gmx_fatal(FARGS,
                  "The -dd or -npme option request a parallel simulation, "
#ifndef GMX_MPI
                  "but mdrun was compiled without threads or MPI enabled"
#else
#ifdef GMX_THREADS
                  "but the number of threads (option -nt) is 1"
#else
                  "but mdrun was not started through mpirun/mpiexec or only one process was requested through mpirun/mpiexec" 
#endif
#endif
            );
    }

    if ((Flags & MD_RERUN) &&
        (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
    }

    if (can_use_allvsall(inputrec,mtop,TRUE,cr,fplog))
    {
        /* All-vs-all loops do not work with domain decomposition */
        Flags |= MD_PARTDEC;
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
    init_disres(fplog,mtop,inputrec,cr,Flags & MD_PARTDEC,fcd,state, repl_ex_nst > 0);

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
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        deform_init_init_step_tpx = inputrec->init_step;
        copy_mat(box,deform_init_box_tpx);
#ifdef GMX_THREADS
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
    }

    if (opt2bSet("-cpi",nfile,fnm)) 
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        if( gmx_fexist_master(opt2fn_master("-cpi",nfile,fnm,cr),cr) )
        {
            load_checkpoint(opt2fn_master("-cpi",nfile,fnm,cr),&fplog,
                            cr,Flags & MD_PARTDEC,ddxyz,
                            inputrec,state,&bReadRNG,&bReadEkin,
                            (Flags & MD_APPENDFILES),
                            (Flags & MD_APPENDFILESSET));
            
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

    if (((MASTER(cr) || (Flags & MD_SEPPOT)) && (Flags & MD_APPENDFILES))
#ifdef GMX_THREADS
        /* With thread MPI only the master node/thread exists in mdrun.c,
         * therefore non-master nodes need to open the "seppot" log file here.
         */
        || (!MASTER(cr) && (Flags & MD_SEPPOT))
#endif
        )
    {
        gmx_log_open(ftp2fn(efLOG,nfile,fnm),cr,!(Flags & MD_SEPPOT),
                             Flags,&fplog);
    }

    if (SIMMASTER(cr)) 
    {
        copy_mat(state->box,box);
    }

    if (PAR(cr)) 
    {
        gmx_bcast(sizeof(box),box,cr);
    }

    /* Essential dynamics */
    if (opt2bSet("-ei",nfile,fnm))
    {
        /* Open input and output files, allocate space for ED data structure */
        ed = ed_open(nfile,fnm,Flags,cr);
    }

    if (bVerbose && SIMMASTER(cr))
    {
        fprintf(stderr,"Loaded with Money\n\n");
    }

    if (PAR(cr) && !((Flags & MD_PARTDEC) ||
                     EI_TPI(inputrec->eI) ||
                     inputrec->eI == eiNM))
    {
        cr->dd = init_domain_decomposition(fplog,cr,Flags,ddxyz,rdd,rconstr,
                                           dddlb_opt,dlb_scale,
                                           ddcsx,ddcsy,ddcsz,
                                           mtop,inputrec,
                                           box,state->x,
                                           &ddbox,&npme_major,&npme_minor);

        make_dd_communicators(fplog,cr,dd_node_order);

        /* Set overallocation to avoid frequent reallocation of arrays */
        set_over_alloc_dd(TRUE);
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty = (DUTY_PP | DUTY_PME);
        npme_major = 1;
        npme_minor = 1;
        if (!EI_TPI(inputrec->eI))
        {
            npme_major = cr->nnodes;
        }
        
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

    wcycle = wallcycle_init(fplog,resetstep,cr);
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
            status = gmx_pme_init(pmedata,cr,npme_major,npme_minor,inputrec,
                                  mtop ? mtop->natoms : 0,nChargePerturbed,
                                  (Flags & MD_REPRODUCIBLE));
            if (status != 0) 
            {
                gmx_fatal(FARGS,"Error %d initializing PME",status);
            }
        }
    }


    if (integrator[inputrec->eI].func == do_md
#ifdef GMX_OPENMM
        ||
        integrator[inputrec->eI].func == do_md_openmm
#endif
        )
    {
        /* Turn on signal handling on all nodes */
        /*
         * (A user signal from the PME nodes (if any)
         * is communicated to the PP nodes.
         */
        signal_handler_install();
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
                                      deviceOptions,
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
        if (SIMMASTER(cr))
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

    rc=(int)gmx_get_stop_condition();

#ifdef GMX_THREADS
    /* we need to join all threads. The sub-threads join when they
       exit this function, but the master thread needs to be told to 
       wait for that. */
    if (PAR(cr) && MASTER(cr))
    {
        tMPI_Finalize();
    }
#endif

    return rc;
}
