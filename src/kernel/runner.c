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
#ifdef __linux
#define _GNU_SOURCE
#include <sched.h>
#include <sys/syscall.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <string.h>

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
#include "gmx_omp_nthreads.h"
#include "pull_rotation.h"
#include "calc_verletbuf.h"
#include "gmx_fatal_collective.h"

#include "md_openmm.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

#ifdef GMX_OPENMM
#include "md_openmm.h"
#endif

#ifdef GMX_OPENMP
#include <omp.h>
#endif

#ifdef GMX_GPU
#include "gpu_utils.h"
#include "nbnxn_cuda_data_mgmt.h"
#endif /* GMX_GPU */


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
#ifdef GMX_THREAD_MPI
tMPI_Thread_mutex_t deform_init_box_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


#ifdef GMX_THREAD_MPI
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
    const char *nbpu_opt;
    int nsteps_cmdline;
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
                      mc.ddcsx, mc.ddcsy, mc.ddcsz,
                      mc.nbpu_opt,
                      mc.nsteps_cmdline, mc.nstepout, mc.resetstep,
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
              const char *nbpu_opt,
              int nsteps_cmdline, int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
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
    mda->nbpu_opt=nbpu_opt;
    mda->nsteps_cmdline=nsteps_cmdline;
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

    fprintf(stderr, "Starting %d tMPI threads\n",nthreads);
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


/* Get the number of threads to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 */
static int get_nthreads_mpi(int nthreads_requested, t_inputrec *inputrec,
                            gmx_mtop_t *mtop)
{
    int nthreads,nthreads_new;
    int min_atoms_per_thread;
    char *env;

    nthreads = nthreads_requested;

    /* determine # of hardware threads. */
    if (nthreads_requested < 1)
    {
        if ((env = getenv("GMX_MAX_MPI_THREADS")) != NULL)
        {
            nthreads = 0;
            sscanf(env,"%d",&nthreads);
            if (nthreads < 1)
            {
                gmx_fatal(FARGS,"GMX_MAX_MPI_THREADS (%d) should be larger than 0",
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
#endif /* GMX_THREAD_MPI */


/* Environment variable for setting nstlist */
#define NSTLIST_ENVVAR      "GMX_NSTLIST"
/* Try to increase nstlist when using a GPU with nstlist less than this */
#define NSTLIST_GPU_ENOUGH  20
/* Increase nstlist until the non-bonded cost increases more than this factor */
#define NBNXN_GPU_LIST_OK_FAC   1.25
/* Don't increase nstlist beyond a non-bonded cost increases of this factor */
#define NBNXN_GPU_LIST_MAX_FAC  1.40

/* Try to increase nstlist when running on a GPU */
static void increase_nstlist(FILE *fp,t_commrec *cr,
                             t_inputrec *ir,gmx_mtop_t *mtop,matrix box)
{
#define NNSTL 4
    int  nstl[NNSTL]={ 20, 25, 40, 50 };
    char *env;
    int  nstlist_orig,nstlist_prev;
    real rlist_inc,rlist_ok,rlist_max,rlist_new,rlist_prev;
    int  i;
    t_state state_tmp;
    gmx_bool bBox,bDD,bCont;
    const char *nstl_fmt="\nFor optimal performace with a GPU nstlist (now %d) should be larger\nThe optimum depends on your CPU and GPU resources\nYou might want to try nstlist values using the env.var. GMX_NSTLIST\n";
    const char *vbd_err="Can not increase nstlist for GPU run because verlet-buffer-drift is not set or used";
    const char *box_err="Can not increase nstlist for GPU run because the box is too small";
    const char *dd_err ="Can not increase nstlist for GPU run because of domain decomposition limitations";
    char buf[STRLEN];

    env = getenv(NSTLIST_ENVVAR);
    if (env == NULL)
    {
        if (MASTER(cr))
        {
            fprintf(stderr,nstl_fmt,ir->nstlist);
        }
        if (fp != NULL)
        {
            fprintf(fp,nstl_fmt,ir->nstlist);
        }
    }

    if (ir->verletbuf_drift == 0)
    {
        gmx_fatal(FARGS,"You are using an old tpr file with a GPU, please generate a new tpr file with an up to date version of grompp");
    }

    if (ir->verletbuf_drift < 0)
    {
        if (MASTER(cr))
        {
            fprintf(stderr,"%s\n",vbd_err);
        }
        if (fp != NULL)
        {
            fprintf(fp,"%s\n",vbd_err);
        }

        return;
    }

    nstlist_orig = ir->nstlist;
    if (env != NULL)
    {
        sprintf(buf,"Getting nstlist from env.var. GMX_NSTLIST=%s",env);
        if (MASTER(cr))
        {
            fprintf(stderr,"%s\n",buf);
        }
        if (fp != NULL)
        {
            fprintf(fp,"%s\n",buf);
        }
        sscanf(env,"%d",&ir->nstlist);
    }

    /* Allow rlist to make the list double the size of the cut-off sphere */
    rlist_inc = nbnxn_rlist_inc(NBNXN_GPU_CLUSTER_SIZE,mtop->natoms/det(box));
    rlist_ok  = (max(ir->rvdw,ir->rcoulomb) + rlist_inc)*pow(NBNXN_GPU_LIST_OK_FAC,1.0/3.0) - rlist_inc;
    rlist_max = (max(ir->rvdw,ir->rcoulomb) + rlist_inc)*pow(NBNXN_GPU_LIST_MAX_FAC,1.0/3.0) - rlist_inc;
    if (debug)
    {
        fprintf(debug,"GPU nstlist tuning: rlist_inc %.3f rlist_max %.3f\n",
                rlist_inc,rlist_max);
    }

    i = 0;
    nstlist_prev = nstlist_orig;
    rlist_prev   = ir->rlist;
    do
    {
        if (env == NULL)
        {
            ir->nstlist = nstl[i];
        }

        /* Set the pair-list buffer size in ir */
        calc_verlet_buffer_size(mtop,det(box),ir,ir->verletbuf_drift,
                                NULL,&rlist_new);

        /* Does rlist fit in the box? */
        bBox = (sqr(rlist_new) < max_cutoff2(ir->ePBC,box));
        bDD  = TRUE;
        if (bBox && DOMAINDECOMP(cr))
        {
            /* Check if rlist fits in the domain decomposition */
            if (inputrec2nboundeddim(ir) < DIM)
            {
                gmx_incons("Changing nstlist with domain decomposition and unbounded dimensions is not implemented yet");
            }
            copy_mat(box,state_tmp.box);
            bDD = change_dd_cutoff(cr,&state_tmp,ir,rlist_new);
        }

        bCont = FALSE;

        if (env == NULL)
        {
            if (bBox && bDD && rlist_new <= rlist_max)
            {
                /* Increase nstlist */
                nstlist_prev = ir->nstlist;
                rlist_prev   = rlist_new;
                bCont = (rlist_new < rlist_ok);
            }
            else
            {
                /* Stick with the previous nstlist */
                ir->nstlist = nstlist_prev;
                rlist_new   = rlist_prev;
                bBox = TRUE;
                bDD  = TRUE;
            }
        }

        i++;
    }
    while (bCont);

    if (!bBox || !bDD)
    {
        gmx_warning(!bBox ? box_err : dd_err);
        if (fp != NULL)
        {
            fprintf(fp,"\n%s\n",bBox ? box_err : dd_err);
        }
        ir->nstlist = nstlist_orig;
    }
    else
    {
        sprintf(buf,"Changing nstlist from %d to %d, rlist from %g to %g",
                nstlist_orig,ir->nstlist,
                ir->rlist,rlist_new);
        if (MASTER(cr))
        {
            fprintf(stderr,"%s\n\n",buf);
        }
        if (fp != NULL)
        {
            fprintf(fp,"%s\n\n",buf);
        }
        ir->rlist     = rlist_new;
        ir->rlistlong = rlist_new;
    }
}

static void convert_to_verlet_scheme(FILE *fplog,
                                     t_inputrec *ir,
                                     gmx_mtop_t *mtop,real box_vol)
{
    char *conv_mesg="Converting input file with group cut-off scheme to the Verlet cut-off scheme";

    if (fplog != NULL)
    {
        fprintf(fplog,"\n%s\n\n",conv_mesg);
    }
    fprintf(stderr,"\n%s\n\n",conv_mesg);

    if (!(ir->vdwtype == evdwCUT &&
          (ir->coulombtype == eelCUT ||
           EEL_RF(ir->coulombtype) ||
           ir->coulombtype == eelPME) &&
          ir->rcoulomb == ir->rvdw))
    {
        gmx_fatal(FARGS,"Can only convert old tpr files to the Verlet cut-off scheme with cut-off LJ interactions and PME, RF or cut-off electrostatics and rcoulomb=rvdw");
    }

    if (inputrec2nboundeddim(ir) != 3)
    {
        gmx_fatal(FARGS,"Can only convert old tpr files to the Verlet cut-off scheme with 3D pbc");
    }

    if (EI_DYNAMICS(ir->eI) && ir->etc == etcNO)
    {
        gmx_fatal(FARGS,"Will not convert old tpr files to the Verlet cut-off scheme without temperature coupling");
    }

    if (ir->efep != efepNO || ir->implicit_solvent != eisNO)
    {
        gmx_fatal(FARGS,"Will not convert old tpr files to the Verlet cut-off scheme with free-energy calculations or implicit solvent");
    }

    ir->cutoff_scheme   = ecutsVERLET;
    ir->verletbuf_drift = 0.005;

    if (EI_DYNAMICS(ir->eI))
    {
        calc_verlet_buffer_size(mtop,box_vol,ir,ir->verletbuf_drift,
                                NULL,&ir->rlist);
    }
    else
    {
        ir->rlist = 1.05*max(ir->rvdw,ir->rcoulomb);
    }

    gmx_mtop_remove_chargegroups(mtop);
}

int mdrunner(int nthreads_requested, FILE *fplog,t_commrec *cr,int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm,
             ivec ddxyz,int dd_node_order,real rdd,real rconstr,
             const char *dddlb_opt,real dlb_scale,
             const char *ddcsx,const char *ddcsy,const char *ddcsz,
             const char *nbpu_opt,
             int nsteps_cmdline, int nstepout,int resetstep,int nmultisim,int repl_ex_nst,
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
    int        nbnxn_kernel=-1;
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
    int         nthreads_mpi=1;
    int         nthreads_pme=1;
    int         nthreads_pp=1;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    snew(inputrec,1);
    snew(mtop,1);
    
    if (Flags & MD_APPENDFILES) 
    {
        fplog = NULL;
    }

    gmx_omp_nthreads_detecthw();

    snew(state,1);
    if (MASTER(cr)) 
    {
        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPX,nfile,fnm),inputrec,state,NULL,mtop);

        if (inputrec->cutoff_scheme != ecutsVERLET &&
            getenv("GMX_VERLET_SCHEME") != NULL)
        {
            convert_to_verlet_scheme(fplog,inputrec,mtop,det(state->box));
        }

        /* NOW the threads will be started: */
#ifdef GMX_THREAD_MPI
        nthreads_mpi = get_nthreads_mpi(nthreads_requested, inputrec, mtop);

        if (nthreads_mpi > 1)
        {
            /* now start the threads. */
            cr=mdrunner_start_threads(nthreads_mpi, fplog, cr_old, nfile, fnm, 
                                      oenv, bVerbose, bCompact, nstglobalcomm, 
                                      ddxyz, dd_node_order, rdd, rconstr, 
                                      dddlb_opt, dlb_scale, ddcsx, ddcsy, ddcsz,
                                      nbpu_opt,
                                      nsteps_cmdline, nstepout, resetstep, nmultisim, 
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
#ifdef GMX_THREAD_MPI
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
        if (cr->npmenodes > 0)
        {
            if (!EEL_PME(inputrec->coulombtype))
            {
                gmx_fatal_collective(FARGS,cr,NULL,
                                     "PME nodes are requested, but the system does not use PME electrostatics");
            }
            if (Flags & MD_PARTDEC)
            {
                gmx_fatal_collective(FARGS,cr,NULL,
                                     "PME nodes are requested, but particle decomposition does not support separate PME nodes");
            }
        }

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
#ifdef GMX_THREAD_MPI
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        deform_init_init_step_tpx = inputrec->init_step;
        copy_mat(box,deform_init_box_tpx);
#ifdef GMX_THREAD_MPI
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
#ifdef GMX_THREAD_MPI
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

    /* Initialize per-node process ID and counters. */
    gmx_init_intra_counters(cr);

    gmx_omp_nthreads_init(fplog, cr, (cr->duty & DUTY_PP) == 0,
                          inputrec->cutoff_scheme == ecutsVERLET);

    /* getting number of PP/PME threads
       PME: env variable should be read only on one node to make sure it is 
       identical everywhere;
     */
    /* TODO nthreads_pp is only used for pinning threads.
     * This is a temporary solution.
     */
    nthreads_pp  = gmx_omp_nthreads_get(emntNonbonded);
    nthreads_pme = gmx_omp_nthreads_get(emntPME);

    wcycle = wallcycle_init(fplog,resetstep,cr,nthreads_pp,nthreads_pme);

    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes 
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast_sim(sizeof(reset_counters),&reset_counters,cr);
        wcycle_set_reset_counters(wcycle, reset_counters);
    }

    /* override nsteps if defined on the command line */
    if (nsteps_cmdline >= -1)
    {
        char stmp[STRLEN];

        inputrec->nsteps = nsteps_cmdline;
        if (EI_DYNAMICS(inputrec->eI))
        {
            sprintf(stmp, "Overriding nsteps with value passed on the command line: %d steps, %.3f ps",
                    nsteps_cmdline, nsteps_cmdline*inputrec->delta_t);
        }
        else
        {
            sprintf(stmp, "Overriding nsteps with value passed on the command line: %d steps",
                    nsteps_cmdline);
        }

        if (SIMMASTER(cr))
        {
            fprintf(stderr, "\n%s\n\n", stmp);
        }
        if (fplog)
        {
            fprintf(fplog, "%s\n\n", stmp);
        }
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

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            gmx_bool useGPU;

            /* With GPU we should check nstlist for performance */
            pick_nbnxn_kernel(fplog, cr,
                              (strncmp(nbpu_opt,"gpu",3) == 0 ||
                               strcmp(nbpu_opt,"auto") == 0),
                              &useGPU,
                              strncmp(nbpu_opt,"gpu",3) == 0,
                              EEL_FULL(inputrec->coulombtype),
                              &nbnxn_kernel);
            if ((EI_DYNAMICS(inputrec->eI) &&
                 (nbnxn_kernel == nbk8x8x8CUDA ||
                  nbnxn_kernel == nbk8x8x8PlainC) &&
                 inputrec->nstlist < NSTLIST_GPU_ENOUGH) ||
                getenv(NSTLIST_ENVVAR) != NULL)
            {
                /* Choose a better nstlist */
                increase_nstlist(fplog,cr,inputrec,mtop,box);
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
                      opt2fn("-tabletf",nfile,fnm),
                      opt2fn("-tablep",nfile,fnm),
                      opt2fn("-tableb",nfile,fnm),
                      nbpu_opt,nbnxn_kernel,
                      FALSE,pforce);

        /* version for PCA_NOT_READ_NODE (see md.c) */
        /*init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
          "nofile","nofile","nofile","nofile",FALSE,pforce);
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

        /* Set CPU affinity. Can be important for performance.
           On some systems (e.g. Cray) CPU Affinity is set by default.
           But default assigning doesn't work (well) with only some ranks
           having threads. This causes very low performance.
           External tools have cumbersome syntax for setting affinity
           in the case that only some ranks have threads.
           Thus it is important that GROMACS sets the affinity internally at
           if only PME is using threads.
        */

#ifdef GMX_OPENMP /* TODO: actually we could do this even without OpenMP! */
#ifdef __linux
    if (getenv("GMX_NO_THREAD_PINNING") == NULL)
    {
        int core, local_nthreads, offset;

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            local_nthreads = gmx_omp_nthreads_get(emntNonbonded);
        }
        else
        {
            /* threads on this node */
            local_nthreads = (cr->duty & DUTY_PME) ? nthreads_pme : 1;
        }

        /* map the current process to cores */
        if (PAR(cr) || MULTISIM(cr))
        {
            core = cr->nodeid_intra*local_nthreads;
        }
        else
        {
            core = 0;
        }

        char *env;
        if ((env = getenv("GMX_THREAD_PINNING_OFFSET")) != NULL)
        {
            char *end;
            offset = strtol(env, &end, 10);
            if (!end || (*end != 0))
            {
                gmx_fatal(FARGS, "Invalid thread pinning offset: %s", env);
            }

            fprintf(stderr, "Applying thread pinning offset %d\n", offset);
            core += offset;
        }

        /* set the per-thread affinity */
#pragma omp parallel firstprivate(core) num_threads(local_nthreads)
        {
            cpu_set_t mask;
            CPU_ZERO(&mask);
            core += omp_get_thread_num();
            CPU_SET(core, &mask);
            sched_setaffinity((pid_t) syscall (SYS_gettid), sizeof(cpu_set_t), &mask);
        }
    }
    else
    {
        char sbuf[STRLEN];
        sprintf(sbuf, "NOTE: Thread pinning turned off by the "
                "GMX_NO_THREAD_PINNING environment variable");
        if (SIMMASTER(cr))
        {
            fprintf(stderr, "%s\n\n", sbuf);
        }
        if (fplog)
        {
            fprintf(fplog, "\n%s\n", sbuf);
        }
    }
#endif /* __linux    */
#endif /* GMX_OPENMP */


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
                                  (Flags & MD_REPRODUCIBLE),nthreads_pme);
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
        
        if (inputrec->bRot)
        {
           /* Initialize enforced rotation code */
           init_rot(fplog,inputrec,nfile,fnm,cr,state->x,state->box,mtop,oenv,
                    bVerbose,Flags);
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
        
        if (inputrec->bRot)
        {
            finish_rot(fplog,inputrec->rot);
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
#ifdef GMX_GPU
               fr != NULL && fr->nbv != NULL && fr->nbv->useGPU ?
                 nbnxn_cuda_get_timings(fr->nbv->cu_nbv) :
#endif
               NULL,
               nthreads_pp, 
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

#ifdef GMX_GPU
    if (cr->duty & DUTY_PP && fr->nbv != NULL && fr->nbv->useGPU)
    {
        int gpu_device_id = cr->nodeid; /* FIXME get dev_id */

        /* free GPU memory and uninitialize GPU */
        nbnxn_cuda_free(fplog, fr->nbv->cu_nbv, DOMAINDECOMP(cr));

        if (uninit_gpu(fplog, gpu_device_id) != 0)
        {
            gmx_warning("Failed to uninitialize GPU.");
        }
    }
#endif

    /* Does what it says */  
    print_date_and_time(fplog,cr->nodeid,"Finished mdrun",&runtime);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_close(fplog);
    }	

    rc=(int)gmx_get_stop_condition();

#ifdef GMX_THREAD_MPI
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
