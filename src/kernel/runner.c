/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "statutil.h"
#include "mdrun.h"
#include "md_logging.h"
#include "md_support.h"
#include "network.h"
#include "pull.h"
#include "names.h"
#include "disre.h"
#include "orires.h"
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
#include "gmx_detect_hardware.h"
#include "gmx_omp_nthreads.h"
#include "pull_rotation.h"
#include "calc_verletbuf.h"
#include "../mdlib/nbnxn_search.h"
#include "../mdlib/nbnxn_consts.h"
#include "gmx_fatal_collective.h"
#include "membed.h"
#include "gmx_omp.h"
#include "gmx_thread_affinity.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

#include "gpu_utils.h"
#include "nbnxn_cuda_data_mgmt.h"

typedef struct {
    gmx_integrator_t *func;
} gmx_intp_t;

/* The array should match the eI array in include/types/enums.h */
const gmx_intp_t    integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md}, {do_md}, {do_md}};

gmx_large_int_t     deform_init_init_step_tpx;
matrix              deform_init_box_tpx;
#ifdef GMX_THREAD_MPI
tMPI_Thread_mutex_t deform_init_box_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif


#ifdef GMX_THREAD_MPI
struct mdrunner_arglist
{
    gmx_hw_opt_t   *hw_opt;
    FILE           *fplog;
    t_commrec      *cr;
    int             nfile;
    const t_filenm *fnm;
    output_env_t    oenv;
    gmx_bool        bVerbose;
    gmx_bool        bCompact;
    int             nstglobalcomm;
    ivec            ddxyz;
    int             dd_node_order;
    real            rdd;
    real            rconstr;
    const char     *dddlb_opt;
    real            dlb_scale;
    const char     *ddcsx;
    const char     *ddcsy;
    const char     *ddcsz;
    const char     *nbpu_opt;
    gmx_large_int_t nsteps_cmdline;
    int             nstepout;
    int             resetstep;
    int             nmultisim;
    int             repl_ex_nst;
    int             repl_ex_nex;
    int             repl_ex_seed;
    real            pforce;
    real            cpt_period;
    real            max_hours;
    const char     *deviceOptions;
    unsigned long   Flags;
    int             ret; /* return value */
};


/* The function used for spawning threads. Extracts the mdrunner()
   arguments from its one argument and calls mdrunner(), after making
   a commrec. */
static void mdrunner_start_fn(void *arg)
{
    struct mdrunner_arglist *mda = (struct mdrunner_arglist*)arg;
    struct mdrunner_arglist  mc  = *mda; /* copy the arg list to make sure
                                            that it's thread-local. This doesn't
                                            copy pointed-to items, of course,
                                            but those are all const. */
    t_commrec *cr;                       /* we need a local version of this */
    FILE      *fplog = NULL;
    t_filenm  *fnm;

    fnm = dup_tfn(mc.nfile, mc.fnm);

    cr = init_par_threads(mc.cr);

    if (MASTER(cr))
    {
        fplog = mc.fplog;
    }

    mda->ret = mdrunner(mc.hw_opt, fplog, cr, mc.nfile, fnm, mc.oenv,
                        mc.bVerbose, mc.bCompact, mc.nstglobalcomm,
                        mc.ddxyz, mc.dd_node_order, mc.rdd,
                        mc.rconstr, mc.dddlb_opt, mc.dlb_scale,
                        mc.ddcsx, mc.ddcsy, mc.ddcsz,
                        mc.nbpu_opt,
                        mc.nsteps_cmdline, mc.nstepout, mc.resetstep,
                        mc.nmultisim, mc.repl_ex_nst, mc.repl_ex_nex, mc.repl_ex_seed, mc.pforce,
                        mc.cpt_period, mc.max_hours, mc.deviceOptions, mc.Flags);
}

/* called by mdrunner() to start a specific number of threads (including
   the main thread) for thread-parallel runs. This in turn calls mdrunner()
   for each thread.
   All options besides nthreads are the same as for mdrunner(). */
static t_commrec *mdrunner_start_threads(gmx_hw_opt_t *hw_opt,
                                         FILE *fplog, t_commrec *cr, int nfile,
                                         const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
                                         gmx_bool bCompact, int nstglobalcomm,
                                         ivec ddxyz, int dd_node_order, real rdd, real rconstr,
                                         const char *dddlb_opt, real dlb_scale,
                                         const char *ddcsx, const char *ddcsy, const char *ddcsz,
                                         const char *nbpu_opt,
                                         gmx_large_int_t nsteps_cmdline,
                                         int nstepout, int resetstep,
                                         int nmultisim, int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
                                         real pforce, real cpt_period, real max_hours,
                                         const char *deviceOptions, unsigned long Flags)
{
    int                      ret;
    struct mdrunner_arglist *mda;
    t_commrec               *crn; /* the new commrec */
    t_filenm                *fnmn;

    /* first check whether we even need to start tMPI */
    if (hw_opt->nthreads_tmpi < 2)
    {
        return cr;
    }

    /* a few small, one-time, almost unavoidable memory leaks: */
    snew(mda, 1);
    fnmn = dup_tfn(nfile, fnm);

    /* fill the data structure to pass as void pointer to thread start fn */
    mda->hw_opt         = hw_opt;
    mda->fplog          = fplog;
    mda->cr             = cr;
    mda->nfile          = nfile;
    mda->fnm            = fnmn;
    mda->oenv           = oenv;
    mda->bVerbose       = bVerbose;
    mda->bCompact       = bCompact;
    mda->nstglobalcomm  = nstglobalcomm;
    mda->ddxyz[XX]      = ddxyz[XX];
    mda->ddxyz[YY]      = ddxyz[YY];
    mda->ddxyz[ZZ]      = ddxyz[ZZ];
    mda->dd_node_order  = dd_node_order;
    mda->rdd            = rdd;
    mda->rconstr        = rconstr;
    mda->dddlb_opt      = dddlb_opt;
    mda->dlb_scale      = dlb_scale;
    mda->ddcsx          = ddcsx;
    mda->ddcsy          = ddcsy;
    mda->ddcsz          = ddcsz;
    mda->nbpu_opt       = nbpu_opt;
    mda->nsteps_cmdline = nsteps_cmdline;
    mda->nstepout       = nstepout;
    mda->resetstep      = resetstep;
    mda->nmultisim      = nmultisim;
    mda->repl_ex_nst    = repl_ex_nst;
    mda->repl_ex_nex    = repl_ex_nex;
    mda->repl_ex_seed   = repl_ex_seed;
    mda->pforce         = pforce;
    mda->cpt_period     = cpt_period;
    mda->max_hours      = max_hours;
    mda->deviceOptions  = deviceOptions;
    mda->Flags          = Flags;

    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns, we set thread affinity later */
    ret = tMPI_Init_fn(TRUE, hw_opt->nthreads_tmpi, TMPI_AFFINITY_NONE,
                       mdrunner_start_fn, (void*)(mda) );
    if (ret != TMPI_SUCCESS)
    {
        return NULL;
    }

    /* make a new comm_rec to reflect the new situation */
    crn = init_par_threads(cr);
    return crn;
}


static int get_tmpi_omp_thread_division(const gmx_hw_info_t *hwinfo,
                                        const gmx_hw_opt_t  *hw_opt,
                                        int                  nthreads_tot,
                                        int                  ngpu)
{
    int nthreads_tmpi;

    /* There are no separate PME nodes here, as we ensured in
     * check_and_update_hw_opt that nthreads_tmpi>0 with PME nodes
     * and a conditional ensures we would not have ended up here.
     * Note that separate PME nodes might be switched on later.
     */
    if (ngpu > 0)
    {
        nthreads_tmpi = ngpu;
        if (nthreads_tot > 0 && nthreads_tot < nthreads_tmpi)
        {
            nthreads_tmpi = nthreads_tot;
        }
    }
    else if (hw_opt->nthreads_omp > 0)
    {
        /* Here we could oversubscribe, when we do, we issue a warning later */
        nthreads_tmpi = max(1, nthreads_tot/hw_opt->nthreads_omp);
    }
    else
    {
        /* TODO choose nthreads_omp based on hardware topology
           when we have a hardware topology detection library */
        /* In general, when running up to 4 threads, OpenMP should be faster.
         * Note: on AMD Bulldozer we should avoid running OpenMP over two dies.
         * On Intel>=Nehalem running OpenMP on a single CPU is always faster,
         * even on two CPUs it's usually faster (but with many OpenMP threads
         * it could be faster not to use HT, currently we always use HT).
         * On Nehalem/Westmere we want to avoid running 16 threads over
         * two CPUs with HT, so we need a limit<16; thus we use 12.
         * A reasonable limit for Intel Sandy and Ivy bridge,
         * not knowing the topology, is 16 threads.
         */
        const int nthreads_omp_always_faster             =  4;
        const int nthreads_omp_always_faster_Nehalem     = 12;
        const int nthreads_omp_always_faster_SandyBridge = 16;
        const int first_model_Nehalem                    = 0x1A;
        const int first_model_SandyBridge                = 0x2A;
        gmx_bool  bIntel_Family6;

        bIntel_Family6 =
            (gmx_cpuid_vendor(hwinfo->cpuid_info) == GMX_CPUID_VENDOR_INTEL &&
             gmx_cpuid_family(hwinfo->cpuid_info) == 6);

        if (nthreads_tot <= nthreads_omp_always_faster ||
            (bIntel_Family6 &&
             ((gmx_cpuid_model(hwinfo->cpuid_info) >= nthreads_omp_always_faster_Nehalem && nthreads_tot <= nthreads_omp_always_faster_Nehalem) ||
              (gmx_cpuid_model(hwinfo->cpuid_info) >= nthreads_omp_always_faster_SandyBridge && nthreads_tot <= nthreads_omp_always_faster_SandyBridge))))
        {
            /* Use pure OpenMP parallelization */
            nthreads_tmpi = 1;
        }
        else
        {
            /* Don't use OpenMP parallelization */
            nthreads_tmpi = nthreads_tot;
        }
    }

    return nthreads_tmpi;
}


/* Get the number of threads to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than #GPU.
 */
static int get_nthreads_mpi(const gmx_hw_info_t *hwinfo,
                            gmx_hw_opt_t *hw_opt,
                            t_inputrec *inputrec, gmx_mtop_t *mtop,
                            const t_commrec *cr,
                            FILE *fplog)
{
    int      nthreads_hw, nthreads_tot_max, nthreads_tmpi, nthreads_new, ngpu;
    int      min_atoms_per_mpi_thread;
    char    *env;
    char     sbuf[STRLEN];
    gmx_bool bCanUseGPU;

    if (hw_opt->nthreads_tmpi > 0)
    {
        /* Trivial, return right away */
        return hw_opt->nthreads_tmpi;
    }

    nthreads_hw = hwinfo->nthreads_hw_avail;

    /* How many total (#tMPI*#OpenMP) threads can we start? */
    if (hw_opt->nthreads_tot > 0)
    {
        nthreads_tot_max = hw_opt->nthreads_tot;
    }
    else
    {
        nthreads_tot_max = nthreads_hw;
    }

    bCanUseGPU = (inputrec->cutoff_scheme == ecutsVERLET && hwinfo->bCanUseGPU);
    if (bCanUseGPU)
    {
        ngpu = hwinfo->gpu_info.ncuda_dev_use;
    }
    else
    {
        ngpu = 0;
    }

    nthreads_tmpi =
        get_tmpi_omp_thread_division(hwinfo, hw_opt, nthreads_tot_max, ngpu);

    if (inputrec->eI == eiNM || EI_TPI(inputrec->eI))
    {
        /* Dims/steps are divided over the nodes iso splitting the atoms */
        min_atoms_per_mpi_thread = 0;
    }
    else
    {
        if (bCanUseGPU)
        {
            min_atoms_per_mpi_thread = MIN_ATOMS_PER_GPU;
        }
        else
        {
            min_atoms_per_mpi_thread = MIN_ATOMS_PER_MPI_THREAD;
        }
    }

    /* Check if an algorithm does not support parallel simulation.  */
    if (nthreads_tmpi != 1 &&
        ( inputrec->eI == eiLBFGS ||
          inputrec->coulombtype == eelEWALD ) )
    {
        nthreads_tmpi = 1;

        md_print_warn(cr, fplog, "The integration or electrostatics algorithm doesn't support parallel runs. Using a single thread-MPI thread.\n");
        if (hw_opt->nthreads_tmpi > nthreads_tmpi)
        {
            gmx_fatal(FARGS, "You asked for more than 1 thread-MPI thread, but an algorithm doesn't support that");
        }
    }
    else if (mtop->natoms/nthreads_tmpi < min_atoms_per_mpi_thread)
    {
        /* the thread number was chosen automatically, but there are too many
           threads (too few atoms per thread) */
        nthreads_new = max(1, mtop->natoms/min_atoms_per_mpi_thread);

        /* Avoid partial use of Hyper-Threading */
        if (gmx_cpuid_x86_smt(hwinfo->cpuid_info) == GMX_CPUID_X86_SMT_ENABLED &&
            nthreads_new > nthreads_hw/2 && nthreads_new < nthreads_hw)
        {
            nthreads_new = nthreads_hw/2;
        }

        /* Avoid large prime numbers in the thread count */
        if (nthreads_new >= 6)
        {
            /* Use only 6,8,10 with additional factors of 2 */
            int fac;

            fac = 2;
            while (3*fac*2 <= nthreads_new)
            {
                fac *= 2;
            }

            nthreads_new = (nthreads_new/fac)*fac;
        }
        else
        {
            /* Avoid 5 */
            if (nthreads_new == 5)
            {
                nthreads_new = 4;
            }
        }

        nthreads_tmpi = nthreads_new;

        fprintf(stderr, "\n");
        fprintf(stderr, "NOTE: Parallelization is limited by the small number of atoms,\n");
        fprintf(stderr, "      only starting %d thread-MPI threads.\n", nthreads_tmpi);
        fprintf(stderr, "      You can use the -nt and/or -ntmpi option to optimize the number of threads.\n\n");
    }

    return nthreads_tmpi;
}
#endif /* GMX_THREAD_MPI */


/* Environment variable for setting nstlist */
static const char*  NSTLIST_ENVVAR          =  "GMX_NSTLIST";
/* Try to increase nstlist when using a GPU with nstlist less than this */
static const int    NSTLIST_GPU_ENOUGH      = 20;
/* Increase nstlist until the non-bonded cost increases more than this factor */
static const float  NBNXN_GPU_LIST_OK_FAC   = 1.25;
/* Don't increase nstlist beyond a non-bonded cost increases of this factor */
static const float  NBNXN_GPU_LIST_MAX_FAC  = 1.40;

/* Try to increase nstlist when running on a GPU */
static void increase_nstlist(FILE *fp, t_commrec *cr,
                             t_inputrec *ir, const gmx_mtop_t *mtop, matrix box)
{
    char                  *env;
    int                    nstlist_orig, nstlist_prev;
    verletbuf_list_setup_t ls;
    real                   rlist_inc, rlist_ok, rlist_max, rlist_new, rlist_prev;
    int                    i;
    t_state                state_tmp;
    gmx_bool               bBox, bDD, bCont;
    const char            *nstl_fmt = "\nFor optimal performance with a GPU nstlist (now %d) should be larger.\nThe optimum depends on your CPU and GPU resources.\nYou might want to try several nstlist values.\n";
    const char            *vbd_err  = "Can not increase nstlist for GPU run because verlet-buffer-drift is not set or used";
    const char            *box_err  = "Can not increase nstlist for GPU run because the box is too small";
    const char            *dd_err   = "Can not increase nstlist for GPU run because of domain decomposition limitations";
    char                   buf[STRLEN];

    /* Number of + nstlist alternative values to try when switching  */
    const int nstl[] = { 20, 25, 40, 50 };
#define NNSTL  sizeof(nstl)/sizeof(nstl[0])

    env = getenv(NSTLIST_ENVVAR);
    if (env == NULL)
    {
        if (fp != NULL)
        {
            fprintf(fp, nstl_fmt, ir->nstlist);
        }
    }

    if (ir->verletbuf_drift == 0)
    {
        gmx_fatal(FARGS, "You are using an old tpr file with a GPU, please generate a new tpr file with an up to date version of grompp");
    }

    if (ir->verletbuf_drift < 0)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n", vbd_err);
        }
        if (fp != NULL)
        {
            fprintf(fp, "%s\n", vbd_err);
        }

        return;
    }

    nstlist_orig = ir->nstlist;
    if (env != NULL)
    {
        sprintf(buf, "Getting nstlist from environment variable GMX_NSTLIST=%s", env);
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n", buf);
        }
        if (fp != NULL)
        {
            fprintf(fp, "%s\n", buf);
        }
        sscanf(env, "%d", &ir->nstlist);
    }

    verletbuf_get_list_setup(TRUE, &ls);

    /* Allow rlist to make the list double the size of the cut-off sphere */
    rlist_inc = nbnxn_get_rlist_effective_inc(NBNXN_GPU_CLUSTER_SIZE, mtop->natoms/det(box));
    rlist_ok  = (max(ir->rvdw, ir->rcoulomb) + rlist_inc)*pow(NBNXN_GPU_LIST_OK_FAC, 1.0/3.0) - rlist_inc;
    rlist_max = (max(ir->rvdw, ir->rcoulomb) + rlist_inc)*pow(NBNXN_GPU_LIST_MAX_FAC, 1.0/3.0) - rlist_inc;
    if (debug)
    {
        fprintf(debug, "GPU nstlist tuning: rlist_inc %.3f rlist_max %.3f\n",
                rlist_inc, rlist_max);
    }

    i            = 0;
    nstlist_prev = nstlist_orig;
    rlist_prev   = ir->rlist;
    do
    {
        if (env == NULL)
        {
            ir->nstlist = nstl[i];
        }

        /* Set the pair-list buffer size in ir */
        calc_verlet_buffer_size(mtop, det(box), ir, ir->verletbuf_drift, &ls,
                                NULL, &rlist_new);

        /* Does rlist fit in the box? */
        bBox = (sqr(rlist_new) < max_cutoff2(ir->ePBC, box));
        bDD  = TRUE;
        if (bBox && DOMAINDECOMP(cr))
        {
            /* Check if rlist fits in the domain decomposition */
            if (inputrec2nboundeddim(ir) < DIM)
            {
                gmx_incons("Changing nstlist with domain decomposition and unbounded dimensions is not implemented yet");
            }
            copy_mat(box, state_tmp.box);
            bDD = change_dd_cutoff(cr, &state_tmp, ir, rlist_new);
        }

        bCont = FALSE;

        if (env == NULL)
        {
            if (bBox && bDD && rlist_new <= rlist_max)
            {
                /* Increase nstlist */
                nstlist_prev = ir->nstlist;
                rlist_prev   = rlist_new;
                bCont        = (i+1 < NNSTL && rlist_new < rlist_ok);
            }
            else
            {
                /* Stick with the previous nstlist */
                ir->nstlist = nstlist_prev;
                rlist_new   = rlist_prev;
                bBox        = TRUE;
                bDD         = TRUE;
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
            fprintf(fp, "\n%s\n", bBox ? box_err : dd_err);
        }
        ir->nstlist = nstlist_orig;
    }
    else if (ir->nstlist != nstlist_orig || rlist_new != ir->rlist)
    {
        sprintf(buf, "Changing nstlist from %d to %d, rlist from %g to %g",
                nstlist_orig, ir->nstlist,
                ir->rlist, rlist_new);
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n\n", buf);
        }
        if (fp != NULL)
        {
            fprintf(fp, "%s\n\n", buf);
        }
        ir->rlist     = rlist_new;
        ir->rlistlong = rlist_new;
    }
}

static void prepare_verlet_scheme(FILE                *fplog,
                                  const gmx_hw_info_t *hwinfo,
                                  t_commrec           *cr,
                                  const char          *nbpu_opt,
                                  t_inputrec          *ir,
                                  const gmx_mtop_t    *mtop,
                                  matrix               box,
                                  gmx_bool            *bUseGPU)
{
    /* Here we only check for GPU usage on the MPI master process,
     * as here we don't know how many GPUs we will use yet.
     * We check for a GPU on all processes later.
     */
    *bUseGPU = hwinfo->bCanUseGPU || (getenv("GMX_EMULATE_GPU") != NULL);

    if (ir->verletbuf_drift > 0)
    {
        /* Update the Verlet buffer size for the current run setup */
        verletbuf_list_setup_t ls;
        real                   rlist_new;

        /* Here we assume CPU acceleration is on. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        verletbuf_get_list_setup(*bUseGPU, &ls);

        calc_verlet_buffer_size(mtop, det(box), ir,
                                ir->verletbuf_drift, &ls,
                                NULL, &rlist_new);
        if (rlist_new != ir->rlist)
        {
            if (fplog != NULL)
            {
                fprintf(fplog, "\nChanging rlist from %g to %g for non-bonded %dx%d atom kernels\n\n",
                        ir->rlist, rlist_new,
                        ls.cluster_size_i, ls.cluster_size_j);
            }
            ir->rlist     = rlist_new;
            ir->rlistlong = rlist_new;
        }
    }

    /* With GPU or emulation we should check nstlist for performance */
    if ((EI_DYNAMICS(ir->eI) &&
         *bUseGPU &&
         ir->nstlist < NSTLIST_GPU_ENOUGH) ||
        getenv(NSTLIST_ENVVAR) != NULL)
    {
        /* Choose a better nstlist */
        increase_nstlist(fplog, cr, ir, mtop, box);
    }
}

static void convert_to_verlet_scheme(FILE *fplog,
                                     t_inputrec *ir,
                                     gmx_mtop_t *mtop, real box_vol)
{
    char *conv_mesg = "Converting input file with group cut-off scheme to the Verlet cut-off scheme";

    md_print_warn(NULL, fplog, "%s\n", conv_mesg);

    ir->cutoff_scheme   = ecutsVERLET;
    ir->verletbuf_drift = 0.005;

    if (ir->rcoulomb != ir->rvdw)
    {
        gmx_fatal(FARGS, "The VdW and Coulomb cut-offs are different, whereas the Verlet scheme only supports equal cut-offs");
    }

    if (ir->vdwtype == evdwUSER || EEL_USER(ir->coulombtype))
    {
        gmx_fatal(FARGS, "User non-bonded potentials are not (yet) supported with the Verlet scheme");
    }
    else if (EVDW_SWITCHED(ir->vdwtype) || EEL_SWITCHED(ir->coulombtype))
    {
        md_print_warn(NULL, fplog, "Converting switched or shifted interactions to a shifted potential (without force shift), this will lead to slightly different interaction potentials");

        if (EVDW_SWITCHED(ir->vdwtype))
        {
            ir->vdwtype = evdwCUT;
        }
        if (EEL_SWITCHED(ir->coulombtype))
        {
            if (EEL_FULL(ir->coulombtype))
            {
                /* With full electrostatic only PME can be switched */
                ir->coulombtype = eelPME;
            }
            else
            {
                md_print_warn(NULL, fplog, "NOTE: Replacing %s electrostatics with reaction-field with epsilon-rf=inf\n", eel_names[ir->coulombtype]);
                ir->coulombtype = eelRF;
                ir->epsilon_rf  = 0.0;
            }
        }

        /* We set the target energy drift to a small number.
         * Note that this is only for testing. For production the user
         * should think about this and set the mdp options.
         */
        ir->verletbuf_drift = 1e-4;
    }

    if (inputrec2nboundeddim(ir) != 3)
    {
        gmx_fatal(FARGS, "Can only convert old tpr files to the Verlet cut-off scheme with 3D pbc");
    }

    if (ir->efep != efepNO || ir->implicit_solvent != eisNO)
    {
        gmx_fatal(FARGS, "Will not convert old tpr files to the Verlet cut-off scheme with free-energy calculations or implicit solvent");
    }

    if (EI_DYNAMICS(ir->eI) && !(EI_MD(ir->eI) && ir->etc == etcNO))
    {
        verletbuf_list_setup_t ls;

        verletbuf_get_list_setup(FALSE, &ls);
        calc_verlet_buffer_size(mtop, box_vol, ir, ir->verletbuf_drift, &ls,
                                NULL, &ir->rlist);
    }
    else
    {
        ir->verletbuf_drift = -1;
        ir->rlist           = 1.05*max(ir->rvdw, ir->rcoulomb);
    }

    gmx_mtop_remove_chargegroups(mtop);
}

static void check_and_update_hw_opt(gmx_hw_opt_t *hw_opt,
                                    int           cutoff_scheme,
                                    gmx_bool      bIsSimMaster)
{
    gmx_omp_nthreads_read_env(&hw_opt->nthreads_omp, bIsSimMaster);

#ifndef GMX_THREAD_MPI
    if (hw_opt->nthreads_tot > 0)
    {
        gmx_fatal(FARGS, "Setting the total number of threads is only supported with thread-MPI and Gromacs was compiled without thread-MPI");
    }
    if (hw_opt->nthreads_tmpi > 0)
    {
        gmx_fatal(FARGS, "Setting the number of thread-MPI threads is only supported with thread-MPI and Gromacs was compiled without thread-MPI");
    }
#endif

#ifndef GMX_OPENMP
    if (hw_opt->nthreads_omp > 1)
    {
        gmx_fatal(FARGS, "More than 1 OpenMP thread requested, but Gromacs was compiled without OpenMP support");
    }
    hw_opt->nthreads_omp = 1;
#endif

    if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp_pme <= 0)
    {
        /* We have the same number of OpenMP threads for PP and PME processes,
         * thus we can perform several consistency checks.
         */
        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_omp > 0 &&
            hw_opt->nthreads_tot != hw_opt->nthreads_tmpi*hw_opt->nthreads_omp)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) does not match the thread-MPI threads (%d) times the OpenMP threads (%d) requested",
                      hw_opt->nthreads_tot, hw_opt->nthreads_tmpi, hw_opt->nthreads_omp);
        }

        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_tot % hw_opt->nthreads_tmpi != 0)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) is not divisible by the number of thread-MPI threads requested (%d)",
                      hw_opt->nthreads_tot, hw_opt->nthreads_tmpi);
        }

        if (hw_opt->nthreads_omp > 0 &&
            hw_opt->nthreads_tot % hw_opt->nthreads_omp != 0)
        {
            gmx_fatal(FARGS, "The total number of threads requested (%d) is not divisible by the number of OpenMP threads requested (%d)",
                      hw_opt->nthreads_tot, hw_opt->nthreads_omp);
        }

        if (hw_opt->nthreads_tmpi > 0 &&
            hw_opt->nthreads_omp <= 0)
        {
            hw_opt->nthreads_omp = hw_opt->nthreads_tot/hw_opt->nthreads_tmpi;
        }
    }

#ifndef GMX_OPENMP
    if (hw_opt->nthreads_omp > 1)
    {
        gmx_fatal(FARGS, "OpenMP threads are requested, but Gromacs was compiled without OpenMP support");
    }
#endif

    if (cutoff_scheme == ecutsGROUP)
    {
        /* We only have OpenMP support for PME only nodes */
        if (hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS, "OpenMP threads have been requested with cut-off scheme %s, but these are only supported with cut-off scheme %s",
                      ecutscheme_names[cutoff_scheme],
                      ecutscheme_names[ecutsVERLET]);
        }
        hw_opt->nthreads_omp = 1;
    }

    if (hw_opt->nthreads_omp_pme > 0 && hw_opt->nthreads_omp <= 0)
    {
        gmx_fatal(FARGS, "You need to specify -ntomp in addition to -ntomp_pme");
    }

    if (hw_opt->nthreads_tot == 1)
    {
        hw_opt->nthreads_tmpi = 1;

        if (hw_opt->nthreads_omp > 1)
        {
            gmx_fatal(FARGS, "You requested %d OpenMP threads with %d total threads",
                      hw_opt->nthreads_tmpi, hw_opt->nthreads_tot);
        }
        hw_opt->nthreads_omp = 1;
    }

    if (hw_opt->nthreads_omp_pme <= 0 && hw_opt->nthreads_omp > 0)
    {
        hw_opt->nthreads_omp_pme = hw_opt->nthreads_omp;
    }

    if (debug)
    {
        fprintf(debug, "hw_opt: nt %d ntmpi %d ntomp %d ntomp_pme %d gpu_id '%s'\n",
                hw_opt->nthreads_tot,
                hw_opt->nthreads_tmpi,
                hw_opt->nthreads_omp,
                hw_opt->nthreads_omp_pme,
                hw_opt->gpu_id != NULL ? hw_opt->gpu_id : "");

    }
}


/* Override the value in inputrec with value passed on the command line (if any) */
static void override_nsteps_cmdline(FILE            *fplog,
                                    gmx_large_int_t  nsteps_cmdline,
                                    t_inputrec      *ir,
                                    const t_commrec *cr)
{
    char sbuf[STEPSTRSIZE];

    assert(ir);
    assert(cr);

    /* override with anything else than the default -2 */
    if (nsteps_cmdline > -2)
    {
        char stmp[STRLEN];

        ir->nsteps = nsteps_cmdline;
        if (EI_DYNAMICS(ir->eI))
        {
            sprintf(stmp, "Overriding nsteps with value passed on the command line: %s steps, %.3f ps",
                    gmx_step_str(nsteps_cmdline, sbuf),
                    nsteps_cmdline*ir->delta_t);
        }
        else
        {
            sprintf(stmp, "Overriding nsteps with value passed on the command line: %s steps",
                    gmx_step_str(nsteps_cmdline, sbuf));
        }

        md_print_warn(cr, fplog, "%s\n", stmp);
    }
}

/* Data structure set by SIMMASTER which needs to be passed to all nodes
 * before the other nodes have read the tpx file and called gmx_detect_hardware.
 */
typedef struct {
    int      cutoff_scheme; /* The cutoff scheme from inputrec_t */
    gmx_bool bUseGPU;       /* Use GPU or GPU emulation          */
} master_inf_t;

int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, t_commrec *cr, int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm,
             ivec ddxyz, int dd_node_order, real rdd, real rconstr,
             const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt,
             gmx_large_int_t nsteps_cmdline, int nstepout, int resetstep,
             int nmultisim, int repl_ex_nst, int repl_ex_nex,
             int repl_ex_seed, real pforce, real cpt_period, real max_hours,
             const char *deviceOptions, unsigned long Flags)
{
    gmx_bool        bForceUseGPU, bTryUseGPU;
    double          nodetime = 0, realtime;
    t_inputrec     *inputrec;
    t_state        *state = NULL;
    matrix          box;
    gmx_ddbox_t     ddbox = {0};
    int             npme_major, npme_minor;
    real            tmpr1, tmpr2;
    t_nrnb         *nrnb;
    gmx_mtop_t     *mtop       = NULL;
    t_mdatoms      *mdatoms    = NULL;
    t_forcerec     *fr         = NULL;
    t_fcdata       *fcd        = NULL;
    real            ewaldcoeff = 0;
    gmx_pme_t      *pmedata    = NULL;
    gmx_vsite_t    *vsite      = NULL;
    gmx_constr_t    constr;
    int             i, m, nChargePerturbed = -1, status, nalloc;
    char           *gro;
    gmx_wallcycle_t wcycle;
    gmx_bool        bReadRNG, bReadEkin;
    int             list;
    gmx_runtime_t   runtime;
    int             rc;
    gmx_large_int_t reset_counters;
    gmx_edsam_t     ed           = NULL;
    t_commrec      *cr_old       = cr;
    int             nthreads_pme = 1;
    int             nthreads_pp  = 1;
    gmx_membed_t    membed       = NULL;
    gmx_hw_info_t  *hwinfo       = NULL;
    master_inf_t    minf         = {-1, FALSE};

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    snew(inputrec, 1);
    snew(mtop, 1);

    if (Flags & MD_APPENDFILES)
    {
        fplog = NULL;
    }

    bForceUseGPU = (strncmp(nbpu_opt, "gpu", 3) == 0);
    bTryUseGPU   = (strncmp(nbpu_opt, "auto", 4) == 0) || bForceUseGPU;

    /* Detect hardware, gather information. This is an operation that is
     * global for this process (MPI rank). */
    hwinfo = gmx_detect_hardware(fplog, cr,
                                 bForceUseGPU, bTryUseGPU, hw_opt->gpu_id);


    snew(state, 1);
    if (SIMMASTER(cr))
    {
        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPX, nfile, fnm), inputrec, state, NULL, mtop);

        if (inputrec->cutoff_scheme != ecutsVERLET &&
            ((Flags & MD_TESTVERLET) || getenv("GMX_VERLET_SCHEME") != NULL))
        {
            convert_to_verlet_scheme(fplog, inputrec, mtop, det(state->box));
        }


        minf.cutoff_scheme = inputrec->cutoff_scheme;
        minf.bUseGPU       = FALSE;

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            prepare_verlet_scheme(fplog, hwinfo, cr, nbpu_opt,
                                  inputrec, mtop, state->box,
                                  &minf.bUseGPU);
        }
        else if (hwinfo->bCanUseGPU)
        {
            md_print_warn(cr, fplog,
                          "NOTE: GPU(s) found, but the current simulation can not use GPUs\n"
                          "      To use a GPU, set the mdp option: cutoff-scheme = Verlet\n"
                          "      (for quick performance testing you can use the -testverlet option)\n");

            if (bForceUseGPU)
            {
                gmx_fatal(FARGS, "GPU requested, but can't be used without cutoff-scheme=Verlet");
            }
        }
#ifdef GMX_IS_BGQ
        else
        {
            md_print_warn(cr, fplog,
                          "NOTE: There is no SIMD implementation of the group scheme kernels on\n"
                          "      BlueGene/Q. You will observe better performance from using the\n"
                          "      Verlet cut-off scheme.\n");
        }
#endif
    }
#ifndef GMX_THREAD_MPI
    if (PAR(cr))
    {
        gmx_bcast_sim(sizeof(minf), &minf, cr);
    }
#endif
    if (minf.bUseGPU && cr->npmenodes == -1)
    {
        /* Don't automatically use PME-only nodes with GPUs */
        cr->npmenodes = 0;
    }

    /* Check for externally set OpenMP affinity and turn off internal
     * pinning if any is found. We need to do this check early to tell
     * thread-MPI whether it should do pinning when spawning threads.
     * TODO: the above no longer holds, we should move these checks down
     */
    gmx_omp_check_thread_affinity(fplog, cr, hw_opt);

#ifdef GMX_THREAD_MPI
    /* With thread-MPI inputrec is only set here on the master thread */
    if (SIMMASTER(cr))
#endif
    {
        check_and_update_hw_opt(hw_opt, minf.cutoff_scheme, SIMMASTER(cr));

#ifdef GMX_THREAD_MPI
        /* Early check for externally set process affinity. Can't do over all
         * MPI processes because hwinfo is not available everywhere, but with
         * thread-MPI it's needed as pinning might get turned off which needs
         * to be known before starting thread-MPI. */
        gmx_check_thread_affinity_set(fplog,
                                      NULL,
                                      hw_opt, hwinfo->nthreads_hw_avail, FALSE);
#endif

#ifdef GMX_THREAD_MPI
        if (cr->npmenodes > 0 && hw_opt->nthreads_tmpi <= 0)
        {
            gmx_fatal(FARGS, "You need to explicitly specify the number of MPI threads (-ntmpi) when using separate PME nodes");
        }
#endif

        if (hw_opt->nthreads_omp_pme != hw_opt->nthreads_omp &&
            cr->npmenodes <= 0)
        {
            gmx_fatal(FARGS, "You need to explicitly specify the number of PME nodes (-npme) when using different number of OpenMP threads for PP and PME nodes");
        }
    }

#ifdef GMX_THREAD_MPI
    if (SIMMASTER(cr))
    {
        /* NOW the threads will be started: */
        hw_opt->nthreads_tmpi = get_nthreads_mpi(hwinfo,
                                                 hw_opt,
                                                 inputrec, mtop,
                                                 cr, fplog);
        if (hw_opt->nthreads_tot > 0 && hw_opt->nthreads_omp <= 0)
        {
            hw_opt->nthreads_omp = hw_opt->nthreads_tot/hw_opt->nthreads_tmpi;
        }

        if (hw_opt->nthreads_tmpi > 1)
        {
            /* now start the threads. */
            cr = mdrunner_start_threads(hw_opt, fplog, cr_old, nfile, fnm,
                                        oenv, bVerbose, bCompact, nstglobalcomm,
                                        ddxyz, dd_node_order, rdd, rconstr,
                                        dddlb_opt, dlb_scale, ddcsx, ddcsy, ddcsz,
                                        nbpu_opt,
                                        nsteps_cmdline, nstepout, resetstep, nmultisim,
                                        repl_ex_nst, repl_ex_nex, repl_ex_seed, pforce,
                                        cpt_period, max_hours, deviceOptions,
                                        Flags);
            /* the main thread continues here with a new cr. We don't deallocate
               the old cr because other threads may still be reading it. */
            if (cr == NULL)
            {
                gmx_comm("Failed to spawn threads");
            }
        }
    }
#endif
    /* END OF CAUTION: cr is now reliable */

    /* g_membed initialisation *
     * Because we change the mtop, init_membed is called before the init_parallel *
     * (in case we ever want to make it run in parallel) */
    if (opt2bSet("-membed", nfile, fnm))
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "Initializing membed");
        }
        membed = init_membed(fplog, nfile, fnm, mtop, inputrec, state, cr, &cpt_period);
    }

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(fplog, cr, inputrec, mtop);

        /* This check needs to happen after get_nthreads_mpi() */
        if (inputrec->cutoff_scheme == ecutsVERLET && (Flags & MD_PARTDEC))
        {
            gmx_fatal_collective(FARGS, cr, NULL,
                                 "The Verlet cut-off scheme is not supported with particle decomposition.\n"
                                 "You can achieve the same effect as particle decomposition by running in parallel using only OpenMP threads.");
        }
    }
    if (fplog != NULL)
    {
        pr_inputrec(fplog, 0, "Input Parameters", inputrec, FALSE);
    }

    /* now make sure the state is initialized and propagated */
    set_state_entries(state, inputrec, cr->nnodes);

    /* A parallel command line option consistency check that we can
       only do after any threads have started. */
    if (!PAR(cr) &&
        (ddxyz[XX] > 1 || ddxyz[YY] > 1 || ddxyz[ZZ] > 1 || cr->npmenodes > 0))
    {
        gmx_fatal(FARGS,
                  "The -dd or -npme option request a parallel simulation, "
#ifndef GMX_MPI
                  "but %s was compiled without threads or MPI enabled"
#else
#ifdef GMX_THREAD_MPI
                  "but the number of threads (option -nt) is 1"
#else
                  "but %s was not started through mpirun/mpiexec or only one process was requested through mpirun/mpiexec"
#endif
#endif
                  , ShortProgram()
                  );
    }

    if ((Flags & MD_RERUN) &&
        (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
    }

    if (can_use_allvsall(inputrec, mtop, TRUE, cr, fplog) && PAR(cr))
    {
        /* Simple neighbour searching and (also?) all-vs-all loops
         * do not work with domain decomposition. */
        Flags |= MD_PARTDEC;
    }

    if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC))
    {
        if (cr->npmenodes > 0)
        {
            if (!EEL_PME(inputrec->coulombtype))
            {
                gmx_fatal_collective(FARGS, cr, NULL,
                                     "PME nodes are requested, but the system does not use PME electrostatics");
            }
            if (Flags & MD_PARTDEC)
            {
                gmx_fatal_collective(FARGS, cr, NULL,
                                     "PME nodes are requested, but particle decomposition does not support separate PME nodes");
            }
        }

        cr->npmenodes = 0;
    }

#ifdef GMX_FAHCORE
    fcRegisterSteps(inputrec->nsteps, inputrec->init_step);
#endif

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */
    snew(fcd, 1);

    /* This needs to be called before read_checkpoint to extend the state */
    init_disres(fplog, mtop, inputrec, cr, Flags & MD_PARTDEC, fcd, state, repl_ex_nst > 0);

    if (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0)
    {
        if (PAR(cr) && !(Flags & MD_PARTDEC))
        {
            gmx_fatal(FARGS, "Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)");
        }
        /* Orientation restraints */
        if (MASTER(cr))
        {
            init_orires(fplog, mtop, state->x, inputrec, cr->ms, &(fcd->orires),
                        state);
        }
    }

    if (DEFORM(*inputrec))
    {
        /* Store the deform reference box before reading the checkpoint */
        if (SIMMASTER(cr))
        {
            copy_mat(state->box, box);
        }
        if (PAR(cr))
        {
            gmx_bcast(sizeof(box), box, cr);
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
        copy_mat(box, deform_init_box_tpx);
#ifdef GMX_THREAD_MPI
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
    }

    if (opt2bSet("-cpi", nfile, fnm))
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        if (gmx_fexist_master(opt2fn_master("-cpi", nfile, fnm, cr), cr) )
        {
            load_checkpoint(opt2fn_master("-cpi", nfile, fnm, cr), &fplog,
                            cr, Flags & MD_PARTDEC, ddxyz,
                            inputrec, state, &bReadRNG, &bReadEkin,
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
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr, !(Flags & MD_SEPPOT),
                     Flags, &fplog);
    }

    /* override nsteps with value from cmdline */
    override_nsteps_cmdline(fplog, nsteps_cmdline, inputrec, cr);

    if (SIMMASTER(cr))
    {
        copy_mat(state->box, box);
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(box), box, cr);
    }

    /* Essential dynamics */
    if (opt2bSet("-ei", nfile, fnm))
    {
        /* Open input and output files, allocate space for ED data structure */
        ed = ed_open(mtop->natoms, &state->edsamstate, nfile, fnm, Flags, oenv, cr);
    }

    if (PAR(cr) && !((Flags & MD_PARTDEC) ||
                     EI_TPI(inputrec->eI) ||
                     inputrec->eI == eiNM))
    {
        cr->dd = init_domain_decomposition(fplog, cr, Flags, ddxyz, rdd, rconstr,
                                           dddlb_opt, dlb_scale,
                                           ddcsx, ddcsy, ddcsz,
                                           mtop, inputrec,
                                           box, state->x,
                                           &ddbox, &npme_major, &npme_minor);

        make_dd_communicators(fplog, cr, dd_node_order);

        /* Set overallocation to avoid frequent reallocation of arrays */
        set_over_alloc_dd(TRUE);
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty      = (DUTY_PP | DUTY_PME);
        npme_major    = 1;
        npme_minor    = 1;
        /* NM and TPI perform single node energy calculations in parallel */
        if (!(inputrec->eI == eiNM || EI_TPI(inputrec->eI)))
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
        gmx_setup_nodecomm(fplog, cr);
    }

    /* Initialize per-physical-node MPI process/thread ID and counters. */
    gmx_init_intranode_counters(cr);

#ifdef GMX_MPI
    md_print_info(cr, fplog, "Using %d MPI %s\n",
                  cr->nnodes,
#ifdef GMX_THREAD_MPI
                  cr->nnodes == 1 ? "thread" : "threads"
#else
                  cr->nnodes == 1 ? "process" : "processes"
#endif
                  );
    fflush(stderr);
#endif

    gmx_omp_nthreads_init(fplog, cr,
                          hwinfo->nthreads_hw_avail,
                          hw_opt->nthreads_omp,
                          hw_opt->nthreads_omp_pme,
                          (cr->duty & DUTY_PP) == 0,
                          inputrec->cutoff_scheme == ecutsVERLET);

    /* check consistency and decide on the number of gpus to use. */
    gmx_check_hw_runconf_consistency(fplog, hwinfo, cr, hw_opt->nthreads_tmpi,
                                     minf.bUseGPU);

    /* getting number of PP/PME threads
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    /* TODO nthreads_pp is only used for pinning threads.
     * This is a temporary solution until we have a hw topology library.
     */
    nthreads_pp  = gmx_omp_nthreads_get(emntNonbonded);
    nthreads_pme = gmx_omp_nthreads_get(emntPME);

    wcycle = wallcycle_init(fplog, resetstep, cr, nthreads_pp, nthreads_pme);

    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast_sim(sizeof(reset_counters), &reset_counters, cr);
        wcycle_set_reset_counters(wcycle, reset_counters);
    }

    snew(nrnb, 1);
    if (cr->duty & DUTY_PP)
    {
        /* For domain decomposition we allocate dynamically
         * in dd_partition_system.
         */
        if (DOMAINDECOMP(cr))
        {
            bcast_state_setup(cr, state);
        }
        else
        {
            if (PAR(cr))
            {
                bcast_state(cr, state, TRUE);
            }
        }

        /* Initiate forcerecord */
        fr         = mk_forcerec();
        fr->hwinfo = hwinfo;
        init_forcerec(fplog, oenv, fr, fcd, inputrec, mtop, cr, box, FALSE,
                      opt2fn("-table", nfile, fnm),
                      opt2fn("-tabletf", nfile, fnm),
                      opt2fn("-tablep", nfile, fnm),
                      opt2fn("-tableb", nfile, fnm),
                      nbpu_opt,
                      FALSE, pforce);

        /* version for PCA_NOT_READ_NODE (see md.c) */
        /*init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
           "nofile","nofile","nofile","nofile",FALSE,pforce);
         */
        fr->bSepDVDL = ((Flags & MD_SEPPOT) == MD_SEPPOT);

        /* Initialize QM-MM */
        if (fr->bQMMM)
        {
            init_QMMMrec(cr, box, mtop, inputrec, fr);
        }

        /* Initialize the mdatoms structure.
         * mdatoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdatoms = init_mdatoms(fplog, mtop, inputrec->efep != efepNO);

        if (mdatoms->nPerturbed > 0 && inputrec->cutoff_scheme == ecutsVERLET)
        {
            gmx_fatal(FARGS, "The Verlet cut-off scheme does not (yet) support free-energy calculations with perturbed atoms, only perturbed interactions. This will be implemented soon. Use the group scheme for now.");
        }

        /* Initialize the virtual site communication */
        vsite = init_vsite(mtop, cr, FALSE);

        calc_shifts(box, fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(cr) &&
            !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->ePBC != epbcNONE)
            {
                do_pbc_first_mtop(fplog, inputrec->ePBC, box, mtop, state->x);
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                construct_vsites_mtop(fplog, vsite, mtop, state->x);
            }
        }

        if (EEL_PME(fr->eeltype))
        {
            ewaldcoeff = fr->ewaldcoeff;
            pmedata    = &fr->pmedata;
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
        snew(pmedata, 1);
    }

    if (hw_opt->thread_affinity != threadaffOFF)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(fplog, cr,
                                      hw_opt, hwinfo->nthreads_hw_avail, TRUE);

        /* Set the CPU affinity */
        gmx_set_thread_affinity(fplog, cr, hw_opt, nthreads_pme, hwinfo,
                                inputrec);
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
            gmx_bcast_sim(sizeof(nChargePerturbed), &nChargePerturbed, cr);
        }

        if (cr->duty & DUTY_PME)
        {
            status = gmx_pme_init(pmedata, cr, npme_major, npme_minor, inputrec,
                                  mtop ? mtop->natoms : 0, nChargePerturbed,
                                  (Flags & MD_REPRODUCIBLE), nthreads_pme);
            if (status != 0)
            {
                gmx_fatal(FARGS, "Error %d initializing PME", status);
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
        signal_handler_install();
    }

    if (cr->duty & DUTY_PP)
    {
        if (inputrec->ePull != epullNO)
        {
            /* Initialize pull code */
            init_pull(fplog, inputrec, nfile, fnm, mtop, cr, oenv, inputrec->fepvals->init_lambda,
                      EI_DYNAMICS(inputrec->eI) && MASTER(cr), Flags);
        }

        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            init_rot(fplog, inputrec, nfile, fnm, cr, state->x, box, mtop, oenv,
                     bVerbose, Flags);
        }

        constr = init_constraints(fplog, mtop, inputrec, ed, state, cr);

        if (DOMAINDECOMP(cr))
        {
            dd_init_bondeds(fplog, cr->dd, mtop, vsite, constr, inputrec,
                            Flags & MD_DDBONDCHECK, fr->cginfo_mb);

            set_dd_parameters(fplog, cr->dd, dlb_scale, inputrec, fr, &ddbox);

            setup_dd_grid(fplog, cr->dd);
        }

        /* Now do whatever the user wants us to do (how flexible...) */
        integrator[inputrec->eI].func(fplog, cr, nfile, fnm,
                                      oenv, bVerbose, bCompact,
                                      nstglobalcomm,
                                      vsite, constr,
                                      nstepout, inputrec, mtop,
                                      fcd, state,
                                      mdatoms, nrnb, wcycle, ed, fr,
                                      repl_ex_nst, repl_ex_nex, repl_ex_seed,
                                      membed,
                                      cpt_period, max_hours,
                                      deviceOptions,
                                      Flags,
                                      &runtime);

        if (inputrec->ePull != epullNO)
        {
            finish_pull(fplog, inputrec->pull);
        }

        if (inputrec->bRot)
        {
            finish_rot(inputrec->rot);
        }

    }
    else
    {
        /* do PME only */
        gmx_pmeonly(*pmedata, cr, nrnb, wcycle, &runtime, ewaldcoeff, FALSE, inputrec);
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

    wallcycle_stop(wcycle, ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog, cr, ftp2fn(efSTO, nfile, fnm),
               inputrec, nrnb, wcycle, &runtime,
               fr != NULL && fr->nbv != NULL && fr->nbv->bUseGPU ?
               nbnxn_cuda_get_timings(fr->nbv->cu_nbv) : NULL,
               nthreads_pp,
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

    if ((cr->duty & DUTY_PP) && fr->nbv != NULL && fr->nbv->bUseGPU)
    {
        char gpu_err_str[STRLEN];

        /* free GPU memory and uninitialize GPU (by destroying the context) */
        nbnxn_cuda_free(fplog, fr->nbv->cu_nbv);

        if (!free_gpu(gpu_err_str))
        {
            gmx_warning("On node %d failed to free GPU #%d: %s",
                        cr->nodeid, get_current_gpu_device_id(), gpu_err_str);
        }
    }

    if (opt2bSet("-membed", nfile, fnm))
    {
        sfree(membed);
    }

    gmx_hardware_info_free(hwinfo);

    /* Does what it says */
    print_date_and_time(fplog, cr->nodeid, "Finished mdrun", &runtime);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_close(fplog);
    }

    rc = (int)gmx_get_stop_condition();

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
