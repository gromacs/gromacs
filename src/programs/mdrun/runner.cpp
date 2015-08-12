/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "config.h"

#include <assert.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/gpu_utils/gpu_utils.h"
#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/disre.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_detect_hardware.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/gmx_thread_affinity.h"
#include "gromacs/legacyheaders/inputrec.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/md_support.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/orires.h"
#include "gromacs/legacyheaders/qmmm.h"
#include "gromacs/legacyheaders/sighandler.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/calculate-ewald-splitting-coefficient.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#include "deform.h"
#include "membed.h"
#include "repl_ex.h"
#include "resource-division.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

typedef struct {
    gmx_integrator_t *func;
} gmx_intp_t;

/* The array should match the eI array in include/types/enums.h */
const gmx_intp_t    integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md}, {do_md}, {do_md}};

gmx_int64_t         deform_init_init_step_tpx;
matrix              deform_init_box_tpx;
tMPI_Thread_mutex_t deform_init_box_mutex = TMPI_THREAD_MUTEX_INITIALIZER;


#ifdef GMX_THREAD_MPI
/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
#define MIN_ATOMS_PER_MPI_THREAD    90
#define MIN_ATOMS_PER_GPU           900

struct mdrunner_arglist
{
    gmx_hw_opt_t    hw_opt;
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
    int             nstlist_cmdline;
    gmx_int64_t     nsteps_cmdline;
    int             nstepout;
    int             resetstep;
    int             nmultisim;
    int             repl_ex_nst;
    int             repl_ex_nex;
    int             repl_ex_seed;
    real            pforce;
    real            cpt_period;
    real            max_hours;
    int             imdport;
    unsigned long   Flags;
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

    cr = reinitialize_commrec_for_this_thread(mc.cr);

    if (MASTER(cr))
    {
        fplog = mc.fplog;
    }

    mdrunner(&mc.hw_opt, fplog, cr, mc.nfile, fnm, mc.oenv,
             mc.bVerbose, mc.bCompact, mc.nstglobalcomm,
             mc.ddxyz, mc.dd_node_order, mc.rdd,
             mc.rconstr, mc.dddlb_opt, mc.dlb_scale,
             mc.ddcsx, mc.ddcsy, mc.ddcsz,
             mc.nbpu_opt, mc.nstlist_cmdline,
             mc.nsteps_cmdline, mc.nstepout, mc.resetstep,
             mc.nmultisim, mc.repl_ex_nst, mc.repl_ex_nex, mc.repl_ex_seed, mc.pforce,
             mc.cpt_period, mc.max_hours, mc.imdport, mc.Flags);
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
                                         const char *nbpu_opt, int nstlist_cmdline,
                                         gmx_int64_t nsteps_cmdline,
                                         int nstepout, int resetstep,
                                         int nmultisim, int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
                                         real pforce, real cpt_period, real max_hours,
                                         unsigned long Flags)
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
    /* hw_opt contains pointers, which should all be NULL at this stage */
    mda->hw_opt          = *hw_opt;
    mda->fplog           = fplog;
    mda->cr              = cr;
    mda->nfile           = nfile;
    mda->fnm             = fnmn;
    mda->oenv            = oenv;
    mda->bVerbose        = bVerbose;
    mda->bCompact        = bCompact;
    mda->nstglobalcomm   = nstglobalcomm;
    mda->ddxyz[XX]       = ddxyz[XX];
    mda->ddxyz[YY]       = ddxyz[YY];
    mda->ddxyz[ZZ]       = ddxyz[ZZ];
    mda->dd_node_order   = dd_node_order;
    mda->rdd             = rdd;
    mda->rconstr         = rconstr;
    mda->dddlb_opt       = dddlb_opt;
    mda->dlb_scale       = dlb_scale;
    mda->ddcsx           = ddcsx;
    mda->ddcsy           = ddcsy;
    mda->ddcsz           = ddcsz;
    mda->nbpu_opt        = nbpu_opt;
    mda->nstlist_cmdline = nstlist_cmdline;
    mda->nsteps_cmdline  = nsteps_cmdline;
    mda->nstepout        = nstepout;
    mda->resetstep       = resetstep;
    mda->nmultisim       = nmultisim;
    mda->repl_ex_nst     = repl_ex_nst;
    mda->repl_ex_nex     = repl_ex_nex;
    mda->repl_ex_seed    = repl_ex_seed;
    mda->pforce          = pforce;
    mda->cpt_period      = cpt_period;
    mda->max_hours       = max_hours;
    mda->Flags           = Flags;

    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns, we set thread affinity later */
    ret = tMPI_Init_fn(TRUE, hw_opt->nthreads_tmpi, TMPI_AFFINITY_NONE,
                       mdrunner_start_fn, (void*)(mda) );
    if (ret != TMPI_SUCCESS)
    {
        return NULL;
    }

    crn = reinitialize_commrec_for_this_thread(cr);
    return crn;
}

#endif /* GMX_THREAD_MPI */


/* We determine the extra cost of the non-bonded kernels compared to
 * a reference nstlist value of 10 (which is the default in grompp).
 */
static const int    nbnxnReferenceNstlist = 10;
/* The values to try when switching  */
const int           nstlist_try[] = { 20, 25, 40 };
#define NNSTL  sizeof(nstlist_try)/sizeof(nstlist_try[0])
/* Increase nstlist until the non-bonded cost increases more than listfac_ok,
 * but never more than listfac_max.
 * A standard (protein+)water system at 300K with PME ewald_rtol=1e-5
 * needs 1.28 at rcoulomb=0.9 and 1.24 at rcoulomb=1.0 to get to nstlist=40.
 * Note that both CPU and GPU factors are conservative. Performance should
 * not go down due to this tuning, except with a relatively slow GPU.
 * On the other hand, at medium/high parallelization or with fast GPUs
 * nstlist will not be increased enough to reach optimal performance.
 */
/* CPU: pair-search is about a factor 1.5 slower than the non-bonded kernel */
static const float  nbnxn_cpu_listfac_ok    = 1.05;
static const float  nbnxn_cpu_listfac_max   = 1.09;
/* GPU: pair-search is a factor 1.5-3 slower than the non-bonded kernel */
static const float  nbnxn_gpu_listfac_ok    = 1.20;
static const float  nbnxn_gpu_listfac_max   = 1.30;

/* Try to increase nstlist when using the Verlet cut-off scheme */
static void increase_nstlist(FILE *fp, t_commrec *cr,
                             t_inputrec *ir, int nstlist_cmdline,
                             const gmx_mtop_t *mtop, matrix box,
                             gmx_bool bGPU)
{
    float                  listfac_ok, listfac_max;
    int                    nstlist_orig, nstlist_prev;
    verletbuf_list_setup_t ls;
    real                   rlistWithReferenceNstlist, rlist_inc, rlist_ok, rlist_max;
    real                   rlist_new, rlist_prev;
    size_t                 nstlist_ind = 0;
    t_state                state_tmp;
    gmx_bool               bBox, bDD, bCont;
    const char            *nstl_gpu = "\nFor optimal performance with a GPU nstlist (now %d) should be larger.\nThe optimum depends on your CPU and GPU resources.\nYou might want to try several nstlist values.\n";
    const char            *nve_err  = "Can not increase nstlist because an NVE ensemble is used";
    const char            *vbd_err  = "Can not increase nstlist because verlet-buffer-tolerance is not set or used";
    const char            *box_err  = "Can not increase nstlist because the box is too small";
    const char            *dd_err   = "Can not increase nstlist because of domain decomposition limitations";
    char                   buf[STRLEN];
    const float            oneThird = 1.0f / 3.0f;

    if (nstlist_cmdline <= 0)
    {
        if (ir->nstlist == 1)
        {
            /* The user probably set nstlist=1 for a reason,
             * don't mess with the settings.
             */
            return;
        }

        if (fp != NULL && bGPU && ir->nstlist < nstlist_try[0])
        {
            fprintf(fp, nstl_gpu, ir->nstlist);
        }
        nstlist_ind = 0;
        while (nstlist_ind < NNSTL && ir->nstlist >= nstlist_try[nstlist_ind])
        {
            nstlist_ind++;
        }
        if (nstlist_ind == NNSTL)
        {
            /* There are no larger nstlist value to try */
            return;
        }
    }

    if (EI_MD(ir->eI) && ir->etc == etcNO)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n", nve_err);
        }
        if (fp != NULL)
        {
            fprintf(fp, "%s\n", nve_err);
        }

        return;
    }

    if (ir->verletbuf_tol == 0 && bGPU)
    {
        gmx_fatal(FARGS, "You are using an old tpr file with a GPU, please generate a new tpr file with an up to date version of grompp");
    }

    if (ir->verletbuf_tol < 0)
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

    if (bGPU)
    {
        listfac_ok  = nbnxn_gpu_listfac_ok;
        listfac_max = nbnxn_gpu_listfac_max;
    }
    else
    {
        listfac_ok  = nbnxn_cpu_listfac_ok;
        listfac_max = nbnxn_cpu_listfac_max;
    }

    nstlist_orig = ir->nstlist;
    if (nstlist_cmdline > 0)
    {
        if (fp)
        {
            sprintf(buf, "Getting nstlist=%d from command line option",
                    nstlist_cmdline);
        }
        ir->nstlist = nstlist_cmdline;
    }

    verletbuf_get_list_setup(TRUE, bGPU, &ls);

    /* Allow rlist to make the list a given factor larger than the list
     * would be with the reference value for nstlist (10).
     */
    nstlist_prev = ir->nstlist;
    ir->nstlist  = nbnxnReferenceNstlist;
    calc_verlet_buffer_size(mtop, det(box), ir, -1, &ls, NULL,
                            &rlistWithReferenceNstlist);
    ir->nstlist  = nstlist_prev;

    /* Determine the pair list size increase due to zero interactions */
    rlist_inc = nbnxn_get_rlist_effective_inc(ls.cluster_size_j,
                                              mtop->natoms/det(box));
    rlist_ok  = (rlistWithReferenceNstlist + rlist_inc)*pow(listfac_ok, oneThird) - rlist_inc;
    rlist_max = (rlistWithReferenceNstlist + rlist_inc)*pow(listfac_max, oneThird) - rlist_inc;
    if (debug)
    {
        fprintf(debug, "nstlist tuning: rlist_inc %.3f rlist_ok %.3f rlist_max %.3f\n",
                rlist_inc, rlist_ok, rlist_max);
    }

    nstlist_prev = nstlist_orig;
    rlist_prev   = ir->rlist;
    do
    {
        if (nstlist_cmdline <= 0)
        {
            ir->nstlist = nstlist_try[nstlist_ind];
        }

        /* Set the pair-list buffer size in ir */
        calc_verlet_buffer_size(mtop, det(box), ir, -1, &ls, NULL, &rlist_new);

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

        if (debug)
        {
            fprintf(debug, "nstlist %d rlist %.3f bBox %d bDD %d\n",
                    ir->nstlist, rlist_new, bBox, bDD);
        }

        bCont = FALSE;

        if (nstlist_cmdline <= 0)
        {
            if (bBox && bDD && rlist_new <= rlist_max)
            {
                /* Increase nstlist */
                nstlist_prev = ir->nstlist;
                rlist_prev   = rlist_new;
                bCont        = (nstlist_ind+1 < NNSTL && rlist_new < rlist_ok);
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

        nstlist_ind++;
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

static void prepare_verlet_scheme(FILE                           *fplog,
                                  t_commrec                      *cr,
                                  t_inputrec                     *ir,
                                  int                             nstlist_cmdline,
                                  const gmx_mtop_t               *mtop,
                                  matrix                          box,
                                  gmx_bool                        bUseGPU)
{
    /* For NVE simulations, we will retain the initial list buffer */
    if (ir->verletbuf_tol > 0 && !(EI_MD(ir->eI) && ir->etc == etcNO))
    {
        /* Update the Verlet buffer size for the current run setup */
        verletbuf_list_setup_t ls;
        real                   rlist_new;

        /* Here we assume SIMD-enabled kernels are being used. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        verletbuf_get_list_setup(TRUE, bUseGPU, &ls);

        calc_verlet_buffer_size(mtop, det(box), ir, -1, &ls, NULL, &rlist_new);

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

    if (nstlist_cmdline > 0 && (!EI_DYNAMICS(ir->eI) || ir->verletbuf_tol <= 0))
    {
        gmx_fatal(FARGS, "Can not set nstlist without %s",
                  !EI_DYNAMICS(ir->eI) ? "dynamics" : "verlet-buffer-tolerance");
    }

    if (EI_DYNAMICS(ir->eI))
    {
        /* Set or try nstlist values */
        increase_nstlist(fplog, cr, ir, nstlist_cmdline, mtop, box, bUseGPU);
    }
}

/* Override the value in inputrec with value passed on the command line (if any) */
static void override_nsteps_cmdline(FILE            *fplog,
                                    gmx_int64_t      nsteps_cmdline,
                                    t_inputrec      *ir,
                                    const t_commrec *cr)
{
    assert(ir);
    assert(cr);

    /* override with anything else than the default -2 */
    if (nsteps_cmdline > -2)
    {
        char sbuf_steps[STEPSTRSIZE];
        char sbuf_msg[STRLEN];

        ir->nsteps = nsteps_cmdline;
        if (EI_DYNAMICS(ir->eI) && nsteps_cmdline != -1)
        {
            sprintf(sbuf_msg, "Overriding nsteps with value passed on the command line: %s steps, %.3g ps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps),
                    fabs(nsteps_cmdline*ir->delta_t));
        }
        else
        {
            sprintf(sbuf_msg, "Overriding nsteps with value passed on the command line: %s steps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps));
        }

        md_print_warn(cr, fplog, "%s\n", sbuf_msg);
    }
    else if (nsteps_cmdline < -2)
    {
        gmx_fatal(FARGS, "Invalid nsteps value passed on the command line: %d",
                  nsteps_cmdline);
    }
    /* Do nothing if nsteps_cmdline == -2 */
}

int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, t_commrec *cr, int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm,
             ivec ddxyz, int dd_node_order, real rdd, real rconstr,
             const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt, int nstlist_cmdline,
             gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
             int gmx_unused nmultisim, int repl_ex_nst, int repl_ex_nex,
             int repl_ex_seed, real pforce, real cpt_period, real max_hours,
             int imdport, unsigned long Flags)
{
    gmx_bool                  bForceUseGPU, bTryUseGPU, bRerunMD, bCantUseGPU;
    t_inputrec               *inputrec;
    t_state                  *state = NULL;
    matrix                    box;
    gmx_ddbox_t               ddbox = {0};
    int                       npme_major, npme_minor;
    t_nrnb                   *nrnb;
    gmx_mtop_t               *mtop          = NULL;
    t_mdatoms                *mdatoms       = NULL;
    t_forcerec               *fr            = NULL;
    t_fcdata                 *fcd           = NULL;
    real                      ewaldcoeff_q  = 0;
    real                      ewaldcoeff_lj = 0;
    struct gmx_pme_t        **pmedata       = NULL;
    gmx_vsite_t              *vsite         = NULL;
    gmx_constr_t              constr;
    int                       nChargePerturbed = -1, nTypePerturbed = 0, status;
    gmx_wallcycle_t           wcycle;
    gmx_bool                  bReadEkin;
    gmx_walltime_accounting_t walltime_accounting = NULL;
    int                       rc;
    gmx_int64_t               reset_counters;
    gmx_edsam_t               ed           = NULL;
    int                       nthreads_pme = 1;
    int                       nthreads_pp  = 1;
    gmx_membed_t              membed       = NULL;
    gmx_hw_info_t            *hwinfo       = NULL;
    /* The master rank decides early on bUseGPU and broadcasts this later */
    gmx_bool                  bUseGPU      = FALSE;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    snew(inputrec, 1);
    snew(mtop, 1);

    if (Flags & MD_APPENDFILES)
    {
        fplog = NULL;
    }

    bRerunMD     = (Flags & MD_RERUN);
    bForceUseGPU = (strncmp(nbpu_opt, "gpu", 3) == 0);
    bTryUseGPU   = (strncmp(nbpu_opt, "auto", 4) == 0) || bForceUseGPU;
    /* Rerun execution time is dominated by I/O and pair search, so
     * GPUs are not very useful, plus they do not support more than
     * one energy group. Don't select them when they can't be used,
     * unless the user requested it, then fatal_error is called later.
     *
     * TODO it would be nice to notify the user that if this check
     * causes GPUs not to be used that this is what is happening, and
     * why, but that will be easier to do after some future
     * cleanup. */
    bCantUseGPU = bRerunMD && (inputrec->opts.ngener > 1);
    bTryUseGPU  = bTryUseGPU && !(bCantUseGPU && !bForceUseGPU);

    /* Detect hardware, gather information. This is an operation that is
     * global for this process (MPI rank). */
    hwinfo = gmx_detect_hardware(fplog, cr, bTryUseGPU);

    gmx_print_detected_hardware(fplog, cr, hwinfo);

    if (fplog != NULL)
    {
        /* Print references after all software/hardware printing */
        please_cite(fplog, "Pall2015");
        please_cite(fplog, "Pronk2013");
        please_cite(fplog, "Hess2008b");
        please_cite(fplog, "Spoel2005a");
        please_cite(fplog, "Lindahl2001a");
        please_cite(fplog, "Berendsen95a");
    }

    snew(state, 1);
    if (SIMMASTER(cr))
    {
        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPR, nfile, fnm), inputrec, state, NULL, mtop);

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            /* Here the master rank decides if all ranks will use GPUs */
            bUseGPU = (hwinfo->gpu_info.n_dev_compatible > 0 ||
                       getenv("GMX_EMULATE_GPU") != NULL);

            /* TODO add GPU kernels for this and replace this check by:
             * (bUseGPU && (ir->vdwtype == evdwPME &&
             *               ir->ljpme_combination_rule == eljpmeLB))
             * update the message text and the content of nbnxn_acceleration_supported.
             */
            if (bUseGPU &&
                !nbnxn_acceleration_supported(fplog, cr, inputrec, bUseGPU))
            {
                /* Fallback message printed by nbnxn_acceleration_supported */
                if (bForceUseGPU)
                {
                    gmx_fatal(FARGS, "GPU acceleration requested, but not supported with the given input settings");
                }
                bUseGPU = FALSE;
            }

            prepare_verlet_scheme(fplog, cr,
                                  inputrec, nstlist_cmdline, mtop, state->box,
                                  bUseGPU);
        }
        else
        {
            if (nstlist_cmdline > 0)
            {
                gmx_fatal(FARGS, "Can not set nstlist with the group cut-off scheme");
            }

            if (hwinfo->gpu_info.n_dev_compatible > 0)
            {
                md_print_warn(cr, fplog,
                              "NOTE: GPU(s) found, but the current simulation can not use GPUs\n"
                              "      To use a GPU, set the mdp option: cutoff-scheme = Verlet\n");
            }

            if (bForceUseGPU)
            {
                gmx_fatal(FARGS, "GPU requested, but can't be used without cutoff-scheme=Verlet");
            }

#ifdef GMX_TARGET_BGQ
            md_print_warn(cr, fplog,
                          "NOTE: There is no SIMD implementation of the group scheme kernels on\n"
                          "      BlueGene/Q. You will observe better performance from using the\n"
                          "      Verlet cut-off scheme.\n");
#endif
        }

        if (inputrec->eI == eiSD2)
        {
            md_print_warn(cr, fplog, "The stochastic dynamics integrator %s is deprecated, since\n"
                          "it is slower than integrator %s and is slightly less accurate\n"
                          "with constraints. Use the %s integrator.",
                          ei_names[inputrec->eI], ei_names[eiSD1], ei_names[eiSD1]);
        }
    }

    /* Check and update the hardware options for internal consistency */
    check_and_update_hw_opt_1(hw_opt, cr);

    /* Early check for externally set process affinity. */
    gmx_check_thread_affinity_set(fplog, cr,
                                  hw_opt, hwinfo->nthreads_hw_avail, FALSE);

#ifdef GMX_THREAD_MPI
    if (SIMMASTER(cr))
    {
        if (cr->npmenodes > 0 && hw_opt->nthreads_tmpi <= 0)
        {
            gmx_fatal(FARGS, "You need to explicitly specify the number of MPI threads (-ntmpi) when using separate PME ranks");
        }

        /* Since the master knows the cut-off scheme, update hw_opt for this.
         * This is done later for normal MPI and also once more with tMPI
         * for all tMPI ranks.
         */
        check_and_update_hw_opt_2(hw_opt, inputrec->cutoff_scheme);

        /* NOW the threads will be started: */
        hw_opt->nthreads_tmpi = get_nthreads_mpi(hwinfo,
                                                 hw_opt,
                                                 inputrec, mtop,
                                                 cr, fplog);

        if (hw_opt->nthreads_tmpi > 1)
        {
            t_commrec *cr_old       = cr;
            /* now start the threads. */
            cr = mdrunner_start_threads(hw_opt, fplog, cr_old, nfile, fnm,
                                        oenv, bVerbose, bCompact, nstglobalcomm,
                                        ddxyz, dd_node_order, rdd, rconstr,
                                        dddlb_opt, dlb_scale, ddcsx, ddcsy, ddcsz,
                                        nbpu_opt, nstlist_cmdline,
                                        nsteps_cmdline, nstepout, resetstep, nmultisim,
                                        repl_ex_nst, repl_ex_nex, repl_ex_seed, pforce,
                                        cpt_period, max_hours,
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
        init_parallel(cr, inputrec, mtop);

        /* The master rank decided on the use of GPUs,
         * broadcast this information to all ranks.
         */
        gmx_bcast_sim(sizeof(bUseGPU), &bUseGPU, cr);
    }

    if (fplog != NULL)
    {
        pr_inputrec(fplog, 0, "Input Parameters", inputrec, FALSE);
        fprintf(fplog, "\n");
    }

    /* now make sure the state is initialized and propagated */
    set_state_entries(state, inputrec);

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
                  "but %s was not started through mpirun/mpiexec or only one rank was requested through mpirun/mpiexec"
#endif
#endif
                  , output_env_get_program_display_name(oenv)
                  );
    }

    if (bRerunMD &&
        (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
    }

    if (can_use_allvsall(inputrec, TRUE, cr, fplog) && DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "All-vs-all loops do not work with domain decomposition, use a single MPI rank");
    }

    if (!(EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype)))
    {
        if (cr->npmenodes > 0)
        {
            gmx_fatal_collective(FARGS, cr, NULL,
                                 "PME-only ranks are requested, but the system does not use PME for electrostatics or LJ");
        }

        cr->npmenodes = 0;
    }

    if (bUseGPU && cr->npmenodes < 0)
    {
        /* With GPUs we don't automatically use PME-only ranks. PME ranks can
         * improve performance with many threads per GPU, since our OpenMP
         * scaling is bad, but it's difficult to automate the setup.
         */
        cr->npmenodes = 0;
    }

#ifdef GMX_FAHCORE
    if (MASTER(cr))
    {
        fcRegisterSteps(inputrec->nsteps, inputrec->init_step);
    }
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
    init_disres(fplog, mtop, inputrec, cr, fcd, state, repl_ex_nst > 0);

    init_orires(fplog, mtop, state->x, inputrec, cr, &(fcd->orires),
                state);

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
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        deform_init_init_step_tpx = inputrec->init_step;
        copy_mat(box, deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    if (opt2bSet("-cpi", nfile, fnm))
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        if (gmx_fexist_master(opt2fn_master("-cpi", nfile, fnm, cr), cr) )
        {
            load_checkpoint(opt2fn_master("-cpi", nfile, fnm, cr), &fplog,
                            cr, ddxyz,
                            inputrec, state, &bReadEkin,
                            (Flags & MD_APPENDFILES),
                            (Flags & MD_APPENDFILESSET));

            if (bReadEkin)
            {
                Flags |= MD_READ_EKIN;
            }
        }
    }

    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr,
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

    if (PAR(cr) && !(EI_TPI(inputrec->eI) ||
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
    if (MULTISIM(cr))
    {
        md_print_info(cr, fplog,
                      "This is simulation %d out of %d running as a composite GROMACS\n"
                      "multi-simulation job. Setup for this simulation:\n\n",
                      cr->ms->sim, cr->ms->nsim);
    }
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

    /* Check and update hw_opt for the cut-off scheme */
    check_and_update_hw_opt_2(hw_opt, inputrec->cutoff_scheme);

    /* Check and update hw_opt for the number of MPI ranks */
    check_and_update_hw_opt_3(hw_opt);

    gmx_omp_nthreads_init(fplog, cr,
                          hwinfo->nthreads_hw_avail,
                          hw_opt->nthreads_omp,
                          hw_opt->nthreads_omp_pme,
                          (cr->duty & DUTY_PP) == 0,
                          inputrec->cutoff_scheme == ecutsVERLET);

#ifndef NDEBUG
    if (integrator[inputrec->eI].func != do_tpi &&
        inputrec->cutoff_scheme == ecutsVERLET)
    {
        gmx_feenableexcept();
    }
#endif

    if (bUseGPU)
    {
        /* Select GPU id's to use */
        gmx_select_gpu_ids(fplog, cr, &hwinfo->gpu_info, bForceUseGPU,
                           &hw_opt->gpu_opt);
    }
    else
    {
        /* Ignore (potentially) manually selected GPUs */
        hw_opt->gpu_opt.n_dev_use = 0;
    }

    /* check consistency across ranks of things like SIMD
     * support and number of GPUs selected */
    gmx_check_hw_runconf_consistency(fplog, hwinfo, cr, hw_opt, bUseGPU);

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo, hw_opt, Flags & MD_NTOMPSET,
                                       cr, fplog);

    if (DOMAINDECOMP(cr))
    {
        /* When we share GPUs over ranks, we need to know this for the DLB */
        dd_setup_dlb_resource_sharing(cr, hwinfo, hw_opt);
    }

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
        bcast_state(cr, state);

        /* Initiate forcerecord */
        fr          = mk_forcerec();
        fr->hwinfo  = hwinfo;
        fr->gpu_opt = &hw_opt->gpu_opt;
        init_forcerec(fplog, oenv, fr, fcd, inputrec, mtop, cr, box,
                      opt2fn("-table", nfile, fnm),
                      opt2fn("-tabletf", nfile, fnm),
                      opt2fn("-tablep", nfile, fnm),
                      opt2fn("-tableb", nfile, fnm),
                      nbpu_opt,
                      FALSE,
                      pforce);

        /* version for PCA_NOT_READ_NODE (see md.c) */
        /*init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
           "nofile","nofile","nofile","nofile",FALSE,pforce);
         */

        /* Initialize QM-MM */
        if (fr->bQMMM)
        {
            init_QMMMrec(cr, mtop, inputrec, fr);
        }

        /* Initialize the mdatoms structure.
         * mdatoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdatoms = init_mdatoms(fplog, mtop, inputrec->efep != efepNO);

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
                construct_vsites_mtop(vsite, mtop, state->x);
            }
        }

        if (EEL_PME(fr->eeltype) || EVDW_PME(fr->vdwtype))
        {
            ewaldcoeff_q  = fr->ewaldcoeff_q;
            ewaldcoeff_lj = fr->ewaldcoeff_lj;
            pmedata       = &fr->pmedata;
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

        ewaldcoeff_q  = calc_ewaldcoeff_q(inputrec->rcoulomb, inputrec->ewald_rtol);
        ewaldcoeff_lj = calc_ewaldcoeff_lj(inputrec->rvdw, inputrec->ewald_rtol_lj);
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
        gmx_set_thread_affinity(fplog, cr, hw_opt, hwinfo);
    }

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype))
    {
        if (mdatoms)
        {
            nChargePerturbed = mdatoms->nChargePerturbed;
            if (EVDW_PME(inputrec->vdwtype))
            {
                nTypePerturbed   = mdatoms->nTypePerturbed;
            }
        }
        if (cr->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed(FEP on Q) and nTypePerturbed(FEP on LJ)*/
            gmx_bcast_sim(sizeof(nChargePerturbed), &nChargePerturbed, cr);
            gmx_bcast_sim(sizeof(nTypePerturbed), &nTypePerturbed, cr);
        }

        if (cr->duty & DUTY_PME)
        {
            status = gmx_pme_init(pmedata, cr, npme_major, npme_minor, inputrec,
                                  mtop ? mtop->natoms : 0, nChargePerturbed, nTypePerturbed,
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
        /* Assumes uniform use of the number of OpenMP threads */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntDefault));

        if (inputrec->bPull)
        {
            /* Initialize pull code */
            inputrec->pull_work =
                init_pull(fplog, inputrec->pull, inputrec, nfile, fnm,
                          mtop, cr, oenv, inputrec->fepvals->init_lambda,
                          EI_DYNAMICS(inputrec->eI) && MASTER(cr), Flags);
        }

        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            init_rot(fplog, inputrec, nfile, fnm, cr, state->x, box, mtop, oenv,
                     bVerbose, Flags);
        }

        if (inputrec->eSwapCoords != eswapNO)
        {
            /* Initialize ion swapping code */
            init_swapcoords(fplog, bVerbose, inputrec, opt2fn_master("-swap", nfile, fnm, cr),
                            mtop, state->x, state->box, &state->swapstate, cr, oenv, Flags);
        }

        constr = init_constraints(fplog, mtop, inputrec, ed, state, cr);

        if (DOMAINDECOMP(cr))
        {
            GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
            dd_init_bondeds(fplog, cr->dd, mtop, vsite, inputrec,
                            Flags & MD_DDBONDCHECK, fr->cginfo_mb);

            set_dd_parameters(fplog, cr->dd, dlb_scale, inputrec, &ddbox);

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
                                      imdport,
                                      Flags,
                                      walltime_accounting);

        if (inputrec->bPull)
        {
            finish_pull(inputrec->pull_work);
        }

        if (inputrec->bRot)
        {
            finish_rot(inputrec->rot);
        }

    }
    else
    {
        GMX_RELEASE_ASSERT(pmedata, "pmedata was NULL while cr->duty was not DUTY_PP");
        /* do PME only */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntPME));
        gmx_pmeonly(*pmedata, cr, nrnb, wcycle, walltime_accounting, ewaldcoeff_q, ewaldcoeff_lj, inputrec);
    }

    wallcycle_stop(wcycle, ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog, cr,
               inputrec, nrnb, wcycle, walltime_accounting,
               fr ? fr->nbv : NULL,
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));


    /* Free GPU memory and context */
    free_gpu_resources(fr, cr, &hwinfo->gpu_info, fr ? fr->gpu_opt : NULL);

    if (opt2bSet("-membed", nfile, fnm))
    {
        sfree(membed);
    }

    gmx_hardware_info_free(hwinfo);

    /* Does what it says */
    print_date_and_time(fplog, cr->nodeid, "Finished mdrun", gmx_gettime());
    walltime_accounting_destroy(walltime_accounting);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && (Flags & MD_APPENDFILES))
    {
        gmx_log_close(fplog);
    }

    rc = (int)gmx_get_stop_condition();

    done_ed(&ed);

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
