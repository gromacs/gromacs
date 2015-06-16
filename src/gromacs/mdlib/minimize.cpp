/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <math.h>
#include <string.h>
#include <time.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/fileio/trajectory_writing.h"
#include "gromacs/imd/imd.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/md_support.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    t_state  s;
    rvec    *f;
    real     epot;
    real     fnorm;
    real     fmax;
    int      a_fmax;
} em_state_t;

static em_state_t *init_em_state()
{
    em_state_t *ems;

    snew(ems, 1);

    /* does this need to be here?  Should the array be declared differently (staticaly)in the state definition? */
    snew(ems->s.lambda, efptNR);

    return ems;
}

static void print_em_start(FILE                     *fplog,
                           t_commrec                *cr,
                           gmx_walltime_accounting_t walltime_accounting,
                           gmx_wallcycle_t           wcycle,
                           const char               *name)
{
    walltime_accounting_start(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, name);
}
static void em_time_end(gmx_walltime_accounting_t walltime_accounting,
                        gmx_wallcycle_t           wcycle)
{
    wallcycle_stop(wcycle, ewcRUN);

    walltime_accounting_end(walltime_accounting);
}

static void sp_header(FILE *out, const char *minimizer, real ftol, int nsteps)
{
    fprintf(out, "\n");
    fprintf(out, "%s:\n", minimizer);
    fprintf(out, "   Tolerance (Fmax)   = %12.5e\n", ftol);
    fprintf(out, "   Number of steps    = %12d\n", nsteps);
}

static void warn_step(FILE *fp, real ftol, gmx_bool bLastStep, gmx_bool bConstrain)
{
    char buffer[2048];
    if (bLastStep)
    {
        sprintf(buffer,
                "\nEnergy minimization reached the maximum number "
                "of steps before the forces reached the requested "
                "precision Fmax < %g.\n", ftol);
    }
    else
    {
        sprintf(buffer,
                "\nEnergy minimization has stopped, but the forces have "
                "not converged to the requested precision Fmax < %g (which "
                "may not be possible for your system). It stopped "
                "because the algorithm tried to make a new step whose size "
                "was too small, or there was no change in the energy since "
                "last step. Either way, we regard the minimization as "
                "converged to within the available machine precision, "
                "given your starting configuration and EM parameters.\n%s%s",
                ftol,
                sizeof(real) < sizeof(double) ?
                "\nDouble precision normally gives you higher accuracy, but "
                "this is often not needed for preparing to run molecular "
                "dynamics.\n" :
                "",
                bConstrain ?
                "You might need to increase your constraint accuracy, or turn\n"
                "off constraints altogether (set constraints = none in mdp file)\n" :
                "");
    }
    fputs(wrap_lines(buffer, 78, 0, FALSE), fp);
}



static void print_converged(FILE *fp, const char *alg, real ftol,
                            gmx_int64_t count, gmx_bool bDone, gmx_int64_t nsteps,
                            real epot, real fmax, int nfmax, real fnorm)
{
    char buf[STEPSTRSIZE];

    if (bDone)
    {
        fprintf(fp, "\n%s converged to Fmax < %g in %s steps\n",
                alg, ftol, gmx_step_str(count, buf));
    }
    else if (count < nsteps)
    {
        fprintf(fp, "\n%s converged to machine precision in %s steps,\n"
                "but did not reach the requested Fmax < %g.\n",
                alg, gmx_step_str(count, buf), ftol);
    }
    else
    {
        fprintf(fp, "\n%s did not converge to Fmax < %g in %s steps.\n",
                alg, ftol, gmx_step_str(count, buf));
    }

#ifdef GMX_DOUBLE
    fprintf(fp, "Potential Energy  = %21.14e\n", epot);
    fprintf(fp, "Maximum force     = %21.14e on atom %d\n", fmax, nfmax+1);
    fprintf(fp, "Norm of force     = %21.14e\n", fnorm);
#else
    fprintf(fp, "Potential Energy  = %14.7e\n", epot);
    fprintf(fp, "Maximum force     = %14.7e on atom %d\n", fmax, nfmax+1);
    fprintf(fp, "Norm of force     = %14.7e\n", fnorm);
#endif
}

static void get_f_norm_max(t_commrec *cr,
                           t_grpopts *opts, t_mdatoms *mdatoms, rvec *f,
                           real *fnorm, real *fmax, int *a_fmax)
{
    double fnorm2, *sum;
    real   fmax2, fam;
    int    la_max, a_max, start, end, i, m, gf;

    /* This routine finds the largest force and returns it.
     * On parallel machines the global max is taken.
     */
    fnorm2 = 0;
    fmax2  = 0;
    la_max = -1;
    start  = 0;
    end    = mdatoms->homenr;
    if (mdatoms->cFREEZE)
    {
        for (i = start; i < end; i++)
        {
            gf  = mdatoms->cFREEZE[i];
            fam = 0;
            for (m = 0; m < DIM; m++)
            {
                if (!opts->nFreeze[gf][m])
                {
                    fam += sqr(f[i][m]);
                }
            }
            fnorm2 += fam;
            if (fam > fmax2)
            {
                fmax2  = fam;
                la_max = i;
            }
        }
    }
    else
    {
        for (i = start; i < end; i++)
        {
            fam     = norm2(f[i]);
            fnorm2 += fam;
            if (fam > fmax2)
            {
                fmax2  = fam;
                la_max = i;
            }
        }
    }

    if (la_max >= 0 && DOMAINDECOMP(cr))
    {
        a_max = cr->dd->gatindex[la_max];
    }
    else
    {
        a_max = la_max;
    }
    if (PAR(cr))
    {
        snew(sum, 2*cr->nnodes+1);
        sum[2*cr->nodeid]   = fmax2;
        sum[2*cr->nodeid+1] = a_max;
        sum[2*cr->nnodes]   = fnorm2;
        gmx_sumd(2*cr->nnodes+1, sum, cr);
        fnorm2 = sum[2*cr->nnodes];
        /* Determine the global maximum */
        for (i = 0; i < cr->nnodes; i++)
        {
            if (sum[2*i] > fmax2)
            {
                fmax2 = sum[2*i];
                a_max = (int)(sum[2*i+1] + 0.5);
            }
        }
        sfree(sum);
    }

    if (fnorm)
    {
        *fnorm = sqrt(fnorm2);
    }
    if (fmax)
    {
        *fmax  = sqrt(fmax2);
    }
    if (a_fmax)
    {
        *a_fmax = a_max;
    }
}

static void get_state_f_norm_max(t_commrec *cr,
                                 t_grpopts *opts, t_mdatoms *mdatoms,
                                 em_state_t *ems)
{
    get_f_norm_max(cr, opts, mdatoms, ems->f, &ems->fnorm, &ems->fmax, &ems->a_fmax);
}

void init_em(FILE *fplog, const char *title,
             t_commrec *cr, t_inputrec *ir,
             t_state *state_global, gmx_mtop_t *top_global,
             em_state_t *ems, gmx_localtop_t **top,
             rvec **f, rvec **f_global,
             t_nrnb *nrnb, rvec mu_tot,
             t_forcerec *fr, gmx_enerdata_t **enerd,
             t_graph **graph, t_mdatoms *mdatoms, gmx_global_stat_t *gstat,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             int imdport, unsigned long gmx_unused Flags,
             gmx_wallcycle_t wcycle)
{
    int  i;
    real dvdl_constr;

    if (fplog)
    {
        fprintf(fplog, "Initiating %s\n", title);
    }

    state_global->ngtc = 0;

    /* Initialize lambda variables */
    initialize_lambdas(fplog, ir, &(state_global->fep_state), state_global->lambda, NULL);

    init_nrnb(nrnb);

    /* Interactive molecular dynamics */
    init_IMD(ir, cr, top_global, fplog, 1, state_global->x,
             nfile, fnm, NULL, imdport, Flags);

    if (DOMAINDECOMP(cr))
    {
        *top = dd_init_local_top(top_global);

        dd_init_local_state(cr->dd, state_global, &ems->s);

        *f = NULL;

        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            &ems->s, &ems->f, mdatoms, *top,
                            fr, vsite, NULL, constr,
                            nrnb, NULL, FALSE);
        dd_store_state(cr->dd, &ems->s);

        if (ir->nstfout)
        {
            snew(*f_global, top_global->natoms);
        }
        else
        {
            *f_global = NULL;
        }
        *graph = NULL;
    }
    else
    {
        snew(*f, top_global->natoms);

        /* Just copy the state */
        ems->s = *state_global;
        snew(ems->s.x, ems->s.nalloc);
        snew(ems->f, ems->s.nalloc);
        for (i = 0; i < state_global->natoms; i++)
        {
            copy_rvec(state_global->x[i], ems->s.x[i]);
        }
        copy_mat(state_global->box, ems->s.box);

        *top      = gmx_mtop_generate_local_top(top_global, ir);
        *f_global = *f;

        forcerec_set_excl_load(fr, *top);

        setup_bonded_threading(fr, &(*top)->idef);

        if (ir->ePBC != epbcNONE && !fr->bMolPBC)
        {
            *graph = mk_graph(fplog, &((*top)->idef), 0, top_global->natoms, FALSE, FALSE);
        }
        else
        {
            *graph = NULL;
        }

        atoms2md(top_global, ir, 0, NULL, top_global->natoms, mdatoms);
        update_mdatoms(mdatoms, state_global->lambda[efptFEP]);

        if (vsite)
        {
            set_vsite_top(vsite, *top, mdatoms, cr);
        }
    }

    if (constr)
    {
        if (ir->eConstrAlg == econtSHAKE &&
            gmx_mtop_ftype_count(top_global, F_CONSTR) > 0)
        {
            gmx_fatal(FARGS, "Can not do energy minimization with %s, use %s\n",
                      econstr_names[econtSHAKE], econstr_names[econtLINCS]);
        }

        if (!DOMAINDECOMP(cr))
        {
            set_constraints(constr, *top, ir, mdatoms, cr);
        }

        if (!ir->bContinuation)
        {
            /* Constrain the starting coordinates */
            dvdl_constr = 0;
            constrain(PAR(cr) ? NULL : fplog, TRUE, TRUE, constr, &(*top)->idef,
                      ir, cr, -1, 0, 1.0, mdatoms,
                      ems->s.x, ems->s.x, NULL, fr->bMolPBC, ems->s.box,
                      ems->s.lambda[efptFEP], &dvdl_constr,
                      NULL, NULL, nrnb, econqCoord);
        }
    }

    if (PAR(cr))
    {
        *gstat = global_stat_init(ir);
    }
    else
    {
        *gstat = NULL;
    }

    *outf = init_mdoutf(fplog, nfile, fnm, 0, cr, ir, top_global, NULL, wcycle);

    snew(*enerd, 1);
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  *enerd);

    if (mdebin != NULL)
    {
        /* Init bin for energy stuff */
        *mdebin = init_mdebin(mdoutf_get_fp_ene(*outf), top_global, ir, NULL);
    }

    clear_rvec(mu_tot);
    calc_shifts(ems->s.box, fr->shift_vec);
}

static void finish_em(t_commrec *cr, gmx_mdoutf_t outf,
                      gmx_walltime_accounting_t walltime_accounting,
                      gmx_wallcycle_t wcycle)
{
    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    done_mdoutf(outf);

    em_time_end(walltime_accounting, wcycle);
}

static void swap_em_state(em_state_t *ems1, em_state_t *ems2)
{
    em_state_t tmp;

    tmp   = *ems1;
    *ems1 = *ems2;
    *ems2 = tmp;
}

static void copy_em_coords(em_state_t *ems, t_state *state)
{
    int i;

    for (i = 0; (i < state->natoms); i++)
    {
        copy_rvec(ems->s.x[i], state->x[i]);
    }
}

static void write_em_traj(FILE *fplog, t_commrec *cr,
                          gmx_mdoutf_t outf,
                          gmx_bool bX, gmx_bool bF, const char *confout,
                          gmx_mtop_t *top_global,
                          t_inputrec *ir, gmx_int64_t step,
                          em_state_t *state,
                          t_state *state_global, rvec *f_global)
{
    int      mdof_flags;
    gmx_bool bIMDout = FALSE;


    /* Shall we do IMD output? */
    if (ir->bIMD)
    {
        bIMDout = do_per_step(step, IMD_get_step(ir->imd->setup));
    }

    if ((bX || bF || bIMDout || confout != NULL) && !DOMAINDECOMP(cr))
    {
        copy_em_coords(state, state_global);
        f_global = state->f;
    }

    mdof_flags = 0;
    if (bX)
    {
        mdof_flags |= MDOF_X;
    }
    if (bF)
    {
        mdof_flags |= MDOF_F;
    }

    /* If we want IMD output, set appropriate MDOF flag */
    if (ir->bIMD)
    {
        mdof_flags |= MDOF_IMD;
    }

    mdoutf_write_to_trajectory_files(fplog, cr, outf, mdof_flags,
                                     top_global, step, (double)step,
                                     &state->s, state_global, state->f, f_global);

    if (confout != NULL && MASTER(cr))
    {
        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols && DOMAINDECOMP(cr))
        {
            /* Make molecules whole only for confout writing */
            do_pbc_mtop(fplog, ir->ePBC, state_global->box, top_global,
                        state_global->x);
        }

        write_sto_conf_mtop(confout,
                            *top_global->name, top_global,
                            state_global->x, NULL, ir->ePBC, state_global->box);
    }
}

static void do_em_step(t_commrec *cr, t_inputrec *ir, t_mdatoms *md,
                       gmx_bool bMolPBC,
                       em_state_t *ems1, real a, rvec *f, em_state_t *ems2,
                       gmx_constr_t constr, gmx_localtop_t *top,
                       t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                       gmx_int64_t count)

{
    t_state *s1, *s2;
    int      i;
    int      start, end;
    rvec    *x1, *x2;
    real     dvdl_constr;
    int      nthreads gmx_unused;

    s1 = &ems1->s;
    s2 = &ems2->s;

    if (DOMAINDECOMP(cr) && s1->ddp_count != cr->dd->ddp_count)
    {
        gmx_incons("state mismatch in do_em_step");
    }

    s2->flags = s1->flags;

    if (s2->nalloc != s1->nalloc)
    {
        s2->nalloc = s1->nalloc;
        srenew(s2->x, s1->nalloc);
        srenew(ems2->f,  s1->nalloc);
        if (s2->flags & (1<<estCGP))
        {
            srenew(s2->cg_p,  s1->nalloc);
        }
    }

    s2->natoms = s1->natoms;
    copy_mat(s1->box, s2->box);
    /* Copy free energy state */
    for (i = 0; i < efptNR; i++)
    {
        s2->lambda[i] = s1->lambda[i];
    }
    copy_mat(s1->box, s2->box);

    start = 0;
    end   = md->homenr;

    x1 = s1->x;
    x2 = s2->x;

    // cppcheck-suppress unreadVariable
    nthreads = gmx_omp_nthreads_get(emntUpdate);
#pragma omp parallel num_threads(nthreads)
    {
        int gf, i, m;

        gf = 0;
#pragma omp for schedule(static) nowait
        for (i = start; i < end; i++)
        {
            if (md->cFREEZE)
            {
                gf = md->cFREEZE[i];
            }
            for (m = 0; m < DIM; m++)
            {
                if (ir->opts.nFreeze[gf][m])
                {
                    x2[i][m] = x1[i][m];
                }
                else
                {
                    x2[i][m] = x1[i][m] + a*f[i][m];
                }
            }
        }

        if (s2->flags & (1<<estCGP))
        {
            /* Copy the CG p vector */
            x1 = s1->cg_p;
            x2 = s2->cg_p;
#pragma omp for schedule(static) nowait
            for (i = start; i < end; i++)
            {
                copy_rvec(x1[i], x2[i]);
            }
        }

        if (DOMAINDECOMP(cr))
        {
            s2->ddp_count = s1->ddp_count;
            if (s2->cg_gl_nalloc < s1->cg_gl_nalloc)
            {
#pragma omp barrier
                s2->cg_gl_nalloc = s1->cg_gl_nalloc;
                srenew(s2->cg_gl, s2->cg_gl_nalloc);
#pragma omp barrier
            }
            s2->ncg_gl = s1->ncg_gl;
#pragma omp for schedule(static) nowait
            for (i = 0; i < s2->ncg_gl; i++)
            {
                s2->cg_gl[i] = s1->cg_gl[i];
            }
            s2->ddp_count_cg_gl = s1->ddp_count_cg_gl;
        }
    }

    if (constr)
    {
        wallcycle_start(wcycle, ewcCONSTR);
        dvdl_constr = 0;
        constrain(NULL, TRUE, TRUE, constr, &top->idef,
                  ir, cr, count, 0, 1.0, md,
                  s1->x, s2->x, NULL, bMolPBC, s2->box,
                  s2->lambda[efptBONDED], &dvdl_constr,
                  NULL, NULL, nrnb, econqCoord);
        wallcycle_stop(wcycle, ewcCONSTR);
    }
}

static void em_dd_partition_system(FILE *fplog, int step, t_commrec *cr,
                                   gmx_mtop_t *top_global, t_inputrec *ir,
                                   em_state_t *ems, gmx_localtop_t *top,
                                   t_mdatoms *mdatoms, t_forcerec *fr,
                                   gmx_vsite_t *vsite, gmx_constr_t constr,
                                   t_nrnb *nrnb, gmx_wallcycle_t wcycle)
{
    /* Repartition the domain decomposition */
    dd_partition_system(fplog, step, cr, FALSE, 1,
                        NULL, top_global, ir,
                        &ems->s, &ems->f,
                        mdatoms, top, fr, vsite, NULL, constr,
                        nrnb, wcycle, FALSE);
    dd_store_state(cr->dd, &ems->s);
}

static void evaluate_energy(FILE *fplog, t_commrec *cr,
                            gmx_mtop_t *top_global,
                            em_state_t *ems, gmx_localtop_t *top,
                            t_inputrec *inputrec,
                            t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                            gmx_global_stat_t gstat,
                            gmx_vsite_t *vsite, gmx_constr_t constr,
                            t_fcdata *fcd,
                            t_graph *graph, t_mdatoms *mdatoms,
                            t_forcerec *fr, rvec mu_tot,
                            gmx_enerdata_t *enerd, tensor vir, tensor pres,
                            gmx_int64_t count, gmx_bool bFirst)
{
    real     t;
    gmx_bool bNS;
    tensor   force_vir, shake_vir, ekin;
    real     dvdl_constr, prescorr, enercorr, dvdlcorr;
    real     terminate = 0;

    /* Set the time to the initial time, the time does not change during EM */
    t = inputrec->init_t;

    if (bFirst ||
        (DOMAINDECOMP(cr) && ems->s.ddp_count < cr->dd->ddp_count))
    {
        /* This is the first state or an old state used before the last ns */
        bNS = TRUE;
    }
    else
    {
        bNS = FALSE;
        if (inputrec->nstlist > 0)
        {
            bNS = TRUE;
        }
    }

    if (vsite)
    {
        construct_vsites(vsite, ems->s.x, 1, NULL,
                         top->idef.iparams, top->idef.il,
                         fr->ePBC, fr->bMolPBC, cr, ems->s.box);
    }

    if (DOMAINDECOMP(cr) && bNS)
    {
        /* Repartition the domain decomposition */
        em_dd_partition_system(fplog, count, cr, top_global, inputrec,
                               ems, top, mdatoms, fr, vsite, constr,
                               nrnb, wcycle);
    }

    /* Calc force & energy on new trial position  */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in congrad.c
     */
    do_force(fplog, cr, inputrec,
             count, nrnb, wcycle, top, &top_global->groups,
             ems->s.box, ems->s.x, &ems->s.hist,
             ems->f, force_vir, mdatoms, enerd, fcd,
             ems->s.lambda, graph, fr, vsite, mu_tot, t, NULL, NULL, TRUE,
             GMX_FORCE_STATECHANGED | GMX_FORCE_ALLFORCES |
             GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY |
             (bNS ? GMX_FORCE_NS | GMX_FORCE_DO_LR : 0));

    /* Clear the unused shake virial and pressure */
    clear_mat(shake_vir);
    clear_mat(pres);

    /* Communicate stuff when parallel */
    if (PAR(cr) && inputrec->eI != eiNM)
    {
        wallcycle_start(wcycle, ewcMoveE);

        global_stat(fplog, gstat, cr, enerd, force_vir, shake_vir, mu_tot,
                    inputrec, NULL, NULL, NULL, 1, &terminate,
                    top_global, &ems->s, FALSE,
                    CGLO_ENERGY |
                    CGLO_PRESSURE |
                    CGLO_CONSTRAINT);

        wallcycle_stop(wcycle, ewcMoveE);
    }

    /* Calculate long range corrections to pressure and energy */
    calc_dispcorr(inputrec, fr, top_global->natoms, ems->s.box, ems->s.lambda[efptVDW],
                  pres, force_vir, &prescorr, &enercorr, &dvdlcorr);
    enerd->term[F_DISPCORR] = enercorr;
    enerd->term[F_EPOT]    += enercorr;
    enerd->term[F_PRES]    += prescorr;
    enerd->term[F_DVDL]    += dvdlcorr;

    ems->epot = enerd->term[F_EPOT];

    if (constr)
    {
        /* Project out the constraint components of the force */
        wallcycle_start(wcycle, ewcCONSTR);
        dvdl_constr = 0;
        constrain(NULL, FALSE, FALSE, constr, &top->idef,
                  inputrec, cr, count, 0, 1.0, mdatoms,
                  ems->s.x, ems->f, ems->f, fr->bMolPBC, ems->s.box,
                  ems->s.lambda[efptBONDED], &dvdl_constr,
                  NULL, &shake_vir, nrnb, econqForceDispl);
        enerd->term[F_DVDL_CONSTR] += dvdl_constr;
        m_add(force_vir, shake_vir, vir);
        wallcycle_stop(wcycle, ewcCONSTR);
    }
    else
    {
        copy_mat(force_vir, vir);
    }

    clear_mat(ekin);
    enerd->term[F_PRES] =
        calc_pres(fr->ePBC, inputrec->nwall, ems->s.box, ekin, vir, pres);

    sum_dhdl(enerd, ems->s.lambda, inputrec->fepvals);

    if (EI_ENERGY_MINIMIZATION(inputrec->eI))
    {
        get_state_f_norm_max(cr, &(inputrec->opts), mdatoms, ems);
    }
}

static double reorder_partsum(t_commrec *cr, t_grpopts *opts, t_mdatoms *mdatoms,
                              gmx_mtop_t *mtop,
                              em_state_t *s_min, em_state_t *s_b)
{
    rvec          *fm, *fb, *fmg;
    t_block       *cgs_gl;
    int            ncg, *cg_gl, *index, c, cg, i, a0, a1, a, gf, m;
    double         partsum;
    unsigned char *grpnrFREEZE;

    if (debug)
    {
        fprintf(debug, "Doing reorder_partsum\n");
    }

    fm = s_min->f;
    fb = s_b->f;

    cgs_gl = dd_charge_groups_global(cr->dd);
    index  = cgs_gl->index;

    /* Collect fm in a global vector fmg.
     * This conflicts with the spirit of domain decomposition,
     * but to fully optimize this a much more complicated algorithm is required.
     */
    snew(fmg, mtop->natoms);

    ncg   = s_min->s.ncg_gl;
    cg_gl = s_min->s.cg_gl;
    i     = 0;
    for (c = 0; c < ncg; c++)
    {
        cg = cg_gl[c];
        a0 = index[cg];
        a1 = index[cg+1];
        for (a = a0; a < a1; a++)
        {
            copy_rvec(fm[i], fmg[a]);
            i++;
        }
    }
    gmx_sum(mtop->natoms*3, fmg[0], cr);

    /* Now we will determine the part of the sum for the cgs in state s_b */
    ncg         = s_b->s.ncg_gl;
    cg_gl       = s_b->s.cg_gl;
    partsum     = 0;
    i           = 0;
    gf          = 0;
    grpnrFREEZE = mtop->groups.grpnr[egcFREEZE];
    for (c = 0; c < ncg; c++)
    {
        cg = cg_gl[c];
        a0 = index[cg];
        a1 = index[cg+1];
        for (a = a0; a < a1; a++)
        {
            if (mdatoms->cFREEZE && grpnrFREEZE)
            {
                gf = grpnrFREEZE[i];
            }
            for (m = 0; m < DIM; m++)
            {
                if (!opts->nFreeze[gf][m])
                {
                    partsum += (fb[i][m] - fmg[a][m])*fb[i][m];
                }
            }
            i++;
        }
    }

    sfree(fmg);

    return partsum;
}

static real pr_beta(t_commrec *cr, t_grpopts *opts, t_mdatoms *mdatoms,
                    gmx_mtop_t *mtop,
                    em_state_t *s_min, em_state_t *s_b)
{
    rvec  *fm, *fb;
    double sum;
    int    gf, i, m;

    /* This is just the classical Polak-Ribiere calculation of beta;
     * it looks a bit complicated since we take freeze groups into account,
     * and might have to sum it in parallel runs.
     */

    if (!DOMAINDECOMP(cr) ||
        (s_min->s.ddp_count == cr->dd->ddp_count &&
         s_b->s.ddp_count   == cr->dd->ddp_count))
    {
        fm  = s_min->f;
        fb  = s_b->f;
        sum = 0;
        gf  = 0;
        /* This part of code can be incorrect with DD,
         * since the atom ordering in s_b and s_min might differ.
         */
        for (i = 0; i < mdatoms->homenr; i++)
        {
            if (mdatoms->cFREEZE)
            {
                gf = mdatoms->cFREEZE[i];
            }
            for (m = 0; m < DIM; m++)
            {
                if (!opts->nFreeze[gf][m])
                {
                    sum += (fb[i][m] - fm[i][m])*fb[i][m];
                }
            }
        }
    }
    else
    {
        /* We need to reorder cgs while summing */
        sum = reorder_partsum(cr, opts, mdatoms, mtop, s_min, s_b);
    }
    if (PAR(cr))
    {
        gmx_sumd(1, &sum, cr);
    }

    return sum/sqr(s_min->fnorm);
}

double do_cg(FILE *fplog, t_commrec *cr,
             int nfile, const t_filenm fnm[],
             const output_env_t gmx_unused oenv, gmx_bool bVerbose, gmx_bool gmx_unused bCompact,
             int gmx_unused nstglobalcomm,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int gmx_unused stepout,
             t_inputrec *inputrec,
             gmx_mtop_t *top_global, t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb, gmx_wallcycle_t wcycle,
             gmx_edsam_t gmx_unused ed,
             t_forcerec *fr,
             int gmx_unused repl_ex_nst, int gmx_unused repl_ex_nex, int gmx_unused repl_ex_seed,
             gmx_membed_t gmx_unused membed,
             real gmx_unused cpt_period, real gmx_unused max_hours,
             int imdport,
             unsigned long gmx_unused Flags,
             gmx_walltime_accounting_t walltime_accounting)
{
    const char       *CG = "Polak-Ribiere Conjugate Gradients";

    em_state_t       *s_min, *s_a, *s_b, *s_c;
    gmx_localtop_t   *top;
    gmx_enerdata_t   *enerd;
    rvec             *f;
    gmx_global_stat_t gstat;
    t_graph          *graph;
    rvec             *f_global, *p, *sf;
    double            gpa, gpb, gpc, tmp, minstep;
    real              fnormn;
    real              stepsize;
    real              a, b, c, beta = 0.0;
    real              epot_repl = 0;
    real              pnorm;
    t_mdebin         *mdebin;
    gmx_bool          converged, foundlower;
    rvec              mu_tot;
    gmx_bool          do_log = FALSE, do_ene = FALSE, do_x, do_f;
    tensor            vir, pres;
    int               number_steps, neval = 0, nstcg = inputrec->nstcgsteep;
    gmx_mdoutf_t      outf;
    int               i, m, gf, step, nminstep;

    step = 0;

    s_min = init_em_state();
    s_a   = init_em_state();
    s_b   = init_em_state();
    s_c   = init_em_state();

    /* Init em and store the local state in s_min */
    init_em(fplog, CG, cr, inputrec,
            state_global, top_global, s_min, &top, &f, &f_global,
            nrnb, mu_tot, fr, &enerd, &graph, mdatoms, &gstat, vsite, constr,
            nfile, fnm, &outf, &mdebin, imdport, Flags, wcycle);

    /* Print to log file */
    print_em_start(fplog, cr, walltime_accounting, wcycle, CG);

    /* Max number of steps */
    number_steps = inputrec->nsteps;

    if (MASTER(cr))
    {
        sp_header(stderr, CG, inputrec->em_tol, number_steps);
    }
    if (fplog)
    {
        sp_header(fplog, CG, inputrec->em_tol, number_steps);
    }

    /* Call the force routine and some auxiliary (neighboursearching etc.) */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in congrad.c
     */
    evaluate_energy(fplog, cr,
                    top_global, s_min, top,
                    inputrec, nrnb, wcycle, gstat,
                    vsite, constr, fcd, graph, mdatoms, fr,
                    mu_tot, enerd, vir, pres, -1, TRUE);
    where();

    if (MASTER(cr))
    {
        /* Copy stuff to the energy bin for easy printing etc. */
        upd_mdebin(mdebin, FALSE, FALSE, (double)step,
                   mdatoms->tmass, enerd, &s_min->s, inputrec->fepvals, inputrec->expandedvals, s_min->s.box,
                   NULL, NULL, vir, pres, NULL, mu_tot, constr);

        print_ebin_header(fplog, step, step, s_min->s.lambda[efptFEP]);
        print_ebin(mdoutf_get_fp_ene(outf), TRUE, FALSE, FALSE, fplog, step, step, eprNORMAL,
                   TRUE, mdebin, fcd, &(top_global->groups), &(inputrec->opts));
    }
    where();

    /* Estimate/guess the initial stepsize */
    stepsize = inputrec->em_stepsize/s_min->fnorm;

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        fprintf(stderr, "   F-max             = %12.5e on atom %d\n",
                s_min->fmax, s_min->a_fmax+1);
        fprintf(stderr, "   F-Norm            = %12.5e\n",
                s_min->fnorm/sqrtNumAtoms);
        fprintf(stderr, "\n");
        /* and copy to the log file too... */
        fprintf(fplog, "   F-max             = %12.5e on atom %d\n",
                s_min->fmax, s_min->a_fmax+1);
        fprintf(fplog, "   F-Norm            = %12.5e\n",
                s_min->fnorm/sqrtNumAtoms);
        fprintf(fplog, "\n");
    }
    /* Start the loop over CG steps.
     * Each successful step is counted, and we continue until
     * we either converge or reach the max number of steps.
     */
    converged = FALSE;
    for (step = 0; (number_steps < 0 || (number_steps >= 0 && step <= number_steps)) && !converged; step++)
    {

        /* start taking steps in a new direction
         * First time we enter the routine, beta=0, and the direction is
         * simply the negative gradient.
         */

        /* Calculate the new direction in p, and the gradient in this direction, gpa */
        p   = s_min->s.cg_p;
        sf  = s_min->f;
        gpa = 0;
        gf  = 0;
        for (i = 0; i < mdatoms->homenr; i++)
        {
            if (mdatoms->cFREEZE)
            {
                gf = mdatoms->cFREEZE[i];
            }
            for (m = 0; m < DIM; m++)
            {
                if (!inputrec->opts.nFreeze[gf][m])
                {
                    p[i][m] = sf[i][m] + beta*p[i][m];
                    gpa    -= p[i][m]*sf[i][m];
                    /* f is negative gradient, thus the sign */
                }
                else
                {
                    p[i][m] = 0;
                }
            }
        }

        /* Sum the gradient along the line across CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &gpa, cr);
        }

        /* Calculate the norm of the search vector */
        get_f_norm_max(cr, &(inputrec->opts), mdatoms, p, &pnorm, NULL, NULL);

        /* Just in case stepsize reaches zero due to numerical precision... */
        if (stepsize <= 0)
        {
            stepsize = inputrec->em_stepsize/pnorm;
        }

        /*
         * Double check the value of the derivative in the search direction.
         * If it is positive it must be due to the old information in the
         * CG formula, so just remove that and start over with beta=0.
         * This corresponds to a steepest descent step.
         */
        if (gpa > 0)
        {
            beta = 0;
            step--;   /* Don't count this step since we are restarting */
            continue; /* Go back to the beginning of the big for-loop */
        }

        /* Calculate minimum allowed stepsize, before the average (norm)
         * relative change in coordinate is smaller than precision
         */
        minstep = 0;
        for (i = 0; i < mdatoms->homenr; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                tmp = fabs(s_min->s.x[i][m]);
                if (tmp < 1.0)
                {
                    tmp = 1.0;
                }
                tmp      = p[i][m]/tmp;
                minstep += tmp*tmp;
            }
        }
        /* Add up from all CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &minstep, cr);
        }

        minstep = GMX_REAL_EPS/sqrt(minstep/(3*state_global->natoms));

        if (stepsize < minstep)
        {
            converged = TRUE;
            break;
        }

        /* Write coordinates if necessary */
        do_x = do_per_step(step, inputrec->nstxout);
        do_f = do_per_step(step, inputrec->nstfout);

        write_em_traj(fplog, cr, outf, do_x, do_f, NULL,
                      top_global, inputrec, step,
                      s_min, state_global, f_global);

        /* Take a step downhill.
         * In theory, we should minimize the function along this direction.
         * That is quite possible, but it turns out to take 5-10 function evaluations
         * for each line. However, we dont really need to find the exact minimum -
         * it is much better to start a new CG step in a modified direction as soon
         * as we are close to it. This will save a lot of energy evaluations.
         *
         * In practice, we just try to take a single step.
         * If it worked (i.e. lowered the energy), we increase the stepsize but
         * the continue straight to the next CG step without trying to find any minimum.
         * If it didn't work (higher energy), there must be a minimum somewhere between
         * the old position and the new one.
         *
         * Due to the finite numerical accuracy, it turns out that it is a good idea
         * to even accept a SMALL increase in energy, if the derivative is still downhill.
         * This leads to lower final energies in the tests I've done. / Erik
         */
        s_a->epot = s_min->epot;
        a         = 0.0;
        c         = a + stepsize; /* reference position along line is zero */

        if (DOMAINDECOMP(cr) && s_min->s.ddp_count < cr->dd->ddp_count)
        {
            em_dd_partition_system(fplog, step, cr, top_global, inputrec,
                                   s_min, top, mdatoms, fr, vsite, constr,
                                   nrnb, wcycle);
        }

        /* Take a trial step (new coords in s_c) */
        do_em_step(cr, inputrec, mdatoms, fr->bMolPBC, s_min, c, s_min->s.cg_p, s_c,
                   constr, top, nrnb, wcycle, -1);

        neval++;
        /* Calculate energy for the trial step */
        evaluate_energy(fplog, cr,
                        top_global, s_c, top,
                        inputrec, nrnb, wcycle, gstat,
                        vsite, constr, fcd, graph, mdatoms, fr,
                        mu_tot, enerd, vir, pres, -1, FALSE);

        /* Calc derivative along line */
        p   = s_c->s.cg_p;
        sf  = s_c->f;
        gpc = 0;
        for (i = 0; i < mdatoms->homenr; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                gpc -= p[i][m]*sf[i][m]; /* f is negative gradient, thus the sign */
            }
        }
        /* Sum the gradient along the line across CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &gpc, cr);
        }

        /* This is the max amount of increase in energy we tolerate */
        tmp = sqrt(GMX_REAL_EPS)*fabs(s_a->epot);

        /* Accept the step if the energy is lower, or if it is not significantly higher
         * and the line derivative is still negative.
         */
        if (s_c->epot < s_a->epot || (gpc < 0 && s_c->epot < (s_a->epot + tmp)))
        {
            foundlower = TRUE;
            /* Great, we found a better energy. Increase step for next iteration
             * if we are still going down, decrease it otherwise
             */
            if (gpc < 0)
            {
                stepsize *= 1.618034; /* The golden section */
            }
            else
            {
                stepsize *= 0.618034; /* 1/golden section */
            }
        }
        else
        {
            /* New energy is the same or higher. We will have to do some work
             * to find a smaller value in the interval. Take smaller step next time!
             */
            foundlower = FALSE;
            stepsize  *= 0.618034;
        }




        /* OK, if we didn't find a lower value we will have to locate one now - there must
         * be one in the interval [a=0,c].
         * The same thing is valid here, though: Don't spend dozens of iterations to find
         * the line minimum. We try to interpolate based on the derivative at the endpoints,
         * and only continue until we find a lower value. In most cases this means 1-2 iterations.
         *
         * I also have a safeguard for potentially really pathological functions so we never
         * take more than 20 steps before we give up ...
         *
         * If we already found a lower value we just skip this step and continue to the update.
         */
        if (!foundlower)
        {
            nminstep = 0;

            do
            {
                /* Select a new trial point.
                 * If the derivatives at points a & c have different sign we interpolate to zero,
                 * otherwise just do a bisection.
                 */
                if (gpa < 0 && gpc > 0)
                {
                    b = a + gpa*(a-c)/(gpc-gpa);
                }
                else
                {
                    b = 0.5*(a+c);
                }

                /* safeguard if interpolation close to machine accuracy causes errors:
                 * never go outside the interval
                 */
                if (b <= a || b >= c)
                {
                    b = 0.5*(a+c);
                }

                if (DOMAINDECOMP(cr) && s_min->s.ddp_count != cr->dd->ddp_count)
                {
                    /* Reload the old state */
                    em_dd_partition_system(fplog, -1, cr, top_global, inputrec,
                                           s_min, top, mdatoms, fr, vsite, constr,
                                           nrnb, wcycle);
                }

                /* Take a trial step to this new point - new coords in s_b */
                do_em_step(cr, inputrec, mdatoms, fr->bMolPBC, s_min, b, s_min->s.cg_p, s_b,
                           constr, top, nrnb, wcycle, -1);

                neval++;
                /* Calculate energy for the trial step */
                evaluate_energy(fplog, cr,
                                top_global, s_b, top,
                                inputrec, nrnb, wcycle, gstat,
                                vsite, constr, fcd, graph, mdatoms, fr,
                                mu_tot, enerd, vir, pres, -1, FALSE);

                /* p does not change within a step, but since the domain decomposition
                 * might change, we have to use cg_p of s_b here.
                 */
                p   = s_b->s.cg_p;
                sf  = s_b->f;
                gpb = 0;
                for (i = 0; i < mdatoms->homenr; i++)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        gpb -= p[i][m]*sf[i][m]; /* f is negative gradient, thus the sign */
                    }
                }
                /* Sum the gradient along the line across CPUs */
                if (PAR(cr))
                {
                    gmx_sumd(1, &gpb, cr);
                }

                if (debug)
                {
                    fprintf(debug, "CGE: EpotA %f EpotB %f EpotC %f gpb %f\n",
                            s_a->epot, s_b->epot, s_c->epot, gpb);
                }

                epot_repl = s_b->epot;

                /* Keep one of the intervals based on the value of the derivative at the new point */
                if (gpb > 0)
                {
                    /* Replace c endpoint with b */
                    swap_em_state(s_b, s_c);
                    c   = b;
                    gpc = gpb;
                }
                else
                {
                    /* Replace a endpoint with b */
                    swap_em_state(s_b, s_a);
                    a   = b;
                    gpa = gpb;
                }

                /*
                 * Stop search as soon as we find a value smaller than the endpoints.
                 * Never run more than 20 steps, no matter what.
                 */
                nminstep++;
            }
            while ((epot_repl > s_a->epot || epot_repl > s_c->epot) &&
                   (nminstep < 20));

            if (fabs(epot_repl - s_min->epot) < fabs(s_min->epot)*GMX_REAL_EPS ||
                nminstep >= 20)
            {
                /* OK. We couldn't find a significantly lower energy.
                 * If beta==0 this was steepest descent, and then we give up.
                 * If not, set beta=0 and restart with steepest descent before quitting.
                 */
                if (beta == 0.0)
                {
                    /* Converged */
                    converged = TRUE;
                    break;
                }
                else
                {
                    /* Reset memory before giving up */
                    beta = 0.0;
                    continue;
                }
            }

            /* Select min energy state of A & C, put the best in B.
             */
            if (s_c->epot < s_a->epot)
            {
                if (debug)
                {
                    fprintf(debug, "CGE: C (%f) is lower than A (%f), moving C to B\n",
                            s_c->epot, s_a->epot);
                }
                swap_em_state(s_b, s_c);
                gpb = gpc;
            }
            else
            {
                if (debug)
                {
                    fprintf(debug, "CGE: A (%f) is lower than C (%f), moving A to B\n",
                            s_a->epot, s_c->epot);
                }
                swap_em_state(s_b, s_a);
                gpb = gpa;
            }

        }
        else
        {
            if (debug)
            {
                fprintf(debug, "CGE: Found a lower energy %f, moving C to B\n",
                        s_c->epot);
            }
            swap_em_state(s_b, s_c);
            gpb = gpc;
        }

        /* new search direction */
        /* beta = 0 means forget all memory and restart with steepest descents. */
        if (nstcg && ((step % nstcg) == 0))
        {
            beta = 0.0;
        }
        else
        {
            /* s_min->fnorm cannot be zero, because then we would have converged
             * and broken out.
             */

            /* Polak-Ribiere update.
             * Change to fnorm2/fnorm2_old for Fletcher-Reeves
             */
            beta = pr_beta(cr, &inputrec->opts, mdatoms, top_global, s_min, s_b);
        }
        /* Limit beta to prevent oscillations */
        if (fabs(beta) > 5.0)
        {
            beta = 0.0;
        }


        /* update positions */
        swap_em_state(s_min, s_b);
        gpa = gpb;

        /* Print it if necessary */
        if (MASTER(cr))
        {
            if (bVerbose)
            {
                double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
                fprintf(stderr, "\rStep %d, Epot=%12.6e, Fnorm=%9.3e, Fmax=%9.3e (atom %d)\n",
                        step, s_min->epot, s_min->fnorm/sqrtNumAtoms,
                        s_min->fmax, s_min->a_fmax+1);
            }
            /* Store the new (lower) energies */
            upd_mdebin(mdebin, FALSE, FALSE, (double)step,
                       mdatoms->tmass, enerd, &s_min->s, inputrec->fepvals, inputrec->expandedvals, s_min->s.box,
                       NULL, NULL, vir, pres, NULL, mu_tot, constr);

            do_log = do_per_step(step, inputrec->nstlog);
            do_ene = do_per_step(step, inputrec->nstenergy);

            /* Prepare IMD energy record, if bIMD is TRUE. */
            IMD_fill_energy_record(inputrec->bIMD, inputrec->imd, enerd, step, TRUE);

            if (do_log)
            {
                print_ebin_header(fplog, step, step, s_min->s.lambda[efptFEP]);
            }
            print_ebin(mdoutf_get_fp_ene(outf), do_ene, FALSE, FALSE,
                       do_log ? fplog : NULL, step, step, eprNORMAL,
                       TRUE, mdebin, fcd, &(top_global->groups), &(inputrec->opts));
        }

        /* Send energies and positions to the IMD client if bIMD is TRUE. */
        if (do_IMD(inputrec->bIMD, step, cr, TRUE, state_global->box, state_global->x, inputrec, 0, wcycle) && MASTER(cr))
        {
            IMD_send_positions(inputrec->imd);
        }

        /* Stop when the maximum force lies below tolerance.
         * If we have reached machine precision, converged is already set to true.
         */
        converged = converged || (s_min->fmax < inputrec->em_tol);

    } /* End of the loop */

    /* IMD cleanup, if bIMD is TRUE. */
    IMD_finalize(inputrec->bIMD, inputrec->imd);

    if (converged)
    {
        step--; /* we never took that last step in this case */

    }
    if (s_min->fmax > inputrec->em_tol)
    {
        if (MASTER(cr))
        {
            warn_step(stderr, inputrec->em_tol, step-1 == number_steps, FALSE);
            warn_step(fplog, inputrec->em_tol, step-1 == number_steps, FALSE);
        }
        converged = FALSE;
    }

    if (MASTER(cr))
    {
        /* If we printed energy and/or logfile last step (which was the last step)
         * we don't have to do it again, but otherwise print the final values.
         */
        if (!do_log)
        {
            /* Write final value to log since we didn't do anything the last step */
            print_ebin_header(fplog, step, step, s_min->s.lambda[efptFEP]);
        }
        if (!do_ene || !do_log)
        {
            /* Write final energy file entries */
            print_ebin(mdoutf_get_fp_ene(outf), !do_ene, FALSE, FALSE,
                       !do_log ? fplog : NULL, step, step, eprNORMAL,
                       TRUE, mdebin, fcd, &(top_global->groups), &(inputrec->opts));
        }
    }

    /* Print some stuff... */
    if (MASTER(cr))
    {
        fprintf(stderr, "\nwriting lowest energy coordinates.\n");
    }

    /* IMPORTANT!
     * For accurate normal mode calculation it is imperative that we
     * store the last conformation into the full precision binary trajectory.
     *
     * However, we should only do it if we did NOT already write this step
     * above (which we did if do_x or do_f was true).
     */
    do_x = !do_per_step(step, inputrec->nstxout);
    do_f = (inputrec->nstfout > 0 && !do_per_step(step, inputrec->nstfout));

    write_em_traj(fplog, cr, outf, do_x, do_f, ftp2fn(efSTO, nfile, fnm),
                  top_global, inputrec, step,
                  s_min, state_global, f_global);


    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        fnormn = s_min->fnorm/sqrtNumAtoms;
        print_converged(stderr, CG, inputrec->em_tol, step, converged, number_steps,
                        s_min->epot, s_min->fmax, s_min->a_fmax, fnormn);
        print_converged(fplog, CG, inputrec->em_tol, step, converged, number_steps,
                        s_min->epot, s_min->fmax, s_min->a_fmax, fnormn);

        fprintf(fplog, "\nPerformed %d energy evaluations in total.\n", neval);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    /* To print the actual number of steps we needed somewhere */
    walltime_accounting_set_nsteps_done(walltime_accounting, step);

    return 0;
} /* That's all folks */


double do_lbfgs(FILE *fplog, t_commrec *cr,
                int nfile, const t_filenm fnm[],
                const output_env_t gmx_unused oenv, gmx_bool bVerbose, gmx_bool gmx_unused bCompact,
                int gmx_unused nstglobalcomm,
                gmx_vsite_t *vsite, gmx_constr_t constr,
                int gmx_unused stepout,
                t_inputrec *inputrec,
                gmx_mtop_t *top_global, t_fcdata *fcd,
                t_state *state,
                t_mdatoms *mdatoms,
                t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                gmx_edsam_t gmx_unused ed,
                t_forcerec *fr,
                int gmx_unused repl_ex_nst, int gmx_unused repl_ex_nex, int gmx_unused repl_ex_seed,
                gmx_membed_t gmx_unused membed,
                real gmx_unused cpt_period, real gmx_unused max_hours,
                int imdport,
                unsigned long gmx_unused Flags,
                gmx_walltime_accounting_t walltime_accounting)
{
    static const char *LBFGS = "Low-Memory BFGS Minimizer";
    em_state_t         ems;
    gmx_localtop_t    *top;
    gmx_enerdata_t    *enerd;
    rvec              *f;
    gmx_global_stat_t  gstat;
    t_graph           *graph;
    rvec              *f_global;
    int                ncorr, nmaxcorr, point, cp, neval, nminstep;
    double             stepsize, step_taken, gpa, gpb, gpc, tmp, minstep;
    real              *rho, *alpha, *ff, *xx, *p, *s, *lastx, *lastf, **dx, **dg;
    real              *xa, *xb, *xc, *fa, *fb, *fc, *xtmp, *ftmp;
    real               a, b, c, maxdelta, delta;
    real               diag, Epot0, Epot, EpotA, EpotB, EpotC;
    real               dgdx, dgdg, sq, yr, beta;
    t_mdebin          *mdebin;
    gmx_bool           converged;
    rvec               mu_tot;
    real               fnorm, fmax;
    gmx_bool           do_log, do_ene, do_x, do_f, foundlower, *frozen;
    tensor             vir, pres;
    int                start, end, number_steps;
    gmx_mdoutf_t       outf;
    int                i, k, m, n, nfmax, gf, step;
    int                mdof_flags;

    if (PAR(cr))
    {
        gmx_fatal(FARGS, "Cannot do parallel L-BFGS Minimization - yet.\n");
    }

    if (NULL != constr)
    {
        gmx_fatal(FARGS, "The combination of constraints and L-BFGS minimization is not implemented. Either do not use constraints, or use another minimizer (e.g. steepest descent).");
    }

    n        = 3*state->natoms;
    nmaxcorr = inputrec->nbfgscorr;

    /* Allocate memory */
    /* Use pointers to real so we dont have to loop over both atoms and
     * dimensions all the time...
     * x/f are allocated as rvec *, so make new x0/f0 pointers-to-real
     * that point to the same memory.
     */
    snew(xa, n);
    snew(xb, n);
    snew(xc, n);
    snew(fa, n);
    snew(fb, n);
    snew(fc, n);
    snew(frozen, n);

    snew(p, n);
    snew(lastx, n);
    snew(lastf, n);
    snew(rho, nmaxcorr);
    snew(alpha, nmaxcorr);

    snew(dx, nmaxcorr);
    for (i = 0; i < nmaxcorr; i++)
    {
        snew(dx[i], n);
    }

    snew(dg, nmaxcorr);
    for (i = 0; i < nmaxcorr; i++)
    {
        snew(dg[i], n);
    }

    step  = 0;
    neval = 0;

    /* Init em */
    init_em(fplog, LBFGS, cr, inputrec,
            state, top_global, &ems, &top, &f, &f_global,
            nrnb, mu_tot, fr, &enerd, &graph, mdatoms, &gstat, vsite, constr,
            nfile, fnm, &outf, &mdebin, imdport, Flags, wcycle);
    /* Do_lbfgs is not completely updated like do_steep and do_cg,
     * so we free some memory again.
     */
    sfree(ems.s.x);
    sfree(ems.f);

    xx = (real *)state->x;
    ff = (real *)f;

    start = 0;
    end   = mdatoms->homenr;

    /* Print to log file */
    print_em_start(fplog, cr, walltime_accounting, wcycle, LBFGS);

    do_log = do_ene = do_x = do_f = TRUE;

    /* Max number of steps */
    number_steps = inputrec->nsteps;

    /* Create a 3*natoms index to tell whether each degree of freedom is frozen */
    gf = 0;
    for (i = start; i < end; i++)
    {
        if (mdatoms->cFREEZE)
        {
            gf = mdatoms->cFREEZE[i];
        }
        for (m = 0; m < DIM; m++)
        {
            frozen[3*i+m] = inputrec->opts.nFreeze[gf][m];
        }
    }
    if (MASTER(cr))
    {
        sp_header(stderr, LBFGS, inputrec->em_tol, number_steps);
    }
    if (fplog)
    {
        sp_header(fplog, LBFGS, inputrec->em_tol, number_steps);
    }

    if (vsite)
    {
        construct_vsites(vsite, state->x, 1, NULL,
                         top->idef.iparams, top->idef.il,
                         fr->ePBC, fr->bMolPBC, cr, state->box);
    }

    /* Call the force routine and some auxiliary (neighboursearching etc.) */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole
     */
    neval++;
    ems.s.x = state->x;
    ems.f   = f;
    evaluate_energy(fplog, cr,
                    top_global, &ems, top,
                    inputrec, nrnb, wcycle, gstat,
                    vsite, constr, fcd, graph, mdatoms, fr,
                    mu_tot, enerd, vir, pres, -1, TRUE);
    where();

    if (MASTER(cr))
    {
        /* Copy stuff to the energy bin for easy printing etc. */
        upd_mdebin(mdebin, FALSE, FALSE, (double)step,
                   mdatoms->tmass, enerd, state, inputrec->fepvals, inputrec->expandedvals, state->box,
                   NULL, NULL, vir, pres, NULL, mu_tot, constr);

        print_ebin_header(fplog, step, step, state->lambda[efptFEP]);
        print_ebin(mdoutf_get_fp_ene(outf), TRUE, FALSE, FALSE, fplog, step, step, eprNORMAL,
                   TRUE, mdebin, fcd, &(top_global->groups), &(inputrec->opts));
    }
    where();

    /* This is the starting energy */
    Epot = enerd->term[F_EPOT];

    fnorm = ems.fnorm;
    fmax  = ems.fmax;
    nfmax = ems.a_fmax;

    /* Set the initial step.
     * since it will be multiplied by the non-normalized search direction
     * vector (force vector the first time), we scale it by the
     * norm of the force.
     */

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state->natoms));
        fprintf(stderr, "Using %d BFGS correction steps.\n\n", nmaxcorr);
        fprintf(stderr, "   F-max             = %12.5e on atom %d\n", fmax, nfmax+1);
        fprintf(stderr, "   F-Norm            = %12.5e\n", fnorm/sqrtNumAtoms);
        fprintf(stderr, "\n");
        /* and copy to the log file too... */
        fprintf(fplog, "Using %d BFGS correction steps.\n\n", nmaxcorr);
        fprintf(fplog, "   F-max             = %12.5e on atom %d\n", fmax, nfmax+1);
        fprintf(fplog, "   F-Norm            = %12.5e\n", fnorm/sqrtNumAtoms);
        fprintf(fplog, "\n");
    }

    // Point is an index to the memory of search directions, where 0 is the first one.
    point = 0;

    // Set initial search direction to the force (-gradient), or 0 for frozen particles.
    for (i = 0; i < n; i++)
    {
        if (!frozen[i])
        {
            dx[point][i] = ff[i]; /* Initial search direction */
        }
        else
        {
            dx[point][i] = 0;
        }
    }

    // Stepsize will be modified during the search, and actually it is not critical
    // (the main efficiency in the algorithm comes from changing directions), but
    // we still need an initial value, so estimate it as the inverse of the norm
    // so we take small steps where the potential fluctuates a lot.
    stepsize  = 1.0/fnorm;

    /* Start the loop over BFGS steps.
     * Each successful step is counted, and we continue until
     * we either converge or reach the max number of steps.
     */

    ncorr = 0;

    /* Set the gradient from the force */
    converged = FALSE;
    for (step = 0; (number_steps < 0 || (number_steps >= 0 && step <= number_steps)) && !converged; step++)
    {

        /* Write coordinates if necessary */
        do_x = do_per_step(step, inputrec->nstxout);
        do_f = do_per_step(step, inputrec->nstfout);

        mdof_flags = 0;
        if (do_x)
        {
            mdof_flags |= MDOF_X;
        }

        if (do_f)
        {
            mdof_flags |= MDOF_F;
        }

        if (inputrec->bIMD)
        {
            mdof_flags |= MDOF_IMD;
        }

        mdoutf_write_to_trajectory_files(fplog, cr, outf, mdof_flags,
                                         top_global, step, (real)step, state, state, f, f);

        /* Do the linesearching in the direction dx[point][0..(n-1)] */

        /* make s a pointer to current search direction - point=0 first time we get here */
        s = dx[point];

        // calculate line gradient in position A
        for (gpa = 0, i = 0; i < n; i++)
        {
            gpa -= s[i]*ff[i];
        }

        /* Calculate minimum allowed stepsize along the line, before the average (norm)
         * relative change in coordinate is smaller than precision
         */
        for (minstep = 0, i = 0; i < n; i++)
        {
            tmp = fabs(xx[i]);
            if (tmp < 1.0)
            {
                tmp = 1.0;
            }
            tmp      = s[i]/tmp;
            minstep += tmp*tmp;
        }
        minstep = GMX_REAL_EPS/sqrt(minstep/n);

        if (stepsize < minstep)
        {
            converged = TRUE;
            break;
        }

        // Before taking any steps along the line, store the old position
        for (i = 0; i < n; i++)
        {
            lastx[i] = xx[i];
            lastf[i] = ff[i];
        }
        Epot0 = Epot;

        for (i = 0; i < n; i++)
        {
            xa[i] = xx[i];
        }

        /* Take a step downhill.
         * In theory, we should find the actual minimum of the function in this
         * direction, somewhere along the line.
         * That is quite possible, but it turns out to take 5-10 function evaluations
         * for each line. However, we dont really need to find the exact minimum -
         * it is much better to start a new BFGS step in a modified direction as soon
         * as we are close to it. This will save a lot of energy evaluations.
         *
         * In practice, we just try to take a single step.
         * If it worked (i.e. lowered the energy), we increase the stepsize but
         * continue straight to the next BFGS step without trying to find any minimum,
         * i.e. we change the search direction too. If the line was smooth, it is
         * likely we are in a smooth region, and then it makes sense to take longer
         * steps in the modified search direction too.
         *
         * If it didn't work (higher energy), there must be a minimum somewhere between
         * the old position and the new one. Then we need to start by finding a lower
         * value before we change search direction. Since the energy was apparently
         * quite rough, we need to decrease the step size.
         *
         * Due to the finite numerical accuracy, it turns out that it is a good idea
         * to accept a SMALL increase in energy, if the derivative is still downhill.
         * This leads to lower final energies in the tests I've done. / Erik
         */

        // State "A" is the first position along the line.
        // reference position along line is initially zero
        EpotA      = Epot0;
        a          = 0.0;

        // Check stepsize first. We do not allow displacements
        // larger than emstep.
        //
        do
        {
            // Pick a new position C by adding stepsize to A.
            c        = a + stepsize;

            // Calculate what the largest change in any individual coordinate
            // would be (translation along line * gradient along line)
            maxdelta = 0;
            for (i = 0; i < n; i++)
            {
                delta = c*s[i];
                if (delta > maxdelta)
                {
                    maxdelta = delta;
                }
            }
            // If any displacement is larger than the stepsize limit, reduce the step
            if (maxdelta > inputrec->em_stepsize)
            {
                stepsize *= 0.1;
            }
        }
        while (maxdelta > inputrec->em_stepsize);

        // Take a trial step and move the coordinate array xc[] to position C
        for (i = 0; i < n; i++)
        {
            xc[i] = lastx[i] + c*s[i];
        }

        neval++;
        // Calculate energy for the trial step in position C
        ems.s.x = (rvec *)xc;
        ems.f   = (rvec *)fc;
        evaluate_energy(fplog, cr,
                        top_global, &ems, top,
                        inputrec, nrnb, wcycle, gstat,
                        vsite, constr, fcd, graph, mdatoms, fr,
                        mu_tot, enerd, vir, pres, step, FALSE);
        EpotC = ems.epot;

        // Calc line gradient in position C
        for (gpc = 0, i = 0; i < n; i++)
        {
            gpc -= s[i]*fc[i]; /* f is negative gradient, thus the sign */
        }
        /* Sum the gradient along the line across CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &gpc, cr);
        }

        // This is the max amount of increase in energy we tolerate.
        // By allowing VERY small changes (close to numerical precision) we
        // frequently find even better (lower) final energies.
        tmp = sqrt(GMX_REAL_EPS)*fabs(EpotA);

        // Accept the step if the energy is lower in the new position C (compared to A),
        // or if it is not significantly higher and the line derivative is still negative.
        if (EpotC < EpotA || (gpc < 0 && EpotC < (EpotA+tmp)))
        {
            // Great, we found a better energy. We no longer try to alter the
            // stepsize, but simply accept this new better position. The we select a new
            // search direction instead, which will be much more efficient than continuing
            // to take smaller steps along a line. Set fnorm based on the new C position,
            // which will be used to update the stepsize to 1/fnorm further down.
            foundlower = TRUE;
            fnorm      = ems.fnorm;
        }
        else
        {
            // If we got here, the energy is NOT lower in point C, i.e. it will be the same
            // or higher than in point A. In this case it is pointless to move to point C,
            // so we will have to do more iterations along the same line to find a smaller
            // value in the interval [A=0.0,C].
            // Here, A is still 0.0, but that will change when we do a search in the interval
            // [0.0,C] below. That search we will do by interpolation or bisection rather
            // than with the stepsize, so no need to modify it. For the next search direction
            // it will be reset to 1/fnorm anyway.
            foundlower = FALSE;
        }

        if (!foundlower)
        {
            // OK, if we didn't find a lower value we will have to locate one now - there must
            // be one in the interval [a,c].
            // The same thing is valid here, though: Don't spend dozens of iterations to find
            // the line minimum. We try to interpolate based on the derivative at the endpoints,
            // and only continue until we find a lower value. In most cases this means 1-2 iterations.
            // I also have a safeguard for potentially really pathological functions so we never
            // take more than 20 steps before we give up.
            // If we already found a lower value we just skip this step and continue to the update.
            nminstep = 0;
            do
            {
                // Select a new trial point B in the interval [A,C].
                // If the derivatives at points a & c have different sign we interpolate to zero,
                // otherwise just do a bisection since there might be multiple minima/maxima
                // inside the interval.
                if (gpa < 0 && gpc > 0)
                {
                    b = a + gpa*(a-c)/(gpc-gpa);
                }
                else
                {
                    b = 0.5*(a+c);
                }

                /* safeguard if interpolation close to machine accuracy causes errors:
                 * never go outside the interval
                 */
                if (b <= a || b >= c)
                {
                    b = 0.5*(a+c);
                }

                // Take a trial step to point B
                for (i = 0; i < n; i++)
                {
                    xb[i] = lastx[i] + b*s[i];
                }

                neval++;
                // Calculate energy for the trial step in point B
                ems.s.x = (rvec *)xb;
                ems.f   = (rvec *)fb;
                evaluate_energy(fplog, cr,
                                top_global, &ems, top,
                                inputrec, nrnb, wcycle, gstat,
                                vsite, constr, fcd, graph, mdatoms, fr,
                                mu_tot, enerd, vir, pres, step, FALSE);
                EpotB = ems.epot;
                fnorm = ems.fnorm;

                // Calculate gradient in point B
                for (gpb = 0, i = 0; i < n; i++)
                {
                    gpb -= s[i]*fb[i]; /* f is negative gradient, thus the sign */

                }
                /* Sum the gradient along the line across CPUs */
                if (PAR(cr))
                {
                    gmx_sumd(1, &gpb, cr);
                }

                // Keep one of the intervals [A,B] or [B,C] based on the value of the derivative
                // at the new point B, and rename the endpoints of this new interval A and C.
                if (gpb > 0)
                {
                    /* Replace c endpoint with b */
                    EpotC = EpotB;
                    c     = b;
                    gpc   = gpb;
                    /* swap coord pointers b/c */
                    xtmp = xb;
                    ftmp = fb;
                    xb   = xc;
                    fb   = fc;
                    xc   = xtmp;
                    fc   = ftmp;
                }
                else
                {
                    /* Replace a endpoint with b */
                    EpotA = EpotB;
                    a     = b;
                    gpa   = gpb;
                    /* swap coord pointers a/b */
                    xtmp = xb;
                    ftmp = fb;
                    xb   = xa;
                    fb   = fa;
                    xa   = xtmp;
                    fa   = ftmp;
                }

                /*
                 * Stop search as soon as we find a value smaller than the endpoints,
                 * or if the tolerance is below machine precision.
                 * Never run more than 20 steps, no matter what.
                 */
                nminstep++;
            }
            while ((EpotB > EpotA || EpotB > EpotC) && (nminstep < 20));

            if (fabs(EpotB-Epot0) < GMX_REAL_EPS || nminstep >= 20)
            {
                /* OK. We couldn't find a significantly lower energy.
                 * If ncorr==0 this was steepest descent, and then we give up.
                 * If not, reset memory to restart as steepest descent before quitting.
                 */
                if (ncorr == 0)
                {
                    /* Converged */
                    converged = TRUE;
                    break;
                }
                else
                {
                    /* Reset memory */
                    ncorr = 0;
                    /* Search in gradient direction */
                    for (i = 0; i < n; i++)
                    {
                        dx[point][i] = ff[i];
                    }
                    /* Reset stepsize */
                    stepsize = 1.0/fnorm;
                    continue;
                }
            }

            /* Select min energy state of A & C, put the best in xx/ff/Epot
             */
            if (EpotC < EpotA)
            {
                Epot = EpotC;
                /* Use state C */
                for (i = 0; i < n; i++)
                {
                    xx[i] = xc[i];
                    ff[i] = fc[i];
                }
                step_taken = c;
            }
            else
            {
                Epot = EpotA;
                /* Use state A */
                for (i = 0; i < n; i++)
                {
                    xx[i] = xa[i];
                    ff[i] = fa[i];
                }
                step_taken = a;
            }

        }
        else
        {
            /* found lower */
            Epot = EpotC;
            /* Use state C */
            for (i = 0; i < n; i++)
            {
                xx[i] = xc[i];
                ff[i] = fc[i];
            }
            step_taken = c;
        }

        /* Update the memory information, and calculate a new
         * approximation of the inverse hessian
         */

        /* Have new data in Epot, xx, ff */
        if (ncorr < nmaxcorr)
        {
            ncorr++;
        }

        for (i = 0; i < n; i++)
        {
            dg[point][i]  = lastf[i]-ff[i];
            dx[point][i] *= step_taken;
        }

        dgdg = 0;
        dgdx = 0;
        for (i = 0; i < n; i++)
        {
            dgdg += dg[point][i]*dg[point][i];
            dgdx += dg[point][i]*dx[point][i];
        }

        diag = dgdx/dgdg;

        rho[point] = 1.0/dgdx;
        point++;

        if (point >= nmaxcorr)
        {
            point = 0;
        }

        /* Update */
        for (i = 0; i < n; i++)
        {
            p[i] = ff[i];
        }

        cp = point;

        /* Recursive update. First go back over the memory points */
        for (k = 0; k < ncorr; k++)
        {
            cp--;
            if (cp < 0)
            {
                cp = ncorr-1;
            }

            sq = 0;
            for (i = 0; i < n; i++)
            {
                sq += dx[cp][i]*p[i];
            }

            alpha[cp] = rho[cp]*sq;

            for (i = 0; i < n; i++)
            {
                p[i] -= alpha[cp]*dg[cp][i];
            }
        }

        for (i = 0; i < n; i++)
        {
            p[i] *= diag;
        }

        /* And then go forward again */
        for (k = 0; k < ncorr; k++)
        {
            yr = 0;
            for (i = 0; i < n; i++)
            {
                yr += p[i]*dg[cp][i];
            }

            beta = rho[cp]*yr;
            beta = alpha[cp]-beta;

            for (i = 0; i < n; i++)
            {
                p[i] += beta*dx[cp][i];
            }

            cp++;
            if (cp >= ncorr)
            {
                cp = 0;
            }
        }

        for (i = 0; i < n; i++)
        {
            if (!frozen[i])
            {
                dx[point][i] = p[i];
            }
            else
            {
                dx[point][i] = 0;
            }
        }

        /* Test whether the convergence criterion is met */
        get_f_norm_max(cr, &(inputrec->opts), mdatoms, f, &fnorm, &fmax, &nfmax);

        /* Print it if necessary */
        if (MASTER(cr))
        {
            if (bVerbose)
            {
                double sqrtNumAtoms = sqrt(static_cast<double>(state->natoms));
                fprintf(stderr, "\rStep %d, Epot=%12.6e, Fnorm=%9.3e, Fmax=%9.3e (atom %d)\n",
                        step, Epot, fnorm/sqrtNumAtoms, fmax, nfmax+1);
            }
            /* Store the new (lower) energies */
            upd_mdebin(mdebin, FALSE, FALSE, (double)step,
                       mdatoms->tmass, enerd, state, inputrec->fepvals, inputrec->expandedvals, state->box,
                       NULL, NULL, vir, pres, NULL, mu_tot, constr);
            do_log = do_per_step(step, inputrec->nstlog);
            do_ene = do_per_step(step, inputrec->nstenergy);
            if (do_log)
            {
                print_ebin_header(fplog, step, step, state->lambda[efptFEP]);
            }
            print_ebin(mdoutf_get_fp_ene(outf), do_ene, FALSE, FALSE,
                       do_log ? fplog : NULL, step, step, eprNORMAL,
                       TRUE, mdebin, fcd, &(top_global->groups), &(inputrec->opts));
        }

        /* Send x and E to IMD client, if bIMD is TRUE. */
        if (do_IMD(inputrec->bIMD, step, cr, TRUE, state->box, state->x, inputrec, 0, wcycle) && MASTER(cr))
        {
            IMD_send_positions(inputrec->imd);
        }

        // Reset stepsize in we are doing more iterations
        stepsize = 1.0/fnorm;

        /* Stop when the maximum force lies below tolerance.
         * If we have reached machine precision, converged is already set to true.
         */
        converged = converged || (fmax < inputrec->em_tol);

    } /* End of the loop */

    /* IMD cleanup, if bIMD is TRUE. */
    IMD_finalize(inputrec->bIMD, inputrec->imd);

    if (converged)
    {
        step--; /* we never took that last step in this case */

    }
    if (fmax > inputrec->em_tol)
    {
        if (MASTER(cr))
        {
            warn_step(stderr, inputrec->em_tol, step-1 == number_steps, FALSE);
            warn_step(fplog, inputrec->em_tol, step-1 == number_steps, FALSE);
        }
        converged = FALSE;
    }

    /* If we printed energy and/or logfile last step (which was the last step)
     * we don't have to do it again, but otherwise print the final values.
     */
    if (!do_log) /* Write final value to log since we didn't do anythin last step */
    {
        print_ebin_header(fplog, step, step, state->lambda[efptFEP]);
    }
    if (!do_ene || !do_log) /* Write final energy file entries */
    {
        print_ebin(mdoutf_get_fp_ene(outf), !do_ene, FALSE, FALSE,
                   !do_log ? fplog : NULL, step, step, eprNORMAL,
                   TRUE, mdebin, fcd, &(top_global->groups), &(inputrec->opts));
    }

    /* Print some stuff... */
    if (MASTER(cr))
    {
        fprintf(stderr, "\nwriting lowest energy coordinates.\n");
    }

    /* IMPORTANT!
     * For accurate normal mode calculation it is imperative that we
     * store the last conformation into the full precision binary trajectory.
     *
     * However, we should only do it if we did NOT already write this step
     * above (which we did if do_x or do_f was true).
     */
    do_x = !do_per_step(step, inputrec->nstxout);
    do_f = !do_per_step(step, inputrec->nstfout);
    write_em_traj(fplog, cr, outf, do_x, do_f, ftp2fn(efSTO, nfile, fnm),
                  top_global, inputrec, step,
                  &ems, state, f);

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state->natoms));
        print_converged(stderr, LBFGS, inputrec->em_tol, step, converged,
                        number_steps, Epot, fmax, nfmax, fnorm/sqrtNumAtoms);
        print_converged(fplog, LBFGS, inputrec->em_tol, step, converged,
                        number_steps, Epot, fmax, nfmax, fnorm/sqrtNumAtoms);

        fprintf(fplog, "\nPerformed %d energy evaluations in total.\n", neval);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    /* To print the actual number of steps we needed somewhere */
    walltime_accounting_set_nsteps_done(walltime_accounting, step);

    return 0;
} /* That's all folks */


double do_steep(FILE *fplog, t_commrec *cr,
                int nfile, const t_filenm fnm[],
                const output_env_t gmx_unused oenv, gmx_bool bVerbose, gmx_bool gmx_unused bCompact,
                int gmx_unused nstglobalcomm,
                gmx_vsite_t *vsite, gmx_constr_t constr,
                int gmx_unused stepout,
                t_inputrec *inputrec,
                gmx_mtop_t *top_global, t_fcdata *fcd,
                t_state *state_global,
                t_mdatoms *mdatoms,
                t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                gmx_edsam_t gmx_unused  ed,
                t_forcerec *fr,
                int gmx_unused repl_ex_nst, int gmx_unused repl_ex_nex, int gmx_unused repl_ex_seed,
                gmx_membed_t gmx_unused membed,
                real gmx_unused cpt_period, real gmx_unused max_hours,
                int imdport,
                unsigned long gmx_unused Flags,
                gmx_walltime_accounting_t walltime_accounting)
{
    const char       *SD = "Steepest Descents";
    em_state_t       *s_min, *s_try;
    rvec             *f_global;
    gmx_localtop_t   *top;
    gmx_enerdata_t   *enerd;
    rvec             *f;
    gmx_global_stat_t gstat;
    t_graph          *graph;
    real              stepsize;
    real              ustep, fnormn;
    gmx_mdoutf_t      outf;
    t_mdebin         *mdebin;
    gmx_bool          bDone, bAbort, do_x, do_f;
    tensor            vir, pres;
    rvec              mu_tot;
    int               nsteps;
    int               count          = 0;
    int               steps_accepted = 0;

    s_min = init_em_state();
    s_try = init_em_state();

    /* Init em and store the local state in s_try */
    init_em(fplog, SD, cr, inputrec,
            state_global, top_global, s_try, &top, &f, &f_global,
            nrnb, mu_tot, fr, &enerd, &graph, mdatoms, &gstat, vsite, constr,
            nfile, fnm, &outf, &mdebin, imdport, Flags, wcycle);

    /* Print to log file  */
    print_em_start(fplog, cr, walltime_accounting, wcycle, SD);

    /* Set variables for stepsize (in nm). This is the largest
     * step that we are going to make in any direction.
     */
    ustep    = inputrec->em_stepsize;
    stepsize = 0;

    /* Max number of steps  */
    nsteps = inputrec->nsteps;

    if (MASTER(cr))
    {
        /* Print to the screen  */
        sp_header(stderr, SD, inputrec->em_tol, nsteps);
    }
    if (fplog)
    {
        sp_header(fplog, SD, inputrec->em_tol, nsteps);
    }

    /**** HERE STARTS THE LOOP ****
     * count is the counter for the number of steps
     * bDone will be TRUE when the minimization has converged
     * bAbort will be TRUE when nsteps steps have been performed or when
     * the stepsize becomes smaller than is reasonable for machine precision
     */
    count  = 0;
    bDone  = FALSE;
    bAbort = FALSE;
    while (!bDone && !bAbort)
    {
        bAbort = (nsteps >= 0) && (count == nsteps);

        /* set new coordinates, except for first step */
        if (count > 0)
        {
            do_em_step(cr, inputrec, mdatoms, fr->bMolPBC,
                       s_min, stepsize, s_min->f, s_try,
                       constr, top, nrnb, wcycle, count);
        }

        evaluate_energy(fplog, cr,
                        top_global, s_try, top,
                        inputrec, nrnb, wcycle, gstat,
                        vsite, constr, fcd, graph, mdatoms, fr,
                        mu_tot, enerd, vir, pres, count, count == 0);

        if (MASTER(cr))
        {
            print_ebin_header(fplog, count, count, s_try->s.lambda[efptFEP]);
        }

        if (count == 0)
        {
            s_min->epot = s_try->epot;
        }

        /* Print it if necessary  */
        if (MASTER(cr))
        {
            if (bVerbose)
            {
                fprintf(stderr, "Step=%5d, Dmax= %6.1e nm, Epot= %12.5e Fmax= %11.5e, atom= %d%c",
                        count, ustep, s_try->epot, s_try->fmax, s_try->a_fmax+1,
                        ( (count == 0) || (s_try->epot < s_min->epot) ) ? '\n' : '\r');
            }

            if ( (count == 0) || (s_try->epot < s_min->epot) )
            {
                /* Store the new (lower) energies  */
                upd_mdebin(mdebin, FALSE, FALSE, (double)count,
                           mdatoms->tmass, enerd, &s_try->s, inputrec->fepvals, inputrec->expandedvals,
                           s_try->s.box, NULL, NULL, vir, pres, NULL, mu_tot, constr);

                /* Prepare IMD energy record, if bIMD is TRUE. */
                IMD_fill_energy_record(inputrec->bIMD, inputrec->imd, enerd, count, TRUE);

                print_ebin(mdoutf_get_fp_ene(outf), TRUE,
                           do_per_step(steps_accepted, inputrec->nstdisreout),
                           do_per_step(steps_accepted, inputrec->nstorireout),
                           fplog, count, count, eprNORMAL, TRUE,
                           mdebin, fcd, &(top_global->groups), &(inputrec->opts));
                fflush(fplog);
            }
        }

        /* Now if the new energy is smaller than the previous...
         * or if this is the first step!
         * or if we did random steps!
         */

        if ( (count == 0) || (s_try->epot < s_min->epot) )
        {
            steps_accepted++;

            /* Test whether the convergence criterion is met...  */
            bDone = (s_try->fmax < inputrec->em_tol);

            /* Copy the arrays for force, positions and energy  */
            /* The 'Min' array always holds the coords and forces of the minimal
               sampled energy  */
            swap_em_state(s_min, s_try);
            if (count > 0)
            {
                ustep *= 1.2;
            }

            /* Write to trn, if necessary */
            do_x = do_per_step(steps_accepted, inputrec->nstxout);
            do_f = do_per_step(steps_accepted, inputrec->nstfout);
            write_em_traj(fplog, cr, outf, do_x, do_f, NULL,
                          top_global, inputrec, count,
                          s_min, state_global, f_global);
        }
        else
        {
            /* If energy is not smaller make the step smaller...  */
            ustep *= 0.5;

            if (DOMAINDECOMP(cr) && s_min->s.ddp_count != cr->dd->ddp_count)
            {
                /* Reload the old state */
                em_dd_partition_system(fplog, count, cr, top_global, inputrec,
                                       s_min, top, mdatoms, fr, vsite, constr,
                                       nrnb, wcycle);
            }
        }

        /* Determine new step  */
        stepsize = ustep/s_min->fmax;

        /* Check if stepsize is too small, with 1 nm as a characteristic length */
#ifdef GMX_DOUBLE
        if (count == nsteps || ustep < 1e-12)
#else
        if (count == nsteps || ustep < 1e-6)
#endif
        {
            if (MASTER(cr))
            {
                warn_step(stderr, inputrec->em_tol, count == nsteps, constr != NULL);
                warn_step(fplog, inputrec->em_tol, count == nsteps, constr != NULL);
            }
            bAbort = TRUE;
        }

        /* Send IMD energies and positions, if bIMD is TRUE. */
        if (do_IMD(inputrec->bIMD, count, cr, TRUE, state_global->box, state_global->x, inputrec, 0, wcycle) && MASTER(cr))
        {
            IMD_send_positions(inputrec->imd);
        }

        count++;
    } /* End of the loop  */

    /* IMD cleanup, if bIMD is TRUE. */
    IMD_finalize(inputrec->bIMD, inputrec->imd);

    /* Print some data...  */
    if (MASTER(cr))
    {
        fprintf(stderr, "\nwriting lowest energy coordinates.\n");
    }
    write_em_traj(fplog, cr, outf, TRUE, inputrec->nstfout, ftp2fn(efSTO, nfile, fnm),
                  top_global, inputrec, count,
                  s_min, state_global, f_global);

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        fnormn = s_min->fnorm/sqrtNumAtoms;

        print_converged(stderr, SD, inputrec->em_tol, count, bDone, nsteps,
                        s_min->epot, s_min->fmax, s_min->a_fmax, fnormn);
        print_converged(fplog, SD, inputrec->em_tol, count, bDone, nsteps,
                        s_min->epot, s_min->fmax, s_min->a_fmax, fnormn);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    /* To print the actual number of steps we needed somewhere */
    inputrec->nsteps = count;

    walltime_accounting_set_nsteps_done(walltime_accounting, count);

    return 0;
} /* That's all folks */


double do_nm(FILE *fplog, t_commrec *cr,
             int nfile, const t_filenm fnm[],
             const output_env_t gmx_unused oenv, gmx_bool bVerbose, gmx_bool gmx_unused  bCompact,
             int gmx_unused nstglobalcomm,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int gmx_unused stepout,
             t_inputrec *inputrec,
             gmx_mtop_t *top_global, t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb, gmx_wallcycle_t wcycle,
             gmx_edsam_t  gmx_unused ed,
             t_forcerec *fr,
             int gmx_unused repl_ex_nst, int gmx_unused repl_ex_nex, int gmx_unused repl_ex_seed,
             gmx_membed_t gmx_unused membed,
             real gmx_unused cpt_period, real gmx_unused max_hours,
             int imdport,
             unsigned long gmx_unused Flags,
             gmx_walltime_accounting_t walltime_accounting)
{
    const char          *NM = "Normal Mode Analysis";
    gmx_mdoutf_t         outf;
    int                  natoms, atom, d;
    int                  nnodes, node;
    rvec                *f_global;
    gmx_localtop_t      *top;
    gmx_enerdata_t      *enerd;
    rvec                *f;
    gmx_global_stat_t    gstat;
    t_graph             *graph;
    tensor               vir, pres;
    rvec                 mu_tot;
    rvec                *fneg, *dfdx;
    gmx_bool             bSparse; /* use sparse matrix storage format */
    size_t               sz = 0;
    gmx_sparsematrix_t * sparse_matrix           = NULL;
    real           *     full_matrix             = NULL;
    em_state_t       *   state_work;

    /* added with respect to mdrun */
    int        i, j, k, row, col;
    real       der_range = 10.0*sqrt(GMX_REAL_EPS);
    real       x_min;
    bool       bIsMaster = MASTER(cr);

    if (constr != NULL)
    {
        gmx_fatal(FARGS, "Constraints present with Normal Mode Analysis, this combination is not supported");
    }

    state_work = init_em_state();

    /* Init em and store the local state in state_minimum */
    init_em(fplog, NM, cr, inputrec,
            state_global, top_global, state_work, &top,
            &f, &f_global,
            nrnb, mu_tot, fr, &enerd, &graph, mdatoms, &gstat, vsite, constr,
            nfile, fnm, &outf, NULL, imdport, Flags, wcycle);

    natoms = top_global->natoms;
    snew(fneg, natoms);
    snew(dfdx, natoms);

#ifndef GMX_DOUBLE
    if (bIsMaster)
    {
        fprintf(stderr,
                "NOTE: This version of GROMACS has been compiled in single precision,\n"
                "      which MIGHT not be accurate enough for normal mode analysis.\n"
                "      GROMACS now uses sparse matrix storage, so the memory requirements\n"
                "      are fairly modest even if you recompile in double precision.\n\n");
    }
#endif

    /* Check if we can/should use sparse storage format.
     *
     * Sparse format is only useful when the Hessian itself is sparse, which it
     * will be when we use a cutoff.
     * For small systems (n<1000) it is easier to always use full matrix format, though.
     */
    if (EEL_FULL(fr->eeltype) || fr->rlist == 0.0)
    {
        md_print_info(cr, fplog, "Non-cutoff electrostatics used, forcing full Hessian format.\n");
        bSparse = FALSE;
    }
    else if (top_global->natoms < 1000)
    {
        md_print_info(cr, fplog, "Small system size (N=%d), using full Hessian format.\n", top_global->natoms);
        bSparse = FALSE;
    }
    else
    {
        md_print_info(cr, fplog, "Using compressed symmetric sparse Hessian format.\n");
        bSparse = TRUE;
    }

    if (bIsMaster)
    {
        sz = DIM*top_global->natoms;

        fprintf(stderr, "Allocating Hessian memory...\n\n");

        if (bSparse)
        {
            sparse_matrix = gmx_sparsematrix_init(sz);
            sparse_matrix->compressed_symmetric = TRUE;
        }
        else
        {
            snew(full_matrix, sz*sz);
        }
    }

    init_nrnb(nrnb);

    where();

    /* Write start time and temperature */
    print_em_start(fplog, cr, walltime_accounting, wcycle, NM);

    /* fudge nr of steps to nr of atoms */
    inputrec->nsteps = natoms*2;

    if (bIsMaster)
    {
        fprintf(stderr, "starting normal mode calculation '%s'\n%d steps.\n\n",
                *(top_global->name), (int)inputrec->nsteps);
    }

    nnodes = cr->nnodes;

    /* Make evaluate_energy do a single node force calculation */
    cr->nnodes = 1;
    evaluate_energy(fplog, cr,
                    top_global, state_work, top,
                    inputrec, nrnb, wcycle, gstat,
                    vsite, constr, fcd, graph, mdatoms, fr,
                    mu_tot, enerd, vir, pres, -1, TRUE);
    cr->nnodes = nnodes;

    /* if forces are not small, warn user */
    get_state_f_norm_max(cr, &(inputrec->opts), mdatoms, state_work);

    md_print_info(cr, fplog, "Maximum force:%12.5e\n", state_work->fmax);
    if (state_work->fmax > 1.0e-3)
    {
        md_print_info(cr, fplog,
                      "The force is probably not small enough to "
                      "ensure that you are at a minimum.\n"
                      "Be aware that negative eigenvalues may occur\n"
                      "when the resulting matrix is diagonalized.\n\n");
    }

    /***********************************************************
     *
     *      Loop over all pairs in matrix
     *
     *      do_force called twice. Once with positive and
     *      once with negative displacement
     *
     ************************************************************/

    /* Steps are divided one by one over the nodes */
    for (atom = cr->nodeid; atom < natoms; atom += nnodes)
    {

        for (d = 0; d < DIM; d++)
        {
            x_min = state_work->s.x[atom][d];

            state_work->s.x[atom][d] = x_min - der_range;

            /* Make evaluate_energy do a single node force calculation */
            cr->nnodes = 1;
            evaluate_energy(fplog, cr,
                            top_global, state_work, top,
                            inputrec, nrnb, wcycle, gstat,
                            vsite, constr, fcd, graph, mdatoms, fr,
                            mu_tot, enerd, vir, pres, atom*2, FALSE);

            for (i = 0; i < natoms; i++)
            {
                copy_rvec(state_work->f[i], fneg[i]);
            }

            state_work->s.x[atom][d] = x_min + der_range;

            evaluate_energy(fplog, cr,
                            top_global, state_work, top,
                            inputrec, nrnb, wcycle, gstat,
                            vsite, constr, fcd, graph, mdatoms, fr,
                            mu_tot, enerd, vir, pres, atom*2+1, FALSE);
            cr->nnodes = nnodes;

            /* x is restored to original */
            state_work->s.x[atom][d] = x_min;

            for (j = 0; j < natoms; j++)
            {
                for (k = 0; (k < DIM); k++)
                {
                    dfdx[j][k] =
                        -(state_work->f[j][k] - fneg[j][k])/(2*der_range);
                }
            }

            if (!bIsMaster)
            {
#ifdef GMX_MPI
#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif
                MPI_Send(dfdx[0], natoms*DIM, mpi_type, MASTERNODE(cr), cr->nodeid,
                         cr->mpi_comm_mygroup);
#endif
            }
            else
            {
                for (node = 0; (node < nnodes && atom+node < natoms); node++)
                {
                    if (node > 0)
                    {
#ifdef GMX_MPI
                        MPI_Status stat;
                        MPI_Recv(dfdx[0], natoms*DIM, mpi_type, node, node,
                                 cr->mpi_comm_mygroup, &stat);
#undef mpi_type
#endif
                    }

                    row = (atom + node)*DIM + d;

                    for (j = 0; j < natoms; j++)
                    {
                        for (k = 0; k < DIM; k++)
                        {
                            col = j*DIM + k;

                            if (bSparse)
                            {
                                if (col >= row && dfdx[j][k] != 0.0)
                                {
                                    gmx_sparsematrix_increment_value(sparse_matrix,
                                                                     row, col, dfdx[j][k]);
                                }
                            }
                            else
                            {
                                full_matrix[row*sz+col] = dfdx[j][k];
                            }
                        }
                    }
                }
            }

            if (bVerbose && fplog)
            {
                fflush(fplog);
            }
        }
        /* write progress */
        if (bIsMaster && bVerbose)
        {
            fprintf(stderr, "\rFinished step %d out of %d",
                    std::min(atom+nnodes, natoms), natoms);
            fflush(stderr);
        }
    }

    if (bIsMaster)
    {
        fprintf(stderr, "\n\nWriting Hessian...\n");
        gmx_mtxio_write(ftp2fn(efMTX, nfile, fnm), sz, sz, full_matrix, sparse_matrix);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    walltime_accounting_set_nsteps_done(walltime_accounting, natoms*2);

    return 0;
}
