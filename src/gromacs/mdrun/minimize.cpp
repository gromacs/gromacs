/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019,2020, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file defines integrators for energy minimization
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "config.h"

#include <cmath>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/imd/imd.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/listed_forces/manage_threading.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"

#include "legacysimulator.h"
#include "shellfc.h"

using gmx::MdrunScheduleWorkload;

//! Utility structure for manipulating states during EM
typedef struct
{
    //! Copy of the global state
    t_state s;
    //! Force array
    PaddedHostVector<gmx::RVec> f;
    //! Potential energy
    real epot;
    //! Norm of the force
    real fnorm;
    //! Maximum force
    real fmax;
    //! Direction
    int a_fmax;
} em_state_t;

//! Print the EM starting conditions
static void print_em_start(FILE*                     fplog,
                           const t_commrec*          cr,
                           gmx_walltime_accounting_t walltime_accounting,
                           gmx_wallcycle_t           wcycle,
                           const char*               name)
{
    walltime_accounting_start_time(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, name);
}

//! Stop counting time for EM
static void em_time_end(gmx_walltime_accounting_t walltime_accounting, gmx_wallcycle_t wcycle)
{
    wallcycle_stop(wcycle, ewcRUN);

    walltime_accounting_end_time(walltime_accounting);
}

//! Printing a log file and console header
static void sp_header(FILE* out, const char* minimizer, real ftol, int nsteps)
{
    fprintf(out, "\n");
    fprintf(out, "%s:\n", minimizer);
    fprintf(out, "   Tolerance (Fmax)   = %12.5e\n", ftol);
    fprintf(out, "   Number of steps    = %12d\n", nsteps);
}

//! Print warning message
static void warn_step(FILE* fp, real ftol, real fmax, gmx_bool bLastStep, gmx_bool bConstrain)
{
    constexpr bool realIsDouble = GMX_DOUBLE;
    char           buffer[2048];

    if (!std::isfinite(fmax))
    {
        sprintf(buffer,
                "\nEnergy minimization has stopped because the force "
                "on at least one atom is not finite. This usually means "
                "atoms are overlapping. Modify the input coordinates to "
                "remove atom overlap or use soft-core potentials with "
                "the free energy code to avoid infinite forces.\n%s",
                !realIsDouble ? "You could also be lucky that switching to double precision "
                                "is sufficient to obtain finite forces.\n"
                              : "");
    }
    else if (bLastStep)
    {
        sprintf(buffer,
                "\nEnergy minimization reached the maximum number "
                "of steps before the forces reached the requested "
                "precision Fmax < %g.\n",
                ftol);
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
                !realIsDouble ? "\nDouble precision normally gives you higher accuracy, but "
                                "this is often not needed for preparing to run molecular "
                                "dynamics.\n"
                              : "",
                bConstrain ? "You might need to increase your constraint accuracy, or turn\n"
                             "off constraints altogether (set constraints = none in mdp file)\n"
                           : "");
    }

    fputs(wrap_lines(buffer, 78, 0, FALSE), stderr);
    fputs(wrap_lines(buffer, 78, 0, FALSE), fp);
}

//! Print message about convergence of the EM
static void print_converged(FILE*             fp,
                            const char*       alg,
                            real              ftol,
                            int64_t           count,
                            gmx_bool          bDone,
                            int64_t           nsteps,
                            const em_state_t* ems,
                            double            sqrtNumAtoms)
{
    char buf[STEPSTRSIZE];

    if (bDone)
    {
        fprintf(fp, "\n%s converged to Fmax < %g in %s steps\n", alg, ftol, gmx_step_str(count, buf));
    }
    else if (count < nsteps)
    {
        fprintf(fp,
                "\n%s converged to machine precision in %s steps,\n"
                "but did not reach the requested Fmax < %g.\n",
                alg, gmx_step_str(count, buf), ftol);
    }
    else
    {
        fprintf(fp, "\n%s did not converge to Fmax < %g in %s steps.\n", alg, ftol,
                gmx_step_str(count, buf));
    }

#if GMX_DOUBLE
    fprintf(fp, "Potential Energy  = %21.14e\n", ems->epot);
    fprintf(fp, "Maximum force     = %21.14e on atom %d\n", ems->fmax, ems->a_fmax + 1);
    fprintf(fp, "Norm of force     = %21.14e\n", ems->fnorm / sqrtNumAtoms);
#else
    fprintf(fp, "Potential Energy  = %14.7e\n", ems->epot);
    fprintf(fp, "Maximum force     = %14.7e on atom %d\n", ems->fmax, ems->a_fmax + 1);
    fprintf(fp, "Norm of force     = %14.7e\n", ems->fnorm / sqrtNumAtoms);
#endif
}

//! Compute the norm and max of the force array in parallel
static void get_f_norm_max(const t_commrec* cr,
                           t_grpopts*       opts,
                           t_mdatoms*       mdatoms,
                           const rvec*      f,
                           real*            fnorm,
                           real*            fmax,
                           int*             a_fmax)
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
                    fam += gmx::square(f[i][m]);
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
            fam = norm2(f[i]);
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
        a_max = cr->dd->globalAtomIndices[la_max];
    }
    else
    {
        a_max = la_max;
    }
    if (PAR(cr))
    {
        snew(sum, 2 * cr->nnodes + 1);
        sum[2 * cr->nodeid]     = fmax2;
        sum[2 * cr->nodeid + 1] = a_max;
        sum[2 * cr->nnodes]     = fnorm2;
        gmx_sumd(2 * cr->nnodes + 1, sum, cr);
        fnorm2 = sum[2 * cr->nnodes];
        /* Determine the global maximum */
        for (i = 0; i < cr->nnodes; i++)
        {
            if (sum[2 * i] > fmax2)
            {
                fmax2 = sum[2 * i];
                a_max = gmx::roundToInt(sum[2 * i + 1]);
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
        *fmax = sqrt(fmax2);
    }
    if (a_fmax)
    {
        *a_fmax = a_max;
    }
}

//! Compute the norm of the force
static void get_state_f_norm_max(const t_commrec* cr, t_grpopts* opts, t_mdatoms* mdatoms, em_state_t* ems)
{
    get_f_norm_max(cr, opts, mdatoms, ems->f.rvec_array(), &ems->fnorm, &ems->fmax, &ems->a_fmax);
}

//! Initialize the energy minimization
static void init_em(FILE*                fplog,
                    const gmx::MDLogger& mdlog,
                    const char*          title,
                    const t_commrec*     cr,
                    t_inputrec*          ir,
                    gmx::ImdSession*     imdSession,
                    pull_t*              pull_work,
                    t_state*             state_global,
                    gmx_mtop_t*          top_global,
                    em_state_t*          ems,
                    gmx_localtop_t*      top,
                    t_nrnb*              nrnb,
                    t_forcerec*          fr,
                    t_graph**            graph,
                    gmx::MDAtoms*        mdAtoms,
                    gmx_global_stat_t*   gstat,
                    gmx_vsite_t*         vsite,
                    gmx::Constraints*    constr,
                    gmx_shellfc_t**      shellfc)
{
    real dvdl_constr;

    if (fplog)
    {
        fprintf(fplog, "Initiating %s\n", title);
    }

    if (MASTER(cr))
    {
        state_global->ngtc = 0;
    }
    initialize_lambdas(fplog, *ir, MASTER(cr), &(state_global->fep_state), state_global->lambda, nullptr);

    if (ir->eI == eiNM)
    {
        GMX_ASSERT(shellfc != nullptr, "With NM we always support shells");

        *shellfc = init_shell_flexcon(stdout, top_global, constr ? constr->numFlexibleConstraints() : 0,
                                      ir->nstcalcenergy, DOMAINDECOMP(cr));
    }
    else
    {
        GMX_ASSERT(EI_ENERGY_MINIMIZATION(ir->eI),
                   "This else currently only handles energy minimizers, consider if your algorithm "
                   "needs shell/flexible-constraint support");

        /* With energy minimization, shells and flexible constraints are
         * automatically minimized when treated like normal DOFS.
         */
        if (shellfc != nullptr)
        {
            *shellfc = nullptr;
        }
    }

    auto mdatoms = mdAtoms->mdatoms();
    if (DOMAINDECOMP(cr))
    {
        top->useInDomainDecomp_ = true;
        dd_init_local_top(*top_global, top);

        dd_init_local_state(cr->dd, state_global, &ems->s);

        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, mdlog, ir->init_step, cr, TRUE, 1, state_global, *top_global, ir,
                            imdSession, pull_work, &ems->s, &ems->f, mdAtoms, top, fr, vsite,
                            constr, nrnb, nullptr, FALSE);
        dd_store_state(cr->dd, &ems->s);

        *graph = nullptr;
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        /* Just copy the state */
        ems->s = *state_global;
        state_change_natoms(&ems->s, ems->s.natoms);
        ems->f.resizeWithPadding(ems->s.natoms);

        mdAlgorithmsSetupAtomData(cr, ir, *top_global, top, fr, graph, mdAtoms, constr, vsite,
                                  shellfc ? *shellfc : nullptr);

        if (vsite)
        {
            set_vsite_top(vsite, top, mdatoms);
        }
    }

    update_mdatoms(mdAtoms->mdatoms(), ems->s.lambda[efptMASS]);

    if (constr)
    {
        // TODO how should this cross-module support dependency be managed?
        if (ir->eConstrAlg == econtSHAKE && gmx_mtop_ftype_count(top_global, F_CONSTR) > 0)
        {
            gmx_fatal(FARGS, "Can not do energy minimization with %s, use %s\n",
                      econstr_names[econtSHAKE], econstr_names[econtLINCS]);
        }

        if (!ir->bContinuation)
        {
            /* Constrain the starting coordinates */
            dvdl_constr = 0;
            constr->apply(TRUE, TRUE, -1, 0, 1.0, ems->s.x.rvec_array(), ems->s.x.rvec_array(),
                          nullptr, ems->s.box, ems->s.lambda[efptFEP], &dvdl_constr, nullptr,
                          nullptr, gmx::ConstraintVariable::Positions);
        }
    }

    if (PAR(cr))
    {
        *gstat = global_stat_init(ir);
    }
    else
    {
        *gstat = nullptr;
    }

    calc_shifts(ems->s.box, fr->shift_vec);
}

//! Finalize the minimization
static void finish_em(const t_commrec*          cr,
                      gmx_mdoutf_t              outf,
                      gmx_walltime_accounting_t walltime_accounting,
                      gmx_wallcycle_t           wcycle)
{
    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    done_mdoutf(outf);

    em_time_end(walltime_accounting, wcycle);
}

//! Swap two different EM states during minimization
static void swap_em_state(em_state_t** ems1, em_state_t** ems2)
{
    em_state_t* tmp;

    tmp   = *ems1;
    *ems1 = *ems2;
    *ems2 = tmp;
}

//! Save the EM trajectory
static void write_em_traj(FILE*               fplog,
                          const t_commrec*    cr,
                          gmx_mdoutf_t        outf,
                          gmx_bool            bX,
                          gmx_bool            bF,
                          const char*         confout,
                          gmx_mtop_t*         top_global,
                          t_inputrec*         ir,
                          int64_t             step,
                          em_state_t*         state,
                          t_state*            state_global,
                          ObservablesHistory* observablesHistory)
{
    int mdof_flags = 0;

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

    mdoutf_write_to_trajectory_files(fplog, cr, outf, mdof_flags, top_global->natoms, step,
                                     static_cast<double>(step), &state->s, state_global,
                                     observablesHistory, state->f);

    if (confout != nullptr)
    {
        if (DOMAINDECOMP(cr))
        {
            /* If bX=true, x was collected to state_global in the call above */
            if (!bX)
            {
                auto globalXRef = MASTER(cr) ? state_global->x : gmx::ArrayRef<gmx::RVec>();
                dd_collect_vec(cr->dd, &state->s, state->s.x, globalXRef);
            }
        }
        else
        {
            /* Copy the local state pointer */
            state_global = &state->s;
        }

        if (MASTER(cr))
        {
            if (ir->ePBC != epbcNONE && !ir->bPeriodicMols && DOMAINDECOMP(cr))
            {
                /* Make molecules whole only for confout writing */
                do_pbc_mtop(ir->ePBC, state->s.box, top_global, state_global->x.rvec_array());
            }

            write_sto_conf_mtop(confout, *top_global->name, top_global,
                                state_global->x.rvec_array(), nullptr, ir->ePBC, state->s.box);
        }
    }
}

//! \brief Do one minimization step
//
// \returns true when the step succeeded, false when a constraint error occurred
static bool do_em_step(const t_commrec*                   cr,
                       t_inputrec*                        ir,
                       t_mdatoms*                         md,
                       em_state_t*                        ems1,
                       real                               a,
                       const PaddedHostVector<gmx::RVec>* force,
                       em_state_t*                        ems2,
                       gmx::Constraints*                  constr,
                       int64_t                            count)

{
    t_state *s1, *s2;
    int      start, end;
    real     dvdl_constr;
    int nthreads gmx_unused;

    bool validStep = true;

    s1 = &ems1->s;
    s2 = &ems2->s;

    if (DOMAINDECOMP(cr) && s1->ddp_count != cr->dd->ddp_count)
    {
        gmx_incons("state mismatch in do_em_step");
    }

    s2->flags = s1->flags;

    if (s2->natoms != s1->natoms)
    {
        state_change_natoms(s2, s1->natoms);
        ems2->f.resizeWithPadding(s2->natoms);
    }
    if (DOMAINDECOMP(cr) && s2->cg_gl.size() != s1->cg_gl.size())
    {
        s2->cg_gl.resize(s1->cg_gl.size());
    }

    copy_mat(s1->box, s2->box);
    /* Copy free energy state */
    s2->lambda = s1->lambda;
    copy_mat(s1->box, s2->box);

    start = 0;
    end   = md->homenr;

    nthreads = gmx_omp_nthreads_get(emntUpdate);
#pragma omp parallel num_threads(nthreads)
    {
        const rvec* x1 = s1->x.rvec_array();
        rvec*       x2 = s2->x.rvec_array();
        const rvec* f  = force->rvec_array();

        int gf = 0;
#pragma omp for schedule(static) nowait
        for (int i = start; i < end; i++)
        {
            try
            {
                if (md->cFREEZE)
                {
                    gf = md->cFREEZE[i];
                }
                for (int m = 0; m < DIM; m++)
                {
                    if (ir->opts.nFreeze[gf][m])
                    {
                        x2[i][m] = x1[i][m];
                    }
                    else
                    {
                        x2[i][m] = x1[i][m] + a * f[i][m];
                    }
                }
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }

        if (s2->flags & (1 << estCGP))
        {
            /* Copy the CG p vector */
            const rvec* p1 = s1->cg_p.rvec_array();
            rvec*       p2 = s2->cg_p.rvec_array();
#pragma omp for schedule(static) nowait
            for (int i = start; i < end; i++)
            {
                // Trivial OpenMP block that does not throw
                copy_rvec(p1[i], p2[i]);
            }
        }

        if (DOMAINDECOMP(cr))
        {
            /* OpenMP does not supported unsigned loop variables */
#pragma omp for schedule(static) nowait
            for (gmx::index i = 0; i < gmx::ssize(s2->cg_gl); i++)
            {
                s2->cg_gl[i] = s1->cg_gl[i];
            }
        }
    }

    if (DOMAINDECOMP(cr))
    {
        s2->ddp_count       = s1->ddp_count;
        s2->ddp_count_cg_gl = s1->ddp_count_cg_gl;
    }

    if (constr)
    {
        dvdl_constr = 0;
        validStep = constr->apply(TRUE, TRUE, count, 0, 1.0, s1->x.rvec_array(), s2->x.rvec_array(),
                                  nullptr, s2->box, s2->lambda[efptBONDED], &dvdl_constr, nullptr,
                                  nullptr, gmx::ConstraintVariable::Positions);

        if (cr->nnodes > 1)
        {
            /* This global reduction will affect performance at high
             * parallelization, but we can not really avoid it.
             * But usually EM is not run at high parallelization.
             */
            int reductionBuffer = static_cast<int>(!validStep);
            gmx_sumi(1, &reductionBuffer, cr);
            validStep = (reductionBuffer == 0);
        }

        // We should move this check to the different minimizers
        if (!validStep && ir->eI != eiSteep)
        {
            gmx_fatal(FARGS,
                      "The coordinates could not be constrained. Minimizer '%s' can not handle "
                      "constraint failures, use minimizer '%s' before using '%s'.",
                      EI(ir->eI), EI(eiSteep), EI(ir->eI));
        }
    }

    return validStep;
}

//! Prepare EM for using domain decomposition parallellization
static void em_dd_partition_system(FILE*                fplog,
                                   const gmx::MDLogger& mdlog,
                                   int                  step,
                                   const t_commrec*     cr,
                                   gmx_mtop_t*          top_global,
                                   t_inputrec*          ir,
                                   gmx::ImdSession*     imdSession,
                                   pull_t*              pull_work,
                                   em_state_t*          ems,
                                   gmx_localtop_t*      top,
                                   gmx::MDAtoms*        mdAtoms,
                                   t_forcerec*          fr,
                                   gmx_vsite_t*         vsite,
                                   gmx::Constraints*    constr,
                                   t_nrnb*              nrnb,
                                   gmx_wallcycle_t      wcycle)
{
    /* Repartition the domain decomposition */
    dd_partition_system(fplog, mdlog, step, cr, FALSE, 1, nullptr, *top_global, ir, imdSession, pull_work,
                        &ems->s, &ems->f, mdAtoms, top, fr, vsite, constr, nrnb, wcycle, FALSE);
    dd_store_state(cr->dd, &ems->s);
}

namespace
{

/*! \brief Class to handle the work of setting and doing an energy evaluation.
 *
 * This class is a mere aggregate of parameters to pass to evaluate an
 * energy, so that future changes to names and types of them consume
 * less time when refactoring other code.
 *
 * Aggregate initialization is used, for which the chief risk is that
 * if a member is added at the end and not all initializer lists are
 * updated, then the member will be value initialized, which will
 * typically mean initialization to zero.
 *
 * Use a braced initializer list to construct one of these. */
class EnergyEvaluator
{
public:
    /*! \brief Evaluates an energy on the state in \c ems.
     *
     * \todo In practice, the same objects mu_tot, vir, and pres
     * are always passed to this function, so we would rather have
     * them as data members. However, their C-array types are
     * unsuited for aggregate initialization. When the types
     * improve, the call signature of this method can be reduced.
     */
    void run(em_state_t* ems, rvec mu_tot, tensor vir, tensor pres, int64_t count, gmx_bool bFirst);
    //! Handles logging (deprecated).
    FILE* fplog;
    //! Handles logging.
    const gmx::MDLogger& mdlog;
    //! Handles communication.
    const t_commrec* cr;
    //! Coordinates multi-simulations.
    const gmx_multisim_t* ms;
    //! Holds the simulation topology.
    gmx_mtop_t* top_global;
    //! Holds the domain topology.
    gmx_localtop_t* top;
    //! User input options.
    t_inputrec* inputrec;
    //! The Interactive Molecular Dynamics session.
    gmx::ImdSession* imdSession;
    //! The pull work object.
    pull_t* pull_work;
    //! Manages flop accounting.
    t_nrnb* nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle_t wcycle;
    //! Coordinates global reduction.
    gmx_global_stat_t gstat;
    //! Handles virtual sites.
    gmx_vsite_t* vsite;
    //! Handles constraints.
    gmx::Constraints* constr;
    //! Handles strange things.
    t_fcdata* fcd;
    //! Molecular graph for SHAKE.
    t_graph* graph;
    //! Per-atom data for this domain.
    gmx::MDAtoms* mdAtoms;
    //! Handles how to calculate the forces.
    t_forcerec* fr;
    //! Schedule of force-calculation work each step for this task.
    MdrunScheduleWorkload* runScheduleWork;
    //! Stores the computed energies.
    gmx_enerdata_t* enerd;
};

void EnergyEvaluator::run(em_state_t* ems, rvec mu_tot, tensor vir, tensor pres, int64_t count, gmx_bool bFirst)
{
    real     t;
    gmx_bool bNS;
    tensor   force_vir, shake_vir, ekin;
    real     dvdl_constr;
    real     terminate = 0;

    /* Set the time to the initial time, the time does not change during EM */
    t = inputrec->init_t;

    if (bFirst || (DOMAINDECOMP(cr) && ems->s.ddp_count < cr->dd->ddp_count))
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
        construct_vsites(vsite, ems->s.x.rvec_array(), 1, nullptr, top->idef.iparams, top->idef.il,
                         fr->ePBC, fr->bMolPBC, cr, ems->s.box);
    }

    if (DOMAINDECOMP(cr) && bNS)
    {
        /* Repartition the domain decomposition */
        em_dd_partition_system(fplog, mdlog, count, cr, top_global, inputrec, imdSession, pull_work,
                               ems, top, mdAtoms, fr, vsite, constr, nrnb, wcycle);
    }

    /* Calc force & energy on new trial position  */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in congrad.c
     */
    do_force(fplog, cr, ms, inputrec, nullptr, nullptr, imdSession, pull_work, count, nrnb, wcycle,
             top, ems->s.box, ems->s.x.arrayRefWithPadding(), &ems->s.hist,
             ems->f.arrayRefWithPadding(), force_vir, mdAtoms->mdatoms(), enerd, fcd, ems->s.lambda,
             graph, fr, runScheduleWork, vsite, mu_tot, t, nullptr,
             GMX_FORCE_STATECHANGED | GMX_FORCE_ALLFORCES | GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY
                     | (bNS ? GMX_FORCE_NS : 0),
             DDBalanceRegionHandler(cr));

    /* Clear the unused shake virial and pressure */
    clear_mat(shake_vir);
    clear_mat(pres);

    /* Communicate stuff when parallel */
    if (PAR(cr) && inputrec->eI != eiNM)
    {
        wallcycle_start(wcycle, ewcMoveE);

        global_stat(gstat, cr, enerd, force_vir, shake_vir, mu_tot, inputrec, nullptr, nullptr, nullptr,
                    1, &terminate, nullptr, FALSE, CGLO_ENERGY | CGLO_PRESSURE | CGLO_CONSTRAINT);

        wallcycle_stop(wcycle, ewcMoveE);
    }

    if (fr->dispersionCorrection)
    {
        /* Calculate long range corrections to pressure and energy */
        const DispersionCorrection::Correction correction =
                fr->dispersionCorrection->calculate(ems->s.box, ems->s.lambda[efptVDW]);

        enerd->term[F_DISPCORR] = correction.energy;
        enerd->term[F_EPOT] += correction.energy;
        enerd->term[F_PRES] += correction.pressure;
        enerd->term[F_DVDL] += correction.dvdl;
    }
    else
    {
        enerd->term[F_DISPCORR] = 0;
    }

    ems->epot = enerd->term[F_EPOT];

    if (constr)
    {
        /* Project out the constraint components of the force */
        dvdl_constr  = 0;
        rvec* f_rvec = ems->f.rvec_array();
        constr->apply(FALSE, FALSE, count, 0, 1.0, ems->s.x.rvec_array(), f_rvec, f_rvec,
                      ems->s.box, ems->s.lambda[efptBONDED], &dvdl_constr, nullptr, &shake_vir,
                      gmx::ConstraintVariable::ForceDispl);
        enerd->term[F_DVDL_CONSTR] += dvdl_constr;
        m_add(force_vir, shake_vir, vir);
    }
    else
    {
        copy_mat(force_vir, vir);
    }

    clear_mat(ekin);
    enerd->term[F_PRES] = calc_pres(fr->ePBC, inputrec->nwall, ems->s.box, ekin, vir, pres);

    sum_dhdl(enerd, ems->s.lambda, *inputrec->fepvals);

    if (EI_ENERGY_MINIMIZATION(inputrec->eI))
    {
        get_state_f_norm_max(cr, &(inputrec->opts), mdAtoms->mdatoms(), ems);
    }
}

} // namespace

//! Parallel utility summing energies and forces
static double reorder_partsum(const t_commrec* cr,
                              t_grpopts*       opts,
                              gmx_mtop_t*      top_global,
                              em_state_t*      s_min,
                              em_state_t*      s_b)
{
    if (debug)
    {
        fprintf(debug, "Doing reorder_partsum\n");
    }

    const rvec* fm = s_min->f.rvec_array();
    const rvec* fb = s_b->f.rvec_array();

    /* Collect fm in a global vector fmg.
     * This conflicts with the spirit of domain decomposition,
     * but to fully optimize this a much more complicated algorithm is required.
     */
    const int natoms = top_global->natoms;
    rvec*     fmg;
    snew(fmg, natoms);

    gmx::ArrayRef<const int> indicesMin = s_min->s.cg_gl;
    int                      i          = 0;
    for (int a : indicesMin)
    {
        copy_rvec(fm[i], fmg[a]);
        i++;
    }
    gmx_sum(top_global->natoms * 3, fmg[0], cr);

    /* Now we will determine the part of the sum for the cgs in state s_b */
    gmx::ArrayRef<const int> indicesB = s_b->s.cg_gl;

    double partsum                  = 0;
    i                               = 0;
    int                          gf = 0;
    gmx::ArrayRef<unsigned char> grpnrFREEZE =
            top_global->groups.groupNumbers[SimulationAtomGroupType::Freeze];
    for (int a : indicesB)
    {
        if (!grpnrFREEZE.empty())
        {
            gf = grpnrFREEZE[i];
        }
        for (int m = 0; m < DIM; m++)
        {
            if (!opts->nFreeze[gf][m])
            {
                partsum += (fb[i][m] - fmg[a][m]) * fb[i][m];
            }
        }
        i++;
    }

    sfree(fmg);

    return partsum;
}

//! Print some stuff, like beta, whatever that means.
static real pr_beta(const t_commrec* cr,
                    t_grpopts*       opts,
                    t_mdatoms*       mdatoms,
                    gmx_mtop_t*      top_global,
                    em_state_t*      s_min,
                    em_state_t*      s_b)
{
    double sum;

    /* This is just the classical Polak-Ribiere calculation of beta;
     * it looks a bit complicated since we take freeze groups into account,
     * and might have to sum it in parallel runs.
     */

    if (!DOMAINDECOMP(cr)
        || (s_min->s.ddp_count == cr->dd->ddp_count && s_b->s.ddp_count == cr->dd->ddp_count))
    {
        const rvec* fm = s_min->f.rvec_array();
        const rvec* fb = s_b->f.rvec_array();
        sum            = 0;
        int gf         = 0;
        /* This part of code can be incorrect with DD,
         * since the atom ordering in s_b and s_min might differ.
         */
        for (int i = 0; i < mdatoms->homenr; i++)
        {
            if (mdatoms->cFREEZE)
            {
                gf = mdatoms->cFREEZE[i];
            }
            for (int m = 0; m < DIM; m++)
            {
                if (!opts->nFreeze[gf][m])
                {
                    sum += (fb[i][m] - fm[i][m]) * fb[i][m];
                }
            }
        }
    }
    else
    {
        /* We need to reorder cgs while summing */
        sum = reorder_partsum(cr, opts, top_global, s_min, s_b);
    }
    if (PAR(cr))
    {
        gmx_sumd(1, &sum, cr);
    }

    return sum / gmx::square(s_min->fnorm);
}

namespace gmx
{

void LegacySimulator::do_cg()
{
    const char* CG = "Polak-Ribiere Conjugate Gradients";

    gmx_localtop_t    top;
    gmx_global_stat_t gstat;
    t_graph*          graph;
    double            tmp, minstep;
    real              stepsize;
    real              a, b, c, beta = 0.0;
    real              epot_repl = 0;
    real              pnorm;
    gmx_bool          converged, foundlower;
    rvec              mu_tot = { 0 };
    gmx_bool          do_log = FALSE, do_ene = FALSE, do_x, do_f;
    tensor            vir, pres;
    int               number_steps, neval = 0, nstcg = inputrec->nstcgsteep;
    int               m, step, nminstep;
    auto              mdatoms = mdAtoms->mdatoms();

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Note that activating conjugate gradient energy minimization via the "
                    "integrator .mdp option and the command gmx mdrun may "
                    "be available in a different form in a future version of GROMACS, "
                    "e.g. gmx minimize and an .mdp option.");

    step = 0;

    if (MASTER(cr))
    {
        // In CG, the state is extended with a search direction
        state_global->flags |= (1 << estCGP);

        // Ensure the extra per-atom state array gets allocated
        state_change_natoms(state_global, state_global->natoms);

        // Initialize the search direction to zero
        for (RVec& cg_p : state_global->cg_p)
        {
            cg_p = { 0, 0, 0 };
        }
    }

    /* Create 4 states on the stack and extract pointers that we will swap */
    em_state_t  s0{}, s1{}, s2{}, s3{};
    em_state_t* s_min = &s0;
    em_state_t* s_a   = &s1;
    em_state_t* s_b   = &s2;
    em_state_t* s_c   = &s3;

    /* Init em and store the local state in s_min */
    init_em(fplog, mdlog, CG, cr, inputrec, imdSession, pull_work, state_global, top_global, s_min,
            &top, nrnb, fr, &graph, mdAtoms, &gstat, vsite, constr, nullptr);
    const bool        simulationsShareState = false;
    gmx_mdoutf*       outf = init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider,
                                   mdModulesNotifier, inputrec, top_global, nullptr, wcycle,
                                   StartingBehavior::NewSimulation, simulationsShareState, ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf), top_global, inputrec, pull_work, nullptr,
                                   false, StartingBehavior::NewSimulation, mdModulesNotifier);

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

    EnergyEvaluator energyEvaluator{
        fplog,      mdlog,     cr,      ms,     top_global,      &top,  inputrec,
        imdSession, pull_work, nrnb,    wcycle, gstat,           vsite, constr,
        fcd,        graph,     mdAtoms, fr,     runScheduleWork, enerd
    };
    /* Call the force routine and some auxiliary (neighboursearching etc.) */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in congrad.c
     */
    energyEvaluator.run(s_min, mu_tot, vir, pres, -1, TRUE);

    if (MASTER(cr))
    {
        /* Copy stuff to the energy bin for easy printing etc. */
        matrix nullBox = {};
        energyOutput.addDataAtEnergyStep(false, false, static_cast<double>(step), mdatoms->tmass,
                                         enerd, nullptr, nullptr, nullptr, nullBox, nullptr,
                                         nullptr, vir, pres, nullptr, mu_tot, constr);

        energyOutput.printHeader(fplog, step, step);
        energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), TRUE, FALSE, FALSE, fplog, step,
                                           step, fcd, nullptr);
    }

    /* Estimate/guess the initial stepsize */
    stepsize = inputrec->em_stepsize / s_min->fnorm;

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        fprintf(stderr, "   F-max             = %12.5e on atom %d\n", s_min->fmax, s_min->a_fmax + 1);
        fprintf(stderr, "   F-Norm            = %12.5e\n", s_min->fnorm / sqrtNumAtoms);
        fprintf(stderr, "\n");
        /* and copy to the log file too... */
        fprintf(fplog, "   F-max             = %12.5e on atom %d\n", s_min->fmax, s_min->a_fmax + 1);
        fprintf(fplog, "   F-Norm            = %12.5e\n", s_min->fnorm / sqrtNumAtoms);
        fprintf(fplog, "\n");
    }
    /* Start the loop over CG steps.
     * Each successful step is counted, and we continue until
     * we either converge or reach the max number of steps.
     */
    converged = FALSE;
    for (step = 0; (number_steps < 0 || step <= number_steps) && !converged; step++)
    {

        /* start taking steps in a new direction
         * First time we enter the routine, beta=0, and the direction is
         * simply the negative gradient.
         */

        /* Calculate the new direction in p, and the gradient in this direction, gpa */
        rvec*       pm  = s_min->s.cg_p.rvec_array();
        const rvec* sfm = s_min->f.rvec_array();
        double      gpa = 0;
        int         gf  = 0;
        for (int i = 0; i < mdatoms->homenr; i++)
        {
            if (mdatoms->cFREEZE)
            {
                gf = mdatoms->cFREEZE[i];
            }
            for (m = 0; m < DIM; m++)
            {
                if (!inputrec->opts.nFreeze[gf][m])
                {
                    pm[i][m] = sfm[i][m] + beta * pm[i][m];
                    gpa -= pm[i][m] * sfm[i][m];
                    /* f is negative gradient, thus the sign */
                }
                else
                {
                    pm[i][m] = 0;
                }
            }
        }

        /* Sum the gradient along the line across CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &gpa, cr);
        }

        /* Calculate the norm of the search vector */
        get_f_norm_max(cr, &(inputrec->opts), mdatoms, pm, &pnorm, nullptr, nullptr);

        /* Just in case stepsize reaches zero due to numerical precision... */
        if (stepsize <= 0)
        {
            stepsize = inputrec->em_stepsize / pnorm;
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
        minstep      = 0;
        auto s_min_x = makeArrayRef(s_min->s.x);
        for (int i = 0; i < mdatoms->homenr; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                tmp = fabs(s_min_x[i][m]);
                if (tmp < 1.0)
                {
                    tmp = 1.0;
                }
                tmp = pm[i][m] / tmp;
                minstep += tmp * tmp;
            }
        }
        /* Add up from all CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &minstep, cr);
        }

        minstep = GMX_REAL_EPS / sqrt(minstep / (3 * top_global->natoms));

        if (stepsize < minstep)
        {
            converged = TRUE;
            break;
        }

        /* Write coordinates if necessary */
        do_x = do_per_step(step, inputrec->nstxout);
        do_f = do_per_step(step, inputrec->nstfout);

        write_em_traj(fplog, cr, outf, do_x, do_f, nullptr, top_global, inputrec, step, s_min,
                      state_global, observablesHistory);

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
            em_dd_partition_system(fplog, mdlog, step, cr, top_global, inputrec, imdSession,
                                   pull_work, s_min, &top, mdAtoms, fr, vsite, constr, nrnb, wcycle);
        }

        /* Take a trial step (new coords in s_c) */
        do_em_step(cr, inputrec, mdatoms, s_min, c, &s_min->s.cg_p, s_c, constr, -1);

        neval++;
        /* Calculate energy for the trial step */
        energyEvaluator.run(s_c, mu_tot, vir, pres, -1, FALSE);

        /* Calc derivative along line */
        const rvec* pc  = s_c->s.cg_p.rvec_array();
        const rvec* sfc = s_c->f.rvec_array();
        double      gpc = 0;
        for (int i = 0; i < mdatoms->homenr; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                gpc -= pc[i][m] * sfc[i][m]; /* f is negative gradient, thus the sign */
            }
        }
        /* Sum the gradient along the line across CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &gpc, cr);
        }

        /* This is the max amount of increase in energy we tolerate */
        tmp = std::sqrt(GMX_REAL_EPS) * fabs(s_a->epot);

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
            stepsize *= 0.618034;
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
        double gpb;
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
                    b = a + gpa * (a - c) / (gpc - gpa);
                }
                else
                {
                    b = 0.5 * (a + c);
                }

                /* safeguard if interpolation close to machine accuracy causes errors:
                 * never go outside the interval
                 */
                if (b <= a || b >= c)
                {
                    b = 0.5 * (a + c);
                }

                if (DOMAINDECOMP(cr) && s_min->s.ddp_count != cr->dd->ddp_count)
                {
                    /* Reload the old state */
                    em_dd_partition_system(fplog, mdlog, -1, cr, top_global, inputrec, imdSession, pull_work,
                                           s_min, &top, mdAtoms, fr, vsite, constr, nrnb, wcycle);
                }

                /* Take a trial step to this new point - new coords in s_b */
                do_em_step(cr, inputrec, mdatoms, s_min, b, &s_min->s.cg_p, s_b, constr, -1);

                neval++;
                /* Calculate energy for the trial step */
                energyEvaluator.run(s_b, mu_tot, vir, pres, -1, FALSE);

                /* p does not change within a step, but since the domain decomposition
                 * might change, we have to use cg_p of s_b here.
                 */
                const rvec* pb  = s_b->s.cg_p.rvec_array();
                const rvec* sfb = s_b->f.rvec_array();
                gpb             = 0;
                for (int i = 0; i < mdatoms->homenr; i++)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        gpb -= pb[i][m] * sfb[i][m]; /* f is negative gradient, thus the sign */
                    }
                }
                /* Sum the gradient along the line across CPUs */
                if (PAR(cr))
                {
                    gmx_sumd(1, &gpb, cr);
                }

                if (debug)
                {
                    fprintf(debug, "CGE: EpotA %f EpotB %f EpotC %f gpb %f\n", s_a->epot, s_b->epot,
                            s_c->epot, gpb);
                }

                epot_repl = s_b->epot;

                /* Keep one of the intervals based on the value of the derivative at the new point */
                if (gpb > 0)
                {
                    /* Replace c endpoint with b */
                    swap_em_state(&s_b, &s_c);
                    c   = b;
                    gpc = gpb;
                }
                else
                {
                    /* Replace a endpoint with b */
                    swap_em_state(&s_b, &s_a);
                    a   = b;
                    gpa = gpb;
                }

                /*
                 * Stop search as soon as we find a value smaller than the endpoints.
                 * Never run more than 20 steps, no matter what.
                 */
                nminstep++;
            } while ((epot_repl > s_a->epot || epot_repl > s_c->epot) && (nminstep < 20));

            if (std::fabs(epot_repl - s_min->epot) < fabs(s_min->epot) * GMX_REAL_EPS || nminstep >= 20)
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
                    fprintf(debug, "CGE: C (%f) is lower than A (%f), moving C to B\n", s_c->epot,
                            s_a->epot);
                }
                swap_em_state(&s_b, &s_c);
                gpb = gpc;
            }
            else
            {
                if (debug)
                {
                    fprintf(debug, "CGE: A (%f) is lower than C (%f), moving A to B\n", s_a->epot,
                            s_c->epot);
                }
                swap_em_state(&s_b, &s_a);
                gpb = gpa;
            }
        }
        else
        {
            if (debug)
            {
                fprintf(debug, "CGE: Found a lower energy %f, moving C to B\n", s_c->epot);
            }
            swap_em_state(&s_b, &s_c);
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
        swap_em_state(&s_min, &s_b);
        gpa = gpb;

        /* Print it if necessary */
        if (MASTER(cr))
        {
            if (mdrunOptions.verbose)
            {
                double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
                fprintf(stderr, "\rStep %d, Epot=%12.6e, Fnorm=%9.3e, Fmax=%9.3e (atom %d)\n", step,
                        s_min->epot, s_min->fnorm / sqrtNumAtoms, s_min->fmax, s_min->a_fmax + 1);
                fflush(stderr);
            }
            /* Store the new (lower) energies */
            matrix nullBox = {};
            energyOutput.addDataAtEnergyStep(false, false, static_cast<double>(step), mdatoms->tmass,
                                             enerd, nullptr, nullptr, nullptr, nullBox, nullptr,
                                             nullptr, vir, pres, nullptr, mu_tot, constr);

            do_log = do_per_step(step, inputrec->nstlog);
            do_ene = do_per_step(step, inputrec->nstenergy);

            imdSession->fillEnergyRecord(step, TRUE);

            if (do_log)
            {
                energyOutput.printHeader(fplog, step, step);
            }
            energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), do_ene, FALSE, FALSE,
                                               do_log ? fplog : nullptr, step, step, fcd, nullptr);
        }

        /* Send energies and positions to the IMD client if bIMD is TRUE. */
        if (MASTER(cr) && imdSession->run(step, TRUE, state_global->box, state_global->x.rvec_array(), 0))
        {
            imdSession->sendPositionsAndEnergies();
        }

        /* Stop when the maximum force lies below tolerance.
         * If we have reached machine precision, converged is already set to true.
         */
        converged = converged || (s_min->fmax < inputrec->em_tol);

    } /* End of the loop */

    if (converged)
    {
        step--; /* we never took that last step in this case */
    }
    if (s_min->fmax > inputrec->em_tol)
    {
        if (MASTER(cr))
        {
            warn_step(fplog, inputrec->em_tol, s_min->fmax, step - 1 == number_steps, FALSE);
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
            energyOutput.printHeader(fplog, step, step);
        }
        if (!do_ene || !do_log)
        {
            /* Write final energy file entries */
            energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), !do_ene, FALSE, FALSE,
                                               !do_log ? fplog : nullptr, step, step, fcd, nullptr);
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
    /* Note that with 0 < nstfout != nstxout we can end up with two frames
     * in the trajectory with the same step number.
     */
    do_x = !do_per_step(step, inputrec->nstxout);
    do_f = (inputrec->nstfout > 0 && !do_per_step(step, inputrec->nstfout));

    write_em_traj(fplog, cr, outf, do_x, do_f, ftp2fn(efSTO, nfile, fnm), top_global, inputrec,
                  step, s_min, state_global, observablesHistory);


    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        print_converged(stderr, CG, inputrec->em_tol, step, converged, number_steps, s_min, sqrtNumAtoms);
        print_converged(fplog, CG, inputrec->em_tol, step, converged, number_steps, s_min, sqrtNumAtoms);

        fprintf(fplog, "\nPerformed %d energy evaluations in total.\n", neval);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    /* To print the actual number of steps we needed somewhere */
    walltime_accounting_set_nsteps_done(walltime_accounting, step);
}


void LegacySimulator::do_lbfgs()
{
    static const char* LBFGS = "Low-Memory BFGS Minimizer";
    em_state_t         ems;
    gmx_localtop_t     top;
    gmx_global_stat_t  gstat;
    t_graph*           graph;
    int                ncorr, nmaxcorr, point, cp, neval, nminstep;
    double             stepsize, step_taken, gpa, gpb, gpc, tmp, minstep;
    real *             rho, *alpha, *p, *s, **dx, **dg;
    real               a, b, c, maxdelta, delta;
    real               diag, Epot0;
    real               dgdx, dgdg, sq, yr, beta;
    gmx_bool           converged;
    rvec               mu_tot = { 0 };
    gmx_bool           do_log, do_ene, do_x, do_f, foundlower, *frozen;
    tensor             vir, pres;
    int                start, end, number_steps;
    int                i, k, m, n, gf, step;
    int                mdof_flags;
    auto               mdatoms = mdAtoms->mdatoms();

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Note that activating L-BFGS energy minimization via the "
                    "integrator .mdp option and the command gmx mdrun may "
                    "be available in a different form in a future version of GROMACS, "
                    "e.g. gmx minimize and an .mdp option.");

    if (PAR(cr))
    {
        gmx_fatal(FARGS, "L-BFGS minimization only supports a single rank");
    }

    if (nullptr != constr)
    {
        gmx_fatal(
                FARGS,
                "The combination of constraints and L-BFGS minimization is not implemented. Either "
                "do not use constraints, or use another minimizer (e.g. steepest descent).");
    }

    n        = 3 * state_global->natoms;
    nmaxcorr = inputrec->nbfgscorr;

    snew(frozen, n);

    snew(p, n);
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
    init_em(fplog, mdlog, LBFGS, cr, inputrec, imdSession, pull_work, state_global, top_global,
            &ems, &top, nrnb, fr, &graph, mdAtoms, &gstat, vsite, constr, nullptr);
    const bool        simulationsShareState = false;
    gmx_mdoutf*       outf = init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider,
                                   mdModulesNotifier, inputrec, top_global, nullptr, wcycle,
                                   StartingBehavior::NewSimulation, simulationsShareState, ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf), top_global, inputrec, pull_work, nullptr,
                                   false, StartingBehavior::NewSimulation, mdModulesNotifier);

    start = 0;
    end   = mdatoms->homenr;

    /* We need 4 working states */
    em_state_t  s0{}, s1{}, s2{}, s3{};
    em_state_t* sa   = &s0;
    em_state_t* sb   = &s1;
    em_state_t* sc   = &s2;
    em_state_t* last = &s3;
    /* Initialize by copying the state from ems (we could skip x and f here) */
    *sa = ems;
    *sb = ems;
    *sc = ems;

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
            frozen[3 * i + m] = (inputrec->opts.nFreeze[gf][m] != 0);
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
        construct_vsites(vsite, state_global->x.rvec_array(), 1, nullptr, top.idef.iparams,
                         top.idef.il, fr->ePBC, fr->bMolPBC, cr, state_global->box);
    }

    /* Call the force routine and some auxiliary (neighboursearching etc.) */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole
     */
    neval++;
    EnergyEvaluator energyEvaluator{
        fplog,      mdlog,     cr,      ms,     top_global,      &top,  inputrec,
        imdSession, pull_work, nrnb,    wcycle, gstat,           vsite, constr,
        fcd,        graph,     mdAtoms, fr,     runScheduleWork, enerd
    };
    energyEvaluator.run(&ems, mu_tot, vir, pres, -1, TRUE);

    if (MASTER(cr))
    {
        /* Copy stuff to the energy bin for easy printing etc. */
        matrix nullBox = {};
        energyOutput.addDataAtEnergyStep(false, false, static_cast<double>(step), mdatoms->tmass,
                                         enerd, nullptr, nullptr, nullptr, nullBox, nullptr,
                                         nullptr, vir, pres, nullptr, mu_tot, constr);

        energyOutput.printHeader(fplog, step, step);
        energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), TRUE, FALSE, FALSE, fplog, step,
                                           step, fcd, nullptr);
    }

    /* Set the initial step.
     * since it will be multiplied by the non-normalized search direction
     * vector (force vector the first time), we scale it by the
     * norm of the force.
     */

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        fprintf(stderr, "Using %d BFGS correction steps.\n\n", nmaxcorr);
        fprintf(stderr, "   F-max             = %12.5e on atom %d\n", ems.fmax, ems.a_fmax + 1);
        fprintf(stderr, "   F-Norm            = %12.5e\n", ems.fnorm / sqrtNumAtoms);
        fprintf(stderr, "\n");
        /* and copy to the log file too... */
        fprintf(fplog, "Using %d BFGS correction steps.\n\n", nmaxcorr);
        fprintf(fplog, "   F-max             = %12.5e on atom %d\n", ems.fmax, ems.a_fmax + 1);
        fprintf(fplog, "   F-Norm            = %12.5e\n", ems.fnorm / sqrtNumAtoms);
        fprintf(fplog, "\n");
    }

    // Point is an index to the memory of search directions, where 0 is the first one.
    point = 0;

    // Set initial search direction to the force (-gradient), or 0 for frozen particles.
    real* fInit = static_cast<real*>(ems.f.rvec_array()[0]);
    for (i = 0; i < n; i++)
    {
        if (!frozen[i])
        {
            dx[point][i] = fInit[i]; /* Initial search direction */
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
    stepsize = 1.0 / ems.fnorm;

    /* Start the loop over BFGS steps.
     * Each successful step is counted, and we continue until
     * we either converge or reach the max number of steps.
     */

    ncorr = 0;

    /* Set the gradient from the force */
    converged = FALSE;
    for (step = 0; (number_steps < 0 || step <= number_steps) && !converged; step++)
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

        mdoutf_write_to_trajectory_files(fplog, cr, outf, mdof_flags, top_global->natoms, step,
                                         static_cast<real>(step), &ems.s, state_global,
                                         observablesHistory, ems.f);

        /* Do the linesearching in the direction dx[point][0..(n-1)] */

        /* make s a pointer to current search direction - point=0 first time we get here */
        s = dx[point];

        real* xx = static_cast<real*>(ems.s.x.rvec_array()[0]);
        real* ff = static_cast<real*>(ems.f.rvec_array()[0]);

        // calculate line gradient in position A
        for (gpa = 0, i = 0; i < n; i++)
        {
            gpa -= s[i] * ff[i];
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
            tmp = s[i] / tmp;
            minstep += tmp * tmp;
        }
        minstep = GMX_REAL_EPS / sqrt(minstep / n);

        if (stepsize < minstep)
        {
            converged = TRUE;
            break;
        }

        // Before taking any steps along the line, store the old position
        *last       = ems;
        real* lastx = static_cast<real*>(last->s.x.data()[0]);
        real* lastf = static_cast<real*>(last->f.data()[0]);
        Epot0       = ems.epot;

        *sa = ems;

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
        a = 0.0;

        // Check stepsize first. We do not allow displacements
        // larger than emstep.
        //
        do
        {
            // Pick a new position C by adding stepsize to A.
            c = a + stepsize;

            // Calculate what the largest change in any individual coordinate
            // would be (translation along line * gradient along line)
            maxdelta = 0;
            for (i = 0; i < n; i++)
            {
                delta = c * s[i];
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
        } while (maxdelta > inputrec->em_stepsize);

        // Take a trial step and move the coordinate array xc[] to position C
        real* xc = static_cast<real*>(sc->s.x.rvec_array()[0]);
        for (i = 0; i < n; i++)
        {
            xc[i] = lastx[i] + c * s[i];
        }

        neval++;
        // Calculate energy for the trial step in position C
        energyEvaluator.run(sc, mu_tot, vir, pres, step, FALSE);

        // Calc line gradient in position C
        real* fc = static_cast<real*>(sc->f.rvec_array()[0]);
        for (gpc = 0, i = 0; i < n; i++)
        {
            gpc -= s[i] * fc[i]; /* f is negative gradient, thus the sign */
        }
        /* Sum the gradient along the line across CPUs */
        if (PAR(cr))
        {
            gmx_sumd(1, &gpc, cr);
        }

        // This is the max amount of increase in energy we tolerate.
        // By allowing VERY small changes (close to numerical precision) we
        // frequently find even better (lower) final energies.
        tmp = std::sqrt(GMX_REAL_EPS) * fabs(sa->epot);

        // Accept the step if the energy is lower in the new position C (compared to A),
        // or if it is not significantly higher and the line derivative is still negative.
        foundlower = sc->epot < sa->epot || (gpc < 0 && sc->epot < (sa->epot + tmp));
        // If true, great, we found a better energy. We no longer try to alter the
        // stepsize, but simply accept this new better position. The we select a new
        // search direction instead, which will be much more efficient than continuing
        // to take smaller steps along a line. Set fnorm based on the new C position,
        // which will be used to update the stepsize to 1/fnorm further down.

        // If false, the energy is NOT lower in point C, i.e. it will be the same
        // or higher than in point A. In this case it is pointless to move to point C,
        // so we will have to do more iterations along the same line to find a smaller
        // value in the interval [A=0.0,C].
        // Here, A is still 0.0, but that will change when we do a search in the interval
        // [0.0,C] below. That search we will do by interpolation or bisection rather
        // than with the stepsize, so no need to modify it. For the next search direction
        // it will be reset to 1/fnorm anyway.

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
            real fnorm = 0;
            nminstep   = 0;
            do
            {
                // Select a new trial point B in the interval [A,C].
                // If the derivatives at points a & c have different sign we interpolate to zero,
                // otherwise just do a bisection since there might be multiple minima/maxima
                // inside the interval.
                if (gpa < 0 && gpc > 0)
                {
                    b = a + gpa * (a - c) / (gpc - gpa);
                }
                else
                {
                    b = 0.5 * (a + c);
                }

                /* safeguard if interpolation close to machine accuracy causes errors:
                 * never go outside the interval
                 */
                if (b <= a || b >= c)
                {
                    b = 0.5 * (a + c);
                }

                // Take a trial step to point B
                real* xb = static_cast<real*>(sb->s.x.rvec_array()[0]);
                for (i = 0; i < n; i++)
                {
                    xb[i] = lastx[i] + b * s[i];
                }

                neval++;
                // Calculate energy for the trial step in point B
                energyEvaluator.run(sb, mu_tot, vir, pres, step, FALSE);
                fnorm = sb->fnorm;

                // Calculate gradient in point B
                real* fb = static_cast<real*>(sb->f.rvec_array()[0]);
                for (gpb = 0, i = 0; i < n; i++)
                {
                    gpb -= s[i] * fb[i]; /* f is negative gradient, thus the sign */
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
                    c = b;
                    /* copy state b to c */
                    *sc = *sb;
                }
                else
                {
                    /* Replace a endpoint with b */
                    a = b;
                    /* copy state b to a */
                    *sa = *sb;
                }

                /*
                 * Stop search as soon as we find a value smaller than the endpoints,
                 * or if the tolerance is below machine precision.
                 * Never run more than 20 steps, no matter what.
                 */
                nminstep++;
            } while ((sb->epot > sa->epot || sb->epot > sc->epot) && (nminstep < 20));

            if (std::fabs(sb->epot - Epot0) < GMX_REAL_EPS || nminstep >= 20)
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
                    stepsize = 1.0 / fnorm;
                    continue;
                }
            }

            /* Select min energy state of A & C, put the best in xx/ff/Epot
             */
            if (sc->epot < sa->epot)
            {
                /* Use state C */
                ems        = *sc;
                step_taken = c;
            }
            else
            {
                /* Use state A */
                ems        = *sa;
                step_taken = a;
            }
        }
        else
        {
            /* found lower */
            /* Use state C */
            ems        = *sc;
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
            dg[point][i] = lastf[i] - ff[i];
            dx[point][i] *= step_taken;
        }

        dgdg = 0;
        dgdx = 0;
        for (i = 0; i < n; i++)
        {
            dgdg += dg[point][i] * dg[point][i];
            dgdx += dg[point][i] * dx[point][i];
        }

        diag = dgdx / dgdg;

        rho[point] = 1.0 / dgdx;
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
                cp = ncorr - 1;
            }

            sq = 0;
            for (i = 0; i < n; i++)
            {
                sq += dx[cp][i] * p[i];
            }

            alpha[cp] = rho[cp] * sq;

            for (i = 0; i < n; i++)
            {
                p[i] -= alpha[cp] * dg[cp][i];
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
                yr += p[i] * dg[cp][i];
            }

            beta = rho[cp] * yr;
            beta = alpha[cp] - beta;

            for (i = 0; i < n; i++)
            {
                p[i] += beta * dx[cp][i];
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

        /* Print it if necessary */
        if (MASTER(cr))
        {
            if (mdrunOptions.verbose)
            {
                double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
                fprintf(stderr, "\rStep %d, Epot=%12.6e, Fnorm=%9.3e, Fmax=%9.3e (atom %d)\n", step,
                        ems.epot, ems.fnorm / sqrtNumAtoms, ems.fmax, ems.a_fmax + 1);
                fflush(stderr);
            }
            /* Store the new (lower) energies */
            matrix nullBox = {};
            energyOutput.addDataAtEnergyStep(false, false, static_cast<double>(step), mdatoms->tmass,
                                             enerd, nullptr, nullptr, nullptr, nullBox, nullptr,
                                             nullptr, vir, pres, nullptr, mu_tot, constr);

            do_log = do_per_step(step, inputrec->nstlog);
            do_ene = do_per_step(step, inputrec->nstenergy);

            imdSession->fillEnergyRecord(step, TRUE);

            if (do_log)
            {
                energyOutput.printHeader(fplog, step, step);
            }
            energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), do_ene, FALSE, FALSE,
                                               do_log ? fplog : nullptr, step, step, fcd, nullptr);
        }

        /* Send x and E to IMD client, if bIMD is TRUE. */
        if (imdSession->run(step, TRUE, state_global->box, state_global->x.rvec_array(), 0) && MASTER(cr))
        {
            imdSession->sendPositionsAndEnergies();
        }

        // Reset stepsize in we are doing more iterations
        stepsize = 1.0;

        /* Stop when the maximum force lies below tolerance.
         * If we have reached machine precision, converged is already set to true.
         */
        converged = converged || (ems.fmax < inputrec->em_tol);

    } /* End of the loop */

    if (converged)
    {
        step--; /* we never took that last step in this case */
    }
    if (ems.fmax > inputrec->em_tol)
    {
        if (MASTER(cr))
        {
            warn_step(fplog, inputrec->em_tol, ems.fmax, step - 1 == number_steps, FALSE);
        }
        converged = FALSE;
    }

    /* If we printed energy and/or logfile last step (which was the last step)
     * we don't have to do it again, but otherwise print the final values.
     */
    if (!do_log) /* Write final value to log since we didn't do anythin last step */
    {
        energyOutput.printHeader(fplog, step, step);
    }
    if (!do_ene || !do_log) /* Write final energy file entries */
    {
        energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), !do_ene, FALSE, FALSE,
                                           !do_log ? fplog : nullptr, step, step, fcd, nullptr);
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
    write_em_traj(fplog, cr, outf, do_x, do_f, ftp2fn(efSTO, nfile, fnm), top_global, inputrec,
                  step, &ems, state_global, observablesHistory);

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));
        print_converged(stderr, LBFGS, inputrec->em_tol, step, converged, number_steps, &ems, sqrtNumAtoms);
        print_converged(fplog, LBFGS, inputrec->em_tol, step, converged, number_steps, &ems, sqrtNumAtoms);

        fprintf(fplog, "\nPerformed %d energy evaluations in total.\n", neval);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    /* To print the actual number of steps we needed somewhere */
    walltime_accounting_set_nsteps_done(walltime_accounting, step);
}

void LegacySimulator::do_steep()
{
    const char*       SD = "Steepest Descents";
    gmx_localtop_t    top;
    gmx_global_stat_t gstat;
    t_graph*          graph;
    real              stepsize;
    real              ustep;
    gmx_bool          bDone, bAbort, do_x, do_f;
    tensor            vir, pres;
    rvec              mu_tot = { 0 };
    int               nsteps;
    int               count          = 0;
    int               steps_accepted = 0;
    auto              mdatoms        = mdAtoms->mdatoms();

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Note that activating steepest-descent energy minimization via the "
                    "integrator .mdp option and the command gmx mdrun may "
                    "be available in a different form in a future version of GROMACS, "
                    "e.g. gmx minimize and an .mdp option.");

    /* Create 2 states on the stack and extract pointers that we will swap */
    em_state_t  s0{}, s1{};
    em_state_t* s_min = &s0;
    em_state_t* s_try = &s1;

    /* Init em and store the local state in s_try */
    init_em(fplog, mdlog, SD, cr, inputrec, imdSession, pull_work, state_global, top_global, s_try,
            &top, nrnb, fr, &graph, mdAtoms, &gstat, vsite, constr, nullptr);
    const bool        simulationsShareState = false;
    gmx_mdoutf*       outf = init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider,
                                   mdModulesNotifier, inputrec, top_global, nullptr, wcycle,
                                   StartingBehavior::NewSimulation, simulationsShareState, ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf), top_global, inputrec, pull_work, nullptr,
                                   false, StartingBehavior::NewSimulation, mdModulesNotifier);

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
    EnergyEvaluator energyEvaluator{
        fplog,      mdlog,     cr,      ms,     top_global,      &top,  inputrec,
        imdSession, pull_work, nrnb,    wcycle, gstat,           vsite, constr,
        fcd,        graph,     mdAtoms, fr,     runScheduleWork, enerd
    };

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
        bool validStep = true;
        if (count > 0)
        {
            validStep = do_em_step(cr, inputrec, mdatoms, s_min, stepsize, &s_min->f, s_try, constr, count);
        }

        if (validStep)
        {
            energyEvaluator.run(s_try, mu_tot, vir, pres, count, count == 0);
        }
        else
        {
            // Signal constraint error during stepping with energy=inf
            s_try->epot = std::numeric_limits<real>::infinity();
        }

        if (MASTER(cr))
        {
            energyOutput.printHeader(fplog, count, count);
        }

        if (count == 0)
        {
            s_min->epot = s_try->epot;
        }

        /* Print it if necessary  */
        if (MASTER(cr))
        {
            if (mdrunOptions.verbose)
            {
                fprintf(stderr, "Step=%5d, Dmax= %6.1e nm, Epot= %12.5e Fmax= %11.5e, atom= %d%c",
                        count, ustep, s_try->epot, s_try->fmax, s_try->a_fmax + 1,
                        ((count == 0) || (s_try->epot < s_min->epot)) ? '\n' : '\r');
                fflush(stderr);
            }

            if ((count == 0) || (s_try->epot < s_min->epot))
            {
                /* Store the new (lower) energies  */
                matrix nullBox = {};
                energyOutput.addDataAtEnergyStep(false, false, static_cast<double>(count), mdatoms->tmass,
                                                 enerd, nullptr, nullptr, nullptr, nullBox, nullptr,
                                                 nullptr, vir, pres, nullptr, mu_tot, constr);

                imdSession->fillEnergyRecord(count, TRUE);

                const bool do_dr = do_per_step(steps_accepted, inputrec->nstdisreout);
                const bool do_or = do_per_step(steps_accepted, inputrec->nstorireout);
                energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), TRUE, do_dr, do_or,
                                                   fplog, count, count, fcd, nullptr);
                fflush(fplog);
            }
        }

        /* Now if the new energy is smaller than the previous...
         * or if this is the first step!
         * or if we did random steps!
         */

        if ((count == 0) || (s_try->epot < s_min->epot))
        {
            steps_accepted++;

            /* Test whether the convergence criterion is met...  */
            bDone = (s_try->fmax < inputrec->em_tol);

            /* Copy the arrays for force, positions and energy  */
            /* The 'Min' array always holds the coords and forces of the minimal
               sampled energy  */
            swap_em_state(&s_min, &s_try);
            if (count > 0)
            {
                ustep *= 1.2;
            }

            /* Write to trn, if necessary */
            do_x = do_per_step(steps_accepted, inputrec->nstxout);
            do_f = do_per_step(steps_accepted, inputrec->nstfout);
            write_em_traj(fplog, cr, outf, do_x, do_f, nullptr, top_global, inputrec, count, s_min,
                          state_global, observablesHistory);
        }
        else
        {
            /* If energy is not smaller make the step smaller...  */
            ustep *= 0.5;

            if (DOMAINDECOMP(cr) && s_min->s.ddp_count != cr->dd->ddp_count)
            {
                /* Reload the old state */
                em_dd_partition_system(fplog, mdlog, count, cr, top_global, inputrec, imdSession,
                                       pull_work, s_min, &top, mdAtoms, fr, vsite, constr, nrnb, wcycle);
            }
        }

        // If the force is very small after finishing minimization,
        // we risk dividing by zero when calculating the step size.
        // So we check first if the minimization has stopped before
        // trying to obtain a new step size.
        if (!bDone)
        {
            /* Determine new step  */
            stepsize = ustep / s_min->fmax;
        }

        /* Check if stepsize is too small, with 1 nm as a characteristic length */
#if GMX_DOUBLE
        if (count == nsteps || ustep < 1e-12)
#else
        if (count == nsteps || ustep < 1e-6)
#endif
        {
            if (MASTER(cr))
            {
                warn_step(fplog, inputrec->em_tol, s_min->fmax, count == nsteps, constr != nullptr);
            }
            bAbort = TRUE;
        }

        /* Send IMD energies and positions, if bIMD is TRUE. */
        if (imdSession->run(count, TRUE, state_global->box,
                            MASTER(cr) ? state_global->x.rvec_array() : nullptr, 0)
            && MASTER(cr))
        {
            imdSession->sendPositionsAndEnergies();
        }

        count++;
    } /* End of the loop  */

    /* Print some data...  */
    if (MASTER(cr))
    {
        fprintf(stderr, "\nwriting lowest energy coordinates.\n");
    }
    write_em_traj(fplog, cr, outf, TRUE, inputrec->nstfout != 0, ftp2fn(efSTO, nfile, fnm),
                  top_global, inputrec, count, s_min, state_global, observablesHistory);

    if (MASTER(cr))
    {
        double sqrtNumAtoms = sqrt(static_cast<double>(state_global->natoms));

        print_converged(stderr, SD, inputrec->em_tol, count, bDone, nsteps, s_min, sqrtNumAtoms);
        print_converged(fplog, SD, inputrec->em_tol, count, bDone, nsteps, s_min, sqrtNumAtoms);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    /* To print the actual number of steps we needed somewhere */
    inputrec->nsteps = count;

    walltime_accounting_set_nsteps_done(walltime_accounting, count);
}

void LegacySimulator::do_nm()
{
    const char*         NM = "Normal Mode Analysis";
    int                 nnodes;
    gmx_localtop_t      top;
    gmx_global_stat_t   gstat;
    t_graph*            graph;
    tensor              vir, pres;
    rvec                mu_tot = { 0 };
    rvec*               dfdx;
    gmx_bool            bSparse; /* use sparse matrix storage format */
    size_t              sz;
    gmx_sparsematrix_t* sparse_matrix = nullptr;
    real*               full_matrix   = nullptr;

    /* added with respect to mdrun */
    int  row, col;
    real der_range = 10.0 * std::sqrt(GMX_REAL_EPS);
    real x_min;
    bool bIsMaster = MASTER(cr);
    auto mdatoms   = mdAtoms->mdatoms();

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(
                    "Note that activating normal-mode analysis via the integrator "
                    ".mdp option and the command gmx mdrun may "
                    "be available in a different form in a future version of GROMACS, "
                    "e.g. gmx normal-modes.");

    if (constr != nullptr)
    {
        gmx_fatal(
                FARGS,
                "Constraints present with Normal Mode Analysis, this combination is not supported");
    }

    gmx_shellfc_t* shellfc;

    em_state_t state_work{};

    /* Init em and store the local state in state_minimum */
    init_em(fplog, mdlog, NM, cr, inputrec, imdSession, pull_work, state_global, top_global,
            &state_work, &top, nrnb, fr, &graph, mdAtoms, &gstat, vsite, constr, &shellfc);
    const bool  simulationsShareState = false;
    gmx_mdoutf* outf = init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider,
                                   mdModulesNotifier, inputrec, top_global, nullptr, wcycle,
                                   StartingBehavior::NewSimulation, simulationsShareState, ms);

    std::vector<int>       atom_index = get_atom_index(top_global);
    std::vector<gmx::RVec> fneg(atom_index.size(), { 0, 0, 0 });
    snew(dfdx, atom_index.size());

#if !GMX_DOUBLE
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
    if (EEL_FULL(fr->ic->eeltype) || fr->rlist == 0.0)
    {
        GMX_LOG(mdlog.warning)
                .appendText("Non-cutoff electrostatics used, forcing full Hessian format.");
        bSparse = FALSE;
    }
    else if (atom_index.size() < 1000)
    {
        GMX_LOG(mdlog.warning)
                .appendTextFormatted("Small system size (N=%zu), using full Hessian format.",
                                     atom_index.size());
        bSparse = FALSE;
    }
    else
    {
        GMX_LOG(mdlog.warning).appendText("Using compressed symmetric sparse Hessian format.");
        bSparse = TRUE;
    }

    /* Number of dimensions, based on real atoms, that is not vsites or shell */
    sz = DIM * atom_index.size();

    fprintf(stderr, "Allocating Hessian memory...\n\n");

    if (bSparse)
    {
        sparse_matrix                       = gmx_sparsematrix_init(sz);
        sparse_matrix->compressed_symmetric = TRUE;
    }
    else
    {
        snew(full_matrix, sz * sz);
    }

    /* Write start time and temperature */
    print_em_start(fplog, cr, walltime_accounting, wcycle, NM);

    /* fudge nr of steps to nr of atoms */
    inputrec->nsteps = atom_index.size() * 2;

    if (bIsMaster)
    {
        fprintf(stderr, "starting normal mode calculation '%s'\n%" PRId64 " steps.\n\n",
                *(top_global->name), inputrec->nsteps);
    }

    nnodes = cr->nnodes;

    /* Make evaluate_energy do a single node force calculation */
    cr->nnodes = 1;
    EnergyEvaluator energyEvaluator{
        fplog,      mdlog,     cr,      ms,     top_global,      &top,  inputrec,
        imdSession, pull_work, nrnb,    wcycle, gstat,           vsite, constr,
        fcd,        graph,     mdAtoms, fr,     runScheduleWork, enerd
    };
    energyEvaluator.run(&state_work, mu_tot, vir, pres, -1, TRUE);
    cr->nnodes = nnodes;

    /* if forces are not small, warn user */
    get_state_f_norm_max(cr, &(inputrec->opts), mdatoms, &state_work);

    GMX_LOG(mdlog.warning).appendTextFormatted("Maximum force:%12.5e", state_work.fmax);
    if (state_work.fmax > 1.0e-3)
    {
        GMX_LOG(mdlog.warning)
                .appendText(
                        "The force is probably not small enough to "
                        "ensure that you are at a minimum.\n"
                        "Be aware that negative eigenvalues may occur\n"
                        "when the resulting matrix is diagonalized.");
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
    bool bNS          = true;
    auto state_work_x = makeArrayRef(state_work.s.x);
    auto state_work_f = makeArrayRef(state_work.f);
    for (index aid = cr->nodeid; aid < ssize(atom_index); aid += nnodes)
    {
        size_t atom = atom_index[aid];
        for (size_t d = 0; d < DIM; d++)
        {
            int64_t step        = 0;
            int     force_flags = GMX_FORCE_STATECHANGED | GMX_FORCE_ALLFORCES;
            double  t           = 0;

            x_min = state_work_x[atom][d];

            for (unsigned int dx = 0; (dx < 2); dx++)
            {
                if (dx == 0)
                {
                    state_work_x[atom][d] = x_min - der_range;
                }
                else
                {
                    state_work_x[atom][d] = x_min + der_range;
                }

                /* Make evaluate_energy do a single node force calculation */
                cr->nnodes = 1;
                if (shellfc)
                {
                    /* Now is the time to relax the shells */
                    relax_shell_flexcon(fplog, cr, ms, mdrunOptions.verbose, nullptr, step, inputrec,
                                        imdSession, pull_work, bNS, force_flags, &top, constr, enerd,
                                        fcd, state_work.s.natoms, state_work.s.x.arrayRefWithPadding(),
                                        state_work.s.v.arrayRefWithPadding(), state_work.s.box,
                                        state_work.s.lambda, &state_work.s.hist,
                                        state_work.f.arrayRefWithPadding(), vir, mdatoms, nrnb,
                                        wcycle, graph, shellfc, fr, runScheduleWork, t, mu_tot,
                                        vsite, DDBalanceRegionHandler(nullptr));
                    bNS = false;
                    step++;
                }
                else
                {
                    energyEvaluator.run(&state_work, mu_tot, vir, pres, aid * 2 + dx, FALSE);
                }

                cr->nnodes = nnodes;

                if (dx == 0)
                {
                    std::copy(state_work_f.begin(), state_work_f.begin() + atom_index.size(),
                              fneg.begin());
                }
            }

            /* x is restored to original */
            state_work_x[atom][d] = x_min;

            for (size_t j = 0; j < atom_index.size(); j++)
            {
                for (size_t k = 0; (k < DIM); k++)
                {
                    dfdx[j][k] = -(state_work_f[atom_index[j]][k] - fneg[j][k]) / (2 * der_range);
                }
            }

            if (!bIsMaster)
            {
#if GMX_MPI
#    define mpi_type GMX_MPI_REAL
                MPI_Send(dfdx[0], atom_index.size() * DIM, mpi_type, MASTER(cr), cr->nodeid,
                         cr->mpi_comm_mygroup);
#endif
            }
            else
            {
                for (index node = 0; (node < nnodes && aid + node < ssize(atom_index)); node++)
                {
                    if (node > 0)
                    {
#if GMX_MPI
                        MPI_Status stat;
                        MPI_Recv(dfdx[0], atom_index.size() * DIM, mpi_type, node, node,
                                 cr->mpi_comm_mygroup, &stat);
#    undef mpi_type
#endif
                    }

                    row = (aid + node) * DIM + d;

                    for (size_t j = 0; j < atom_index.size(); j++)
                    {
                        for (size_t k = 0; k < DIM; k++)
                        {
                            col = j * DIM + k;

                            if (bSparse)
                            {
                                if (col >= row && dfdx[j][k] != 0.0)
                                {
                                    gmx_sparsematrix_increment_value(sparse_matrix, row, col, dfdx[j][k]);
                                }
                            }
                            else
                            {
                                full_matrix[row * sz + col] = dfdx[j][k];
                            }
                        }
                    }
                }
            }

            if (mdrunOptions.verbose && fplog)
            {
                fflush(fplog);
            }
        }
        /* write progress */
        if (bIsMaster && mdrunOptions.verbose)
        {
            fprintf(stderr, "\rFinished step %d out of %td",
                    std::min<int>(atom + nnodes, atom_index.size()), ssize(atom_index));
            fflush(stderr);
        }
    }

    if (bIsMaster)
    {
        fprintf(stderr, "\n\nWriting Hessian...\n");
        gmx_mtxio_write(ftp2fn(efMTX, nfile, fnm), sz, sz, full_matrix, sparse_matrix);
    }

    finish_em(cr, outf, walltime_accounting, wcycle);

    walltime_accounting_set_nsteps_done(walltime_accounting, atom_index.size() * 2);
}

} // namespace gmx
