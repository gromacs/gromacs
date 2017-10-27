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

#include "sim_util.h"

#include "config.h"

#include <stdio.h>
#include <string.h>

#include <cmath>
#include <cstdint>

#include <array>

#include "gromacs/awh/awh.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/chargegroup.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/listed-forces/disre.h"
#include "gromacs/listed-forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/calcmu.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/genborn.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_grid.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_gpu_ref.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/wallcyclereporting.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"

#include "nbnxn_gpu.h"
#include "nbnxn_kernels/nbnxn_kernel_cpu.h"
#include "nbnxn_kernels/nbnxn_kernel_prune.h"

// TODO: this environment variable allows us to verify before release
// that on less common architectures the total cost of polling is not larger than
// a blocking wait (so polling does not introduce overhead when the static
// PME-first ordering would suffice).
static const bool c_disableAlternatingWait = (getenv("GMX_DISABLE_ALTERNATING_GPU_WAIT") != nullptr);


void print_time(FILE                     *out,
                gmx_walltime_accounting_t walltime_accounting,
                gmx_int64_t               step,
                t_inputrec               *ir,
                t_commrec gmx_unused     *cr)
{
    time_t finish;
    char   timebuf[STRLEN];
    double dt, elapsed_seconds, time_per_step;
    char   buf[48];

#if !GMX_THREAD_MPI
    if (!PAR(cr))
#endif
    {
        fprintf(out, "\r");
    }
    fprintf(out, "step %s", gmx_step_str(step, buf));
    fflush(out);

    if ((step >= ir->nstlist))
    {
        double seconds_since_epoch = gmx_gettime();
        elapsed_seconds = seconds_since_epoch - walltime_accounting_get_start_time_stamp(walltime_accounting);
        time_per_step   = elapsed_seconds/(step - ir->init_step + 1);
        dt              = (ir->nsteps + ir->init_step - step) * time_per_step;

        if (ir->nsteps >= 0)
        {
            if (dt >= 300)
            {
                finish = (time_t) (seconds_since_epoch + dt);
                gmx_ctime_r(&finish, timebuf, STRLEN);
                sprintf(buf, "%s", timebuf);
                buf[strlen(buf)-1] = '\0';
                fprintf(out, ", will finish %s", buf);
            }
            else
            {
                fprintf(out, ", remaining wall clock time: %5d s          ", (int)dt);
            }
        }
        else
        {
            fprintf(out, " performance: %.1f ns/day    ",
                    ir->delta_t/1000*24*60*60/time_per_step);
        }
    }
#if !GMX_THREAD_MPI
    if (PAR(cr))
    {
        fprintf(out, "\n");
    }
#endif

    fflush(out);
}

void print_date_and_time(FILE *fplog, int nodeid, const char *title,
                         double the_time)
{
    char   time_string[STRLEN];

    if (!fplog)
    {
        return;
    }

    {
        int    i;
        char   timebuf[STRLEN];
        time_t temp_time = (time_t) the_time;

        gmx_ctime_r(&temp_time, timebuf, STRLEN);
        for (i = 0; timebuf[i] >= ' '; i++)
        {
            time_string[i] = timebuf[i];
        }
        time_string[i] = '\0';
    }

    fprintf(fplog, "%s on rank %d %s\n", title, nodeid, time_string);
}

void print_start(FILE *fplog, t_commrec *cr,
                 gmx_walltime_accounting_t walltime_accounting,
                 const char *name)
{
    char buf[STRLEN];

    sprintf(buf, "Started %s", name);
    print_date_and_time(fplog, cr->nodeid, buf,
                        walltime_accounting_get_start_time_stamp(walltime_accounting));
}

static void sum_forces(rvec f[], gmx::ArrayRef<const gmx::RVec> forceToAdd)
{
    const int end = forceToAdd.size();

    // cppcheck-suppress unreadVariable
    int gmx_unused nt = gmx_omp_nthreads_get(emntDefault);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (int i = 0; i < end; i++)
    {
        rvec_inc(f[i], forceToAdd[i]);
    }
}

static void pme_gpu_reduce_outputs(gmx_wallcycle_t            wcycle,
                                   ForceWithVirial           *forceWithVirial,
                                   ArrayRef<const gmx::RVec>  pmeForces,
                                   gmx_enerdata_t            *enerd,
                                   const tensor               vir_Q,
                                   real                       Vlr_q)
{
    wallcycle_start(wcycle, ewcPME_GPU_F_REDUCTION);
    GMX_ASSERT(forceWithVirial, "Invalid force pointer");
    forceWithVirial->addVirialContribution(vir_Q);
    enerd->term[F_COUL_RECIP] += Vlr_q;
    sum_forces(as_rvec_array(forceWithVirial->force_.data()), pmeForces);
    wallcycle_stop(wcycle, ewcPME_GPU_F_REDUCTION);
}

static void calc_virial(int start, int homenr, rvec x[], rvec f[],
                        tensor vir_part, t_graph *graph, matrix box,
                        t_nrnb *nrnb, const t_forcerec *fr, int ePBC)
{
    int    i;

    /* The short-range virial from surrounding boxes */
    calc_vir(SHIFTS, fr->shift_vec, fr->fshift, vir_part, ePBC == epbcSCREW, box);
    inc_nrnb(nrnb, eNR_VIRIAL, SHIFTS);

    /* Calculate partial virial, for local atoms only, based on short range.
     * Total virial is computed in global_stat, called from do_md
     */
    f_calc_vir(start, start+homenr, x, f, vir_part, graph, box);
    inc_nrnb(nrnb, eNR_VIRIAL, homenr);

    /* Add wall contribution */
    for (i = 0; i < DIM; i++)
    {
        vir_part[i][ZZ] += fr->vir_wall_z[i];
    }

    if (debug)
    {
        pr_rvecs(debug, 0, "vir_part", vir_part, DIM);
    }
}

static void pull_potential_wrapper(t_commrec *cr,
                                   t_inputrec *ir,
                                   matrix box, rvec x[],
                                   ForceWithVirial *force,
                                   t_mdatoms *mdatoms,
                                   gmx_enerdata_t *enerd,
                                   real *lambda,
                                   double t,
                                   gmx_wallcycle_t wcycle)
{
    t_pbc  pbc;
    real   dvdl;

    /* Calculate the center of mass forces, this requires communication,
     * which is why pull_potential is called close to other communication.
     */
    wallcycle_start(wcycle, ewcPULLPOT);
    set_pbc(&pbc, ir->ePBC, box);
    dvdl                     = 0;
    enerd->term[F_COM_PULL] +=
        pull_potential(ir->pull_work, mdatoms, &pbc,
                       cr, t, lambda[efptRESTRAINT], x, force, &dvdl);
    enerd->dvdl_lin[efptRESTRAINT] += dvdl;
    wallcycle_stop(wcycle, ewcPULLPOT);
}

static void pme_receive_force_ener(t_commrec       *cr,
                                   ForceWithVirial *forceWithVirial,
                                   gmx_enerdata_t  *enerd,
                                   gmx_wallcycle_t  wcycle)
{
    real   e_q, e_lj, dvdl_q, dvdl_lj;
    float  cycles_ppdpme, cycles_seppme;

    cycles_ppdpme = wallcycle_stop(wcycle, ewcPPDURINGPME);
    dd_cycles_add(cr->dd, cycles_ppdpme, ddCyclPPduringPME);

    /* In case of node-splitting, the PP nodes receive the long-range
     * forces, virial and energy from the PME nodes here.
     */
    wallcycle_start(wcycle, ewcPP_PMEWAITRECVF);
    dvdl_q  = 0;
    dvdl_lj = 0;
    gmx_pme_receive_f(cr, forceWithVirial, &e_q, &e_lj, &dvdl_q, &dvdl_lj,
                      &cycles_seppme);
    enerd->term[F_COUL_RECIP] += e_q;
    enerd->term[F_LJ_RECIP]   += e_lj;
    enerd->dvdl_lin[efptCOUL] += dvdl_q;
    enerd->dvdl_lin[efptVDW]  += dvdl_lj;

    if (wcycle)
    {
        dd_cycles_add(cr->dd, cycles_seppme, ddCyclPME);
    }
    wallcycle_stop(wcycle, ewcPP_PMEWAITRECVF);
}

static void print_large_forces(FILE *fp, t_mdatoms *md, t_commrec *cr,
                               gmx_int64_t step, real forceTolerance,
                               const rvec *x, const rvec *f)
{
    real           force2Tolerance = gmx::square(forceTolerance);
    std::uintmax_t numNonFinite    = 0;
    for (int i = 0; i < md->homenr; i++)
    {
        real force2    = norm2(f[i]);
        bool nonFinite = !std::isfinite(force2);
        if (force2 >= force2Tolerance || nonFinite)
        {
            fprintf(fp, "step %" GMX_PRId64 " atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
                    step,
                    ddglatnr(cr->dd, i), x[i][XX], x[i][YY], x[i][ZZ], std::sqrt(force2));
        }
        if (nonFinite)
        {
            numNonFinite++;
        }
    }
    if (numNonFinite > 0)
    {
        /* Note that with MPI this fatal call on one rank might interrupt
         * the printing on other ranks. But we can only avoid that with
         * an expensive MPI barrier that we would need at each step.
         */
        gmx_fatal(FARGS, "At step %" GMX_PRId64 " detected non-finite forces on %ju atoms", step, numNonFinite);
    }
}

static void post_process_forces(t_commrec *cr,
                                gmx_int64_t step,
                                t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                                gmx_localtop_t *top,
                                matrix box, rvec x[],
                                rvec f[],
                                ForceWithVirial *forceWithVirial,
                                tensor vir_force,
                                t_mdatoms *mdatoms,
                                t_graph *graph,
                                t_forcerec *fr, gmx_vsite_t *vsite,
                                int flags)
{
    if (fr->haveDirectVirialContributions)
    {
        rvec *fDirectVir = as_rvec_array(forceWithVirial->force_.data());

        if (vsite)
        {
            /* Spread the mesh force on virtual sites to the other particles...
             * This is parallellized. MPI communication is performed
             * if the constructing atoms aren't local.
             */
            wallcycle_start(wcycle, ewcVSITESPREAD);
            matrix virial = { { 0 } };
            spread_vsite_f(vsite, x, fDirectVir, nullptr,
                           (flags & GMX_FORCE_VIRIAL), virial,
                           nrnb,
                           &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
            forceWithVirial->addVirialContribution(virial);
            wallcycle_stop(wcycle, ewcVSITESPREAD);
        }

        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Now add the forces, this is local */
            sum_forces(f, forceWithVirial->force_);

            /* Add the direct virial contributions */
            GMX_ASSERT(forceWithVirial->computeVirial_, "forceWithVirial should request virial computation when we request the virial");
            m_add(vir_force, forceWithVirial->getVirial(), vir_force);

            if (debug)
            {
                pr_rvecs(debug, 0, "vir_force", vir_force, DIM);
            }
        }
    }

    if (fr->print_force >= 0)
    {
        print_large_forces(stderr, mdatoms, cr, step, fr->print_force, x, f);
    }
}

static void do_nb_verlet(t_forcerec *fr,
                         interaction_const_t *ic,
                         gmx_enerdata_t *enerd,
                         int flags, int ilocality,
                         int clearF,
                         gmx_int64_t step,
                         t_nrnb *nrnb,
                         gmx_wallcycle_t wcycle)
{
    if (!(flags & GMX_FORCE_NONBONDED))
    {
        /* skip non-bonded calculation */
        return;
    }

    nonbonded_verlet_t       *nbv  = fr->nbv;
    nonbonded_verlet_group_t *nbvg = &nbv->grp[ilocality];

    /* GPU kernel launch overhead is already timed separately */
    if (fr->cutoff_scheme != ecutsVERLET)
    {
        gmx_incons("Invalid cut-off scheme passed!");
    }

    bool bUsingGpuKernels = (nbvg->kernel_type == nbnxnk8x8x8_GPU);

    if (!bUsingGpuKernels)
    {
        /* When dynamic pair-list  pruning is requested, we need to prune
         * at nstlistPrune steps.
         */
        if (nbv->listParams->useDynamicPruning &&
            (step - nbvg->nbl_lists.outerListCreationStep) % nbv->listParams->nstlistPrune == 0)
        {
            /* Prune the pair-list beyond fr->ic->rlistPrune using
             * the current coordinates of the atoms.
             */
            wallcycle_sub_start(wcycle, ewcsNONBONDED_PRUNING);
            nbnxn_kernel_cpu_prune(nbvg, nbv->nbat, fr->shift_vec, nbv->listParams->rlistInner);
            wallcycle_sub_stop(wcycle, ewcsNONBONDED_PRUNING);
        }

        wallcycle_sub_start(wcycle, ewcsNONBONDED);
    }

    switch (nbvg->kernel_type)
    {
        case nbnxnk4x4_PlainC:
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
            nbnxn_kernel_cpu(nbvg,
                             nbv->nbat,
                             ic,
                             fr->shift_vec,
                             flags,
                             clearF,
                             fr->fshift[0],
                             enerd->grpp.ener[egCOULSR],
                             fr->bBHAM ?
                             enerd->grpp.ener[egBHAMSR] :
                             enerd->grpp.ener[egLJSR]
#ifdef BUILD_WITH_FDA
                             ,
                             fr->fda,
                             fr->nbv->nbs->a
#endif
                            );
            break;

        case nbnxnk8x8x8_GPU:
            nbnxn_gpu_launch_kernel(nbv->gpu_nbv, nbv->nbat, flags, ilocality);
            break;

        case nbnxnk8x8x8_PlainC:
            nbnxn_kernel_gpu_ref(nbvg->nbl_lists.nbl[0],
                                 nbv->nbat, ic,
                                 fr->shift_vec,
                                 flags,
                                 clearF,
                                 nbv->nbat->out[0].f,
                                 fr->fshift[0],
                                 enerd->grpp.ener[egCOULSR],
                                 fr->bBHAM ?
                                 enerd->grpp.ener[egBHAMSR] :
                                 enerd->grpp.ener[egLJSR]
#ifdef BUILD_WITH_FDA
                                 ,
                                 fr->fda,
                                 fr->nbv->nbs->a
#endif
                                );
            break;

        default:
            GMX_RELEASE_ASSERT(false, "Invalid nonbonded kernel type passed!");

    }
    if (!bUsingGpuKernels)
    {
        wallcycle_sub_stop(wcycle, ewcsNONBONDED);
    }

    int enr_nbnxn_kernel_ljc, enr_nbnxn_kernel_lj;
    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_RF;
    }
    else if ((!bUsingGpuKernels && nbvg->ewald_excl == ewaldexclAnalytical) ||
             (bUsingGpuKernels && nbnxn_gpu_is_kernel_ewald_analytical(nbv->gpu_nbv)))
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_EWALD;
    }
    else
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_TAB;
    }
    enr_nbnxn_kernel_lj = eNR_NBNXN_LJ;
    if (flags & GMX_FORCE_ENERGY)
    {
        /* In eNR_??? the nbnxn F+E kernels are always the F kernel + 1 */
        enr_nbnxn_kernel_ljc += 1;
        enr_nbnxn_kernel_lj  += 1;
    }

    inc_nrnb(nrnb, enr_nbnxn_kernel_ljc,
             nbvg->nbl_lists.natpair_ljq);
    inc_nrnb(nrnb, enr_nbnxn_kernel_lj,
             nbvg->nbl_lists.natpair_lj);
    /* The Coulomb-only kernels are offset -eNR_NBNXN_LJ_RF+eNR_NBNXN_RF */
    inc_nrnb(nrnb, enr_nbnxn_kernel_ljc-eNR_NBNXN_LJ_RF+eNR_NBNXN_RF,
             nbvg->nbl_lists.natpair_q);

    if (ic->vdw_modifier == eintmodFORCESWITCH)
    {
        /* We add up the switch cost separately */
        inc_nrnb(nrnb, eNR_NBNXN_ADD_LJ_FSW+((flags & GMX_FORCE_ENERGY) ? 1 : 0),
                 nbvg->nbl_lists.natpair_ljq + nbvg->nbl_lists.natpair_lj);
    }
    if (ic->vdw_modifier == eintmodPOTSWITCH)
    {
        /* We add up the switch cost separately */
        inc_nrnb(nrnb, eNR_NBNXN_ADD_LJ_PSW+((flags & GMX_FORCE_ENERGY) ? 1 : 0),
                 nbvg->nbl_lists.natpair_ljq + nbvg->nbl_lists.natpair_lj);
    }
    if (ic->vdwtype == evdwPME)
    {
        /* We add up the LJ Ewald cost separately */
        inc_nrnb(nrnb, eNR_NBNXN_ADD_LJ_EWALD+((flags & GMX_FORCE_ENERGY) ? 1 : 0),
                 nbvg->nbl_lists.natpair_ljq + nbvg->nbl_lists.natpair_lj);
    }
}

static void do_nb_verlet_fep(nbnxn_pairlist_set_t *nbl_lists,
                             t_forcerec           *fr,
                             rvec                  x[],
                             rvec                  f[],
                             t_mdatoms            *mdatoms,
                             t_lambda             *fepvals,
                             real                 *lambda,
                             gmx_enerdata_t       *enerd,
                             int                   flags,
                             t_nrnb               *nrnb,
                             gmx_wallcycle_t       wcycle)
{
    int              donb_flags;
    nb_kernel_data_t kernel_data;
    real             lam_i[efptNR];
    real             dvdl_nb[efptNR];
    int              th;
    int              i, j;

    donb_flags = 0;
    /* Add short-range interactions */
    donb_flags |= GMX_NONBONDED_DO_SR;

    /* Currently all group scheme kernels always calculate (shift-)forces */
    if (flags & GMX_FORCE_FORCES)
    {
        donb_flags |= GMX_NONBONDED_DO_FORCE;
    }
    if (flags & GMX_FORCE_VIRIAL)
    {
        donb_flags |= GMX_NONBONDED_DO_SHIFTFORCE;
    }
    if (flags & GMX_FORCE_ENERGY)
    {
        donb_flags |= GMX_NONBONDED_DO_POTENTIAL;
    }

    kernel_data.flags  = donb_flags;
    kernel_data.lambda = lambda;
    kernel_data.dvdl   = dvdl_nb;

    kernel_data.energygrp_elec = enerd->grpp.ener[egCOULSR];
    kernel_data.energygrp_vdw  = enerd->grpp.ener[egLJSR];

    /* reset free energy components */
    for (i = 0; i < efptNR; i++)
    {
        dvdl_nb[i]  = 0;
    }

    GMX_ASSERT(gmx_omp_nthreads_get(emntNonbonded) == nbl_lists->nnbl, "Number of lists should be same as number of NB threads");

    wallcycle_sub_start(wcycle, ewcsNONBONDED);
#pragma omp parallel for schedule(static) num_threads(nbl_lists->nnbl)
    for (th = 0; th < nbl_lists->nnbl; th++)
    {
        try
        {
            gmx_nb_free_energy_kernel(nbl_lists->nbl_fep[th],
                                      x, f, fr, mdatoms, &kernel_data, nrnb);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    if (fepvals->sc_alpha != 0)
    {
        enerd->dvdl_nonlin[efptVDW]  += dvdl_nb[efptVDW];
        enerd->dvdl_nonlin[efptCOUL] += dvdl_nb[efptCOUL];
    }
    else
    {
        enerd->dvdl_lin[efptVDW]  += dvdl_nb[efptVDW];
        enerd->dvdl_lin[efptCOUL] += dvdl_nb[efptCOUL];
    }

    /* If we do foreign lambda and we have soft-core interactions
     * we have to recalculate the (non-linear) energies contributions.
     */
    if (fepvals->n_lambda > 0 && (flags & GMX_FORCE_DHDL) && fepvals->sc_alpha != 0)
    {
        kernel_data.flags          = (donb_flags & ~(GMX_NONBONDED_DO_FORCE | GMX_NONBONDED_DO_SHIFTFORCE)) | GMX_NONBONDED_DO_FOREIGNLAMBDA;
        kernel_data.lambda         = lam_i;
        kernel_data.energygrp_elec = enerd->foreign_grpp.ener[egCOULSR];
        kernel_data.energygrp_vdw  = enerd->foreign_grpp.ener[egLJSR];
        /* Note that we add to kernel_data.dvdl, but ignore the result */

        for (i = 0; i < enerd->n_lambda; i++)
        {
            for (j = 0; j < efptNR; j++)
            {
                lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
            }
            reset_foreign_enerdata(enerd);
#pragma omp parallel for schedule(static) num_threads(nbl_lists->nnbl)
            for (th = 0; th < nbl_lists->nnbl; th++)
            {
                try
                {
                    gmx_nb_free_energy_kernel(nbl_lists->nbl_fep[th],
                                              x, f, fr, mdatoms, &kernel_data, nrnb);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

            sum_epot(&(enerd->foreign_grpp), enerd->foreign_term);
            enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
        }
    }

    wallcycle_sub_stop(wcycle, ewcsNONBONDED);
}

gmx_bool use_GPU(const nonbonded_verlet_t *nbv)
{
    return nbv != nullptr && nbv->bUseGPU;
}

static gmx_inline void clear_rvecs_omp(int n, rvec v[])
{
    int nth = gmx_omp_nthreads_get_simple_rvec_task(emntDefault, n);

    /* Note that we would like to avoid this conditional by putting it
     * into the omp pragma instead, but then we still take the full
     * omp parallel for overhead (at least with gcc5).
     */
    if (nth == 1)
    {
        for (int i = 0; i < n; i++)
        {
            clear_rvec(v[i]);
        }
    }
    else
    {
#pragma omp parallel for num_threads(nth) schedule(static)
        for (int i = 0; i < n; i++)
        {
            clear_rvec(v[i]);
        }
    }
}

/*! \brief Return an estimate of the average kinetic energy or 0 when unreliable
 *
 * \param groupOptions  Group options, containing T-coupling options
 */
static real averageKineticEnergyEstimate(const t_grpopts &groupOptions)
{
    real nrdfCoupled   = 0;
    real nrdfUncoupled = 0;
    real kineticEnergy = 0;
    for (int g = 0; g < groupOptions.ngtc; g++)
    {
        if (groupOptions.tau_t[g] >= 0)
        {
            nrdfCoupled   += groupOptions.nrdf[g];
            kineticEnergy += groupOptions.nrdf[g]*0.5*groupOptions.ref_t[g]*BOLTZ;
        }
        else
        {
            nrdfUncoupled += groupOptions.nrdf[g];
        }
    }

    /* This conditional with > also catches nrdf=0 */
    if (nrdfCoupled > nrdfUncoupled)
    {
        return kineticEnergy*(nrdfCoupled + nrdfUncoupled)/nrdfCoupled;
    }
    else
    {
        return 0;
    }
}

/*! \brief This routine checks that the potential energy is finite.
 *
 * Always checks that the potential energy is finite. If step equals
 * inputrec.init_step also checks that the magnitude of the potential energy
 * is reasonable. Terminates with a fatal error when a check fails.
 * Note that passing this check does not guarantee finite forces,
 * since those use slightly different arithmetics. But in most cases
 * there is just a narrow coordinate range where forces are not finite
 * and energies are finite.
 *
 * \param[in] step      The step number, used for checking and printing
 * \param[in] enerd     The energy data; the non-bonded group energies need to be added to enerd.term[F_EPOT] before calling this routine
 * \param[in] inputrec  The input record
 */
static void checkPotentialEnergyValidity(gmx_int64_t           step,
                                         const gmx_enerdata_t &enerd,
                                         const t_inputrec     &inputrec)
{
    /* Threshold valid for comparing absolute potential energy against
     * the kinetic energy. Normally one should not consider absolute
     * potential energy values, but with a factor of one million
     * we should never get false positives.
     */
    constexpr real c_thresholdFactor = 1e6;

    bool           energyIsNotFinite    = !std::isfinite(enerd.term[F_EPOT]);
    real           averageKineticEnergy = 0;
    /* We only check for large potential energy at the initial step,
     * because that is by far the most likely step for this too occur
     * and because computing the average kinetic energy is not free.
     * Note: nstcalcenergy >> 1 often does not allow to catch large energies
     * before they become NaN.
     */
    if (step == inputrec.init_step && EI_DYNAMICS(inputrec.eI))
    {
        averageKineticEnergy = averageKineticEnergyEstimate(inputrec.opts);
    }

    if (energyIsNotFinite || (averageKineticEnergy > 0 &&
                              enerd.term[F_EPOT] > c_thresholdFactor*averageKineticEnergy))
    {
        gmx_fatal(FARGS, "Step %" GMX_PRId64 ": The total potential energy is %g, which is %s. The LJ and electrostatic contributions to the energy are %g and %g, respectively. A %s potential energy can be caused by overlapping interactions in bonded interactions or very large%s coordinate values. Usually this is caused by a badly- or non-equilibrated initial configuration, incorrect interactions or parameters in the topology.",
                  step,
                  enerd.term[F_EPOT],
                  energyIsNotFinite ? "not finite" : "extremely high",
                  enerd.term[F_LJ],
                  enerd.term[F_COUL_SR],
                  energyIsNotFinite ? "non-finite" : "very high",
                  energyIsNotFinite ? " or Nan" : "");
    }
}

/*! \brief Compute forces and/or energies for special algorithms
 *
 * The intention is to collect all calls to algorithms that compute
 * forces on local atoms only and that do not contribute to the local
 * virial sum (but add their virial contribution separately).
 * Eventually these should likely all become ForceProviders.
 * Within this function the intention is to have algorithms that do
 * global communication at the end, so global barriers within the MD loop
 * are as close together as possible.
 *
 * \param[in]     fplog            The log file
 * \param[in]     cr               The communication record
 * \param[in]     inputrec         The input record
 * \param[in]     step             The current MD step
 * \param[in]     t                The current time
 * \param[in,out] wcycle           Wallcycle accounting struct
 * \param[in,out] forceProviders   Pointer to a list of force providers
 * \param[in]     box              The unit cell
 * \param[in]     x                The coordinates
 * \param[in]     mdatoms          Per atom properties
 * \param[in]     lambda           Array of free-energy lambda values
 * \param[in]     forceFlags       Flags that tell whether we should compute forces/energies/virial
 * \param[in,out] forceWithVirial  Force and virial buffers
 * \param[in,out] enerd            Energy buffer
 * \param[in,out] ed               Essential dynamics pointer
 * \param[in]     bNS              Tells if we did neighbor searching this step, used for ED sampling
 *
 * \todo Remove bNS, which is used incorrectly.
 * \todo Convert all other algorithms called here to ForceProviders.
 */
static void
computeSpecialForces(FILE             *fplog,
                     t_commrec        *cr,
                     t_inputrec       *inputrec,
                     gmx_int64_t       step,
                     double            t,
                     gmx_wallcycle_t   wcycle,
                     ForceProviders   *forceProviders,
                     matrix            box,
                     rvec             *x,
                     t_mdatoms        *mdatoms,
                     real             *lambda,
                     int               forceFlags,
                     ForceWithVirial  *forceWithVirial,
                     gmx_enerdata_t   *enerd,
                     gmx_edsam_t       ed,
                     gmx_bool          bNS)
{
    const bool computeForces = (forceFlags & GMX_FORCE_FORCES);

    /* NOTE: Currently all ForceProviders only provide forces.
     *       When they also provide energies, remove this conditional.
     */
    if (computeForces)
    {
        /* Collect forces from modules */
        forceProviders->calculateForces(cr, mdatoms, box, t, x, forceWithVirial);
    }

    if (inputrec->bPull && pull_have_potential(inputrec->pull_work))
    {
        pull_potential_wrapper(cr, inputrec, box, x,
                               forceWithVirial,
                               mdatoms, enerd, lambda, t,
                               wcycle);

        if (inputrec->bDoAwh)
        {
            Awh &awh = *inputrec->awh;
            enerd->term[F_COM_PULL] +=
                awh.applyBiasForcesAndUpdateBias(inputrec->ePBC, *mdatoms, box,
                                                 forceWithVirial,
                                                 t, step, wcycle, fplog);
        }
    }

    rvec *f = as_rvec_array(forceWithVirial->force_.data());

    /* Add the forces from enforced rotation potentials (if any) */
    if (inputrec->bRot)
    {
        wallcycle_start(wcycle, ewcROTadd);
        enerd->term[F_COM_PULL] += add_rot_forces(inputrec->rot, f, cr, step, t);
        wallcycle_stop(wcycle, ewcROTadd);
    }

    if (ed)
    {
        /* Note that since init_edsam() is called after the initialization
         * of forcerec, edsam doesn't request the noVirSum force buffer.
         * Thus if no other algorithm (e.g. PME) requires it, the forces
         * here will contribute to the virial.
         */
        do_flood(cr, inputrec, x, f, ed, box, step, bNS);
    }

    /* Add forces from interactive molecular dynamics (IMD), if bIMD == TRUE. */
    if (inputrec->bIMD && computeForces)
    {
        IMD_apply_forces(inputrec->bIMD, inputrec->imd, cr, f, wcycle);
    }
}

/*! \brief Launch the prepare_step and spread stages of PME GPU.
 *
 * \param[in]  pmedata       The PME structure
 * \param[in]  box           The box matrix
 * \param[in]  x             Coordinate array
 * \param[in]  flags         Force flags
 * \param[in]  wcycle        The wallcycle structure
 */
static inline void launchPmeGpuSpread(gmx_pme_t      *pmedata,
                                      matrix          box,
                                      rvec            x[],
                                      int             flags,
                                      gmx_wallcycle_t wcycle)
{
    int pmeFlags = GMX_PME_SPREAD | GMX_PME_SOLVE;
    pmeFlags |= (flags & GMX_FORCE_FORCES) ? GMX_PME_CALC_F : 0;
    pmeFlags |= (flags & GMX_FORCE_VIRIAL) ? GMX_PME_CALC_ENER_VIR : 0;

    pme_gpu_prepare_computation(pmedata, flags & GMX_FORCE_DYNAMICBOX, box, wcycle, pmeFlags);
    pme_gpu_launch_spread(pmedata, x, wcycle);
}

/*! \brief Launch the FFT and gather stages of PME GPU
 *
 * This function only implements setting the output forces (no accumulation).
 *
 * \param[in]  pmedata        The PME structure
 * \param[in]  wcycle         The wallcycle structure
 */
static void launchPmeGpuFftAndGather(gmx_pme_t        *pmedata,
                                     gmx_wallcycle_t   wcycle)
{
    pme_gpu_launch_complex_transforms(pmedata, wcycle);
    pme_gpu_launch_gather(pmedata, wcycle, PmeForceOutputHandling::Set);
}

/*! \brief
 *  Polling wait for either of the PME or nonbonded GPU tasks.
 *
 * Instead of a static order in waiting for GPU tasks, this function
 * polls checking which of the two tasks completes first, and does the
 * associated force buffer reduction overlapped with the other task.
 * By doing that, unlike static scheduling order, it can always overlap
 * one of the reductions, regardless of the GPU task completion order.
 *
 * \param[in]     nbv              Nonbonded verlet structure
 * \param[in]     pmedata          PME module data
 * \param[in,out] force            Force array to reduce task outputs into.
 * \param[in,out] forceWithVirial  Force and virial buffers
 * \param[in,out] fshift           Shift force output vector results are reduced into
 * \param[in,out] enerd            Energy data structure results are reduced into
 * \param[in]     flags            Force flags
 * \param[in]     wcycle           The wallcycle structure
 */
static void alternatePmeNbGpuWaitReduce(nonbonded_verlet_t             *nbv,
                                        const gmx_pme_t                *pmedata,
                                        gmx::PaddedArrayRef<gmx::RVec> *force,
                                        ForceWithVirial                *forceWithVirial,
                                        rvec                            fshift[],
                                        gmx_enerdata_t                 *enerd,
                                        int                             flags,
                                        gmx_wallcycle_t                 wcycle)
{
    bool isPmeGpuDone = false;
    bool isNbGpuDone  = false;


    gmx::ArrayRef<const gmx::RVec> pmeGpuForces;

    while (!isPmeGpuDone || !isNbGpuDone)
    {
        if (!isPmeGpuDone)
        {
            matrix            vir_Q;
            real              Vlr_q;

            GpuTaskCompletion completionType = (isNbGpuDone) ? GpuTaskCompletion::Wait : GpuTaskCompletion::Check;
            isPmeGpuDone = pme_gpu_try_finish_task(pmedata, wcycle, &pmeGpuForces,
                                                   vir_Q, &Vlr_q, completionType);

            if (isPmeGpuDone)
            {
                pme_gpu_reduce_outputs(wcycle, forceWithVirial, pmeGpuForces,
                                       enerd, vir_Q, Vlr_q);
            }
        }

        if (!isNbGpuDone)
        {
            GpuTaskCompletion completionType = (isPmeGpuDone) ? GpuTaskCompletion::Wait : GpuTaskCompletion::Check;
            wallcycle_start_nocount(wcycle, ewcWAIT_GPU_NB_L);
            isNbGpuDone = nbnxn_gpu_try_finish_task(nbv->gpu_nbv,
                                                    flags, eatLocal,
                                                    enerd->grpp.ener[egLJSR], enerd->grpp.ener[egCOULSR],
                                                    fshift, completionType);
            wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L);
            // To get the call count right, when the task finished we
            // issue a start/stop.
            // TODO: move the ewcWAIT_GPU_NB_L cycle counting into nbnxn_gpu_try_finish_task()
            // and ewcNB_XF_BUF_OPS counting into nbnxn_atomdata_add_nbat_f_to_f().
            if (isNbGpuDone)
            {
                wallcycle_start(wcycle, ewcWAIT_GPU_NB_L);
                wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L);

                wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
                wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
                nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatLocal,
                                               nbv->nbat, as_rvec_array(force->data()));
                wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
                wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
            }
        }
    }
}

/*! \brief
 *  Launch the dynamic rolling pruning GPU task.
 *
 *  We currently alternate local/non-local list pruning in odd-even steps
 *  (only pruning every second step without DD).
 *
 * \param[in]     cr               The communication record
 * \param[in]     nbv              Nonbonded verlet structure
 * \param[in]     inputrec         The input record
 * \param[in]     step             The current MD step
 */
static inline void launchGpuRollingPruning(const t_commrec          *cr,
                                           const nonbonded_verlet_t *nbv,
                                           const t_inputrec         *inputrec,
                                           const gmx_int64_t         step)
{
    /* We should not launch the rolling pruning kernel at a search
     * step or just before search steps, since that's useless.
     * Without domain decomposition we prune at even steps.
     * With domain decomposition we alternate local and non-local
     * pruning at even and odd steps.
     */
    int  numRollingParts     = nbv->listParams->numRollingParts;
    GMX_ASSERT(numRollingParts == nbv->listParams->nstlistPrune/2, "Since we alternate local/non-local at even/odd steps, we need numRollingParts<=nstlistPrune/2 for correctness and == for efficiency");
    int  stepWithCurrentList = step - nbv->grp[eintLocal].nbl_lists.outerListCreationStep;
    bool stepIsEven          = ((stepWithCurrentList & 1) == 0);
    if (stepWithCurrentList > 0 &&
        stepWithCurrentList < inputrec->nstlist - 1 &&
        (stepIsEven || DOMAINDECOMP(cr)))
    {
        nbnxn_gpu_launch_kernel_pruneonly(nbv->gpu_nbv,
                                          stepIsEven ? eintLocal : eintNonlocal,
                                          numRollingParts);
    }
}

static void do_force_cutsVERLET(FILE *fplog, t_commrec *cr,
                                t_inputrec *inputrec,
                                gmx_int64_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                                gmx_localtop_t *top,
                                gmx_groups_t gmx_unused *groups,
                                matrix box, gmx::PaddedArrayRef<gmx::RVec> x, history_t *hist,
                                gmx::PaddedArrayRef<gmx::RVec> force,
                                tensor vir_force,
                                t_mdatoms *mdatoms,
                                gmx_enerdata_t *enerd, t_fcdata *fcd,
                                real *lambda, t_graph *graph,
                                t_forcerec *fr, interaction_const_t *ic,
                                gmx_vsite_t *vsite, rvec mu_tot,
                                double t, gmx_edsam_t ed,
                                gmx_bool bBornRadii,
                                int flags,
                                DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
                                DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion)
{
    int                 cg1, i, j;
    double              mu[2*DIM];
    gmx_bool            bStateChanged, bNS, bFillGrid, bCalcCGCM;
    gmx_bool            bDoForces, bUseGPU, bUseOrEmulGPU;
    rvec                vzero, box_diag;
    float               cycles_pme, cycles_wait_gpu;
    nonbonded_verlet_t *nbv = fr->nbv;

    bStateChanged = (flags & GMX_FORCE_STATECHANGED);
    bNS           = (flags & GMX_FORCE_NS) && (fr->bAllvsAll == FALSE);
    bFillGrid     = (bNS && bStateChanged);
    bCalcCGCM     = (bFillGrid && !DOMAINDECOMP(cr));
    bDoForces     = (flags & GMX_FORCE_FORCES);
    bUseGPU       = fr->nbv->bUseGPU;
    bUseOrEmulGPU = bUseGPU || (fr->nbv->emulateGpu == EmulateGpuNonbonded::Yes);

    const auto pmeRunMode = fr->pmedata ? pme_run_mode(fr->pmedata) : PmeRunMode::CPU;
    // TODO slim this conditional down - inputrec and duty checks should mean the same in proper code!
    const bool useGpuPme  = EEL_PME(fr->ic->eeltype) && thisRankHasDuty(cr, DUTY_PME) &&
        ((pmeRunMode == PmeRunMode::GPU) || (pmeRunMode == PmeRunMode::Mixed));

    /* At a search step we need to start the first balancing region
     * somewhere early inside the step after communication during domain
     * decomposition (and not during the previous step as usual).
     */
    if (bNS &&
        ddOpenBalanceRegion == DdOpenBalanceRegionBeforeForceComputation::yes)
    {
        ddOpenBalanceRegionCpu(cr->dd, DdAllowBalanceRegionReopen::yes);
    }

    cycles_wait_gpu = 0;

    const int start  = 0;
    const int homenr = mdatoms->homenr;

    clear_mat(vir_force);

    if (DOMAINDECOMP(cr))
    {
        cg1 = cr->dd->ncg_tot;
    }
    else
    {
        cg1 = top->cgs.nr;
    }
    if (fr->n_tpi > 0)
    {
        cg1--;
    }

    if (bStateChanged)
    {
        update_forcerec(fr, box);

        if (inputrecNeedMutot(inputrec))
        {
            /* Calculate total (local) dipole moment in a temporary common array.
             * This makes it possible to sum them over nodes faster.
             */
            calc_mu(start, homenr,
                    x, mdatoms->chargeA, mdatoms->chargeB, mdatoms->nChargePerturbed,
                    mu, mu+DIM);
        }
    }

    if (fr->ePBC != epbcNONE)
    {
        /* Compute shift vectors every step,
         * because of pressure coupling or box deformation!
         */
        if ((flags & GMX_FORCE_DYNAMICBOX) && bStateChanged)
        {
            calc_shifts(box, fr->shift_vec);
        }

        if (bCalcCGCM)
        {
            put_atoms_in_box_omp(fr->ePBC, box, x.subArray(0, homenr));
            inc_nrnb(nrnb, eNR_SHIFTX, homenr);
        }
        else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph)
        {
            unshift_self(graph, box, as_rvec_array(x.data()));
        }
    }

    nbnxn_atomdata_copy_shiftvec(flags & GMX_FORCE_DYNAMICBOX,
                                 fr->shift_vec, nbv->nbat);

#if GMX_MPI
    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Send particle coordinates to the pme nodes.
         * Since this is only implemented for domain decomposition
         * and domain decomposition does not use the graph,
         * we do not need to worry about shifting.
         */
        wallcycle_start(wcycle, ewcPP_PMESENDX);
        gmx_pme_send_coordinates(cr, box, as_rvec_array(x.data()),
                                 lambda[efptCOUL], lambda[efptVDW],
                                 (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY)),
                                 step);
        wallcycle_stop(wcycle, ewcPP_PMESENDX);
    }
#endif /* GMX_MPI */

    if (useGpuPme)
    {
        launchPmeGpuSpread(fr->pmedata, box, as_rvec_array(x.data()), flags, wcycle);
    }

    /* do gridding for pair search */
    if (bNS)
    {
        if (graph && bStateChanged)
        {
            /* Calculate intramolecular shift vectors to make molecules whole */
            mk_mshift(fplog, graph, fr->ePBC, box, as_rvec_array(x.data()));
        }

        clear_rvec(vzero);
        box_diag[XX] = box[XX][XX];
        box_diag[YY] = box[YY][YY];
        box_diag[ZZ] = box[ZZ][ZZ];

        wallcycle_start(wcycle, ewcNS);
        if (!fr->bDomDec)
        {
            wallcycle_sub_start(wcycle, ewcsNBS_GRID_LOCAL);
            nbnxn_put_on_grid(nbv->nbs, fr->ePBC, box,
                              0, vzero, box_diag,
                              0, mdatoms->homenr, -1, fr->cginfo, as_rvec_array(x.data()),
                              0, nullptr,
                              nbv->grp[eintLocal].kernel_type,
                              nbv->nbat);
            wallcycle_sub_stop(wcycle, ewcsNBS_GRID_LOCAL);
        }
        else
        {
            wallcycle_sub_start(wcycle, ewcsNBS_GRID_NONLOCAL);
            nbnxn_put_on_grid_nonlocal(nbv->nbs, domdec_zones(cr->dd),
                                       fr->cginfo, as_rvec_array(x.data()),
                                       nbv->grp[eintNonlocal].kernel_type,
                                       nbv->nbat);
            wallcycle_sub_stop(wcycle, ewcsNBS_GRID_NONLOCAL);
        }

        nbnxn_atomdata_set(nbv->nbat, nbv->nbs, mdatoms, fr->cginfo);
        wallcycle_stop(wcycle, ewcNS);
    }

    /* initialize the GPU atom data and copy shift vector */
    if (bUseGPU)
    {
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_NONBONDED);

        if (bNS)
        {
            nbnxn_gpu_init_atomdata(nbv->gpu_nbv, nbv->nbat);
        }

        nbnxn_gpu_upload_shiftvec(nbv->gpu_nbv, nbv->nbat);

        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);
    }

    /* do local pair search */
    if (bNS)
    {
        wallcycle_start_nocount(wcycle, ewcNS);
        wallcycle_sub_start(wcycle, ewcsNBS_SEARCH_LOCAL);
        nbnxn_make_pairlist(nbv->nbs, nbv->nbat,
                            &top->excls,
                            nbv->listParams->rlistOuter,
                            nbv->min_ci_balanced,
                            &nbv->grp[eintLocal].nbl_lists,
                            eintLocal,
                            nbv->grp[eintLocal].kernel_type,
                            nrnb);
        nbv->grp[eintLocal].nbl_lists.outerListCreationStep = step;
        if (nbv->listParams->useDynamicPruning && !bUseGPU)
        {
            nbnxnPrepareListForDynamicPruning(&nbv->grp[eintLocal].nbl_lists);
        }
        wallcycle_sub_stop(wcycle, ewcsNBS_SEARCH_LOCAL);

        if (bUseGPU)
        {
            /* initialize local pair-list on the GPU */
            nbnxn_gpu_init_pairlist(nbv->gpu_nbv,
                                    nbv->grp[eintLocal].nbl_lists.nbl[0],
                                    eintLocal);
        }
        wallcycle_stop(wcycle, ewcNS);
    }
    else
    {
        wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_sub_start(wcycle, ewcsNB_X_BUF_OPS);
        nbnxn_atomdata_copy_x_to_nbat_x(nbv->nbs, eatLocal, FALSE, as_rvec_array(x.data()),
                                        nbv->nbat);
        wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
        wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
    }

    if (bUseGPU)
    {
        if (DOMAINDECOMP(cr))
        {
            ddOpenBalanceRegionGpu(cr->dd);
        }

        wallcycle_start(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        /* launch local nonbonded work on GPU */
        do_nb_verlet(fr, ic, enerd, flags, eintLocal, enbvClearFNo,
                     step, nrnb, wcycle);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);
    }

    if (useGpuPme)
    {
        // In PME GPU and mixed mode we launch FFT / gather after the
        // X copy/transform to allow overlap as well as after the GPU NB
        // launch to avoid FFT launch overhead hijacking the CPU and delaying
        // the nonbonded kernel.
        launchPmeGpuFftAndGather(fr->pmedata, wcycle);
    }

    /* Communicate coordinates and sum dipole if necessary +
       do non-local pair search */
    if (DOMAINDECOMP(cr))
    {
        if (bNS)
        {
            wallcycle_start_nocount(wcycle, ewcNS);
            wallcycle_sub_start(wcycle, ewcsNBS_SEARCH_NONLOCAL);

            nbnxn_make_pairlist(nbv->nbs, nbv->nbat,
                                &top->excls,
                                nbv->listParams->rlistOuter,
                                nbv->min_ci_balanced,
                                &nbv->grp[eintNonlocal].nbl_lists,
                                eintNonlocal,
                                nbv->grp[eintNonlocal].kernel_type,
                                nrnb);
            nbv->grp[eintNonlocal].nbl_lists.outerListCreationStep = step;
            if (nbv->listParams->useDynamicPruning && !bUseGPU)
            {
                nbnxnPrepareListForDynamicPruning(&nbv->grp[eintNonlocal].nbl_lists);
            }
            wallcycle_sub_stop(wcycle, ewcsNBS_SEARCH_NONLOCAL);

            if (nbv->grp[eintNonlocal].kernel_type == nbnxnk8x8x8_GPU)
            {
                /* initialize non-local pair-list on the GPU */
                nbnxn_gpu_init_pairlist(nbv->gpu_nbv,
                                        nbv->grp[eintNonlocal].nbl_lists.nbl[0],
                                        eintNonlocal);
            }
            wallcycle_stop(wcycle, ewcNS);
        }
        else
        {
            wallcycle_start(wcycle, ewcMOVEX);
            dd_move_x(cr->dd, box, as_rvec_array(x.data()));
            wallcycle_stop(wcycle, ewcMOVEX);

            wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
            wallcycle_sub_start(wcycle, ewcsNB_X_BUF_OPS);
            nbnxn_atomdata_copy_x_to_nbat_x(nbv->nbs, eatNonlocal, FALSE, as_rvec_array(x.data()),
                                            nbv->nbat);
            wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
            wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        }

        if (bUseGPU)
        {
            wallcycle_start(wcycle, ewcLAUNCH_GPU);
            wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_NONBONDED);
            /* launch non-local nonbonded tasks on GPU */
            do_nb_verlet(fr, ic, enerd, flags, eintNonlocal, enbvClearFNo,
                         step, nrnb, wcycle);
            wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_NONBONDED);
            wallcycle_stop(wcycle, ewcLAUNCH_GPU);
        }
    }

    if (bUseGPU)
    {
        /* launch D2H copy-back F */
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        if (DOMAINDECOMP(cr))
        {
            nbnxn_gpu_launch_cpyback(nbv->gpu_nbv, nbv->nbat,
                                     flags, eatNonlocal);
        }
        nbnxn_gpu_launch_cpyback(nbv->gpu_nbv, nbv->nbat,
                                 flags, eatLocal);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);
    }

    if (bStateChanged && inputrecNeedMutot(inputrec))
    {
        if (PAR(cr))
        {
            gmx_sumd(2*DIM, mu, cr);
            ddReopenBalanceRegionCpu(cr->dd);
        }

        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                fr->mu_tot[i][j] = mu[i*DIM + j];
            }
        }
    }
    if (fr->efep == efepNO)
    {
        copy_rvec(fr->mu_tot[0], mu_tot);
    }
    else
    {
        for (j = 0; j < DIM; j++)
        {
            mu_tot[j] =
                (1.0 - lambda[efptCOUL])*fr->mu_tot[0][j] +
                lambda[efptCOUL]*fr->mu_tot[1][j];
        }
    }

    /* Reset energies */
    reset_enerdata(enerd);
    clear_rvecs(SHIFTS, fr->fshift);

    if (DOMAINDECOMP(cr) && !thisRankHasDuty(cr, DUTY_PME))
    {
        wallcycle_start(wcycle, ewcPPDURINGPME);
        dd_force_flop_start(cr->dd, nrnb);
    }

    if (inputrec->bRot)
    {
        wallcycle_start(wcycle, ewcROT);
        do_rotation(cr, inputrec, box, as_rvec_array(x.data()), t, step, bNS);
        wallcycle_stop(wcycle, ewcROT);
    }

    /* Temporary solution until all routines take PaddedRVecVector */
    rvec *f = as_rvec_array(force.data());

    /* Start the force cycle counter.
     * Note that a different counter is used for dynamic load balancing.
     */
    wallcycle_start(wcycle, ewcFORCE);

    gmx::ArrayRef<gmx::RVec> forceRef = force;
    if (bDoForces)
    {
        /* If we need to compute the virial, we might need a separate
         * force buffer for algorithms for which the virial is calculated
         * directly, such as PME.
         */
        if ((flags & GMX_FORCE_VIRIAL) && fr->haveDirectVirialContributions)
        {
            forceRef = *fr->forceBufferForDirectVirialContributions;

            /* We only compute forces on local atoms. Note that vsites can
             * spread to non-local atoms, but that part of the buffer is
             * cleared separately in the vsite spreading code.
             */
            clear_rvecs_omp(mdatoms->homenr, as_rvec_array(forceRef.data()));
        }
        /* Clear the short- and long-range forces */
        clear_rvecs_omp(fr->natoms_force_constr, f);
    }

    /* forceWithVirial uses the local atom range only */
    ForceWithVirial forceWithVirial(forceRef, flags & GMX_FORCE_VIRIAL);

    if (inputrec->bPull && pull_have_constraint(inputrec->pull_work))
    {
        clear_pull_forces(inputrec->pull_work);
    }

    /* We calculate the non-bonded forces, when done on the CPU, here.
     * We do this before calling do_force_lowlevel, because in that
     * function, the listed forces are calculated before PME, which
     * does communication.  With this order, non-bonded and listed
     * force calculation imbalance can be balanced out by the domain
     * decomposition load balancing.
     */

    if (!bUseOrEmulGPU)
    {
        do_nb_verlet(fr, ic, enerd, flags, eintLocal, enbvClearFYes,
                     step, nrnb, wcycle);
    }

    if (fr->efep != efepNO)
    {
        /* Calculate the local and non-local free energy interactions here.
         * Happens here on the CPU both with and without GPU.
         */
        if (fr->nbv->grp[eintLocal].nbl_lists.nbl_fep[0]->nrj > 0)
        {
            do_nb_verlet_fep(&fr->nbv->grp[eintLocal].nbl_lists,
                             fr, as_rvec_array(x.data()), f, mdatoms,
                             inputrec->fepvals, lambda,
                             enerd, flags, nrnb, wcycle);
        }

        if (DOMAINDECOMP(cr) &&
            fr->nbv->grp[eintNonlocal].nbl_lists.nbl_fep[0]->nrj > 0)
        {
            do_nb_verlet_fep(&fr->nbv->grp[eintNonlocal].nbl_lists,
                             fr, as_rvec_array(x.data()), f, mdatoms,
                             inputrec->fepvals, lambda,
                             enerd, flags, nrnb, wcycle);
        }
    }

    if (!bUseOrEmulGPU)
    {
        int aloc;

        if (DOMAINDECOMP(cr))
        {
            do_nb_verlet(fr, ic, enerd, flags, eintNonlocal, enbvClearFNo,
                         step, nrnb, wcycle);
        }

        if (!bUseOrEmulGPU)
        {
            aloc = eintLocal;
        }
        else
        {
            aloc = eintNonlocal;
        }

        /* Add all the non-bonded force to the normal force array.
         * This can be split into a local and a non-local part when overlapping
         * communication with calculation with domain decomposition.
         */
        wallcycle_stop(wcycle, ewcFORCE);
        wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
        nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatAll, nbv->nbat, f);
        wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
        wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_start_nocount(wcycle, ewcFORCE);

        /* if there are multiple fshift output buffers reduce them */
        if ((flags & GMX_FORCE_VIRIAL) &&
            nbv->grp[aloc].nbl_lists.nnbl > 1)
        {
            /* This is not in a subcounter because it takes a
               negligible and constant-sized amount of time */
            nbnxn_atomdata_add_nbat_fshift_to_fshift(nbv->nbat,
                                                     fr->fshift);
        }
    }

    /* update QMMMrec, if necessary */
    if (fr->bQMMM)
    {
        update_QMMMrec(cr, fr, as_rvec_array(x.data()), mdatoms, box);
    }

    /* Compute the bonded and non-bonded energies and optionally forces */
    do_force_lowlevel(fr, inputrec, &(top->idef),
                      cr, nrnb, wcycle, mdatoms,
                      as_rvec_array(x.data()), hist, f, &forceWithVirial, enerd, fcd, top, fr->born,
                      bBornRadii, box,
                      inputrec->fepvals, lambda, graph, &(top->excls), fr->mu_tot,
                      flags, &cycles_pme);

    wallcycle_stop(wcycle, ewcFORCE);

    computeSpecialForces(fplog, cr, inputrec, step, t, wcycle,
                         fr->forceProviders, box, as_rvec_array(x.data()), mdatoms, lambda,
                         flags, &forceWithVirial, enerd,
                         ed, bNS);

    if (bUseOrEmulGPU)
    {
        /* wait for non-local forces (or calculate in emulation mode) */
        if (DOMAINDECOMP(cr))
        {
            if (bUseGPU)
            {
                wallcycle_start(wcycle, ewcWAIT_GPU_NB_NL);
                nbnxn_gpu_wait_finish_task(nbv->gpu_nbv,
                                           flags, eatNonlocal,
                                           enerd->grpp.ener[egLJSR], enerd->grpp.ener[egCOULSR],
                                           fr->fshift);
                cycles_wait_gpu += wallcycle_stop(wcycle, ewcWAIT_GPU_NB_NL);
            }
            else
            {
                wallcycle_start_nocount(wcycle, ewcFORCE);
                do_nb_verlet(fr, ic, enerd, flags, eintNonlocal, enbvClearFYes,
                             step, nrnb, wcycle);
                wallcycle_stop(wcycle, ewcFORCE);
            }
            wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
            wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
            /* skip the reduction if there was no non-local work to do */
            if (nbv->grp[eintNonlocal].nbl_lists.nbl[0]->nsci > 0)
            {
                nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatNonlocal,
                                               nbv->nbat, f);
            }
            wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
            wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* We are done with the CPU compute.
         * We will now communicate the non-local forces.
         * If we use a GPU this will overlap with GPU work, so in that case
         * we do not close the DD force balancing region here.
         */
        if (ddCloseBalanceRegion == DdCloseBalanceRegionAfterForceComputation::yes)
        {
            ddCloseBalanceRegionCpu(cr->dd);
        }
        if (bDoForces)
        {
            wallcycle_start(wcycle, ewcMOVEF);
            dd_move_f(cr->dd, f, fr->fshift);
            wallcycle_stop(wcycle, ewcMOVEF);
        }
    }

    // With both nonbonded and PME offloaded a GPU on the same rank, we use
    // an alternating wait/reduction scheme.
    bool alternateGpuWait = (!c_disableAlternatingWait && useGpuPme && bUseGPU && !DOMAINDECOMP(cr));
    if (alternateGpuWait)
    {
        alternatePmeNbGpuWaitReduce(fr->nbv, fr->pmedata, &force, &forceWithVirial, fr->fshift, enerd, flags, wcycle);
    }

    if (!alternateGpuWait && useGpuPme)
    {
        gmx::ArrayRef<const gmx::RVec> pmeGpuForces;
        matrix vir_Q;
        real   Vlr_q;
        pme_gpu_wait_finish_task(fr->pmedata, wcycle, &pmeGpuForces, vir_Q, &Vlr_q);
        pme_gpu_reduce_outputs(wcycle, &forceWithVirial, pmeGpuForces, enerd, vir_Q, Vlr_q);
    }

    /* Wait for local GPU NB outputs on the non-alternating wait path */
    if (!alternateGpuWait && bUseGPU)
    {
        /* Measured overhead on CUDA and OpenCL with(out) GPU sharing
         * is between 0.5 and 1.5 Mcycles. So 2 MCycles is an overestimate,
         * but even with a step of 0.1 ms the difference is less than 1%
         * of the step time.
         */
        const float gpuWaitApiOverheadMargin = 2e6f; /* cycles */

        wallcycle_start(wcycle, ewcWAIT_GPU_NB_L);
        nbnxn_gpu_wait_finish_task(nbv->gpu_nbv,
                                   flags, eatLocal,
                                   enerd->grpp.ener[egLJSR], enerd->grpp.ener[egCOULSR],
                                   fr->fshift);
        float cycles_tmp = wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L);

        if (ddCloseBalanceRegion == DdCloseBalanceRegionAfterForceComputation::yes)
        {
            DdBalanceRegionWaitedForGpu waitedForGpu = DdBalanceRegionWaitedForGpu::yes;
            if (bDoForces && cycles_tmp <= gpuWaitApiOverheadMargin)
            {
                /* We measured few cycles, it could be that the kernel
                 * and transfer finished earlier and there was no actual
                 * wait time, only API call overhead.
                 * Then the actual time could be anywhere between 0 and
                 * cycles_wait_est. We will use half of cycles_wait_est.
                 */
                waitedForGpu = DdBalanceRegionWaitedForGpu::no;
            }
            ddCloseBalanceRegionGpu(cr->dd, cycles_wait_gpu, waitedForGpu);
        }
    }

    if (fr->nbv->emulateGpu == EmulateGpuNonbonded::Yes)
    {
        // NOTE: emulation kernel is not included in the balancing region,
        // but emulation mode does not target performance anyway
        wallcycle_start_nocount(wcycle, ewcFORCE);
        do_nb_verlet(fr, ic, enerd, flags, eintLocal,
                     DOMAINDECOMP(cr) ? enbvClearFNo : enbvClearFYes,
                     step, nrnb, wcycle);
        wallcycle_stop(wcycle, ewcFORCE);
    }

    if (bUseGPU)
    {
        /* now clear the GPU outputs while we finish the step on the CPU */
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        nbnxn_gpu_clear_outputs(nbv->gpu_nbv, flags);

        /* Is dynamic pair-list pruning activated? */
        if (nbv->listParams->useDynamicPruning)
        {
            launchGpuRollingPruning(cr, nbv, inputrec, step);
        }
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_NONBONDED);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);

        // TODO: move here the PME buffer clearing call pme_gpu_reinit_computation()
    }

    /* Do the nonbonded GPU (or emulation) force buffer reduction
     * on the non-alternating path. */
    if (bUseOrEmulGPU && !alternateGpuWait)
    {
        wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
        nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatLocal,
                                       nbv->nbat, f);
        wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
        wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
    }

    if (DOMAINDECOMP(cr))
    {
        dd_force_flop_stop(cr->dd, nrnb);
    }

    if (bDoForces)
    {
        /* If we have NoVirSum forces, but we do not calculate the virial,
         * we sum fr->f_novirsum=f later.
         */
        if (vsite && !(fr->haveDirectVirialContributions && !(flags & GMX_FORCE_VIRIAL)))
        {
            wallcycle_start(wcycle, ewcVSITESPREAD);
            spread_vsite_f(vsite, as_rvec_array(x.data()), f, fr->fshift, FALSE, nullptr, nrnb,
                           &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
            wallcycle_stop(wcycle, ewcVSITESPREAD);
        }

        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Calculation of the virial must be done after vsites! */
            calc_virial(0, mdatoms->homenr, as_rvec_array(x.data()), f,
                        vir_force, graph, box, nrnb, fr, inputrec->ePBC);
        }
    }

    if (PAR(cr) && !thisRankHasDuty(cr, DUTY_PME))
    {
        /* In case of node-splitting, the PP nodes receive the long-range
         * forces, virial and energy from the PME nodes here.
         */
        pme_receive_force_ener(cr, &forceWithVirial, enerd, wcycle);
    }

    if (bDoForces)
    {
        post_process_forces(cr, step, nrnb, wcycle,
                            top, box, as_rvec_array(x.data()), f, &forceWithVirial,
                            vir_force, mdatoms, graph, fr, vsite,
                            flags);
    }

    if (flags & GMX_FORCE_ENERGY)
    {
        /* Sum the potential energy terms from group contributions */
        sum_epot(&(enerd->grpp), enerd->term);

        if (!EI_TPI(inputrec->eI))
        {
            checkPotentialEnergyValidity(step, *enerd, *inputrec);
        }
    }
}

static void do_force_cutsGROUP(FILE *fplog, t_commrec *cr,
                               t_inputrec *inputrec,
                               gmx_int64_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                               gmx_localtop_t *top,
                               gmx_groups_t *groups,
                               matrix box, gmx::PaddedArrayRef<gmx::RVec> x, history_t *hist,
                               gmx::PaddedArrayRef<gmx::RVec> force,
                               tensor vir_force,
                               t_mdatoms *mdatoms,
                               gmx_enerdata_t *enerd, t_fcdata *fcd,
                               real *lambda, t_graph *graph,
                               t_forcerec *fr, gmx_vsite_t *vsite, rvec mu_tot,
                               double t, gmx_edsam_t ed,
                               gmx_bool bBornRadii,
                               int flags,
                               DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
                               DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion)
{
    int        cg0, cg1, i, j;
    double     mu[2*DIM];
    gmx_bool   bStateChanged, bNS, bFillGrid, bCalcCGCM;
    gmx_bool   bDoForces;
    float      cycles_pme;

    const int  start  = 0;
    const int  homenr = mdatoms->homenr;

    clear_mat(vir_force);

    cg0 = 0;
    if (DOMAINDECOMP(cr))
    {
        cg1 = cr->dd->ncg_tot;
    }
    else
    {
        cg1 = top->cgs.nr;
    }
    if (fr->n_tpi > 0)
    {
        cg1--;
    }

    bStateChanged  = (flags & GMX_FORCE_STATECHANGED);
    bNS            = (flags & GMX_FORCE_NS) && (fr->bAllvsAll == FALSE);
    /* Should we perform the long-range nonbonded evaluation inside the neighborsearching? */
    bFillGrid      = (bNS && bStateChanged);
    bCalcCGCM      = (bFillGrid && !DOMAINDECOMP(cr));
    bDoForces      = (flags & GMX_FORCE_FORCES);

    if (bStateChanged)
    {
        update_forcerec(fr, box);

        if (inputrecNeedMutot(inputrec))
        {
            /* Calculate total (local) dipole moment in a temporary common array.
             * This makes it possible to sum them over nodes faster.
             */
            calc_mu(start, homenr,
                    x, mdatoms->chargeA, mdatoms->chargeB, mdatoms->nChargePerturbed,
                    mu, mu+DIM);
        }
    }

    if (fr->ePBC != epbcNONE)
    {
        /* Compute shift vectors every step,
         * because of pressure coupling or box deformation!
         */
        if ((flags & GMX_FORCE_DYNAMICBOX) && bStateChanged)
        {
            calc_shifts(box, fr->shift_vec);
        }

        if (bCalcCGCM)
        {
            put_charge_groups_in_box(fplog, cg0, cg1, fr->ePBC, box,
                                     &(top->cgs), as_rvec_array(x.data()), fr->cg_cm);
            inc_nrnb(nrnb, eNR_CGCM, homenr);
            inc_nrnb(nrnb, eNR_RESETX, cg1-cg0);
        }
        else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph)
        {
            unshift_self(graph, box, as_rvec_array(x.data()));
        }
    }
    else if (bCalcCGCM)
    {
        calc_cgcm(fplog, cg0, cg1, &(top->cgs), as_rvec_array(x.data()), fr->cg_cm);
        inc_nrnb(nrnb, eNR_CGCM, homenr);
    }

    if (bCalcCGCM && gmx_debug_at)
    {
        pr_rvecs(debug, 0, "cgcm", fr->cg_cm, top->cgs.nr);
    }

#if GMX_MPI
    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Send particle coordinates to the pme nodes.
         * Since this is only implemented for domain decomposition
         * and domain decomposition does not use the graph,
         * we do not need to worry about shifting.
         */
        wallcycle_start(wcycle, ewcPP_PMESENDX);
        gmx_pme_send_coordinates(cr, box, as_rvec_array(x.data()),
                                 lambda[efptCOUL], lambda[efptVDW],
                                 (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY)),
                                 step);
        wallcycle_stop(wcycle, ewcPP_PMESENDX);
    }
#endif /* GMX_MPI */

    /* Communicate coordinates and sum dipole if necessary */
    if (DOMAINDECOMP(cr))
    {
        wallcycle_start(wcycle, ewcMOVEX);
        dd_move_x(cr->dd, box, as_rvec_array(x.data()));
        wallcycle_stop(wcycle, ewcMOVEX);
        /* No GPU support, no move_x overlap, so reopen the balance region here */
        if (ddOpenBalanceRegion == DdOpenBalanceRegionBeforeForceComputation::yes)
        {
            ddReopenBalanceRegionCpu(cr->dd);
        }
    }

    if (inputrecNeedMutot(inputrec))
    {
        if (bStateChanged)
        {
            if (PAR(cr))
            {
                gmx_sumd(2*DIM, mu, cr);
                ddReopenBalanceRegionCpu(cr->dd);
            }
            for (i = 0; i < 2; i++)
            {
                for (j = 0; j < DIM; j++)
                {
                    fr->mu_tot[i][j] = mu[i*DIM + j];
                }
            }
        }
        if (fr->efep == efepNO)
        {
            copy_rvec(fr->mu_tot[0], mu_tot);
        }
        else
        {
            for (j = 0; j < DIM; j++)
            {
                mu_tot[j] =
                    (1.0 - lambda[efptCOUL])*fr->mu_tot[0][j] + lambda[efptCOUL]*fr->mu_tot[1][j];
            }
        }
    }

    /* Reset energies */
    reset_enerdata(enerd);
    clear_rvecs(SHIFTS, fr->fshift);

    if (bNS)
    {
        wallcycle_start(wcycle, ewcNS);

        if (graph && bStateChanged)
        {
            /* Calculate intramolecular shift vectors to make molecules whole */
            mk_mshift(fplog, graph, fr->ePBC, box, as_rvec_array(x.data()));
        }

        /* Do the actual neighbour searching */
        ns(fplog, fr, box,
           groups, top, mdatoms,
           cr, nrnb, bFillGrid);

        wallcycle_stop(wcycle, ewcNS);
    }

    if (inputrec->implicit_solvent && bNS)
    {
        make_gb_nblist(cr, inputrec->gb_algorithm,
                       as_rvec_array(x.data()), box, fr, &top->idef, graph, fr->born);
    }

    if (DOMAINDECOMP(cr) && !thisRankHasDuty(cr, DUTY_PME))
    {
        wallcycle_start(wcycle, ewcPPDURINGPME);
        dd_force_flop_start(cr->dd, nrnb);
    }

    if (inputrec->bRot)
    {
        wallcycle_start(wcycle, ewcROT);
        do_rotation(cr, inputrec, box, as_rvec_array(x.data()), t, step, bNS);
        wallcycle_stop(wcycle, ewcROT);
    }

    /* Temporary solution until all routines take PaddedRVecVector */
    rvec *f = as_rvec_array(force.data());

    /* Start the force cycle counter.
     * Note that a different counter is used for dynamic load balancing.
     */
    wallcycle_start(wcycle, ewcFORCE);

    gmx::ArrayRef<gmx::RVec> forceRef = force;
    if (bDoForces)
    {
        /* If we need to compute the virial, we might need a separate
         * force buffer for algorithms for which the virial is calculated
         * separately, such as PME.
         */
        if ((flags & GMX_FORCE_VIRIAL) && fr->haveDirectVirialContributions)
        {
            forceRef = *fr->forceBufferForDirectVirialContributions;

            clear_rvecs_omp(forceRef.size(), as_rvec_array(forceRef.data()));
        }

        /* Clear the short- and long-range forces */
        clear_rvecs(fr->natoms_force_constr, f);
    }

    /* forceWithVirial might need the full force atom range */
    ForceWithVirial forceWithVirial(forceRef, flags & GMX_FORCE_VIRIAL);

    if (inputrec->bPull && pull_have_constraint(inputrec->pull_work))
    {
        clear_pull_forces(inputrec->pull_work);
    }

    /* update QMMMrec, if necessary */
    if (fr->bQMMM)
    {
        update_QMMMrec(cr, fr, as_rvec_array(x.data()), mdatoms, box);
    }

    /* Compute the bonded and non-bonded energies and optionally forces */
    do_force_lowlevel(fr, inputrec, &(top->idef),
                      cr, nrnb, wcycle, mdatoms,
                      as_rvec_array(x.data()), hist, f, &forceWithVirial, enerd, fcd, top, fr->born,
                      bBornRadii, box,
                      inputrec->fepvals, lambda,
                      graph, &(top->excls), fr->mu_tot,
                      flags,
                      &cycles_pme);

    wallcycle_stop(wcycle, ewcFORCE);

    if (DOMAINDECOMP(cr))
    {
        dd_force_flop_stop(cr->dd, nrnb);

        if (ddCloseBalanceRegion == DdCloseBalanceRegionAfterForceComputation::yes)
        {
            ddCloseBalanceRegionCpu(cr->dd);
        }
    }

    computeSpecialForces(fplog, cr, inputrec, step, t, wcycle,
                         fr->forceProviders, box, as_rvec_array(x.data()), mdatoms, lambda,
                         flags, &forceWithVirial, enerd,
                         ed, bNS);

    if (bDoForces)
    {
        /* Communicate the forces */
        if (DOMAINDECOMP(cr))
        {
            wallcycle_start(wcycle, ewcMOVEF);
            dd_move_f(cr->dd, f, fr->fshift);
            /* Do we need to communicate the separate force array
             * for terms that do not contribute to the single sum virial?
             * Position restraints and electric fields do not introduce
             * inter-cg forces, only full electrostatics methods do.
             * When we do not calculate the virial, fr->f_novirsum = f,
             * so we have already communicated these forces.
             */
            if (EEL_FULL(fr->ic->eeltype) && cr->dd->n_intercg_excl &&
                (flags & GMX_FORCE_VIRIAL))
            {
                dd_move_f(cr->dd, as_rvec_array(forceWithVirial.force_.data()), nullptr);
            }
            wallcycle_stop(wcycle, ewcMOVEF);
        }

        /* If we have NoVirSum forces, but we do not calculate the virial,
         * we sum fr->f_novirsum=f later.
         */
        if (vsite && !(fr->haveDirectVirialContributions && !(flags & GMX_FORCE_VIRIAL)))
        {
            wallcycle_start(wcycle, ewcVSITESPREAD);
            spread_vsite_f(vsite, as_rvec_array(x.data()), f, fr->fshift, FALSE, nullptr, nrnb,
                           &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
            wallcycle_stop(wcycle, ewcVSITESPREAD);
        }

        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Calculation of the virial must be done after vsites! */
            calc_virial(0, mdatoms->homenr, as_rvec_array(x.data()), f,
                        vir_force, graph, box, nrnb, fr, inputrec->ePBC);
        }
    }

    if (PAR(cr) && !thisRankHasDuty(cr, DUTY_PME))
    {
        /* In case of node-splitting, the PP nodes receive the long-range
         * forces, virial and energy from the PME nodes here.
         */
        pme_receive_force_ener(cr, &forceWithVirial, enerd, wcycle);
    }

    if (bDoForces)
    {
        post_process_forces(cr, step, nrnb, wcycle,
                            top, box, as_rvec_array(x.data()), f, &forceWithVirial,
                            vir_force, mdatoms, graph, fr, vsite,
                            flags);
    }

    if (flags & GMX_FORCE_ENERGY)
    {
        /* Sum the potential energy terms from group contributions */
        sum_epot(&(enerd->grpp), enerd->term);

        if (!EI_TPI(inputrec->eI))
        {
            checkPotentialEnergyValidity(step, *enerd, *inputrec);
        }
    }

}

void do_force(FILE *fplog, t_commrec *cr,
              t_inputrec *inputrec,
              gmx_int64_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
              gmx_localtop_t *top,
              gmx_groups_t *groups,
              matrix box, gmx::PaddedArrayRef<gmx::RVec> x, history_t *hist,
              gmx::PaddedArrayRef<gmx::RVec> force,
              tensor vir_force,
              t_mdatoms *mdatoms,
              gmx_enerdata_t *enerd, t_fcdata *fcd,
              gmx::ArrayRef<real> lambda, t_graph *graph,
              t_forcerec *fr,
              gmx_vsite_t *vsite, rvec mu_tot,
              double t, gmx_edsam_t ed,
              gmx_bool bBornRadii,
              int flags,
              DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
              DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion)
{
    /* modify force flag if not doing nonbonded */
    if (!fr->bNonbonded)
    {
        flags &= ~GMX_FORCE_NONBONDED;
    }

    GMX_ASSERT(x.size()     >= gmx::paddedRVecVectorSize(fr->natoms_force), "coordinates should be padded");
    GMX_ASSERT(force.size() >= gmx::paddedRVecVectorSize(fr->natoms_force), "force should be padded");

    switch (inputrec->cutoff_scheme)
    {
        case ecutsVERLET:
            do_force_cutsVERLET(fplog, cr, inputrec,
                                step, nrnb, wcycle,
                                top,
                                groups,
                                box, x, hist,
                                force, vir_force,
                                mdatoms,
                                enerd, fcd,
                                lambda.data(), graph,
                                fr, fr->ic,
                                vsite, mu_tot,
                                t, ed,
                                bBornRadii,
                                flags,
                                ddOpenBalanceRegion,
                                ddCloseBalanceRegion);
            break;
        case ecutsGROUP:
            do_force_cutsGROUP(fplog, cr, inputrec,
                               step, nrnb, wcycle,
                               top,
                               groups,
                               box, x, hist,
                               force, vir_force,
                               mdatoms,
                               enerd, fcd,
                               lambda.data(), graph,
                               fr, vsite, mu_tot,
                               t, ed,
                               bBornRadii,
                               flags,
                               ddOpenBalanceRegion,
                               ddCloseBalanceRegion);
            break;
        default:
            gmx_incons("Invalid cut-off scheme passed!");
    }

    /* In case we don't have constraints and are using GPUs, the next balancing
     * region starts here.
     * Some "special" work at the end of do_force_cuts?, such as vsite spread,
     * virial calculation and COM pulling, is not thus not included in
     * the balance timing, which is ok as most tasks do communication.
     */
    if (ddOpenBalanceRegion == DdOpenBalanceRegionBeforeForceComputation::yes)
    {
        ddOpenBalanceRegionCpu(cr->dd, DdAllowBalanceRegionReopen::no);
    }
}


void do_constrain_first(FILE *fplog, gmx_constr_t constr,
                        t_inputrec *ir, t_mdatoms *md,
                        t_state *state, t_commrec *cr, t_nrnb *nrnb,
                        t_forcerec *fr, gmx_localtop_t *top)
{
    int             i, m, start, end;
    gmx_int64_t     step;
    real            dt = ir->delta_t;
    real            dvdl_dum;
    rvec           *savex;

    /* We need to allocate one element extra, since we might use
     * (unaligned) 4-wide SIMD loads to access rvec entries.
     */
    snew(savex, state->natoms + 1);

    start = 0;
    end   = md->homenr;

    if (debug)
    {
        fprintf(debug, "vcm: start=%d, homenr=%d, end=%d\n",
                start, md->homenr, end);
    }
    /* Do a first constrain to reset particles... */
    step = ir->init_step;
    if (fplog)
    {
        char buf[STEPSTRSIZE];
        fprintf(fplog, "\nConstraining the starting coordinates (step %s)\n",
                gmx_step_str(step, buf));
    }
    dvdl_dum = 0;

    /* constrain the current position */
    constrain(nullptr, TRUE, FALSE, constr, &(top->idef),
              ir, cr, step, 0, 1.0, md,
              as_rvec_array(state->x.data()), as_rvec_array(state->x.data()), nullptr,
              fr->bMolPBC, state->box,
              state->lambda[efptBONDED], &dvdl_dum,
              nullptr, nullptr, nrnb, econqCoord);
    if (EI_VV(ir->eI))
    {
        /* constrain the inital velocity, and save it */
        /* also may be useful if we need the ekin from the halfstep for velocity verlet */
        constrain(nullptr, TRUE, FALSE, constr, &(top->idef),
                  ir, cr, step, 0, 1.0, md,
                  as_rvec_array(state->x.data()), as_rvec_array(state->v.data()), as_rvec_array(state->v.data()),
                  fr->bMolPBC, state->box,
                  state->lambda[efptBONDED], &dvdl_dum,
                  nullptr, nullptr, nrnb, econqVeloc);
    }
    /* constrain the inital velocities at t-dt/2 */
    if (EI_STATE_VELOCITY(ir->eI) && ir->eI != eiVV)
    {
        for (i = start; (i < end); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                /* Reverse the velocity */
                state->v[i][m] = -state->v[i][m];
                /* Store the position at t-dt in buf */
                savex[i][m] = state->x[i][m] + dt*state->v[i][m];
            }
        }
        /* Shake the positions at t=-dt with the positions at t=0
         * as reference coordinates.
         */
        if (fplog)
        {
            char buf[STEPSTRSIZE];
            fprintf(fplog, "\nConstraining the coordinates at t0-dt (step %s)\n",
                    gmx_step_str(step, buf));
        }
        dvdl_dum = 0;
        constrain(nullptr, TRUE, FALSE, constr, &(top->idef),
                  ir, cr, step, -1, 1.0, md,
                  as_rvec_array(state->x.data()), savex, nullptr,
                  fr->bMolPBC, state->box,
                  state->lambda[efptBONDED], &dvdl_dum,
                  as_rvec_array(state->v.data()), nullptr, nrnb, econqCoord);

        for (i = start; i < end; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                /* Re-reverse the velocities */
                state->v[i][m] = -state->v[i][m];
            }
        }
    }
    sfree(savex);
}


static void
integrate_table(real vdwtab[], real scale, int offstart, int rstart, int rend,
                double *enerout, double *virout)
{
    double enersum, virsum;
    double invscale, invscale2, invscale3;
    double r, ea, eb, ec, pa, pb, pc, pd;
    double y0, f, g, h;
    int    ri, offset;
    double tabfactor;

    invscale  = 1.0/scale;
    invscale2 = invscale*invscale;
    invscale3 = invscale*invscale2;

    /* Following summation derived from cubic spline definition,
     * Numerical Recipies in C, second edition, p. 113-116.  Exact for
     * the cubic spline.  We first calculate the negative of the
     * energy from rvdw to rvdw_switch, assuming that g(r)=1, and then
     * add the more standard, abrupt cutoff correction to that result,
     * yielding the long-range correction for a switched function.  We
     * perform both the pressure and energy loops at the same time for
     * simplicity, as the computational cost is low. */

    if (offstart == 0)
    {
        /* Since the dispersion table has been scaled down a factor
         * 6.0 and the repulsion a factor 12.0 to compensate for the
         * c6/c12 parameters inside nbfp[] being scaled up (to save
         * flops in kernels), we need to correct for this.
         */
        tabfactor = 6.0;
    }
    else
    {
        tabfactor = 12.0;
    }

    enersum = 0.0;
    virsum  = 0.0;
    for (ri = rstart; ri < rend; ++ri)
    {
        r  = ri*invscale;
        ea = invscale3;
        eb = 2.0*invscale2*r;
        ec = invscale*r*r;

        pa = invscale3;
        pb = 3.0*invscale2*r;
        pc = 3.0*invscale*r*r;
        pd = r*r*r;

        /* this "8" is from the packing in the vdwtab array - perhaps
           should be defined? */

        offset = 8*ri + offstart;
        y0     = vdwtab[offset];
        f      = vdwtab[offset+1];
        g      = vdwtab[offset+2];
        h      = vdwtab[offset+3];

        enersum += y0*(ea/3 + eb/2 + ec) + f*(ea/4 + eb/3 + ec/2) + g*(ea/5 + eb/4 + ec/3) + h*(ea/6 + eb/5 + ec/4);
        virsum  +=  f*(pa/4 + pb/3 + pc/2 + pd) + 2*g*(pa/5 + pb/4 + pc/3 + pd/2) + 3*h*(pa/6 + pb/5 + pc/4 + pd/3);
    }
    *enerout = 4.0*M_PI*enersum*tabfactor;
    *virout  = 4.0*M_PI*virsum*tabfactor;
}

void calc_enervirdiff(FILE *fplog, int eDispCorr, t_forcerec *fr)
{
    double   eners[2], virs[2], enersum, virsum;
    double   r0, rc3, rc9;
    int      ri0, ri1, i;
    real     scale, *vdwtab;

    fr->enershiftsix    = 0;
    fr->enershifttwelve = 0;
    fr->enerdiffsix     = 0;
    fr->enerdifftwelve  = 0;
    fr->virdiffsix      = 0;
    fr->virdifftwelve   = 0;

    const interaction_const_t *ic = fr->ic;

    if (eDispCorr != edispcNO)
    {
        for (i = 0; i < 2; i++)
        {
            eners[i] = 0;
            virs[i]  = 0;
        }
        if ((ic->vdw_modifier == eintmodPOTSHIFT) ||
            (ic->vdw_modifier == eintmodPOTSWITCH) ||
            (ic->vdw_modifier == eintmodFORCESWITCH) ||
            (ic->vdwtype == evdwSHIFT) ||
            (ic->vdwtype == evdwSWITCH))
        {
            if (((ic->vdw_modifier == eintmodPOTSWITCH) ||
                 (ic->vdw_modifier == eintmodFORCESWITCH) ||
                 (ic->vdwtype == evdwSWITCH)) && ic->rvdw_switch == 0)
            {
                gmx_fatal(FARGS,
                          "With dispersion correction rvdw-switch can not be zero "
                          "for vdw-type = %s", evdw_names[ic->vdwtype]);
            }

            /* TODO This code depends on the logic in tables.c that
               constructs the table layout, which should be made
               explicit in future cleanup. */
            GMX_ASSERT(fr->dispersionCorrectionTable->interaction == GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
                       "Dispersion-correction code needs a table with both repulsion and dispersion terms");
            scale  = fr->dispersionCorrectionTable->scale;
            vdwtab = fr->dispersionCorrectionTable->data;

            /* Round the cut-offs to exact table values for precision */
            ri0  = static_cast<int>(floor(ic->rvdw_switch*scale));
            ri1  = static_cast<int>(ceil(ic->rvdw*scale));

            /* The code below has some support for handling force-switching, i.e.
             * when the force (instead of potential) is switched over a limited
             * region. This leads to a constant shift in the potential inside the
             * switching region, which we can handle by adding a constant energy
             * term in the force-switch case just like when we do potential-shift.
             *
             * For now this is not enabled, but to keep the functionality in the
             * code we check separately for switch and shift. When we do force-switch
             * the shifting point is rvdw_switch, while it is the cutoff when we
             * have a classical potential-shift.
             *
             * For a pure potential-shift the potential has a constant shift
             * all the way out to the cutoff, and that is it. For other forms
             * we need to calculate the constant shift up to the point where we
             * start modifying the potential.
             */
            ri0  = (ic->vdw_modifier == eintmodPOTSHIFT) ? ri1 : ri0;

            r0   = ri0/scale;
            rc3  = r0*r0*r0;
            rc9  = rc3*rc3*rc3;

            if ((ic->vdw_modifier == eintmodFORCESWITCH) ||
                (ic->vdwtype == evdwSHIFT))
            {
                /* Determine the constant energy shift below rvdw_switch.
                 * Table has a scale factor since we have scaled it down to compensate
                 * for scaling-up c6/c12 with the derivative factors to save flops in analytical kernels.
                 */
                fr->enershiftsix    = (real)(-1.0/(rc3*rc3)) - 6.0*vdwtab[8*ri0];
                fr->enershifttwelve = (real)( 1.0/(rc9*rc3)) - 12.0*vdwtab[8*ri0 + 4];
            }
            else if (ic->vdw_modifier == eintmodPOTSHIFT)
            {
                fr->enershiftsix    = (real)(-1.0/(rc3*rc3));
                fr->enershifttwelve = (real)( 1.0/(rc9*rc3));
            }

            /* Add the constant part from 0 to rvdw_switch.
             * This integration from 0 to rvdw_switch overcounts the number
             * of interactions by 1, as it also counts the self interaction.
             * We will correct for this later.
             */
            eners[0] += 4.0*M_PI*fr->enershiftsix*rc3/3.0;
            eners[1] += 4.0*M_PI*fr->enershifttwelve*rc3/3.0;

            /* Calculate the contribution in the range [r0,r1] where we
             * modify the potential. For a pure potential-shift modifier we will
             * have ri0==ri1, and there will not be any contribution here.
             */
            for (i = 0; i < 2; i++)
            {
                enersum = 0;
                virsum  = 0;
                integrate_table(vdwtab, scale, (i == 0 ? 0 : 4), ri0, ri1, &enersum, &virsum);
                eners[i] -= enersum;
                virs[i]  -= virsum;
            }

            /* Alright: Above we compensated by REMOVING the parts outside r0
             * corresponding to the ideal VdW 1/r6 and /r12 potentials.
             *
             * Regardless of whether r0 is the point where we start switching,
             * or the cutoff where we calculated the constant shift, we include
             * all the parts we are missing out to infinity from r0 by
             * calculating the analytical dispersion correction.
             */
            eners[0] += -4.0*M_PI/(3.0*rc3);
            eners[1] +=  4.0*M_PI/(9.0*rc9);
            virs[0]  +=  8.0*M_PI/rc3;
            virs[1]  += -16.0*M_PI/(3.0*rc9);
        }
        else if (ic->vdwtype == evdwCUT ||
                 EVDW_PME(ic->vdwtype) ||
                 ic->vdwtype == evdwUSER)
        {
            if (ic->vdwtype == evdwUSER && fplog)
            {
                fprintf(fplog,
                        "WARNING: using dispersion correction with user tables\n");
            }

            /* Note that with LJ-PME, the dispersion correction is multiplied
             * by the difference between the actual C6 and the value of C6
             * that would produce the combination rule.
             * This means the normal energy and virial difference formulas
             * can be used here.
             */

            rc3  = ic->rvdw*ic->rvdw*ic->rvdw;
            rc9  = rc3*rc3*rc3;
            /* Contribution beyond the cut-off */
            eners[0] += -4.0*M_PI/(3.0*rc3);
            eners[1] +=  4.0*M_PI/(9.0*rc9);
            if (ic->vdw_modifier == eintmodPOTSHIFT)
            {
                /* Contribution within the cut-off */
                eners[0] += -4.0*M_PI/(3.0*rc3);
                eners[1] +=  4.0*M_PI/(3.0*rc9);
            }
            /* Contribution beyond the cut-off */
            virs[0]  +=  8.0*M_PI/rc3;
            virs[1]  += -16.0*M_PI/(3.0*rc9);
        }
        else
        {
            gmx_fatal(FARGS,
                      "Dispersion correction is not implemented for vdw-type = %s",
                      evdw_names[ic->vdwtype]);
        }

        /* When we deprecate the group kernels the code below can go too */
        if (ic->vdwtype == evdwPME && fr->cutoff_scheme == ecutsGROUP)
        {
            /* Calculate self-interaction coefficient (assuming that
             * the reciprocal-space contribution is constant in the
             * region that contributes to the self-interaction).
             */
            fr->enershiftsix = gmx::power6(ic->ewaldcoeff_lj) / 6.0;

            eners[0] += -gmx::power3(std::sqrt(M_PI)*ic->ewaldcoeff_lj)/3.0;
            virs[0]  +=  gmx::power3(std::sqrt(M_PI)*ic->ewaldcoeff_lj);
        }

        fr->enerdiffsix    = eners[0];
        fr->enerdifftwelve = eners[1];
        /* The 0.5 is due to the Gromacs definition of the virial */
        fr->virdiffsix     = 0.5*virs[0];
        fr->virdifftwelve  = 0.5*virs[1];
    }
}

void calc_dispcorr(t_inputrec *ir, t_forcerec *fr,
                   matrix box, real lambda, tensor pres, tensor virial,
                   real *prescorr, real *enercorr, real *dvdlcorr)
{
    gmx_bool bCorrAll, bCorrPres;
    real     dvdlambda, invvol, dens, ninter, avcsix, avctwelve, enerdiff, svir = 0, spres = 0;
    int      m;

    *prescorr = 0;
    *enercorr = 0;
    *dvdlcorr = 0;

    clear_mat(virial);
    clear_mat(pres);

    if (ir->eDispCorr != edispcNO)
    {
        bCorrAll  = (ir->eDispCorr == edispcAllEner ||
                     ir->eDispCorr == edispcAllEnerPres);
        bCorrPres = (ir->eDispCorr == edispcEnerPres ||
                     ir->eDispCorr == edispcAllEnerPres);

        invvol = 1/det(box);
        if (fr->n_tpi)
        {
            /* Only correct for the interactions with the inserted molecule */
            dens   = (fr->numAtomsForDispersionCorrection - fr->n_tpi)*invvol;
            ninter = fr->n_tpi;
        }
        else
        {
            dens   = fr->numAtomsForDispersionCorrection*invvol;
            ninter = 0.5*fr->numAtomsForDispersionCorrection;
        }

        if (ir->efep == efepNO)
        {
            avcsix    = fr->avcsix[0];
            avctwelve = fr->avctwelve[0];
        }
        else
        {
            avcsix    = (1 - lambda)*fr->avcsix[0]    + lambda*fr->avcsix[1];
            avctwelve = (1 - lambda)*fr->avctwelve[0] + lambda*fr->avctwelve[1];
        }

        enerdiff   = ninter*(dens*fr->enerdiffsix - fr->enershiftsix);
        *enercorr += avcsix*enerdiff;
        dvdlambda  = 0.0;
        if (ir->efep != efepNO)
        {
            dvdlambda += (fr->avcsix[1] - fr->avcsix[0])*enerdiff;
        }
        if (bCorrAll)
        {
            enerdiff   = ninter*(dens*fr->enerdifftwelve - fr->enershifttwelve);
            *enercorr += avctwelve*enerdiff;
            if (fr->efep != efepNO)
            {
                dvdlambda += (fr->avctwelve[1] - fr->avctwelve[0])*enerdiff;
            }
        }

        if (bCorrPres)
        {
            svir = ninter*dens*avcsix*fr->virdiffsix/3.0;
            if (ir->eDispCorr == edispcAllEnerPres)
            {
                svir += ninter*dens*avctwelve*fr->virdifftwelve/3.0;
            }
            /* The factor 2 is because of the Gromacs virial definition */
            spres = -2.0*invvol*svir*PRESFAC;

            for (m = 0; m < DIM; m++)
            {
                virial[m][m] += svir;
                pres[m][m]   += spres;
            }
            *prescorr += spres;
        }

        /* Can't currently control when it prints, for now, just print when degugging */
        if (debug)
        {
            if (bCorrAll)
            {
                fprintf(debug, "Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
                        avcsix, avctwelve);
            }
            if (bCorrPres)
            {
                fprintf(debug,
                        "Long Range LJ corr.: Epot %10g, Pres: %10g, Vir: %10g\n",
                        *enercorr, spres, svir);
            }
            else
            {
                fprintf(debug, "Long Range LJ corr.: Epot %10g\n", *enercorr);
            }
        }

        if (fr->efep != efepNO)
        {
            *dvdlcorr += dvdlambda;
        }
    }
}

static void low_do_pbc_mtop(FILE *fplog, int ePBC, const matrix box,
                            const gmx_mtop_t *mtop, rvec x[],
                            gmx_bool bFirst)
{
    t_graph        *graph;
    int             mb, as, mol;
    gmx_molblock_t *molb;

    if (bFirst && fplog)
    {
        fprintf(fplog, "Removing pbc first time\n");
    }

    snew(graph, 1);
    as = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        if (molb->natoms_mol == 1 ||
            (!bFirst && mtop->moltype[molb->type].cgs.nr == 1))
        {
            /* Just one atom or charge group in the molecule, no PBC required */
            as += molb->nmol*molb->natoms_mol;
        }
        else
        {
            /* Pass NULL iso fplog to avoid graph prints for each molecule type */
            mk_graph_ilist(nullptr, mtop->moltype[molb->type].ilist,
                           0, molb->natoms_mol, FALSE, FALSE, graph);

            for (mol = 0; mol < molb->nmol; mol++)
            {
                mk_mshift(fplog, graph, ePBC, box, x+as);

                shift_self(graph, box, x+as);
                /* The molecule is whole now.
                 * We don't need the second mk_mshift call as in do_pbc_first,
                 * since we no longer need this graph.
                 */

                as += molb->natoms_mol;
            }
            done_graph(graph);
        }
    }
    sfree(graph);
}

void do_pbc_first_mtop(FILE *fplog, int ePBC, const matrix box,
                       const gmx_mtop_t *mtop, rvec x[])
{
    low_do_pbc_mtop(fplog, ePBC, box, mtop, x, TRUE);
}

void do_pbc_mtop(FILE *fplog, int ePBC, const matrix box,
                 const gmx_mtop_t *mtop, rvec x[])
{
    low_do_pbc_mtop(fplog, ePBC, box, mtop, x, FALSE);
}

void put_atoms_in_box_omp(int ePBC, const matrix box, gmx::ArrayRef<gmx::RVec> x)
{
    int t, nth;
    nth = gmx_omp_nthreads_get(emntDefault);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (t = 0; t < nth; t++)
    {
        try
        {
            size_t natoms = x.size();
            size_t offset = (natoms*t    )/nth;
            size_t len    = (natoms*(t + 1))/nth - offset;
            put_atoms_in_box(ePBC, box, x.subArray(offset, len));
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

// TODO This can be cleaned up a lot, and move back to runner.cpp
void finish_run(FILE *fplog, const gmx::MDLogger &mdlog, t_commrec *cr,
                const t_inputrec *inputrec,
                t_nrnb nrnb[], gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                nonbonded_verlet_t *nbv,
                gmx_pme_t *pme,
                gmx_bool bWriteStat)
{
    t_nrnb *nrnb_tot = nullptr;
    double  delta_t  = 0;
    double  nbfs     = 0, mflop = 0;
    double  elapsed_time,
            elapsed_time_over_all_ranks,
            elapsed_time_over_all_threads,
            elapsed_time_over_all_threads_over_all_ranks;
    /* Control whether it is valid to print a report. Only the
       simulation master may print, but it should not do so if the run
       terminated e.g. before a scheduled reset step. This is
       complicated by the fact that PME ranks are unaware of the
       reason why they were sent a pmerecvqxFINISH. To avoid
       communication deadlocks, we always do the communication for the
       report, even if we've decided not to write the report, because
       how long it takes to finish the run is not important when we've
       decided not to report on the simulation performance. */
    bool    printReport = SIMMASTER(cr);

    if (!walltime_accounting_get_valid_finish(walltime_accounting))
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText("Simulation ended prematurely, no performance report will be written.");
        printReport = false;
    }

    if (cr->nnodes > 1)
    {
        snew(nrnb_tot, 1);
#if GMX_MPI
        MPI_Allreduce(nrnb->n, nrnb_tot->n, eNRNB, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mysim);
#endif
    }
    else
    {
        nrnb_tot = nrnb;
    }

    elapsed_time                                 = walltime_accounting_get_elapsed_time(walltime_accounting);
    elapsed_time_over_all_ranks                  = elapsed_time;
    elapsed_time_over_all_threads                = walltime_accounting_get_elapsed_time_over_all_threads(walltime_accounting);
    elapsed_time_over_all_threads_over_all_ranks = elapsed_time_over_all_threads;
#if GMX_MPI
    if (cr->nnodes > 1)
    {
        /* reduce elapsed_time over all MPI ranks in the current simulation */
        MPI_Allreduce(&elapsed_time,
                      &elapsed_time_over_all_ranks,
                      1, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mysim);
        elapsed_time_over_all_ranks /= cr->nnodes;
        /* Reduce elapsed_time_over_all_threads over all MPI ranks in the
         * current simulation. */
        MPI_Allreduce(&elapsed_time_over_all_threads,
                      &elapsed_time_over_all_threads_over_all_ranks,
                      1, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mysim);
    }
#endif

    if (printReport)
    {
        print_flop(fplog, nrnb_tot, &nbfs, &mflop);
    }
    if (cr->nnodes > 1)
    {
        sfree(nrnb_tot);
    }

    if (thisRankHasDuty(cr, DUTY_PP) && DOMAINDECOMP(cr))
    {
        print_dd_statistics(cr, inputrec, fplog);
    }

    /* TODO Move the responsibility for any scaling by thread counts
     * to the code that handled the thread region, so that there's a
     * mechanism to keep cycle counting working during the transition
     * to task parallelism. */
    int nthreads_pp  = gmx_omp_nthreads_get(emntNonbonded);
    int nthreads_pme = gmx_omp_nthreads_get(emntPME);
    wallcycle_scale_by_num_threads(wcycle, thisRankHasDuty(cr, DUTY_PME) && !thisRankHasDuty(cr, DUTY_PP), nthreads_pp, nthreads_pme);
    auto cycle_sum(wallcycle_sum(cr, wcycle));

    if (printReport)
    {
        auto                    nbnxn_gpu_timings = use_GPU(nbv) ? nbnxn_gpu_get_timings(nbv->gpu_nbv) : nullptr;
        gmx_wallclock_gpu_pme_t pme_gpu_timings   = {};
        if (pme_gpu_task_enabled(pme))
        {
            pme_gpu_get_timings(pme, &pme_gpu_timings);
        }
        wallcycle_print(fplog, mdlog, cr->nnodes, cr->npmenodes, nthreads_pp, nthreads_pme,
                        elapsed_time_over_all_ranks,
                        wcycle, cycle_sum,
                        nbnxn_gpu_timings,
                        &pme_gpu_timings);

        if (EI_DYNAMICS(inputrec->eI))
        {
            delta_t = inputrec->delta_t;
        }

        if (fplog)
        {
            print_perf(fplog, elapsed_time_over_all_threads_over_all_ranks,
                       elapsed_time_over_all_ranks,
                       walltime_accounting_get_nsteps_done(walltime_accounting),
                       delta_t, nbfs, mflop);
        }
        if (bWriteStat)
        {
            print_perf(stderr, elapsed_time_over_all_threads_over_all_ranks,
                       elapsed_time_over_all_ranks,
                       walltime_accounting_get_nsteps_done(walltime_accounting),
                       delta_t, nbfs, mflop);
        }
    }
}

extern void initialize_lambdas(FILE *fplog, t_inputrec *ir, int *fep_state, gmx::ArrayRef<real> lambda, double *lam0)
{
    /* this function works, but could probably use a logic rewrite to keep all the different
       types of efep straight. */

    if ((ir->efep == efepNO) && (ir->bSimTemp == FALSE))
    {
        return;
    }

    t_lambda *fep = ir->fepvals;
    *fep_state    = fep->init_fep_state; /* this might overwrite the checkpoint
                                            if checkpoint is set -- a kludge is in for now
                                            to prevent this.*/

    for (int i = 0; i < efptNR; i++)
    {
        /* overwrite lambda state with init_lambda for now for backwards compatibility */
        if (fep->init_lambda >= 0) /* if it's -1, it was never initializd */
        {
            lambda[i] = fep->init_lambda;
            if (lam0)
            {
                lam0[i] = lambda[i];
            }
        }
        else
        {
            lambda[i] = fep->all_lambda[i][*fep_state];
            if (lam0)
            {
                lam0[i] = lambda[i];
            }
        }
    }
    if (ir->bSimTemp)
    {
        /* need to rescale control temperatures to match current state */
        for (int i = 0; i < ir->opts.ngtc; i++)
        {
            if (ir->opts.ref_t[i] > 0)
            {
                ir->opts.ref_t[i] = ir->simtempvals->temperatures[*fep_state];
            }
        }
    }

    /* Send to the log the information on the current lambdas */
    if (fplog != nullptr)
    {
        fprintf(fplog, "Initial vector of lambda components:[ ");
        for (const auto &l : lambda)
        {
            fprintf(fplog, "%10.4f ", l);
        }
        fprintf(fplog, "]\n");
    }
    return;
}


void init_md(FILE *fplog,
             t_commrec *cr, gmx::IMDOutputProvider *outputProvider,
             t_inputrec *ir, const gmx_output_env_t *oenv,
             const MdrunOptions &mdrunOptions,
             double *t, double *t0,
             t_state *globalState, double *lam0,
             t_nrnb *nrnb, gmx_mtop_t *mtop,
             gmx_update_t **upd,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             tensor force_vir, tensor shake_vir, rvec mu_tot,
             gmx_bool *bSimAnn, t_vcm **vcm,
             gmx_wallcycle_t wcycle)
{
    int  i;

    /* Initial values */
    *t = *t0       = ir->init_t;

    *bSimAnn = FALSE;
    for (i = 0; i < ir->opts.ngtc; i++)
    {
        /* set bSimAnn if any group is being annealed */
        if (ir->opts.annealing[i] != eannNO)
        {
            *bSimAnn = TRUE;
        }
    }

    /* Initialize lambda variables */
    /* TODO: Clean up initialization of fep_state and lambda in t_state.
     * We currently need to call initialize_lambdas on non-master ranks
     * to initialize lam0.
     */
    if (MASTER(cr))
    {
        initialize_lambdas(fplog, ir, &globalState->fep_state, globalState->lambda, lam0);
    }
    else
    {
        int                      tmpFepState;
        std::array<real, efptNR> tmpLambda;
        initialize_lambdas(fplog, ir, &tmpFepState, tmpLambda, lam0);
    }

    // TODO upd is never NULL in practice, but the analysers don't know that
    if (upd)
    {
        *upd = init_update(ir);
    }
    if (*bSimAnn)
    {
        update_annealing_target_temp(ir, ir->init_t, upd ? *upd : nullptr);
    }

    if (vcm != nullptr)
    {
        *vcm = init_vcm(fplog, &mtop->groups, ir);
    }

    if (EI_DYNAMICS(ir->eI) && !mdrunOptions.continuationOptions.appendFiles)
    {
        if (ir->etc == etcBERENDSEN)
        {
            please_cite(fplog, "Berendsen84a");
        }
        if (ir->etc == etcVRESCALE)
        {
            please_cite(fplog, "Bussi2007a");
        }
        if (ir->eI == eiSD1)
        {
            please_cite(fplog, "Goga2012");
        }
    }
    init_nrnb(nrnb);

    if (nfile != -1)
    {
        *outf = init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider, ir, mtop, oenv, wcycle);

        *mdebin = init_mdebin(mdrunOptions.continuationOptions.appendFiles ? nullptr : mdoutf_get_fp_ene(*outf),
                              mtop, ir, mdoutf_get_fp_dhdl(*outf));
    }

    /* Initiate variables */
    clear_mat(force_vir);
    clear_mat(shake_vir);
    clear_rvec(mu_tot);
}
