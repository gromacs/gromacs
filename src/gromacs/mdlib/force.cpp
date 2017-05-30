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

#include "force.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/ewald.h"
#include "gromacs/ewald/long-range-correction.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/listed-forces/listed-forces.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/forcerec-threading.h"
#include "gromacs/mdlib/genborn.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void ns(FILE              *fp,
        t_forcerec        *fr,
        matrix             box,
        gmx_groups_t      *groups,
        gmx_localtop_t    *top,
        t_mdatoms         *md,
        t_commrec         *cr,
        t_nrnb            *nrnb,
        gmx_bool           bFillGrid)
{
    int     nsearch;


    if (!fr->ns->nblist_initialized)
    {
        init_neighbor_list(fp, fr, md->homenr);
    }

    nsearch = search_neighbours(fp, fr, box, top, groups, cr, nrnb, md,
                                bFillGrid);
    if (debug)
    {
        fprintf(debug, "nsearch = %d\n", nsearch);
    }

    /* Check whether we have to do dynamic load balancing */
    /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
       count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
       &(top->idef),opts->ngener);
     */
    if (fr->ns->dump_nl > 0)
    {
        dump_nblist(fp, cr, fr, fr->ns->dump_nl);
    }
}

static void clearEwaldThreadOutput(ewald_corr_thread_t *ewc_t)
{
    ewc_t->Vcorr_q        = 0;
    ewc_t->Vcorr_lj       = 0;
    ewc_t->dvdl[efptCOUL] = 0;
    ewc_t->dvdl[efptVDW]  = 0;
    clear_mat(ewc_t->vir_q);
    clear_mat(ewc_t->vir_lj);
}

static void reduceEwaldThreadOuput(int nthreads, ewald_corr_thread_t *ewc_t)
{
    ewald_corr_thread_t &dest = ewc_t[0];

    for (int t = 1; t < nthreads; t++)
    {
        dest.Vcorr_q        += ewc_t[t].Vcorr_q;
        dest.Vcorr_lj       += ewc_t[t].Vcorr_lj;
        dest.dvdl[efptCOUL] += ewc_t[t].dvdl[efptCOUL];
        dest.dvdl[efptVDW]  += ewc_t[t].dvdl[efptVDW];
        m_add(dest.vir_q,  ewc_t[t].vir_q,  dest.vir_q);
        m_add(dest.vir_lj, ewc_t[t].vir_lj, dest.vir_lj);
    }
}

void do_force_lowlevel(t_forcerec *fr,      t_inputrec *ir,
                       t_idef     *idef,    t_commrec  *cr,
                       t_nrnb     *nrnb,    gmx_wallcycle_t wcycle,
                       t_mdatoms  *md,
                       rvec       x[],      history_t  *hist,
                       rvec      *forceForUseWithShiftForces,
                       gmx::ForceWithVirial *forceWithVirial,
                       gmx_enerdata_t *enerd,
                       t_fcdata   *fcd,
                       gmx_localtop_t *top,
                       gmx_genborn_t *born,
                       gmx_bool       bBornRadii,
                       matrix     box,
                       t_lambda   *fepvals,
                       real       *lambda,
                       t_graph    *graph,
                       t_blocka   *excl,
                       rvec       mu_tot[],
                       int        flags,
                       float     *cycles_pme)
{
    int         i, j;
    int         donb_flags;
    int         pme_flags;
    t_pbc       pbc;
    real        dvdl_dum[efptNR], dvdl_nb[efptNR];

#if GMX_MPI
    double  t0 = 0.0, t1, t2, t3; /* time measurement for coarse load balancing */
#endif

    set_pbc(&pbc, fr->ePBC, box);

    /* reset free energy components */
    for (i = 0; i < efptNR; i++)
    {
        dvdl_nb[i]  = 0;
        dvdl_dum[i] = 0;
    }

    /* do QMMM first if requested */
    if (fr->bQMMM)
    {
        enerd->term[F_EQM] = calculate_QMMM(cr, forceForUseWithShiftForces, fr);
    }

    /* Call the short range functions all in one go. */

#if GMX_MPI
    /*#define TAKETIME ((cr->npmenodes) && (fr->timesteps < 12))*/
#define TAKETIME FALSE
    if (TAKETIME)
    {
        MPI_Barrier(cr->mpi_comm_mygroup);
        t0 = MPI_Wtime();
    }
#endif

    if (ir->nwall)
    {
        /* foreign lambda component for walls */
        real dvdl_walls = do_walls(ir, fr, box, md, x, forceForUseWithShiftForces, lambda[efptVDW],
                                   enerd->grpp.ener[egLJSR], nrnb);
        enerd->dvdl_lin[efptVDW] += dvdl_walls;
    }

    /* If doing GB, reset dvda and calculate the Born radii */
    if (ir->implicit_solvent)
    {
        wallcycle_sub_start(wcycle, ewcsNONBONDED);

        for (i = 0; i < born->nr; i++)
        {
            fr->dvda[i] = 0;
        }

        if (bBornRadii)
        {
            calc_gb_rad(cr, fr, ir, top, x, fr->gblist, born, md, nrnb);
        }

        wallcycle_sub_stop(wcycle, ewcsNONBONDED);
    }

    where();
    /* We only do non-bonded calculation with group scheme here, the verlet
     * calls are done from do_force_cutsVERLET(). */
    if (fr->cutoff_scheme == ecutsGROUP && (flags & GMX_FORCE_NONBONDED))
    {
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

        wallcycle_sub_start(wcycle, ewcsNONBONDED);
        do_nonbonded(fr, x, forceForUseWithShiftForces, md, excl,
                     &enerd->grpp, nrnb,
                     lambda, dvdl_nb, -1, -1, donb_flags);

        /* If we do foreign lambda and we have soft-core interactions
         * we have to recalculate the (non-linear) energies contributions.
         */
        if (fepvals->n_lambda > 0 && (flags & GMX_FORCE_DHDL) && fepvals->sc_alpha != 0)
        {
            for (i = 0; i < enerd->n_lambda; i++)
            {
                real lam_i[efptNR];

                for (j = 0; j < efptNR; j++)
                {
                    lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
                }
                reset_foreign_enerdata(enerd);
                do_nonbonded(fr, x, forceForUseWithShiftForces, md, excl,
                             &(enerd->foreign_grpp), nrnb,
                             lam_i, dvdl_dum, -1, -1,
                             (donb_flags & ~GMX_NONBONDED_DO_FORCE) | GMX_NONBONDED_DO_FOREIGNLAMBDA);
                sum_epot(&(enerd->foreign_grpp), enerd->foreign_term);
                enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
            }
        }
        wallcycle_sub_stop(wcycle, ewcsNONBONDED);
        where();
    }

    /* If we are doing GB, calculate bonded forces and apply corrections
     * to the solvation forces */
    /* MRS: Eventually, many need to include free energy contribution here! */
    if (ir->implicit_solvent)
    {
        wallcycle_sub_start(wcycle, ewcsLISTED);
        calc_gb_forces(cr, md, born, top, x, forceForUseWithShiftForces, fr, idef,
                       ir->gb_algorithm, ir->sa_algorithm, nrnb, &pbc, graph, enerd);
        wallcycle_sub_stop(wcycle, ewcsLISTED);
    }

#if GMX_MPI
    if (TAKETIME)
    {
        t1          = MPI_Wtime();
        fr->t_fnbf += t1-t0;
    }
#endif

    if (fepvals->sc_alpha != 0)
    {
        enerd->dvdl_nonlin[efptVDW] += dvdl_nb[efptVDW];
    }
    else
    {
        enerd->dvdl_lin[efptVDW] += dvdl_nb[efptVDW];
    }

    if (fepvals->sc_alpha != 0)

    /* even though coulomb part is linear, we already added it, beacuse we
       need to go through the vdw calculation anyway */
    {
        enerd->dvdl_nonlin[efptCOUL] += dvdl_nb[efptCOUL];
    }
    else
    {
        enerd->dvdl_lin[efptCOUL] += dvdl_nb[efptCOUL];
    }

    if (debug)
    {
        pr_rvecs(debug, 0, "fshift after SR", fr->fshift, SHIFTS);
    }

    /* Shift the coordinates. Must be done before listed forces and PPPM,
     * but is also necessary for SHAKE and update, therefore it can NOT
     * go when no listed forces have to be evaluated.
     *
     * The shifting and PBC code is deliberately not timed, since with
     * the Verlet scheme it only takes non-zero time with triclinic
     * boxes, and even then the time is around a factor of 100 less
     * than the next smallest counter.
     */


    /* Here sometimes we would not need to shift with NBFonly,
     * but we do so anyhow for consistency of the returned coordinates.
     */
    if (graph)
    {
        shift_self(graph, box, x);
        if (TRICLINIC(box))
        {
            inc_nrnb(nrnb, eNR_SHIFTX, 2*graph->nnodes);
        }
        else
        {
            inc_nrnb(nrnb, eNR_SHIFTX, graph->nnodes);
        }
    }
    /* Check whether we need to do listed interactions or correct for exclusions */
    if (fr->bMolPBC &&
        ((flags & GMX_FORCE_LISTED)
         || EEL_RF(fr->ic->eeltype) || EEL_FULL(fr->ic->eeltype) || EVDW_PME(fr->ic->vdwtype)))
    {
        /* TODO There are no electrostatics methods that require this
           transformation, when using the Verlet scheme, so update the
           above conditional. */
        /* Since all atoms are in the rectangular or triclinic unit-cell,
         * only single box vector shifts (2 in x) are required.
         */
        set_pbc_dd(&pbc, fr->ePBC, DOMAINDECOMP(cr) ? cr->dd->nc : nullptr,
                   TRUE, box);
    }

    do_force_listed(wcycle, box, ir->fepvals, cr,
                    idef, (const rvec *) x, hist,
                    forceForUseWithShiftForces, forceWithVirial,
                    fr, &pbc, graph, enerd, nrnb, lambda, md, fcd,
                    DOMAINDECOMP(cr) ? cr->dd->gatindex : nullptr,
                    flags);

    where();

    *cycles_pme = 0;

    /* Do long-range electrostatics and/or LJ-PME, including related short-range
     * corrections.
     */
    if (EEL_FULL(fr->ic->eeltype) || EVDW_PME(fr->ic->vdwtype))
    {
        int  status            = 0;
        real Vlr_q             = 0, Vlr_lj = 0;

        /* We reduce all virial, dV/dlambda and energy contributions, except
         * for the reciprocal energies (Vlr_q, Vlr_lj) into the same struct.
         */
        ewald_corr_thread_t &ewaldOutput = fr->ewc_t[0];
        clearEwaldThreadOutput(&ewaldOutput);

        if (EEL_PME_EWALD(fr->ic->eeltype) || EVDW_PME(fr->ic->vdwtype))
        {
            /* With the Verlet scheme exclusion forces are calculated
             * in the non-bonded kernel.
             */
            /* The TPI molecule does not have exclusions with the rest
             * of the system and no intra-molecular PME grid
             * contributions will be calculated in
             * gmx_pme_calc_energy.
             */
            if ((ir->cutoff_scheme == ecutsGROUP && fr->n_tpi == 0) ||
                ir->ewald_geometry != eewg3D ||
                ir->epsilon_surface != 0)
            {
                int nthreads, t;

                wallcycle_sub_start(wcycle, ewcsEWALD_CORRECTION);

                if (fr->n_tpi > 0)
                {
                    gmx_fatal(FARGS, "TPI with PME currently only works in a 3D geometry with tin-foil boundary conditions");
                }

                nthreads = fr->nthread_ewc;
#pragma omp parallel for num_threads(nthreads) schedule(static)
                for (t = 0; t < nthreads; t++)
                {
                    try
                    {
                        ewald_corr_thread_t &ewc_t = fr->ewc_t[t];
                        if (t > 0)
                        {
                            clearEwaldThreadOutput(&ewc_t);
                        }

                        /* Threading is only supported with the Verlet cut-off
                         * scheme and then only single particle forces (no
                         * exclusion forces) are calculated, so we can store
                         * the forces in the normal, single forceWithVirial->force_ array.
                         */
                        ewald_LRcorrection(md->homenr, cr, nthreads, t, fr, ir,
                                           md->chargeA, md->chargeB,
                                           md->sqrt_c6A, md->sqrt_c6B,
                                           md->sigmaA, md->sigmaB,
                                           md->sigma3A, md->sigma3B,
                                           md->nChargePerturbed || md->nTypePerturbed,
                                           ir->cutoff_scheme != ecutsVERLET,
                                           excl, x, box, mu_tot,
                                           ir->ewald_geometry,
                                           ir->epsilon_surface,
                                           as_rvec_array(forceWithVirial->force_.data()),
                                           ewc_t.vir_q, ewc_t.vir_lj,
                                           &ewc_t.Vcorr_q, &ewc_t.Vcorr_lj,
                                           lambda[efptCOUL], lambda[efptVDW],
                                           &ewc_t.dvdl[efptCOUL], &ewc_t.dvdl[efptVDW]);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                }
                if (nthreads > 1)
                {
                    reduceEwaldThreadOuput(nthreads, fr->ewc_t);
                }
                wallcycle_sub_stop(wcycle, ewcsEWALD_CORRECTION);
            }

            if (EEL_PME_EWALD(fr->ic->eeltype) && fr->n_tpi == 0)
            {
                /* This is not in a subcounter because it takes a
                   negligible and constant-sized amount of time */
                ewaldOutput.Vcorr_q +=
                    ewald_charge_correction(cr, fr, lambda[efptCOUL], box,
                                            &ewaldOutput.dvdl[efptCOUL],
                                            ewaldOutput.vir_q);
            }

            if ((EEL_PME(fr->ic->eeltype) || EVDW_PME(fr->ic->vdwtype)) && (cr->duty & DUTY_PME) && (pme_run_mode(fr->pmedata) == PmeRunMode::CPU))
            {
                /* Do reciprocal PME for Coulomb and/or LJ. */
                assert(fr->n_tpi >= 0);
                if (fr->n_tpi == 0 || (flags & GMX_FORCE_STATECHANGED))
                {
                    pme_flags = GMX_PME_SPREAD | GMX_PME_SOLVE;

                    if (flags & GMX_FORCE_FORCES)
                    {
                        pme_flags |= GMX_PME_CALC_F;
                    }
                    if (flags & GMX_FORCE_VIRIAL)
                    {
                        pme_flags |= GMX_PME_CALC_ENER_VIR;
                    }
                    if (fr->n_tpi > 0)
                    {
                        /* We don't calculate f, but we do want the potential */
                        pme_flags |= GMX_PME_CALC_POT;
                    }

                    /* With domain decomposition we close the CPU side load
                     * balancing region here, because PME does global
                     * communication that acts as a global barrier.
                     */
                    if (DOMAINDECOMP(cr))
                    {
                        ddCloseBalanceRegionCpu(cr->dd);
                    }

                    wallcycle_start(wcycle, ewcPMEMESH);
                    status = gmx_pme_do(fr->pmedata,
                                        0, md->homenr - fr->n_tpi,
                                        x,
                                        as_rvec_array(forceWithVirial->force_.data()),
                                        md->chargeA, md->chargeB,
                                        md->sqrt_c6A, md->sqrt_c6B,
                                        md->sigmaA, md->sigmaB,
                                        box, cr,
                                        DOMAINDECOMP(cr) ? dd_pme_maxshift_x(cr->dd) : 0,
                                        DOMAINDECOMP(cr) ? dd_pme_maxshift_y(cr->dd) : 0,
                                        nrnb, wcycle,
                                        ewaldOutput.vir_q, ewaldOutput.vir_lj,
                                        &Vlr_q, &Vlr_lj,
                                        lambda[efptCOUL], lambda[efptVDW],
                                        &ewaldOutput.dvdl[efptCOUL],
                                        &ewaldOutput.dvdl[efptVDW],
                                        pme_flags);
                    *cycles_pme = wallcycle_stop(wcycle, ewcPMEMESH);
                    if (status != 0)
                    {
                        gmx_fatal(FARGS, "Error %d in reciprocal PME routine", status);
                    }

                    /* We should try to do as little computation after
                     * this as possible, because parallel PME synchronizes
                     * the nodes, so we want all load imbalance of the
                     * rest of the force calculation to be before the PME
                     * call.  DD load balancing is done on the whole time
                     * of the force call (without PME).
                     */
                }
                if (fr->n_tpi > 0)
                {
                    if (EVDW_PME(ir->vdwtype))
                    {

                        gmx_fatal(FARGS, "Test particle insertion not implemented with LJ-PME");
                    }
                    /* Determine the PME grid energy of the test molecule
                     * with the PME grid potential of the other charges.
                     */
                    gmx_pme_calc_energy(fr->pmedata, fr->n_tpi,
                                        x + md->homenr - fr->n_tpi,
                                        md->chargeA + md->homenr - fr->n_tpi,
                                        &Vlr_q);
                }
            }
        }

        if (!EEL_PME(fr->ic->eeltype) && EEL_PME_EWALD(fr->ic->eeltype))
        {
            Vlr_q = do_ewald(ir, x, as_rvec_array(forceWithVirial->force_.data()),
                             md->chargeA, md->chargeB,
                             box, cr, md->homenr,
                             ewaldOutput.vir_q, fr->ic->ewaldcoeff_q,
                             lambda[efptCOUL], &ewaldOutput.dvdl[efptCOUL],
                             fr->ewald_table);
        }

        /* Note that with separate PME nodes we get the real energies later */
        forceWithVirial->addVirialContribution(ewaldOutput.vir_q);
        forceWithVirial->addVirialContribution(ewaldOutput.vir_lj);
        enerd->dvdl_lin[efptCOUL] += ewaldOutput.dvdl[efptCOUL];
        enerd->dvdl_lin[efptVDW]  += ewaldOutput.dvdl[efptVDW];
        enerd->term[F_COUL_RECIP]  = Vlr_q + ewaldOutput.Vcorr_q;
        enerd->term[F_LJ_RECIP]    = Vlr_lj + ewaldOutput.Vcorr_lj;

        if (debug)
        {
            fprintf(debug, "Vlr_q = %g, Vcorr_q = %g, Vlr_corr_q = %g\n",
                    Vlr_q, ewaldOutput.Vcorr_q, enerd->term[F_COUL_RECIP]);
            pr_rvecs(debug, 0, "vir_el_recip after corr", ewaldOutput.vir_q, DIM);
            pr_rvecs(debug, 0, "fshift after LR Corrections", fr->fshift, SHIFTS);
            fprintf(debug, "Vlr_lj: %g, Vcorr_lj = %g, Vlr_corr_lj = %g\n",
                    Vlr_lj, ewaldOutput.Vcorr_lj, enerd->term[F_LJ_RECIP]);
            pr_rvecs(debug, 0, "vir_lj_recip after corr", ewaldOutput.vir_lj, DIM);
        }
    }
    else
    {
        /* Is there a reaction-field exclusion correction needed?
         * With the Verlet scheme, exclusion forces are calculated
         * in the non-bonded kernel.
         */
        if (ir->cutoff_scheme != ecutsVERLET && EEL_RF(fr->ic->eeltype))
        {
            real dvdl_rf_excl      = 0;
            enerd->term[F_RF_EXCL] =
                RF_excl_correction(fr, graph, md, excl, x, forceForUseWithShiftForces,
                                   fr->fshift, &pbc, lambda[efptCOUL], &dvdl_rf_excl);

            enerd->dvdl_lin[efptCOUL] += dvdl_rf_excl;
        }
    }
    where();

    if (debug)
    {
        print_nrnb(debug, nrnb);
    }

#if GMX_MPI
    if (TAKETIME)
    {
        t2 = MPI_Wtime();
        MPI_Barrier(cr->mpi_comm_mygroup);
        t3          = MPI_Wtime();
        fr->t_wait += t3-t2;
        if (fr->timesteps == 11)
        {
            char buf[22];
            fprintf(stderr, "* PP load balancing info: rank %d, step %s, rel wait time=%3.0f%% , load string value: %7.2f\n",
                    cr->nodeid, gmx_step_str(fr->timesteps, buf),
                    100*fr->t_wait/(fr->t_wait+fr->t_fnbf),
                    (fr->t_fnbf+fr->t_wait)/fr->t_fnbf);
        }
        fr->timesteps++;
    }
#endif

    if (debug)
    {
        pr_rvecs(debug, 0, "fshift after bondeds", fr->fshift, SHIFTS);
    }

}

void init_enerdata(int ngener, int n_lambda, gmx_enerdata_t *enerd)
{
    int i, n2;

    for (i = 0; i < F_NRE; i++)
    {
        enerd->term[i]         = 0;
        enerd->foreign_term[i] = 0;
    }


    for (i = 0; i < efptNR; i++)
    {
        enerd->dvdl_lin[i]     = 0;
        enerd->dvdl_nonlin[i]  = 0;
    }

    n2 = ngener*ngener;
    if (debug)
    {
        fprintf(debug, "Creating %d sized group matrix for energies\n", n2);
    }
    enerd->grpp.nener         = n2;
    enerd->foreign_grpp.nener = n2;
    for (i = 0; (i < egNR); i++)
    {
        snew(enerd->grpp.ener[i], n2);
        snew(enerd->foreign_grpp.ener[i], n2);
    }

    if (n_lambda)
    {
        enerd->n_lambda = 1 + n_lambda;
        snew(enerd->enerpart_lambda, enerd->n_lambda);
    }
    else
    {
        enerd->n_lambda = 0;
    }
}

void destroy_enerdata(gmx_enerdata_t *enerd)
{
    int i;

    for (i = 0; (i < egNR); i++)
    {
        sfree(enerd->grpp.ener[i]);
    }

    for (i = 0; (i < egNR); i++)
    {
        sfree(enerd->foreign_grpp.ener[i]);
    }

    if (enerd->n_lambda)
    {
        sfree(enerd->enerpart_lambda);
    }
}

static real sum_v(int n, real v[])
{
    real t;
    int  i;

    t = 0.0;
    for (i = 0; (i < n); i++)
    {
        t = t + v[i];
    }

    return t;
}

void sum_epot(gmx_grppairener_t *grpp, real *epot)
{
    int i;

    /* Accumulate energies */
    epot[F_COUL_SR]  = sum_v(grpp->nener, grpp->ener[egCOULSR]);
    epot[F_LJ]       = sum_v(grpp->nener, grpp->ener[egLJSR]);
    epot[F_LJ14]     = sum_v(grpp->nener, grpp->ener[egLJ14]);
    epot[F_COUL14]   = sum_v(grpp->nener, grpp->ener[egCOUL14]);
    /* We have already added 1-2,1-3, and 1-4 terms to F_GBPOL */
    epot[F_GBPOL]   += sum_v(grpp->nener, grpp->ener[egGB]);

/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
    epot[F_BHAM]     = sum_v(grpp->nener, grpp->ener[egBHAMSR]);

    epot[F_EPOT] = 0;
    for (i = 0; (i < F_EPOT); i++)
    {
        if (i != F_DISRESVIOL && i != F_ORIRESDEV)
        {
            epot[F_EPOT] += epot[i];
        }
    }
}

void sum_dhdl(gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, t_lambda *fepvals)
{
    int    index;
    double dlam;

    enerd->dvdl_lin[efptVDW] += enerd->term[F_DVDL_VDW];  /* include dispersion correction */
    enerd->term[F_DVDL]       = 0.0;
    for (int i = 0; i < efptNR; i++)
    {
        if (fepvals->separate_dvdl[i])
        {
            /* could this be done more readably/compactly? */
            switch (i)
            {
                case (efptMASS):
                    index = F_DKDL;
                    break;
                case (efptCOUL):
                    index = F_DVDL_COUL;
                    break;
                case (efptVDW):
                    index = F_DVDL_VDW;
                    break;
                case (efptBONDED):
                    index = F_DVDL_BONDED;
                    break;
                case (efptRESTRAINT):
                    index = F_DVDL_RESTRAINT;
                    break;
                default:
                    index = F_DVDL;
                    break;
            }
            enerd->term[index] = enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvdl-%s[%2d]: %f: non-linear %f + linear %f\n",
                        efpt_names[i], i, enerd->term[index], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
        else
        {
            enerd->term[F_DVDL] += enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvd-%sl[%2d]: %f: non-linear %f + linear %f\n",
                        efpt_names[0], i, enerd->term[F_DVDL], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
    }

    /* Notes on the foreign lambda free energy difference evaluation:
     * Adding the potential and ekin terms that depend linearly on lambda
     * as delta lam * dvdl to the energy differences is exact.
     * For the constraints this is not exact, but we have no other option
     * without literally changing the lengths and reevaluating the energies at each step.
     * (try to remedy this post 4.6 - MRS)
     */
    if (fepvals->separate_dvdl[efptBONDED])
    {
        enerd->term[F_DVDL_BONDED] += enerd->term[F_DVDL_CONSTR];
    }
    else
    {
        enerd->term[F_DVDL] += enerd->term[F_DVDL_CONSTR];
    }
    enerd->term[F_DVDL_CONSTR] = 0;

    for (int i = 0; i < fepvals->n_lambda; i++)
    {
        /* note we are iterating over fepvals here!
           For the current lam, dlam = 0 automatically,
           so we don't need to add anything to the
           enerd->enerpart_lambda[0] */

        /* we don't need to worry about dvdl_lin contributions to dE at
           current lambda, because the contributions to the current
           lambda are automatically zeroed */

        for (size_t j = 0; j < lambda.size(); j++)
        {
            /* Note that this loop is over all dhdl components, not just the separated ones */
            dlam = (fepvals->all_lambda[j][i] - lambda[j]);
            enerd->enerpart_lambda[i+1] += dlam*enerd->dvdl_lin[j];
            if (debug)
            {
                fprintf(debug, "enerdiff lam %g: (%15s), non-linear %f linear %f*%f\n",
                        fepvals->all_lambda[j][i], efpt_names[j],
                        (enerd->enerpart_lambda[i+1] - enerd->enerpart_lambda[0]),
                        dlam, enerd->dvdl_lin[j]);
            }
        }
    }
}


void reset_foreign_enerdata(gmx_enerdata_t *enerd)
{
    int  i, j;

    /* First reset all foreign energy components.  Foreign energies always called on
       neighbor search steps */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->foreign_grpp.ener[i][j] = 0.0;
        }
    }

    /* potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->foreign_term[i] = 0.0;
    }
}

void reset_enerdata(gmx_enerdata_t *enerd)
{
    int      i, j;

    /* First reset all energy components. */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->grpp.ener[i][j] = 0.0;
        }
    }
    for (i = 0; i < efptNR; i++)
    {
        enerd->dvdl_lin[i]    = 0.0;
        enerd->dvdl_nonlin[i] = 0.0;
    }

    /* Normal potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->term[i] = 0.0;
    }
    enerd->term[F_DVDL]            = 0.0;
    enerd->term[F_DVDL_COUL]       = 0.0;
    enerd->term[F_DVDL_VDW]        = 0.0;
    enerd->term[F_DVDL_BONDED]     = 0.0;
    enerd->term[F_DVDL_RESTRAINT]  = 0.0;
    enerd->term[F_DKDL]            = 0.0;
    if (enerd->n_lambda > 0)
    {
        for (i = 0; i < enerd->n_lambda; i++)
        {
            enerd->enerpart_lambda[i] = 0.0;
        }
    }
    /* reset foreign energy data - separate function since we also call it elsewhere */
    reset_foreign_enerdata(enerd);
}
