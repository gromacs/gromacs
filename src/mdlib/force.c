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

#include <math.h>
#include <string.h>
#include <assert.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "macros.h"
#include "physics.h"
#include "force.h"
#include "nonbonded.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "coulomb.h"
#include "pme.h"
#include "mdrun.h"
#include "domdec.h"
#include "partdec.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "gmx_omp_nthreads.h"

void ns(FILE              *fp,
        t_forcerec        *fr,
        rvec               x[],
        matrix             box,
        gmx_groups_t      *groups,
        t_grpopts         *opts,
        gmx_localtop_t    *top,
        t_mdatoms         *md,
        t_commrec         *cr,
        t_nrnb            *nrnb,
        real              *lambda,
        real              *dvdlambda,
        gmx_grppairener_t *grppener,
        gmx_bool           bFillGrid,
        gmx_bool           bDoLongRangeNS)
{
    char   *ptr;
    int     nsearch;

    GMX_MPE_LOG(ev_ns_start);
    if (!fr->ns.nblist_initialized)
    {
        init_neighbor_list(fp, fr, md->homenr);
    }

    if (fr->bTwinRange)
    {
        fr->nlr = 0;
    }

    nsearch = search_neighbours(fp, fr, x, box, top, groups, cr, nrnb, md,
                                lambda, dvdlambda, grppener,
                                bFillGrid, bDoLongRangeNS, TRUE);
    if (debug)
    {
        fprintf(debug, "nsearch = %d\n", nsearch);
    }

    /* Check whether we have to do dynamic load balancing */
    /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
       count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
       &(top->idef),opts->ngener);
     */
    if (fr->ns.dump_nl > 0)
    {
        dump_nblist(fp, cr, fr, fr->ns.dump_nl);
    }

    GMX_MPE_LOG(ev_ns_finish);
}

static void reduce_thread_forces(int n, rvec *f,
                                 tensor vir,
                                 real *Vcorr,
                                 int efpt_ind, real *dvdl,
                                 int nthreads, f_thread_t *f_t)
{
    int t, i;

    /* This reduction can run over any number of threads */
#pragma omp parallel for num_threads(gmx_omp_nthreads_get(emntBonded)) private(t) schedule(static)
    for (i = 0; i < n; i++)
    {
        for (t = 1; t < nthreads; t++)
        {
            rvec_inc(f[i], f_t[t].f[i]);
        }
    }
    for (t = 1; t < nthreads; t++)
    {
        *Vcorr += f_t[t].Vcorr;
        *dvdl  += f_t[t].dvdl[efpt_ind];
        m_add(vir, f_t[t].vir, vir);
    }
}

void do_force_lowlevel(FILE       *fplog,   gmx_large_int_t step,
                       t_forcerec *fr,      t_inputrec *ir,
                       t_idef     *idef,    t_commrec  *cr,
                       t_nrnb     *nrnb,    gmx_wallcycle_t wcycle,
                       t_mdatoms  *md,
                       t_grpopts  *opts,
                       rvec       x[],      history_t  *hist,
                       rvec       f[],
                       rvec       f_longrange[],
                       gmx_enerdata_t *enerd,
                       t_fcdata   *fcd,
                       gmx_mtop_t     *mtop,
                       gmx_localtop_t *top,
                       gmx_genborn_t *born,
                       t_atomtypes *atype,
                       gmx_bool       bBornRadii,
                       matrix     box,
                       t_lambda   *fepvals,
                       real       *lambda,
                       t_graph    *graph,
                       t_blocka   *excl,
                       rvec       mu_tot[],
                       int        flags,
                       float      *cycles_pme)
{
    int         i, j, status;
    int         donb_flags;
    gmx_bool    bDoEpot, bSepDVDL, bSB;
    int         pme_flags;
    matrix      boxs;
    rvec        box_size;
    real        Vsr, Vlr, Vcorr = 0;
    t_pbc       pbc;
    real        dvdgb;
    char        buf[22];
    double      clam_i, vlam_i;
    real        dvdl_dum[efptNR], dvdl, dvdl_nb[efptNR], lam_i[efptNR];
    real        dvdlsum;

#ifdef GMX_MPI
    double  t0 = 0.0, t1, t2, t3; /* time measurement for coarse load balancing */
#endif

#define PRINT_SEPDVDL(s, v, dvdlambda) if (bSepDVDL) {fprintf(fplog, sepdvdlformat, s, v, dvdlambda); }

    GMX_MPE_LOG(ev_force_start);
    set_pbc(&pbc, fr->ePBC, box);

    /* reset free energy components */
    for (i = 0; i < efptNR; i++)
    {
        dvdl_nb[i]  = 0;
        dvdl_dum[i] = 0;
    }

    /* Reset box */
    for (i = 0; (i < DIM); i++)
    {
        box_size[i] = box[i][i];
    }

    bSepDVDL = (fr->bSepDVDL && do_per_step(step, ir->nstlog));
    debug_gmx();

    /* do QMMM first if requested */
    if (fr->bQMMM)
    {
        enerd->term[F_EQM] = calculate_QMMM(cr, x, f, fr, md);
    }

    if (bSepDVDL)
    {
        fprintf(fplog, "Step %s: non-bonded V and dVdl for node %d:\n",
                gmx_step_str(step, buf), cr->nodeid);
    }

    /* Call the short range functions all in one go. */
    GMX_MPE_LOG(ev_do_fnbf_start);

#ifdef GMX_MPI
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
        dvdl = do_walls(ir, fr, box, md, x, f, lambda[efptVDW],
                        enerd->grpp.ener[egLJSR], nrnb);
        PRINT_SEPDVDL("Walls", 0.0, dvdl);
        enerd->dvdl_lin[efptVDW] += dvdl;
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
            calc_gb_rad(cr, fr, ir, top, atype, x, &(fr->gblist), born, md, nrnb);
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

        if (flags & GMX_FORCE_FORCES)
        {
            donb_flags |= GMX_NONBONDED_DO_FORCE;
        }
        if (flags & GMX_FORCE_ENERGY)
        {
            donb_flags |= GMX_NONBONDED_DO_POTENTIAL;
        }
        if (flags & GMX_FORCE_DO_LR)
        {
            donb_flags |= GMX_NONBONDED_DO_LR;
        }

        wallcycle_sub_start(wcycle, ewcsNONBONDED);
        do_nonbonded(cr, fr, x, f, f_longrange, md, excl,
                     &enerd->grpp, box_size, nrnb,
                     lambda, dvdl_nb, -1, -1, donb_flags);

        /* If we do foreign lambda and we have soft-core interactions
         * we have to recalculate the (non-linear) energies contributions.
         */
        if (fepvals->n_lambda > 0 && (flags & GMX_FORCE_DHDL) && fepvals->sc_alpha != 0)
        {
            for (i = 0; i < enerd->n_lambda; i++)
            {
                for (j = 0; j < efptNR; j++)
                {
                    lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
                }
                reset_foreign_enerdata(enerd);
                do_nonbonded(cr, fr, x, f, f_longrange, md, excl,
                             &(enerd->foreign_grpp), box_size, nrnb,
                             lam_i, dvdl_dum, -1, -1,
                             (donb_flags & ~GMX_NONBONDED_DO_FORCE) | GMX_NONBONDED_DO_FOREIGNLAMBDA);
                sum_epot(&ir->opts, &(enerd->foreign_grpp), enerd->foreign_term);
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
        wallcycle_sub_start(wcycle, ewcsBONDED);
        calc_gb_forces(cr, md, born, top, atype, x, f, fr, idef,
                       ir->gb_algorithm, ir->sa_algorithm, nrnb, bBornRadii, &pbc, graph, enerd);
        wallcycle_sub_stop(wcycle, ewcsBONDED);
    }

#ifdef GMX_MPI
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

    Vsr = 0;
    if (bSepDVDL)
    {
        for (i = 0; i < enerd->grpp.nener; i++)
        {
            Vsr +=
                (fr->bBHAM ?
                 enerd->grpp.ener[egBHAMSR][i] :
                 enerd->grpp.ener[egLJSR][i])
                + enerd->grpp.ener[egCOULSR][i] + enerd->grpp.ener[egGB][i];
        }
        dvdlsum = dvdl_nb[efptVDW] + dvdl_nb[efptCOUL];
        PRINT_SEPDVDL("VdW and Coulomb SR particle-p.", Vsr, dvdlsum);
    }
    debug_gmx();

    GMX_MPE_LOG(ev_do_fnbf_finish);

    if (debug)
    {
        pr_rvecs(debug, 0, "fshift after SR", fr->fshift, SHIFTS);
    }

    /* Shift the coordinates. Must be done before bonded forces and PPPM,
     * but is also necessary for SHAKE and update, therefore it can NOT
     * go when no bonded forces have to be evaluated.
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
    /* Check whether we need to do bondeds or correct for exclusions */
    if (fr->bMolPBC &&
        ((flags & GMX_FORCE_BONDED)
         || EEL_RF(fr->eeltype) || EEL_FULL(fr->eeltype)))
    {
        /* Since all atoms are in the rectangular or triclinic unit-cell,
         * only single box vector shifts (2 in x) are required.
         */
        set_pbc_dd(&pbc, fr->ePBC, cr->dd, TRUE, box);
    }
    debug_gmx();

    if (flags & GMX_FORCE_BONDED)
    {
        GMX_MPE_LOG(ev_calc_bonds_start);

        wallcycle_sub_start(wcycle, ewcsBONDED);
        calc_bonds(fplog, cr->ms,
                   idef, x, hist, f, fr, &pbc, graph, enerd, nrnb, lambda, md, fcd,
                   DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL, atype, born,
                   flags,
                   fr->bSepDVDL && do_per_step(step, ir->nstlog), step);

        /* Check if we have to determine energy differences
         * at foreign lambda's.
         */
        if (fepvals->n_lambda > 0 && (flags & GMX_FORCE_DHDL) &&
            idef->ilsort != ilsortNO_FE)
        {
            if (idef->ilsort != ilsortFE_SORTED)
            {
                gmx_incons("The bonded interactions are not sorted for free energy");
            }
            for (i = 0; i < enerd->n_lambda; i++)
            {
                reset_foreign_enerdata(enerd);
                for (j = 0; j < efptNR; j++)
                {
                    lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
                }
                calc_bonds_lambda(fplog, idef, x, fr, &pbc, graph, &(enerd->foreign_grpp), enerd->foreign_term, nrnb, lam_i, md,
                                  fcd, DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL);
                sum_epot(&ir->opts, &(enerd->foreign_grpp), enerd->foreign_term);
                enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
            }
        }
        debug_gmx();
        GMX_MPE_LOG(ev_calc_bonds_finish);
        wallcycle_sub_stop(wcycle, ewcsBONDED);
    }

    where();

    *cycles_pme = 0;
    if (EEL_FULL(fr->eeltype))
    {
        bSB = (ir->nwall == 2);
        if (bSB)
        {
            copy_mat(box, boxs);
            svmul(ir->wall_ewald_zfac, boxs[ZZ], boxs[ZZ]);
            box_size[ZZ] *= ir->wall_ewald_zfac;
        }

        clear_mat(fr->vir_el_recip);

        if (fr->bEwald)
        {
            Vcorr = 0;
            dvdl  = 0;

            /* With the Verlet scheme exclusion forces are calculated
             * in the non-bonded kernel.
             */
            /* The TPI molecule does not have exclusions with the rest
             * of the system and no intra-molecular PME grid contributions
             * will be calculated in gmx_pme_calc_energy.
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

                nthreads = gmx_omp_nthreads_get(emntBonded);
#pragma omp parallel for num_threads(nthreads) schedule(static)
                for (t = 0; t < nthreads; t++)
                {
                    int     s, e, i;
                    rvec   *fnv;
                    tensor *vir;
                    real   *Vcorrt, *dvdlt;
                    if (t == 0)
                    {
                        fnv    = fr->f_novirsum;
                        vir    = &fr->vir_el_recip;
                        Vcorrt = &Vcorr;
                        dvdlt  = &dvdl;
                    }
                    else
                    {
                        fnv    = fr->f_t[t].f;
                        vir    = &fr->f_t[t].vir;
                        Vcorrt = &fr->f_t[t].Vcorr;
                        dvdlt  = &fr->f_t[t].dvdl[efptCOUL];
                        for (i = 0; i < fr->natoms_force; i++)
                        {
                            clear_rvec(fnv[i]);
                        }
                        clear_mat(*vir);
                    }
                    *dvdlt  = 0;
                    *Vcorrt =
                        ewald_LRcorrection(fplog,
                                           fr->excl_load[t], fr->excl_load[t+1],
                                           cr, t, fr,
                                           md->chargeA,
                                           md->nChargePerturbed ? md->chargeB : NULL,
                                           ir->cutoff_scheme != ecutsVERLET,
                                           excl, x, bSB ? boxs : box, mu_tot,
                                           ir->ewald_geometry,
                                           ir->epsilon_surface,
                                           fnv, *vir,
                                           lambda[efptCOUL], dvdlt);
                }
                if (nthreads > 1)
                {
                    reduce_thread_forces(fr->natoms_force, fr->f_novirsum,
                                         fr->vir_el_recip,
                                         &Vcorr, efptCOUL, &dvdl,
                                         nthreads, fr->f_t);
                }

                wallcycle_sub_stop(wcycle, ewcsEWALD_CORRECTION);
            }

            if (fr->n_tpi == 0)
            {
                Vcorr += ewald_charge_correction(cr, fr, lambda[efptCOUL], box,
                                                 &dvdl, fr->vir_el_recip);
            }

            PRINT_SEPDVDL("Ewald excl./charge/dip. corr.", Vcorr, dvdl);
            enerd->dvdl_lin[efptCOUL] += dvdl;
        }

        status = 0;
        Vlr    = 0;
        dvdl   = 0;
        switch (fr->eeltype)
        {
            case eelPME:
            case eelPMESWITCH:
            case eelPMEUSER:
            case eelPMEUSERSWITCH:
            case eelP3M_AD:
                if (cr->duty & DUTY_PME)
                {
                    assert(fr->n_tpi >= 0);
                    if (fr->n_tpi == 0 || (flags & GMX_FORCE_STATECHANGED))
                    {
                        pme_flags = GMX_PME_SPREAD_Q | GMX_PME_SOLVE;
                        if (flags & GMX_FORCE_FORCES)
                        {
                            pme_flags |= GMX_PME_CALC_F;
                        }
                        if (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY))
                        {
                            pme_flags |= GMX_PME_CALC_ENER_VIR;
                        }
                        if (fr->n_tpi > 0)
                        {
                            /* We don't calculate f, but we do want the potential */
                            pme_flags |= GMX_PME_CALC_POT;
                        }
                        wallcycle_start(wcycle, ewcPMEMESH);
                        status = gmx_pme_do(fr->pmedata,
                                            md->start, md->homenr - fr->n_tpi,
                                            x, fr->f_novirsum,
                                            md->chargeA, md->chargeB,
                                            bSB ? boxs : box, cr,
                                            DOMAINDECOMP(cr) ? dd_pme_maxshift_x(cr->dd) : 0,
                                            DOMAINDECOMP(cr) ? dd_pme_maxshift_y(cr->dd) : 0,
                                            nrnb, wcycle,
                                            fr->vir_el_recip, fr->ewaldcoeff,
                                            &Vlr, lambda[efptCOUL], &dvdl,
                                            pme_flags);
                        *cycles_pme = wallcycle_stop(wcycle, ewcPMEMESH);

                        /* We should try to do as little computation after
                         * this as possible, because parallel PME synchronizes
                         * the nodes, so we want all load imbalance of the rest
                         * of the force calculation to be before the PME call.
                         * DD load balancing is done on the whole time of
                         * the force call (without PME).
                         */
                    }
                    if (fr->n_tpi > 0)
                    {
                        /* Determine the PME grid energy of the test molecule
                         * with the PME grid potential of the other charges.
                         */
                        gmx_pme_calc_energy(fr->pmedata, fr->n_tpi,
                                            x + md->homenr - fr->n_tpi,
                                            md->chargeA + md->homenr - fr->n_tpi,
                                            &Vlr);
                    }
                    PRINT_SEPDVDL("PME mesh", Vlr, dvdl);
                }
                break;
            case eelEWALD:
                Vlr = do_ewald(fplog, FALSE, ir, x, fr->f_novirsum,
                               md->chargeA, md->chargeB,
                               box_size, cr, md->homenr,
                               fr->vir_el_recip, fr->ewaldcoeff,
                               lambda[efptCOUL], &dvdl, fr->ewald_table);
                PRINT_SEPDVDL("Ewald long-range", Vlr, dvdl);
                break;
            default:
                gmx_fatal(FARGS, "No such electrostatics method implemented %s",
                          eel_names[fr->eeltype]);
        }
        if (status != 0)
        {
            gmx_fatal(FARGS, "Error %d in long range electrostatics routine %s",
                      status, EELTYPE(fr->eeltype));
        }
        /* Note that with separate PME nodes we get the real energies later */
        enerd->dvdl_lin[efptCOUL] += dvdl;
        enerd->term[F_COUL_RECIP]  = Vlr + Vcorr;
        if (debug)
        {
            fprintf(debug, "Vlr = %g, Vcorr = %g, Vlr_corr = %g\n",
                    Vlr, Vcorr, enerd->term[F_COUL_RECIP]);
            pr_rvecs(debug, 0, "vir_el_recip after corr", fr->vir_el_recip, DIM);
            pr_rvecs(debug, 0, "fshift after LR Corrections", fr->fshift, SHIFTS);
        }
    }
    else
    {
        if (EEL_RF(fr->eeltype))
        {
            /* With the Verlet scheme exclusion forces are calculated
             * in the non-bonded kernel.
             */
            if (ir->cutoff_scheme != ecutsVERLET && fr->eeltype != eelRF_NEC)
            {
                dvdl                   = 0;
                enerd->term[F_RF_EXCL] =
                    RF_excl_correction(fplog, fr, graph, md, excl, x, f,
                                       fr->fshift, &pbc, lambda[efptCOUL], &dvdl);
            }

            enerd->dvdl_lin[efptCOUL] += dvdl;
            PRINT_SEPDVDL("RF exclusion correction",
                          enerd->term[F_RF_EXCL], dvdl);
        }
    }
    where();
    debug_gmx();

    if (debug)
    {
        print_nrnb(debug, nrnb);
    }
    debug_gmx();

#ifdef GMX_MPI
    if (TAKETIME)
    {
        t2 = MPI_Wtime();
        MPI_Barrier(cr->mpi_comm_mygroup);
        t3          = MPI_Wtime();
        fr->t_wait += t3-t2;
        if (fr->timesteps == 11)
        {
            fprintf(stderr, "* PP load balancing info: node %d, step %s, rel wait time=%3.0f%% , load string value: %7.2f\n",
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

    GMX_MPE_LOG(ev_force_finish);

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

void sum_epot(t_grpopts *opts, gmx_grppairener_t *grpp, real *epot)
{
    int i;

    /* Accumulate energies */
    epot[F_COUL_SR]  = sum_v(grpp->nener, grpp->ener[egCOULSR]);
    epot[F_LJ]       = sum_v(grpp->nener, grpp->ener[egLJSR]);
    epot[F_LJ14]     = sum_v(grpp->nener, grpp->ener[egLJ14]);
    epot[F_COUL14]   = sum_v(grpp->nener, grpp->ener[egCOUL14]);
    epot[F_COUL_LR]  = sum_v(grpp->nener, grpp->ener[egCOULLR]);
    epot[F_LJ_LR]    = sum_v(grpp->nener, grpp->ener[egLJLR]);
    /* We have already added 1-2,1-3, and 1-4 terms to F_GBPOL */
    epot[F_GBPOL]   += sum_v(grpp->nener, grpp->ener[egGB]);

/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
    epot[F_BHAM]     = sum_v(grpp->nener, grpp->ener[egBHAMSR]);
    epot[F_BHAM_LR]  = sum_v(grpp->nener, grpp->ener[egBHAMLR]);

    epot[F_EPOT] = 0;
    for (i = 0; (i < F_EPOT); i++)
    {
        if (i != F_DISRESVIOL && i != F_ORIRESDEV)
        {
            epot[F_EPOT] += epot[i];
        }
    }
}

void sum_dhdl(gmx_enerdata_t *enerd, real *lambda, t_lambda *fepvals)
{
    int    i, j, index;
    double dlam;

    enerd->dvdl_lin[efptVDW] += enerd->term[F_DVDL_VDW];  /* include dispersion correction */
    enerd->term[F_DVDL]       = 0.0;
    for (i = 0; i < efptNR; i++)
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
     * For the non-bonded LR term we assume that the soft-core (if present)
     * no longer affects the energy beyond the short-range cut-off,
     * which is a very good approximation (except for exotic settings).
     * (investigate how to overcome this post 4.6 - MRS)
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

    for (i = 0; i < fepvals->n_lambda; i++)
    {                                         /* note we are iterating over fepvals here!
                                                 For the current lam, dlam = 0 automatically,
                                                 so we don't need to add anything to the
                                                 enerd->enerpart_lambda[0] */

        /* we don't need to worry about dvdl_lin contributions to dE at
           current lambda, because the contributions to the current
           lambda are automatically zeroed */

        for (j = 0; j < efptNR; j++)
        {
            /* Note that this loop is over all dhdl components, not just the separated ones */
            dlam = (fepvals->all_lambda[j][i]-lambda[j]);
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

void reset_enerdata(t_grpopts *opts,
                    t_forcerec *fr, gmx_bool bNS,
                    gmx_enerdata_t *enerd,
                    gmx_bool bMaster)
{
    gmx_bool bKeepLR;
    int      i, j;

    /* First reset all energy components, except for the long range terms
     * on the master at non neighbor search steps, since the long range
     * terms have already been summed at the last neighbor search step.
     */
    bKeepLR = (fr->bTwinRange && !bNS);
    for (i = 0; (i < egNR); i++)
    {
        if (!(bKeepLR && bMaster && (i == egCOULLR || i == egLJLR)))
        {
            for (j = 0; (j < enerd->grpp.nener); j++)
            {
                enerd->grpp.ener[i][j] = 0.0;
            }
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
    /* Initialize the dVdlambda term with the long range contribution */
    /* Initialize the dvdl term with the long range contribution */
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
