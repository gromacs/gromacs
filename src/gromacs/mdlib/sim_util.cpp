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

#include "gromacs/legacyheaders/sim_util.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/imd/imd.h"
#include "gromacs/legacyheaders/calcmu.h"
#include "gromacs/legacyheaders/chargegroup.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/disre.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/genborn.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nonbonded.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/orires.h"
#include "gromacs/legacyheaders/qmmm.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_gpu_ref.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref.h"
#include "gromacs/mdlib/nbnxn_kernels/simd_2xnn/nbnxn_kernel_simd_2xnn.h"
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"

#include "adress.h"
#include "nbnxn_gpu.h"

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

#ifndef GMX_THREAD_MPI
    if (!PAR(cr))
#endif
    {
        fprintf(out, "\r");
    }
    fprintf(out, "step %s", gmx_step_str(step, buf));
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
#ifndef GMX_THREAD_MPI
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

static void sum_forces(int start, int end, rvec f[], rvec flr[])
{
    int i;

    if (gmx_debug_at)
    {
        pr_rvecs(debug, 0, "fsr", f+start, end-start);
        pr_rvecs(debug, 0, "flr", flr+start, end-start);
    }
    for (i = start; (i < end); i++)
    {
        rvec_inc(f[i], flr[i]);
    }
}

/*
 * calc_f_el calculates forces due to an electric field.
 *
 * force is kJ mol^-1 nm^-1 = e * kJ mol^-1 nm^-1 / e
 *
 * Et[] contains the parameters for the time dependent
 * part of the field.
 * Ex[] contains the parameters for
 * the spatial dependent part of the field.
 * The function should return the energy due to the electric field
 * (if any) but for now returns 0.
 *
 * WARNING:
 * There can be problems with the virial.
 * Since the field is not self-consistent this is unavoidable.
 * For neutral molecules the virial is correct within this approximation.
 * For neutral systems with many charged molecules the error is small.
 * But for systems with a net charge or a few charged molecules
 * the error can be significant when the field is high.
 * Solution: implement a self-consistent electric field into PME.
 */
static void calc_f_el(FILE *fp, int  start, int homenr,
                      real charge[], rvec f[],
                      t_cosines Ex[], t_cosines Et[], double t)
{
    rvec Ext;
    real t0;
    int  i, m;

    for (m = 0; (m < DIM); m++)
    {
        if (Et[m].n > 0)
        {
            if (Et[m].n == 3)
            {
                t0     = Et[m].a[1];
                Ext[m] = cos(Et[m].a[0]*(t-t0))*exp(-sqr(t-t0)/(2.0*sqr(Et[m].a[2])));
            }
            else
            {
                Ext[m] = cos(Et[m].a[0]*t);
            }
        }
        else
        {
            Ext[m] = 1.0;
        }
        if (Ex[m].n > 0)
        {
            /* Convert the field strength from V/nm to MD-units */
            Ext[m] *= Ex[m].a[0]*FIELDFAC;
            for (i = start; (i < start+homenr); i++)
            {
                f[i][m] += charge[i]*Ext[m];
            }
        }
        else
        {
            Ext[m] = 0;
        }
    }
    if (fp != NULL)
    {
        fprintf(fp, "%10g  %10g  %10g  %10g #FIELD\n", t,
                Ext[XX]/FIELDFAC, Ext[YY]/FIELDFAC, Ext[ZZ]/FIELDFAC);
    }
}

static void calc_virial(int start, int homenr, rvec x[], rvec f[],
                        tensor vir_part, t_graph *graph, matrix box,
                        t_nrnb *nrnb, const t_forcerec *fr, int ePBC)
{
    int    i;

    /* The short-range virial from surrounding boxes */
    clear_mat(vir_part);
    calc_vir(SHIFTS, fr->shift_vec, fr->fshift, vir_part, ePBC == epbcSCREW, box);
    inc_nrnb(nrnb, eNR_VIRIAL, SHIFTS);

    /* Calculate partial virial, for local atoms only, based on short range.
     * Total virial is computed in global_stat, called from do_md
     */
    f_calc_vir(start, start+homenr, x, f, vir_part, graph, box);
    inc_nrnb(nrnb, eNR_VIRIAL, homenr);

    /* Add position restraint contribution */
    for (i = 0; i < DIM; i++)
    {
        vir_part[i][i] += fr->vir_diag_posres[i];
    }

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
                                   rvec f[],
                                   tensor vir_force,
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
     * The virial contribution is calculated directly,
     * which is why we call pull_potential after calc_virial.
     */
    wallcycle_start(wcycle, ewcPULLPOT);
    set_pbc(&pbc, ir->ePBC, box);
    dvdl                     = 0;
    enerd->term[F_COM_PULL] +=
        pull_potential(ir->pull_work, mdatoms, &pbc,
                       cr, t, lambda[efptRESTRAINT], x, f, vir_force, &dvdl);
    enerd->dvdl_lin[efptRESTRAINT] += dvdl;
    wallcycle_stop(wcycle, ewcPULLPOT);
}

static void pme_receive_force_ener(t_commrec      *cr,
                                   gmx_wallcycle_t wcycle,
                                   gmx_enerdata_t *enerd,
                                   t_forcerec     *fr)
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
    gmx_pme_receive_f(cr, fr->f_novirsum, fr->vir_el_recip, &e_q,
                      fr->vir_lj_recip, &e_lj, &dvdl_q, &dvdl_lj,
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
                               gmx_int64_t step, real pforce, rvec *x, rvec *f)
{
    int  i;
    real pf2, fn2;
    char buf[STEPSTRSIZE];

    pf2 = sqr(pforce);
    for (i = 0; i < md->homenr; i++)
    {
        fn2 = norm2(f[i]);
        /* We also catch NAN, if the compiler does not optimize this away. */
        if (fn2 >= pf2 || fn2 != fn2)
        {
            fprintf(fp, "step %s  atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
                    gmx_step_str(step, buf),
                    ddglatnr(cr->dd, i), x[i][XX], x[i][YY], x[i][ZZ], sqrt(fn2));
        }
    }
}

static void post_process_forces(t_commrec *cr,
                                gmx_int64_t step,
                                t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                                gmx_localtop_t *top,
                                matrix box, rvec x[],
                                rvec f[],
                                tensor vir_force,
                                t_mdatoms *mdatoms,
                                t_graph *graph,
                                t_forcerec *fr, gmx_vsite_t *vsite,
                                int flags)
{
    if (fr->bF_NoVirSum)
    {
        if (vsite)
        {
            /* Spread the mesh force on virtual sites to the other particles...
             * This is parallellized. MPI communication is performed
             * if the constructing atoms aren't local.
             */
            wallcycle_start(wcycle, ewcVSITESPREAD);
            spread_vsite_f(vsite, x, fr->f_novirsum, NULL,
                           (flags & GMX_FORCE_VIRIAL), fr->vir_el_recip,
                           nrnb,
                           &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
            wallcycle_stop(wcycle, ewcVSITESPREAD);
        }
        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Now add the forces, this is local */
            if (fr->bDomDec)
            {
                sum_forces(0, fr->f_novirsum_n, f, fr->f_novirsum);
            }
            else
            {
                sum_forces(0, mdatoms->homenr,
                           f, fr->f_novirsum);
            }
            if (EEL_FULL(fr->eeltype))
            {
                /* Add the mesh contribution to the virial */
                m_add(vir_force, fr->vir_el_recip, vir_force);
            }
            if (EVDW_PME(fr->vdwtype))
            {
                /* Add the mesh contribution to the virial */
                m_add(vir_force, fr->vir_lj_recip, vir_force);
            }
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
                         t_nrnb *nrnb,
                         gmx_wallcycle_t wcycle)
{
    int                        enr_nbnxn_kernel_ljc, enr_nbnxn_kernel_lj;
    nonbonded_verlet_group_t  *nbvg;
    gmx_bool                   bUsingGpuKernels;

    if (!(flags & GMX_FORCE_NONBONDED))
    {
        /* skip non-bonded calculation */
        return;
    }

    nbvg = &fr->nbv->grp[ilocality];

    /* GPU kernel launch overhead is already timed separately */
    if (fr->cutoff_scheme != ecutsVERLET)
    {
        gmx_incons("Invalid cut-off scheme passed!");
    }

    bUsingGpuKernels = (nbvg->kernel_type == nbnxnk8x8x8_GPU);

    if (!bUsingGpuKernels)
    {
        wallcycle_sub_start(wcycle, ewcsNONBONDED);
    }
    switch (nbvg->kernel_type)
    {
        case nbnxnk4x4_PlainC:
            nbnxn_kernel_ref(&nbvg->nbl_lists,
                             nbvg->nbat, ic,
                             fr->shift_vec,
                             flags,
                             clearF,
                             fr->fshift[0],
                             enerd->grpp.ener[egCOULSR],
                             fr->bBHAM ?
                             enerd->grpp.ener[egBHAMSR] :
                             enerd->grpp.ener[egLJSR]);
            break;

        case nbnxnk4xN_SIMD_4xN:
            nbnxn_kernel_simd_4xn(&nbvg->nbl_lists,
                                  nbvg->nbat, ic,
                                  nbvg->ewald_excl,
                                  fr->shift_vec,
                                  flags,
                                  clearF,
                                  fr->fshift[0],
                                  enerd->grpp.ener[egCOULSR],
                                  fr->bBHAM ?
                                  enerd->grpp.ener[egBHAMSR] :
                                  enerd->grpp.ener[egLJSR]);
            break;
        case nbnxnk4xN_SIMD_2xNN:
            nbnxn_kernel_simd_2xnn(&nbvg->nbl_lists,
                                   nbvg->nbat, ic,
                                   nbvg->ewald_excl,
                                   fr->shift_vec,
                                   flags,
                                   clearF,
                                   fr->fshift[0],
                                   enerd->grpp.ener[egCOULSR],
                                   fr->bBHAM ?
                                   enerd->grpp.ener[egBHAMSR] :
                                   enerd->grpp.ener[egLJSR]);
            break;

        case nbnxnk8x8x8_GPU:
            nbnxn_gpu_launch_kernel(fr->nbv->gpu_nbv, nbvg->nbat, flags, ilocality);
            break;

        case nbnxnk8x8x8_PlainC:
            nbnxn_kernel_gpu_ref(nbvg->nbl_lists.nbl[0],
                                 nbvg->nbat, ic,
                                 fr->shift_vec,
                                 flags,
                                 clearF,
                                 nbvg->nbat->out[0].f,
                                 fr->fshift[0],
                                 enerd->grpp.ener[egCOULSR],
                                 fr->bBHAM ?
                                 enerd->grpp.ener[egBHAMSR] :
                                 enerd->grpp.ener[egLJSR]);
            break;

        default:
            gmx_incons("Invalid nonbonded kernel type passed!");

    }
    if (!bUsingGpuKernels)
    {
        wallcycle_sub_stop(wcycle, ewcsNONBONDED);
    }

    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        enr_nbnxn_kernel_ljc = eNR_NBNXN_LJ_RF;
    }
    else if ((!bUsingGpuKernels && nbvg->ewald_excl == ewaldexclAnalytical) ||
             (bUsingGpuKernels && nbnxn_gpu_is_kernel_ewald_analytical(fr->nbv->gpu_nbv)))
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
    if (flags & GMX_FORCE_DO_LR)
    {
        donb_flags |= GMX_NONBONDED_DO_LR;
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

    assert(gmx_omp_nthreads_get(emntNonbonded) == nbl_lists->nnbl);

    wallcycle_sub_start(wcycle, ewcsNONBONDED);
#pragma omp parallel for schedule(static) num_threads(nbl_lists->nnbl)
    for (th = 0; th < nbl_lists->nnbl; th++)
    {
        gmx_nb_free_energy_kernel(nbl_lists->nbl_fep[th],
                                  x, f, fr, mdatoms, &kernel_data, nrnb);
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
                gmx_nb_free_energy_kernel(nbl_lists->nbl_fep[th],
                                          x, f, fr, mdatoms, &kernel_data, nrnb);
            }

            sum_epot(&(enerd->foreign_grpp), enerd->foreign_term);
            enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
        }
    }

    wallcycle_sub_stop(wcycle, ewcsNONBONDED);
}

gmx_bool use_GPU(const nonbonded_verlet_t *nbv)
{
    return nbv != NULL && nbv->bUseGPU;
}

void do_force_cutsVERLET(FILE *fplog, t_commrec *cr,
                         t_inputrec *inputrec,
                         gmx_int64_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                         gmx_localtop_t *top,
                         gmx_groups_t gmx_unused *groups,
                         matrix box, rvec x[], history_t *hist,
                         rvec f[],
                         tensor vir_force,
                         t_mdatoms *mdatoms,
                         gmx_enerdata_t *enerd, t_fcdata *fcd,
                         real *lambda, t_graph *graph,
                         t_forcerec *fr, interaction_const_t *ic,
                         gmx_vsite_t *vsite, rvec mu_tot,
                         double t, FILE *field, gmx_edsam_t ed,
                         gmx_bool bBornRadii,
                         int flags)
{
    int                 cg1, i, j;
    int                 start, homenr;
    double              mu[2*DIM];
    gmx_bool            bStateChanged, bNS, bFillGrid, bCalcCGCM;
    gmx_bool            bDoLongRange, bDoForces, bSepLRF, bUseGPU, bUseOrEmulGPU;
    gmx_bool            bDiffKernels = FALSE;
    rvec                vzero, box_diag;
    float               cycles_pme, cycles_force, cycles_wait_gpu;
    nonbonded_verlet_t *nbv;

    cycles_force    = 0;
    cycles_wait_gpu = 0;
    nbv             = fr->nbv;

    start  = 0;
    homenr = mdatoms->homenr;

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

    bStateChanged = (flags & GMX_FORCE_STATECHANGED);
    bNS           = (flags & GMX_FORCE_NS) && (fr->bAllvsAll == FALSE);
    bFillGrid     = (bNS && bStateChanged);
    bCalcCGCM     = (bFillGrid && !DOMAINDECOMP(cr));
    bDoLongRange  = (fr->bTwinRange && bNS && (flags & GMX_FORCE_DO_LR));
    bDoForces     = (flags & GMX_FORCE_FORCES);
    bSepLRF       = (bDoLongRange && bDoForces && (flags & GMX_FORCE_SEPLRF));
    bUseGPU       = fr->nbv->bUseGPU;
    bUseOrEmulGPU = bUseGPU || (nbv->grp[0].kernel_type == nbnxnk8x8x8_PlainC);

    if (bStateChanged)
    {
        update_forcerec(fr, box);

        if (NEED_MUTOT(*inputrec))
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
            put_atoms_in_box_omp(fr->ePBC, box, homenr, x);
            inc_nrnb(nrnb, eNR_SHIFTX, homenr);
        }
        else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph)
        {
            unshift_self(graph, box, x);
        }
    }

    nbnxn_atomdata_copy_shiftvec(flags & GMX_FORCE_DYNAMICBOX,
                                 fr->shift_vec, nbv->grp[0].nbat);

#ifdef GMX_MPI
    if (!(cr->duty & DUTY_PME))
    {
        gmx_bool bBS;
        matrix   boxs;

        /* Send particle coordinates to the pme nodes.
         * Since this is only implemented for domain decomposition
         * and domain decomposition does not use the graph,
         * we do not need to worry about shifting.
         */

        int pme_flags = 0;

        wallcycle_start(wcycle, ewcPP_PMESENDX);

        bBS = (inputrec->nwall == 2);
        if (bBS)
        {
            copy_mat(box, boxs);
            svmul(inputrec->wall_ewald_zfac, boxs[ZZ], boxs[ZZ]);
        }

        if (EEL_PME(fr->eeltype))
        {
            pme_flags |= GMX_PME_DO_COULOMB;
        }

        if (EVDW_PME(fr->vdwtype))
        {
            pme_flags |= GMX_PME_DO_LJ;
        }

        gmx_pme_send_coordinates(cr, bBS ? boxs : box, x,
                                 mdatoms->nChargePerturbed, mdatoms->nTypePerturbed, lambda[efptCOUL], lambda[efptVDW],
                                 (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY)),
                                 pme_flags, step);

        wallcycle_stop(wcycle, ewcPP_PMESENDX);
    }
#endif /* GMX_MPI */

    /* do gridding for pair search */
    if (bNS)
    {
        if (graph && bStateChanged)
        {
            /* Calculate intramolecular shift vectors to make molecules whole */
            mk_mshift(fplog, graph, fr->ePBC, box, x);
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
                              0, mdatoms->homenr, -1, fr->cginfo, x,
                              0, NULL,
                              nbv->grp[eintLocal].kernel_type,
                              nbv->grp[eintLocal].nbat);
            wallcycle_sub_stop(wcycle, ewcsNBS_GRID_LOCAL);
        }
        else
        {
            wallcycle_sub_start(wcycle, ewcsNBS_GRID_NONLOCAL);
            nbnxn_put_on_grid_nonlocal(nbv->nbs, domdec_zones(cr->dd),
                                       fr->cginfo, x,
                                       nbv->grp[eintNonlocal].kernel_type,
                                       nbv->grp[eintNonlocal].nbat);
            wallcycle_sub_stop(wcycle, ewcsNBS_GRID_NONLOCAL);
        }

        if (nbv->ngrp == 1 ||
            nbv->grp[eintNonlocal].nbat == nbv->grp[eintLocal].nbat)
        {
            nbnxn_atomdata_set(nbv->grp[eintLocal].nbat, eatAll,
                               nbv->nbs, mdatoms, fr->cginfo);
        }
        else
        {
            nbnxn_atomdata_set(nbv->grp[eintLocal].nbat, eatLocal,
                               nbv->nbs, mdatoms, fr->cginfo);
            nbnxn_atomdata_set(nbv->grp[eintNonlocal].nbat, eatAll,
                               nbv->nbs, mdatoms, fr->cginfo);
        }
        wallcycle_stop(wcycle, ewcNS);
    }

    /* initialize the GPU atom data and copy shift vector */
    if (bUseGPU)
    {
        if (bNS)
        {
            wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU_NB);
            nbnxn_gpu_init_atomdata(nbv->gpu_nbv, nbv->grp[eintLocal].nbat);
            wallcycle_stop(wcycle, ewcLAUNCH_GPU_NB);
        }

        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU_NB);
        nbnxn_gpu_upload_shiftvec(nbv->gpu_nbv, nbv->grp[eintLocal].nbat);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU_NB);
    }

    /* do local pair search */
    if (bNS)
    {
        wallcycle_start_nocount(wcycle, ewcNS);
        wallcycle_sub_start(wcycle, ewcsNBS_SEARCH_LOCAL);
        nbnxn_make_pairlist(nbv->nbs, nbv->grp[eintLocal].nbat,
                            &top->excls,
                            ic->rlist,
                            nbv->min_ci_balanced,
                            &nbv->grp[eintLocal].nbl_lists,
                            eintLocal,
                            nbv->grp[eintLocal].kernel_type,
                            nrnb);
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
        nbnxn_atomdata_copy_x_to_nbat_x(nbv->nbs, eatLocal, FALSE, x,
                                        nbv->grp[eintLocal].nbat);
        wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
        wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
    }

    if (bUseGPU)
    {
        wallcycle_start(wcycle, ewcLAUNCH_GPU_NB);
        /* launch local nonbonded F on GPU */
        do_nb_verlet(fr, ic, enerd, flags, eintLocal, enbvClearFNo,
                     nrnb, wcycle);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU_NB);
    }

    /* Communicate coordinates and sum dipole if necessary +
       do non-local pair search */
    if (DOMAINDECOMP(cr))
    {
        bDiffKernels = (nbv->grp[eintNonlocal].kernel_type !=
                        nbv->grp[eintLocal].kernel_type);

        if (bDiffKernels)
        {
            /* With GPU+CPU non-bonded calculations we need to copy
             * the local coordinates to the non-local nbat struct
             * (in CPU format) as the non-local kernel call also
             * calculates the local - non-local interactions.
             */
            wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
            wallcycle_sub_start(wcycle, ewcsNB_X_BUF_OPS);
            nbnxn_atomdata_copy_x_to_nbat_x(nbv->nbs, eatLocal, TRUE, x,
                                            nbv->grp[eintNonlocal].nbat);
            wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
            wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        }

        if (bNS)
        {
            wallcycle_start_nocount(wcycle, ewcNS);
            wallcycle_sub_start(wcycle, ewcsNBS_SEARCH_NONLOCAL);

            if (bDiffKernels)
            {
                nbnxn_grid_add_simple(nbv->nbs, nbv->grp[eintNonlocal].nbat);
            }

            nbnxn_make_pairlist(nbv->nbs, nbv->grp[eintNonlocal].nbat,
                                &top->excls,
                                ic->rlist,
                                nbv->min_ci_balanced,
                                &nbv->grp[eintNonlocal].nbl_lists,
                                eintNonlocal,
                                nbv->grp[eintNonlocal].kernel_type,
                                nrnb);

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
            dd_move_x(cr->dd, box, x);

            /* When we don't need the total dipole we sum it in global_stat */
            if (bStateChanged && NEED_MUTOT(*inputrec))
            {
                gmx_sumd(2*DIM, mu, cr);
            }
            wallcycle_stop(wcycle, ewcMOVEX);

            wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
            wallcycle_sub_start(wcycle, ewcsNB_X_BUF_OPS);
            nbnxn_atomdata_copy_x_to_nbat_x(nbv->nbs, eatNonlocal, FALSE, x,
                                            nbv->grp[eintNonlocal].nbat);
            wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
            cycles_force += wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        }

        if (bUseGPU && !bDiffKernels)
        {
            wallcycle_start(wcycle, ewcLAUNCH_GPU_NB);
            /* launch non-local nonbonded F on GPU */
            do_nb_verlet(fr, ic, enerd, flags, eintNonlocal, enbvClearFNo,
                         nrnb, wcycle);
            cycles_force += wallcycle_stop(wcycle, ewcLAUNCH_GPU_NB);
        }
    }

    if (bUseGPU)
    {
        /* launch D2H copy-back F */
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU_NB);
        if (DOMAINDECOMP(cr) && !bDiffKernels)
        {
            nbnxn_gpu_launch_cpyback(nbv->gpu_nbv, nbv->grp[eintNonlocal].nbat,
                                     flags, eatNonlocal);
        }
        nbnxn_gpu_launch_cpyback(nbv->gpu_nbv, nbv->grp[eintLocal].nbat,
                                 flags, eatLocal);
        cycles_force += wallcycle_stop(wcycle, ewcLAUNCH_GPU_NB);
    }

    if (bStateChanged && NEED_MUTOT(*inputrec))
    {
        if (PAR(cr))
        {
            gmx_sumd(2*DIM, mu, cr);
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
    reset_enerdata(fr, bNS, enerd, MASTER(cr));
    clear_rvecs(SHIFTS, fr->fshift);

    if (DOMAINDECOMP(cr) && !(cr->duty & DUTY_PME))
    {
        wallcycle_start(wcycle, ewcPPDURINGPME);
        dd_force_flop_start(cr->dd, nrnb);
    }

    if (inputrec->bRot)
    {
        /* Enforced rotation has its own cycle counter that starts after the collective
         * coordinates have been communicated. It is added to ddCyclF to allow
         * for proper load-balancing */
        wallcycle_start(wcycle, ewcROT);
        do_rotation(cr, inputrec, box, x, t, step, wcycle, bNS);
        wallcycle_stop(wcycle, ewcROT);
    }

    /* Start the force cycle counter.
     * This counter is stopped after do_force_lowlevel.
     * No parallel communication should occur while this counter is running,
     * since that will interfere with the dynamic load balancing.
     */
    wallcycle_start(wcycle, ewcFORCE);
    if (bDoForces)
    {
        /* Reset forces for which the virial is calculated separately:
         * PME/Ewald forces if necessary */
        if (fr->bF_NoVirSum)
        {
            if (flags & GMX_FORCE_VIRIAL)
            {
                fr->f_novirsum = fr->f_novirsum_alloc;
                if (fr->bDomDec)
                {
                    clear_rvecs(fr->f_novirsum_n, fr->f_novirsum);
                }
                else
                {
                    clear_rvecs(homenr, fr->f_novirsum+start);
                }
            }
            else
            {
                /* We are not calculating the pressure so we do not need
                 * a separate array for forces that do not contribute
                 * to the pressure.
                 */
                fr->f_novirsum = f;
            }
        }

        /* Clear the short- and long-range forces */
        clear_rvecs(fr->natoms_force_constr, f);
        if (bSepLRF && do_per_step(step, inputrec->nstcalclr))
        {
            clear_rvecs(fr->natoms_force_constr, fr->f_twin);
        }

        clear_rvec(fr->vir_diag_posres);
    }

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
        /* Maybe we should move this into do_force_lowlevel */
        do_nb_verlet(fr, ic, enerd, flags, eintLocal, enbvClearFYes,
                     nrnb, wcycle);
    }

    if (fr->efep != efepNO)
    {
        /* Calculate the local and non-local free energy interactions here.
         * Happens here on the CPU both with and without GPU.
         */
        if (fr->nbv->grp[eintLocal].nbl_lists.nbl_fep[0]->nrj > 0)
        {
            do_nb_verlet_fep(&fr->nbv->grp[eintLocal].nbl_lists,
                             fr, x, f, mdatoms,
                             inputrec->fepvals, lambda,
                             enerd, flags, nrnb, wcycle);
        }

        if (DOMAINDECOMP(cr) &&
            fr->nbv->grp[eintNonlocal].nbl_lists.nbl_fep[0]->nrj > 0)
        {
            do_nb_verlet_fep(&fr->nbv->grp[eintNonlocal].nbl_lists,
                             fr, x, f, mdatoms,
                             inputrec->fepvals, lambda,
                             enerd, flags, nrnb, wcycle);
        }
    }

    if (!bUseOrEmulGPU || bDiffKernels)
    {
        int aloc;

        if (DOMAINDECOMP(cr))
        {
            do_nb_verlet(fr, ic, enerd, flags, eintNonlocal,
                         bDiffKernels ? enbvClearFYes : enbvClearFNo,
                         nrnb, wcycle);
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
        cycles_force += wallcycle_stop(wcycle, ewcFORCE);
        wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
        nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatAll, nbv->grp[aloc].nbat, f);
        wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
        cycles_force += wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_start_nocount(wcycle, ewcFORCE);

        /* if there are multiple fshift output buffers reduce them */
        if ((flags & GMX_FORCE_VIRIAL) &&
            nbv->grp[aloc].nbl_lists.nnbl > 1)
        {
            /* This is not in a subcounter because it takes a
               negligible and constant-sized amount of time */
            nbnxn_atomdata_add_nbat_fshift_to_fshift(nbv->grp[aloc].nbat,
                                                     fr->fshift);
        }
    }

    /* update QMMMrec, if necessary */
    if (fr->bQMMM)
    {
        update_QMMMrec(cr, fr, x, mdatoms, box, top);
    }

    /* Compute the bonded and non-bonded energies and optionally forces */
    do_force_lowlevel(fr, inputrec, &(top->idef),
                      cr, nrnb, wcycle, mdatoms,
                      x, hist, f, bSepLRF ? fr->f_twin : f, enerd, fcd, top, fr->born,
                      bBornRadii, box,
                      inputrec->fepvals, lambda, graph, &(top->excls), fr->mu_tot,
                      flags, &cycles_pme);

    if (bSepLRF)
    {
        if (do_per_step(step, inputrec->nstcalclr))
        {
            /* Add the long range forces to the short range forces */
            for (i = 0; i < fr->natoms_force_constr; i++)
            {
                rvec_add(fr->f_twin[i], f[i], f[i]);
            }
        }
    }

    cycles_force += wallcycle_stop(wcycle, ewcFORCE);

    if (ed)
    {
        do_flood(cr, inputrec, x, f, ed, box, step, bNS);
    }

    if (bUseOrEmulGPU && !bDiffKernels)
    {
        /* wait for non-local forces (or calculate in emulation mode) */
        if (DOMAINDECOMP(cr))
        {
            if (bUseGPU)
            {
                float cycles_tmp;

                wallcycle_start(wcycle, ewcWAIT_GPU_NB_NL);
                nbnxn_gpu_wait_for_gpu(nbv->gpu_nbv,
                                       nbv->grp[eintNonlocal].nbat,
                                       flags, eatNonlocal,
                                       enerd->grpp.ener[egLJSR], enerd->grpp.ener[egCOULSR],
                                       fr->fshift);
                cycles_tmp       = wallcycle_stop(wcycle, ewcWAIT_GPU_NB_NL);
                cycles_wait_gpu += cycles_tmp;
                cycles_force    += cycles_tmp;
            }
            else
            {
                wallcycle_start_nocount(wcycle, ewcFORCE);
                do_nb_verlet(fr, ic, enerd, flags, eintNonlocal, enbvClearFYes,
                             nrnb, wcycle);
                cycles_force += wallcycle_stop(wcycle, ewcFORCE);
            }
            wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
            wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
            /* skip the reduction if there was no non-local work to do */
            if (nbv->grp[eintNonlocal].nbl_lists.nbl[0]->nsci > 0)
            {
                nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatNonlocal,
                                               nbv->grp[eintNonlocal].nbat, f);
            }
            wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
            cycles_force += wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
        }
    }

    if (bDoForces && DOMAINDECOMP(cr))
    {
        if (bUseGPU)
        {
            /* We are done with the CPU compute, but the GPU local non-bonded
             * kernel can still be running while we communicate the forces.
             * We start a counter here, so we can, hopefully, time the rest
             * of the GPU kernel execution and data transfer.
             */
            wallcycle_start(wcycle, ewcWAIT_GPU_NB_L_EST);
        }

        /* Communicate the forces */
        wallcycle_start(wcycle, ewcMOVEF);
        dd_move_f(cr->dd, f, fr->fshift);
        if (bSepLRF)
        {
            /* We should not update the shift forces here,
             * since f_twin is already included in f.
             */
            dd_move_f(cr->dd, fr->f_twin, NULL);
        }
        wallcycle_stop(wcycle, ewcMOVEF);
    }

    if (bUseOrEmulGPU)
    {
        /* wait for local forces (or calculate in emulation mode) */
        if (bUseGPU)
        {
#if defined(GMX_GPU) && !defined(GMX_USE_OPENCL)
            float       cycles_tmp, cycles_wait_est;
            const float cuda_api_overhead_margin = 50000.0f; /* cycles */

            wallcycle_start(wcycle, ewcWAIT_GPU_NB_L);
            nbnxn_gpu_wait_for_gpu(nbv->gpu_nbv,
                                   nbv->grp[eintLocal].nbat,
                                   flags, eatLocal,
                                   enerd->grpp.ener[egLJSR], enerd->grpp.ener[egCOULSR],
                                   fr->fshift);
            cycles_tmp      = wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L);

            if (bDoForces && DOMAINDECOMP(cr))
            {
                cycles_wait_est = wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L_EST);

                if (cycles_tmp < cuda_api_overhead_margin)
                {
                    /* We measured few cycles, it could be that the kernel
                     * and transfer finished earlier and there was no actual
                     * wait time, only API call overhead.
                     * Then the actual time could be anywhere between 0 and
                     * cycles_wait_est. As a compromise, we use half the time.
                     */
                    cycles_wait_est *= 0.5f;
                }
            }
            else
            {
                /* No force communication so we actually timed the wait */
                cycles_wait_est = cycles_tmp;
            }
            /* Even though this is after dd_move_f, the actual task we are
             * waiting for runs asynchronously with dd_move_f and we usually
             * have nothing to balance it with, so we can and should add
             * the time to the force time for load balancing.
             */
            cycles_force    += cycles_wait_est;
            cycles_wait_gpu += cycles_wait_est;

#elif defined(GMX_GPU) && defined(GMX_USE_OPENCL)

            wallcycle_start(wcycle, ewcWAIT_GPU_NB_L);
            nbnxn_gpu_wait_for_gpu(nbv->gpu_nbv,
                                   nbv->grp[eintLocal].nbat,
                                   flags, eatLocal,
                                   enerd->grpp.ener[egLJSR], enerd->grpp.ener[egCOULSR],
                                   fr->fshift);
            cycles_wait_gpu += wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L);
#endif

            /* now clear the GPU outputs while we finish the step on the CPU */
            wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU_NB);
            nbnxn_gpu_clear_outputs(nbv->gpu_nbv, flags);
            wallcycle_stop(wcycle, ewcLAUNCH_GPU_NB);
        }
        else
        {
            wallcycle_start_nocount(wcycle, ewcFORCE);
            do_nb_verlet(fr, ic, enerd, flags, eintLocal,
                         DOMAINDECOMP(cr) ? enbvClearFNo : enbvClearFYes,
                         nrnb, wcycle);
            wallcycle_stop(wcycle, ewcFORCE);
        }
        wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
        wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);
        nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs, eatLocal,
                                       nbv->grp[eintLocal].nbat, f);
        wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
        wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
    }

    if (DOMAINDECOMP(cr))
    {
        dd_force_flop_stop(cr->dd, nrnb);
        if (wcycle)
        {
            dd_cycles_add(cr->dd, cycles_force-cycles_pme, ddCyclF);
            if (bUseGPU)
            {
                dd_cycles_add(cr->dd, cycles_wait_gpu, ddCyclWaitGPU);
            }
        }
    }

    if (bDoForces)
    {
        if (IR_ELEC_FIELD(*inputrec))
        {
            /* Compute forces due to electric field */
            calc_f_el(MASTER(cr) ? field : NULL,
                      start, homenr, mdatoms->chargeA, fr->f_novirsum,
                      inputrec->ex, inputrec->et, t);
        }

        /* If we have NoVirSum forces, but we do not calculate the virial,
         * we sum fr->f_novirsum=f later.
         */
        if (vsite && !(fr->bF_NoVirSum && !(flags & GMX_FORCE_VIRIAL)))
        {
            wallcycle_start(wcycle, ewcVSITESPREAD);
            spread_vsite_f(vsite, x, f, fr->fshift, FALSE, NULL, nrnb,
                           &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
            wallcycle_stop(wcycle, ewcVSITESPREAD);

            if (bSepLRF)
            {
                wallcycle_start(wcycle, ewcVSITESPREAD);
                spread_vsite_f(vsite, x, fr->f_twin, NULL, FALSE, NULL,
                               nrnb,
                               &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
                wallcycle_stop(wcycle, ewcVSITESPREAD);
            }
        }

        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Calculation of the virial must be done after vsites! */
            calc_virial(0, mdatoms->homenr, x, f,
                        vir_force, graph, box, nrnb, fr, inputrec->ePBC);
        }
    }

    if (inputrec->bPull && pull_have_potential(inputrec->pull_work))
    {
        /* Since the COM pulling is always done mass-weighted, no forces are
         * applied to vsites and this call can be done after vsite spreading.
         */
        pull_potential_wrapper(cr, inputrec, box, x,
                               f, vir_force, mdatoms, enerd, lambda, t,
                               wcycle);
    }

    /* Add the forces from enforced rotation potentials (if any) */
    if (inputrec->bRot)
    {
        wallcycle_start(wcycle, ewcROTadd);
        enerd->term[F_COM_PULL] += add_rot_forces(inputrec->rot, f, cr, step, t);
        wallcycle_stop(wcycle, ewcROTadd);
    }

    /* Add forces from interactive molecular dynamics (IMD), if bIMD == TRUE. */
    IMD_apply_forces(inputrec->bIMD, inputrec->imd, cr, f, wcycle);

    if (PAR(cr) && !(cr->duty & DUTY_PME))
    {
        /* In case of node-splitting, the PP nodes receive the long-range
         * forces, virial and energy from the PME nodes here.
         */
        pme_receive_force_ener(cr, wcycle, enerd, fr);
    }

    if (bDoForces)
    {
        post_process_forces(cr, step, nrnb, wcycle,
                            top, box, x, f, vir_force, mdatoms, graph, fr, vsite,
                            flags);
    }

    /* Sum the potential energy terms from group contributions */
    sum_epot(&(enerd->grpp), enerd->term);
}

void do_force_cutsGROUP(FILE *fplog, t_commrec *cr,
                        t_inputrec *inputrec,
                        gmx_int64_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                        gmx_localtop_t *top,
                        gmx_groups_t *groups,
                        matrix box, rvec x[], history_t *hist,
                        rvec f[],
                        tensor vir_force,
                        t_mdatoms *mdatoms,
                        gmx_enerdata_t *enerd, t_fcdata *fcd,
                        real *lambda, t_graph *graph,
                        t_forcerec *fr, gmx_vsite_t *vsite, rvec mu_tot,
                        double t, FILE *field, gmx_edsam_t ed,
                        gmx_bool bBornRadii,
                        int flags)
{
    int        cg0, cg1, i, j;
    int        start, homenr;
    double     mu[2*DIM];
    gmx_bool   bStateChanged, bNS, bFillGrid, bCalcCGCM;
    gmx_bool   bDoLongRangeNS, bDoForces, bSepLRF;
    gmx_bool   bDoAdressWF;
    t_pbc      pbc;
    float      cycles_pme, cycles_force;

    start  = 0;
    homenr = mdatoms->homenr;

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
    /* Should we update the long-range neighborlists at this step? */
    bDoLongRangeNS = fr->bTwinRange && bNS;
    /* Should we perform the long-range nonbonded evaluation inside the neighborsearching? */
    bFillGrid      = (bNS && bStateChanged);
    bCalcCGCM      = (bFillGrid && !DOMAINDECOMP(cr));
    bDoForces      = (flags & GMX_FORCE_FORCES);
    bSepLRF        = ((inputrec->nstcalclr > 1) && bDoForces &&
                      (flags & GMX_FORCE_SEPLRF) && (flags & GMX_FORCE_DO_LR));

    /* should probably move this to the forcerec since it doesn't change */
    bDoAdressWF   = ((fr->adress_type != eAdressOff));

    if (bStateChanged)
    {
        update_forcerec(fr, box);

        if (NEED_MUTOT(*inputrec))
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
                                     &(top->cgs), x, fr->cg_cm);
            inc_nrnb(nrnb, eNR_CGCM, homenr);
            inc_nrnb(nrnb, eNR_RESETX, cg1-cg0);
        }
        else if (EI_ENERGY_MINIMIZATION(inputrec->eI) && graph)
        {
            unshift_self(graph, box, x);
        }
    }
    else if (bCalcCGCM)
    {
        calc_cgcm(fplog, cg0, cg1, &(top->cgs), x, fr->cg_cm);
        inc_nrnb(nrnb, eNR_CGCM, homenr);
    }

    if (bCalcCGCM && gmx_debug_at)
    {
        pr_rvecs(debug, 0, "cgcm", fr->cg_cm, top->cgs.nr);
    }

#ifdef GMX_MPI
    if (!(cr->duty & DUTY_PME))
    {
        gmx_bool bBS;
        matrix   boxs;

        /* Send particle coordinates to the pme nodes.
         * Since this is only implemented for domain decomposition
         * and domain decomposition does not use the graph,
         * we do not need to worry about shifting.
         */

        int pme_flags = 0;

        wallcycle_start(wcycle, ewcPP_PMESENDX);

        bBS = (inputrec->nwall == 2);
        if (bBS)
        {
            copy_mat(box, boxs);
            svmul(inputrec->wall_ewald_zfac, boxs[ZZ], boxs[ZZ]);
        }

        if (EEL_PME(fr->eeltype))
        {
            pme_flags |= GMX_PME_DO_COULOMB;
        }

        if (EVDW_PME(fr->vdwtype))
        {
            pme_flags |= GMX_PME_DO_LJ;
        }

        gmx_pme_send_coordinates(cr, bBS ? boxs : box, x,
                                 mdatoms->nChargePerturbed, mdatoms->nTypePerturbed, lambda[efptCOUL], lambda[efptVDW],
                                 (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY)),
                                 pme_flags, step);

        wallcycle_stop(wcycle, ewcPP_PMESENDX);
    }
#endif /* GMX_MPI */

    /* Communicate coordinates and sum dipole if necessary */
    if (DOMAINDECOMP(cr))
    {
        wallcycle_start(wcycle, ewcMOVEX);
        dd_move_x(cr->dd, box, x);
        wallcycle_stop(wcycle, ewcMOVEX);
    }

    /* update adress weight beforehand */
    if (bStateChanged && bDoAdressWF)
    {
        /* need pbc for adress weight calculation with pbc_dx */
        set_pbc(&pbc, inputrec->ePBC, box);
        if (fr->adress_site == eAdressSITEcog)
        {
            update_adress_weights_cog(top->idef.iparams, top->idef.il, x, fr, mdatoms,
                                      inputrec->ePBC == epbcNONE ? NULL : &pbc);
        }
        else if (fr->adress_site == eAdressSITEcom)
        {
            update_adress_weights_com(fplog, cg0, cg1, &(top->cgs), x, fr, mdatoms,
                                      inputrec->ePBC == epbcNONE ? NULL : &pbc);
        }
        else if (fr->adress_site == eAdressSITEatomatom)
        {
            update_adress_weights_atom_per_atom(cg0, cg1, &(top->cgs), x, fr, mdatoms,
                                                inputrec->ePBC == epbcNONE ? NULL : &pbc);
        }
        else
        {
            update_adress_weights_atom(cg0, cg1, &(top->cgs), x, fr, mdatoms,
                                       inputrec->ePBC == epbcNONE ? NULL : &pbc);
        }
    }

    if (NEED_MUTOT(*inputrec))
    {

        if (bStateChanged)
        {
            if (PAR(cr))
            {
                gmx_sumd(2*DIM, mu, cr);
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
    reset_enerdata(fr, bNS, enerd, MASTER(cr));
    clear_rvecs(SHIFTS, fr->fshift);

    if (bNS)
    {
        wallcycle_start(wcycle, ewcNS);

        if (graph && bStateChanged)
        {
            /* Calculate intramolecular shift vectors to make molecules whole */
            mk_mshift(fplog, graph, fr->ePBC, box, x);
        }

        /* Do the actual neighbour searching */
        ns(fplog, fr, box,
           groups, top, mdatoms,
           cr, nrnb, bFillGrid,
           bDoLongRangeNS);

        wallcycle_stop(wcycle, ewcNS);
    }

    if (inputrec->implicit_solvent && bNS)
    {
        make_gb_nblist(cr, inputrec->gb_algorithm,
                       x, box, fr, &top->idef, graph, fr->born);
    }

    if (DOMAINDECOMP(cr) && !(cr->duty & DUTY_PME))
    {
        wallcycle_start(wcycle, ewcPPDURINGPME);
        dd_force_flop_start(cr->dd, nrnb);
    }

    if (inputrec->bRot)
    {
        /* Enforced rotation has its own cycle counter that starts after the collective
         * coordinates have been communicated. It is added to ddCyclF to allow
         * for proper load-balancing */
        wallcycle_start(wcycle, ewcROT);
        do_rotation(cr, inputrec, box, x, t, step, wcycle, bNS);
        wallcycle_stop(wcycle, ewcROT);
    }

    /* Start the force cycle counter.
     * This counter is stopped after do_force_lowlevel.
     * No parallel communication should occur while this counter is running,
     * since that will interfere with the dynamic load balancing.
     */
    wallcycle_start(wcycle, ewcFORCE);

    if (bDoForces)
    {
        /* Reset forces for which the virial is calculated separately:
         * PME/Ewald forces if necessary */
        if (fr->bF_NoVirSum)
        {
            if (flags & GMX_FORCE_VIRIAL)
            {
                fr->f_novirsum = fr->f_novirsum_alloc;
                if (fr->bDomDec)
                {
                    clear_rvecs(fr->f_novirsum_n, fr->f_novirsum);
                }
                else
                {
                    clear_rvecs(homenr, fr->f_novirsum+start);
                }
            }
            else
            {
                /* We are not calculating the pressure so we do not need
                 * a separate array for forces that do not contribute
                 * to the pressure.
                 */
                fr->f_novirsum = f;
            }
        }

        /* Clear the short- and long-range forces */
        clear_rvecs(fr->natoms_force_constr, f);
        if (bSepLRF && do_per_step(step, inputrec->nstcalclr))
        {
            clear_rvecs(fr->natoms_force_constr, fr->f_twin);
        }

        clear_rvec(fr->vir_diag_posres);
    }
    if (inputrec->bPull && pull_have_constraint(inputrec->pull_work))
    {
        clear_pull_forces(inputrec->pull_work);
    }

    /* update QMMMrec, if necessary */
    if (fr->bQMMM)
    {
        update_QMMMrec(cr, fr, x, mdatoms, box, top);
    }

    /* Compute the bonded and non-bonded energies and optionally forces */
    do_force_lowlevel(fr, inputrec, &(top->idef),
                      cr, nrnb, wcycle, mdatoms,
                      x, hist, f, bSepLRF ? fr->f_twin : f, enerd, fcd, top, fr->born,
                      bBornRadii, box,
                      inputrec->fepvals, lambda,
                      graph, &(top->excls), fr->mu_tot,
                      flags,
                      &cycles_pme);

    if (bSepLRF)
    {
        if (do_per_step(step, inputrec->nstcalclr))
        {
            /* Add the long range forces to the short range forces */
            for (i = 0; i < fr->natoms_force_constr; i++)
            {
                rvec_add(fr->f_twin[i], f[i], f[i]);
            }
        }
    }

    cycles_force = wallcycle_stop(wcycle, ewcFORCE);

    if (ed)
    {
        do_flood(cr, inputrec, x, f, ed, box, step, bNS);
    }

    if (DOMAINDECOMP(cr))
    {
        dd_force_flop_stop(cr->dd, nrnb);
        if (wcycle)
        {
            dd_cycles_add(cr->dd, cycles_force-cycles_pme, ddCyclF);
        }
    }

    if (bDoForces)
    {
        if (IR_ELEC_FIELD(*inputrec))
        {
            /* Compute forces due to electric field */
            calc_f_el(MASTER(cr) ? field : NULL,
                      start, homenr, mdatoms->chargeA, fr->f_novirsum,
                      inputrec->ex, inputrec->et, t);
        }

        if (bDoAdressWF && fr->adress_icor == eAdressICThermoForce)
        {
            /* Compute thermodynamic force in hybrid AdResS region */
            adress_thermo_force(start, homenr, &(top->cgs), x, fr->f_novirsum, fr, mdatoms,
                                inputrec->ePBC == epbcNONE ? NULL : &pbc);
        }

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
            if (EEL_FULL(fr->eeltype) && cr->dd->n_intercg_excl &&
                (flags & GMX_FORCE_VIRIAL))
            {
                dd_move_f(cr->dd, fr->f_novirsum, NULL);
            }
            if (bSepLRF)
            {
                /* We should not update the shift forces here,
                 * since f_twin is already included in f.
                 */
                dd_move_f(cr->dd, fr->f_twin, NULL);
            }
            wallcycle_stop(wcycle, ewcMOVEF);
        }

        /* If we have NoVirSum forces, but we do not calculate the virial,
         * we sum fr->f_novirsum=f later.
         */
        if (vsite && !(fr->bF_NoVirSum && !(flags & GMX_FORCE_VIRIAL)))
        {
            wallcycle_start(wcycle, ewcVSITESPREAD);
            spread_vsite_f(vsite, x, f, fr->fshift, FALSE, NULL, nrnb,
                           &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
            wallcycle_stop(wcycle, ewcVSITESPREAD);

            if (bSepLRF)
            {
                wallcycle_start(wcycle, ewcVSITESPREAD);
                spread_vsite_f(vsite, x, fr->f_twin, NULL, FALSE, NULL,
                               nrnb,
                               &top->idef, fr->ePBC, fr->bMolPBC, graph, box, cr);
                wallcycle_stop(wcycle, ewcVSITESPREAD);
            }
        }

        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Calculation of the virial must be done after vsites! */
            calc_virial(0, mdatoms->homenr, x, f,
                        vir_force, graph, box, nrnb, fr, inputrec->ePBC);
        }
    }

    if (inputrec->bPull && pull_have_potential(inputrec->pull_work))
    {
        pull_potential_wrapper(cr, inputrec, box, x,
                               f, vir_force, mdatoms, enerd, lambda, t,
                               wcycle);
    }

    /* Add the forces from enforced rotation potentials (if any) */
    if (inputrec->bRot)
    {
        wallcycle_start(wcycle, ewcROTadd);
        enerd->term[F_COM_PULL] += add_rot_forces(inputrec->rot, f, cr, step, t);
        wallcycle_stop(wcycle, ewcROTadd);
    }

    /* Add forces from interactive molecular dynamics (IMD), if bIMD == TRUE. */
    IMD_apply_forces(inputrec->bIMD, inputrec->imd, cr, f, wcycle);

    if (PAR(cr) && !(cr->duty & DUTY_PME))
    {
        /* In case of node-splitting, the PP nodes receive the long-range
         * forces, virial and energy from the PME nodes here.
         */
        pme_receive_force_ener(cr, wcycle, enerd, fr);
    }

    if (bDoForces)
    {
        post_process_forces(cr, step, nrnb, wcycle,
                            top, box, x, f, vir_force, mdatoms, graph, fr, vsite,
                            flags);
    }

    /* Sum the potential energy terms from group contributions */
    sum_epot(&(enerd->grpp), enerd->term);
}

void do_force(FILE *fplog, t_commrec *cr,
              t_inputrec *inputrec,
              gmx_int64_t step, t_nrnb *nrnb, gmx_wallcycle_t wcycle,
              gmx_localtop_t *top,
              gmx_groups_t *groups,
              matrix box, rvec x[], history_t *hist,
              rvec f[],
              tensor vir_force,
              t_mdatoms *mdatoms,
              gmx_enerdata_t *enerd, t_fcdata *fcd,
              real *lambda, t_graph *graph,
              t_forcerec *fr,
              gmx_vsite_t *vsite, rvec mu_tot,
              double t, FILE *field, gmx_edsam_t ed,
              gmx_bool bBornRadii,
              int flags)
{
    /* modify force flag if not doing nonbonded */
    if (!fr->bNonbonded)
    {
        flags &= ~GMX_FORCE_NONBONDED;
    }

    switch (inputrec->cutoff_scheme)
    {
        case ecutsVERLET:
            do_force_cutsVERLET(fplog, cr, inputrec,
                                step, nrnb, wcycle,
                                top,
                                groups,
                                box, x, hist,
                                f, vir_force,
                                mdatoms,
                                enerd, fcd,
                                lambda, graph,
                                fr, fr->ic,
                                vsite, mu_tot,
                                t, field, ed,
                                bBornRadii,
                                flags);
            break;
        case ecutsGROUP:
            do_force_cutsGROUP(fplog, cr, inputrec,
                               step, nrnb, wcycle,
                               top,
                               groups,
                               box, x, hist,
                               f, vir_force,
                               mdatoms,
                               enerd, fcd,
                               lambda, graph,
                               fr, vsite, mu_tot,
                               t, field, ed,
                               bBornRadii,
                               flags);
            break;
        default:
            gmx_incons("Invalid cut-off scheme passed!");
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

    snew(savex, state->natoms);

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
    constrain(NULL, TRUE, FALSE, constr, &(top->idef),
              ir, cr, step, 0, 1.0, md,
              state->x, state->x, NULL,
              fr->bMolPBC, state->box,
              state->lambda[efptBONDED], &dvdl_dum,
              NULL, NULL, nrnb, econqCoord);
    if (EI_VV(ir->eI))
    {
        /* constrain the inital velocity, and save it */
        /* also may be useful if we need the ekin from the halfstep for velocity verlet */
        constrain(NULL, TRUE, FALSE, constr, &(top->idef),
                  ir, cr, step, 0, 1.0, md,
                  state->x, state->v, state->v,
                  fr->bMolPBC, state->box,
                  state->lambda[efptBONDED], &dvdl_dum,
                  NULL, NULL, nrnb, econqVeloc);
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
        constrain(NULL, TRUE, FALSE, constr, &(top->idef),
                  ir, cr, step, -1, 1.0, md,
                  state->x, savex, NULL,
                  fr->bMolPBC, state->box,
                  state->lambda[efptBONDED], &dvdl_dum,
                  state->v, NULL, nrnb, econqCoord);

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

    if (eDispCorr != edispcNO)
    {
        for (i = 0; i < 2; i++)
        {
            eners[i] = 0;
            virs[i]  = 0;
        }
        if ((fr->vdw_modifier == eintmodPOTSHIFT) ||
            (fr->vdw_modifier == eintmodPOTSWITCH) ||
            (fr->vdw_modifier == eintmodFORCESWITCH) ||
            (fr->vdwtype == evdwSHIFT) ||
            (fr->vdwtype == evdwSWITCH))
        {
            if (((fr->vdw_modifier == eintmodPOTSWITCH) ||
                 (fr->vdw_modifier == eintmodFORCESWITCH) ||
                 (fr->vdwtype == evdwSWITCH)) && fr->rvdw_switch == 0)
            {
                gmx_fatal(FARGS,
                          "With dispersion correction rvdw-switch can not be zero "
                          "for vdw-type = %s", evdw_names[fr->vdwtype]);
            }

            scale  = fr->nblists[0].table_vdw.scale;
            vdwtab = fr->nblists[0].table_vdw.data;

            /* Round the cut-offs to exact table values for precision */
            ri0  = static_cast<int>(floor(fr->rvdw_switch*scale));
            ri1  = static_cast<int>(ceil(fr->rvdw*scale));

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
            ri0  = (fr->vdw_modifier == eintmodPOTSHIFT) ? ri1 : ri0;

            r0   = ri0/scale;
            rc3  = r0*r0*r0;
            rc9  = rc3*rc3*rc3;

            if ((fr->vdw_modifier == eintmodFORCESWITCH) ||
                (fr->vdwtype == evdwSHIFT))
            {
                /* Determine the constant energy shift below rvdw_switch.
                 * Table has a scale factor since we have scaled it down to compensate
                 * for scaling-up c6/c12 with the derivative factors to save flops in analytical kernels.
                 */
                fr->enershiftsix    = (real)(-1.0/(rc3*rc3)) - 6.0*vdwtab[8*ri0];
                fr->enershifttwelve = (real)( 1.0/(rc9*rc3)) - 12.0*vdwtab[8*ri0 + 4];
            }
            else if (fr->vdw_modifier == eintmodPOTSHIFT)
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
        else if (fr->vdwtype == evdwCUT ||
                 EVDW_PME(fr->vdwtype) ||
                 fr->vdwtype == evdwUSER)
        {
            if (fr->vdwtype == evdwUSER && fplog)
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

            rc3  = fr->rvdw*fr->rvdw*fr->rvdw;
            rc9  = rc3*rc3*rc3;
            /* Contribution beyond the cut-off */
            eners[0] += -4.0*M_PI/(3.0*rc3);
            eners[1] +=  4.0*M_PI/(9.0*rc9);
            if (fr->vdw_modifier == eintmodPOTSHIFT)
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
                      evdw_names[fr->vdwtype]);
        }

        /* When we deprecate the group kernels the code below can go too */
        if (fr->vdwtype == evdwPME && fr->cutoff_scheme == ecutsGROUP)
        {
            /* Calculate self-interaction coefficient (assuming that
             * the reciprocal-space contribution is constant in the
             * region that contributes to the self-interaction).
             */
            fr->enershiftsix = pow(fr->ewaldcoeff_lj, 6) / 6.0;

            eners[0] += -pow(sqrt(M_PI)*fr->ewaldcoeff_lj, 3)/3.0;
            virs[0]  +=  pow(sqrt(M_PI)*fr->ewaldcoeff_lj, 3);
        }

        fr->enerdiffsix    = eners[0];
        fr->enerdifftwelve = eners[1];
        /* The 0.5 is due to the Gromacs definition of the virial */
        fr->virdiffsix     = 0.5*virs[0];
        fr->virdifftwelve  = 0.5*virs[1];
    }
}

void calc_dispcorr(t_inputrec *ir, t_forcerec *fr,
                   int natoms,
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
            dens   = (natoms - fr->n_tpi)*invvol;
            ninter = fr->n_tpi;
        }
        else
        {
            dens   = natoms*invvol;
            ninter = 0.5*natoms;
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

void do_pbc_first(FILE *fplog, matrix box, t_forcerec *fr,
                  t_graph *graph, rvec x[])
{
    if (fplog)
    {
        fprintf(fplog, "Removing pbc first time\n");
    }
    calc_shifts(box, fr->shift_vec);
    if (graph)
    {
        mk_mshift(fplog, graph, fr->ePBC, box, x);
        if (gmx_debug_at)
        {
            p_graph(debug, "do_pbc_first 1", graph);
        }
        shift_self(graph, box, x);
        /* By doing an extra mk_mshift the molecules that are broken
         * because they were e.g. imported from another software
         * will be made whole again. Such are the healing powers
         * of GROMACS.
         */
        mk_mshift(fplog, graph, fr->ePBC, box, x);
        if (gmx_debug_at)
        {
            p_graph(debug, "do_pbc_first 2", graph);
        }
    }
    if (fplog)
    {
        fprintf(fplog, "Done rmpbc\n");
    }
}

static void low_do_pbc_mtop(FILE *fplog, int ePBC, matrix box,
                            gmx_mtop_t *mtop, rvec x[],
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
            mk_graph_ilist(NULL, mtop->moltype[molb->type].ilist,
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

void do_pbc_first_mtop(FILE *fplog, int ePBC, matrix box,
                       gmx_mtop_t *mtop, rvec x[])
{
    low_do_pbc_mtop(fplog, ePBC, box, mtop, x, TRUE);
}

void do_pbc_mtop(FILE *fplog, int ePBC, matrix box,
                 gmx_mtop_t *mtop, rvec x[])
{
    low_do_pbc_mtop(fplog, ePBC, box, mtop, x, FALSE);
}

void finish_run(FILE *fplog, t_commrec *cr,
                t_inputrec *inputrec,
                t_nrnb nrnb[], gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                nonbonded_verlet_t *nbv,
                gmx_bool bWriteStat)
{
    t_nrnb *nrnb_tot = NULL;
    double  delta_t  = 0;
    double  nbfs     = 0, mflop = 0;
    double  elapsed_time,
            elapsed_time_over_all_ranks,
            elapsed_time_over_all_threads,
            elapsed_time_over_all_threads_over_all_ranks;
    wallcycle_sum(cr, wcycle);

    if (cr->nnodes > 1)
    {
        snew(nrnb_tot, 1);
#ifdef GMX_MPI
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
#ifdef GMX_MPI
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

    if (SIMMASTER(cr))
    {
        print_flop(fplog, nrnb_tot, &nbfs, &mflop);
    }
    if (cr->nnodes > 1)
    {
        sfree(nrnb_tot);
    }

    if ((cr->duty & DUTY_PP) && DOMAINDECOMP(cr))
    {
        print_dd_statistics(cr, inputrec, fplog);
    }

    if (SIMMASTER(cr))
    {
        struct gmx_wallclock_gpu_t* gputimes = use_GPU(nbv) ? nbnxn_gpu_get_timings(nbv->gpu_nbv) : NULL;

        wallcycle_print(fplog, cr->nnodes, cr->npmenodes,
                        elapsed_time_over_all_ranks,
                        wcycle, gputimes);

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

extern void initialize_lambdas(FILE *fplog, t_inputrec *ir, int *fep_state, real *lambda, double *lam0)
{
    /* this function works, but could probably use a logic rewrite to keep all the different
       types of efep straight. */

    int       i;
    t_lambda *fep = ir->fepvals;

    if ((ir->efep == efepNO) && (ir->bSimTemp == FALSE))
    {
        for (i = 0; i < efptNR; i++)
        {
            lambda[i] = 0.0;
            if (lam0)
            {
                lam0[i] = 0.0;
            }
        }
        return;
    }
    else
    {
        *fep_state = fep->init_fep_state; /* this might overwrite the checkpoint
                                             if checkpoint is set -- a kludge is in for now
                                             to prevent this.*/
        for (i = 0; i < efptNR; i++)
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
            for (i = 0; i < ir->opts.ngtc; i++)
            {
                if (ir->opts.ref_t[i] > 0)
                {
                    ir->opts.ref_t[i] = ir->simtempvals->temperatures[*fep_state];
                }
            }
        }
    }

    /* Send to the log the information on the current lambdas */
    if (fplog != NULL)
    {
        fprintf(fplog, "Initial vector of lambda components:[ ");
        for (i = 0; i < efptNR; i++)
        {
            fprintf(fplog, "%10.4f ", lambda[i]);
        }
        fprintf(fplog, "]\n");
    }
    return;
}


void init_md(FILE *fplog,
             t_commrec *cr, t_inputrec *ir, const output_env_t oenv,
             double *t, double *t0,
             real *lambda, int *fep_state, double *lam0,
             t_nrnb *nrnb, gmx_mtop_t *mtop,
             gmx_update_t *upd,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             tensor force_vir, tensor shake_vir, rvec mu_tot,
             gmx_bool *bSimAnn, t_vcm **vcm, unsigned long Flags,
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
    if (*bSimAnn)
    {
        update_annealing_target_temp(&(ir->opts), ir->init_t);
    }

    /* Initialize lambda variables */
    initialize_lambdas(fplog, ir, fep_state, lambda, lam0);

    if (upd)
    {
        *upd = init_update(ir);
    }


    if (vcm != NULL)
    {
        *vcm = init_vcm(fplog, &mtop->groups, ir);
    }

    if (EI_DYNAMICS(ir->eI) && !(Flags & MD_APPENDFILES))
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
    if ((ir->et[XX].n > 0) || (ir->et[YY].n > 0) || (ir->et[ZZ].n > 0))
    {
        please_cite(fplog, "Caleman2008a");
    }
    init_nrnb(nrnb);

    if (nfile != -1)
    {
        *outf = init_mdoutf(fplog, nfile, fnm, Flags, cr, ir, mtop, oenv, wcycle);

        *mdebin = init_mdebin((Flags & MD_APPENDFILES) ? NULL : mdoutf_get_fp_ene(*outf),
                              mtop, ir, mdoutf_get_fp_dhdl(*outf));
    }

    if (ir->bAdress)
    {
        please_cite(fplog, "Fritsch12");
        please_cite(fplog, "Junghans10");
    }
    /* Initiate variables */
    clear_mat(force_vir);
    clear_mat(shake_vir);
    clear_rvec(mu_tot);

    debug_gmx();
}
