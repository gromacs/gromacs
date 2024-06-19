/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "disre.h"

#include "config.h"

#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

void init_disres(FILE*                 fplog,
                 const gmx_mtop_t&     mtop,
                 t_inputrec*           ir,
                 DisResRunMode         disResRunMode,
                 DDRole                ddRole,
                 NumRanks              numRanks,
                 MPI_Comm              communicator,
                 const gmx_multisim_t* ms,
                 t_disresdata*         dd,
                 t_state*              state,
                 gmx_bool              bIsREMD)
{
    history_t* hist;
    char*      ptr;
    int        type_min, type_max;

    if (gmx_mtop_ftype_count(mtop, F_DISRES) == 0)
    {
        dd->nres = 0;

        return;
    }

    if (fplog)
    {
        fprintf(fplog, "Initializing the distance restraints\n");
    }

    dd->dr_weighting = ir->eDisreWeighting;
    dd->dr_fc        = ir->dr_fc;
    if (EI_DYNAMICS(ir->eI))
    {
        dd->dr_tau = ir->dr_tau;
    }
    else
    {
        dd->dr_tau = 0.0;
    }
    if (dd->dr_tau == 0.0)
    {
        dd->dr_bMixed = FALSE;
        dd->ETerm     = 0.0;
    }
    else
    {
        /* We store the r^-6 time averages in an array that is indexed
         * with the local disres iatom index, so this doesn't work with DD.
         * Note that DD is not initialized yet here, so we check that we are on multiple ranks,
         * but there are probably also issues with e.g. NM MPI parallelization.
         */
        if ((disResRunMode == DisResRunMode::MDRun) && (numRanks == NumRanks::Multiple))
        {
            gmx_fatal(FARGS,
                      "Time-averaged distance restraints are not supported with MPI "
                      "parallelization. You can use OpenMP parallelization on a single node.");
        }

        dd->dr_bMixed = ir->bDisreMixed;
        dd->ETerm     = std::exp(-(ir->delta_t / ir->dr_tau));
    }
    dd->ETerm1 = 1.0 - dd->ETerm;

    dd->nres  = 0;
    dd->npair = 0;
    type_min  = INT_MAX;
    type_max  = 0;
    for (const auto il : IListRange(mtop))
    {
        if (il.nmol() > 1 && !il.list()[F_DISRES].empty()
            && ir->eDisre != DistanceRestraintRefinement::Ensemble)
        {
            gmx_fatal(FARGS,
                      "NMR distance restraints with multiple copies of the same molecule are "
                      "currently only supported with ensemble averaging. If you just want to "
                      "restrain distances between atom pairs using a flat-bottomed potential, use "
                      "a restraint potential (bonds type 10) instead.");
        }

        int np = 0;
        for (int fa = 0; fa < il.list()[F_DISRES].size(); fa += 3)
        {
            int type = il.list()[F_DISRES].iatoms[fa];

            np++;
            int npair = mtop.ffparams.iparams[type].disres.npair;
            if (np == npair)
            {
                dd->nres += (ir->eDisre == DistanceRestraintRefinement::Ensemble ? 1 : il.nmol());
                dd->npair += il.nmol() * npair;
                np = 0;

                type_min = std::min(type_min, type);
                type_max = std::max(type_max, type);
            }
        }
    }

    if ((disResRunMode == DisResRunMode::MDRun) && (numRanks == NumRanks::Multiple) && ir->nstdisreout > 0)
    {
        /* With DD we currently only have local pair information available */
        gmx_fatal(FARGS,
                  "With MPI parallelization distance-restraint pair output is not supported. Use "
                  "nstdisreout=0 or use OpenMP parallelization on a single node.");
    }

    /* For communicating and/or reducing (sums of) r^-6 for pairs over threads
     * we use multiple arrays in t_disresdata. We need to have unique indices
     * for each restraint that work over threads and MPI ranks. To this end
     * we use the type index. These should all be in one block and there should
     * be dd->nres types, but we check for this here.
     * This setup currently does not allow for multiple copies of the same
     * molecule without ensemble averaging, this is check for above.
     */
    GMX_RELEASE_ASSERT(
            type_max - type_min + 1 == dd->nres,
            "All distance restraint parameter entries in the topology should be consecutive");

    dd->type_min = type_min;

    snew(dd->rt, dd->npair);

    if (dd->dr_tau != 0.0)
    {
        GMX_RELEASE_ASSERT(state != nullptr,
                           "We need a valid state when using time-averaged distance restraints");

        hist = &state->hist;
        /* Set the "history lack" factor to 1 */
        state->addEntry(StateEntry::DisreInitF);
        hist->disre_initf = 1.0;
        /* Allocate space for the r^-3 time averages */
        state->addEntry(StateEntry::DisreRm3Tav);
        hist->disre_rm3tav.resize(dd->npair);
    }
    /* Allocate space for a copy of rm3tav,
     * so we can call do_force without modifying the state.
     */
    snew(dd->rm3tav, dd->npair);

    /* Allocate Rt_6 and Rtav_6 consecutively in memory so they can be
     * averaged over the processors in one call (in calc_disre_R_6)
     */
    snew(dd->Rt_6, 2 * dd->nres);
    dd->Rtav_6 = &(dd->Rt_6[dd->nres]);

    ptr = getenv("GMX_DISRE_ENSEMBLE_SIZE");
    if ((disResRunMode == DisResRunMode::MDRun) && ms != nullptr && ptr != nullptr && !bIsREMD)
    {
#if GMX_MPI
        dd->nsystems = 0;
        sscanf(ptr, "%d", &dd->nsystems);
        if (fplog)
        {
            fprintf(fplog, "Found GMX_DISRE_ENSEMBLE_SIZE set to %d systems per ensemble\n", dd->nsystems);
        }
        /* This check is only valid on MAIN(cr), so probably
         * ensemble-averaged distance restraints are broken on more
         * than one processor per simulation system. */
        if (ddRole == DDRole::Main)
        {
            check_multi_int(fplog, ms, dd->nsystems, "the number of systems per ensemble", FALSE);
        }
        gmx_bcast(sizeof(int), &dd->nsystems, communicator);

        /* We use to allow any value of nsystems which was a divisor
         * of ms->numSimulations_. But this required an extra communicator which
         * pulled in mpi.h in nearly all C files.
         */
        if (!(ms->numSimulations_ == 1 || ms->numSimulations_ == dd->nsystems))
        {
            gmx_fatal(FARGS,
                      "GMX_DISRE_ENSEMBLE_SIZE (%d) is not equal to 1 or the number of systems "
                      "(option -multidir) %d",
                      dd->nsystems,
                      ms->numSimulations_);
        }
        if (fplog)
        {
            fprintf(fplog, "Our ensemble consists of systems:");
            for (int i = 0; i < dd->nsystems; i++)
            {
                fprintf(fplog, " %d", (ms->simulationIndex_ / dd->nsystems) * dd->nsystems + i);
            }
            fprintf(fplog, "\n");
        }
#else
        GMX_UNUSED_VALUE(communicator);
#endif
    }
    else
    {
        dd->nsystems = 1;
    }

    if (dd->nsystems == 1)
    {
        dd->Rtl_6 = dd->Rt_6;
    }
    else
    {
        snew(dd->Rtl_6, dd->nres);
    }

    if (dd->npair > 0)
    {
        if (fplog)
        {
            fprintf(fplog, "There are %d distance restraints involving %d atom pairs\n", dd->nres, dd->npair);
        }
        /* Have to avoid g_disre de-referencing cr blindly, mdrun not
         * doing consistency checks for ensemble-averaged distance
         * restraints when that's not happening, and only doing those
         * checks from appropriate processes (since check_multi_int is
         * too broken to check whether the communication will
         * succeed...) */
        if ((disResRunMode == DisResRunMode::MDRun) && ms && dd->nsystems > 1 && (ddRole == DDRole::Main))
        {
            check_multi_int(fplog, ms, dd->nres, "the number of distance restraints", FALSE);
        }
        please_cite(fplog, "Tropp80a");
        please_cite(fplog, "Torda89a");
    }
}

void calc_disres_R_6(const t_commrec*      cr,
                     const gmx_multisim_t* ms,
                     int                   nfa,
                     const t_iatom         forceatoms[],
                     const rvec            x[],
                     const t_pbc*          pbc,
                     t_disresdata*         dd,
                     const history_t*      hist)
{
    rvec     dx;
    real *   rt, *rm3tav, *Rtl_6, *Rt_6, *Rtav_6;
    real     ETerm, ETerm1, cf1 = 0, cf2 = 0;
    gmx_bool bTav;

    bTav   = (dd->dr_tau != 0);
    ETerm  = dd->ETerm;
    ETerm1 = dd->ETerm1;
    rt     = dd->rt;
    rm3tav = dd->rm3tav;
    Rtl_6  = dd->Rtl_6;
    Rt_6   = dd->Rt_6;
    Rtav_6 = dd->Rtav_6;

    if (bTav)
    {
        /* scaling factor to smoothly turn on the restraint forces *
         * when using time averaging                               */
        dd->exp_min_t_tau = hist->disre_initf * ETerm;

        cf1 = dd->exp_min_t_tau;
        cf2 = 1.0 / (1.0 - dd->exp_min_t_tau);
    }

    for (int res = 0; res < dd->nres; res++)
    {
        Rtav_6[res] = 0.0;
        Rt_6[res]   = 0.0;
    }

    /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
     * the total number of atoms pairs is nfa/3                          */
    for (int fa = 0; fa < nfa; fa += 3)
    {
        int type = forceatoms[fa];
        int res  = type - dd->type_min;
        int pair = fa / 3;
        int ai   = forceatoms[fa + 1];
        int aj   = forceatoms[fa + 2];

        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            rvec_sub(x[ai], x[aj], dx);
        }
        real rt2  = iprod(dx, dx);
        real rt_1 = gmx::invsqrt(rt2);
        real rt_3 = rt_1 * rt_1 * rt_1;

        rt[pair] = rt2 * rt_1;
        if (bTav)
        {
            /* Here we update rm3tav in t_disresdata using the data
             * in history_t.
             * Thus the results stay correct when this routine
             * is called multiple times.
             */
            rm3tav[pair] = cf2 * ((ETerm - cf1) * hist->disre_rm3tav[pair] + ETerm1 * rt_3);
        }
        else
        {
            rm3tav[pair] = rt_3;
        }

        /* Currently divide_bondeds_over_threads() ensures that pairs within
         * the same restraint get assigned to the same thread, so we could
         * run this loop thread-parallel.
         */
        Rt_6[res] += rt_3 * rt_3;
        Rtav_6[res] += rm3tav[pair] * rm3tav[pair];
    }

    /* NOTE: Rt_6 and Rtav_6 are stored consecutively in memory */
    if (cr && haveDDAtomOrdering(*cr))
    {
        gmx_sum(2 * dd->nres, dd->Rt_6, cr);
    }

    if (dd->nsystems > 1)
    {
        real invn = 1.0 / dd->nsystems;

        for (int res = 0; res < dd->nres; res++)
        {
            Rtl_6[res] = Rt_6[res];
            Rt_6[res] *= invn;
            Rtav_6[res] *= invn;
        }

        GMX_ASSERT(cr != nullptr && ms != nullptr, "We need multisim with nsystems>1");
        gmx_sum_sim(2 * dd->nres, dd->Rt_6, ms);

        if (haveDDAtomOrdering(*cr))
        {
            gmx_bcast(2 * dd->nres, dd->Rt_6, cr->mpi_comm_mygroup);
        }
    }

    /* Store the base forceatoms pointer, so we can re-calculate the pair
     * index in ta_disres() for indexing pair data in t_disresdata when
     * using thread parallelization.
     */
    dd->forceatomsStart = forceatoms;

    dd->sumviol = 0;
}

real ta_disres(int              nfa,
               const t_iatom*   forceatoms,
               const t_iparams* ip,
               const rvec*      x,
               rvec4*           f,
               rvec*            fshift,
               const t_pbc*     pbc,
               real gmx_unused  lambda,
               real gmx_unused* dvdlambda,
               gmx::ArrayRef<const real> /*charge*/,
               t_fcdata gmx_unused* fcd,
               t_disresdata*        disresdata,
               t_oriresdata gmx_unused* oriresdata,
               int gmx_unused* global_atom_index)
{
    const real seven_three = 7.0 / 3.0;

    rvec     dx;
    real     weight_rt_1;
    real     smooth_fc, Rt, Rtav, rt2, *Rtl_6, *Rt_6, *Rtav_6;
    real     k0, f_scal = 0, fmax_scal, fk_scal, fij;
    real     tav_viol, instant_viol, mixed_viol, violtot, vtot;
    real     tav_viol_Rtav7, instant_viol_Rtav7;
    real     up1, up2, low;
    gmx_bool bConservative, bMixed, bViolation;
    gmx_bool dr_bMixed;

    DistanceRestraintWeighting dr_weighting = disresdata->dr_weighting;
    dr_bMixed                               = disresdata->dr_bMixed;
    Rtl_6                                   = disresdata->Rtl_6;
    Rt_6                                    = disresdata->Rt_6;
    Rtav_6                                  = disresdata->Rtav_6;

    tav_viol = instant_viol = mixed_viol = tav_viol_Rtav7 = instant_viol_Rtav7 = 0;

    smooth_fc = disresdata->dr_fc;
    if (disresdata->dr_tau != 0)
    {
        /* scaling factor to smoothly turn on the restraint forces *
         * when using time averaging                               */
        smooth_fc *= (1.0 - disresdata->exp_min_t_tau);
    }

    violtot = 0;
    vtot    = 0;

    /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
     * the total number of atoms pairs is nfa/3                          */
    int faOffset = static_cast<int>(forceatoms - disresdata->forceatomsStart);
    for (int fa = 0; fa < nfa; fa += 3)
    {
        int type  = forceatoms[fa];
        int npair = ip[type].disres.npair;
        up1       = ip[type].disres.up1;
        up2       = ip[type].disres.up2;
        low       = ip[type].disres.low;
        k0        = smooth_fc * ip[type].disres.kfac;

        int res = type - disresdata->type_min;

        /* save some flops when there is only one pair */
        if (ip[type].disres.type != 2)
        {
            bConservative = (dr_weighting == DistanceRestraintWeighting::Conservative) && (npair > 1);
            bMixed        = dr_bMixed;
            Rt            = gmx::invsixthroot(Rt_6[res]);
            Rtav          = gmx::invsixthroot(Rtav_6[res]);
        }
        else
        {
            /* When rtype=2 use instantaneous not ensemble averaged distance */
            bConservative = (npair > 1);
            bMixed        = FALSE;
            Rt            = gmx::invsixthroot(Rtl_6[res]);
            Rtav          = Rt;
        }

        if (Rtav > up1)
        {
            bViolation = TRUE;
            tav_viol   = Rtav - up1;
        }
        else if (Rtav < low)
        {
            bViolation = TRUE;
            tav_viol   = Rtav - low;
        }
        else
        {
            bViolation = FALSE;
        }

        if (bViolation)
        {
            /* Add 1/npair energy and violation for each of the npair pairs */
            real pairFac = 1 / static_cast<real>(npair);

            /* NOTE:
             * there is no real potential when time averaging is applied
             */
            if (tav_viol > up2 - up1)
            {
                vtot += 0.5 * k0 * (up2 - up1) * (2 * tav_viol + up1 - up2) * pairFac;
            }
            else
            {
                vtot += 0.5 * k0 * gmx::square(tav_viol) * pairFac;
            }
            if (!bMixed)
            {
                f_scal = -k0 * tav_viol;
                violtot += std::fabs(tav_viol) * pairFac;
            }
            else
            {
                if (Rt > up1)
                {
                    if (tav_viol > 0)
                    {
                        instant_viol = Rt - up1;
                    }
                    else
                    {
                        bViolation = FALSE;
                    }
                }
                else if (Rt < low)
                {
                    if (tav_viol < 0)
                    {
                        instant_viol = Rt - low;
                    }
                    else
                    {
                        bViolation = FALSE;
                    }
                }
                else
                {
                    bViolation = FALSE;
                }
                if (bViolation)
                {
                    mixed_viol = std::sqrt(tav_viol * instant_viol);
                    f_scal     = -k0 * mixed_viol;
                    violtot += mixed_viol * pairFac;
                }
            }
        }

        if (bViolation)
        {
            fmax_scal = -k0 * (up2 - up1);
            /* Correct the force for the number of restraints */
            if (bConservative)
            {
                if ((k0 != 0) && (-f_scal / k0 > up2 - up1))
                {
                    f_scal = fmax_scal;
                }
                if (!bMixed)
                {
                    f_scal *= Rtav / Rtav_6[res];
                }
                else
                {
                    f_scal /= 2 * mixed_viol;
                    tav_viol_Rtav7     = tav_viol * Rtav / Rtav_6[res];
                    instant_viol_Rtav7 = instant_viol * Rt / Rt_6[res];
                }
            }
            else
            {
                f_scal /= npair;
                if ((k0 != 0) && (-f_scal / k0 > up2 - up1))
                {
                    f_scal = fmax_scal;
                }
            }

            /* Exert the force ... */

            int pair = (faOffset + fa) / 3;
            int ai   = forceatoms[fa + 1];
            int aj   = forceatoms[fa + 2];
            int ki   = gmx::c_centralShiftIndex;
            if (pbc)
            {
                ki = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
            }
            else
            {
                rvec_sub(x[ai], x[aj], dx);
            }
            rt2 = iprod(dx, dx);

            weight_rt_1 = gmx::invsqrt(rt2);

            if (bConservative)
            {
                if (!dr_bMixed)
                {
                    weight_rt_1 *= std::pow(disresdata->rm3tav[pair], seven_three);
                }
                else
                {
                    weight_rt_1 *=
                            tav_viol_Rtav7 * std::pow(disresdata->rm3tav[pair], seven_three)
                            + instant_viol_Rtav7
                                      / (disresdata->rt[pair] * gmx::power6(disresdata->rt[pair]));
                }
            }

            fk_scal = f_scal * weight_rt_1;

            for (int m = 0; m < DIM; m++)
            {
                fij = fk_scal * dx[m];

                f[ai][m] += fij;
                f[aj][m] -= fij;
                if (fshift)
                {
                    fshift[ki][m] += fij;
                    fshift[gmx::c_centralShiftIndex][m] -= fij;
                }
            }
        }
    }

#pragma omp atomic
    disresdata->sumviol += violtot;

    /* Return energy */
    return vtot;
}

void update_disres_history(const t_disresdata& dd, history_t* hist)
{
    if (dd.dr_tau != 0)
    {
        /* Copy the new time averages that have been calculated
         * in calc_disres_R_6.
         */
        hist->disre_initf = dd.exp_min_t_tau;
        for (int pair = 0; pair < dd.npair; pair++)
        {
            hist->disre_rm3tav[pair] = dd.rm3tav[pair];
        }
    }
}
