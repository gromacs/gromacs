/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "gromacs/legacyheaders/disre.h"

#include "config.h"

#include <math.h>
#include <stdlib.h>

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

void init_disres(FILE *fplog, const gmx_mtop_t *mtop,
                 t_inputrec *ir, const t_commrec *cr,
                 t_fcdata *fcd, t_state *state, gmx_bool bIsREMD)
{
    int                  fa, nmol, i, npair, np;
    t_iparams           *ip;
    t_disresdata        *dd;
    history_t           *hist;
    gmx_mtop_ilistloop_t iloop;
    t_ilist             *il;
    char                *ptr;

    dd = &(fcd->disres);

    if (gmx_mtop_ftype_count(mtop, F_DISRES) == 0)
    {
        dd->nres = 0;

        return;
    }

    if (fplog)
    {
        fprintf(fplog, "Initializing the distance restraints\n");
    }


    if (ir->eDisre == edrEnsemble)
    {
        gmx_fatal(FARGS, "Sorry, distance restraints with ensemble averaging over multiple molecules in one system are not functional in this version of GROMACS");
    }

    dd->dr_weighting = ir->eDisreWeighting;
    dd->dr_fc        = ir->dr_fc;
    if (EI_DYNAMICS(ir->eI))
    {
        dd->dr_tau   = ir->dr_tau;
    }
    else
    {
        dd->dr_tau   = 0.0;
    }
    if (dd->dr_tau == 0.0)
    {
        dd->dr_bMixed = FALSE;
        dd->ETerm     = 0.0;
    }
    else
    {
        dd->dr_bMixed = ir->bDisreMixed;
        dd->ETerm     = exp(-(ir->delta_t/ir->dr_tau));
    }
    dd->ETerm1        = 1.0 - dd->ETerm;

    ip = mtop->ffparams.iparams;

    dd->nres  = 0;
    dd->npair = 0;
    iloop     = gmx_mtop_ilistloop_init(mtop);
    while (gmx_mtop_ilistloop_next(iloop, &il, &nmol))
    {
        np = 0;
        for (fa = 0; fa < il[F_DISRES].nr; fa += 3)
        {
            np++;
            npair = mtop->ffparams.iparams[il[F_DISRES].iatoms[fa]].disres.npair;
            if (np == npair)
            {
                dd->nres  += (ir->eDisre == edrEnsemble ? 1 : nmol)*npair;
                dd->npair += nmol*npair;
                np         = 0;
            }
        }
    }

    if (cr && PAR(cr))
    {
        /* Temporary check, will be removed when disre is implemented with DD */
        const char *notestr = "NOTE: atoms involved in distance restraints should be within the same domain. If this is not the case mdrun generates a fatal error. If you encounter this, use a single MPI rank (Verlet+OpenMP+GPUs work fine).";

        if (MASTER(cr))
        {
            fprintf(stderr, "\n%s\n\n", notestr);
        }
        if (fplog)
        {
            fprintf(fplog, "%s\n", notestr);
        }

        if (dd->dr_tau != 0 || ir->eDisre == edrEnsemble || cr->ms != NULL ||
            dd->nres != dd->npair)
        {
            gmx_fatal(FARGS, "Time or ensemble averaged or multiple pair distance restraints do not work (yet) with domain decomposition, use a single MPI rank%s", cr->ms ? " per simulation" : "");
        }
        if (ir->nstdisreout != 0)
        {
            if (fplog)
            {
                fprintf(fplog, "\nWARNING: Can not write distance restraint data to energy file with domain decomposition\n\n");
            }
            if (MASTER(cr))
            {
                fprintf(stderr, "\nWARNING: Can not write distance restraint data to energy file with domain decomposition\n");
            }
            ir->nstdisreout = 0;
        }
    }

    snew(dd->rt, dd->npair);

    if (dd->dr_tau != 0.0)
    {
        hist = &state->hist;
        /* Set the "history lack" factor to 1 */
        state->flags     |= (1<<estDISRE_INITF);
        hist->disre_initf = 1.0;
        /* Allocate space for the r^-3 time averages */
        state->flags     |= (1<<estDISRE_RM3TAV);
        hist->ndisrepairs = dd->npair;
        snew(hist->disre_rm3tav, hist->ndisrepairs);
    }
    /* Allocate space for a copy of rm3tav,
     * so we can call do_force without modifying the state.
     */
    snew(dd->rm3tav, dd->npair);

    /* Allocate Rt_6 and Rtav_6 consecutively in memory so they can be
     * averaged over the processors in one call (in calc_disre_R_6)
     */
    snew(dd->Rt_6, 2*dd->nres);
    dd->Rtav_6 = &(dd->Rt_6[dd->nres]);

    ptr = getenv("GMX_DISRE_ENSEMBLE_SIZE");
    if (cr && cr->ms != NULL && ptr != NULL && !bIsREMD)
    {
#ifdef GMX_MPI
        dd->nsystems = 0;
        sscanf(ptr, "%d", &dd->nsystems);
        if (fplog)
        {
            fprintf(fplog, "Found GMX_DISRE_ENSEMBLE_SIZE set to %d systems per ensemble\n", dd->nsystems);
        }
        /* This check is only valid on MASTER(cr), so probably
         * ensemble-averaged distance restraints are broken on more
         * than one processor per simulation system. */
        if (MASTER(cr))
        {
            check_multi_int(fplog, cr->ms, dd->nsystems,
                            "the number of systems per ensemble",
                            FALSE);
        }
        gmx_bcast_sim(sizeof(int), &dd->nsystems, cr);

        /* We use to allow any value of nsystems which was a divisor
         * of ms->nsim. But this required an extra communicator which
         * was stored in t_fcdata. This pulled in mpi.h in nearly all C files.
         */
        if (!(cr->ms->nsim == 1 || cr->ms->nsim == dd->nsystems))
        {
            gmx_fatal(FARGS, "GMX_DISRE_ENSEMBLE_SIZE (%d) is not equal to 1 or the number of systems (option -multi) %d", dd->nsystems, cr->ms->nsim);
        }
        if (fplog)
        {
            fprintf(fplog, "Our ensemble consists of systems:");
            for (i = 0; i < dd->nsystems; i++)
            {
                fprintf(fplog, " %d",
                        (cr->ms->sim/dd->nsystems)*dd->nsystems+i);
            }
            fprintf(fplog, "\n");
        }
        snew(dd->Rtl_6, dd->nres);
#endif
    }
    else
    {
        dd->nsystems = 1;
        dd->Rtl_6    = dd->Rt_6;
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
        if (cr && cr->ms && dd->nsystems > 1 && MASTER(cr))
        {
            check_multi_int(fplog, cr->ms, fcd->disres.nres,
                            "the number of distance restraints",
                            FALSE);
        }
        please_cite(fplog, "Tropp80a");
        please_cite(fplog, "Torda89a");
    }
}

void calc_disres_R_6(int nfa, const t_iatom forceatoms[], const t_iparams ip[],
                     const rvec x[], const t_pbc *pbc,
                     t_fcdata *fcd, history_t *hist)
{
    atom_id         ai, aj;
    int             fa, res, i, pair, ki, kj, m;
    int             type, npair, np;
    rvec            dx;
    real           *rt, *rm3tav, *Rtl_6, *Rt_6, *Rtav_6;
    real            rt_1, rt_3, rt2;
    ivec            it, jt, dt;
    t_disresdata   *dd;
    real            ETerm, ETerm1, cf1 = 0, cf2 = 0, invn = 0;
    gmx_bool        bTav;

    dd           = &(fcd->disres);
    bTav         = (dd->dr_tau != 0);
    ETerm        = dd->ETerm;
    ETerm1       = dd->ETerm1;
    rt           = dd->rt;
    rm3tav       = dd->rm3tav;
    Rtl_6        = dd->Rtl_6;
    Rt_6         = dd->Rt_6;
    Rtav_6       = dd->Rtav_6;

    if (bTav)
    {
        /* scaling factor to smoothly turn on the restraint forces *
         * when using time averaging                               */
        dd->exp_min_t_tau = hist->disre_initf*ETerm;

        cf1 = dd->exp_min_t_tau;
        cf2 = 1.0/(1.0 - dd->exp_min_t_tau);
    }

    if (dd->nsystems > 1)
    {
        invn = 1.0/dd->nsystems;
    }

    /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
     * the total number of atoms pairs is nfa/3                          */
    res = 0;
    fa  = 0;
    while (fa < nfa)
    {
        type  = forceatoms[fa];
        npair = ip[type].disres.npair;

        Rtav_6[res] = 0.0;
        Rt_6[res]   = 0.0;

        /* Loop over the atom pairs of 'this' restraint */
        np = 0;
        while (fa < nfa && np < npair)
        {
            pair = fa/3;
            ai   = forceatoms[fa+1];
            aj   = forceatoms[fa+2];

            if (pbc)
            {
                pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
            }
            else
            {
                rvec_sub(x[ai], x[aj], dx);
            }
            rt2  = iprod(dx, dx);
            rt_1 = gmx_invsqrt(rt2);
            rt_3 = rt_1*rt_1*rt_1;

            rt[pair]         = sqrt(rt2);
            if (bTav)
            {
                /* Here we update rm3tav in t_fcdata using the data
                 * in history_t.
                 * Thus the results stay correct when this routine
                 * is called multiple times.
                 */
                rm3tav[pair] = cf2*((ETerm - cf1)*hist->disre_rm3tav[pair] +
                                    ETerm1*rt_3);
            }
            else
            {
                rm3tav[pair] = rt_3;
            }

            Rt_6[res]       += rt_3*rt_3;
            Rtav_6[res]     += rm3tav[pair]*rm3tav[pair];

            fa += 3;
            np++;
        }
        if (dd->nsystems > 1)
        {
            Rtl_6[res]   = Rt_6[res];
            Rt_6[res]   *= invn;
            Rtav_6[res] *= invn;
        }

        res++;
    }
}

real ta_disres(int nfa, const t_iatom forceatoms[], const t_iparams ip[],
               const rvec x[], rvec f[], rvec fshift[],
               const t_pbc *pbc, const t_graph *g,
               real gmx_unused lambda, real gmx_unused *dvdlambda,
               const t_mdatoms gmx_unused *md, t_fcdata *fcd,
               int gmx_unused *global_atom_index)
{
    const real      sixth       = 1.0/6.0;
    const real      seven_three = 7.0/3.0;

    atom_id         ai, aj;
    int             fa, res, npair, p, pair, ki = CENTRAL, m;
    int             type;
    rvec            dx;
    real            weight_rt_1;
    real            smooth_fc, Rt, Rtav, rt2, *Rtl_6, *Rt_6, *Rtav_6;
    real            k0, f_scal = 0, fmax_scal, fk_scal, fij;
    real            tav_viol, instant_viol, mixed_viol, violtot, vtot;
    real            tav_viol_Rtav7, instant_viol_Rtav7;
    real            up1, up2, low;
    gmx_bool        bConservative, bMixed, bViolation;
    ivec            it, jt, dt;
    t_disresdata   *dd;
    int             dr_weighting;
    gmx_bool        dr_bMixed;

    dd           = &(fcd->disres);
    dr_weighting = dd->dr_weighting;
    dr_bMixed    = dd->dr_bMixed;
    Rtl_6        = dd->Rtl_6;
    Rt_6         = dd->Rt_6;
    Rtav_6       = dd->Rtav_6;

    tav_viol = instant_viol = mixed_viol = tav_viol_Rtav7 = instant_viol_Rtav7 = 0;

    smooth_fc = dd->dr_fc;
    if (dd->dr_tau != 0)
    {
        /* scaling factor to smoothly turn on the restraint forces *
         * when using time averaging                               */
        smooth_fc *= (1.0 - dd->exp_min_t_tau);
    }

    violtot = 0;
    vtot    = 0;

    /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
     * the total number of atoms pairs is nfa/3                          */
    res  = 0;
    fa   = 0;
    while (fa < nfa)
    {
        type  = forceatoms[fa];
        /* Take action depending on restraint, calculate scalar force */
        npair = ip[type].disres.npair;
        up1   = ip[type].disres.up1;
        up2   = ip[type].disres.up2;
        low   = ip[type].disres.low;
        k0    = smooth_fc*ip[type].disres.kfac;

        /* save some flops when there is only one pair */
        if (ip[type].disres.type != 2)
        {
            bConservative = (dr_weighting == edrwConservative) && (npair > 1);
            bMixed        = dr_bMixed;
            Rt            = pow(Rt_6[res], -sixth);
            Rtav          = pow(Rtav_6[res], -sixth);
        }
        else
        {
            /* When rtype=2 use instantaneous not ensemble avereged distance */
            bConservative = (npair > 1);
            bMixed        = FALSE;
            Rt            = pow(Rtl_6[res], -sixth);
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
            /* NOTE:
             * there is no real potential when time averaging is applied
             */
            vtot += 0.5*k0*sqr(tav_viol);
            if (1/vtot == 0)
            {
                printf("vtot is inf: %f\n", vtot);
            }
            if (!bMixed)
            {
                f_scal   = -k0*tav_viol;
                violtot += fabs(tav_viol);
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
                    mixed_viol = sqrt(tav_viol*instant_viol);
                    f_scal     = -k0*mixed_viol;
                    violtot   += mixed_viol;
                }
            }
        }

        if (bViolation)
        {
            fmax_scal = -k0*(up2-up1);
            /* Correct the force for the number of restraints */
            if (bConservative)
            {
                f_scal  = max(f_scal, fmax_scal);
                if (!bMixed)
                {
                    f_scal *= Rtav/Rtav_6[res];
                }
                else
                {
                    f_scal            /= 2*mixed_viol;
                    tav_viol_Rtav7     = tav_viol*Rtav/Rtav_6[res];
                    instant_viol_Rtav7 = instant_viol*Rt/Rt_6[res];
                }
            }
            else
            {
                f_scal /= (real)npair;
                f_scal  = max(f_scal, fmax_scal);
            }

            /* Exert the force ... */

            /* Loop over the atom pairs of 'this' restraint */
            for (p = 0; p < npair; p++)
            {
                pair = fa/3;
                ai   = forceatoms[fa+1];
                aj   = forceatoms[fa+2];

                if (pbc)
                {
                    ki = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
                }
                else
                {
                    rvec_sub(x[ai], x[aj], dx);
                }
                rt2 = iprod(dx, dx);

                weight_rt_1 = gmx_invsqrt(rt2);

                if (bConservative)
                {
                    if (!dr_bMixed)
                    {
                        weight_rt_1 *= pow(dd->rm3tav[pair], seven_three);
                    }
                    else
                    {
                        weight_rt_1 *= tav_viol_Rtav7*pow(dd->rm3tav[pair], seven_three)+
                            instant_viol_Rtav7*pow(dd->rt[pair], -7);
                    }
                }

                fk_scal  = f_scal*weight_rt_1;

                if (g)
                {
                    ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
                    ki = IVEC2IS(dt);
                }

                for (m = 0; m < DIM; m++)
                {
                    fij            = fk_scal*dx[m];

                    f[ai][m]           += fij;
                    f[aj][m]           -= fij;
                    fshift[ki][m]      += fij;
                    fshift[CENTRAL][m] -= fij;
                }
                fa += 3;
            }
        }
        else
        {
            /* No violation so force and potential contributions */
            fa += 3*npair;
        }
        res++;
    }

    dd->sumviol = violtot;

    /* Return energy */
    return vtot;
}

void update_disres_history(t_fcdata *fcd, history_t *hist)
{
    t_disresdata *dd;
    int           pair;

    dd = &(fcd->disres);
    if (dd->dr_tau != 0)
    {
        /* Copy the new time averages that have been calculated
         * in calc_disres_R_6.
         */
        hist->disre_initf = dd->exp_min_t_tau;
        for (pair = 0; pair < dd->npair; pair++)
        {
            hist->disre_rm3tav[pair] = dd->rm3tav[pair];
        }
    }
}
