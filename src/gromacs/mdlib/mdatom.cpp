/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include <math.h>

#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/qmmm.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/smalloc.h"

#define ALMOST_ZERO 1e-30

t_mdatoms *init_mdatoms(FILE *fp, gmx_mtop_t *mtop, gmx_bool bFreeEnergy)
{
    int                     a;
    double                  tmA, tmB;
    t_atom                 *atom;
    t_mdatoms              *md;
    gmx_mtop_atomloop_all_t aloop;

    snew(md, 1);

    md->nenergrp = mtop->groups.grps[egcENER].nr;
    md->bVCMgrps = FALSE;
    tmA          = 0.0;
    tmB          = 0.0;

    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &a, &atom))
    {
        if (ggrpnr(&mtop->groups, egcVCM, a) > 0)
        {
            md->bVCMgrps = TRUE;
        }

        if (bFreeEnergy && PERTURBED(*atom))
        {
            md->nPerturbed++;
            if (atom->mB != atom->m)
            {
                md->nMassPerturbed++;
            }
            if (atom->qB != atom->q)
            {
                md->nChargePerturbed++;
            }
            if (atom->typeB != atom->type)
            {
                md->nTypePerturbed++;
            }
        }

        tmA += atom->m;
        tmB += atom->mB;
    }

    md->tmassA = tmA;
    md->tmassB = tmB;

    if (bFreeEnergy && fp)
    {
        fprintf(fp,
                "There are %d atoms and %d charges for free energy perturbation\n",
                md->nPerturbed, md->nChargePerturbed);
    }

    md->bOrires = gmx_mtop_ftype_count(mtop, F_ORIRES);

    return md;
}

void atoms2md(gmx_mtop_t *mtop, t_inputrec *ir,
              int nindex, int *index,
              int homenr,
              t_mdatoms *md)
{
    gmx_bool              bLJPME;
    gmx_mtop_atomlookup_t alook;
    int                   i;
    t_grpopts            *opts;
    gmx_groups_t         *groups;
    int                   nthreads gmx_unused;
    const real            oneOverSix = 1.0 / 6.0;

    bLJPME = EVDW_PME(ir->vdwtype);

    opts = &ir->opts;

    groups = &mtop->groups;

    /* Index==NULL indicates no DD (unless we have a DD node with no
     * atoms), so also check for homenr. This should be
     * signaled properly with an extra parameter or nindex==-1.
     */
    if (index == NULL && (homenr > 0))
    {
        md->nr = mtop->natoms;
    }
    else
    {
        md->nr = nindex;
    }

    if (md->nr > md->nalloc)
    {
        md->nalloc = over_alloc_dd(md->nr);

        if (md->nMassPerturbed)
        {
            srenew(md->massA, md->nalloc);
            srenew(md->massB, md->nalloc);
        }
        srenew(md->massT, md->nalloc);
        srenew(md->invmass, md->nalloc);
        srenew(md->chargeA, md->nalloc);
        srenew(md->typeA, md->nalloc);
        if (md->nPerturbed)
        {
            srenew(md->chargeB, md->nalloc);
            srenew(md->typeB, md->nalloc);
        }
        if (bLJPME)
        {
            srenew(md->sqrt_c6A, md->nalloc);
            srenew(md->sigmaA, md->nalloc);
            srenew(md->sigma3A, md->nalloc);
            if (md->nPerturbed)
            {
                srenew(md->sqrt_c6B, md->nalloc);
                srenew(md->sigmaB, md->nalloc);
                srenew(md->sigma3B, md->nalloc);
            }
        }
        srenew(md->ptype, md->nalloc);
        if (opts->ngtc > 1)
        {
            srenew(md->cTC, md->nalloc);
            /* We always copy cTC with domain decomposition */
        }
        srenew(md->cENER, md->nalloc);
        if (opts->ngacc > 1)
        {
            srenew(md->cACC, md->nalloc);
        }
        if (opts->nFreeze &&
            (opts->ngfrz > 1 ||
             opts->nFreeze[0][XX] || opts->nFreeze[0][YY] || opts->nFreeze[0][ZZ]))
        {
            srenew(md->cFREEZE, md->nalloc);
        }
        if (md->bVCMgrps)
        {
            srenew(md->cVCM, md->nalloc);
        }
        if (md->bOrires)
        {
            srenew(md->cORF, md->nalloc);
        }
        if (md->nPerturbed)
        {
            srenew(md->bPerturbed, md->nalloc);
        }

        /* Note that these user t_mdatoms array pointers are NULL
         * when there is only one group present.
         * Therefore, when adding code, the user should use something like:
         * gprnrU1 = (md->cU1==NULL ? 0 : md->cU1[localatindex])
         */
        if (mtop->groups.grpnr[egcUser1] != NULL)
        {
            srenew(md->cU1, md->nalloc);
        }
        if (mtop->groups.grpnr[egcUser2] != NULL)
        {
            srenew(md->cU2, md->nalloc);
        }

        if (ir->bQMMM)
        {
            srenew(md->bQM, md->nalloc);
        }
        if (ir->bAdress)
        {
            srenew(md->wf, md->nalloc);
            srenew(md->tf_table_index, md->nalloc);
        }
    }

    alook = gmx_mtop_atomlookup_init(mtop);

    // cppcheck-suppress unreadVariable
    nthreads = gmx_omp_nthreads_get(emntDefault);
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (i = 0; i < md->nr; i++)
    {
        int      g, ag;
        real     mA, mB, fac;
        real     c6, c12;
        t_atom  *atom;

        if (index == NULL)
        {
            ag = i;
        }
        else
        {
            ag   = index[i];
        }
        gmx_mtop_atomnr_to_atom(alook, ag, &atom);

        if (md->cFREEZE)
        {
            md->cFREEZE[i] = ggrpnr(groups, egcFREEZE, ag);
        }
        if (EI_ENERGY_MINIMIZATION(ir->eI))
        {
            /* Displacement is proportional to F, masses used for constraints */
            mA = 1.0;
            mB = 1.0;
        }
        else if (ir->eI == eiBD)
        {
            /* With BD the physical masses are irrelevant.
             * To keep the code simple we use most of the normal MD code path
             * for BD. Thus for constraining the masses should be proportional
             * to the friction coefficient. We set the absolute value such that
             * m/2<(dx/dt)^2> = m/2*2kT/fric*dt = kT/2 => m=fric*dt/2
             * Then if we set the (meaningless) velocity to v=dx/dt, we get the
             * correct kinetic energy and temperature using the usual code path.
             * Thus with BD v*dt will give the displacement and the reported
             * temperature can signal bad integration (too large time step).
             */
            if (ir->bd_fric > 0)
            {
                mA = 0.5*ir->bd_fric*ir->delta_t;
                mB = 0.5*ir->bd_fric*ir->delta_t;
            }
            else
            {
                /* The friction coefficient is mass/tau_t */
                fac = ir->delta_t/opts->tau_t[md->cTC ? groups->grpnr[egcTC][ag] : 0];
                mA  = 0.5*atom->m*fac;
                mB  = 0.5*atom->mB*fac;
            }
        }
        else
        {
            mA = atom->m;
            mB = atom->mB;
        }
        if (md->nMassPerturbed)
        {
            md->massA[i]  = mA;
            md->massB[i]  = mB;
        }
        md->massT[i]    = mA;
        if (mA == 0.0)
        {
            md->invmass[i]    = 0;
        }
        else if (md->cFREEZE)
        {
            g = md->cFREEZE[i];
            if (opts->nFreeze[g][XX] && opts->nFreeze[g][YY] && opts->nFreeze[g][ZZ])
            {
                /* Set the mass of completely frozen particles to ALMOST_ZERO iso 0
                 * to avoid div by zero in lincs or shake.
                 * Note that constraints can still move a partially frozen particle.
                 */
                md->invmass[i]  = ALMOST_ZERO;
            }
            else
            {
                md->invmass[i]  = 1.0/mA;
            }
        }
        else
        {
            md->invmass[i]    = 1.0/mA;
        }
        md->chargeA[i]      = atom->q;
        md->typeA[i]        = atom->type;
        if (bLJPME)
        {
            c6                = mtop->ffparams.iparams[atom->type*(mtop->ffparams.atnr+1)].lj.c6;
            c12               = mtop->ffparams.iparams[atom->type*(mtop->ffparams.atnr+1)].lj.c12;
            md->sqrt_c6A[i]   = sqrt(c6);
            if (c6 == 0.0 || c12 == 0)
            {
                md->sigmaA[i] = 1.0;
            }
            else
            {
                md->sigmaA[i] = pow(c12/c6, oneOverSix);
            }
            md->sigma3A[i]    = 1/(md->sigmaA[i]*md->sigmaA[i]*md->sigmaA[i]);
        }
        if (md->nPerturbed)
        {
            md->bPerturbed[i] = PERTURBED(*atom);
            md->chargeB[i]    = atom->qB;
            md->typeB[i]      = atom->typeB;
            if (bLJPME)
            {
                c6                = mtop->ffparams.iparams[atom->typeB*(mtop->ffparams.atnr+1)].lj.c6;
                c12               = mtop->ffparams.iparams[atom->typeB*(mtop->ffparams.atnr+1)].lj.c12;
                md->sqrt_c6B[i]   = sqrt(c6);
                if (c6 == 0.0 || c12 == 0)
                {
                    md->sigmaB[i] = 1.0;
                }
                else
                {
                    md->sigmaB[i] = pow(c12/c6, oneOverSix);
                }
                md->sigma3B[i]    = 1/(md->sigmaB[i]*md->sigmaB[i]*md->sigmaB[i]);
            }
        }
        md->ptype[i]    = atom->ptype;
        if (md->cTC)
        {
            md->cTC[i]    = groups->grpnr[egcTC][ag];
        }
        md->cENER[i]    =
            (groups->grpnr[egcENER] ? groups->grpnr[egcENER][ag] : 0);
        if (md->cACC)
        {
            md->cACC[i]   = groups->grpnr[egcACC][ag];
        }
        if (md->cVCM)
        {
            md->cVCM[i]       = groups->grpnr[egcVCM][ag];
        }
        if (md->cORF)
        {
            md->cORF[i]       = groups->grpnr[egcORFIT][ag];
        }

        if (md->cU1)
        {
            md->cU1[i]        = groups->grpnr[egcUser1][ag];
        }
        if (md->cU2)
        {
            md->cU2[i]        = groups->grpnr[egcUser2][ag];
        }

        if (ir->bQMMM)
        {
            if (groups->grpnr[egcQMMM] == 0 ||
                groups->grpnr[egcQMMM][ag] < groups->grps[egcQMMM].nr-1)
            {
                md->bQM[i]      = TRUE;
            }
            else
            {
                md->bQM[i]      = FALSE;
            }
        }
        /* Initialize AdResS weighting functions to adressw */
        if (ir->bAdress)
        {
            md->wf[i]           = 1.0;
            /* if no tf table groups specified, use default table */
            md->tf_table_index[i] = DEFAULT_TF_TABLE;
            if (ir->adress->n_tf_grps > 0)
            {
                /* if tf table groups specified, tf is only applied to thoose energy groups*/
                md->tf_table_index[i] = NO_TF_TABLE;
                /* check wether atom is in one of the relevant energy groups and assign a table index */
                for (g = 0; g < ir->adress->n_tf_grps; g++)
                {
                    if (md->cENER[i] == ir->adress->tf_table_index[g])
                    {
                        md->tf_table_index[i] = g;
                    }
                }
            }
        }
    }

    gmx_mtop_atomlookup_destroy(alook);

    md->homenr = homenr;
    md->lambda = 0;
}

void update_mdatoms(t_mdatoms *md, real lambda)
{
    int    al, end;
    real   L1 = 1.0-lambda;

    end = md->nr;

    if (md->nMassPerturbed)
    {
        for (al = 0; (al < end); al++)
        {
            if (md->bPerturbed[al])
            {
                md->massT[al] = L1*md->massA[al]+ lambda*md->massB[al];
                if (md->invmass[al] > 1.1*ALMOST_ZERO)
                {
                    md->invmass[al] = 1.0/md->massT[al];
                }
            }
        }
        md->tmass = L1*md->tmassA + lambda*md->tmassB;
    }
    else
    {
        md->tmass = md->tmassA;
    }
    md->lambda = lambda;
}
