/*
 * This file is part of the GROMACS molecular simulation package.
 *
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

#include "nbnxn_kernel_gpu_ref.h"

#include "config.h"

#include <math.h>

#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
#include "gromacs/pbcutil/ishift.h"

#define NCL_PER_SUPERCL         (NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER)
#define CL_SIZE                 (NBNXN_GPU_CLUSTER_SIZE)

void
nbnxn_kernel_gpu_ref(const nbnxn_pairlist_t     *nbl,
                     const nbnxn_atomdata_t     *nbat,
                     const interaction_const_t  *iconst,
                     rvec                       *shift_vec,
                     int                         force_flags,
                     int                         clearF,
                     real  *                     f,
                     real  *                     fshift,
                     real  *                     Vc,
                     real  *                     Vvdw)
{
    const nbnxn_sci_t  *nbln;
    const real         *x;
    gmx_bool            bEner;
    gmx_bool            bEwald;
    const real         *Ftab = NULL;
    real                rcut2, rvdw2, rlist2;
    int                 ntype;
    real                facel;
    int                 n;
    int                 ish3;
    int                 sci;
    int                 cj4_ind0, cj4_ind1, cj4_ind;
    int                 ci, cj;
    int                 ic, jc, ia, ja, is, ifs, js, jfs, im, jm;
    int                 n0;
    int                 ggid;
    real                shX, shY, shZ;
    real                fscal, tx, ty, tz;
    real                rinvsq;
    real                iq;
    real                qq, vcoul = 0, krsq, vctot;
    int                 nti;
    int                 tj;
    real                rt, r, eps;
    real                rinvsix;
    real                Vvdwtot;
    real                Vvdw_rep, Vvdw_disp;
    real                ix, iy, iz, fix, fiy, fiz;
    real                jx, jy, jz;
    real                dx, dy, dz, rsq, rinv;
    int                 int_bit;
    real                fexcl;
    real                c6, c12, cexp1, cexp2, br;
    const real       *  shiftvec;
    real       *        vdwparam;
    int       *         shift;
    int       *         type;
    const nbnxn_excl_t *excl[2];

    int                 npair_tot, npair;
    int                 nhwu, nhwu_pruned;

    if (nbl->na_ci != CL_SIZE)
    {
        gmx_fatal(FARGS, "The neighborlist cluster size in the GPU reference kernel is %d, expected it to be %d", nbl->na_ci, CL_SIZE);
    }

    if (clearF == enbvClearFYes)
    {
        clear_f(nbat, 0, f);
    }

    bEner = (force_flags & GMX_FORCE_ENERGY);

    bEwald = EEL_FULL(iconst->eeltype);
    if (bEwald)
    {
        Ftab = iconst->tabq_coul_F;
    }

    rcut2               = iconst->rcoulomb*iconst->rcoulomb;
    rvdw2               = iconst->rvdw*iconst->rvdw;

    rlist2              = nbl->rlist*nbl->rlist;

    type                = nbat->type;
    facel               = iconst->epsfac;
    shiftvec            = shift_vec[0];
    vdwparam            = nbat->nbfp;
    ntype               = nbat->ntype;

    x = nbat->x;

    npair_tot   = 0;
    nhwu        = 0;
    nhwu_pruned = 0;

    for (n = 0; n < nbl->nsci; n++)
    {
        nbln = &nbl->sci[n];

        ish3             = 3*nbln->shift;
        shX              = shiftvec[ish3];
        shY              = shiftvec[ish3+1];
        shZ              = shiftvec[ish3+2];
        cj4_ind0         = nbln->cj4_ind_start;
        cj4_ind1         = nbln->cj4_ind_end;
        sci              = nbln->sci;
        vctot            = 0;
        Vvdwtot          = 0;

        if (nbln->shift == CENTRAL &&
            nbl->cj4[cj4_ind0].cj[0] == sci*NCL_PER_SUPERCL)
        {
            /* we have the diagonal:
             * add the charge self interaction energy term
             */
            for (im = 0; im < NCL_PER_SUPERCL; im++)
            {
                ci = sci*NCL_PER_SUPERCL + im;
                for (ic = 0; ic < CL_SIZE; ic++)
                {
                    ia     = ci*CL_SIZE + ic;
                    iq     = x[ia*nbat->xstride+3];
                    vctot += iq*iq;
                }
            }
            if (!bEwald)
            {
                vctot *= -facel*0.5*iconst->c_rf;
            }
            else
            {
                /* last factor 1/sqrt(pi) */
                vctot *= -facel*iconst->ewaldcoeff_q*M_1_SQRTPI;
            }
        }

        for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++)
        {
            excl[0]           = &nbl->excl[nbl->cj4[cj4_ind].imei[0].excl_ind];
            excl[1]           = &nbl->excl[nbl->cj4[cj4_ind].imei[1].excl_ind];

            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++)
            {
                cj               = nbl->cj4[cj4_ind].cj[jm];

                for (im = 0; im < NCL_PER_SUPERCL; im++)
                {
                    /* We're only using the first imask,
                     * but here imei[1].imask is identical.
                     */
                    if ((nbl->cj4[cj4_ind].imei[0].imask >> (jm*NCL_PER_SUPERCL+im)) & 1)
                    {
                        gmx_bool within_rlist;

                        ci               = sci*NCL_PER_SUPERCL + im;

                        within_rlist     = FALSE;
                        npair            = 0;
                        for (ic = 0; ic < CL_SIZE; ic++)
                        {
                            ia               = ci*CL_SIZE + ic;

                            is               = ia*nbat->xstride;
                            ifs              = ia*nbat->fstride;
                            ix               = shX + x[is+0];
                            iy               = shY + x[is+1];
                            iz               = shZ + x[is+2];
                            iq               = facel*x[is+3];
                            nti              = ntype*2*type[ia];

                            fix              = 0;
                            fiy              = 0;
                            fiz              = 0;

                            for (jc = 0; jc < CL_SIZE; jc++)
                            {
                                ja               = cj*CL_SIZE + jc;

                                if (nbln->shift == CENTRAL &&
                                    ci == cj && ja <= ia)
                                {
                                    continue;
                                }

                                int_bit = ((excl[jc>>2]->pair[(jc & 3)*CL_SIZE+ic] >> (jm*NCL_PER_SUPERCL+im)) & 1);

                                js               = ja*nbat->xstride;
                                jfs              = ja*nbat->fstride;
                                jx               = x[js+0];
                                jy               = x[js+1];
                                jz               = x[js+2];
                                dx               = ix - jx;
                                dy               = iy - jy;
                                dz               = iz - jz;
                                rsq              = dx*dx + dy*dy + dz*dz;
                                if (rsq < rlist2)
                                {
                                    within_rlist = TRUE;
                                }
                                if (rsq >= rcut2)
                                {
                                    continue;
                                }

                                if (type[ia] != ntype-1 && type[ja] != ntype-1)
                                {
                                    npair++;
                                }

                                /* avoid NaN for excluded pairs at r=0 */
                                rsq             += (1.0 - int_bit)*NBNXN_AVOID_SING_R2_INC;

                                rinv             = gmx_invsqrt(rsq);
                                rinvsq           = rinv*rinv;
                                fscal            = 0;

                                qq               = iq*x[js+3];
                                if (!bEwald)
                                {
                                    /* Reaction-field */
                                    krsq  = iconst->k_rf*rsq;
                                    fscal = qq*(int_bit*rinv - 2*krsq)*rinvsq;
                                    if (bEner)
                                    {
                                        vcoul = qq*(int_bit*rinv + krsq - iconst->c_rf);
                                    }
                                }
                                else
                                {
                                    r     = rsq*rinv;
                                    rt    = r*iconst->tabq_scale;
                                    n0    = rt;
                                    eps   = rt - n0;

                                    fexcl = (1 - eps)*Ftab[n0] + eps*Ftab[n0+1];

                                    fscal = qq*(int_bit*rinvsq - fexcl)*rinv;

                                    if (bEner)
                                    {
                                        vcoul = qq*((int_bit - gmx_erf(iconst->ewaldcoeff_q*r))*rinv - int_bit*iconst->sh_ewald);
                                    }
                                }

                                if (rsq < rvdw2)
                                {
                                    tj        = nti + 2*type[ja];

                                    /* Vanilla Lennard-Jones cutoff */
                                    c6        = vdwparam[tj];
                                    c12       = vdwparam[tj+1];

                                    rinvsix   = int_bit*rinvsq*rinvsq*rinvsq;
                                    Vvdw_disp = c6*rinvsix;
                                    Vvdw_rep  = c12*rinvsix*rinvsix;
                                    fscal    += (Vvdw_rep - Vvdw_disp)*rinvsq;

                                    if (bEner)
                                    {
                                        vctot   += vcoul;

                                        Vvdwtot +=
                                            (Vvdw_rep - int_bit*c12*iconst->sh_invrc6*iconst->sh_invrc6)/12 -
                                            (Vvdw_disp - int_bit*c6*iconst->sh_invrc6)/6;
                                    }
                                }

                                tx        = fscal*dx;
                                ty        = fscal*dy;
                                tz        = fscal*dz;
                                fix       = fix + tx;
                                fiy       = fiy + ty;
                                fiz       = fiz + tz;
                                f[jfs+0] -= tx;
                                f[jfs+1] -= ty;
                                f[jfs+2] -= tz;
                            }

                            f[ifs+0]        += fix;
                            f[ifs+1]        += fiy;
                            f[ifs+2]        += fiz;
                            fshift[ish3]     = fshift[ish3]   + fix;
                            fshift[ish3+1]   = fshift[ish3+1] + fiy;
                            fshift[ish3+2]   = fshift[ish3+2] + fiz;

                            /* Count in half work-units.
                             * In CUDA one work-unit is 2 warps.
                             */
                            if ((ic+1) % (CL_SIZE/2) == 0)
                            {
                                npair_tot += npair;

                                nhwu++;
                                if (within_rlist)
                                {
                                    nhwu_pruned++;
                                }

                                within_rlist = FALSE;
                                npair        = 0;
                            }
                        }
                    }
                }
            }
        }

        if (bEner)
        {
            ggid             = 0;
            Vc[ggid]         = Vc[ggid]   + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
        }
    }

    if (debug)
    {
        fprintf(debug, "number of half %dx%d atom pairs: %d after pruning: %d fraction %4.2f\n",
                nbl->na_ci, nbl->na_ci,
                nhwu, nhwu_pruned, nhwu_pruned/(double)nhwu);
        fprintf(debug, "generic kernel pair interactions:            %d\n",
                nhwu*nbl->na_ci/2*nbl->na_ci);
        fprintf(debug, "generic kernel post-prune pair interactions: %d\n",
                nhwu_pruned*nbl->na_ci/2*nbl->na_ci);
        fprintf(debug, "generic kernel non-zero pair interactions:   %d\n",
                npair_tot);
        fprintf(debug, "ratio non-zero/post-prune pair interactions: %4.2f\n",
                npair_tot/(double)(nhwu_pruned*nbl->na_ci/2*nbl->na_ci));
    }
}
