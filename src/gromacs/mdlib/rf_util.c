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
#include "gmxpre.h"

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"

real RF_excl_correction(const t_forcerec *fr, t_graph *g,
                        const t_mdatoms *mdatoms, const t_blocka *excl,
                        rvec x[], rvec f[], rvec *fshift, const t_pbc *pbc,
                        real lambda, real *dvdlambda)
{
    /* Calculate the reaction-field energy correction for this node:
     * epsfac q_i q_j (k_rf r_ij^2 - c_rf)
     * and force correction for all excluded pairs, including self pairs.
     */
    int         top, i, j, j1, j2, k, ki;
    double      q2sumA, q2sumB, ener;
    const real *chargeA, *chargeB;
    real        ek, ec, L1, qiA, qiB, qqA, qqB, qqL, v;
    rvec        dx, df;
    atom_id    *AA;
    ivec        dt;
    int         start = 0;
    int         end   = mdatoms->homenr;
    int         niat;
    gmx_bool    bMolPBC = fr->bMolPBC;

    if (fr->n_tpi)
    {
        /* For test particle insertion we only correct for the test molecule */
        start = mdatoms->nr - fr->n_tpi;
    }

    ek      = fr->epsfac*fr->k_rf;
    ec      = fr->epsfac*fr->c_rf;
    chargeA = mdatoms->chargeA;
    chargeB = mdatoms->chargeB;
    AA      = excl->a;
    ki      = CENTRAL;

    if (fr->bDomDec)
    {
        niat = excl->nr;
    }
    else
    {
        niat = end;
    }

    q2sumA = 0;
    q2sumB = 0;
    ener   = 0;
    if (mdatoms->nChargePerturbed == 0)
    {
        for (i = start; i < niat; i++)
        {
            qiA = chargeA[i];
            if (i < end)
            {
                q2sumA += qiA*qiA;
            }
            /* Do the exclusions */
            j1  = excl->index[i];
            j2  = excl->index[i+1];
            for (j = j1; j < j2; j++)
            {
                k = AA[j];
                if (k > i)
                {
                    qqA = qiA*chargeA[k];
                    if (qqA != 0)
                    {
                        if (g)
                        {
                            rvec_sub(x[i], x[k], dx);
                            ivec_sub(SHIFT_IVEC(g, i), SHIFT_IVEC(g, k), dt);
                            ki = IVEC2IS(dt);
                        }
                        else if (bMolPBC)
                        {
                            ki = pbc_dx_aiuc(pbc, x[i], x[k], dx);
                        }
                        else
                        {
                            rvec_sub(x[i], x[k], dx);
                        }
                        ener += qqA*(ek*norm2(dx) - ec);
                        svmul(-2*qqA*ek, dx, df);
                        rvec_inc(f[i], df);
                        rvec_dec(f[k], df);
                        rvec_inc(fshift[ki], df);
                        rvec_dec(fshift[CENTRAL], df);
                    }
                }
            }
        }
        ener += -0.5*ec*q2sumA;
    }
    else
    {
        L1 = 1.0 - lambda;
        for (i = start; i < niat; i++)
        {
            qiA = chargeA[i];
            qiB = chargeB[i];
            if (i < end)
            {
                q2sumA += qiA*qiA;
                q2sumB += qiB*qiB;
            }
            /* Do the exclusions */
            j1  = excl->index[i];
            j2  = excl->index[i+1];
            for (j = j1; j < j2; j++)
            {
                k = AA[j];
                if (k > i)
                {
                    qqA = qiA*chargeA[k];
                    qqB = qiB*chargeB[k];
                    if (qqA != 0 || qqB != 0)
                    {
                        qqL = L1*qqA + lambda*qqB;
                        if (g)
                        {
                            rvec_sub(x[i], x[k], dx);
                            ivec_sub(SHIFT_IVEC(g, i), SHIFT_IVEC(g, k), dt);
                            ki = IVEC2IS(dt);
                        }
                        else if (bMolPBC)
                        {
                            ki = pbc_dx_aiuc(pbc, x[i], x[k], dx);
                        }
                        else
                        {
                            rvec_sub(x[i], x[k], dx);
                        }
                        v     = ek*norm2(dx) - ec;
                        ener += qqL*v;
                        svmul(-2*qqL*ek, dx, df);
                        rvec_inc(f[i], df);
                        rvec_dec(f[k], df);
                        rvec_inc(fshift[ki], df);
                        rvec_dec(fshift[CENTRAL], df);
                        *dvdlambda += (qqB - qqA)*v;
                    }
                }
            }
        }
        ener       += -0.5*ec*(L1*q2sumA + lambda*q2sumB);
        *dvdlambda += -0.5*ec*(q2sumB - q2sumA);
    }

    if (debug)
    {
        fprintf(debug, "RF exclusion energy: %g\n", ener);
    }

    return ener;
}

void calc_rffac(FILE *fplog, int eel, real eps_r, real eps_rf, real Rc, real Temp,
                real zsq, matrix box,
                real *kappa, real *krf, real *crf)
{
    /* Compute constants for Generalized reaction field */
    real   k1, k2, I, vol, rmin;

    if (EEL_RF(eel))
    {
        vol     = det(box);
        if (eel == eelGRF)
        {
            /* Consistency check */
            if (Temp <= 0.0)
            {
                gmx_fatal(FARGS, "Temperature is %f while using"
                          " Generalized Reaction Field\n", Temp);
            }
            /* Ionic strength (only needed for eelGRF */
            I       = 0.5*zsq/vol;
            *kappa  = sqrt(2*I/(EPSILON0*eps_rf*BOLTZ*Temp));
        }
        else
        {
            I      = 0;
            *kappa = 0;
        }

        /* eps == 0 signals infinite dielectric */
        if (eps_rf == 0)
        {
            *krf = 1/(2*Rc*Rc*Rc);
        }
        else
        {
            k1   = 1 + *kappa*Rc;
            k2   = eps_rf*sqr((real)(*kappa*Rc));

            *krf = ((eps_rf - eps_r)*k1 + 0.5*k2)/((2*eps_rf + eps_r)*k1 + k2)/(Rc*Rc*Rc);
        }
        *crf   = 1/Rc + *krf*Rc*Rc;
        rmin   = pow(*krf*2.0, -1.0/3.0);

        if (fplog)
        {
            if (eel == eelGRF)
            {
                please_cite(fplog, "Tironi95a");
                fprintf(fplog, "%s:\n"
                        "epsRF = %10g, I   = %10g, volume = %10g, kappa  = %10g\n"
                        "rc    = %10g, krf = %10g, crf    = %10g, epsfac = %10g\n",
                        eel_names[eel], eps_rf, I, vol, *kappa, Rc, *krf, *crf,
                        ONE_4PI_EPS0/eps_r);
            }
            else
            {
                fprintf(fplog, "%s:\n"
                        "epsRF = %g, rc = %g, krf = %g, crf = %g, epsfac = %g\n",
                        eel_names[eel], eps_rf, Rc, *krf, *crf, ONE_4PI_EPS0/eps_r);
            }
            fprintf(fplog,
                    "The electrostatics potential has its minimum at r = %g\n",
                    rmin);
        }
    }
}

void init_generalized_rf(FILE *fplog,
                         const gmx_mtop_t *mtop, const t_inputrec *ir,
                         t_forcerec *fr)
{
    int                  mb, i, j;
    real                 q, zsq, nrdf, T;
    const gmx_moltype_t *molt;
    const t_block       *cgs;

    if (ir->efep != efepNO && fplog)
    {
        fprintf(fplog, "\nWARNING: the generalized reaction field constants are determined from topology A only\n\n");
    }
    zsq = 0.0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        cgs  = &molt->cgs;
        for (i = 0; (i < cgs->nr); i++)
        {
            q = 0;
            for (j = cgs->index[i]; (j < cgs->index[i+1]); j++)
            {
                q += molt->atoms.atom[j].q;
            }
            zsq += mtop->molblock[mb].nmol*q*q;
        }
        fr->zsquare = zsq;
    }

    T    = 0.0;
    nrdf = 0.0;
    for (i = 0; (i < ir->opts.ngtc); i++)
    {
        nrdf += ir->opts.nrdf[i];
        T    += (ir->opts.nrdf[i] * ir->opts.ref_t[i]);
    }
    if (nrdf == 0)
    {
        gmx_fatal(FARGS, "No degrees of freedom!");
    }
    fr->temp   = T/nrdf;
}
