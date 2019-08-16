/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "rf_util.h"

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

void calc_rffac(FILE *fplog, int eel, real eps_r, real eps_rf, real Rc, real Temp,
                real zsq, matrix box,
                real *krf, real *crf)
{
    /* Compute constants for Generalized reaction field */
    real   kappa, k1, k2, I, rmin;
    real   vol = 0;

    if (EEL_RF(eel))
    {
        if (eel == eelGRF)
        {
            /* Consistency check */
            if (Temp <= 0.0)
            {
                gmx_fatal(FARGS, "Temperature is %f while using"
                          " Generalized Reaction Field\n", Temp);
            }
            /* Ionic strength (only needed for eelGRF */
            vol     = det(box);
            I       = 0.5*zsq/vol;
            kappa   = std::sqrt(2*I/(EPSILON0*eps_rf*BOLTZ*Temp));
        }
        else
        {
            I      = 0;
            kappa  = 0;
        }

        /* eps == 0 signals infinite dielectric */
        if (eps_rf == 0)
        {
            *krf = 1/(2*Rc*Rc*Rc);
        }
        else
        {
            k1   = 1 + kappa*Rc;
            k2   = eps_rf*gmx::square(static_cast<real>(kappa*Rc));

            *krf = ((eps_rf - eps_r)*k1 + 0.5*k2)/((2*eps_rf + eps_r)*k1 + k2)/(Rc*Rc*Rc);
        }
        *crf   = 1/Rc + *krf*Rc*Rc;
        // Make sure we don't lose resolution in pow() by casting real arg to double
        rmin   = gmx::invcbrt(static_cast<double>(*krf*2.0));

        if (fplog)
        {
            if (eel == eelGRF)
            {
                please_cite(fplog, "Tironi95a");
                fprintf(fplog, "%s:\n"
                        "epsRF = %10g, I   = %10g, volume = %10g, kappa  = %10g\n"
                        "rc    = %10g, krf = %10g, crf    = %10g, epsfac = %10g\n",
                        eel_names[eel], eps_rf, I, vol, kappa, Rc, *krf, *crf,
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
    int                  i, j;
    real                 q, zsq, nrdf, T;
    const gmx_moltype_t *molt;
    const t_block       *cgs;

    if (ir->efep != efepNO && fplog)
    {
        fprintf(fplog, "\nWARNING: the generalized reaction field constants are determined from topology A only\n\n");
    }
    zsq = 0.0;
    for (const gmx_molblock_t &molb : mtop->molblock)
    {
        molt = &mtop->moltype[molb.type];
        cgs  = &molt->cgs;
        for (i = 0; (i < cgs->nr); i++)
        {
            q = 0;
            for (j = cgs->index[i]; (j < cgs->index[i+1]); j++)
            {
                q += molt->atoms.atom[j].q;
            }
            zsq += molb.nmol*q*q;
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
