/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/calcmu.h"

#include <stdio.h>
#include <stdlib.h>

#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"

void calc_mu(int start, int homenr, rvec x[], real q[], real qB[],
             int nChargePerturbed,
             dvec mu, dvec mu_B)
{
    int    i, end, m;
    double mu_x, mu_y, mu_z;

    end   = start + homenr;

    mu_x = mu_y = mu_z = 0.0;
#pragma omp parallel for reduction(+: mu_x, mu_y, mu_z) schedule(static) \
    num_threads(gmx_omp_nthreads_get(emntDefault))
    for (i = start; i < end; i++)
    {
        mu_x += q[i]*x[i][XX];
        mu_y += q[i]*x[i][YY];
        mu_z += q[i]*x[i][ZZ];
    }
    mu[XX] = mu_x;
    mu[YY] = mu_y;
    mu[ZZ] = mu_z;

    for (m = 0; (m < DIM); m++)
    {
        mu[m] *= ENM2DEBYE;
    }

    if (nChargePerturbed)
    {
        mu_x = mu_y = mu_z = 0.0;
#pragma omp parallel for reduction(+: mu_x, mu_y, mu_z) schedule(static) \
        num_threads(gmx_omp_nthreads_get(emntDefault))
        for (i = start; i < end; i++)
        {
            mu_x += qB[i]*x[i][XX];
            mu_y += qB[i]*x[i][YY];
            mu_z += qB[i]*x[i][ZZ];
        }
        mu_B[XX] = mu_x * ENM2DEBYE;
        mu_B[YY] = mu_y * ENM2DEBYE;
        mu_B[ZZ] = mu_z * ENM2DEBYE;
    }
    else
    {
        copy_dvec(mu, mu_B);
    }
}

gmx_bool read_mu(FILE *fp, rvec mu, real *vol)
{
    /* For backward compatibility */
    real mmm[4];

    if (fread(mmm, (size_t)(4*sizeof(real)), 1, fp) != 1)
    {
        return FALSE;
    }

    copy_rvec(mmm, mu);
    *vol = mmm[3];

    return TRUE;
}
