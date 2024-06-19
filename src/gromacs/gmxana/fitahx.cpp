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
#include "gmxpre.h"

#include "fitahx.h"

#include <cmath>

#include <filesystem>
#include <vector>

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"

static void my_calc_xcm(int nbb, const int bbind[], rvec x[], rvec xcm)
{
    int i, m, ai;

    clear_rvec(xcm);
    for (i = 0; (i < nbb); i++)
    {
        ai = bbind[i];
        rvec_inc(xcm, x[ai]);
    }
    for (m = 0; (m < DIM); m++)
    {
        xcm[m] /= (nbb);
    }
}

static void my_sub_xcm(int nbb, const int bbind[], rvec x[], rvec xcm)
{
    int i, ai;

    for (i = 0; (i < nbb); i++)
    {
        ai = bbind[i];
        rvec_dec(x[ai], xcm);
    }
}

real fit_ahx(int nres, t_bb bb[], int natoms, int nall, int allindex[], rvec x[], int nca, int caindex[], gmx_bool bFit)
{
    static std::vector<gmx::RVec> xref;
    static std::vector<real>      mass;
    const real                    d   = 0.15;  /* Rise per residue (nm)    */
    const real                    tw  = 1.745; /* Twist per residue (rad)  */
    const real                    rad = 0.23;  /* Radius of the helix (nm) */
    real                          phi0, trms, rms;
    rvec                          dx, xcm;
    int                           ai, i, nmass;

    if (nca < 3)
    {
        gmx_fatal(FARGS, "Need at least 3 Calphas to fit to, (now %d)...\n", nca);
    }

    if (xref.empty())
    {
        xref.resize(natoms);
        mass.resize(natoms);
    }
    phi0 = 0;
    for (i = 0; (i < nca); i++)
    {
        ai           = caindex[i];
        xref[ai][XX] = rad * std::cos(phi0);
        xref[ai][YY] = rad * std::sin(phi0);
        xref[ai][ZZ] = d * i;

        /* Set the mass to select that this Calpha contributes to fitting */
        mass[ai] = 10.0;
/*#define DEBUG*/
#ifdef DEBUG
        fprintf(stderr,
                "%5d  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
                ai,
                x[ai][XX],
                x[ai][YY],
                x[ai][ZZ],
                xref[ai][XX],
                xref[ai][YY],
                xref[ai][ZZ]);
#endif
        phi0 += tw;
    }

    /* Center the referece around the origin */
    my_calc_xcm(nca, caindex, as_rvec_array(xref.data()), xcm);
    my_sub_xcm(nca, caindex, as_rvec_array(xref.data()), xcm);

    if (bFit)
    {
        /* Now center the to-be-fitted coords around the origin */
        my_calc_xcm(nca, caindex, x, xcm);
        my_sub_xcm(nall, allindex, x, xcm);
    }
#ifdef DEBUG
    dump_ahx(nres, bb, xref, box, 0);
#endif

    nmass = 0;
    for (i = 0; (i < natoms); i++)
    {
        if (mass[i] > 0)
        {
            nmass++;
        }
    }
    if (nmass != nca)
    {
        gmx_fatal(FARGS, "nmass=%d, nca=%d\n", nmass, nca);
    }

    /* Now call the fitting routine */
    if (bFit)
    {
        do_fit(natoms, mass.data(), as_rvec_array(xref.data()), x);
    }

    /* Reset masses and calc rms */
    trms = 0.0;
    for (i = 0; (i < nres); i++)
    {
        ai = bb[i].CA;

        if (mass[ai] > 0.0)
        {
            rvec_sub(x[ai], xref[ai], dx);
            rms = iprod(dx, dx);
            bb[i].rmsa += std::sqrt(rms);
            bb[i].nrms++;
            trms += rms;
            mass[ai] = 0.0;
        }
    }
    return std::sqrt(trms / nca);
}
