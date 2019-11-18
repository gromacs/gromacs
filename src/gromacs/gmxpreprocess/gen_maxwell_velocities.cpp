/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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

#include "gen_maxwell_velocities.h"

#include <cmath>

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static void low_mspeed(real tempi, gmx_mtop_t* mtop, rvec v[], gmx::ThreeFry2x64<>* rng)
{
    int                                    nrdf;
    real                                   boltz;
    real                                   ekin, temp;
    gmx::TabulatedNormalDistribution<real> normalDist;

    boltz = BOLTZ * tempi;
    ekin  = 0.0;
    nrdf  = 0;
    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        real          mass  = local.m;
        if (mass > 0)
        {
            rng->restart(i, 0);
            real sd = std::sqrt(boltz / mass);
            for (int m = 0; (m < DIM); m++)
            {
                v[i][m] = sd * normalDist(*rng);
                ekin += 0.5 * mass * v[i][m] * v[i][m];
            }
            nrdf += DIM;
        }
    }
    temp = (2.0 * ekin) / (nrdf * BOLTZ);
    if (temp > 0)
    {
        real scal = std::sqrt(tempi / temp);
        for (int i = 0; (i < mtop->natoms); i++)
        {
            for (int m = 0; (m < DIM); m++)
            {
                v[i][m] *= scal;
            }
        }
    }
    fprintf(stderr, "Velocities were taken from a Maxwell distribution at %g K\n", tempi);
    if (debug)
    {
        fprintf(debug,
                "Velocities were taken from a Maxwell distribution\n"
                "Initial generated temperature: %12.5e (scaled to: %12.5e)\n",
                temp, tempi);
    }
}

void maxwell_speed(real tempi, unsigned int seed, gmx_mtop_t* mtop, rvec v[])
{

    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
        fprintf(stderr, "Using random seed %u for generating velocities\n", seed);
    }
    gmx::ThreeFry2x64<> rng(seed, gmx::RandomDomain::MaxwellVelocities);

    low_mspeed(tempi, mtop, v, &rng);
}

static real calc_cm(int natoms, const real mass[], rvec x[], rvec v[], rvec xcm, rvec vcm, rvec acm, matrix L)
{
    rvec dx, a0;
    real tm, m0;
    int  i, m;

    clear_rvec(xcm);
    clear_rvec(vcm);
    clear_rvec(acm);
    tm = 0.0;
    for (i = 0; (i < natoms); i++)
    {
        m0 = mass[i];
        tm += m0;
        cprod(x[i], v[i], a0);
        for (m = 0; (m < DIM); m++)
        {
            xcm[m] += m0 * x[i][m]; /* c.o.m. position */
            vcm[m] += m0 * v[i][m]; /* c.o.m. velocity */
            acm[m] += m0 * a0[m];   /* rotational velocity around c.o.m. */
        }
    }
    cprod(xcm, vcm, a0);
    for (m = 0; (m < DIM); m++)
    {
        xcm[m] /= tm;
        vcm[m] /= tm;
        acm[m] -= a0[m] / tm;
    }

#define PVEC(str, v) \
    fprintf(log, "%s[X]: %10.5e  %s[Y]: %10.5e  %s[Z]: %10.5e\n", str, (v)[0], str, (v)[1], str, (v)[2])
#ifdef DEBUG
    PVEC("xcm", xcm);
    PVEC("acm", acm);
    PVEC("vcm", vcm);
#endif

    clear_mat(L);
    for (i = 0; (i < natoms); i++)
    {
        m0 = mass[i];
        for (m = 0; (m < DIM); m++)
        {
            dx[m] = x[i][m] - xcm[m];
        }
        L[XX][XX] += dx[XX] * dx[XX] * m0;
        L[XX][YY] += dx[XX] * dx[YY] * m0;
        L[XX][ZZ] += dx[XX] * dx[ZZ] * m0;
        L[YY][YY] += dx[YY] * dx[YY] * m0;
        L[YY][ZZ] += dx[YY] * dx[ZZ] * m0;
        L[ZZ][ZZ] += dx[ZZ] * dx[ZZ] * m0;
    }
#ifdef DEBUG
    PVEC("L-x", L[XX]);
    PVEC("L-y", L[YY]);
    PVEC("L-z", L[ZZ]);
#endif

    return tm;
}

void stop_cm(FILE gmx_unused* log, int natoms, real mass[], rvec x[], rvec v[])
{
    rvec   xcm, vcm, acm;
    tensor L;
    int    i, m;

#ifdef DEBUG
    fprintf(log, "stopping center of mass motion...\n");
#endif
    (void)calc_cm(natoms, mass, x, v, xcm, vcm, acm, L);

    /* Subtract center of mass velocity */
    for (i = 0; (i < natoms); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            v[i][m] -= vcm[m];
        }
    }

#ifdef DEBUG
    (void)calc_cm(log, natoms, mass, x, v, xcm, vcm, acm, L);
#endif
}
