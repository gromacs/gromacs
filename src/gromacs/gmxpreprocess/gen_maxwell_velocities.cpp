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

#include "gen_maxwell_velocities.h"

#include <cmath>

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"

static void low_mspeed(real tempi, gmx_mtop_t* mtop, rvec v[], gmx::ThreeFry2x64<>* rng, const gmx::MDLogger& logger)
{
    int                                    nrdf;
    real                                   ekin, temp;
    gmx::TabulatedNormalDistribution<real> normalDist;

    ekin = 0.0;
    nrdf = 0;
    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        real          mass  = local.m;
        if (mass > 0)
        {
            rng->restart(i, 0);
            real sd = std::sqrt(gmx::c_boltz * tempi / mass);
            for (int m = 0; (m < DIM); m++)
            {
                v[i][m] = sd * normalDist(*rng);
                ekin += 0.5 * mass * v[i][m] * v[i][m];
            }
            nrdf += DIM;
        }
    }
    temp = (2.0 * ekin) / (nrdf * gmx::c_boltz);
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
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted("Velocities were taken from a Maxwell distribution at %g K", tempi);
    if (debug)
    {
        fprintf(debug,
                "Velocities were taken from a Maxwell distribution\n"
                "Initial generated temperature: %12.5e (scaled to: %12.5e)\n",
                temp,
                tempi);
    }
}

void maxwell_speed(real tempi, int seed, gmx_mtop_t* mtop, rvec v[], const gmx::MDLogger& logger)
{

    if (seed == -1)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Using random seed %d for generating velocities", seed);
    }
    gmx::ThreeFry2x64<> rng(seed, gmx::RandomDomain::MaxwellVelocities);

    low_mspeed(tempi, mtop, v, &rng, logger);
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

    return tm;
}

void stop_cm(const gmx::MDLogger gmx_unused& logger, int natoms, real mass[], rvec x[], rvec v[])
{
    rvec   xcm, vcm, acm;
    tensor L;
    int    i, m;

#ifdef DEBUG
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("stopping center of mass motion...");
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
}
