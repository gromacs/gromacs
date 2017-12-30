/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "entropy.h"

#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"

static real calc_entropy_qh(int  n,
                            real eigval[],
                            real temperature,
                            int  nskip_eigval)
{
    double S    = 0;
    double hbar = PLANCK1/(2*M_PI);
    for (int i = 0; (i < n-nskip_eigval); i++)
    {
        if (eigval[i] > 0)
        {
            double lambda = eigval[i]*AMU;
            double w      = std::sqrt(BOLTZMANN*temperature/lambda)/NANO;
            double hwkT   = (hbar*w)/(BOLTZMANN*temperature);
            double dS     = (hwkT/std::expm1(hwkT) - std::log1p(-std::exp(-hwkT)));
            S            += dS;
            if (debug)
            {
                fprintf(debug, "i = %5d w = %10g lam = %10g hwkT = %10g dS = %10g\n",
                        i, w, lambda, hwkT, dS);
            }
        }
        else
        {
            fprintf(stderr, "eigval[%d] = %g\n", i, eigval[i]);
        }
    }
    return S*RGAS;
}

static real calc_entropy_schlitter(int  n,
                                   real eigval[],
                                   real temperature,
                                   int  nskip_eigval)
{
    double hbar = PLANCK1/(2*M_PI);
    double kt   = BOLTZMANN*temperature;
    double kteh = kt*std::exp(2.0)/(hbar*hbar)*AMU*(NANO*NANO);
    if (debug)
    {
        fprintf(debug, "n = %d, nskip_eigval = %d kteh = %g\n", n, nskip_eigval, kteh);
    }

    double deter = 0;
    for (int i = 0; (i < n-nskip_eigval); i++)
    {
        double dd = 1+kteh*eigval[i];
        deter    += std::log(dd);
    }
    return 0.5*RGAS*deter;
}

real calc_entropy(EntropyAlgorithm algorithm,
                  int              n,
                  real             eigval[],
                  real             temperature,
                  int              nskip_eigval)
{
    switch (algorithm)
    {
        case eaQuasiHarmonic:
            return calc_entropy_qh(n, eigval, temperature, nskip_eigval);
            break;
        case eaSchlitter:
            return calc_entropy_schlitter(n, eigval, temperature, nskip_eigval);
            break;
    }
}
