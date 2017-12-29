/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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

#include "entropy.h"

#include <cmath>
#include <cstdio>

#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"

real S_Translate(real mass,
                 real temperature,
                 real pressure)
{
    // Following the equations in J. W. Ochterski, 
    // Thermochemistry in Gaussian, Gaussian, Inc., 2000
    // Pitssburg PA
    double kT = BOLTZ*temperature;
    
    GMX_RELEASE_ASSERT(mass > 0, "Molecular mass should be larger than zero");
    GMX_RELEASE_ASSERT(pressure > 0, "Pressure should be larger than zero");
    GMX_RELEASE_ASSERT(temperature > 0, "Temperature should be larger than zero");
    // Convert bar to Pascal
    double P  = pressure*1e5;
    double qT = (std::pow(2*M_PI*mass*kT/gmx::square(PLANCK), 1.5) * 
                 (kT/P) * (1e30/AVOGADRO));
    return RGAS*(std::log(qT) + 2.5);
}

real S_Rotate(real temperature,
              int natom,
              bool linear,
              rvec theta,
              real sigma_r)
{
    GMX_RELEASE_ASSERT(sigma_r > 0, "Symmetry factor should be larger than zero");
    GMX_RELEASE_ASSERT(temperature > 0, "Temperature should be larger than zero");

    double sR = 0;
    if (natom > 1)
    {
        if (linear)
        {
            GMX_RELEASE_ASSERT(theta[0] > 0, "Theta should be larger than zero");
            double qR = temperature/(sigma_r * theta[0]);
            sR        = RGAS * (std::log(qR) + 1);
        }
        else
        {
            double Q  = theta[XX]*theta[YY]*theta[ZZ];
            GMX_RELEASE_ASSERT(Q > 0, "Q should be larger than zero");
            double qR = std::sqrt(M_PI * std::pow(temperature, 3)/Q)/sigma_r;
            sR        = RGAS * (std::log(qR) + 1.5);
        }
    }
    return sR;
}

static real eigval_to_frequency(real eigval)
{
    double factor_gmx_to_omega2       = 1.0E21/(AVOGADRO*AMU);
    return std::sqrt(eigval*factor_gmx_to_omega2);
}

real calc_entropy_quasi_harmonic(int      n,
                                 real     eigval[],
                                 real     temperature,
                                 gmx_bool bLinear)
{
    int    nskip = bLinear ? 5 : 6;
    double S     = 0;
    double hbar  = PLANCK1/(2*M_PI);
    for (int i = nskip; (i < n); i++)
    {
        if (eigval[i] > 0)
        {
            double omega  = eigval_to_frequency(eigval[i]);
            double hwkT   = (hbar*omega)/(BOLTZMANN*temperature);
            double dS     = (hwkT/std::expm1(hwkT) - std::log1p(-std::exp(-hwkT)));
            S            += dS;
            if (debug)
            {
                fprintf(debug, "i = %5d eigval = %10g w = %10g hwkT = %10g dS = %10g\n",
                        i, eigval[i], omega, hwkT, dS);
            }
        }
        else
        {
            fprintf(stderr, "eigval[%d] = %g\n", i, eigval[i]);
        }
    }
    return S*RGAS;
}

real calc_entropy_schlitter(int      n,
                            real     eigval[],
                            real     temperature,
                            gmx_bool bLinear)
{
    int    nskip  = bLinear ? 5 : 6;
    double hbar   = PLANCK1/(2*M_PI);             // J s
    double kt     = BOLTZMANN*temperature;        // J
    double kteh   = kt*std::exp(2.0)/(hbar*hbar); // 1/(J s^2) = 1/(kg m^2)
    double evcorr = NANO*NANO*AMU;
    if (debug)
    {
        fprintf(debug, "n = %d, kteh = %g evcorr = %g\n", n, kteh, evcorr);
    }
    printf("n = %d, kteh = %g evcorr = %g\n", n, kteh, evcorr);
    double deter = 0;
    for (int i = nskip; (i < n); i++)
    {
        double dd    = 1+kteh*eigval[i]*evcorr;
        deter       += std::log(dd);
    }
    return 0.5*RGAS*deter;
}
