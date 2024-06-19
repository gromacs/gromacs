/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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

#include "thermochemistry.h"

#include <cmath>
#include <cstdio>

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

static double eigval_to_frequency(double eigval)
{
    double factor_gmx_to_omega2 = 1.0E21 / (gmx::c_avogadro * gmx::c_amu);
    return std::sqrt(std::max(0.0, eigval) * factor_gmx_to_omega2);
}

double calcZeroPointEnergy(gmx::ArrayRef<const real> eigval, real scale_factor)
{
    // Convert frequency (ps^-1) to energy (kJ/mol)
    double factor = gmx::c_planck * gmx::c_pico / (2.0 * M_PI);
    double zpe    = 0;
    for (const auto& r : eigval)
    {
        double omega = eigval_to_frequency(r);
        zpe += 0.5 * factor * scale_factor * omega;
    }
    return zpe;
}

double calcVibrationalInternalEnergy(gmx::ArrayRef<const real> eigval, real temperature, gmx_bool linear, real scale_factor)
{
    size_t nskip = linear ? 5 : 6;
    double Evib  = 0;
    double hbar  = gmx::c_planck1 / (2 * M_PI);
    for (gmx::Index i = nskip; i < eigval.ssize(); i++)
    {
        if (eigval[i] > 0)
        {
            double omega = scale_factor * eigval_to_frequency(eigval[i]);
            double hwkT  = (hbar * omega) / (gmx::c_boltzmann * temperature);
            // Prevent overflow by checking for unreasonably large numbers.
            if (hwkT < 100)
            {
                double dEvib = hwkT * (0.5 + 1.0 / (std::expm1(hwkT)));
                if (debug)
                {
                    fprintf(debug,
                            "i %d eigval %g omega %g hwkT %g dEvib %g\n",
                            static_cast<int>(i + 1),
                            static_cast<double>(eigval[i]),
                            omega,
                            hwkT,
                            dEvib);
                }
                Evib += dEvib;
            }
        }
    }
    return temperature * gmx::c_boltz * Evib;
}

double calcVibrationalHeatCapacity(gmx::ArrayRef<const real> eigval, real temperature, gmx_bool linear, real scale_factor)
{
    size_t nskip = linear ? 5 : 6;
    double cv    = 0;
    double hbar  = gmx::c_planck1 / (2 * M_PI);
    for (gmx::Index i = nskip; i < eigval.ssize(); i++)
    {
        if (eigval[i] > 0)
        {
            double omega = scale_factor * eigval_to_frequency(eigval[i]);
            double hwkT  = (hbar * omega) / (gmx::c_boltzmann * temperature);
            // Prevent overflow by checking for unreasonably large numbers.
            if (hwkT < 100)
            {
                double dcv = std::exp(hwkT) * gmx::square(hwkT / std::expm1(hwkT));
                if (debug)
                {
                    fprintf(debug,
                            "i %d eigval %g omega %g hwkT %g dcv %g\n",
                            static_cast<int>(i + 1),
                            static_cast<double>(eigval[i]),
                            omega,
                            hwkT,
                            dcv);
                }
                cv += dcv;
            }
        }
    }
    return gmx::c_universalGasConstant * cv;
}

double calcTranslationalEntropy(real mass, real temperature, real pressure)
{
    double kT = gmx::c_boltz * temperature;

    GMX_RELEASE_ASSERT(mass > 0, "Molecular mass should be larger than zero");
    GMX_RELEASE_ASSERT(pressure > 0, "Pressure should be larger than zero");
    GMX_RELEASE_ASSERT(temperature > 0, "Temperature should be larger than zero");
    // Convert bar to Pascal
    double P  = pressure * 1e5;
    double qT = (std::pow(2 * M_PI * mass * kT / gmx::square(gmx::c_planck), 1.5) * (kT / P)
                 * (1e30 / gmx::c_avogadro));
    return gmx::c_universalGasConstant * (std::log(qT) + 2.5);
}

double calcRotationalEntropy(real temperature, int natom, gmx_bool linear, const rvec theta, real sigma_r)
{
    GMX_RELEASE_ASSERT(sigma_r > 0, "Symmetry factor should be larger than zero");
    GMX_RELEASE_ASSERT(temperature > 0, "Temperature should be larger than zero");

    double sR = 0;
    if (natom > 1)
    {
        if (linear)
        {
            GMX_RELEASE_ASSERT(theta[0] > 0, "Theta should be larger than zero");
            double qR = temperature / (sigma_r * theta[0]);
            sR        = gmx::c_universalGasConstant * (std::log(qR) + 1);
        }
        else
        {
            double Q = theta[XX] * theta[YY] * theta[ZZ];
            GMX_RELEASE_ASSERT(Q > 0, "Q should be larger than zero");
            double qR = std::sqrt(M_PI * std::pow(temperature, 3) / Q) / sigma_r;
            sR        = gmx::c_universalGasConstant * (std::log(qR) + 1.5);
        }
    }
    return sR;
}

double calcQuasiHarmonicEntropy(gmx::ArrayRef<const real> eigval, real temperature, gmx_bool bLinear, real scale_factor)
{
    size_t nskip = bLinear ? 5 : 6;
    double S     = 0;
    double hbar  = gmx::c_planck1 / (2 * M_PI);
    for (gmx::Index i = nskip; (i < eigval.ssize()); i++)
    {
        if (eigval[i] > 0)
        {
            double omega = scale_factor * eigval_to_frequency(eigval[i]);
            double hwkT  = (hbar * omega) / (gmx::c_boltzmann * temperature);
            double dS    = (hwkT / std::expm1(hwkT) - std::log1p(-std::exp(-hwkT)));
            S += dS;
            if (debug)
            {
                fprintf(debug,
                        "i = %5d eigval = %10g w = %10g hwkT = %10g dS = %10g\n",
                        static_cast<int>(i + 1),
                        static_cast<double>(eigval[i]),
                        omega,
                        hwkT,
                        dS);
            }
        }
        else if (debug)
        {
            fprintf(debug, "eigval[%d] = %g\n", static_cast<int>(i + 1), static_cast<double>(eigval[i]));
        }
    }
    return S * gmx::c_universalGasConstant;
}

double calcSchlitterEntropy(gmx::ArrayRef<const real> eigval, real temperature, gmx_bool bLinear)
{
    size_t nskip  = bLinear ? 5 : 6;
    double hbar   = gmx::c_planck1 / (2 * M_PI);        // J s
    double kt     = gmx::c_boltzmann * temperature;     // J
    double kteh   = kt * std::exp(2.0) / (hbar * hbar); // 1/(J s^2) = 1/(kg m^2)
    double evcorr = gmx::c_nano * gmx::c_nano * gmx::c_amu;
    if (debug)
    {
        fprintf(debug, "n = %td, kteh = %g evcorr = %g\n", ssize(eigval), kteh, evcorr);
    }
    double deter = 0;
    for (gmx::Index i = nskip; i < eigval.ssize(); i++)
    {
        double dd = 1 + kteh * eigval[i] * evcorr;
        deter += std::log(dd);
    }
    return 0.5 * gmx::c_universalGasConstant * deter;
}
