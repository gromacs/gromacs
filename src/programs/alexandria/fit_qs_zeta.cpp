/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <ctype.h>
#include <stdlib.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/random.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "poldata.h"
#include "poldata_low.h"
#include "poldata_xml.h"

static double chi2_coulomb(double zeta, double qs, double alpha, 
                           double delta_q, double rmax,
                           int row, 
                           alexandria::ChargeDistributionModel iChargeDistributionModel)
{
    double chi2 = 0;
    int    imax = 10;
    for(int i = 1; i<=imax; i++)
    {
        double r        = i*(rmax/imax);
        double fpol     = qs*r/alpha;
        double fcoulomb = 0;
        
        switch (iChargeDistributionModel)
        {
        case alexandria::eqdAXpg:
            fcoulomb = -DNuclear_GG(r, zeta);
            break;
        case alexandria::eqdAXps:
            fcoulomb = -DNuclear_SS(r, row, zeta);
            break;
        default:
            fcoulomb = 0;
        }
        double fg = (delta_q-qs)*fcoulomb;
        chi2 += gmx::square(fpol - fg);
    }
    return chi2;
}

static void fit_polarization(double alpha, double delta_q, double rmax, int row,
                             alexandria::ChargeDistributionModel iChargeDistributionModel,
                             int maxiter, double tolerance,
                             double *zeta_opt, double *qs_opt)
{
    double zeta = 0;
    double qs   = -4;
            
    if (alpha > 0)
    {
        // Guess for delta_q == 0 and ChargeDistributionModel == eqdAXpg
        zeta = std::pow(3.0*std::sqrt(M_PI)/(4.0*alpha), 1.0/3.0);
        printf("zeta0 %g qs %g delta_q %g\n", zeta, qs, delta_q);
        if (delta_q != 0.0)
        {
            double qs_min    = 2*qs;
            double qs_max    = qs/2;
            double zeta_max  = 2*zeta;
            double zeta_min  = zeta/2;
            double chi2      = chi2_coulomb(zeta, qs, alpha, delta_q, rmax, row, 
                                            iChargeDistributionModel);
            double chi2_opt  = chi2;
            double kT        = 1;// 100*tolerance;
            bool   bZeta     = true;
            int    iter      = 0;
            gmx::ThreeFry2x64<64>         rng(123457, gmx::RandomDomain::Other);
            gmx::NormalDistribution<real> dist(1.0, 0.02);
            gmx::NormalDistribution<real> mc(0.0, 1.0);
            printf("Iter: %5d  chi2: %8g  qs: %7.3f  zeta: %7.3f\n", iter, chi2, qs, zeta);
            while (chi2 > tolerance && iter < maxiter)
            {
                double zeta_new = zeta, qs_new = qs;
                real factor = dist(rng);
                if (bZeta)
                {
                    zeta_new = std::max(zeta_min, std::min(zeta_max, zeta*factor));
                }
                else
                {
                    qs_new = std::min(qs_max, std::max(qs_min, qs*factor));
                }
                double chi2_new = chi2_coulomb(zeta_new, qs_new, alpha, delta_q, rmax, row, 
                                               iChargeDistributionModel);
                                           
                if (chi2_new < chi2 || std::exp(-(chi2_new - chi2)/chi2) < mc(rng))
                {
                    if (chi2_new < chi2_opt)
                    {
                        *qs_opt   = qs_new;
                        *zeta_opt = zeta_new;
                        chi2_opt  = chi2_new;
                    }
                    qs   = qs_new;
                    zeta = zeta_new;
                    chi2 = chi2_new;
                    printf("Iter: %5d  chi2: %8g  qs: %7.3f  zeta: %7.3f\n", iter, chi2_opt, *qs_opt, *zeta_opt);
                }
                bZeta = !bZeta;
                iter += 1;
            }
            if (chi2 > tolerance)
            {
                zeta = 0;
            }
        }
    }
}

int alex_fit_qs_zeta(int argc, char *argv[])
{
    static const char               *desc[] = {
        "Determine distributed charge parameters based on polarizabilities."
    };
    gmx_output_env_t                *oenv;
    gmx_atomprop_t                   aps = nullptr;

    t_filenm                         fnm[] = {
        { efDAT, "-d",  "gentop",     ffOPTRD },
        { efDAT, "-do", "gentop_out", ffOPTWR },
        { efTEX, "-t",  "beta",       ffOPTWR }
    };

    const  int         NFILE          = asize(fnm);
    static const char *cqdist[]       = { nullptr, "AXpg", "AXps", nullptr };
    gmx_bool           bVerbose       = TRUE;
    real               alpha          = 0;
    real               delta_q        = 0;
    real               rmax           = 0.005;
    real               tolerance      = 0.1;
    int                maxiter        = 10000;
    int                row            = 1;
    alexandria::ChargeDistributionModel iChargeDistributionModel;

    t_pargs                  pa[]     = {
        { "-v",       FALSE, etBOOL, {&bVerbose},
          "Generate verbose output on the terminal." },
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Charge distribution used" },
        { "-alpha",   FALSE, etREAL, {&alpha},
          "Polarizability (nm^3)" },
        { "-delta_q", FALSE, etREAL, {&delta_q},
          "Charge on the particle" },
        { "-rmax",    FALSE, etREAL, {&rmax},
          "Maximum core-shell distance to use in fitting" },
        { "-toler",    FALSE, etREAL, {&tolerance},
          "Tolerance for optimization" },
        { "-maxiter",     FALSE, etINT,  {&maxiter},
          "Maximum number of iterations" },
        { "-row",     FALSE, etINT,  {&row},
          "Row in the periodic table" }
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    if ((iChargeDistributionModel = alexandria::name2eemtype(cqdist[0])) == 
        alexandria::eqdNR)
    {
        gmx_fatal(FARGS, "Invalid Charge Distribution model %s.\n",
                  cqdist[0]);
    }

    alexandria::Poldata       pd;
    const char *gentop_fnm = opt2fn_null("-d", NFILE, fnm);
    if (gentop_fnm)
    {
        try
        {
            const char *gentop_fnm = opt2fn_null("-d", NFILE, fnm);
            alexandria::readPoldata(gentop_fnm ? gentop_fnm : "", pd, aps);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    else
    {
        double zeta, qs;
        fit_polarization(alpha, delta_q, rmax, row,
                         iChargeDistributionModel, maxiter, tolerance,
                         &zeta, &qs);
        if (zeta > 0)
        {
            printf("qdist = %s alpha = %g delta_q = %g zeta = %g qs = %g\n",
                   cqdist[0], alpha, delta_q, zeta, qs);
        }
        else
        {
            printf("Do not know how to deal with qdist = %s, alpha = %g delta_q = %g\n",
                   cqdist[0], alpha, delta_q);
        }
    }


    return 0;
}
