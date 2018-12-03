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
        double fpol     = -qs*r/alpha;
        double fcoulomb = 0;
        
        switch (iChargeDistributionModel)
        {
        case alexandria::eqdAXpg:
            fcoulomb = DNuclear_GG(r, zeta);
            break;
        case alexandria::eqdAXps:
            fcoulomb = DNuclear_SS(r, row, zeta);
            break;
        default:
            fcoulomb = 0;
        }
        double fg = (delta_q-qs)*fcoulomb;
        chi2 += gmx::square(fpol - fg);
    }
    return chi2;
}

static void print_stats(bool verbose, int iter,
                        double chi2, double qs, double zeta)
{
    if (verbose)
    {
        printf("Iter: %5d  chi2: %8g  qs: %7.3f  zeta: %7.3f\n", iter, chi2, qs, zeta);
    }
}

static double zeta0(alexandria::ChargeDistributionModel iChargeDistributionModel,
                    double alpha)
{
    double zeta = 0;
    if (alpha > 0)
    {
        zeta = std::pow(3.0*std::sqrt(M_PI)/(4.0*alpha), 1.0/3.0);
    }
    if (iChargeDistributionModel == alexandria::eqdAXps)
    {
        zeta *= 2.0;
    }
    return zeta;
}

static void fit_polarization(double alpha, double delta_q, double rmax, int row,
                             alexandria::ChargeDistributionModel iChargeDistributionModel,
                             int maxiter, double tolerance,
                             double *zeta_opt, double *qs_opt,
                             bool verbose)
{
    *zeta_opt       = 0;
    *qs_opt         = -4;
    
    if (alpha <= 0)
    {
        return;
    }
     
    double zeta[2] = {   0,   0 };
    double qs[2]   = {  -4,  -4 };
    double chi2[2] = { 1e8, 1e8 };
    int    Min     = 0;
#define Try (1-Min)
     
    gmx::ThreeFry2x64<64>              rng(123457, gmx::RandomDomain::Other);
    gmx::NormalDistribution<real>      dist(1.0, 0.02);
    gmx::UniformRealDistribution<real> mc(0.0, 1.0);
    double qs_min    = -8;
    double qs_max    = -1;
    double zeta_max  = 40;
    double zeta_min  = 2;
    
    int    iter = 0;
    do
    {
        if (iter % 100 == 0)
        {
            zeta[Try] = zeta_min + (zeta_max-zeta_min)*mc(rng);
            qs[Try]   = qs_min   + (qs_max-qs_min)*mc(rng);
        }
        else
        {
            double factor   = dist(rng);
            if (iter % 2)
            {
                zeta[Try] = std::max(zeta_min, std::min(zeta_max, zeta[Min]*factor));
            }
            else
            {
                qs[Try] = std::min(qs_max, std::max(qs_min, qs[Min]*factor));
            }
        }
        chi2[Try] = chi2_coulomb(zeta[Try], qs[Try],
                                 alpha, delta_q, rmax, row, 
                                 iChargeDistributionModel);
        if (chi2[Try] < chi2[Min])
        {
            Min = Try;
            print_stats(verbose, iter, chi2[Min], qs[Min], zeta[Min]);
        }
        iter += 1;
    } 
    while (chi2[Min] > tolerance && iter <= maxiter);

    if (chi2[Min] <= tolerance)
    {
        *zeta_opt = zeta[Min];
        *qs_opt   = qs[Min];
    }
    else
    {
        *zeta_opt = 0;
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
        for(auto ai = pd.getAtypeBegin(); ai < pd.getAtypeEnd(); ++ai)
        {
            auto eem = pd.findEem(iChargeDistributionModel, ai->getType());
            if (pd.EndEemprops() != eem)
            {
                auto pname = ai->getPtype();
                auto zname = ai->getZtype();
                if (true || pname == zname)
                {
                    double zeta, qs;
                    auto ptype = pd.findPtype(pname);
                    fit_polarization(ptype->getPolarizability()/1000, 
                                     0.5, rmax, 
                                     eem->getRow(1),
                                     iChargeDistributionModel, 
                                     maxiter, tolerance,
                                     &zeta, &qs, true);
                    printf("atype %s zeta_old %7g qs_old %3g alpha %7g zeta_new %7g qs_new %7g\n",
                           ai->getType().c_str(), eem->getZeta(1), eem->getQ(1),  
                           ptype->getPolarizability(), zeta, qs);
                    eem->setZeta(1, zeta);
                    eem->setQ(1, qs);
                }
                else
                {
                    fprintf(stderr, "Ignoring atomtype %s because ptype %s differs from ztype %s\n", ai->getType().c_str(), pname.c_str(), zname.c_str());
                }
            }
        }
        writePoldata(opt2fn("-do", NFILE, fnm), pd, true);
    }
    else
    {
        double zeta, qs;
        fit_polarization(alpha, delta_q, rmax, row,
                         iChargeDistributionModel, maxiter, tolerance,
                         &zeta, &qs, true);
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
