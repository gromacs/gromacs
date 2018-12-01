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
                           double qtot, double rmax,
                           int row, 
                           alexandria::ChargeDistributionModel iChargeDistributionModel)
{
    double chi2 = 0;
    int    imax = 10;
    double sqrtpi = std::sqrt(M_PI);
    for(int i = 1; i<=imax; i++)
    {
        double r        = i*(rmax/imax);
        double fpol     = qs*qs*r/alpha;
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
        double fg = qs*(qtot-qs)*fcoulomb;
        chi2 += gmx::square(fpol - fg);
    }
    return chi2;
}

static double fit_polarization(double alpha, double delta_q, double rmax, 
                               alexandria::ChargeDistributionModel iChargeDistributionModel)
{
    double zeta = 0;
    
    if (alpha > 0)
    {
        double zeta0 = std::pow(3.0*std::sqrt(M_PI)/(4.0*alpha), 1.0/3.0);
        if (delta_q == 0)
        {
            zeta = zeta0;
        }
        else
        {
            
        }
    }
    return zeta;
}

int alex_fit_qs_zeta(int argc, char *argv[])
{
    static const char               *desc[] = {
        "Determine distributed charge parameters based on polarizability."
    };
    gmx_output_env_t                *oenv;
    gmx_atomprop_t                   aps = nullptr;
    char                             forcefield[STRLEN];
    char                             ffdir[STRLEN];
    char                             ffname[STRLEN];

    t_filenm                         fnm[] = {
        { efDAT, "-d",        "gentop",    ffOPTRD },
        { efDAT, "-do",       "gentop_out",    ffOPTWR },
        { efTEX, "-t",        "beta",      ffOPTWR }
    };

    const  int         NFILE          = asize(fnm);
    static const char *cqdist[]       = {nullptr, "AXpg", "AXps", nullptr};
    gmx_bool           bVerbose       = TRUE;
    real               alpha          = 0;
    real               delta_q        = 0;
    real               rmax           = 0.005;
    int                row            = 1;
    alexandria::ChargeDistributionModel iChargeDistributionModel;

    t_pargs                          pa[]     = {
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
        double zeta = fit_polarization(alpha, delta_q, rmax, 
                                       iChargeDistributionModel);
        if (zeta > 0)
        {
            printf("qdist = %s alpha = %g delta_q = %g zeta = %g\n",
                   cqdist[0], alpha, delta_q, zeta);
        }
        else
        {
            printf("Do not know how to deal with qdist = %s, alpha = %g delta_q = %g\n",
                   cqdist[0], alpha, delta_q);
        }
    }


    return 0;
}
