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

#include <cmath>
#include <cstdio>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;

static real pot(real x, real qq, real c6, real cn, int npow)
{
    return cn * std::pow(x, -npow) - c6 / gmx::power6(x) + qq * gmx::c_one4PiEps0 / x;
}

static real bhpot(real x, real A, real B, real C)
{
    return A * std::exp(-B * x) - C / gmx::power6(x);
}

static real dpot(real x, real qq, real c6, real cn, int npow)
{
    return -(npow * cn * std::pow(x, -npow - 1) - 6 * c6 / (x * gmx::power6(x))
             + qq * gmx::c_one4PiEps0 / gmx::square(x));
}

int gmx_sigeps(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] is a simple utility that converts C6/C12 or C6/Cn combinations",
        "to [GRK]sigma[grk] and [GRK]epsilon[grk], or vice versa. It can also plot the potential",
        "in  file. In addition, it makes an approximation of a Buckingham potential",
        "to a Lennard-Jones potential."
    };
    static real       c6 = 1.0e-3, cn = 1.0e-6, qi = 0, qj = 0, sig = 0.3, eps = 1, sigfac = 0.7;
    static real       Abh = 1e5, Bbh = 32, Cbh = 1e-3;
    static int        npow  = 12;
    t_pargs           pa[]  = { { "-c6", FALSE, etREAL, { &c6 }, "C6" },
                     { "-cn", FALSE, etREAL, { &cn }, "Constant for repulsion" },
                     { "-pow", FALSE, etINT, { &npow }, "Power of the repulsion term" },
                     { "-sig", FALSE, etREAL, { &sig }, "[GRK]sigma[grk]" },
                     { "-eps", FALSE, etREAL, { &eps }, "[GRK]epsilon[grk]" },
                     { "-A", FALSE, etREAL, { &Abh }, "Buckingham A" },
                     { "-B", FALSE, etREAL, { &Bbh }, "Buckingham B" },
                     { "-C", FALSE, etREAL, { &Cbh }, "Buckingham C" },
                     { "-qi", FALSE, etREAL, { &qi }, "qi" },
                     { "-qj", FALSE, etREAL, { &qj }, "qj" },
                     { "-sigfac",
                       FALSE,
                       etREAL,
                       { &sigfac },
                       "Factor in front of [GRK]sigma[grk] for starting the plot" } };
    t_filenm          fnm[] = { { efXVG, "-o", "potje", ffWRITE } };
    gmx_output_env_t* oenv;
#define NFILE asize(fnm)
    std::array<std::string, 2> legend = { "Lennard-Jones", "Buckingham" };
    FILE*                      fp;
    int                        i;
    gmx_bool                   bBham;
    real                       qq, x, oldx, minimum, mval, dp[2];
    int                        cur = 0;
#define next (1 - cur)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    bBham = (opt2parg_bSet("-A", asize(pa), pa) || opt2parg_bSet("-B", asize(pa), pa)
             || opt2parg_bSet("-C", asize(pa), pa));

    if (bBham)
    {
        c6  = Cbh;
        sig = std::pow((6.0 / npow) * std::pow(npow / Bbh, npow - 6), 1.0 / (npow - 6));
        eps = c6 / (4 * gmx::power6(sig));
        cn  = 4 * eps * std::pow(sig, npow);
    }
    else
    {
        if (opt2parg_bSet("-sig", asize(pa), pa) || opt2parg_bSet("-eps", asize(pa), pa))
        {
            c6 = 4 * eps * gmx::power6(sig);
            cn = 4 * eps * std::pow(sig, npow);
        }
        else if (opt2parg_bSet("-c6", asize(pa), pa) || opt2parg_bSet("-cn", asize(pa), pa)
                 || opt2parg_bSet("-pow", asize(pa), pa))
        {
            sig = std::pow(cn / c6, static_cast<real>(1.0 / (npow - 6)));
            eps = 0.25 * c6 / gmx::power6(sig);
        }
        else
        {
            sig = eps = 0;
        }
        printf("c6    = %12.5e, c%d    = %12.5e\n", c6, npow, cn);
        printf("sigma = %12.5f, epsilon = %12.5f\n", sig, eps);

        minimum = std::pow(npow / 6.0 * std::pow(sig, npow - 6), 1.0 / (npow - 6));
        printf("Van der Waals minimum at %g, V = %g\n\n", minimum, pot(minimum, 0, c6, cn, npow));
        printf("Fit of Lennard Jones (%d-6) to Buckingham:\n", npow);
        Bbh = npow / minimum;
        Cbh = c6;
        Abh = 4 * eps * std::pow(sig / minimum, static_cast<real>(npow))
              * std::exp(static_cast<real>(npow));
        printf("A = %g, B = %g, C = %g\n", Abh, Bbh, Cbh);
    }
    qq = qi * qj;

    fp = xvgropen(ftp2fn(efXVG, NFILE, fnm), "Potential", "r (nm)", "E (kJ/mol)", oenv);
    xvgrLegend(fp, legend, oenv);
    if (sig == 0)
    {
        sig = 0.25;
    }
    oldx = 0;
    for (i = 0; (i < 100); i++)
    {
        x        = sigfac * sig + sig * i * 0.02;
        dp[next] = dpot(x, qq, c6, cn, npow);
        fprintf(fp, "%10g  %10g  %10g\n", x, pot(x, qq, c6, cn, npow), bhpot(x, Abh, Bbh, Cbh));
        if (qq != 0)
        {
            if ((i > 0) && (dp[cur] * dp[next] < 0))
            {
                minimum = oldx + dp[cur] * (x - oldx) / (dp[cur] - dp[next]);
                mval    = pot(minimum, qq, c6, cn, npow);
                printf("Van der Waals + Coulomb minimum at r = %g (nm). Value = %g (kJ/mol)\n", minimum, mval);
            }
        }
        cur  = next;
        oldx = x;
    }
    xvgrclose(fp);

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), nullptr);

    return 0;
}
