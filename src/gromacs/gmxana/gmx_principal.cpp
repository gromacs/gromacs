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

#include <cstdio>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

enum class PbcType : int;
struct gmx_output_env_t;


static void calc_principal_axes(const t_topology* top, rvec* x, int* index, int n, matrix axes, rvec inertia)
{
    rvec xcm;

    sub_xcm(x, n, index, top->atoms.atom, xcm, FALSE);
    principal_comp(n, index, top->atoms.atom, x, axes, inertia);
}

int gmx_principal(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] calculates the three principal axes of inertia for a group",
        "of atoms. NOTE: Old versions of GROMACS wrote the output data in a",
        "strange transposed way. As of GROMACS 5.0, the output file paxis1.dat",
        "contains the x/y/z components of the first (major) principal axis for",
        "each frame, and similarly for the middle and minor axes in paxis2.dat",
        "and paxis3.dat."
    };
    static gmx_bool foo = FALSE;

    t_pargs pa[] = { { "-foo", FALSE, etBOOL, { &foo }, "Dummy option to avoid empty array" } };
    t_trxstatus* status;
    t_topology   top;
    PbcType      pbcType;
    real         t;
    rvec*        x;

    int                      natoms;
    char*                    grpname;
    int                      i, gnx;
    int*                     index;
    rvec                     moi;
    FILE*                    axis1;
    FILE*                    axis2;
    FILE*                    axis3;
    FILE*                    fmoi;
    matrix                   axes, box;
    gmx_output_env_t*        oenv;
    gmx_rmpbc_t              gpbc = nullptr;
    std::vector<std::string> legend;

    t_filenm fnm[] = { { efTRX, "-f", nullptr, ffREAD },     { efTPS, nullptr, nullptr, ffREAD },
                       { efNDX, nullptr, nullptr, ffOPTRD }, { efXVG, "-a1", "paxis1", ffWRITE },
                       { efXVG, "-a2", "paxis2", ffWRITE },  { efXVG, "-a3", "paxis3", ffWRITE },
                       { efXVG, "-om", "moi", ffWRITE } };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc,
                           argv,
                           PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW,
                           NFILE,
                           fnm,
                           asize(pa),
                           pa,
                           asize(desc),
                           desc,
                           0,
                           nullptr,
                           &oenv))
    {
        return 0;
    }

    for (i = 0; i < DIM; i++)
    {
        legend.emplace_back(gmx::formatString("%c component", 'X' + i));
    }

    axis1 = xvgropen(opt2fn("-a1", NFILE, fnm),
                     "Principal axis 1 (major axis)",
                     output_env_get_xvgr_tlabel(oenv),
                     "Component (nm)",
                     oenv);
    xvgrLegend(axis1, legend, oenv);

    axis2 = xvgropen(opt2fn("-a2", NFILE, fnm),
                     "Principal axis 2 (middle axis)",
                     output_env_get_xvgr_tlabel(oenv),
                     "Component (nm)",
                     oenv);
    xvgrLegend(axis2, legend, oenv);

    axis3 = xvgropen(opt2fn("-a3", NFILE, fnm),
                     "Principal axis 3 (minor axis)",
                     output_env_get_xvgr_tlabel(oenv),
                     "Component (nm)",
                     oenv);
    xvgrLegend(axis3, legend, oenv);

    legend[XX] = "Axis 1 (major)";
    legend[YY] = "Axis 2 (middle)";
    legend[ZZ] = "Axis 3 (minor)";

    fmoi = xvgropen(opt2fn("-om", NFILE, fnm),
                    "Moments of inertia around inertial axes",
                    output_env_get_xvgr_tlabel(oenv),
                    "I (au nm\\S2\\N)",
                    oenv);
    xvgrLegend(fmoi, legend, oenv);

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, nullptr, nullptr, box, TRUE);

    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &gnx, &index, &grpname);

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    gpbc = gmx_rmpbc_init(&top.idef, pbcType, natoms);

    do
    {
        gmx_rmpbc_apply(gpbc, natoms, box, x);

        calc_principal_axes(&top, x, index, gnx, axes, moi);

        fprintf(axis1, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, axes[XX][XX], axes[XX][YY], axes[XX][ZZ]);
        fprintf(axis2, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, axes[YY][XX], axes[YY][YY], axes[YY][ZZ]);
        fprintf(axis3, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, axes[ZZ][XX], axes[ZZ][YY], axes[ZZ][ZZ]);
        fprintf(fmoi, "%15.10f     %15.10f  %15.10f  %15.10f\n", t, moi[XX], moi[YY], moi[ZZ]);
    } while (read_next_x(oenv, status, &t, x, box));

    gmx_rmpbc_done(gpbc);

    close_trx(status);

    xvgrclose(axis1);
    xvgrclose(axis2);
    xvgrclose(axis3);
    xvgrclose(fmoi);

    return 0;
}
