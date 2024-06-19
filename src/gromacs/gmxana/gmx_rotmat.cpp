/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
#include <cstring>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;
struct gmx_output_env_t;

static void get_refx(gmx_output_env_t* oenv,
                     const char*       trxfn,
                     int               nfitdim,
                     int               skip,
                     int               gnx,
                     int*              index,
                     gmx_bool          bMW,
                     const t_topology* top,
                     PbcType           pbcType,
                     rvec*             x_ref)
{
    int          natoms, nfr_all, nfr, i, j, a, r, c, min_fr;
    t_trxstatus* status;
    real *       ti, min_t;
    double       tot_mass, msd, *srmsd, min_srmsd, srmsd_tot;
    rvec *       x, **xi;
    real         xf;
    matrix       box, R;
    real*        w_rls;
    gmx_rmpbc_t  gpbc = nullptr;


    nfr_all = 0;
    nfr     = 0;
    snew(ti, 100);
    snew(xi, 100);
    natoms = read_first_x(oenv, &status, trxfn, &ti[nfr], &x, box);

    snew(w_rls, gnx);
    tot_mass = 0;
    for (a = 0; a < gnx; a++)
    {
        if (index[a] >= natoms)
        {
            gmx_fatal(FARGS,
                      "Atom index (%d) is larger than the number of atoms in the trajecory (%d)",
                      index[a] + 1,
                      natoms);
        }
        w_rls[a] = (bMW ? top->atoms.atom[index[a]].m : 1.0);
        tot_mass += w_rls[a];
    }
    gpbc = gmx_rmpbc_init(&top->idef, pbcType, natoms);

    do
    {
        if (nfr_all % skip == 0)
        {
            gmx_rmpbc_apply(gpbc, natoms, box, x);
            snew(xi[nfr], gnx);
            for (i = 0; i < gnx; i++)
            {
                copy_rvec(x[index[i]], xi[nfr][i]);
            }
            reset_x(gnx, nullptr, gnx, nullptr, xi[nfr], w_rls);
            nfr++;
            if (nfr % 100 == 0)
            {
                srenew(ti, nfr + 100);
                srenew(xi, nfr + 100);
            }
        }
        nfr_all++;
    } while (read_next_x(oenv, status, &ti[nfr], x, box));
    close_trx(status);
    sfree(x);

    gmx_rmpbc_done(gpbc);

    snew(srmsd, nfr);
    for (i = 0; i < nfr; i++)
    {
        fprintf(stdout, "\rProcessing frame %d of %d", i, nfr);
        fflush(stdout);
        for (j = i + 1; j < nfr; j++)
        {
            calc_fit_R(nfitdim, gnx, w_rls, xi[i], xi[j], R);

            msd = 0;
            for (a = 0; a < gnx; a++)
            {
                for (r = 0; r < DIM; r++)
                {
                    xf = 0;
                    for (c = 0; c < DIM; c++)
                    {
                        xf += R[r][c] * xi[j][a][c];
                    }
                    msd += w_rls[a] * gmx::square(xi[i][a][r] - xf);
                }
            }
            msd /= tot_mass;
            srmsd[i] += std::sqrt(msd);
            srmsd[j] += std::sqrt(msd);
        }
        sfree(xi[i]);
    }
    printf("\n");
    sfree(w_rls);

    min_srmsd = GMX_REAL_MAX;
    min_fr    = -1;
    min_t     = -1;
    srmsd_tot = 0;
    for (i = 0; i < nfr; i++)
    {
        srmsd[i] /= (nfr - 1);
        if (srmsd[i] < min_srmsd)
        {
            min_srmsd = srmsd[i];
            min_fr    = i;
            min_t     = ti[i];
        }
        srmsd_tot += srmsd[i];
    }
    sfree(srmsd);

    printf("Average RMSD between all structures: %.3f\n", srmsd_tot / nfr);
    printf("Structure with lowest RMSD to all others: time %g, av. RMSD %.3f\n", min_t, min_srmsd);

    for (a = 0; a < gnx; a++)
    {
        copy_rvec(xi[min_fr][a], x_ref[index[a]]);
    }

    sfree(xi);
}

int gmx_rotmat(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] plots the rotation matrix required for least squares fitting",
        "a conformation onto the reference conformation provided with",
        "[TT]-s[tt]. Translation is removed before fitting.",
        "The output are the three vectors that give the new directions",
        "of the x, y and z directions of the reference conformation,",
        "for example: (zx,zy,zz) is the orientation of the reference",
        "z-axis in the trajectory frame.",
        "[PAR]",
        "This tool is useful for, for instance,",
        "determining the orientation of a molecule",
        "at an interface, possibly on a trajectory produced with",
        "[TT]gmx trjconv -fit rotxy+transxy[tt] to remove the rotation",
        "in the [IT]x-y[it] plane.",
        "[PAR]",
        "Option [TT]-ref[tt] determines a reference structure for fitting,",
        "instead of using the structure from [TT]-s[tt]. The structure with",
        "the lowest sum of RMSD's to all other structures is used.",
        "Since the computational cost of this procedure grows with",
        "the square of the number of frames, the [TT]-skip[tt] option",
        "can be useful. A full fit or only a fit in the [IT]x-y[it] plane can",
        "be performed.",
        "[PAR]",
        "Option [TT]-fitxy[tt] fits in the [IT]x-y[it] plane before determining",
        "the rotation matrix."
    };
    const char*     reffit[] = { nullptr, "none", "xyz", "xy", nullptr };
    static int      skip     = 1;
    static gmx_bool bFitXY = FALSE, bMW = TRUE;
    t_pargs         pa[] = {
        { "-ref", FALSE, etENUM, { reffit }, "Determine the optimal reference structure" },
        { "-skip", FALSE, etINT, { &skip }, "Use every nr-th frame for [TT]-ref[tt]" },
        { "-fitxy",
          FALSE,
          etBOOL,
          { &bFitXY },
          "Fit the x/y rotation before determining the rotation" },
        { "-mw", FALSE, etBOOL, { &bMW }, "Use mass weighted fitting" }
    };
    FILE*                      out;
    t_trxstatus*               status;
    t_topology                 top;
    PbcType                    pbcType;
    rvec *                     x_ref, *x;
    matrix                     box, R;
    real                       t;
    int                        natoms, i;
    char*                      grpname;
    int                        gnx;
    gmx_rmpbc_t                gpbc = nullptr;
    int*                       index;
    gmx_output_env_t*          oenv;
    real*                      w_rls;
    std::array<std::string, 9> leg = { "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz" };
#define NLEG asize(leg)
    t_filenm fnm[] = { { efTRX, "-f", nullptr, ffREAD },
                       { efTPS, nullptr, nullptr, ffREAD },
                       { efNDX, nullptr, nullptr, ffOPTRD },
                       { efXVG, nullptr, "rotmat", ffWRITE } };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &x_ref, nullptr, box, bMW);

    gpbc = gmx_rmpbc_init(&top.idef, pbcType, top.atoms.nr);

    gmx_rmpbc_apply(gpbc, top.atoms.nr, box, x_ref);

    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &gnx, &index, &grpname);

    GMX_RELEASE_ASSERT(reffit[0] != nullptr, "Options inconsistency; reffit[0] is NULL");
    if (reffit[0][0] != 'n')
    {
        get_refx(oenv, ftp2fn(efTRX, NFILE, fnm), reffit[0][2] == 'z' ? 3 : 2, skip, gnx, index, bMW, &top, pbcType, x_ref);
    }

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    snew(w_rls, natoms);
    for (i = 0; i < gnx; i++)
    {
        if (index[i] >= natoms)
        {
            gmx_fatal(FARGS,
                      "Atom index (%d) is larger than the number of atoms in the trajecory (%d)",
                      index[i] + 1,
                      natoms);
        }
        w_rls[index[i]] = (bMW ? top.atoms.atom[index[i]].m : 1.0);
    }

    if (reffit[0][0] == 'n')
    {
        reset_x(gnx, index, natoms, nullptr, x_ref, w_rls);
    }

    out = xvgropen(ftp2fn(efXVG, NFILE, fnm), "Fit matrix", "Time (ps)", "", oenv);
    xvgrLegend(out, leg, oenv);

    do
    {
        gmx_rmpbc_apply(gpbc, natoms, box, x);

        reset_x(gnx, index, natoms, nullptr, x, w_rls);

        if (bFitXY)
        {
            do_fit_ndim(2, natoms, w_rls, x_ref, x);
        }

        calc_fit_R(DIM, natoms, w_rls, x_ref, x, R);

        fprintf(out,
                "%7g %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                t,
                R[XX][XX],
                R[XX][YY],
                R[XX][ZZ],
                R[YY][XX],
                R[YY][YY],
                R[YY][ZZ],
                R[ZZ][XX],
                R[ZZ][YY],
                R[ZZ][ZZ]);
    } while (read_next_x(oenv, status, &t, x, box));

    gmx_rmpbc_done(gpbc);

    close_trx(status);

    xvgrclose(out);

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    return 0;
}
