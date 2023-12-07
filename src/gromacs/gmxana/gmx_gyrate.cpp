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
#include <cstring>

#include <array>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static real calc_gyro(rvec     x[],
                      int      gnx,
                      int      index[],
                      t_atom   atom[],
                      real     tm,
                      rvec     gvec,
                      rvec     d,
                      gmx_bool bQ,
                      gmx_bool bRot,
                      gmx_bool bMOI,
                      matrix   trans)
{
    int  i, ii, m;
    real gyro, dx2, m0, Itot;
    rvec comp;

    if (bRot)
    {
        principal_comp(gnx, index, atom, x, trans, d);
        Itot = norm(d);
        if (bMOI)
        {
            return Itot;
        }
        for (m = 0; (m < DIM); m++)
        {
            d[m] = std::sqrt(d[m] / tm);
        }
        /* rotate_atoms(gnx,index,x,trans); */
    }
    clear_rvec(comp);
    for (i = 0; (i < gnx); i++)
    {
        ii = index[i];
        if (bQ)
        {
            m0 = std::abs(atom[ii].q);
        }
        else
        {
            m0 = atom[ii].m;
        }
        for (m = 0; (m < DIM); m++)
        {
            dx2 = x[ii][m] * x[ii][m];
            comp[m] += dx2 * m0;
        }
    }
    gyro = comp[XX] + comp[YY] + comp[ZZ];

    for (m = 0; (m < DIM); m++)
    {
        gvec[m] = std::sqrt((gyro - comp[m]) / tm);
    }

    return std::sqrt(gyro / tm);
}

static void calc_gyro_z(rvec x[], matrix box, int gnx, const int index[], t_atom atom[], int nz, real time, FILE* out)
{
    static std::vector<gmx::DVec> inertia;
    static std::vector<double>    tm;
    int                           i, ii, j, zi;
    real                          zf, w, sdet, e1, e2;

    if (inertia.empty())
    {
        inertia.resize(nz);
        tm.resize(nz);
    }

    for (i = 0; i < nz; i++)
    {
        clear_dvec(inertia[i]);
        tm[i] = 0;
    }

    for (i = 0; (i < gnx); i++)
    {
        ii = index[i];
        zf = nz * x[ii][ZZ] / box[ZZ][ZZ];
        if (zf >= nz)
        {
            zf -= nz;
        }
        if (zf < 0)
        {
            zf += nz;
        }
        for (j = 0; j < 2; j++)
        {
            zi = static_cast<int>(zf + j);
            if (zi == nz)
            {
                zi = 0;
            }
            w = atom[ii].m * (1 + std::cos(M_PI * (zf - zi)));
            inertia[zi][0] += w * gmx::square(x[ii][YY]);
            inertia[zi][1] += w * gmx::square(x[ii][XX]);
            inertia[zi][2] -= w * x[ii][XX] * x[ii][YY];
            tm[zi] += w;
        }
    }
    fprintf(out, "%10g", time);
    for (j = 0; j < nz; j++)
    {
        for (i = 0; i < 3; i++)
        {
            inertia[j][i] /= tm[j];
        }
        sdet = std::sqrt(gmx::square(inertia[j][0] - inertia[j][1]) + 4 * gmx::square(inertia[j][2]));
        e1   = std::sqrt(0.5 * (inertia[j][0] + inertia[j][1] + sdet));
        e2   = std::sqrt(0.5 * (inertia[j][0] + inertia[j][1] - sdet));
        fprintf(out, " %5.3f %5.3f", e1, e2);
    }
    fprintf(out, "\n");
}


int gmx_gyrate(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes the radius of gyration of a molecule",
        "and the radii of gyration about the [IT]x[it]-, [IT]y[it]- and [IT]z[it]-axes,",
        "as a function of time. The atoms are explicitly mass weighted.[PAR]",
        "The axis components corresponds to the mass-weighted root-mean-square",
        "of the radii components orthogonal to each axis, for example:[PAR]",
        "Rg(x) = sqrt((sum_i m_i (R_i(y)^2 + R_i(z)^2))/(sum_i m_i)).[PAR]",
        "With the [TT]-nmol[tt] option the radius of gyration will be calculated",
        "for multiple molecules by splitting the analysis group in equally",
        "sized parts.[PAR]",
        "With the option [TT]-nz[tt] 2D radii of gyration in the [IT]x-y[it] plane",
        "of slices along the [IT]z[it]-axis are calculated."
    };
    static int      nmol = 1, nz = 0;
    static gmx_bool bQ = FALSE, bRot = FALSE, bMOI = FALSE;
    t_pargs         pa[] = {
        { "-nmol", FALSE, etINT, { &nmol }, "The number of molecules to analyze" },
        { "-q",
          FALSE,
          etBOOL,
          { &bQ },
          "Use absolute value of the charge of an atom as weighting factor instead of mass" },
        { "-p",
          FALSE,
          etBOOL,
          { &bRot },
          "Calculate the radii of gyration about the principal axes." },
        { "-moi",
          FALSE,
          etBOOL,
          { &bMOI },
          "Calculate the moments of inertia (defined by the principal axes)." },
        { "-nz",
          FALSE,
          etINT,
          { &nz },
          "Calculate the 2D radii of gyration of this number of slices along the z-axis" },
    };
    FILE*                      out;
    t_trxstatus*               status;
    t_topology                 top;
    PbcType                    pbcType;
    rvec *                     x, *x_s;
    rvec                       xcm, gvec, gvec1;
    matrix                     box, trans;
    gmx_bool                   bACF;
    real**                     moi_trans = nullptr;
    int                        max_moi = 0, delta_moi = 100;
    rvec                       d, d1; /* eigenvalues of inertia tensor */
    real                       t, t0, tm, gyro;
    int                        natoms;
    char*                      grpname;
    int                        j, m, gnx, nam, mol;
    int*                       index;
    gmx_output_env_t*          oenv;
    gmx_rmpbc_t                gpbc = nullptr;
    std::array<std::string, 4> leg  = { "Rg", "Rg\\sX\\N", "Rg\\sY\\N", "Rg\\sZ\\N" };
    std::array<std::string, 4> legI = { "Itot", "I1", "I2", "I3" };
#define NLEG asize(leg)
    t_filenm fnm[] = {
        { efTRX, "-f", nullptr, ffREAD },      { efTPS, nullptr, nullptr, ffREAD },
        { efNDX, nullptr, nullptr, ffOPTRD },  { efXVG, nullptr, "gyrate", ffWRITE },
        { efXVG, "-acf", "moi-acf", ffOPTWR },
    };
#define NFILE asize(fnm)
    int      npargs;
    t_pargs* ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);

    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW, NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    std::fprintf(stdout,
                 "You are going to use a deprecated gmx tool. Please migrate to the new one, gmx "
                 "gyrate");

    bACF = opt2bSet("-acf", NFILE, fnm);
    if (bACF && nmol != 1)
    {
        gmx_fatal(FARGS, "Can only do acf with nmol=1");
    }
    bRot = bRot || bMOI || bACF;
    /*
       if (nz > 0)
       bMOI = TRUE;
     */
    if (bRot)
    {
        printf("Will rotate system along principal axes\n");
        snew(moi_trans, DIM);
    }
    if (bMOI)
    {
        printf("Will print moments of inertia\n");
        bQ = FALSE;
    }
    if (bQ)
    {
        printf("Will print radius normalised by charge\n");
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &x, nullptr, box, TRUE);
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &gnx, &index, &grpname);

    if (nmol > gnx || gnx % nmol != 0)
    {
        gmx_fatal(FARGS, "The number of atoms in the group (%d) is not a multiple of nmol (%d)", gnx, nmol);
    }
    nam = gnx / nmol;

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    snew(x_s, natoms);

    j  = 0;
    t0 = t;
    if (bQ)
    {
        out = xvgropen(ftp2fn(efXVG, NFILE, fnm),
                       "Radius of Charge (total and around axes)",
                       "Time (ps)",
                       "Rg (nm)",
                       oenv);
    }
    else if (bMOI)
    {
        out = xvgropen(ftp2fn(efXVG, NFILE, fnm),
                       "Moments of inertia (total and around axes)",
                       "Time (ps)",
                       "I (a.m.u. nm\\S2\\N)",
                       oenv);
    }
    else
    {
        out = xvgropen(ftp2fn(efXVG, NFILE, fnm),
                       "Radius of gyration (total and around axes)",
                       "Time (ps)",
                       "Rg (nm)",
                       oenv);
    }
    if (bMOI)
    {
        xvgrLegend(out, legI, oenv);
    }
    else
    {
        if (bRot)
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(out, "@ subtitle \"Axes are principal component axes\"\n");
            }
        }
        xvgrLegend(out, leg, oenv);
    }
    if (nz == 0)
    {
        gpbc = gmx_rmpbc_init(&top.idef, pbcType, natoms);
    }
    do
    {
        if (nz == 0)
        {
            gmx_rmpbc_copy(gpbc, natoms, box, x, x_s);
        }
        gyro = 0;
        clear_rvec(gvec);
        clear_rvec(gvec1);
        clear_rvec(d);
        clear_rvec(d1);
        for (mol = 0; mol < nmol; mol++)
        {
            tm = sub_xcm(nz == 0 ? x_s : x, nam, index + mol * nam, top.atoms.atom, xcm, bQ);
            if (nz == 0)
            {
                gyro += calc_gyro(
                        x_s, nam, index + mol * nam, top.atoms.atom, tm, gvec1, d1, bQ, bRot, bMOI, trans);
            }
            else
            {
                calc_gyro_z(x, box, nam, index + mol * nam, top.atoms.atom, nz, t, out);
            }
            rvec_inc(gvec, gvec1);
            rvec_inc(d, d1);
        }
        if (nmol > 0)
        {
            gyro /= nmol;
            svmul(1.0 / nmol, gvec, gvec);
            svmul(1.0 / nmol, d, d);
        }

        if (nz == 0)
        {
            if (bRot)
            {
                if (j >= max_moi)
                {
                    max_moi += delta_moi;
                    for (m = 0; (m < DIM); m++)
                    {
                        srenew(moi_trans[m], max_moi * DIM);
                    }
                }
                for (m = 0; (m < DIM); m++)
                {
                    copy_rvec(trans[m], moi_trans[m] + DIM * j);
                }
                fprintf(out, "%10g  %10g  %10g  %10g  %10g\n", t, gyro, d[XX], d[YY], d[ZZ]);
            }
            else
            {
                fprintf(out, "%10g  %10g  %10g  %10g  %10g\n", t, gyro, gvec[XX], gvec[YY], gvec[ZZ]);
            }
        }
        j++;
    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);
    if (nz == 0)
    {
        gmx_rmpbc_done(gpbc);
    }

    xvgrclose(out);

    if (bACF)
    {
        int mode = eacVector;

        do_autocorr(opt2fn("-acf", NFILE, fnm),
                    oenv,
                    "Moment of inertia vector ACF",
                    j,
                    3,
                    moi_trans,
                    (t - t0) / j,
                    mode,
                    FALSE);
        do_view(oenv, opt2fn("-acf", NFILE, fnm), "-nxy");
    }

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    return 0;
}
