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

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/rgb.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/eigio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "thermochemistry.h"

struct gmx_output_env_t;

static real tick_spacing(real range, int minticks)
{
    real sp;

    if (range <= 0)
    {
        return 1.0;
    }

    sp = 0.2 * std::exp(std::log(10.0) * std::ceil(std::log(range) / std::log(10.0)));
    while (range / sp < minticks - 1)
    {
        sp = sp / 2;
    }

    return sp;
}

static void write_xvgr_graphs(const char*             file,
                              int                     ngraphs,
                              int                     nsetspergraph,
                              const char*             title,
                              const char*             subtitle,
                              const std::string&      xlabel,
                              const char**            ylabel,
                              int                     n,
                              real*                   x,
                              real**                  y,
                              real***                 sy,
                              real                    scale_x,
                              gmx_bool                bZero,
                              gmx_bool                bSplit,
                              const gmx_output_env_t* oenv)
{
    FILE* out;
    int   g, s, i;
    real  ymin, ymax, xsp, ysp;

    out = gmx_ffopen(file, "w");
    if (output_env_get_xvg_format(oenv) == XvgFormat::Xmgrace)
    {
        fprintf(out, "@ autoscale onread none\n");
    }
    for (g = 0; g < ngraphs; g++)
    {
        if (y)
        {
            ymin = y[g][0];
            ymax = y[g][0];
            for (i = 0; i < n; i++)
            {
                if (y[g][i] < ymin)
                {
                    ymin = y[g][i];
                }
                if (y[g][i] > ymax)
                {
                    ymax = y[g][i];
                }
            }
        }
        else
        {
            assert(sy);
            ymin = sy[g][0][0];
            ymax = sy[g][0][0];
            for (s = 0; s < nsetspergraph; s++)
            {
                for (i = 0; i < n; i++)
                {
                    if (sy[g][s][i] < ymin)
                    {
                        ymin = sy[g][s][i];
                    }
                    if (sy[g][s][i] > ymax)
                    {
                        ymax = sy[g][s][i];
                    }
                }
            }
        }
        if (bZero)
        {
            ymin = 0;
        }
        else
        {
            ymin = ymin - 0.1 * (ymax - ymin);
        }
        ymax = ymax + 0.1 * (ymax - ymin);
        xsp  = tick_spacing((x[n - 1] - x[0]) * scale_x, 4);
        ysp  = tick_spacing(ymax - ymin, 3);
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(out, "@ with g%d\n@ g%d on\n", g, g);
            if (g == 0)
            {
                fprintf(out, "@ title \"%s\"\n", title);
                if (subtitle)
                {
                    fprintf(out, "@ subtitle \"%s\"\n", subtitle);
                }
            }
            if (g == ngraphs - 1)
            {
                fprintf(out, "@ xaxis  label \"%s\"\n", xlabel.c_str());
            }
            else
            {
                fprintf(out, "@ xaxis  ticklabel off\n");
            }
            if (n > 1)
            {
                fprintf(out, "@ world xmin %g\n", x[0] * scale_x);
                fprintf(out, "@ world xmax %g\n", x[n - 1] * scale_x);
                fprintf(out, "@ world ymin %g\n", ymin);
                fprintf(out, "@ world ymax %g\n", ymax);
            }
            fprintf(out, "@ view xmin 0.15\n");
            fprintf(out, "@ view xmax 0.85\n");
            fprintf(out, "@ view ymin %g\n", 0.15 + (ngraphs - 1 - g) * 0.7 / ngraphs);
            fprintf(out, "@ view ymax %g\n", 0.15 + (ngraphs - g) * 0.7 / ngraphs);
            fprintf(out, "@ yaxis  label \"%s\"\n", ylabel[g]);
            fprintf(out, "@ xaxis tick major %g\n", xsp);
            fprintf(out, "@ xaxis tick minor %g\n", xsp / 2);
            fprintf(out, "@ xaxis ticklabel start type spec\n");
            fprintf(out, "@ xaxis ticklabel start %g\n", std::ceil(ymin / xsp) * xsp);
            fprintf(out, "@ yaxis tick major %g\n", ysp);
            fprintf(out, "@ yaxis tick minor %g\n", ysp / 2);
            fprintf(out, "@ yaxis ticklabel start type spec\n");
            fprintf(out, "@ yaxis ticklabel start %g\n", std::ceil(ymin / ysp) * ysp);
            if ((ymin < 0) && (ymax > 0))
            {
                fprintf(out, "@ zeroxaxis bar on\n");
                fprintf(out, "@ zeroxaxis bar linestyle 3\n");
            }
        }
        for (s = 0; s < nsetspergraph; s++)
        {
            for (i = 0; i < n; i++)
            {
                if (bSplit && i > 0 && std::abs(x[i]) < 1e-5)
                {
                    fprintf(out, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
                }
                fprintf(out, "%10.4f %10.5f\n", x[i] * scale_x, y ? y[g][i] : sy[g][s][i]);
            }
            fprintf(out, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
        }
    }
    gmx_ffclose(out);
}

static void
compare(int natoms, int n1, rvec** eigvec1, int n2, rvec** eigvec2, real* eigval1, int neig1, real* eigval2, int neig2)
{
    int    n;
    int    i, j, k;
    double sum1, sum2, trace1, trace2, sab, samsb2, tmp, ip;

    n = std::min(n1, n2);

    n = std::min(n, std::min(neig1, neig2));
    fprintf(stdout, "Will compare the covariance matrices using %d dimensions\n", n);

    sum1 = 0;
    for (i = 0; i < n; i++)
    {
        if (eigval1[i] < 0)
        {
            eigval1[i] = 0;
        }
        sum1 += eigval1[i];
        eigval1[i] = std::sqrt(eigval1[i]);
    }
    trace1 = sum1;
    for (i = n; i < neig1; i++)
    {
        trace1 += eigval1[i];
    }
    sum2 = 0;
    for (i = 0; i < n; i++)
    {
        if (eigval2[i] < 0)
        {
            eigval2[i] = 0;
        }
        sum2 += eigval2[i];
        eigval2[i] = std::sqrt(eigval2[i]);
    }
    trace2 = sum2;


    // The static analyzer appears to be confused by the fact that the loop below
    // starts from n instead of 0. However, given all the complex code it's
    // better to be safe than sorry, so we check it with an assert.
    // If we are in this comparison routine in the first place, neig2 should not be 0,
    // so eigval2 should always be a valid pointer.
    GMX_RELEASE_ASSERT(eigval2 != nullptr, "NULL pointer provided for eigval2");

    for (i = n; i < neig2; i++)
    {
        trace2 += eigval2[i];
    }

    fprintf(stdout, "Trace of the two matrices: %g and %g\n", sum1, sum2);
    if (neig1 != n || neig2 != n)
    {
        fprintf(stdout,
                "this is %d%% and %d%% of the total trace\n",
                gmx::roundToInt(100 * sum1 / trace1),
                gmx::roundToInt(100 * sum2 / trace2));
    }
    fprintf(stdout, "Square root of the traces: %g and %g\n", std::sqrt(sum1), std::sqrt(sum2));

    sab = 0;
    for (i = 0; i < n; i++)
    {
        tmp = 0;
        for (j = 0; j < n; j++)
        {
            ip = 0;
            for (k = 0; k < natoms; k++)
            {
                ip += iprod(eigvec1[i][k], eigvec2[j][k]);
            }
            tmp += eigval2[j] * ip * ip;
        }
        sab += eigval1[i] * tmp;
    }

    samsb2 = sum1 + sum2 - 2 * sab;
    if (samsb2 < 0)
    {
        samsb2 = 0;
    }

    fprintf(stdout, "The overlap of the covariance matrices:\n");
    fprintf(stdout, "  normalized:  %.3f\n", 1 - std::sqrt(samsb2 / (sum1 + sum2)));
    tmp = 1 - sab / std::sqrt(sum1 * sum2);
    if (tmp < 0)
    {
        tmp = 0;
    }
    fprintf(stdout, "       shape:  %.3f\n", 1 - std::sqrt(tmp));
}


static void inprod_matrix(const char* matfile,
                          int         natoms,
                          int         nvec1,
                          int*        eignr1,
                          rvec**      eigvec1,
                          int         nvec2,
                          const int*  eignr2,
                          rvec**      eigvec2,
                          gmx_bool    bSelect,
                          int         noutvec,
                          const int*  outvec)
{
    FILE*  out;
    real** mat;
    int    i, x1, y1, x, y, nlevels;
    int    nx, ny;
    real   inp, *t_x, *t_y, maxval;
    t_rgb  rlo, rhi;

    snew(t_y, nvec2);
    if (bSelect)
    {
        nx = noutvec;
        ny = 0;
        for (y1 = 0; y1 < nx; y1++)
        {
            if (outvec[y1] < nvec2)
            {
                t_y[ny] = eignr2[outvec[y1]] + 1;
                ny++;
            }
        }
    }
    else
    {
        nx = nvec1;
        ny = nvec2;
        for (y = 0; y < ny; y++)
        {
            t_y[y] = eignr2[y] + 1;
        }
    }

    fprintf(stderr, "Calculating inner-product matrix of %dx%d eigenvectors\n", nx, nvec2);

    snew(mat, nx);
    snew(t_x, nx);
    maxval = 0;
    for (x1 = 0; x1 < nx; x1++)
    {
        snew(mat[x1], ny);
        if (bSelect)
        {
            x = outvec[x1];
        }
        else
        {
            x = x1;
        }
        t_x[x1] = eignr1[x] + 1;
        fprintf(stderr, " %d", eignr1[x] + 1);
        for (y1 = 0; y1 < ny; y1++)
        {
            inp = 0;
            if (bSelect)
            {
                while (outvec[y1] >= nvec2)
                {
                    y1++;
                }
                y = outvec[y1];
            }
            else
            {
                y = y1;
            }
            for (i = 0; i < natoms; i++)
            {
                inp += iprod(eigvec1[x][i], eigvec2[y][i]);
            }
            mat[x1][y1] = std::abs(inp);
            if (mat[x1][y1] > maxval)
            {
                maxval = mat[x1][y1];
            }
        }
    }
    fprintf(stderr, "\n");
    rlo.r   = 1;
    rlo.g   = 1;
    rlo.b   = 1;
    rhi.r   = 0;
    rhi.g   = 0;
    rhi.b   = 0;
    nlevels = 41;
    out     = gmx_ffopen(matfile, "w");
    write_xpm(out, 0, "Eigenvector inner-products", "in.prod.", "run 1", "run 2", nx, ny, t_x, t_y, mat, 0.0, maxval, rlo, rhi, &nlevels);
    gmx_ffclose(out);
}

static void overlap(const char*             outfile,
                    int                     natoms,
                    rvec**                  eigvec1,
                    int                     nvec2,
                    int*                    eignr2,
                    rvec**                  eigvec2,
                    int                     noutvec,
                    int*                    outvec,
                    const gmx_output_env_t* oenv)
{
    FILE* out;
    int   i, v, vec, x;
    real  overlap, inp;

    fprintf(stderr, "Calculating overlap between eigenvectors of set 2 with eigenvectors\n");
    for (i = 0; i < noutvec; i++)
    {
        fprintf(stderr, "%d ", outvec[i] + 1);
    }
    fprintf(stderr, "\n");

    out = xvgropen(outfile, "Subspace overlap", "Eigenvectors of trajectory 2", "Overlap", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(out, "@ subtitle \"using %d eigenvectors of trajectory 1\"\n", noutvec);
    }
    overlap = 0;
    for (x = 0; x < nvec2; x++)
    {
        for (v = 0; v < noutvec; v++)
        {
            vec = outvec[v];
            inp = 0;
            for (i = 0; i < natoms; i++)
            {
                inp += iprod(eigvec1[vec][i], eigvec2[x][i]);
            }
            overlap += gmx::square(inp);
        }
        fprintf(out, "%5d  %5.3f\n", eignr2[x] + 1, overlap / noutvec);
    }

    xvgrclose(out);
}

static void project(const char*             trajfile,
                    const t_topology*       top,
                    PbcType                 pbcType,
                    matrix                  topbox,
                    const char*             projfile,
                    const char*             twodplotfile,
                    const char*             threedplotfile,
                    const char*             filterfile,
                    const char*             projUnit,
                    int                     skip,
                    const char*             extremefile,
                    gmx_bool                bExtrAll,
                    real                    extreme,
                    int                     nextr,
                    const t_atoms*          atoms,
                    int                     natoms,
                    int*                    index,
                    gmx_bool                bFit,
                    rvec*                   xref,
                    int                     nfit,
                    int*                    ifit,
                    real*                   w_rls,
                    const real*             sqrtm,
                    rvec*                   xav,
                    int*                    eignr,
                    rvec**                  eigvec,
                    int                     noutvec,
                    int*                    outvec,
                    gmx_bool                bSplit,
                    const gmx_output_env_t* oenv)
{
    FILE*        xvgrout = nullptr;
    int          nat, i, j, d, v, vec, nfr, nframes = 0, snew_size, frame;
    t_trxstatus* out = nullptr;
    t_trxstatus* status;
    int          noutvec_extr, imin, imax;
    real *       pmin, *pmax;
    int*         all_at;
    matrix       box;
    rvec *       xread, *x;
    real         t, inp, **inprod = nullptr;
    char         str[STRLEN], str2[STRLEN];
    const char** ylabel;
    real         fact;
    gmx_rmpbc_t  gpbc = nullptr;

    snew(x, natoms);

    if (bExtrAll)
    {
        noutvec_extr = noutvec;
    }
    else
    {
        noutvec_extr = 1;
    }


    if (trajfile)
    {
        snew(inprod, noutvec + 1);

        if (filterfile)
        {
            fprintf(stderr, "Writing a filtered trajectory to %s using eigenvectors\n", filterfile);
            for (i = 0; i < noutvec; i++)
            {
                fprintf(stderr, "%d ", outvec[i] + 1);
            }
            fprintf(stderr, "\n");
            out = open_trx(filterfile, "w");
        }
        snew_size = 0;
        nfr       = 0;
        nframes   = 0;
        nat       = read_first_x(oenv, &status, trajfile, &t, &xread, box);
        if (nat > atoms->nr)
        {
            gmx_fatal(FARGS,
                      "the number of atoms in your trajectory (%d) is larger than the number of "
                      "atoms in your structure file (%d)",
                      nat,
                      atoms->nr);
        }
        snew(all_at, nat);

        if (top)
        {
            gpbc = gmx_rmpbc_init(&top->idef, pbcType, nat);
        }

        for (i = 0; i < nat; i++)
        {
            all_at[i] = i;
        }
        do
        {
            if (nfr % skip == 0)
            {
                if (top)
                {
                    gmx_rmpbc_apply(gpbc, nat, box, xread);
                }
                if (nframes >= snew_size)
                {
                    snew_size += 100;
                    for (i = 0; i < noutvec + 1; i++)
                    {
                        srenew(inprod[i], snew_size);
                    }
                }
                inprod[noutvec][nframes] = t;
                /* calculate x: a fitted struture of the selected atoms */
                if (bFit)
                {
                    reset_x(nfit, ifit, nat, nullptr, xread, w_rls);
                    do_fit(nat, w_rls, xref, xread);
                }
                for (i = 0; i < natoms; i++)
                {
                    copy_rvec(xread[index[i]], x[i]);
                }

                for (v = 0; v < noutvec; v++)
                {
                    vec = outvec[v];
                    /* calculate (mass-weighted) projection */
                    inp = 0;
                    for (i = 0; i < natoms; i++)
                    {
                        inp += (eigvec[vec][i][0] * (x[i][0] - xav[i][0])
                                + eigvec[vec][i][1] * (x[i][1] - xav[i][1])
                                + eigvec[vec][i][2] * (x[i][2] - xav[i][2]))
                               * sqrtm[i];
                    }
                    inprod[v][nframes] = inp;
                }
                if (filterfile)
                {
                    for (i = 0; i < natoms; i++)
                    {
                        for (d = 0; d < DIM; d++)
                        {
                            /* misuse xread for output */
                            xread[index[i]][d] = xav[i][d];
                            for (v = 0; v < noutvec; v++)
                            {
                                xread[index[i]][d] +=
                                        inprod[v][nframes] * eigvec[outvec[v]][i][d] / sqrtm[i];
                            }
                        }
                    }
                    write_trx(out, natoms, index, atoms, 0, t, box, xread, nullptr, nullptr);
                }
                nframes++;
            }
            nfr++;
        } while (read_next_x(oenv, status, &t, xread, box));
        close_trx(status);
        sfree(x);
        if (filterfile)
        {
            close_trx(out);
        }
    }
    else
    {
        snew(xread, atoms->nr);
    }

    if (top)
    {
        gmx_rmpbc_done(gpbc);
    }


    if (projfile)
    {
        GMX_RELEASE_ASSERT(inprod != nullptr, "inprod must be non-NULL if projfile is non-NULL");
        snew(ylabel, noutvec);
        for (v = 0; v < noutvec; v++)
        {
            sprintf(str, "vec %d", eignr[outvec[v]] + 1);
            ylabel[v] = gmx_strdup(str);
        }
        sprintf(str, "projection on eigenvectors (%s)", projUnit);
        write_xvgr_graphs(projfile,
                          noutvec,
                          1,
                          str,
                          nullptr,
                          output_env_get_xvgr_tlabel(oenv),
                          ylabel,
                          nframes,
                          inprod[noutvec],
                          inprod,
                          nullptr,
                          output_env_get_time_factor(oenv),
                          FALSE,
                          bSplit,
                          oenv);
    }

    if (twodplotfile)
    {
        sprintf(str, "projection on eigenvector %d (%s)", eignr[outvec[0]] + 1, projUnit);
        sprintf(str2, "projection on eigenvector %d (%s)", eignr[outvec[noutvec - 1]] + 1, projUnit);
        xvgrout = xvgropen(twodplotfile, "2D projection of trajectory", str, str2, oenv);
        for (i = 0; i < nframes; i++)
        {
            if (bSplit && i > 0 && std::abs(inprod[noutvec][i]) < 1e-5)
            {
                fprintf(xvgrout, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
            }
            fprintf(xvgrout, "%10.5f %10.5f\n", inprod[0][i], inprod[noutvec - 1][i]);
        }
        xvgrclose(xvgrout);
    }

    if (threedplotfile)
    {
        t_atoms  atoms3D;
        rvec*    x;
        real*    b = nullptr;
        matrix   box;
        char *   resnm, *atnm;
        gmx_bool bPDB, b4D;
        FILE*    out;

        if (noutvec < 3)
        {
            gmx_fatal(FARGS, "You have selected less than 3 eigenvectors");
        }

        /* initialize */
        bPDB = fn2ftp(threedplotfile) == efPDB;
        clear_mat(box);
        box[XX][XX] = box[YY][YY] = box[ZZ][ZZ] = 1;

        b4D = bPDB && (noutvec >= 4);
        if (b4D)
        {
            fprintf(stderr,
                    "You have selected four or more eigenvectors:\n"
                    "fourth eigenvector will be plotted "
                    "in bfactor field of pdb file\n");
            sprintf(str,
                    "4D proj. of traj. on eigenv. %d, %d, %d and %d",
                    eignr[outvec[0]] + 1,
                    eignr[outvec[1]] + 1,
                    eignr[outvec[2]] + 1,
                    eignr[outvec[3]] + 1);
        }
        else
        {
            sprintf(str,
                    "3D proj. of traj. on eigenv. %d, %d and %d",
                    eignr[outvec[0]] + 1,
                    eignr[outvec[1]] + 1,
                    eignr[outvec[2]] + 1);
        }
        init_t_atoms(&atoms3D, nframes, FALSE);
        snew(x, nframes);
        snew(b, nframes);
        atnm  = gmx_strdup("C");
        resnm = gmx_strdup("PRJ");

        if (nframes > 10000)
        {
            fact = 10000.0 / nframes;
        }
        else
        {
            fact = 1.0;
        }

        for (i = 0; i < nframes; i++)
        {
            atoms3D.atomname[i]     = &atnm;
            atoms3D.atom[i].resind  = i;
            atoms3D.resinfo[i].name = &resnm;
            atoms3D.resinfo[i].nr   = static_cast<int>(std::ceil(i * fact));
            atoms3D.resinfo[i].ic   = ' ';
            x[i][XX]                = inprod[0][i];
            x[i][YY]                = inprod[1][i];
            x[i][ZZ]                = inprod[2][i];
            if (b4D)
            {
                b[i] = inprod[3][i];
            }
        }
        if ((b4D || bSplit) && bPDB)
        {
            GMX_RELEASE_ASSERT(inprod != nullptr,
                               "inprod must be non-NULL with 4D or split PDB output options");

            out = gmx_ffopen(threedplotfile, "w");
            fprintf(out, "HEADER    %s\n", str);
            if (b4D)
            {
                fprintf(out, "REMARK    %s\n", "fourth dimension plotted as B-factor");
            }
            j = 0;
            for (i = 0; i < atoms3D.nr; i++)
            {
                if (j > 0 && bSplit && std::abs(inprod[noutvec][i]) < 1e-5)
                {
                    fprintf(out, "TER\n");
                    j = 0;
                }
                gmx_fprintf_pdb_atomline(out,
                                         PdbRecordType::Atom,
                                         i + 1,
                                         "C",
                                         ' ',
                                         "PRJ",
                                         ' ',
                                         j + 1,
                                         ' ',
                                         10 * x[i][XX],
                                         10 * x[i][YY],
                                         10 * x[i][ZZ],
                                         1.0,
                                         10 * b[i],
                                         "");
                if (j > 0)
                {
                    fprintf(out, "CONECT%5d%5d\n", i, i + 1);
                }
                j++;
            }
            fprintf(out, "TER\n");
            gmx_ffclose(out);
        }
        else
        {
            write_sto_conf(threedplotfile, str, &atoms3D, x, nullptr, pbcType, box);
        }
        done_atom(&atoms3D);
    }

    if (extremefile)
    {
        snew(pmin, noutvec_extr);
        snew(pmax, noutvec_extr);
        if (extreme == 0)
        {
            GMX_RELEASE_ASSERT(inprod != nullptr, "inprod must be non-NULL");
            fprintf(stderr, "%11s %17s %17s\n", "eigenvector", "Minimum", "Maximum");
            fprintf(stderr, "%11s %10s %10s %10s %10s\n", "", "value", "frame", "value", "frame");
            imin = 0;
            imax = 0;
            for (v = 0; v < noutvec_extr; v++)
            {
                for (i = 0; i < nframes; i++)
                {
                    if (inprod[v][i] < inprod[v][imin])
                    {
                        imin = i;
                    }
                    if (inprod[v][i] > inprod[v][imax])
                    {
                        imax = i;
                    }
                }
                pmin[v] = inprod[v][imin];
                pmax[v] = inprod[v][imax];
                fprintf(stderr, "%7d     %10.6f %10d %10.6f %10d\n", eignr[outvec[v]] + 1, pmin[v], imin, pmax[v], imax);
            }
        }
        else
        {
            pmin[0] = -extreme;
            pmax[0] = extreme;
        }
        std::filesystem::path extremeFileToUse = extremefile;
        for (v = 0; v < noutvec_extr; v++)
        {
            /* make filename using format string */
            if (noutvec_extr != 1)
            {
                extremeFileToUse = gmx::concatenateBeforeExtension(
                        extremefile, std::to_string(eignr[outvec[v]] + 1));
            }
            fprintf(stderr,
                    "Writing %d frames along eigenvector %d to %s\n",
                    nextr,
                    outvec[v] + 1,
                    extremeFileToUse.string().c_str());
            out = open_trx(extremeFileToUse, "w");
            for (frame = 0; frame < nextr; frame++)
            {
                if ((extreme == 0) && (nextr <= 3))
                {
                    for (i = 0; i < natoms; i++)
                    {
                        atoms->resinfo[atoms->atom[index[i]].resind].chainid = 'A' + frame;
                    }
                }
                for (i = 0; i < natoms; i++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        xread[index[i]][d] =
                                (xav[i][d]
                                 + (pmin[v] * (nextr - frame - 1) + pmax[v] * frame) / (nextr - 1)
                                           * eigvec[outvec[v]][i][d] / sqrtm[i]);
                    }
                }
                write_trx(out, natoms, index, atoms, 0, frame, topbox, xread, nullptr, nullptr);
            }
            close_trx(out);
        }
        sfree(pmin);
        sfree(pmax);
    }
    fprintf(stderr, "\n");
}

static void components(const char*             outfile,
                       int                     natoms,
                       int*                    eignr,
                       rvec**                  eigvec,
                       int                     noutvec,
                       const int*              outvec,
                       const gmx_output_env_t* oenv)
{
    int          g, s, v, i;
    real *       x, ***y;
    char         str[STRLEN];
    const char** ylabel;

    fprintf(stderr, "Writing eigenvector components to %s\n", outfile);

    snew(ylabel, noutvec);
    snew(y, noutvec);
    snew(x, natoms);
    for (i = 0; i < natoms; i++)
    {
        x[i] = i + 1;
    }
    for (g = 0; g < noutvec; g++)
    {
        v = outvec[g];
        sprintf(str, "vec %d", eignr[v] + 1);
        ylabel[g] = gmx_strdup(str);
        snew(y[g], 4);
        for (s = 0; s < 4; s++)
        {
            snew(y[g][s], natoms);
        }
        for (i = 0; i < natoms; i++)
        {
            y[g][0][i] = norm(eigvec[v][i]);
            for (s = 0; s < 3; s++)
            {
                y[g][s + 1][i] = eigvec[v][i][s];
            }
        }
    }
    write_xvgr_graphs(outfile,
                      noutvec,
                      4,
                      "Eigenvector components",
                      "black: total, red: x, green: y, blue: z",
                      "Atom number",
                      ylabel,
                      natoms,
                      x,
                      nullptr,
                      y,
                      1,
                      FALSE,
                      FALSE,
                      oenv);
    fprintf(stderr, "\n");
}

static void rmsf(const char*             outfile,
                 int                     natoms,
                 const real*             sqrtm,
                 int*                    eignr,
                 rvec**                  eigvec,
                 int                     noutvec,
                 const int*              outvec,
                 real*                   eigval,
                 int                     neig,
                 const gmx_output_env_t* oenv)
{
    int          g, v, i;
    real *       x, **y;
    char         str[STRLEN];
    const char** ylabel;

    for (i = 0; i < neig; i++)
    {
        if (eigval[i] < 0)
        {
            eigval[i] = 0;
        }
    }

    fprintf(stderr, "Writing rmsf to %s\n", outfile);

    snew(ylabel, noutvec);
    snew(y, noutvec);
    snew(x, natoms);
    for (i = 0; i < natoms; i++)
    {
        x[i] = i + 1;
    }
    for (g = 0; g < noutvec; g++)
    {
        v = outvec[g];
        if (eignr[v] >= neig)
        {
            gmx_fatal(FARGS,
                      "Selected vector %d is larger than the number of eigenvalues (%d)",
                      eignr[v] + 1,
                      neig);
        }
        sprintf(str, "vec %d", eignr[v] + 1);
        ylabel[g] = gmx_strdup(str);
        snew(y[g], natoms);
        for (i = 0; i < natoms; i++)
        {
            y[g][i] = std::sqrt(eigval[eignr[v]] * norm2(eigvec[v][i])) / sqrtm[i];
        }
    }
    write_xvgr_graphs(
            outfile, noutvec, 1, "RMS fluctuation (nm) ", nullptr, "Atom number", ylabel, natoms, x, y, nullptr, 1, TRUE, FALSE, oenv);
    fprintf(stderr, "\n");
}

int gmx_anaeig(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] analyzes eigenvectors. The eigenvectors can be of a",
        "covariance matrix ([gmx-covar]) or of a Normal Modes analysis",
        "([gmx-nmeig]).[PAR]",

        "When a trajectory is projected on eigenvectors, all structures are",
        "fitted to the structure in the eigenvector file, if present, otherwise",
        "to the structure in the structure file. When no run input file is",
        "supplied, periodicity will not be taken into account. Most analyses",
        "are performed on eigenvectors [TT]-first[tt] to [TT]-last[tt], but when",
        "[TT]-first[tt] is set to -1 you will be prompted for a selection.[PAR]",

        "[TT]-comp[tt]: plot the vector components per atom of eigenvectors",
        "[TT]-first[tt] to [TT]-last[tt].[PAR]",

        "[TT]-rmsf[tt]: plot the RMS fluctuation per atom of eigenvectors",
        "[TT]-first[tt] to [TT]-last[tt] (requires [TT]-eig[tt]).[PAR]",

        "[TT]-proj[tt]: calculate projections of a trajectory on eigenvectors",
        "[TT]-first[tt] to [TT]-last[tt].",
        "The projections of a trajectory on the eigenvectors of its",
        "covariance matrix are called principal components (pc's).",
        "It is often useful to check the cosine content of the pc's,",
        "since the pc's of random diffusion are cosines with the number",
        "of periods equal to half the pc index.",
        "The cosine content of the pc's can be calculated with the program",
        "[gmx-analyze].[PAR]",

        "[TT]-2d[tt]: calculate a 2d projection of a trajectory on eigenvectors",
        "[TT]-first[tt] and [TT]-last[tt].[PAR]",

        "[TT]-3d[tt]: calculate a 3d projection of a trajectory on the first",
        "three selected eigenvectors.[PAR]",

        "[TT]-filt[tt]: filter the trajectory to show only the motion along",
        "eigenvectors [TT]-first[tt] to [TT]-last[tt].[PAR]",

        "[TT]-extr[tt]: calculate the two extreme projections along a trajectory",
        "on the average structure and interpolate [TT]-nframes[tt] frames",
        "between them, or set your own extremes with [TT]-max[tt]. The",
        "eigenvector [TT]-first[tt] will be written unless [TT]-first[tt] and",
        "[TT]-last[tt] have been set explicitly, in which case all eigenvectors",
        "will be written to separate files. Chain identifiers will be added",
        "when writing a [REF].pdb[ref] file with two or three structures (you",
        "can use [TT]rasmol -nmrpdb[tt] to view such a [REF].pdb[ref] file).[PAR]",

        "Overlap calculations between covariance analysis",
        "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",
        "",
        "[BB]Note:[bb] the analysis should use the same fitting structure",
        "",
        "[TT]-over[tt]: calculate the subspace overlap of the eigenvectors in",
        "file [TT]-v2[tt] with eigenvectors [TT]-first[tt] to [TT]-last[tt]",
        "in file [TT]-v[tt].[PAR]",

        "[TT]-inpr[tt]: calculate a matrix of inner-products between",
        "eigenvectors in files [TT]-v[tt] and [TT]-v2[tt]. All eigenvectors",
        "of both files will be used unless [TT]-first[tt] and [TT]-last[tt]",
        "have been set explicitly.[PAR]",

        "When [TT]-v[tt] and [TT]-v2[tt] are given, a single number for the",
        "overlap between the covariance matrices is generated. Note that the",
        "eigenvalues are by default read from the timestamp field in the",
        "eigenvector input files, but when [TT]-eig[tt], or [TT]-eig2[tt] are",
        "given, the corresponding eigenvalues are used instead. The formulas are::",
        "",
        "         difference = sqrt(tr((sqrt(M1) - sqrt(M2))^2))",
        " normalized overlap = 1 - difference/sqrt(tr(M1) + tr(M2))",
        "      shape overlap = 1 - sqrt(tr((sqrt(M1/tr(M1)) - sqrt(M2/tr(M2)))^2))",
        "",
        "where M1 and M2 are the two covariance matrices and tr is the trace",
        "of a matrix. The numbers are proportional to the overlap of the square",
        "root of the fluctuations. The normalized overlap is the most useful",
        "number, it is 1 for identical matrices and 0 when the sampled",
        "subspaces are orthogonal.[PAR]",
        "When the [TT]-entropy[tt] flag is given an entropy estimate will be",
        "computed based on the Quasiharmonic approach and based on",
        "Schlitter's formula."
    };
    static int      first = 1, last = -1, skip = 1, nextr = 2, nskip = 6;
    static real     max = 0.0, temp = 298.15;
    static gmx_bool bSplit = FALSE, bEntropy = FALSE;
    t_pargs         pa[] = {
        { "-first", FALSE, etINT, { &first }, "First eigenvector for analysis (-1 is select)" },
        { "-last", FALSE, etINT, { &last }, "Last eigenvector for analysis (-1 is till the last)" },
        { "-skip", FALSE, etINT, { &skip }, "Only analyse every nr-th frame" },
        { "-max",
          FALSE,
          etREAL,
          { &max },
          "Maximum for projection of the eigenvector on the average structure, "
          "max=0 gives the extremes" },
        { "-nframes", FALSE, etINT, { &nextr }, "Number of frames for the extremes output" },
        { "-split", FALSE, etBOOL, { &bSplit }, "Split eigenvector projections where time is zero" },
        { "-entropy",
          FALSE,
          etBOOL,
          { &bEntropy },
          "Compute entropy according to the Quasiharmonic formula or Schlitter's method." },
        { "-temp", FALSE, etREAL, { &temp }, "Temperature for entropy calculations" },
        { "-nevskip",
          FALSE,
          etINT,
          { &nskip },
          "Number of eigenvalues to skip when computing the entropy due to the quasi harmonic "
          "approximation. When you do a rotational and/or translational fit prior to the "
          "covariance analysis, you get 3 or 6 eigenvalues that are very close to zero, and which "
          "should not be taken into account when computing the entropy." }
    };
#define NPA asize(pa)

    t_topology        top;
    PbcType           pbcType = PbcType::Unset;
    const t_atoms*    atoms   = nullptr;
    rvec *            xtop, *xref1, *xref2, *xrefp = nullptr;
    gmx_bool          bDMR1, bDMA1, bDMR2, bDMA2;
    int               nvec1, nvec2, *eignr1 = nullptr, *eignr2 = nullptr;
    rvec *            xav1, *xav2, **eigvec1 = nullptr, **eigvec2 = nullptr;
    matrix            topbox;
    real              totmass, *sqrtm, *w_rls, t;
    int               natoms;
    char*             grpname;
    const char*       indexfile;
    int               i, j, d;
    int               nout, *iout, noutvec, *outvec, nfit;
    int *             index = nullptr, *ifit = nullptr;
    const char *      VecFile, *Vec2File, *topfile;
    const char *      EigFile, *Eig2File;
    const char *      CompFile, *RmsfFile, *ProjOnVecFile;
    const char *      TwoDPlotFile, *ThreeDPlotFile;
    const char *      FilterFile, *ExtremeFile;
    const char *      OverlapFile, *InpMatFile;
    gmx_bool          bFit1, bFit2, bM, bIndex, bTPS, bTop, bVec2, bProj;
    gmx_bool          bFirstToLast, bFirstLastSet, bTraj, bCompare, bPDB3D;
    real *            eigval1 = nullptr, *eigval2 = nullptr;
    int               neig1, neig2;
    double**          xvgdata;
    gmx_output_env_t* oenv;
    gmx_rmpbc_t       gpbc;

    t_filenm fnm[] = {
        { efTRN, "-v", "eigenvec", ffREAD },      { efTRN, "-v2", "eigenvec2", ffOPTRD },
        { efTRX, "-f", nullptr, ffOPTRD },        { efTPS, nullptr, nullptr, ffOPTRD },
        { efNDX, nullptr, nullptr, ffOPTRD },     { efXVG, "-eig", "eigenval", ffOPTRD },
        { efXVG, "-eig2", "eigenval2", ffOPTRD }, { efXVG, "-comp", "eigcomp", ffOPTWR },
        { efXVG, "-rmsf", "eigrmsf", ffOPTWR },   { efXVG, "-proj", "proj", ffOPTWR },
        { efXVG, "-2d", "2dproj", ffOPTWR },      { efSTO, "-3d", "3dproj.pdb", ffOPTWR },
        { efTRX, "-filt", "filtered", ffOPTWR },  { efTRX, "-extr", "extreme.pdb", ffOPTWR },
        { efXVG, "-over", "overlap", ffOPTWR },   { efXPM, "-inpr", "inprod", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW, NFILE, fnm, NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    indexfile = ftp2fn_null(efNDX, NFILE, fnm);

    VecFile        = opt2fn("-v", NFILE, fnm);
    Vec2File       = opt2fn_null("-v2", NFILE, fnm);
    topfile        = ftp2fn(efTPS, NFILE, fnm);
    EigFile        = opt2fn_null("-eig", NFILE, fnm);
    Eig2File       = opt2fn_null("-eig2", NFILE, fnm);
    CompFile       = opt2fn_null("-comp", NFILE, fnm);
    RmsfFile       = opt2fn_null("-rmsf", NFILE, fnm);
    ProjOnVecFile  = opt2fn_null("-proj", NFILE, fnm);
    TwoDPlotFile   = opt2fn_null("-2d", NFILE, fnm);
    ThreeDPlotFile = opt2fn_null("-3d", NFILE, fnm);
    FilterFile     = opt2fn_null("-filt", NFILE, fnm);
    ExtremeFile    = opt2fn_null("-extr", NFILE, fnm);
    OverlapFile    = opt2fn_null("-over", NFILE, fnm);
    InpMatFile     = ftp2fn_null(efXPM, NFILE, fnm);

    bProj = (ProjOnVecFile != nullptr) || (TwoDPlotFile != nullptr) || (ThreeDPlotFile != nullptr)
            || (FilterFile != nullptr) || (ExtremeFile != nullptr);
    bFirstLastSet = opt2parg_bSet("-first", NPA, pa) && opt2parg_bSet("-last", NPA, pa);
    bFirstToLast  = (CompFile != nullptr) || (RmsfFile != nullptr) || (ProjOnVecFile != nullptr)
                   || (FilterFile != nullptr) || (OverlapFile != nullptr)
                   || (((ExtremeFile != nullptr) || (InpMatFile != nullptr)) && bFirstLastSet);
    bVec2 = (Vec2File != nullptr) || (OverlapFile != nullptr) || (InpMatFile != nullptr);
    bM    = (RmsfFile != nullptr) || bProj;
    bTraj = (ProjOnVecFile != nullptr) || (FilterFile != nullptr)
            || ((ExtremeFile != nullptr) && (max == 0)) || (TwoDPlotFile != nullptr)
            || (ThreeDPlotFile != nullptr);
    bIndex = bM || bProj;
    bTPS   = ftp2bSet(efTPS, NFILE, fnm) || bM || bTraj || (FilterFile != nullptr)
           || (bIndex && (indexfile != nullptr));
    bCompare = (Vec2File != nullptr) || (Eig2File != nullptr);
    bPDB3D   = fn2ftp(ThreeDPlotFile) == efPDB;

    read_eigenvectors(
            VecFile, &natoms, &bFit1, &xref1, &bDMR1, &xav1, &bDMA1, &nvec1, &eignr1, &eigvec1, &eigval1);
    neig1 = std::min(nvec1, DIM * natoms);
    if (nvec1 != DIM * natoms)
    {
        fprintf(stderr,
                "Warning: number of eigenvectors %d does not match three times\n"
                "the number of atoms %d in %s. Using %d eigenvectors.\n\n",
                nvec1,
                natoms,
                VecFile,
                neig1);
    }

    /* Overwrite eigenvalues from separate files if the user provides them */
    if (EigFile != nullptr)
    {
        int neig_tmp = read_xvg(EigFile, &xvgdata, &i);
        if (neig_tmp != neig1)
        {
            fprintf(stderr,
                    "Warning: number of eigenvalues in xvg file (%d) does not mtch trr file (%d)\n",
                    neig1,
                    natoms);
        }
        neig1 = neig_tmp;
        srenew(eigval1, neig1);
        for (j = 0; j < neig1; j++)
        {
            real tmp   = eigval1[j];
            eigval1[j] = xvgdata[1][j];
            if (debug && (eigval1[j] != tmp))
            {
                fprintf(debug, "Replacing eigenvalue %d. From trr: %10g, from xvg: %10g\n", j, tmp, eigval1[j]);
            }
        }
        for (j = 0; j < i; j++)
        {
            sfree(xvgdata[j]);
        }
        sfree(xvgdata);
        fprintf(stderr, "Read %d eigenvalues from %s\n", neig1, EigFile);
    }

    if (bEntropy)
    {
        if (bDMA1)
        {
            gmx_fatal(FARGS,
                      "Can not calculate entropies from mass-weighted eigenvalues, redo the "
                      "analysis without mass-weighting");
        }
        printf("The Entropy due to the Schlitter formula is %g J/mol K\n",
               calcSchlitterEntropy(gmx::arrayRefFromArray(eigval1, neig1), temp, FALSE));
        printf("The Entropy due to the Quasiharmonic analysis is %g J/mol K\n",
               calcQuasiHarmonicEntropy(gmx::arrayRefFromArray(eigval1, neig1), temp, FALSE, 1.0));
    }

    if (bVec2)
    {
        int natoms2;
        if (Vec2File)
        {
            read_eigenvectors(
                    Vec2File, &natoms2, &bFit2, &xref2, &bDMR2, &xav2, &bDMA2, &nvec2, &eignr2, &eigvec2, &eigval2);
        }
        else
        {
            gmx_fatal(FARGS, "Need a second eigenvector file to do this analysis.");
        }

        neig2 = std::min(nvec2, DIM * natoms2);
        if (neig2 != neig1)
        {
            gmx_fatal(FARGS, "Dimensions in the eigenvector files don't match");
        }
    }
    else
    {
        nvec2 = 0;
        neig2 = 0;
    }

    if (Eig2File != nullptr)
    {
        neig2 = read_xvg(Eig2File, &xvgdata, &i);
        srenew(eigval2, neig2);
        for (j = 0; j < neig2; j++)
        {
            eigval2[j] = xvgdata[1][j];
        }
        for (j = 0; j < i; j++)
        {
            sfree(xvgdata[j]);
        }
        sfree(xvgdata);
        fprintf(stderr, "Read %d eigenvalues from %s\n", neig2, Eig2File);
    }


    if ((!bFit1 || xref1) && !bDMR1 && !bDMA1)
    {
        bM = FALSE;
    }
    if ((xref1 == nullptr) && (bM || bTraj))
    {
        bTPS = TRUE;
    }

    xtop  = nullptr;
    nfit  = 0;
    ifit  = nullptr;
    w_rls = nullptr;

    if (!bTPS)
    {
        bTop = FALSE;
    }
    else
    {
        bTop = read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, topbox, bM);
        atoms = &top.atoms;
        gpbc  = gmx_rmpbc_init(&top.idef, pbcType, atoms->nr);
        gmx_rmpbc_apply(gpbc, atoms->nr, topbox, xtop);
        /* Fitting is only required for the projection */
        if (bProj && bFit1)
        {
            if (xref1 == nullptr)
            {
                printf("\nNote: the structure in %s should be the same\n"
                       "      as the one used for the fit in gmx covar\n",
                       topfile);
            }
            printf("\nSelect the index group that was used for the least squares fit in gmx "
                   "covar\n");
            get_index(atoms, indexfile, 1, &nfit, &ifit, &grpname);

            snew(w_rls, atoms->nr);
            for (i = 0; (i < nfit); i++)
            {
                if (bDMR1)
                {
                    w_rls[ifit[i]] = atoms->atom[ifit[i]].m;
                }
                else
                {
                    w_rls[ifit[i]] = 1.0;
                }
            }

            snew(xrefp, atoms->nr);
            if (xref1 != nullptr)
            {
                /* Safety check between selected fit-group and reference structure read from the eigenvector file */
                if (natoms != nfit)
                {
                    gmx_fatal(FARGS,
                              "you selected a group with %d elements instead of %d, your selection "
                              "does not fit the reference structure in the eigenvector file.",
                              nfit,
                              natoms);
                }
                for (i = 0; (i < nfit); i++)
                {
                    copy_rvec(xref1[i], xrefp[ifit[i]]);
                }
            }
            else
            {
                /* The top coordinates are the fitting reference */
                for (i = 0; (i < nfit); i++)
                {
                    copy_rvec(xtop[ifit[i]], xrefp[ifit[i]]);
                }
                reset_x(nfit, ifit, atoms->nr, nullptr, xrefp, w_rls);
            }
        }
        gmx_rmpbc_done(gpbc);
    }

    if (bIndex)
    {
        printf("\nSelect an index group of %d elements that corresponds to the eigenvectors\n", natoms);
        get_index(atoms, indexfile, 1, &i, &index, &grpname);
        if (i != natoms)
        {
            gmx_fatal(FARGS, "you selected a group with %d elements instead of %d", i, natoms);
        }
        printf("\n");
    }

    snew(sqrtm, natoms);
    std::string projUnit;
    if (bM && bDMA1)
    {
        projUnit = "u\\S1/2\\Nnm";
        for (i = 0; (i < natoms); i++)
        {
            sqrtm[i] = std::sqrt(atoms->atom[index[i]].m);
        }
    }
    else
    {
        projUnit = "nm";
        for (i = 0; (i < natoms); i++)
        {
            sqrtm[i] = 1.0;
        }
    }

    if (bVec2)
    {
        t       = 0;
        totmass = 0;
        for (i = 0; (i < natoms); i++)
        {
            for (d = 0; (d < DIM); d++)
            {
                t += gmx::square((xav1[i][d] - xav2[i][d]) * sqrtm[i]);
                totmass += gmx::square(sqrtm[i]);
            }
        }
        fprintf(stdout,
                "RMSD (without fit) between the two average structures:"
                " %.3f (nm)\n\n",
                std::sqrt(t / totmass));
    }

    if (last == -1)
    {
        last = natoms * DIM;
    }
    if (first > -1)
    {
        if (bFirstToLast)
        {
            /* make an index from first to last */
            nout = last - first + 1;
            snew(iout, nout);
            for (i = 0; i < nout; i++)
            {
                iout[i] = first - 1 + i;
            }
        }
        else if (ThreeDPlotFile)
        {
            /* make an index of first+(0,1,2) and last */
            nout = bPDB3D ? 4 : 3;
            nout = std::min(last - first + 1, nout);
            snew(iout, nout);
            iout[0] = first - 1;
            iout[1] = first;
            if (nout > 3)
            {
                iout[2] = first + 1;
            }
            iout[nout - 1] = last - 1;
        }
        else
        {
            /* make an index of first and last */
            nout = 2;
            snew(iout, nout);
            iout[0] = first - 1;
            iout[1] = last - 1;
        }
    }
    else
    {
        printf("Select eigenvectors for output, end your selection with 0\n");
        nout = -1;
        iout = nullptr;

        do
        {
            nout++;
            srenew(iout, nout + 1);
            if (1 != scanf("%d", &iout[nout]))
            {
                gmx_fatal(FARGS, "Error reading user input");
            }
            iout[nout]--;
        } while (iout[nout] >= 0);

        printf("\n");
    }
    /* make an index of the eigenvectors which are present */
    snew(outvec, nout);
    noutvec = 0;
    for (i = 0; i < nout; i++)
    {
        j = 0;
        while ((j < nvec1) && (eignr1[j] != iout[i]))
        {
            j++;
        }
        if ((j < nvec1) && (eignr1[j] == iout[i]))
        {
            outvec[noutvec] = j;
            noutvec++;
        }
    }
    fprintf(stderr, "%d eigenvectors selected for output", noutvec);
    if (noutvec <= 100)
    {
        fprintf(stderr, ":");
        for (j = 0; j < noutvec; j++)
        {
            fprintf(stderr, " %d", eignr1[outvec[j]] + 1);
        }
    }
    fprintf(stderr, "\n");

    if (CompFile)
    {
        components(CompFile, natoms, eignr1, eigvec1, noutvec, outvec, oenv);
    }

    if (RmsfFile)
    {
        rmsf(RmsfFile, natoms, sqrtm, eignr1, eigvec1, noutvec, outvec, eigval1, neig1, oenv);
    }

    if (bProj)
    {
        project(bTraj ? opt2fn("-f", NFILE, fnm) : nullptr,
                bTop ? &top : nullptr,
                pbcType,
                topbox,
                ProjOnVecFile,
                TwoDPlotFile,
                ThreeDPlotFile,
                FilterFile,
                projUnit.c_str(),
                skip,
                ExtremeFile,
                bFirstLastSet,
                max,
                nextr,
                atoms,
                natoms,
                index,
                bFit1,
                xrefp,
                nfit,
                ifit,
                w_rls,
                sqrtm,
                xav1,
                eignr1,
                eigvec1,
                noutvec,
                outvec,
                bSplit,
                oenv);
    }

    if (OverlapFile)
    {
        overlap(OverlapFile, natoms, eigvec1, nvec2, eignr2, eigvec2, noutvec, outvec, oenv);
    }

    if (InpMatFile)
    {
        inprod_matrix(
                InpMatFile, natoms, nvec1, eignr1, eigvec1, nvec2, eignr2, eigvec2, bFirstLastSet, noutvec, outvec);
    }

    if (bCompare)
    {
        compare(natoms, nvec1, eigvec1, nvec2, eigvec2, eigval1, neig1, eigval2, neig2);
    }


    if (!CompFile && !bProj && !OverlapFile && !InpMatFile && !bCompare && !bEntropy)
    {
        fprintf(stderr,
                "\nIf you want some output,"
                " set one (or two or ...) of the output file options\n");
    }


    view_all(oenv, NFILE, fnm);

    return 0;
}
