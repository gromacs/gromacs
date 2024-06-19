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

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/rgb.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

typedef struct
{
    int    nx, ny;
    real   dx, dy;
    real** data;
} t_shiftdata;

static real interpolate(real phi, real psi, t_shiftdata* sd)
{
    int  iphi, ipsi, iphi1, ipsi1;
    real fphi, fpsi, wx0, wx1, wy0, wy1;

    /*phi  += M_PI;
       if (phi > 2*M_PI) phi -= 2*M_PI;
       psi  += M_PI;
       if (psi > 2*M_PI) psi -= 2*M_PI;
     */
    while (phi < 0)
    {
        phi += 2 * M_PI;
    }
    while (psi < 0)
    {
        psi += 2 * M_PI;
    }
    phi = 2 * M_PI - phi;

    fphi = phi * sd->dx;
    fpsi = psi * sd->dy;

    iphi = static_cast<int>(fphi);
    ipsi = static_cast<int>(fpsi);
    fphi -= iphi; /* Fraction (offset from gridpoint) */
    fpsi -= ipsi;

    wx0   = 1.0 - fphi;
    wx1   = fphi;
    wy0   = 1.0 - fpsi;
    wy1   = fpsi;
    iphi  = iphi % sd->nx;
    ipsi  = ipsi % sd->ny;
    iphi1 = (iphi + 1) % sd->nx;
    ipsi1 = (ipsi + 1) % sd->ny;

    return (sd->data[iphi][ipsi] * wx0 * wy0 + sd->data[iphi1][ipsi] * wx1 * wy0
            + sd->data[iphi][ipsi1] * wx0 * wy1 + sd->data[iphi1][ipsi1] * wx1 * wy1);
}

static void dump_sd(const char* fn, t_shiftdata* sd)
{
    FILE* fp;
    int   i, j;
    char  buf[256];
    int   nnx, nny, nfac = 4, nlevels = 20;
    real  phi, psi, *x_phi, *y_psi, **newdata;
    real  lo, hi;
    t_rgb rlo = { 1, 0, 0 }, rhi = { 0, 0, 1 };

    nnx = sd->nx * nfac + 1;
    nny = sd->ny * nfac + 1;
    snew(x_phi, nnx);
    snew(y_psi, nny);
    snew(newdata, nnx);
    lo = 100000;
    hi = -100000;
    for (i = 0; (i < nnx); i++)
    {
        snew(newdata[i], nny);
        phi      = i * 2 * M_PI / (nnx - 1);
        x_phi[i] = phi * gmx::c_rad2Deg - 180;
        for (j = 0; (j < nny); j++)
        {
            psi = j * 2 * M_PI / (nny - 1);
            if (i == 0)
            {
                y_psi[j] = psi * gmx::c_rad2Deg - 180;
            }
            /*if (((i % nfac) == 0) && ((j % nfac) == 0))
               newdata[i][j] = sd->data[i/nfac][j/nfac];
               else*/
            newdata[i][j] = interpolate(phi, psi, sd);
            lo            = std::min(lo, newdata[i][j]);
            hi            = std::max(hi, newdata[i][j]);
        }
    }
    sprintf(buf, "%s.xpm", fn);
    fp = gmx_ffopen(buf, "w");
    write_xpm(fp, 0, fn, fn, "Phi", "Psi", nnx, nny, x_phi, y_psi, newdata, lo, hi, rlo, rhi, &nlevels);
    for (i = 0; (i < nnx); i++)
    {
        sfree(newdata[i]);
    }
    sfree(newdata);
    sfree(x_phi);
    sfree(y_psi);
}

static t_shiftdata* read_shifts(const char* fn)
{
    double       xx;
    int          i, j, nx, ny;
    t_shiftdata* sd;

    snew(sd, 1);
    gmx::FilePtr fp = gmx::openLibraryFile(fn);
    if (2 != fscanf(fp.get(), "%d%d", &nx, &ny))
    {
        gmx_fatal(FARGS, "Error reading from file %s", fn);
    }
    GMX_ASSERT(nx > 0, "");
    sd->nx = nx;
    sd->ny = ny;
    sd->dx = nx / (2 * M_PI);
    sd->dy = ny / (2 * M_PI);
    snew(sd->data, nx + 1);
    for (i = 0; (i <= nx); i++)
    {
        snew(sd->data[i], ny + 1);
        for (j = 0; (j < ny); j++)
        {
            if (i == nx)
            {
                sd->data[i][j] = sd->data[0][j];
            }
            else
            {
                if (1 != fscanf(fp.get(), "%lf", &xx))
                {
                    gmx_fatal(FARGS, "Error reading from file %s", fn);
                }
                sd->data[i][j] = xx;
            }
        }
        sd->data[i][j] = sd->data[i][0];
    }

    if (bDebugMode())
    {
        dump_sd(fn, sd);
    }

    return sd;
}


static void done_shifts(t_shiftdata* sd)
{
    int i;

    for (i = 0; (i <= sd->nx); i++)
    {
        sfree(sd->data[i]);
    }
    sfree(sd->data);
    sfree(sd);
}

void do_pp2shifts(FILE* fp, int nf, gmx::ArrayRef<const t_dlist> dlist, real** dih)
{
    /* Read the shift files */
    t_shiftdata* ca_sd = read_shifts("ca-shift.dat");
    t_shiftdata* cb_sd = read_shifts("cb-shift.dat");
    t_shiftdata* ha_sd = read_shifts("ha-shift.dat");
    t_shiftdata* co_sd = read_shifts("co-shift.dat");

    fprintf(fp, "\n *** Chemical shifts from the chemical shift index ***\n");
    please_cite(fp, "Wishart98a");
    fprintf(fp,
            "%12s  %10s  %10s  %10s  %10s\n",
            "Residue",
            "delta Ca",
            "delta Ha",
            "delta CO",
            "delta Cb");
    for (const auto& dihedral : dlist)
    {
        if ((has_dihedral(edPhi, dihedral)) && (has_dihedral(edPsi, dihedral)))
        {
            int  Phi = dihedral.j0[edPhi];
            int  Psi = dihedral.j0[edPsi];
            real ca = 0, cb = 0, co = 0, ha = 0;
            for (int j = 0; (j < nf); j++)
            {
                real phi = dih[Phi][j];
                real psi = dih[Psi][j];

                ca += interpolate(phi, psi, ca_sd);
                cb += interpolate(phi, psi, cb_sd);
                co += interpolate(phi, psi, co_sd);
                ha += interpolate(phi, psi, ha_sd);
            }
            fprintf(fp, "%12s  %10g  %10g  %10g  %10g\n", dihedral.name, ca / nf, ha / nf, co / nf, cb / nf);
        }
    }
    fprintf(fp, "\n");

    /* Free memory */
    done_shifts(ca_sd);
    done_shifts(cb_sd);
    done_shifts(co_sd);
    done_shifts(ha_sd);
}
