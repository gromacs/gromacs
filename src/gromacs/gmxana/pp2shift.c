/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include <math.h>
#include <stdlib.h>

#include "gromacs/fileio/matio.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    int    nx, ny;
    real   dx, dy;
    real **data;
} t_shiftdata;

static real interpolate(real phi, real psi, t_shiftdata *sd)
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
        phi += 2*M_PI;
    }
    while (psi < 0)
    {
        psi += 2*M_PI;
    }
    phi    = 2*M_PI-phi;

    fphi  = phi*sd->dx;
    fpsi  = psi*sd->dy;

    iphi  = (int)fphi;
    ipsi  = (int)fpsi;
    fphi -= iphi; /* Fraction (offset from gridpoint) */
    fpsi -= ipsi;

    wx0   = 1.0-fphi;
    wx1   = fphi;
    wy0   = 1.0-fpsi;
    wy1   = fpsi;
    iphi  = iphi % sd->nx;
    ipsi  = ipsi % sd->ny;
    iphi1 = (iphi+1) % sd->nx;
    ipsi1 = (ipsi+1) % sd->ny;

    return (sd->data[iphi]  [ipsi]  * wx0*wy0 +
            sd->data[iphi1] [ipsi]  * wx1*wy0 +
            sd->data[iphi]  [ipsi1] * wx0*wy1 +
            sd->data[iphi1] [ipsi1] * wx1*wy1);
}

static void dump_sd(const char *fn, t_shiftdata *sd)
{
    FILE  *fp;
    int    i, j;
    char   buf[256];
    int    nnx, nny, nfac = 4, nlevels = 20;
    real   phi, psi, *x_phi, *y_psi, **newdata;
    real   lo, hi;
    t_rgb  rlo = { 1, 0, 0 }, rhi = { 0, 0, 1 };

    nnx = sd->nx*nfac+1;
    nny = sd->ny*nfac+1;
    snew(x_phi, nnx);
    snew(y_psi, nny);
    snew(newdata, nnx);
    lo =  100000;
    hi = -100000;
    for (i = 0; (i < nnx); i++)
    {
        snew(newdata[i], nny);
        phi      = i*2*M_PI/(nnx-1);
        x_phi[i] = phi*RAD2DEG-180;
        for (j = 0; (j < nny); j++)
        {
            psi = j*2*M_PI/(nny-1);
            if (i == 0)
            {
                y_psi[j] = psi*RAD2DEG-180;
            }
            /*if (((i % nfac) == 0) && ((j % nfac) == 0))
               newdata[i][j] = sd->data[i/nfac][j/nfac];
               else*/
            newdata[i][j] = interpolate(phi, psi, sd);
            lo            = min(lo, newdata[i][j]);
            hi            = max(hi, newdata[i][j]);
        }
    }
    sprintf(buf, "%s.xpm", fn);
    fp = gmx_ffopen(buf, "w");
    write_xpm(fp, 0, fn, fn, "Phi", "Psi", nnx, nny,
              x_phi, y_psi, newdata, lo, hi, rlo, rhi, &nlevels);
    for (i = 0; (i < nnx); i++)
    {
        sfree(newdata[i]);
    }
    sfree(newdata);
    sfree(x_phi);
    sfree(y_psi);
}

static t_shiftdata *read_shifts(const char *fn)
{
    FILE        *fp;
    double       xx;
    int          i, j, nx, ny;
    t_shiftdata *sd;

    snew(sd, 1);
    fp = libopen(fn);
    if (2 != fscanf(fp, "%d%d", &nx, &ny))
    {
        gmx_fatal(FARGS, "Error reading from file %s", fn);
    }

    sd->nx = nx;
    sd->ny = ny;
    sd->dx = nx/(2*M_PI);
    sd->dy = ny/(2*M_PI);
    snew(sd->data, nx+1);
    for (i = 0; (i <= nx); i++)
    {
        snew(sd->data[i], ny+1);
        for (j = 0; (j < ny); j++)
        {
            if (i == nx)
            {
                sd->data[i][j] = sd->data[0][j];
            }
            else
            {
                if (1 != fscanf(fp, "%lf", &xx))
                {
                    gmx_fatal(FARGS, "Error reading from file %s", fn);
                }
                sd->data[i][j] = xx;
            }
        }
        sd->data[i][j] = sd->data[i][0];
    }
    gmx_ffclose(fp);

    if (bDebugMode())
    {
        dump_sd(fn, sd);
    }

    return sd;
}


static void done_shifts(t_shiftdata *sd)
{
    int i;

    for (i = 0; (i <= sd->nx); i++)
    {
        sfree(sd->data[i]);
    }
    sfree(sd->data);
    sfree(sd);
}

void do_pp2shifts(FILE *fp, int nf, int nlist, t_dlist dlist[], real **dih)
{
    t_shiftdata *ca_sd, *co_sd, *ha_sd, *cb_sd;
    int          i, j, Phi, Psi;
    real         phi, psi;
    real         ca, co, ha, cb;

    /* Read the shift files */
    ca_sd = read_shifts("ca-shift.dat");
    cb_sd = read_shifts("cb-shift.dat");
    ha_sd = read_shifts("ha-shift.dat");
    co_sd = read_shifts("co-shift.dat");

    fprintf(fp, "\n *** Chemical shifts from the chemical shift index ***\n");
    please_cite(fp, "Wishart98a");
    fprintf(fp, "%12s  %10s  %10s  %10s  %10s\n",
            "Residue", "delta Ca", "delta Ha", "delta CO", "delta Cb");
    for (i = 0; (i < nlist); i++)
    {
        if ((has_dihedral(edPhi, &(dlist[i]))) &&
            (has_dihedral(edPsi, &(dlist[i]))))
        {
            Phi = dlist[i].j0[edPhi];
            Psi = dlist[i].j0[edPsi];
            ca  = cb = co = ha = 0;
            for (j = 0; (j < nf); j++)
            {
                phi = dih[Phi][j];
                psi = dih[Psi][j];

                ca += interpolate(phi, psi, ca_sd);
                cb += interpolate(phi, psi, cb_sd);
                co += interpolate(phi, psi, co_sd);
                ha += interpolate(phi, psi, ha_sd);
            }
            fprintf(fp, "%12s  %10g  %10g  %10g  %10g\n",
                    dlist[i].name, ca/nf, ha/nf, co/nf, cb/nf);
        }
    }
    fprintf(fp, "\n");

    /* Free memory */
    done_shifts(ca_sd);
    done_shifts(cb_sd);
    done_shifts(co_sd);
    done_shifts(ha_sd);
}
