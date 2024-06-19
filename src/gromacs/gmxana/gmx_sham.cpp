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
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/rgb.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;


static int index2(const int* ibox, int x, int y)
{
    return (ibox[1] * x + y);
}

static int index3(const int* ibox, int x, int y, int z)
{
    return (ibox[2] * (ibox[1] * x + y) + z);
}

static int64_t indexn(int ndim, const int* ibox, const int* nxyz)
{
    int64_t d, dd;
    int     k, kk;

    /* Compute index in 1-D array */
    d = 0;
    for (k = 0; (k < ndim); k++)
    {
        dd = nxyz[k];
        for (kk = k + 1; (kk < ndim); kk++)
        {
            dd = dd * ibox[kk];
        }
        d += dd;
    }
    return d;
}

typedef struct
{
    int   Nx;      /* x grid points in unit cell */
    int   Ny;      /* y grid points in unit cell */
    int   Nz;      /* z grid points in unit cell */
    int   dmin[3]; /* starting point x,y,z*/
    int   dmax[3]; /* ending point x,y,z */
    real  cell[6]; /* usual cell parameters */
    real* ed;      /* data */
} XplorMap;

static void lo_write_xplor(XplorMap* map, const char* file)
{
    FILE* fp;
    int   z, i, j, n;

    fp = gmx_ffopen(file, "w");
    /* The REMARKS part is the worst part of the XPLOR format
     * and may cause problems with some programs
     */
    fprintf(fp, "\n       2 !NTITLE\n");
    fprintf(fp, " REMARKS Energy Landscape from GROMACS\n");
    fprintf(fp, " REMARKS DATE: 2004-12-21 \n");
    fprintf(fp,
            " %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",
            map->Nx,
            map->dmin[0],
            map->dmax[0],
            map->Ny,
            map->dmin[1],
            map->dmax[1],
            map->Nz,
            map->dmin[2],
            map->dmax[2]);
    fprintf(fp,
            "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
            map->cell[0],
            map->cell[1],
            map->cell[2],
            map->cell[3],
            map->cell[4],
            map->cell[5]);
    fprintf(fp, "ZYX\n");

    z = map->dmin[2];
    for (n = 0; n < map->Nz; n++, z++)
    {
        fprintf(fp, "%8d\n", z);
        for (i = 0; i < map->Nx * map->Ny; i += 6)
        {
            for (j = 0; j < 6; j++)
            {
                if (i + j < map->Nx * map->Ny)
                {
                    fprintf(fp, "%12.5E", map->ed[n * map->Nx * map->Ny + i + j]);
                }
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "   -9999\n");
    gmx_ffclose(fp);
}

static void write_xplor(const char* file, const real* data, int* ibox, const real dmin[], const real dmax[])
{
    XplorMap* xm;
    int       i, j, k, n;

    snew(xm, 1);
    xm->Nx = ibox[XX];
    xm->Ny = ibox[YY];
    xm->Nz = ibox[ZZ];
    snew(xm->ed, xm->Nx * xm->Ny * xm->Nz);
    n = 0;
    for (k = 0; (k < xm->Nz); k++)
    {
        for (j = 0; (j < xm->Ny); j++)
        {
            for (i = 0; (i < xm->Nx); i++)
            {
                xm->ed[n++] = data[index3(ibox, i, j, k)];
            }
        }
    }
    xm->cell[0] = dmax[XX] - dmin[XX];
    xm->cell[1] = dmax[YY] - dmin[YY];
    xm->cell[2] = dmax[ZZ] - dmin[ZZ];
    xm->cell[3] = xm->cell[4] = xm->cell[5] = 90;

    clear_ivec(xm->dmin);
    xm->dmax[XX] = ibox[XX] - 1;
    xm->dmax[YY] = ibox[YY] - 1;
    xm->dmax[ZZ] = ibox[ZZ] - 1;

    lo_write_xplor(xm, file);

    sfree(xm->ed);
    sfree(xm);
}

static void normalize_p_e(int len, double* P, const int* nbin, real* E, real pmin)
{
    int    i;
    double Ptot = 0;

    for (i = 0; (i < len); i++)
    {
        Ptot += P[i];
        if (nbin[i] > 0)
        {
            E[i] = E[i] / nbin[i];
        }
    }
    printf("Ptot = %g\n", Ptot);
    for (i = 0; (i < len); i++)
    {
        P[i] = P[i] / Ptot;
        /* Have to check for pmin after normalizing to prevent "stretching"
         * the energies.
         */
        if (P[i] < pmin)
        {
            P[i] = 0;
        }
    }
}

typedef struct
{
    int64_t index;
    real    ener;
} t_minimum;

static int comp_minima(const void* a, const void* b)
{
    const t_minimum* ma = reinterpret_cast<const t_minimum*>(a);
    const t_minimum* mb = reinterpret_cast<const t_minimum*>(b);


    if (ma->ener < mb->ener)
    {
        return -1;
    }
    else if (ma->ener > mb->ener)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

static inline void print_minimum(FILE* fp, int num, const t_minimum* min)
{
    fprintf(fp,
            "Minimum %d at index "
            "%" PRId64 " energy %10.3f\n",
            num,
            min->index,
            min->ener);
}

static inline void add_minimum(FILE* fp, int num, const t_minimum* min, t_minimum* mm)
{
    print_minimum(fp, num, min);
    mm[num].index = min->index;
    mm[num].ener  = min->ener;
}

static inline gmx_bool is_local_minimum_from_below(const t_minimum* this_min,
                                                   int              dimension_index,
                                                   int              dimension_min,
                                                   int              neighbour_index,
                                                   const real*      W)
{
    return ((dimension_index == dimension_min)
            || ((dimension_index > dimension_min) && (this_min->ener < W[neighbour_index])));
    /* Note over/underflow within W cannot occur. */
}

static inline gmx_bool is_local_minimum_from_above(const t_minimum* this_min,
                                                   int              dimension_index,
                                                   int              dimension_max,
                                                   int              neighbour_index,
                                                   const real*      W)
{
    return ((dimension_index == dimension_max)
            || ((dimension_index < dimension_max) && (this_min->ener < W[neighbour_index])));
    /* Note over/underflow within W cannot occur. */
}

static void pick_minima(const char* logfile, int* ibox, int ndim, int len, real W[])
{
    FILE*      fp;
    int        i, j, k, nmin;
    t_minimum *mm, this_min;
    int*       this_point;
    int        loopmax, loopcounter;

    snew(mm, len);
    nmin = 0;
    fp   = gmx_ffopen(logfile, "w");
    /* Loop over each element in the array of dimenion ndim seeking
     * minima with respect to every dimension. Specialized loops for
     * speed with ndim == 2 and ndim == 3. */
    switch (ndim)
    {
        case 0:
            /* This is probably impossible to reach anyway. */
            break;
        case 2:
            for (i = 0; (i < ibox[0]); i++)
            {
                for (j = 0; (j < ibox[1]); j++)
                {
                    /* Get the index of this point in the flat array */
                    this_min.index = index2(ibox, i, j);
                    this_min.ener  = W[this_min.index];
                    if (is_local_minimum_from_below(&this_min, i, 0, index2(ibox, i - 1, j), W)
                        && is_local_minimum_from_above(&this_min, i, ibox[0] - 1, index2(ibox, i + 1, j), W)
                        && is_local_minimum_from_below(&this_min, j, 0, index2(ibox, i, j - 1), W)
                        && is_local_minimum_from_above(&this_min, j, ibox[1] - 1, index2(ibox, i, j + 1), W))
                    {
                        add_minimum(fp, nmin, &this_min, mm);
                        nmin++;
                    }
                }
            }
            break;
        case 3:
            for (i = 0; (i < ibox[0]); i++)
            {
                for (j = 0; (j < ibox[1]); j++)
                {
                    for (k = 0; (k < ibox[2]); k++)
                    {
                        /* Get the index of this point in the flat array */
                        this_min.index = index3(ibox, i, j, k);
                        this_min.ener  = W[this_min.index];
                        if (is_local_minimum_from_below(&this_min, i, 0, index3(ibox, i - 1, j, k), W)
                            && is_local_minimum_from_above(
                                    &this_min, i, ibox[0] - 1, index3(ibox, i + 1, j, k), W)
                            && is_local_minimum_from_below(&this_min, j, 0, index3(ibox, i, j - 1, k), W)
                            && is_local_minimum_from_above(
                                    &this_min, j, ibox[1] - 1, index3(ibox, i, j + 1, k), W)
                            && is_local_minimum_from_below(&this_min, k, 0, index3(ibox, i, j, k - 1), W)
                            && is_local_minimum_from_above(
                                    &this_min, k, ibox[2] - 1, index3(ibox, i, j, k + 1), W))
                        {
                            add_minimum(fp, nmin, &this_min, mm);
                            nmin++;
                        }
                    }
                }
            }
            break;
        default:
            /* Note this treats ndim == 1 and ndim > 3 */

            /* Set up an ndim-dimensional vector to loop over the points
             * on the grid. (0,0,0, ... 0) is an acceptable place to
             * start. */
            snew(this_point, ndim);

            /* Determine the number of points of the ndim-dimensional
             * grid. */
            loopmax = ibox[0];
            for (i = 1; i < ndim; i++)
            {
                loopmax *= ibox[i];
            }

            loopcounter = 0;
            while (loopmax > loopcounter)
            {
                gmx_bool bMin = TRUE;

                /* Get the index of this_point in the flat array */
                this_min.index = indexn(ndim, ibox, this_point);
                this_min.ener  = W[this_min.index];

                /* Is this_point a minimum from above and below in each
                 * dimension? */
                for (i = 0; bMin && (i < ndim); i++)
                {
                    /* Save the index of this_point within the curent
                     * dimension so we can change that index in the
                     * this_point array for use with indexn(). */
                    int index = this_point[i];
                    this_point[i]--;
                    bMin = bMin
                           && is_local_minimum_from_below(
                                   &this_min, index, 0, indexn(ndim, ibox, this_point), W);
                    this_point[i] += 2;
                    bMin = bMin
                           && is_local_minimum_from_above(
                                   &this_min, index, ibox[i] - 1, indexn(ndim, ibox, this_point), W);
                    this_point[i]--;
                }
                if (bMin)
                {
                    add_minimum(fp, nmin, &this_min, mm);
                    nmin++;
                }

                /* update global loop counter */
                loopcounter++;

                /* Avoid underflow of this_point[i] */
                if (loopmax > loopcounter)
                {
                    /* update this_point non-recursively */
                    i = ndim - 1;
                    this_point[i]++;
                    while (ibox[i] == this_point[i])
                    {
                        this_point[i] = 0;
                        i--;
                        /* this_point[i] cannot underflow because
                         * loopmax > loopcounter. */
                        this_point[i]++;
                    }
                }
            }

            sfree(this_point);
            break;
    }
    qsort(mm, nmin, sizeof(mm[0]), comp_minima);
    fprintf(fp, "Minima sorted after energy\n");
    for (i = 0; (i < nmin); i++)
    {
        print_minimum(fp, i, &mm[i]);
    }
    gmx_ffclose(fp);
    sfree(mm);
}

static void do_sham(const char* fn,
                    const char* ndx,
                    const char* xpmP,
                    const char* xpm,
                    const char* xpm2,
                    const char* xpm3,
                    const char* pdb,
                    const char* logf,
                    int         n,
                    int         neig,
                    real**      eig,
                    gmx_bool    bGE,
                    int         nenerT,
                    real**      enerT,
                    real        Tref,
                    real        pmax,
                    real        gmax,
                    const real* emin,
                    const real* emax,
                    int         nlevels,
                    real        pmin,
                    const int*  idim,
                    int*        ibox,
                    gmx_bool    bXmin,
                    real*       xmin,
                    gmx_bool    bXmax,
                    real*       xmax)
{
    FILE*        fp;
    real *       min_eig, *max_eig;
    real *       axis_x, *axis_y, *axis_z, *axis = nullptr;
    double*      P;
    real **      PP, *W, *E, **WW, **EE, *S, **SS, *M, *bE;
    rvec         xxx;
    char*        buf;
    double *     bfac, efac, bref, Pmax, Wmin, Wmax, Winf, Emin, Emax, Einf, Smin, Smax, Sinf;
    real*        delta;
    int          i, j, k, imin, len, index, *nbin, *bindex, bi;
    int *        nxyz, maxbox;
    t_blocka*    b;
    gmx_bool     bOutside;
    unsigned int flags;
    t_rgb        rlo = { 0, 0, 0 };
    t_rgb        rhi = { 1, 1, 1 };

    /* Determine extremes for the eigenvectors */
    snew(min_eig, neig);
    snew(max_eig, neig);
    snew(nxyz, neig);
    snew(bfac, neig);
    snew(delta, neig);

    for (i = 0; (i < neig); i++)
    {
        /* Check for input constraints */
        min_eig[i] = max_eig[i] = eig[i][0];
        for (j = 0; (j < n); j++)
        {
            min_eig[i] = std::min(min_eig[i], eig[i][j]);
            max_eig[i] = std::max(max_eig[i], eig[i][j]);
            delta[i]   = (max_eig[i] - min_eig[i]) / (2.0 * ibox[i]);
        }
        /* Add some extra space, half a bin on each side, unless the
         * user has set the limits.
         */
        if (bXmax)
        {
            if (max_eig[i] > xmax[i])
            {
                gmx_warning("Your xmax[%d] value %f is smaller than the largest data point %f",
                            i,
                            xmax[i],
                            max_eig[i]);
            }
            max_eig[i] = xmax[i];
        }
        else
        {
            max_eig[i] += delta[i];
        }

        if (bXmin)
        {
            if (min_eig[i] < xmin[i])
            {
                gmx_warning("Your xmin[%d] value %f is larger than the smallest data point %f",
                            i,
                            xmin[i],
                            min_eig[i]);
            }
            min_eig[i] = xmin[i];
        }
        else
        {
            min_eig[i] -= delta[i];
        }
        bfac[i] = ibox[i] / (max_eig[i] - min_eig[i]);
    }
    /* Do the binning */
    bref = 1 / (gmx::c_boltz * Tref);
    snew(bE, n);
    if (bGE || nenerT == 2)
    {
        assert(enerT);
        Emin = 1e8;
        for (j = 0; (j < n); j++)
        {
            if (bGE)
            {
                bE[j] = bref * enerT[0][j];
            }
            else
            {
                bE[j] = (bref - 1 / (gmx::c_boltz * enerT[1][j])) * enerT[0][j];
            }
            Emin = std::min(Emin, static_cast<double>(bE[j]));
        }
    }
    else
    {
        Emin = 0;
    }
    len = 1;
    for (i = 0; (i < neig); i++)
    {
        len = len * ibox[i];
    }
    printf("There are %d bins in the %d-dimensional histogram. Beta-Emin = %g\n", len, neig, Emin);
    snew(P, len);
    snew(W, len);
    snew(E, len);
    snew(S, len);
    snew(M, len);
    snew(nbin, len);
    snew(bindex, n);


    /* Loop over projections */
    for (j = 0; (j < n); j++)
    {
        /* Loop over dimensions */
        bOutside = FALSE;
        for (i = 0; (i < neig); i++)
        {
            nxyz[i] = static_cast<int>(bfac[i] * (eig[i][j] - min_eig[i]));
            if (nxyz[i] < 0 || nxyz[i] >= ibox[i])
            {
                bOutside = TRUE;
            }
        }
        if (!bOutside)
        {
            index = indexn(neig, ibox, nxyz);
            range_check(index, 0, len);
            /* Compute the exponential factor */
            if (enerT)
            {
                efac = std::exp(-bE[j] + Emin);
            }
            else
            {
                efac = 1;
            }
            /* Apply the bin volume correction for a multi-dimensional distance */
            for (i = 0; i < neig; i++)
            {
                if (idim[i] == 2)
                {
                    efac /= eig[i][j];
                }
                else if (idim[i] == 3)
                {
                    efac /= gmx::square(eig[i][j]);
                }
                else if (idim[i] == -1)
                {
                    efac /= std::sin(gmx::c_deg2Rad * eig[i][j]);
                }
            }
            /* Update the probability */
            P[index] += efac;
            /* Update the energy */
            if (enerT)
            {
                E[index] += enerT[0][j];
            }
            /* Statistics: which "structure" in which bin */
            nbin[index]++;
            bindex[j] = index;
        }
    }
    /* Normalize probability */
    normalize_p_e(len, P, nbin, E, pmin);
    Pmax = 0;
    /* Compute boundaries for the Free energy */
    Wmin = 1e8;
    imin = -1;
    Wmax = -1e8;
    /* Recompute Emin: it may have changed due to averaging */
    Emin = 1e8;
    Emax = -1e8;
    for (i = 0; (i < len); i++)
    {
        if (P[i] != 0)
        {
            Pmax = std::max(P[i], Pmax);
            W[i] = -gmx::c_boltz * Tref * std::log(P[i]);
            if (W[i] < Wmin)
            {
                Wmin = W[i];
                imin = i;
            }
            Emin = std::min(static_cast<double>(E[i]), Emin);
            Emax = std::max(static_cast<double>(E[i]), Emax);
            Wmax = std::max(static_cast<double>(W[i]), Wmax);
        }
    }
    if (pmax > 0)
    {
        Pmax = pmax;
    }
    if (gmax > 0)
    {
        Wmax = gmax;
    }
    else
    {
        Wmax -= Wmin;
    }
    Winf = Wmax + 1;
    Einf = Emax + 1;
    Smin = Emin - Wmax;
    Smax = Emax - Smin;
    Sinf = Smax + 1;
    /* Write out the free energy as a function of bin index */
    fp = gmx_ffopen(fn, "w");
    for (i = 0; (i < len); i++)
    {
        if (P[i] != 0)
        {
            W[i] -= Wmin;
            S[i] = E[i] - W[i] - Smin;
            fprintf(fp, "%5d  %10.5e  %10.5e  %10.5e\n", i, W[i], E[i], S[i]);
        }
        else
        {
            W[i] = Winf;
            E[i] = Einf;
            S[i] = Sinf;
        }
    }
    gmx_ffclose(fp);
    /* Organize the structures in the bins */
    snew(b, 1);
    snew(b->index, len + 1);
    snew(b->a, n);
    b->index[0] = 0;
    for (i = 0; (i < len); i++)
    {
        b->index[i + 1] = b->index[i] + nbin[i];
        nbin[i]         = 0;
    }
    for (i = 0; (i < n); i++)
    {
        bi                            = bindex[i];
        b->a[b->index[bi] + nbin[bi]] = i;
        nbin[bi]++;
    }
    /* Consistency check */
    /* This no longer applies when we allow the plot to be smaller
       than the sampled space.
       for(i=0; (i<len); i++) {
       if (nbin[i] != (b->index[i+1] - b->index[i]))
        gmx_fatal(FARGS,"nbin[%d] = %d, should be %d",i,nbin[i],
          b->index[i+1] - b->index[i]);
       }
     */
    /* Write the index file */
    fp = gmx_ffopen(ndx, "w");
    for (i = 0; (i < len); i++)
    {
        if (nbin[i] > 0)
        {
            fprintf(fp, "[ %d ]\n", i);
            for (j = b->index[i]; (j < b->index[i + 1]); j++)
            {
                fprintf(fp, "%d\n", b->a[j] + 1);
            }
        }
    }
    gmx_ffclose(fp);
    snew(axis_x, ibox[0] + 1);
    snew(axis_y, ibox[1] + 1);
    snew(axis_z, ibox[2] + 1);
    maxbox = std::max(ibox[0], std::max(ibox[1], ibox[2]));
    snew(PP, maxbox * maxbox);
    snew(WW, maxbox * maxbox);
    snew(EE, maxbox * maxbox);
    snew(SS, maxbox * maxbox);
    for (i = 0; (i < std::min(neig, 3)); i++)
    {
        switch (i)
        {
            case 0: axis = axis_x; break;
            case 1: axis = axis_y; break;
            case 2: axis = axis_z; break;
            default: break;
        }
        for (j = 0; j <= ibox[i]; j++)
        {
            axis[j] = min_eig[i] + j / bfac[i];
        }
    }

    pick_minima(logf, ibox, neig, len, W);
    if (gmax <= 0)
    {
        gmax = Winf;
    }
    flags = MAT_SPATIAL_X | MAT_SPATIAL_Y;
    if (neig == 2)
    {
        /* Dump to XPM file */
        snew(PP, ibox[0]);
        for (i = 0; (i < ibox[0]); i++)
        {
            snew(PP[i], ibox[1]);
            for (j = 0; j < ibox[1]; j++)
            {
                PP[i][j] = P[i * ibox[1] + j];
            }
            WW[i] = &(W[i * ibox[1]]);
            EE[i] = &(E[i * ibox[1]]);
            SS[i] = &(S[i * ibox[1]]);
        }
        fp = gmx_ffopen(xpmP, "w");
        write_xpm(fp,
                  flags,
                  "Probability Distribution",
                  "",
                  "PC1",
                  "PC2",
                  ibox[0],
                  ibox[1],
                  axis_x,
                  axis_y,
                  PP,
                  0,
                  Pmax,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
        fp = gmx_ffopen(xpm, "w");
        write_xpm(fp,
                  flags,
                  "Gibbs Energy Landscape",
                  "G (kJ/mol)",
                  "PC1",
                  "PC2",
                  ibox[0],
                  ibox[1],
                  axis_x,
                  axis_y,
                  WW,
                  0,
                  gmax,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
        fp = gmx_ffopen(xpm2, "w");
        write_xpm(fp,
                  flags,
                  "Enthalpy Landscape",
                  "H (kJ/mol)",
                  "PC1",
                  "PC2",
                  ibox[0],
                  ibox[1],
                  axis_x,
                  axis_y,
                  EE,
                  emin ? *emin : Emin,
                  emax ? *emax : Einf,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
        fp = gmx_ffopen(xpm3, "w");
        write_xpm(fp,
                  flags,
                  "Entropy Landscape",
                  "TDS (kJ/mol)",
                  "PC1",
                  "PC2",
                  ibox[0],
                  ibox[1],
                  axis_x,
                  axis_y,
                  SS,
                  0,
                  Sinf,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
    }
    else if (neig == 3)
    {
        /* Dump to PDB file */
        fp = gmx_ffopen(pdb, "w");
        for (i = 0; (i < ibox[0]); i++)
        {
            xxx[XX] = 3 * i + 1.5 * (1 - ibox[0]);
            for (j = 0; (j < ibox[1]); j++)
            {
                xxx[YY] = 3 * j + 1.5 * (1 - ibox[1]);
                for (k = 0; (k < ibox[2]); k++)
                {
                    xxx[ZZ] = 3 * k + 1.5 * (1 - ibox[2]);
                    index   = index3(ibox, i, j, k);
                    if (P[index] > 0)
                    {
                        fprintf(fp,
                                "%-6s%5d  %-4.4s%3.3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                "ATOM",
                                (index + 1) % 10000,
                                "H",
                                "H",
                                (index + 1) % 10000,
                                xxx[XX],
                                xxx[YY],
                                xxx[ZZ],
                                1.0,
                                W[index]);
                    }
                }
            }
        }
        gmx_ffclose(fp);
        write_xplor("out.xplor", W, ibox, min_eig, max_eig);
        nxyz[XX] = imin / (ibox[1] * ibox[2]);
        nxyz[YY] = (imin - nxyz[XX] * ibox[1] * ibox[2]) / ibox[2];
        nxyz[ZZ] = imin % ibox[2];
        for (i = 0; (i < ibox[0]); i++)
        {
            snew(WW[i], maxbox);
            for (j = 0; (j < ibox[1]); j++)
            {
                WW[i][j] = W[index3(ibox, i, j, nxyz[ZZ])];
            }
        }
        snew(buf, std::strlen(xpm) + 4);
        sprintf(buf, "%s", xpm);
        sprintf(&buf[std::strlen(xpm) - 4], "12.xpm");
        fp = gmx_ffopen(buf, "w");
        write_xpm(fp,
                  flags,
                  "Gibbs Energy Landscape",
                  "W (kJ/mol)",
                  "PC1",
                  "PC2",
                  ibox[0],
                  ibox[1],
                  axis_x,
                  axis_y,
                  WW,
                  0,
                  gmax,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
        for (i = 0; (i < ibox[0]); i++)
        {
            for (j = 0; (j < ibox[2]); j++)
            {
                WW[i][j] = W[index3(ibox, i, nxyz[YY], j)];
            }
        }
        sprintf(&buf[std::strlen(xpm) - 4], "13.xpm");
        fp = gmx_ffopen(buf, "w");
        write_xpm(fp,
                  flags,
                  "SHAM Energy Landscape",
                  "kJ/mol",
                  "PC1",
                  "PC3",
                  ibox[0],
                  ibox[2],
                  axis_x,
                  axis_z,
                  WW,
                  0,
                  gmax,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
        for (i = 0; (i < ibox[1]); i++)
        {
            for (j = 0; (j < ibox[2]); j++)
            {
                WW[i][j] = W[index3(ibox, nxyz[XX], i, j)];
            }
        }
        sprintf(&buf[std::strlen(xpm) - 4], "23.xpm");
        fp = gmx_ffopen(buf, "w");
        write_xpm(fp,
                  flags,
                  "SHAM Energy Landscape",
                  "kJ/mol",
                  "PC2",
                  "PC3",
                  ibox[1],
                  ibox[2],
                  axis_y,
                  axis_z,
                  WW,
                  0,
                  gmax,
                  rlo,
                  rhi,
                  &nlevels);
        gmx_ffclose(fp);
        sfree(buf);
    }
}

static void ehisto(const char* fh, int n, real** enerT, const gmx_output_env_t* oenv)
{
    FILE* fp;
    int   i, j, k, nbin, blength;
    int*  bindex;
    real *T, bmin, bmax, bwidth;
    int** histo;

    bmin = 1e8;
    bmax = -1e8;
    snew(bindex, n);
    snew(T, n);
    nbin = 0;
    for (j = 1; (j < n); j++)
    {
        for (k = 0; (k < nbin); k++)
        {
            if (T[k] == enerT[1][j])
            {
                bindex[j] = k;
                break;
            }
        }
        if (k == nbin)
        {
            bindex[j] = nbin;
            T[nbin]   = enerT[1][j];
            nbin++;
        }
        bmin = std::min(enerT[0][j], bmin);
        bmax = std::max(enerT[0][j], bmax);
    }
    bwidth  = 1.0;
    blength = static_cast<int>((bmax - bmin) / bwidth + 2);
    snew(histo, nbin);
    for (i = 0; (i < nbin); i++)
    {
        snew(histo[i], blength);
    }
    for (j = 0; (j < n); j++)
    {
        k = static_cast<int>((enerT[0][j] - bmin) / bwidth);
        histo[bindex[j]][k]++;
    }
    fp = xvgropen(fh, "Energy distribution", "E (kJ/mol)", "", oenv);
    for (j = 0; (j < blength); j++)
    {
        fprintf(fp, "%8.3f", bmin + j * bwidth);
        for (k = 0; (k < nbin); k++)
        {
            fprintf(fp, "  %6d", histo[k][j]);
        }
        fprintf(fp, "\n");
    }
    xvgrclose(fp);
}

int gmx_sham(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] makes multi-dimensional free-energy, enthalpy and entropy plots.",
        "[THISMODULE] reads one or more [REF].xvg[ref] files and analyzes data sets.",
        "The basic purpose of [THISMODULE] is to plot Gibbs free energy landscapes",
        "(option [TT]-ls[tt])",
        "by Bolzmann inverting multi-dimensional histograms (option [TT]-lp[tt]),",
        "but it can also",
        "make enthalpy (option [TT]-lsh[tt]) and entropy (option [TT]-lss[tt])",
        "plots. The histograms can be made for any quantities the user supplies.",
        "A line in the input file may start with a time",
        "(see option [TT]-time[tt]) and any number of [IT]y[it]-values may follow.",
        "Multiple sets can also be",
        "read when they are separated by & (option [TT]-n[tt]),",
        "in this case only one [IT]y[it]-value is read from each line.",
        "All lines starting with # and @ are skipped.",
        "[PAR]",
        "Option [TT]-ge[tt] can be used to supply a file with free energies",
        "when the ensemble is not a Boltzmann ensemble, but needs to be biased",
        "by this free energy. One free energy value is required for each",
        "(multi-dimensional) data point in the [TT]-f[tt] input.",
        "[PAR]",
        "Option [TT]-ene[tt] can be used to supply a file with energies.",
        "These energies are used as a weighting function in the single",
        "histogram analysis method by Kumar et al. When temperatures",
        "are supplied (as a second column in the file), an experimental",
        "weighting scheme is applied. In addition the vales",
        "are used for making enthalpy and entropy plots.",
        "[PAR]",
        "With option [TT]-dim[tt], dimensions can be gives for distances.",
        "When a distance is 2- or 3-dimensional, the circumference or surface",
        "sampled by two particles increases with increasing distance.",
        "Depending on what one would like to show, one can choose to correct",
        "the histogram and free-energy for this volume effect.",
        "The probability is normalized by r and r^2 for dimensions of 2 and 3, ",
        "respectively.",
        "A value of -1 is used to indicate an angle in degrees between two",
        "vectors: a sin(angle) normalization will be applied.",
        "[BB]Note[bb] that for angles between vectors the inner-product or cosine",
        "is the natural quantity to use, as it will produce bins of the same",
        "volume."
    };
    static real     tb = -1, te = -1;
    static gmx_bool bHaveT = TRUE, bDer = FALSE;
    static gmx_bool bSham = TRUE;
    static real     Tref = 298.15, pmin = 0, ttol = 0, pmax = 0, gmax = 0, emin = 0, emax = 0;
    static rvec     nrdim = { 1, 1, 1 };
    static rvec     nrbox = { 32, 32, 32 };
    static rvec     xmin = { 0, 0, 0 }, xmax = { 1, 1, 1 };
    static int      nsets_in = 1, nlevels = 25;
    t_pargs         pa[] = {
        { "-time", FALSE, etBOOL, { &bHaveT }, "Expect a time in the input" },
        { "-b", FALSE, etREAL, { &tb }, "First time to read from set" },
        { "-e", FALSE, etREAL, { &te }, "Last time to read from set" },
        { "-ttol", FALSE, etREAL, { &ttol }, "Tolerance on time in appropriate units (usually ps)" },
        { "-n",
          FALSE,
          etINT,
          { &nsets_in },
          "Read this number of sets separated by lines containing only an ampersand" },
        { "-d", FALSE, etBOOL, { &bDer }, "Use the derivative" },
        { "-sham",
          FALSE,
          etBOOL,
          { &bSham },
          "Turn off energy weighting even if energies are given" },
        { "-tsham", FALSE, etREAL, { &Tref }, "Temperature for single histogram analysis" },
        { "-pmin",
          FALSE,
          etREAL,
          { &pmin },
          "Minimum probability. Anything lower than this will be set to zero" },
        { "-dim",
          FALSE,
          etRVEC,
          { nrdim },
          "Dimensions for distances, used for volume correction (max 3 values, dimensions > 3 will "
          "get the same value as the last)" },
        { "-ngrid",
          FALSE,
          etRVEC,
          { nrbox },
          "Number of bins for energy landscapes (max 3 values, dimensions > 3 will get the same "
          "value as the last)" },
        { "-xmin",
          FALSE,
          etRVEC,
          { xmin },
          "Minimum for the axes in energy landscape (see above for > 3 dimensions)" },
        { "-xmax",
          FALSE,
          etRVEC,
          { xmax },
          "Maximum for the axes in energy landscape (see above for > 3 dimensions)" },
        { "-pmax", FALSE, etREAL, { &pmax }, "Maximum probability in output, default is calculate" },
        { "-gmax", FALSE, etREAL, { &gmax }, "Maximum free energy in output, default is calculate" },
        { "-emin", FALSE, etREAL, { &emin }, "Minimum enthalpy in output, default is calculate" },
        { "-emax", FALSE, etREAL, { &emax }, "Maximum enthalpy in output, default is calculate" },
        { "-nlevels", FALSE, etINT, { &nlevels }, "Number of levels for energy landscape" },
    };
#define NPA asize(pa)

    int               n, e_n, nset, e_nset = 0, i, *idim, *ibox;
    real **           val, **et_val, *t, *e_t, e_dt, dt;
    real *            rmin, *rmax;
    const char *      fn_ge, *fn_ene;
    gmx_output_env_t* oenv;
    int64_t           num_grid_points;

    t_filenm fnm[] = {
        { efXVG, "-f", "graph", ffREAD },       { efXVG, "-ge", "gibbs", ffOPTRD },
        { efXVG, "-ene", "esham", ffOPTRD },    { efXVG, "-dist", "ener", ffOPTWR },
        { efXVG, "-histo", "edist", ffOPTWR },  { efNDX, "-bin", "bindex", ffOPTWR },
        { efXPM, "-lp", "prob", ffOPTWR },      { efXPM, "-ls", "gibbs", ffOPTWR },
        { efXPM, "-lsh", "enthalpy", ffOPTWR }, { efXPM, "-lss", "entropy", ffOPTWR },
        { efPDB, "-ls3", "gibbs3", ffOPTWR },   { efLOG, "-g", "shamlog", ffOPTWR }
    };
#define NFILE asize(fnm)

    int npargs;

    npargs = asize(pa);
    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW, NFILE, fnm, npargs, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    val = read_xvg_time(opt2fn("-f", NFILE, fnm),
                        bHaveT,
                        opt2parg_bSet("-b", npargs, pa),
                        tb - ttol,
                        opt2parg_bSet("-e", npargs, pa),
                        te + ttol,
                        nsets_in,
                        &nset,
                        &n,
                        &dt,
                        &t);
    printf("Read %d sets of %d points, dt = %g\n\n", nset, n, dt);

    fn_ge  = opt2fn_null("-ge", NFILE, fnm);
    fn_ene = opt2fn_null("-ene", NFILE, fnm);

    if (fn_ge && fn_ene)
    {
        gmx_fatal(FARGS, "Can not do free energy and energy corrections at the same time");
    }

    if (fn_ge || fn_ene)
    {
        et_val = read_xvg_time(fn_ge ? fn_ge : fn_ene,
                               bHaveT,
                               opt2parg_bSet("-b", npargs, pa),
                               tb - ttol,
                               opt2parg_bSet("-e", npargs, pa),
                               te + ttol,
                               1,
                               &e_nset,
                               &e_n,
                               &e_dt,
                               &e_t);
        if (fn_ge)
        {
            if (e_nset != 1)
            {
                gmx_fatal(FARGS, "Can only handle one free energy component in %s", fn_ge);
            }
        }
        else
        {
            if (e_nset != 1 && e_nset != 2)
            {
                gmx_fatal(FARGS,
                          "Can only handle one energy component or one energy and one T in %s",
                          fn_ene);
            }
        }
        if (e_n != n)
        {
            gmx_fatal(FARGS,
                      "Number of energies (%d) does not match number of entries (%d) in %s",
                      e_n,
                      n,
                      opt2fn("-f", NFILE, fnm));
        }
    }
    else
    {
        et_val = nullptr;
    }

    if (fn_ene && et_val)
    {
        ehisto(opt2fn("-histo", NFILE, fnm), e_n, et_val, oenv);
    }

    snew(idim, std::max(3, nset));
    snew(ibox, std::max(3, nset));
    snew(rmin, std::max(3, nset));
    snew(rmax, std::max(3, nset));
    for (i = 0; (i < std::min(3, nset)); i++)
    {
        idim[i] = static_cast<int>(nrdim[i]);
        ibox[i] = static_cast<int>(nrbox[i]);
        rmin[i] = xmin[i];
        rmax[i] = xmax[i];
    }
    for (; (i < nset); i++)
    {
        idim[i] = static_cast<int>(nrdim[2]);
        ibox[i] = static_cast<int>(nrbox[2]);
        rmin[i] = xmin[2];
        rmax[i] = xmax[2];
    }

    /* Check that the grid size is manageable. */
    num_grid_points = ibox[0];
    for (i = 1; i < nset; i++)
    {
        int64_t result;
        if (!check_int_multiply_for_overflow(num_grid_points, ibox[i], &result))
        {
            gmx_fatal(FARGS,
                      "The number of dimensions and grid points is too large for this tool.\n");
        }
        num_grid_points = result;
    }
    /* The number of grid points fits in a int64_t. */

    do_sham(opt2fn("-dist", NFILE, fnm),
            opt2fn("-bin", NFILE, fnm),
            opt2fn("-lp", NFILE, fnm),
            opt2fn("-ls", NFILE, fnm),
            opt2fn("-lsh", NFILE, fnm),
            opt2fn("-lss", NFILE, fnm),
            opt2fn("-ls3", NFILE, fnm),
            opt2fn("-g", NFILE, fnm),
            n,
            nset,
            val,
            fn_ge != nullptr,
            e_nset,
            et_val,
            Tref,
            pmax,
            gmax,
            opt2parg_bSet("-emin", NPA, pa) ? &emin : nullptr,
            opt2parg_bSet("-emax", NPA, pa) ? &emax : nullptr,
            nlevels,
            pmin,
            idim,
            ibox,
            opt2parg_bSet("-xmin", NPA, pa),
            rmin,
            opt2parg_bSet("-xmax", NPA, pa),
            rmax);

    return 0;
}
