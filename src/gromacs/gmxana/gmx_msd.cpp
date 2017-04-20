/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <cmath>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#define FACTOR  1000.0  /* Convert nm^2/ps to 10e-5 cm^2/s */
/* NORMAL = total diffusion coefficient (default). X,Y,Z is diffusion
   coefficient in X,Y,Z direction. LATERAL is diffusion coefficient in
   plane perpendicular to axis
 */
typedef enum {
    NOT_USED, NORMAL, X, Y, Z, LATERAL
} msd_type;

typedef struct {
    real          t0;         /* start time and time increment between  */
    real          delta_t;    /* time between restart points */
    real          beginfit,   /* the begin/end time for fits as reals between */
                  endfit;     /* 0 and 1 */
    real          dim_factor; /* the dimensionality factor for the diffusion
                                 constant */
    real        **data;       /* the displacement data. First index is the group
                                 number, second is frame number */
    real         *time;       /* frame time */
    real         *mass;       /* masses for mass-weighted msd */
    matrix      **datam;
    rvec        **x0;         /* original positions */
    rvec         *com;        /* center of mass correction for each frame */
    gmx_stats_t **lsq;        /* fitting stats for individual molecule msds */
    msd_type      type;       /* the type of msd to calculate (lateral, etc.)*/
    int           axis;       /* the axis along which to calculate */
    int           ncoords;
    int           nrestart;   /* number of restart points */
    int           nmol;       /* number of molecules (for bMol) */
    int           nframes;    /* number of frames */
    int           nlast;
    int           ngrp;       /* number of groups to use for msd calculation */
    int          *n_offs;
    int         **ndata;      /* the number of msds (particles/mols) per data
                                 point. */
} t_corr;

typedef real t_calc_func (t_corr *curr, int nx, int index[], int nx0, rvec xc[],
                          rvec dcom, gmx_bool bTen, matrix mat);

static real thistime(t_corr *curr)
{
    return curr->time[curr->nframes];
}

static gmx_bool in_data(t_corr *curr, int nx00)
{
    return curr->nframes-curr->n_offs[nx00];
}

t_corr *init_corr(int nrgrp, int type, int axis, real dim_factor,
                  int nmol, gmx_bool bTen, gmx_bool bMass, real dt, const t_topology *top,
                  real beginfit, real endfit)
{
    t_corr  *curr;
    int      i;

    snew(curr, 1);
    curr->type       = (msd_type)type;
    curr->axis       = axis;
    curr->ngrp       = nrgrp;
    curr->nrestart   = 0;
    curr->delta_t    = dt;
    curr->beginfit   = (1 - 2*GMX_REAL_EPS)*beginfit;
    curr->endfit     = (1 + 2*GMX_REAL_EPS)*endfit;
    curr->x0         = nullptr;
    curr->n_offs     = nullptr;
    curr->nframes    = 0;
    curr->nlast      = 0;
    curr->dim_factor = dim_factor;

    snew(curr->ndata, nrgrp);
    snew(curr->data, nrgrp);
    if (bTen)
    {
        snew(curr->datam, nrgrp);
    }
    for (i = 0; (i < nrgrp); i++)
    {
        curr->ndata[i] = nullptr;
        curr->data[i]  = nullptr;
        if (bTen)
        {
            curr->datam[i] = nullptr;
        }
    }
    curr->time = nullptr;
    curr->lsq  = nullptr;
    curr->nmol = nmol;
    if (curr->nmol > 0)
    {
        snew(curr->mass, curr->nmol);
        for (i = 0; i < curr->nmol; i++)
        {
            curr->mass[i] = 1;
        }
    }
    else
    {
        if (bMass)
        {
            const t_atoms *atoms = &top->atoms;
            snew(curr->mass, atoms->nr);
            for (i = 0; (i < atoms->nr); i++)
            {
                curr->mass[i] = atoms->atom[i].m;
            }
        }
    }

    return curr;
}

static void corr_print(t_corr *curr, gmx_bool bTen, const char *fn, const char *title,
                       const char *yaxis,
                       real msdtime, real beginfit, real endfit,
                       real *DD, real *SigmaD, char *grpname[],
                       const gmx_output_env_t *oenv)
{
    FILE *out;
    int   i, j;

    out = xvgropen(fn, title, output_env_get_xvgr_tlabel(oenv), yaxis, oenv);
    if (DD)
    {
        fprintf(out, "# MSD gathered over %g %s with %d restarts\n",
                msdtime, output_env_get_time_unit(oenv).c_str(), curr->nrestart);
        fprintf(out, "# Diffusion constants fitted from time %g to %g %s\n",
                beginfit, endfit, output_env_get_time_unit(oenv).c_str());
        for (i = 0; i < curr->ngrp; i++)
        {
            fprintf(out, "# D[%10s] = %.4f (+/- %.4f) (1e-5 cm^2/s)\n",
                    grpname[i], DD[i], SigmaD[i]);
        }
    }
    for (i = 0; i < curr->nframes; i++)
    {
        fprintf(out, "%10g", output_env_conv_time(oenv, curr->time[i]));
        for (j = 0; j < curr->ngrp; j++)
        {
            fprintf(out, "  %10g", curr->data[j][i]);
            if (bTen)
            {
                fprintf(out, " %10g %10g %10g %10g %10g %10g",
                        curr->datam[j][i][XX][XX],
                        curr->datam[j][i][YY][YY],
                        curr->datam[j][i][ZZ][ZZ],
                        curr->datam[j][i][YY][XX],
                        curr->datam[j][i][ZZ][XX],
                        curr->datam[j][i][ZZ][YY]);
            }
        }
        fprintf(out, "\n");
    }
    xvgrclose(out);
}

/* called from corr_loop, to do the main calculations */
static void calc_corr(t_corr *curr, int nr, int nx, int index[], rvec xc[],
                      gmx_bool bRmCOMM, rvec com, t_calc_func *calc1, gmx_bool bTen)
{
    int    nx0;
    real   g;
    matrix mat;
    rvec   dcom;

    /* Check for new starting point */
    if (curr->nlast < curr->nrestart)
    {
        if ((thistime(curr) >= (curr->nlast*curr->delta_t)) && (nr == 0))
        {
            std::memcpy(curr->x0[curr->nlast], xc, curr->ncoords*sizeof(xc[0]));
            curr->n_offs[curr->nlast] = curr->nframes;
            copy_rvec(com, curr->com[curr->nlast]);
            curr->nlast++;
        }
    }

    /* nx0 appears to be the number of new starting points,
     * so for all starting points, call calc1.
     */
    for (nx0 = 0; (nx0 < curr->nlast); nx0++)
    {
        if (bRmCOMM)
        {
            rvec_sub(com, curr->com[nx0], dcom);
        }
        else
        {
            clear_rvec(dcom);
        }
        g = calc1(curr, nx, index, nx0, xc, dcom, bTen, mat);
#ifdef DEBUG2
        printf("g[%d]=%g\n", nx0, g);
#endif
        curr->data[nr][in_data(curr, nx0)] += g;
        if (bTen)
        {
            m_add(curr->datam[nr][in_data(curr, nx0)], mat,
                  curr->datam[nr][in_data(curr, nx0)]);
        }
        curr->ndata[nr][in_data(curr, nx0)]++;
    }
}

/* the non-mass-weighted mean-squared displacement calcuation */
static real calc1_norm(t_corr *curr, int nx, int index[], int nx0, rvec xc[],
                       rvec dcom, gmx_bool bTen, matrix mat)
{
    int  i, ix, m, m2;
    real g, r, r2;
    rvec rv;

    g = 0.0;
    clear_mat(mat);

    for (i = 0; (i < nx); i++)
    {
        ix = index[i];
        r2 = 0.0;
        switch (curr->type)
        {
            case NORMAL:
                for (m = 0; (m < DIM); m++)
                {
                    rv[m] = xc[ix][m] - curr->x0[nx0][ix][m] - dcom[m];
                    r2   += rv[m]*rv[m];
                    if (bTen)
                    {
                        for (m2 = 0; m2 <= m; m2++)
                        {
                            mat[m][m2] += rv[m]*rv[m2];
                        }
                    }
                }
                break;
            case X:
            case Y:
            case Z:
                r = xc[ix][curr->type-X] - curr->x0[nx0][ix][curr->type-X] -
                    dcom[curr->type-X];
                r2 += r*r;
                break;
            case LATERAL:
                for (m = 0; (m < DIM); m++)
                {
                    if (m != curr->axis)
                    {
                        r   = xc[ix][m] - curr->x0[nx0][ix][m] - dcom[m];
                        r2 += r*r;
                    }
                }
                break;
            default:
                gmx_fatal(FARGS, "Error: did not expect option value %d", curr->type);
        }
        g += r2;
    }
    g /= nx;
    msmul(mat, 1.0/nx, mat);

    return g;
}

/* calculate the com of molecules in x and put it into xa */
static void calc_mol_com(int nmol, int *molindex, const t_block *mols, const t_atoms *atoms,
                         rvec *x, rvec *xa)
{
    int  m, mol, i, d;
    rvec xm;
    real mass, mtot;

    for (m = 0; m < nmol; m++)
    {
        mol = molindex[m];
        clear_rvec(xm);
        mtot = 0;
        for (i = mols->index[mol]; i < mols->index[mol+1]; i++)
        {
            mass = atoms->atom[i].m;
            for (d = 0; d < DIM; d++)
            {
                xm[d] += mass*x[i][d];
            }
            mtot += mass;
        }
        svmul(1/mtot, xm, xa[m]);
    }
}

static real calc_one_mw(t_corr *curr, int ix, int nx0, rvec xc[], real *tm,
                        rvec dcom, gmx_bool bTen, matrix mat)
{
    real r2, r, mm;
    rvec rv;
    int  m, m2;

    mm = curr->mass[ix];
    if (mm == 0)
    {
        return 0;
    }
    (*tm) += mm;
    r2     = 0.0;
    switch (curr->type)
    {
        case NORMAL:
            for (m = 0; (m < DIM); m++)
            {
                rv[m] = xc[ix][m] - curr->x0[nx0][ix][m] - dcom[m];
                r2   += mm*rv[m]*rv[m];
                if (bTen)
                {
                    for (m2 = 0; m2 <= m; m2++)
                    {
                        mat[m][m2] += mm*rv[m]*rv[m2];
                    }
                }
            }
            break;
        case X:
        case Y:
        case Z:
            r  = xc[ix][curr->type-X] - curr->x0[nx0][ix][curr->type-X] -
                dcom[curr->type-X];
            r2 = mm*r*r;
            break;
        case LATERAL:
            for (m = 0; (m < DIM); m++)
            {
                if (m != curr->axis)
                {
                    r   = xc[ix][m] - curr->x0[nx0][ix][m] - dcom[m];
                    r2 += mm*r*r;
                }
            }
            break;
        default:
            gmx_fatal(FARGS, "Options got screwed. Did not expect value %d\n", curr->type);
    } /* end switch */
    return r2;
}

/* the normal, mass-weighted mean-squared displacement calcuation */
static real calc1_mw(t_corr *curr, int nx, int index[], int nx0, rvec xc[],
                     rvec dcom, gmx_bool bTen, matrix mat)
{
    int  i;
    real g, tm;

    g = tm = 0.0;
    clear_mat(mat);
    for (i = 0; (i < nx); i++)
    {
        g += calc_one_mw(curr, index[i], nx0, xc, &tm, dcom, bTen, mat);
    }

    g /= tm;
    if (bTen)
    {
        msmul(mat, 1/tm, mat);
    }

    return g;
}

/* prepare the coordinates by removing periodic boundary crossings.
   gnx = the number of atoms/molecules
   index = the indices
   xcur = the current coordinates
   xprev = the previous coordinates
   box = the box matrix */
static void prep_data(gmx_bool bMol, int gnx, int index[],
                      rvec xcur[], rvec xprev[], matrix box)
{
    int  i, m, ind;
    rvec hbox;

    /* Remove periodicity */
    for (m = 0; (m < DIM); m++)
    {
        hbox[m] = 0.5*box[m][m];
    }

    for (i = 0; (i < gnx); i++)
    {
        if (bMol)
        {
            ind = i;
        }
        else
        {
            ind = index[i];
        }

        for (m = DIM-1; m >= 0; m--)
        {
            if (hbox[m] == 0)
            {
                continue;
            }
            while (xcur[ind][m]-xprev[ind][m] <= -hbox[m])
            {
                rvec_inc(xcur[ind], box[m]);
            }
            while (xcur[ind][m]-xprev[ind][m] >  hbox[m])
            {
                rvec_dec(xcur[ind], box[m]);
            }
        }
    }
}

/* calculate the center of mass for a group
   gnx = the number of atoms/molecules
   index = the indices
   xcur = the current coordinates
   xprev = the previous coordinates
   box = the box matrix
   atoms = atom data (for mass)
   com(output) = center of mass  */
static void calc_com(gmx_bool bMol, int gnx, int index[],
                     rvec xcur[], rvec xprev[], matrix box, const t_atoms *atoms,
                     rvec com)
{
    int    i, m, ind;
    real   mass;
    double tmass;
    dvec   sx;

    clear_dvec(sx);
    tmass = 0;

    prep_data(bMol, gnx, index, xcur, xprev, box);
    for (i = 0; (i < gnx); i++)
    {
        if (bMol)
        {
            ind = i;
        }
        else
        {
            ind = index[i];
        }


        mass = atoms->atom[ind].m;
        for (m = 0; m < DIM; m++)
        {
            sx[m] += mass*xcur[ind][m];
        }
        tmass += mass;
    }
    for (m = 0; m < DIM; m++)
    {
        com[m] = sx[m]/tmass;
    }
}


static real calc1_mol(t_corr *curr, int nx, int gmx_unused index[], int nx0, rvec xc[],
                      rvec dcom, gmx_bool bTen, matrix mat)
{
    int  i;
    real g, tm, gtot, tt;

    tt   = curr->time[in_data(curr, nx0)];
    gtot = 0;
    tm   = 0;
    clear_mat(mat);
    for (i = 0; (i < nx); i++)
    {
        g = calc_one_mw(curr, i, nx0, xc, &tm, dcom, bTen, mat);
        /* We don't need to normalize as the mass was set to 1 */
        gtot += g;
        if (tt >= curr->beginfit && (curr->endfit < 0 || tt <= curr->endfit))
        {
            gmx_stats_add_point(curr->lsq[nx0][i], tt, g, 0, 0);
        }
    }
    msmul(mat, 1.0/nx, mat);

    return gtot/nx;
}

static void printmol(t_corr *curr, const char *fn,
                     const char *fn_pdb, int *molindex, const t_topology *top,
                     rvec *x, int ePBC, matrix box, const gmx_output_env_t *oenv)
{
#define NDIST 100
    FILE       *out;
    gmx_stats_t lsq1;
    int         i, j;
    real        a, b, D, Dav, D2av, VarD, sqrtD, sqrtD_max, scale;
    t_pdbinfo  *pdbinfo = nullptr;
    const int  *mol2a   = nullptr;

    out = xvgropen(fn, "Diffusion Coefficients / Molecule", "Molecule", "D", oenv);

    if (fn_pdb)
    {
        pdbinfo = top->atoms.pdbinfo;
        mol2a   = top->mols.index;
    }

    Dav       = D2av = 0;
    sqrtD_max = 0;
    for (i = 0; (i < curr->nmol); i++)
    {
        lsq1 = gmx_stats_init();
        for (j = 0; (j < curr->nrestart); j++)
        {
            real xx, yy, dx, dy;

            while (gmx_stats_get_point(curr->lsq[j][i], &xx, &yy, &dx, &dy, 0) == estatsOK)
            {
                gmx_stats_add_point(lsq1, xx, yy, dx, dy);
            }
        }
        gmx_stats_get_ab(lsq1, elsqWEIGHT_NONE, &a, &b, nullptr, nullptr, nullptr, nullptr);
        gmx_stats_free(lsq1);
        D     = a*FACTOR/curr->dim_factor;
        if (D < 0)
        {
            D   = 0;
        }
        Dav  += D;
        D2av += gmx::square(D);
        fprintf(out, "%10d  %10g\n", i, D);
        if (pdbinfo)
        {
            sqrtD = std::sqrt(D);
            if (sqrtD > sqrtD_max)
            {
                sqrtD_max = sqrtD;
            }
            for (j = mol2a[molindex[i]]; j < mol2a[molindex[i]+1]; j++)
            {
                pdbinfo[j].bfac = sqrtD;
            }
        }
    }
    xvgrclose(out);
    do_view(oenv, fn, "-graphtype bar");

    /* Compute variance, stddev and error */
    Dav  /= curr->nmol;
    D2av /= curr->nmol;
    VarD  = D2av - gmx::square(Dav);
    printf("<D> = %.4f Std. Dev. = %.4f Error = %.4f\n",
           Dav, std::sqrt(VarD), std::sqrt(VarD/curr->nmol));

    if (fn_pdb && x)
    {
        scale = 1;
        while (scale*sqrtD_max > 10)
        {
            scale *= 0.1;
        }
        while (scale*sqrtD_max < 0.1)
        {
            scale *= 10;
        }
        GMX_RELEASE_ASSERT(pdbinfo != nullptr, "Internal error - pdbinfo not set for PDB input");
        for (i = 0; i < top->atoms.nr; i++)
        {
            pdbinfo[i].bfac *= scale;
        }
        write_sto_conf(fn_pdb, "molecular MSD", &top->atoms, x, nullptr, ePBC, box);
    }
}

/* this is the main loop for the correlation type functions
 * fx and nx are file pointers to things like read_first_x and
 * read_next_x
 */
int corr_loop(t_corr *curr, const char *fn, const t_topology *top, int ePBC,
              gmx_bool bMol, int gnx[], int *index[],
              t_calc_func *calc1, gmx_bool bTen, int *gnx_com, int *index_com[],
              real dt, real t_pdb, rvec **x_pdb, matrix box_pdb,
              const gmx_output_env_t *oenv)
{
    rvec            *x[2];  /* the coordinates to read */
    rvec            *xa[2]; /* the coordinates to calculate displacements for */
    rvec             com = {0};
    real             t, t_prev = 0;
    int              natoms, i, j, cur = 0, maxframes = 0;
    t_trxstatus     *status;
#define        prev (1-cur)
    matrix           box;
    gmx_bool         bFirst;
    gmx_rmpbc_t      gpbc = nullptr;

    natoms = read_first_x(oenv, &status, fn, &curr->t0, &(x[cur]), box);
#ifdef DEBUG
    fprintf(stderr, "Read %d atoms for first frame\n", natoms);
#endif
    if ((gnx_com != nullptr) && natoms < top->atoms.nr)
    {
        fprintf(stderr, "WARNING: The trajectory only contains part of the system (%d of %d atoms) and therefore the COM motion of only this part of the system will be removed\n", natoms, top->atoms.nr);
    }

    snew(x[prev], natoms);

    if (bMol)
    {
        curr->ncoords = curr->nmol;
        snew(xa[0], curr->ncoords);
        snew(xa[1], curr->ncoords);
    }
    else
    {
        curr->ncoords = natoms;
        xa[0]         = x[0];
        xa[1]         = x[1];
    }

    bFirst = TRUE;
    t      = curr->t0;
    if (x_pdb)
    {
        *x_pdb = nullptr;
    }

    if (bMol)
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }

    /* the loop over all frames */
    do
    {
        if (x_pdb && ((bFirst && t_pdb < t) ||
                      (!bFirst &&
                       t_pdb > t - 0.5*(t - t_prev) &&
                       t_pdb < t + 0.5*(t - t_prev))))
        {
            if (*x_pdb == nullptr)
            {
                snew(*x_pdb, natoms);
            }
            for (i = 0; i < natoms; i++)
            {
                copy_rvec(x[cur][i], (*x_pdb)[i]);
            }
            copy_mat(box, box_pdb);
        }


        /* check whether we've reached a restart point */
        if (bRmod(t, curr->t0, dt))
        {
            curr->nrestart++;

            srenew(curr->x0, curr->nrestart);
            snew(curr->x0[curr->nrestart-1], curr->ncoords);
            srenew(curr->com, curr->nrestart);
            srenew(curr->n_offs, curr->nrestart);
            srenew(curr->lsq, curr->nrestart);
            snew(curr->lsq[curr->nrestart-1], curr->nmol);
            for (i = 0; i < curr->nmol; i++)
            {
                curr->lsq[curr->nrestart-1][i]  = gmx_stats_init();
            }

            if (debug)
            {
                fprintf(debug, "Extended data structures because of new restart %d\n",
                        curr->nrestart);
            }
        }
        /* create or extend the frame-based arrays */
        if (curr->nframes >= maxframes-1)
        {
            if (maxframes == 0)
            {
                for (i = 0; (i < curr->ngrp); i++)
                {
                    curr->ndata[i] = nullptr;
                    curr->data[i]  = nullptr;
                    if (bTen)
                    {
                        curr->datam[i] = nullptr;
                    }
                }
                curr->time = nullptr;
            }
            maxframes += 10;
            for (i = 0; (i < curr->ngrp); i++)
            {
                srenew(curr->ndata[i], maxframes);
                srenew(curr->data[i], maxframes);
                if (bTen)
                {
                    srenew(curr->datam[i], maxframes);
                }
                for (j = maxframes-10; j < maxframes; j++)
                {
                    curr->ndata[i][j] = 0;
                    curr->data[i][j]  = 0;
                    if (bTen)
                    {
                        clear_mat(curr->datam[i][j]);
                    }
                }
            }
            srenew(curr->time, maxframes);
        }

        /* set the time */
        curr->time[curr->nframes] = t - curr->t0;

        /* for the first frame, the previous frame is a copy of the first frame */
        if (bFirst)
        {
            std::memcpy(xa[prev], xa[cur], curr->ncoords*sizeof(xa[prev][0]));
            bFirst = FALSE;
        }

        /* make the molecules whole */
        if (bMol)
        {
            gmx_rmpbc(gpbc, natoms, box, x[cur]);
        }

        /* calculate the molecules' centers of masses and put them into xa */
        if (bMol)
        {
            calc_mol_com(gnx[0], index[0], &top->mols, &top->atoms, x[cur], xa[cur]);
        }

        /* first remove the periodic boundary condition crossings */
        for (i = 0; i < curr->ngrp; i++)
        {
            prep_data(bMol, gnx[i], index[i], xa[cur], xa[prev], box);
        }

        /* calculate the center of mass */
        if (gnx_com)
        {
            prep_data(bMol, gnx_com[0], index_com[0], xa[cur], xa[prev], box);
            calc_com(bMol, gnx_com[0], index_com[0], xa[cur], xa[prev], box,
                     &top->atoms, com);
        }

        /* loop over all groups in index file */
        for (i = 0; (i < curr->ngrp); i++)
        {
            /* calculate something useful, like mean square displacements */
            calc_corr(curr, i, gnx[i], index[i], xa[cur], (gnx_com != nullptr), com,
                      calc1, bTen);
        }
        cur    = prev;
        t_prev = t;

        curr->nframes++;
    }
    while (read_next_x(oenv, status, &t, x[cur], box));
    fprintf(stderr, "\nUsed %d restart points spaced %g %s over %g %s\n\n",
            curr->nrestart,
            output_env_conv_time(oenv, dt), output_env_get_time_unit(oenv).c_str(),
            output_env_conv_time(oenv, curr->time[curr->nframes-1]),
            output_env_get_time_unit(oenv).c_str() );

    if (bMol)
    {
        gmx_rmpbc_done(gpbc);
    }

    close_trj(status);

    return natoms;
}

static void index_atom2mol(int *n, int *index, const t_block *mols)
{
    int nat, i, nmol, mol, j;

    nat  = *n;
    i    = 0;
    nmol = 0;
    mol  = 0;
    while (i < nat)
    {
        while (index[i] > mols->index[mol])
        {
            mol++;
            if (mol >= mols->nr)
            {
                gmx_fatal(FARGS, "Atom index out of range: %d", index[i]+1);
            }
        }
        for (j = mols->index[mol]; j < mols->index[mol+1]; j++)
        {
            if (i >= nat || index[i] != j)
            {
                gmx_fatal(FARGS, "The index group does not consist of whole molecules");
            }
            i++;
        }
        index[nmol++] = mol;
    }

    fprintf(stderr, "Split group of %d atoms into %d molecules\n", nat, nmol);

    *n = nmol;
}

void do_corr(const char *trx_file, const char *ndx_file, const char *msd_file,
             const char *mol_file, const char *pdb_file, real t_pdb,
             int nrgrp, t_topology *top, int ePBC,
             gmx_bool bTen, gmx_bool bMW, gmx_bool bRmCOMM,
             int type, real dim_factor, int axis,
             real dt, real beginfit, real endfit, const gmx_output_env_t *oenv)
{
    t_corr        *msd;
    int           *gnx;   /* the selected groups' sizes */
    int          **index; /* selected groups' indices */
    char         **grpname;
    int            i, i0, i1, j, N, nat_trx;
    real          *DD, *SigmaD, a, a2, b, r, chi2;
    rvec          *x;
    matrix         box;
    int           *gnx_com     = nullptr; /* the COM removal group size  */
    int          **index_com   = nullptr; /* the COM removal group atom indices */
    char         **grpname_com = nullptr; /* the COM removal group name */

    snew(gnx, nrgrp);
    snew(index, nrgrp);
    snew(grpname, nrgrp);

    fprintf(stderr, "\nSelect a group to calculate mean squared displacement for:\n");
    get_index(&top->atoms, ndx_file, nrgrp, gnx, index, grpname);

    if (bRmCOMM)
    {
        snew(gnx_com, 1);
        snew(index_com, 1);
        snew(grpname_com, 1);

        fprintf(stderr, "\nNow select a group for center of mass removal:\n");
        get_index(&top->atoms, ndx_file, 1, gnx_com, index_com, grpname_com);
    }

    if (mol_file)
    {
        index_atom2mol(&gnx[0], index[0], &top->mols);
    }

    msd = init_corr(nrgrp, type, axis, dim_factor,
                    mol_file == nullptr ? 0 : gnx[0], bTen, bMW, dt, top,
                    beginfit, endfit);

    nat_trx =
        corr_loop(msd, trx_file, top, ePBC, mol_file ? gnx[0] : 0, gnx, index,
                  (mol_file != nullptr) ? calc1_mol : (bMW ? calc1_mw : calc1_norm),
                  bTen, gnx_com, index_com, dt, t_pdb,
                  pdb_file ? &x : nullptr, box, oenv);

    /* Correct for the number of points */
    for (j = 0; (j < msd->ngrp); j++)
    {
        for (i = 0; (i < msd->nframes); i++)
        {
            msd->data[j][i] /= msd->ndata[j][i];
            if (bTen)
            {
                msmul(msd->datam[j][i], 1.0/msd->ndata[j][i], msd->datam[j][i]);
            }
        }
    }

    if (mol_file)
    {
        if (pdb_file && x == nullptr)
        {
            fprintf(stderr, "\nNo frame found need time tpdb = %g ps\n"
                    "Can not write %s\n\n", t_pdb, pdb_file);
        }
        i             = top->atoms.nr;
        top->atoms.nr = nat_trx;
        if (pdb_file && top->atoms.pdbinfo == nullptr)
        {
            snew(top->atoms.pdbinfo, top->atoms.nr);
        }
        printmol(msd, mol_file, pdb_file, index[0], top, x, ePBC, box, oenv);
        top->atoms.nr = i;
    }

    DD     = nullptr;
    SigmaD = nullptr;

    if (beginfit == -1)
    {
        i0       = static_cast<int>(0.1*(msd->nframes - 1) + 0.5);
        beginfit = msd->time[i0];
    }
    else
    {
        for (i0 = 0; i0 < msd->nframes && msd->time[i0] < beginfit; i0++)
        {
            ;
        }
    }

    if (endfit == -1)
    {
        i1     = static_cast<int>(0.9*(msd->nframes - 1) + 0.5) + 1;
        endfit = msd->time[i1-1];
    }
    else
    {
        for (i1 = i0; i1 < msd->nframes && msd->time[i1] <= endfit; i1++)
        {
            ;
        }
    }
    fprintf(stdout, "Fitting from %g to %g %s\n\n", beginfit, endfit,
            output_env_get_time_unit(oenv).c_str());

    N = i1-i0;
    if (N <= 2)
    {
        fprintf(stdout, "Not enough points for fitting (%d).\n"
                "Can not determine the diffusion constant.\n", N);
    }
    else
    {
        snew(DD, msd->ngrp);
        snew(SigmaD, msd->ngrp);
        for (j = 0; j < msd->ngrp; j++)
        {
            if (N >= 4)
            {
                lsq_y_ax_b(N/2, &(msd->time[i0]), &(msd->data[j][i0]), &a, &b, &r, &chi2);
                lsq_y_ax_b(N/2, &(msd->time[i0+N/2]), &(msd->data[j][i0+N/2]), &a2, &b, &r, &chi2);
                SigmaD[j] = std::abs(a-a2);
            }
            else
            {
                SigmaD[j] = 0;
            }
            lsq_y_ax_b(N, &(msd->time[i0]), &(msd->data[j][i0]), &(DD[j]), &b, &r, &chi2);
            DD[j]     *= FACTOR/msd->dim_factor;
            SigmaD[j] *= FACTOR/msd->dim_factor;
            if (DD[j] > 0.01 && DD[j] < 1e4)
            {
                fprintf(stdout, "D[%10s] %.4f (+/- %.4f) 1e-5 cm^2/s\n",
                        grpname[j], DD[j], SigmaD[j]);
            }
            else
            {
                fprintf(stdout, "D[%10s] %.4g (+/- %.4g) 1e-5 cm^2/s\n",
                        grpname[j], DD[j], SigmaD[j]);
            }
        }
    }
    /* Print mean square displacement */
    corr_print(msd, bTen, msd_file,
               "Mean Square Displacement",
               "MSD (nm\\S2\\N)",
               msd->time[msd->nframes-1], beginfit, endfit, DD, SigmaD, grpname, oenv);
}

int gmx_msd(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] computes the mean square displacement (MSD) of atoms from",
        "a set of initial positions. This provides an easy way to compute",
        "the diffusion constant using the Einstein relation.",
        "The time between the reference points for the MSD calculation",
        "is set with [TT]-trestart[tt].",
        "The diffusion constant is calculated by least squares fitting a",
        "straight line (D*t + c) through the MSD(t) from [TT]-beginfit[tt] to",
        "[TT]-endfit[tt] (note that t is time from the reference positions,",
        "not simulation time). An error estimate given, which is the difference",
        "of the diffusion coefficients obtained from fits over the two halves",
        "of the fit interval.[PAR]",
        "There are three, mutually exclusive, options to determine different",
        "types of mean square displacement: [TT]-type[tt], [TT]-lateral[tt]",
        "and [TT]-ten[tt]. Option [TT]-ten[tt] writes the full MSD tensor for",
        "each group, the order in the output is: trace xx yy zz yx zx zy.[PAR]",
        "If [TT]-mol[tt] is set, [THISMODULE] plots the MSD for individual molecules",
        "(including making molecules whole across periodic boundaries): ",
        "for each individual molecule a diffusion constant is computed for ",
        "its center of mass. The chosen index group will be split into ",
        "molecules.[PAR]",
        "The default way to calculate a MSD is by using mass-weighted averages.",
        "This can be turned off with [TT]-nomw[tt].[PAR]",
        "With the option [TT]-rmcomm[tt], the center of mass motion of a ",
        "specific group can be removed. For trajectories produced with ",
        "GROMACS this is usually not necessary, ",
        "as [gmx-mdrun] usually already removes the center of mass motion.",
        "When you use this option be sure that the whole system is stored",
        "in the trajectory file.[PAR]",
        "The diffusion coefficient is determined by linear regression of the MSD,",
        "where, unlike for the normal output of D, the times are weighted",
        "according to the number of reference points, i.e. short times have",
        "a higher weight. Also when [TT]-beginfit[tt]=-1,fitting starts at 10%",
        "and when [TT]-endfit[tt]=-1, fitting goes to 90%.",
        "Using this option one also gets an accurate error estimate",
        "based on the statistics between individual molecules.",
        "Note that this diffusion coefficient and error estimate are only",
        "accurate when the MSD is completely linear between",
        "[TT]-beginfit[tt] and [TT]-endfit[tt].[PAR]",
        "Option [TT]-pdb[tt] writes a [REF].pdb[ref] file with the coordinates of the frame",
        "at time [TT]-tpdb[tt] with in the B-factor field the square root of",
        "the diffusion coefficient of the molecule.",
        "This option implies option [TT]-mol[tt]."
    };
    static const char *normtype[] = { nullptr, "no", "x", "y", "z", nullptr };
    static const char *axtitle[]  = { nullptr, "no", "x", "y", "z", nullptr };
    static int         ngroup     = 1;
    static real        dt         = 10;
    static real        t_pdb      = 0;
    static real        beginfit   = -1;
    static real        endfit     = -1;
    static gmx_bool    bTen       = FALSE;
    static gmx_bool    bMW        = TRUE;
    static gmx_bool    bRmCOMM    = FALSE;
    t_pargs            pa[]       = {
        { "-type",    FALSE, etENUM, {normtype},
          "Compute diffusion coefficient in one direction" },
        { "-lateral", FALSE, etENUM, {axtitle},
          "Calculate the lateral diffusion in a plane perpendicular to" },
        { "-ten",      FALSE, etBOOL, {&bTen},
          "Calculate the full tensor" },
        { "-ngroup",  FALSE, etINT,  {&ngroup},
          "Number of groups to calculate MSD for" },
        { "-mw",      FALSE, etBOOL, {&bMW},
          "Mass weighted MSD" },
        { "-rmcomm",      FALSE, etBOOL, {&bRmCOMM},
          "Remove center of mass motion" },
        { "-tpdb", FALSE, etTIME, {&t_pdb},
          "The frame to use for option [TT]-pdb[tt] (%t)" },
        { "-trestart", FALSE, etTIME, {&dt},
          "Time between restarting points in trajectory (%t)" },
        { "-beginfit", FALSE, etTIME, {&beginfit},
          "Start time for fitting the MSD (%t), -1 is 10%" },
        { "-endfit", FALSE, etTIME, {&endfit},
          "End time for fitting the MSD (%t), -1 is 90%" }
    };

    t_filenm           fnm[] = {
        { efTRX, nullptr, nullptr,  ffREAD },
        { efTPS, nullptr, nullptr,  ffREAD },
        { efNDX, nullptr, nullptr,  ffOPTRD },
        { efXVG, nullptr, "msd", ffWRITE },
        { efXVG, "-mol", "diff_mol", ffOPTWR },
        { efPDB, "-pdb", "diff_mol", ffOPTWR }
    };
#define NFILE asize(fnm)

    t_topology        top;
    int               ePBC;
    matrix            box;
    const char       *trx_file, *tps_file, *ndx_file, *msd_file, *mol_file, *pdb_file;
    rvec             *xdum;
    gmx_bool          bTop;
    int               axis, type;
    real              dim_factor;
    gmx_output_env_t *oenv;

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    trx_file = ftp2fn_null(efTRX, NFILE, fnm);
    tps_file = ftp2fn_null(efTPS, NFILE, fnm);
    ndx_file = ftp2fn_null(efNDX, NFILE, fnm);
    msd_file = ftp2fn_null(efXVG, NFILE, fnm);
    pdb_file = opt2fn_null("-pdb", NFILE, fnm);
    if (pdb_file)
    {
        mol_file = opt2fn("-mol", NFILE, fnm);
    }
    else
    {
        mol_file = opt2fn_null("-mol", NFILE, fnm);
    }

    if (ngroup < 1)
    {
        gmx_fatal(FARGS, "Must have at least 1 group (now %d)", ngroup);
    }
    if (mol_file && ngroup > 1)
    {
        gmx_fatal(FARGS, "With molecular msd can only have 1 group (now %d)",
                  ngroup);
    }


    if (mol_file)
    {
        bMW  = TRUE;
        fprintf(stderr, "Calculating diffusion coefficients for molecules.\n");
    }

    GMX_RELEASE_ASSERT(normtype[0] != nullptr, "Options inconsistency; normtype[0] is NULL");
    GMX_RELEASE_ASSERT(axtitle[0] != nullptr, "Options inconsistency; axtitle[0] is NULL");

    if (normtype[0][0] != 'n')
    {
        type       = normtype[0][0] - 'x' + X; /* See defines above */
        dim_factor = 2.0;
    }
    else
    {
        type       = NORMAL;
        dim_factor = 6.0;
    }
    if ((type == NORMAL) && (axtitle[0][0] != 'n'))
    {
        type       = LATERAL;
        dim_factor = 4.0;
        axis       = (axtitle[0][0] - 'x'); /* See defines above */
    }
    else
    {
        axis = 0;
    }

    if (bTen && type != NORMAL)
    {
        gmx_fatal(FARGS, "Can only calculate the full tensor for 3D msd");
    }

    bTop = read_tps_conf(tps_file, &top, &ePBC, &xdum, nullptr, box, bMW || bRmCOMM);
    if (mol_file && !bTop)
    {
        gmx_fatal(FARGS,
                  "Could not read a topology from %s. Try a tpr file instead.",
                  tps_file);
    }

    do_corr(trx_file, ndx_file, msd_file, mol_file, pdb_file, t_pdb, ngroup,
            &top, ePBC, bTen, bMW, bRmCOMM, type, dim_factor, axis, dt, beginfit, endfit,
            oenv);

    view_all(oenv, NFILE, fnm);

    return 0;
}
