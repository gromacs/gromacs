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

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define e2d(x) ENM2DEBYE*(x)
#define EANG2CM  E_CHARGE*1.0e-10       /* e Angstrom to Coulomb meter */
#define CM2D  SPEED_OF_LIGHT*1.0e+24    /* Coulomb meter to Debye */

typedef struct {
    int      nelem;
    real     spacing, radius;
    real    *elem;
    int     *count;
    gmx_bool bPhi;
    int      nx, ny;
    real   **cmap;
} t_gkrbin;

static t_gkrbin *mk_gkrbin(real radius, real rcmax, gmx_bool bPhi, int ndegrees)
{
    t_gkrbin *gb;
    char     *ptr;
    int       i;

    snew(gb, 1);

    if ((ptr = getenv("GMX_DIPOLE_SPACING")) != nullptr)
    {
        double bw = strtod(ptr, nullptr);
        gb->spacing = bw;
    }
    else
    {
        gb->spacing = 0.01; /* nm */
    }
    gb->nelem   = 1 + static_cast<int>(radius/gb->spacing);
    if (rcmax == 0)
    {
        gb->nx = gb->nelem;
    }
    else
    {
        gb->nx = 1 + static_cast<int>(rcmax/gb->spacing);
    }
    gb->radius  = radius;
    snew(gb->elem, gb->nelem);
    snew(gb->count, gb->nelem);

    snew(gb->cmap, gb->nx);
    gb->ny = std::max(2, ndegrees);
    for (i = 0; (i < gb->nx); i++)
    {
        snew(gb->cmap[i], gb->ny);
    }
    gb->bPhi = bPhi;

    return gb;
}

static void done_gkrbin(t_gkrbin **gb)
{
    sfree((*gb)->elem);
    sfree((*gb)->count);
    sfree((*gb));
    *gb = nullptr;
}

static void add2gkr(t_gkrbin *gb, real r, real cosa, real phi)
{
    int  cy, index = std::round(r/gb->spacing);
    real alpha;

    if (index < gb->nelem)
    {
        gb->elem[index]  += cosa;
        gb->count[index]++;
    }
    if (index < gb->nx)
    {
        alpha = std::acos(cosa);
        if (gb->bPhi)
        {
            cy = static_cast<int>((M_PI+phi)*gb->ny/(2*M_PI));
        }
        else
        {
            cy = static_cast<int>((alpha*gb->ny)/M_PI); /*((1+cosa)*0.5*(gb->ny));*/
        }
        cy = std::min(gb->ny-1, std::max(0, cy));
        if (debug)
        {
            fprintf(debug, "CY: %10f  %5d\n", alpha, cy);
        }
        gb->cmap[index][cy] += 1;
    }
}

static void rvec2sprvec(rvec dipcart, rvec dipsp)
{
    clear_rvec(dipsp);
    dipsp[0] = std::sqrt(dipcart[XX]*dipcart[XX]+dipcart[YY]*dipcart[YY]+dipcart[ZZ]*dipcart[ZZ]);  /* R */
    dipsp[1] = std::atan2(dipcart[YY], dipcart[XX]);                                                /* Theta */
    dipsp[2] = std::atan2(std::sqrt(dipcart[XX]*dipcart[XX]+dipcart[YY]*dipcart[YY]), dipcart[ZZ]); /* Phi */
}



static void do_gkr(t_gkrbin *gb, int ncos, int *ngrp, int *molindex[],
                   int mindex[], rvec x[], rvec mu[],
                   int ePBC, const matrix box, const t_atom *atom, const int *nAtom)
{
    static rvec *xcm[2] = { nullptr, nullptr};
    int          gi, gj, j0, j1, i, j, k, n, grp0, grp1;
    real         qtot, q, cosa, rr, phi;
    rvec         dx;
    t_pbc        pbc;

    for (n = 0; (n < ncos); n++)
    {
        if (!xcm[n])
        {
            snew(xcm[n], ngrp[n]);
        }
        for (i = 0; (i < ngrp[n]); i++)
        {
            /* Calculate center of mass of molecule */
            gi = molindex[n][i];
            j0 = mindex[gi];

            if (nAtom[n] > 0)
            {
                copy_rvec(x[j0+nAtom[n]-1], xcm[n][i]);
            }
            else
            {
                j1 = mindex[gi+1];
                clear_rvec(xcm[n][i]);
                qtot = 0;
                for (j = j0; j < j1; j++)
                {
                    q     = std::abs(atom[j].q);
                    qtot += q;
                    for (k = 0; k < DIM; k++)
                    {
                        xcm[n][i][k] += q*x[j][k];
                    }
                }
                svmul(1/qtot, xcm[n][i], xcm[n][i]);
            }
        }
    }
    set_pbc(&pbc, ePBC, box);
    grp0 = 0;
    grp1 = ncos-1;
    for (i = 0; i < ngrp[grp0]; i++)
    {
        gi = molindex[grp0][i];
        j0 = (ncos == 2) ? 0 : i+1;
        for (j = j0; j < ngrp[grp1]; j++)
        {
            gj   = molindex[grp1][j];
            if ((iprod(mu[gi], mu[gi]) > 0) && (iprod(mu[gj], mu[gj]) > 0))
            {
                /* Compute distance between molecules including PBC */
                pbc_dx(&pbc, xcm[grp0][i], xcm[grp1][j], dx);
                rr = norm(dx);

                if (gb->bPhi)
                {
                    rvec xi, xj, xk, xl;
                    rvec r_ij, r_kj, r_kl, mm, nn;
                    int  t1, t2, t3;

                    copy_rvec(xcm[grp0][i], xj);
                    copy_rvec(xcm[grp1][j], xk);
                    rvec_add(xj, mu[gi], xi);
                    rvec_add(xk, mu[gj], xl);
                    phi = dih_angle(xi, xj, xk, xl, &pbc,
                                    r_ij, r_kj, r_kl, mm, nn, /* out */
                                    &t1, &t2, &t3);
                    cosa = std::cos(phi);
                }
                else
                {
                    cosa = cos_angle(mu[gi], mu[gj]);
                    phi  = 0;
                }
                if (debug || std::isnan(cosa))
                {
                    fprintf(debug ? debug : stderr,
                            "mu[%d] = %5.2f %5.2f %5.2f |mi| = %5.2f, mu[%d] = %5.2f %5.2f %5.2f |mj| = %5.2f rr = %5.2f cosa = %5.2f\n",
                            gi, mu[gi][XX], mu[gi][YY], mu[gi][ZZ], norm(mu[gi]),
                            gj, mu[gj][XX], mu[gj][YY], mu[gj][ZZ], norm(mu[gj]),
                            rr, cosa);
                }

                add2gkr(gb, rr, cosa, phi);
            }
        }
    }
}

static real normalize_cmap(t_gkrbin *gb)
{
    int    i, j;
    real   hi;
    double vol;

    hi = 0;
    for (i = 0; (i < gb->nx); i++)
    {
        vol = 4*M_PI*gmx::square(gb->spacing*i)*gb->spacing;
        for (j = 0; (j < gb->ny); j++)
        {
            gb->cmap[i][j] /= vol;
            hi              = std::max(hi, gb->cmap[i][j]);
        }
    }
    if (hi <= 0)
    {
        gmx_fatal(FARGS, "No data in the cmap");
    }
    return hi;
}

static void print_cmap(const char *cmap, t_gkrbin *gb, int *nlevels)
{
    FILE   *out;
    int     i, j;
    real    hi;

    real   *xaxis, *yaxis;
    t_rgb   rlo = { 1, 1, 1 };
    t_rgb   rhi = { 0, 0, 0 };

    hi = normalize_cmap(gb);
    snew(xaxis, gb->nx+1);
    for (i = 0; (i < gb->nx+1); i++)
    {
        xaxis[i] = i*gb->spacing;
    }
    snew(yaxis, gb->ny);
    for (j = 0; (j < gb->ny); j++)
    {
        if (gb->bPhi)
        {
            yaxis[j] = (360.0*j)/(gb->ny-1.0)-180;
        }
        else
        {
            yaxis[j] = (180.0*j)/(gb->ny-1.0);
        }
        /*2.0*j/(gb->ny-1.0)-1.0;*/
    }
    out = gmx_ffopen(cmap, "w");
    write_xpm(out, 0,
              "Dipole Orientation Distribution", "Fraction", "r (nm)",
              gb->bPhi ? "Phi" : "Alpha",
              gb->nx, gb->ny, xaxis, yaxis,
              gb->cmap, 0, hi, rlo, rhi, nlevels);
    gmx_ffclose(out);
    sfree(xaxis);
    sfree(yaxis);
}

static void print_gkrbin(const char *fn, t_gkrbin *gb,
                         int ngrp, int nframes, real volume,
                         const gmx_output_env_t *oenv)
{
    /* We compute Gk(r), gOO and hOO according to
     * Nymand & Linse, JCP 112 (2000) pp 6386-6395.
     * In this implementation the angle between dipoles is stored
     * rather than their inner product. This allows to take polarizible
     * models into account. The RDF is calculated as well, almost for free!
     */
    FILE       *fp;
    const char *leg[] = { "G\\sk\\N(r)", "< cos >", "h\\sOO\\N", "g\\sOO\\N", "Energy" };
    int         i, last;
    real        x0, x1, ggg, Gkr, vol_s, rho, gOO, hOO, cosav, ener;
    double      fac;

    fp = xvgropen(fn, "Distance dependent Gk", "r (nm)", "G\\sk\\N(r)", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);

    Gkr = 1;  /* Self-dipole inproduct = 1 */
    rho = ngrp/volume;

    if (debug)
    {
        fprintf(debug, "Number density is %g molecules / nm^3\n", rho);
        fprintf(debug, "ngrp = %d, nframes = %d\n", ngrp, nframes);
    }

    last = gb->nelem-1;
    while (last > 1 && gb->elem[last-1] == 0)
    {
        last--;
    }

    /* Divide by dipole squared, by number of frames, by number of origins.
     * Multiply by 2 because we only take half the matrix of interactions
     * into account.
     */
    fac  = 2.0/(ngrp * nframes);

    x0 = 0;
    for (i = 0; i < last; i++)
    {
        /* Centre of the coordinate in the spherical layer */
        x1    = x0+gb->spacing;

        /* Volume of the layer */
        vol_s = (4.0/3.0)*M_PI*(x1*x1*x1-x0*x0*x0);

        /* gOO */
        gOO   = gb->count[i]*fac/(rho*vol_s);

        /* Dipole correlation hOO, normalized by the relative number density, like
         * in a Radial distribution function.
         */
        ggg  = gb->elem[i]*fac;
        hOO  = 3.0*ggg/(rho*vol_s);
        Gkr += ggg;
        if (gb->count[i])
        {
            cosav = gb->elem[i]/gb->count[i];
        }
        else
        {
            cosav = 0;
        }
        ener = -0.5*cosav*ONE_4PI_EPS0/(x1*x1*x1);

        fprintf(fp, "%10.5e %12.5e %12.5e %12.5e %12.5e  %12.5e\n",
                x1, Gkr, cosav, hOO, gOO, ener);

        /* Swap x0 and x1 */
        x0 = x1;
    }
    xvgrclose(fp);
}

static gmx_bool read_mu_from_enx(ener_file_t fmu, int Vol, ivec iMu, rvec mu, real *vol,
                                 real *t, int nre, t_enxframe *fr)
{
    int          i;
    gmx_bool     bCont;
    char         buf[22];

    bCont = do_enx(fmu, fr);
    if (fr->nre != nre)
    {
        fprintf(stderr, "Something strange: expected %d entries in energy file at step %s\n(time %g) but found %d entries\n",
                nre, gmx_step_str(fr->step, buf), fr->t, fr->nre);
    }

    if (bCont)
    {
        if (Vol != -1)          /* we've got Volume in the energy file */
        {
            *vol = fr->ener[Vol].e;
        }
        for (i = 0; i < DIM; i++)
        {
            mu[i] = fr->ener[iMu[i]].e;
        }
        *t = fr->t;
    }

    return bCont;
}

static void neutralize_mols(int n, int *index, const t_block *mols, t_atom *atom)
{
    double mtot, qtot;
    int    ncharged, m, a0, a1, a;

    ncharged = 0;
    for (m = 0; m < n; m++)
    {
        a0   = mols->index[index[m]];
        a1   = mols->index[index[m]+1];
        mtot = 0;
        qtot = 0;
        for (a = a0; a < a1; a++)
        {
            mtot += atom[a].m;
            qtot += atom[a].q;
        }
        /* This check is only for the count print */
        if (std::abs(qtot) > 0.01)
        {
            ncharged++;
        }
        if (mtot > 0)
        {
            /* Remove the net charge at the center of mass */
            for (a = a0; a < a1; a++)
            {
                atom[a].q -= qtot*atom[a].m/mtot;
            }
        }
    }

    if (ncharged)
    {
        printf("There are %d charged molecules in the selection,\n"
               "will subtract their charge at their center of mass\n", ncharged);
    }
}

static void mol_dip(int k0, int k1, rvec x[], const t_atom atom[], rvec mu)
{
    int  k, m;
    real q;

    clear_rvec(mu);
    for (k = k0; (k < k1); k++)
    {
        q  = e2d(atom[k].q);
        for (m = 0; (m < DIM); m++)
        {
            mu[m] += q*x[k][m];
        }
    }
}

static void mol_quad(int k0, int k1, rvec x[], const t_atom atom[], rvec quad)
{
    int      i, k, m, n, niter;
    real     q, r2, mass, masstot;
    rvec     com;        /* center of mass */
    rvec     r;          /* distance of atoms to center of mass */
    double **inten;
    double   dd[DIM], **ev, tmp;

    snew(inten, DIM);
    snew(ev, DIM);
    for (i = 0; (i < DIM); i++)
    {
        snew(inten[i], DIM);
        snew(ev[i], DIM);
        dd[i] = 0.0;
    }

    /* Compute center of mass */
    clear_rvec(com);
    masstot = 0.0;
    for (k = k0; (k < k1); k++)
    {
        mass     = atom[k].m;
        masstot += mass;
        for (i = 0; (i < DIM); i++)
        {
            com[i] += mass*x[k][i];
        }
    }
    svmul((1.0/masstot), com, com);

    /* We want traceless quadrupole moments, so let us calculate the complete
     * quadrupole moment tensor and diagonalize this tensor to get
     * the individual components on the diagonal.
     */

#define delta(a, b) (( a == b ) ? 1.0 : 0.0)

    for (m = 0; (m < DIM); m++)
    {
        for (n = 0; (n < DIM); n++)
        {
            inten[m][n] = 0;
        }
    }
    for (k = k0; (k < k1); k++)         /* loop over atoms in a molecule */
    {
        q  = (atom[k].q)*100.0;
        rvec_sub(x[k], com, r);
        r2 = iprod(r, r);
        for (m = 0; (m < DIM); m++)
        {
            for (n = 0; (n < DIM); n++)
            {
                inten[m][n] += 0.5*q*(3.0*r[m]*r[n] - r2*delta(m, n))*EANG2CM*CM2D;
            }
        }
    }
    if (debug)
    {
        for (i = 0; (i < DIM); i++)
        {
            fprintf(debug, "Q[%d] = %8.3f  %8.3f  %8.3f\n",
                    i, inten[i][XX], inten[i][YY], inten[i][ZZ]);
        }
    }

    /* We've got the quadrupole tensor, now diagonalize the sucker */
    jacobi(inten, 3, dd, ev, &niter);

    if (debug)
    {
        for (i = 0; (i < DIM); i++)
        {
            fprintf(debug, "ev[%d] = %8.3f  %8.3f  %8.3f\n",
                    i, ev[i][XX], ev[i][YY], ev[i][ZZ]);
        }
        for (i = 0; (i < DIM); i++)
        {
            fprintf(debug, "Q'[%d] = %8.3f  %8.3f  %8.3f\n",
                    i, inten[i][XX], inten[i][YY], inten[i][ZZ]);
        }
    }
    /* Sort the eigenvalues, for water we know that the order is as follows:
     *
     * Q_yy, Q_zz, Q_xx
     *
     * At the moment I have no idea how this will work out for other molecules...
     */

#define SWAP(i)                                 \
    if (dd[i+1] > dd[i]) {                      \
        tmp     = dd[i];                              \
        dd[i]   = dd[i+1];                          \
        dd[i+1] = tmp;                            \
    }
    SWAP(0);
    SWAP(1);
    SWAP(0);

    for (m = 0; (m < DIM); m++)
    {
        quad[0] = dd[2];  /* yy */
        quad[1] = dd[0];  /* zz */
        quad[2] = dd[1];  /* xx */
    }

    if (debug)
    {
        pr_rvec(debug, 0, "Quadrupole", quad, DIM, TRUE);
    }

    /* clean-up */
    for (i = 0; (i < DIM); i++)
    {
        sfree(inten[i]);
        sfree(ev[i]);
    }
    sfree(inten);
    sfree(ev);
}

/*
 * Calculates epsilon according to M. Neumann, Mol. Phys. 50, 841 (1983)
 */
static real calc_eps(double M_diff, double volume, double epsRF, double temp)
{
    double eps, A, teller, noemer;
    double eps_0 = 8.854187817e-12;   /* epsilon_0 in C^2 J^-1 m^-1 */
    double fac   = 1.112650021e-59;   /* converts Debye^2 to C^2 m^2 */

    A = M_diff*fac/(3*eps_0*volume*NANO*NANO*NANO*BOLTZMANN*temp);

    if (epsRF == 0.0)
    {
        teller = 1 + A;
        noemer = 1;
    }
    else
    {
        teller = 1 + (A*2*epsRF/(2*epsRF+1));
        noemer = 1 - (A/(2*epsRF+1));
    }
    eps = teller / noemer;

    return eps;
}

static void update_slab_dipoles(int k0, int k1, rvec x[], rvec mu,
                                int idim, int nslice, rvec slab_dipole[],
                                matrix box)
{
    int  k;
    real xdim = 0;

    for (k = k0; (k < k1); k++)
    {
        xdim += x[k][idim];
    }
    xdim /= (k1-k0);
    k     = (static_cast<int>(xdim*nslice/box[idim][idim] + nslice)) % nslice;
    rvec_inc(slab_dipole[k], mu);
}

static void dump_slab_dipoles(const char *fn, int idim, int nslice,
                              rvec slab_dipole[], matrix box, int nframes,
                              const gmx_output_env_t *oenv)
{
    FILE       *fp;
    char        buf[STRLEN];
    int         i;
    real        mutot;
    const char *leg_dim[4] = {
        "\\f{12}m\\f{4}\\sX\\N",
        "\\f{12}m\\f{4}\\sY\\N",
        "\\f{12}m\\f{4}\\sZ\\N",
        "\\f{12}m\\f{4}\\stot\\N"
    };

    sprintf(buf, "Box-%c (nm)", 'X'+idim);
    fp = xvgropen(fn, "Average dipole moment per slab", buf, "\\f{12}m\\f{4} (D)",
                  oenv);
    xvgr_legend(fp, DIM, leg_dim, oenv);
    for (i = 0; (i < nslice); i++)
    {
        mutot = norm(slab_dipole[i])/nframes;
        fprintf(fp, "%10.3f  %10.3f  %10.3f  %10.3f  %10.3f\n",
                ((i+0.5)*box[idim][idim])/nslice,
                slab_dipole[i][XX]/nframes,
                slab_dipole[i][YY]/nframes,
                slab_dipole[i][ZZ]/nframes,
                mutot);
    }
    xvgrclose(fp);
    do_view(oenv, fn, "-autoscale xy -nxy");
}

static void compute_avercos(int n, rvec dip[], real *dd, rvec axis, gmx_bool bPairs)
{
    int    i, j, k;
    double dc, d, ddc1, ddc2, ddc3;
    rvec   xxx = { 1, 0, 0 };
    rvec   yyy = { 0, 1, 0 };
    rvec   zzz = { 0, 0, 1 };

    d    = 0;
    ddc1 = ddc2 = ddc3 = 0;
    for (i = k = 0; (i < n); i++)
    {
        ddc1 += std::abs(cos_angle(dip[i], xxx));
        ddc2 += std::abs(cos_angle(dip[i], yyy));
        ddc3 += std::abs(cos_angle(dip[i], zzz));
        if (bPairs)
        {
            for (j = i+1; (j < n); j++, k++)
            {
                dc  = cos_angle(dip[i], dip[j]);
                d  += std::abs(dc);
            }
        }
    }
    *dd      = d/k;
    axis[XX] = ddc1/n;
    axis[YY] = ddc2/n;
    axis[ZZ] = ddc3/n;
}

static void do_dip(const t_topology *top, int ePBC, real volume,
                   const char *fn,
                   const char *out_mtot, const char *out_eps,
                   const char *out_aver, const char *dipdist,
                   const char *cosaver, const char *fndip3d,
                   const char *fnadip,  gmx_bool bPairs,
                   const char *corrtype, const char *corf,
                   gmx_bool bGkr,     const char *gkrfn,
                   gmx_bool bPhi,     int  *nlevels,  int ndegrees,
                   int  ncos,
                   const char *cmap,    real rcmax,
                   gmx_bool bQuad,
                   gmx_bool bMU,      const char *mufn,
                   int  *gnx,     int  *molindex[],
                   real mu_max,   real mu_aver,
                   real epsilonRF, real temp,
                   int  *gkatom,  int skip,
                   gmx_bool bSlab,    int nslices,
                   const char *axtitle, const char *slabfn,
                   const gmx_output_env_t *oenv)
{
    const char *leg_mtot[] = {
        "M\\sx \\N",
        "M\\sy \\N",
        "M\\sz \\N",
        "|M\\stot \\N|"
    };
#define NLEGMTOT asize(leg_mtot)
    const char *leg_eps[] = {
        "epsilon",
        "G\\sk",
        "g\\sk"
    };
#define NLEGEPS asize(leg_eps)
    const char *leg_aver[] = {
        "< |M|\\S2\\N >",
        "< |M| >\\S2\\N",
        "< |M|\\S2\\N > - < |M| >\\S2\\N",
        "< |M| >\\S2\\N / < |M|\\S2\\N >"
    };
#define NLEGAVER asize(leg_aver)
    const char *leg_cosaver[] = {
        "\\f{4}<|cos\\f{12}q\\f{4}\\sij\\N|>",
        "RMSD cos",
        "\\f{4}<|cos\\f{12}q\\f{4}\\siX\\N|>",
        "\\f{4}<|cos\\f{12}q\\f{4}\\siY\\N|>",
        "\\f{4}<|cos\\f{12}q\\f{4}\\siZ\\N|>"
    };
#define NLEGCOSAVER asize(leg_cosaver)
    const char *leg_adip[] = {
        "<mu>",
        "Std. Dev.",
        "Error"
    };
#define NLEGADIP asize(leg_adip)

    FILE          *outdd, *outmtot, *outaver, *outeps, *caver = nullptr;
    FILE          *dip3d = nullptr, *adip = nullptr;
    rvec          *x, *dipole = nullptr, mu_t, quad, *dipsp = nullptr;
    t_gkrbin      *gkrbin = nullptr;
    gmx_enxnm_t   *enm    = nullptr;
    t_enxframe    *fr;
    int            nframes = 1000, nre, timecheck = 0, ncolour = 0;
    ener_file_t    fmu     = nullptr;
    int            i, n, m, natom = 0, gnx_tot, teller, tel3;
    t_trxstatus   *status;
    int           *dipole_bin, ndipbin, ibin, iVol, idim = -1;
    unsigned long  mode;
    real           rcut = 0, t, t0, t1, dt, dd, rms_cos;
    rvec           dipaxis;
    matrix         box;
    gmx_bool       bCorr, bTotal, bCont;
    double         M_diff = 0, epsilon, invtel, vol_aver;
    double         mu_ave, mu_mol, M2_ave = 0, M_ave2 = 0, M_av[DIM], M_av2[DIM];
    double         M[3], M2[3], M4[3], Gk = 0, g_k = 0;
    gmx_stats_t   *Qlsq, mulsq, muframelsq = nullptr;
    ivec           iMu;
    real         **muall        = nullptr;
    rvec          *slab_dipoles = nullptr;
    const t_atom  *atom         = nullptr;
    const t_block *mols         = nullptr;
    gmx_rmpbc_t    gpbc         = nullptr;

    gnx_tot = gnx[0];
    if (ncos > 1)
    {
        gnx_tot += gnx[1];
    }

    vol_aver = 0.0;

    iVol = -1;

    /* initialize to a negative value so we can see that it wasn't picked up */
    iMu[XX] = iMu[YY] = iMu[ZZ] = -1;
    if (bMU)
    {
        fmu = open_enx(mufn, "r");
        do_enxnms(fmu, &nre, &enm);

        /* Determine the indexes of the energy grps we need */
        for (i = 0; (i < nre); i++)
        {
            if (std::strstr(enm[i].name, "Volume"))
            {
                iVol = i;
            }
            else if (std::strstr(enm[i].name, "Mu-X"))
            {
                iMu[XX] = i;
            }
            else if (std::strstr(enm[i].name, "Mu-Y"))
            {
                iMu[YY] = i;
            }
            else if (std::strstr(enm[i].name, "Mu-Z"))
            {
                iMu[ZZ] = i;
            }
        }
        if (iMu[XX] < 0 || iMu[YY] < 0 || iMu[ZZ] < 0)
        {
            gmx_fatal(FARGS, "No index for Mu-X, Mu-Y or Mu-Z energy group.");
        }
    }
    else
    {
        atom = top->atoms.atom;
        mols = &(top->mols);
    }
    if ((iVol == -1) && bMU)
    {
        printf("Using Volume from topology: %g nm^3\n", volume);
    }

    /* Correlation stuff */
    bCorr  = (corrtype[0] != 'n');
    bTotal = (corrtype[0] == 't');
    if (bCorr)
    {
        if (bTotal)
        {
            snew(muall, 1);
            snew(muall[0], nframes*DIM);
        }
        else
        {
            snew(muall, gnx[0]);
            for (i = 0; (i < gnx[0]); i++)
            {
                snew(muall[i], nframes*DIM);
            }
        }
    }

    /* Allocate array which contains for every molecule in a frame the
     * dipole moment.
     */
    if (!bMU)
    {
        snew(dipole, gnx_tot);
    }

    /* Statistics */
    snew(Qlsq, DIM);
    for (i = 0; (i < DIM); i++)
    {
        Qlsq[i] = gmx_stats_init();
    }
    mulsq = gmx_stats_init();

    /* Open all the files */
    outmtot = xvgropen(out_mtot,
                       "Total dipole moment of the simulation box vs. time",
                       "Time (ps)", "Total Dipole Moment (Debye)", oenv);
    outeps  = xvgropen(out_eps, "Epsilon and Kirkwood factors",
                       "Time (ps)", "", oenv);
    outaver = xvgropen(out_aver, "Total dipole moment",
                       "Time (ps)", "D", oenv);
    if (bSlab)
    {
        idim = axtitle[0] - 'X';
        if ((idim < 0) || (idim >= DIM))
        {
            idim = axtitle[0] - 'x';
        }
        if ((idim < 0) || (idim >= DIM))
        {
            bSlab = FALSE;
        }
        if (nslices < 2)
        {
            bSlab = FALSE;
        }
        fprintf(stderr, "axtitle = %s, nslices = %d, idim = %d\n",
                axtitle, nslices, idim);
        if (bSlab)
        {
            snew(slab_dipoles, nslices);
            fprintf(stderr, "Doing slab analysis\n");
        }
    }

    if (fnadip)
    {
        adip = xvgropen(fnadip, "Average molecular dipole", "Dipole (D)", "", oenv);
        xvgr_legend(adip, NLEGADIP, leg_adip, oenv);

    }
    if (cosaver)
    {
        caver = xvgropen(cosaver, bPairs ? "Average pair orientation" :
                         "Average absolute dipole orientation", "Time (ps)", "", oenv);
        xvgr_legend(caver, NLEGCOSAVER, bPairs ? leg_cosaver : &(leg_cosaver[1]),
                    oenv);
    }

    if (fndip3d)
    {
        snew(dipsp, gnx_tot);

        /* we need a dummy file for gnuplot */
        dip3d = (FILE *)gmx_ffopen("dummy.dat", "w");
        fprintf(dip3d, "%f %f %f", 0.0, 0.0, 0.0);
        gmx_ffclose(dip3d);

        dip3d = (FILE *)gmx_ffopen(fndip3d, "w");
        try
        {
            gmx::BinaryInformationSettings settings;
            settings.generatedByHeader(true);
            settings.linePrefix("# ");
            gmx::printBinaryInformation(dip3d, output_env_get_program_context(oenv),
                                        settings);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    /* Write legends to all the files */
    xvgr_legend(outmtot, NLEGMTOT, leg_mtot, oenv);
    xvgr_legend(outaver, NLEGAVER, leg_aver, oenv);

    if (bMU && (mu_aver == -1))
    {
        xvgr_legend(outeps, NLEGEPS-2, leg_eps, oenv);
    }
    else
    {
        xvgr_legend(outeps, NLEGEPS, leg_eps, oenv);
    }

    snew(fr, 1);
    clear_rvec(mu_t);
    teller = 0;
    /* Read the first frame from energy or traj file */
    if (bMU)
    {
        do
        {
            bCont = read_mu_from_enx(fmu, iVol, iMu, mu_t, &volume, &t, nre, fr);
            if (bCont)
            {
                timecheck = check_times(t);
                if (timecheck < 0)
                {
                    teller++;
                }
                if ((teller % 10) == 0)
                {
                    fprintf(stderr, "\r Skipping Frame %6d, time: %8.3f", teller, t);
                    fflush(stderr);
                }
            }
            else
            {
                printf("End of %s reached\n", mufn);
                break;
            }
        }
        while (bCont && (timecheck < 0));
    }
    else
    {
        natom  = read_first_x(oenv, &status, fn, &t, &x, box);
    }

    /* Calculate spacing for dipole bin (simple histogram) */
    ndipbin = 1 + static_cast<int>(mu_max/0.01);
    snew(dipole_bin, ndipbin);
    mu_ave     = 0.0;
    for (m = 0; (m < DIM); m++)
    {
        M[m] = M2[m] = M4[m] = 0.0;
    }

    if (bGkr)
    {
        /* Use 0.7 iso 0.5 to account for pressure scaling */
        /*  rcut   = 0.7*sqrt(max_cutoff2(box)); */
        rcut   = 0.7*std::sqrt(gmx::square(box[XX][XX])+gmx::square(box[YY][YY])+gmx::square(box[ZZ][ZZ]));

        gkrbin = mk_gkrbin(rcut, rcmax, bPhi, ndegrees);
    }
    gpbc = gmx_rmpbc_init(&top->idef, ePBC, natom);

    /* Start while loop over frames */
    t0     = t;
    teller = 0;
    do
    {
        if (bCorr && (teller >= nframes))
        {
            nframes += 1000;
            if (bTotal)
            {
                srenew(muall[0], nframes*DIM);
            }
            else
            {
                for (i = 0; (i < gnx_tot); i++)
                {
                    srenew(muall[i], nframes*DIM);
                }
            }
        }
        t1 = t;

        muframelsq = gmx_stats_init();

        /* Initialise */
        for (m = 0; (m < DIM); m++)
        {
            M_av2[m] = 0;
        }

        if (bMU)
        {
            /* Copy rvec into double precision local variable */
            for (m = 0; (m < DIM); m++)
            {
                M_av[m]  = mu_t[m];
            }
        }
        else
        {
            /* Initialise */
            for (m = 0; (m < DIM); m++)
            {
                M_av[m] = 0;
            }

            gmx_rmpbc(gpbc, natom, box, x);

            /* Begin loop of all molecules in frame */
            for (n = 0; (n < ncos); n++)
            {
                for (i = 0; (i < gnx[n]); i++)
                {
                    int ind0, ind1;

                    ind0  = mols->index[molindex[n][i]];
                    ind1  = mols->index[molindex[n][i]+1];

                    mol_dip(ind0, ind1, x, atom, dipole[i]);
                    gmx_stats_add_point(mulsq, 0, norm(dipole[i]), 0, 0);
                    gmx_stats_add_point(muframelsq, 0, norm(dipole[i]), 0, 0);
                    if (bSlab)
                    {
                        update_slab_dipoles(ind0, ind1, x,
                                            dipole[i], idim, nslices, slab_dipoles, box);
                    }
                    if (bQuad)
                    {
                        mol_quad(ind0, ind1, x, atom, quad);
                        for (m = 0; (m < DIM); m++)
                        {
                            gmx_stats_add_point(Qlsq[m], 0, quad[m], 0, 0);
                        }
                    }
                    if (bCorr && !bTotal)
                    {
                        tel3              = DIM*teller;
                        muall[i][tel3+XX] = dipole[i][XX];
                        muall[i][tel3+YY] = dipole[i][YY];
                        muall[i][tel3+ZZ] = dipole[i][ZZ];
                    }
                    mu_mol = 0.0;
                    for (m = 0; (m < DIM); m++)
                    {
                        M_av[m]  += dipole[i][m];               /* M per frame */
                        mu_mol   += dipole[i][m]*dipole[i][m];  /* calc. mu for distribution */
                    }
                    mu_mol = std::sqrt(mu_mol);

                    mu_ave += mu_mol;                         /* calc. the average mu */

                    /* Update the dipole distribution */
                    ibin = static_cast<int>(ndipbin*mu_mol/mu_max + 0.5);
                    if (ibin < ndipbin)
                    {
                        dipole_bin[ibin]++;
                    }

                    if (fndip3d)
                    {
                        rvec2sprvec(dipole[i], dipsp[i]);

                        if (dipsp[i][YY] > -M_PI && dipsp[i][YY] < -0.5*M_PI)
                        {
                            if (dipsp[i][ZZ] < 0.5 * M_PI)
                            {
                                ncolour = 1;
                            }
                            else
                            {
                                ncolour = 2;
                            }
                        }
                        else if (dipsp[i][YY] > -0.5*M_PI && dipsp[i][YY] < 0.0*M_PI)
                        {
                            if (dipsp[i][ZZ] < 0.5 * M_PI)
                            {
                                ncolour = 3;
                            }
                            else
                            {
                                ncolour = 4;
                            }
                        }
                        else if (dipsp[i][YY] > 0.0 && dipsp[i][YY] < 0.5*M_PI)
                        {
                            if (dipsp[i][ZZ] < 0.5 * M_PI)
                            {
                                ncolour = 5;
                            }
                            else
                            {
                                ncolour = 6;
                            }
                        }
                        else if (dipsp[i][YY] > 0.5*M_PI && dipsp[i][YY] < M_PI)
                        {
                            if (dipsp[i][ZZ] < 0.5 * M_PI)
                            {
                                ncolour = 7;
                            }
                            else
                            {
                                ncolour = 8;
                            }
                        }
                        if (dip3d)
                        {
                            fprintf(dip3d, "set arrow %d from %f, %f, %f to %f, %f, %f lt %d  # %d %d\n",
                                    i+1,
                                    x[ind0][XX],
                                    x[ind0][YY],
                                    x[ind0][ZZ],
                                    x[ind0][XX]+dipole[i][XX]/25,
                                    x[ind0][YY]+dipole[i][YY]/25,
                                    x[ind0][ZZ]+dipole[i][ZZ]/25,
                                    ncolour, ind0, i);
                        }
                    }
                } /* End loop of all molecules in frame */

                if (dip3d)
                {
                    fprintf(dip3d, "set title \"t = %4.3f\"\n", t);
                    fprintf(dip3d, "set xrange [0.0:%4.2f]\n", box[XX][XX]);
                    fprintf(dip3d, "set yrange [0.0:%4.2f]\n", box[YY][YY]);
                    fprintf(dip3d, "set zrange [0.0:%4.2f]\n\n", box[ZZ][ZZ]);
                    fprintf(dip3d, "splot 'dummy.dat' using 1:2:3 w vec\n");
                    fprintf(dip3d, "pause -1 'Hit return to continue'\n");
                }
            }
        }
        /* Compute square of total dipole */
        for (m = 0; (m < DIM); m++)
        {
            M_av2[m] = M_av[m]*M_av[m];
        }

        if (cosaver)
        {
            compute_avercos(gnx_tot, dipole, &dd, dipaxis, bPairs);
            rms_cos = std::sqrt(gmx::square(dipaxis[XX]-0.5)+
                                gmx::square(dipaxis[YY]-0.5)+
                                gmx::square(dipaxis[ZZ]-0.5));
            if (bPairs)
            {
                fprintf(caver, "%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n",
                        t, dd, rms_cos, dipaxis[XX], dipaxis[YY], dipaxis[ZZ]);
            }
            else
            {
                fprintf(caver, "%10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n",
                        t, rms_cos, dipaxis[XX], dipaxis[YY], dipaxis[ZZ]);
            }
        }

        if (bGkr)
        {
            do_gkr(gkrbin, ncos, gnx, molindex, mols->index, x, dipole, ePBC, box,
                   atom, gkatom);
        }

        if (bTotal)
        {
            tel3              = DIM*teller;
            muall[0][tel3+XX] = M_av[XX];
            muall[0][tel3+YY] = M_av[YY];
            muall[0][tel3+ZZ] = M_av[ZZ];
        }

        /* Write to file the total dipole moment of the box, and its components
         * for this frame.
         */
        if ((skip == 0) || ((teller % skip) == 0))
        {
            fprintf(outmtot, "%10g  %12.8e %12.8e %12.8e %12.8e\n",
                    t, M_av[XX], M_av[YY], M_av[ZZ],
                    std::sqrt(M_av2[XX]+M_av2[YY]+M_av2[ZZ]));
        }

        for (m = 0; (m < DIM); m++)
        {
            M[m]  += M_av[m];
            M2[m] += M_av2[m];
            M4[m] += gmx::square(M_av2[m]);
        }
        /* Increment loop counter */
        teller++;

        /* Calculate for output the running averages */
        invtel  = 1.0/teller;
        M2_ave  = (M2[XX]+M2[YY]+M2[ZZ])*invtel;
        M_ave2  = invtel*(invtel*(M[XX]*M[XX] + M[YY]*M[YY] + M[ZZ]*M[ZZ]));
        M_diff  = M2_ave - M_ave2;

        /* Compute volume from box in traj, else we use the one from above */
        if (!bMU)
        {
            volume  = det(box);
        }
        vol_aver += volume;

        epsilon = calc_eps(M_diff, (vol_aver/teller), epsilonRF, temp);

        /* Calculate running average for dipole */
        if (mu_ave != 0)
        {
            mu_aver = (mu_ave/gnx_tot)*invtel;
        }

        if ((skip == 0) || ((teller % skip) == 0))
        {
            /* Write to file < |M|^2 >, |< M >|^2. And the difference between
             * the two. Here M is sum mu_i. Further write the finite system
             * Kirkwood G factor and epsilon.
             */
            fprintf(outaver, "%10g  %10.3e %10.3e %10.3e %10.3e\n",
                    t, M2_ave, M_ave2, M_diff, M_ave2/M2_ave);

            if (fnadip)
            {
                real aver;
                gmx_stats_get_average(muframelsq, &aver);
                fprintf(adip, "%10g %f \n", t, aver);
            }
            /*if (dipole)
               printf("%f %f\n", norm(dipole[0]), norm(dipole[1]));
             */
            if (!bMU || (mu_aver != -1))
            {
                /* Finite system Kirkwood G-factor */
                Gk = M_diff/(gnx_tot*mu_aver*mu_aver);
                /* Infinite system Kirkwood G-factor */
                if (epsilonRF == 0.0)
                {
                    g_k = ((2*epsilon+1)*Gk/(3*epsilon));
                }
                else
                {
                    g_k = ((2*epsilonRF+epsilon)*(2*epsilon+1)*
                           Gk/(3*epsilon*(2*epsilonRF+1)));
                }

                fprintf(outeps, "%10g  %10.3e %10.3e %10.3e\n", t, epsilon, Gk, g_k);

            }
            else
            {
                fprintf(outeps, "%10g  %12.8e\n", t, epsilon);
            }
        }
        gmx_stats_free(muframelsq);

        if (bMU)
        {
            bCont = read_mu_from_enx(fmu, iVol, iMu, mu_t, &volume, &t, nre, fr);
        }
        else
        {
            bCont = read_next_x(oenv, status, &t, x, box);
        }
        timecheck = check_times(t);
    }
    while (bCont && (timecheck == 0) );

    gmx_rmpbc_done(gpbc);

    if (!bMU)
    {
        close_trx(status);
    }

    xvgrclose(outmtot);
    xvgrclose(outaver);
    xvgrclose(outeps);

    if (fnadip)
    {
        xvgrclose(adip);
    }

    if (cosaver)
    {
        xvgrclose(caver);
    }

    if (dip3d)
    {
        fprintf(dip3d, "set xrange [0.0:%4.2f]\n", box[XX][XX]);
        fprintf(dip3d, "set yrange [0.0:%4.2f]\n", box[YY][YY]);
        fprintf(dip3d, "set zrange [0.0:%4.2f]\n\n", box[ZZ][ZZ]);
        fprintf(dip3d, "splot 'dummy.dat' using 1:2:3 w vec\n");
        fprintf(dip3d, "pause -1 'Hit return to continue'\n");
        gmx_ffclose(dip3d);
    }

    if (bSlab)
    {
        dump_slab_dipoles(slabfn, idim, nslices, slab_dipoles, box, teller, oenv);
        sfree(slab_dipoles);
    }

    vol_aver /= teller;
    printf("Average volume over run is %g\n", vol_aver);
    if (bGkr)
    {
        print_gkrbin(gkrfn, gkrbin, gnx[0], teller, vol_aver, oenv);
        print_cmap(cmap, gkrbin, nlevels);
    }
    /* Autocorrelation function */
    if (bCorr)
    {
        if (teller < 2)
        {
            printf("Not enough frames for autocorrelation\n");
        }
        else
        {
            dt = (t1 - t0)/(teller-1);
            printf("t0 %g, t %g, teller %d\n", t0, t, teller);

            mode = eacVector;

            if (bTotal)
            {
                do_autocorr(corf, oenv, "Autocorrelation Function of Total Dipole",
                            teller, 1, muall, dt, mode, TRUE);
            }
            else
            {
                do_autocorr(corf, oenv, "Dipole Autocorrelation Function",
                            teller, gnx_tot, muall, dt,
                            mode, std::strcmp(corrtype, "molsep"));
            }
        }
    }
    if (!bMU)
    {
        real aver, sigma, error;

        gmx_stats_get_ase(mulsq, &aver, &sigma, &error);
        printf("\nDipole moment (Debye)\n");
        printf("---------------------\n");
        printf("Average  = %8.4f  Std. Dev. = %8.4f  Error = %8.4f\n",
               aver, sigma, error);
        if (bQuad)
        {
            rvec a, s, e;
            for (m = 0; (m < DIM); m++)
            {
                gmx_stats_get_ase(mulsq, &(a[m]), &(s[m]), &(e[m]));
            }

            printf("\nQuadrupole moment (Debye-Ang)\n");
            printf("-----------------------------\n");
            printf("Averages  = %8.4f  %8.4f  %8.4f\n", a[XX], a[YY], a[ZZ]);
            printf("Std. Dev. = %8.4f  %8.4f  %8.4f\n", s[XX], s[YY], s[ZZ]);
            printf("Error     = %8.4f  %8.4f  %8.4f\n", e[XX], e[YY], e[ZZ]);
        }
        printf("\n");
    }
    printf("The following averages for the complete trajectory have been calculated:\n\n");
    printf(" Total < M_x > = %g Debye\n", M[XX]/teller);
    printf(" Total < M_y > = %g Debye\n", M[YY]/teller);
    printf(" Total < M_z > = %g Debye\n\n", M[ZZ]/teller);

    printf(" Total < M_x^2 > = %g Debye^2\n", M2[XX]/teller);
    printf(" Total < M_y^2 > = %g Debye^2\n", M2[YY]/teller);
    printf(" Total < M_z^2 > = %g Debye^2\n\n", M2[ZZ]/teller);

    printf(" Total < |M|^2 > = %g Debye^2\n", M2_ave);
    printf(" Total |< M >|^2 = %g Debye^2\n\n", M_ave2);

    printf(" < |M|^2 > - |< M >|^2 = %g Debye^2\n\n", M_diff);

    if (!bMU || (mu_aver != -1))
    {
        printf("Finite system Kirkwood g factor G_k = %g\n", Gk);
        printf("Infinite system Kirkwood g factor g_k = %g\n\n", g_k);
    }
    printf("Epsilon = %g\n", epsilon);

    if (!bMU)
    {
        /* Write to file the dipole moment distibution during the simulation.
         */
        outdd = xvgropen(dipdist, "Dipole Moment Distribution", "mu (Debye)", "", oenv);
        for (i = 0; (i < ndipbin); i++)
        {
            fprintf(outdd, "%10g  %10f\n",
                    (i*mu_max)/ndipbin, static_cast<real>(dipole_bin[i])/teller);
        }
        xvgrclose(outdd);
        sfree(dipole_bin);
    }
    if (bGkr)
    {
        done_gkrbin(&gkrbin);
    }
}

static void dipole_atom2molindex(int *n, int *index, const t_block *mols)
{
    int nmol, i, j, m;

    nmol = 0;
    i    = 0;
    while (i < *n)
    {
        m = 0;
        while (m < mols->nr && index[i] != mols->index[m])
        {
            m++;
        }
        if (m == mols->nr)
        {
            gmx_fatal(FARGS, "index[%d]=%d does not correspond to the first atom of a molecule", i+1, index[i]+1);
        }
        for (j = mols->index[m]; j < mols->index[m+1]; j++)
        {
            if (i >= *n || index[i] != j)
            {
                gmx_fatal(FARGS, "The index group is not a set of whole molecules");
            }
            i++;
        }
        /* Modify the index in place */
        index[nmol++] = m;
    }
    printf("There are %d molecules in the selection\n", nmol);

    *n = nmol;
}
int gmx_dipoles(int argc, char *argv[])
{
    const char       *desc[] = {
        "[THISMODULE] computes the total dipole plus fluctuations of a simulation",
        "system. From this you can compute e.g. the dielectric constant for",
        "low-dielectric media.",
        "For molecules with a net charge, the net charge is subtracted at",
        "center of mass of the molecule.[PAR]",
        "The file [TT]Mtot.xvg[tt] contains the total dipole moment of a frame, the",
        "components as well as the norm of the vector.",
        "The file [TT]aver.xvg[tt] contains [CHEVRON][MAG][GRK]mu[grk][mag]^2[chevron] and [MAG][CHEVRON][GRK]mu[grk][chevron][mag]^2 during the",
        "simulation.",
        "The file [TT]dipdist.xvg[tt] contains the distribution of dipole moments during",
        "the simulation",
        "The value of [TT]-mumax[tt] is used as the highest value in the distribution graph.[PAR]",
        "Furthermore, the dipole autocorrelation function will be computed when",
        "option [TT]-corr[tt] is used. The output file name is given with the [TT]-c[tt]",
        "option.",
        "The correlation functions can be averaged over all molecules",
        "([TT]mol[tt]), plotted per molecule separately ([TT]molsep[tt])",
        "or it can be computed over the total dipole moment of the simulation box",
        "([TT]total[tt]).[PAR]",
        "Option [TT]-g[tt] produces a plot of the distance dependent Kirkwood",
        "G-factor, as well as the average cosine of the angle between the dipoles",
        "as a function of the distance. The plot also includes gOO and hOO",
        "according to Nymand & Linse, J. Chem. Phys. 112 (2000) pp 6386-6395. In the same plot, ",
        "we also include the energy per scale computed by taking the inner product of",
        "the dipoles divided by the distance to the third power.[PAR]",
        "[PAR]",
        "EXAMPLES[PAR]",
        "[TT]gmx dipoles -corr mol -P 1 -o dip_sqr -mu 2.273 -mumax 5.0[tt][PAR]",
        "This will calculate the autocorrelation function of the molecular",
        "dipoles using a first order Legendre polynomial of the angle of the",
        "dipole vector and itself a time t later. For this calculation 1001",
        "frames will be used. Further, the dielectric constant will be calculated",
        "using an [TT]-epsilonRF[tt] of infinity (default), temperature of 300 K (default) and",
        "an average dipole moment of the molecule of 2.273 (SPC). For the",
        "distribution function a maximum of 5.0 will be used."
    };
    real              mu_max     = 5, mu_aver = -1, rcmax = 0;
    real              epsilonRF  = 0.0, temp = 300;
    gmx_bool          bPairs     = TRUE, bPhi = FALSE, bQuad = FALSE;
    const char       *corrtype[] = {nullptr, "none", "mol", "molsep", "total", nullptr};
    const char       *axtitle    = "Z";
    int               nslices    = 10; /* nr of slices defined       */
    int               skip       = 0, nFA = 0, nFB = 0, ncos = 1;
    int               nlevels    = 20, ndegrees = 90;
    gmx_output_env_t *oenv;
    t_pargs           pa[] = {
        { "-mu",       FALSE, etREAL, {&mu_aver},
          "dipole of a single molecule (in Debye)" },
        { "-mumax",    FALSE, etREAL, {&mu_max},
          "max dipole in Debye (for histogram)" },
        { "-epsilonRF", FALSE, etREAL, {&epsilonRF},
          "[GRK]epsilon[grk] of the reaction field used during the simulation, needed for dielectric constant calculation. WARNING: 0.0 means infinity (default)" },
        { "-skip",     FALSE, etINT, {&skip},
          "Skip steps in the output (but not in the computations)" },
        { "-temp",     FALSE, etREAL, {&temp},
          "Average temperature of the simulation (needed for dielectric constant calculation)" },
        { "-corr",     FALSE, etENUM, {corrtype},
          "Correlation function to calculate" },
        { "-pairs",    FALSE, etBOOL, {&bPairs},
          "Calculate [MAG][COS][GRK]theta[grk][cos][mag] between all pairs of molecules. May be slow" },
        { "-quad",     FALSE, etBOOL, {&bQuad},
          "Take quadrupole into account"},
        { "-ncos",     FALSE, etINT, {&ncos},
          "Must be 1 or 2. Determines whether the [CHEVRON][COS][GRK]theta[grk][cos][chevron] is computed between all molecules in one group, or between molecules in two different groups. This turns on the [TT]-g[tt] flag." },
        { "-axis",     FALSE, etSTR, {&axtitle},
          "Take the normal on the computational box in direction X, Y or Z." },
        { "-sl",       FALSE, etINT, {&nslices},
          "Divide the box into this number of slices." },
        { "-gkratom",  FALSE, etINT, {&nFA},
          "Use the n-th atom of a molecule (starting from 1) to calculate the distance between molecules rather than the center of charge (when 0) in the calculation of distance dependent Kirkwood factors" },
        { "-gkratom2", FALSE, etINT, {&nFB},
          "Same as previous option in case ncos = 2, i.e. dipole interaction between two groups of molecules" },
        { "-rcmax",    FALSE, etREAL, {&rcmax},
          "Maximum distance to use in the dipole orientation distribution (with ncos == 2). If zero, a criterion based on the box length will be used." },
        { "-phi",      FALSE, etBOOL, {&bPhi},
          "Plot the 'torsion angle' defined as the rotation of the two dipole vectors around the distance vector between the two molecules in the [REF].xpm[ref] file from the [TT]-cmap[tt] option. By default the cosine of the angle between the dipoles is plotted." },
        { "-nlevels",  FALSE, etINT, {&nlevels},
          "Number of colors in the cmap output" },
        { "-ndegrees", FALSE, etINT, {&ndegrees},
          "Number of divisions on the [IT]y[it]-axis in the cmap output (for 180 degrees)" }
    };
    int              *gnx;
    int               nFF[2];
    int             **grpindex;
    char            **grpname = nullptr;
    gmx_bool          bGkr, bMU, bSlab;
    t_filenm          fnm[] = {
        { efEDR, "-en", nullptr,         ffOPTRD },
        { efTRX, "-f", nullptr,           ffREAD },
        { efTPR, nullptr, nullptr,           ffREAD },
        { efNDX, nullptr, nullptr,           ffOPTRD },
        { efXVG, "-o",   "Mtot",       ffWRITE },
        { efXVG, "-eps", "epsilon",    ffWRITE },
        { efXVG, "-a",   "aver",       ffWRITE },
        { efXVG, "-d",   "dipdist",    ffWRITE },
        { efXVG, "-c",   "dipcorr",    ffOPTWR },
        { efXVG, "-g",   "gkr",        ffOPTWR },
        { efXVG, "-adip", "adip",       ffOPTWR },
        { efXVG, "-dip3d", "dip3d",    ffOPTWR },
        { efXVG, "-cos", "cosaver",    ffOPTWR },
        { efXPM, "-cmap", "cmap",       ffOPTWR },
        { efXVG, "-slab", "slab",       ffOPTWR }
    };
#define NFILE asize(fnm)
    int               npargs;
    t_pargs          *ppa;
    t_topology       *top;
    int               ePBC;
    int               k, natoms;
    matrix            box;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW,
                           NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    printf("Using %g as mu_max and %g as the dipole moment.\n",
           mu_max, mu_aver);
    if (epsilonRF == 0.0)
    {
        printf("WARNING: EpsilonRF = 0.0, this really means EpsilonRF = infinity\n");
    }

    bMU   = opt2bSet("-en", NFILE, fnm);
    if (bMU)
    {
        gmx_fatal(FARGS, "Due to new ways of treating molecules in GROMACS the total dipole in the energy file may be incorrect, because molecules can be split over periodic boundary conditions before computing the dipole. Please use your trajectory file.");
    }
    bGkr  = opt2bSet("-g", NFILE, fnm);
    if (opt2parg_bSet("-ncos", asize(pa), pa))
    {
        if ((ncos != 1) && (ncos != 2))
        {
            gmx_fatal(FARGS, "ncos has to be either 1 or 2");
        }
        bGkr = TRUE;
    }
    bSlab = (opt2bSet("-slab", NFILE, fnm) || opt2parg_bSet("-sl", asize(pa), pa) ||
             opt2parg_bSet("-axis", asize(pa), pa));
    if (bMU)
    {
        if (bQuad)
        {
            printf("WARNING: Can not determine quadrupoles from energy file\n");
            bQuad = FALSE;
        }
        if (bGkr)
        {
            printf("WARNING: Can not determine Gk(r) from energy file\n");
            bGkr  = FALSE;
            ncos  = 1;
        }
        if (mu_aver == -1)
        {
            printf("WARNING: Can not calculate Gk and gk, since you did\n"
                   "         not enter a valid dipole for the molecules\n");
        }
    }

    snew(top, 1);
    ePBC = read_tpx_top(ftp2fn(efTPR, NFILE, fnm), nullptr, box,
                        &natoms, nullptr, nullptr, top);

    snew(gnx, ncos);
    snew(grpname, ncos);
    snew(grpindex, ncos);
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm),
              ncos, gnx, grpindex, grpname);
    for (k = 0; (k < ncos); k++)
    {
        dipole_atom2molindex(&gnx[k], grpindex[k], &(top->mols));
        neutralize_mols(gnx[k], grpindex[k], &(top->mols), top->atoms.atom);
    }
    nFF[0] = nFA;
    nFF[1] = nFB;
    do_dip(top, ePBC, det(box), ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn("-eps", NFILE, fnm),
           opt2fn("-a", NFILE, fnm), opt2fn("-d", NFILE, fnm),
           opt2fn_null("-cos", NFILE, fnm),
           opt2fn_null("-dip3d", NFILE, fnm), opt2fn_null("-adip", NFILE, fnm),
           bPairs, corrtype[0],
           opt2fn("-c", NFILE, fnm),
           bGkr,    opt2fn("-g", NFILE, fnm),
           bPhi,    &nlevels,  ndegrees,
           ncos,
           opt2fn("-cmap", NFILE, fnm), rcmax,
           bQuad, bMU,     opt2fn("-en", NFILE, fnm),
           gnx, grpindex, mu_max, mu_aver, epsilonRF, temp, nFF, skip,
           bSlab, nslices, axtitle, opt2fn("-slab", NFILE, fnm), oenv);

    do_view(oenv, opt2fn("-o", NFILE, fnm), "-autoscale xy -nxy");
    do_view(oenv, opt2fn("-eps", NFILE, fnm), "-autoscale xy -nxy");
    do_view(oenv, opt2fn("-a", NFILE, fnm), "-autoscale xy -nxy");
    do_view(oenv, opt2fn("-d", NFILE, fnm), "-autoscale xy");
    do_view(oenv, opt2fn("-c", NFILE, fnm), "-autoscale xy");

    return 0;
}
