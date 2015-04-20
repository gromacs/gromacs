/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/tables.h"

#include <math.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/* All the possible (implemented) table functions */
enum {
    etabLJ6,
    etabLJ12,
    etabLJ6Shift,
    etabLJ12Shift,
    etabShift,
    etabRF,
    etabRF_ZERO,
    etabCOUL,
    etabEwald,
    etabEwaldSwitch,
    etabEwaldUser,
    etabEwaldUserSwitch,
    etabLJ6Ewald,
    etabLJ6Switch,
    etabLJ12Switch,
    etabCOULSwitch,
    etabLJ6Encad,
    etabLJ12Encad,
    etabCOULEncad,
    etabEXPMIN,
    etabUSER,
    etabNR
};

/** Evaluates to true if the table type contains user data. */
#define ETAB_USER(e)  ((e) == etabUSER || \
                       (e) == etabEwaldUser || (e) == etabEwaldUserSwitch)

typedef struct {
    const char *name;
    gmx_bool    bCoulomb;
} t_tab_props;

/* This structure holds name and a flag that tells whether
   this is a Coulomb type funtion */
static const t_tab_props tprops[etabNR] = {
    { "LJ6",  FALSE },
    { "LJ12", FALSE },
    { "LJ6Shift", FALSE },
    { "LJ12Shift", FALSE },
    { "Shift", TRUE },
    { "RF", TRUE },
    { "RF-zero", TRUE },
    { "COUL", TRUE },
    { "Ewald", TRUE },
    { "Ewald-Switch", TRUE },
    { "Ewald-User", TRUE },
    { "Ewald-User-Switch", TRUE },
    { "LJ6Ewald", FALSE },
    { "LJ6Switch", FALSE },
    { "LJ12Switch", FALSE },
    { "COULSwitch", TRUE },
    { "LJ6-Encad shift", FALSE },
    { "LJ12-Encad shift", FALSE },
    { "COUL-Encad shift",  TRUE },
    { "EXPMIN", FALSE },
    { "USER", FALSE },
};

/* Index in the table that says which function to use */
enum {
    etiCOUL, etiLJ6, etiLJ12, etiNR
};

typedef struct {
    int     nx, nx0;
    double  tabscale;
    double *x, *v, *f;
} t_tabledata;

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))

double v_q_ewald_lr(double beta, double r)
{
    if (r == 0)
    {
        return beta*2/sqrt(M_PI);
    }
    else
    {
        return gmx_erfd(beta*r)/r;
    }
}

double v_lj_ewald_lr(double beta, double r)
{
    double br, br2, br4, r6, factor;
    if (r == 0)
    {
        return pow(beta, 6)/6;
    }
    else
    {
        br     = beta*r;
        br2    = br*br;
        br4    = br2*br2;
        r6     = pow(r, 6.0);
        factor = (1.0 - exp(-br2)*(1 + br2 + 0.5*br4))/r6;
        return factor;
    }
}

void table_spline3_fill_ewald_lr(real                                 *table_f,
                                 real                                 *table_v,
                                 real                                 *table_fdv0,
                                 int                                   ntab,
                                 double                                dx,
                                 real                                  beta,
                                 real_space_grid_contribution_computer v_lr)
{
    real     tab_max;
    int      i, i_inrange;
    double   dc, dc_new;
    gmx_bool bOutOfRange;
    double   v_r0, v_r1, v_inrange, vi, a0, a1, a2dx;
    double   x_r0;

    /* This function is called using either v_ewald_lr or v_lj_ewald_lr as a function argument
     * depending on wether we should create electrostatic or Lennard-Jones Ewald tables.
     */

    if (ntab < 2)
    {
        gmx_fatal(FARGS, "Can not make a spline table with less than 2 points");
    }

    /* We need some margin to be able to divide table values by r
     * in the kernel and also to do the integration arithmetics
     * without going out of range. Furthemore, we divide by dx below.
     */
    tab_max = GMX_REAL_MAX*0.0001;

    /* This function produces a table with:
     * maximum energy error: V'''/(6*12*sqrt(3))*dx^3
     * maximum force error:  V'''/(6*4)*dx^2
     * The rms force error is the max error times 1/sqrt(5)=0.45.
     */

    bOutOfRange = FALSE;
    i_inrange   = ntab;
    v_inrange   = 0;
    dc          = 0;
    for (i = ntab-1; i >= 0; i--)
    {
        x_r0 = i*dx;

        v_r0 = (*v_lr)(beta, x_r0);

        if (!bOutOfRange)
        {
            i_inrange = i;
            v_inrange = v_r0;

            vi = v_r0;
        }
        else
        {
            /* Linear continuation for the last point in range */
            vi = v_inrange - dc*(i - i_inrange)*dx;
        }

        if (table_v != NULL)
        {
            table_v[i] = vi;
        }

        if (i == 0)
        {
            continue;
        }

        /* Get the potential at table point i-1 */
        v_r1 = (*v_lr)(beta, (i-1)*dx);

        if (v_r1 != v_r1 || v_r1 < -tab_max || v_r1 > tab_max)
        {
            bOutOfRange = TRUE;
        }

        if (!bOutOfRange)
        {
            /* Calculate the average second derivative times dx over interval i-1 to i.
             * Using the function values at the end points and in the middle.
             */
            a2dx = (v_r0 + v_r1 - 2*(*v_lr)(beta, x_r0-0.5*dx))/(0.25*dx);
            /* Set the derivative of the spline to match the difference in potential
             * over the interval plus the average effect of the quadratic term.
             * This is the essential step for minimizing the error in the force.
             */
            dc = (v_r0 - v_r1)/dx + 0.5*a2dx;
        }

        if (i == ntab - 1)
        {
            /* Fill the table with the force, minus the derivative of the spline */
            table_f[i] = -dc;
        }
        else
        {
            /* tab[i] will contain the average of the splines over the two intervals */
            table_f[i] += -0.5*dc;
        }

        if (!bOutOfRange)
        {
            /* Make spline s(x) = a0 + a1*(x - xr) + 0.5*a2*(x - xr)^2
             * matching the potential at the two end points
             * and the derivative dc at the end point xr.
             */
            a0   = v_r0;
            a1   = dc;
            a2dx = (a1*dx + v_r1 - a0)*2/dx;

            /* Set dc to the derivative at the next point */
            dc_new = a1 - a2dx;

            if (dc_new != dc_new || dc_new < -tab_max || dc_new > tab_max)
            {
                bOutOfRange = TRUE;
            }
            else
            {
                dc = dc_new;
            }
        }

        table_f[(i-1)] = -0.5*dc;
    }
    /* Currently the last value only contains half the force: double it */
    table_f[0] *= 2;

    if (table_v != NULL && table_fdv0 != NULL)
    {
        /* Copy to FDV0 table too. Allocation occurs in forcerec.c,
         * init_ewald_f_table().
         */
        for (i = 0; i < ntab-1; i++)
        {
            table_fdv0[4*i]     = table_f[i];
            table_fdv0[4*i+1]   = table_f[i+1]-table_f[i];
            table_fdv0[4*i+2]   = table_v[i];
            table_fdv0[4*i+3]   = 0.0;
        }
        table_fdv0[4*(ntab-1)]    = table_f[(ntab-1)];
        table_fdv0[4*(ntab-1)+1]  = -table_f[(ntab-1)];
        table_fdv0[4*(ntab-1)+2]  = table_v[(ntab-1)];
        table_fdv0[4*(ntab-1)+3]  = 0.0;
    }
}

/* Returns the spacing for a function using the maximum of
 * the third derivative, x_scale (unit 1/length)
 * and function tolerance.
 */
static double spline3_table_scale(double third_deriv_max,
                                  double x_scale,
                                  double func_tol)
{
    double deriv_tol;
    double sc_deriv, sc_func;

    /* Force tolerance: single precision accuracy */
    deriv_tol = GMX_FLOAT_EPS;
    sc_deriv  = sqrt(third_deriv_max/(6*4*deriv_tol*x_scale))*x_scale;

    /* Don't try to be more accurate on energy than the precision */
    func_tol  = max(func_tol, GMX_REAL_EPS);
    sc_func   = pow(third_deriv_max/(6*12*sqrt(3)*func_tol), 1.0/3.0)*x_scale;

    return max(sc_deriv, sc_func);
}

/* The scale (1/spacing) for third order spline interpolation
 * of the Ewald mesh contribution which needs to be subtracted
 * from the non-bonded interactions.
 * Since there is currently only one spacing for Coulomb and LJ,
 * the finest spacing is used if both Ewald types are present.
 *
 * Note that we could also implement completely separate tables
 * for Coulomb and LJ Ewald, each with their own spacing.
 * The current setup with the same spacing can provide slightly
 * faster kernels with both Coulomb and LJ Ewald, especially
 * when interleaving both tables (currently not implemented).
 */
real ewald_spline3_table_scale(const interaction_const_t *ic)
{
    real sc;

    sc = 0;

    if (ic->eeltype == eelEWALD || EEL_PME(ic->eeltype))
    {
        double erf_x_d3 = 1.0522; /* max of (erf(x)/x)''' */
        double etol;
        real   sc_q;

        /* Energy tolerance: 0.1 times the cut-off jump */
        etol  = 0.1*gmx_erfc(ic->ewaldcoeff_q*ic->rcoulomb);

        sc_q  = spline3_table_scale(erf_x_d3, ic->ewaldcoeff_q, etol);

        if (debug)
        {
            fprintf(debug, "Ewald Coulomb quadratic spline table spacing: %f 1/nm\n", 1/sc_q);
        }

        sc    = max(sc, sc_q);
    }

    if (EVDW_PME(ic->vdwtype))
    {
        double func_d3 = 0.42888; /* max of (x^-6 (1 - exp(-x^2)(1+x^2+x^4/2)))''' */
        double xrc2, etol;
        real   sc_lj;

        /* Energy tolerance: 0.1 times the cut-off jump */
        xrc2  = sqr(ic->ewaldcoeff_lj*ic->rvdw);
        etol  = 0.1*exp(-xrc2)*(1 + xrc2 + xrc2*xrc2/2.0);

        sc_lj = spline3_table_scale(func_d3, ic->ewaldcoeff_lj, etol);

        if (debug)
        {
            fprintf(debug, "Ewald LJ quadratic spline table spacing: %f 1/nm\n", 1/sc_lj);
        }

        sc = max(sc, sc_lj);
    }

    return sc;
}

/* Calculate the potential and force for an r value
 * in exactly the same way it is done in the inner loop.
 * VFtab is a pointer to the table data, offset is
 * the point where we should begin and stride is
 * 4 if we have a buckingham table, 3 otherwise.
 * If you want to evaluate table no N, set offset to 4*N.
 *
 * We use normal precision here, since that is what we
 * will use in the inner loops.
 */
static void evaluate_table(real VFtab[], int offset, int stride,
                           real tabscale, real r, real *y, real *yp)
{
    int  n;
    real rt, eps, eps2;
    real Y, F, Geps, Heps2, Fp;

    rt       =  r*tabscale;
    n        =  (int)rt;
    eps      =  rt - n;
    eps2     =  eps*eps;
    n        =  offset+stride*n;
    Y        =  VFtab[n];
    F        =  VFtab[n+1];
    Geps     =  eps*VFtab[n+2];
    Heps2    =  eps2*VFtab[n+3];
    Fp       =  F+Geps+Heps2;
    *y       =  Y+eps*Fp;
    *yp      =  (Fp+Geps+2.0*Heps2)*tabscale;
}

static void copy2table(int n, int offset, int stride,
                       double x[], double Vtab[], double Ftab[], real scalefactor,
                       real dest[])
{
/* Use double prec. for the intermediary variables
 * and temporary x/vtab/vtab2 data to avoid unnecessary
 * loss of precision.
 */
    int    i, nn0;
    double F, G, H, h;

    h = 0;
    for (i = 0; (i < n); i++)
    {
        if (i < n-1)
        {
            h   = x[i+1] - x[i];
            F   = -Ftab[i]*h;
            G   =  3*(Vtab[i+1] - Vtab[i]) + (Ftab[i+1] + 2*Ftab[i])*h;
            H   = -2*(Vtab[i+1] - Vtab[i]) - (Ftab[i+1] +   Ftab[i])*h;
        }
        else
        {
            /* Fill the last entry with a linear potential,
             * this is mainly for rounding issues with angle and dihedral potentials.
             */
            F   = -Ftab[i]*h;
            G   = 0;
            H   = 0;
        }
        nn0         = offset + i*stride;
        dest[nn0]   = scalefactor*Vtab[i];
        dest[nn0+1] = scalefactor*F;
        dest[nn0+2] = scalefactor*G;
        dest[nn0+3] = scalefactor*H;
    }
}

static void init_table(int n, int nx0,
                       double tabscale, t_tabledata *td, gmx_bool bAlloc)
{
    int i;

    td->nx       = n;
    td->nx0      = nx0;
    td->tabscale = tabscale;
    if (bAlloc)
    {
        snew(td->x, td->nx);
        snew(td->v, td->nx);
        snew(td->f, td->nx);
    }
    for (i = 0; (i < td->nx); i++)
    {
        td->x[i] = i/tabscale;
    }
}

static void spline_forces(int nx, double h, double v[], gmx_bool bS3, gmx_bool bE3,
                          double f[])
{
    int    start, end, i;
    double v3, b_s, b_e, b;
    double beta, *gamma;

    /* Formulas can be found in:
     * H.J.C. Berendsen, Simulating the Physical World, Cambridge 2007
     */

    if (nx < 4 && (bS3 || bE3))
    {
        gmx_fatal(FARGS, "Can not generate splines with third derivative boundary conditions with less than 4 (%d) points", nx);
    }

    /* To make life easy we initially set the spacing to 1
     * and correct for this at the end.
     */
    beta = 2;
    if (bS3)
    {
        /* Fit V''' at the start */
        v3  = v[3] - 3*v[2] + 3*v[1] - v[0];
        if (debug)
        {
            fprintf(debug, "The left third derivative is %g\n", v3/(h*h*h));
        }
        b_s   = 2*(v[1] - v[0]) + v3/6;
        start = 0;

        if (FALSE)
        {
            /* Fit V'' at the start */
            real v2;

            v2  = -v[3] + 4*v[2] - 5*v[1] + 2*v[0];
            /* v2  = v[2] - 2*v[1] + v[0]; */
            if (debug)
            {
                fprintf(debug, "The left second derivative is %g\n", v2/(h*h));
            }
            b_s   = 3*(v[1] - v[0]) - v2/2;
            start = 0;
        }
    }
    else
    {
        b_s   = 3*(v[2] - v[0]) + f[0]*h;
        start = 1;
    }
    if (bE3)
    {
        /* Fit V''' at the end */
        v3  = v[nx-1] - 3*v[nx-2] + 3*v[nx-3] - v[nx-4];
        if (debug)
        {
            fprintf(debug, "The right third derivative is %g\n", v3/(h*h*h));
        }
        b_e = 2*(v[nx-1] - v[nx-2]) + v3/6;
        end = nx;
    }
    else
    {
        /* V'=0 at the end */
        b_e = 3*(v[nx-1] - v[nx-3]) + f[nx-1]*h;
        end = nx - 1;
    }

    snew(gamma, nx);
    beta = (bS3 ? 1 : 4);

    /* For V'' fitting */
    /* beta = (bS3 ? 2 : 4); */

    f[start] = b_s/beta;
    for (i = start+1; i < end; i++)
    {
        gamma[i] = 1/beta;
        beta     = 4 - gamma[i];
        b        =  3*(v[i+1] - v[i-1]);
        f[i]     = (b - f[i-1])/beta;
    }
    gamma[end-1] = 1/beta;
    beta         = (bE3 ? 1 : 4) - gamma[end-1];
    f[end-1]     = (b_e - f[end-2])/beta;

    for (i = end-2; i >= start; i--)
    {
        f[i] -= gamma[i+1]*f[i+1];
    }
    sfree(gamma);

    /* Correct for the minus sign and the spacing */
    for (i = start; i < end; i++)
    {
        f[i] = -f[i]/h;
    }
}

static void set_forces(FILE *fp, int angle,
                       int nx, double h, double v[], double f[],
                       int table)
{
    int start, end;

    if (angle == 2)
    {
        gmx_fatal(FARGS,
                  "Force generation for dihedral tables is not (yet) implemented");
    }

    start = 0;
    while (v[start] == 0)
    {
        start++;
    }

    end = nx;
    while (v[end-1] == 0)
    {
        end--;
    }
    if (end > nx - 2)
    {
        end = nx;
    }
    else
    {
        end++;
    }

    if (fp)
    {
        fprintf(fp, "Generating forces for table %d, boundary conditions: V''' at %g, %s at %g\n",
                table+1, start*h, end == nx ? "V'''" : "V'=0", (end-1)*h);
    }
    spline_forces(end-start, h, v+start, TRUE, end == nx, f+start);
}

static void read_tables(FILE *fp, const char *fn,
                        int ntab, int angle, t_tabledata td[])
{
    char    *libfn;
    char     buf[STRLEN];
    double **yy = NULL, start, end, dx0, dx1, ssd, vm, vp, f, numf;
    int      k, i, nx, nx0 = 0, ny, nny, ns;
    gmx_bool bAllZero, bZeroV, bZeroF;
    double   tabscale;

    nny   = 2*ntab+1;
    libfn = gmxlibfn(fn);
    nx    = read_xvg(libfn, &yy, &ny);
    if (ny != nny)
    {
        gmx_fatal(FARGS, "Trying to read file %s, but nr columns = %d, should be %d",
                  libfn, ny, nny);
    }
    if (angle == 0)
    {
        if (yy[0][0] != 0.0)
        {
            gmx_fatal(FARGS,
                      "The first distance in file %s is %f nm instead of %f nm",
                      libfn, yy[0][0], 0.0);
        }
    }
    else
    {
        if (angle == 1)
        {
            start = 0.0;
        }
        else
        {
            start = -180.0;
        }
        end = 180.0;
        if (yy[0][0] != start || yy[0][nx-1] != end)
        {
            gmx_fatal(FARGS, "The angles in file %s should go from %f to %f instead of %f to %f\n",
                      libfn, start, end, yy[0][0], yy[0][nx-1]);
        }
    }

    tabscale = (nx-1)/(yy[0][nx-1] - yy[0][0]);

    if (fp)
    {
        fprintf(fp, "Read user tables from %s with %d data points.\n", libfn, nx);
        if (angle == 0)
        {
            fprintf(fp, "Tabscale = %g points/nm\n", tabscale);
        }
    }

    bAllZero = TRUE;
    for (k = 0; k < ntab; k++)
    {
        bZeroV = TRUE;
        bZeroF = TRUE;
        for (i = 0; (i < nx); i++)
        {
            if (i >= 2)
            {
                dx0 = yy[0][i-1] - yy[0][i-2];
                dx1 = yy[0][i]   - yy[0][i-1];
                /* Check for 1% deviation in spacing */
                if (fabs(dx1 - dx0) >= 0.005*(fabs(dx0) + fabs(dx1)))
                {
                    gmx_fatal(FARGS, "In table file '%s' the x values are not equally spaced: %f %f %f", fn, yy[0][i-2], yy[0][i-1], yy[0][i]);
                }
            }
            if (yy[1+k*2][i] != 0)
            {
                bZeroV = FALSE;
                if (bAllZero)
                {
                    bAllZero = FALSE;
                    nx0      = i;
                }
                if (yy[1+k*2][i] >  0.01*GMX_REAL_MAX ||
                    yy[1+k*2][i] < -0.01*GMX_REAL_MAX)
                {
                    gmx_fatal(FARGS, "Out of range potential value %g in file '%s'",
                              yy[1+k*2][i], fn);
                }
            }
            if (yy[1+k*2+1][i] != 0)
            {
                bZeroF = FALSE;
                if (bAllZero)
                {
                    bAllZero = FALSE;
                    nx0      = i;
                }
                if (yy[1+k*2+1][i] >  0.01*GMX_REAL_MAX ||
                    yy[1+k*2+1][i] < -0.01*GMX_REAL_MAX)
                {
                    gmx_fatal(FARGS, "Out of range force value %g in file '%s'",
                              yy[1+k*2+1][i], fn);
                }
            }
        }

        if (!bZeroV && bZeroF)
        {
            set_forces(fp, angle, nx, 1/tabscale, yy[1+k*2], yy[1+k*2+1], k);
        }
        else
        {
            /* Check if the second column is close to minus the numerical
             * derivative of the first column.
             */
            ssd = 0;
            ns  = 0;
            for (i = 1; (i < nx-1); i++)
            {
                vm = yy[1+2*k][i-1];
                vp = yy[1+2*k][i+1];
                f  = yy[1+2*k+1][i];
                if (vm != 0 && vp != 0 && f != 0)
                {
                    /* Take the centered difference */
                    numf = -(vp - vm)*0.5*tabscale;
                    ssd += fabs(2*(f - numf)/(f + numf));
                    ns++;
                }
            }
            if (ns > 0)
            {
                ssd /= ns;
                sprintf(buf, "For the %d non-zero entries for table %d in %s the forces deviate on average %d%% from minus the numerical derivative of the potential\n", ns, k, libfn, (int)(100*ssd+0.5));
                if (debug)
                {
                    fprintf(debug, "%s", buf);
                }
                if (ssd > 0.2)
                {
                    if (fp)
                    {
                        fprintf(fp, "\nWARNING: %s\n", buf);
                    }
                    fprintf(stderr, "\nWARNING: %s\n", buf);
                }
            }
        }
    }
    if (bAllZero && fp)
    {
        fprintf(fp, "\nNOTE: All elements in table %s are zero\n\n", libfn);
    }

    for (k = 0; (k < ntab); k++)
    {
        init_table(nx, nx0, tabscale, &(td[k]), TRUE);
        for (i = 0; (i < nx); i++)
        {
            td[k].x[i] = yy[0][i];
            td[k].v[i] = yy[2*k+1][i];
            td[k].f[i] = yy[2*k+2][i];
        }
    }
    for (i = 0; (i < ny); i++)
    {
        sfree(yy[i]);
    }
    sfree(yy);
    sfree(libfn);
}

static void done_tabledata(t_tabledata *td)
{
    int i;

    if (!td)
    {
        return;
    }

    sfree(td->x);
    sfree(td->v);
    sfree(td->f);
}

static void fill_table(t_tabledata *td, int tp, const t_forcerec *fr,
                       gmx_bool b14only)
{
    /* Fill the table according to the formulas in the manual.
     * In principle, we only need the potential and the second
     * derivative, but then we would have to do lots of calculations
     * in the inner loop. By precalculating some terms (see manual)
     * we get better eventual performance, despite a larger table.
     *
     * Since some of these higher-order terms are very small,
     * we always use double precision to calculate them here, in order
     * to avoid unnecessary loss of precision.
     */
#ifdef DEBUG_SWITCH
    FILE    *fp;
#endif
    int      i;
    double   reppow, p;
    double   r1, rc, r12, r13;
    double   r, r2, r6, rc2, rc6, rc12;
    double   expr, Vtab, Ftab;
    /* Parameters for David's function */
    double   A = 0, B = 0, C = 0, A_3 = 0, B_4 = 0;
    /* Parameters for the switching function */
    double   ksw, swi, swi1;
    /* Temporary parameters */
    gmx_bool bPotentialSwitch, bForceSwitch, bPotentialShift;
    double   ewc   = fr->ewaldcoeff_q;
    double   ewclj = fr->ewaldcoeff_lj;
    double   Vcut  = 0;

    if (b14only)
    {
        bPotentialSwitch = FALSE;
        bForceSwitch     = FALSE;
        bPotentialShift  = FALSE;
    }
    else
    {
        bPotentialSwitch = ((tp == etabLJ6Switch) || (tp == etabLJ12Switch) ||
                            (tp == etabCOULSwitch) ||
                            (tp == etabEwaldSwitch) || (tp == etabEwaldUserSwitch) ||
                            (tprops[tp].bCoulomb && (fr->coulomb_modifier == eintmodPOTSWITCH)) ||
                            (!tprops[tp].bCoulomb && (fr->vdw_modifier == eintmodPOTSWITCH)));
        bForceSwitch  = ((tp == etabLJ6Shift) || (tp == etabLJ12Shift) ||
                         (tp == etabShift) ||
                         (tprops[tp].bCoulomb && (fr->coulomb_modifier == eintmodFORCESWITCH)) ||
                         (!tprops[tp].bCoulomb && (fr->vdw_modifier == eintmodFORCESWITCH)));
        bPotentialShift = ((tprops[tp].bCoulomb && (fr->coulomb_modifier == eintmodPOTSHIFT)) ||
                           (!tprops[tp].bCoulomb && (fr->vdw_modifier == eintmodPOTSHIFT)));
    }

    reppow = fr->reppow;

    if (tprops[tp].bCoulomb)
    {
        r1 = fr->rcoulomb_switch;
        rc = fr->rcoulomb;
    }
    else
    {
        r1 = fr->rvdw_switch;
        rc = fr->rvdw;
    }
    if (bPotentialSwitch)
    {
        ksw  = 1.0/(pow5(rc-r1));
    }
    else
    {
        ksw  = 0.0;
    }
    if (bForceSwitch)
    {
        if (tp == etabShift)
        {
            p = 1;
        }
        else if (tp == etabLJ6Shift)
        {
            p = 6;
        }
        else
        {
            p = reppow;
        }

        A = p * ((p+1)*r1-(p+4)*rc)/(pow(rc, p+2)*pow2(rc-r1));
        B = -p * ((p+1)*r1-(p+3)*rc)/(pow(rc, p+2)*pow3(rc-r1));
        C = 1.0/pow(rc, p)-A/3.0*pow3(rc-r1)-B/4.0*pow4(rc-r1);
        if (tp == etabLJ6Shift)
        {
            A = -A;
            B = -B;
            C = -C;
        }
        A_3 = A/3.0;
        B_4 = B/4.0;
    }
    if (debug)
    {
        fprintf(debug, "Setting up tables\n"); fflush(debug);
    }

#ifdef DEBUG_SWITCH
    fp = xvgropen("switch.xvg", "switch", "r", "s");
#endif

    if (bPotentialShift)
    {
        rc2   = rc*rc;
        rc6   = 1.0/(rc2*rc2*rc2);
        if (gmx_within_tol(reppow, 12.0, 10*GMX_DOUBLE_EPS))
        {
            rc12 = rc6*rc6;
        }
        else
        {
            rc12 = pow(rc, -reppow);
        }

        switch (tp)
        {
            case etabLJ6:
                /* Dispersion */
                Vcut = -rc6;
                break;
            case etabLJ6Ewald:
                Vcut  = -rc6*exp(-ewclj*ewclj*rc2)*(1 + ewclj*ewclj*rc2 + pow4(ewclj)*rc2*rc2/2);
                break;
            case etabLJ12:
                /* Repulsion */
                Vcut  = rc12;
                break;
            case etabCOUL:
                Vcut  = 1.0/rc;
                break;
            case etabEwald:
            case etabEwaldSwitch:
                Vtab  = gmx_erfc(ewc*rc)/rc;
                break;
            case etabEwaldUser:
                /* Only calculate minus the reciprocal space contribution */
                Vtab  = -gmx_erf(ewc*rc)/rc;
                break;
            case etabRF:
            case etabRF_ZERO:
                /* No need for preventing the usage of modifiers with RF */
                Vcut  = 0.0;
                break;
            case etabEXPMIN:
                Vcut  = exp(-rc);
                break;
            default:
                gmx_fatal(FARGS, "Cannot apply new potential-shift modifier to interaction type '%s' yet. (%s,%d)",
                          tprops[tp].name, __FILE__, __LINE__);
        }
    }

    for (i = td->nx0; (i < td->nx); i++)
    {
        r     = td->x[i];
        r2    = r*r;
        r6    = 1.0/(r2*r2*r2);
        if (gmx_within_tol(reppow, 12.0, 10*GMX_DOUBLE_EPS))
        {
            r12 = r6*r6;
        }
        else
        {
            r12 = pow(r, -reppow);
        }
        Vtab  = 0.0;
        Ftab  = 0.0;
        if (bPotentialSwitch)
        {
            /* swi is function, swi1 1st derivative and swi2 2nd derivative */
            /* The switch function is 1 for r<r1, 0 for r>rc, and smooth for
             * r1<=r<=rc. The 1st and 2nd derivatives are both zero at
             * r1 and rc.
             * ksw is just the constant 1/(rc-r1)^5, to save some calculations...
             */
            if (r <= r1)
            {
                swi  = 1.0;
                swi1 = 0.0;
            }
            else if (r >= rc)
            {
                swi  = 0.0;
                swi1 = 0.0;
            }
            else
            {
                swi      = 1 - 10*pow3(r-r1)*ksw*pow2(rc-r1)
                    + 15*pow4(r-r1)*ksw*(rc-r1) - 6*pow5(r-r1)*ksw;
                swi1     = -30*pow2(r-r1)*ksw*pow2(rc-r1)
                    + 60*pow3(r-r1)*ksw*(rc-r1) - 30*pow4(r-r1)*ksw;
            }
        }
        else /* not really needed, but avoids compiler warnings... */
        {
            swi  = 1.0;
            swi1 = 0.0;
        }
#ifdef DEBUG_SWITCH
        fprintf(fp, "%10g  %10g  %10g  %10g\n", r, swi, swi1, swi2);
#endif

        rc6 = rc*rc*rc;
        rc6 = 1.0/(rc6*rc6);

        switch (tp)
        {
            case etabLJ6:
                /* Dispersion */
                Vtab = -r6;
                Ftab = 6.0*Vtab/r;
                break;
            case etabLJ6Switch:
            case etabLJ6Shift:
                /* Dispersion */
                if (r < rc)
                {
                    Vtab = -r6;
                    Ftab = 6.0*Vtab/r;
                    break;
                }
                break;
            case etabLJ12:
                /* Repulsion */
                Vtab  = r12;
                Ftab  = reppow*Vtab/r;
                break;
            case etabLJ12Switch:
            case etabLJ12Shift:
                /* Repulsion */
                if (r < rc)
                {
                    Vtab  = r12;
                    Ftab  = reppow*Vtab/r;
                }
                break;
            case etabLJ6Encad:
                if (r < rc)
                {
                    Vtab  = -(r6-6.0*(rc-r)*rc6/rc-rc6);
                    Ftab  = -(6.0*r6/r-6.0*rc6/rc);
                }
                else /* r>rc */
                {
                    Vtab  = 0;
                    Ftab  = 0;
                }
                break;
            case etabLJ12Encad:
                if (r < rc)
                {
                    Vtab  = -(r6-6.0*(rc-r)*rc6/rc-rc6);
                    Ftab  = -(6.0*r6/r-6.0*rc6/rc);
                }
                else /* r>rc */
                {
                    Vtab  = 0;
                    Ftab  = 0;
                }
                break;
            case etabCOUL:
                Vtab  = 1.0/r;
                Ftab  = 1.0/r2;
                break;
            case etabCOULSwitch:
            case etabShift:
                if (r < rc)
                {
                    Vtab  = 1.0/r;
                    Ftab  = 1.0/r2;
                }
                break;
            case etabEwald:
            case etabEwaldSwitch:
                Vtab  = gmx_erfc(ewc*r)/r;
                Ftab  = gmx_erfc(ewc*r)/r2+exp(-(ewc*ewc*r2))*ewc*M_2_SQRTPI/r;
                break;
            case etabEwaldUser:
            case etabEwaldUserSwitch:
                /* Only calculate the negative of the reciprocal space contribution */
                Vtab  = -gmx_erf(ewc*r)/r;
                Ftab  = -gmx_erf(ewc*r)/r2+exp(-(ewc*ewc*r2))*ewc*M_2_SQRTPI/r;
                break;
            case etabLJ6Ewald:
                Vtab  = -r6*exp(-ewclj*ewclj*r2)*(1 + ewclj*ewclj*r2 + pow4(ewclj)*r2*r2/2);
                Ftab  = 6.0*Vtab/r - r6*exp(-ewclj*ewclj*r2)*pow5(ewclj)*ewclj*r2*r2*r;
                break;
            case etabRF:
            case etabRF_ZERO:
                Vtab  = 1.0/r      +   fr->k_rf*r2 - fr->c_rf;
                Ftab  = 1.0/r2     - 2*fr->k_rf*r;
                if (tp == etabRF_ZERO && r >= rc)
                {
                    Vtab = 0;
                    Ftab = 0;
                }
                break;
            case etabEXPMIN:
                expr  = exp(-r);
                Vtab  = expr;
                Ftab  = expr;
                break;
            case etabCOULEncad:
                if (r < rc)
                {
                    Vtab  = 1.0/r-(rc-r)/(rc*rc)-1.0/rc;
                    Ftab  = 1.0/r2-1.0/(rc*rc);
                }
                else /* r>rc */
                {
                    Vtab  = 0;
                    Ftab  = 0;
                }
                break;
            default:
                gmx_fatal(FARGS, "Table type %d not implemented yet. (%s,%d)",
                          tp, __FILE__, __LINE__);
        }
        if (bForceSwitch)
        {
            /* Normal coulomb with cut-off correction for potential */
            if (r < rc)
            {
                Vtab -= C;
                /* If in Shifting range add something to it */
                if (r > r1)
                {
                    r12    = (r-r1)*(r-r1);
                    r13    = (r-r1)*r12;
                    Vtab  += -A_3*r13 - B_4*r12*r12;
                    Ftab  +=   A*r12 + B*r13;
                }
            }
            else
            {
                /* Make sure interactions are zero outside cutoff with modifiers */
                Vtab = 0;
                Ftab = 0;
            }
        }
        if (bPotentialShift)
        {
            if (r < rc)
            {
                Vtab -= Vcut;
            }
            else
            {
                /* Make sure interactions are zero outside cutoff with modifiers */
                Vtab = 0;
                Ftab = 0;
            }
        }

        if (ETAB_USER(tp))
        {
            Vtab += td->v[i];
            Ftab += td->f[i];
        }

        if (bPotentialSwitch)
        {
            if (r >= rc)
            {
                /* Make sure interactions are zero outside cutoff with modifiers */
                Vtab = 0;
                Ftab = 0;
            }
            else if (r > r1)
            {
                Ftab = Ftab*swi - Vtab*swi1;
                Vtab = Vtab*swi;
            }
        }
        /* Convert to single precision when we store to mem */
        td->v[i]  = Vtab;
        td->f[i]  = Ftab;
    }

    /* Continue the table linearly from nx0 to 0.
     * These values are only required for energy minimization with overlap or TPI.
     */
    for (i = td->nx0-1; i >= 0; i--)
    {
        td->v[i] = td->v[i+1] + td->f[i+1]*(td->x[i+1] - td->x[i]);
        td->f[i] = td->f[i+1];
    }

#ifdef DEBUG_SWITCH
    xvgrclose(fp);
#endif
}

static void set_table_type(int tabsel[], const t_forcerec *fr, gmx_bool b14only)
{
    int eltype, vdwtype;

    /* Set the different table indices.
     * Coulomb first.
     */


    if (b14only)
    {
        switch (fr->eeltype)
        {
            case eelRF_NEC:
                eltype = eelRF;
                break;
            case eelUSER:
            case eelPMEUSER:
            case eelPMEUSERSWITCH:
                eltype = eelUSER;
                break;
            default:
                eltype = eelCUT;
        }
    }
    else
    {
        eltype = fr->eeltype;
    }

    switch (eltype)
    {
        case eelCUT:
            tabsel[etiCOUL] = etabCOUL;
            break;
        case eelPOISSON:
            tabsel[etiCOUL] = etabShift;
            break;
        case eelSHIFT:
            if (fr->rcoulomb > fr->rcoulomb_switch)
            {
                tabsel[etiCOUL] = etabShift;
            }
            else
            {
                tabsel[etiCOUL] = etabCOUL;
            }
            break;
        case eelEWALD:
        case eelPME:
        case eelP3M_AD:
            tabsel[etiCOUL] = etabEwald;
            break;
        case eelPMESWITCH:
            tabsel[etiCOUL] = etabEwaldSwitch;
            break;
        case eelPMEUSER:
            tabsel[etiCOUL] = etabEwaldUser;
            break;
        case eelPMEUSERSWITCH:
            tabsel[etiCOUL] = etabEwaldUserSwitch;
            break;
        case eelRF:
        case eelGRF:
        case eelRF_NEC:
            tabsel[etiCOUL] = etabRF;
            break;
        case eelRF_ZERO:
            tabsel[etiCOUL] = etabRF_ZERO;
            break;
        case eelSWITCH:
            tabsel[etiCOUL] = etabCOULSwitch;
            break;
        case eelUSER:
            tabsel[etiCOUL] = etabUSER;
            break;
        case eelENCADSHIFT:
            tabsel[etiCOUL] = etabCOULEncad;
            break;
        default:
            gmx_fatal(FARGS, "Invalid eeltype %d", eltype);
    }

    /* Van der Waals time */
    if (fr->bBHAM && !b14only)
    {
        tabsel[etiLJ6]  = etabLJ6;
        tabsel[etiLJ12] = etabEXPMIN;
    }
    else
    {
        if (b14only && fr->vdwtype != evdwUSER)
        {
            vdwtype = evdwCUT;
        }
        else
        {
            vdwtype = fr->vdwtype;
        }

        switch (vdwtype)
        {
            case evdwSWITCH:
                tabsel[etiLJ6]  = etabLJ6Switch;
                tabsel[etiLJ12] = etabLJ12Switch;
                break;
            case evdwSHIFT:
                tabsel[etiLJ6]  = etabLJ6Shift;
                tabsel[etiLJ12] = etabLJ12Shift;
                break;
            case evdwUSER:
                tabsel[etiLJ6]  = etabUSER;
                tabsel[etiLJ12] = etabUSER;
                break;
            case evdwCUT:
                tabsel[etiLJ6]  = etabLJ6;
                tabsel[etiLJ12] = etabLJ12;
                break;
            case evdwENCADSHIFT:
                tabsel[etiLJ6]  = etabLJ6Encad;
                tabsel[etiLJ12] = etabLJ12Encad;
                break;
            case evdwPME:
                tabsel[etiLJ6]  = etabLJ6Ewald;
                tabsel[etiLJ12] = etabLJ12;
                break;
            default:
                gmx_fatal(FARGS, "Invalid vdwtype %d in %s line %d", vdwtype,
                          __FILE__, __LINE__);
        }

        if (!b14only && fr->vdw_modifier != eintmodNONE)
        {
            if (fr->vdw_modifier != eintmodPOTSHIFT &&
                fr->vdwtype != evdwCUT)
            {
                gmx_incons("Potential modifiers other than potential-shift are only implemented for LJ cut-off");
            }

            /* LJ-PME and other (shift-only) modifiers are handled by applying the modifiers
             * to the original interaction forms when we fill the table, so we only check cutoffs here.
             */
            if (fr->vdwtype == evdwCUT)
            {
                switch (fr->vdw_modifier)
                {
                    case eintmodNONE:
                    case eintmodPOTSHIFT:
                    case eintmodEXACTCUTOFF:
                        /* No modification */
                        break;
                    case eintmodPOTSWITCH:
                        tabsel[etiLJ6]  = etabLJ6Switch;
                        tabsel[etiLJ12] = etabLJ12Switch;
                        break;
                    case eintmodFORCESWITCH:
                        tabsel[etiLJ6]  = etabLJ6Shift;
                        tabsel[etiLJ12] = etabLJ12Shift;
                        break;
                    default:
                        gmx_incons("Unsupported vdw_modifier");
                }
            }
        }
    }
}

t_forcetable make_tables(FILE *out, const output_env_t oenv,
                         const t_forcerec *fr,
                         gmx_bool bVerbose, const char *fn,
                         real rtab, int flags)
{
    const char     *fns[3]   = { "ctab.xvg", "dtab.xvg", "rtab.xvg" };
    const char     *fns14[3] = { "ctab14.xvg", "dtab14.xvg", "rtab14.xvg" };
    FILE           *fp;
    t_tabledata    *td;
    gmx_bool        b14only, bReadTab, bGenTab;
    real            x0, y0, yp;
    int             i, j, k, nx, nx0, tabsel[etiNR];
    real            scalefactor;

    t_forcetable    table;

    b14only = (flags & GMX_MAKETABLES_14ONLY);

    if (flags & GMX_MAKETABLES_FORCEUSER)
    {
        tabsel[etiCOUL] = etabUSER;
        tabsel[etiLJ6]  = etabUSER;
        tabsel[etiLJ12] = etabUSER;
    }
    else
    {
        set_table_type(tabsel, fr, b14only);
    }
    snew(td, etiNR);
    table.r         = rtab;
    table.scale     = 0;
    table.n         = 0;
    table.scale_exp = 0;
    nx0             = 10;
    nx              = 0;

    table.interaction   = GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP;
    table.format        = GMX_TABLE_FORMAT_CUBICSPLINE_YFGH;
    table.formatsize    = 4;
    table.ninteractions = 3;
    table.stride        = table.formatsize*table.ninteractions;

    /* Check whether we have to read or generate */
    bReadTab = FALSE;
    bGenTab  = FALSE;
    for (i = 0; (i < etiNR); i++)
    {
        if (ETAB_USER(tabsel[i]))
        {
            bReadTab = TRUE;
        }
        if (tabsel[i] != etabUSER)
        {
            bGenTab  = TRUE;
        }
    }
    if (bReadTab)
    {
        read_tables(out, fn, etiNR, 0, td);
        if (rtab == 0 || (flags & GMX_MAKETABLES_14ONLY))
        {
            rtab      = td[0].x[td[0].nx-1];
            table.n   = td[0].nx;
            nx        = table.n;
        }
        else
        {
            if (td[0].x[td[0].nx-1] < rtab)
            {
                gmx_fatal(FARGS, "Tables in file %s not long enough for cut-off:\n"
                          "\tshould be at least %f nm\n", fn, rtab);
            }
            nx        = table.n = (int)(rtab*td[0].tabscale + 0.5);
        }
        table.scale = td[0].tabscale;
        nx0         = td[0].nx0;
    }
    if (bGenTab)
    {
        if (!bReadTab)
        {
#ifdef GMX_DOUBLE
            table.scale = 2000.0;
#else
            table.scale = 500.0;
#endif
            nx = table.n = rtab*table.scale;
        }
    }
    if (fr->bBHAM)
    {
        if (fr->bham_b_max != 0)
        {
            table.scale_exp = table.scale/fr->bham_b_max;
        }
        else
        {
            table.scale_exp = table.scale;
        }
    }

    /* Each table type (e.g. coul,lj6,lj12) requires four
     * numbers per nx+1 data points. For performance reasons we want
     * the table data to be aligned to 16-byte.
     */
    snew_aligned(table.data, 12*(nx+1)*sizeof(real), 32);

    for (k = 0; (k < etiNR); k++)
    {
        if (tabsel[k] != etabUSER)
        {
            init_table(nx, nx0,
                       (tabsel[k] == etabEXPMIN) ? table.scale_exp : table.scale,
                       &(td[k]), !bReadTab);
            fill_table(&(td[k]), tabsel[k], fr, b14only);
            if (out)
            {
                fprintf(out, "%s table with %d data points for %s%s.\n"
                        "Tabscale = %g points/nm\n",
                        ETAB_USER(tabsel[k]) ? "Modified" : "Generated",
                        td[k].nx, b14only ? "1-4 " : "", tprops[tabsel[k]].name,
                        td[k].tabscale);
            }
        }

        /* Set scalefactor for c6/c12 tables. This is because we save flops in the non-table kernels
         * by including the derivative constants (6.0 or 12.0) in the parameters, since
         * we no longer calculate force in most steps. This means the c6/c12 parameters
         * have been scaled up, so we need to scale down the table interactions too.
         * It comes here since we need to scale user tables too.
         */
        if (k == etiLJ6)
        {
            scalefactor = 1.0/6.0;
        }
        else if (k == etiLJ12 && tabsel[k] != etabEXPMIN)
        {
            scalefactor = 1.0/12.0;
        }
        else
        {
            scalefactor = 1.0;
        }

        copy2table(table.n, k*4, 12, td[k].x, td[k].v, td[k].f, scalefactor, table.data);

        if (bDebugMode() && bVerbose)
        {
            if (b14only)
            {
                fp = xvgropen(fns14[k], fns14[k], "r", "V", oenv);
            }
            else
            {
                fp = xvgropen(fns[k], fns[k], "r", "V", oenv);
            }
            /* plot the output 5 times denser than the table data */
            for (i = 5*((nx0+1)/2); i < 5*table.n; i++)
            {
                x0 = i*table.r/(5*(table.n-1));
                evaluate_table(table.data, 4*k, 12, table.scale, x0, &y0, &yp);
                fprintf(fp, "%15.10e  %15.10e  %15.10e\n", x0, y0, yp);
            }
            xvgrclose(fp);
        }
        done_tabledata(&(td[k]));
    }
    sfree(td);

    return table;
}

t_forcetable make_gb_table(const output_env_t oenv,
                           const t_forcerec  *fr)
{
    const char     *fns[3]   = { "gbctab.xvg", "gbdtab.xvg", "gbrtab.xvg" };
    const char     *fns14[3] = { "gbctab14.xvg", "gbdtab14.xvg", "gbrtab14.xvg" };
    FILE           *fp;
    t_tabledata    *td;
    gmx_bool        bReadTab, bGenTab;
    real            x0, y0, yp;
    int             i, j, k, nx, nx0, tabsel[etiNR];
    double          r, r2, Vtab, Ftab, expterm;

    t_forcetable    table;

    double          abs_error_r, abs_error_r2;
    double          rel_error_r, rel_error_r2;
    double          rel_error_r_old = 0, rel_error_r2_old = 0;
    double          x0_r_error, x0_r2_error;


    /* Only set a Coulomb table for GB */
    /*
       tabsel[0]=etabGB;
       tabsel[1]=-1;
       tabsel[2]=-1;
     */

    /* Set the table dimensions for GB, not really necessary to
     * use etiNR (since we only have one table, but ...)
     */
    snew(td, 1);
    table.interaction   = GMX_TABLE_INTERACTION_ELEC;
    table.format        = GMX_TABLE_FORMAT_CUBICSPLINE_YFGH;
    table.r             = fr->gbtabr;
    table.scale         = fr->gbtabscale;
    table.scale_exp     = 0;
    table.n             = table.scale*table.r;
    table.formatsize    = 4;
    table.ninteractions = 1;
    table.stride        = table.formatsize*table.ninteractions;
    nx0                 = 0;
    nx                  = table.scale*table.r;

    /* Check whether we have to read or generate
     * We will always generate a table, so remove the read code
     * (Compare with original make_table function
     */
    bReadTab = FALSE;
    bGenTab  = TRUE;

    /* Each table type (e.g. coul,lj6,lj12) requires four
     * numbers per datapoint. For performance reasons we want
     * the table data to be aligned to 16-byte. This is accomplished
     * by allocating 16 bytes extra to a temporary pointer, and then
     * calculating an aligned pointer. This new pointer must not be
     * used in a free() call, but thankfully we're sloppy enough not
     * to do this :-)
     */

    snew_aligned(table.data, 4*nx, 32);

    init_table(nx, nx0, table.scale, &(td[0]), !bReadTab);

    /* Local implementation so we don't have to use the etabGB
     * enum above, which will cause problems later when
     * making the other tables (right now even though we are using
     * GB, the normal Coulomb tables will be created, but this
     * will cause a problem since fr->eeltype==etabGB which will not
     * be defined in fill_table and set_table_type
     */

    for (i = nx0; i < nx; i++)
    {
        r       = td->x[i];
        r2      = r*r;
        expterm = exp(-0.25*r2);

        Vtab = 1/sqrt(r2+expterm);
        Ftab = (r-0.25*r*expterm)/((r2+expterm)*sqrt(r2+expterm));

        /* Convert to single precision when we store to mem */
        td->v[i]  = Vtab;
        td->f[i]  = Ftab;

    }

    copy2table(table.n, 0, 4, td[0].x, td[0].v, td[0].f, 1.0, table.data);

    if (bDebugMode())
    {
        fp = xvgropen(fns[0], fns[0], "r", "V", oenv);
        /* plot the output 5 times denser than the table data */
        /* for(i=5*nx0;i<5*table.n;i++) */
        for (i = nx0; i < table.n; i++)
        {
            /* x0=i*table.r/(5*table.n); */
            x0 = i*table.r/table.n;
            evaluate_table(table.data, 0, 4, table.scale, x0, &y0, &yp);
            fprintf(fp, "%15.10e  %15.10e  %15.10e\n", x0, y0, yp);

        }
        xvgrclose(fp);
    }

    /*
       for(i=100*nx0;i<99.81*table.n;i++)
       {
       r = i*table.r/(100*table.n);
       r2      = r*r;
       expterm = exp(-0.25*r2);

       Vtab = 1/sqrt(r2+expterm);
       Ftab = (r-0.25*r*expterm)/((r2+expterm)*sqrt(r2+expterm));


       evaluate_table(table.data,0,4,table.scale,r,&y0,&yp);
       printf("gb: i=%d, x0=%g, y0=%15.15f, Vtab=%15.15f, yp=%15.15f, Ftab=%15.15f\n",i,r, y0, Vtab, yp, Ftab);

       abs_error_r=fabs(y0-Vtab);
       abs_error_r2=fabs(yp-(-1)*Ftab);

       rel_error_r=abs_error_r/y0;
       rel_error_r2=fabs(abs_error_r2/yp);


       if(rel_error_r>rel_error_r_old)
       {
       rel_error_r_old=rel_error_r;
       x0_r_error=x0;
       }

       if(rel_error_r2>rel_error_r2_old)
       {
       rel_error_r2_old=rel_error_r2;
       x0_r2_error=x0;
       }
       }

       printf("gb: MAX REL ERROR IN R=%15.15f, MAX REL ERROR IN R2=%15.15f\n",rel_error_r_old, rel_error_r2_old);
       printf("gb: XO_R=%g, X0_R2=%g\n",x0_r_error, x0_r2_error);

       exit(1); */
    done_tabledata(&(td[0]));
    sfree(td);

    return table;


}

t_forcetable make_atf_table(FILE *out, const output_env_t oenv,
                            const t_forcerec *fr,
                            const char *fn,
                            matrix box)
{
    const char  *fns[3] = { "tf_tab.xvg", "atfdtab.xvg", "atfrtab.xvg" };
    FILE        *fp;
    t_tabledata *td;
    real         x0, y0, yp, rtab;
    int          i, nx, nx0;
    real         rx, ry, rz, box_r;

    t_forcetable table;


    /* Set the table dimensions for ATF, not really necessary to
     * use etiNR (since we only have one table, but ...)
     */
    snew(td, 1);

    if (fr->adress_type == eAdressSphere)
    {
        /* take half box diagonal direction as tab range */
        rx    = 0.5*box[0][0]+0.5*box[1][0]+0.5*box[2][0];
        ry    = 0.5*box[0][1]+0.5*box[1][1]+0.5*box[2][1];
        rz    = 0.5*box[0][2]+0.5*box[1][2]+0.5*box[2][2];
        box_r = sqrt(rx*rx+ry*ry+rz*rz);

    }
    else
    {
        /* xsplit: take half box x direction as tab range */
        box_r        = box[0][0]/2;
    }
    table.r         = box_r;
    table.scale     = 0;
    table.n         = 0;
    table.scale_exp = 0;
    nx0             = 10;
    nx              = 0;

    read_tables(out, fn, 1, 0, td);
    rtab      = td[0].x[td[0].nx-1];

    if (fr->adress_type == eAdressXSplit && (rtab < box[0][0]/2))
    {
        gmx_fatal(FARGS, "AdResS full box therm force table in file %s extends to %f:\n"
                  "\tshould extend to at least half the length of the box in x-direction"
                  "%f\n", fn, rtab, box[0][0]/2);
    }
    if (rtab < box_r)
    {
        gmx_fatal(FARGS, "AdResS full box therm force table in file %s extends to %f:\n"
                  "\tshould extend to at least for spherical adress"
                  "%f (=distance from center to furthermost point in box \n", fn, rtab, box_r);
    }


    table.n     = td[0].nx;
    nx          = table.n;
    table.scale = td[0].tabscale;
    nx0         = td[0].nx0;

    /* Each table type (e.g. coul,lj6,lj12) requires four
     * numbers per datapoint. For performance reasons we want
     * the table data to be aligned to 16-byte. This is accomplished
     * by allocating 16 bytes extra to a temporary pointer, and then
     * calculating an aligned pointer. This new pointer must not be
     * used in a free() call, but thankfully we're sloppy enough not
     * to do this :-)
     */

    snew_aligned(table.data, 4*nx, 32);

    copy2table(table.n, 0, 4, td[0].x, td[0].v, td[0].f, 1.0, table.data);

    if (bDebugMode())
    {
        fp = xvgropen(fns[0], fns[0], "r", "V", oenv);
        /* plot the output 5 times denser than the table data */
        /* for(i=5*nx0;i<5*table.n;i++) */

        for (i = 5*((nx0+1)/2); i < 5*table.n; i++)
        {
            /* x0=i*table.r/(5*table.n); */
            x0 = i*table.r/(5*(table.n-1));
            evaluate_table(table.data, 0, 4, table.scale, x0, &y0, &yp);
            fprintf(fp, "%15.10e  %15.10e  %15.10e\n", x0, y0, yp);

        }
        xvgrclose(fp);
    }

    done_tabledata(&(td[0]));
    sfree(td);

    table.interaction   = GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP;
    table.format        = GMX_TABLE_FORMAT_CUBICSPLINE_YFGH;
    table.formatsize    = 4;
    table.ninteractions = 3;
    table.stride        = table.formatsize*table.ninteractions;


    return table;
}

bondedtable_t make_bonded_table(FILE *fplog, char *fn, int angle)
{
    t_tabledata   td;
    double        start;
    int           i;
    bondedtable_t tab;

    if (angle < 2)
    {
        start = 0;
    }
    else
    {
        start = -180.0;
    }
    read_tables(fplog, fn, 1, angle, &td);
    if (angle > 0)
    {
        /* Convert the table from degrees to radians */
        for (i = 0; i < td.nx; i++)
        {
            td.x[i] *= DEG2RAD;
            td.f[i] *= RAD2DEG;
        }
        td.tabscale *= RAD2DEG;
    }
    tab.n     = td.nx;
    tab.scale = td.tabscale;
    snew(tab.data, tab.n*4);
    copy2table(tab.n, 0, 4, td.x, td.v, td.f, 1.0, tab.data);
    done_tabledata(&td);

    return tab;
}
