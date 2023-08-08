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

#include "forcetable.h"

#include <cmath>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

/* All the possible (implemented) table functions */
enum
{
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
    etabEXPMIN,
    etabUSER,
    etabNR
};

/** Evaluates to true if the table type contains user data. */
#define ETAB_USER(e) ((e) == etabUSER || (e) == etabEwaldUser || (e) == etabEwaldUserSwitch)

typedef struct
{
    const char* name;
    gmx_bool    bCoulomb;
} t_tab_props;

/* This structure holds name and a flag that tells whether
   this is a Coulomb type funtion */
static const t_tab_props tprops[etabNR] = {
    { "LJ6", FALSE },         { "LJ12", FALSE },      { "LJ6Shift", FALSE },
    { "LJ12Shift", FALSE },   { "Shift", TRUE },      { "RF", TRUE },
    { "RF-zero", TRUE },      { "COUL", TRUE },       { "Ewald", TRUE },
    { "Ewald-Switch", TRUE }, { "Ewald-User", TRUE }, { "Ewald-User-Switch", TRUE },
    { "LJ6Ewald", FALSE },    { "LJ6Switch", FALSE }, { "LJ12Switch", FALSE },
    { "COULSwitch", TRUE },   { "EXPMIN", FALSE },    { "USER", FALSE },
};

struct t_tabledata
{
    t_tabledata() = default;
    t_tabledata(int n, int firstTableNum, double scale, bool bAlloc);
    int                 nx;
    int                 nx0;
    double              tabscale;
    std::vector<double> x;
    std::vector<double> v;
    std::vector<double> f;
};

double v_q_ewald_lr(double beta, double r)
{
    if (r == 0)
    {
        return beta * 2 / std::sqrt(M_PI);
    }
    else
    {
        return std::erf(beta * r) / r;
    }
}

double v_lj_ewald_lr(double beta, double r)
{
    double br, br2, br4, r6, factor;
    if (r == 0)
    {
        return gmx::power6(beta) / 6;
    }
    else
    {
        br     = beta * r;
        br2    = br * br;
        br4    = br2 * br2;
        r6     = gmx::power6(r);
        factor = (1.0 - std::exp(-br2) * (1 + br2 + 0.5 * br4)) / r6;
        return factor;
    }
}

EwaldCorrectionTables generateEwaldCorrectionTables(const int    numPoints,
                                                    const double tableScaling,
                                                    const real   beta,
                                                    real_space_grid_contribution_computer v_lr)
{
    real     tab_max;
    int      i, i_inrange;
    double   dc, dc_new;
    gmx_bool bOutOfRange;
    double   v_r0, v_r1, v_inrange, vi, a0, a1, a2dx;
    double   x_r0;

    /* This function is called using either v_ewald_lr or v_lj_ewald_lr as a function argument
     * depending on whether we should create electrostatic or Lennard-Jones Ewald tables.
     */

    if (numPoints < 2)
    {
        gmx_fatal(FARGS, "Can not make a spline table with less than 2 points");
    }

    const double dx = 1 / tableScaling;

    EwaldCorrectionTables tables;
    tables.scale = tableScaling;
    tables.tableF.resize(numPoints);
    tables.tableV.resize(numPoints);
    tables.tableFDV0.resize(numPoints * 4);
    gmx::ArrayRef<real> table_f    = tables.tableF;
    gmx::ArrayRef<real> table_v    = tables.tableV;
    gmx::ArrayRef<real> table_fdv0 = tables.tableFDV0;

    /* We need some margin to be able to divide table values by r
     * in the kernel and also to do the integration arithmetics
     * without going out of range. Furthemore, we divide by dx below.
     */
    tab_max = GMX_REAL_MAX * 0.0001;

    /* This function produces a table with:
     * maximum energy error: V'''/(6*12*sqrt(3))*dx^3
     * maximum force error:  V'''/(6*4)*dx^2
     * The rms force error is the max error times 1/sqrt(5)=0.45.
     */

    bOutOfRange = FALSE;
    i_inrange   = numPoints;
    v_inrange   = 0;
    dc          = 0;
    for (i = numPoints - 1; i >= 0; i--)
    {
        x_r0 = i * dx;

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
            vi = v_inrange - dc * (i - i_inrange) * dx;
        }

        table_v[i] = vi;

        if (i == 0)
        {
            continue;
        }

        /* Get the potential at table point i-1 */
        v_r1 = (*v_lr)(beta, (i - 1) * dx);

        if (v_r1 != v_r1 || v_r1 < -tab_max || v_r1 > tab_max)
        {
            bOutOfRange = TRUE;
        }

        if (!bOutOfRange)
        {
            /* Calculate the average second derivative times dx over interval i-1 to i.
             * Using the function values at the end points and in the middle.
             */
            a2dx = (v_r0 + v_r1 - 2 * (*v_lr)(beta, x_r0 - 0.5 * dx)) / (0.25 * dx);
            /* Set the derivative of the spline to match the difference in potential
             * over the interval plus the average effect of the quadratic term.
             * This is the essential step for minimizing the error in the force.
             */
            dc = (v_r0 - v_r1) / dx + 0.5 * a2dx;
        }

        if (i == numPoints - 1)
        {
            /* Fill the table with the force, minus the derivative of the spline */
            table_f[i] = -dc;
        }
        else
        {
            /* tab[i] will contain the average of the splines over the two intervals */
            table_f[i] += -0.5 * dc;
        }

        if (!bOutOfRange)
        {
            /* Make spline s(x) = a0 + a1*(x - xr) + 0.5*a2*(x - xr)^2
             * matching the potential at the two end points
             * and the derivative dc at the end point xr.
             */
            a0   = v_r0;
            a1   = dc;
            a2dx = (a1 * dx + v_r1 - a0) * 2 / dx;

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

        table_f[(i - 1)] = -0.5 * dc;
    }
    /* Currently the last value only contains half the force: double it */
    table_f[0] *= 2;

    if (!table_fdv0.empty())
    {
        /* Copy to FDV0 table too. Allocation occurs in forcerec.c,
         * init_ewald_f_table().
         */
        for (i = 0; i < numPoints - 1; i++)
        {
            table_fdv0[4 * i]     = table_f[i];
            table_fdv0[4 * i + 1] = table_f[i + 1] - table_f[i];
            table_fdv0[4 * i + 2] = table_v[i];
            table_fdv0[4 * i + 3] = 0.0;
        }
        const int lastPoint           = numPoints - 1;
        table_fdv0[4 * lastPoint]     = table_f[lastPoint];
        table_fdv0[4 * lastPoint + 1] = -table_f[lastPoint];
        table_fdv0[4 * lastPoint + 2] = table_v[lastPoint];
        table_fdv0[4 * lastPoint + 3] = 0.0;
    }

    return tables;
}

/* Returns the spacing for a function using the maximum of
 * the third derivative, x_scale (unit 1/length)
 * and function tolerance.
 */
static double spline3_table_scale(double third_deriv_max, double x_scale, double func_tol)
{
    double deriv_tol;
    double sc_deriv, sc_func;

    /* Force tolerance: single precision accuracy */
    deriv_tol = GMX_FLOAT_EPS;
    sc_deriv  = std::sqrt(third_deriv_max / (6 * 4 * deriv_tol * x_scale)) * x_scale;

    /* Don't try to be more accurate on energy than the precision */
    func_tol = std::max(func_tol, static_cast<double>(GMX_REAL_EPS));
    sc_func  = std::cbrt(third_deriv_max / (6 * 12 * std::sqrt(3.0) * func_tol)) * x_scale;

    return std::max(sc_deriv, sc_func);
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
real ewald_spline3_table_scale(const interaction_const_t& ic,
                               const bool                 generateCoulombTables,
                               const bool                 generateVdwTables)
{
    GMX_RELEASE_ASSERT(!generateCoulombTables || usingPmeOrEwald(ic.eeltype),
                       "Can only use tables with Ewald");
    GMX_RELEASE_ASSERT(!generateVdwTables || usingLJPme(ic.vdwtype),
                       "Can only use tables with Ewald");

    real sc = 0;

    if (generateCoulombTables)
    {
        GMX_RELEASE_ASSERT(ic.ewaldcoeff_q > 0, "The Ewald coefficient should be positive");

        double erf_x_d3 = 1.0522; /* max of (erf(x)/x)''' */
        double etol;
        real   sc_q;

        /* Energy tolerance: 0.1 times the cut-off jump */
        etol = 0.1 * std::erfc(ic.ewaldcoeff_q * ic.rcoulomb);

        sc_q = spline3_table_scale(erf_x_d3, ic.ewaldcoeff_q, etol);

        if (debug)
        {
            fprintf(debug, "Ewald Coulomb quadratic spline table spacing: %f nm\n", 1 / sc_q);
        }

        sc = std::max(sc, sc_q);
    }

    if (generateVdwTables)
    {
        GMX_RELEASE_ASSERT(ic.ewaldcoeff_lj > 0, "The Ewald coefficient should be positive");

        double func_d3 = 0.42888; /* max of (x^-6 (1 - exp(-x^2)(1+x^2+x^4/2)))''' */
        double xrc2, etol;
        real   sc_lj;

        /* Energy tolerance: 0.1 times the cut-off jump */
        xrc2 = gmx::square(ic.ewaldcoeff_lj * ic.rvdw);
        etol = 0.1 * std::exp(-xrc2) * (1 + xrc2 + xrc2 * xrc2 / 2.0);

        sc_lj = spline3_table_scale(func_d3, ic.ewaldcoeff_lj, etol);

        if (debug)
        {
            fprintf(debug, "Ewald LJ quadratic spline table spacing: %f nm\n", 1 / sc_lj);
        }

        sc = std::max(sc, sc_lj);
    }

    return sc;
}

static void copy2table(int                         n,
                       int                         offset,
                       int                         stride,
                       gmx::ArrayRef<const double> x,
                       gmx::ArrayRef<const double> Vtab,
                       gmx::ArrayRef<const double> Ftab,
                       real                        scalefactor,
                       gmx::ArrayRef<real>         dest)
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
        if (i < n - 1)
        {
            h = x[i + 1] - x[i];
            F = -Ftab[i] * h;
            G = 3 * (Vtab[i + 1] - Vtab[i]) + (Ftab[i + 1] + 2 * Ftab[i]) * h;
            H = -2 * (Vtab[i + 1] - Vtab[i]) - (Ftab[i + 1] + Ftab[i]) * h;
        }
        else
        {
            /* Fill the last entry with a linear potential,
             * this is mainly for rounding issues with angle and dihedral potentials.
             */
            F = -Ftab[i] * h;
            G = 0;
            H = 0;
        }
        nn0           = offset + i * stride;
        dest[nn0]     = scalefactor * Vtab[i];
        dest[nn0 + 1] = scalefactor * F;
        dest[nn0 + 2] = scalefactor * G;
        dest[nn0 + 3] = scalefactor * H;
    }
}

t_tabledata::t_tabledata(int n, int firstTableNum, double scale, bool bAlloc) :
    nx(n), nx0(firstTableNum), tabscale(scale)
{
    if (bAlloc)
    {
        x.resize(nx);
        v.resize(nx);
        f.resize(nx);
    }
}

static void spline_forces(int nx, double h, const double v[], gmx_bool bS3, gmx_bool bE3, double f[])
{
    int    start, end, i;
    double v3, b_s, b_e, b;
    double beta;

    /* Formulas can be found in:
     * H.J.C. Berendsen, Simulating the Physical World, Cambridge 2007
     */

    if (nx < 4 && (bS3 || bE3))
    {
        gmx_fatal(FARGS,
                  "Can not generate splines with third derivative boundary conditions with less "
                  "than 4 (%d) points",
                  nx);
    }

    /* To make life easy we initially set the spacing to 1
     * and correct for this at the end.
     */
    if (bS3)
    {
        /* Fit V''' at the start */
        v3 = v[3] - 3 * v[2] + 3 * v[1] - v[0];
        if (debug)
        {
            fprintf(debug, "The left third derivative is %g\n", v3 / (h * h * h));
        }
        b_s   = 2 * (v[1] - v[0]) + v3 / 6;
        start = 0;
    }
    else
    {
        b_s   = 3 * (v[2] - v[0]) + f[0] * h;
        start = 1;
    }
    if (bE3)
    {
        /* Fit V''' at the end */
        v3 = v[nx - 1] - 3 * v[nx - 2] + 3 * v[nx - 3] - v[nx - 4];
        if (debug)
        {
            fprintf(debug, "The right third derivative is %g\n", v3 / (h * h * h));
        }
        b_e = 2 * (v[nx - 1] - v[nx - 2]) + v3 / 6;
        end = nx;
    }
    else
    {
        /* V'=0 at the end */
        b_e = 3 * (v[nx - 1] - v[nx - 3]) + f[nx - 1] * h;
        end = nx - 1;
    }

    std::vector<double> gamma(nx);
    beta = (bS3 ? 1 : 4);

    /* For V'' fitting */
    /* beta = (bS3 ? 2 : 4); */

    f[start] = b_s / beta;
    for (i = start + 1; i < end; i++)
    {
        gamma[i] = 1 / beta;
        beta     = 4 - gamma[i];
        b        = 3 * (v[i + 1] - v[i - 1]);
        f[i]     = (b - f[i - 1]) / beta;
    }
    gamma[end - 1] = 1 / beta;
    beta           = (bE3 ? 1 : 4) - gamma[end - 1];
    f[end - 1]     = (b_e - f[end - 2]) / beta;

    for (i = end - 2; i >= start; i--)
    {
        f[i] -= gamma[i + 1] * f[i + 1];
    }

    /* Correct for the minus sign and the spacing */
    for (i = start; i < end; i++)
    {
        f[i] = -f[i] / h;
    }
}

static void set_forces(FILE* fp, int angle, int nx, double h, double v[], double f[], int table)
{
    int start, end;

    if (angle == 2)
    {
        gmx_fatal(FARGS, "Force generation for dihedral tables is not (yet) implemented");
    }

    start = 0;
    while (v[start] == 0)
    {
        start++;
    }

    end = nx;
    while (v[end - 1] == 0)
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
        fprintf(fp,
                "Generating forces for table %d, boundary conditions: V''' at %g, %s at %g\n",
                table + 1,
                start * h,
                end == nx ? "V'''" : "V'=0",
                (end - 1) * h);
    }
    spline_forces(end - start, h, v + start, TRUE, end == nx, f + start);
}

static std::vector<t_tabledata> read_tables(FILE* fp, const char* filename, int ntab, int angle)
{
    char     buf[STRLEN];
    double   start, end, dx0, dx1, ssd, vm, vp, f, numf;
    int      k, i, nx0 = 0, nny, ns;
    gmx_bool bAllZero, bZeroV, bZeroF;
    double   tabscale;

    nny                               = 2 * ntab + 1;
    const std::filesystem::path libfn = gmx::findLibraryFile(filename);
    gmx::MultiDimArray<std::vector<double>, gmx::dynamicExtents2D> xvgData    = readXvgData(libfn);
    int                                                            numColumns = xvgData.extent(0);
    if (numColumns != nny)
    {
        gmx_fatal(FARGS,
                  "Trying to read file %s, but nr columns = %d, should be %d",
                  libfn.string().c_str(),
                  numColumns,
                  nny);
    }
    int numRows = xvgData.extent(1);

    const auto& yy = xvgData.asView();
    if (angle == 0)
    {
        if (yy[0][0] != 0.0)
        {
            gmx_fatal(FARGS,
                      "The first distance in file %s is %f nm instead of %f nm",
                      libfn.string().c_str(),
                      yy[0][0],
                      0.0);
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
        if (yy[0][0] != start || yy[0][numRows - 1] != end)
        {
            gmx_fatal(FARGS,
                      "The angles in file %s should go from %f to %f instead of %f to %f\n",
                      libfn.string().c_str(),
                      start,
                      end,
                      yy[0][0],
                      yy[0][numRows - 1]);
        }
    }

    tabscale = (numRows - 1) / (yy[0][numRows - 1] - yy[0][0]);

    if (fp)
    {
        fprintf(fp, "Read user tables from %s with %d data points.\n", libfn.string().c_str(), numRows);
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
        for (i = 0; (i < numRows); i++)
        {
            if (i >= 2)
            {
                dx0 = yy[0][i - 1] - yy[0][i - 2];
                dx1 = yy[0][i] - yy[0][i - 1];
                /* Check for 1% deviation in spacing */
                if (std::fabs(dx1 - dx0) >= 0.005 * (std::fabs(dx0) + std::fabs(dx1)))
                {
                    gmx_fatal(FARGS,
                              "In table file '%s' the x values are not equally spaced: %f %f %f",
                              filename,
                              yy[0][i - 2],
                              yy[0][i - 1],
                              yy[0][i]);
                }
            }
            if (yy[1 + k * 2][i] != 0)
            {
                bZeroV = FALSE;
                if (bAllZero)
                {
                    bAllZero = FALSE;
                    nx0      = i;
                }
                if (yy[1 + k * 2][i] > 0.01 * GMX_REAL_MAX || yy[1 + k * 2][i] < -0.01 * GMX_REAL_MAX)
                {
                    gmx_fatal(FARGS, "Out of range potential value %g in file '%s'", yy[1 + k * 2][i], filename);
                }
            }
            if (yy[1 + k * 2 + 1][i] != 0)
            {
                bZeroF = FALSE;
                if (bAllZero)
                {
                    bAllZero = FALSE;
                    nx0      = i;
                }
                if (yy[1 + k * 2 + 1][i] > 0.01 * GMX_REAL_MAX || yy[1 + k * 2 + 1][i] < -0.01 * GMX_REAL_MAX)
                {
                    gmx_fatal(FARGS, "Out of range force value %g in file '%s'", yy[1 + k * 2 + 1][i], filename);
                }
            }
        }

        if (!bZeroV && bZeroF)
        {
            set_forces(fp, angle, numRows, 1 / tabscale, yy[1 + k * 2].data(), yy[1 + k * 2 + 1].data(), k);
        }
        else
        {
            /* Check if the second column is close to minus the numerical
             * derivative of the first column.
             */
            ssd = 0;
            ns  = 0;
            for (i = 1; (i < numRows - 1); i++)
            {
                vm = yy[1 + 2 * k][i - 1];
                vp = yy[1 + 2 * k][i + 1];
                f  = yy[1 + 2 * k + 1][i];
                if (vm != 0 && vp != 0 && f != 0)
                {
                    /* Take the centered difference */
                    numf = -(vp - vm) * 0.5 * tabscale;
                    if (f + numf != 0)
                    {
                        ssd += std::fabs(2 * (f - numf) / (f + numf));
                    }
                    ns++;
                }
            }
            if (ns > 0)
            {
                ssd /= ns;
                sprintf(buf,
                        "For the %d non-zero entries for table %d in %s the forces deviate on "
                        "average %" PRId64
                        "%% from minus the numerical derivative of the potential\n",
                        ns,
                        k,
                        libfn.string().c_str(),
                        gmx::roundToInt64(100 * ssd));
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
        fprintf(fp, "\nNOTE: All elements in table %s are zero\n\n", libfn.string().c_str());
    }

    std::vector<t_tabledata> td;
    for (k = 0; (k < ntab); k++)
    {
        td.emplace_back(numRows, nx0, tabscale, true);
        for (i = 0; (i < numRows); i++)
        {
            td[k].x[i] = yy[0][i];
            td[k].v[i] = yy[2 * k + 1][i];
            td[k].f[i] = yy[2 * k + 2][i];
        }
    }
    return td;
}

static void fill_table(t_tabledata* td, int tp, const interaction_const_t* ic, gmx_bool b14only)
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
    int    i;
    double reppow, p;
    double r1, rc, r12, r13;
    double r, r2, r6, rc2;
    double expr, Vtab, Ftab;
    /* Parameters for David's function */
    double A = 0, B = 0, C = 0, A_3 = 0, B_4 = 0;
    /* Parameters for the switching function */
    double ksw, swi, swi1;
    /* Temporary parameters */
    gmx_bool bPotentialSwitch, bForceSwitch, bPotentialShift;
    double   ewc   = ic->ewaldcoeff_q;
    double   ewclj = ic->ewaldcoeff_lj;
    double   Vcut  = 0;

    if (b14only)
    {
        bPotentialSwitch = FALSE;
        bForceSwitch     = FALSE;
        bPotentialShift  = FALSE;
    }
    else
    {
        bPotentialSwitch =
                ((tp == etabLJ6Switch) || (tp == etabLJ12Switch) || (tp == etabCOULSwitch)
                 || (tp == etabEwaldSwitch) || (tp == etabEwaldUserSwitch)
                 || (tprops[tp].bCoulomb && (ic->coulomb_modifier == InteractionModifiers::PotSwitch))
                 || (!tprops[tp].bCoulomb && (ic->vdw_modifier == InteractionModifiers::PotSwitch)));
        bForceSwitch =
                ((tp == etabLJ6Shift) || (tp == etabLJ12Shift) || (tp == etabShift)
                 || (tprops[tp].bCoulomb && (ic->coulomb_modifier == InteractionModifiers::ForceSwitch))
                 || (!tprops[tp].bCoulomb && (ic->vdw_modifier == InteractionModifiers::ForceSwitch)));
        bPotentialShift =
                ((tprops[tp].bCoulomb && (ic->coulomb_modifier == InteractionModifiers::PotShift))
                 || (!tprops[tp].bCoulomb && (ic->vdw_modifier == InteractionModifiers::PotShift)));
    }

    reppow = ic->reppow;

    if (tprops[tp].bCoulomb)
    {
        r1 = ic->rcoulomb_switch;
        rc = ic->rcoulomb;
    }
    else
    {
        r1 = ic->rvdw_switch;
        rc = ic->rvdw;
    }
    if (bPotentialSwitch)
    {
        ksw = 1.0 / (gmx::power5(rc - r1));
    }
    else
    {
        ksw = 0.0;
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

        A = p * ((p + 1) * r1 - (p + 4) * rc) / (std::pow(rc, p + 2) * gmx::square(rc - r1));
        B = -p * ((p + 1) * r1 - (p + 3) * rc) / (std::pow(rc, p + 2) * gmx::power3(rc - r1));
        C = 1.0 / std::pow(rc, p) - A / 3.0 * gmx::power3(rc - r1) - B / 4.0 * gmx::power4(rc - r1);
        if (tp == etabLJ6Shift)
        {
            A = -A;
            B = -B;
            C = -C;
        }
        A_3 = A / 3.0;
        B_4 = B / 4.0;
    }
    if (debug)
    {
        fprintf(debug, "Setting up tables\n");
        fflush(debug);
    }

    if (bPotentialShift)
    {
        rc2        = rc * rc;
        double rc6 = 1.0 / (rc2 * rc2 * rc2);
        double rc12;
        if (gmx_within_tol(reppow, 12.0, 10 * GMX_DOUBLE_EPS))
        {
            rc12 = rc6 * rc6;
        }
        else
        {
            rc12 = std::pow(rc, -reppow);
        }

        switch (tp)
        {
            case etabLJ6:
                /* Dispersion */
                Vcut = -rc6;
                break;
            case etabLJ6Ewald:
                Vcut = -rc6 * std::exp(-ewclj * ewclj * rc2)
                       * (1 + ewclj * ewclj * rc2 + gmx::power4(ewclj) * rc2 * rc2 / 2);
                break;
            case etabLJ12:
                /* Repulsion */
                Vcut = rc12;
                break;
            case etabCOUL: Vcut = 1.0 / rc; break;
            case etabEwald:
            case etabEwaldSwitch: Vcut = std::erfc(ewc * rc) / rc; break;
            case etabEwaldUser:
                /* Only calculate minus the reciprocal space contribution */
                Vcut = -std::erf(ewc * rc) / rc;
                break;
            case etabRF:
            case etabRF_ZERO:
                /* No need for preventing the usage of modifiers with RF */
                Vcut = 0.0;
                break;
            case etabEXPMIN: Vcut = std::exp(-rc); break;
            default:
                gmx_fatal(FARGS,
                          "Cannot apply new potential-shift modifier to interaction type '%s' yet. "
                          "(%s,%d)",
                          tprops[tp].name,
                          __FILE__,
                          __LINE__);
        }
    }

    for (i = 0; (i < td->nx); i++)
    {
        td->x[i] = i / td->tabscale;
    }
    for (i = td->nx0; (i < td->nx); i++)
    {
        r  = td->x[i];
        r2 = r * r;
        r6 = 1.0 / (r2 * r2 * r2);
        if (gmx_within_tol(reppow, 12.0, 10 * GMX_DOUBLE_EPS))
        {
            r12 = r6 * r6;
        }
        else
        {
            r12 = std::pow(r, -reppow);
        }
        Vtab = 0.0;
        Ftab = 0.0;
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
                swi = 1 - 10 * gmx::power3(r - r1) * ksw * gmx::square(rc - r1)
                      + 15 * gmx::power4(r - r1) * ksw * (rc - r1) - 6 * gmx::power5(r - r1) * ksw;
                swi1 = -30 * gmx::square(r - r1) * ksw * gmx::square(rc - r1)
                       + 60 * gmx::power3(r - r1) * ksw * (rc - r1) - 30 * gmx::power4(r - r1) * ksw;
            }
        }
        else /* not really needed, but avoids compiler warnings... */
        {
            swi  = 1.0;
            swi1 = 0.0;
        }

        switch (tp)
        {
            case etabLJ6:
                /* Dispersion */
                Vtab = -r6;
                Ftab = 6.0 * Vtab / r;
                break;
            case etabLJ6Switch:
            case etabLJ6Shift:
                /* Dispersion */
                if (r < rc)
                {
                    Vtab = -r6;
                    Ftab = 6.0 * Vtab / r;
                    break;
                }
                break;
            case etabLJ12:
                /* Repulsion */
                Vtab = r12;
                Ftab = reppow * Vtab / r;
                break;
            case etabLJ12Switch:
            case etabLJ12Shift:
                /* Repulsion */
                if (r < rc)
                {
                    Vtab = r12;
                    Ftab = reppow * Vtab / r;
                }
                break;
            case etabCOUL:
                Vtab = 1.0 / r;
                Ftab = 1.0 / r2;
                break;
            case etabCOULSwitch:
            case etabShift:
                if (r < rc)
                {
                    Vtab = 1.0 / r;
                    Ftab = 1.0 / r2;
                }
                break;
            case etabEwald:
            case etabEwaldSwitch:
                Vtab = std::erfc(ewc * r) / r;
                Ftab = std::erfc(ewc * r) / r2 + std::exp(-(ewc * ewc * r2)) * ewc * M_2_SQRTPI / r;
                break;
            case etabEwaldUser:
            case etabEwaldUserSwitch:
                /* Only calculate the negative of the reciprocal space contribution */
                Vtab = -std::erf(ewc * r) / r;
                Ftab = -std::erf(ewc * r) / r2 + std::exp(-(ewc * ewc * r2)) * ewc * M_2_SQRTPI / r;
                break;
            case etabLJ6Ewald:
                Vtab = -r6 * std::exp(-ewclj * ewclj * r2)
                       * (1 + ewclj * ewclj * r2 + gmx::power4(ewclj) * r2 * r2 / 2);
                Ftab = 6.0 * Vtab / r
                       - r6 * std::exp(-ewclj * ewclj * r2) * gmx::power5(ewclj) * ewclj * r2 * r2 * r;
                break;
            case etabRF:
            case etabRF_ZERO:
                Vtab = 1.0 / r + ic->reactionFieldCoefficient * r2 - ic->reactionFieldShift;
                Ftab = 1.0 / r2 - 2 * ic->reactionFieldCoefficient * r;
                if (tp == etabRF_ZERO && r >= rc)
                {
                    Vtab = 0;
                    Ftab = 0;
                }
                break;
            case etabEXPMIN:
                expr = std::exp(-r);
                Vtab = expr;
                Ftab = expr;
                break;
            default:
                gmx_fatal(FARGS, "Table type %d not implemented yet. (%s,%d)", tp, __FILE__, __LINE__);
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
                    r12 = (r - r1) * (r - r1);
                    r13 = (r - r1) * r12;
                    Vtab += -A_3 * r13 - B_4 * r12 * r12;
                    Ftab += A * r12 + B * r13;
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
                Ftab = Ftab * swi - Vtab * swi1;
                Vtab = Vtab * swi;
            }
        }
        /* Convert to single precision when we store to mem */
        td->v[i] = Vtab;
        td->f[i] = Ftab;
    }

    /* Continue the table linearly from nx0 to 0.
     * These values are only required for energy minimization with overlap or TPI.
     */
    for (i = td->nx0 - 1; i >= 0; i--)
    {
        td->v[i] = td->v[i + 1] + td->f[i + 1] * (td->x[i + 1] - td->x[i]);
        td->f[i] = td->f[i + 1];
    }
}

static void set_table_type(int tabsel[], const interaction_const_t* ic, gmx_bool b14only)
{
    /* Set the different table indices.
     * Coulomb first.
     */

    CoulombInteractionType eltype;
    VanDerWaalsType        vdwtype;

    if (b14only)
    {
        switch (ic->eeltype)
        {
            case CoulombInteractionType::User:
            case CoulombInteractionType::PmeUser:
            case CoulombInteractionType::PmeUserSwitch:
                eltype = CoulombInteractionType::User;
                break;
            default: eltype = CoulombInteractionType::Cut;
        }
    }
    else
    {
        eltype = ic->eeltype;
    }

    switch (eltype)
    {
        case CoulombInteractionType::Cut: tabsel[etiCOUL] = etabCOUL; break;
        case CoulombInteractionType::Poisson: tabsel[etiCOUL] = etabShift; break;
        case CoulombInteractionType::Shift:
            if (ic->rcoulomb > ic->rcoulomb_switch)
            {
                tabsel[etiCOUL] = etabShift;
            }
            else
            {
                tabsel[etiCOUL] = etabCOUL;
            }
            break;
        case CoulombInteractionType::Ewald:
        case CoulombInteractionType::Pme:
        case CoulombInteractionType::P3mAD: tabsel[etiCOUL] = etabEwald; break;
        case CoulombInteractionType::PmeSwitch: tabsel[etiCOUL] = etabEwaldSwitch; break;
        case CoulombInteractionType::PmeUser: tabsel[etiCOUL] = etabEwaldUser; break;
        case CoulombInteractionType::PmeUserSwitch: tabsel[etiCOUL] = etabEwaldUserSwitch; break;
        case CoulombInteractionType::RF:
        case CoulombInteractionType::RFZero: tabsel[etiCOUL] = etabRF_ZERO; break;
        case CoulombInteractionType::Switch: tabsel[etiCOUL] = etabCOULSwitch; break;
        case CoulombInteractionType::User: tabsel[etiCOUL] = etabUSER; break;
        default: gmx_fatal(FARGS, "Invalid eeltype %s", enumValueToString(eltype));
    }

    /* Van der Waals time */
    if (ic->useBuckingham && !b14only)
    {
        tabsel[etiLJ6]  = etabLJ6;
        tabsel[etiLJ12] = etabEXPMIN;
    }
    else
    {
        if (b14only && ic->vdwtype != VanDerWaalsType::User)
        {
            vdwtype = VanDerWaalsType::Cut;
        }
        else
        {
            vdwtype = ic->vdwtype;
        }

        switch (vdwtype)
        {
            case VanDerWaalsType::Switch:
                tabsel[etiLJ6]  = etabLJ6Switch;
                tabsel[etiLJ12] = etabLJ12Switch;
                break;
            case VanDerWaalsType::Shift:
                tabsel[etiLJ6]  = etabLJ6Shift;
                tabsel[etiLJ12] = etabLJ12Shift;
                break;
            case VanDerWaalsType::User:
                tabsel[etiLJ6]  = etabUSER;
                tabsel[etiLJ12] = etabUSER;
                break;
            case VanDerWaalsType::Cut:
                tabsel[etiLJ6]  = etabLJ6;
                tabsel[etiLJ12] = etabLJ12;
                break;
            case VanDerWaalsType::Pme:
                tabsel[etiLJ6]  = etabLJ6Ewald;
                tabsel[etiLJ12] = etabLJ12;
                break;
            default:
                gmx_fatal(FARGS, "Invalid vdwtype %s in %s line %d", enumValueToString(vdwtype), __FILE__, __LINE__);
        }

        if (!b14only && ic->vdw_modifier != InteractionModifiers::None)
        {
            if (ic->vdw_modifier != InteractionModifiers::PotShift && ic->vdwtype != VanDerWaalsType::Cut)
            {
                gmx_incons(
                        "Potential modifiers other than potential-shift are only implemented for "
                        "LJ cut-off");
            }

            /* LJ-PME and other (shift-only) modifiers are handled by applying the modifiers
             * to the original interaction forms when we fill the table, so we only check cutoffs here.
             */
            if (ic->vdwtype == VanDerWaalsType::Cut)
            {
                switch (ic->vdw_modifier)
                {
                    case InteractionModifiers::None:
                    case InteractionModifiers::PotShift:
                    case InteractionModifiers::ExactCutoff:
                        /* No modification */
                        break;
                    case InteractionModifiers::PotSwitch:
                        tabsel[etiLJ6]  = etabLJ6Switch;
                        tabsel[etiLJ12] = etabLJ12Switch;
                        break;
                    case InteractionModifiers::ForceSwitch:
                        tabsel[etiLJ6]  = etabLJ6Shift;
                        tabsel[etiLJ12] = etabLJ12Shift;
                        break;
                    default: gmx_incons("Unsupported vdw_modifier");
                }
            }
        }
    }
}

std::unique_ptr<t_forcetable>
make_tables(FILE* fp, const interaction_const_t* ic, const char* fn, real rtab, int flags)
{
    gmx_bool b14only, useUserTable;
    int      nx0, tabsel[etiNR];
    real     scalefactor;

    auto table = std::make_unique<t_forcetable>(
            TableInteraction::ElectrostaticVdwRepulsionVdwDispersion, TableFormat::CubicsplineYfgh);

    b14only = ((flags & GMX_MAKETABLES_14ONLY) != 0);

    if (flags & GMX_MAKETABLES_FORCEUSER)
    {
        tabsel[etiCOUL] = etabUSER;
        tabsel[etiLJ6]  = etabUSER;
        tabsel[etiLJ12] = etabUSER;
    }
    else
    {
        set_table_type(tabsel, ic, b14only);
    }
    std::vector<t_tabledata> td;
    table->interactionRange = rtab;
    table->scale            = 0;
    table->numTablePoints   = 0;

    table->numInteractions = etiNR;
    table->stride          = table->formatsize * table->numInteractions;

    /* Check whether we have to read or generate */
    useUserTable = FALSE;
    for (unsigned int i = 0; (i < etiNR); i++)
    {
        if (ETAB_USER(tabsel[i]))
        {
            useUserTable = TRUE;
        }
    }
    if (useUserTable)
    {
        td = read_tables(fp, fn, etiNR, 0);
        if (rtab == 0 || (flags & GMX_MAKETABLES_14ONLY))
        {
            table->numTablePoints = td[0].nx;
        }
        else
        {
            if (td[0].x[td[0].nx - 1] < rtab)
            {
                gmx_fatal(FARGS,
                          "Tables in file %s not long enough for cut-off:\n"
                          "\tshould be at least %f nm\n",
                          fn,
                          rtab);
            }
            table->numTablePoints = gmx::roundToInt(rtab * td[0].tabscale);
        }
        table->scale = td[0].tabscale;
        nx0          = td[0].nx0;
    }
    else
    {
        td.resize(etiNR);
        // No tables are read
#if GMX_DOUBLE
        table->scale = 2000.0;
#else
        table->scale = 500.0;
#endif
        table->numTablePoints = static_cast<int>(rtab * table->scale);
        nx0                   = 10;
    }

    /* Each table type (e.g. coul,lj6,lj12) requires four
     * numbers per table->numTablePoints+1 data points. For performance reasons we want
     * the table data to be aligned to (at least) a 32-byte boundary.
     */
    table->data.resize(table->stride * (table->numTablePoints + 1) * sizeof(real));

    for (int k = 0; (k < etiNR); k++)
    {
        /* Now fill data for tables that have not been read
         * or add the Ewald long-range correction for Ewald user tables.
         */
        if (tabsel[k] != etabUSER)
        {
            real scale = table->scale;
            if (ic->useBuckingham && (ic->buckinghamBMax != 0) && tabsel[k] == etabEXPMIN)
            {
                scale /= ic->buckinghamBMax;
            }
            td[k] = t_tabledata(table->numTablePoints, nx0, scale, !useUserTable);

            fill_table(&(td[k]), tabsel[k], ic, b14only);
            if (fp)
            {
                fprintf(fp,
                        "Generated table with %d data points for %s%s.\n"
                        "Tabscale = %g points/nm\n",
                        td[k].nx,
                        b14only ? "1-4 " : "",
                        tprops[tabsel[k]].name,
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
            scalefactor = 1.0 / 6.0;
        }
        else if (k == etiLJ12 && tabsel[k] != etabEXPMIN)
        {
            scalefactor = 1.0 / 12.0;
        }
        else
        {
            scalefactor = 1.0;
        }

        copy2table(table->numTablePoints,
                   k * table->formatsize,
                   table->stride,
                   td[k].x,
                   td[k].v,
                   td[k].f,
                   scalefactor,
                   table->data);
    }

    return table;
}

bondedtable_t make_bonded_table(FILE* fplog, const char* fn, int angle)
{
    int           i;
    bondedtable_t tab;
    int           stride = 4;

    t_tabledata td = read_tables(fplog, fn, 1, angle)[0];
    if (angle > 0)
    {
        /* Convert the table from degrees to radians */
        for (i = 0; i < td.nx; i++)
        {
            td.x[i] *= gmx::c_deg2Rad;
            td.f[i] *= gmx::c_rad2Deg;
        }
        td.tabscale *= gmx::c_rad2Deg;
    }
    tab.n     = td.nx;
    tab.scale = td.tabscale;
    tab.data.resize(tab.n * stride);
    copy2table(tab.n, 0, stride, td.x, td.v, td.f, 1.0, tab.data);

    return tab;
}

std::unique_ptr<t_forcetable>
makeDispersionCorrectionTable(FILE* fp, const interaction_const_t* ic, real rtab, const char* tabfn)
{
    GMX_RELEASE_ASSERT(ic->vdwtype != VanDerWaalsType::User || tabfn,
                       "With VdW user tables we need a table file name");

    std::unique_ptr<t_forcetable> fullTable = make_tables(fp, ic, tabfn, rtab, 0);
    /* Copy the contents of the table to one that has just dispersion
     * and repulsion, to improve cache performance. We want the table
     * data to be aligned to 32-byte boundaries.
     */
    std::unique_ptr<t_forcetable> dispersionCorrectionTable = std::make_unique<t_forcetable>(
            TableInteraction::VdwRepulsionVdwDispersion, fullTable->format_);
    dispersionCorrectionTable->interactionRange = fullTable->interactionRange;
    dispersionCorrectionTable->numTablePoints   = fullTable->numTablePoints;
    dispersionCorrectionTable->scale            = fullTable->scale;
    dispersionCorrectionTable->numInteractions  = 2;
    dispersionCorrectionTable->stride =
            dispersionCorrectionTable->formatsize * dispersionCorrectionTable->numInteractions;
    dispersionCorrectionTable->data.resize(dispersionCorrectionTable->stride
                                           * (dispersionCorrectionTable->numTablePoints + 1));

    for (int i = 0; i <= fullTable->numTablePoints; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            dispersionCorrectionTable->data[8 * i + j] = fullTable->data[12 * i + 4 + j];
        }
    }

    return dispersionCorrectionTable;
}

t_forcetable::t_forcetable(enum TableInteraction interaction, enum TableFormat format) :
    interaction_(interaction),
    format_(format),
    interactionRange(0),
    numTablePoints(0),
    scale(0),
    numInteractions(0),
    stride(0)
{
}

t_forcetable::~t_forcetable() = default;
