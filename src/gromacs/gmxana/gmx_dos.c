/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/correlationfunctions/manyautocorrelation.h"
#include "gromacs/utility/cstringutil.h"
#include "gmx_ana.h"

/* TYPE OF VELOCITIES: trans, angular, vibrational, rotational, total */
enum Velocity {
    V_T, V_A, V_V, V_R, V_ATOM, V_SUM, V_NR
};

static const char *velocityName[V_NR] = {
    "Trans", "Angul", "Vibr", "Rotat", "Atomic", "Sum"
};

static const char *velocityLong[V_NR] = {
    "Translation", "Angular motion", "Vibration", "Rotational motion",
    "Atomic motion", "Sum of it all"
};

enum Component {
    C_SOLID, C_GAS, C_TOTAL, C_NR
};

static const char *cName[C_NR] = {
    "Solid", "Gas", "Total"
};

enum ThermoMethod {
    TM_1PT, TM_2PT, TM_MD, TM_NR
};

static const char *tmName[TM_NR] = {
    "1PT", "2PT", "MD"
};

enum LevelOfTheory {
    LOT_CLASSICAL, LOT_QUANTUM, LOT_QUANTUM_CORRECTION, LOT_NR
};

static const char *lotName[LOT_NR] = {
    "Classical", "Quantum", "Quant Corr"
};

enum Density {
    DOS_CV, DOS_S, DOS_A, DOS_E, DOS_NR
};

static const char *dosName[DOS_NR] = {
    "Heat capacity",
    "Entropy",
    "Helmholtz energy",
    "Internal energy"
};

static const char *dosShort[DOS_NR] = {
    "cV",
    "S",
    "A",
    "E"
};

static const char *dosUnit[DOS_NR] = {
    "(J/mol K)",
    "(J/mol K)",
    "(kJ/mol)",
    "(kJ/mol)"
};

typedef struct {
    double  X[DOS_NR][C_NR];
    real   *dos[C_NR];
    real   *wdos[DOS_NR][C_NR];
} t_dos_data;

typedef struct {
    int         npart;
    t_dos_data  dd[LOT_NR][TM_NR];
    real       *alldos;
    double      wDiff[DOS_NR];
    double      Emd, E0;
    double      f, y, fy, ry, z, nrho, Delta, beta, bfac, partMass;
    double      hsdf, ttdf;
    double      sigHS, Shs, Sig;
    double      DiffACF, DiffDoS;
    double      dostot, dof, DoS0;
    double      Tdiff, T;
    /* Velocity components. Dimension DIM*npart x nframes */
    real      **vec;
} t_dosprops;

static void init_dos_data(t_dosprops *dp, int N_tot)
{
    int i, j, k, l;

    snew(dp->alldos, N_tot);
    for (k = 0; (k < LOT_NR); k++)
    {
        for (l = 0; (l < TM_NR); l++)
        {
            for (i = 0; (i < C_NR); i++)
            {
                snew_aligned(dp->dd[k][l].dos[i], N_tot, 16);
                for (j = 0; (j < DOS_NR); j++)
                {
                    dp->dd[k][l].X[j][i] = 0;
                    snew_aligned(dp->dd[k][l].wdos[j][i], N_tot, 16);
                }
            }
        }
    }
}

/* double precision dsub_xcm */
static void dsub_xcm(rvec x[], int gnx, atom_id *index, t_atom atom[], rvec xcm)
{
    int    i, m;
    double m0, tm, dxcm[DIM], xx[DIM], a, b, c;

    dxcm[XX] = dxcm[YY] = dxcm[ZZ] = 0.0;

    tm = 0.0;
    for (i = 0; (i < gnx); i++)
    {
        int ii     = index ? index[i] : i;
        xx[XX] = x[ii][XX];
        xx[YY] = x[ii][YY];
        xx[ZZ] = x[ii][ZZ];
        m0     = atom[ii].m;
        tm    += m0;
        for (m = 0; (m < DIM); m++)
        {
            dxcm[m] += m0*xx[m];
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        dxcm[m] /= tm;
    }


    for (i = 0; (i < gnx); i++)
    {
        int ii        = index ? index[i] : i;
        a         = (double)x[ii][XX] - dxcm[XX];
        b         = (double)x[ii][YY] - dxcm[YY];
        c         = (double)x[ii][ZZ] - dxcm[ZZ];
        x[ii][XX] = a;
        x[ii][YY] = b;
        x[ii][ZZ] = c;
    }

    xcm[XX] = dxcm[XX];
    xcm[YY] = dxcm[YY];
    xcm[ZZ] = dxcm[ZZ];
}

/* ********************************************************** */

static double FD(double Delta, double f) /* eq 34 JCP 119, 11792 (2003) */
{
    if (0)
    {
        return (2*pow(Delta, -4.5)*pow(f, 7.5) -
                6*pow(Delta, -3)*pow(f, 5) -
                pow(Delta, -1.5)*pow(f, 3.5) +
                6*pow(Delta, -1.5)*pow(f, 2.5) +
                2*f - 2);
        // ~14+80 = 94 flops
    }
    else
    {
        double D_15f25 = f*f*sqrt(f/(Delta*Delta*Delta));
        return (2.0*D_15f25*D_15f25*(D_15f25 - 3.0) -
                D_15f25*f +
                6.0*D_15f25 +
                2.0*f - 2.0);
        // ~29 flops
    }
}

/* fluidicity eq. for hard spheres packing, eq 31 JCP 119, 11792 (2003) */
static double YYY(double f, double y)
{
    //return (2*pow(y*f, 3) - sqr(f)*y*(1+6*y) + (2+6*y)*f - 2);
    double yf = y*f;
    return 2.0*yf*yf*yf - f*yf*(1.0+6.0*y) + (2.0+6.0*y)*f - 2.0;
}

/* compressibility of hard spheres z=(1+y+y^3-y^3)/(1-y)^3 */
static double calc_compress(double y)
{
    //    return ((1+y+sqr(y)-pow(y, 3))/(pow((1-y), 3)));
    return (1 + y + y*y - y*y*y)/((1-y) * (1-y) * (1-y));
}

static double bisector(double Delta, double tol,
                       double ff0, double ff1,
                       double ff(double, double))
{
    double fd0, fd, fd1, f, f0, f1;
    double tolmin = 1e-8;

    f0 = ff0;
    f1 = ff1;
    if (tol < tolmin)
    {
        fprintf(stderr, "Unrealistic tolerance %g for bisector. Setting it to %g\n",
                tol, tolmin);
        tol = tolmin;
    }

    do
    {
        fd0 = ff(Delta, f0);
        fd1 = ff(Delta, f1);
        f   = (f0+f1)*0.5;
        fd  = ff(Delta, f);
        if (fd < 0)
        {
            f0 = f;
        }
        else if (fd > 0)
        {
            f1 = f;
        }
        else
        {
            return f;
        }
    }
    while ((f1-f0) > tol);

    return f;
}

/* calculate fluidicity f */
static double calc_fluidicity(double Delta, double tol)
{
    /* solution for FD */
    return bisector(Delta, tol, 0, 1, FD);
}

/* hard sphere packing fraction, eq 32 JCP 119, 11792 (2003) */
static double calc_y(double f, double Delta, double toler)
{
    double y1, y2;

    y1 = pow(f/Delta, 1.5);
    y2 = bisector(f, toler, 0, 10000, YYY);
    if (fabs((y1-y2)/(y1+y2)) > 100*toler)
    {
        fprintf(stderr, "Inconsistency computing y: y1 = %f, y2 = %f, using y1.\n",
                y1, y2);
    }

    return y1;
}

static double calc_Shs(double fy)
/* entropy for hard spheres */
{
    double z = calc_compress(fy);

    if (0 == z)
    {
        fprintf(stderr, "Zero compressibility found in calc_Shs\n");
        return 0;
    }
    else if (1 == fy)
    {
        fprintf(stderr, "Hard sphere packing fraction is 1 in calc_Shs\n");
        return 0;
    }
    else
    {
        return BOLTZ*(log(z) + fy*(3*fy-4)/sqr(1-fy));
    }
}

static real weightE(int LevelOfTheory, int ThermoMethod,
                    int Component, real nu, real beta)
{
    real u = beta*PLANCK*nu;

    switch (LevelOfTheory)
    {
        case LOT_CLASSICAL:
            /* Berens, eq. 3.33 */
            return 1;
        case LOT_QUANTUM:
            if (u == 0)
            {
                return 1.0;
            }
            else
            {
                /* Berens, eq. 3.40 */
                return u/2 + u/(exp(u)-1);
            }
        case LOT_QUANTUM_CORRECTION:
            return (weightE(LOT_QUANTUM, ThermoMethod, Component, nu, beta) -
                    weightE(LOT_CLASSICAL, ThermoMethod, Component, nu, beta));
        default:
            gmx_fatal(FARGS, "Incorrect LevelOfTheory %d", LevelOfTheory);
    }
    return 0;
}

static real weightS(int LevelOfTheory, int ThermoMethod, int Component,
                    real nu, real beta)
{
    real u = beta*PLANCK*nu;

    if (u == 0)
    {
        /* in the limit u -> 0, eq 3.43 -> + Infinity.
         * set to 0 because S(v)  =0 for Solid
         */
        return 0;
    }
    switch (LevelOfTheory)
    {
        case LOT_CLASSICAL:
            /* Berens, eq. 3.35 */
            return (1.0 - log(u));
        case LOT_QUANTUM:
            /* Berens, eq. 3.43 */
            return u/(exp(u)-1) - log(1-exp(-u));
        case LOT_QUANTUM_CORRECTION:
            return (weightS(LOT_QUANTUM, ThermoMethod, Component, nu, beta) -
                    weightS(LOT_CLASSICAL, ThermoMethod, Component, nu, beta));
        default:
            gmx_fatal(FARGS, "Incorrect LevelOfTheory %d", LevelOfTheory);
    }
    return 0;
}

static real weightA(int LevelOfTheory, int ThermoMethod, int Component,
                    real nu, real beta)
{
    real u = beta*PLANCK*nu;

    if (u == 0)
    {
        /* In the limit u -> 0, both Berens eq 3.35 and eq 3.42 -> - Infinity
         * Return 0 because S(v) = 0 for solids anyway.
         */
        return 0;
    }

    switch (LevelOfTheory)
    {
        case LOT_CLASSICAL:
            /* Berens, eq 3.35 */
            return log(u);
        case LOT_QUANTUM:
            /* Berense eq 3.42 */
            return log((1-exp(-u))/(exp(-u/2)));
        case LOT_QUANTUM_CORRECTION:
            return (weightA(LOT_QUANTUM, ThermoMethod, Component, nu, beta) -
                    weightA(LOT_CLASSICAL, ThermoMethod, Component, nu, beta));
        default:
            gmx_fatal(FARGS, "Incorrect LevelOfTheory %d", LevelOfTheory);
    }
    return 0;
}

static real weightCv(int LevelOfTheory, int ThermoMethod, int Component,
                     real nu, real beta)
{
    switch (LevelOfTheory)
    {
        case LOT_CLASSICAL:
            /* Berens, eq 3.34 */
            return 1.0;
        case LOT_QUANTUM:
        {
            real u = beta*PLANCK*nu;
            real ebn, koko;

            if (u == 0)
            {
                /* in the limit u -> 0, eq 3.41 -> 1 */
                return 1.0;

            }
            if (u > 44.3437) /* empirical limit */
            {
                return 0.0;
            }
            else
            {
                /* Berens, eq. 3.41 */
                ebn  = exp(u);
                koko = sqr(1-ebn);
                return sqr(u)*ebn/koko;
            }
        }
        case LOT_QUANTUM_CORRECTION:
            return (weightCv(LOT_QUANTUM, ThermoMethod, Component, nu, beta) -
                    weightCv(LOT_CLASSICAL, ThermoMethod, Component, nu, beta));
        default:
            gmx_fatal(FARGS, "Incorrect LevelOfTheory %d", LevelOfTheory);
    }
    return 0;
}

static real weight(int Density, int LevelOfTheory,
                   int ThermoMethod, int Component,
                   real nu, real beta)
{
    real (*w[DOS_NR]) (int, int, int, real, real) = {
        weightCv, weightS, weightA, weightE
    };
    range_check(Density, 0, DOS_NR);
    range_check(LevelOfTheory, 0, LOT_NR);
    range_check(ThermoMethod, 0, TM_NR);
    range_check(Component, 0, C_NR);

    return w[Density](LevelOfTheory, ThermoMethod, Component, nu, beta);
}

static void dump_fy(output_env_t oenv, real toler)
{
    FILE       *fp;
    double      Delta, f, y, DD;
    const char *leg[] = { "f", "fy", "y" };

    DD = pow(10.0, 0.125);
    fp = xvgropen("fy.xvg", "Fig. 2, Lin2003a", "Delta", "y or fy", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@    world 1e-05, 0, 1000, 1\n");
        fprintf(fp, "@    xaxes scale Logarithmic\n");
    }
    for (Delta = 1e-5; (Delta <= 1000); Delta *= DD)
    {
        f = calc_fluidicity(Delta, toler);
        y = calc_y(f, Delta, toler);
        fprintf(fp, "%10g  %10g  %10g  %10g\n", Delta, f, f*y, y);
    }
    xvgrclose(fp);
}

static void dump_w0(output_env_t oenv, real beta, char *title, int lot)
{
    FILE       *fp;
    double      nu;
    const char *leg[]    = { "wCv", "wS", "wA", "wE" };
    char        buf[128];

    range_check(lot, 0, LOT_NR);
    snprintf(buf, 128, "w%s.xvg", lotName[lot]);
    fp = xvgropen(buf, title, "\\f{12}b\\f{4}h\\f{12}n",
                  "w", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);
    for (nu = 0; (nu < 100); nu += 0.05)
    {
        fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n", beta*PLANCK*nu,
                weight(DOS_CV, lot, TM_2PT, C_TOTAL, nu, beta),
                weight(DOS_S,  lot, TM_2PT, C_TOTAL, nu, beta),
                weight(DOS_A,  lot, TM_2PT, C_TOTAL, nu, beta),
                weight(DOS_E,  lot, TM_2PT, C_TOTAL, nu, beta));
    }
    fclose(fp);
}

static void dump_w(output_env_t oenv, real beta)
{
    dump_w0(oenv, beta, "Fig. 1, Berens1983a", LOT_QUANTUM);
    dump_w0(oenv, beta, "Classical weighting", LOT_CLASSICAL);
}

#define NDIM 4
void SWAPPER(double dd[DIM], double **ev, int i)
{
    double temp;
    double tvec[NDIM];
    int    j;

    if (fabs(dd[i+1]) > fabs(dd[i]))
    {
        temp = dd[i];
        for (j = 0; (j < NDIM); j++)
        {
            tvec[j] = ev[j][i];
        }
        dd[i] = dd[i+1];
        for (j = 0; (j < NDIM); j++)
        {
            ev[j][i] = ev[j][i+1];
        }
        dd[i+1] = temp;
        for (j = 0; (j < NDIM); j++)
        {
            ev[j][i+1] = tvec[j];
        }
    }
}

static void principal(int n, atom_id index[], t_atom atom[], rvec x[],
                      matrix trans, rvec d, double **inten, double **ev)
{
#define NDIM 4
    int      i, j, ai, m, nrot;
    real     mm, rx, ry, rz;
    double   dd[NDIM];

    /*for (i = 0; (i < NDIM); i++)
       {
        dd[i] = 0.0;
        for (j = 0; (j < NDIM); j++)
        {
            ev[i][j] = 0;
        }
        }*/

    for (i = 0; (i < NDIM); i++)
    {
        for (m = 0; (m < NDIM); m++)
        {
            inten[i][m] = 0;
        }
    }
    for (i = 0; (i < n); i++)
    {
        ai           = index[i];
        mm           = atom[ai].m;
        rx           = x[ai][XX];
        ry           = x[ai][YY];
        rz           = x[ai][ZZ];
        inten[0][0] += mm*(sqr(ry)+sqr(rz));
        inten[1][1] += mm*(sqr(rx)+sqr(rz));
        inten[2][2] += mm*(sqr(rx)+sqr(ry));
        inten[1][0] -= mm*(ry*rx);
        inten[2][0] -= mm*(rx*rz);
        inten[2][1] -= mm*(rz*ry);
    }
    inten[0][1] = inten[1][0];
    inten[0][2] = inten[2][0];
    inten[1][2] = inten[2][1];
#ifdef DEBUG
    ptrans("initial", inten, dd, e);
#endif

    /* Call numerical recipe routines */
    jacobi(inten, 3, dd, ev, &nrot);

    /* Sort eigenvalues in ascending order */
    SWAPPER(dd, ev, 0);
    SWAPPER(dd, ev, 1);
    SWAPPER(dd, ev, 0);

    for (i = 0; (i < DIM); i++)
    {
        d[i] = dd[i];
        for (m = 0; (m < DIM); m++)
        {
            trans[i][m] = ev[i][m];
        }
    }
#undef NDIM
}

static void print_dp(FILE *fplog, int t, t_dosprops *dp)
{
    int i, j;

    if (dp->npart == 0)
    {
        fprintf(fplog, "\n No particles in %s\n", velocityLong[t]);
        return;
    }
    fprintf(fplog, "\n+++ Analyzing %s +++\n", velocityLong[t]);
    fprintf(fplog, "Npart                           = %d\n", dp->npart);
    fprintf(fplog, "T                               = %g\n", dp->T);
    fprintf(fplog, "Number density                  = %g particles/nm^3\n", dp->nrho);
    fprintf(fplog, "Delta                           = %g\n", dp->Delta);
    fprintf(fplog, "fluidicity f                    = %g\n", dp->f);
    fprintf(fplog, "ry                              = %g\n", dp->ry);
    fprintf(fplog, "hsdf                            = %g\n", dp->hsdf);
    fprintf(fplog, "ttdf                            = %g\n", dp->ttdf);
    fprintf(fplog, "fy                              = %g nm\n", dp->fy);
    fprintf(fplog, "hard sphere packing fraction y  = %g\n", dp->y);
    fprintf(fplog, "hard sphere compressibility z   = %g\n", dp->z);
    fprintf(fplog, "Particle mass                   = %g amu\n", dp->partMass);
    fprintf(fplog, "ideal gas entropy Sig           = %g 1/K\n", dp->Sig / BOLTZ);
    fprintf(fplog, "hard sphere entropy Shs         = %g 1/K\n", dp->Shs / BOLTZ);
    fprintf(fplog, "sigma_HS                        = %g nm\n", dp->sigHS);
    fprintf(fplog, "DoS0                            = %g\n", dp->DoS0);
    fprintf(fplog, "DoSTot                          = %g\n", dp->dostot);
    if (V_V != t)
    {
        fprintf(fplog, "Diffusion coefficient from dos[VACF] %g 10^-5 cm^2/s\n",
                dp->DiffACF);
        fprintf(fplog, "Diffusion coefficient from DoS0 %g 10^-5 cm^2/s\n",
                dp->DiffDoS);
    }
}

static void print_legend(FILE *fp, t_dosprops dp[], output_env_t oenv)
{
    const char *leg[V_NR];
    int         t, nleg = 0;

    for (t = 0; (t < V_SUM); t++)
    {
        if (dp[t].npart > 0)
        {
            leg[nleg++] = velocityName[t];
        }
    }
    xvgr_legend(fp, nleg, leg, oenv);
}

static void print_stuff(int n, real tt[], t_dosprops dp[], real dos_vacf[],
                        const char *vacf, const char *mvacf, output_env_t oenv)
{
    FILE *fp;
    int   j, t;

    if ((NULL != vacf) && (NULL != dos_vacf))
    {
        fp = xvgropen(vacf, "Velocity ACF",
                      "Time (ps)", "C(t)", oenv);
        for (j = 0; (j < n); j++)
        {
            fprintf(fp, "%10g %10g\n", tt[j], dos_vacf[j]);
        }
        fclose(fp);
    }
    if (NULL != mvacf)
    {
        fp = xvgropen(mvacf, "Mass-weighted velocity ACF",
                      "Time (ps)", "C(t)", oenv);
        print_legend(fp, dp, oenv);
        for (j = 0; (j < n); j++)
        {
            fprintf(fp, "%10g", tt[j]);
            for (t = 0; (t < V_SUM); t++)
            {
                if ((dp[t].npart > 0) && (dp[t].dof > 0))
                {
                    fprintf(fp, "  %10g", dp[t].alldos[j]);
                }
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

static void print_dos(gmx_bool bRecip, double recip_fac,
                      int n, real nu[], t_dosprops dp[],
                      output_env_t oenv)
{
    int   j, t;
    FILE *fp = xvgropen("DoS.xvg", "DoS",
                        bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)", "DoS(v)", oenv);

    print_legend(fp, dp, oenv);
    for (j = 0; (j < n); j++)
    {
        fprintf(fp, "%10g", recip_fac*nu[j]);
        for (t = 0; (t < V_SUM); t++)
        {
            if ((dp[t].npart > 0) && (dp[t].dof > 0))
            {
                fprintf(fp, " %10g", dp[t].alldos[j]/recip_fac);
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

static void print_dos_component(const char *fn, gmx_bool bRecip, double recip_fac,
                                int n,
                                real nu[], int t, t_dosprops dp[],
                                output_env_t oenv)
{
    FILE       *fp;
    int         j, k;
    char        buf[STRLEN];

    snprintf(buf, STRLEN, "%s-%s", velocityName[t], fn);

    fp = xvgropen(buf, "Density of states",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                  "\\f{4}S(\\f{12}n\\f{4})", oenv);
    xvgr_legend(fp, asize(cName), cName, oenv);
    for (j = 0; (j < n); j++)
    {
        fprintf(fp, "%10g", recip_fac*nu[j]);
        for (k = 0; (k < C_NR); k++)
        {
            fprintf(fp, "  %10g", dp[t].dd[LOT_QUANTUM][TM_2PT].dos[k][j]/recip_fac);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

typedef struct {
    int      ftype, fff, nesum;
    gmx_bool bCheckDrift;
    double   energy, stddev;
    double   esum, e2sum, esum0, e2sum0;
} t_energy_term;

void get_energy_terms(const char *enx_fn, int net, t_energy_term et[])
{
    int          i, nre, fff;
    gmx_enxnm_t *enm = NULL;
    gmx_bool     bCont;
    ener_file_t  fp = open_enx(enx_fn, "r");
    t_enxframe   frame;
    int          timecheck = 0;
    gmx_bool     b0set     = FALSE;
    gmx_int64_t  step0;
    int          nesum     = 0;

    for (i = 0; (i < net); i++)
    {
        et[i].nesum  = 0;
        et[i].energy = 0;
        et[i].stddev = 0;
        et[i].fff    = -1;
        et[i].esum0  = et[i].e2sum0 = 0;
        et[i].esum   = et[i].e2sum  = 0;
    }
    do_enxnms(fp, &nre, &enm);
    for (fff = 0; (fff < nre); fff++)
    {
        for (i = 0; (i < net); i++)
        {
            if (strcmp(interaction_function[et[i].ftype].longname, enm[fff].name) == 0)
            {
                et[i].fff = fff;
            }
        }
    }

    for (i = 0; (i < net); i++)
    {
        if (et[i].fff == nre)
        {
            fprintf(stderr, "No energy %s found in energy file %s\n",
                    interaction_function[et[i].ftype].longname, enx_fn);
            return;
        }
    }
    init_enxframe(&frame);

    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(fp, &frame);
            if (bCont)
            {
                timecheck = check_times(frame.t);
            }
        }
        while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            /* We read a valid frame, so we can use it */
            if (frame.nre > 0)
            {
                for (i = 0; (i < net); i++)
                {
                    if (et[i].fff < frame.nre)
                    {
                        if (!b0set)
                        {
                            step0  = frame.step;
                            b0set  = TRUE;
                        }
                        et[i].esum0  += frame.ener[et[i].fff].esum;
                        et[i].e2sum0 += frame.ener[et[i].fff].eav;
                        et[i].esum   += frame.ener[et[i].fff].e;
                        et[i].e2sum  += frame.ener[et[i].fff].e*frame.ener[et[i].fff].e;
                        et[i].nesum++;
                    }
                }
            }
        }
    }
    while (bCont);
    close_enx(fp);
    for (i = 0; (i < net); i++)
    {
        if (0 && b0set && (frame.step > step0))
        {
            et[i].stddev = (et[i].e2sum0)/(frame.step-step0);
            et[i].energy = (et[i].esum0)/(frame.step-step0);
        }
        else if (et[i].nesum > 0)
        {
            et[i].energy = et[i].esum/et[i].nesum;
            et[i].stddev = sqrt(et[i].e2sum/et[i].nesum - pow(et[i].energy, 2));
        }
    }
    free_enxnms(nre, enm);
    /* Printing a new line, just because the gromacs library prints step info while reading. */
    fprintf(stderr, "\n");
}

void get_energy_term(const char *enx_fn, int ftype, double *e, double *stddev)
{
    t_energy_term *et;

    snew(et, 1);
    et->ftype  = ftype;
    get_energy_terms(enx_fn, 1, et);
    *e      = et->energy;
    *stddev = et->stddev;
    sfree(et);
}

int gmx_dos(int argc, char *argv[])
{
    const char         *desc[] = {
        "[THISMODULE] computes the Density of States from a constant volume simulation.",
        "In order for this to be meaningful the coordinates and velocities must be saved",
        "in the trajectory with sufficiently high frequency such as to cover",
        "all vibrations, typically at most 5 fs.",
        "between saving. Properties based on the DoS are printed in the",
        "log file."
    };
    const char         *bugs[] = {
        "This program needs a lot of memory: total usage equals the number of atoms times 3 times number of frames times 4 (or 8 when run in double precision)."
    };
    /* Command line options */
    static     gmx_bool bVerbose    = TRUE, bAtomic = FALSE;
    static     gmx_bool bRecip      = FALSE, bDump = FALSE;
    static     real     toler       = 1e-6, Emd = 0;
    static     real     cPclassical = 0;
    static     real     cVclassical = 0;
    static     real     rot_symm    = 2;
    static     int      nthreads    = 0;
    static     int      nframes0    = 0;
    /*static     int      split       = 1;*/

    t_pargs             pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy." },
        { "-atomic", FALSE, etBOOL, {&bAtomic},
          "Do the Density of State calculation treating the system as independent atoms as well. This requires a lot of extra memory and computer time and is therefore by default turned off." },
        { "-recip", FALSE, etBOOL, {&bRecip},
          "Use cm^-1 on X-axis instead of 1/ps for DoS plots." },
        { "-Emd", FALSE, etREAL, {&Emd},
          "The average total energy in the MD simulation" },
        { "-cp",  FALSE, etREAL, {&cPclassical},
          "Reference heat capacity at constant pressure. This will be corrected using the quantum corrections to cV due to Berens et al."},
        { "-cv",  FALSE, etREAL, {&cVclassical},
          "Reference heat capacity at constant volume. This will be used instead of the one computed from the energy file."},
        { "-rot_symm", FALSE, etREAL, {&rot_symm},
          "Rotational symmetry in the molecule" },
        { "-nt",      FALSE, etINT, {&nthreads},
          "Number of threads to run the analysis (<= 0 means determine automatically)" },
        { "-nframes", FALSE, etINT, {&nframes0},
          "Number of frames in the trajectory (used to optimize memory allocation). <= 0 indicates, go figure it out by yourself." },
        /*{ "-split", FALSE, etINT, {&split},
           "Split the calculation in this many bits in order to limit the memory use. The trajectory will be read multiple times so this will be much slower." },*/
        { "-toler", FALSE, etREAL, {&toler},
          "[HIDDEN]Tolerance when computing the fluidicity using bisection algorithm" },
        { "-dump", FALSE, etBOOL, {&bDump},
          "[HIDDEN]Dump the y/fy plot corresponding to Fig. 2 inLin2003a and the and the weighting functions corresponding to Fig. 1 in Berens1983a." },
    };

    t_filenm            fnm[] = {
        { efTRN, "-f",    NULL,    ffREAD  },
        { efTPR, "-s",    NULL,    ffREAD  },
        { efNDX, NULL,    NULL,    ffOPTRD },
        { efEDR, "-enx",    NULL,      ffOPTRD  },
        { efXVG, "-vacf", "molvacf",  ffWRITE },
        { efXVG, "-mvacf", "molmvacf", ffWRITE },
        { efXVG, "-dos",  "moldos",   ffWRITE },
        { efLOG, "-g",    "moldos",   ffWRITE },
    };
#define NFILE asize(fnm)

    /* Regular variables */
    FILE               *fplog;
    gmx_mtop_t          mtop;
    t_topology          top;
    int                 Natom, Nmol, Natmxmol, NatomTotal;
    int                 ePBC;
    t_trxframe          fr;
    matrix              box;
    t_trxstatus        *status;
    int                 nV, nframes, n_alloc, fftcode;
    double              rho, V2sum, Vsum, Volume, VolError, massSystem;
    double           ***intensity, ***eigenvector;
    real               *dos_vacf, *nu, *tt, stddev;
    output_env_t        oenv;
    gmx_fft_t           fft2;
    t_dosprops         *dp;
    double              recip_fac, dt, t0, t1;
    int                 id, rdof;
    real                crT;
    rvec               *r1, *v1, *vecI;
    char                title[STRLEN];
    /* rotational temperatures */
    real                vib;
    int                 t, ai, nconstr, count;
    atom_id            *particle_index, *dummy_index;
    t_atom             *atom;
    real               *massMol;
    gmx_rmpbc_t         gpbc = NULL;
    const char         *enx_fn;
    double              fudge = 1.5;
    int                 N_tot, nfreq;
    double              recip0 = (1e7/SPEED_OF_LIGHT);

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc,
                           asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    if (bDump)
    {
        real beta = 1/(298.15*BOLTZ);
        printf("Dumping reference figures. Thanks for your patience.\n");
        dump_fy(oenv, toler);
        dump_w(oenv, beta);

        return 0;
    }

    fplog = gmx_fio_fopen(ftp2fn(efLOG, NFILE, fnm), "w");
    fprintf(fplog, "Doing density of states analysis based on trajectory.\n");
    please_cite(fplog, "Berens1983a");
    please_cite(fplog, "Pascal2011a");
    please_cite(fplog, "Caleman2011b");

    read_tps_conf(ftp2fn(efTPR, NFILE, fnm), title, &top, &ePBC, NULL, NULL, box,
                  TRUE);
    Volume    = det(box);
    recip_fac = bRecip ? recip0 : 1.0;

    /* Allocate stuff */
    snew(dp, V_NR);

    {
        /* Read topology, allocate memory and such. */
        t_tpxheader tpx;
        t_inputrec  ir;
        int         version, generation;

        read_tpxheader(ftp2fn(efTPR, NFILE, fnm), &tpx, TRUE, &version, &generation);
        snew(v1, tpx.natoms);
        snew(r1, tpx.natoms);

        ePBC = read_tpx(ftp2fn(efTPR, NFILE, fnm), &ir, box,
                        &Natom, r1, v1, NULL, &mtop);
        Nmol = mtop.molblock[0].nmol;
        if (mtop.nmoltype != 1)
        {
            fprintf(stderr, "WARNING: the system contains more than 1 molecule type.\n");
            fprintf(stderr, "Name: %s, nmoltype: %d\n", *mtop.name, mtop.nmoltype);
            fprintf(stderr, "Will only analyze the first molecule type (%d molecules).\n",
                    Nmol);
        }
    }

    enx_fn      = opt2fn_null("-enx", NFILE, fnm);
    if (NULL != enx_fn)
    {
        /* Get the energy from the energy file */
        t_energy_term *et;

        snew(et, 3);
        et[0].ftype       = F_PRES;
        et[1].ftype       = F_TEMP;
        et[2].ftype       = F_ETOT;
        et[1].bCheckDrift = et[2].bCheckDrift = TRUE;
        get_energy_terms(enx_fn, 3, et);

        if (cVclassical == 0)
        {
            cVclassical = KILO*(pow(et[2].stddev, 2))/(Nmol*BOLTZ*pow(et[1].energy, 2));
        }
        fprintf(fplog, "Read energy file %s.\n", enx_fn);
        fprintf(fplog, "Etot        = %10g +/- %10g\n", et[2].energy, et[2].stddev);
        fprintf(fplog, "T           = %10g +/- %10g\n", et[1].energy, et[1].stddev);
        fprintf(fplog, "P           = %10g +/- %10g\n", et[0].energy, et[0].stddev);
        fprintf(fplog, "cVclassical = %10g\n", cVclassical);
        if (opt2parg_bSet("-Emd", asize(pa), pa))
        {
            fprintf(stderr, "WARNING: You set both -Emd %g and -enx %s.\n", Emd, enx_fn);
            fprintf(stderr, "Using energy = %g from energy file.\n", et[2].energy);
        }
        Emd = et[2].energy;
        sfree(et);
    }
    else if (!opt2parg_bSet("-Emd", asize(pa), pa))
    {
        gmx_fatal(FARGS, "Neither the -Emd nor the -enx (preferred) option were set.\n"
                  "You need to supply at least one of these to get correct energies.");
    }
    else
    {
        fprintf(stderr, "WARNING: without energy file no correct cV can be computed.\n");
    }

    /* allocation of memory */
    snew(vecI,  Nmol);
    snew(particle_index, Natom);
    snew(dummy_index, Natom);
    snew(atom,  Natom);
    massSystem = 0;

    snew(massMol, mtop.nmoltype);
    {
        int k;
        for (k = 0; (k < mtop.nmoltype); k++)
        {
            int i;
            for (i = 0; (i < mtop.moltype[k].atoms.nr); i++)
            {
                massMol[k] += mtop.moltype[k].atoms.atom[i].m;
            }
        }
    }
    /* Re-count number of *real* atoms */
    Natom = 0;
    count = 0;
    {
        int j;
        for (j = 0; (j < mtop.nmolblock); j++)
        {
            int k, type = mtop.molblock[j].type;
            for (k = 0; (k < mtop.molblock[j].nmol); k++)
            {
                int i;
                for (i = 0; (i < mtop.moltype[type].atoms.nr); i++)
                {
                    if (eptAtom == mtop.moltype[type].atoms.atom[i].ptype)
                    {
                        /* we need to store the propertie of each atom */
                        particle_index[Natom] = count;
                        dummy_index[Natom]    = Natom;
                        atom[Natom]           = mtop.moltype[type].atoms.atom[i];
                        massSystem           += atom[Natom].m;
                        Natom++;
                    }
                    count++;
                }
            }
        }
    }
    Natmxmol = Natom/Nmol;

    /* Compute degrees of freedom */
    /* Should be MORE molecule specific! */
    nconstr = (gmx_mtop_ftype_count(&mtop, F_CONSTR) +
               3*gmx_mtop_ftype_count(&mtop, F_SETTLE))/Nmol;
    fprintf(fplog, "Number of constraints in molecule type 0 = %d\n", nconstr);
    dp[V_ATOM].dof  = max(0, 3*Natmxmol - nconstr);
    dp[V_T].dof     = min(3, dp[V_ATOM].dof);
    dp[V_A].dof     = min(3, dp[V_ATOM].dof - dp[V_T].dof);
    if (Natmxmol == 1)
    {
        dp[V_A].dof = 0;
    }
    else if (Natmxmol == 2)
    {
        dp[V_A].dof = min(2, dp[V_A].dof);
    }
    dp[V_R].dof = dp[V_A].dof;
    dp[V_V].dof = dp[V_ATOM].dof - dp[V_T].dof - dp[V_A].dof;

    /* allocation of memory */
    dp[V_T].npart = dp[V_A].npart = Nmol;
    dp[V_V].npart = Natom;
    /* To avoid allocating memory */
    dp[V_SUM].npart = 0;
    if (bAtomic)
    {
        dp[V_R].npart = dp[V_ATOM].npart = Natom;
    }

    {
        /* Pbc corrections.
         * Note that this routine breaks the mtop structure so it can not be
         * used afterwards!
         */
        top        = gmx_mtop_t_to_t_topology(&mtop);
        NatomTotal = top.atoms.nr;
        gpbc       = gmx_rmpbc_init(&top.idef, ePBC, NatomTotal);

        fprintf(fplog, "Natom  %d, Nmol %d, NatmXmol %d, NatomTotal %d.\n",
                Natom, Nmol, Natmxmol, NatomTotal);
    }

    read_first_frame(oenv, &status, ftp2fn(efTRN, NFILE, fnm),
                     &fr, TRX_NEED_V | TRX_NEED_X);
    t0 = fr.time;

    n_alloc = 0;
    nframes = 0;
    Vsum    = V2sum = VolError = 0;
    nV      = 0;

    for (t = 0; (t < V_NR); t++)
    {
        if (dp[t].npart > 0)
        {
            snew(dp[t].vec, dp[t].npart*DIM);
        }
    }
    if (nthreads > 1)
    {
        gmx_omp_set_num_threads(nthreads);
    }
    else
    {
        nthreads  = gmx_omp_get_max_threads();
    }
    fprintf(fplog, "Running the analysis on %d threads\n", nthreads);
    /* Allocate temp variables for principal component calc */
    {
        int i, j;
        snew(intensity, nthreads);
        snew(eigenvector, nthreads);
        for (j = 0; (j < nthreads); j++)
        {
            snew(intensity[j], 4);
            snew(eigenvector[j], 4);
            for (i = 0; (i < 4); i++)
            {
                snew(intensity[j][i], 4);
                snew(eigenvector[j][i], 4);
            }
        }
    }

    /* Read each frame, save coordinates and velocities,
     * do calculation and save results.
     */
    do
    {
        if (fr.bBox)
        {
            Volume = det(fr.box);
            V2sum += Volume*Volume;
            Vsum  += Volume;
            nV++;
        }

        if (nframes >= n_alloc)
        {
            if ((n_alloc == 0) && (nframes0 > 0))
            {
                n_alloc += fudge*(nframes0 + 16);
            }
            else
            {
                n_alloc += 12800;
            }
            for (t = 0; (t < V_NR); t++)
            {
#pragma omp parallel
                {
                    int i, i0, i1, thread_id, NN, NT;

                    NT        = nthreads;
                    NN        = dp[t].npart*DIM;
                    thread_id = gmx_omp_get_thread_num();
                    i0        = thread_id*NN/NT;
                    i1        = min(NN, (thread_id+1)*NN/NT);

                    for (i = i0; i < i1; i++)
                    {
                        srenew(dp[t].vec[i], n_alloc);
                    }
                }
                // End parallel section.
            }
        }
        /* Remove periodicity */
        gmx_rmpbc_copy(gpbc, NatomTotal, box, fr.x, fr.x);

        /* Read velocities and coordinates of all atoms:
         * each atoms is pointed to by j
         */
#pragma omp parallel
        {
            int j, j0, j1, thread_id, NN, NT;

            NT        = nthreads;
            NN        = Natom;
            thread_id = gmx_omp_get_thread_num();
            j0        = thread_id*NN/NT;
            j1        = min(NN, (thread_id+1)*NN/NT);

            for (j = j0; j < j1; j++)
            {
                int  aj  = particle_index[j];

                /* load atomic positions and velocities */
                copy_rvec(fr.v[aj], v1[j]);
                copy_rvec(fr.x[aj], r1[j]);

                if (bAtomic)
                {
                    int k;
                    for (k = 0; (k < DIM); k++)
                    {
                        /* save total velocities */
                        dp[V_ATOM].vec[DIM*j+k][nframes] = v1[j][k];
                    }
                }
            }
        }
        // End parallel section
        /* Dividing atoms in molecules, id. Do it in parallel using OpenMP */
#pragma omp parallel
        {
            int    id, i0, i1, thread_id, NN, NT;
            matrix trans;
            rvec   rT, tmp, I, angmom, pomega, omega, vcm, xcm, vr;

            NT        = nthreads;
            NN        = Nmol;
            thread_id = gmx_omp_get_thread_num();
            i0        = thread_id*NN/NT;
            i1        = min(NN, (thread_id+1)*NN/NT);

            for (id = i0; id < i1; id++) /* Nmol: number of molecules */
            {
                int  k, lid  = id*DIM;
                int *d_index = &dummy_index[id*Natmxmol];

                /* Find center of mass and shift all position to this new origin:
                 * it must be done for each atom in the molecule.
                 * Returns total mass of molecule
                 */
                dsub_xcm(r1, Natmxmol, d_index, atom, xcm);
                /* find the velocity of center of mass and shift all position to
                 * this new origin: must be done for each molecule
                 */
                dsub_xcm(v1, Natmxmol, d_index, atom, vcm);
                /* save translational velocities for the molecule id */
                for (k = 0; (k < DIM); k++)
                {
                    dp[V_T].vec[lid+k][nframes] = vcm[k];
                }

                if (Natmxmol > 1.0)
                {
                    int ai, j;
                    /* compute principal moment of inertia
                     * return trans eigenvector tensor and I the eigenvalues
                     * (principal moment of inertia)
                     */
                    principal(Natmxmol, d_index, atom, r1, trans, I,
                              intensity[thread_id], eigenvector[thread_id]);
                    /* Save the principal moment: it will be used to compute the
                     * rotational temperatures.
                     */
                    rvec_inc(vecI[id], I);

                    /* we run through the j-th atom of id-th molecule and we
                     * save velocities for i-th atoms from d_index
                     * reset the needed variables
                     */
                    clear_rvec(angmom);
                    clear_rvec(pomega);
                    clear_rvec(omega);

                    for (j = 0; j < Natmxmol; j++)
                    {
                        int ai = d_index[j];
                        /* compute the angular momentum angmom */
                        clear_rvec(tmp);
                        cprod(r1[ai], v1[ai], tmp);
                        /* mass weigth it for each atom in the molecule */
                        angmom[XX] += atom[ai].m*tmp[XX];
                        angmom[YY] += atom[ai].m*tmp[YY];
                        angmom[ZZ] += atom[ai].m*tmp[ZZ];
                    }
                    /* compute angular velocity along principle axis, first step */
                    for (j = 0; (j < DIM); j++)
                    {
                        if (I[j] > 0.0)
                        {
                            for (k = 0; k < 3; k++)
                            {
                                pomega[j] += angmom[k]*trans[k][j];
                            }
                            pomega[j] /= I[j];
                        }
                    }
                    /* calculate angular velocities. Here we use the transpose of the trans
                     * matrix by swapping the indexing
                     */
                    for (j = 0; (j < 3); j++)
                    {
                        for (k = 0; k < 3; k++)
                        {
                            omega[j] += pomega[k]*trans[j][k];
                        }
                    }
                    /* calculate inertia weighted angular velocities and save */
                    for (j = 0; (j < DIM); j++)
                    {
                        real vvv = 0;
                        for (k = 0; k < 3; k++)
                        {
                            if (I[k] > 0)
                            {
                                vvv += pomega[k]*trans[j][k]*sqrt(I[k]);
                            }
                        }
                        /* save angular velocities of molecule id,
                         * weighted by sqrt(I)
                         */
                        dp[V_A].vec[id*DIM+j][nframes] = vvv;
                    }

                    for (j = 0; j < Natmxmol; j++)
                    {
                        ai = d_index[j];
                        /* calculate velocity due to rotation vr = w x r */
                        cprod(omega, r1[ai], vr);
                        for (k = 0; (k < DIM); k++)
                        {
                            if (bAtomic)
                            {
                                /* save rotational velocities */
                                dp[V_R].vec[ai*DIM+k][nframes] = vr[k];
                            }
                            /* calculate vibrational velocities and save
                             * v1 is the relative velocity versus center of mass
                             */
                            dp[V_V].vec[ai*DIM+k][nframes] = v1[ai][k] - vr[k];
                        }
                    }
                }
            }
        }
        // End parallel section

        t1 = fr.time;
        nframes++;
    }
    while (read_next_frame(oenv, status, &fr));
    N_tot = fudge*nframes;
    nfreq = 1+N_tot/2;
    printf("N_tot = %d n_alloc = %d\n", N_tot, n_alloc);

    gmx_rmpbc_done(gpbc);

    sfree(r1);   /* empty r1, we need memory */
    sfree(v1);   /* empty v1, we need memory */
    close_trj(status);

    {
        /* Normalize moment of inertia */
        real ddd = 1.0/nframes;
        rvec sumVec;
        clear_rvec(sumVec);
        for (id = 0; id < Nmol; id++)
        {
            svmul(ddd, vecI[id], vecI[id]);
            rvec_inc(sumVec, vecI[id]);
        }
        if (NULL != debug)
        {
            pr_rvecs(debug, 0, "Average vecI", vecI, Nmol);
            svmul(1.0/Nmol, sumVec, sumVec);
            pr_rvec(debug, 0, "Overall Average vecI", sumVec, DIM, FALSE);
        }
    }

    dt = (t1-t0)/(nframes-1);
    if (nV > 0)
    {
        double V2aver = V2sum/nV;
        double Vol2   = Volume*Volume;
        Volume   = Vsum/nV;
        VolError = 0;
        if (V2aver > Vol2)
        {
            VolError = sqrt(V2aver - Vol2);
        }
        fprintf(fplog, "Vsum = %g V2sum = %g Volume = %g VolError = %g\n",
                Vsum, V2sum, Volume, VolError);
    }

    snew(dos_vacf, N_tot);
    {
        int nfft = 0;
        for (t = 0; (t < V_NR); t++)
        {
            int i;
            /* Set the particle mass */
            if (dp[t].npart > 0)
            {
                dp[t].partMass = massSystem/dp[t].npart;
            }

            /* Padding the vec array with zeroes */
#pragma omp parallel
            {
                int i, i0, i1, thread_id, NN, NT;

                NT        = nthreads;
                NN        = DIM*dp[t].npart;
                thread_id = gmx_omp_get_thread_num();
                i0        = thread_id*NN/NT;
                i1        = min(NN, (thread_id+1)*NN/NT);
                for (i = i0; (i < i1); i++)
                {
                    int j;
                    if (N_tot > n_alloc)
                    {
                        srenew(dp[t].vec[i], N_tot);
                    }
                    for (j = nframes; (j < N_tot); j++)
                    {
                        dp[t].vec[i][j] = 0;
                    }
                }
            }
            // End parallel section
            nfft += DIM*dp[t].npart;
        }
        if (bVerbose)
        {
            printf("Going to determine %d correlation functions of length %d. Hang on.\n",
                   nfft, N_tot);
            printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");
        }
    }
    /* Allocate arrays */
    for (t = 0; (t < V_NR); t++)
    {
        init_dos_data(&dp[t], N_tot);
    }

    for (t = 0; (t < V_NR); t++)
    {
        if ((dp[t].dof > 0) && (dp[t].npart > 0))
        {
            int i, fftcode;
            fftcode = many_auto_correl(DIM*dp[t].npart, nframes, N_tot, dp[t].vec);
            if (0 != fftcode)
            {
                return -1;
            }
#pragma omp parallel
            //TODO: Swap i and j loops for parallellization.
            {
                int i, j, j0, j1, thread_id, NN, NT;

                NT        = nthreads;
                NN        = N_tot;
                thread_id = gmx_omp_get_thread_num();
                j0        = thread_id*NN/NT;
                j1        = min(NN, (thread_id+1)*NN/NT);

                for (j = j0; (j < j1); j++)
                {
                    for (i = 0; (i < dp[t].npart); i++)
                    {
                        int  ai    = i*DIM;
                        real v1j, mass  = atom[i].m;
                        if (V_T == t)
                        {
                            /* NOTE: should be molecule specific */
                            mass = massMol[0];
                        }
                        else if (V_A == t)
                        {
                            mass = 1;
                        }

                        v1j = (dp[t].vec[ai+XX][j] +
                               dp[t].vec[ai+YY][j] +
                               dp[t].vec[ai+ZZ][j]);
                        dp[t].alldos[j] += v1j*mass;
                        if (V_ATOM == t)
                        {
                            dos_vacf[j] += v1j/dp[t].npart;
                        }
                    }
                }
            }
            // End parallel section
        }
    }

    /* Compute the temperatures. The number of degrees of freedom
     * are per molecule, therefore we divide the DoS by the number
     * of molecules.
     * From here, the temperature is taken from simulation.
     */
    for (t = 0; (t < V_NR); t++)
    {
        if ((dp[t].dof > 0) && (dp[t].npart > 0))
        {
            dp[t].T    = dp[t].alldos[0]/(dp[t].dof*BOLTZ*Nmol);
            dp[t].beta = 1.0/(dp[t].T*BOLTZ);
            dp[t].bfac = 2.0*dt*dp[t].beta;
        }
    }
    /* Make time- and frequency arrays */
    snew(tt, N_tot);
    snew(nu, N_tot);
    {
        int    j;
        double pwrfreq = 1/(dt*(N_tot+1));
        for (j = 0; (j < N_tot); j++)
        {
            tt[j] = j*dt;
            nu[j] = j*pwrfreq;
        }
    }
    print_stuff(N_tot, tt, dp, dos_vacf, opt2fn_null("-vacf", NFILE, fnm),
                opt2fn_null("-mvacf", NFILE, fnm), oenv);
    if ((fftcode = gmx_fft_init_1d(&fft2, N_tot,
                                   GMX_FFT_FLAG_CONSERVATIVE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d returned %d", fftcode);
    }

    /* compute density of state for each kind of velocities,
     * return in same vector
     */
    for (t = 0; t < V_NR; t++)
    {
        if ((dp[t].dof > 0) && (dp[t].npart > 0))
        {
            typedef real complex[2];
            complex *in, *out;
            int      j;
            snew(in, N_tot);
            snew(out, N_tot);
            for (j = 0; (j < N_tot); j++)
            {
                in[j][0] = dp[t].bfac*dp[t].alldos[j];
                in[j][1] = 0;
            }
            for (; (j < N_tot); j++)
            {
                in[j][0] = in[j][1] = 0;
            }
            /* The 2pt code uses a backword transform, but the articles
             * write a forward transform. However for real data it
             * does not affect the real component of the output.
             */
            if ((fftcode = gmx_fft_1d(fft2, GMX_FFT_BACKWARD,
                                      (void *)in, (void *)out)) != 0)
            {
                gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
            }
            for (j = 0; (j < nfreq); j++)
            {
                dp[t].alldos[j] = sqrt(sqr(out[j][0])+sqr(out[j][1]));
            }
            sfree(in);
            sfree(out);
        }
    }
    gmx_fft_destroy(fft2);

    print_dos(bRecip, recip_fac, nfreq, nu, dp, oenv);

    /* Density */
    rho   = (massSystem*AMU)/(Volume*NANO*NANO*NANO);

    /* Analyze the DoS[t] */
    for (t = 0; (t < V_NR); t++)
    {
        if ((dp[t].dof == 0) || (dp[t].npart == 0))
        {
            continue;
        }
        dp[t].DoS0 = dp[t].alldos[0];

        if (V_V != t)
        {
            int j;
            /* eq 13 JPC B, 114, 8191 (2010) */
            dp[t].nrho  = dp[t].npart/Volume;
            dp[t].Delta = ((2*dp[t].DoS0/(9*dp[t].npart))*sqrt(M_PI*BOLTZ*dp[t].T/dp[t].partMass)*
                           pow(dp[t].nrho, 1.0/3.0)*pow(6/M_PI, 2.0/3.0));
            dp[t].f     = calc_fluidicity(dp[t].Delta, toler);
            dp[t].ry    = calc_y(dp[t].f, dp[t].Delta, toler);

            for (j = 0; (j < nfreq); j++)
            {
                int  mm;
                real dj = (dp[t].DoS0/(1+sqr(dp[t].DoS0*M_PI*nu[j]/
                                             (6*dp[t].f*dp[t].npart))));
                dp[t].dd[LOT_CLASSICAL][TM_2PT].dos[C_GAS][j] = dj;
                dp[t].dd[LOT_QUANTUM][TM_2PT].dos[C_GAS][j]   = dj;
                /* Difference between the two */
                dp[t].dd[LOT_QUANTUM_CORRECTION][TM_2PT].dos[C_GAS][j] = 0;
            }

            /* Determine hard sphere degrees of freedom hsdf.
             * Since all the 2PT variants have the same diffusivity, we just take the first.
             */
            dp[t].hsdf  = evaluate_integral(nfreq, nu, dp[t].dd[LOT_QUANTUM][TM_2PT].dos[C_GAS],
                                            NULL, 0, &stddev);
            /* Determine total degrees of freedom ttdf */
            dp[t].ttdf  = evaluate_integral(nfreq, nu, dp[t].alldos, NULL, 0, &stddev);

            dp[t].y     = dp[t].ry*(dp[t].hsdf/dp[t].ttdf);

            dp[t].fy    = dp[t].f * dp[t].y;

            dp[t].z     = calc_compress(dp[t].y);
            /* Lin 2010, eq. 16, Sackur-Tetrode equation */
            {
                double xlog = (pow(2*M_PI*dp[t].partMass*BOLTZ*dp[t].T/
                                   (sqr(PLANCK)), 1.5)*Volume/(dp[t].f * dp[t].npart));
                dp[t].Sig   = BOLTZ*(2.5+log(xlog));
                fprintf(fplog, "m %g k %g T %g h %g V %g N %d xlog %g\n",
                        dp[t].partMass, BOLTZ, dp[t].T, PLANCK, Volume, dp[t].npart, xlog);
            }
            dp[t].Shs   = dp[t].Sig+calc_Shs(dp[t].y);
            dp[t].sigHS = pow(6*dp[t].ry*Volume/(M_PI*dp[t].npart), 1.0/3.0);
        }
        else /* t == V_V */
        {
            int j;
            /* All other parameters are 0 */
            dp[t].y     = 1;
            for (j = 0; (j < nfreq); j++)
            {
                int lot;
                for (lot = 0; (lot < LOT_NR); lot++)
                {
                    dp[t].dd[lot][TM_2PT].dos[C_GAS][j] = 0;
                }
            }
        }
        /* Now compute solid (2) component, for 2PT only.
         * Fill the total arrays for the other methods.
         */
        {
            int j;
            for (j = 0; (j < nfreq); j++)
            {
                int lot, tm;
                for (lot = 0; (lot < LOT_NR); lot++)
                {
                    for (tm = 0; (tm < TM_NR); tm++)
                    {
                        dp[t].dd[lot][tm].dos[C_TOTAL][j] = dp[t].alldos[j];
                        if (TM_2PT == tm)
                        {
                            dp[t].dd[lot][tm].dos[C_SOLID][j] = (dp[t].dd[lot][tm].dos[C_TOTAL][j]-
                                                                 dp[t].dd[lot][tm].dos[C_GAS][j]);
                        }
                    }
                }
            }
        }
        dp[t].dostot = evaluate_integral(nfreq, nu, dp[t].alldos, NULL, 0, &stddev);

        print_dos_component(opt2fn("-dos", NFILE, fnm), bRecip, recip_fac,
                            nfreq, nu, t, dp, oenv);

        /* Finally analyze the results! */
        {
            int i;
            for (i = 0; (i < DOS_NR); i++)
            {
                dp[t].wDiff[i] = 0;
            }
        }
        if (V_A == t)
        {
            int          i;
            const double Q = PLANCK*PLANCK/(8.0*M_PI*M_PI*BOLTZ);
            rvec         I, rT;
            clear_rvec(rT);
            clear_rvec(I);

            for (i = 0; i < dp[t].npart; i++)   /* run for each molecule */
            {
                int k;
                for (k = 0; (k < DIM); k++)
                {
                    if (vecI[i][k] > toler)
                    {
                        rT[k] += Q/vecI[i][k];
                        I[k]  += vecI[i][k];
                    }
                }
            }
            svmul(1.0/dp[t].npart, I,  I);
            svmul(1.0/dp[t].npart, rT, rT);

            fprintf(fplog, "Momenta of Inertia Ix %g Iy %g Iz %g amu nm^2\n",
                    I[XX], I[YY], I[ZZ]);
            fprintf(fplog, "Characteristic Rotational T x %g y %g z %g K\n",
                    rT[XX], rT[YY], rT[ZZ]);

            crT  = 1.0; /* characteristic rotational temperature */
            rdof = 0;   /* rotational degrees of freedom */

            for (i = 0; i < DIM; i++)
            {
                if (rT[i] > 0.0)
                {
                    crT *= rT[i];
                    rdof++;
                }
            }
            if (rdof == 0)
            {
                crT = -1.0;
            }
            {
                real SR  = (3.0/2.0 + log(sqrt(M_PI*pow(dp[V_A].T, rdof)/crT)/rot_symm))/3;
                dp[t].wDiff[DOS_S]  = SR;
                fprintf(fplog, "SR = %g, rdof = %d, crT = %g rot_symm = %g\n", SR, rdof, crT, rot_symm);
            }
        }
        else
        {
            dp[t].wDiff[DOS_S] = dp[t].Shs/(3*BOLTZ);
        }
        if (V_V != t)
        {
            dp[t].wDiff[DOS_CV] = 0.5;
            dp[t].wDiff[DOS_E]  = 0.5;
            dp[t].wDiff[DOS_A]  = dp[t].wDiff[DOS_E] - dp[t].wDiff[DOS_S];
        }
        {
            real   fff[DOS_NR];
            int    dos;
            fff[DOS_CV] = fff[DOS_S] = (BOLTZ*KILO/Nmol);
            fff[DOS_A]  = fff[DOS_E] = (BOLTZ*dp[t].T/Nmol);
            for (dos = 0; (dos < DOS_NR); dos++)
            {
                int j, lot;
                for (j = 0; (j < nfreq); j++)
                {
                    for (lot = 0; (lot < LOT_QUANTUM_CORRECTION); lot++)
                    {
                        int tm;
                        for (tm = 0; (tm < TM_NR); tm++)
                        {
                            if (TM_2PT == tm)
                            {
                                dp[t].dd[lot][tm].wdos[dos][C_SOLID][j] =
                                    (dp[t].dd[lot][tm].dos[C_SOLID][j] *
                                     weight(dos, lot, tm, C_SOLID, nu[j], dp[t].beta));
                                dp[t].dd[lot][tm].wdos[dos][C_GAS][j]   =
                                    (dp[t].dd[lot][tm].dos[C_GAS][j] *
                                     dp[t].wDiff[dos]);
                                dp[t].dd[lot][tm].wdos[dos][C_TOTAL][j] =
                                    (dp[t].dd[lot][tm].wdos[dos][C_SOLID][j] +
                                     dp[t].dd[lot][tm].wdos[dos][C_GAS][j]);
                            }
                            else
                            {
                                dp[t].dd[lot][tm].wdos[dos][C_TOTAL][j] =
                                    (dp[t].dd[lot][tm].dos[C_TOTAL][j] *
                                     weight(dos, lot, tm, C_TOTAL, nu[j], dp[t].beta));
                            }
                        }
                        for (tm = 0; (tm < TM_NR); tm++)
                        {
                            int cc;
                            for (cc = 0; (cc < C_NR); cc++)
                            {
                                dp[t].dd[LOT_QUANTUM_CORRECTION][tm].wdos[dos][cc][j] =
                                    (dp[t].dd[LOT_QUANTUM][tm].wdos[dos][cc][j] -
                                     dp[t].dd[LOT_CLASSICAL][tm].wdos[dos][cc][j]);
                            }
                        }
                    }
                }
                for (lot = 0; (lot < LOT_NR); lot++)
                {
                    int tm;
                    for (tm = 0; (tm < TM_NR); tm++)
                    {
                        int cc;
                        for (cc = 0; (cc < C_NR); cc++)
                        {
                            dp[t].dd[lot][tm].X[dos][cc] =
                                fff[dos]*evaluate_integral(nfreq, nu,
                                                           dp[t].dd[lot][tm].wdos[dos][cc],
                                                           NULL, 0, &stddev);
                        }
                    }
                }
            }
        }

        if (V_V != t)
        {
            dp[t].DiffACF = (1000/3.0)*evaluate_integral(nfreq, tt, dos_vacf, NULL,
                                                         0, &stddev);
            dp[t].DiffDoS = 1000*dp[t].DoS0/(12*massSystem*dp[t].beta);
        }
    }
    for (t = 0; (t <= V_ATOM); t++)
    {
        int    lot, tm;
        double E0[C_NR];
        double dof_frac;

        dof_frac  = (dp[t].dof/dp[V_ATOM].dof);
        dp[t].Emd = Emd*dof_frac/Nmol;
        dp[t].E0  = dp[t].Emd - dp[t].dd[LOT_CLASSICAL][TM_1PT].X[DOS_E][C_TOTAL];

        E0[C_SOLID] = (1-dp[t].f)*dp[t].E0;
        E0[C_GAS]   = dp[t].f*dp[t].E0;
        E0[C_TOTAL] = dp[t].E0;
        /* Since we add E0 to both QM and Classical
         * the difference to be added to the quantum corrections is 0.
         */
        for (lot = 0; (lot < LOT_QUANTUM_CORRECTION); lot++)
        {
            for (tm = 0; (tm < TM_NR); tm++)
            {
                int cc;
                for (cc = 0; (cc < C_NR); cc++)
                {
                    dp[t].dd[lot][tm].X[DOS_A][cc] += E0[cc];
                    dp[t].dd[lot][tm].X[DOS_E][cc] += E0[cc];
                }
            }
        }

        print_dp(fplog, t, &dp[t]);
    }

    /* Sum the energies etc. */
    for (t = 0; (t <= V_V); t++)
    {
        int dos;
        dp[V_SUM].dostot += dp[t].dostot;
        dp[V_SUM].DoS0   += dp[t].DoS0;
        for (dos = 0; (dos < DOS_NR); dos++)
        {
            int lot;
            for (lot = 0; (lot < LOT_NR); lot++)
            {
                int tm;
                for (tm = 0; (tm < TM_NR); tm++)
                {
                    int cc;
                    for (cc = 0; (cc < C_NR); cc++)
                    {
                        dp[V_SUM].dd[lot][tm].X[dos][cc] += dp[t].dd[lot][tm].X[dos][cc];
                    }
                }
            }
        }
        dp[V_SUM].Emd    += dp[t].Emd;
        dp[V_SUM].E0     += dp[t].E0;
        dp[V_SUM].dof    += dp[t].dof;
        dp[V_SUM].T      += dp[t].T*dp[t].dof;
    }
    dp[V_SUM].T /= dp[V_SUM].dof;
    /* ---------------------- */
    fprintf(fplog, "\n\t\t Final resume \n");
    fprintf(fplog, "System  = %s\n", *top.name);
    fprintf(fplog, "Nmol    = %d\n", Nmol);
    fprintf(fplog, "Natom   = %d\n", Natom);
    fprintf(fplog, "dt      = %g ps\n", dt);
    fprintf(fplog, "tmass   = %g amu\n", massSystem);
    if (VolError/Volume > toler)
    {
        fprintf(fplog, "Volume  = %g +/- %g nm^3\n", Volume, VolError);
        fprintf(fplog, "WARNING: The two phase thermodynamics method may not work well\n");
        fprintf(fplog, "in constant pressure simulations.\n");
    }
    fprintf(fplog, "Volume  = %g nm^3\n", Volume);
    fprintf(fplog, "Density = %g g/l\n", rho);
    fprintf(fplog, "Emd     = %g kJ/mol/system\n", Emd);
    fprintf(fplog, "bRecip  = %s\n", bool_names[bRecip]);

    {
        int         k;
        const char *items[] = {
            "", "Temperature", "DOF", "Dos integral", "DoS0",
            "E0", "Emd", "fluidicity",
        };

        for (k = 0; (k < asize(items)); k++)
        {
            fprintf(fplog, "%-29s", items[k]);
            for (t = 0; (t <= V_SUM); t++)
            {
                if ((t != V_SUM) &&
                    ((dp[t].dof == 0) || (dp[t].npart == 0)))
                {
                    continue;
                }
                switch (k)
                {
                    case 0:
                        fprintf(fplog, "  %10s", velocityName[t]);
                        break;
                    case 1:
                        fprintf(fplog, "  %10g", dp[t].T);
                        break;
                    case 2:
                        fprintf(fplog, "  %10g", dp[t].dof);
                        break;
                    case 3:
                        fprintf(fplog, "  %10g", dp[t].dostot);
                        break;
                    case 4:
                        fprintf(fplog, "  %10g", dp[t].DoS0);
                        break;
                    case 5:
                        fprintf(fplog, "  %10g", dp[t].E0);
                        break;
                    case 6:
                        fprintf(fplog, "  %10g", dp[t].Emd);
                        break;
                    case 7:
                        fprintf(fplog, "  %10g", dp[t].f);
                        break;
                }
            }
            fprintf(fplog, "\n");
        }
        {
            int dos;
            for (dos = 0; (dos < DOS_NR); dos++)
            {
                int lot, tm;
                for (tm = 0; (tm < TM_NR); tm++)
                {
                    for (lot = 0; (lot < LOT_NR); lot++)
                    {
                        int j, j0 = C_TOTAL;
                        if (TM_2PT == tm)
                        {
                            j0 = 0;
                        }
                        for (j = j0; (j < C_NR); j++)
                        {

                            fprintf(fplog, "%-2s  %-10s  %-6s %-6s",
                                    dosShort[dos], lotName[lot], tmName[tm], cName[j]);
                            for (t = 0; (t <= V_SUM); t++)
                            {
                                if ((t != V_SUM) &&
                                    ((dp[t].dof == 0) || (dp[t].npart == 0)))
                                {
                                    continue;
                                }
                                fprintf(fplog, "  %10g", dp[t].dd[lot][tm].X[dos][j]);
                            }
                            fprintf(fplog, " %s\n", dosUnit[dos]);
                        }
                    }
                }
            }
        }
    }
    sfree(atom);
    fprintf(fplog, "\nRecommended thermodynamics values:\n");
    fprintf(fplog, "----------------------------------\n");
    fprintf(fplog, "%-16s %2s = %10.3f %s (%s/%s method)\n",
            dosName[DOS_S], dosShort[DOS_S],
            dp[V_SUM].dd[LOT_QUANTUM][TM_2PT].X[DOS_S][C_TOTAL],
            dosUnit[DOS_A],
            tmName[TM_2PT], lotName[LOT_QUANTUM]);

    if (cVclassical > 0)
    {
        fprintf(fplog, "%-16s %2s = %10.3f %s (%s/%s method)\n",
                dosName[DOS_CV], dosShort[DOS_CV],
                cVclassical+dp[V_SUM].dd[LOT_QUANTUM_CORRECTION][TM_1PT].X[DOS_CV][C_TOTAL],
                dosUnit[DOS_CV],
                tmName[TM_MD], lotName[LOT_QUANTUM_CORRECTION]);
        fprintf(fplog, "Note: that heat capacities in MD simulations converge poorly.\n");
        fprintf(fplog, "Please check the convergence using gmx fluctprops.\n");
    }
    else
    {
        fprintf(fplog, "No cV based on the MD simulation. Please rerun with the -enx option.\n");
    }
    if (cPclassical > 0)
    {
        fprintf(fplog, "%-16s %2s = %10.3f %s (%s/%s method)\n",
                "Heat capacity", "cP",
                cPclassical+dp[V_SUM].dd[LOT_QUANTUM_CORRECTION][TM_1PT].X[DOS_CV][C_TOTAL],
                dosUnit[DOS_CV],
                tmName[TM_MD], lotName[LOT_QUANTUM_CORRECTION]);
    }


    fprintf(fplog, "\nNote that the Helmholtz energy and the internal energy\n");
    fprintf(fplog, "include the MD energy which usually is on an unknown scale\n");
    fprintf(fplog, "meaning the energies can not be compared to experiment.\n");
    fprintf(fplog, "%-16s %2s = %10.3f %s (%s/%s method)\n",
            dosName[DOS_A], dosShort[DOS_A],
            dp[V_SUM].dd[LOT_QUANTUM][TM_2PT].X[DOS_A][C_TOTAL],
            dosUnit[DOS_A],
            tmName[TM_2PT], lotName[LOT_QUANTUM]);

    fprintf(fplog, "%-16s %2s = %10.3f %s (%s/%s method)\n",
            dosName[DOS_E], dosShort[DOS_E],
            dp[V_SUM].Emd + dp[V_SUM].dd[LOT_QUANTUM_CORRECTION][TM_1PT].X[DOS_E][C_TOTAL],
            dosUnit[DOS_CV],
            tmName[TM_MD], lotName[LOT_QUANTUM_CORRECTION]);
    fprintf(fplog, "\nArrivederci!\n");
    fclose(fplog);
    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    return 0;
}
