/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include "dosutils.h"

#include <cmath>
#include <cstring>

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

const char *lotName(LevelOfTheory lot)
{
    static const char *lotname[LOT_NR] = {
        "Classical", "Quantum", "Quant Corr"
    };

    switch (lot)
    {
        case LOT_CLASSICAL:
            return lotname[0];
        case LOT_QUANTUM:
            return lotname[1];
        case LOT_QUANTUM_CORRECTION:
            return lotname[2];
        case LOT_NR:
            gmx_fatal(FARGS, "lotName called with wrong argument");
    }
    return nullptr;
};


/* Utility to do double precision subtract center of mass */
void dsub_xcm(rvec x[], int gnx, int *index, t_atom atom[], rvec xcm)
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

double FD(double Delta, double f) /* eq 34 JCP 119, 11792 (2003) */
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
double YYY(double f, double y)
{
    //return (2*pow(y*f, 3) - gmx::square(f)*y*(1+6*y) + (2+6*y)*f - 2);
    double yf = y*f;
    return 2.0*yf*yf*yf - f*yf*(1.0+6.0*y) + (2.0+6.0*y)*f - 2.0;
}

/* compressibility of hard spheres z=(1+y+y^2-y^3)/(1-y)^3 */
double calc_compress(double y)
{
    //    return ((1+y+gmx::square(y)-pow(y, 3))/(pow((1-y), 3)));
    if (y == 0)
    {
        return 1;
    }
    else
    {
        return (1 + y + y*y - y*y*y)/((1-y) * (1-y) * (1-y));
    }
}

double bisector(double Delta, double tol,
                double ff0, double ff1,
                double ff(double, double))
{
    double fd, f;
    double tolmin = 1e-8;

    printf("bisector: Delta %g tol %g ff0 %g ff1 %g\n", Delta, tol, ff0, ff1);
    if (tol < tolmin)
    {
        fprintf(stderr, "Unrealistic tolerance %g for bisector. Setting it to %g\n",
                tol, tolmin);
        tol = tolmin;
    }

    do
    {
        f   = (ff0+ff1)*0.5;
        fd  = ff(Delta, f);
        if (fd < 0)
        {
            ff0 = f;
        }
        else if (fd > 0)
        {
            ff1 = f;
        }
        else
        {
            return f;
        }
    }
    while ((ff1-ff0) > tol);

    return f;
}

/* calculate fluidicity f */
double calc_fluidicity(double Delta, double tol)
{
    /* solution for FD */
    return bisector(Delta, tol, 0, 1, FD);
}

/* hard sphere packing fraction, eq 32 JCP 119, 11792 (2003) */
double calc_y(double f, double Delta, double toler)
{
    double y1, y2;

    y1 = pow(f/Delta, 1.5);
    y2 = bisector(f, toler, 0, 10000, YYY);
    if (fabs((y1-y2)/(y1+y2)) > 100*toler)
    {
        fprintf(stderr, "Inconsistency computing y: y1 = %f, y2 = %f, using y1.\n",
                y1, y2);
    }

    return y2;
}

/* entropy for hard spheres */
double calc_Shs(double fy)
{
    double z = calc_compress(fy);

    if (0 == z)
    {
        fprintf(stderr, "Zero compressibility found in calc_Shs\n");
        return 0;
    }
    else if (1 <= fy)
    {
        fprintf(stderr, "Hard sphere packing fraction is %g (> 1!) in calc_Shs\n", fy);
        return 0;
    }
    else
    {
        return BOLTZ*(log(z) + fy*(3*fy-4)/gmx::square(1-fy));
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
                return u/2 + u/(std::expm1(u));
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
            return u/(std::expm1(u)) - log(1-exp(-u));
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
                koko = gmx::square(1-ebn);
                return gmx::square(u)*ebn/koko;
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

real weight(int Density, int LevelOfTheory,
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

static void dump_fy(gmx_output_env_t *oenv, real toler)
{
    FILE       *fp;
    double      Delta, f, y, DD;
#define NLEG 3
    const char *leg[NLEG] = { "f", "fy", "y" };

    DD = pow(10.0, 0.125);
    fp = xvgropen("fy.xvg", "Fig. 2, Lin2003a", "Delta", "y or fy", oenv);
    xvgr_legend(fp, NLEG, leg, oenv);
#undef NLEG
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

static void dump_w0(gmx_output_env_t *oenv,
                    real              beta,
                    const char       *title,
                    LevelOfTheory     lot)
{
    FILE       *fp;
    double      nu;
#define NLEG 4
    const char *leg[NLEG]    = { "wCv", "wS", "wA", "wE" };
    char        buf[128];

    snprintf(buf, 128, "w%s.xvg", lotName(lot));
    fp = xvgropen(buf, title, "\\f{12}b\\f{4}h\\f{12}n",
                  "w", oenv);
    xvgr_legend(fp, NLEG, leg, oenv);
#undef NLEG
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

void dump_w(gmx_output_env_t *oenv, real beta)
{
    dump_w0(oenv, beta, "Fig. 1, Berens1983a", LOT_QUANTUM);
    dump_w0(oenv, beta, "Classical weighting", LOT_CLASSICAL);
}

#define NDIM 4
static void SWAPPER(double dd[DIM], double **ev, int i)
{
    double temp;
    double tvec[NDIM];

    if (fabs(dd[i+1]) > fabs(dd[i]))
    {
        temp = dd[i];
        for (int j = 0; (j < NDIM); j++)
        {
            tvec[j] = ev[j][i];
        }
        dd[i] = dd[i+1];
        for (int j = 0; (j < NDIM); j++)
        {
            ev[j][i] = ev[j][i+1];
        }
        dd[i+1] = temp;
        for (int j = 0; (j < NDIM); j++)
        {
            ev[j][i+1] = tvec[j];
        }
    }
}

void principal(int n, int index[], t_atom atom[], rvec x[],
               matrix trans, rvec d, double **inten, double **ev)
{
    int      i, ai, m, nrot;
    real     mm, rx, ry, rz;
    double   dd[NDIM];

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
        inten[0][0] += mm*(gmx::square(ry)+gmx::square(rz));
        inten[1][1] += mm*(gmx::square(rx)+gmx::square(rz));
        inten[2][2] += mm*(gmx::square(rx)+gmx::square(ry));
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

void get_energy_terms2(const char *enx_fn, int net, t_energy_term et[])
{
    int          i, nre, fff;
    gmx_enxnm_t *enm = NULL;
    gmx_bool     bCont;
    ener_file_t  fp = open_enx(enx_fn, "r");
    t_enxframe   frame;
    int          timecheck = 0;
    gmx_bool     b0set     = FALSE;
    gmx_int64_t  step0;

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
        if (b0set && (frame.step > step0))
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
    for (i = 0; (i < net); i++)
    {
        if (et[i].nesum > 0)
        {
            et[i].energy = et[i].esum/et[i].nesum;
            et[i].stddev = sqrt(et[i].e2sum/et[i].nesum - gmx::square(et[i].energy));
        }
    }
}
