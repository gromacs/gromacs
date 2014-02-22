/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/utilities.h"
#include "gromacs/fft/fft.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "txtdump.h"
#include "gstat.h"
#include "macros.h"
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "txtdump.h"
#include "typedefs.h"
#include "names.h"
#include "vec.h"
#include "xvgr.h"
#include "correl.h"
#include "gmx_ana.h"
#include "nrjac.h"

#define SPEED_OF_LIGHT 2.99792458E5
// m/s
// TYPE OF VELOCITIES: trans, angular, vibrational, rotational, total
enum Velocity {
    V_T, V_A, V_V, V_R, VEL_TOT, V_SUM, VEL_NR
};

static const char *velocityName[VEL_NR] = {
    "Trans", "Angul", "Vibr", "Rotat", "Total", "Sum"
};

static const char *velocityLong[VEL_NR] = {
    "Translation", "Angular motion", "Vibration", "Rotational motion", "Total motion", "Sum of it all"
};

typedef real complex[2];

typedef struct {
    int    npart;
    double cP, A, S, E, Ec, Emd, E0;
    double f, y, fy, ry, z, nrho, Delta, beta, bfac, partMass;
    double hsdf, ttdf;
    double sigHS, Shs, Sig;
    double DiffACF, DiffDoS;
    double dostot, dof, DoS0;
    double wCdiff, wSdiff, wEdiff, wAdiff, wEcdiff;
    double Tdiff, T;
} t_dosprops;

enum Density {
    VACF = VEL_NR, DOS_SOLID, DOS_DIFF, DOS_CP, DOS_S, DOS_A, DOS_E, DOS_EC, DOS_NR
};

/* double precision dsub_xcm */

real dsub_xcm(rvec x[], int gnx, atom_id *index, t_atom atom[], rvec xcm,
              gmx_bool bQ)
{
    int    i, m;
    double m0, tm, dxcm[3], xx[3], a, b, c;

    dxcm[XX] = dxcm[YY] = dxcm[ZZ] = 0.0;

    tm = 0.0;
    for (i = 0; (i < gnx); i++)
    {
        int ii     = index ? index[i] : i;
        xx[XX] = x[ii][XX];
        xx[YY] = x[ii][YY];
        xx[ZZ] = x[ii][ZZ];
        if (bQ)
        {
            m0 = fabs(atom[ii].q);
        }
        else
        {
            m0 = atom[ii].m;
        }
        tm += m0;
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

    return ((real) tm);
}

/* ********************************************************** */

static double FD(double Delta, double f) // eq 34 JCP 119, 11792 (2003)
{
    return (2*pow(Delta, -4.5)*pow(f, 7.5) -
            6*pow(Delta, -3)*pow(f, 5) -
            pow(Delta, -1.5)*pow(f, 3.5) +
            6*pow(Delta, -1.5)*pow(f, 2.5) +
            2*f - 2);
}

//! fluidicity eq. for hard spheres packing, eq 31 JCP 119, 11792 (2003)
static double YYY(double f, double y)
{
    return (2*pow(y*f, 3) - sqr(f)*y*(1+6*y) +
            (2+6*y)*f - 2);
}

//! compressibility of hard spheres z=(1+y+y^3-y^3)/(1-y)^3
static double calc_compress(double y)
{
    return ((1+y+sqr(y)-pow(y, 3))/(pow((1-y), 3)));
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

// calculate fluidicity f
static double calc_fluidicity(double Delta, double tol)
{
    // solution for FD
    return bisector(Delta, tol, 0, 1, FD);
}

// hard sphere packing fraction, eq 32 JCP 119, 11792 (2003)
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
// entropy for hard spheres
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

static real old_calc_fluidicity(real Delta, real tol)
{
    real fd0, fd, fd1, f, f0 = 0, f1 = 1;
    real tolmin = 1e-6;

    /* Since the fluidity is monotonous according to Fig. 2 in Lin2003a,
       J. Chem. Phys. 112 (2003) p. 11792 we can use a bisection method
       to get it. */
    if (tol < tolmin)
    {
        fprintf(stderr, "Unrealistic tolerance %g for calc_fluidity. Setting it to %g\n", tol, tolmin);
        tol = 1e-6;
    }

    do
    {
        fd0 = FD(Delta, f0);
        fd1 = FD(Delta, f1);
        f   = (f0+f1)*0.5;
        fd  = FD(Delta, f);
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

static real wCsolid(real nu, real beta) // weight function heat capacity constant volume for solid, eq 20 (2013)
{
    real bhn = beta*PLANCK*nu;
    real ebn, koko;

    if (bhn == 0)
    {
        return 1.0;    // in the limit bhn -> 0, eq 3.41 -> 1
    }
    if (bhn > 44.3437) // empirical limit
    {
        return 0.0;
    }
    else
    {
        ebn  = exp(bhn);
        koko = sqr(1-ebn);
        return sqr(bhn)*ebn/koko; // eq 3.41 Berens
    }
}

static real wSsolid(real nu, real beta) // weight function entropy for solid, eq 19 (2013)
{
    real bhn = beta*PLANCK*nu;

    if (bhn == 0)
    {
        return 0; // in the limit bhn -> 0, eq 3.43 -> +inf, set to 0 because S(v)=0 for Solid
    }
    else
    {
        return bhn/(exp(bhn)-1) - log(1-exp(-bhn)); // eq 3.43 Berens
    }
}

static real wAsolid(real nu, real beta) // weight function hentalpy for solid, eq 18 (2013)?
{
    real bhn = beta*PLANCK*nu;

    if (bhn == 0)
    {
        return 0; // in the limit bhn -> 0, eq 3.42 -> -inf, set to 0 because S(v)=0 for Solid
    }
    else
    {
        return log((1-exp(-bhn))/(exp(-bhn/2))); // eq 3.42 Berens
    }
}

static real wEcsolid() // weight function energy for solid, eq 17 (2013)
{
    return 1;
}

static real wEsolid(real nu, real beta) // weight function energy for solid, eq 17 (2013)
{
    real bhn = beta*PLANCK*nu;

    if (bhn == 0)
    {
        return 1; // bhn -> 0, eq 3.40 -> 1
    }
    else
    {
        return bhn/2 + bhn/(exp(bhn)-1); // eq 3.40 Berens
    }
}

double search2PT(double K)
{
    double P, fold, fnew, dPdf, tol;
    int    count;

    fold = 0.0;
    fnew = 0.7293*pow(K, 0.5727); /*initial guess*/
    if (fnew > 0.5)
    {
        fnew = 0.5;
    }
    tol   = 1e-10;
    count = 0;

    while (fabs(fnew-fold) > tol && count < 999)
    {
        dPdf = 0.0;
        fold = fnew;
        P    = 2.0*pow(K, -4.5)*pow(fnew, 7.5)-6.0*pow(K, -3.0)*pow(fnew, 5.0)-pow(K, -1.5)*pow(fnew, 3.5)
            + 6.0*pow(K, -1.5)*pow(fnew, 2.5)+2.0*fnew-2;
        dPdf = 15.0*pow(K, -4.5)*pow(fnew, 6.5)-30.0*pow(K, -3.0)*pow(fnew, 4.0)-3.5*pow(K, -1.5)*pow(fnew, 2.5)
            + 15.0*pow(K, -1.5)*pow(fnew, 1.5)+2.0;
        fnew = fold-(P)/dPdf;
        count++;
    }
    return fnew;
}

static void dump_fy(output_env_t oenv, real toler)
{
    FILE       *fp;
    double      Delta, f, y, DD;
    const char *leg[] = { "f", "fy", "y" };

    DD = pow(10.0, 0.125);
    fp = xvgropen("fy.xvg", "Fig. 2, Lin2003a", "Delta", "y or fy", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);
    fprintf(fp, "@    world 1e-05, 0, 1000, 1\n");
    fprintf(fp, "@    xaxes scale Logarithmic\n");
    for (Delta = 1e-5; (Delta <= 1000); Delta *= DD)
    {
        f = calc_fluidicity(Delta, toler);
        y = calc_y(f, Delta, toler);
        fprintf(fp, "%10g  %10g  %10g  %10g\n", Delta, f, f*y, y);
    }
    fclose(fp);
}

static void dump_w(output_env_t oenv, real beta)
{
    FILE       *fp;
    double      nu;
    const char *leg[] = { "wCv", "wS", "wA", "wE" };

    fp = xvgropen("w.xvg", "Fig. 1, Berens1983a", "\\f{12}b\\f{4}h\\f{12}n",
                  "w", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);
    for (nu = 0; (nu < 100); nu += 0.05)
    {
        fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n", beta*PLANCK*nu,
                wCsolid(nu, beta), wSsolid(nu, beta),
                wAsolid(nu, beta), wEsolid(nu, beta));
    }
    fclose(fp);
}

static void principal(int n, atom_id index[], t_atom atom[], rvec x[],
                      matrix trans, rvec d)
{
#define NDIM 4
    int      i, j, ai, m, nrot;
    real     mm, rx, ry, rz;
    double **inten, dd[NDIM], tvec[NDIM], **ev;
#ifdef DEBUG
    real     e[NDIM];
#endif
    real     temp;

    snew(inten, NDIM);
    snew(ev, NDIM);
    for (i = 0; (i < NDIM); i++)
    {
        snew(inten[i], NDIM);
        snew(ev[i], NDIM);
        dd[i] = 0.0;
#ifdef DEBUG
        e[i] = 0.0;
#endif
    }

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
#define SWAPPER(i)          \
    if (fabs(dd[i+1]) > fabs(dd[i])) {    \
        temp = dd[i];         \
        for (j = 0; (j < NDIM); j++) { tvec[j] = ev[j][i]; } \
        dd[i] = dd[i+1];          \
        for (j = 0; (j < NDIM); j++) { ev[j][i] = ev[j][i+1]; }        \
        dd[i+1] = temp;           \
        for (j = 0; (j < NDIM); j++) { ev[j][i+1] = tvec[j]; }         \
    }
    SWAPPER(0)
    SWAPPER(1)
    SWAPPER(0)

    for (i = 0; (i < DIM); i++)
    {
        d[i] = dd[i];
        for (m = 0; (m < DIM); m++)
        {
            trans[i][m] = ev[i][m];
        }
    }

    for (i = 0; (i < NDIM); i++)
    {
        sfree(inten[i]);
        sfree(ev[i]);
    }
    sfree(inten);
    sfree(ev);
}

static void twoPT(double *tmpg,double *tmps,double s0,double sv,double v,
                  int nmol,double fract_f)
{
  if(fract_f==0) 
  {
      *tmps=sv; *tmpg=0; 
  }
  else if(v==0) 
  {
      *tmpg=s0; *tmps=0; 
  }
  else 
  {
    *tmpg=s0/(1.0+pow(M_PI*s0*v/(6.0*nmol*fract_f),2.0));
    //*tmpg=3*fract_f*nmol*(s0/(1.0+pow(PI*s0*v/2,2.0)));
    if(*tmpg>sv) 
    {
        *tmpg=sv;
    }
    *tmps=sv-*tmpg;
  }
}

static double HSDF(real *pwr,int lenth,double fmin,int nmol,double fract_f)
{
  double hsdf,tmpg,tmps;
  int j;
  hsdf=tmpg=tmps=0;
  for(j=0;j<lenth;j++) 
  {
    twoPT(&tmpg,&tmps,pwr[0],pwr[j],fmin*j,nmol,fract_f);
    if(j==0 || j==lenth-1) 
    { 
        hsdf += tmpg*fmin*0.5;
    }
    else 
    { 
        hsdf += tmpg*fmin; 
    }
  }
  return hsdf;
}

/* ---------------------------------------------------------------------- */

static double TTDF(real *pwr,int lenth,double fmin)
{
  double ttdf;
  int j;
  ttdf=0;
  for(j=0;j<lenth;j++) {
    if(j==0 || j==lenth-1) { ttdf += pwr[j]*fmin*0.5;}
    else { ttdf += pwr[j]*fmin; }
  }
  return ttdf;
}
static void print_dp(FILE *fplog, int t, t_dosprops *dp)
{
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
    fprintf(fplog, "Heat capacity cV                = %g J/mol K\n", dp->cP);
    fprintf(fplog, "Entropy S                       = %g J/mol K\n", dp->S);
    fprintf(fplog, "Helmholtz energy A              = %g kJ/mol\n", dp->A);
    fprintf(fplog, "Internal energy E               = %g kJ/mol\n", dp->E);
    if (V_V != t)
    {
        fprintf(fplog, "Diffusion coefficient from dos[VACF] %g 10^-5 cm^2/s\n",
                dp->DiffACF);
        fprintf(fplog, "Diffusion coefficient from DoS0 %g 10^-5 cm^2/s\n",
                dp->DiffDoS);
    }
}

static void many_auto_correl(int nfunc, int nframes, int nfft, real **c)
{
    gmx_fft_t    fft1;
    complex     *in, *out;
    int          i, j, fftcode;
    
    snew(in, nfft);
    snew(out, nfft);
    if ((fftcode = gmx_fft_init_1d(&fft1, nfft,
                                   GMX_FFT_FLAG_CONSERVATIVE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d returned %d", fftcode);
    }
    for(i = 0; (i<nfunc); i++)
    {
        for(j = 0; j < nframes; j++ ) 
        {
            in[j][0] = c[i][j];
            in[j][1] = 0;
        }
        for( ; (j<nfft); j++)
        {
            in[j][0] = in[j][1] = 0;
        }
        if ((fftcode = gmx_fft_1d(fft1, GMX_FFT_BACKWARD,
                                  (void *)in, (void *)out)) != 0)
        {
            gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
        }
        
        for(j = 0; j < nfft; j++) 
        {
            in[j][0] = (sqr(out[j][0]) + sqr(out[j][1]))/nfft;
            in[j][1] = 0;
        }
        for( ; (j<nfft); j++)
        {
            in[j][0] = in[j][1] = 0;
        }
        
        if ((fftcode = gmx_fft_1d(fft1, GMX_FFT_FORWARD,
                                  (void *)in, (void *)out)) != 0)
        {
            gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
        }
        
        for (j = 0; (j < nfft); j++)
        {
            c[i][j] = out[j][0]/nframes;
        }
    }
    sfree(in);
    sfree(out);
    gmx_fft_destroy(fft1);
}

static void print_stuff(int n, real tt[], real **dos,
                        const char *vacf, const char *mvacf, output_env_t oenv)
{
    FILE *fp;
    int   j;
    
    if (NULL != vacf)
    {
        fp = xvgropen(vacf, "Velocity ACF",
                      "Time (ps)", "C(t)", oenv);
        for (j = 0; (j < n); j++)
        {
            fprintf(fp, "%10g %10g\n", tt[j], dos[VACF][j]);
        }
        fclose(fp);
    }
    if (NULL != mvacf)
    {
        fp = xvgropen(mvacf, "Mass-weighted velocity ACF",
                      "Time (ps)", "C(t)", oenv);
        xvgr_legend(fp, asize(velocityName), velocityName, oenv);
        for (j = 0; (j < n); j++)
        {
            fprintf(fp, "%10g %10g %10g %10g %10g %10g\n",
                    tt[j], dos[V_T][j], dos[V_A][j], dos[V_V][j],
                    dos[V_R][j], dos[VEL_TOT][j]);
        }
        fclose(fp);
    }
}

int gmx_dos(int argc, char *argv[])
{
    const char         *desc[] = {
        "[THISMODULE] computes the Density of States from a simulations.",
        "In order for this to be meaningful the velocities must be saved",
        "in the trajecotry with sufficiently high frequency such as to cover",
        "all vibrations. For flexible systems that would be around a few fs",
        "between saving. Properties based on the DoS are printed in the",
        "log file."
    };
    const char         *bugs[] = {
        "This program needs a lot of memory: total usage equals the number of atoms times 3 times number of frames times 4 (or 8 when run in double precision)."
    };
    FILE               *fp, *fplog;
    t_topology          top;
    int                 ePBC = -1;
    t_trxframe          fr;
    matrix              box;
    char                title[256];
    real                t0, t1;
    t_trxstatus        *status;
    int                 nV, nframes, n_alloc, i, j, k, l, fftcode;
    double              rho, V2sum, Vsum, V, tmass, dosabs;
    real              **dos, mi, *nu, *tt, stddev, c1j;
    output_env_t        oenv;
    gmx_fft_t           fft2;
    t_dosprops          dp[VEL_NR];
    double              recip_fac, dt;

    static     gmx_bool bVerbose = TRUE;
    static     gmx_bool bRecip   = FALSE, bDump = FALSE, bQ = FALSE;
    static     real     toler = 1e-6, Emd = 0;

    /* New variables */
    static     gmx_bool axmol = FALSE;
    int                 Natmxmol, tot_coor, mol_id, id, id_end, rdof;
    real             ***vec, crT;
    /* vec: vector for velocities vec[Velocity][ATOM_ID + POS][FRAME] */
    rvec               *r1, *v1, *vecI;

    /* rotational temperatures */
    rvec                vr, I, xcm, vcm, tmp, angmom, pomega, omega, rT;
    matrix              trans;
    real                v1j, tm, vib;
    int                 t, ai, lid;
    atom_id            *index, *Index;
    t_atom             *atom;
    gmx_rmpbc_t         gpbc = NULL;

    /* Set the memory to zero */
    memset(dp, 0, sizeof(dp));

    /* used for line debugging */
    gmx_bool    D = TRUE, F = TRUE;

    t_pargs     pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy." },
        { "-recip", FALSE, etBOOL, {&bRecip},
          "Use cm^-1 on X-axis instead of 1/ps for DoS plots." },
        { "-Emd", FALSE, etREAL, {&Emd},
          "The average total energy in the MD simulation" },
        { "-toler", FALSE, etREAL, {&toler},
          "[HIDDEN]Tolerance when computing the fluidicity using bisection algorithm" },
        { "-dump", FALSE, etBOOL, {&bDump},
          "[HIDDEN]Dump the y/fy plot corresponding to Fig. 2 inLin2003a and the and the weighting functions corresponding to Fig. 1 in Berens1983a." },
    };

    t_filenm    fnm[] = {
        { efTRN, "-f",    NULL,    ffREAD  },
        { efTPX, "-s",    NULL,    ffREAD  },
        { efNDX, NULL,    NULL,    ffOPTRD },
        { efXVG, "-vacf", "molvacf",  ffWRITE },
        { efXVG, "-mvacf", "molmvacf", ffWRITE },
        { efXVG, "-dos",  "moldos",   ffWRITE },
        { efLOG, "-g",    "moldos",   ffWRITE },
    };
#define NFILE asize(fnm)
    int         npargs;
    t_pargs    *ppa;
    const char *DoSlegend[] = {
        "DoS(v)", "DoS(v)[Solid]", "DoS(v)[Diff]"
    };
    double fudge = 1.5;
    int N_tot, nf2, nfreq;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                           NFILE, fnm, npargs, ppa, asize(desc), desc,
                           asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    /* Allocate beta */
    dp[VEL_TOT].beta = 1/(298.15*BOLTZ);
    if (bDump)
    {
        printf("Dumping reference figures. Thanks for your patience.\n");
        dump_fy(oenv, toler);
        dump_w(oenv, dp[VEL_TOT].beta);
        exit(0);
    }

    fplog = gmx_fio_fopen(ftp2fn(efLOG, NFILE, fnm), "w");
    fprintf(fplog, "Doing density of states analysis based on trajectory.\n");
    please_cite(fplog, "Pascal2011a");
    please_cite(fplog, "Caleman2011b");

    read_tps_conf(ftp2fn(efTPX, NFILE, fnm), title, &top, &ePBC, NULL, NULL, box,
                  TRUE);
    V     = det(box);
    tmass = 0;
    // allocate memory for the mass
    snew(atom, top.atoms.nr);
    for (i = 0; (i < top.atoms.nr); i++)
    {
        // we need to know the mass of each atom
        atom[i] = top.atoms.atom[i];
        tmass  += top.atoms.atom[i].m;
    }

    Natmxmol = top.atoms.nr/top.mols.nr;
    fprintf(stdout, "Natom  %d, Nmol %d, NatmXmol %d\n",
            top.atoms.nr, top.mols.nr, Natmxmol);

    /* pbc corrections */
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, top.atoms.nr);

    read_first_frame(oenv, &status, ftp2fn(efTRN, NFILE, fnm),
                     &fr, TRX_NEED_V | TRX_NEED_X);
    t0 = fr.time;

    n_alloc = 0;
    nframes = 0;
    Vsum    = V2sum = 0;
    nV      = 0;

    /* allocation of memory */
    snew(vec, VEL_NR);
    for (j = 0; (j < VEL_NR); j++)
    {
        snew(vec[j], top.atoms.nr*DIM);
    }

    snew(v1, top.atoms.nr);
    snew(r1, top.atoms.nr);
    snew(index, top.atoms.nr);

    for (i = 0; i < top.atoms.nr; i++)
    {
        index[i] = i;
    }

    /* allocation of memory */
    snew(vecI, top.mols.nr);

    // read each frame, save coordinates and velocities,
    // do calculation and save results
    do
    {

        if (fr.bBox)
        {
            V      = det(fr.box);
            V2sum += V*V;
            Vsum  += V;
            nV++;
        }

        if (nframes >= n_alloc)
        {
            n_alloc += 100;
            for (i = 0; i < top.atoms.nr*DIM; i++)
            {
                for (t = 0; (t < VEL_NR); t++)
                {
                    srenew(vec[t][i], n_alloc);
                }
            }
        }

        /* read velocities and coordinated of all atoms:
         * each atoms is pointed by j
         */
        for (j = i = 0; j < top.atoms.nr; j++, i += DIM)
        {
            /* load velocities */
            v1[j][XX] = fr.v[j][XX];
            v1[j][YY] = fr.v[j][YY];
            v1[j][ZZ] = fr.v[j][ZZ];
            // save total velocities
            vec[VEL_TOT][i+XX][nframes] = fr.v[j][XX];
            vec[VEL_TOT][i+YY][nframes] = fr.v[j][YY];
            vec[VEL_TOT][i+ZZ][nframes] = fr.v[j][ZZ];
            /* load atomic positions */
            r1[j][XX] = fr.x[j][XX];
            r1[j][YY] = fr.x[j][YY];
            r1[j][ZZ] = fr.x[j][ZZ];
        }
        /* Remove periodicity */
        gmx_rmpbc_copy(gpbc, top.atoms.nr, box, r1, r1);

        /* dividing atoms in molecules, id */
        for (id = 0; id < top.mols.nr; id++) // Nmol: number of molecules
        {
            lid   = id*DIM;
            Index = &index[id*Natmxmol];
            // find center of mass and shift all position to this new origin: 
            // it must be done for each
            // atom in the molecule
            tm = dsub_xcm(r1, Natmxmol, Index, atom, xcm, bQ); 
            // return total mass of molecule
            /* find the velocity of center of mass and shift all position to
             * this new origin: must be done for each molecule
             */
            tm = dsub_xcm(v1, Natmxmol, Index, atom, vcm, bQ);
            // save translational velocities for the molecule id
            vec[V_T][lid+XX][nframes] = vcm[XX];
            vec[V_T][lid+YY][nframes] = vcm[YY];
            vec[V_T][lid+ZZ][nframes] = vcm[ZZ];
            
            if (Natmxmol > 1.0)
            {
                // compute principal moment of inertia
                // return trans eigenvector tensor and I the eigenvalues
                // (principal moment of inertia)
                principal(Natmxmol, Index, atom, r1, trans, I);
                /* Save the principal moment: it will be used to compute the 
                 * rotational temperatures. */
                rvec_inc(vecI[id], I);
                
                // we run through the j-th atom of id-th molecule and we
                // save velocities for i-th atoms from Index
                // reset the needed variables
                clear_rvec(angmom);
                clear_rvec(pomega);
                clear_rvec(omega);

                // Index = &index[id*Natmxmol];
                for (j = 0; j < Natmxmol; j++)
                {
                    ai = Index[j];
                    // compute the angular momentum angmom
                    clear_rvec(tmp);
                    cprod(r1[ai], v1[ai], tmp);
                    // mass weigth it for each atom in the molecule
                    angmom[XX] += atom[ai].m*tmp[XX];
                    angmom[YY] += atom[ai].m*tmp[YY];
                    angmom[ZZ] += atom[ai].m*tmp[ZZ];
                }
                // compute angular velocity along principle axis, first step
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
                // calculate angular velocities. Here we use the transpose of the trans
                // matrix by swapping the indexing
                for (j = 0; (j < 3); j++)
                {
                    for (k = 0; k < 3; k++)
                    {
                        omega[j] += pomega[k]*trans[j][k];
                    }
                }
                // calculate inertia weighted angular velocities and save
                for (j = 0; (j < DIM); j++)
                {
                    real vvv = 0;
                    for (k = 0; k < 3; k++)
                    {
                        vvv += pomega[k]*trans[j][k]*sqrt(I[k]);
                    }
                    // save angular velocities of molecule id, weighted by sqrt(I)
                    vec[V_A][id*DIM+j][nframes] = vvv;
                }

                // calculate velocity due to rotation vr = w x r
                for (j = 0; j < Natmxmol; j++)
                {
                    ai = Index[j];
                    // calculate velocity due to rotation vr = w x r
                    cprod(omega, r1[ai], vr);
                    // save rotational velocities
                    vec[V_R][ai*DIM+XX][nframes] = vr[XX];
                    vec[V_R][ai*DIM+YY][nframes] = vr[YY];
                    vec[V_R][ai*DIM+ZZ][nframes] = vr[ZZ];
                    // calculate vibrational velocities and save
                    // v1 is the relative velocity versus center of mass
                    vec[V_V][ai*DIM+XX][nframes] = v1[ai][XX] - vr[XX];
                    vec[V_V][ai*DIM+YY][nframes] = v1[ai][YY] - vr[YY];
                    vec[V_V][ai*DIM+ZZ][nframes] = v1[ai][ZZ] - vr[ZZ];
                }
            } //(Natmxmol>1.0){
        }     // for (id=0; id<Nmol; id++) {

        t1 = fr.time;
        nframes++;
    }
    while (read_next_frame(oenv, status, &fr));
    N_tot = fudge*nframes;
    nf2   = 1+N_tot/2;
    nfreq = nf2;
    printf("N_tot = %d\n", N_tot);

    gmx_rmpbc_done(gpbc);

    sfree(r1);   // empty r1, we need memory
    sfree(v1);   // empty v1, we need memory
    close_trj(status);

    /* Normalize moment of inertia */
    real ddd = 1.0/nframes;
    rvec sumVec;
    clear_rvec(sumVec);
    for (id = 0; id < top.mols.nr; id++)
    {
        svmul(ddd, vecI[id], vecI[id]);
        rvec_inc(sumVec, vecI[id]);
    }
    if (NULL != debug)
    {
        pr_rvecs(debug, 0, "Average vecI", vecI, top.mols.nr);
        svmul(1.0/top.mols.nr, sumVec, sumVec);
        pr_rvec(debug, 0, "Overall Average vecI", sumVec, DIM, FALSE);
    }
        
    dt = (t1-t0)/(nframes-1);
    if (nV > 0)
    {
        V = Vsum/nV;
    }

    /* Compute temperatures, degrees of freedom must be changed to the right one */
    dp[VEL_TOT].dof = 3*Natmxmol;
    dp[V_T].dof     = 3;
    dp[V_A].dof     = 3;
    if (Natmxmol == 1)
    {
        dp[V_A].dof = 0;
    }
    else if (Natmxmol == 2)
    {
        dp[V_A].dof = 2;
    }
    dp[V_R].dof = dp[V_A].dof;
    dp[V_V].dof = dp[VEL_TOT].dof - dp[V_T].dof - dp[V_A].dof;

    if (bVerbose)
    {
        printf("Going to do %d fourier transforms of length %d. Hang on.\n",
               top.atoms.nr*DIM*VEL_NR, N_tot);
    }

    // Allocate temp arrays
    snew(dos, DOS_NR);
    for (j = 0; (j < DOS_NR); j++)
    {
        snew_aligned(dos[j], N_tot, 16);
    }
    dp[V_T].npart = dp[V_A].npart = top.mols.nr;
    dp[V_V].npart = dp[V_R].npart = dp[VEL_TOT].npart = top.atoms.nr;

    /* Set the particle mass */
    for (t = 0; (t < VEL_NR); t++)
    {
        dp[t].partMass = tmass/dp[t].npart;
    }
    /* Padding the vec array with zeroes */
    for (t = 0; (t < VEL_NR); t++)
    {
        for(i=0; (i<DIM*dp[t].npart); i++)
        {
            srenew(vec[t][i], N_tot);
            for(j = nframes; (j < N_tot); j++)
            {
                vec[t][i][j] = 0;
            }
        }
    }
    
    if (bVerbose)
    {
        printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");
    }

    for (t = 0; (t < VEL_NR); t++)
    {
        if (dp[t].dof > 0)
        {
            many_auto_correl(DIM*dp[t].npart, nframes, N_tot, vec[t]);
            for (i = 0; (i < dp[t].npart); i++)
            {
                int ai    = i*DIM;
                real mass = atom[i].m;
                if (V_T == t)
                {
                    mass = tm;
                }
                else if (V_A == t)
                {
                    mass = 1;
                }
                
                for (j = 0; (j < N_tot); j++)
                {
                    v1j        = (vec[t][ai+XX][j] + 
                                  vec[t][ai+YY][j] +
                                  vec[t][ai+ZZ][j]);
                    dos[t][j] += v1j*mass;
                    if (VEL_TOT == t)
                    {
                        dos[VACF][j] += v1j/dp[t].npart;
                    }
                }
            }
        }
    }

    /* Compute the temperatures. The number of degrees of freedom
     * are per molecule, therefore we divide the DoS by the number
     * of molecules.
     * From here, the temperature is taken from simulation.
     */
    for (t = 0; (t < VEL_NR); t++)
    {
        if (dp[t].dof > 0)
        {
            dp[t].T    = dos[t][0]/(dp[t].dof*BOLTZ*top.mols.nr);
            dp[t].beta = 1.0/(dp[t].T*BOLTZ);
            dp[t].bfac = 2.0*dt*dp[t].beta;
        }
    }
    /* Make time- and frequency arrays */
    snew(tt, N_tot);
    snew(nu, N_tot);
    double pwrfreq = 1/(dt*(N_tot+1));
    for (j = 0; (j < N_tot); j++)
    {
        tt[j] = j*dt;
        nu[j] = j*pwrfreq;
    }
    print_stuff(N_tot, tt, dos, opt2fn_null("-vacf", NFILE, fnm), 
                opt2fn_null("-mvacf", NFILE, fnm), oenv);
    if ((fftcode = gmx_fft_init_1d(&fft2, N_tot,
                                   GMX_FFT_FLAG_CONSERVATIVE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d returned %d", fftcode);
    }

    //compute density of state for each kind of velocities, return in same vector
    for (t = 0; t < VEL_NR; t++)
    {
        if (dp[t].dof > 0)
        {
            complex *in, *out;
            snew(in, N_tot);
            snew(out, N_tot);
            for (j = 0; (j < N_tot); j++)
            {
                in[j][0] = dp[t].bfac*dos[t][j];
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
                dos[t][j] = sqrt(sqr(out[j][0])+sqr(out[j][1]));
            }
            sfree(in);
            sfree(out);
        }
    }
    gmx_fft_destroy(fft2);
    
    fp = xvgropen("DoS.xvg", "DoS",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)", "DoS(v)", oenv);
    double recip0 = (1e7/SPEED_OF_LIGHT);
    recip_fac = bRecip ? recip0 : 1.0;
    xvgr_legend(fp, asize(velocityName), velocityName, oenv);
    for (j = 0; (j < nfreq); j++)
    {
        fprintf(fp, "%10g %10g %10g %10g %10g %10g %10g\n",
                recip_fac*nu[j], dos[V_T][j]/recip_fac, dos[V_A][j]/recip_fac, dos[V_V][j]/recip_fac, dos[V_R][j]/recip_fac, dos[VEL_TOT][j]/recip_fac, wEsolid(nu[j], dp[V_V].beta));
    }
    fclose(fp);

    /* Density */
    rho   = (tmass*AMU)/(V*NANO*NANO*NANO);

    if (0)
    {
        char **strings;
        int nline = get_file("test.pwr", &strings);
        for(j = 0; (j<nline); j++)
    {
        double nnu, ddos;
        if (2 == sscanf(strings[j],"%lf %lf", &nnu, &ddos)) 
        {
            nnu = nnu/recip0;
            dos[V_V][j] = ddos*recip0;
            if (fabs(nnu - nu[j]) > 0.01) 
            {
                printf("Strange nu %g (expected %g) dos = %g line %d\n", 
                       nnu, nu[j], ddos, j);
            }
        }
        else
        {
            printf("Couldn't read two arguments from line %d in test.pwr\n", j);
        }
        sfree(strings[j]);
    }
    sfree(strings);
    }
    /* Analyze the DoS[t] */
    for (t = 0; (t < VEL_NR); t++)
    {
        if (dp[t].dof == 0)
        {
            continue;
        }
        dp[t].DoS0 = dos[t][0];

        if (V_V != t)
        {
            // eq 13 JPC B, 114, 8191 (2010)
            dp[t].nrho  = dp[t].npart/V;
            dp[t].Delta = ((2*dp[t].DoS0/(9*dp[t].npart))*sqrt(M_PI*BOLTZ*dp[t].T/dp[t].partMass)*
                           pow(dp[t].nrho, 1.0/3.0)*pow(6/M_PI, 2.0/3.0));
            dp[t].f     = search2PT(dp[t].Delta);
            //calc_fluidicity(dp[t].Delta, toler);
            dp[t].ry    = calc_y(dp[t].f, dp[t].Delta, toler);
            
            dp[t].hsdf  = HSDF(dos[t], nfreq, nu[1], top.mols.nr, dp[t].f);
            /* determine hard sphere degrees of freedom hsdf */
            dp[t].ttdf  = TTDF(dos[t], nfreq, nu[1]); 
            /* determine total degrees of freedom ttdf */
            dp[t].y     = dp[t].ry*(dp[t].hsdf/dp[t].ttdf);

            dp[t].fy    = dp[t].f * dp[t].y;
            
            dp[t].z     = calc_compress(dp[t].y);
            // eq 16, Lin 2010, Sackur-Tetrode eq.
            double xlog = (pow(2*M_PI*dp[t].partMass*BOLTZ*dp[t].T/
                               (sqr(PLANCK)), 1.5)*V/(dp[t].f * dp[t].npart));
            dp[t].Sig   = BOLTZ*(2.5+log(xlog));
            fprintf(fplog, "m %g k %g T %g h %g V %g N %d xlog %g\n",
                    dp[t].partMass,BOLTZ,dp[t].T,PLANCK,V,dp[t].npart, xlog);
                    
            dp[t].Shs   = dp[t].Sig+calc_Shs(dp[t].y);
            dp[t].sigHS = pow(6*dp[t].ry*V/(M_PI*dp[t].npart), 1.0/3.0);
        }
        else
        {
            // All other parameters are 0
            dp[t].y     = 1;
        }

        dp[t].dostot = evaluate_integral(nfreq, nu, dos[t], NULL, nfreq, &stddev);

        /* Now compute solid (2) and diffusive (3) components */
        char buf[STRLEN];
        snprintf(buf, STRLEN, "%s-%s", velocityName[t],
                 opt2fn("-dos", NFILE, fnm));
        fp = xvgropen(buf, "Density of states",
                      bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                      "\\f{4}S(\\f{12}n\\f{4})", oenv);
        xvgr_legend(fp, asize(DoSlegend), DoSlegend, oenv);
        for (j = 0; (j < nfreq); j++)
        {
            dos[DOS_DIFF][j] = 0;
            if (V_V != t)
            {
                dos[DOS_DIFF][j]  = dp[t].DoS0/(1+sqr(dp[t].DoS0*M_PI*nu[j]/(6*dp[t].f*dp[t].npart)));
            }
            dos[DOS_SOLID][j] = dos[t][j]-dos[DOS_DIFF][j];

            fprintf(fp, "%10g  %10g  %10g  %10g\n",
                    recip_fac*nu[j],
                    dos[t][j]/recip_fac,
                    dos[DOS_SOLID][j]/recip_fac,
                    dos[DOS_DIFF][j]/recip_fac);
        }
        fclose(fp);

        /* Finally analyze the results! */
        dp[t].wCdiff = dp[t].wEdiff = dp[t].wEcdiff = dp[t].wAdiff = dp[t].wSdiff = 0;
        if (V_A == t)
        {
            clear_rvec(rT);
            clear_rvec(I);

            const double Q = PLANCK*PLANCK/(8.0*M_PI*M_PI*BOLTZ);
            for (i = 0; i < dp[t].npart; i++)   // run for each molecule
            {
                for (k = 0; (k < DIM); k++)
                {
                    if (vecI[i][k] > 0.0)
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

            crT  = 1.0; // characteristic rotational temperature
            rdof = 0;   // rotational degrees of freedom

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
            real rot = 2.0;
            real SR  = (3.0/2.0 + log(sqrt(M_PI*pow(dp[V_A].T, rdof)/crT)/rot))/3;
            dp[t].wSdiff  = SR;
            fprintf(fplog, "SR = %g, rdof = %d, crT = %g rot = %g\n", SR, rdof, crT, rot);
        }
        else
        {
            dp[t].wSdiff = dp[t].Shs/(3*BOLTZ);
        }
        if (V_V != t)
        {
            dp[t].wCdiff  = 0.5;
            dp[t].wEdiff  = 0.5;
            dp[t].wEcdiff = 0.5;
            dp[t].wAdiff  = dp[t].wEdiff - dp[t].wSdiff;
        }
        for (j = 0; (j < nfreq); j++)
        {
            dos[DOS_CP][j]  = (dos[DOS_DIFF][j]*dp[t].wCdiff +
                               dos[DOS_SOLID][j]*wCsolid(nu[j], dp[t].beta));
            dos[DOS_S][j]   = (dos[DOS_DIFF][j]*dp[t].wSdiff +
                               dos[DOS_SOLID][j]*wSsolid(nu[j], dp[t].beta));
            dos[DOS_A][j]   = (dos[DOS_DIFF][j]*dp[t].wAdiff +
                               dos[DOS_SOLID][j]*wAsolid(nu[j], dp[t].beta));
            dos[DOS_E][j]   = (dos[DOS_DIFF][j]*dp[t].wEdiff +
                               dos[DOS_SOLID][j]*wEsolid(nu[j], dp[t].beta));
            dos[DOS_EC][j]  = (dos[DOS_DIFF][j]*dp[t].wEcdiff +
                               dos[DOS_SOLID][j]*wEcsolid());
        }

        if (V_V != t)
        {
            dp[t].DiffACF = (1000/3.0)*evaluate_integral(nfreq, tt, dos[VACF], NULL,
                                                         0, &stddev);
            dp[t].DiffDoS = 1000*dp[t].DoS0/(12*tmass*dp[t].beta);
        }
        double fac = BOLTZ/top.mols.nr;
        dp[t].cP = KILO * fac * evaluate_integral(nfreq, nu, dos[DOS_CP], NULL,
                                                  0, &stddev);
        dp[t].S  = KILO * fac * evaluate_integral(nfreq, nu, dos[DOS_S], NULL,
                                                  0, &stddev);

        dp[t].A  = fac * dp[t].T * evaluate_integral(nfreq, nu, dos[DOS_A], NULL,
                                                     0, &stddev);
        dp[t].E  = fac * dp[t].T * evaluate_integral(nfreq, nu, dos[DOS_E], NULL,
                                                     0, &stddev);
        dp[t].Ec = fac * dp[t].T * evaluate_integral(nfreq, nu, dos[DOS_EC], NULL,
                                                     0, &stddev);
    }
    for (t = 0; (t <= VEL_TOT); t++)
    {    
        double Ec2 = (dp[t].dof/dp[VEL_TOT].dof)*DIM*top.atoms.nr*BOLTZ*dp[t].T*(1 - dp[V_T].wEdiff*dp[V_T].f - dp[V_A].wEdiff*dp[V_A].f);
        fprintf(fplog, "Ec %10g  Ec2 %10g\n", dp[t].Ec, Ec2);
        dp[t].Emd = Emd*dp[t].dof/dp[VEL_TOT].dof;
        dp[t].E0  = dp[t].Emd - dp[t].Ec*top.mols.nr;
        
        print_dp(fplog, t, &dp[t]);
    }
    /* Sum the energies etc. */
    //double E0 = Emd/3 - DIM*BOLTZ*dp[VEL_TOT].T*top.atoms.nr*(1 - dp[V_T].wEdiff*dp[V_T].f - dp[V_A].wEdiff*dp[V_A].f);
    //E0 = Emd - DIM*BOLTZ*dp[VEL_TOT].T*top.atoms.nr;
    
    for (t = 0; (t <= V_V); t++)
    {    
        dp[V_SUM].dostot += dp[t].dostot;
        dp[V_SUM].DoS0 += dp[t].DoS0;
        dp[V_SUM].cP   += dp[t].cP;
        dp[V_SUM].A    += dp[t].A;
        dp[V_SUM].S    += dp[t].S;
        dp[V_SUM].E    += dp[t].E;
        dp[V_SUM].Ec   += dp[t].Ec;
        dp[V_SUM].E0   += dp[t].E0;
        dp[V_SUM].dof  += dp[t].dof;
        dp[V_SUM].T    += dp[t].T*dp[t].dof;
    }
    dp[V_SUM].T /= dp[V_SUM].dof;
    /* ---------------------- */
    fprintf(fplog, "\n\t\t Final resume \n");
    fprintf(fplog, "System  = \"%s\"\n", title);
    fprintf(fplog, "Nmol    = %d\n", top.mols.nr);
    fprintf(fplog, "Natom   = %d\n", top.atoms.nr);
    fprintf(fplog, "dt      = %g ps\n", dt);
    fprintf(fplog, "tmass   = %g amu\n", tmass);
    fprintf(fplog, "Volume  = %g nm^3\n", V);
    fprintf(fplog, "Density = %g g/l\n", rho);
    fprintf(fplog, "Emd     = %g kJ/mol/system\n", Emd);
    fprintf(fplog, "bRecip  = %s\n",bool_names[bRecip]);
    //    fprintf(fplog, "E0      = %g kJ/mol/system\n", E0);

    {
        const char *items[] = {
            "", "Temperature", "DOF", "Dos integral", "DoS0",
            "cV", "A", "S", "E", "Ec", "fluidicity", 
            "wCdiff", "wAdiff", "wSdiff", "wEdiff", "wEcdiff", "E0"
        };

        for (k = 0; (k < asize(items)); k++)
        {
            fprintf(fplog, "%-12s", items[k]);
            for (t = 0; (t <= V_SUM); t++)
            {
                if (dp[t].dof == 0)
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
                    fprintf(fplog, "  %10g", dos[t][0]);
                    break;
                case 5:
                    fprintf(fplog, "  %10g", dp[t].cP*top.mols.nr);
                    break;
                case 6:
                    fprintf(fplog, "  %10g", dp[t].E0+dp[t].A*top.mols.nr);
                    break;
                case 7:
                    fprintf(fplog, "  %10g", dp[t].S*top.mols.nr);
                    break;
                case 8:
                    fprintf(fplog, "  %10g", dp[t].E0+dp[t].E*top.mols.nr);
                    break;
                case 9:
                    fprintf(fplog, "  %10g", dp[t].Ec*top.mols.nr);
                    break;
                case 10:
                    fprintf(fplog, "  %10g", dp[t].f);
                    break;
                case 11:
                    fprintf(fplog, "  %10g", dp[t].wCdiff);
                    break;
                case 12:
                    fprintf(fplog, "  %10g", dp[t].wAdiff);
                    break;
                case 13:
                    fprintf(fplog, "  %10g", dp[t].wSdiff);
                    break;
                case 14:
                    fprintf(fplog, "  %10g", dp[t].wEdiff);
                    break;
                case 15:
                    fprintf(fplog, "  %10g", dp[t].wEcdiff);
                    break;
                case 16:
                    fprintf(fplog, "  %10g", dp[t].E0);
                    break;
                }
            }
            fprintf(fplog, "\n");
        }
    }
    sfree(atom);

    fprintf(fplog, "\nArrivederci!\n");
    fclose(fplog);
    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    return 0;
}
