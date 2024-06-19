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
/*! \internal \file
 * \brief
 * Implements function to compute many autocorrelation functions
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "autocorr.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/correlationfunctions/manyautocorrelation.h"
#include "gromacs/correlationfunctions/polynomials.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

/*! \brief Shortcut macro to select modes. */
#define MODE(x) ((mode & (x)) == (x))

typedef struct
{
    unsigned long mode;
    int           nrestart, nout, P, fitfn;
    gmx_bool      bFour, bNormalize;
    real          tbeginfit, tendfit;
} t_acf;

/*! \brief Global variable set true if initialization routines are called. */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx_bool bACFinit = FALSE;

/*! \brief Data structure for storing command line variables. */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static t_acf acf;

enum
{
    enNorm,
    enCos,
    enSin
};

/*! \brief Routine to compute ACF using FFT. */
static void low_do_four_core(int nframes, real c1[], real cfour[], int nCos)
{
    int                            i = 0;
    std::vector<std::vector<real>> data;
    data.resize(1);
    data[0].resize(nframes, 0);
    switch (nCos)
    {
        case enNorm:
            for (i = 0; (i < nframes); i++)
            {
                data[0][i] = c1[i];
            }
            break;
        case enCos:
            for (i = 0; (i < nframes); i++)
            {
                data[0][i] = std::cos(c1[i]);
            }
            break;
        case enSin:
            for (i = 0; (i < nframes); i++)
            {
                data[0][i] = std::sin(c1[i]);
            }
            break;
        default: gmx_fatal(FARGS, "nCos = %d, %s %d", nCos, __FILE__, __LINE__);
    }

    many_auto_correl(&data);
    for (i = 0; (i < nframes); i++)
    {
        cfour[i] = data[0][i];
    }
}

/*! \brief Routine to comput ACF without FFT. */
static void do_ac_core(int nframes, int nout, real corr[], real c1[], int nrestart, unsigned long mode)
{
    int  j, k, j3, jk3, m, n;
    real ccc, cth;
    rvec xj, xk;

    if (nrestart < 1)
    {
        printf("WARNING: setting number of restarts to 1\n");
        nrestart = 1;
    }
    if (debug)
    {
        fprintf(debug, "Starting do_ac_core: nframes=%d, nout=%d, nrestart=%d,mode=%lu\n", nframes, nout, nrestart, mode);
    }

    for (j = 0; (j < nout); j++)
    {
        corr[j] = 0;
    }

    /* Loop over starting points. */
    for (j = 0; (j < nframes); j += nrestart)
    {
        j3 = DIM * j;

        /* Loop over the correlation length for this starting point */
        for (k = 0; (k < nout) && (j + k < nframes); k++)
        {
            jk3 = DIM * (j + k);

            /* Switch over possible ACF types.
             * It might be more efficient to put the loops inside the switch,
             * but this is more clear, and save development time!
             */
            if (MODE(eacNormal))
            {
                corr[k] += c1[j] * c1[j + k];
            }
            else if (MODE(eacCos))
            {
                /* Compute the cos (phi(t)-phi(t+dt)) */
                corr[k] += std::cos(c1[j] - c1[j + k]);
            }
            else if (MODE(eacIden))
            {
                /* Check equality (phi(t)==phi(t+dt)) */
                corr[k] += (c1[j] == c1[j + k]) ? 1 : 0;
            }
            else if (MODE(eacP1) || MODE(eacP2) || MODE(eacP3))
            {
                unsigned int mmm;

                for (m = 0; (m < DIM); m++)
                {
                    xj[m] = c1[j3 + m];
                    xk[m] = c1[jk3 + m];
                }
                cth = cos_angle(xj, xk);

                if (cth - 1.0 > 1.0e-15)
                {
                    printf("j: %d, k: %d, xj:(%g,%g,%g), xk:(%g,%g,%g)\n",
                           j,
                           k,
                           xj[XX],
                           xj[YY],
                           xj[ZZ],
                           xk[XX],
                           xk[YY],
                           xk[ZZ]);
                }
                mmm = 1;
                if (MODE(eacP2))
                {
                    mmm = 2;
                }
                else if (MODE(eacP3))
                {
                    mmm = 3;
                }
                corr[k] += LegendreP(cth, mmm); /* 1.5*cth*cth-0.5; */
            }
            else if (MODE(eacRcross))
            {
                rvec xj, xk, rr;
                for (m = 0; (m < DIM); m++)
                {
                    xj[m] = c1[j3 + m];
                    xk[m] = c1[jk3 + m];
                }
                cprod(xj, xk, rr);

                corr[k] += iprod(rr, rr);
            }
            else if (MODE(eacVector))
            {
                for (m = 0; (m < DIM); m++)
                {
                    xj[m] = c1[j3 + m];
                    xk[m] = c1[jk3 + m];
                }
                ccc = iprod(xj, xk);

                corr[k] += ccc;
            }
            else
            {
                gmx_fatal(FARGS, "\nInvalid mode (%lu) in do_ac_core", mode);
            }
        }
    }
    /* Correct for the number of points and copy results to the data array */
    for (j = 0; (j < nout); j++)
    {
        n     = gmx::divideRoundUp(nframes - j, nrestart);
        c1[j] = corr[j] / n;
    }
}

/*! \brief Routine to normalize ACF, dividing by corr[0]. */
static void normalize_acf(int nout, real corr[])
{
    int    j;
    double c0;

    if (debug)
    {
        fprintf(debug, "Before normalization\n");
        for (j = 0; (j < nout); j++)
        {
            fprintf(debug, "%5d  %10f\n", j, corr[j]);
        }
    }

    /* Normalisation makes that c[0] = 1.0 and that other points are scaled
     * accordingly.
     */
    if (std::fabs(corr[0]) < 1e-5)
    {
        c0 = 1.0;
    }
    else
    {
        c0 = 1.0 / corr[0];
    }
    for (j = 0; (j < nout); j++)
    {
        corr[j] *= c0;
    }

    if (debug)
    {
        fprintf(debug, "After normalization\n");
        for (j = 0; (j < nout); j++)
        {
            fprintf(debug, "%5d  %10f\n", j, corr[j]);
        }
    }
}

/*! \brief Routine that averages ACFs. */
static void average_acf(gmx_bool bVerbose, int n, int nitem, real** c1)
{
    real c0;
    int  i, j;

    if (bVerbose)
    {
        printf("Averaging correlation functions\n");
    }

    for (j = 0; (j < n); j++)
    {
        c0 = 0;
        for (i = 0; (i < nitem); i++)
        {
            c0 += c1[i][j];
        }
        c1[0][j] = c0 / nitem;
    }
}

/*! \brief Normalize ACFs. */
static void norm_and_scale_vectors(int nframes, real c1[], real scale)
{
    int   j, m;
    real* rij;

    for (j = 0; (j < nframes); j++)
    {
        rij = &(c1[j * DIM]);
        unitv(rij, rij);
        for (m = 0; (m < DIM); m++)
        {
            rij[m] *= scale;
        }
    }
}

/*! \brief Debugging */
static void dump_tmp(char* s, int n, real c[])
{
    FILE* fp;
    int   i;

    fp = gmx_ffopen(s, "w");
    for (i = 0; (i < n); i++)
    {
        fprintf(fp, "%10d  %10g\n", i, c[i]);
    }
    gmx_ffclose(fp);
}

/*! \brief High level ACF routine. */
static void do_four_core(unsigned long mode, int nframes, real c1[], real csum[], real ctmp[])
{
    real* cfour;
    char  buf[32];
    real  fac;
    int   j, m, m1;

    snew(cfour, nframes);

    if (MODE(eacNormal))
    {
        /********************************************
         *  N O R M A L
         ********************************************/
        low_do_four_core(nframes, c1, csum, enNorm);
    }
    else if (MODE(eacCos))
    {
        /***************************************************
         * C O S I N E
         ***************************************************/
        /* Copy the data to temp array. Since we need it twice
         * we can't overwrite original.
         */
        for (j = 0; (j < nframes); j++)
        {
            ctmp[j] = c1[j];
        }

        /* Cosine term of AC function */
        low_do_four_core(nframes, ctmp, cfour, enCos);
        for (j = 0; (j < nframes); j++)
        {
            c1[j] = cfour[j];
        }

        /* Sine term of AC function */
        low_do_four_core(nframes, ctmp, cfour, enSin);
        for (j = 0; (j < nframes); j++)
        {
            c1[j] += cfour[j];
            csum[j] = c1[j];
        }
    }
    else if (MODE(eacP2))
    {
        /***************************************************
         * Legendre polynomials
         ***************************************************/
        /* First normalize the vectors */
        norm_and_scale_vectors(nframes, c1, 1.0);

        /* For P2 thingies we have to do six FFT based correls
         * First for XX^2, then for YY^2, then for ZZ^2
         * Then we have to do XY, YZ and XZ (counting these twice)
         * After that we sum them and normalise
         * P2(x) = (3 * cos^2 (x) - 1)/2
         * for unit vectors u and v we compute the cosine as the inner product
         * cos(u,v) = uX vX + uY vY + uZ vZ
         *
         *        oo
         *        /
         * C(t) = |  (3 cos^2(u(t'),u(t'+t)) - 1)/2 dt'
         *        /
         *        0
         *
         * For ACF we need:
         * P2(u(0),u(t)) = [3 * (uX(0) uX(t) +
         *                       uY(0) uY(t) +
         *                       uZ(0) uZ(t))^2 - 1]/2
         *               = [3 * ((uX(0) uX(t))^2 +
         *                       (uY(0) uY(t))^2 +
         *                       (uZ(0) uZ(t))^2 +
         *                 2(uX(0) uY(0) uX(t) uY(t)) +
         *                 2(uX(0) uZ(0) uX(t) uZ(t)) +
         *                 2(uY(0) uZ(0) uY(t) uZ(t))) - 1]/2
         *
         *               = [(3/2) * (<uX^2> + <uY^2> + <uZ^2> +
         *                         2<uXuY> + 2<uXuZ> + 2<uYuZ>) - 0.5]
         *
         */

        /* Because of normalization the number of -0.5 to subtract
         * depends on the number of data points!
         */
        for (j = 0; (j < nframes); j++)
        {
            csum[j] = -0.5 * (nframes - j);
        }

        /***** DIAGONAL ELEMENTS ************/
        for (m = 0; (m < DIM); m++)
        {
            /* Copy the vector data in a linear array */
            for (j = 0; (j < nframes); j++)
            {
                ctmp[j] = gmx::square(c1[DIM * j + m]);
            }
            if (debug)
            {
                sprintf(buf, "c1diag%d.xvg", m);
                dump_tmp(buf, nframes, ctmp);
            }

            low_do_four_core(nframes, ctmp, cfour, enNorm);

            if (debug)
            {
                sprintf(buf, "c1dfout%d.xvg", m);
                dump_tmp(buf, nframes, cfour);
            }
            fac = 1.5;
            for (j = 0; (j < nframes); j++)
            {
                csum[j] += fac * (cfour[j]);
            }
        }
        /******* OFF-DIAGONAL ELEMENTS **********/
        for (m = 0; (m < DIM); m++)
        {
            /* Copy the vector data in a linear array */
            m1 = (m + 1) % DIM;
            for (j = 0; (j < nframes); j++)
            {
                ctmp[j] = c1[DIM * j + m] * c1[DIM * j + m1];
            }

            if (debug)
            {
                sprintf(buf, "c1off%d.xvg", m);
                dump_tmp(buf, nframes, ctmp);
            }
            low_do_four_core(nframes, ctmp, cfour, enNorm);
            if (debug)
            {
                sprintf(buf, "c1ofout%d.xvg", m);
                dump_tmp(buf, nframes, cfour);
            }
            fac = 3.0;
            for (j = 0; (j < nframes); j++)
            {
                csum[j] += fac * cfour[j];
            }
        }
    }
    else if (MODE(eacP1) || MODE(eacVector))
    {
        /***************************************************
         * V E C T O R & P1
         ***************************************************/
        if (MODE(eacP1))
        {
            /* First normalize the vectors */
            norm_and_scale_vectors(nframes, c1, 1.0);
        }

        /* For vector thingies we have to do three FFT based correls
         * First for XX, then for YY, then for ZZ
         * After that we sum them and normalise
         */
        for (j = 0; (j < nframes); j++)
        {
            csum[j] = 0.0;
        }
        for (m = 0; (m < DIM); m++)
        {
            /* Copy the vector data in a linear array */
            for (j = 0; (j < nframes); j++)
            {
                ctmp[j] = c1[DIM * j + m];
            }
            low_do_four_core(nframes, ctmp, cfour, enNorm);
            for (j = 0; (j < nframes); j++)
            {
                csum[j] += cfour[j];
            }
        }
    }
    else
    {
        gmx_fatal(FARGS, "\nUnknown mode in do_autocorr (%lu)", mode);
    }

    sfree(cfour);
    for (j = 0; (j < nframes); j++)
    {
        c1[j] = csum[j] / static_cast<real>(nframes - j);
    }
}

void low_do_autocorr(const char*             fn,
                     const gmx_output_env_t* oenv,
                     const char*             title,
                     int                     nframes,
                     int                     nitem,
                     int                     nout,
                     real**                  c1,
                     real                    dt,
                     unsigned long           mode,
                     int                     nrestart,
                     gmx_bool                bAver,
                     gmx_bool                bNormalize,
                     gmx_bool                bVerbose,
                     real                    tbeginfit,
                     real                    tendfit,
                     int                     eFitFn)
{
    FILE *   fp, *gp = nullptr;
    int      i;
    real*    csum;
    real *   ctmp, *fit;
    real     sum, Ct2av, Ctav;
    gmx_bool bFour = acf.bFour;

    /* Check flags and parameters */
    nout = get_acfnout();
    if (nout == -1)
    {
        nout = acf.nout = (nframes + 1) / 2;
    }
    else if (nout > nframes)
    {
        nout = nframes;
    }

    if (MODE(eacCos) && MODE(eacVector))
    {
        gmx_fatal(FARGS, "Incompatible options bCos && bVector (%s, %d)", __FILE__, __LINE__);
    }
    if ((MODE(eacP3) || MODE(eacRcross)) && bFour)
    {
        if (bVerbose)
        {
            fprintf(stderr, "Can't combine mode %lu with FFT, turning off FFT\n", mode);
        }
        bFour = FALSE;
    }
    if (MODE(eacNormal) && MODE(eacVector))
    {
        gmx_fatal(FARGS, "Incompatible mode bits: normal and vector (or Legendre)");
    }

    /* Print flags and parameters */
    if (bVerbose)
    {
        printf("Will calculate %s of %d thingies for %d frames\n", title ? title : "autocorrelation", nitem, nframes);
        printf("bAver = %s, bFour = %s bNormalize= %s\n",
               gmx::boolToString(bAver),
               gmx::boolToString(bFour),
               gmx::boolToString(bNormalize));
        printf("mode = %lu, dt = %g, nrestart = %d\n", mode, dt, nrestart);
    }
    /* Allocate temp arrays */
    snew(csum, nframes);
    snew(ctmp, nframes);

    /* Loop over items (e.g. molecules or dihedrals)
     * In this loop the actual correlation functions are computed, but without
     * normalizing them.
     */
    for (int i = 0; i < nitem; i++)
    {
        if (bVerbose && (((i % 100) == 0) || (i == nitem - 1)))
        {
            fprintf(stderr, "\rThingie %d", i + 1);
            fflush(stderr);
        }

        if (bFour)
        {
            do_four_core(mode, nframes, c1[i], csum, ctmp);
        }
        else
        {
            do_ac_core(nframes, nout, ctmp, c1[i], nrestart, mode);
        }
    }
    if (bVerbose)
    {
        fprintf(stderr, "\n");
    }
    sfree(ctmp);
    sfree(csum);

    if (fn)
    {
        snew(fit, nout);
        fp = xvgropen(fn, title, "Time (ps)", "C(t)", oenv);
    }
    else
    {
        fit = nullptr;
        fp  = nullptr;
    }
    if (bAver)
    {
        if (nitem > 1)
        {
            average_acf(bVerbose, nframes, nitem, c1);
        }

        if (bNormalize)
        {
            normalize_acf(nout, c1[0]);
        }

        if (eFitFn != effnNONE)
        {
            fit_acf(nout, eFitFn, oenv, fn != nullptr, tbeginfit, tendfit, dt, c1[0], fit);
            sum = print_and_integrate(fp, nout, dt, c1[0], fit, 1);
        }
        else
        {
            sum = print_and_integrate(fp, nout, dt, c1[0], nullptr, 1);
        }
        if (bVerbose)
        {
            printf("Correlation time (integral over corrfn): %g (ps)\n", sum);
        }
    }
    else
    {
        /* Not averaging. Normalize individual ACFs */
        Ctav = Ct2av = 0;
        if (debug)
        {
            gp = xvgropen("ct-distr.xvg", "Correlation times", "item", "time (ps)", oenv);
        }
        for (i = 0; i < nitem; i++)
        {
            if (bNormalize)
            {
                normalize_acf(nout, c1[i]);
            }
            if (eFitFn != effnNONE)
            {
                fit_acf(nout, eFitFn, oenv, fn != nullptr, tbeginfit, tendfit, dt, c1[i], fit);
                sum = print_and_integrate(fp, nout, dt, c1[i], fit, 1);
            }
            else
            {
                sum = print_and_integrate(fp, nout, dt, c1[i], nullptr, 1);
                if (debug)
                {
                    fprintf(debug, "CORRelation time (integral over corrfn %d): %g (ps)\n", i, sum);
                }
            }
            Ctav += sum;
            Ct2av += sum * sum;
            if (debug)
            {
                fprintf(gp, "%5d  %.3f\n", i, sum);
            }
        }
        if (debug)
        {
            xvgrclose(gp);
        }
        if (nitem > 1)
        {
            Ctav /= nitem;
            Ct2av /= nitem;
            printf("Average correlation time %.3f Std. Dev. %.3f Error %.3f (ps)\n",
                   Ctav,
                   std::sqrt((Ct2av - gmx::square(Ctav))),
                   std::sqrt((Ct2av - gmx::square(Ctav)) / (nitem - 1)));
        }
    }
    if (fp)
    {
        xvgrclose(fp);
    }
    sfree(fit);
}

/*! \brief Legend for selecting Legendre polynomials. */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static const char* Leg[] = { nullptr, "0", "1", "2", "3", nullptr };

t_pargs* add_acf_pargs(int* npargs, t_pargs* pa)
{
    t_pargs acfpa[] = {
        { "-acflen",
          FALSE,
          etINT,
          { &acf.nout },
          "Length of the ACF, default is half the number of frames" },
        { "-normalize", FALSE, etBOOL, { &acf.bNormalize }, "Normalize ACF" },
        { "-fftcorr",
          FALSE,
          etBOOL,
          { &acf.bFour },
          "HIDDENUse fast fourier transform for correlation function" },
        { "-nrestart",
          FALSE,
          etINT,
          { &acf.nrestart },
          "HIDDENNumber of frames between time origins for ACF when no FFT is used" },
        { "-P", FALSE, etENUM, { Leg }, "Order of Legendre polynomial for ACF (0 indicates none)" },
        { "-fitfn", FALSE, etENUM, { s_ffn }, "Fit function" },
        { "-beginfit",
          FALSE,
          etREAL,
          { &acf.tbeginfit },
          "Time where to begin the exponential fit of the correlation function" },
        { "-endfit",
          FALSE,
          etREAL,
          { &acf.tendfit },
          "Time where to end the exponential fit of the correlation function, -1 is until the "
          "end" },
    };
    t_pargs* ppa;
    int      i, npa;

    npa = asize(acfpa);
    snew(ppa, *npargs + npa);
    for (i = 0; (i < *npargs); i++)
    {
        ppa[i] = pa[i];
    }
    for (i = 0; (i < npa); i++)
    {
        ppa[*npargs + i] = acfpa[i];
    }
    (*npargs) += npa;

    acf.mode       = 0;
    acf.nrestart   = 1;
    acf.nout       = -1;
    acf.P          = 0;
    acf.fitfn      = effnEXP1;
    acf.bFour      = TRUE;
    acf.bNormalize = TRUE;
    acf.tbeginfit  = 0.0;
    acf.tendfit    = -1;

    bACFinit = TRUE;

    return ppa;
}

void do_autocorr(const char*             fn,
                 const gmx_output_env_t* oenv,
                 const char*             title,
                 int                     nframes,
                 int                     nitem,
                 real**                  c1,
                 real                    dt,
                 unsigned long           mode,
                 gmx_bool                bAver)
{
    if (!bACFinit)
    {
        printf("ACF data structures have not been initialised. Call add_acf_pargs\n");
    }

    /* Handle enumerated types */
    sscanf(Leg[0], "%d", &acf.P);
    acf.fitfn = sffn2effn(s_ffn);

    switch (acf.P)
    {
        case 1: mode = mode | eacP1; break;
        case 2: mode = mode | eacP2; break;
        case 3: mode = mode | eacP3; break;
        default: break;
    }

    low_do_autocorr(fn,
                    oenv,
                    title,
                    nframes,
                    nitem,
                    acf.nout,
                    c1,
                    dt,
                    mode,
                    acf.nrestart,
                    bAver,
                    acf.bNormalize,
                    bDebugMode(),
                    acf.tbeginfit,
                    acf.tendfit,
                    acf.fitfn);
}

int get_acfnout()
{
    if (!bACFinit)
    {
        gmx_fatal(FARGS, "ACF data not initialized yet");
    }

    return acf.nout;
}

int get_acffitfn()
{
    if (!bACFinit)
    {
        gmx_fatal(FARGS, "ACF data not initialized yet");
    }

    return sffn2effn(s_ffn);
}
