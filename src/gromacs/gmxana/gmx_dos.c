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

#include "gromacs/fileio/confio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "gromacs/fileio/futil.h"
#include "gstat.h"
#include "macros.h"
#include "gromacs/math/utilities.h"
#include "physics.h"
#include "index.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/commandline/pargs.h"
#include <string.h>
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "xvgr.h"
#include "correl.h"
#include "gmx_ana.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fileio/trxio.h"

enum {
    VACF, MVACF, DOS, DOS_SOLID, DOS_DIFF, DOS_CP, DOS_S, DOS_A, DOS_E, DOS_NR
};

static double FD(double Delta, double f)
{
    return (2*pow(Delta, -4.5)*pow(f, 7.5) -
            6*pow(Delta, -3)*pow(f, 5) -
            pow(Delta, -1.5)*pow(f, 3.5) +
            6*pow(Delta, -1.5)*pow(f, 2.5) +
            2*f - 2);
}

static double YYY(double f, double y)
{
    return (2*pow(y*f, 3) - sqr(f)*y*(1+6*y) +
            (2+6*y)*f - 2);
}

static double calc_compress(double y)
{
    if (y == 1)
    {
        return 0;
    }
    return ((1+y+sqr(y)-pow(y, 3))/(pow(1-y, 3)));
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
        fprintf(stderr, "Unrealistic tolerance %g for bisector. Setting it to %g\n", tol, tolmin);
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

static double calc_fluidicity(double Delta, double tol)
{
    return bisector(Delta, tol, 0, 1, FD);
}

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

static double calc_Shs(double f, double y)
{
    double fy  = f*y;

    return BOLTZ*(log(calc_compress(fy)) + fy*(3*fy-4)/sqr(1-fy));
}

static real wCsolid(real nu, real beta)
{
    real bhn = beta*PLANCK*nu;
    real ebn, koko;

    if (bhn == 0)
    {
        return 1.0;
    }
    else
    {
        ebn  = exp(bhn);
        koko = sqr(1-ebn);
        return sqr(bhn)*ebn/koko;
    }
}

static real wSsolid(real nu, real beta)
{
    real bhn = beta*PLANCK*nu;

    if (bhn == 0)
    {
        return 1;
    }
    else
    {
        return bhn/(exp(bhn)-1) - log(1-exp(-bhn));
    }
}

static real wAsolid(real nu, real beta)
{
    real bhn = beta*PLANCK*nu;

    if (bhn == 0)
    {
        return 0;
    }
    else
    {
        return log((1-exp(-bhn))/(exp(-bhn/2))) - log(bhn);
    }
}

static real wEsolid(real nu, real beta)
{
    real bhn = beta*PLANCK*nu;

    if (bhn == 0)
    {
        return 1;
    }
    else
    {
        return bhn/2 + bhn/(exp(bhn)-1)-1;
    }
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

static void dump_w(output_env_t oenv, real beta)
{
    FILE       *fp;
    double      nu;
    const char *leg[] = { "wCv", "wS", "wA", "wE" };

    fp = xvgropen("w.xvg", "Fig. 1, Berens1983a", "\\f{12}b\\f{4}h\\f{12}n",
                  "w", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);
    for (nu = 1; (nu < 100); nu += 0.05)
    {
        fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n", beta*PLANCK*nu,
                wCsolid(nu, beta), wSsolid(nu, beta),
                wAsolid(nu, beta), wEsolid(nu, beta));
    }
    xvgrclose(fp);
}

int gmx_dos(int argc, char *argv[])
{
    const char         *desc[] = {
        "[THISMODULE] computes the Density of States from a simulations.",
        "In order for this to be meaningful the velocities must be saved",
        "in the trajecotry with sufficiently high frequency such as to cover",
        "all vibrations. For flexible systems that would be around a few fs",
        "between saving. Properties based on the DoS are printed on the",
        "standard output."
    };
    const char         *bugs[] = {
        "This program needs a lot of memory: total usage equals the number of atoms times 3 times number of frames times 4 (or 8 when run in double precision)."
    };
    FILE               *fp, *fplog;
    t_topology          top;
    int                 ePBC = -1;
    t_trxframe          fr;
    matrix              box;
    int                 gnx;
    char                title[256];
    real                t0, t1, m;
    t_trxstatus        *status;
    int                 nV, nframes, n_alloc, i, j, k, l, fftcode, Nmol, Natom;
    double              rho, dt, V2sum, Vsum, V, tmass, dostot, dos2, dosabs;
    real              **c1, **dos, mi, beta, bfac, *nu, *tt, stddev, c1j;
    output_env_t        oenv;
    gmx_fft_t           fft;
    double              cP, S, A, E, DiffCoeff, Delta, f, y, z, sigHS, Shs, Sig, DoS0, recip_fac;
    double              wCdiff, wSdiff, wAdiff, wEdiff;

    static     gmx_bool bVerbose = TRUE, bAbsolute = FALSE, bNormalize = FALSE;
    static     gmx_bool bRecip   = FALSE, bDump = FALSE;
    static     real     Temp     = 298.15, toler = 1e-6;
    t_pargs             pa[]     = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy." },
        { "-recip", FALSE, etBOOL, {&bRecip},
          "Use cm^-1 on X-axis instead of 1/ps for DoS plots." },
        { "-abs", FALSE, etBOOL, {&bAbsolute},
          "Use the absolute value of the Fourier transform of the VACF as the Density of States. Default is to use the real component only" },
        { "-normdos", FALSE, etBOOL, {&bNormalize},
          "Normalize the DoS such that it adds up to 3N. This is a hack that should not be necessary." },
        { "-T", FALSE, etREAL, {&Temp},
          "Temperature in the simulation" },
        { "-toler", FALSE, etREAL, {&toler},
          "[HIDDEN]Tolerance when computing the fluidicity using bisection algorithm" },
        { "-dump", FALSE, etBOOL, {&bDump},
          "[HIDDEN]Dump the y/fy plot corresponding to Fig. 2 inLin2003a and the and the weighting functions corresponding to Fig. 1 in Berens1983a." }
    };

    t_filenm            fnm[] = {
        { efTRN, "-f",    NULL,    ffREAD  },
        { efTPX, "-s",    NULL,    ffREAD  },
        { efNDX, NULL,    NULL,    ffOPTRD },
        { efXVG, "-vacf", "vacf",  ffWRITE },
        { efXVG, "-mvacf", "mvacf", ffWRITE },
        { efXVG, "-dos",  "dos",   ffWRITE },
        { efLOG, "-g",    "dos",   ffWRITE },
    };
#define NFILE asize(fnm)
    int                 npargs;
    t_pargs            *ppa;
    const char         *DoSlegend[] = {
        "DoS(v)", "DoS(v)[Solid]", "DoS(v)[Diff]"
    };

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                           NFILE, fnm, npargs, ppa, asize(desc), desc,
                           asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    beta = 1/(Temp*BOLTZ);
    if (bDump)
    {
        printf("Dumping reference figures. Thanks for your patience.\n");
        dump_fy(oenv, toler);
        dump_w(oenv, beta);
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
    for (i = 0; (i < top.atoms.nr); i++)
    {
        tmass += top.atoms.atom[i].m;
    }

    Natom = top.atoms.nr;
    Nmol  = top.mols.nr;
    gnx   = Natom*DIM;

    /* Correlation stuff */
    snew(c1, gnx);
    for (i = 0; (i < gnx); i++)
    {
        c1[i] = NULL;
    }

    read_first_frame(oenv, &status, ftp2fn(efTRN, NFILE, fnm), &fr, TRX_NEED_V);
    t0 = fr.time;

    n_alloc = 0;
    nframes = 0;
    Vsum    = V2sum = 0;
    nV      = 0;
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
            for (i = 0; i < gnx; i++)
            {
                srenew(c1[i], n_alloc);
            }
        }
        for (i = 0; i < gnx; i += DIM)
        {
            c1[i+XX][nframes] = fr.v[i/DIM][XX];
            c1[i+YY][nframes] = fr.v[i/DIM][YY];
            c1[i+ZZ][nframes] = fr.v[i/DIM][ZZ];
        }

        t1 = fr.time;

        nframes++;
    }
    while (read_next_frame(oenv, status, &fr));

    close_trj(status);

    dt = (t1-t0)/(nframes-1);
    if (nV > 0)
    {
        V = Vsum/nV;
    }
    if (bVerbose)
    {
        printf("Going to do %d fourier transforms of length %d. Hang on.\n",
               gnx, nframes);
    }
    low_do_autocorr(NULL, oenv, NULL, nframes, gnx, nframes, c1, dt, eacNormal, 0, FALSE,
                    FALSE, FALSE, -1, -1, 0);
    snew(dos, DOS_NR);
    for (j = 0; (j < DOS_NR); j++)
    {
        snew(dos[j], nframes+4);
    }

    if (bVerbose)
    {
        printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");
    }
    for (i = 0; (i < gnx); i += DIM)
    {
        mi = top.atoms.atom[i/DIM].m;
        for (j = 0; (j < nframes/2); j++)
        {
            c1j            = (c1[i+XX][j] + c1[i+YY][j] + c1[i+ZZ][j]);
            dos[VACF][j]  += c1j/Natom;
            dos[MVACF][j] += mi*c1j;
        }
    }
    fp = xvgropen(opt2fn("-vacf", NFILE, fnm), "Velocity ACF",
                  "Time (ps)", "C(t)", oenv);
    snew(tt, nframes/2);
    for (j = 0; (j < nframes/2); j++)
    {
        tt[j] = j*dt;
        fprintf(fp, "%10g  %10g\n", tt[j], dos[VACF][j]);
    }
    xvgrclose(fp);
    fp = xvgropen(opt2fn("-mvacf", NFILE, fnm), "Mass-weighted velocity ACF",
                  "Time (ps)", "C(t)", oenv);
    for (j = 0; (j < nframes/2); j++)
    {
        fprintf(fp, "%10g  %10g\n", tt[j], dos[MVACF][j]);
    }
    xvgrclose(fp);

    if ((fftcode = gmx_fft_init_1d_real(&fft, nframes/2,
                                        GMX_FFT_FLAG_NONE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d_real returned %d", fftcode);
    }
    if ((fftcode = gmx_fft_1d_real(fft, GMX_FFT_REAL_TO_COMPLEX,
                                   (void *)dos[MVACF], (void *)dos[DOS])) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
    }

    /* First compute the DoS */
    /* Magic factor of 8 included now. */
    bfac = 8*dt*beta/2;
    dos2 = 0;
    snew(nu, nframes/4);
    for (j = 0; (j < nframes/4); j++)
    {
        nu[j] = 2*j/(t1-t0);
        dos2 += sqr(dos[DOS][2*j]) + sqr(dos[DOS][2*j+1]);
        if (bAbsolute)
        {
            dos[DOS][j] = bfac*sqrt(sqr(dos[DOS][2*j]) + sqr(dos[DOS][2*j+1]));
        }
        else
        {
            dos[DOS][j] = bfac*dos[DOS][2*j];
        }
    }
    /* Normalize it */
    dostot = evaluate_integral(nframes/4, nu, dos[DOS], NULL, nframes/4, &stddev);
    if (bNormalize)
    {
        for (j = 0; (j < nframes/4); j++)
        {
            dos[DOS][j] *= 3*Natom/dostot;
        }
    }

    /* Now analyze it */
    DoS0 = dos[DOS][0];

    /* Note this eqn. is incorrect in Pascal2011a! */
    Delta = ((2*DoS0/(9*Natom))*sqrt(M_PI*BOLTZ*Temp*Natom/tmass)*
             pow((Natom/V), 1.0/3.0)*pow(6/M_PI, 2.0/3.0));
    f     = calc_fluidicity(Delta, toler);
    y     = calc_y(f, Delta, toler);
    z     = calc_compress(y);
    Sig   = BOLTZ*(5.0/2.0+log(2*M_PI*BOLTZ*Temp/(sqr(PLANCK))*V/(f*Natom)));
    Shs   = Sig+calc_Shs(f, y);
    rho   = (tmass*AMU)/(V*NANO*NANO*NANO);
    sigHS = pow(6*y*V/(M_PI*Natom), 1.0/3.0);

    fprintf(fplog, "System = \"%s\"\n", title);
    fprintf(fplog, "Nmol = %d\n", Nmol);
    fprintf(fplog, "Natom = %d\n", Natom);
    fprintf(fplog, "dt = %g ps\n", dt);
    fprintf(fplog, "tmass = %g amu\n", tmass);
    fprintf(fplog, "V = %g nm^3\n", V);
    fprintf(fplog, "rho = %g g/l\n", rho);
    fprintf(fplog, "T = %g K\n", Temp);
    fprintf(fplog, "beta = %g mol/kJ\n", beta);

    fprintf(fplog, "\nDoS parameters\n");
    fprintf(fplog, "Delta = %g\n", Delta);
    fprintf(fplog, "fluidicity = %g\n", f);
    fprintf(fplog, "hard sphere packing fraction = %g\n", y);
    fprintf(fplog, "hard sphere compressibility = %g\n", z);
    fprintf(fplog, "ideal gas entropy = %g\n", Sig);
    fprintf(fplog, "hard sphere entropy = %g\n", Shs);
    fprintf(fplog, "sigma_HS = %g nm\n", sigHS);
    fprintf(fplog, "DoS0 = %g\n", DoS0);
    fprintf(fplog, "Dos2 = %g\n", dos2);
    fprintf(fplog, "DoSTot = %g\n", dostot);

    /* Now compute solid (2) and diffusive (3) components */
    fp = xvgropen(opt2fn("-dos", NFILE, fnm), "Density of states",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                  "\\f{4}S(\\f{12}n\\f{4})", oenv);
    xvgr_legend(fp, asize(DoSlegend), DoSlegend, oenv);
    recip_fac = bRecip ? (1e7/SPEED_OF_LIGHT) : 1.0;
    for (j = 0; (j < nframes/4); j++)
    {
        dos[DOS_DIFF][j]  = DoS0/(1+sqr(DoS0*M_PI*nu[j]/(6*f*Natom)));
        dos[DOS_SOLID][j] = dos[DOS][j]-dos[DOS_DIFF][j];
        fprintf(fp, "%10g  %10g  %10g  %10g\n",
                recip_fac*nu[j],
                dos[DOS][j]/recip_fac,
                dos[DOS_SOLID][j]/recip_fac,
                dos[DOS_DIFF][j]/recip_fac);
    }
    xvgrclose(fp);

    /* Finally analyze the results! */
    wCdiff = 0.5;
    wSdiff = Shs/(3*BOLTZ); /* Is this correct? */
    wEdiff = 0.5;
    wAdiff = wEdiff-wSdiff;
    for (j = 0; (j < nframes/4); j++)
    {
        dos[DOS_CP][j] = (dos[DOS_DIFF][j]*wCdiff +
                          dos[DOS_SOLID][j]*wCsolid(nu[j], beta));
        dos[DOS_S][j]  = (dos[DOS_DIFF][j]*wSdiff +
                          dos[DOS_SOLID][j]*wSsolid(nu[j], beta));
        dos[DOS_A][j]  = (dos[DOS_DIFF][j]*wAdiff +
                          dos[DOS_SOLID][j]*wAsolid(nu[j], beta));
        dos[DOS_E][j]  = (dos[DOS_DIFF][j]*wEdiff +
                          dos[DOS_SOLID][j]*wEsolid(nu[j], beta));
    }
    DiffCoeff = evaluate_integral(nframes/2, tt, dos[VACF], NULL, nframes/2, &stddev);
    DiffCoeff = 1000*DiffCoeff/3.0;
    fprintf(fplog, "Diffusion coefficient from VACF %g 10^-5 cm^2/s\n",
            DiffCoeff);
    fprintf(fplog, "Diffusion coefficient from DoS %g 10^-5 cm^2/s\n",
            1000*DoS0/(12*tmass*beta));

    cP = BOLTZ * evaluate_integral(nframes/4, nu, dos[DOS_CP], NULL,
                                   nframes/4, &stddev);
    fprintf(fplog, "Heat capacity %g J/mol K\n", 1000*cP/Nmol);

    /*
       S  = BOLTZ * evaluate_integral(nframes/4,nu,dos[DOS_S],NULL,
                                   nframes/4,&stddev);
       fprintf(fplog,"Entropy %g J/mol K\n",1000*S/Nmol);
       A  = BOLTZ * evaluate_integral(nframes/4,nu,dos[DOS_A],NULL,
                                   nframes/4,&stddev);
       fprintf(fplog,"Helmholtz energy %g kJ/mol\n",A/Nmol);
       E  = BOLTZ * evaluate_integral(nframes/4,nu,dos[DOS_E],NULL,
                                   nframes/4,&stddev);
       fprintf(fplog,"Internal energy %g kJ/mol\n",E/Nmol);
     */
    fprintf(fplog, "\nArrivederci!\n");
    gmx_fio_fclose(fplog);

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    return 0;
}
