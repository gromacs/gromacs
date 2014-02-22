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
#include "smalloc.h"
#include "string.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "gromacs/fileio/strdb.h"
#include "xvgr.h"
#include "correl.h"
#include "gmx_ana.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fileio/trxio.h"
#include "nrjac.h"

// TYPE OF VELOCITIES: trans, angular, vibrational, rotational, total
enum Velocity {
    V_T, V_A, V_V, V_R, VEL_TOT, VEL_NR
};
enum Density {
    VACF = VEL_NR, DOS_SOLID, DOS_DIFF, DOS_CP, DOS_S, DOS_A, DOS_E, DOS_NR
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
    if (y == 1)
    {
        return 0;
    }
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

static double calc_fluidicity(double Delta, double tol) // calculate fluidicity f
{
    return bisector(Delta, tol, 0, 1, FD);              // solution for FD
}

static double calc_y(double f, double Delta, double toler) // hard sphere packing fraction, eq 32 JCP 119, 11792 (2003)
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

static double calc_Shs(double f, double y) // entropy for hard spheres
{
    double fy  = f*y;

    return BOLTZ*(log(calc_compress(fy)) + fy*(3*fy-4)/sqr(1-fy));
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
    for (nu = 1; (nu < 100); nu += 0.05)
    {
        fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n", beta*PLANCK*nu,
                wCsolid(nu, beta), wSsolid(nu, beta),
                wAsolid(nu, beta), wEsolid(nu, beta));
    }
    fclose(fp);
}
/*
   enum Velocity { V_T, V_A, V_V, V_R, VEL_TOT, VEL_NR }; // TYPE OF VELOCITIES: trans, angular, vibrational, rotational, total

   enum Density { VACF=VEL_NR, DOS_SOLID, DOS_DIFF, DOS_CP, DOS_S, DOS_A, DOS_E, DOS_NR
   };
 */

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

    for (i = 0; (i < DIM); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            trans[i][m] = inten[i][m];
        }
    }

    /* Call numerical recipe routines */
    jacobi(inten, 3, dd, ev, &nrot);
#ifdef DEBUG
    ptrans("jacobi", ev, dd, e);
#endif

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
#ifdef DEBUG
    ptrans("swap", ev, dd, e);
    t_trans(trans, dd, ev);
#endif

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
    int                 gnx;
    char                title[256];
    real                t0, t1;
    t_trxstatus        *status;
    int                 nV, nframes, n_alloc, i, j, k, l, fftcode, Nmol, Natom, NNN;
    double              rho, dt, V2sum, Vsum, V, tmass, dostot[VEL_NR], *dos2, dosabs, mmm;
    real              **dos, mi, *beta, *bfac, *nu, *tt, stddev, c1j;
    output_env_t        oenv;
    gmx_fft_t           fft;
    double              cP[VEL_NR], S[VEL_NR], A[VEL_NR], E[VEL_NR], f[VEL_NR];
    double              DiffCoeff, Delta, y, z, sigHS, Shs, Sig, DoS0, recip_fac;
    double              wCdiff, wSdiff, wAdiff, wEdiff, V0;

    static     gmx_bool bVerbose = TRUE, bAbsolute = FALSE, bNormalize = FALSE;
    static     gmx_bool bRecip   = FALSE, bDump = FALSE, bQ = FALSE;
    static     real     Temp     = 298.15, toler = 1e-6, Emd = 0;

    /* New variables */
    static     gmx_bool axmol = FALSE;
    int                 Natmxmol, tot_coor, mol_id, id, id_end, rdof;
    real             ***vec, *m, *vecT, **vecI, crT;
    /* vec: vector for velocities vec[Velocity][ATOM_ID + POS][FRAME] */
    /* vecT: vector for temperature from velocities vecT[Velocity] */
    rvec           *r1, *v1;
    rvec            vr, I, xcm, vcm, tmp, angmom, pomega, omega, rT /*rotational temperatures*/;
    matrix          trans;
    real            v1j, tm /*mass of molecule*/, vib, dof[VEL_NR], vax, vay, vaz;
    real            dofT, T_f, T_diff; /* total quantities*/
    int             t, ai, lid;
    atom_id        *index, *Index;
    t_atom         *atom;
    gmx_rmpbc_t     gpbc = NULL;

    // reset variables
    angmom[XX] = 0.0; angmom[YY] = 0.0; angmom[ZZ] = 0.0;
    vax        = 0.0;        vay = 0.0;        vaz = 0.0;
    pomega[XX] = 0.0; pomega[YY] = 0.0; pomega[ZZ] = 0.0;
    omega[XX]  = 0.0;  omega[YY] = 0.0;  omega[ZZ] = 0.0;
    for(t = V_T; (t < VEL_NR); t++)
    {
        cP[t] = S[t] = A[t] = E[t] = 0;
    }
    T_f = 0.0;        T_diff = 0.0;
    // used for line debugging
    gmx_bool   D = TRUE, F = TRUE;


    t_pargs     pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy." },
        { "-recip", FALSE, etBOOL, {&bRecip},
          "Use cm^-1 on X-axis instead of 1/ps for DoS plots." },
        { "-abs", FALSE, etBOOL, {&bAbsolute},
          "[HIDDEN]Use the absolute value of the Fourier transform of the VACF as the Density of States. Default is to use the real component only" },
        { "-normdos", FALSE, etBOOL, {&bNormalize},
          "[HIDDEN]Normalize the DoS such that it adds up to 3N. This is a hack that should not be necessary." },
        { "-T", FALSE, etREAL, {&Temp},
          "Temperature in the simulation" },
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

    const char *Mvacflegend[VEL_NR] = {
        "V_T", "V_A", "V_V", "V_R", "VEL_TOT"
    };
    //    CopyRight(stderr,argv[0]); give error in my compilation
    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE, fnm, npargs, ppa, asize(desc), desc,
                      asize(bugs), bugs, &oenv);

    /* Allocate beta */
    snew(beta, VEL_NR);
    beta[V_T] = 1/(Temp*BOLTZ);
    if (bDump)
    {
        printf("Dumping reference figures. Thanks for your patience.\n");
        dump_fy(oenv, toler);
        dump_w(oenv, beta[0]);
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
    snew(m, top.atoms.nr);
    snew(atom, top.atoms.nr);
    for (i = 0; (i < top.atoms.nr); i++)
    {
        tmass    += top.atoms.atom[i].m;
        m[i]      = top.atoms.atom[i].m; // we need to know the mass of each atom
        atom[i].m = top.atoms.atom[i].m; // we need to know the mass of each atom
    }

    Natom    = top.atoms.nr;
    Nmol     = top.mols.nr;
    Natmxmol = Natom/Nmol;
    fprintf(stdout, "Natom  %d, Nmol %d, NatmXmol %d\n", Natom, Nmol, Natmxmol);
    gnx = Natom*DIM;

    /* pbc corrections */
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, Natom);

    read_first_frame(oenv, &status, ftp2fn(efTRN, NFILE, fnm), &fr, TRX_NEED_V | TRX_NEED_X);
    t0 = fr.time;

    n_alloc = 0;
    nframes = 0;
    Vsum    = V2sum = 0;
    nV      = 0;

    /* allocation of memory */
    snew(vec, VEL_NR);
    for (j = 0; (j < VEL_NR); j++)
    {
        snew(vec[j], gnx);
    }

    snew(vecT, VEL_NR);
    snew(v1, Natom);
    snew(r1, Natom);
    snew(index, Natom);

    for (i = 0; i < Natom; i++)
    {
        index[i] = i;
    }

    /* allocation of memory */
    snew(vecI, Nmol);
    for (j = 0; (j < Nmol); j++)
    {
        snew(vecI[j], DIM);
    }

// read each frame, save coordinates and velocities, do calculation and save results
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
                srenew(vec[V_T][i], n_alloc);
                srenew(vec[V_A][i], n_alloc);
                srenew(vec[V_V][i], n_alloc);
                srenew(vec[V_R][i], n_alloc);
                srenew(vec[VEL_TOT][i], n_alloc);
            }
        }

/* read velocities and coordinated of all atoms: each atoms is pointed by j */
        for (j = 0, i = 0; j < Natom; j++, i += DIM)
        {
            /* load velocities */
            v1[j][XX] = fr.v[j][XX];
            v1[j][YY] = fr.v[j][YY];
            v1[j][ZZ] = fr.v[j][ZZ];
            // save mass weighted total velocities
            vec[VEL_TOT][i+XX][nframes] = fr.v[j][XX];
            vec[VEL_TOT][i+YY][nframes] = fr.v[j][YY];
            vec[VEL_TOT][i+ZZ][nframes] = fr.v[j][ZZ];
            /* load atomic positions */
            r1[j][XX] = fr.x[j][XX];
            r1[j][YY] = fr.x[j][YY];
            r1[j][ZZ] = fr.x[j][ZZ];
        }
        /* Remove periodicity */
        gmx_rmpbc_copy(gpbc, Natom, box, r1, r1);

        /* dividing atoms in molecules, id */
        for (id = 0; id < Nmol; id++) // Nmol: number of molecules
        {
            lid   = id*DIM;
            Index = &index[id*Natmxmol];
#if 0
            if (D && (nframes == 1))
            {
                fprintf(stderr, "Index ");
                for (k = 0; k < Natmxmol; k++)
                {
                    fprintf(stderr, "  %d  %g", Index[k], m[k]);
                }
                fprintf(stderr, "\n");
            }
#endif
            // find center of mass and shift all position to this new origin: it must be done for each
            // atom in the molecule
            // real sub_xcm(rvec x[], int gnx, atom_id *index, t_atom atom[], rvec xcm, gmx_bool bQ)
            tm = dsub_xcm(r1, Natmxmol, Index, atom, xcm, bQ); // return total mass of molecule
#if 0
            if (!D)
            {
                fprintf(stderr, "dxc %g %g %g ", xcm[XX], xcm[YY], xcm[ZZ]);
            }
            if (!D)
            {
                fprintf(stderr, "sub %g %g %g \n", r1[Index[0]][XX], r1[Index[0]][YY], r1[Index[0]][ZZ]);
            }
#endif
            // find the velocity of center of mass and shift all position to this new origin: it must
            //  be done for each molecule
            tm = dsub_xcm(v1, Natmxmol, Index, atom, vcm, bQ);
#if 0
            if (nframes == 1)
            {
                fprintf(stderr, "vcm %g %g %g ", vcm[XX], vcm[YY], vcm[ZZ]);
            }
            if (nframes == 1)
            {
                fprintf(stderr, "sub %g %g %g \n", v1[Index[0]][XX], v1[Index[0]][YY], v1[Index[0]][ZZ]);
            }
#endif
            // save translational velocities for the molecule id
            vec[V_T][lid+XX][nframes] = vcm[XX];
            vec[V_T][lid+YY][nframes] = vcm[YY];
            vec[V_T][lid+ZZ][nframes] = vcm[ZZ];
#if 0
            for (j = 0; j < Natmxmol; j++)
            {
                // save atomic translational velocities
                vec[V_T][Index[j]*3+XX][nframes] = vcm[XX];
                vec[V_T][Index[j]*3+YY][nframes] = vcm[YY];
                vec[V_T][Index[j]*3+ZZ][nframes] = vcm[ZZ];

            }
#endif
            if (Natmxmol > 1.0)
            {
                // compute principal moment of inertia
                // return trans eigenvector and I the eigenvalue (principal moment of inertia)
                principal(Natmxmol, Index, atom, r1, trans, I);
#if 0
                if (D)
                {
                    fprintf(stderr, "pI\t%d %g %g %g \n", id, I[XX], I[YY], I[ZZ]);
                }

                if (Natmxmol == 2)   // Goddard set so
                {
                    I[ZZ] = 0.0;
                }

                if (D)
                {
                    fprintf(stderr, "trans %g %g %g \n", trans[XX][XX], trans[XX][YY], trans[XX][ZZ]);
                    fprintf(stderr, "trans %g %g %g \n", trans[YY][XX], trans[YY][YY], trans[YY][ZZ]);
                    fprintf(stderr, "trans %g %g %g \n", trans[ZZ][XX], trans[ZZ][YY], trans[ZZ][ZZ]);
                }
#endif
                /* save the principal moment: it will be used to compute the rotational temperatures */
                /* if the I is close to zero then the following gives wrong solution due to the
                   propagation of the division. Threshold 5e-5 comes from benchmark on O2
                   so setting a high treshold */
                if (I[XX] < 5.0e-5)
                {
                    I[XX] = 0.0;
                }
                else
                {
                    vecI[id][XX] += I[XX];
                }
                if (I[YY] < 5.0e-5)
                {
                    I[YY] = 0.0;
                }
                else
                {
                    vecI[id][YY] += I[YY];
                }
                if (I[ZZ] < 5.0e-5)
                {
                    I[ZZ] = 0.0;
                }
                else
                {
                    vecI[id][ZZ] += I[ZZ];
                }

                // we run through the j-th atom of id-th molecule and we save velocities for i-th atoms from Index
                // reset the needed variables
                for (i = 0; i < DIM; i++)
                {
                    angmom[i] = 0.0;
                    pomega[i] = 0.0;
                    omega[i]  = 0.0;
                }

                vax = 0.0;
                vay = 0.0;
                vaz = 0.0;
                // Index = &index[id*Natmxmol];
                for (j = 0; j < Natmxmol; j++)
                {
                    ai = Index[j];
                    // compute the angular momentum angmom
                    clear_rvec(tmp);
                    cprod(r1[ai], v1[ai], tmp);
                    // mass weigth it for each atom in the molecule
                    angmom[XX] += m[ai]*tmp[XX];
                    angmom[YY] += m[ai]*tmp[YY];
                    angmom[ZZ] += m[ai]*tmp[ZZ];
                }
                // compute angular velocity along principle axis, first step
                for(j = 0; (j < DIM); j++)
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
                for(j = 0; (j< 3); j++) 
                {
                    for (k = 0; k < 3; k++)
                    {
                        omega[j] += pomega[k]*trans[j][k];
                    }
                }
                // calculate inertia weighted angular velocities and save
                for(j = 0; (j<DIM); j++)
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

#if 0
                    if (D)
                    {
                        fprintf(stderr, "I %d\tpI %g %g %g\n\tan %g %g %g\n\tvc %g %g %g\n\tpo %g %g %g\n\tom %g %g %g\n\tvr %g %g %g\n\tvv %g %g %g\n", i, I[XX], I[YY], I[ZZ],
                                angmom[XX], angmom[YY], angmom[ZZ],
                                v1[ai][XX], v1[ai][YY], v1[ai][ZZ],
                                pomega[XX], pomega[YY], pomega[ZZ],
                                omega[XX], omega[YY], omega[ZZ],
                                vr[XX], vr[YY], vr[ZZ],
                                vec[V_V][ai*DIM+XX][nframes], vec[V_V][ai*DIM+YY][nframes], vec[V_V][ai*DIM+ZZ][nframes]);
                    }
#endif
                }
            } //(Natmxmol>1.0){
        }     // for (id=0; id<Nmol; id++) {

        t1 = fr.time;
        nframes++;
    }
    while (read_next_frame(oenv, status, &fr));

    gmx_rmpbc_done(gpbc);

    sfree(r1);   // empty r1, we need memory
    sfree(v1);   // empty v1, we need memory
    sfree(atom); // empty atom, we need memory
    close_trj(status);

    dt = (t1-t0)/(nframes-1);
    if (nV > 0)
    {
        V = Vsum/nV;
    }

    if (bVerbose)
    {
        printf("Going to do %d fourier transforms of length %d. Hang on.\n",
               gnx*VEL_NR, nframes);
    }

    // V_T is counted on number of molecules
    low_do_autocorr(NULL, oenv, NULL, nframes, Nmol*DIM, nframes, vec[V_T], dt, eacNormal, 0, FALSE,
                    FALSE, FALSE, -1, -1, 0);
    low_do_autocorr(NULL, oenv, NULL, nframes, gnx, nframes, vec[VEL_TOT], dt, eacNormal, 0, FALSE,
                    FALSE, FALSE, -1, -1, 0);
    if (Natmxmol > 1)
    {
        // V_A will be counted on number of molecules
        low_do_autocorr(NULL, oenv, NULL, nframes, Nmol*DIM, nframes, vec[V_A], dt, eacNormal, 0, FALSE,
                        FALSE, FALSE, -1, -1, 0);
        low_do_autocorr(NULL, oenv, NULL, nframes, gnx, nframes, vec[V_V], dt, eacNormal, 0, FALSE,
                        FALSE, FALSE, -1, -1, 0);
        low_do_autocorr(NULL, oenv, NULL, nframes, gnx, nframes, vec[V_R], dt, eacNormal, 0, FALSE,
                        FALSE, FALSE, -1, -1, 0);
    }

/*
   enum Velocity { V_T, V_A, V_V, V_R, VEL_TOT, VEL_NR }; trans, angular, vibrational, rotational, total, STOP
   enum Density { VACF=VEL_NR, DOS_SOLID, DOS_DIFF, DOS_CP, DOS_S, DOS_A, DOS_E, DOS_NR };
 */
    snew(dos, DOS_NR);
    for (j = 0; (j < DOS_NR); j++)
    {
        /* snew(dos[j],nframes+4); */
        snew(dos[j], nframes); //(nframes+1)/2);// in the next we sum up to nframes/2

    }
    if (bVerbose)
    {
        printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");
    }
    for (t = V_T; t < VEL_NR; t++)
    {
        switch (t)
        {
            case V_T: // center of mass -> number of molecules
                for (i = 0; (i < 3*Nmol); i += DIM)
                {     /* originally on nframes/2 */
                    for (j = 0; (j < nframes); j++)
                    {
                        v1j          = vec[V_T][i+XX][j] + vec[V_T][i+YY][j] + vec[V_T][i+ZZ][j];
                        dos[V_T][j] += v1j*tm;
                    }
                }
                break;
            case V_A:
                if (Natmxmol > 1)  // angular velocity of the molecules -> number of molecules
                {
                    for (i = 0; (i < 3*Nmol); i += DIM)
                    { /* originally on nframes/2 */
                        for (j = 0; (j < nframes); j++)
                        {
                            v1j          = vec[V_A][i+XX][j] + vec[V_A][i+YY][j] + vec[V_A][i+ZZ][j]; // already I weighted
                            dos[V_A][j] += v1j;
                        }
                    }
                }
                break;
            case V_V:
                if (Natmxmol > 1)
                {
                    real am;
                    for (i = 0; (i < gnx); i += DIM)
                    { /* originally on nframes/2 */
                        am = m[i/DIM];
                        for (j = 0; (j < nframes); j++)
                        {
                            v1j          = vec[V_V][i+XX][j] + vec[V_V][i+YY][j] + vec[V_V][i+ZZ][j];
                            dos[V_V][j] += v1j*am;
                        }
                    }
                }
                break;
            case V_R:
                if (Natmxmol > 1)
                {
                    real am;
                    for (i = 0; (i < gnx); i += DIM)
                    { /* originally on nframes/2 */
                        am = m[i/DIM];
                        for (j = 0; (j < nframes); j++)
                        {
                            v1j          = vec[V_R][i+XX][j] + vec[V_R][i+YY][j] + vec[V_R][i+ZZ][j];
                            dos[V_R][j] += v1j*am;
                        }
                    }
                }
                break;
            case VEL_TOT:
            {
                real iNatom = 1.0/Natom;
                real am;
                for (i = 0; (i < gnx); i += DIM)
                {   /* originally on nframes/2 */
                    am = m[i/DIM];
                    for (j = 0; (j < nframes); j++)
                    {
                        v1j              = vec[VEL_TOT][i+XX][j] + vec[VEL_TOT][i+YY][j] + vec[VEL_TOT][i+ZZ][j];
                        dos[VACF][j]    += v1j*iNatom;
                        dos[VEL_TOT][j] += v1j*am;
                    }
                }
            }
            break;
            default:
                break;
        }
    }
    /* Compute temperatures, degrees of freedom must be changed to the right one */
    if (Natmxmol == 1)
    {
        dof[V_T]      = 3.0;
        dof[V_A]      = 0.0;
        dof[V_V]      = 0.0;
        dof[V_R]      = 0.0;
        dof[VEL_TOT]  = 3.0;
        vecT[V_T]     = dos[V_T][0]/(dof[V_T]*BOLTZ*Natom);
        vecT[VEL_TOT] = dos[VEL_TOT][0]/(dof[VEL_TOT]*BOLTZ*Natom);

    }
    else if (Natmxmol == 2)   /* Corrected degrees of freedom as in Goddard */
    {
        real corr = 5.0;
        dofT          = 6.0;
        dof[V_T]      = 3.0*(1-corr/dofT);
        dof[V_A]      = 2.0*(1-corr/dofT);
        dof[V_V]      = 1.0*(1-corr/dofT);
        dof[V_R]      = 2.0*(1-corr/dofT);
        dof[VEL_TOT]  = 3.0;
        vecT[V_T]     = dos[V_T][0]/(dof[V_T]*BOLTZ*Nmol);
        vecT[V_A]     = dos[V_A][0]/(dof[V_A]*BOLTZ*Nmol);
        vecT[V_V]     = dos[V_V][0]/(dof[V_V]*BOLTZ*Natom);
        vecT[V_R]     = dos[V_R][0]/(dof[V_R]*BOLTZ*Natom);
        vecT[VEL_TOT] = dos[VEL_TOT][0]/(dof[VEL_TOT]*BOLTZ*Natom);
    }
    else
    {
        dof[VEL_TOT]  = 3*Natmxmol;
        dof[V_T]      = 3;
        dof[V_A]      = 3;
        dof[V_V]      = dof[VEL_TOT] - dof[V_T] - dof[V_A];
        dof[V_R]      = dof[V_A];
        dofT          = dof[V_T]+dof[V_A]+dof[V_V];
        vecT[V_T]     = dos[V_T][0]/(dof[V_T]*BOLTZ*Nmol);
        vecT[V_A]     = dos[V_A][0]/(dof[V_A]*BOLTZ*Nmol);
        vecT[V_V]     = dos[V_V][0]/(dof[V_V]*BOLTZ*Nmol);
        vecT[V_R]     = dos[V_R][0]/(dof[V_R]*BOLTZ*Nmol);
        vecT[VEL_TOT] = dos[VEL_TOT][0]/(dof[VEL_TOT]*BOLTZ*Nmol);

    }
    //real corr = 6.0;
    //printf("\t\t\t%g %g %g %g %g\n",dof[V_T]*(1-corr/dofT),dof[V_A]*(1-corr/dofT),dof[V_V]*(1-corr/dofT),dof[VEL_TOT],dofT);
    /* from here, the temperature is taken from simulation */
#if 1
    fprintf(stdout, "vecT[V_T] %g vecT[V_A] %g vecT[V_V] %g vecT[V_R] %g vecT[VEL_TOT] %g\n",
            vecT[V_T], vecT[V_A], vecT[V_V], vecT[V_R], vecT[VEL_TOT]);
#endif
    fp = xvgropen(opt2fn("-vacf", NFILE, fnm), "Velocity ACF",
                  "Time (ps)", "C(t)", oenv);
    xvgr_legend(fp, asize(Mvacflegend), Mvacflegend, oenv);
    snew(tt, nframes/2);
    for (j = 0; (j < nframes/2); j++)
    {
        tt[j] = j*dt;
        fprintf(fp, "%10g %10g\n", tt[j], dos[VACF][j]);
    }
    fclose(fp);

    fp = xvgropen(opt2fn("-mvacf", NFILE, fnm), "Mass-weighted velocity ACF",
                  "Time (ps)", "C(t)", oenv);
    xvgr_legend(fp, asize(Mvacflegend), Mvacflegend, oenv);
    for (j = 0; (j < nframes/2); j++)
    {
        fprintf(fp, "%10g %10g %10g %10g %10g %10g\n", tt[j], dos[V_T][j], dos[V_A][j], dos[V_V][j], dos[V_R][j], dos[VEL_TOT][j]);
    }
    fclose(fp);

    if ((fftcode = gmx_fft_init_1d_real(&fft, nframes,
                                        GMX_FFT_FLAG_NONE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d_real returned %d", fftcode);
    }
    //compute density of state for each kind of velocities, return in same vector
    for (t = V_T; t < VEL_NR; t++)
    {
        real *dost;
        snew(dost, nframes);
        for(j=0; (j<nframes/2); j++)
        {
            dost[2*j]   = dos[t][j];
            dost[2*j+1] = 0;
        }
        if ((fftcode = gmx_fft_1d_real(fft, GMX_FFT_COMPLEX_TO_REAL,
                                       (void *)dost, (void *)dos[t])) != 0)
        {
            gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
        }
        sfree(dost);
    }

    /* First compute the DoS */
    snew(bfac, VEL_NR);
    /* Magic factor of 8 included now. WHY MAGIC?  Do we need now? */
    /* Now a factor of 4 only DvdS 20140222 */
    for (i = V_T; i < VEL_NR; i++)
    {
        beta[i] = 1.0/(vecT[i]*BOLTZ);
        bfac[i] = 2.0*dt*beta[i];
        /* fprintf(stderr,"beta %g bfac %g \n", beta[i], bfac[i]); */
    }

    snew(dos2, VEL_NR);
    int nfreq = nframes/4;
    snew(nu, nfreq);

    for (j = 0; (j < nfreq); j++)
    {
        nu[j] = 2*j/(t1-t0);
        for (t = V_T; t < VEL_NR; t++)
        {
            real tmp = sqr(dos[t][2*j]) + sqr(dos[t][2*j+1]);
            dos2[t] += tmp; //sqr(dos[t][2*j]) + sqr(dos[t][2*j+1]);
            if (bAbsolute)
            {
                dos[t][j] = bfac[t]*sqrt(tmp); //(sqr(dos[t][2*j]) + sqr(dos[t][2*j+1]));
            }
            else
            {
                dos[t][j] = bfac[t]*dos[t][2*j];
            }
        }
    }
    fp = xvgropen("DoS.xvg", "DoS",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)", "DoS(v)", oenv);
    recip_fac = bRecip ? (1e7/SPEED_OF_LIGHT) : 1.0;
    xvgr_legend(fp, asize(Mvacflegend), Mvacflegend, oenv);
    for (j = 0; (j < nfreq); j++)
    {
        fprintf(fp, "%10g %10g %10g %10g %10g %10g\n",
                recip_fac*nu[j], dos[V_T][j]/recip_fac, dos[V_A][j]/recip_fac, dos[V_V][j]/recip_fac, dos[V_R][j]/recip_fac, dos[VEL_TOT][j]/recip_fac);
    }
    fclose(fp);
    /* Normalize it */
    dostot[VEL_TOT] = evaluate_integral(nfreq, nu, dos[VEL_TOT], NULL, nframes/4, &stddev);
    if (bNormalize)
    {
        for (j = 0; (j < nfreq); j++)
        {
            dos[VEL_TOT][j] *= 3*Natom/dostot[VEL_TOT];
        }
    }

    /* Analyze the translational DoS[V_T] */
    fprintf(stderr, "\nAnalyzing translational terms\n");
    DoS0 = dos[V_T][0];

    mmm   = tmass/Nmol;
    Delta = ((2*DoS0/(9*Nmol))*sqrt(M_PI*BOLTZ*vecT[V_T]*Nmol/tmass)*pow((Nmol/V), 1.0/3.0)*pow(6/M_PI, 2.0/3.0)); // eq 13 JPC B, 114, 8191 (2010)
    f[V_T]     = calc_fluidicity(Delta, toler);
    y     = calc_y(f[V_T], Delta, toler);
    z     = calc_compress(y);
    T_f  += f[V_T];
    Sig   = BOLTZ*(5.0/2.0+log(pow(2*mmm*M_PI*BOLTZ*vecT[V_T]/(sqr(PLANCK)), 1.5)*V/(f[V_T]*Natom))); // eq 16, Lin 2010, Sackur-Tetrode eq.
    Shs   = Sig+calc_Shs(f[V_T], y);
    rho   = (tmass*AMU)/(V*NANO*NANO*NANO);
    sigHS = pow(6*y*V/(M_PI*Natom), 1.0/3.0);
    /* Base level energy: that from the MD simulation minus the energy
     * of 3N harmonic oscillators. (Lin2003a, eq. 17)
     */
    V0 = (Emd - 3*Natom*BOLTZ*vecT[VEL_TOT])/Nmol;

    fprintf(fplog, "System = \"%s\"\n", title);
    fprintf(fplog, "Nmol = %d\n", Nmol);
    fprintf(fplog, "Natom = %d\n", Natom);
    fprintf(fplog, "dt = %g ps\n", dt);
    fprintf(fplog, "tmass = %g amu\n", tmass);
    fprintf(fplog, "V = %g nm^3\n", V);
    fprintf(fplog, "rho = %g g/l\n", rho);
    fprintf(fplog, "T = %g K\n", Temp);
    fprintf(fplog, "beta = %g mol/kJ\n", beta[V_T]);
    fprintf(fplog, "V0 = %g kJ/mol\n", V0);

    //fprintf(fplog,"\nDoS parameters\n");
    printf("\t\tAnalyzing translational terms\n");
    fprintf(fplog, "\t\tTransllational quantities\n");
    fprintf(fplog, "Delta = %g\n", Delta);
    fprintf(fplog, "fluidicity f = %g\n", f[V_T]);
    fprintf(fplog, "hard sphere packing fraction fy = %g\n", f[V_T]*y);
    fprintf(fplog, "hard sphere compressibility z = %g\n", z);
    fprintf(fplog, "ideal gas entropy Sig = %g\n", Sig);
    fprintf(fplog, "hard sphere entropy Shs = %g\n", Shs);
    fprintf(fplog, "sigma_HS = %g nm\n", sigHS);
    fprintf(fplog, "DoS0 = %g\n", DoS0);
    fprintf(fplog, "Dos2[V_T] = %g\n", dos2[V_T]);
    fprintf(fplog, "DoSTot[VEL_TOT] = %g (should be 3*Natom = %d?)\n", dostot[VEL_TOT], 3*Natom);
    fprintf(stdout, "DoSTot[VEL_TOT] = %g (should be 3*Natom = %d?)\n", dostot[VEL_TOT], 3*Natom);

    dostot[V_T] = evaluate_integral(nfreq, nu, dos[V_T], NULL, nfreq, &stddev);
    fprintf(fplog, "DoSTot[V_T] = %g\n", dostot[V_T]);
    /* Now compute solid (2) and diffusive (3) components */
    fp = xvgropen(opt2fn("-dos", NFILE, fnm), "Moldos Density of states",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                  "\\f{4}S(\\f{12}n\\f{4})", oenv);
    xvgr_legend(fp, asize(DoSlegend), DoSlegend, oenv);
    recip_fac = bRecip ? (1e7/SPEED_OF_LIGHT) : 1.0;
    for (j = 0; (j < nfreq); j++)
    {
        dos[DOS_DIFF][j]  = DoS0/(1+sqr(DoS0*M_PI*nu[j]/(6*f[V_T]*Natom)));
        dos[DOS_SOLID][j] = dos[V_T][j]-dos[DOS_DIFF][j];

        fprintf(fp, "%10g  %10g  %10g  %10g\n",
                recip_fac*nu[j],
                dos[V_T][j]/recip_fac,
                dos[DOS_SOLID][j]/recip_fac,
                dos[DOS_DIFF][j]/recip_fac);
    }
    fclose(fp);

    /* Finally analyze the results! */
    wCdiff = 0.5;
    wSdiff = Shs/(3*BOLTZ);
    wEdiff = 0.5;
    wAdiff = wEdiff-wSdiff;
    for (j = 0; (j < nfreq); j++)
    {
        dos[DOS_CP][j] = (dos[DOS_DIFF][j]*wCdiff +
                          dos[DOS_SOLID][j]*wCsolid(nu[j], beta[V_T]));
        dos[DOS_S][j]  = (dos[DOS_DIFF][j]*wSdiff +
                          dos[DOS_SOLID][j]*wSsolid(nu[j], beta[V_T]));
        dos[DOS_A][j]  = (dos[DOS_DIFF][j]*wAdiff +
                          dos[DOS_SOLID][j]*wAsolid(nu[j], beta[V_T]));
        dos[DOS_E][j]  = (dos[DOS_DIFF][j]*wEdiff +
                          dos[DOS_SOLID][j]*wEsolid(nu[j], beta[V_T]));
    }
#if 0
    for (j = 0; (j < nfreq); j++)
    {
        fprintf(stderr, "%g %g \n", nu[j], dos[DOS_CP][j]);
    }
#endif
    DiffCoeff = evaluate_integral(nframes/2, tt, dos[VACF], NULL, nframes/2, &stddev);
    DiffCoeff = 1000*DiffCoeff/3.0;
    fprintf(fplog, "Diffusion coefficient from dos[VACF] %g 10^-5 cm^2/s\n",
            DiffCoeff);
    fprintf(fplog, "Diffusion coefficient from DoS0 %g 10^-5 cm^2/s\n",
            1000*DoS0/(12*tmass*beta[V_T]));

    fprintf(stdout, "Diffusion coefficient from dos[VACF] %g 10^-5 cm^2/s\n", DiffCoeff);
    fprintf(stdout, "Diffusion coefficient from DoS0 %g 10^-5 cm^2/s\n", 1000*DoS0/(12*tmass*beta[V_T]));
    T_diff += 1000*DoS0/(12*tmass*beta[V_T]);

    cP[V_T] = BOLTZ * evaluate_integral(nfreq, nu, dos[DOS_CP], NULL, nfreq, &stddev);
    fprintf(fplog, "Heat capacity %g J/mol K\n", 1000*cP[V_T]/Nmol);

    S[V_T]  = BOLTZ * evaluate_integral(nfreq, nu, dos[DOS_S], NULL, nfreq, &stddev);
    fprintf(fplog, "Entropy %g J/mol K\n", 1000*S[V_T]/Nmol);

    A[V_T]  = BOLTZ * Temp * evaluate_integral(nfreq, nu, dos[DOS_A], NULL, nfreq, &stddev);
    fprintf(fplog, "Helmholtz energy %g kJ/mol\n", V0+A[V_T]/Nmol);

    E[V_T]  = BOLTZ * Temp * evaluate_integral(nfreq, nu, dos[DOS_E], NULL, nfreq, &stddev);
    fprintf(fplog, "Internal energy %g kJ/mol\n", V0+E[V_T]/Nmol);

    if (D)
    {
        fprintf(stdout, "Heat capacity %g J/mol K\n", 1000*cP[V_T]/Nmol);
        fprintf(stdout, "Entropy %g J/mol K\n", 1000*S[V_T]/Nmol);
        fprintf(stdout, "Helmholtz energy %g kJ/mol\n", V0+A[V_T]/Nmol);
        fprintf(stdout, "Internal energy %g kJ/mol\n", V0+E[V_T]/Nmol);
    }
    double T_CP = 0, T_S = 0, T_E = 0, T_A = 0;
    T_CP += cP[V_T];
    T_S  += S[V_T];
    T_E  += E[V_T];
    T_A  += A[V_T];

    /*                                 */
    /* Analyze the angular dos[V_A] */
    /*                                 */

    /* add a check if we really need to compute it, for example if Natom!=Nmol */
    /* Compute the characteristic rotational temperature http://en.wikipedia.org/wiki/Rotational_temperature */

    /* reset rT */

    rT[XX] = rT[YY] = rT[ZZ] = 0.0;
    I[XX]  = I[YY] = I[ZZ] = 0.0;
    for (id = 0; id < Nmol; id++)
    {
        vecI[id][XX] /= nframes;
        vecI[id][YY] /= nframes;
        vecI[id][ZZ] /= nframes;
    }

    for (i = 0; i < Nmol; i++)   // run for each molecule
    {
        real       tmp = 0.0;
        const real Q   = PLANCK*PLANCK/(8.0*M_PI*M_PI*BOLTZ);

        if (vecI[i][XX] > 0.0)
        {
            tmp     = Q/vecI[i][XX];
            rT[XX] += tmp;
            I[XX]  += vecI[i][XX];
        }
        if (vecI[i][YY] > 0.0)
        {
            tmp     = Q/vecI[i][YY];
            rT[YY] += tmp;
            I[YY]  += vecI[i][YY];
        }
        if (vecI[i][ZZ] > 0.0)
        {
            tmp     = Q/vecI[i][ZZ];
            rT[ZZ] += tmp;
            I[ZZ]  += vecI[i][ZZ];
        }
    }

    I[XX] = I[XX]/Nmol;
    I[YY] = I[YY]/Nmol;
    I[ZZ] = I[ZZ]/Nmol;

    fprintf(stdout, "Momenta of Inertia Ix %g Iy %g Iz %g\n", I[XX], I[YY], I[ZZ]);
    rT[XX] = rT[XX]/Nmol;
    rT[YY] = rT[YY]/Nmol;
    rT[ZZ] = rT[ZZ]/Nmol;
    if (!D)
    {
        fprintf(stdout, "Oxygen rT = 2.08K\nCharacteristic Rotational T x %g y %g z %g\n", rT[XX], rT[YY], rT[ZZ]);
    }
    fprintf(stdout, "Characteristic Rotational T x %g y %g z %g\n", rT[XX], rT[YY], rT[ZZ]);

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

    if (crT > 0.0)
    {
        printf("\n\t\tAnalyzing angular terms\n");
        fprintf(fplog, "\t\tAngular quantities\n");

        /* rot= rotation simmetry */
        real rot = 1.0, SR;
        DoS0 = dos[V_A][0];

        // Goddard use number of molecule
        Delta = ((2*DoS0/(9*Nmol))*sqrt(M_PI*BOLTZ*vecT[V_A]*Nmol/tmass)*pow((Nmol/V), 1.0/3.0)*pow(6/M_PI, 2.0/3.0)); // eq 13 JPC B, 114, 8191 (2010)
        f[V_A]     = calc_fluidicity(Delta, toler);
        y     = calc_y(f[V_A], Delta, toler);
        z     = calc_compress(y);
        T_f  += f[V_A];
        fp    = xvgropen("moldos_ang.xvg", "Angular Component",
                         bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                         "\\f{4}S(\\f{12}n\\f{4})", oenv);
        xvgr_legend(fp, asize(DoSlegend), DoSlegend, oenv);

        for (j = 0; (j < nfreq); j++)
        {
            dos[DOS_DIFF][j]  = DoS0/(1+sqr(DoS0*M_PI*nu[j]/(6*f[V_A]*Natom)));
            dos[DOS_SOLID][j] = dos[V_A][j]-dos[DOS_DIFF][j];
            fprintf(fp, "%10g  %10g  %10g  %10g\n",
                    recip_fac*nu[j],
                    dos[V_A][j]/recip_fac,
                    dos[DOS_SOLID][j]/recip_fac,
                    dos[DOS_DIFF][j]/recip_fac);
        }
        fclose(fp);

        //fprintf(fplog,"\nDoS parameters\n");
        fprintf(fplog, "Delta = %g\n", Delta);
        fprintf(fplog, "fluidicity f = %g\n", f[V_A]);
        fprintf(fplog, "hard sphere packing fraction fy = %g\n", f[V_A]*y);
        fprintf(fplog, "hard sphere compressibility z = %g\n", z);
        fprintf(fplog, "ideal gas entropy Sig = %g\n", Sig);
        fprintf(fplog, "hard sphere entropy Shs = %g\n", Shs);
        fprintf(fplog, "sigma_HS = %g nm\n", sigHS);
        fprintf(fplog, "DoS0 = %g\n", DoS0);
        fprintf(fplog, "Dos2[V_A] = %g\n", dos2[V_A]);

        dostot[V_A] = evaluate_integral(nfreq, nu, dos[V_A], NULL, nfreq, &stddev);
        fprintf(fplog, "DoSTot[V_A] = %g\n", dostot[V_A]);

        fprintf(stdout, "Diffusion coefficient from DoS0[V_A] %g 10^-5 cm^2/s\n", 1000*DoS0/(12*tmass*beta[V_A]));
        /* SR/BOLTZ Rotational entropy of rigid body solid with rotational temperature rT[XX, YY, ZZ] */
        T_diff += 1000*DoS0/(12*tmass*beta[V_A]);
        SR      = log(((pow(M_PI, 0.5)*exp(3.0/2.0)/rot))*pow(vecT[V_A], rdof)/crT);

        /* Finally analyze the results! */
        wCdiff = 0.5;
        wSdiff = SR/(3*BOLTZ);
        wEdiff = 0.5;
        wAdiff = wEdiff-wSdiff;
        for (j = 0; (j < nfreq); j++)
        {
            dos[DOS_CP][j] = (dos[DOS_DIFF][j]*wCdiff +
                              dos[DOS_SOLID][j]*wCsolid(nu[j], beta[V_A]));
            dos[DOS_S][j]  = (dos[DOS_DIFF][j]*wSdiff +
                              dos[DOS_SOLID][j]*wSsolid(nu[j], beta[V_A]));
            dos[DOS_A][j]  = (dos[DOS_DIFF][j]*wAdiff +
                              dos[DOS_SOLID][j]*wAsolid(nu[j], beta[V_A]));
            dos[DOS_E][j]  = (dos[DOS_DIFF][j]*wEdiff +
                              dos[DOS_SOLID][j]*wEsolid(nu[j], beta[V_A]));
        }

        cP[V_A] = BOLTZ * evaluate_integral(nfreq, nu, dos[DOS_CP], NULL, nfreq, &stddev);
        fprintf(fplog, "Heat capacity %g J/mol K\n", 1000*cP[V_A]/Nmol);

        S[V_A]  = BOLTZ * evaluate_integral(nfreq, nu, dos[DOS_S], NULL, nfreq, &stddev);
        fprintf(fplog, "Entropy %g J/mol K\n", 1000*S[V_A]/Nmol);

        A[V_A]  = BOLTZ * Temp * evaluate_integral(nfreq, nu, dos[DOS_A], NULL, nfreq, &stddev);
        fprintf(fplog, "Helmholtz energy %g kJ/mol\n", V0+A[V_A]/Nmol);

        E[V_A]  = BOLTZ * Temp * evaluate_integral(nfreq, nu, dos[DOS_E], NULL, nfreq, &stddev);
        fprintf(fplog, "Internal energy %g kJ/mol\n", V0+E[V_A]/Nmol);

        if (D)
        {
            fprintf(stdout, "Heat capacity %g J/mol K\n", 1000*cP[V_A]/Nmol);
            fprintf(stdout, "Entropy %g J/mol K\n", 1000*S[V_A]/Nmol);
            fprintf(stdout, "Helmholtz energy %g kJ/mol\n", V0+A[V_A]/Nmol);
            fprintf(stdout, "Internal energy %g kJ/mol\n", V0+E[V_A]/Nmol);
        }
        T_CP += cP[V_A];
        T_S  += S[V_A];
        T_E  += E[V_A];
        T_A  += A[V_A];
    }

    /* Analyze the vibrational dos[V_V] */
    fprintf(fplog, "\n\t\tVibrational quantities\n");

    if (vecT[V_V] > 0.0)
    {

        dostot[V_V] = evaluate_integral(nfreq, nu, dos[V_V], NULL, nfreq, &stddev);
        fprintf(fplog, "DoSTot[V_V] = %g\n", dostot[V_V]);

        /* add a check if we really need to compute it, for example if Natom!=Nmol */
        printf("\n\t\tAnalyzing vibrational terms\n");
        for (j = 0; (j < nfreq); j++)
        {
            dos[DOS_CP][j] = dos[V_V][j]*wCsolid(nu[j], beta[V_V]);
            dos[DOS_S][j]  = dos[V_V][j]*wSsolid(nu[j], beta[V_V]);
            dos[DOS_A][j]  = dos[V_V][j]*wAsolid(nu[j], beta[V_V]);
            dos[DOS_E][j]  = dos[V_V][j]*wEsolid(nu[j], beta[V_V]);
        }

        cP[V_V] = BOLTZ * evaluate_integral(nfreq, nu, dos[DOS_CP], NULL, nfreq, &stddev);
        fprintf(fplog, "Heat capacity %g J/mol K\n", 1000*cP[V_V]/Nmol);

        S[V_V]  = BOLTZ * evaluate_integral(nfreq, nu, dos[DOS_S], NULL, nfreq, &stddev);
        fprintf(fplog, "Entropy %g J/mol K\n", 1000*S[V_V]/Nmol);

        A[V_V]  = BOLTZ * Temp * evaluate_integral(nfreq, nu, dos[DOS_A], NULL, nfreq, &stddev);
        fprintf(fplog, "Helmholtz energy %g kJ/mol\n", V0+A[V_V]/Nmol);

        E[V_V]  = BOLTZ * Temp * evaluate_integral(nfreq, nu, dos[DOS_E], NULL, nfreq, &stddev);
        fprintf(fplog, "Internal energy %g kJ/mol\n", V0+E[V_V]/Nmol);
        if (D)
        {
            fprintf(stdout, "Heat capacity %g J/mol K\n", 1000*cP[V_V]/Nmol);
            fprintf(stdout, "Entropy %g J/mol K\n", 1000*S[V_V]/Nmol);
            fprintf(stdout, "Helmholtz energy %g kJ/mol\n", V0+A[V_V]/Nmol);
            fprintf(stdout, "Internal energy %g kJ/mol\n", V0+E[V_V]/Nmol);
        }
        T_CP += cP[V_V];
        T_S  += S[V_V];
        T_E  += E[V_V];
        T_A  += A[V_V];
    }

    /* ---------------------- */
    fprintf(fplog, "\n\t\t Final resume \n");
    fprintf(fplog, "Temperature:\nT[V_T] %g T[V_A] %g T[V_V] %g T[V_R] %g T[VEL_TOT] %g\n", 
            vecT[V_T], vecT[V_A], vecT[V_V], vecT[V_R], vecT[VEL_TOT]);
    fprintf(fplog, "Degrees of freedom:\n V_T %g V_A %g V_V %g VEL_TOT %g\n", 
            dof[V_T], dof[V_A], dof[V_V], dof[VEL_TOT]);
    fprintf(fplog, "Dos integral:\n V_T %g V_A %g V_V %g VEL_TOT %g theoretical %d\n", 
            dostot[V_T], dostot[V_A], dostot[V_V], dostot[VEL_TOT], gnx);
    fprintf(fplog, "Dos0:\n V_T %g V_A %g V_V %g VEL_TOT %g\n", 
            dos[V_T][0], dos[V_A][0], dos[V_V][0], dos[VEL_TOT][0]);
    
    {
        const char *items[] = { "", "Temperature", "DOF", "Dos integral", "DoS0",
                                "cV", "A", "S", "E", "fluidicity" }; 
        
        for(k = 0; (k < asize(items)); k++) 
        {
            fprintf(fplog, "%-12s", items[k]);
            for(t = V_T; (t <= VEL_TOT); t++)
            {
                switch (k) 
                {
                case 0:
                    fprintf(fplog, "  %10s", Mvacflegend[t]);
                    break;
                case 1:
                    fprintf(fplog, "  %10g", vecT[t]);
                    break;
                case 2:
                    fprintf(fplog, "  %10g", dof[t]);
                    break;
                case 3:
                    fprintf(fplog, "  %10g", dostot[t]);
                    break;
                case 4:
                    fprintf(fplog, "  %10g", dos[t][0]);
                    break;
                case 5:
                    fprintf(fplog, "  %10g", cP[t]);
                    break;
                case 6:
                    fprintf(fplog, "  %10g", A[t]);
                    break;
                case 7:
                    fprintf(fplog, "  %10g", S[t]);
                    break;
                case 8:
                    fprintf(fplog, "  %10g", E[t]);
                    break;
                case 9:
                    fprintf(fplog, "  %10g", f[t]);
                    break;
                }
            }
            fprintf(fplog, "\n");
        }
    }
    fprintf(fplog, "Total Fluidicity %g \n", T_f);
    fprintf(fplog, "Total diffusivity %g 10^-5 cm^2/s\n", T_diff);
    fprintf(fplog, "Total Heat capacity %g J/mol K\n", 1000*T_CP/Nmol);
    fprintf(fplog, "Total Entropy %g J/mol K\n", 1000*T_S/Nmol);
    fprintf(fplog, "Total Helmholtz energy %g kJ/mol (wrong due to the unknown V0)\n", T_A/Nmol);
    fprintf(fplog, "Total Internal energy %g kJ/mol (wrong due to the unknown V0)\n", T_E/Nmol);

    fprintf(fplog, "\nArrivederci!\n");
    fclose(fplog);
    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    //    thanx(stderr);  give error in my compilation

    return 0;
}
