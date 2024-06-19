/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

enum
{
    VACF,
    MVACF,
    DOS,
    DOS_SOLID,
    DOS_DIFF,
    DOS_CP,
    DOS_S,
    DOS_A,
    DOS_E,
    DOS_NR
};

static int calcMoleculesInIndexGroup(const t_block* mols, int natoms, const int* index, int nindex)
{
    int i    = 0;
    int mol  = 0;
    int nMol = 0;
    int j;

    while (i < nindex)
    {
        while (index[i] > mols->index[mol])
        {
            mol++;
            if (mol >= mols->nr)
            {
                gmx_fatal(FARGS, "Atom index out of range: %d", index[i] + 1);
            }
        }
        for (j = mols->index[mol]; j < mols->index[mol + 1]; j++)
        {
            if (index[i] != j)
            {
                gmx_fatal(FARGS, "The index group does not consist of whole molecules");
            }
            i++;
            if (i > natoms)
            {
                gmx_fatal(FARGS, "Index contains atom numbers larger than the topology");
            }
        }
        nMol++;
    }
    return nMol;
}

static double FD(double Delta, double f)
{
    return (2 * std::pow(Delta, -4.5) * std::pow(f, 7.5)
            - 6 * std::pow(Delta, -3.0) * std::pow(f, 5.0) - std::pow(Delta, -1.5) * std::pow(f, 3.5)
            + 6 * std::pow(Delta, -1.5) * std::pow(f, 2.5) + 2 * f - 2);
}

static double YYY(double f, double y)
{
    return (2 * gmx::power3(y * f) - gmx::square(f) * y * (1 + 6 * y) + (2 + 6 * y) * f - 2);
}

static double calc_compress(double y)
{
    if (y == 1)
    {
        return 0;
    }
    return ((1 + y + gmx::square(y) - gmx::power3(y)) / (gmx::power3(1 - y)));
}

static double bisector(double Delta, double tol, double ff0, double ff1, double ff(double, double))
{
    double fd, f, f0, f1;
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
        f  = (f0 + f1) * 0.5;
        fd = ff(Delta, f);
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
    } while ((f1 - f0) > tol);

    return f;
}

static double calc_fluidicity(double Delta, double tol)
{
    return bisector(Delta, tol, 0, 1, FD);
}

static double calc_y(double f, double Delta, double toler)
{
    double y1, y2;

    y1 = std::pow(f / Delta, 1.5);
    y2 = bisector(f, toler, 0, 10000, YYY);
    if (std::abs((y1 - y2) / (y1 + y2)) > 100 * toler)
    {
        fprintf(stderr, "Inconsistency computing y: y1 = %f, y2 = %f, using y1.\n", y1, y2);
    }

    return y1;
}

static double calc_Shs(double f, double y)
{
    double fy = f * y;

    return gmx::c_boltz * (std::log(calc_compress(fy)) + fy * (3 * fy - 4) / gmx::square(1 - fy));
}

static real wCsolid(real nu, real beta)
{
    real bhn = beta * gmx::c_planck * nu;
    real ebn, koko;

    if (bhn == 0)
    {
        return 1.0;
    }
    else
    {
        ebn  = std::exp(bhn);
        koko = gmx::square(1 - ebn);
        return gmx::square(bhn) * ebn / koko;
    }
}

static real wSsolid(real nu, real beta)
{
    real bhn = beta * gmx::c_planck * nu;

    if (bhn == 0)
    {
        return 1;
    }
    else
    {
        return bhn / std::expm1(bhn) - std::log1p(-std::exp(-bhn));
    }
}

static real wAsolid(real nu, real beta)
{
    real bhn = beta * gmx::c_planck * nu;

    if (bhn == 0)
    {
        return 0;
    }
    else
    {
        return std::log((1 - std::exp(-bhn)) / (std::exp(-bhn / 2))) - std::log(bhn);
    }
}

static real wEsolid(real nu, real beta)
{
    real bhn = beta * gmx::c_planck * nu;

    if (bhn == 0)
    {
        return 1;
    }
    else
    {
        return bhn / 2 + bhn / std::expm1(bhn) - 1;
    }
}

int gmx_dos(int argc, char* argv[])
{
    const char* desc[] = { "[THISMODULE] computes the Density of States from a simulations.",
                           "In order for this to be meaningful the velocities must be saved",
                           "in the trajecotry with sufficiently high frequency such as to cover",
                           "all vibrations. For flexible systems that would be around a few fs",
                           "between saving. Properties based on the DoS are printed on the",
                           "standard output.",
                           "Note that the density of states is calculated from the mass-weighted",
                           "autocorrelation, and by default only from the square of the real",
                           "component rather than absolute value. This means the shape can differ",
                           "substantially from the plain vibrational power spectrum you can",
                           "calculate with gmx velacc." };
    const char* bugs[] = {
        "This program needs a lot of memory: total usage equals the number of atoms times "
        "3 times number of frames times 4 (or 8 when run in double precision)."
    };
    FILE *            fp, *fplog;
    t_topology        top;
    PbcType           pbcType = PbcType::Unset;
    t_trxframe        fr;
    matrix            box;
    int               gnx;
    real              t0, t1;
    t_trxstatus*      status;
    int               nV, nframes, n_alloc, i, j, fftcode, Nmol, Natom;
    double            rho, dt, Vsum, V, tmass, dostot, dos2;
    real **           c1, **dos, mi, beta, bfac, *nu, *tt, stddev, c1j;
    gmx_output_env_t* oenv;
    gmx_fft_t         fft;
    double            cP, DiffCoeff, Delta, f, y, z, sigHS, Shs, Sig, DoS0, recip_fac;
    double            wCdiff, wSdiff, wAdiff, wEdiff;
    int               grpNatoms;
    int*              index;
    char*             grpname;
    double            invNormalize;
    gmx_bool          normalizeAutocorrelation;

    static gmx_bool bVerbose = TRUE, bAbsolute = FALSE, bNormalizeDos = FALSE;
    static gmx_bool bRecip = FALSE;
    static real     Temp = 298.15, toler = 1e-6;
    int             min_frames = 100;

    t_pargs pa[] = {
        { "-v", FALSE, etBOOL, { &bVerbose }, "Be loud and noisy." },
        { "-recip",
          FALSE,
          etBOOL,
          { &bRecip },
          "Use cm^-1 on X-axis instead of 1/ps for DoS plots." },
        { "-abs",
          FALSE,
          etBOOL,
          { &bAbsolute },
          "Use the absolute value of the Fourier transform of the VACF as the Density of States. "
          "Default is to use the real component only" },
        { "-normdos",
          FALSE,
          etBOOL,
          { &bNormalizeDos },
          "Normalize the DoS such that it adds up to 3N. This should usually not be necessary." },
        { "-T", FALSE, etREAL, { &Temp }, "Temperature in the simulation" },
        { "-toler",
          FALSE,
          etREAL,
          { &toler },
          "HIDDENTolerance when computing the fluidicity using bisection algorithm" }
    };

    t_filenm fnm[] = {
        { efTRN, "-f", nullptr, ffREAD },      { efTPR, "-s", nullptr, ffREAD },
        { efNDX, nullptr, nullptr, ffOPTRD },  { efXVG, "-vacf", "vacf", ffWRITE },
        { efXVG, "-mvacf", "mvacf", ffWRITE }, { efXVG, "-dos", "dos", ffWRITE },
        { efLOG, "-g", "dos", ffWRITE },
    };
#define NFILE asize(fnm)
    int                        npargs;
    t_pargs*                   ppa;
    std::array<std::string, 3> DoSlegend = { "DoS(v)", "DoS(v)[Solid]", "DoS(v)[Diff]" };

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, npargs, ppa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    beta = 1 / (Temp * gmx::c_boltz);

    fplog = gmx_fio_fopen(ftp2fn(efLOG, NFILE, fnm), "w");
    fprintf(fplog, "Doing density of states analysis based on trajectory.\n");
    please_cite(fplog, "Pascal2011a");
    please_cite(fplog, "Caleman2011b");

    read_tps_conf(ftp2fn(efTPR, NFILE, fnm), &top, &pbcType, nullptr, nullptr, box, TRUE);

    /* Handle index groups */
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &grpNatoms, &index, &grpname);

    V     = det(box);
    tmass = 0;
    for (i = 0; i < grpNatoms; i++)
    {
        tmass += top.atoms.atom[index[i]].m;
    }

    Natom = grpNatoms;
    Nmol  = calcMoleculesInIndexGroup(&top.mols, top.atoms.nr, index, grpNatoms);
    gnx   = Natom * DIM;

    /* Correlation stuff */
    snew(c1, gnx);
    for (i = 0; (i < gnx); i++)
    {
        c1[i] = nullptr;
    }

    read_first_frame(oenv, &status, ftp2fn(efTRN, NFILE, fnm), &fr, TRX_NEED_V);
    t0 = fr.time;

    n_alloc = 0;
    nframes = 0;
    Vsum    = 0;
    nV      = 0;
    do
    {
        if (fr.bBox)
        {
            V = det(fr.box);
            Vsum += V;
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
            c1[i + XX][nframes] = fr.v[index[i / DIM]][XX];
            c1[i + YY][nframes] = fr.v[index[i / DIM]][YY];
            c1[i + ZZ][nframes] = fr.v[index[i / DIM]][ZZ];
        }

        t1 = fr.time;

        nframes++;
    } while (read_next_frame(oenv, status, &fr));

    close_trx(status);

    if (nframes < min_frames)
    {
        gmx_fatal(FARGS, "You need at least %d frames in the trajectory and you only have %d.", min_frames, nframes);
    }
    dt = (t1 - t0) / (nframes - 1);
    if (nV > 0)
    {
        V = Vsum / nV;
    }
    if (bVerbose)
    {
        printf("Going to do %d fourier transforms of length %d. Hang on.\n", gnx, nframes);
    }
    /* Unfortunately the -normalize program option for the autocorrelation
     * function calculation is added as a hack with a static variable in the
     * autocorrelation.c source. That would work if we called the normal
     * do_autocorr(), but this routine overrides that by directly calling
     * the low-level functionality. That unfortunately leads to ignoring the
     * default value for the option (which is to normalize).
     * Since the absolute value seems to be important for the subsequent
     * analysis below, we detect the value directly from the option, calculate
     * the autocorrelation without normalization, and then apply the
     * normalization just to the autocorrelation output
     * (or not, if the user asked for a non-normalized autocorrelation).
     */
    normalizeAutocorrelation = opt2parg_bool("-normalize", npargs, ppa);

    /* Note that we always disable normalization here, regardless of user settings */
    low_do_autocorr(
            nullptr, oenv, nullptr, nframes, gnx, nframes, c1, dt, eacNormal, 0, FALSE, FALSE, FALSE, -1, -1, 0);
    snew(dos, DOS_NR);
    for (j = 0; (j < DOS_NR); j++)
    {
        snew(dos[j], nframes + 4);
    }

    if (bVerbose)
    {
        printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");
    }
    for (i = 0; (i < gnx); i += DIM)
    {
        mi = top.atoms.atom[index[i / DIM]].m;
        for (j = 0; (j < nframes / 2); j++)
        {
            c1j = (c1[i + XX][j] + c1[i + YY][j] + c1[i + ZZ][j]);
            dos[VACF][j] += c1j / Natom;
            dos[MVACF][j] += mi * c1j;
        }
    }

    fp = xvgropen(
            opt2fn("-vacf", NFILE, fnm), "Velocity autocorrelation function", "Time (ps)", "C(t)", oenv);
    snew(tt, nframes / 2);

    invNormalize = normalizeAutocorrelation ? 1.0 / dos[VACF][0] : 1.0;

    for (j = 0; (j < nframes / 2); j++)
    {
        tt[j] = j * dt;
        fprintf(fp, "%10g  %10g\n", tt[j], dos[VACF][j] * invNormalize);
    }
    xvgrclose(fp);

    fp = xvgropen(opt2fn("-mvacf", NFILE, fnm),
                  "Mass-weighted velocity autocorrelation function",
                  "Time (ps)",
                  "C(t)",
                  oenv);

    invNormalize = normalizeAutocorrelation ? 1.0 / dos[VACF][0] : 1.0;

    for (j = 0; (j < nframes / 2); j++)
    {
        fprintf(fp, "%10g  %10g\n", tt[j], dos[MVACF][j] * invNormalize);
    }
    xvgrclose(fp);

    if ((fftcode = gmx_fft_init_1d_real(&fft, nframes / 2, GMX_FFT_FLAG_NONE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d_real returned %d", fftcode);
    }
    if ((fftcode = gmx_fft_1d_real(fft, GMX_FFT_REAL_TO_COMPLEX, dos[MVACF], dos[DOS])) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
    }

    /* First compute the DoS */
    /* Magic factor of 8 included now. */
    bfac = 8 * dt * beta / 2;
    dos2 = 0;
    snew(nu, nframes / 4);
    for (j = 0; (j < nframes / 4); j++)
    {
        nu[j] = 2 * j / (t1 - t0);
        dos2 += gmx::square(dos[DOS][2 * j]) + gmx::square(dos[DOS][2 * j + 1]);
        if (bAbsolute)
        {
            dos[DOS][j] = bfac * std::hypot(dos[DOS][2 * j], dos[DOS][2 * j + 1]);
        }
        else
        {
            dos[DOS][j] = bfac * dos[DOS][2 * j];
        }
    }
    /* Normalize it */
    dostot = evaluate_integral(nframes / 4, nu, dos[DOS], nullptr, int{ nframes / 4 }, &stddev);
    if (bNormalizeDos)
    {
        for (j = 0; (j < nframes / 4); j++)
        {
            dos[DOS][j] *= 3 * Natom / dostot;
        }
    }

    /* Now analyze it */
    DoS0 = dos[DOS][0];

    /* Note this eqn. is incorrect in Pascal2011a! */
    Delta = ((2 * DoS0 / (9 * Natom)) * std::sqrt(M_PI * gmx::c_boltz * Temp * Natom / tmass)
             * std::pow((Natom / V), 1.0 / 3.0) * std::pow(6.0 / M_PI, 2.0 / 3.0));
    f     = calc_fluidicity(Delta, toler);
    y     = calc_y(f, Delta, toler);
    z     = calc_compress(y);
    Sig   = gmx::c_boltz
          * (5.0 / 2.0
             + std::log(2 * M_PI * gmx::c_boltz * Temp / (gmx::square(gmx::c_planck)) * V / (f * Natom)));
    Shs   = Sig + calc_Shs(f, y);
    rho   = (tmass * gmx::c_amu) / (V * gmx::c_nano * gmx::c_nano * gmx::c_nano);
    sigHS = std::cbrt(6 * y * V / (M_PI * Natom));

    fprintf(fplog, "System = \"%s\"\n", *top.name);
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
    fp = xvgropen(opt2fn("-dos", NFILE, fnm),
                  "Density of states",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                  "\\f{4}S(\\f{12}n\\f{4})",
                  oenv);
    xvgrLegend(fp, DoSlegend, oenv);
    recip_fac = bRecip ? (1e7 / gmx::c_speedOfLight) : 1.0;
    for (j = 0; (j < nframes / 4); j++)
    {
        dos[DOS_DIFF][j]  = DoS0 / (1 + gmx::square(DoS0 * M_PI * nu[j] / (6 * f * Natom)));
        dos[DOS_SOLID][j] = dos[DOS][j] - dos[DOS_DIFF][j];
        fprintf(fp,
                "%10g  %10g  %10g  %10g\n",
                recip_fac * nu[j],
                dos[DOS][j] / recip_fac,
                dos[DOS_SOLID][j] / recip_fac,
                dos[DOS_DIFF][j] / recip_fac);
    }
    xvgrclose(fp);

    /* Finally analyze the results! */
    wCdiff = 0.5;
    wSdiff = Shs / (3 * gmx::c_boltz); /* Is this correct? */
    wEdiff = 0.5;
    wAdiff = wEdiff - wSdiff;
    for (j = 0; (j < nframes / 4); j++)
    {
        dos[DOS_CP][j] = (dos[DOS_DIFF][j] * wCdiff + dos[DOS_SOLID][j] * wCsolid(nu[j], beta));
        dos[DOS_S][j]  = (dos[DOS_DIFF][j] * wSdiff + dos[DOS_SOLID][j] * wSsolid(nu[j], beta));
        dos[DOS_A][j]  = (dos[DOS_DIFF][j] * wAdiff + dos[DOS_SOLID][j] * wAsolid(nu[j], beta));
        dos[DOS_E][j]  = (dos[DOS_DIFF][j] * wEdiff + dos[DOS_SOLID][j] * wEsolid(nu[j], beta));
    }
    DiffCoeff = evaluate_integral(nframes / 2, tt, dos[VACF], nullptr, nframes / 2., &stddev);
    DiffCoeff = 1000 * DiffCoeff / 3.0;
    fprintf(fplog, "Diffusion coefficient from VACF %g 10^-5 cm^2/s\n", DiffCoeff);
    fprintf(fplog, "Diffusion coefficient from DoS %g 10^-5 cm^2/s\n", 1000 * DoS0 / (12 * tmass * beta));

    cP = gmx::c_boltz
         * evaluate_integral(nframes / 4, nu, dos[DOS_CP], nullptr, int{ nframes / 4 }, &stddev);
    fprintf(fplog, "Heat capacity %g J/mol K\n", 1000 * cP / Nmol);
    fprintf(fplog, "\nArrivederci!\n");
    gmx_fio_fclose(fplog);

    do_view(oenv, ftp2fn(efXVG, NFILE, fnm), "-nxy");

    return 0;
}
