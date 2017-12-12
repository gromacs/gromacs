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
#include <cstdio>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static void index_atom2mol(int *n, int *index, const t_block *mols)
{
    int nat, i, nmol, mol, j;

    nat  = *n;
    i    = 0;
    nmol = 0;
    mol  = 0;
    while (i < nat)
    {
        while (index[i] > mols->index[mol])
        {
            mol++;
            if (mol >= mols->nr)
            {
                gmx_fatal(FARGS, "Atom index out of range: %d", index[i]+1);
            }
        }
        for (j = mols->index[mol]; j < mols->index[mol+1]; j++)
        {
            if (i >= nat || index[i] != j)
            {
                gmx_fatal(FARGS, "The index group does not consist of whole molecules");
            }
            i++;
        }
        index[nmol++] = mol;
    }

    fprintf(stderr, "\nSplit group of %d atoms into %d molecules\n", nat, nmol);

    *n = nmol;
}

static void precalc(const t_topology &top, real normm[])
{

    real mtot;
    int  i, j, k, l;

    for (i = 0; i < top.mols.nr; i++)
    {
        k    = top.mols.index[i];
        l    = top.mols.index[i+1];
        mtot = 0.0;

        for (j = k; j < l; j++)
        {
            mtot += top.atoms.atom[j].m;
        }

        for (j = k; j < l; j++)
        {
            normm[j] = top.atoms.atom[j].m/mtot;
        }

    }

}

static void calc_spectrum(int n, real c[], real dt, const char *fn,
                          gmx_output_env_t *oenv, gmx_bool bRecip)
{
    FILE     *fp;
    gmx_fft_t fft;
    int       i, status;
    real     *data;
    real      nu, omega, recip_fac;

    snew(data, n*2);
    for (i = 0; (i < n); i++)
    {
        data[i] = c[i];
    }

    if ((status = gmx_fft_init_1d_real(&fft, n, GMX_FFT_FLAG_NONE)) != 0)
    {
        gmx_fatal(FARGS, "Invalid fft return status %d", status);
    }
    if ((status = gmx_fft_1d_real(fft, GMX_FFT_REAL_TO_COMPLEX, data, data)) != 0)
    {
        gmx_fatal(FARGS, "Invalid fft return status %d", status);
    }
    fp = xvgropen(fn, "Vibrational Power Spectrum",
                  bRecip ? "\\f{12}w\\f{4} (cm\\S-1\\N)" :
                  "\\f{12}n\\f{4} (ps\\S-1\\N)",
                  "a.u.", oenv);
    /* This is difficult.
     * The length of the ACF is dt (as passed to this routine).
     * We pass the vacf with N time steps from 0 to dt.
     * That means that after FFT we have lowest frequency = 1/dt
     * then 1/(2 dt) etc. (this is the X-axis of the data after FFT).
     * To convert to 1/cm we need to have to realize that
     * E = hbar w = h nu = h c/lambda. We want to have reciprokal cm
     * on the x-axis, that is 1/lambda, so we then have
     * 1/lambda = nu/c. Since nu has units of 1/ps and c has gromacs units
     * of nm/ps, we need to multiply by 1e7.
     * The timestep between saving the trajectory is
     * 1e7 is to convert nanometer to cm
     */
    recip_fac = bRecip ? (1e7/SPEED_OF_LIGHT) : 1.0;
    for (i = 0; (i < n); i += 2)
    {
        nu    = i/(2*dt);
        omega = nu*recip_fac;
        /* Computing the square magnitude of a complex number, since this is a power
         * spectrum.
         */
        fprintf(fp, "%10g  %10g\n", omega, gmx::square(data[i])+gmx::square(data[i+1]));
    }
    xvgrclose(fp);
    gmx_fft_destroy(fft);
    sfree(data);
}

int gmx_velacc(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] computes the velocity autocorrelation function.",
        "When the [TT]-m[tt] option is used, the momentum autocorrelation",
        "function is calculated.[PAR]",
        "With option [TT]-mol[tt] the velocity autocorrelation function of",
        "molecules is calculated. In this case the index group should consist",
        "of molecule numbers instead of atom numbers.[PAR]",
        "By using option [TT]-os[tt] you can also extract the estimated",
        "(vibrational) power spectrum, which is the Fourier transform of the",
        "velocity autocorrelation function."
        "Be sure that your trajectory contains frames with velocity information",
        "(i.e. [TT]nstvout[tt] was set in your original [REF].mdp[ref] file),",
        "and that the time interval between data collection points is",
        "much shorter than the time scale of the autocorrelation."
    };

    static gmx_bool bMass = FALSE, bMol = FALSE, bRecip = TRUE;
    t_pargs         pa[]  = {
        { "-m", FALSE, etBOOL, {&bMass},
          "Calculate the momentum autocorrelation function" },
        { "-recip", FALSE, etBOOL, {&bRecip},
          "Use cm^-1 on X-axis instead of 1/ps for spectra." },
        { "-mol", FALSE, etBOOL, {&bMol},
          "Calculate the velocity acf of molecules" }
    };

    t_topology      top;
    int             ePBC = -1;
    t_trxframe      fr;
    matrix          box;
    gmx_bool        bTPS = FALSE, bTop = FALSE;
    int             gnx;
    int            *index;
    char           *grpname;
    /* t0, t1 are the beginning and end time respectively.
     * dt is the time step, mass is temp variable for atomic mass.
     */
    real              t0, t1, dt, mass;
    t_trxstatus      *status;
    int               counter, n_alloc, i, j, counter_dim, k, l;
    rvec              mv_mol;
    /* Array for the correlation function */
    real            **c1;
    real             *normm = nullptr;
    gmx_output_env_t *oenv;

#define NHISTO 360

    t_filenm  fnm[] = {
        { efTRN, "-f",    nullptr,   ffREAD  },
        { efTPS, nullptr,    nullptr,   ffOPTRD },
        { efNDX, nullptr,    nullptr,   ffOPTRD },
        { efXVG, "-o",    "vac",  ffWRITE },
        { efXVG, "-os",   "spectrum", ffOPTWR }
    };
#define NFILE asize(fnm)
    int       npargs;
    t_pargs  *ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    if (bMol || bMass)
    {
        bTPS = ftp2bSet(efTPS, NFILE, fnm) || !ftp2bSet(efNDX, NFILE, fnm);
    }

    if (bTPS)
    {
        bTop = read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &ePBC, nullptr, nullptr, box,
                             TRUE);
        get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &gnx, &index, &grpname);
    }
    else
    {
        rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &gnx, &index, &grpname);
    }

    if (bMol)
    {
        if (!bTop)
        {
            gmx_fatal(FARGS, "Need a topology to determine the molecules");
        }
        snew(normm, top.atoms.nr);
        precalc(top, normm);
        index_atom2mol(&gnx, index, &top.mols);
    }

    /* Correlation stuff */
    snew(c1, gnx);
    for (i = 0; (i < gnx); i++)
    {
        c1[i] = nullptr;
    }

    read_first_frame(oenv, &status, ftp2fn(efTRN, NFILE, fnm), &fr, TRX_NEED_V);
    t0 = fr.time;

    n_alloc = 0;
    counter = 0;
    do
    {
        if (counter >= n_alloc)
        {
            n_alloc += 100;
            for (i = 0; i < gnx; i++)
            {
                srenew(c1[i], DIM*n_alloc);
            }
        }
        counter_dim = DIM*counter;
        if (bMol)
        {
            for (i = 0; i < gnx; i++)
            {
                clear_rvec(mv_mol);
                k = top.mols.index[index[i]];
                l = top.mols.index[index[i]+1];
                for (j = k; j < l; j++)
                {
                    if (bMass)
                    {
                        mass = top.atoms.atom[j].m;
                    }
                    else
                    {
                        mass = normm[j];
                    }
                    mv_mol[XX] += mass*fr.v[j][XX];
                    mv_mol[YY] += mass*fr.v[j][YY];
                    mv_mol[ZZ] += mass*fr.v[j][ZZ];
                }
                c1[i][counter_dim+XX] = mv_mol[XX];
                c1[i][counter_dim+YY] = mv_mol[YY];
                c1[i][counter_dim+ZZ] = mv_mol[ZZ];
            }
        }
        else
        {
            for (i = 0; i < gnx; i++)
            {
                if (bMass)
                {
                    mass = top.atoms.atom[index[i]].m;
                }
                else
                {
                    mass = 1;
                }
                c1[i][counter_dim+XX] = mass*fr.v[index[i]][XX];
                c1[i][counter_dim+YY] = mass*fr.v[index[i]][YY];
                c1[i][counter_dim+ZZ] = mass*fr.v[index[i]][ZZ];
            }
        }

        t1 = fr.time;

        counter++;
    }
    while (read_next_frame(oenv, status, &fr));

    close_trx(status);

    if (counter >= 4)
    {
        /* Compute time step between frames */
        dt = (t1-t0)/(counter-1);
        do_autocorr(opt2fn("-o", NFILE, fnm), oenv,
                    bMass ?
                    "Momentum Autocorrelation Function" :
                    "Velocity Autocorrelation Function",
                    counter, gnx, c1, dt, eacVector, TRUE);

        do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy");

        if (opt2bSet("-os", NFILE, fnm))
        {
            calc_spectrum(counter/2, (real *) (c1[0]), (t1-t0)/2, opt2fn("-os", NFILE, fnm),
                          oenv, bRecip);
            do_view(oenv, opt2fn("-os", NFILE, fnm), "-nxy");
        }
    }
    else
    {
        fprintf(stderr, "Not enough frames in trajectory - no output generated.\n");
    }

    return 0;
}
