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

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/eigio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

enum class PbcType : int;
struct gmx_output_env_t;

int gmx_nmtraj(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] generates an virtual trajectory from an eigenvector, ",
        "corresponding to a harmonic Cartesian oscillation around the average ",
        "structure. The eigenvectors should normally be mass-weighted, but you can ",
        "use non-weighted eigenvectors to generate orthogonal motions. ",
        "The output frames are written as a trajectory file covering an entire period, and ",
        "the first frame is the average structure. If you write the trajectory in (or convert to) ",
        "PDB format you can view it directly in PyMol and also render a photorealistic movie. ",
        "Motion amplitudes are calculated from the eigenvalues and a preset temperature, ",
        "assuming equipartition of the energy over all modes. To make the motion clearly visible ",
        "in PyMol you might want to amplify it by setting an unrealistically high temperature. ",
        "However, be aware that both the linear Cartesian displacements and mass weighting will ",
        "lead to serious structure deformation for high amplitudes - this is is simply a ",
        "limitation of the Cartesian normal mode model. By default the selected eigenvector ",
        "is set to 7, since ",
        "the first six normal modes are the translational and rotational degrees of freedom."
    };

    static real        refamplitude = 0.25;
    static int         nframes      = 30;
    static real        temp         = 300.0;
    static const char* eignrvec     = "7";
    static const char* phasevec     = "0.0";

    t_pargs pa[] = {
        { "-eignr", FALSE, etSTR, { &eignrvec }, "String of eigenvectors to use (first is 1)" },
        { "-phases", FALSE, etSTR, { &phasevec }, "String of phases (default is 0.0)" },
        { "-temp", FALSE, etREAL, { &temp }, "Temperature (K)" },
        { "-amplitude", FALSE, etREAL, { &refamplitude }, "Amplitude for modes with eigenvalue<=0" },
        { "-nframes", FALSE, etINT, { &nframes }, "Number of frames to generate" }
    };

#define NPA asize(pa)

    t_trxstatus* out;
    t_topology   top;
    PbcType      pbcType;
    t_atoms*     atoms;
    rvec *       xtop, *xref, *xav, *xout;
    int          nvec, *eignr = nullptr;
    rvec**       eigvec = nullptr;
    matrix       box;
    int          natoms;
    int          i, j, k, kmode, d;
    gmx_bool     bDMR, bDMA, bFit;

    real*             eigval;
    int*              dummy;
    real*             invsqrtm;
    int*              out_eigidx;
    rvec*             this_eigvec;
    real              omega, Ekin, m, vel;
    real*             amplitude;
    gmx_output_env_t* oenv;

    t_filenm fnm[] = { { efTPS, nullptr, nullptr, ffREAD },
                       { efTRN, "-v", "eigenvec", ffREAD },
                       { efTRO, "-o", "nmtraj", ffWRITE } };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    read_eigenvectors(
            opt2fn("-v", NFILE, fnm), &natoms, &bFit, &xref, &bDMR, &xav, &bDMA, &nvec, &eignr, &eigvec, &eigval);

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, box, bDMA);

    /* Find vectors and phases */

    std::vector<int> imodes;
    for (const auto& imodeString : gmx::splitString(eignrvec))
    {
        imodes.emplace_back(gmx::fromStdString<int>(imodeString));
    }
    int nmodes = gmx::ssize(imodes);

    std::vector<real> phases;
    phases.reserve(nmodes);
    for (const auto& phaseString : gmx::splitString(phasevec))
    {
        phases.emplace_back(gmx::fromStdString<real>(phaseString));
    }
    int nphases = gmx::ssize(phases);

    if (nphases > nmodes)
    {
        gmx_fatal(FARGS, "More phases than eigenvector indices specified.\n");
    }

    if (nmodes > nphases)
    {
        printf("Warning: Setting phase of last %d modes to zero...\n", nmodes - nphases);
        phases.resize(nmodes, 0);
    }

    atoms = &top.atoms;

    if (atoms->nr != natoms)
    {
        gmx_fatal(FARGS, "Different number of atoms in topology and eigenvectors.\n");
    }

    snew(dummy, natoms);
    for (i = 0; i < natoms; i++)
    {
        dummy[i] = i;
    }

    /* Find the eigenvalue/vector to match our select one */
    snew(out_eigidx, nmodes);
    for (i = 0; i < nmodes; i++)
    {
        out_eigidx[i] = -1;
    }

    for (i = 0; i < nvec; i++)
    {
        for (j = 0; j < nmodes; j++)
        {
            if (imodes[j] == eignr[i])
            {
                out_eigidx[j] = i;
            }
        }
    }
    for (i = 0; i < nmodes; i++)
    {
        if (out_eigidx[i] == -1)
        {
            gmx_fatal(FARGS, "Could not find mode %d in eigenvector file.\n", imodes[i]);
        }
    }


    snew(invsqrtm, natoms);

    if (bDMA)
    {
        for (i = 0; (i < natoms); i++)
        {
            invsqrtm[i] = gmx::invsqrt(atoms->atom[i].m);
        }
    }
    else
    {
        for (i = 0; (i < natoms); i++)
        {
            invsqrtm[i] = 1.0;
        }
    }

    snew(xout, natoms);
    snew(amplitude, nmodes);

    printf("mode phases: %g %g\n", phases[0], phases[1]);

    for (i = 0; i < nmodes; i++)
    {
        kmode       = out_eigidx[i];
        this_eigvec = eigvec[kmode];

        if ((kmode >= 6) && (eigval[kmode] > 0))
        {
            /* Derive amplitude from temperature and eigenvalue if we can */

            /* Convert eigenvalue to angular frequency, in units s^(-1) */
            omega = std::sqrt(eigval[kmode] * 1.0E21 / (gmx::c_avogadro * gmx::c_amu));
            /* Harmonic motion will be x=x0 + A*sin(omega*t)*eigenvec.
             * The velocity is thus:
             *
             * v = A*omega*cos(omega*t)*eigenvec.
             *
             * And the average kinetic energy the integral of mass*v*v/2 over a
             * period:
             *
             * (1/4)*mass*A*omega*eigenvec
             *
             * For t =2*pi*n, all energy will be kinetic, and v=A*omega*eigenvec.
             * The kinetic energy will be sum(0.5*mass*v*v) if we temporarily set A to 1,
             * and the average over a period half of this.
             */

            Ekin = 0;
            for (k = 0; k < natoms; k++)
            {
                m = atoms->atom[k].m;
                for (d = 0; d < DIM; d++)
                {
                    vel = omega * this_eigvec[k][d];
                    Ekin += 0.5 * 0.5 * m * vel * vel;
                }
            }

            /* Convert Ekin from amu*(nm/s)^2 to J, i.e., kg*(m/s)^2
             * This will also be proportional to A^2
             */
            Ekin *= gmx::c_amu * 1E-18;

            /* Set the amplitude so the energy is kT/2 */
            amplitude[i] = std::sqrt(0.5 * gmx::c_boltzmann * temp / Ekin);
        }
        else
        {
            amplitude[i] = refamplitude;
        }
    }

    out = open_trx(ftp2fn(efTRO, NFILE, fnm), "w");

    /* Write a sine oscillation around the average structure,
     * modulated by the eigenvector with selected amplitude.
     */

    for (i = 0; i < nframes; i++)
    {
        real fraction = static_cast<real>(i) / static_cast<real>(nframes);
        for (j = 0; j < natoms; j++)
        {
            copy_rvec(xav[j], xout[j]);
        }

        for (k = 0; k < nmodes; k++)
        {
            kmode       = out_eigidx[k];
            this_eigvec = eigvec[kmode];

            for (j = 0; j < natoms; j++)
            {
                for (d = 0; d < DIM; d++)
                {
                    xout[j][d] += amplitude[k] * std::sin(2 * M_PI * (fraction + phases[k] / 360.0))
                                  * this_eigvec[j][d];
                }
            }
        }
        write_trx(out, natoms, dummy, atoms, i, fraction, box, xout, nullptr, nullptr);
    }

    fprintf(stderr, "\n");
    close_trx(out);

    return 0;
}
