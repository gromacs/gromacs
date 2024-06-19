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

#include <cmath>
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <string>

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
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;
struct gmx_output_env_t;


int gmx_nmens(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] generates an ensemble around an average structure",
        "in a subspace that is defined by a set of normal modes (eigenvectors).",
        "The eigenvectors are assumed to be mass-weighted.",
        "The position along each eigenvector is randomly taken from a Gaussian",
        "distribution with variance kT/eigenvalue.[PAR]",
        "By default the starting eigenvector is set to 7, since the first six",
        "normal modes are the translational and rotational degrees of freedom."
    };
    static int  nstruct = 100, first = 7, last = -1, seed = 0;
    static real temp = 300.0;
    t_pargs     pa[] = {
        { "-temp", FALSE, etREAL, { &temp }, "Temperature in Kelvin" },
        { "-seed", FALSE, etINT, { &seed }, "Random seed (0 means generate)" },
        { "-num", FALSE, etINT, { &nstruct }, "Number of structures to generate" },
        { "-first", FALSE, etINT, { &first }, "First eigenvector to use (-1 is select)" },
        { "-last", FALSE, etINT, { &last }, "Last eigenvector to use (-1 is till the last)" }
    };
#define NPA asize(pa)

    t_trxstatus*        out;
    t_topology          top;
    PbcType             pbcType;
    t_atoms*            atoms;
    rvec *              xtop, *xref, *xav, *xout1, *xout2;
    gmx_bool            bDMR, bDMA, bFit;
    int                 nvec, *eignr = nullptr;
    rvec**              eigvec = nullptr;
    matrix              box;
    real *              eigval, *invsqrtm, t, disp;
    int                 natoms;
    char*               grpname;
    const char*         indexfile;
    int                 i, j, d, s, v;
    int                 nout, *iout, noutvec, *outvec;
    int*                index;
    real                rfac, rhalf, jr;
    gmx_output_env_t*   oenv;
    int                 jran;
    const unsigned long im = 0xffff;
    const unsigned long ia = 1093;
    const unsigned long ic = 18257;


    t_filenm fnm[] = { { efTRN, "-v", "eigenvec", ffREAD },
                       { efXVG, "-e", "eigenval", ffREAD },
                       { efTPS, nullptr, nullptr, ffREAD },
                       { efNDX, nullptr, nullptr, ffOPTRD },
                       { efTRO, "-o", "ensemble", ffWRITE } };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    indexfile = ftp2fn_null(efNDX, NFILE, fnm);

    read_eigenvectors(
            opt2fn("-v", NFILE, fnm), &natoms, &bFit, &xref, &bDMR, &xav, &bDMA, &nvec, &eignr, &eigvec, &eigval);

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &pbcType, &xtop, nullptr, box, bDMA);
    atoms = &top.atoms;

    printf("\nSelect an index group of %d elements that corresponds to the eigenvectors\n", natoms);
    get_index(atoms, indexfile, 1, &i, &index, &grpname);
    if (i != natoms)
    {
        gmx_fatal(FARGS, "you selected a group with %d elements instead of %d", i, natoms);
    }
    printf("\n");

    snew(invsqrtm, natoms);
    if (bDMA)
    {
        for (i = 0; (i < natoms); i++)
        {
            invsqrtm[i] = gmx::invsqrt(atoms->atom[index[i]].m);
        }
    }
    else
    {
        for (i = 0; (i < natoms); i++)
        {
            invsqrtm[i] = 1.0;
        }
    }

    if (last == -1)
    {
        last = natoms * DIM;
    }
    if (first > -1)
    {
        /* make an index from first to last */
        nout = last - first + 1;
        snew(iout, nout);
        for (i = 0; i < nout; i++)
        {
            iout[i] = first - 1 + i;
        }
    }
    else
    {
        printf("Select eigenvectors for output, end your selection with 0\n");
        nout = -1;
        iout = nullptr;
        do
        {
            nout++;
            srenew(iout, nout + 1);
            if (1 != scanf("%d", &iout[nout]))
            {
                gmx_fatal(FARGS, "Error reading user input");
            }
            iout[nout]--;
        } while (iout[nout] >= 0);
        printf("\n");
    }

    /* make an index of the eigenvectors which are present */
    snew(outvec, nout);
    noutvec = 0;
    for (i = 0; i < nout; i++)
    {
        j = 0;
        while ((j < nvec) && (eignr[j] != iout[i]))
        {
            j++;
        }
        if ((j < nvec) && (eignr[j] == iout[i]))
        {
            outvec[noutvec] = j;
            iout[noutvec]   = iout[i];
            noutvec++;
        }
    }

    fprintf(stderr, "%d eigenvectors selected for output\n", noutvec);


    if (seed == 0)
    {
        // Make do with 32 bits for now to avoid changing user input to hex
        seed = static_cast<int>(gmx::makeRandomSeed());
    }

    gmx::DefaultRandomEngine rng(seed);

    fprintf(stderr, "Using random seed %d and a temperature of %g K.\n", seed, temp);

    gmx::UniformIntDistribution<int> dist(0, im - 1);
    jran = dist(rng);

    snew(xout1, natoms);
    snew(xout2, atoms->nr);
    out = open_trx(ftp2fn(efTRO, NFILE, fnm), "w");

    for (s = 0; s < nstruct; s++)
    {
        for (i = 0; i < natoms; i++)
        {
            copy_rvec(xav[i], xout1[i]);
        }
        for (j = 0; j < noutvec; j++)
        {
            v = outvec[j];
            /* (r-0.5) n times:  var_n = n * var_1 = n/12
               n=4:  var_n = 1/3, so multiply with 3 */

            rfac  = std::sqrt(3.0 * gmx::c_boltz * temp / eigval[iout[j]]);
            rhalf = 2.0 * rfac;
            rfac  = rfac / im;

            jran = (jran * ia + ic) & im;
            jr   = jran;
            jran = (jran * ia + ic) & im;
            jr += jran;
            jran = (jran * ia + ic) & im;
            jr += jran;
            jran = (jran * ia + ic) & im;
            jr += jran;
            disp = rfac * jr - rhalf;

            for (i = 0; i < natoms; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    xout1[i][d] += disp * eigvec[v][i][d] * invsqrtm[i];
                }
            }
        }
        for (i = 0; i < natoms; i++)
        {
            copy_rvec(xout1[i], xout2[index[i]]);
        }
        t = s + 1;
        write_trx(out, natoms, index, atoms, 0, t, box, xout2, nullptr, nullptr);
        fprintf(stderr, "\rGenerated %d structures", s + 1);
        fflush(stderr);
    }
    fprintf(stderr, "\n");
    close_trx(out);

    return 0;
}
