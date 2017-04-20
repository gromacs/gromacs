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
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/eigio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"

int gmx_covar(int argc, char *argv[])
{
    const char       *desc[] = {
        "[THISMODULE] calculates and diagonalizes the (mass-weighted)",
        "covariance matrix.",
        "All structures are fitted to the structure in the structure file.",
        "When this is not a run input file periodicity will not be taken into",
        "account. When the fit and analysis groups are identical and the analysis",
        "is non mass-weighted, the fit will also be non mass-weighted.",
        "[PAR]",
        "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
        "When the same atoms are used for the fit and the covariance analysis,",
        "the reference structure for the fit is written first with t=-1.",
        "The average (or reference when [TT]-ref[tt] is used) structure is",
        "written with t=0, the eigenvectors",
        "are written as frames with the eigenvector number and eigenvalue",
        "as step number and timestamp, respectively.",
        "[PAR]",
        "The eigenvectors can be analyzed with [gmx-anaeig].",
        "[PAR]",
        "Option [TT]-ascii[tt] writes the whole covariance matrix to",
        "an ASCII file. The order of the elements is: x1x1, x1y1, x1z1, x1x2, ...",
        "[PAR]",
        "Option [TT]-xpm[tt] writes the whole covariance matrix to an [REF].xpm[ref] file.",
        "[PAR]",
        "Option [TT]-xpma[tt] writes the atomic covariance matrix to an [REF].xpm[ref] file,",
        "i.e. for each atom pair the sum of the xx, yy and zz covariances is",
        "written.",
        "[PAR]",
        "Note that the diagonalization of a matrix requires memory and time",
        "that will increase at least as fast as than the square of the number",
        "of atoms involved. It is easy to run out of memory, in which",
        "case this tool will probably exit with a 'Segmentation fault'. You",
        "should consider carefully whether a reduced set of atoms will meet",
        "your needs for lower costs."
    };
    static gmx_bool   bFit = TRUE, bRef = FALSE, bM = FALSE, bPBC = TRUE;
    static int        end  = -1;
    t_pargs           pa[] = {
        { "-fit",  FALSE, etBOOL, {&bFit},
          "Fit to a reference structure"},
        { "-ref",  FALSE, etBOOL, {&bRef},
          "Use the deviation from the conformation in the structure file instead of from the average" },
        { "-mwa",  FALSE, etBOOL, {&bM},
          "Mass-weighted covariance analysis"},
        { "-last",  FALSE, etINT, {&end},
          "Last eigenvector to write away (-1 is till the last)" },
        { "-pbc",  FALSE,  etBOOL, {&bPBC},
          "Apply corrections for periodic boundary conditions" }
    };
    FILE             *out = nullptr; /* initialization makes all compilers happy */
    t_trxstatus      *status;
    t_topology        top;
    int               ePBC;
    t_atoms          *atoms;
    rvec             *x, *xread, *xref, *xav, *xproj;
    matrix            box, zerobox;
    real             *sqrtm, *mat, *eigenvalues, sum, trace, inv_nframes;
    real              t, tstart, tend, **mat2;
    real              xj, *w_rls = nullptr;
    real              min, max, *axis;
    int               natoms, nat, nframes0, nframes, nlevels;
    gmx_int64_t       ndim, i, j, k, l;
    int               WriteXref;
    const char       *fitfile, *trxfile, *ndxfile;
    const char       *eigvalfile, *eigvecfile, *averfile, *logfile;
    const char       *asciifile, *xpmfile, *xpmafile;
    char              str[STRLEN], *fitname, *ananame;
    int               d, dj, nfit;
    int              *index, *ifit;
    gmx_bool          bDiffMass1, bDiffMass2;
    char              timebuf[STRLEN];
    t_rgb             rlo, rmi, rhi;
    real             *eigenvectors;
    gmx_output_env_t *oenv;
    gmx_rmpbc_t       gpbc = nullptr;

    t_filenm          fnm[] = {
        { efTRX, "-f",  nullptr, ffREAD },
        { efTPS, nullptr,  nullptr, ffREAD },
        { efNDX, nullptr,  nullptr, ffOPTRD },
        { efXVG, nullptr,  "eigenval", ffWRITE },
        { efTRN, "-v",  "eigenvec", ffWRITE },
        { efSTO, "-av", "average.pdb", ffWRITE },
        { efLOG, nullptr,  "covar", ffWRITE },
        { efDAT, "-ascii", "covar", ffOPTWR },
        { efXPM, "-xpm", "covar", ffOPTWR },
        { efXPM, "-xpma", "covara", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    clear_mat(zerobox);

    fitfile    = ftp2fn(efTPS, NFILE, fnm);
    trxfile    = ftp2fn(efTRX, NFILE, fnm);
    ndxfile    = ftp2fn_null(efNDX, NFILE, fnm);
    eigvalfile = ftp2fn(efXVG, NFILE, fnm);
    eigvecfile = ftp2fn(efTRN, NFILE, fnm);
    averfile   = ftp2fn(efSTO, NFILE, fnm);
    logfile    = ftp2fn(efLOG, NFILE, fnm);
    asciifile  = opt2fn_null("-ascii", NFILE, fnm);
    xpmfile    = opt2fn_null("-xpm", NFILE, fnm);
    xpmafile   = opt2fn_null("-xpma", NFILE, fnm);

    read_tps_conf(fitfile, &top, &ePBC, &xref, nullptr, box, TRUE);
    atoms = &top.atoms;

    if (bFit)
    {
        printf("\nChoose a group for the least squares fit\n");
        get_index(atoms, ndxfile, 1, &nfit, &ifit, &fitname);
        if (nfit < 3)
        {
            gmx_fatal(FARGS, "Need >= 3 points to fit!\n");
        }
    }
    else
    {
        nfit = 0;
    }
    printf("\nChoose a group for the covariance analysis\n");
    get_index(atoms, ndxfile, 1, &natoms, &index, &ananame);

    bDiffMass1 = FALSE;
    if (bFit)
    {
        snew(w_rls, atoms->nr);
        for (i = 0; (i < nfit); i++)
        {
            w_rls[ifit[i]] = atoms->atom[ifit[i]].m;
            if (i)
            {
                bDiffMass1 = bDiffMass1 || (w_rls[ifit[i]] != w_rls[ifit[i-1]]);
            }
        }
    }
    bDiffMass2 = FALSE;
    snew(sqrtm, natoms);
    for (i = 0; (i < natoms); i++)
    {
        if (bM)
        {
            sqrtm[i] = std::sqrt(atoms->atom[index[i]].m);
            if (i)
            {
                bDiffMass2 = bDiffMass2 || (sqrtm[i] != sqrtm[i-1]);
            }
        }
        else
        {
            sqrtm[i] = 1.0;
        }
    }

    if (bFit && bDiffMass1 && !bDiffMass2)
    {
        bDiffMass1 = natoms != nfit;
        for (i = 0; (i < natoms) && !bDiffMass1; i++)
        {
            bDiffMass1 = index[i] != ifit[i];
        }
        if (!bDiffMass1)
        {
            fprintf(stderr, "\n"
                    "Note: the fit and analysis group are identical,\n"
                    "      while the fit is mass weighted and the analysis is not.\n"
                    "      Making the fit non mass weighted.\n\n");
            for (i = 0; (i < nfit); i++)
            {
                w_rls[ifit[i]] = 1.0;
            }
        }
    }

    /* Prepare reference frame */
    if (bPBC)
    {
        gpbc = gmx_rmpbc_init(&top.idef, ePBC, atoms->nr);
        gmx_rmpbc(gpbc, atoms->nr, box, xref);
    }
    if (bFit)
    {
        reset_x(nfit, ifit, atoms->nr, nullptr, xref, w_rls);
    }

    snew(x, natoms);
    snew(xav, natoms);
    ndim = natoms*DIM;
    if (std::sqrt(static_cast<real>(GMX_INT64_MAX)) < static_cast<real>(ndim))
    {
        gmx_fatal(FARGS, "Number of degrees of freedoms to large for matrix.\n");
    }
    snew(mat, ndim*ndim);

    fprintf(stderr, "Calculating the average structure ...\n");
    nframes0 = 0;
    nat      = read_first_x(oenv, &status, trxfile, &t, &xread, box);
    if (nat != atoms->nr)
    {
        fprintf(stderr, "\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n", natoms, nat);
    }
    do
    {
        nframes0++;
        /* calculate x: a fitted struture of the selected atoms */
        if (bPBC)
        {
            gmx_rmpbc(gpbc, nat, box, xread);
        }
        if (bFit)
        {
            reset_x(nfit, ifit, nat, nullptr, xread, w_rls);
            do_fit(nat, w_rls, xref, xread);
        }
        for (i = 0; i < natoms; i++)
        {
            rvec_inc(xav[i], xread[index[i]]);
        }
    }
    while (read_next_x(oenv, status, &t, xread, box));
    close_trj(status);

    inv_nframes = 1.0/nframes0;
    for (i = 0; i < natoms; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            xav[i][d]         *= inv_nframes;
            xread[index[i]][d] = xav[i][d];
        }
    }
    write_sto_conf_indexed(opt2fn("-av", NFILE, fnm), "Average structure",
                           atoms, xread, nullptr, epbcNONE, zerobox, natoms, index);
    sfree(xread);

    fprintf(stderr, "Constructing covariance matrix (%dx%d) ...\n", static_cast<int>(ndim), static_cast<int>(ndim));
    nframes = 0;
    nat     = read_first_x(oenv, &status, trxfile, &t, &xread, box);
    tstart  = t;
    do
    {
        nframes++;
        tend = t;
        /* calculate x: a (fitted) structure of the selected atoms */
        if (bPBC)
        {
            gmx_rmpbc(gpbc, nat, box, xread);
        }
        if (bFit)
        {
            reset_x(nfit, ifit, nat, nullptr, xread, w_rls);
            do_fit(nat, w_rls, xref, xread);
        }
        if (bRef)
        {
            for (i = 0; i < natoms; i++)
            {
                rvec_sub(xread[index[i]], xref[index[i]], x[i]);
            }
        }
        else
        {
            for (i = 0; i < natoms; i++)
            {
                rvec_sub(xread[index[i]], xav[i], x[i]);
            }
        }

        for (j = 0; j < natoms; j++)
        {
            for (dj = 0; dj < DIM; dj++)
            {
                k  = ndim*(DIM*j+dj);
                xj = x[j][dj];
                for (i = j; i < natoms; i++)
                {
                    l = k+DIM*i;
                    for (d = 0; d < DIM; d++)
                    {
                        mat[l+d] += x[i][d]*xj;
                    }
                }
            }
        }
    }
    while (read_next_x(oenv, status, &t, xread, box) &&
           (bRef || nframes < nframes0));
    close_trj(status);
    gmx_rmpbc_done(gpbc);

    fprintf(stderr, "Read %d frames\n", nframes);

    if (bRef)
    {
        /* copy the reference structure to the ouput array x */
        snew(xproj, natoms);
        for (i = 0; i < natoms; i++)
        {
            copy_rvec(xref[index[i]], xproj[i]);
        }
    }
    else
    {
        xproj = xav;
    }

    /* correct the covariance matrix for the mass */
    inv_nframes = 1.0/nframes;
    for (j = 0; j < natoms; j++)
    {
        for (dj = 0; dj < DIM; dj++)
        {
            for (i = j; i < natoms; i++)
            {
                k = ndim*(DIM*j+dj)+DIM*i;
                for (d = 0; d < DIM; d++)
                {
                    mat[k+d] = mat[k+d]*inv_nframes*sqrtm[i]*sqrtm[j];
                }
            }
        }
    }

    /* symmetrize the matrix */
    for (j = 0; j < ndim; j++)
    {
        for (i = j; i < ndim; i++)
        {
            mat[ndim*i+j] = mat[ndim*j+i];
        }
    }

    trace = 0;
    for (i = 0; i < ndim; i++)
    {
        trace += mat[i*ndim+i];
    }
    fprintf(stderr, "\nTrace of the covariance matrix: %g (%snm^2)\n",
            trace, bM ? "u " : "");

    if (asciifile)
    {
        out = gmx_ffopen(asciifile, "w");
        for (j = 0; j < ndim; j++)
        {
            for (i = 0; i < ndim; i += 3)
            {
                fprintf(out, "%g %g %g\n",
                        mat[ndim*j+i], mat[ndim*j+i+1], mat[ndim*j+i+2]);
            }
        }
        gmx_ffclose(out);
    }

    if (xpmfile)
    {
        min = 0;
        max = 0;
        snew(mat2, ndim);
        for (j = 0; j < ndim; j++)
        {
            mat2[j] = &(mat[ndim*j]);
            for (i = 0; i <= j; i++)
            {
                if (mat2[j][i] < min)
                {
                    min = mat2[j][i];
                }
                if (mat2[j][j] > max)
                {
                    max = mat2[j][i];
                }
            }
        }
        snew(axis, ndim);
        for (i = 0; i < ndim; i++)
        {
            axis[i] = i+1;
        }
        rlo.r   = 0; rlo.g = 0; rlo.b = 1;
        rmi.r   = 1; rmi.g = 1; rmi.b = 1;
        rhi.r   = 1; rhi.g = 0; rhi.b = 0;
        out     = gmx_ffopen(xpmfile, "w");
        nlevels = 80;
        write_xpm3(out, 0, "Covariance", bM ? "u nm^2" : "nm^2",
                   "dim", "dim", ndim, ndim, axis, axis,
                   mat2, min, 0.0, max, rlo, rmi, rhi, &nlevels);
        gmx_ffclose(out);
        sfree(axis);
        sfree(mat2);
    }

    if (xpmafile)
    {
        min = 0;
        max = 0;
        snew(mat2, ndim/DIM);
        for (i = 0; i < ndim/DIM; i++)
        {
            snew(mat2[i], ndim/DIM);
        }
        for (j = 0; j < ndim/DIM; j++)
        {
            for (i = 0; i <= j; i++)
            {
                mat2[j][i] = 0;
                for (d = 0; d < DIM; d++)
                {
                    mat2[j][i] += mat[ndim*(DIM*j+d)+DIM*i+d];
                }
                if (mat2[j][i] < min)
                {
                    min = mat2[j][i];
                }
                if (mat2[j][j] > max)
                {
                    max = mat2[j][i];
                }
                mat2[i][j] = mat2[j][i];
            }
        }
        snew(axis, ndim/DIM);
        for (i = 0; i < ndim/DIM; i++)
        {
            axis[i] = i+1;
        }
        rlo.r   = 0; rlo.g = 0; rlo.b = 1;
        rmi.r   = 1; rmi.g = 1; rmi.b = 1;
        rhi.r   = 1; rhi.g = 0; rhi.b = 0;
        out     = gmx_ffopen(xpmafile, "w");
        nlevels = 80;
        write_xpm3(out, 0, "Covariance", bM ? "u nm^2" : "nm^2",
                   "atom", "atom", ndim/DIM, ndim/DIM, axis, axis,
                   mat2, min, 0.0, max, rlo, rmi, rhi, &nlevels);
        gmx_ffclose(out);
        sfree(axis);
        for (i = 0; i < ndim/DIM; i++)
        {
            sfree(mat2[i]);
        }
        sfree(mat2);
    }


    /* call diagonalization routine */

    snew(eigenvalues, ndim);
    snew(eigenvectors, ndim*ndim);

    std::memcpy(eigenvectors, mat, ndim*ndim*sizeof(real));
    fprintf(stderr, "\nDiagonalizing ...\n");
    fflush(stderr);
    eigensolver(eigenvectors, ndim, 0, ndim, eigenvalues, mat);
    sfree(eigenvectors);

    /* now write the output */

    sum = 0;
    for (i = 0; i < ndim; i++)
    {
        sum += eigenvalues[i];
    }
    fprintf(stderr, "\nSum of the eigenvalues: %g (%snm^2)\n",
            sum, bM ? "u " : "");
    if (std::abs(trace-sum) > 0.01*trace)
    {
        fprintf(stderr, "\nWARNING: eigenvalue sum deviates from the trace of the covariance matrix\n");
    }

    /* Set 'end', the maximum eigenvector and -value index used for output */
    if (end == -1)
    {
        if (nframes-1 < ndim)
        {
            end = nframes-1;
            fprintf(stderr, "\nWARNING: there are fewer frames in your trajectory than there are\n");
            fprintf(stderr, "degrees of freedom in your system. Only generating the first\n");
            fprintf(stderr, "%d out of %d eigenvectors and eigenvalues.\n", end, static_cast<int>(ndim));
        }
        else
        {
            end = ndim;
        }
    }

    fprintf(stderr, "\nWriting eigenvalues to %s\n", eigvalfile);

    sprintf(str, "(%snm\\S2\\N)", bM ? "u " : "");
    out = xvgropen(eigvalfile,
                   "Eigenvalues of the covariance matrix",
                   "Eigenvector index", str, oenv);
    for (i = 0; (i < end); i++)
    {
        fprintf (out, "%10d %g\n", static_cast<int>(i+1), eigenvalues[ndim-1-i]);
    }
    xvgrclose(out);

    if (bFit)
    {
        /* misuse lambda: 0/1 mass weighted analysis no/yes */
        if (nfit == natoms)
        {
            WriteXref = eWXR_YES;
            for (i = 0; i < nfit; i++)
            {
                copy_rvec(xref[ifit[i]], x[i]);
            }
        }
        else
        {
            WriteXref = eWXR_NO;
        }
    }
    else
    {
        /* misuse lambda: -1 for no fit */
        WriteXref = eWXR_NOFIT;
    }

    write_eigenvectors(eigvecfile, natoms, mat, TRUE, 1, end,
                       WriteXref, x, bDiffMass1, xproj, bM, eigenvalues);

    out = gmx_ffopen(logfile, "w");

    gmx_format_current_time(timebuf, STRLEN);
    fprintf(out, "Covariance analysis log, written %s\n", timebuf);

    fprintf(out, "Program: %s\n", argv[0]);
    gmx_getcwd(str, STRLEN);

    fprintf(out, "Working directory: %s\n\n", str);

    fprintf(out, "Read %d frames from %s (time %g to %g %s)\n", nframes, trxfile,
            output_env_conv_time(oenv, tstart), output_env_conv_time(oenv, tend), output_env_get_time_unit(oenv).c_str());
    if (bFit)
    {
        fprintf(out, "Read reference structure for fit from %s\n", fitfile);
    }
    if (ndxfile)
    {
        fprintf(out, "Read index groups from %s\n", ndxfile);
    }
    fprintf(out, "\n");

    fprintf(out, "Analysis group is '%s' (%d atoms)\n", ananame, natoms);
    if (bFit)
    {
        fprintf(out, "Fit group is '%s' (%d atoms)\n", fitname, nfit);
    }
    else
    {
        fprintf(out, "No fit was used\n");
    }
    fprintf(out, "Analysis is %smass weighted\n", bDiffMass2 ? "" : "non-");
    if (bFit)
    {
        fprintf(out, "Fit is %smass weighted\n", bDiffMass1 ? "" : "non-");
    }
    fprintf(out, "Diagonalized the %dx%d covariance matrix\n", static_cast<int>(ndim), static_cast<int>(ndim));
    fprintf(out, "Trace of the covariance matrix before diagonalizing: %g\n",
            trace);
    fprintf(out, "Trace of the covariance matrix after diagonalizing: %g\n\n",
            sum);

    fprintf(out, "Wrote %d eigenvalues to %s\n", static_cast<int>(end), eigvalfile);
    if (WriteXref == eWXR_YES)
    {
        fprintf(out, "Wrote reference structure to %s\n", eigvecfile);
    }
    fprintf(out, "Wrote average structure to %s and %s\n", averfile, eigvecfile);
    fprintf(out, "Wrote eigenvectors %d to %d to %s\n", 1, end, eigvecfile);

    gmx_ffclose(out);

    fprintf(stderr, "Wrote the log to %s\n", logfile);

    return 0;
}
