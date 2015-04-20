/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <math.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/eigio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static double cv_corr(double nu, double T)
{
    double x  = PLANCK*nu/(BOLTZ*T);
    double ex = exp(x);

    if (nu <= 0)
    {
        return BOLTZ*KILO;
    }
    else
    {
        return BOLTZ*KILO*(ex*sqr(x)/sqr(ex-1) - 1);
    }
}

static double u_corr(double nu, double T)
{
    double x  = PLANCK*nu/(BOLTZ*T);
    double ex = exp(x);

    if (nu <= 0)
    {
        return BOLTZ*T;
    }
    else
    {
        return BOLTZ*T*(0.5*x - 1 + x/(ex-1));
    }
}

static int get_nharm_mt(gmx_moltype_t *mt)
{
    static int harm_func[] = { F_BONDS };
    int        i, ft, nh;

    nh = 0;
    for (i = 0; (i < asize(harm_func)); i++)
    {
        ft  = harm_func[i];
        nh += mt->ilist[ft].nr/(interaction_function[ft].nratoms+1);
    }
    return nh;
}

static int get_nvsite_mt(gmx_moltype_t *mt)
{
    static int vs_func[] = {
        F_VSITE2, F_VSITE3, F_VSITE3FD, F_VSITE3FAD,
        F_VSITE3OUT, F_VSITE4FD, F_VSITE4FDN, F_VSITEN
    };
    int        i, ft, nh;

    nh = 0;
    for (i = 0; (i < asize(vs_func)); i++)
    {
        ft  = vs_func[i];
        nh += mt->ilist[ft].nr/(interaction_function[ft].nratoms+1);
    }
    return nh;
}

static int get_nharm(gmx_mtop_t *mtop, int *nvsites)
{
    int j, mt, nh, nv;

    nh = 0;
    nv = 0;
    for (j = 0; (j < mtop->nmolblock); j++)
    {
        mt  = mtop->molblock[j].type;
        nh += mtop->molblock[j].nmol * get_nharm_mt(&(mtop->moltype[mt]));
        nv += mtop->molblock[j].nmol * get_nvsite_mt(&(mtop->moltype[mt]));
    }
    *nvsites = nv;
    return nh;
}

static void
nma_full_hessian(real     *           hess,
                 int                  ndim,
                 gmx_bool             bM,
                 t_topology     *     top,
                 int                  begin,
                 int                  end,
                 real     *           eigenvalues,
                 real     *           eigenvectors)
{
    int  i, j, k, l;
    real mass_fac, rdum;
    int  natoms;

    natoms = top->atoms.nr;

    /* divide elements hess[i][j] by sqrt(mas[i])*sqrt(mas[j]) when required */

    if (bM)
    {
        for (i = 0; (i < natoms); i++)
        {
            for (j = 0; (j < DIM); j++)
            {
                for (k = 0; (k < natoms); k++)
                {
                    mass_fac = gmx_invsqrt(top->atoms.atom[i].m*top->atoms.atom[k].m);
                    for (l = 0; (l < DIM); l++)
                    {
                        hess[(i*DIM+j)*ndim+k*DIM+l] *= mass_fac;
                    }
                }
            }
        }
    }

    /* call diagonalization routine. */

    fprintf(stderr, "\nDiagonalizing to find vectors %d through %d...\n", begin, end);
    fflush(stderr);

    eigensolver(hess, ndim, begin-1, end-1, eigenvalues, eigenvectors);

    /* And scale the output eigenvectors */
    if (bM && eigenvectors != NULL)
    {
        for (i = 0; i < (end-begin+1); i++)
        {
            for (j = 0; j < natoms; j++)
            {
                mass_fac = gmx_invsqrt(top->atoms.atom[j].m);
                for (k = 0; (k < DIM); k++)
                {
                    eigenvectors[i*ndim+j*DIM+k] *= mass_fac;
                }
            }
        }
    }
}



static void
nma_sparse_hessian(gmx_sparsematrix_t     *     sparse_hessian,
                   gmx_bool                     bM,
                   t_topology     *             top,
                   int                          neig,
                   real     *                   eigenvalues,
                   real     *                   eigenvectors)
{
    int  i, j, k;
    int  row, col;
    real mass_fac;
    int  iatom, katom;
    int  natoms;
    int  ndim;

    natoms = top->atoms.nr;
    ndim   = DIM*natoms;

    /* Cannot check symmetry since we only store half matrix */
    /* divide elements hess[i][j] by sqrt(mas[i])*sqrt(mas[j]) when required */

    if (bM)
    {
        for (iatom = 0; (iatom < natoms); iatom++)
        {
            for (j = 0; (j < DIM); j++)
            {
                row = DIM*iatom+j;
                for (k = 0; k < sparse_hessian->ndata[row]; k++)
                {
                    col      = sparse_hessian->data[row][k].col;
                    katom    = col/3;
                    mass_fac = gmx_invsqrt(top->atoms.atom[iatom].m*top->atoms.atom[katom].m);
                    sparse_hessian->data[row][k].value *= mass_fac;
                }
            }
        }
    }
    fprintf(stderr, "\nDiagonalizing to find eigenvectors 1 through %d...\n", neig);
    fflush(stderr);

    sparse_eigensolver(sparse_hessian, neig, eigenvalues, eigenvectors, 10000000);

    /* Scale output eigenvectors */
    if (bM && eigenvectors != NULL)
    {
        for (i = 0; i < neig; i++)
        {
            for (j = 0; j < natoms; j++)
            {
                mass_fac = gmx_invsqrt(top->atoms.atom[j].m);
                for (k = 0; (k < DIM); k++)
                {
                    eigenvectors[i*ndim+j*DIM+k] *= mass_fac;
                }
            }
        }
    }
}



int gmx_nmeig(int argc, char *argv[])
{
    const char            *desc[] = {
        "[THISMODULE] calculates the eigenvectors/values of a (Hessian) matrix,",
        "which can be calculated with [gmx-mdrun].",
        "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
        "The structure is written first with t=0. The eigenvectors",
        "are written as frames with the eigenvector number as timestamp.",
        "The eigenvectors can be analyzed with [gmx-anaeig].",
        "An ensemble of structures can be generated from the eigenvectors with",
        "[gmx-nmens]. When mass weighting is used, the generated eigenvectors",
        "will be scaled back to plain Cartesian coordinates before generating the",
        "output. In this case, they will no longer be exactly orthogonal in the",
        "standard Cartesian norm, but in the mass-weighted norm they would be.[PAR]",
        "This program can be optionally used to compute quantum corrections to heat capacity",
        "and enthalpy by providing an extra file argument [TT]-qcorr[tt]. See the GROMACS",
        "manual, Chapter 1, for details. The result includes subtracting a harmonic",
        "degree of freedom at the given temperature.",
        "The total correction is printed on the terminal screen.",
        "The recommended way of getting the corrections out is:[PAR]",
        "[TT]gmx nmeig -s topol.tpr -f nm.mtx -first 7 -last 10000 -T 300 -qc [-constr][tt][PAR]",
        "The [TT]-constr[tt] option should be used when bond constraints were used during the",
        "simulation [BB]for all the covalent bonds[bb]. If this is not the case, ",
        "you need to analyze the [TT]quant_corr.xvg[tt] file yourself.[PAR]",
        "To make things more flexible, the program can also take virtual sites into account",
        "when computing quantum corrections. When selecting [TT]-constr[tt] and",
        "[TT]-qc[tt], the [TT]-begin[tt] and [TT]-end[tt] options will be set automatically as well.",
        "Again, if you think you know it better, please check the [TT]eigenfreq.xvg[tt]",
        "output."
    };

    static gmx_bool        bM    = TRUE, bCons = FALSE;
    static int             begin = 1, end = 50, maxspec = 4000;
    static real            T     = 298.15, width = 1;
    t_pargs                pa[]  =
    {
        { "-m",  FALSE, etBOOL, {&bM},
          "Divide elements of Hessian by product of sqrt(mass) of involved "
          "atoms prior to diagonalization. This should be used for 'Normal Modes' "
          "analysis" },
        { "-first", FALSE, etINT, {&begin},
          "First eigenvector to write away" },
        { "-last",  FALSE, etINT, {&end},
          "Last eigenvector to write away" },
        { "-maxspec", FALSE, etINT, {&maxspec},
          "Highest frequency (1/cm) to consider in the spectrum" },
        { "-T",     FALSE, etREAL, {&T},
          "Temperature for computing quantum heat capacity and enthalpy when using normal mode calculations to correct classical simulations" },
        { "-constr", FALSE, etBOOL, {&bCons},
          "If constraints were used in the simulation but not in the normal mode analysis (this is the recommended way of doing it) you will need to set this for computing the quantum corrections." },
        { "-width",  FALSE, etREAL, {&width},
          "Width (sigma) of the gaussian peaks (1/cm) when generating a spectrum" }
    };
    FILE                  *out, *qc, *spec;
    int                    status, trjout;
    t_topology             top;
    gmx_mtop_t             mtop;
    int                    ePBC;
    rvec                  *top_x;
    matrix                 box;
    real                  *eigenvalues;
    real                  *eigenvectors;
    real                   rdum, mass_fac, qcvtot, qutot, qcv, qu;
    int                    natoms, ndim, nrow, ncol, count, nharm, nvsite;
    char                  *grpname;
    int                    i, j, k, l, d, gnx;
    gmx_bool               bSuck;
    atom_id               *index;
    t_tpxheader            tpx;
    int                    version, generation;
    real                   value, omega, nu;
    real                   factor_gmx_to_omega2;
    real                   factor_omega_to_wavenumber;
    real                  *spectrum = NULL;
    real                   wfac;
    output_env_t           oenv;
    const char            *qcleg[] = {
        "Heat Capacity cV (J/mol K)",
        "Enthalpy H (kJ/mol)"
    };
    real *                 full_hessian   = NULL;
    gmx_sparsematrix_t *   sparse_hessian = NULL;

    t_filenm               fnm[] = {
        { efMTX, "-f", "hessian",    ffREAD  },
        { efTPR, NULL, NULL,         ffREAD  },
        { efXVG, "-of", "eigenfreq", ffWRITE },
        { efXVG, "-ol", "eigenval",  ffWRITE },
        { efXVG, "-os", "spectrum",  ffOPTWR },
        { efXVG, "-qc", "quant_corr",  ffOPTWR },
        { efTRN, "-v", "eigenvec",  ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    /* Read tpr file for volume and number of harmonic terms */
    read_tpxheader(ftp2fn(efTPR, NFILE, fnm), &tpx, TRUE, &version, &generation);
    snew(top_x, tpx.natoms);

    read_tpx(ftp2fn(efTPR, NFILE, fnm), NULL, box, &natoms,
             top_x, NULL, NULL, &mtop);
    if (bCons)
    {
        nharm = get_nharm(&mtop, &nvsite);
    }
    else
    {
        nharm  = 0;
        nvsite = 0;
    }
    top = gmx_mtop_t_to_t_topology(&mtop);

    bM   = TRUE;
    ndim = DIM*natoms;

    if (opt2bSet("-qc", NFILE, fnm))
    {
        begin = 7+DIM*nvsite;
        end   = DIM*natoms;
    }
    if (begin < 1)
    {
        begin = 1;
    }
    if (end > ndim)
    {
        end = ndim;
    }
    printf("Using begin = %d and end = %d\n", begin, end);

    /*open Hessian matrix */
    gmx_mtxio_read(ftp2fn(efMTX, NFILE, fnm), &nrow, &ncol, &full_hessian, &sparse_hessian);

    /* Memory for eigenvalues and eigenvectors (begin..end) */
    snew(eigenvalues, nrow);
    snew(eigenvectors, nrow*(end-begin+1));

    /* If the Hessian is in sparse format we can calculate max (ndim-1) eigenvectors,
     * and they must start at the first one. If this is not valid we convert to full matrix
     * storage, but warn the user that we might run out of memory...
     */
    if ((sparse_hessian != NULL) && (begin != 1 || end == ndim))
    {
        if (begin != 1)
        {
            fprintf(stderr, "Cannot use sparse Hessian with first eigenvector != 1.\n");
        }
        else if (end == ndim)
        {
            fprintf(stderr, "Cannot use sparse Hessian to calculate all eigenvectors.\n");
        }

        fprintf(stderr, "Will try to allocate memory and convert to full matrix representation...\n");

        snew(full_hessian, nrow*ncol);
        for (i = 0; i < nrow*ncol; i++)
        {
            full_hessian[i] = 0;
        }

        for (i = 0; i < sparse_hessian->nrow; i++)
        {
            for (j = 0; j < sparse_hessian->ndata[i]; j++)
            {
                k                      = sparse_hessian->data[i][j].col;
                value                  = sparse_hessian->data[i][j].value;
                full_hessian[i*ndim+k] = value;
                full_hessian[k*ndim+i] = value;
            }
        }
        gmx_sparsematrix_destroy(sparse_hessian);
        sparse_hessian = NULL;
        fprintf(stderr, "Converted sparse to full matrix storage.\n");
    }

    if (full_hessian != NULL)
    {
        /* Using full matrix storage */
        nma_full_hessian(full_hessian, nrow, bM, &top, begin, end,
                         eigenvalues, eigenvectors);
    }
    else
    {
        /* Sparse memory storage, allocate memory for eigenvectors */
        snew(eigenvectors, ncol*end);
        nma_sparse_hessian(sparse_hessian, bM, &top, end, eigenvalues, eigenvectors);
    }

    /* check the output, first 6 eigenvalues should be reasonably small */
    bSuck = FALSE;
    for (i = begin-1; (i < 6); i++)
    {
        if (fabs(eigenvalues[i]) > 1.0e-3)
        {
            bSuck = TRUE;
        }
    }
    if (bSuck)
    {
        fprintf(stderr, "\nOne of the lowest 6 eigenvalues has a non-zero value.\n");
        fprintf(stderr, "This could mean that the reference structure was not\n");
        fprintf(stderr, "properly energy minimized.\n");
    }

    /* now write the output */
    fprintf (stderr, "Writing eigenvalues...\n");
    out = xvgropen(opt2fn("-ol", NFILE, fnm),
                   "Eigenvalues", "Eigenvalue index", "Eigenvalue [Gromacs units]",
                   oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        if (bM)
        {
            fprintf(out, "@ subtitle \"mass weighted\"\n");
        }
        else
        {
            fprintf(out, "@ subtitle \"not mass weighted\"\n");
        }
    }

    for (i = 0; i <= (end-begin); i++)
    {
        fprintf (out, "%6d %15g\n", begin+i, eigenvalues[i]);
    }
    xvgrclose(out);


    if (opt2bSet("-qc", NFILE, fnm))
    {
        qc = xvgropen(opt2fn("-qc", NFILE, fnm), "Quantum Corrections", "Eigenvector index", "", oenv);
        xvgr_legend(qc, asize(qcleg), qcleg, oenv);
        qcvtot = qutot = 0;
    }
    else
    {
        qc = NULL;
    }
    printf("Writing eigenfrequencies - negative eigenvalues will be set to zero.\n");

    out = xvgropen(opt2fn("-of", NFILE, fnm),
                   "Eigenfrequencies", "Eigenvector index", "Wavenumber [cm\\S-1\\N]",
                   oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        if (bM)
        {
            fprintf(out, "@ subtitle \"mass weighted\"\n");
        }
        else
        {
            fprintf(out, "@ subtitle \"not mass weighted\"\n");
        }
    }
    /* Spectrum ? */
    spec = NULL;
    if (opt2bSet("-os", NFILE, fnm) && (maxspec > 0))
    {
        snew(spectrum, maxspec);
        spec = xvgropen(opt2fn("-os", NFILE, fnm),
                        "Vibrational spectrum based on harmonic approximation",
                        "\\f{12}w\\f{4} (cm\\S-1\\N)",
                        "Intensity [Gromacs units]",
                        oenv);
        for (i = 0; (i < maxspec); i++)
        {
            spectrum[i] = 0;
        }
    }

    /* Gromacs units are kJ/(mol*nm*nm*amu),
     * where amu is the atomic mass unit.
     *
     * For the eigenfrequencies we want to convert this to spectroscopic absorption
     * wavenumbers given in cm^(-1), which is the frequency divided by the speed of
     * light. Do this by first converting to omega^2 (units 1/s), take the square
     * root, and finally divide by the speed of light (nm/ps in gromacs).
     */
    factor_gmx_to_omega2       = 1.0E21/(AVOGADRO*AMU);
    factor_omega_to_wavenumber = 1.0E-5/(2.0*M_PI*SPEED_OF_LIGHT);

    for (i = begin; (i <= end); i++)
    {
        value = eigenvalues[i-begin];
        if (value < 0)
        {
            value = 0;
        }
        omega = sqrt(value*factor_gmx_to_omega2);
        nu    = 1e-12*omega/(2*M_PI);
        value = omega*factor_omega_to_wavenumber;
        fprintf (out, "%6d %15g\n", i, value);
        if (NULL != spec)
        {
            wfac = eigenvalues[i-begin]/(width*sqrt(2*M_PI));
            for (j = 0; (j < maxspec); j++)
            {
                spectrum[j] += wfac*exp(-sqr(j-value)/(2*sqr(width)));
            }
        }
        if (NULL != qc)
        {
            qcv = cv_corr(nu, T);
            qu  = u_corr(nu, T);
            if (i > end-nharm)
            {
                qcv += BOLTZ*KILO;
                qu  += BOLTZ*T;
            }
            fprintf (qc, "%6d %15g %15g\n", i, qcv, qu);
            qcvtot += qcv;
            qutot  += qu;
        }
    }
    xvgrclose(out);
    if (NULL != spec)
    {
        for (j = 0; (j < maxspec); j++)
        {
            fprintf(spec, "%10g  %10g\n", 1.0*j, spectrum[j]);
        }
        xvgrclose(spec);
    }
    if (NULL != qc)
    {
        printf("Quantum corrections for harmonic degrees of freedom\n");
        printf("Use appropriate -first and -last options to get reliable results.\n");
        printf("There were %d constraints and %d vsites in the simulation\n",
               nharm, nvsite);
        printf("Total correction to cV = %g J/mol K\n", qcvtot);
        printf("Total correction to  H = %g kJ/mol\n", qutot);
        xvgrclose(qc);
        please_cite(stdout, "Caleman2011b");
    }
    /* Writing eigenvectors. Note that if mass scaling was used, the eigenvectors
     * were scaled back from mass weighted cartesian to plain cartesian in the
     * nma_full_hessian() or nma_sparse_hessian() routines. Mass scaled vectors
     * will not be strictly orthogonal in plain cartesian scalar products.
     */
    write_eigenvectors(opt2fn("-v", NFILE, fnm), natoms, eigenvectors, FALSE, begin, end,
                       eWXR_NO, NULL, FALSE, top_x, bM, eigenvalues);

    return 0;
}
