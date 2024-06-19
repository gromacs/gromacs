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

#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/eigio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "thermochemistry.h"

struct gmx_output_env_t;

static double cv_corr(double nu, double T)
{
    double x  = gmx::c_planck * nu / (gmx::c_boltz * T);
    double ex = std::exp(x);

    if (nu <= 0)
    {
        return gmx::c_boltz * gmx::c_kilo;
    }
    else
    {
        return gmx::c_boltz * gmx::c_kilo * (ex * gmx::square(x) / gmx::square(ex - 1) - 1);
    }
}

static double u_corr(double nu, double T)
{
    double x  = gmx::c_planck * nu / (gmx::c_boltz * T);
    double ex = std::exp(x);

    if (nu <= 0)
    {
        return gmx::c_boltz * T;
    }
    else
    {
        return gmx::c_boltz * T * (0.5 * x - 1 + x / (ex - 1));
    }
}

static size_t get_nharm_mt(const gmx_moltype_t* mt)
{
    static int harm_func[] = { F_BONDS };
    int        i, ft;
    size_t     nh = 0;

    for (i = 0; (i < asize(harm_func)); i++)
    {
        ft = harm_func[i];
        nh += mt->ilist[ft].size() / (interaction_function[ft].nratoms + 1);
    }
    return nh;
}

static int get_nharm(const gmx_mtop_t* mtop)
{
    int nh = 0;

    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        nh += molb.nmol * get_nharm_mt(&(mtop->moltype[molb.type]));
    }
    return nh;
}

static void nma_full_hessian(real*                    hess,
                             int                      ndim,
                             gmx_bool                 bM,
                             const t_topology*        top,
                             gmx::ArrayRef<const int> atom_index,
                             int                      begin,
                             int                      end,
                             real*                    eigenvalues,
                             real*                    eigenvectors)
{
    real mass_fac;

    /* divide elements hess[i][j] by sqrt(mas[i])*sqrt(mas[j]) when required */

    if (bM)
    {
        for (int i = 0; (i < atom_index.ssize()); i++)
        {
            size_t ai = atom_index[i];
            for (size_t j = 0; (j < DIM); j++)
            {
                for (int k = 0; (k < atom_index.ssize()); k++)
                {
                    size_t ak = atom_index[k];
                    mass_fac  = gmx::invsqrt(top->atoms.atom[ai].m * top->atoms.atom[ak].m);
                    for (size_t l = 0; (l < DIM); l++)
                    {
                        hess[(i * DIM + j) * ndim + k * DIM + l] *= mass_fac;
                    }
                }
            }
        }
    }

    /* call diagonalization routine. */

    fprintf(stderr, "\nDiagonalizing to find vectors %d through %d...\n", begin, end);
    fflush(stderr);

    eigensolver(hess, ndim, begin - 1, end - 1, eigenvalues, eigenvectors);

    /* And scale the output eigenvectors */
    if (bM && eigenvectors != nullptr)
    {
        for (int i = 0; i < (end - begin + 1); i++)
        {
            for (gmx::Index j = 0; j < atom_index.ssize(); j++)
            {
                size_t aj = atom_index[j];
                mass_fac  = gmx::invsqrt(top->atoms.atom[aj].m);
                for (size_t k = 0; (k < DIM); k++)
                {
                    eigenvectors[i * ndim + j * DIM + k] *= mass_fac;
                }
            }
        }
    }
}


static void nma_sparse_hessian(gmx_sparsematrix_t*      sparse_hessian,
                               gmx_bool                 bM,
                               const t_topology*        top,
                               gmx::ArrayRef<const int> atom_index,
                               int                      neig,
                               real*                    eigenvalues,
                               real*                    eigenvectors)
{
    int    i, k;
    int    row, col;
    real   mass_fac;
    int    katom;
    size_t ndim;

    ndim = DIM * atom_index.size();

    /* Cannot check symmetry since we only store half matrix */
    /* divide elements hess[i][j] by sqrt(mas[i])*sqrt(mas[j]) when required */

    GMX_RELEASE_ASSERT(sparse_hessian != nullptr,
                       "NULL matrix pointer provided to nma_sparse_hessian");

    if (bM)
    {
        for (int iatom = 0; (iatom < atom_index.ssize()); iatom++)
        {
            size_t ai = atom_index[iatom];
            for (size_t j = 0; (j < DIM); j++)
            {
                row = DIM * iatom + j;
                for (k = 0; k < sparse_hessian->ndata[row]; k++)
                {
                    col       = sparse_hessian->data[row][k].col;
                    katom     = col / 3;
                    size_t ak = atom_index[katom];
                    mass_fac  = gmx::invsqrt(top->atoms.atom[ai].m * top->atoms.atom[ak].m);
                    sparse_hessian->data[row][k].value *= mass_fac;
                }
            }
        }
    }
    fprintf(stderr, "\nDiagonalizing to find eigenvectors 1 through %d...\n", neig);
    fflush(stderr);

    sparse_eigensolver(sparse_hessian, neig, eigenvalues, eigenvectors, 10000000);

    /* Scale output eigenvectors */
    if (bM && eigenvectors != nullptr)
    {
        for (i = 0; i < neig; i++)
        {
            for (gmx::Index j = 0; j < atom_index.ssize(); j++)
            {
                size_t aj = atom_index[j];
                mass_fac  = gmx::invsqrt(top->atoms.atom[aj].m);
                for (k = 0; (k < DIM); k++)
                {
                    eigenvectors[i * ndim + j * DIM + k] *= mass_fac;
                }
            }
        }
    }
}


/* Returns a pointer for eigenvector storage */
static real* allocateEigenvectors(int nrow, int first, int last, bool ignoreBegin)
{
    int numVector;
    if (ignoreBegin)
    {
        numVector = last;
    }
    else
    {
        numVector = last - first + 1;
    }
    size_t vectorsSize = static_cast<size_t>(nrow) * static_cast<size_t>(numVector);
    /* We can't have more than INT_MAX elements.
     * Relaxing this restriction probably requires changing lots of loop
     * variable types in the linear algebra code.
     */
    if (vectorsSize > INT_MAX)
    {
        gmx_fatal(FARGS,
                  "You asked to store %d eigenvectors of size %d, which requires more than the "
                  "supported %d elements; %sdecrease -last",
                  numVector,
                  nrow,
                  INT_MAX,
                  ignoreBegin ? "" : "increase -first and/or ");
    }

    real* eigenvectors;
    snew(eigenvectors, vectorsSize);

    return eigenvectors;
}

/*! \brief Compute heat capacity due to translational motion
 */
static double calcTranslationalHeatCapacity()
{
    return gmx::c_universalGasConstant * 1.5;
}

/*! \brief Compute internal energy due to translational motion
 */
static double calcTranslationalInternalEnergy(double T)
{
    return gmx::c_boltz * T * 1.5;
}

/*! \brief Compute heat capacity due to rotational motion
 *
 * \param[in] linear Should be TRUE if this is a linear molecule
 * \param[in] T      Temperature
 * \return The heat capacity at constant volume
 */
static double calcRotationalInternalEnergy(gmx_bool linear, double T)
{
    if (linear)
    {
        return gmx::c_boltz * T;
    }
    else
    {
        return gmx::c_boltz * T * 1.5;
    }
}

/*! \brief Compute heat capacity due to rotational motion
 *
 * \param[in] linear Should be TRUE if this is a linear molecule
 * \return The heat capacity at constant volume
 */
static double calcRotationalHeatCapacity(gmx_bool linear)
{
    if (linear)
    {
        return gmx::c_universalGasConstant;
    }
    else
    {
        return gmx::c_universalGasConstant * 1.5;
    }
}

static void analyzeThermochemistry(FILE*                    fp,
                                   const t_topology&        top,
                                   rvec                     top_x[],
                                   gmx::ArrayRef<const int> atom_index,
                                   real                     eigfreq[],
                                   real                     T,
                                   real                     P,
                                   int                      sigma_r,
                                   real                     scale_factor,
                                   real                     linear_toler)
{
    std::vector<int> index(atom_index.begin(), atom_index.end());

    rvec   xcm;
    double tmass  = calc_xcm(top_x, index.size(), index.data(), top.atoms.atom, xcm, FALSE);
    double Strans = calcTranslationalEntropy(tmass, T, P);
    std::vector<gmx::RVec> x_com;
    x_com.resize(top.atoms.nr);
    for (int i = 0; i < top.atoms.nr; i++)
    {
        copy_rvec(top_x[i], x_com[i]);
    }
    (void)sub_xcm(as_rvec_array(x_com.data()), index.size(), index.data(), top.atoms.atom, xcm, FALSE);

    rvec   inertia;
    matrix trans;
    principal_comp(index.size(), index.data(), top.atoms.atom, as_rvec_array(x_com.data()), trans, inertia);
    bool linear = (inertia[XX] / inertia[YY] < linear_toler && inertia[XX] / inertia[ZZ] < linear_toler);
    // (kJ/mol ps)^2/(Dalton nm^2 kJ/mol K) =
    // c_kilo kg m^2 ps^2/(s^2 mol g/mol nm^2 K) =
    // c_kilo^2 10^18 / 10^24 K = 1/K
    double rot_const = gmx::square(gmx::c_planck) / (8 * gmx::square(M_PI) * gmx::c_boltz);
    // Rotational temperature (1/K)
    rvec theta = { 0, 0, 0 };
    if (linear)
    {
        // For linear molecules the first element of the inertia
        // vector is zero.
        theta[0] = rot_const / inertia[1];
    }
    else
    {
        for (int m = 0; m < DIM; m++)
        {
            theta[m] = rot_const / inertia[m];
        }
    }
    if (debug)
    {
        pr_rvec(debug, 0, "inertia", inertia, DIM, TRUE);
        pr_rvec(debug, 0, "theta", theta, DIM, TRUE);
        pr_rvecs(debug, 0, "trans", trans, DIM);
        fprintf(debug, "linear molecule = %s\n", linear ? "true" : "no");
    }
    size_t nFreq = index.size() * DIM;
    auto   eFreq = gmx::arrayRefFromArray(eigfreq, nFreq);
    double Svib  = calcQuasiHarmonicEntropy(eFreq, T, linear, scale_factor);

    double Srot = calcRotationalEntropy(T, top.atoms.nr, linear, theta, sigma_r);
    fprintf(fp, "Translational entropy %g J/mol K\n", Strans);
    fprintf(fp, "Rotational entropy    %g J/mol K\n", Srot);
    fprintf(fp, "Vibrational entropy   %g J/mol K\n", Svib);
    fprintf(fp, "Total entropy         %g J/mol K\n", Svib + Strans + Srot);
    double HeatCapacity = (calcTranslationalHeatCapacity() + calcRotationalHeatCapacity(linear)
                           + calcVibrationalHeatCapacity(eFreq, T, linear, scale_factor));
    fprintf(fp, "Heat capacity         %g J/mol K\n", HeatCapacity);
    double Evib = (calcTranslationalInternalEnergy(T) + calcRotationalInternalEnergy(linear, T)
                   + calcVibrationalInternalEnergy(eFreq, T, linear, scale_factor));
    fprintf(fp, "Internal energy       %g kJ/mol\n", Evib);
    double ZPE = calcZeroPointEnergy(eFreq, scale_factor);
    fprintf(fp, "Zero-point energy     %g kJ/mol\n", ZPE);
}


int gmx_nmeig(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] calculates the eigenvectors/values of a (Hessian) matrix,",
        "which can be calculated with [gmx-mdrun].",
        "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
        "The structure is written first with t=0. The eigenvectors",
        "are written as frames with the eigenvector number and eigenvalue",
        "written as step number and timestamp, respectively.",
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
        "[TT]-qc[tt], the [TT]-begin[tt] and [TT]-end[tt] options will be set automatically as ",
        "well.[PAR]",
        "Based on a harmonic analysis of the normal mode frequencies,",
        "thermochemical properties S0 (Standard Entropy),",
        "Cv (Heat capacity at constant volume), Zero-point energy and the internal energy are",
        "computed, much in the same manner as popular quantum chemistry",
        "programs."
    };

    static gmx_bool bM = TRUE, bCons = FALSE;
    static int      begin = 1, end = 50, maxspec = 4000, sigma_r = 1;
    static real     T = 298.15, width = 1, P = 1, scale_factor = 1;
    static real     linear_toler = 1e-5;

    t_pargs pa[] = {
        { "-m",
          FALSE,
          etBOOL,
          { &bM },
          "Divide elements of Hessian by product of sqrt(mass) of involved "
          "atoms prior to diagonalization. This should be used for 'Normal Modes' "
          "analysis" },
        { "-first", FALSE, etINT, { &begin }, "First eigenvector to write away" },
        { "-last",
          FALSE,
          etINT,
          { &end },
          "Last eigenvector to write away. -1 is use all dimensions." },
        { "-maxspec",
          FALSE,
          etINT,
          { &maxspec },
          "Highest frequency (1/cm) to consider in the spectrum" },
        { "-T",
          FALSE,
          etREAL,
          { &T },
          "Temperature for computing entropy, quantum heat capacity and enthalpy when "
          "using normal mode calculations to correct classical simulations" },
        { "-P", FALSE, etREAL, { &P }, "Pressure (bar) when computing entropy" },
        { "-sigma",
          FALSE,
          etINT,
          { &sigma_r },
          "Number of symmetric copies used when computing entropy. E.g. for water the "
          "number is 2, for NH3 it is 3 and for methane it is 12." },
        { "-scale",
          FALSE,
          etREAL,
          { &scale_factor },
          "Factor to scale frequencies before computing thermochemistry values" },
        { "-linear_toler",
          FALSE,
          etREAL,
          { &linear_toler },
          "Tolerance for determining whether a compound is linear as determined from the "
          "ration of the moments inertion Ix/Iy and Ix/Iz." },
        { "-constr",
          FALSE,
          etBOOL,
          { &bCons },
          "If constraints were used in the simulation but not in the normal mode analysis "
          "you will need to set this for computing the quantum corrections." },
        { "-width",
          FALSE,
          etREAL,
          { &width },
          "Width (sigma) of the gaussian peaks (1/cm) when generating a spectrum" }
    };
    FILE *                     out, *qc, *spec;
    t_topology                 top;
    gmx_mtop_t                 mtop;
    rvec*                      top_x;
    matrix                     box;
    real*                      eigenvalues;
    real*                      eigenvectors;
    real                       qcvtot, qutot, qcv, qu;
    int                        i, j, k;
    real                       value, omega, nu;
    real                       factor_gmx_to_omega2;
    real                       factor_omega_to_wavenumber;
    real*                      spectrum = nullptr;
    real                       wfac;
    gmx_output_env_t*          oenv;
    std::array<std::string, 2> qcleg = { "Heat Capacity cV (J/mol K)", "Enthalpy H (kJ/mol)" };
    real*                      full_hessian   = nullptr;
    gmx_sparsematrix_t*        sparse_hessian = nullptr;

    t_filenm fnm[] = {
        { efMTX, "-f", "hessian", ffREAD },     { efTPR, nullptr, nullptr, ffREAD },
        { efXVG, "-of", "eigenfreq", ffWRITE }, { efXVG, "-ol", "eigenval", ffWRITE },
        { efXVG, "-os", "spectrum", ffOPTWR },  { efXVG, "-qc", "quant_corr", ffOPTWR },
        { efTRN, "-v", "eigenvec", ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    /* Read tpr file for volume and number of harmonic terms */
    TpxFileHeader tpx = readTpxHeader(ftp2fn(efTPR, NFILE, fnm), true);
    snew(top_x, tpx.natoms);

    int natoms_tpx;
    read_tpx(ftp2fn(efTPR, NFILE, fnm), nullptr, box, &natoms_tpx, top_x, nullptr, &mtop);
    int nharm = 0;
    if (bCons)
    {
        nharm = get_nharm(&mtop);
    }
    std::vector<int> atom_index = get_atom_index(mtop);

    top = gmx_mtop_t_to_t_topology(&mtop, true);

    bM       = TRUE;
    int ndim = DIM * atom_index.size();

    if (opt2bSet("-qc", NFILE, fnm))
    {
        begin = 7;
        end   = ndim;
    }
    if (begin < 1)
    {
        begin = 1;
    }
    if (end == -1 || end > ndim)
    {
        end = ndim;
    }
    printf("Using begin = %d and end = %d\n", begin, end);

    /*open Hessian matrix */
    int nrow, ncol;
    gmx_mtxio_read(ftp2fn(efMTX, NFILE, fnm), &nrow, &ncol, &full_hessian, &sparse_hessian);

    /* If the Hessian is in sparse format we can calculate max (ndim-1) eigenvectors,
     * If this is not valid we convert to full matrix storage,
     * but warn the user that we might run out of memory...
     */
    if ((sparse_hessian != nullptr) && (end == ndim))
    {
        fprintf(stderr, "Cannot use sparse Hessian to calculate all eigenvectors.\n");

        fprintf(stderr,
                "Will try to allocate memory and convert to full matrix representation...\n");

        size_t hessianSize = static_cast<size_t>(nrow) * static_cast<size_t>(ncol);
        /* Allowing  Hessians larger than INT_MAX probably only makes sense
         * with (OpenMP) parallel diagonalization routines, since with a single
         * thread it will takes months.
         */
        if (hessianSize > INT_MAX)
        {
            gmx_fatal(FARGS,
                      "Hessian size is %d x %d, which is larger than the maximum allowed %d "
                      "elements.",
                      nrow,
                      ncol,
                      INT_MAX);
        }
        snew(full_hessian, hessianSize);
        for (i = 0; i < nrow * ncol; i++)
        {
            full_hessian[i] = 0;
        }

        for (i = 0; i < sparse_hessian->nrow; i++)
        {
            for (j = 0; j < sparse_hessian->ndata[i]; j++)
            {
                k                          = sparse_hessian->data[i][j].col;
                value                      = sparse_hessian->data[i][j].value;
                full_hessian[i * ndim + k] = value;
                full_hessian[k * ndim + i] = value;
            }
        }
        gmx_sparsematrix_destroy(sparse_hessian);
        sparse_hessian = nullptr;
        fprintf(stderr, "Converted sparse to full matrix storage.\n");
    }

    snew(eigenvalues, nrow);

    if (full_hessian != nullptr)
    {
        /* Using full matrix storage */
        eigenvectors = allocateEigenvectors(nrow, begin, end, false);

        nma_full_hessian(full_hessian, nrow, bM, &top, atom_index, begin, end, eigenvalues, eigenvectors);
    }
    else
    {
        assert(sparse_hessian);
        /* Sparse memory storage, allocate memory for eigenvectors */
        eigenvectors = allocateEigenvectors(nrow, begin, end, true);

        nma_sparse_hessian(sparse_hessian, bM, &top, atom_index, end, eigenvalues, eigenvectors);
    }

    /* check the output, first 6 eigenvalues should be reasonably small */
    gmx_bool bSuck = FALSE;
    for (i = begin - 1; (i < 6); i++)
    {
        if (std::abs(eigenvalues[i]) > 1.0e-3)
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
    fprintf(stderr, "Writing eigenvalues...\n");
    out = xvgropen(opt2fn("-ol", NFILE, fnm),
                   "Eigenvalues",
                   "Eigenvalue index",
                   "Eigenvalue [Gromacs units]",
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

    for (i = 0; i <= (end - begin); i++)
    {
        fprintf(out, "%6d %15g\n", begin + i, eigenvalues[i]);
    }
    xvgrclose(out);


    if (opt2bSet("-qc", NFILE, fnm))
    {
        qc = xvgropen(opt2fn("-qc", NFILE, fnm), "Quantum Corrections", "Eigenvector index", "", oenv);
        xvgrLegend(qc, qcleg, oenv);
        qcvtot = qutot = 0;
    }
    else
    {
        qc = nullptr;
    }
    printf("Writing eigenfrequencies - negative eigenvalues will be set to zero.\n");

    out = xvgropen(opt2fn("-of", NFILE, fnm),
                   "Eigenfrequencies",
                   "Eigenvector index",
                   "Wavenumber [cm\\S-1\\N]",
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
    spec = nullptr;
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
    factor_gmx_to_omega2       = 1.0E21 / (gmx::c_avogadro * gmx::c_amu);
    factor_omega_to_wavenumber = 1.0E-5 / (2.0 * M_PI * gmx::c_speedOfLight);

    value = 0;
    for (i = begin; (i <= end); i++)
    {
        value = eigenvalues[i - begin];
        if (value < 0)
        {
            value = 0;
        }
        omega = std::sqrt(value * factor_gmx_to_omega2);
        nu    = 1e-12 * omega / (2 * M_PI);
        value = omega * factor_omega_to_wavenumber;
        fprintf(out, "%6d %15g\n", i, value);
        if (nullptr != spec)
        {
            wfac = eigenvalues[i - begin] / (width * std::sqrt(2 * M_PI));
            for (j = 0; (j < maxspec); j++)
            {
                spectrum[j] += wfac * std::exp(-gmx::square(j - value) / (2 * gmx::square(width)));
            }
        }
        if (nullptr != qc)
        {
            qcv = cv_corr(nu, T);
            qu  = u_corr(nu, T);
            if (i > end - nharm)
            {
                qcv += gmx::c_boltz * gmx::c_kilo;
                qu += gmx::c_boltz * T;
            }
            fprintf(qc, "%6d %15g %15g\n", i, qcv, qu);
            qcvtot += qcv;
            qutot += qu;
        }
    }
    xvgrclose(out);
    if (value >= maxspec)
    {
        printf("WARNING: high frequencies encountered (%g cm^-1).\n", value);
        printf("Your calculations may be incorrect due to e.g. improper minimization of\n");
        printf("your starting structure or due to issues in your topology.\n");
    }
    if (nullptr != spec)
    {
        for (j = 0; (j < maxspec); j++)
        {
            fprintf(spec, "%10g  %10g\n", 1.0 * j, spectrum[j]);
        }
        xvgrclose(spec);
    }
    if (nullptr != qc)
    {
        printf("Quantum corrections for harmonic degrees of freedom\n");
        printf("Use appropriate -first and -last options to get reliable results.\n");
        printf("There were %d constraints in the simulation\n", nharm);
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
    const real* eigenvectorPtr;
    if (full_hessian != nullptr)
    {
        eigenvectorPtr = eigenvectors;
    }
    else
    {
        /* The sparse matrix diagonalization store all eigenvectors up to end */
        eigenvectorPtr = eigenvectors + (begin - 1) * atom_index.size();
    }
    write_eigenvectors(opt2fn("-v", NFILE, fnm),
                       atom_index.size(),
                       eigenvectorPtr,
                       FALSE,
                       begin,
                       end,
                       eWXR_NO,
                       nullptr,
                       FALSE,
                       top_x,
                       bM,
                       eigenvalues);

    if (begin == 1)
    {
        analyzeThermochemistry(
                stdout, top, top_x, atom_index, eigenvalues, T, P, sigma_r, scale_factor, linear_toler);
        please_cite(stdout, "Spoel2018a");
    }
    else
    {
        printf("Cannot compute entropy when -first = %d\n", begin);
    }

    return 0;
}
