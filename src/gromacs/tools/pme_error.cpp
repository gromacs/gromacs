/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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

#include "pme_error.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fft/calcgrid.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

struct gmx_output_env_t;

/* #define TAKETIME */
/* #define DEBUG  */

/* Enum for situations that can occur during log file parsing */
enum
{
    eParselogOK,
    eParselogNotFound,
    eParselogNoPerfData,
    eParselogTerm,
    eParselogResetProblem,
    eParselogNr
};


struct PmeErrorInputs
{
    int64_t orig_sim_steps;  /* Number of steps to be done in the real simulation  */
    int     n_entries;       /* Number of entries in arrays                        */
    real    volume;          /* The volume of the box                              */
    matrix  recipbox;        /* The reciprocal box                                 */
    int     natoms;          /* The number of atoms in the MD system               */
    real*   fac;             /* The scaling factor                                 */
    real*   rcoulomb;        /* The coulomb radii [0...nr_inputfiles]              */
    real*   rvdw;            /* The vdW radii                                      */
    int *   nkx, *nky, *nkz; /* Number of k vectors in each spatial dimension      */
    real*   fourier_sp;      /* Fourierspacing                                     */
    real*   ewald_rtol;      /* Real space tolerance for Ewald, determines         */
                             /* the real/reciprocal space relative weight          */
    real*    ewald_beta;     /* Splitting parameter [1/nm]                         */
    real     fracself;       /* fraction of particles for SI error                 */
    real     q2all;          /* sum ( q ^2 )                                       */
    real     q2allnr;        /* nr of charges                                      */
    int*     pme_order;      /* Interpolation order for PME (bsplines)             */
    char**   fn_out;         /* Name of the output tpr file                        */
    real*    e_dir;          /* Direct space part of PME error with these settings */
    real*    e_rec;          /* Reciprocal space part of PME error                 */
    gmx_bool bTUNE;          /* flag for tuning */
};


/* Returns TRUE when atom is charged */
static gmx_bool is_charge(real charge)
{
    return charge * charge > GMX_REAL_EPS;
}


/* calculate charge density */
static void calc_q2all(const gmx_mtop_t* mtop, /* molecular topology */
                       real*             q2all,
                       real*             q2allnr)
{
    real q2_all = 0; /* Sum of squared charges */
    int  nrq_mol;    /* Number of charges in a single molecule */
    int  nrq_all;    /* Total number of charges in the MD system */
    real qi, q2_mol;

#ifdef DEBUG
    fprintf(stderr, "\nCharge density:\n");
#endif
    q2_all  = 0.0; /* total q squared */
    nrq_all = 0;   /* total number of charges in the system */
    for (const gmx_molblock_t& molblock : mtop->molblock) /* Loop over molecule types */
    {
        q2_mol                        = 0.0; /* q squared value of this molecule */
        nrq_mol                       = 0;   /* number of charges this molecule carries */
        const gmx_moltype_t& molecule = mtop->moltype[molblock.type];
        for (int i = 0; i < molecule.atoms.nr; i++)
        {
            qi = molecule.atoms.atom[i].q;
            /* Is this charge worth to be considered? */
            if (is_charge(qi))
            {
                q2_mol += qi * qi;
                nrq_mol++;
            }
        }
        /* Multiply with the number of molecules present of this type and add */
        q2_all += q2_mol * molblock.nmol;
        nrq_all += nrq_mol * molblock.nmol;
#ifdef DEBUG
        fprintf(stderr,
                "Molecule %2d (%5d atoms) q2_mol=%10.3e nr.mol.charges=%5d (%6dx)  q2_all=%10.3e  "
                "tot.charges=%d\n",
                imol,
                molecule.atoms.nr,
                q2_mol,
                nrq_mol,
                molblock.nmol,
                q2_all,
                nrq_all);
#endif
    }

    *q2all   = q2_all;
    *q2allnr = nrq_all;
}


/* Estimate the direct space part error of the SPME Ewald sum */
static real estimate_direct(PmeErrorInputs* info)
{
    real e_dir     = 0; /* Error estimate */
    real beta      = 0; /* Splitting parameter (1/nm) */
    real r_coulomb = 0; /* Cut-off in direct space */


    beta      = info->ewald_beta[0];
    r_coulomb = info->rcoulomb[0];

    e_dir = 2.0 * info->q2all * gmx::invsqrt(info->q2allnr * r_coulomb * info->volume);
    e_dir *= std::exp(-beta * beta * r_coulomb * r_coulomb);

    return gmx::c_one4PiEps0 * e_dir;
}

#define SUMORDER 6

/* the following 4 functions determine polynomials required for the reciprocal error estimate */

static inline real eps_poly1(real m, /* grid coordinate in certain direction */
                             real K, /* grid size in corresponding direction */
                             real n) /* spline interpolation order of the SPME */
{
    int  i;
    real nom   = 0; /* nominator */
    real denom = 0; /* denominator */
    real tmp   = 0;

    if (m == 0.0)
    {
        return 0.0;
    }

    for (i = -SUMORDER; i < 0; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += std::pow(tmp, -n);
    }

    for (i = SUMORDER; i > 0; i--)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += std::pow(tmp, -n);
    }

    tmp = m / K;
    tmp *= 2.0 * M_PI;
    denom = std::pow(tmp, -n) + nom;

    return -nom / denom;
}

static inline real eps_poly2(real m, /* grid coordinate in certain direction */
                             real K, /* grid size in corresponding direction */
                             real n) /* spline interpolation order of the SPME */
{
    int  i;
    real nom   = 0; /* nominator */
    real denom = 0; /* denominator */
    real tmp   = 0;

    if (m == 0.0)
    {
        return 0.0;
    }

    for (i = -SUMORDER; i < 0; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += std::pow(tmp, -2 * n);
    }

    for (i = SUMORDER; i > 0; i--)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += std::pow(tmp, -2 * n);
    }

    for (i = -SUMORDER; i < SUMORDER + 1; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        denom += std::pow(tmp, -n);
    }
    tmp = eps_poly1(m, K, n);
    return nom / denom / denom + tmp * tmp;
}

static inline real eps_poly3(real m, /* grid coordinate in certain direction */
                             real K, /* grid size in corresponding direction */
                             real n) /* spline interpolation order of the SPME */
{
    int  i;
    real nom   = 0; /* nominator */
    real denom = 0; /* denominator */
    real tmp   = 0;

    if (m == 0.0)
    {
        return 0.0;
    }

    for (i = -SUMORDER; i < 0; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += i * std::pow(tmp, -2 * n);
    }

    for (i = SUMORDER; i > 0; i--)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += i * std::pow(tmp, -2 * n);
    }

    for (i = -SUMORDER; i < SUMORDER + 1; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        denom += std::pow(tmp, -n);
    }

    return 2.0 * M_PI * nom / denom / denom;
}

static inline real eps_poly4(real m, /* grid coordinate in certain direction */
                             real K, /* grid size in corresponding direction */
                             real n) /* spline interpolation order of the SPME */
{
    int  i;
    real nom   = 0; /* nominator */
    real denom = 0; /* denominator */
    real tmp   = 0;

    if (m == 0.0)
    {
        return 0.0;
    }

    for (i = -SUMORDER; i < 0; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += i * i * std::pow(tmp, -2 * n);
    }

    for (i = SUMORDER; i > 0; i--)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        nom += i * i * std::pow(tmp, -2 * n);
    }

    for (i = -SUMORDER; i < SUMORDER + 1; i++)
    {
        tmp = m / K + i;
        tmp *= 2.0 * M_PI;
        denom += std::pow(tmp, -n);
    }

    return 4.0 * M_PI * M_PI * nom / denom / denom;
}

static inline real eps_self(real m,     /* grid coordinate in certain direction */
                            real K,     /* grid size in corresponding direction */
                            rvec rboxv, /* reciprocal box vector */
                            real n,     /* spline interpolation order of the SPME */
                            rvec x)     /* coordinate of charge */
{
    int  i;
    real tmp    = 0; /* temporary variables for computations */
    real tmp1   = 0; /* temporary variables for computations */
    real tmp2   = 0; /* temporary variables for computations */
    real rcoord = 0; /* coordinate in certain reciprocal space direction */
    real nom    = 0; /* nominator */
    real denom  = 0; /* denominator */


    if (m == 0.0)
    {
        return 0.0;
    }

    rcoord = iprod(rboxv, x);


    for (i = -SUMORDER; i < 0; i++)
    {
        tmp  = -std::sin(2.0 * M_PI * i * K * rcoord);
        tmp1 = 2.0 * M_PI * m / K + 2.0 * M_PI * i;
        tmp2 = std::pow(tmp1, -n);
        nom += tmp * tmp2 * i;
        denom += tmp2;
    }

    for (i = SUMORDER; i > 0; i--)
    {
        tmp  = -std::sin(2.0 * M_PI * i * K * rcoord);
        tmp1 = 2.0 * M_PI * m / K + 2.0 * M_PI * i;
        tmp2 = std::pow(tmp1, -n);
        nom += tmp * tmp2 * i;
        denom += tmp2;
    }


    tmp  = 2.0 * M_PI * m / K;
    tmp1 = std::pow(tmp, -n);
    denom += tmp1;

    return 2.0 * M_PI * nom / denom * K;
}

#undef SUMORDER

/* The following routine is just a copy from pme.c */

static void calc_recipbox(matrix box, matrix recipbox)
{
    /* Save some time by assuming upper right part is zero */

    real tmp = 1.0 / (box[XX][XX] * box[YY][YY] * box[ZZ][ZZ]);

    recipbox[XX][XX] = box[YY][YY] * box[ZZ][ZZ] * tmp;
    recipbox[XX][YY] = 0;
    recipbox[XX][ZZ] = 0;
    recipbox[YY][XX] = -box[YY][XX] * box[ZZ][ZZ] * tmp;
    recipbox[YY][YY] = box[XX][XX] * box[ZZ][ZZ] * tmp;
    recipbox[YY][ZZ] = 0;
    recipbox[ZZ][XX] = (box[YY][XX] * box[ZZ][YY] - box[YY][YY] * box[ZZ][XX]) * tmp;
    recipbox[ZZ][YY] = -box[ZZ][YY] * box[XX][XX] * tmp;
    recipbox[ZZ][ZZ] = box[XX][XX] * box[YY][YY] * tmp;
}


/* Estimate the reciprocal space part error of the SPME Ewald sum. */
static real estimate_reciprocal(PmeErrorInputs* info,
                                rvec            x[], /* array of particles */
                                const real      q[], /* array of charges */
                                int nr, /* number of charges = size of the charge array */
                                FILE gmx_unused* fp_out,
                                gmx_bool         bVerbose,
                                int  seed,     /* The seed for the random number generator */
                                int* nsamples, /* Return the number of samples used if Monte Carlo
                                                * algorithm is used for self energy error estimate */
                                const gmx::MpiComm& mpiComm)
{
    real     e_rec   = 0; /* reciprocal error estimate */
    real     e_rec1  = 0; /* Error estimate term 1*/
    real     e_rec2  = 0; /* Error estimate term 2*/
    real     e_rec3  = 0; /* Error estimate term 3 */
    real     e_rec3x = 0; /* part of Error estimate term 3 in x */
    real     e_rec3y = 0; /* part of Error estimate term 3 in y */
    real     e_rec3z = 0; /* part of Error estimate term 3 in z */
    int      i, ci;
    int      nx, ny, nz; /* grid coordinates */
    real     q2_all = 0; /* sum of squared charges */
    rvec     gridpx;     /* reciprocal grid point in x direction*/
    rvec     gridpxy;    /* reciprocal grid point in x and y direction*/
    rvec     gridp;      /* complete reciprocal grid point in 3 directions*/
    rvec     tmpvec;     /* template to create points from basis vectors */
    rvec     tmpvec2;    /* template to create points from basis vectors */
    real     coeff  = 0; /* variable to compute coefficients of the error estimate */
    real     coeff2 = 0; /* variable to compute coefficients of the error estimate */
    real     tmp    = 0; /* variables to compute different factors from vectors */
    real     tmp1   = 0;
    real     tmp2   = 0;
    gmx_bool bFraction;

    int* numbers = nullptr;

    /* Index variables for parallel work distribution */
    int startglobal, stopglobal;
    int startlocal, stoplocal;
    int x_per_core;
    int xtot;

#ifdef TAKETIME
    double t0 = 0.0;
    double t1 = 0.0;
#endif

    GMX_RELEASE_ASSERT(q != nullptr, "Must have charges");

    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
    }
    fprintf(stderr, "Using random seed %d.\n", seed);

    gmx::DefaultRandomEngine         rng(seed);
    gmx::UniformIntDistribution<int> dist(0, nr - 1);

    clear_rvec(gridpx);
    clear_rvec(gridpxy);
    clear_rvec(gridp);
    clear_rvec(tmpvec);
    clear_rvec(tmpvec2);

    for (i = 0; i < nr; i++)
    {
        q2_all += q[i] * q[i];
    }

    /* Calculate indices for work distribution */
    startglobal = -info->nkx[0] / 2;
    stopglobal  = info->nkx[0] / 2;
    xtot        = stopglobal * 2 + 1;
    if (mpiComm.size() > 1)
    {
        x_per_core = static_cast<int>(std::ceil(static_cast<real>(xtot) / mpiComm.size()));
        startlocal = startglobal + x_per_core * mpiComm.rank();
        stoplocal  = startlocal + x_per_core - 1;
        if (stoplocal > stopglobal)
        {
            stoplocal = stopglobal;
        }
    }
    else
    {
        startlocal = startglobal;
        stoplocal  = stopglobal;
        x_per_core = xtot;
    }

#if GMX_LIB_MPI
#    ifdef TAKETIME
    if (mpiComm.isMainRank())
    {
        t0 = MPI_Wtime();
    }
#    endif
#endif

    if (mpiComm.isMainRank())
    {

        fprintf(stderr, "Calculating reciprocal error part 1 ...");
    }

    for (nx = startlocal; nx <= stoplocal; nx++)
    {
        svmul(nx, info->recipbox[XX], gridpx);
        for (ny = -info->nky[0] / 2; ny < info->nky[0] / 2 + 1; ny++)
        {
            svmul(ny, info->recipbox[YY], tmpvec);
            rvec_add(gridpx, tmpvec, gridpxy);
            for (nz = -info->nkz[0] / 2; nz < info->nkz[0] / 2 + 1; nz++)
            {
                if (0 == nx && 0 == ny && 0 == nz)
                {
                    continue;
                }
                svmul(nz, info->recipbox[ZZ], tmpvec);
                rvec_add(gridpxy, tmpvec, gridp);
                tmp = norm2(gridp);
                coeff = std::exp(-1.0 * M_PI * M_PI * tmp / info->ewald_beta[0] / info->ewald_beta[0]);
                coeff /= 2.0 * M_PI * info->volume * tmp;
                coeff2 = tmp;


                tmp = eps_poly2(nx, info->nkx[0], info->pme_order[0]);
                tmp += eps_poly2(ny, info->nkx[0], info->pme_order[0]);
                tmp += eps_poly2(nz, info->nkx[0], info->pme_order[0]);

                tmp1 = eps_poly1(nx, info->nkx[0], info->pme_order[0]);
                tmp2 = eps_poly1(ny, info->nky[0], info->pme_order[0]);

                tmp += 2.0 * tmp1 * tmp2;

                tmp1 = eps_poly1(nz, info->nkz[0], info->pme_order[0]);
                tmp2 = eps_poly1(ny, info->nky[0], info->pme_order[0]);

                tmp += 2.0 * tmp1 * tmp2;

                tmp1 = eps_poly1(nz, info->nkz[0], info->pme_order[0]);
                tmp2 = eps_poly1(nx, info->nkx[0], info->pme_order[0]);

                tmp += 2.0 * tmp1 * tmp2;

                tmp1 = eps_poly1(nx, info->nkx[0], info->pme_order[0]);
                tmp1 += eps_poly1(ny, info->nky[0], info->pme_order[0]);
                tmp1 += eps_poly1(nz, info->nkz[0], info->pme_order[0]);

                tmp += tmp1 * tmp1;

                e_rec1 += 32.0 * M_PI * M_PI * coeff * coeff * coeff2 * tmp * q2_all * q2_all / nr;

                tmp1 = eps_poly3(nx, info->nkx[0], info->pme_order[0]);
                tmp1 *= info->nkx[0];
                tmp2 = iprod(gridp, info->recipbox[XX]);

                tmp = tmp1 * tmp2;

                tmp1 = eps_poly3(ny, info->nky[0], info->pme_order[0]);
                tmp1 *= info->nky[0];
                tmp2 = iprod(gridp, info->recipbox[YY]);

                tmp += tmp1 * tmp2;

                tmp1 = eps_poly3(nz, info->nkz[0], info->pme_order[0]);
                tmp1 *= info->nkz[0];
                tmp2 = iprod(gridp, info->recipbox[ZZ]);

                tmp += tmp1 * tmp2;

                tmp *= 4.0 * M_PI;

                tmp1 = eps_poly4(nx, info->nkx[0], info->pme_order[0]);
                tmp1 *= norm2(info->recipbox[XX]);
                tmp1 *= info->nkx[0] * info->nkx[0];

                tmp += tmp1;

                tmp1 = eps_poly4(ny, info->nky[0], info->pme_order[0]);
                tmp1 *= norm2(info->recipbox[YY]);
                tmp1 *= info->nky[0] * info->nky[0];

                tmp += tmp1;

                tmp1 = eps_poly4(nz, info->nkz[0], info->pme_order[0]);
                tmp1 *= norm2(info->recipbox[ZZ]);
                tmp1 *= info->nkz[0] * info->nkz[0];

                tmp += tmp1;

                e_rec2 += 4.0 * coeff * coeff * tmp * q2_all * q2_all / nr;
            }
        }
        if (mpiComm.isMainRank())
        {
            fprintf(stderr,
                    "\rCalculating reciprocal error part 1 ... %3.0f%%",
                    100.0 * (nx - startlocal + 1) / (x_per_core));
            std::fflush(stderr);
        }
    }

    if (mpiComm.isMainRank())
    {
        fprintf(stderr, "\n");
    }

    /* Use just a fraction of all charges to estimate the self energy error term? */
    bFraction = (info->fracself > 0.0) && (info->fracself < 1.0);

    if (bFraction)
    {
        /* Here xtot is the number of samples taken for the Monte Carlo calculation
         * of the average of term IV of equation 35 in Wang2010. Round up to a
         * number of samples that is divisible by the number of nodes */
        x_per_core = static_cast<int>(std::ceil(info->fracself * nr / mpiComm.size()));
        xtot       = x_per_core * mpiComm.size();
    }
    else
    {
        /* In this case we use all nr particle positions */
        xtot       = nr;
        x_per_core = static_cast<int>(std::ceil(static_cast<real>(xtot) / mpiComm.size()));
    }

    startlocal = x_per_core * mpiComm.rank();
    stoplocal  = std::min(startlocal + x_per_core, xtot); /* min needed if xtot == nr */

    if (bFraction)
    {
        /* Make shure we get identical results in serial and parallel. Therefore,
         * take the sample indices from a single, global random number array that
         * is constructed on the main node and that only depends on the seed */
        snew(numbers, xtot);
        if (mpiComm.isMainRank())
        {
            for (i = 0; i < xtot; i++)
            {
                numbers[i] = dist(rng); // [0,nr-1]
            }
        }
        /* Broadcast the random number array to the other nodes */
        if (mpiComm.size() > 1)
        {
            nblock_bc(mpiComm.comm(), xtot, numbers);
        }

        if (bVerbose && mpiComm.isMainRank())
        {
            fprintf(stdout,
                    "Using %d sample%s to approximate the self interaction error term",
                    xtot,
                    xtot == 1 ? "" : "s");
            if (mpiComm.size() > 1)
            {
                fprintf(stdout, " (%d sample%s per rank)", x_per_core, x_per_core == 1 ? "" : "s");
            }
            fprintf(stdout, ".\n");
        }
    }

    /* Return the number of positions used for the Monte Carlo algorithm */
    *nsamples = xtot;

    for (i = startlocal; i < stoplocal; i++)
    {
        e_rec3x = 0;
        e_rec3y = 0;
        e_rec3z = 0;

        if (bFraction)
        {
            /* Randomly pick a charge */
            ci = numbers[i];
        }
        else
        {
            /* Use all charges */
            ci = i;
        }

        /* for(nx=startlocal; nx<=stoplocal; nx++)*/
        for (nx = -info->nkx[0] / 2; nx < info->nkx[0] / 2 + 1; nx++)
        {
            svmul(nx, info->recipbox[XX], gridpx);
            for (ny = -info->nky[0] / 2; ny < info->nky[0] / 2 + 1; ny++)
            {
                svmul(ny, info->recipbox[YY], tmpvec);
                rvec_add(gridpx, tmpvec, gridpxy);
                for (nz = -info->nkz[0] / 2; nz < info->nkz[0] / 2 + 1; nz++)
                {

                    if (0 == nx && 0 == ny && 0 == nz)
                    {
                        continue;
                    }

                    svmul(nz, info->recipbox[ZZ], tmpvec);
                    rvec_add(gridpxy, tmpvec, gridp);
                    tmp = norm2(gridp);
                    coeff = std::exp(-1.0 * M_PI * M_PI * tmp / info->ewald_beta[0] / info->ewald_beta[0]);
                    coeff /= tmp;
                    e_rec3x += coeff
                               * eps_self(nx, info->nkx[0], info->recipbox[XX], info->pme_order[0], x[ci]);
                    e_rec3y += coeff
                               * eps_self(ny, info->nky[0], info->recipbox[YY], info->pme_order[0], x[ci]);
                    e_rec3z += coeff
                               * eps_self(nz, info->nkz[0], info->recipbox[ZZ], info->pme_order[0], x[ci]);
                }
            }
        }

        clear_rvec(tmpvec2);

        svmul(e_rec3x, info->recipbox[XX], tmpvec);
        rvec_inc(tmpvec2, tmpvec);
        svmul(e_rec3y, info->recipbox[YY], tmpvec);
        rvec_inc(tmpvec2, tmpvec);
        svmul(e_rec3z, info->recipbox[ZZ], tmpvec);
        rvec_inc(tmpvec2, tmpvec);

        e_rec3 += q[ci] * q[ci] * q[ci] * q[ci] * norm2(tmpvec2)
                  / (xtot * M_PI * info->volume * M_PI * info->volume);
        if (mpiComm.isMainRank())
        {
            fprintf(stderr, "\rCalculating reciprocal error part 2 ... %3.0f%%", 100.0 * (i + 1) / stoplocal);
            std::fflush(stderr);
        }
    }

    if (mpiComm.isMainRank())
    {
        fprintf(stderr, "\n");
    }

#if GMX_LIB_MPI
#    ifdef TAKETIME
    if (mpiComm.isMainRank())
    {
        t1 = MPI_Wtime() - t0;
        fprintf(fp_out, "Recip. err. est. took   : %lf s\n", t1);
    }
#    endif
#endif

#ifdef DEBUG
    if (mpiComm.size() > 1)
    {
        fprintf(stderr, "Rank %3d: nx=[%3d...%3d]  e_rec3=%e\n", mpiComm.rank(), startlocal, stoplocal, e_rec3);
    }
#endif

    if (mpiComm.size() > 1)
    {
        mpiComm.sumReduce(1, &e_rec1);
        mpiComm.sumReduce(1, &e_rec2);
        mpiComm.sumReduce(1, &e_rec3);
    }

    /* e_rec1*=8.0 * q2_all / info->volume / info->volume / nr ;
       e_rec2*=  q2_all / M_PI / M_PI / info->volume / info->volume / nr ;
       e_rec3/= M_PI * M_PI * info->volume * info->volume * nr ;
     */
    e_rec = std::sqrt(e_rec1 + e_rec2 + e_rec3);


    return gmx::c_one4PiEps0 * e_rec;
}


/* Allocate memory for the PmeErrorInputs struct: */
static void create_info(PmeErrorInputs* info)
{
    snew(info->fac, info->n_entries);
    snew(info->rcoulomb, info->n_entries);
    snew(info->rvdw, info->n_entries);
    snew(info->nkx, info->n_entries);
    snew(info->nky, info->n_entries);
    snew(info->nkz, info->n_entries);
    snew(info->fourier_sp, info->n_entries);
    snew(info->ewald_rtol, info->n_entries);
    snew(info->ewald_beta, info->n_entries);
    snew(info->pme_order, info->n_entries);
    snew(info->fn_out, info->n_entries);
    snew(info->e_dir, info->n_entries);
    snew(info->e_rec, info->n_entries);

    // Keep the static-analyzer happy
    info->volume  = 0;
    info->q2all   = 0;
    info->q2allnr = 0;
}


/* Allocate and fill an array with coordinates and charges,
 * returns the number of charges found
 */
static int prepare_x_q(real* q[], rvec* x[], const gmx_mtop_t* mtop, const rvec x_orig[], const gmx::MpiComm& mpiComm)
{
    int nq = 0; /* number of charged particles, keep static-analyzer happy by zeroing here */


    if (mpiComm.isMainRank())
    {
        snew(*q, mtop->natoms);
        snew(*x, mtop->natoms);
        nq = 0;

        for (const AtomProxy atomP : AtomRange(*mtop))
        {
            const t_atom& local = atomP.atom();
            int           i     = atomP.globalAtomNumber();
            if (is_charge(local.q))
            {
                (*q)[nq]     = local.q;
                (*x)[nq][XX] = x_orig[i][XX];
                (*x)[nq][YY] = x_orig[i][YY];
                (*x)[nq][ZZ] = x_orig[i][ZZ];
                nq++;
            }
        }
        /* Give back some unneeded memory */
        srenew(*q, nq);
        srenew(*x, nq);
    }
    /* Broadcast x and q in the parallel case */
    if (mpiComm.size() > 1)
    {
        /* Transfer the number of charges */
        block_bc(mpiComm.comm(), nq);
        snew_bc(mpiComm.isMainRank(), *x, nq);
        snew_bc(mpiComm.isMainRank(), *q, nq);
        nblock_bc(mpiComm.comm(), nq, *x);
        nblock_bc(mpiComm.comm(), nq, *q);
    }

    return nq;
}


/* Read in the tpr file and save information we need later in info */
static void read_tpr_file(const char*     fn_sim_tpr,
                          PmeErrorInputs* info,
                          t_state*        state,
                          gmx_mtop_t*     mtop,
                          t_inputrec*     ir,
                          real            user_beta,
                          real            fracself)
{
    read_tpx_state(fn_sim_tpr, ir, state, mtop);

    /* The values of the original tpr input file are save in the first
     * place [0] of the arrays */
    info->orig_sim_steps = ir->nsteps;
    info->pme_order[0]   = ir->pme_order;
    info->rcoulomb[0]    = ir->rcoulomb;
    info->rvdw[0]        = ir->rvdw;
    info->nkx[0]         = ir->nkx;
    info->nky[0]         = ir->nky;
    info->nkz[0]         = ir->nkz;
    info->ewald_rtol[0]  = ir->ewald_rtol;
    info->fracself       = fracself;
    if (user_beta > 0)
    {
        info->ewald_beta[0] = user_beta;
    }
    else
    {
        info->ewald_beta[0] = calc_ewaldcoeff_q(info->rcoulomb[0], info->ewald_rtol[0]);
    }

    /* Check if PME was chosen */
    if (!usingPme(ir->coulombtype))
    {
        gmx_fatal(FARGS, "Can only do optimizations for simulations with PME");
    }

    /* Check if rcoulomb == rlist, which is necessary for PME */
    if (!(ir->rcoulomb == ir->rlist))
    {
        gmx_fatal(FARGS, "PME requires rcoulomb (%f) to be equal to rlist (%f).", ir->rcoulomb, ir->rlist);
    }
}


/* Transfer what we need for parallelizing the reciprocal error estimate */
static void bcast_info(PmeErrorInputs* info, const gmx::MpiComm& mpiComm)
{
    nblock_bc(mpiComm.comm(), info->n_entries, info->nkx);
    nblock_bc(mpiComm.comm(), info->n_entries, info->nky);
    nblock_bc(mpiComm.comm(), info->n_entries, info->nkz);
    nblock_bc(mpiComm.comm(), info->n_entries, info->ewald_beta);
    nblock_bc(mpiComm.comm(), info->n_entries, info->pme_order);
    nblock_bc(mpiComm.comm(), info->n_entries, info->e_dir);
    nblock_bc(mpiComm.comm(), info->n_entries, info->e_rec);
    block_bc(mpiComm.comm(), info->volume);
    block_bc(mpiComm.comm(), info->recipbox);
    block_bc(mpiComm.comm(), info->natoms);
    block_bc(mpiComm.comm(), info->fracself);
    block_bc(mpiComm.comm(), info->bTUNE);
    block_bc(mpiComm.comm(), info->q2all);
    block_bc(mpiComm.comm(), info->q2allnr);
}


/* Estimate the error of the SPME Ewald sum. This estimate is based upon
 * a) a homogeneous distribution of the charges
 * b) a total charge of zero.
 */
static void estimate_PME_error(PmeErrorInputs*     info,
                               const t_state*      state,
                               const gmx_mtop_t*   mtop,
                               FILE*               fp_out,
                               gmx_bool            bVerbose,
                               unsigned int        seed,
                               const gmx::MpiComm& mpiComm)
{
    rvec* x     = nullptr; /* The coordinates */
    real* q     = nullptr; /* The charges     */
    real  edir  = 0.0;     /* real space error */
    real  erec  = 0.0;     /* reciprocal space error */
    real  derr  = 0.0;     /* difference of real and reciprocal space error */
    real  derr0 = 0.0;     /* difference of real and reciprocal space error */
    real  beta  = 0.0;     /* splitting parameter beta */
    real  beta0 = 0.0;     /* splitting parameter beta */
    int   ncharges;        /* The number of atoms with charges */
    int   nsamples;        /* The number of samples used for the calculation of the
                            * self-energy error term */
    int i = 0;

    // Help humans and analyzers understand that the main rank has
    // a valid pointer
    GMX_RELEASE_ASSERT(mpiComm.isMainRank() == (fp_out != nullptr), "Inconsistent file pointer");
    if (mpiComm.isMainRank())
    {
        fprintf(fp_out, "\n--- PME ERROR ESTIMATE ---\n");
    }

    /* Prepare an x and q array with only the charged atoms */
    ncharges = prepare_x_q(&q, &x, mtop, state->x.rvec_array(), mpiComm);
    if (mpiComm.isMainRank())
    {
        calc_q2all(mtop, &(info->q2all), &(info->q2allnr));
        info->ewald_rtol[0] = std::erfc(info->rcoulomb[0] * info->ewald_beta[0]);
        /* Write some info to log file */
        fprintf(fp_out, "Box volume              : %g nm^3\n", info->volume);
        fprintf(fp_out, "Number of charged atoms : %d (total atoms %d)\n", ncharges, info->natoms);
        fprintf(fp_out, "Coulomb radius          : %g nm\n", info->rcoulomb[0]);
        fprintf(fp_out, "Ewald_rtol              : %g\n", info->ewald_rtol[0]);
        fprintf(fp_out, "Ewald parameter beta    : %g\n", info->ewald_beta[0]);
        fprintf(fp_out, "Interpolation order     : %d\n", info->pme_order[0]);
        fprintf(fp_out, "Fourier grid (nx,ny,nz) : %d x %d x %d\n", info->nkx[0], info->nky[0], info->nkz[0]);
        std::fflush(fp_out);
    }

    if (mpiComm.size() > 1)
    {
        bcast_info(info, mpiComm);
    }


    /* Calculate direct space error */
    info->e_dir[0] = estimate_direct(info);

    /* Calculate reciprocal space error */
    info->e_rec[0] = estimate_reciprocal(info, x, q, ncharges, fp_out, bVerbose, seed, &nsamples, mpiComm);

    if (mpiComm.size() > 1)
    {
        bcast_info(info, mpiComm);
    }

    // Help humans and analyzers understand that the main rank has
    // a valid pointer
    GMX_RELEASE_ASSERT((mpiComm.isMainRank()) == (fp_out != nullptr), "Inconsistent file pointer");
    if (mpiComm.isMainRank())
    {
        fprintf(fp_out, "Direct space error est. : %10.3e kJ/(mol*nm)\n", info->e_dir[0]);
        fprintf(fp_out, "Reciprocal sp. err. est.: %10.3e kJ/(mol*nm)\n", info->e_rec[0]);
        fprintf(fp_out, "Self-energy error term was estimated using %d samples\n", nsamples);
        std::fflush(fp_out);
        fprintf(stderr, "Direct space error est. : %10.3e kJ/(mol*nm)\n", info->e_dir[0]);
        fprintf(stderr, "Reciprocal sp. err. est.: %10.3e kJ/(mol*nm)\n", info->e_rec[0]);
    }

    i = 0;

    if (info->bTUNE)
    {
        if (mpiComm.isMainRank())
        {
            fprintf(stderr, "Starting tuning ...\n");
        }
        edir  = info->e_dir[0];
        erec  = info->e_rec[0];
        derr0 = edir - erec;
        beta0 = info->ewald_beta[0];
        if (derr > 0.0)
        {
            info->ewald_beta[0] += 0.1;
        }
        else
        {
            info->ewald_beta[0] -= 0.1;
        }
        info->e_dir[0] = estimate_direct(info);
        info->e_rec[0] =
                estimate_reciprocal(info, x, q, ncharges, fp_out, bVerbose, seed, &nsamples, mpiComm);

        if (mpiComm.size() > 1)
        {
            bcast_info(info, mpiComm);
        }


        edir = info->e_dir[0];
        erec = info->e_rec[0];
        derr = edir - erec;
        while (std::abs(derr / std::min(erec, edir)) > 1e-4)
        {

            beta = info->ewald_beta[0];
            beta -= derr * (info->ewald_beta[0] - beta0) / (derr - derr0);
            beta0               = info->ewald_beta[0];
            info->ewald_beta[0] = beta;
            derr0               = derr;

            info->e_dir[0] = estimate_direct(info);
            info->e_rec[0] = estimate_reciprocal(
                    info, x, q, ncharges, fp_out, bVerbose, seed, &nsamples, mpiComm);

            if (mpiComm.size() > 1)
            {
                bcast_info(info, mpiComm);
            }

            edir = info->e_dir[0];
            erec = info->e_rec[0];
            derr = edir - erec;

            if (mpiComm.isMainRank())
            {
                i++;
                fprintf(stderr,
                        "difference between real and rec. space error (step %d): %g\n",
                        i,
                        std::abs(derr));
                fprintf(stderr, "old beta: %f\n", beta0);
                fprintf(stderr, "new beta: %f\n", beta);
            }
        }

        info->ewald_rtol[0] = std::erfc(info->rcoulomb[0] * info->ewald_beta[0]);

        if (mpiComm.isMainRank())
        {
            /* Write some info to log file */
            std::fflush(fp_out);
            fprintf(fp_out, "=========  After tuning ========\n");
            fprintf(fp_out, "Direct space error est. : %10.3e kJ/(mol*nm)\n", info->e_dir[0]);
            fprintf(fp_out, "Reciprocal sp. err. est.: %10.3e kJ/(mol*nm)\n", info->e_rec[0]);
            fprintf(stderr, "Direct space error est. : %10.3e kJ/(mol*nm)\n", info->e_dir[0]);
            fprintf(stderr, "Reciprocal sp. err. est.: %10.3e kJ/(mol*nm)\n", info->e_rec[0]);
            fprintf(fp_out, "Ewald_rtol              : %g\n", info->ewald_rtol[0]);
            fprintf(fp_out, "Ewald parameter beta    : %g\n", info->ewald_beta[0]);
            std::fflush(fp_out);
        }
    }
}


int gmx_pme_error(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] estimates the error of the electrostatic forces",
        "if using the sPME algorithm. The flag [TT]-tune[tt] will determine",
        "the splitting parameter such that the error is equally",
        "distributed over the real and reciprocal space part.",
        "The part of the error that stems from self interaction of the particles ",
        "is computationally demanding. However, a good a approximation is to",
        "just use a fraction of the particles for this term which can be",
        "indicated by the flag [TT]-self[tt].[PAR]",
    };

    real           fs        = 0.0; /* 0 indicates: not set by the user */
    real           user_beta = -1.0;
    real           fracself  = 1.0;
    PmeErrorInputs info;
    t_state        state; /* The state from the tpr input file */
    gmx_mtop_t     mtop;  /* The topology from the tpr input file */
    FILE*          fp = nullptr;
    unsigned long  PCA_Flags;
    gmx_bool       bTUNE    = FALSE;
    gmx_bool       bVerbose = FALSE;
    int            seed     = 0;


    static t_filenm fnm[] = { { efTPR, "-s", nullptr, ffREAD },
                              { efOUT, "-o", "error", ffWRITE },
                              { efTPR, "-so", "tuned", ffOPTWR } };

    gmx_output_env_t* oenv = nullptr;

    t_pargs pa[] = {
        { "-beta",
          FALSE,
          etREAL,
          { &user_beta },
          "If positive, overwrite ewald_beta from [REF].tpr[ref] file with this value" },
        { "-tune",
          FALSE,
          etBOOL,
          { &bTUNE },
          "Tune the splitting parameter such that the error is equally distributed between "
          "real and reciprocal space" },
        { "-self",
          FALSE,
          etREAL,
          { &fracself },
          "If between 0.0 and 1.0, determine self interaction error from just this "
          "fraction of the charged particles" },
        { "-seed",
          FALSE,
          etINT,
          { &seed },
          "Random number seed used for Monte Carlo algorithm when [TT]-self[tt] is set to "
          "a value between 0.0 and 1.0" },
        { "-v", FALSE, etBOOL, { &bVerbose }, "Be loud and noisy" }
    };


#define NFILE asize(fnm)

    const gmx::MpiComm mpiComm(MPI_COMM_WORLD);
    PCA_Flags = PCA_NOEXIT_ON_ARGS;

    if (!parse_common_args(
                &argc, argv, PCA_Flags, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    if (!bTUNE)
    {
        bTUNE = opt2bSet("-so", NFILE, fnm);
    }

    info.n_entries = 1;

    /* Allocate memory for the inputinfo struct: */
    create_info(&info);
    info.fourier_sp[0] = fs;

    t_inputrec ir;
    if (mpiComm.isMainRank())
    {
        read_tpr_file(opt2fn("-s", NFILE, fnm), &info, &state, &mtop, &ir, user_beta, fracself);
        /* Open logfile for reading */
        fp = std::fopen(opt2fn("-o", NFILE, fnm), "w");

        /* Determine the volume of the simulation box */
        info.volume = det(state.box);
        calc_recipbox(state.box, info.recipbox);
        info.natoms = mtop.natoms;
        info.bTUNE  = bTUNE;
    }

    /* Check consistency if the user provided fourierspacing */
    if (fs > 0 && mpiComm.isMainRank())
    {
        /* Recalculate the grid dimensions using fourierspacing from user input */
        info.nkx[0] = 0;
        info.nky[0] = 0;
        info.nkz[0] = 0;
        calcFftGrid(stdout,
                    state.box,
                    info.fourier_sp[0],
                    minimalPmeGridSize(info.pme_order[0]),
                    &(info.nkx[0]),
                    &(info.nky[0]),
                    &(info.nkz[0]));
        if ((ir.nkx != info.nkx[0]) || (ir.nky != info.nky[0]) || (ir.nkz != info.nkz[0]))
        {
            gmx_fatal(FARGS,
                      "Wrong fourierspacing %f nm, input file grid = %d x %d x %d, computed grid = "
                      "%d x %d x %d",
                      fs,
                      ir.nkx,
                      ir.nky,
                      ir.nkz,
                      info.nkx[0],
                      info.nky[0],
                      info.nkz[0]);
        }
    }

    /* Estimate (S)PME force error */

    if (mpiComm.size() > 1)
    {
        bcast_info(&info, mpiComm);
    }

    /* Get an error estimate of the input tpr file and do some tuning if requested */
    estimate_PME_error(&info, &state, &mtop, fp, bVerbose, seed, mpiComm);

    if (mpiComm.isMainRank())
    {
        /* Write out optimized tpr file if requested */
        if (opt2bSet("-so", NFILE, fnm) || bTUNE)
        {
            ir.ewald_rtol = info.ewald_rtol[0];
            write_tpx_state(opt2fn("-so", NFILE, fnm), &ir, &state, mtop);
        }
    }
    if (fp)
    {
        please_cite(fp, "Wang2010");
        std::fclose(fp);
    }

    return 0;
}
