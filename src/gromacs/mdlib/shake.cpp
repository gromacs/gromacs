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
/*! \internal \file
 * \brief Defines SHAKE code.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "shake.h"

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/splitter.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"

namespace gmx
{

typedef struct
{
    int iatom[3];
    int blocknr;
} t_sortblock;

//! Compares sort blocks.
static int pcomp(const void* p1, const void* p2)
{
    int                db;
    int                min1, min2, max1, max2;
    const t_sortblock* a1 = reinterpret_cast<const t_sortblock*>(p1);
    const t_sortblock* a2 = reinterpret_cast<const t_sortblock*>(p2);

    db = a1->blocknr - a2->blocknr;

    if (db != 0)
    {
        return db;
    }

    min1 = std::min(a1->iatom[1], a1->iatom[2]);
    max1 = std::max(a1->iatom[1], a1->iatom[2]);
    min2 = std::min(a2->iatom[1], a2->iatom[2]);
    max2 = std::max(a2->iatom[1], a2->iatom[2]);

    if (min1 == min2)
    {
        return max1 - max2;
    }
    else
    {
        return min1 - min2;
    }
}

//! Prints sortblocks
static void pr_sortblock(FILE* fp, const char* title, gmx::ArrayRef<const t_sortblock> sb)
{
    fprintf(fp, "%s\n", title);
    for (gmx::Index i = 0; i < sb.ssize(); i++)
    {
        fprintf(fp,
                "i: %5td, iatom: (%5d %5d %5d), blocknr: %5d\n",
                i,
                sb[i].iatom[0],
                sb[i].iatom[1],
                sb[i].iatom[2],
                sb[i].blocknr);
    }
}

//! Reallocates a vector.
static void resizeLagrangianData(shakedata* shaked, int ncons)
{
    shaked->scaled_lagrange_multiplier.resize(ncons);
}

void make_shake_sblock_serial(shakedata* shaked, InteractionDefinitions* idef, const int numAtoms)
{
    int bstart, bnr;

    /* Since we are processing the local topology,
     * the F_CONSTRNC ilist has been concatenated to the F_CONSTR ilist.
     */
    const int ncons = idef->il[F_CONSTR].size() / 3;

    gmx::ListOfLists<int> sblocks = gen_sblocks(nullptr, numAtoms, *idef, false);

    /*
       bstart=(idef->nodeid > 0) ? blocks->multinr[idef->nodeid-1] : 0;
       nblocks=blocks->multinr[idef->nodeid] - bstart;
     */
    bstart = 0;
    if (debug)
    {
        fprintf(debug, "ncons: %d, bstart: %d, nblocks: %td\n", ncons, bstart, sblocks.ssize());
    }

    /* Calculate block number for each atom */
    std::vector<int> inv_sblock = make_invblock(sblocks, numAtoms);

    /* Store the block number in temp array and
     * sort the constraints in order of the sblock number
     * and the atom numbers, really sorting a segment of the array!
     */
    gmx::ArrayRef<int>       iatom = idef->il[F_CONSTR].iatoms;
    std::vector<t_sortblock> sb(ncons);
    for (int i = 0; i < ncons; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            sb[i].iatom[m] = iatom[i * 3 + m];
        }
        sb[i].blocknr = inv_sblock[iatom[i * 3 + 1]];
    }

    /* Now sort the blocks */
    if (debug)
    {
        pr_sortblock(debug, "Before sorting", sb);
        fprintf(debug, "Going to sort constraints\n");
    }

    std::qsort(sb.data(), gmx::ssize(sb), sizeof(sb[0]), pcomp);

    if (debug)
    {
        pr_sortblock(debug, "After sorting", sb);
    }

    for (int i = 0; i < ncons; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            iatom[i * 3 + m] = sb[i].iatom[m];
        }
    }

    shaked->sblock.clear();
    bnr = -2;
    for (int i = 0; i < ncons; i++)
    {
        if (sb[i].blocknr != bnr)
        {
            bnr = sb[i].blocknr;
            shaked->sblock.push_back(3 * i);
        }
    }
    /* Last block... */
    shaked->sblock.push_back(3 * ncons);

    resizeLagrangianData(shaked, ncons);
}

void make_shake_sblock_dd(shakedata* shaked, const InteractionList& ilcon)
{
    int ncons, c, cg;

    ncons            = ilcon.size() / 3;
    const int* iatom = ilcon.iatoms.data();
    shaked->sblock.clear();
    cg = 0;
    for (c = 0; c < ncons; c++)
    {
        if (c == 0 || iatom[1] >= cg + 1)
        {
            shaked->sblock.push_back(3 * c);
            while (iatom[1] >= cg + 1)
            {
                cg++;
            }
        }
        iatom += 3;
    }
    shaked->sblock.push_back(3 * ncons);
    resizeLagrangianData(shaked, ncons);
}

/*! \brief Inner kernel for SHAKE constraints
 *
 * Original implementation from R.C. van Schaik and W.F. van Gunsteren
 * (ETH Zuerich, June 1992), adapted for GROMACS by David van der
 * Spoel November 1992.
 *
 * The algorithm here is based section five of Ryckaert, Ciccotti and
 * Berendsen, J Comp Phys, 23, 327, 1977.
 *
 * \param[in]    iatom                         Mini-topology of triplets of constraint type (unused
 *                                             in this function) and indices of two atoms involved
 * \param[in]    ncon                          Number of constraints
 * \param[out]   nnit                          Number of iterations performed
 * \param[in]    maxnit                        Maximum number of iterations permitted
 * \param[in]    constraint_distance_squared   The objective value for each constraint
 * \param[inout] positions                     The initial (and final) values of the positions
 *                                             of all atoms
 * \param[in]    pbc                           PBC information
 * \param[in]    initial_displacements         The initial displacements of each constraint
 * \param[in]    half_of_reduced_mass          Half of the reduced mass for each constraint
 * \param[in]    omega                         SHAKE over-relaxation factor (set non-1.0 by
 *                                             using shake-sor=yes in the .mdp,
 *                                             but there is no documentation anywhere)
 * \param[in]    invmass                       Inverse mass of each atom
 * \param[in]    distance_squared_tolerance    Multiplicative tolerance on the difference in the
 *                                             square of the constrained distance (see code)
 * \param[out]   scaled_lagrange_multiplier    Scaled Lagrange multiplier for each constraint
 *                                             (-2 * eta from p. 336 of the paper, divided by
 *                                             the constraint distance)
 * \param[out]   nerror                        Zero upon success, returns one more than the index of
 *                                             the problematic constraint if the input was malformed
 *
 * \todo Make SHAKE use better data structures, in particular for iatom. */
void cshake(const int            iatom[],
            int                  ncon,
            int*                 nnit,
            int                  maxnit,
            ArrayRef<const real> constraint_distance_squared,
            ArrayRef<RVec>       positions,
            const t_pbc*         pbc,
            ArrayRef<const RVec> initial_displacements,
            ArrayRef<const real> half_of_reduced_mass,
            real                 omega,
            ArrayRef<const real> invmass,
            ArrayRef<const real> distance_squared_tolerance,
            ArrayRef<real>       scaled_lagrange_multiplier,
            int*                 nerror)
{
    /* default should be increased! MRS 8/4/2009 */
    const real mytol = 1e-10;

    // TODO nconv is used solely as a boolean, so we should write the
    // code like that
    int error = 0;
    int nconv = 1;
    int nit;
    for (nit = 0; (nit < maxnit) && (nconv != 0) && (error == 0); nit++)
    {
        nconv = 0;
        for (int ll = 0; (ll < ncon) && (error == 0); ll++)
        {
            const int  l3   = 3 * ll;
            const real rijx = initial_displacements[ll][XX];
            const real rijy = initial_displacements[ll][YY];
            const real rijz = initial_displacements[ll][ZZ];
            const int  i    = iatom[l3 + 1];
            const int  j    = iatom[l3 + 2];

            /* Compute r prime between atoms i and j, which is the
               displacement *before* this update stage */
            rvec r_prime;
            if (pbc)
            {
                pbc_dx(pbc, positions[i], positions[j], r_prime);
            }
            else
            {
                rvec_sub(positions[i], positions[j], r_prime);
            }
            const real r_prime_squared                = norm2(r_prime);
            const real constraint_distance_squared_ll = constraint_distance_squared[ll];
            const real diff = constraint_distance_squared_ll - r_prime_squared;

            /* iconvf is less than 1 when the error is smaller than a bound */
            const real iconvf = std::abs(diff) * distance_squared_tolerance[ll];

            if (iconvf > 1.0_real)
            {
                nconv = static_cast<int>(iconvf);
                const real r_dot_r_prime =
                        (rijx * r_prime[XX] + rijy * r_prime[YY] + rijz * r_prime[ZZ]);

                if (r_dot_r_prime < constraint_distance_squared_ll * mytol)
                {
                    error = ll + 1;
                }
                else
                {
                    /* The next line solves equation 5.6 (neglecting
                       the term in g^2), for g */
                    real scaled_lagrange_multiplier_ll =
                            omega * diff * half_of_reduced_mass[ll] / r_dot_r_prime;
                    scaled_lagrange_multiplier[ll] += scaled_lagrange_multiplier_ll;
                    const real xh = rijx * scaled_lagrange_multiplier_ll;
                    const real yh = rijy * scaled_lagrange_multiplier_ll;
                    const real zh = rijz * scaled_lagrange_multiplier_ll;
                    const real im = invmass[i];
                    const real jm = invmass[j];
                    positions[i][XX] += xh * im;
                    positions[i][YY] += yh * im;
                    positions[i][ZZ] += zh * im;
                    positions[j][XX] -= xh * jm;
                    positions[j][YY] -= yh * jm;
                    positions[j][ZZ] -= zh * jm;
                }
            }
        }
    }
    *nnit   = nit;
    *nerror = error;
}

//! Implements RATTLE (ie. SHAKE for velocity verlet integrators)
static void crattle(const int            iatom[],
                    int                  ncon,
                    int*                 nnit,
                    int                  maxnit,
                    ArrayRef<const real> constraint_distance_squared,
                    ArrayRef<RVec>       vp,
                    ArrayRef<const RVec> rij,
                    ArrayRef<const real> m2,
                    real                 omega,
                    ArrayRef<const real> invmass,
                    ArrayRef<const real> distance_squared_tolerance,
                    ArrayRef<real>       scaled_lagrange_multiplier,
                    int*                 nerror,
                    real                 invdt)
{
    /*
     *     r.c. van schaik and w.f. van gunsteren
     *     eth zuerich
     *     june 1992
     *     Adapted for use with Gromacs by David van der Spoel november 92 and later.
     *     rattle added by M.R. Shirts, April 2004, from code written by Jay Ponder in TINKER
     *     second part of rattle algorithm
     */

    // TODO nconv is used solely as a boolean, so we should write the
    // code like that
    int error = 0;
    int nconv = 1;
    int nit;
    for (nit = 0; (nit < maxnit) && (nconv != 0) && (error == 0); nit++)
    {
        nconv = 0;
        for (int ll = 0; (ll < ncon) && (error == 0); ll++)
        {
            const int  l3   = 3 * ll;
            const real rijx = rij[ll][XX];
            const real rijy = rij[ll][YY];
            const real rijz = rij[ll][ZZ];
            const int  i    = iatom[l3 + 1];
            const int  j    = iatom[l3 + 2];
            rvec       v;
            rvec_sub(vp[i], vp[j], v);

            const real vpijd                          = v[XX] * rijx + v[YY] * rijy + v[ZZ] * rijz;
            const real constraint_distance_squared_ll = constraint_distance_squared[ll];

            /* iconv is zero when the error is smaller than a bound */
            const real iconvf = std::fabs(vpijd) * (distance_squared_tolerance[ll] / invdt);

            if (iconvf > 1.0_real)
            {
                nconv           = static_cast<int>(iconvf);
                const real fac  = omega * 2.0_real * m2[ll] / constraint_distance_squared_ll;
                const real acor = -fac * vpijd;
                scaled_lagrange_multiplier[ll] += acor;
                const real xh = rijx * acor;
                const real yh = rijy * acor;
                const real zh = rijz * acor;

                const real im = invmass[i];
                const real jm = invmass[j];

                vp[i][XX] += xh * im;
                vp[i][YY] += yh * im;
                vp[i][ZZ] += zh * im;
                vp[j][XX] -= xh * jm;
                vp[j][YY] -= yh * jm;
                vp[j][ZZ] -= zh * jm;
            }
        }
    }
    *nnit   = nit;
    *nerror = error;
}

//! Applies SHAKE
static int vec_shakef(FILE*                     fplog,
                      shakedata*                shaked,
                      ArrayRef<const real>      invmass,
                      int                       ncon,
                      ArrayRef<const t_iparams> ip,
                      const int*                iatom,
                      real                      tol,
                      ArrayRef<const RVec>      x,
                      ArrayRef<RVec>            prime,
                      const t_pbc*              pbc,
                      real                      omega,
                      bool                      bFEP,
                      real                      lambda,
                      ArrayRef<real>            scaled_lagrange_multiplier,
                      real                      invdt,
                      ArrayRef<RVec>            v,
                      bool                      bCalcVir,
                      tensor                    vir_r_m_dr,
                      ConstraintVariable        econq)
{
    int  maxnit = 1000;
    int  nit    = 0, ll, i, j, d, d2, type;
    real L1;
    int  error = 0;
    real constraint_distance;

    shaked->rij.resize(ncon);
    shaked->half_of_reduced_mass.resize(ncon);
    shaked->distance_squared_tolerance.resize(ncon);
    shaked->constraint_distance_squared.resize(ncon);

    ArrayRef<RVec> rij                         = shaked->rij;
    ArrayRef<real> half_of_reduced_mass        = shaked->half_of_reduced_mass;
    ArrayRef<real> distance_squared_tolerance  = shaked->distance_squared_tolerance;
    ArrayRef<real> constraint_distance_squared = shaked->constraint_distance_squared;

    L1            = 1.0_real - lambda;
    const int* ia = iatom;
    for (ll = 0; (ll < ncon); ll++, ia += 3)
    {
        type = ia[0];
        i    = ia[1];
        j    = ia[2];

        if (pbc)
        {
            pbc_dx(pbc, x[i], x[j], rij[ll]);
        }
        else
        {
            rvec_sub(x[i], x[j], rij[ll]);
        }
        const real mm            = 2.0_real * (invmass[i] + invmass[j]);
        half_of_reduced_mass[ll] = 1.0_real / mm;
        if (bFEP)
        {
            constraint_distance = L1 * ip[type].constr.dA + lambda * ip[type].constr.dB;
        }
        else
        {
            constraint_distance = ip[type].constr.dA;
        }
        constraint_distance_squared[ll] = gmx::square(constraint_distance);
        distance_squared_tolerance[ll]  = 0.5 / (constraint_distance_squared[ll] * tol);
    }

    switch (econq)
    {
        case ConstraintVariable::Positions:
            cshake(iatom,
                   ncon,
                   &nit,
                   maxnit,
                   constraint_distance_squared,
                   prime,
                   pbc,
                   rij,
                   half_of_reduced_mass,
                   omega,
                   invmass,
                   distance_squared_tolerance,
                   scaled_lagrange_multiplier,
                   &error);
            break;
        case ConstraintVariable::Velocities:
            crattle(iatom,
                    ncon,
                    &nit,
                    maxnit,
                    constraint_distance_squared,
                    prime,
                    rij,
                    half_of_reduced_mass,
                    omega,
                    invmass,
                    distance_squared_tolerance,
                    scaled_lagrange_multiplier,
                    &error,
                    invdt);
            break;
        default: gmx_incons("Unknown constraint quantity for SHAKE");
    }

    if (nit >= maxnit)
    {
        if (fplog)
        {
            fprintf(fplog, "Shake did not converge in %d steps\n", maxnit);
        }
        fprintf(stderr, "Shake did not converge in %d steps\n", maxnit);
        nit = 0;
    }
    else if (error != 0)
    {
        if (fplog)
        {
            fprintf(fplog,
                    "Inner product between old and new vector <= 0.0!\n"
                    "constraint #%d atoms %d and %d\n",
                    error - 1,
                    iatom[3 * (error - 1) + 1] + 1,
                    iatom[3 * (error - 1) + 2] + 1);
        }
        fprintf(stderr,
                "Inner product between old and new vector <= 0.0!\n"
                "constraint #%d atoms %d and %d\n",
                error - 1,
                iatom[3 * (error - 1) + 1] + 1,
                iatom[3 * (error - 1) + 2] + 1);
        nit = 0;
    }

    /* Constraint virial and correct the Lagrange multipliers for the length */

    ia = iatom;

    for (ll = 0; (ll < ncon); ll++, ia += 3)
    {
        type = ia[0];
        i    = ia[1];
        j    = ia[2];

        if ((econq == ConstraintVariable::Positions) && !v.empty())
        {
            /* Correct the velocities */
            real mm = scaled_lagrange_multiplier[ll] * invmass[i] * invdt;
            for (d = 0; d < DIM; d++)
            {
                v[ia[1]][d] += mm * rij[ll][d];
            }
            mm = scaled_lagrange_multiplier[ll] * invmass[j] * invdt;
            for (d = 0; d < DIM; d++)
            {
                v[ia[2]][d] -= mm * rij[ll][d];
            }
            /* 16 flops */
        }

        /* constraint virial */
        if (bCalcVir)
        {
            const real mm = scaled_lagrange_multiplier[ll];
            for (d = 0; d < DIM; d++)
            {
                const real tmp = mm * rij[ll][d];
                for (d2 = 0; d2 < DIM; d2++)
                {
                    vir_r_m_dr[d][d2] -= tmp * rij[ll][d2];
                }
            }
            /* 21 flops */
        }

        /* cshake and crattle produce Lagrange multipliers scaled by
           the reciprocal of the constraint length, so fix that */
        if (bFEP)
        {
            constraint_distance = L1 * ip[type].constr.dA + lambda * ip[type].constr.dB;
        }
        else
        {
            constraint_distance = ip[type].constr.dA;
        }
        scaled_lagrange_multiplier[ll] *= constraint_distance;
    }

    return nit;
}

//! Check that constraints are satisfied.
static void check_cons(FILE*                     log,
                       int                       nc,
                       ArrayRef<const RVec>      x,
                       ArrayRef<const RVec>      prime,
                       ArrayRef<const RVec>      v,
                       const t_pbc*              pbc,
                       ArrayRef<const t_iparams> ip,
                       const int*                iatom,
                       ArrayRef<const real>      invmass,
                       ConstraintVariable        econq)
{
    int  ai, aj;
    int  i;
    real d, dp;
    rvec dx, dv;

    GMX_ASSERT(!v.empty(), "Input has to be non-null");
    fprintf(log, "    i     mi      j     mj      before       after   should be\n");
    const int* ia = iatom;
    for (i = 0; (i < nc); i++, ia += 3)
    {
        ai = ia[1];
        aj = ia[2];
        rvec_sub(x[ai], x[aj], dx);
        d = norm(dx);

        switch (econq)
        {
            case ConstraintVariable::Positions:
                if (pbc)
                {
                    pbc_dx(pbc, prime[ai], prime[aj], dx);
                }
                else
                {
                    rvec_sub(prime[ai], prime[aj], dx);
                }
                dp = norm(dx);
                fprintf(log,
                        "%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n",
                        ai + 1,
                        1.0 / invmass[ai],
                        aj + 1,
                        1.0 / invmass[aj],
                        d,
                        dp,
                        ip[ia[0]].constr.dA);
                break;
            case ConstraintVariable::Velocities:
                rvec_sub(v[ai], v[aj], dv);
                d = iprod(dx, dv);
                rvec_sub(prime[ai], prime[aj], dv);
                dp = iprod(dx, dv);
                fprintf(log,
                        "%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n",
                        ai + 1,
                        1.0 / invmass[ai],
                        aj + 1,
                        1.0 / invmass[aj],
                        d,
                        dp,
                        0.);
                break;
            default: gmx_incons("Unknown constraint quantity for SHAKE");
        }
    }
}

//! Applies SHAKE.
static bool bshakef(FILE*                         log,
                    shakedata*                    shaked,
                    ArrayRef<const real>          invmass,
                    const InteractionDefinitions& idef,
                    const t_inputrec&             ir,
                    ArrayRef<const RVec>          x_s,
                    ArrayRef<RVec>                prime,
                    const t_pbc*                  pbc,
                    t_nrnb*                       nrnb,
                    real                          lambda,
                    real*                         dvdlambda,
                    real                          invdt,
                    ArrayRef<RVec>                v,
                    bool                          bCalcVir,
                    tensor                        vir_r_m_dr,
                    bool                          bDumpOnError,
                    ConstraintVariable            econq)
{
    real dt_2, dvdl;
    int  i, n0, ncon, blen, type, ll;
    int  tnit = 0, trij = 0;

    ncon = idef.il[F_CONSTR].size() / 3;

    for (ll = 0; ll < ncon; ll++)
    {
        shaked->scaled_lagrange_multiplier[ll] = 0;
    }

    // TODO Rewrite this block so that it is obvious that i, iatoms
    // and lam are all iteration variables. Is this easier if the
    // sblock data structure is organized differently?
    const int*     iatoms = &(idef.il[F_CONSTR].iatoms[shaked->sblock[0]]);
    ArrayRef<real> lam    = shaked->scaled_lagrange_multiplier;
    for (i = 0; (i < shaked->numShakeBlocks());)
    {
        blen = (shaked->sblock[i + 1] - shaked->sblock[i]);
        blen /= 3;
        n0 = vec_shakef(log,
                        shaked,
                        invmass,
                        blen,
                        idef.iparams,
                        iatoms,
                        ir.shake_tol,
                        x_s,
                        prime,
                        pbc,
                        shaked->omega,
                        ir.efep != FreeEnergyPerturbationType::No,
                        lambda,
                        lam,
                        invdt,
                        v,
                        bCalcVir,
                        vir_r_m_dr,
                        econq);

        if (n0 == 0)
        {
            if (bDumpOnError && log)
            {
                {
                    check_cons(log, blen, x_s, prime, v, pbc, idef.iparams, iatoms, invmass, econq);
                }
            }
            return FALSE;
        }
        tnit += n0 * blen;
        trij += blen;
        iatoms += 3 * blen; /* Increment pointer! */
        lam = lam.subArray(blen, lam.ssize() - blen);
        i++;
    }
    /* only for position part? */
    if (econq == ConstraintVariable::Positions)
    {
        if (ir.efep != FreeEnergyPerturbationType::No)
        {
            ArrayRef<const t_iparams> iparams = idef.iparams;

            /* TODO This should probably use invdt, so that sd integrator scaling works properly */
            dt_2 = 1 / gmx::square(ir.delta_t);
            dvdl = 0;
            for (ll = 0; ll < ncon; ll++)
            {
                type = idef.il[F_CONSTR].iatoms[3 * ll];

                /* Per equations in the manual, dv/dl = -2 \sum_ll lagrangian_ll * r_ll * (d_B - d_A) */
                /* The vector scaled_lagrange_multiplier[ll] contains the value -2 r_ll eta_ll
                   (eta_ll is the estimate of the Langrangian, definition on page 336 of Ryckaert et
                   al 1977), so the pre-factors are already present. */
                const real bondA = iparams[type].constr.dA;
                const real bondB = iparams[type].constr.dB;
                dvdl += shaked->scaled_lagrange_multiplier[ll] * dt_2 * (bondB - bondA);
            }
            *dvdlambda += dvdl;
        }
    }
    if (ir.bShakeSOR)
    {
        if (tnit > shaked->gamma)
        {
            shaked->delta *= -0.5;
        }
        shaked->omega += shaked->delta;
        shaked->gamma = tnit;
    }
    inc_nrnb(nrnb, eNR_SHAKE, tnit);
    inc_nrnb(nrnb, eNR_SHAKE_RIJ, trij);
    if (!v.empty())
    {
        inc_nrnb(nrnb, eNR_CONSTR_V, trij * 2);
    }
    if (bCalcVir)
    {
        inc_nrnb(nrnb, eNR_CONSTR_VIR, trij);
    }

    return TRUE;
}

bool constrain_shake(FILE*                         log,
                     shakedata*                    shaked,
                     ArrayRef<const real>          invmass,
                     const InteractionDefinitions& idef,
                     const t_inputrec&             ir,
                     ArrayRef<const RVec>          x_s,
                     ArrayRef<RVec>                xprime,
                     ArrayRef<RVec>                vprime,
                     const t_pbc*                  pbc,
                     t_nrnb*                       nrnb,
                     real                          lambda,
                     real*                         dvdlambda,
                     real                          invdt,
                     ArrayRef<RVec>                v,
                     bool                          bCalcVir,
                     tensor                        vir_r_m_dr,
                     bool                          bDumpOnError,
                     ConstraintVariable            econq)
{
    if (shaked->numShakeBlocks() == 0)
    {
        return true;
    }
    bool bOK;
    switch (econq)
    {
        case (ConstraintVariable::Positions):
            bOK = bshakef(
                    log, shaked, invmass, idef, ir, x_s, xprime, pbc, nrnb, lambda, dvdlambda, invdt, v, bCalcVir, vir_r_m_dr, bDumpOnError, econq);
            break;
        case (ConstraintVariable::Velocities):
            bOK = bshakef(log,
                          shaked,
                          invmass,
                          idef,
                          ir,
                          x_s,
                          vprime,
                          pbc,
                          nrnb,
                          lambda,
                          dvdlambda,
                          invdt,
                          {},
                          bCalcVir,
                          vir_r_m_dr,
                          bDumpOnError,
                          econq);
            break;
        default:
            gmx_fatal(FARGS,
                      "Internal error, SHAKE called for constraining something else than "
                      "coordinates");
    }
    return bOK;
}

} // namespace gmx
