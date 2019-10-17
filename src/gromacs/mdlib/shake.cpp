/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/splitter.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

struct shakedata
{
    rvec* rij;
    real* half_of_reduced_mass;
    real* distance_squared_tolerance;
    real* constraint_distance_squared;
    int   nalloc;
    /* SOR stuff */
    real delta;
    real omega;
    real gamma;
    int  nblocks;       /* The number of SHAKE blocks         */
    int* sblock;        /* The SHAKE blocks                   */
    int  sblock_nalloc; /* The allocation size of sblock      */
    /*! \brief Scaled Lagrange multiplier for each constraint.
     *
     * Value is -2 * eta from p. 336 of the paper, divided by the
     * constraint distance. */
    real* scaled_lagrange_multiplier;
    int   lagr_nalloc; /* The allocation size of scaled_lagrange_multiplier */
};

shakedata* shake_init()
{
    shakedata* d;

    snew(d, 1);

    d->nalloc                      = 0;
    d->rij                         = nullptr;
    d->half_of_reduced_mass        = nullptr;
    d->distance_squared_tolerance  = nullptr;
    d->constraint_distance_squared = nullptr;

    /* SOR initialization */
    d->delta = 0.1;
    d->omega = 1.0;
    d->gamma = 1000000;

    return d;
}

void done_shake(shakedata* d)
{
    sfree(d->rij);
    sfree(d->half_of_reduced_mass);
    sfree(d->distance_squared_tolerance);
    sfree(d->constraint_distance_squared);
    sfree(d->sblock);
    sfree(d->scaled_lagrange_multiplier);
    sfree(d);
}

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
static void pr_sortblock(FILE* fp, const char* title, int nsb, t_sortblock sb[])
{
    int i;

    fprintf(fp, "%s\n", title);
    for (i = 0; (i < nsb); i++)
    {
        fprintf(fp, "i: %5d, iatom: (%5d %5d %5d), blocknr: %5d\n", i, sb[i].iatom[0],
                sb[i].iatom[1], sb[i].iatom[2], sb[i].blocknr);
    }
}

//! Reallocates a vector.
static void resizeLagrangianData(shakedata* shaked, int ncons)
{
    if (ncons > shaked->lagr_nalloc)
    {
        shaked->lagr_nalloc = over_alloc_dd(ncons);
        srenew(shaked->scaled_lagrange_multiplier, shaked->lagr_nalloc);
    }
}

void make_shake_sblock_serial(shakedata* shaked, const t_idef* idef, const t_mdatoms& md)
{
    int          i, j, m, ncons;
    int          bstart, bnr;
    t_blocka     sblocks;
    t_sortblock* sb;
    t_iatom*     iatom;
    int*         inv_sblock;

    /* Since we are processing the local topology,
     * the F_CONSTRNC ilist has been concatenated to the F_CONSTR ilist.
     */
    ncons = idef->il[F_CONSTR].nr / 3;

    init_blocka(&sblocks);
    sfree(sblocks.index); // To solve memory leak
    gen_sblocks(nullptr, 0, md.homenr, idef, &sblocks, FALSE);

    /*
       bstart=(idef->nodeid > 0) ? blocks->multinr[idef->nodeid-1] : 0;
       nblocks=blocks->multinr[idef->nodeid] - bstart;
     */
    bstart          = 0;
    shaked->nblocks = sblocks.nr;
    if (debug)
    {
        fprintf(debug, "ncons: %d, bstart: %d, nblocks: %d\n", ncons, bstart, shaked->nblocks);
    }

    /* Calculate block number for each atom */
    inv_sblock = make_invblocka(&sblocks, md.nr);

    done_blocka(&sblocks);

    /* Store the block number in temp array and
     * sort the constraints in order of the sblock number
     * and the atom numbers, really sorting a segment of the array!
     */
    iatom = idef->il[F_CONSTR].iatoms;
    snew(sb, ncons);
    for (i = 0; (i < ncons); i++, iatom += 3)
    {
        for (m = 0; (m < 3); m++)
        {
            sb[i].iatom[m] = iatom[m];
        }
        sb[i].blocknr = inv_sblock[iatom[1]];
    }

    /* Now sort the blocks */
    if (debug)
    {
        pr_sortblock(debug, "Before sorting", ncons, sb);
        fprintf(debug, "Going to sort constraints\n");
    }

    std::qsort(sb, ncons, sizeof(*sb), pcomp);

    if (debug)
    {
        pr_sortblock(debug, "After sorting", ncons, sb);
    }

    iatom = idef->il[F_CONSTR].iatoms;
    for (i = 0; (i < ncons); i++, iatom += 3)
    {
        for (m = 0; (m < 3); m++)
        {
            iatom[m] = sb[i].iatom[m];
        }
    }

    j = 0;
    snew(shaked->sblock, shaked->nblocks + 1);
    bnr = -2;
    for (i = 0; (i < ncons); i++)
    {
        if (sb[i].blocknr != bnr)
        {
            bnr                 = sb[i].blocknr;
            shaked->sblock[j++] = 3 * i;
        }
    }
    /* Last block... */
    shaked->sblock[j++] = 3 * ncons;

    if (j != (shaked->nblocks + 1))
    {
        fprintf(stderr, "bstart: %d\n", bstart);
        fprintf(stderr, "j: %d, nblocks: %d, ncons: %d\n", j, shaked->nblocks, ncons);
        for (i = 0; (i < ncons); i++)
        {
            fprintf(stderr, "i: %5d  sb[i].blocknr: %5d\n", i, sb[i].blocknr);
        }
        for (j = 0; (j <= shaked->nblocks); j++)
        {
            fprintf(stderr, "sblock[%3d]=%5d\n", j, shaked->sblock[j]);
        }
        gmx_fatal(FARGS,
                  "DEATH HORROR: "
                  "sblocks does not match idef->il[F_CONSTR]");
    }
    sfree(sb);
    sfree(inv_sblock);
    resizeLagrangianData(shaked, ncons);
}

// TODO: Check if this code is useful. It might never be called.
void make_shake_sblock_dd(shakedata* shaked, const t_ilist* ilcon, const gmx_domdec_t* dd)
{
    int      ncons, c, cg;
    t_iatom* iatom;

    if (dd->ncg_home + 1 > shaked->sblock_nalloc)
    {
        shaked->sblock_nalloc = over_alloc_dd(dd->ncg_home + 1);
        srenew(shaked->sblock, shaked->sblock_nalloc);
    }

    ncons           = ilcon->nr / 3;
    iatom           = ilcon->iatoms;
    shaked->nblocks = 0;
    cg              = 0;
    for (c = 0; c < ncons; c++)
    {
        if (c == 0 || iatom[1] >= cg + 1)
        {
            shaked->sblock[shaked->nblocks++] = 3 * c;
            while (iatom[1] >= cg + 1)
            {
                cg++;
            }
        }
        iatom += 3;
    }
    shaked->sblock[shaked->nblocks] = 3 * ncons;
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
 * \param[in]    iatom                         Mini-topology of triples of constraint type (unused in this
 *                                             function) and indices of the two atoms involved
 * \param[in]    ncon                          Number of constraints
 * \param[out]   nnit                          Number of iterations performed
 * \param[in]    maxnit                        Maximum number of iterations permitted
 * \param[in]    constraint_distance_squared   The objective value for each constraint
 * \param[inout] positions                     The initial (and final) values of the positions of all atoms
 * \param[in]    initial_displacements         The initial displacements of each constraint
 * \param[in]    half_of_reduced_mass          Half of the reduced mass for each constraint
 * \param[in]    omega                         SHAKE over-relaxation factor (set non-1.0 by
 *                                             using shake-sor=yes in the .mdp, but there is no documentation anywhere)
 * \param[in]    invmass                       Inverse mass of each atom
 * \param[in]    distance_squared_tolerance    Multiplicative tolerance on the difference in the
 *                                             square of the constrained distance (see code)
 * \param[out]   scaled_lagrange_multiplier    Scaled Lagrange multiplier for each constraint (-2 * eta from p. 336
 *                                             of the paper, divided by the constraint distance)
 * \param[out]   nerror                        Zero upon success, returns one more than the index of the
 *                                             problematic constraint if the input was malformed
 *
 * \todo Make SHAKE use better data structures, in particular for iatom. */
void cshake(const int  iatom[],
            int        ncon,
            int*       nnit,
            int        maxnit,
            const real constraint_distance_squared[],
            real       positions[],
            const real initial_displacements[],
            const real half_of_reduced_mass[],
            real       omega,
            const real invmass[],
            const real distance_squared_tolerance[],
            real       scaled_lagrange_multiplier[],
            int*       nerror)
{
    /* default should be increased! MRS 8/4/2009 */
    const real mytol = 1e-10;

    int  ll, i, j, i3, j3, l3;
    int  ix, iy, iz, jx, jy, jz;
    real r_dot_r_prime;
    real constraint_distance_squared_ll;
    real r_prime_squared;
    real scaled_lagrange_multiplier_ll;
    real r_prime_x, r_prime_y, r_prime_z, diff, im, jm;
    real xh, yh, zh, rijx, rijy, rijz;
    int  nit, error, nconv;
    real iconvf;

    // TODO nconv is used solely as a boolean, so we should write the
    // code like that
    error = 0;
    nconv = 1;
    for (nit = 0; (nit < maxnit) && (nconv != 0) && (error == 0); nit++)
    {
        nconv = 0;
        for (ll = 0; (ll < ncon) && (error == 0); ll++)
        {
            l3   = 3 * ll;
            rijx = initial_displacements[l3 + XX];
            rijy = initial_displacements[l3 + YY];
            rijz = initial_displacements[l3 + ZZ];
            i    = iatom[l3 + 1];
            j    = iatom[l3 + 2];
            i3   = 3 * i;
            j3   = 3 * j;
            ix   = i3 + XX;
            iy   = i3 + YY;
            iz   = i3 + ZZ;
            jx   = j3 + XX;
            jy   = j3 + YY;
            jz   = j3 + ZZ;

            /* Compute r prime between atoms i and j, which is the
               displacement *before* this update stage */
            r_prime_x       = positions[ix] - positions[jx];
            r_prime_y       = positions[iy] - positions[jy];
            r_prime_z       = positions[iz] - positions[jz];
            r_prime_squared = (r_prime_x * r_prime_x + r_prime_y * r_prime_y + r_prime_z * r_prime_z);
            constraint_distance_squared_ll = constraint_distance_squared[ll];
            diff                           = constraint_distance_squared_ll - r_prime_squared;

            /* iconvf is less than 1 when the error is smaller than a bound */
            iconvf = fabs(diff) * distance_squared_tolerance[ll];

            if (iconvf > 1.0)
            {
                nconv         = static_cast<int>(iconvf);
                r_dot_r_prime = (rijx * r_prime_x + rijy * r_prime_y + rijz * r_prime_z);

                if (r_dot_r_prime < constraint_distance_squared_ll * mytol)
                {
                    error = ll + 1;
                }
                else
                {
                    /* The next line solves equation 5.6 (neglecting
                       the term in g^2), for g */
                    scaled_lagrange_multiplier_ll = omega * diff * half_of_reduced_mass[ll] / r_dot_r_prime;
                    scaled_lagrange_multiplier[ll] += scaled_lagrange_multiplier_ll;
                    xh = rijx * scaled_lagrange_multiplier_ll;
                    yh = rijy * scaled_lagrange_multiplier_ll;
                    zh = rijz * scaled_lagrange_multiplier_ll;
                    im = invmass[i];
                    jm = invmass[j];
                    positions[ix] += xh * im;
                    positions[iy] += yh * im;
                    positions[iz] += zh * im;
                    positions[jx] -= xh * jm;
                    positions[jy] -= yh * jm;
                    positions[jz] -= zh * jm;
                }
            }
        }
    }
    *nnit   = nit;
    *nerror = error;
}

//! Implements RATTLE (ie. SHAKE for velocity verlet integrators)
static void crattle(const int  iatom[],
                    int        ncon,
                    int*       nnit,
                    int        maxnit,
                    const real constraint_distance_squared[],
                    real       vp[],
                    const real rij[],
                    const real m2[],
                    real       omega,
                    const real invmass[],
                    const real distance_squared_tolerance[],
                    real       scaled_lagrange_multiplier[],
                    int*       nerror,
                    real       invdt)
{
    /*
     *     r.c. van schaik and w.f. van gunsteren
     *     eth zuerich
     *     june 1992
     *     Adapted for use with Gromacs by David van der Spoel november 92 and later.
     *     rattle added by M.R. Shirts, April 2004, from code written by Jay Ponder in TINKER
     *     second part of rattle algorithm
     */

    int  ll, i, j, i3, j3, l3;
    int  ix, iy, iz, jx, jy, jz;
    real constraint_distance_squared_ll;
    real vpijd, vx, vy, vz, acor, fac, im, jm;
    real xh, yh, zh, rijx, rijy, rijz;
    int  nit, error, nconv;
    real iconvf;

    // TODO nconv is used solely as a boolean, so we should write the
    // code like that
    error = 0;
    nconv = 1;
    for (nit = 0; (nit < maxnit) && (nconv != 0) && (error == 0); nit++)
    {
        nconv = 0;
        for (ll = 0; (ll < ncon) && (error == 0); ll++)
        {
            l3   = 3 * ll;
            rijx = rij[l3 + XX];
            rijy = rij[l3 + YY];
            rijz = rij[l3 + ZZ];
            i    = iatom[l3 + 1];
            j    = iatom[l3 + 2];
            i3   = 3 * i;
            j3   = 3 * j;
            ix   = i3 + XX;
            iy   = i3 + YY;
            iz   = i3 + ZZ;
            jx   = j3 + XX;
            jy   = j3 + YY;
            jz   = j3 + ZZ;
            vx   = vp[ix] - vp[jx];
            vy   = vp[iy] - vp[jy];
            vz   = vp[iz] - vp[jz];

            vpijd                          = vx * rijx + vy * rijy + vz * rijz;
            constraint_distance_squared_ll = constraint_distance_squared[ll];

            /* iconv is zero when the error is smaller than a bound */
            iconvf = fabs(vpijd) * (distance_squared_tolerance[ll] / invdt);

            if (iconvf > 1)
            {
                nconv = static_cast<int>(iconvf);
                fac   = omega * 2.0 * m2[ll] / constraint_distance_squared_ll;
                acor  = -fac * vpijd;
                scaled_lagrange_multiplier[ll] += acor;
                xh = rijx * acor;
                yh = rijy * acor;
                zh = rijz * acor;

                im = invmass[i];
                jm = invmass[j];

                vp[ix] += xh * im;
                vp[iy] += yh * im;
                vp[iz] += zh * im;
                vp[jx] -= xh * jm;
                vp[jy] -= yh * jm;
                vp[jz] -= zh * jm;
            }
        }
    }
    *nnit   = nit;
    *nerror = error;
}

//! Applies SHAKE
static int vec_shakef(FILE*              fplog,
                      shakedata*         shaked,
                      const real         invmass[],
                      int                ncon,
                      t_iparams          ip[],
                      t_iatom*           iatom,
                      real               tol,
                      const rvec         x[],
                      rvec               prime[],
                      real               omega,
                      bool               bFEP,
                      real               lambda,
                      real               scaled_lagrange_multiplier[],
                      real               invdt,
                      rvec*              v,
                      bool               bCalcVir,
                      tensor             vir_r_m_dr,
                      ConstraintVariable econq)
{
    rvec*    rij;
    real *   half_of_reduced_mass, *distance_squared_tolerance, *constraint_distance_squared;
    int      maxnit = 1000;
    int      nit    = 0, ll, i, j, d, d2, type;
    t_iatom* ia;
    real     L1;
    real     mm    = 0., tmp;
    int      error = 0;
    real     constraint_distance;

    if (ncon > shaked->nalloc)
    {
        shaked->nalloc = over_alloc_dd(ncon);
        srenew(shaked->rij, shaked->nalloc);
        srenew(shaked->half_of_reduced_mass, shaked->nalloc);
        srenew(shaked->distance_squared_tolerance, shaked->nalloc);
        srenew(shaked->constraint_distance_squared, shaked->nalloc);
    }
    rij                         = shaked->rij;
    half_of_reduced_mass        = shaked->half_of_reduced_mass;
    distance_squared_tolerance  = shaked->distance_squared_tolerance;
    constraint_distance_squared = shaked->constraint_distance_squared;

    L1 = 1.0 - lambda;
    ia = iatom;
    for (ll = 0; (ll < ncon); ll++, ia += 3)
    {
        type = ia[0];
        i    = ia[1];
        j    = ia[2];

        mm                       = 2.0 * (invmass[i] + invmass[j]);
        rij[ll][XX]              = x[i][XX] - x[j][XX];
        rij[ll][YY]              = x[i][YY] - x[j][YY];
        rij[ll][ZZ]              = x[i][ZZ] - x[j][ZZ];
        half_of_reduced_mass[ll] = 1.0 / mm;
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
            cshake(iatom, ncon, &nit, maxnit, constraint_distance_squared, prime[0], rij[0],
                   half_of_reduced_mass, omega, invmass, distance_squared_tolerance,
                   scaled_lagrange_multiplier, &error);
            break;
        case ConstraintVariable::Velocities:
            crattle(iatom, ncon, &nit, maxnit, constraint_distance_squared, prime[0], rij[0],
                    half_of_reduced_mass, omega, invmass, distance_squared_tolerance,
                    scaled_lagrange_multiplier, &error, invdt);
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
                    error - 1, iatom[3 * (error - 1) + 1] + 1, iatom[3 * (error - 1) + 2] + 1);
        }
        fprintf(stderr,
                "Inner product between old and new vector <= 0.0!\n"
                "constraint #%d atoms %d and %d\n",
                error - 1, iatom[3 * (error - 1) + 1] + 1, iatom[3 * (error - 1) + 2] + 1);
        nit = 0;
    }

    /* Constraint virial and correct the Lagrange multipliers for the length */

    ia = iatom;

    for (ll = 0; (ll < ncon); ll++, ia += 3)
    {
        type = ia[0];
        i    = ia[1];
        j    = ia[2];

        if ((econq == ConstraintVariable::Positions) && v != nullptr)
        {
            /* Correct the velocities */
            mm = scaled_lagrange_multiplier[ll] * invmass[i] * invdt;
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
            mm = scaled_lagrange_multiplier[ll];
            for (d = 0; d < DIM; d++)
            {
                tmp = mm * rij[ll][d];
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
static void check_cons(FILE*              log,
                       int                nc,
                       const rvec         x[],
                       rvec               prime[],
                       rvec               v[],
                       t_iparams          ip[],
                       t_iatom*           iatom,
                       const real         invmass[],
                       ConstraintVariable econq)
{
    t_iatom* ia;
    int      ai, aj;
    int      i;
    real     d, dp;
    rvec     dx, dv;

    GMX_ASSERT(v, "Input has to be non-null");
    fprintf(log, "    i     mi      j     mj      before       after   should be\n");
    ia = iatom;
    for (i = 0; (i < nc); i++, ia += 3)
    {
        ai = ia[1];
        aj = ia[2];
        rvec_sub(x[ai], x[aj], dx);
        d = norm(dx);

        switch (econq)
        {
            case ConstraintVariable::Positions:
                rvec_sub(prime[ai], prime[aj], dx);
                dp = norm(dx);
                fprintf(log, "%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n", ai + 1,
                        1.0 / invmass[ai], aj + 1, 1.0 / invmass[aj], d, dp, ip[ia[0]].constr.dA);
                break;
            case ConstraintVariable::Velocities:
                rvec_sub(v[ai], v[aj], dv);
                d = iprod(dx, dv);
                rvec_sub(prime[ai], prime[aj], dv);
                dp = iprod(dx, dv);
                fprintf(log, "%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n", ai + 1,
                        1.0 / invmass[ai], aj + 1, 1.0 / invmass[aj], d, dp, 0.);
                break;
            default: gmx_incons("Unknown constraint quantity for SHAKE");
        }
    }
}

//! Applies SHAKE.
static bool bshakef(FILE*              log,
                    shakedata*         shaked,
                    const real         invmass[],
                    const t_idef&      idef,
                    const t_inputrec&  ir,
                    const rvec         x_s[],
                    rvec               prime[],
                    t_nrnb*            nrnb,
                    real               lambda,
                    real*              dvdlambda,
                    real               invdt,
                    rvec*              v,
                    bool               bCalcVir,
                    tensor             vir_r_m_dr,
                    bool               bDumpOnError,
                    ConstraintVariable econq)
{
    t_iatom* iatoms;
    real *   lam, dt_2, dvdl;
    int      i, n0, ncon, blen, type, ll;
    int      tnit = 0, trij = 0;

    ncon = idef.il[F_CONSTR].nr / 3;

    for (ll = 0; ll < ncon; ll++)
    {
        shaked->scaled_lagrange_multiplier[ll] = 0;
    }

    // TODO Rewrite this block so that it is obvious that i, iatoms
    // and lam are all iteration variables. Is this easier if the
    // sblock data structure is organized differently?
    iatoms = &(idef.il[F_CONSTR].iatoms[shaked->sblock[0]]);
    lam    = shaked->scaled_lagrange_multiplier;
    for (i = 0; (i < shaked->nblocks);)
    {
        blen = (shaked->sblock[i + 1] - shaked->sblock[i]);
        blen /= 3;
        n0 = vec_shakef(log, shaked, invmass, blen, idef.iparams, iatoms, ir.shake_tol, x_s, prime,
                        shaked->omega, ir.efep != efepNO, lambda, lam, invdt, v, bCalcVir,
                        vir_r_m_dr, econq);

        if (n0 == 0)
        {
            if (bDumpOnError && log)
            {
                {
                    check_cons(log, blen, x_s, prime, v, idef.iparams, iatoms, invmass, econq);
                }
            }
            return FALSE;
        }
        tnit += n0 * blen;
        trij += blen;
        iatoms += 3 * blen; /* Increment pointer! */
        lam += blen;
        i++;
    }
    /* only for position part? */
    if (econq == ConstraintVariable::Positions)
    {
        if (ir.efep != efepNO)
        {
            real bondA, bondB;
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
                bondA = idef.iparams[type].constr.dA;
                bondB = idef.iparams[type].constr.dB;
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
    if (v)
    {
        inc_nrnb(nrnb, eNR_CONSTR_V, trij * 2);
    }
    if (bCalcVir)
    {
        inc_nrnb(nrnb, eNR_CONSTR_VIR, trij);
    }

    return TRUE;
}

bool constrain_shake(FILE*              log,
                     shakedata*         shaked,
                     const real         invmass[],
                     const t_idef&      idef,
                     const t_inputrec&  ir,
                     const rvec         x_s[],
                     rvec               xprime[],
                     rvec               vprime[],
                     t_nrnb*            nrnb,
                     real               lambda,
                     real*              dvdlambda,
                     real               invdt,
                     rvec*              v,
                     bool               bCalcVir,
                     tensor             vir_r_m_dr,
                     bool               bDumpOnError,
                     ConstraintVariable econq)
{
    if (shaked->nblocks == 0)
    {
        return true;
    }
    bool bOK;
    switch (econq)
    {
        case (ConstraintVariable::Positions):
            bOK = bshakef(log, shaked, invmass, idef, ir, x_s, xprime, nrnb, lambda, dvdlambda,
                          invdt, v, bCalcVir, vir_r_m_dr, bDumpOnError, econq);
            break;
        case (ConstraintVariable::Velocities):
            bOK = bshakef(log, shaked, invmass, idef, ir, x_s, vprime, nrnb, lambda, dvdlambda,
                          invdt, nullptr, bCalcVir, vir_r_m_dr, bDumpOnError, econq);
            break;
        default:
            gmx_fatal(FARGS,
                      "Internal error, SHAKE called for constraining something else than "
                      "coordinates");
    }
    return bOK;
}

} // namespace gmx
