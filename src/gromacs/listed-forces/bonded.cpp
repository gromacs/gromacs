/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 *
 * \brief This file defines low-level functions necessary for
 * computing energies and forces for listed interactions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed-forces
 */
#include "gmxpre.h"

#include "bonded.h"

#include "config.h"

#include <cassert>
#include <cmath>

#include <algorithm>

#include "gromacs/listed-forces/pairs.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc-simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "listed-internal.h"
#include "restcbt.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

/*! \brief Mysterious CMAP coefficient matrix */
const int cmap_coeff_matrix[] = {
    1, 0, -3,  2, 0, 0,  0,  0, -3,  0,  9, -6,  2,  0, -6,  4,
    0, 0,  0,  0, 0, 0,  0,  0,  3,  0, -9,  6, -2,  0,  6, -4,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  9, -6,  0,  0, -6,  4,
    0, 0,  3, -2, 0, 0,  0,  0,  0,  0, -9,  6,  0,  0,  6, -4,
    0, 0,  0,  0, 1, 0, -3,  2, -2,  0,  6, -4,  1,  0, -3,  2,
    0, 0,  0,  0, 0, 0,  0,  0, -1,  0,  3, -2,  1,  0, -3,  2,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  2,  0,  0,  3, -2,
    0, 0,  0,  0, 0, 0,  3, -2,  0,  0, -6,  4,  0,  0,  3, -2,
    0, 1, -2,  1, 0, 0,  0,  0,  0, -3,  6, -3,  0,  2, -4,  2,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  3, -6,  3,  0, -2,  4, -2,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  3,  0,  0,  2, -2,
    0, 0, -1,  1, 0, 0,  0,  0,  0,  0,  3, -3,  0,  0, -2,  2,
    0, 0,  0,  0, 0, 1, -2,  1,  0, -2,  4, -2,  0,  1, -2,  1,
    0, 0,  0,  0, 0, 0,  0,  0,  0, -1,  2, -1,  0,  1, -2,  1,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1,
    0, 0,  0,  0, 0, 0, -1,  1,  0,  0,  2, -2,  0,  0, -1,  1
};


/*! \brief Compute dx = xi - xj, modulo PBC if non-NULL
 *
 * \todo This kind of code appears in many places. Consolidate it */
static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

/*! \brief Morse potential bond
 *
 * By Frank Everdij. Three parameters needed:
 *
 * b0 = equilibrium distance in nm
 * be = beta in nm^-1 (actually, it's nu_e*Sqrt(2*pi*pi*mu/D_e))
 * cb = well depth in kJ/mol
 *
 * Note: the potential is referenced to be +cb at infinite separation
 *       and zero at the equilibrium distance!
 */
real morse_bonds(int nbonds,
                 const t_iatom forceatoms[], const t_iparams forceparams[],
                 const rvec x[], rvec4 f[], rvec fshift[],
                 const t_pbc *pbc, const t_graph *g,
                 real lambda, real *dvdlambda,
                 const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                 int gmx_unused *global_atom_index)
{
    const real one = 1.0;
    const real two = 2.0;
    real       dr, dr2, temp, omtemp, cbomtemp, fbond, vbond, fij, vtot;
    real       b0, be, cb, b0A, beA, cbA, b0B, beB, cbB, L1;
    rvec       dx;
    int        i, m, ki, type, ai, aj;
    ivec       dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        b0A   = forceparams[type].morse.b0A;
        beA   = forceparams[type].morse.betaA;
        cbA   = forceparams[type].morse.cbA;

        b0B   = forceparams[type].morse.b0B;
        beB   = forceparams[type].morse.betaB;
        cbB   = forceparams[type].morse.cbB;

        L1 = one-lambda;                            /* 1 */
        b0 = L1*b0A + lambda*b0B;                   /* 3 */
        be = L1*beA + lambda*beB;                   /* 3 */
        cb = L1*cbA + lambda*cbB;                   /* 3 */

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx); /*   3          */
        dr2  = iprod(dx, dx);                       /*   5          */
        dr   = dr2*gmx::invsqrt(dr2);               /*  10          */
        temp = std::exp(-be*(dr-b0));               /*  12          */

        if (temp == one)
        {
            /* bonds are constrainted. This may _not_ include bond constraints if they are lambda dependent */
            *dvdlambda += cbB-cbA;
            continue;
        }

        omtemp    = one-temp;                                                                                        /*   1          */
        cbomtemp  = cb*omtemp;                                                                                       /*   1          */
        vbond     = cbomtemp*omtemp;                                                                                 /*   1          */
        fbond     = -two*be*temp*cbomtemp*gmx::invsqrt(dr2);                                                         /*   9          */
        vtot     += vbond;                                                                                           /*   1          */

        *dvdlambda += (cbB - cbA) * omtemp * omtemp - (2-2*omtemp)*omtemp * cb * ((b0B-b0A)*be - (beB-beA)*(dr-b0)); /* 15 */

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }

        for (m = 0; (m < DIM); m++)                    /*  15          */
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }                                         /*  83 TOTAL    */
    return vtot;
}

//! \cond
real cubic_bonds(int nbonds,
                 const t_iatom forceatoms[], const t_iparams forceparams[],
                 const rvec x[], rvec4 f[], rvec fshift[],
                 const t_pbc *pbc, const t_graph *g,
                 real gmx_unused lambda, real gmx_unused *dvdlambda,
                 const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                 int gmx_unused *global_atom_index)
{
    const real three = 3.0;
    const real two   = 2.0;
    real       kb, b0, kcub;
    real       dr, dr2, dist, kdist, kdist2, fbond, vbond, fij, vtot;
    rvec       dx;
    int        i, m, ki, type, ai, aj;
    ivec       dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        b0   = forceparams[type].cubic.b0;
        kb   = forceparams[type].cubic.kb;
        kcub = forceparams[type].cubic.kcub;

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx);     /*   3          */
        dr2  = iprod(dx, dx);                           /*   5          */

        if (dr2 == 0.0)
        {
            continue;
        }

        dr         = dr2*gmx::invsqrt(dr2);                  /*  10          */
        dist       = dr-b0;
        kdist      = kb*dist;
        kdist2     = kdist*dist;

        vbond      = kdist2 + kcub*kdist2*dist;
        fbond      = -(two*kdist + three*kdist2*kcub)/dr;

        vtot      += vbond;   /* 21 */

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)                    /*  15          */
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }                                         /*  54 TOTAL    */
    return vtot;
}

real FENE_bonds(int nbonds,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const rvec x[], rvec4 f[], rvec fshift[],
                const t_pbc *pbc, const t_graph *g,
                real gmx_unused lambda, real gmx_unused *dvdlambda,
                const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                int *global_atom_index)
{
    const real half = 0.5;
    const real one  = 1.0;
    real       bm, kb;
    real       dr2, bm2, omdr2obm2, fbond, vbond, fij, vtot;
    rvec       dx;
    int        i, m, ki, type, ai, aj;
    ivec       dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        bm   = forceparams[type].fene.bm;
        kb   = forceparams[type].fene.kb;

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx);     /*   3          */
        dr2  = iprod(dx, dx);                           /*   5          */

        if (dr2 == 0.0)
        {
            continue;
        }

        bm2 = bm*bm;

        if (dr2 >= bm2)
        {
            gmx_fatal(FARGS,
                      "r^2 (%f) >= bm^2 (%f) in FENE bond between atoms %d and %d",
                      dr2, bm2,
                      glatnr(global_atom_index, ai),
                      glatnr(global_atom_index, aj));
        }

        omdr2obm2  = one - dr2/bm2;

        vbond      = -half*kb*bm2*std::log(omdr2obm2);
        fbond      = -kb/omdr2obm2;

        vtot      += vbond;   /* 35 */

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)                    /*  15          */
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }                                         /*  58 TOTAL    */
    return vtot;
}

static real harmonic(real kA, real kB, real xA, real xB, real x, real lambda,
                     real *V, real *F)
{
    const real half = 0.5;
    real       L1, kk, x0, dx, dx2;
    real       v, f, dvdlambda;

    L1    = 1.0-lambda;
    kk    = L1*kA+lambda*kB;
    x0    = L1*xA+lambda*xB;

    dx    = x-x0;
    dx2   = dx*dx;

    f          = -kk*dx;
    v          = half*kk*dx2;
    dvdlambda  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;

    *F    = f;
    *V    = v;

    return dvdlambda;

    /* That was 19 flops */
}


real bonds(int nbonds,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec4 f[], rvec fshift[],
           const t_pbc *pbc, const t_graph *g,
           real lambda, real *dvdlambda,
           const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
           int gmx_unused *global_atom_index)
{
    int  i, m, ki, ai, aj, type;
    real dr, dr2, fbond, vbond, fij, vtot;
    rvec dx;
    ivec dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx); /*   3      */
        dr2  = iprod(dx, dx);                       /*   5		*/
        dr   = std::sqrt(dr2);                      /*  10		*/

        *dvdlambda += harmonic(forceparams[type].harmonic.krA,
                               forceparams[type].harmonic.krB,
                               forceparams[type].harmonic.rA,
                               forceparams[type].harmonic.rB,
                               dr, lambda, &vbond, &fbond); /*  19  */

        if (dr2 == 0.0)
        {
            continue;
        }


        vtot  += vbond;             /* 1*/
        fbond *= gmx::invsqrt(dr2); /*   6		*/
        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)     /*  15		*/
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }               /* 59 TOTAL	*/
    return vtot;
}

real restraint_bonds(int nbonds,
                     const t_iatom forceatoms[], const t_iparams forceparams[],
                     const rvec x[], rvec4 f[], rvec fshift[],
                     const t_pbc *pbc, const t_graph *g,
                     real lambda, real *dvdlambda,
                     const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                     int gmx_unused *global_atom_index)
{
    int  i, m, ki, ai, aj, type;
    real dr, dr2, fbond, vbond, fij, vtot;
    real L1;
    real low, dlow, up1, dup1, up2, dup2, k, dk;
    real drh, drh2;
    rvec dx;
    ivec dt;

    L1   = 1.0 - lambda;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx); /*   3      */
        dr2  = iprod(dx, dx);                       /*   5		*/
        dr   = dr2*gmx::invsqrt(dr2);               /*  10		*/

        low  = L1*forceparams[type].restraint.lowA + lambda*forceparams[type].restraint.lowB;
        dlow =   -forceparams[type].restraint.lowA +        forceparams[type].restraint.lowB;
        up1  = L1*forceparams[type].restraint.up1A + lambda*forceparams[type].restraint.up1B;
        dup1 =   -forceparams[type].restraint.up1A +        forceparams[type].restraint.up1B;
        up2  = L1*forceparams[type].restraint.up2A + lambda*forceparams[type].restraint.up2B;
        dup2 =   -forceparams[type].restraint.up2A +        forceparams[type].restraint.up2B;
        k    = L1*forceparams[type].restraint.kA   + lambda*forceparams[type].restraint.kB;
        dk   =   -forceparams[type].restraint.kA   +        forceparams[type].restraint.kB;
        /* 24 */

        if (dr < low)
        {
            drh         = dr - low;
            drh2        = drh*drh;
            vbond       = 0.5*k*drh2;
            fbond       = -k*drh;
            *dvdlambda += 0.5*dk*drh2 - k*dlow*drh;
        } /* 11 */
        else if (dr <= up1)
        {
            vbond = 0;
            fbond = 0;
        }
        else if (dr <= up2)
        {
            drh         = dr - up1;
            drh2        = drh*drh;
            vbond       = 0.5*k*drh2;
            fbond       = -k*drh;
            *dvdlambda += 0.5*dk*drh2 - k*dup1*drh;
        } /* 11	*/
        else
        {
            drh         = dr - up2;
            vbond       = k*(up2 - up1)*(0.5*(up2 - up1) + drh);
            fbond       = -k*(up2 - up1);
            *dvdlambda += dk*(up2 - up1)*(0.5*(up2 - up1) + drh)
                + k*(dup2 - dup1)*(up2 - up1 + drh)
                - k*(up2 - up1)*dup2;
        }

        if (dr2 == 0.0)
        {
            continue;
        }

        vtot  += vbond;             /* 1*/
        fbond *= gmx::invsqrt(dr2); /*   6		*/
        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)             /*  15		*/
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }                   /* 59 TOTAL	*/

    return vtot;
}

real polarize(int nbonds,
              const t_iatom forceatoms[], const t_iparams forceparams[],
              const rvec x[], rvec4 f[], rvec fshift[],
              const t_pbc *pbc, const t_graph *g,
              real lambda, real *dvdlambda,
              const t_mdatoms *md, t_fcdata gmx_unused *fcd,
              int gmx_unused *global_atom_index)
{
    int  i, m, ki, ai, aj, type;
    real dr, dr2, fbond, vbond, fij, vtot, ksh;
    rvec dx;
    ivec dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ksh  = gmx::square(md->chargeA[aj])*ONE_4PI_EPS0/forceparams[type].polarize.alpha;

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx);                         /*   3      */
        dr2  = iprod(dx, dx);                                               /*   5		*/
        dr   = std::sqrt(dr2);                                              /*  10		*/

        *dvdlambda += harmonic(ksh, ksh, 0, 0, dr, lambda, &vbond, &fbond); /*  19  */

        if (dr2 == 0.0)
        {
            continue;
        }

        vtot  += vbond;             /* 1*/
        fbond *= gmx::invsqrt(dr2); /*   6		*/

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)     /*  15		*/
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }               /* 59 TOTAL	*/
    return vtot;
}

real anharm_polarize(int nbonds,
                     const t_iatom forceatoms[], const t_iparams forceparams[],
                     const rvec x[], rvec4 f[], rvec fshift[],
                     const t_pbc *pbc, const t_graph *g,
                     real lambda, real *dvdlambda,
                     const t_mdatoms *md, t_fcdata gmx_unused *fcd,
                     int gmx_unused *global_atom_index)
{
    int  i, m, ki, ai, aj, type;
    real dr, dr2, fbond, vbond, fij, vtot, ksh, khyp, drcut, ddr, ddr3;
    rvec dx;
    ivec dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type  = forceatoms[i++];
        ai    = forceatoms[i++];
        aj    = forceatoms[i++];
        ksh   = gmx::square(md->chargeA[aj])*ONE_4PI_EPS0/forceparams[type].anharm_polarize.alpha; /* 7*/
        khyp  = forceparams[type].anharm_polarize.khyp;
        drcut = forceparams[type].anharm_polarize.drcut;

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx);                         /*   3      */
        dr2  = iprod(dx, dx);                                               /*   5		*/
        dr   = dr2*gmx::invsqrt(dr2);                                       /*  10		*/

        *dvdlambda += harmonic(ksh, ksh, 0, 0, dr, lambda, &vbond, &fbond); /*  19  */

        if (dr2 == 0.0)
        {
            continue;
        }

        if (dr > drcut)
        {
            ddr    = dr-drcut;
            ddr3   = ddr*ddr*ddr;
            vbond += khyp*ddr*ddr3;
            fbond -= 4*khyp*ddr3;
        }
        fbond *= gmx::invsqrt(dr2); /*   6		*/
        vtot  += vbond;             /* 1*/

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)     /*  15		*/
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }               /* 72 TOTAL	*/
    return vtot;
}

real water_pol(int nbonds,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec4 f[], rvec gmx_unused fshift[],
               const t_pbc gmx_unused *pbc, const t_graph gmx_unused *g,
               real gmx_unused lambda, real gmx_unused *dvdlambda,
               const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
               int gmx_unused *global_atom_index)
{
    /* This routine implements anisotropic polarizibility for water, through
     * a shell connected to a dummy with spring constant that differ in the
     * three spatial dimensions in the molecular frame.
     */
    int  i, m, aO, aH1, aH2, aD, aS, type, type0, ki;
    ivec dt;
    rvec dOH1, dOH2, dHH, dOD, dDS, nW, kk, dx, kdx, proj;
    real vtot, fij, r_HH, r_OD, r_nW, tx, ty, tz, qS;

    vtot = 0.0;
    if (nbonds > 0)
    {
        type0  = forceatoms[0];
        aS     = forceatoms[5];
        qS     = md->chargeA[aS];
        kk[XX] = gmx::square(qS)*ONE_4PI_EPS0/forceparams[type0].wpol.al_x;
        kk[YY] = gmx::square(qS)*ONE_4PI_EPS0/forceparams[type0].wpol.al_y;
        kk[ZZ] = gmx::square(qS)*ONE_4PI_EPS0/forceparams[type0].wpol.al_z;
        r_HH   = 1.0/forceparams[type0].wpol.rHH;
        for (i = 0; (i < nbonds); i += 6)
        {
            type = forceatoms[i];
            if (type != type0)
            {
                gmx_fatal(FARGS, "Sorry, type = %d, type0 = %d, file = %s, line = %d",
                          type, type0, __FILE__, __LINE__);
            }
            aO   = forceatoms[i+1];
            aH1  = forceatoms[i+2];
            aH2  = forceatoms[i+3];
            aD   = forceatoms[i+4];
            aS   = forceatoms[i+5];

            /* Compute vectors describing the water frame */
            pbc_rvec_sub(pbc, x[aH1], x[aO], dOH1);
            pbc_rvec_sub(pbc, x[aH2], x[aO], dOH2);
            pbc_rvec_sub(pbc, x[aH2], x[aH1], dHH);
            pbc_rvec_sub(pbc, x[aD], x[aO], dOD);
            ki = pbc_rvec_sub(pbc, x[aS], x[aD], dDS);
            cprod(dOH1, dOH2, nW);

            /* Compute inverse length of normal vector
             * (this one could be precomputed, but I'm too lazy now)
             */
            r_nW = gmx::invsqrt(iprod(nW, nW));
            /* This is for precision, but does not make a big difference,
             * it can go later.
             */
            r_OD = gmx::invsqrt(iprod(dOD, dOD));

            /* Normalize the vectors in the water frame */
            svmul(r_nW, nW, nW);
            svmul(r_HH, dHH, dHH);
            svmul(r_OD, dOD, dOD);

            /* Compute displacement of shell along components of the vector */
            dx[ZZ] = iprod(dDS, dOD);
            /* Compute projection on the XY plane: dDS - dx[ZZ]*dOD */
            for (m = 0; (m < DIM); m++)
            {
                proj[m] = dDS[m]-dx[ZZ]*dOD[m];
            }

            /*dx[XX] = iprod(dDS,nW);
               dx[YY] = iprod(dDS,dHH);*/
            dx[XX] = iprod(proj, nW);
            for (m = 0; (m < DIM); m++)
            {
                proj[m] -= dx[XX]*nW[m];
            }
            dx[YY] = iprod(proj, dHH);
            /* Now compute the forces and energy */
            kdx[XX] = kk[XX]*dx[XX];
            kdx[YY] = kk[YY]*dx[YY];
            kdx[ZZ] = kk[ZZ]*dx[ZZ];
            vtot   += iprod(dx, kdx);

            if (g)
            {
                ivec_sub(SHIFT_IVEC(g, aS), SHIFT_IVEC(g, aD), dt);
                ki = IVEC2IS(dt);
            }

            for (m = 0; (m < DIM); m++)
            {
                /* This is a tensor operation but written out for speed */
                tx                  =  nW[m]*kdx[XX];
                ty                  = dHH[m]*kdx[YY];
                tz                  = dOD[m]*kdx[ZZ];
                fij                 = -tx-ty-tz;
                f[aS][m]           += fij;
                f[aD][m]           -= fij;
                fshift[ki][m]      += fij;
                fshift[CENTRAL][m] -= fij;
            }
        }
    }
    return 0.5*vtot;
}

static real do_1_thole(const rvec xi, const rvec xj, rvec fi, rvec fj,
                       const t_pbc *pbc, real qq,
                       rvec fshift[], real afac)
{
    rvec r12;
    real r12sq, r12_1, r12bar, v0, v1, fscal, ebar, fff;
    int  m, t;

    t      = pbc_rvec_sub(pbc, xi, xj, r12);                      /*  3 */

    r12sq  = iprod(r12, r12);                                     /*  5 */
    r12_1  = gmx::invsqrt(r12sq);                                 /*  5 */
    r12bar = afac/r12_1;                                          /*  5 */
    v0     = qq*ONE_4PI_EPS0*r12_1;                               /*  2 */
    ebar   = std::exp(-r12bar);                                   /*  5 */
    v1     = (1-(1+0.5*r12bar)*ebar);                             /*  4 */
    fscal  = ((v0*r12_1)*v1 - v0*0.5*afac*ebar*(r12bar+1))*r12_1; /* 9 */

    for (m = 0; (m < DIM); m++)
    {
        fff                 = fscal*r12[m];
        fi[m]              += fff;
        fj[m]              -= fff;
        fshift[t][m]       += fff;
        fshift[CENTRAL][m] -= fff;
    }             /* 15 */

    return v0*v1; /* 1 */
    /* 54 */
}

real thole_pol(int nbonds,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec4 f[], rvec fshift[],
               const t_pbc *pbc, const t_graph gmx_unused *g,
               real gmx_unused lambda, real gmx_unused *dvdlambda,
               const t_mdatoms *md, t_fcdata gmx_unused *fcd,
               int gmx_unused *global_atom_index)
{
    /* Interaction between two pairs of particles with opposite charge */
    int        i, type, a1, da1, a2, da2;
    real       q1, q2, qq, a, al1, al2, afac;
    real       V             = 0;

    for (i = 0; (i < nbonds); )
    {
        type  = forceatoms[i++];
        a1    = forceatoms[i++];
        da1   = forceatoms[i++];
        a2    = forceatoms[i++];
        da2   = forceatoms[i++];
        q1    = md->chargeA[da1];
        q2    = md->chargeA[da2];
        a     = forceparams[type].thole.a;
        al1   = forceparams[type].thole.alpha1;
        al2   = forceparams[type].thole.alpha2;
        qq    = q1*q2;
        afac  = a*gmx::invsixthroot(al1*al2);
        V    += do_1_thole(x[a1], x[a2], f[a1], f[a2], pbc, qq, fshift, afac);
        V    += do_1_thole(x[da1], x[a2], f[da1], f[a2], pbc, -qq, fshift, afac);
        V    += do_1_thole(x[a1], x[da2], f[a1], f[da2], pbc, -qq, fshift, afac);
        V    += do_1_thole(x[da1], x[da2], f[da1], f[da2], pbc, qq, fshift, afac);
    }
    /* 290 flops */
    return V;
}

real bond_angle(const rvec xi, const rvec xj, const rvec xk, const t_pbc *pbc,
                rvec r_ij, rvec r_kj, real *costh,
                int *t1, int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
    /* 41 FLOPS */
    real th;

    *t1 = pbc_rvec_sub(pbc, xi, xj, r_ij); /*  3		*/
    *t2 = pbc_rvec_sub(pbc, xk, xj, r_kj); /*  3		*/

    *costh = cos_angle(r_ij, r_kj);        /* 25		*/
    th     = std::acos(*costh);            /* 10		*/
    /* 41 TOTAL	*/
    return th;
}

real angles(int nbonds,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec4 f[], rvec fshift[],
            const t_pbc *pbc, const t_graph *g,
            real lambda, real *dvdlambda,
            const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
            int gmx_unused *global_atom_index)
{
    int  i, ai, aj, ak, t1, t2, type;
    rvec r_ij, r_kj;
    real cos_theta, cos_theta2, theta, dVdt, va, vtot;
    ivec jt, dt_ij, dt_kj;

    vtot = 0.0;
    for (i = 0; i < nbonds; )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];

        theta  = bond_angle(x[ai], x[aj], x[ak], pbc,
                            r_ij, r_kj, &cos_theta, &t1, &t2);  /*  41		*/

        *dvdlambda += harmonic(forceparams[type].harmonic.krA,
                               forceparams[type].harmonic.krB,
                               forceparams[type].harmonic.rA*DEG2RAD,
                               forceparams[type].harmonic.rB*DEG2RAD,
                               theta, lambda, &va, &dVdt);  /*  21  */
        vtot += va;

        cos_theta2 = gmx::square(cos_theta);
        if (cos_theta2 < 1)
        {
            int  m;
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            real nrkj_1, nrij_1;
            rvec f_i, f_j, f_k;

            st    = dVdt*gmx::invsqrt(1 - cos_theta2); /*  12		*/
            sth   = st*cos_theta;                      /*   1		*/
            nrij2 = iprod(r_ij, r_ij);                 /*   5		*/
            nrkj2 = iprod(r_kj, r_kj);                 /*   5		*/

            nrij_1 = gmx::invsqrt(nrij2);              /*  10		*/
            nrkj_1 = gmx::invsqrt(nrkj2);              /*  10		*/

            cik = st*nrij_1*nrkj_1;                    /*   2		*/
            cii = sth*nrij_1*nrij_1;                   /*   2		*/
            ckk = sth*nrkj_1*nrkj_1;                   /*   2		*/

            for (m = 0; m < DIM; m++)
            {           /*  39		*/
                f_i[m]    = -(cik*r_kj[m] - cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m] - ckk*r_kj[m]);
                f_j[m]    = -f_i[m] - f_k[m];
                f[ai][m] += f_i[m];
                f[aj][m] += f_j[m];
                f[ak][m] += f_k[m];
            }
            if (g != nullptr)
            {
                copy_ivec(SHIFT_IVEC(g, aj), jt);

                ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
                ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
                t1 = IVEC2IS(dt_ij);
                t2 = IVEC2IS(dt_kj);
            }
            rvec_inc(fshift[t1], f_i);
            rvec_inc(fshift[CENTRAL], f_j);
            rvec_inc(fshift[t2], f_k);
        }                                           /* 161 TOTAL	*/
    }

    return vtot;
}

#if GMX_SIMD_HAVE_REAL

/* As angles, but using SIMD to calculate many angles at once.
 * This routines does not calculate energies and shift forces.
 */
void
angles_noener_simd(int nbonds,
                   const t_iatom forceatoms[], const t_iparams forceparams[],
                   const rvec x[], rvec4 f[],
                   const t_pbc *pbc, const t_graph gmx_unused *g,
                   real gmx_unused lambda,
                   const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                   int gmx_unused *global_atom_index)
{
    const int            nfa1 = 4;
    int                  i, iu, s;
    int                  type;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    ai[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    aj[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    ak[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) real            coeff[2*GMX_SIMD_REAL_WIDTH];
    SimdReal             deg2rad_S(DEG2RAD);
    SimdReal             xi_S, yi_S, zi_S;
    SimdReal             xj_S, yj_S, zj_S;
    SimdReal             xk_S, yk_S, zk_S;
    SimdReal             k_S, theta0_S;
    SimdReal             rijx_S, rijy_S, rijz_S;
    SimdReal             rkjx_S, rkjy_S, rkjz_S;
    SimdReal             one_S(1.0);
    SimdReal             min_one_plus_eps_S(-1.0 + 2.0*GMX_REAL_EPS); // Smallest number > -1

    SimdReal             rij_rkj_S;
    SimdReal             nrij2_S, nrij_1_S;
    SimdReal             nrkj2_S, nrkj_1_S;
    SimdReal             cos_S, invsin_S;
    SimdReal             theta_S;
    SimdReal             st_S, sth_S;
    SimdReal             cik_S, cii_S, ckk_S;
    SimdReal             f_ix_S, f_iy_S, f_iz_S;
    SimdReal             f_kx_S, f_ky_S, f_kz_S;
    alignas(GMX_SIMD_ALIGNMENT) real    pbc_simd[9*GMX_SIMD_REAL_WIDTH];

    set_pbc_simd(pbc, pbc_simd);

    /* nbonds is the number of angles times nfa1, here we step GMX_SIMD_REAL_WIDTH angles */
    for (i = 0; (i < nbonds); i += GMX_SIMD_REAL_WIDTH*nfa1)
    {
        /* Collect atoms for GMX_SIMD_REAL_WIDTH angles.
         * iu indexes into forceatoms, we should not let iu go beyond nbonds.
         */
        iu = i;
        for (s = 0; s < GMX_SIMD_REAL_WIDTH; s++)
        {
            type  = forceatoms[iu];
            ai[s] = forceatoms[iu+1];
            aj[s] = forceatoms[iu+2];
            ak[s] = forceatoms[iu+3];

            /* At the end fill the arrays with the last atoms and 0 params */
            if (i + s*nfa1 < nbonds)
            {
                coeff[s]                     = forceparams[type].harmonic.krA;
                coeff[GMX_SIMD_REAL_WIDTH+s] = forceparams[type].harmonic.rA;

                if (iu + nfa1 < nbonds)
                {
                    iu += nfa1;
                }
            }
            else
            {
                coeff[s]                     = 0;
                coeff[GMX_SIMD_REAL_WIDTH+s] = 0;
            }
        }

        /* Store the non PBC corrected distances packed and aligned */
        gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), ai, &xi_S, &yi_S, &zi_S);
        gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), aj, &xj_S, &yj_S, &zj_S);
        gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), ak, &xk_S, &yk_S, &zk_S);
        rijx_S = xi_S - xj_S;
        rijy_S = yi_S - yj_S;
        rijz_S = zi_S - zj_S;
        rkjx_S = xk_S - xj_S;
        rkjy_S = yk_S - yj_S;
        rkjz_S = zk_S - zj_S;

        k_S       = load<SimdReal>(coeff);
        theta0_S  = load<SimdReal>(coeff+GMX_SIMD_REAL_WIDTH) * deg2rad_S;

        pbc_correct_dx_simd(&rijx_S, &rijy_S, &rijz_S, pbc_simd);
        pbc_correct_dx_simd(&rkjx_S, &rkjy_S, &rkjz_S, pbc_simd);

        rij_rkj_S = iprod(rijx_S, rijy_S, rijz_S,
                          rkjx_S, rkjy_S, rkjz_S);

        nrij2_S   = norm2(rijx_S, rijy_S, rijz_S);
        nrkj2_S   = norm2(rkjx_S, rkjy_S, rkjz_S);

        nrij_1_S  = invsqrt(nrij2_S);
        nrkj_1_S  = invsqrt(nrkj2_S);

        cos_S     = rij_rkj_S * nrij_1_S * nrkj_1_S;

        /* To allow for 180 degrees, we take the max of cos and -1 + 1bit,
         * so we can safely get the 1/sin from 1/sqrt(1 - cos^2).
         * This also ensures that rounding errors would cause the argument
         * of simdAcos to be < -1.
         * Note that we do not take precautions for cos(0)=1, so the outer
         * atoms in an angle should not be on top of each other.
         */
        cos_S     = max(cos_S, min_one_plus_eps_S);

        theta_S   = acos(cos_S);

        invsin_S  = invsqrt( one_S - cos_S * cos_S );

        st_S      = k_S * (theta0_S - theta_S) * invsin_S;
        sth_S     = st_S * cos_S;

        cik_S     = st_S  * nrij_1_S * nrkj_1_S;
        cii_S     = sth_S * nrij_1_S * nrij_1_S;
        ckk_S     = sth_S * nrkj_1_S * nrkj_1_S;

        f_ix_S    = cii_S * rijx_S;
        f_ix_S    = fnma(cik_S, rkjx_S, f_ix_S);
        f_iy_S    = cii_S * rijy_S;
        f_iy_S    = fnma(cik_S, rkjy_S, f_iy_S);
        f_iz_S    = cii_S * rijz_S;
        f_iz_S    = fnma(cik_S, rkjz_S, f_iz_S);
        f_kx_S    = ckk_S * rkjx_S;
        f_kx_S    = fnma(cik_S, rijx_S, f_kx_S);
        f_ky_S    = ckk_S * rkjy_S;
        f_ky_S    = fnma(cik_S, rijy_S, f_ky_S);
        f_kz_S    = ckk_S * rkjz_S;
        f_kz_S    = fnma(cik_S, rijz_S, f_kz_S);

        transposeScatterIncrU<4>(reinterpret_cast<real *>(f), ai, f_ix_S, f_iy_S, f_iz_S);
        transposeScatterDecrU<4>(reinterpret_cast<real *>(f), aj, f_ix_S + f_kx_S, f_iy_S + f_ky_S, f_iz_S + f_kz_S);
        transposeScatterIncrU<4>(reinterpret_cast<real *>(f), ak, f_kx_S, f_ky_S, f_kz_S);
    }
}

#endif // GMX_SIMD_HAVE_REAL

real linear_angles(int nbonds,
                   const t_iatom forceatoms[], const t_iparams forceparams[],
                   const rvec x[], rvec4 f[], rvec fshift[],
                   const t_pbc *pbc, const t_graph *g,
                   real lambda, real *dvdlambda,
                   const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                   int gmx_unused *global_atom_index)
{
    int  i, m, ai, aj, ak, t1, t2, type;
    rvec f_i, f_j, f_k;
    real L1, kA, kB, aA, aB, dr, dr2, va, vtot, a, b, klin;
    ivec jt, dt_ij, dt_kj;
    rvec r_ij, r_kj, r_ik, dx;

    L1   = 1-lambda;
    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];

        kA   = forceparams[type].linangle.klinA;
        kB   = forceparams[type].linangle.klinB;
        klin = L1*kA + lambda*kB;

        aA   = forceparams[type].linangle.aA;
        aB   = forceparams[type].linangle.aB;
        a    = L1*aA+lambda*aB;
        b    = 1-a;

        t1 = pbc_rvec_sub(pbc, x[ai], x[aj], r_ij);
        t2 = pbc_rvec_sub(pbc, x[ak], x[aj], r_kj);
        rvec_sub(r_ij, r_kj, r_ik);

        dr2 = 0;
        for (m = 0; (m < DIM); m++)
        {
            dr        = -a * r_ij[m] - b * r_kj[m];
            dr2      += dr*dr;
            dx[m]     = dr;
            f_i[m]    = a*klin*dr;
            f_k[m]    = b*klin*dr;
            f_j[m]    = -(f_i[m]+f_k[m]);
            f[ai][m] += f_i[m];
            f[aj][m] += f_j[m];
            f[ak][m] += f_k[m];
        }
        va          = 0.5*klin*dr2;
        *dvdlambda += 0.5*(kB-kA)*dr2 + klin*(aB-aA)*iprod(dx, r_ik);

        vtot += va;

        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);

            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
        }
        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);
    }                                         /* 57 TOTAL	*/
    return vtot;
}

real urey_bradley(int nbonds,
                  const t_iatom forceatoms[], const t_iparams forceparams[],
                  const rvec x[], rvec4 f[], rvec fshift[],
                  const t_pbc *pbc, const t_graph *g,
                  real lambda, real *dvdlambda,
                  const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                  int gmx_unused *global_atom_index)
{
    int  i, m, ai, aj, ak, t1, t2, type, ki;
    rvec r_ij, r_kj, r_ik;
    real cos_theta, cos_theta2, theta;
    real dVdt, va, vtot, dr, dr2, vbond, fbond, fik;
    real kthA, th0A, kUBA, r13A, kthB, th0B, kUBB, r13B;
    ivec jt, dt_ij, dt_kj, dt_ik;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type  = forceatoms[i++];
        ai    = forceatoms[i++];
        aj    = forceatoms[i++];
        ak    = forceatoms[i++];
        th0A  = forceparams[type].u_b.thetaA*DEG2RAD;
        kthA  = forceparams[type].u_b.kthetaA;
        r13A  = forceparams[type].u_b.r13A;
        kUBA  = forceparams[type].u_b.kUBA;
        th0B  = forceparams[type].u_b.thetaB*DEG2RAD;
        kthB  = forceparams[type].u_b.kthetaB;
        r13B  = forceparams[type].u_b.r13B;
        kUBB  = forceparams[type].u_b.kUBB;

        theta  = bond_angle(x[ai], x[aj], x[ak], pbc,
                            r_ij, r_kj, &cos_theta, &t1, &t2);                     /*  41		*/

        *dvdlambda += harmonic(kthA, kthB, th0A, th0B, theta, lambda, &va, &dVdt); /*  21  */
        vtot       += va;

        ki   = pbc_rvec_sub(pbc, x[ai], x[ak], r_ik);                               /*   3      */
        dr2  = iprod(r_ik, r_ik);                                                   /*   5		*/
        dr   = dr2*gmx::invsqrt(dr2);                                               /*  10		*/

        *dvdlambda += harmonic(kUBA, kUBB, r13A, r13B, dr, lambda, &vbond, &fbond); /*  19  */

        cos_theta2 = gmx::square(cos_theta);                                        /*   1		*/
        if (cos_theta2 < 1)
        {
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            rvec f_i, f_j, f_k;

            st    = dVdt*gmx::invsqrt(1 - cos_theta2); /*  12		*/
            sth   = st*cos_theta;                      /*   1		*/
            nrkj2 = iprod(r_kj, r_kj);                 /*   5		*/
            nrij2 = iprod(r_ij, r_ij);

            cik = st*gmx::invsqrt(nrkj2*nrij2); /*  12		*/
            cii = sth/nrij2;                    /*  10		*/
            ckk = sth/nrkj2;                    /*  10		*/

            for (m = 0; (m < DIM); m++)         /*  39		*/
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                f[ai][m] += f_i[m];
                f[aj][m] += f_j[m];
                f[ak][m] += f_k[m];
            }
            if (g)
            {
                copy_ivec(SHIFT_IVEC(g, aj), jt);

                ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
                ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
                t1 = IVEC2IS(dt_ij);
                t2 = IVEC2IS(dt_kj);
            }
            rvec_inc(fshift[t1], f_i);
            rvec_inc(fshift[CENTRAL], f_j);
            rvec_inc(fshift[t2], f_k);
        }                                       /* 161 TOTAL	*/
        /* Time for the bond calculations */
        if (dr2 == 0.0)
        {
            continue;
        }

        vtot  += vbond;             /* 1*/
        fbond *= gmx::invsqrt(dr2); /*   6		*/

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, ak), dt_ik);
            ki = IVEC2IS(dt_ik);
        }
        for (m = 0; (m < DIM); m++)     /*  15		*/
        {
            fik                 = fbond*r_ik[m];
            f[ai][m]           += fik;
            f[ak][m]           -= fik;
            fshift[ki][m]      += fik;
            fshift[CENTRAL][m] -= fik;
        }
    }
    return vtot;
}

#if GMX_SIMD_HAVE_REAL

/* As urey_bradley, but using SIMD to calculate many potentials at once.
 * This routines does not calculate energies and shift forces.
 */
void urey_bradley_noener_simd(int nbonds,
                              const t_iatom forceatoms[], const t_iparams forceparams[],
                              const rvec x[], rvec4 f[],
                              const t_pbc *pbc, const t_graph gmx_unused *g,
                              real gmx_unused lambda,
                              const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                              int gmx_unused *global_atom_index)
{
    constexpr int            nfa1 = 4;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    ai[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    aj[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    ak[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) real            coeff[4*GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) real            pbc_simd[9*GMX_SIMD_REAL_WIDTH];

    set_pbc_simd(pbc, pbc_simd);

    /* nbonds is the number of angles times nfa1, here we step GMX_SIMD_REAL_WIDTH angles */
    for (int i = 0; i < nbonds; i += GMX_SIMD_REAL_WIDTH*nfa1)
    {
        /* Collect atoms for GMX_SIMD_REAL_WIDTH angles.
         * iu indexes into forceatoms, we should not let iu go beyond nbonds.
         */
        int iu = i;
        for (int s = 0; s < GMX_SIMD_REAL_WIDTH; s++)
        {
            const int type                     = forceatoms[iu];
            ai[s] = forceatoms[iu+1];
            aj[s] = forceatoms[iu+2];
            ak[s] = forceatoms[iu+3];

            /* At the end fill the arrays with the last atoms and 0 params */
            if (i + s*nfa1 < nbonds)
            {
                coeff[s]                       = forceparams[type].u_b.kthetaA;
                coeff[GMX_SIMD_REAL_WIDTH+s]   = forceparams[type].u_b.thetaA;
                coeff[GMX_SIMD_REAL_WIDTH*2+s] = forceparams[type].u_b.kUBA;
                coeff[GMX_SIMD_REAL_WIDTH*3+s] = forceparams[type].u_b.r13A;

                if (iu + nfa1 < nbonds)
                {
                    iu += nfa1;
                }
            }
            else
            {
                coeff[s]                       = 0;
                coeff[GMX_SIMD_REAL_WIDTH+s]   = 0;
                coeff[GMX_SIMD_REAL_WIDTH*2+s] = 0;
                coeff[GMX_SIMD_REAL_WIDTH*3+s] = 0;
            }
        }

        SimdReal xi_S, yi_S, zi_S;
        SimdReal xj_S, yj_S, zj_S;
        SimdReal xk_S, yk_S, zk_S;

        /* Store the non PBC corrected distances packed and aligned */
        gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), ai, &xi_S, &yi_S, &zi_S);
        gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), aj, &xj_S, &yj_S, &zj_S);
        gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), ak, &xk_S, &yk_S, &zk_S);
        SimdReal       rijx_S = xi_S - xj_S;
        SimdReal       rijy_S = yi_S - yj_S;
        SimdReal       rijz_S = zi_S - zj_S;
        SimdReal       rkjx_S = xk_S - xj_S;
        SimdReal       rkjy_S = yk_S - yj_S;
        SimdReal       rkjz_S = zk_S - zj_S;
        SimdReal       rikx_S = xi_S - xk_S;
        SimdReal       riky_S = yi_S - yk_S;
        SimdReal       rikz_S = zi_S - zk_S;

        const SimdReal ktheta_S = load<SimdReal>(coeff);
        const SimdReal theta0_S = load<SimdReal>(coeff+GMX_SIMD_REAL_WIDTH) * DEG2RAD;
        const SimdReal kUB_S    = load<SimdReal>(coeff+2*GMX_SIMD_REAL_WIDTH);
        const SimdReal r13_S    = load<SimdReal>(coeff+3*GMX_SIMD_REAL_WIDTH);

        pbc_correct_dx_simd(&rijx_S, &rijy_S, &rijz_S, pbc_simd);
        pbc_correct_dx_simd(&rkjx_S, &rkjy_S, &rkjz_S, pbc_simd);
        pbc_correct_dx_simd(&rikx_S, &riky_S, &rikz_S, pbc_simd);

        const SimdReal rij_rkj_S = iprod(rijx_S, rijy_S, rijz_S,
                                         rkjx_S, rkjy_S, rkjz_S);

        const SimdReal dr2_S     = iprod(rikx_S, riky_S, rikz_S,
                                         rikx_S, riky_S, rikz_S);

        const SimdReal nrij2_S   = norm2(rijx_S, rijy_S, rijz_S);
        const SimdReal nrkj2_S   = norm2(rkjx_S, rkjy_S, rkjz_S);

        const SimdReal nrij_1_S  = invsqrt(nrij2_S);
        const SimdReal nrkj_1_S  = invsqrt(nrkj2_S);
        const SimdReal invdr2_S  = invsqrt(dr2_S);
        const SimdReal dr_S      = dr2_S*invdr2_S;

        constexpr real min_one_plus_eps = -1.0 + 2.0*GMX_REAL_EPS; // Smallest number > -1

        /* To allow for 180 degrees, we take the max of cos and -1 + 1bit,
         * so we can safely get the 1/sin from 1/sqrt(1 - cos^2).
         * This also ensures that rounding errors would cause the argument
         * of simdAcos to be < -1.
         * Note that we do not take precautions for cos(0)=1, so the bonds
         * in an angle should not align at an angle of 0 degrees.
         */
        const SimdReal cos_S     = max(rij_rkj_S * nrij_1_S * nrkj_1_S, min_one_plus_eps);

        const SimdReal theta_S   = acos(cos_S);
        const SimdReal invsin_S  = invsqrt( 1.0 - cos_S * cos_S );
        const SimdReal st_S      = ktheta_S * (theta0_S - theta_S) * invsin_S;
        const SimdReal sth_S     = st_S * cos_S;

        const SimdReal cik_S     = st_S  * nrij_1_S * nrkj_1_S;
        const SimdReal cii_S     = sth_S * nrij_1_S * nrij_1_S;
        const SimdReal ckk_S     = sth_S * nrkj_1_S * nrkj_1_S;

        const SimdReal sUB_S     = kUB_S * (r13_S - dr_S) * invdr2_S;

        const SimdReal f_ikx_S   = sUB_S * rikx_S;
        const SimdReal f_iky_S   = sUB_S * riky_S;
        const SimdReal f_ikz_S   = sUB_S * rikz_S;

        const SimdReal f_ix_S    = fnma(cik_S, rkjx_S, cii_S * rijx_S) + f_ikx_S;
        const SimdReal f_iy_S    = fnma(cik_S, rkjy_S, cii_S * rijy_S) + f_iky_S;
        const SimdReal f_iz_S    = fnma(cik_S, rkjz_S, cii_S * rijz_S) + f_ikz_S;
        const SimdReal f_kx_S    = fnma(cik_S, rijx_S, ckk_S * rkjx_S) - f_ikx_S;
        const SimdReal f_ky_S    = fnma(cik_S, rijy_S, ckk_S * rkjy_S) - f_iky_S;
        const SimdReal f_kz_S    = fnma(cik_S, rijz_S, ckk_S * rkjz_S) - f_ikz_S;

        transposeScatterIncrU<4>(reinterpret_cast<real *>(f), ai, f_ix_S, f_iy_S, f_iz_S);
        transposeScatterDecrU<4>(reinterpret_cast<real *>(f), aj, f_ix_S + f_kx_S, f_iy_S + f_ky_S, f_iz_S + f_kz_S);
        transposeScatterIncrU<4>(reinterpret_cast<real *>(f), ak, f_kx_S, f_ky_S, f_kz_S);
    }
}

#endif // GMX_SIMD_HAVE_REAL

real quartic_angles(int nbonds,
                    const t_iatom forceatoms[], const t_iparams forceparams[],
                    const rvec x[], rvec4 f[], rvec fshift[],
                    const t_pbc *pbc, const t_graph *g,
                    real gmx_unused lambda, real gmx_unused *dvdlambda,
                    const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                    int gmx_unused *global_atom_index)
{
    int  i, j, ai, aj, ak, t1, t2, type;
    rvec r_ij, r_kj;
    real cos_theta, cos_theta2, theta, dt, dVdt, va, dtp, c, vtot;
    ivec jt, dt_ij, dt_kj;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];

        theta  = bond_angle(x[ai], x[aj], x[ak], pbc,
                            r_ij, r_kj, &cos_theta, &t1, &t2); /*  41		*/

        dt = theta - forceparams[type].qangle.theta*DEG2RAD;   /* 2          */

        dVdt = 0;
        va   = forceparams[type].qangle.c[0];
        dtp  = 1.0;
        for (j = 1; j <= 4; j++)
        {
            c     = forceparams[type].qangle.c[j];
            dVdt -= j*c*dtp;
            dtp  *= dt;
            va   += c*dtp;
        }
        /* 20 */

        vtot += va;

        cos_theta2 = gmx::square(cos_theta);            /*   1		*/
        if (cos_theta2 < 1)
        {
            int  m;
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            rvec f_i, f_j, f_k;

            st    = dVdt*gmx::invsqrt(1 - cos_theta2); /*  12		*/
            sth   = st*cos_theta;                      /*   1		*/
            nrkj2 = iprod(r_kj, r_kj);                 /*   5		*/
            nrij2 = iprod(r_ij, r_ij);

            cik = st*gmx::invsqrt(nrkj2*nrij2); /*  12		*/
            cii = sth/nrij2;                    /*  10		*/
            ckk = sth/nrkj2;                    /*  10		*/

            for (m = 0; (m < DIM); m++)         /*  39		*/
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                f[ai][m] += f_i[m];
                f[aj][m] += f_j[m];
                f[ak][m] += f_k[m];
            }
            if (g)
            {
                copy_ivec(SHIFT_IVEC(g, aj), jt);

                ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
                ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
                t1 = IVEC2IS(dt_ij);
                t2 = IVEC2IS(dt_kj);
            }
            rvec_inc(fshift[t1], f_i);
            rvec_inc(fshift[CENTRAL], f_j);
            rvec_inc(fshift[t2], f_k);
        }                                       /* 153 TOTAL	*/
    }
    return vtot;
}

real dih_angle(const rvec xi, const rvec xj, const rvec xk, const rvec xl,
               const t_pbc *pbc,
               rvec r_ij, rvec r_kj, rvec r_kl, rvec m, rvec n,
               int *t1, int *t2, int *t3)
{
    *t1 = pbc_rvec_sub(pbc, xi, xj, r_ij); /*  3        */
    *t2 = pbc_rvec_sub(pbc, xk, xj, r_kj); /*  3		*/
    *t3 = pbc_rvec_sub(pbc, xk, xl, r_kl); /*  3		*/

    cprod(r_ij, r_kj, m);                  /*  9        */
    cprod(r_kj, r_kl, n);                  /*  9		*/
    real phi  = gmx_angle(m, n);           /* 49 (assuming 25 for atan2) */
    real ipr  = iprod(r_ij, n);            /*  5        */
    real sign = (ipr < 0.0) ? -1.0 : 1.0;
    phi       = sign*phi;                  /*  1		*/
    /* 82 TOTAL	*/
    return phi;
}


#if GMX_SIMD_HAVE_REAL

/* As dih_angle above, but calculates 4 dihedral angles at once using SIMD,
 * also calculates the pre-factor required for the dihedral force update.
 * Note that bv and buf should be register aligned.
 */
static inline void
dih_angle_simd(const rvec *x,
               const int *ai, const int *aj, const int *ak, const int *al,
               const real *pbc_simd,
               SimdReal *phi_S,
               SimdReal *mx_S, SimdReal *my_S, SimdReal *mz_S,
               SimdReal *nx_S, SimdReal *ny_S, SimdReal *nz_S,
               SimdReal *nrkj_m2_S,
               SimdReal *nrkj_n2_S,
               SimdReal *p_S,
               SimdReal *q_S)
{
    SimdReal xi_S, yi_S, zi_S;
    SimdReal xj_S, yj_S, zj_S;
    SimdReal xk_S, yk_S, zk_S;
    SimdReal xl_S, yl_S, zl_S;
    SimdReal rijx_S, rijy_S, rijz_S;
    SimdReal rkjx_S, rkjy_S, rkjz_S;
    SimdReal rklx_S, rkly_S, rklz_S;
    SimdReal cx_S, cy_S, cz_S;
    SimdReal cn_S;
    SimdReal s_S;
    SimdReal ipr_S;
    SimdReal iprm_S, iprn_S;
    SimdReal nrkj2_S, nrkj_1_S, nrkj_2_S, nrkj_S;
    SimdReal toler_S;
    SimdReal nrkj2_min_S;
    SimdReal real_eps_S;

    /* Used to avoid division by zero.
     * We take into acount that we multiply the result by real_eps_S.
     */
    nrkj2_min_S = SimdReal(GMX_REAL_MIN/(2*GMX_REAL_EPS));

    /* The value of the last significant bit (GMX_REAL_EPS is half of that) */
    real_eps_S  = SimdReal(2*GMX_REAL_EPS);

    /* Store the non PBC corrected distances packed and aligned */
    gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), ai, &xi_S, &yi_S, &zi_S);
    gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), aj, &xj_S, &yj_S, &zj_S);
    gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), ak, &xk_S, &yk_S, &zk_S);
    gatherLoadUTranspose<3>(reinterpret_cast<const real *>(x), al, &xl_S, &yl_S, &zl_S);
    rijx_S = xi_S - xj_S;
    rijy_S = yi_S - yj_S;
    rijz_S = zi_S - zj_S;
    rkjx_S = xk_S - xj_S;
    rkjy_S = yk_S - yj_S;
    rkjz_S = zk_S - zj_S;
    rklx_S = xk_S - xl_S;
    rkly_S = yk_S - yl_S;
    rklz_S = zk_S - zl_S;

    pbc_correct_dx_simd(&rijx_S, &rijy_S, &rijz_S, pbc_simd);
    pbc_correct_dx_simd(&rkjx_S, &rkjy_S, &rkjz_S, pbc_simd);
    pbc_correct_dx_simd(&rklx_S, &rkly_S, &rklz_S, pbc_simd);

    cprod(rijx_S, rijy_S, rijz_S,
          rkjx_S, rkjy_S, rkjz_S,
          mx_S, my_S, mz_S);

    cprod(rkjx_S, rkjy_S, rkjz_S,
          rklx_S, rkly_S, rklz_S,
          nx_S, ny_S, nz_S);

    cprod(*mx_S, *my_S, *mz_S,
          *nx_S, *ny_S, *nz_S,
          &cx_S, &cy_S, &cz_S);

    cn_S       = sqrt(norm2(cx_S, cy_S, cz_S));

    s_S        = iprod(*mx_S, *my_S, *mz_S, *nx_S, *ny_S, *nz_S);

    /* Determine the dihedral angle, the sign might need correction */
    *phi_S     = atan2(cn_S, s_S);

    ipr_S      = iprod(rijx_S, rijy_S, rijz_S,
                       *nx_S, *ny_S, *nz_S);

    iprm_S     = norm2(*mx_S, *my_S, *mz_S);
    iprn_S     = norm2(*nx_S, *ny_S, *nz_S);

    nrkj2_S    = norm2(rkjx_S, rkjy_S, rkjz_S);

    /* Avoid division by zero. When zero, the result is multiplied by 0
     * anyhow, so the 3 max below do not affect the final result.
     */
    nrkj2_S    = max(nrkj2_S, nrkj2_min_S);
    nrkj_1_S   = invsqrt(nrkj2_S);
    nrkj_2_S   = nrkj_1_S * nrkj_1_S;
    nrkj_S     = nrkj2_S * nrkj_1_S;

    toler_S    = nrkj2_S * real_eps_S;

    /* Here the plain-C code uses a conditional, but we can't do that in SIMD.
     * So we take a max with the tolerance instead. Since we multiply with
     * m or n later, the max does not affect the results.
     */
    iprm_S     = max(iprm_S, toler_S);
    iprn_S     = max(iprn_S, toler_S);
    *nrkj_m2_S = nrkj_S * inv(iprm_S);
    *nrkj_n2_S = nrkj_S * inv(iprn_S);

    /* Set sign of phi_S with the sign of ipr_S; phi_S is currently positive */
    *phi_S     = copysign(*phi_S, ipr_S);
    *p_S       = iprod(rijx_S, rijy_S, rijz_S, rkjx_S, rkjy_S, rkjz_S);
    *p_S       = *p_S * nrkj_2_S;

    *q_S       = iprod(rklx_S, rkly_S, rklz_S, rkjx_S, rkjy_S, rkjz_S);
    *q_S       = *q_S * nrkj_2_S;
}

#endif // GMX_SIMD_HAVE_REAL

void do_dih_fup(int i, int j, int k, int l, real ddphi,
                rvec r_ij, rvec r_kj, rvec r_kl,
                rvec m, rvec n, rvec4 f[], rvec fshift[],
                const t_pbc *pbc, const t_graph *g,
                const rvec x[], int t1, int t2, int t3)
{
    /* 143 FLOPS */
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec, dx_jl;
    real iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    real a, b, p, q, toler;
    ivec jt, dt_ij, dt_kj, dt_lj;

    iprm  = iprod(m, m);       /*  5    */
    iprn  = iprod(n, n);       /*  5	*/
    nrkj2 = iprod(r_kj, r_kj); /*  5	*/
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = gmx::invsqrt(nrkj2); /* 10	*/
        nrkj_2 = nrkj_1*nrkj_1;       /*  1	*/
        nrkj   = nrkj2*nrkj_1;        /*  1	*/
        a      = -ddphi*nrkj/iprm;    /* 11	*/
        svmul(a, m, f_i);             /*  3	*/
        b     = ddphi*nrkj/iprn;      /* 11	*/
        svmul(b, n, f_l);             /*  3  */
        p     = iprod(r_ij, r_kj);    /*  5	*/
        p    *= nrkj_2;               /*  1	*/
        q     = iprod(r_kl, r_kj);    /*  5	*/
        q    *= nrkj_2;               /*  1	*/
        svmul(p, f_i, uvec);          /*  3	*/
        svmul(q, f_l, vvec);          /*  3	*/
        rvec_sub(uvec, vvec, svec);   /*  3	*/
        rvec_sub(f_i, svec, f_j);     /*  3	*/
        rvec_add(f_l, svec, f_k);     /*  3	*/
        rvec_inc(f[i], f_i);          /*  3	*/
        rvec_dec(f[j], f_j);          /*  3	*/
        rvec_dec(f[k], f_k);          /*  3	*/
        rvec_inc(f[l], f_l);          /*  3	*/

        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, j), jt);
            ivec_sub(SHIFT_IVEC(g, i), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, k), jt, dt_kj);
            ivec_sub(SHIFT_IVEC(g, l), jt, dt_lj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
            t3 = IVEC2IS(dt_lj);
        }
        else if (pbc)
        {
            t3 = pbc_rvec_sub(pbc, x[l], x[j], dx_jl);
        }
        else
        {
            t3 = CENTRAL;
        }

        rvec_inc(fshift[t1], f_i);
        rvec_dec(fshift[CENTRAL], f_j);
        rvec_dec(fshift[t2], f_k);
        rvec_inc(fshift[t3], f_l);
    }
    /* 112 TOTAL    */
}

/* As do_dih_fup above, but without shift forces */
static void
do_dih_fup_noshiftf(int i, int j, int k, int l, real ddphi,
                    rvec r_ij, rvec r_kj, rvec r_kl,
                    rvec m, rvec n, rvec4 f[])
{
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec;
    real iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    real a, b, p, q, toler;

    iprm  = iprod(m, m);       /*  5    */
    iprn  = iprod(n, n);       /*  5	*/
    nrkj2 = iprod(r_kj, r_kj); /*  5	*/
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = gmx::invsqrt(nrkj2); /* 10	*/
        nrkj_2 = nrkj_1*nrkj_1;       /*  1	*/
        nrkj   = nrkj2*nrkj_1;        /*  1	*/
        a      = -ddphi*nrkj/iprm;    /* 11	*/
        svmul(a, m, f_i);             /*  3	*/
        b     = ddphi*nrkj/iprn;      /* 11	*/
        svmul(b, n, f_l);             /*  3  */
        p     = iprod(r_ij, r_kj);    /*  5	*/
        p    *= nrkj_2;               /*  1	*/
        q     = iprod(r_kl, r_kj);    /*  5	*/
        q    *= nrkj_2;               /*  1	*/
        svmul(p, f_i, uvec);          /*  3	*/
        svmul(q, f_l, vvec);          /*  3	*/
        rvec_sub(uvec, vvec, svec);   /*  3	*/
        rvec_sub(f_i, svec, f_j);     /*  3	*/
        rvec_add(f_l, svec, f_k);     /*  3	*/
        rvec_inc(f[i], f_i);          /*  3	*/
        rvec_dec(f[j], f_j);          /*  3	*/
        rvec_dec(f[k], f_k);          /*  3	*/
        rvec_inc(f[l], f_l);          /*  3	*/
    }
}

#if GMX_SIMD_HAVE_REAL
/* As do_dih_fup_noshiftf above, but with SIMD and pre-calculated pre-factors */
static inline void gmx_simdcall
do_dih_fup_noshiftf_simd(const int *ai, const int *aj, const int *ak, const int *al,
                         SimdReal p, SimdReal q,
                         SimdReal f_i_x,  SimdReal f_i_y,  SimdReal f_i_z,
                         SimdReal mf_l_x, SimdReal mf_l_y, SimdReal mf_l_z,
                         rvec4 f[])
{
    SimdReal sx    = p * f_i_x + q * mf_l_x;
    SimdReal sy    = p * f_i_y + q * mf_l_y;
    SimdReal sz    = p * f_i_z + q * mf_l_z;
    SimdReal f_j_x = f_i_x - sx;
    SimdReal f_j_y = f_i_y - sy;
    SimdReal f_j_z = f_i_z - sz;
    SimdReal f_k_x = mf_l_x - sx;
    SimdReal f_k_y = mf_l_y - sy;
    SimdReal f_k_z = mf_l_z - sz;
    transposeScatterIncrU<4>(reinterpret_cast<real *>(f), ai, f_i_x, f_i_y, f_i_z);
    transposeScatterDecrU<4>(reinterpret_cast<real *>(f), aj, f_j_x, f_j_y, f_j_z);
    transposeScatterIncrU<4>(reinterpret_cast<real *>(f), ak, f_k_x, f_k_y, f_k_z);
    transposeScatterDecrU<4>(reinterpret_cast<real *>(f), al, mf_l_x, mf_l_y, mf_l_z);
}
#endif // GMX_SIMD_HAVE_REAL

static real dopdihs(real cpA, real cpB, real phiA, real phiB, int mult,
                    real phi, real lambda, real *V, real *F)
{
    real v, dvdlambda, mdphi, v1, sdphi, ddphi;
    real L1   = 1.0 - lambda;
    real ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
    real dph0 = (phiB - phiA)*DEG2RAD;
    real cp   = L1*cpA + lambda*cpB;

    mdphi =  mult*phi - ph0;
    sdphi = std::sin(mdphi);
    ddphi = -cp*mult*sdphi;
    v1    = 1.0 + std::cos(mdphi);
    v     = cp*v1;

    dvdlambda  = (cpB - cpA)*v1 + cp*dph0*sdphi;

    *V = v;
    *F = ddphi;

    return dvdlambda;

    /* That was 40 flops */
}

static void
dopdihs_noener(real cpA, real cpB, real phiA, real phiB, int mult,
               real phi, real lambda, real *F)
{
    real mdphi, sdphi, ddphi;
    real L1   = 1.0 - lambda;
    real ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
    real cp   = L1*cpA + lambda*cpB;

    mdphi = mult*phi - ph0;
    sdphi = std::sin(mdphi);
    ddphi = -cp*mult*sdphi;

    *F = ddphi;

    /* That was 20 flops */
}

static real dopdihs_min(real cpA, real cpB, real phiA, real phiB, int mult,
                        real phi, real lambda, real *V, real *F)
/* similar to dopdihs, except for a minus sign  *
 * and a different treatment of mult/phi0       */
{
    real v, dvdlambda, mdphi, v1, sdphi, ddphi;
    real L1   = 1.0 - lambda;
    real ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
    real dph0 = (phiB - phiA)*DEG2RAD;
    real cp   = L1*cpA + lambda*cpB;

    mdphi = mult*(phi-ph0);
    sdphi = std::sin(mdphi);
    ddphi = cp*mult*sdphi;
    v1    = 1.0-std::cos(mdphi);
    v     = cp*v1;

    dvdlambda  = (cpB-cpA)*v1 + cp*dph0*sdphi;

    *V = v;
    *F = ddphi;

    return dvdlambda;

    /* That was 40 flops */
}

real pdihs(int nbonds,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec4 f[], rvec fshift[],
           const t_pbc *pbc, const t_graph *g,
           real lambda, real *dvdlambda,
           const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
           int gmx_unused *global_atom_index)
{
    int  i, type, ai, aj, ak, al;
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi, ddphi, vpd, vtot;

    vtot = 0.0;

    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];

        phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                        &t1, &t2, &t3);  /*  84      */
        *dvdlambda += dopdihs(forceparams[type].pdihs.cpA,
                              forceparams[type].pdihs.cpB,
                              forceparams[type].pdihs.phiA,
                              forceparams[type].pdihs.phiB,
                              forceparams[type].pdihs.mult,
                              phi, lambda, &vpd, &ddphi);

        vtot += vpd;
        do_dih_fup(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc, g, x, t1, t2, t3); /* 112		*/

    }                                                 /* 223 TOTAL  */

    return vtot;
}

void make_dp_periodic(real *dp)  /* 1 flop? */
{
    /* dp cannot be outside (-pi,pi) */
    if (*dp >= M_PI)
    {
        *dp -= 2*M_PI;
    }
    else if (*dp < -M_PI)
    {
        *dp += 2*M_PI;
    }
}

/* As pdihs above, but without calculating energies and shift forces */
void
pdihs_noener(int nbonds,
             const t_iatom forceatoms[], const t_iparams forceparams[],
             const rvec x[], rvec4 f[],
             const t_pbc gmx_unused *pbc, const t_graph gmx_unused *g,
             real lambda,
             const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
             int gmx_unused *global_atom_index)
{
    int  i, type, ai, aj, ak, al;
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi, ddphi_tot, ddphi;

    for (i = 0; (i < nbonds); )
    {
        ai   = forceatoms[i+1];
        aj   = forceatoms[i+2];
        ak   = forceatoms[i+3];
        al   = forceatoms[i+4];

        phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                        &t1, &t2, &t3);

        ddphi_tot = 0;

        /* Loop over dihedrals working on the same atoms,
         * so we avoid recalculating angles and force distributions.
         */
        do
        {
            type = forceatoms[i];
            dopdihs_noener(forceparams[type].pdihs.cpA,
                           forceparams[type].pdihs.cpB,
                           forceparams[type].pdihs.phiA,
                           forceparams[type].pdihs.phiB,
                           forceparams[type].pdihs.mult,
                           phi, lambda, &ddphi);
            ddphi_tot += ddphi;

            i += 5;
        }
        while (i < nbonds &&
               forceatoms[i+1] == ai &&
               forceatoms[i+2] == aj &&
               forceatoms[i+3] == ak &&
               forceatoms[i+4] == al);

        do_dih_fup_noshiftf(ai, aj, ak, al, ddphi_tot, r_ij, r_kj, r_kl, m, n, f);
    }
}


#if GMX_SIMD_HAVE_REAL

/* As pdihs_noner above, but using SIMD to calculate many dihedrals at once */
void
pdihs_noener_simd(int nbonds,
                  const t_iatom forceatoms[], const t_iparams forceparams[],
                  const rvec x[], rvec4 f[],
                  const t_pbc *pbc, const t_graph gmx_unused *g,
                  real gmx_unused lambda,
                  const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                  int gmx_unused *global_atom_index)
{
    const int             nfa1 = 5;
    int                   i, iu, s;
    int                   type;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    ai[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    aj[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    ak[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    al[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) real            buf[3*GMX_SIMD_REAL_WIDTH];
    real                 *cp, *phi0, *mult;
    SimdReal              deg2rad_S(DEG2RAD);
    SimdReal              p_S, q_S;
    SimdReal              phi0_S, phi_S;
    SimdReal              mx_S, my_S, mz_S;
    SimdReal              nx_S, ny_S, nz_S;
    SimdReal              nrkj_m2_S, nrkj_n2_S;
    SimdReal              cp_S, mdphi_S, mult_S;
    SimdReal              sin_S, cos_S;
    SimdReal              mddphi_S;
    SimdReal              sf_i_S, msf_l_S;
    alignas(GMX_SIMD_ALIGNMENT) real            pbc_simd[9*GMX_SIMD_REAL_WIDTH];

    /* Extract aligned pointer for parameters and variables */
    cp    = buf + 0*GMX_SIMD_REAL_WIDTH;
    phi0  = buf + 1*GMX_SIMD_REAL_WIDTH;
    mult  = buf + 2*GMX_SIMD_REAL_WIDTH;

    set_pbc_simd(pbc, pbc_simd);

    /* nbonds is the number of dihedrals times nfa1, here we step GMX_SIMD_REAL_WIDTH dihs */
    for (i = 0; (i < nbonds); i += GMX_SIMD_REAL_WIDTH*nfa1)
    {
        /* Collect atoms quadruplets for GMX_SIMD_REAL_WIDTH dihedrals.
         * iu indexes into forceatoms, we should not let iu go beyond nbonds.
         */
        iu = i;
        for (s = 0; s < GMX_SIMD_REAL_WIDTH; s++)
        {
            type  = forceatoms[iu];
            ai[s] = forceatoms[iu+1];
            aj[s] = forceatoms[iu+2];
            ak[s] = forceatoms[iu+3];
            al[s] = forceatoms[iu+4];

            /* At the end fill the arrays with the last atoms and 0 params */
            if (i + s*nfa1 < nbonds)
            {
                cp[s]   = forceparams[type].pdihs.cpA;
                phi0[s] = forceparams[type].pdihs.phiA;
                mult[s] = forceparams[type].pdihs.mult;

                if (iu + nfa1 < nbonds)
                {
                    iu += nfa1;
                }
            }
            else
            {
                cp[s]   = 0;
                phi0[s] = 0;
                mult[s] = 0;
            }
        }

        /* Caclulate GMX_SIMD_REAL_WIDTH dihedral angles at once */
        dih_angle_simd(x, ai, aj, ak, al, pbc_simd,
                       &phi_S,
                       &mx_S, &my_S, &mz_S,
                       &nx_S, &ny_S, &nz_S,
                       &nrkj_m2_S,
                       &nrkj_n2_S,
                       &p_S, &q_S);

        cp_S     = load<SimdReal>(cp);
        phi0_S   = load<SimdReal>(phi0) * deg2rad_S;
        mult_S   = load<SimdReal>(mult);

        mdphi_S  = fms(mult_S, phi_S, phi0_S);

        /* Calculate GMX_SIMD_REAL_WIDTH sines at once */
        sincos(mdphi_S, &sin_S, &cos_S);
        mddphi_S = cp_S * mult_S * sin_S;
        sf_i_S   = mddphi_S * nrkj_m2_S;
        msf_l_S  = mddphi_S * nrkj_n2_S;

        /* After this m?_S will contain f[i] */
        mx_S     = sf_i_S * mx_S;
        my_S     = sf_i_S * my_S;
        mz_S     = sf_i_S * mz_S;

        /* After this m?_S will contain -f[l] */
        nx_S     = msf_l_S * nx_S;
        ny_S     = msf_l_S * ny_S;
        nz_S     = msf_l_S * nz_S;

        do_dih_fup_noshiftf_simd(ai, aj, ak, al,
                                 p_S, q_S,
                                 mx_S, my_S, mz_S,
                                 nx_S, ny_S, nz_S,
                                 f);
    }
}

/* This is mostly a copy of pdihs_noener_simd above, but with using
 * the RB potential instead of a harmonic potential.
 * This function can replace rbdihs() when no energy and virial are needed.
 */
void
rbdihs_noener_simd(int nbonds,
                   const t_iatom forceatoms[], const t_iparams forceparams[],
                   const rvec x[], rvec4 f[],
                   const t_pbc *pbc, const t_graph gmx_unused *g,
                   real gmx_unused lambda,
                   const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                   int gmx_unused *global_atom_index)
{
    const int             nfa1 = 5;
    int                   i, iu, s, j;
    int                   type;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t  ai[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t  aj[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t  ak[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t  al[GMX_SIMD_REAL_WIDTH];
    alignas(GMX_SIMD_ALIGNMENT) real          parm[NR_RBDIHS*GMX_SIMD_REAL_WIDTH];

    SimdReal              p_S, q_S;
    SimdReal              phi_S;
    SimdReal              ddphi_S, cosfac_S;
    SimdReal              mx_S, my_S, mz_S;
    SimdReal              nx_S, ny_S, nz_S;
    SimdReal              nrkj_m2_S, nrkj_n2_S;
    SimdReal              parm_S, c_S;
    SimdReal              sin_S, cos_S;
    SimdReal              sf_i_S, msf_l_S;
    alignas(GMX_SIMD_ALIGNMENT) real          pbc_simd[9*GMX_SIMD_REAL_WIDTH];

    SimdReal              pi_S(M_PI);
    SimdReal              one_S(1.0);

    set_pbc_simd(pbc, pbc_simd);

    /* nbonds is the number of dihedrals times nfa1, here we step GMX_SIMD_REAL_WIDTH dihs */
    for (i = 0; (i < nbonds); i += GMX_SIMD_REAL_WIDTH*nfa1)
    {
        /* Collect atoms quadruplets for GMX_SIMD_REAL_WIDTH dihedrals.
         * iu indexes into forceatoms, we should not let iu go beyond nbonds.
         */
        iu = i;
        for (s = 0; s < GMX_SIMD_REAL_WIDTH; s++)
        {
            type  = forceatoms[iu];
            ai[s] = forceatoms[iu+1];
            aj[s] = forceatoms[iu+2];
            ak[s] = forceatoms[iu+3];
            al[s] = forceatoms[iu+4];

            /* At the end fill the arrays with the last atoms and 0 params */
            if (i + s*nfa1 < nbonds)
            {
                /* We don't need the first parameter, since that's a constant
                 * which only affects the energies, not the forces.
                 */
                for (j = 1; j < NR_RBDIHS; j++)
                {
                    parm[j*GMX_SIMD_REAL_WIDTH + s] =
                        forceparams[type].rbdihs.rbcA[j];
                }

                if (iu + nfa1 < nbonds)
                {
                    iu += nfa1;
                }
            }
            else
            {
                for (j = 1; j < NR_RBDIHS; j++)
                {
                    parm[j*GMX_SIMD_REAL_WIDTH + s] = 0;
                }
            }
        }

        /* Caclulate GMX_SIMD_REAL_WIDTH dihedral angles at once */
        dih_angle_simd(x, ai, aj, ak, al, pbc_simd,
                       &phi_S,
                       &mx_S, &my_S, &mz_S,
                       &nx_S, &ny_S, &nz_S,
                       &nrkj_m2_S,
                       &nrkj_n2_S,
                       &p_S, &q_S);

        /* Change to polymer convention */
        phi_S = phi_S - pi_S;

        sincos(phi_S, &sin_S, &cos_S);

        ddphi_S   = setZero();
        c_S       = one_S;
        cosfac_S  = one_S;
        for (j = 1; j < NR_RBDIHS; j++)
        {
            parm_S   = load<SimdReal>(parm + j*GMX_SIMD_REAL_WIDTH);
            ddphi_S  = fma(c_S * parm_S, cosfac_S, ddphi_S);
            cosfac_S = cosfac_S * cos_S;
            c_S      = c_S + one_S;
        }

        /* Note that here we do not use the minus sign which is present
         * in the normal RB code. This is corrected for through (m)sf below.
         */
        ddphi_S  = ddphi_S * sin_S;

        sf_i_S   = ddphi_S * nrkj_m2_S;
        msf_l_S  = ddphi_S * nrkj_n2_S;

        /* After this m?_S will contain f[i] */
        mx_S     = sf_i_S * mx_S;
        my_S     = sf_i_S * my_S;
        mz_S     = sf_i_S * mz_S;

        /* After this m?_S will contain -f[l] */
        nx_S     = msf_l_S * nx_S;
        ny_S     = msf_l_S * ny_S;
        nz_S     = msf_l_S * nz_S;

        do_dih_fup_noshiftf_simd(ai, aj, ak, al,
                                 p_S, q_S,
                                 mx_S, my_S, mz_S,
                                 nx_S, ny_S, nz_S,
                                 f);
    }
}

#endif // GMX_SIMD_HAVE_REAL


real idihs(int nbonds,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec4 f[], rvec fshift[],
           const t_pbc *pbc, const t_graph *g,
           real lambda, real *dvdlambda,
           const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
           int gmx_unused *global_atom_index)
{
    int  i, type, ai, aj, ak, al;
    int  t1, t2, t3;
    real phi, phi0, dphi0, ddphi, vtot;
    rvec r_ij, r_kj, r_kl, m, n;
    real L1, kk, dp, dp2, kA, kB, pA, pB, dvdl_term;

    L1        = 1.0-lambda;
    dvdl_term = 0;
    vtot      = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];

        phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                        &t1, &t2, &t3);  /*  84		*/

        /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
         * force changes if we just apply a normal harmonic.
         * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
         * This means we will never have the periodicity problem, unless
         * the dihedral is Pi away from phiO, which is very unlikely due to
         * the potential.
         */
        kA = forceparams[type].harmonic.krA;
        kB = forceparams[type].harmonic.krB;
        pA = forceparams[type].harmonic.rA;
        pB = forceparams[type].harmonic.rB;

        kk    = L1*kA + lambda*kB;
        phi0  = (L1*pA + lambda*pB)*DEG2RAD;
        dphi0 = (pB - pA)*DEG2RAD;

        dp = phi-phi0;

        make_dp_periodic(&dp);

        dp2 = dp*dp;

        vtot += 0.5*kk*dp2;
        ddphi = -kk*dp;

        dvdl_term += 0.5*(kB - kA)*dp2 - kk*dphi0*dp;

        do_dih_fup(ai, aj, ak, al, -ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc, g, x, t1, t2, t3); /* 112		*/
        /* 218 TOTAL	*/
    }

    *dvdlambda += dvdl_term;
    return vtot;
}

static real low_angres(int nbonds,
                       const t_iatom forceatoms[], const t_iparams forceparams[],
                       const rvec x[], rvec4 f[], rvec fshift[],
                       const t_pbc *pbc, const t_graph *g,
                       real lambda, real *dvdlambda,
                       gmx_bool bZAxis)
{
    int  i, m, type, ai, aj, ak, al;
    int  t1, t2;
    real phi, cos_phi, cos_phi2, vid, vtot, dVdphi;
    rvec r_ij, r_kl, f_i, f_k = {0, 0, 0};
    real st, sth, nrij2, nrkl2, c, cij, ckl;

    ivec dt;
    t2 = 0; /* avoid warning with gcc-3.3. It is never used uninitialized */

    vtot = 0.0;
    ak   = al = 0; /* to avoid warnings */
    for (i = 0; i < nbonds; )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        t1   = pbc_rvec_sub(pbc, x[aj], x[ai], r_ij);       /*  3		*/
        if (!bZAxis)
        {
            ak   = forceatoms[i++];
            al   = forceatoms[i++];
            t2   = pbc_rvec_sub(pbc, x[al], x[ak], r_kl);  /*  3		*/
        }
        else
        {
            r_kl[XX] = 0;
            r_kl[YY] = 0;
            r_kl[ZZ] = 1;
        }

        cos_phi = cos_angle(r_ij, r_kl); /* 25		*/
        phi     = std::acos(cos_phi);    /* 10           */

        *dvdlambda += dopdihs_min(forceparams[type].pdihs.cpA,
                                  forceparams[type].pdihs.cpB,
                                  forceparams[type].pdihs.phiA,
                                  forceparams[type].pdihs.phiB,
                                  forceparams[type].pdihs.mult,
                                  phi, lambda, &vid, &dVdphi); /*  40  */

        vtot += vid;

        cos_phi2 = gmx::square(cos_phi);                /*   1		*/
        if (cos_phi2 < 1)
        {
            st    = -dVdphi*gmx::invsqrt(1 - cos_phi2); /*  12		*/
            sth   = st*cos_phi;                         /*   1		*/
            nrij2 = iprod(r_ij, r_ij);                  /*   5		*/
            nrkl2 = iprod(r_kl, r_kl);                  /*   5          */

            c   = st*gmx::invsqrt(nrij2*nrkl2);         /*  11		*/
            cij = sth/nrij2;                            /*  10		*/
            ckl = sth/nrkl2;                            /*  10		*/

            for (m = 0; m < DIM; m++)                   /*  18+18       */
            {
                f_i[m]    = (c*r_kl[m]-cij*r_ij[m]);
                f[ai][m] += f_i[m];
                f[aj][m] -= f_i[m];
                if (!bZAxis)
                {
                    f_k[m]    = (c*r_ij[m]-ckl*r_kl[m]);
                    f[ak][m] += f_k[m];
                    f[al][m] -= f_k[m];
                }
            }

            if (g)
            {
                ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
                t1 = IVEC2IS(dt);
            }
            rvec_inc(fshift[t1], f_i);
            rvec_dec(fshift[CENTRAL], f_i);
            if (!bZAxis)
            {
                if (g)
                {
                    ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, al), dt);
                    t2 = IVEC2IS(dt);
                }
                rvec_inc(fshift[t2], f_k);
                rvec_dec(fshift[CENTRAL], f_k);
            }
        }
    }

    return vtot; /*  184 / 157 (bZAxis)  total  */
}

real angres(int nbonds,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec4 f[], rvec fshift[],
            const t_pbc *pbc, const t_graph *g,
            real lambda, real *dvdlambda,
            const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
            int gmx_unused *global_atom_index)
{
    return low_angres(nbonds, forceatoms, forceparams, x, f, fshift, pbc, g,
                      lambda, dvdlambda, FALSE);
}

real angresz(int nbonds,
             const t_iatom forceatoms[], const t_iparams forceparams[],
             const rvec x[], rvec4 f[], rvec fshift[],
             const t_pbc *pbc, const t_graph *g,
             real lambda, real *dvdlambda,
             const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
             int gmx_unused *global_atom_index)
{
    return low_angres(nbonds, forceatoms, forceparams, x, f, fshift, pbc, g,
                      lambda, dvdlambda, TRUE);
}

real dihres(int nbonds,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec4 f[], rvec fshift[],
            const t_pbc *pbc, const t_graph *g,
            real lambda, real *dvdlambda,
            const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
            int gmx_unused  *global_atom_index)
{
    real vtot = 0;
    int  ai, aj, ak, al, i, type, t1, t2, t3;
    real phi0A, phi0B, dphiA, dphiB, kfacA, kfacB, phi0, dphi, kfac;
    real phi, ddphi, ddp, ddp2, dp, d2r, L1;
    rvec r_ij, r_kj, r_kl, m, n;

    L1 = 1.0-lambda;

    d2r = DEG2RAD;

    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];

        phi0A  = forceparams[type].dihres.phiA*d2r;
        dphiA  = forceparams[type].dihres.dphiA*d2r;
        kfacA  = forceparams[type].dihres.kfacA;

        phi0B  = forceparams[type].dihres.phiB*d2r;
        dphiB  = forceparams[type].dihres.dphiB*d2r;
        kfacB  = forceparams[type].dihres.kfacB;

        phi0  = L1*phi0A + lambda*phi0B;
        dphi  = L1*dphiA + lambda*dphiB;
        kfac  = L1*kfacA + lambda*kfacB;

        phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                        &t1, &t2, &t3);
        /* 84 flops */

        /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
         * force changes if we just apply a normal harmonic.
         * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
         * This means we will never have the periodicity problem, unless
         * the dihedral is Pi away from phiO, which is very unlikely due to
         * the potential.
         */
        dp = phi-phi0;
        make_dp_periodic(&dp);

        if (dp > dphi)
        {
            ddp = dp-dphi;
        }
        else if (dp < -dphi)
        {
            ddp = dp+dphi;
        }
        else
        {
            ddp = 0;
        }

        if (ddp != 0.0)
        {
            ddp2  = ddp*ddp;
            vtot += 0.5*kfac*ddp2;
            ddphi = kfac*ddp;

            *dvdlambda += 0.5*(kfacB - kfacA)*ddp2;
            /* lambda dependence from changing restraint distances */
            if (ddp > 0)
            {
                *dvdlambda -= kfac*ddp*((dphiB - dphiA)+(phi0B - phi0A));
            }
            else if (ddp < 0)
            {
                *dvdlambda += kfac*ddp*((dphiB - dphiA)-(phi0B - phi0A));
            }
            do_dih_fup(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                       f, fshift, pbc, g, x, t1, t2, t3);      /* 112		*/
        }
    }
    return vtot;
}


real unimplemented(int gmx_unused nbonds,
                   const t_iatom gmx_unused forceatoms[], const t_iparams gmx_unused forceparams[],
                   const rvec gmx_unused x[], rvec4 gmx_unused f[], rvec gmx_unused fshift[],
                   const t_pbc gmx_unused *pbc, const t_graph  gmx_unused *g,
                   real gmx_unused lambda, real gmx_unused *dvdlambda,
                   const t_mdatoms  gmx_unused *md, t_fcdata gmx_unused *fcd,
                   int gmx_unused *global_atom_index)
{
    gmx_impl("*** you are using a not implemented function");
}

real restrangles(int nbonds,
                 const t_iatom forceatoms[], const t_iparams forceparams[],
                 const rvec x[], rvec4 f[], rvec fshift[],
                 const t_pbc *pbc, const t_graph *g,
                 real gmx_unused lambda, real gmx_unused *dvdlambda,
                 const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                 int gmx_unused *global_atom_index)
{
    int  i, d, ai, aj, ak, type, m;
    int  t1, t2;
    real v, vtot;
    ivec jt, dt_ij, dt_kj;
    rvec f_i, f_j, f_k;
    real prefactor, ratio_ante, ratio_post;
    rvec delta_ante, delta_post, vec_temp;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];

        t1 = pbc_rvec_sub(pbc, x[ai], x[aj], vec_temp);
        pbc_rvec_sub(pbc, x[aj], x[ai], delta_ante);
        t2 = pbc_rvec_sub(pbc, x[ak], x[aj], delta_post);


        /* This function computes factors needed for restricted angle potential.
         * The restricted angle potential is used in coarse-grained simulations to avoid singularities
         * when three particles align and the dihedral angle and dihedral potential
         * cannot be calculated. This potential is calculated using the formula:
           real restrangles(int nbonds,
            const t_iatom forceatoms[],const t_iparams forceparams[],
            const rvec x[],rvec4 f[],rvec fshift[],
            const t_pbc *pbc,const t_graph *g,
            real gmx_unused lambda,real gmx_unused *dvdlambda,
            const t_mdatoms gmx_unused *md,t_fcdata gmx_unused *fcd,
            int gmx_unused *global_atom_index)
           {
           int  i, d, ai, aj, ak, type, m;
           int t1, t2;
           rvec r_ij,r_kj;
           real v, vtot;
           ivec jt,dt_ij,dt_kj;
           rvec f_i, f_j, f_k;
           real prefactor, ratio_ante, ratio_post;
           rvec delta_ante, delta_post, vec_temp;

           vtot = 0.0;
           for(i=0; (i<nbonds); )
           {
           type = forceatoms[i++];
           ai   = forceatoms[i++];
           aj   = forceatoms[i++];
           ak   = forceatoms[i++];

         * \f[V_{\rm ReB}(\theta_i) = \frac{1}{2} k_{\theta} \frac{(\cos\theta_i - \cos\theta_0)^2}
         * {\sin^2\theta_i}\f] ({eq:ReB} and ref \cite{MonicaGoga2013} from the manual).
         * For more explanations see comments file "restcbt.h". */

        compute_factors_restangles(type, forceparams,  delta_ante, delta_post,
                                   &prefactor, &ratio_ante, &ratio_post, &v);

        /*   Forces are computed per component */
        for (d = 0; d < DIM; d++)
        {
            f_i[d] = prefactor * (ratio_ante * delta_ante[d] - delta_post[d]);
            f_j[d] = prefactor * ((ratio_post + 1.0) * delta_post[d] - (ratio_ante + 1.0) * delta_ante[d]);
            f_k[d] = prefactor * (delta_ante[d] - ratio_post * delta_post[d]);
        }

        /*   Computation of potential energy   */

        vtot += v;

        /*   Update forces */

        for (m = 0; (m < DIM); m++)
        {
            f[ai][m] += f_i[m];
            f[aj][m] += f_j[m];
            f[ak][m] += f_k[m];
        }

        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);
            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
        }

        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);
    }
    return vtot;
}


real restrdihs(int nbonds,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec4 f[], rvec fshift[],
               const t_pbc *pbc, const t_graph *g,
               real gmx_unused lambda, real gmx_unused *dvlambda,
               const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
               int gmx_unused *global_atom_index)
{
    int  i, d, type, ai, aj, ak, al;
    rvec f_i, f_j, f_k, f_l;
    rvec dx_jl;
    ivec jt, dt_ij, dt_kj, dt_lj;
    int  t1, t2, t3;
    real v, vtot;
    rvec delta_ante,  delta_crnt, delta_post, vec_temp;
    real factor_phi_ai_ante, factor_phi_ai_crnt, factor_phi_ai_post;
    real factor_phi_aj_ante, factor_phi_aj_crnt, factor_phi_aj_post;
    real factor_phi_ak_ante, factor_phi_ak_crnt, factor_phi_ak_post;
    real factor_phi_al_ante, factor_phi_al_crnt, factor_phi_al_post;
    real prefactor_phi;


    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];

        t1 = pbc_rvec_sub(pbc, x[ai], x[aj], vec_temp);
        pbc_rvec_sub(pbc, x[aj], x[ai], delta_ante);
        t2 = pbc_rvec_sub(pbc, x[ak], x[aj], delta_crnt);
        pbc_rvec_sub(pbc, x[ak], x[al], vec_temp);
        pbc_rvec_sub(pbc, x[al], x[ak], delta_post);

        /* This function computes factors needed for restricted angle potential.
         * The restricted angle potential is used in coarse-grained simulations to avoid singularities
         * when three particles align and the dihedral angle and dihedral potential cannot be calculated.
         * This potential is calculated using the formula:
         * \f[V_{\rm ReB}(\theta_i) = \frac{1}{2} k_{\theta}
         * \frac{(\cos\theta_i - \cos\theta_0)^2}{\sin^2\theta_i}\f]
         * ({eq:ReB} and ref \cite{MonicaGoga2013} from the manual).
         * For more explanations see comments file "restcbt.h" */

        compute_factors_restrdihs(type, forceparams,
                                  delta_ante, delta_crnt, delta_post,
                                  &factor_phi_ai_ante, &factor_phi_ai_crnt, &factor_phi_ai_post,
                                  &factor_phi_aj_ante, &factor_phi_aj_crnt, &factor_phi_aj_post,
                                  &factor_phi_ak_ante, &factor_phi_ak_crnt, &factor_phi_ak_post,
                                  &factor_phi_al_ante, &factor_phi_al_crnt, &factor_phi_al_post,
                                  &prefactor_phi, &v);


        /*      Computation of forces per component */
        for (d = 0; d < DIM; d++)
        {
            f_i[d] = prefactor_phi * (factor_phi_ai_ante * delta_ante[d] + factor_phi_ai_crnt * delta_crnt[d] + factor_phi_ai_post * delta_post[d]);
            f_j[d] = prefactor_phi * (factor_phi_aj_ante * delta_ante[d] + factor_phi_aj_crnt * delta_crnt[d] + factor_phi_aj_post * delta_post[d]);
            f_k[d] = prefactor_phi * (factor_phi_ak_ante * delta_ante[d] + factor_phi_ak_crnt * delta_crnt[d] + factor_phi_ak_post * delta_post[d]);
            f_l[d] = prefactor_phi * (factor_phi_al_ante * delta_ante[d] + factor_phi_al_crnt * delta_crnt[d] + factor_phi_al_post * delta_post[d]);
        }
        /*      Computation of the energy */

        vtot += v;



        /*    Updating the forces */

        rvec_inc(f[ai], f_i);
        rvec_inc(f[aj], f_j);
        rvec_inc(f[ak], f_k);
        rvec_inc(f[al], f_l);


        /* Updating the fshift forces for the pressure coupling */
        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);
            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            ivec_sub(SHIFT_IVEC(g, al), jt, dt_lj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
            t3 = IVEC2IS(dt_lj);
        }
        else if (pbc)
        {
            t3 = pbc_rvec_sub(pbc, x[al], x[aj], dx_jl);
        }
        else
        {
            t3 = CENTRAL;
        }

        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);
        rvec_inc(fshift[t3], f_l);

    }

    return vtot;
}


real cbtdihs(int nbonds,
             const t_iatom forceatoms[], const t_iparams forceparams[],
             const rvec x[], rvec4 f[], rvec fshift[],
             const t_pbc *pbc, const t_graph *g,
             real gmx_unused lambda, real gmx_unused *dvdlambda,
             const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
             int gmx_unused *global_atom_index)
{
    int  type, ai, aj, ak, al, i, d;
    int  t1, t2, t3;
    real v, vtot;
    rvec vec_temp;
    rvec f_i, f_j, f_k, f_l;
    ivec jt, dt_ij, dt_kj, dt_lj;
    rvec dx_jl;
    rvec delta_ante, delta_crnt, delta_post;
    rvec f_phi_ai, f_phi_aj, f_phi_ak, f_phi_al;
    rvec f_theta_ante_ai, f_theta_ante_aj, f_theta_ante_ak;
    rvec f_theta_post_aj, f_theta_post_ak, f_theta_post_al;




    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];


        t1 = pbc_rvec_sub(pbc, x[ai], x[aj], vec_temp);
        pbc_rvec_sub(pbc, x[aj], x[ai], delta_ante);
        t2 = pbc_rvec_sub(pbc, x[ak], x[aj], vec_temp);
        pbc_rvec_sub(pbc, x[ak], x[aj], delta_crnt);
        pbc_rvec_sub(pbc, x[ak], x[al], vec_temp);
        pbc_rvec_sub(pbc, x[al], x[ak], delta_post);

        /* \brief Compute factors for CBT potential
         * The combined bending-torsion potential goes to zero in a very smooth manner, eliminating the numerical
         * instabilities, when three coarse-grained particles align and the dihedral angle and standard
         * dihedral potentials cannot be calculated. The CBT potential is calculated using the formula:
         * \f[V_{\rm CBT}(\theta_{i-1}, \theta_i, \phi_i) = k_{\phi} \sin^3\theta_{i-1} \sin^3\theta_{i}
         * \sum_{n=0}^4 { a_n \cos^n\phi_i}\f] ({eq:CBT} and ref \cite{MonicaGoga2013} from the manual).
         * It contains in its expression not only the dihedral angle \f$\phi\f$
         * but also \f[\theta_{i-1}\f] (theta_ante bellow) and \f[\theta_{i}\f] (theta_post bellow)
         * --- the adjacent bending angles.
         * For more explanations see comments file "restcbt.h". */

        compute_factors_cbtdihs(type, forceparams, delta_ante, delta_crnt, delta_post,
                                f_phi_ai, f_phi_aj, f_phi_ak, f_phi_al,
                                f_theta_ante_ai, f_theta_ante_aj, f_theta_ante_ak,
                                f_theta_post_aj, f_theta_post_ak, f_theta_post_al,
                                &v);


        /*      Acumulate the resuts per beads */
        for (d = 0; d < DIM; d++)
        {
            f_i[d] = f_phi_ai[d] + f_theta_ante_ai[d];
            f_j[d] = f_phi_aj[d] + f_theta_ante_aj[d] + f_theta_post_aj[d];
            f_k[d] = f_phi_ak[d] + f_theta_ante_ak[d] + f_theta_post_ak[d];
            f_l[d] = f_phi_al[d] + f_theta_post_al[d];
        }

        /*      Compute the potential energy */

        vtot += v;


        /*  Updating the forces */
        rvec_inc(f[ai], f_i);
        rvec_inc(f[aj], f_j);
        rvec_inc(f[ak], f_k);
        rvec_inc(f[al], f_l);


        /* Updating the fshift forces for the pressure coupling */
        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);
            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            ivec_sub(SHIFT_IVEC(g, al), jt, dt_lj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
            t3 = IVEC2IS(dt_lj);
        }
        else if (pbc)
        {
            t3 = pbc_rvec_sub(pbc, x[al], x[aj], dx_jl);
        }
        else
        {
            t3 = CENTRAL;
        }

        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);
        rvec_inc(fshift[t3], f_l);
    }

    return vtot;
}

real rbdihs(int nbonds,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec4 f[], rvec fshift[],
            const t_pbc *pbc, const t_graph *g,
            real lambda, real *dvdlambda,
            const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
            int gmx_unused *global_atom_index)
{
    const real c0 = 0.0, c1 = 1.0, c2 = 2.0, c3 = 3.0, c4 = 4.0, c5 = 5.0;
    int        type, ai, aj, ak, al, i, j;
    int        t1, t2, t3;
    rvec       r_ij, r_kj, r_kl, m, n;
    real       parmA[NR_RBDIHS];
    real       parmB[NR_RBDIHS];
    real       parm[NR_RBDIHS];
    real       cos_phi, phi, rbp, rbpBA;
    real       v, ddphi, sin_phi;
    real       cosfac, vtot;
    real       L1        = 1.0-lambda;
    real       dvdl_term = 0;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];

        phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                        &t1, &t2, &t3);  /*  84		*/

        /* Change to polymer convention */
        if (phi < c0)
        {
            phi += M_PI;
        }
        else
        {
            phi -= M_PI;    /*   1		*/

        }
        cos_phi = std::cos(phi);
        /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
        sin_phi = std::sin(phi);

        for (j = 0; (j < NR_RBDIHS); j++)
        {
            parmA[j] = forceparams[type].rbdihs.rbcA[j];
            parmB[j] = forceparams[type].rbdihs.rbcB[j];
            parm[j]  = L1*parmA[j]+lambda*parmB[j];
        }
        /* Calculate cosine powers */
        /* Calculate the energy */
        /* Calculate the derivative */

        v            = parm[0];
        dvdl_term   += (parmB[0]-parmA[0]);
        ddphi        = c0;
        cosfac       = c1;

        rbp          = parm[1];
        rbpBA        = parmB[1]-parmA[1];
        ddphi       += rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[2];
        rbpBA        = parmB[2]-parmA[2];
        ddphi       += c2*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[3];
        rbpBA        = parmB[3]-parmA[3];
        ddphi       += c3*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[4];
        rbpBA        = parmB[4]-parmA[4];
        ddphi       += c4*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[5];
        rbpBA        = parmB[5]-parmA[5];
        ddphi       += c5*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;

        ddphi = -ddphi*sin_phi;         /*  11		*/

        do_dih_fup(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc, g, x, t1, t2, t3); /* 112		*/
        vtot += v;
    }
    *dvdlambda += dvdl_term;

    return vtot;
}

//! \endcond

/*! \brief Mysterious undocumented function */
static int
cmap_setup_grid_index(int ip, int grid_spacing, int *ipm1, int *ipp1, int *ipp2)
{
    int im1, ip1, ip2;

    if (ip < 0)
    {
        ip = ip + grid_spacing - 1;
    }
    else if (ip > grid_spacing)
    {
        ip = ip - grid_spacing - 1;
    }

    im1 = ip - 1;
    ip1 = ip + 1;
    ip2 = ip + 2;

    if (ip == 0)
    {
        im1 = grid_spacing - 1;
    }
    else if (ip == grid_spacing-2)
    {
        ip2 = 0;
    }
    else if (ip == grid_spacing-1)
    {
        ip1 = 0;
        ip2 = 1;
    }

    *ipm1 = im1;
    *ipp1 = ip1;
    *ipp2 = ip2;

    return ip;

}

real
cmap_dihs(int nbonds,
          const t_iatom forceatoms[], const t_iparams forceparams[],
          const gmx_cmap_t *cmap_grid,
          const rvec x[], rvec4 f[], rvec fshift[],
          const struct t_pbc *pbc, const struct t_graph *g,
          real gmx_unused lambda, real gmx_unused *dvdlambda,
          const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
          int  gmx_unused *global_atom_index)
{
    int         i, n;
    int         ai, aj, ak, al, am;
    int         a1i, a1j, a1k, a1l, a2i, a2j, a2k, a2l;
    int         type;
    int         t11, t21, t31, t12, t22, t32;
    int         iphi1, ip1m1, ip1p1, ip1p2;
    int         iphi2, ip2m1, ip2p1, ip2p2;
    int         l1, l2, l3;
    int         pos1, pos2, pos3, pos4;

    real        ty[4], ty1[4], ty2[4], ty12[4], tx[16];
    real        phi1, cos_phi1, sin_phi1, xphi1;
    real        phi2, cos_phi2, sin_phi2, xphi2;
    real        dx, tt, tu, e, df1, df2, vtot;
    real        ra21, rb21, rg21, rg1, rgr1, ra2r1, rb2r1, rabr1;
    real        ra22, rb22, rg22, rg2, rgr2, ra2r2, rb2r2, rabr2;
    real        fg1, hg1, fga1, hgb1, gaa1, gbb1;
    real        fg2, hg2, fga2, hgb2, gaa2, gbb2;
    real        fac;

    rvec        r1_ij, r1_kj, r1_kl, m1, n1;
    rvec        r2_ij, r2_kj, r2_kl, m2, n2;
    rvec        f1_i, f1_j, f1_k, f1_l;
    rvec        f2_i, f2_j, f2_k, f2_l;
    rvec        a1, b1, a2, b2;
    rvec        f1, g1, h1, f2, g2, h2;
    rvec        dtf1, dtg1, dth1, dtf2, dtg2, dth2;
    ivec        jt1, dt1_ij, dt1_kj, dt1_lj;
    ivec        jt2, dt2_ij, dt2_kj, dt2_lj;

    int         loop_index[4][4] = {
        {0, 4, 8, 12},
        {1, 5, 9, 13},
        {2, 6, 10, 14},
        {3, 7, 11, 15}
    };

    /* Total CMAP energy */
    vtot = 0;

    for (n = 0; n < nbonds; )
    {
        /* Five atoms are involved in the two torsions */
        type   = forceatoms[n++];
        ai     = forceatoms[n++];
        aj     = forceatoms[n++];
        ak     = forceatoms[n++];
        al     = forceatoms[n++];
        am     = forceatoms[n++];

        /* Which CMAP type is this */
        const int   cmapA = forceparams[type].cmap.cmapA;
        const real *cmapd = cmap_grid->cmapdata[cmapA].cmap.data();

        /* First torsion */
        a1i   = ai;
        a1j   = aj;
        a1k   = ak;
        a1l   = al;

        phi1  = dih_angle(x[a1i], x[a1j], x[a1k], x[a1l], pbc, r1_ij, r1_kj, r1_kl, m1, n1,
                          &t11, &t21, &t31);  /* 84 */

        cos_phi1 = std::cos(phi1);

        a1[0] = r1_ij[1]*r1_kj[2]-r1_ij[2]*r1_kj[1];
        a1[1] = r1_ij[2]*r1_kj[0]-r1_ij[0]*r1_kj[2];
        a1[2] = r1_ij[0]*r1_kj[1]-r1_ij[1]*r1_kj[0]; /* 9 */

        b1[0] = r1_kl[1]*r1_kj[2]-r1_kl[2]*r1_kj[1];
        b1[1] = r1_kl[2]*r1_kj[0]-r1_kl[0]*r1_kj[2];
        b1[2] = r1_kl[0]*r1_kj[1]-r1_kl[1]*r1_kj[0]; /* 9 */

        pbc_rvec_sub(pbc, x[a1l], x[a1k], h1);

        ra21  = iprod(a1, a1);       /* 5 */
        rb21  = iprod(b1, b1);       /* 5 */
        rg21  = iprod(r1_kj, r1_kj); /* 5 */
        rg1   = sqrt(rg21);

        rgr1  = 1.0/rg1;
        ra2r1 = 1.0/ra21;
        rb2r1 = 1.0/rb21;
        rabr1 = sqrt(ra2r1*rb2r1);

        sin_phi1 = rg1 * rabr1 * iprod(a1, h1) * (-1);

        if (cos_phi1 < -0.5 || cos_phi1 > 0.5)
        {
            phi1 = std::asin(sin_phi1);

            if (cos_phi1 < 0)
            {
                if (phi1 > 0)
                {
                    phi1 = M_PI - phi1;
                }
                else
                {
                    phi1 = -M_PI - phi1;
                }
            }
        }
        else
        {
            phi1 = std::acos(cos_phi1);

            if (sin_phi1 < 0)
            {
                phi1 = -phi1;
            }
        }

        xphi1 = phi1 + M_PI; /* 1 */

        /* Second torsion */
        a2i   = aj;
        a2j   = ak;
        a2k   = al;
        a2l   = am;

        phi2  = dih_angle(x[a2i], x[a2j], x[a2k], x[a2l], pbc, r2_ij, r2_kj, r2_kl, m2, n2,
                          &t12, &t22, &t32); /* 84 */

        cos_phi2 = std::cos(phi2);

        a2[0] = r2_ij[1]*r2_kj[2]-r2_ij[2]*r2_kj[1];
        a2[1] = r2_ij[2]*r2_kj[0]-r2_ij[0]*r2_kj[2];
        a2[2] = r2_ij[0]*r2_kj[1]-r2_ij[1]*r2_kj[0]; /* 9 */

        b2[0] = r2_kl[1]*r2_kj[2]-r2_kl[2]*r2_kj[1];
        b2[1] = r2_kl[2]*r2_kj[0]-r2_kl[0]*r2_kj[2];
        b2[2] = r2_kl[0]*r2_kj[1]-r2_kl[1]*r2_kj[0]; /* 9 */

        pbc_rvec_sub(pbc, x[a2l], x[a2k], h2);

        ra22  = iprod(a2, a2);         /* 5 */
        rb22  = iprod(b2, b2);         /* 5 */
        rg22  = iprod(r2_kj, r2_kj);   /* 5 */
        rg2   = sqrt(rg22);

        rgr2  = 1.0/rg2;
        ra2r2 = 1.0/ra22;
        rb2r2 = 1.0/rb22;
        rabr2 = sqrt(ra2r2*rb2r2);

        sin_phi2 = rg2 * rabr2 * iprod(a2, h2) * (-1);

        if (cos_phi2 < -0.5 || cos_phi2 > 0.5)
        {
            phi2 = std::asin(sin_phi2);

            if (cos_phi2 < 0)
            {
                if (phi2 > 0)
                {
                    phi2 = M_PI - phi2;
                }
                else
                {
                    phi2 = -M_PI - phi2;
                }
            }
        }
        else
        {
            phi2 = std::acos(cos_phi2);

            if (sin_phi2 < 0)
            {
                phi2 = -phi2;
            }
        }

        xphi2 = phi2 + M_PI; /* 1 */

        /* Range mangling */
        if (xphi1 < 0)
        {
            xphi1 = xphi1 + 2*M_PI;
        }
        else if (xphi1 >= 2*M_PI)
        {
            xphi1 = xphi1 - 2*M_PI;
        }

        if (xphi2 < 0)
        {
            xphi2 = xphi2 + 2*M_PI;
        }
        else if (xphi2 >= 2*M_PI)
        {
            xphi2 = xphi2 - 2*M_PI;
        }

        /* Number of grid points */
        dx = 2*M_PI / cmap_grid->grid_spacing;

        /* Where on the grid are we */
        iphi1 = static_cast<int>(xphi1/dx);
        iphi2 = static_cast<int>(xphi2/dx);

        iphi1 = cmap_setup_grid_index(iphi1, cmap_grid->grid_spacing, &ip1m1, &ip1p1, &ip1p2);
        iphi2 = cmap_setup_grid_index(iphi2, cmap_grid->grid_spacing, &ip2m1, &ip2p1, &ip2p2);

        pos1    = iphi1*cmap_grid->grid_spacing+iphi2;
        pos2    = ip1p1*cmap_grid->grid_spacing+iphi2;
        pos3    = ip1p1*cmap_grid->grid_spacing+ip2p1;
        pos4    = iphi1*cmap_grid->grid_spacing+ip2p1;

        ty[0]   = cmapd[pos1*4];
        ty[1]   = cmapd[pos2*4];
        ty[2]   = cmapd[pos3*4];
        ty[3]   = cmapd[pos4*4];

        ty1[0]   = cmapd[pos1*4+1];
        ty1[1]   = cmapd[pos2*4+1];
        ty1[2]   = cmapd[pos3*4+1];
        ty1[3]   = cmapd[pos4*4+1];

        ty2[0]   = cmapd[pos1*4+2];
        ty2[1]   = cmapd[pos2*4+2];
        ty2[2]   = cmapd[pos3*4+2];
        ty2[3]   = cmapd[pos4*4+2];

        ty12[0]   = cmapd[pos1*4+3];
        ty12[1]   = cmapd[pos2*4+3];
        ty12[2]   = cmapd[pos3*4+3];
        ty12[3]   = cmapd[pos4*4+3];

        /* Switch to degrees */
        dx    = 360.0 / cmap_grid->grid_spacing;
        xphi1 = xphi1 * RAD2DEG;
        xphi2 = xphi2 * RAD2DEG;

        for (i = 0; i < 4; i++) /* 16 */
        {
            tx[i]    = ty[i];
            tx[i+4]  = ty1[i]*dx;
            tx[i+8]  = ty2[i]*dx;
            tx[i+12] = ty12[i]*dx*dx;
        }

        real tc[16] = {0};
        for (int idx = 0; idx < 16; idx++) /* 1056 */
        {
            for (int k = 0; k < 16; k++)
            {
                tc[idx] += cmap_coeff_matrix[k*16+idx]*tx[k];
            }
        }

        tt    = (xphi1-iphi1*dx)/dx;
        tu    = (xphi2-iphi2*dx)/dx;

        e     = 0;
        df1   = 0;
        df2   = 0;

        for (i = 3; i >= 0; i--)
        {
            l1 = loop_index[i][3];
            l2 = loop_index[i][2];
            l3 = loop_index[i][1];

            e     = tt * e    + ((tc[i*4+3]*tu+tc[i*4+2])*tu + tc[i*4+1])*tu+tc[i*4];
            df1   = tu * df1  + (3.0*tc[l1]*tt+2.0*tc[l2])*tt+tc[l3];
            df2   = tt * df2  + (3.0*tc[i*4+3]*tu+2.0*tc[i*4+2])*tu+tc[i*4+1];
        }

        fac     = RAD2DEG/dx;
        df1     = df1   * fac;
        df2     = df2   * fac;

        /* CMAP energy */
        vtot += e;

        /* Do forces - first torsion */
        fg1       = iprod(r1_ij, r1_kj);
        hg1       = iprod(r1_kl, r1_kj);
        fga1      = fg1*ra2r1*rgr1;
        hgb1      = hg1*rb2r1*rgr1;
        gaa1      = -ra2r1*rg1;
        gbb1      = rb2r1*rg1;

        for (i = 0; i < DIM; i++)
        {
            dtf1[i]   = gaa1 * a1[i];
            dtg1[i]   = fga1 * a1[i] - hgb1 * b1[i];
            dth1[i]   = gbb1 * b1[i];

            f1[i]     = df1  * dtf1[i];
            g1[i]     = df1  * dtg1[i];
            h1[i]     = df1  * dth1[i];

            f1_i[i]   =  f1[i];
            f1_j[i]   = -f1[i] - g1[i];
            f1_k[i]   =  h1[i] + g1[i];
            f1_l[i]   = -h1[i];

            f[a1i][i] = f[a1i][i] + f1_i[i];
            f[a1j][i] = f[a1j][i] + f1_j[i]; /* - f1[i] - g1[i] */
            f[a1k][i] = f[a1k][i] + f1_k[i]; /* h1[i] + g1[i] */
            f[a1l][i] = f[a1l][i] + f1_l[i]; /* h1[i] */
        }

        /* Do forces - second torsion */
        fg2       = iprod(r2_ij, r2_kj);
        hg2       = iprod(r2_kl, r2_kj);
        fga2      = fg2*ra2r2*rgr2;
        hgb2      = hg2*rb2r2*rgr2;
        gaa2      = -ra2r2*rg2;
        gbb2      = rb2r2*rg2;

        for (i = 0; i < DIM; i++)
        {
            dtf2[i]   = gaa2 * a2[i];
            dtg2[i]   = fga2 * a2[i] - hgb2 * b2[i];
            dth2[i]   = gbb2 * b2[i];

            f2[i]     = df2  * dtf2[i];
            g2[i]     = df2  * dtg2[i];
            h2[i]     = df2  * dth2[i];

            f2_i[i]   =  f2[i];
            f2_j[i]   = -f2[i] - g2[i];
            f2_k[i]   =  h2[i] + g2[i];
            f2_l[i]   = -h2[i];

            f[a2i][i] = f[a2i][i] + f2_i[i]; /* f2[i] */
            f[a2j][i] = f[a2j][i] + f2_j[i]; /* - f2[i] - g2[i] */
            f[a2k][i] = f[a2k][i] + f2_k[i]; /* h2[i] + g2[i] */
            f[a2l][i] = f[a2l][i] + f2_l[i]; /* - h2[i] */
        }

        /* Shift forces */
        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, a1j), jt1);
            ivec_sub(SHIFT_IVEC(g, a1i),  jt1, dt1_ij);
            ivec_sub(SHIFT_IVEC(g, a1k),  jt1, dt1_kj);
            ivec_sub(SHIFT_IVEC(g, a1l),  jt1, dt1_lj);
            t11 = IVEC2IS(dt1_ij);
            t21 = IVEC2IS(dt1_kj);
            t31 = IVEC2IS(dt1_lj);

            copy_ivec(SHIFT_IVEC(g, a2j), jt2);
            ivec_sub(SHIFT_IVEC(g, a2i),  jt2, dt2_ij);
            ivec_sub(SHIFT_IVEC(g, a2k),  jt2, dt2_kj);
            ivec_sub(SHIFT_IVEC(g, a2l),  jt2, dt2_lj);
            t12 = IVEC2IS(dt2_ij);
            t22 = IVEC2IS(dt2_kj);
            t32 = IVEC2IS(dt2_lj);
        }
        else if (pbc)
        {
            t31 = pbc_rvec_sub(pbc, x[a1l], x[a1j], h1);
            t32 = pbc_rvec_sub(pbc, x[a2l], x[a2j], h2);
        }
        else
        {
            t31 = CENTRAL;
            t32 = CENTRAL;
        }

        rvec_inc(fshift[t11], f1_i);
        rvec_inc(fshift[CENTRAL], f1_j);
        rvec_inc(fshift[t21], f1_k);
        rvec_inc(fshift[t31], f1_l);

        rvec_inc(fshift[t21], f2_i);
        rvec_inc(fshift[CENTRAL], f2_j);
        rvec_inc(fshift[t22], f2_k);
        rvec_inc(fshift[t32], f2_l);
    }
    return vtot;
}


//! \cond
/***********************************************************
 *
 *   G R O M O S  9 6   F U N C T I O N S
 *
 ***********************************************************/
static real g96harmonic(real kA, real kB, real xA, real xB, real x, real lambda,
                        real *V, real *F)
{
    const real half = 0.5;
    real       L1, kk, x0, dx, dx2;
    real       v, f, dvdlambda;

    L1    = 1.0-lambda;
    kk    = L1*kA+lambda*kB;
    x0    = L1*xA+lambda*xB;

    dx    = x-x0;
    dx2   = dx*dx;

    f          = -kk*dx;
    v          = half*kk*dx2;
    dvdlambda  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;

    *F    = f;
    *V    = v;

    return dvdlambda;

    /* That was 21 flops */
}

real g96bonds(int nbonds,
              const t_iatom forceatoms[], const t_iparams forceparams[],
              const rvec x[], rvec4 f[], rvec fshift[],
              const t_pbc *pbc, const t_graph *g,
              real lambda, real *dvdlambda,
              const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
              int gmx_unused *global_atom_index)
{
    int  i, m, ki, ai, aj, type;
    real dr2, fbond, vbond, fij, vtot;
    rvec dx;
    ivec dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx); /*   3      */
        dr2  = iprod(dx, dx);                       /*   5		*/

        *dvdlambda += g96harmonic(forceparams[type].harmonic.krA,
                                  forceparams[type].harmonic.krB,
                                  forceparams[type].harmonic.rA,
                                  forceparams[type].harmonic.rB,
                                  dr2, lambda, &vbond, &fbond);

        vtot  += 0.5*vbond;                         /* 1*/

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)     /*  15		*/
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }               /* 44 TOTAL	*/
    return vtot;
}

static real g96bond_angle(const rvec xi, const rvec xj, const rvec xk, const t_pbc *pbc,
                          rvec r_ij, rvec r_kj,
                          int *t1, int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
    real costh;

    *t1 = pbc_rvec_sub(pbc, xi, xj, r_ij); /*  3		*/
    *t2 = pbc_rvec_sub(pbc, xk, xj, r_kj); /*  3		*/

    costh = cos_angle(r_ij, r_kj);         /* 25		*/
    /* 41 TOTAL	*/
    return costh;
}

real g96angles(int nbonds,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec4 f[], rvec fshift[],
               const t_pbc *pbc, const t_graph *g,
               real lambda, real *dvdlambda,
               const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
               int gmx_unused *global_atom_index)
{
    int  i, ai, aj, ak, type, m, t1, t2;
    rvec r_ij, r_kj;
    real cos_theta, dVdt, va, vtot;
    real rij_1, rij_2, rkj_1, rkj_2, rijrkj_1;
    rvec f_i, f_j, f_k;
    ivec jt, dt_ij, dt_kj;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];

        cos_theta  = g96bond_angle(x[ai], x[aj], x[ak], pbc, r_ij, r_kj, &t1, &t2);

        *dvdlambda += g96harmonic(forceparams[type].harmonic.krA,
                                  forceparams[type].harmonic.krB,
                                  forceparams[type].harmonic.rA,
                                  forceparams[type].harmonic.rB,
                                  cos_theta, lambda, &va, &dVdt);
        vtot    += va;

        rij_1    = gmx::invsqrt(iprod(r_ij, r_ij));
        rkj_1    = gmx::invsqrt(iprod(r_kj, r_kj));
        rij_2    = rij_1*rij_1;
        rkj_2    = rkj_1*rkj_1;
        rijrkj_1 = rij_1*rkj_1;         /* 23 */

        for (m = 0; (m < DIM); m++)     /*  42	*/
        {
            f_i[m]    = dVdt*(r_kj[m]*rijrkj_1 - r_ij[m]*rij_2*cos_theta);
            f_k[m]    = dVdt*(r_ij[m]*rijrkj_1 - r_kj[m]*rkj_2*cos_theta);
            f_j[m]    = -f_i[m]-f_k[m];
            f[ai][m] += f_i[m];
            f[aj][m] += f_j[m];
            f[ak][m] += f_k[m];
        }

        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);

            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
        }
        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);          /* 9 */
        /* 163 TOTAL	*/
    }
    return vtot;
}

real cross_bond_bond(int nbonds,
                     const t_iatom forceatoms[], const t_iparams forceparams[],
                     const rvec x[], rvec4 f[], rvec fshift[],
                     const t_pbc *pbc, const t_graph *g,
                     real gmx_unused lambda, real gmx_unused *dvdlambda,
                     const t_mdatoms gmx_unused *md, t_fcdata gmx_unused  *fcd,
                     int gmx_unused *global_atom_index)
{
    /* Potential from Lawrence and Skimmer, Chem. Phys. Lett. 372 (2003)
     * pp. 842-847
     */
    int  i, ai, aj, ak, type, m, t1, t2;
    rvec r_ij, r_kj;
    real vtot, vrr, s1, s2, r1, r2, r1e, r2e, krr;
    rvec f_i, f_j, f_k;
    ivec jt, dt_ij, dt_kj;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        r1e  = forceparams[type].cross_bb.r1e;
        r2e  = forceparams[type].cross_bb.r2e;
        krr  = forceparams[type].cross_bb.krr;

        /* Compute distance vectors ... */
        t1 = pbc_rvec_sub(pbc, x[ai], x[aj], r_ij);
        t2 = pbc_rvec_sub(pbc, x[ak], x[aj], r_kj);

        /* ... and their lengths */
        r1 = norm(r_ij);
        r2 = norm(r_kj);

        /* Deviations from ideality */
        s1 = r1-r1e;
        s2 = r2-r2e;

        /* Energy (can be negative!) */
        vrr   = krr*s1*s2;
        vtot += vrr;

        /* Forces */
        svmul(-krr*s2/r1, r_ij, f_i);
        svmul(-krr*s1/r2, r_kj, f_k);

        for (m = 0; (m < DIM); m++)     /*  12	*/
        {
            f_j[m]    = -f_i[m] - f_k[m];
            f[ai][m] += f_i[m];
            f[aj][m] += f_j[m];
            f[ak][m] += f_k[m];
        }

        /* Virial stuff */
        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);

            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
        }
        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);          /* 9 */
        /* 163 TOTAL	*/
    }
    return vtot;
}

real cross_bond_angle(int nbonds,
                      const t_iatom forceatoms[], const t_iparams forceparams[],
                      const rvec x[], rvec4 f[], rvec fshift[],
                      const t_pbc *pbc, const t_graph *g,
                      real gmx_unused lambda, real gmx_unused *dvdlambda,
                      const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                      int gmx_unused *global_atom_index)
{
    /* Potential from Lawrence and Skimmer, Chem. Phys. Lett. 372 (2003)
     * pp. 842-847
     */
    int  i, ai, aj, ak, type, m, t1, t2;
    rvec r_ij, r_kj, r_ik;
    real vtot, vrt, s1, s2, s3, r1, r2, r3, r1e, r2e, r3e, krt, k1, k2, k3;
    rvec f_i, f_j, f_k;
    ivec jt, dt_ij, dt_kj;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        r1e  = forceparams[type].cross_ba.r1e;
        r2e  = forceparams[type].cross_ba.r2e;
        r3e  = forceparams[type].cross_ba.r3e;
        krt  = forceparams[type].cross_ba.krt;

        /* Compute distance vectors ... */
        t1 = pbc_rvec_sub(pbc, x[ai], x[aj], r_ij);
        t2 = pbc_rvec_sub(pbc, x[ak], x[aj], r_kj);
        pbc_rvec_sub(pbc, x[ai], x[ak], r_ik);

        /* ... and their lengths */
        r1 = norm(r_ij);
        r2 = norm(r_kj);
        r3 = norm(r_ik);

        /* Deviations from ideality */
        s1 = r1-r1e;
        s2 = r2-r2e;
        s3 = r3-r3e;

        /* Energy (can be negative!) */
        vrt   = krt*s3*(s1+s2);
        vtot += vrt;

        /* Forces */
        k1 = -krt*(s3/r1);
        k2 = -krt*(s3/r2);
        k3 = -krt*(s1+s2)/r3;
        for (m = 0; (m < DIM); m++)
        {
            f_i[m] = k1*r_ij[m] + k3*r_ik[m];
            f_k[m] = k2*r_kj[m] - k3*r_ik[m];
            f_j[m] = -f_i[m] - f_k[m];
        }

        for (m = 0; (m < DIM); m++)     /*  12	*/
        {
            f[ai][m] += f_i[m];
            f[aj][m] += f_j[m];
            f[ak][m] += f_k[m];
        }

        /* Virial stuff */
        if (g)
        {
            copy_ivec(SHIFT_IVEC(g, aj), jt);

            ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
            ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
            t1 = IVEC2IS(dt_ij);
            t2 = IVEC2IS(dt_kj);
        }
        rvec_inc(fshift[t1], f_i);
        rvec_inc(fshift[CENTRAL], f_j);
        rvec_inc(fshift[t2], f_k);          /* 9 */
        /* 163 TOTAL	*/
    }
    return vtot;
}

static real bonded_tab(const char *type, int table_nr,
                       const bondedtable_t *table, real kA, real kB, real r,
                       real lambda, real *V, real *F)
{
    real k, tabscale, *VFtab, rt, eps, eps2, Yt, Ft, Geps, Heps2, Fp, VV, FF;
    int  n0, nnn;
    real dvdlambda;

    k = (1.0 - lambda)*kA + lambda*kB;

    tabscale = table->scale;
    VFtab    = table->data;

    rt    = r*tabscale;
    n0    = static_cast<int>(rt);
    if (n0 >= table->n)
    {
        gmx_fatal(FARGS, "A tabulated %s interaction table number %d is out of the table range: r %f, between table indices %d and %d, table length %d",
                  type, table_nr, r, n0, n0+1, table->n);
    }
    eps   = rt - n0;
    eps2  = eps*eps;
    nnn   = 4*n0;
    Yt    = VFtab[nnn];
    Ft    = VFtab[nnn+1];
    Geps  = VFtab[nnn+2]*eps;
    Heps2 = VFtab[nnn+3]*eps2;
    Fp    = Ft + Geps + Heps2;
    VV    = Yt + Fp*eps;
    FF    = Fp + Geps + 2.0*Heps2;

    *F         = -k*FF*tabscale;
    *V         = k*VV;
    dvdlambda  = (kB - kA)*VV;

    return dvdlambda;

    /* That was 22 flops */
}

real tab_bonds(int nbonds,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec4 f[], rvec fshift[],
               const t_pbc *pbc, const t_graph *g,
               real lambda, real *dvdlambda,
               const t_mdatoms gmx_unused *md, t_fcdata *fcd,
               int gmx_unused  *global_atom_index)
{
    int  i, m, ki, ai, aj, type, table;
    real dr, dr2, fbond, vbond, fij, vtot;
    rvec dx;
    ivec dt;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];

        ki   = pbc_rvec_sub(pbc, x[ai], x[aj], dx); /*   3      */
        dr2  = iprod(dx, dx);                       /*   5		*/
        dr   = dr2*gmx::invsqrt(dr2);               /*  10		*/

        table = forceparams[type].tab.table;

        *dvdlambda += bonded_tab("bond", table,
                                 &fcd->bondtab[table],
                                 forceparams[type].tab.kA,
                                 forceparams[type].tab.kB,
                                 dr, lambda, &vbond, &fbond); /*  22 */

        if (dr2 == 0.0)
        {
            continue;
        }


        vtot  += vbond;             /* 1*/
        fbond *= gmx::invsqrt(dr2); /*   6		*/
        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }
        for (m = 0; (m < DIM); m++)     /*  15		*/
        {
            fij                 = fbond*dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }
    }               /* 62 TOTAL	*/
    return vtot;
}

real tab_angles(int nbonds,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const rvec x[], rvec4 f[], rvec fshift[],
                const t_pbc *pbc, const t_graph *g,
                real lambda, real *dvdlambda,
                const t_mdatoms gmx_unused  *md, t_fcdata *fcd,
                int gmx_unused *global_atom_index)
{
    int  i, ai, aj, ak, t1, t2, type, table;
    rvec r_ij, r_kj;
    real cos_theta, cos_theta2, theta, dVdt, va, vtot;
    ivec jt, dt_ij, dt_kj;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];

        theta  = bond_angle(x[ai], x[aj], x[ak], pbc,
                            r_ij, r_kj, &cos_theta, &t1, &t2); /*  41		*/

        table = forceparams[type].tab.table;

        *dvdlambda += bonded_tab("angle", table,
                                 &fcd->angletab[table],
                                 forceparams[type].tab.kA,
                                 forceparams[type].tab.kB,
                                 theta, lambda, &va, &dVdt); /*  22  */
        vtot += va;

        cos_theta2 = gmx::square(cos_theta);            /*   1		*/
        if (cos_theta2 < 1)
        {
            int  m;
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            rvec f_i, f_j, f_k;

            st    = dVdt*gmx::invsqrt(1 - cos_theta2); /*  12		*/
            sth   = st*cos_theta;                      /*   1		*/
            nrkj2 = iprod(r_kj, r_kj);                 /*   5		*/
            nrij2 = iprod(r_ij, r_ij);

            cik = st*gmx::invsqrt(nrkj2*nrij2); /*  12		*/
            cii = sth/nrij2;                    /*  10		*/
            ckk = sth/nrkj2;                    /*  10		*/

            for (m = 0; (m < DIM); m++)         /*  39		*/
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                f[ai][m] += f_i[m];
                f[aj][m] += f_j[m];
                f[ak][m] += f_k[m];
            }
            if (g)
            {
                copy_ivec(SHIFT_IVEC(g, aj), jt);

                ivec_sub(SHIFT_IVEC(g, ai), jt, dt_ij);
                ivec_sub(SHIFT_IVEC(g, ak), jt, dt_kj);
                t1 = IVEC2IS(dt_ij);
                t2 = IVEC2IS(dt_kj);
            }
            rvec_inc(fshift[t1], f_i);
            rvec_inc(fshift[CENTRAL], f_j);
            rvec_inc(fshift[t2], f_k);
        }                                       /* 169 TOTAL	*/
    }
    return vtot;
}

real tab_dihs(int nbonds,
              const t_iatom forceatoms[], const t_iparams forceparams[],
              const rvec x[], rvec4 f[], rvec fshift[],
              const t_pbc *pbc, const t_graph *g,
              real lambda, real *dvdlambda,
              const t_mdatoms gmx_unused *md, t_fcdata *fcd,
              int gmx_unused *global_atom_index)
{
    int  i, type, ai, aj, ak, al, table;
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi, ddphi, vpd, vtot;

    vtot = 0.0;
    for (i = 0; (i < nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];

        phi = dih_angle(x[ai], x[aj], x[ak], x[al], pbc, r_ij, r_kj, r_kl, m, n,
                        &t1, &t2, &t3);  /*  84  */

        table = forceparams[type].tab.table;

        /* Hopefully phi+M_PI never results in values < 0 */
        *dvdlambda += bonded_tab("dihedral", table,
                                 &fcd->dihtab[table],
                                 forceparams[type].tab.kA,
                                 forceparams[type].tab.kB,
                                 phi+M_PI, lambda, &vpd, &ddphi);

        vtot += vpd;
        do_dih_fup(ai, aj, ak, al, -ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc, g, x, t1, t2, t3); /* 112	*/

    }                                                 /* 227 TOTAL  */

    return vtot;
}

//! \endcond
