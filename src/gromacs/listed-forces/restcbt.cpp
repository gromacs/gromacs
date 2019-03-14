/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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
 * \brief
 * This file contains function definitions necessary
 * for computations of forces due to restricted angle, restricted dihedral and
 * combined bending-torsion potentials.
 *
 * \author Nicolae Goga
 *
 * \ingroup module_listed-forces
 */
#include "gmxpre.h"

#include "restcbt.h"

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/idef.h"

/* This function computes factors needed for restricted angle potential.
 * For explanations on formula used see file "restcbt.h" */

void compute_factors_restangles(int type, const t_iparams forceparams[],
                                rvec delta_ante,  rvec delta_post,
                                real *prefactor, real *ratio_ante, real *ratio_post, real *v)
{
    real theta_equil, k_bending;
    real cosine_theta_equil;
    real c_ante, c_cros, c_post;
    real norm;
    real delta_cosine, cosine_theta;
    real sine_theta_sq;
    real term_theta_theta_equil;

    k_bending          = forceparams[type].harmonic.krA;
    theta_equil        =  forceparams[type].harmonic.rA*DEG2RAD;
    theta_equil        = M_PI - theta_equil;
    cosine_theta_equil = cos(theta_equil);

    c_ante = iprod(delta_ante, delta_ante);
    c_cros = iprod(delta_ante, delta_post);
    c_post = iprod(delta_post, delta_post);

    norm          = gmx::invsqrt(c_ante * c_post);
    cosine_theta  = c_cros * norm;
    sine_theta_sq = 1 - cosine_theta * cosine_theta;

    *ratio_ante = c_cros / c_ante;
    *ratio_post = c_cros / c_post;

    delta_cosine           = cosine_theta - cosine_theta_equil;
    term_theta_theta_equil = 1 - cosine_theta * cosine_theta_equil;
    *prefactor             = -(k_bending) * delta_cosine * norm * term_theta_theta_equil / (sine_theta_sq * sine_theta_sq);

    *v = k_bending * 0.5 * delta_cosine * delta_cosine / sine_theta_sq;

}


/* Compute factors for restricted dihedral potential
 * For explanations on formula used see file "restcbt.h" */
void compute_factors_restrdihs(int type,  const t_iparams forceparams[],
                               rvec delta_ante, rvec delta_crnt, rvec delta_post,
                               real *factor_phi_ai_ante, real *factor_phi_ai_crnt, real *factor_phi_ai_post,
                               real *factor_phi_aj_ante, real *factor_phi_aj_crnt, real *factor_phi_aj_post,
                               real *factor_phi_ak_ante, real *factor_phi_ak_crnt, real *factor_phi_ak_post,
                               real *factor_phi_al_ante, real *factor_phi_al_crnt, real *factor_phi_al_post,
                               real *prefactor_phi, real *v)
{

    real phi0, cosine_phi0;
    real k_torsion;
    real c_self_ante, c_self_crnt, c_self_post;
    real c_cros_ante, c_cros_acrs, c_cros_post;
    real c_prod, d_post, d_ante;
    real sine_phi_sq, cosine_phi;
    real delta_cosine, term_phi_phi0;
    real ratio_phi_ante, ratio_phi_post;
    real norm_phi;

    /* Read parameters phi0 and k_torsion */
    phi0        = forceparams[type].pdihs.phiA * DEG2RAD;
    cosine_phi0 = cos(phi0);
    k_torsion   = forceparams[type].pdihs.cpA;

    /* Computation of the cosine of the dihedral angle. The scalar ("dot") product  method
     * is used. c_*_* cummulate the scalar products of the differences of particles
     * positions while c_prod, d_ante and d_post are differences of products of scalar
     * terms that are parts of the derivatives of forces */
    c_self_ante = iprod(delta_ante, delta_ante);
    c_self_crnt = iprod(delta_crnt, delta_crnt);
    c_self_post = iprod(delta_post, delta_post);
    c_cros_ante = iprod(delta_ante, delta_crnt);
    c_cros_acrs = iprod(delta_ante, delta_post);
    c_cros_post = iprod(delta_crnt, delta_post);
    c_prod      = c_cros_ante * c_cros_post - c_self_crnt * c_cros_acrs;
    d_ante      = c_self_ante * c_self_crnt - c_cros_ante * c_cros_ante;
    d_post      = c_self_post * c_self_crnt - c_cros_post * c_cros_post;

    /*  When three consecutive beads align, we obtain values close to zero.
     *	Here we avoid small values to prevent round-off errors. */
    if (d_ante < GMX_REAL_EPS)
    {
        d_ante = GMX_REAL_EPS;
    }
    if (d_post < GMX_REAL_EPS)
    {
        d_post = GMX_REAL_EPS;
    }

    /* Computes the square of the sinus of phi  in sine_phi_sq */
    norm_phi    = gmx::invsqrt(d_ante * d_post);
    cosine_phi  = c_prod * norm_phi;
    sine_phi_sq = 1.0 - cosine_phi * cosine_phi;

    /*	It is possible that cosine_phi is slightly bigger than 1.0 due to round-off errors. */
    if (sine_phi_sq < 0.0)
    {
        sine_phi_sq = 0.0;
    }

    /* Computation of the differences of cosines (delta_cosine) and a term (term_phi_phi0)
     * that is part of the common prefactor_phi */

    delta_cosine  = cosine_phi - cosine_phi0;
    term_phi_phi0 = 1 - cosine_phi * cosine_phi0;


    /*      Computation of ratios */
    ratio_phi_ante = c_prod / d_ante;
    ratio_phi_post = c_prod / d_post;

    /*      Computation of the prefactor - common term for all forces */
    *prefactor_phi = -(k_torsion) * delta_cosine * norm_phi * term_phi_phi0 / (sine_phi_sq * sine_phi_sq);

    /* Computation of force factors.  Factors factor_phi_*  are coming from the
     * derivatives of the torsion angle (phi) with respect to the beads ai, aj, al, ak,
     * (four) coordinates and they are multiplied in the force computations with the
     * differences of the particles positions stored in parameters delta_ante,
     * delta_crnt, delta_post. For formulas see file "restcbt.h" */

    *factor_phi_ai_ante = ratio_phi_ante * c_self_crnt;
    *factor_phi_ai_crnt = -c_cros_post - ratio_phi_ante * c_cros_ante;
    *factor_phi_ai_post = c_self_crnt;
    *factor_phi_aj_ante = -c_cros_post - ratio_phi_ante * (c_self_crnt + c_cros_ante);
    *factor_phi_aj_crnt = c_cros_post + c_cros_acrs * 2.0 + ratio_phi_ante * (c_self_ante + c_cros_ante) + ratio_phi_post * c_self_post;
    *factor_phi_aj_post = -(c_cros_ante + c_self_crnt) - ratio_phi_post * c_cros_post;
    *factor_phi_ak_ante = c_cros_post + c_self_crnt + ratio_phi_ante * c_cros_ante;
    *factor_phi_ak_crnt = -(c_cros_ante + c_cros_acrs * 2.0)- ratio_phi_ante * c_self_ante - ratio_phi_post * (c_self_post + c_cros_post);
    *factor_phi_ak_post = c_cros_ante + ratio_phi_post * (c_self_crnt + c_cros_post);
    *factor_phi_al_ante = -c_self_crnt;
    *factor_phi_al_crnt = c_cros_ante + ratio_phi_post * c_cros_post;
    *factor_phi_al_post = -ratio_phi_post * c_self_crnt;

    /* Contribution to energy  - see formula in file "restcbt.h"*/
    *v = k_torsion * 0.5 * delta_cosine * delta_cosine / sine_phi_sq;
}



/* Compute factors for CBT potential
 * For explanations on formula used see file "restcbt.h" */

void compute_factors_cbtdihs(int type,  const t_iparams forceparams[],
                             rvec delta_ante, rvec delta_crnt, rvec delta_post,
                             rvec f_phi_ai, rvec f_phi_aj, rvec f_phi_ak, rvec f_phi_al,
                             rvec f_theta_ante_ai, rvec f_theta_ante_aj, rvec f_theta_ante_ak,
                             rvec f_theta_post_aj, rvec f_theta_post_ak, rvec f_theta_post_al,
                             real * v)
{
    int  j, d;
    real torsion_coef[NR_CBTDIHS];
    real c_self_ante, c_self_crnt, c_self_post;
    real c_cros_ante, c_cros_acrs, c_cros_post;
    real c_prod, d_ante,  d_post;
    real norm_phi, norm_theta_ante, norm_theta_post;
    real cosine_phi, cosine_theta_ante, cosine_theta_post;
    real sine_theta_ante_sq, sine_theta_post_sq;
    real sine_theta_ante, sine_theta_post;
    real prefactor_phi;
    real ratio_phi_ante, ratio_phi_post;
    real r1, r2;
    real factor_phi_ai_ante, factor_phi_ai_crnt, factor_phi_ai_post;
    real factor_phi_aj_ante, factor_phi_aj_crnt, factor_phi_aj_post;
    real factor_phi_ak_ante, factor_phi_ak_crnt, factor_phi_ak_post;
    real factor_phi_al_ante, factor_phi_al_crnt, factor_phi_al_post;
    real prefactor_theta_ante, ratio_theta_ante_ante, ratio_theta_ante_crnt;
    real prefactor_theta_post, ratio_theta_post_crnt, ratio_theta_post_post;

    /* The formula for combined bending-torsion potential (see file "restcbt.h") contains
     * in its expression not only the dihedral angle \f[\phi\f] but also \f[\theta_{i-1}\f]
     * (theta_ante bellow) and \f[\theta_{i}\f] (theta_post bellow)--- the adjacent bending
     * angles. The forces for the particles ai, aj, ak, al have components coming from the
     * derivatives of the potential with respect to all three angles.
     * This function is organised in 4 parts
     * PART 1 - Computes force factors common to all the derivatives for the four particles
     * PART 2 - Computes the force components due to the derivatives of dihedral angle Phi
     * PART 3 - Computes the force components due to the derivatives of bending angle Theta_Ante
     * PART 4 - Computes the force components due to the derivatives of bending angle Theta_Post
     * Bellow we will respct thuis structure */


    /* PART 1 - COMPUTES FORCE FACTORS COMMON TO ALL DERIVATIVES FOR THE FOUR PARTICLES */


    for (j = 0; (j < NR_CBTDIHS); j++)
    {
        torsion_coef[j]  = forceparams[type].cbtdihs.cbtcA[j];
    }

    /* Computation of the cosine of the dihedral angle. The scalar ("dot") product  method
     * is used. c_*_* cummulate the scalar products of the differences of particles
     * positions while c_prod, d_ante and d_post are differences of products of scalar
     * terms that are parts of the derivatives of forces */

    c_self_ante = iprod(delta_ante, delta_ante);
    c_self_crnt = iprod(delta_crnt, delta_crnt);
    c_self_post = iprod(delta_post, delta_post);
    c_cros_ante = iprod(delta_ante, delta_crnt);
    c_cros_acrs = iprod(delta_ante, delta_post);
    c_cros_post = iprod(delta_crnt, delta_post);
    c_prod      = c_cros_ante * c_cros_post - c_self_crnt * c_cros_acrs;
    d_ante      = c_self_ante * c_self_crnt - c_cros_ante * c_cros_ante;
    d_post      = c_self_post * c_self_crnt - c_cros_post * c_cros_post;

    /*  When three consecutive beads align, we obtain values close to zero.
       Here we avoid small values to prevent round-off errors. */
    if (d_ante < GMX_REAL_EPS)
    {
        d_ante = GMX_REAL_EPS;
    }
    if (d_post < GMX_REAL_EPS)
    {
        d_post = GMX_REAL_EPS;
    }

    /* Computations of cosines */
    norm_phi           = gmx::invsqrt(d_ante * d_post);
    norm_theta_ante    = gmx::invsqrt(c_self_ante * c_self_crnt);
    norm_theta_post    = gmx::invsqrt(c_self_crnt * c_self_post);
    cosine_phi         = c_prod * norm_phi;
    cosine_theta_ante  = c_cros_ante * norm_theta_ante;
    cosine_theta_post  = c_cros_post * norm_theta_post;
    sine_theta_ante_sq = 1 - cosine_theta_ante * cosine_theta_ante;
    sine_theta_post_sq = 1 - cosine_theta_post * cosine_theta_post;

    /*	It is possible that cosine_theta is slightly bigger than 1.0 due to round-off errors. */
    if (sine_theta_ante_sq < 0.0)
    {
        sine_theta_ante_sq = 0.0;
    }
    if (sine_theta_post_sq < 0.0)
    {
        sine_theta_post_sq = 0.0;
    }

    sine_theta_ante = std::sqrt(sine_theta_ante_sq);
    sine_theta_post = std::sqrt(sine_theta_post_sq);

    /* PART 2 - COMPUTES FORCE COMPONENTS DUE TO DERIVATIVES TO DIHEDRAL ANGLE PHI */

    /*      Computation of ratios */
    ratio_phi_ante = c_prod / d_ante;
    ratio_phi_post = c_prod / d_post;

    /*       Computation of the prefactor */
    /*      Computing 2nd power */
    r1 = cosine_phi;

    prefactor_phi = -torsion_coef[0] * norm_phi * (torsion_coef[2] + torsion_coef[3] * 2.0 * cosine_phi + torsion_coef[4] * 3.0 * (r1 * r1) + 4*torsion_coef[5]*r1*r1*r1) *
        sine_theta_ante_sq * sine_theta_ante * sine_theta_post_sq * sine_theta_post;

    /* Computation of factors (important for gaining speed). Factors factor_phi_*  are coming from the
     * derivatives of the torsion angle (phi) with respect to the beads ai, aj, al, ak,
     * (four) coordinates and they are multiplied in the force computations with the
     * differences of the particles positions stored in parameters delta_ante,
     * delta_crnt, delta_post. For formulas see file "restcbt.h" */

    factor_phi_ai_ante = ratio_phi_ante * c_self_crnt;
    factor_phi_ai_crnt = -c_cros_post - ratio_phi_ante * c_cros_ante;
    factor_phi_ai_post = c_self_crnt;
    factor_phi_aj_ante = -c_cros_post - ratio_phi_ante * (c_self_crnt + c_cros_ante);
    factor_phi_aj_crnt = c_cros_post + c_cros_acrs * 2.0 + ratio_phi_ante * (c_self_ante + c_cros_ante) +  ratio_phi_post * c_self_post;
    factor_phi_aj_post = -(c_cros_ante + c_self_crnt) - ratio_phi_post * c_cros_post;
    factor_phi_ak_ante = c_cros_post + c_self_crnt + ratio_phi_ante * c_cros_ante;
    factor_phi_ak_crnt = -(c_cros_ante + c_cros_acrs * 2.0) - ratio_phi_ante * c_self_ante - ratio_phi_post * (c_self_post + c_cros_post);
    factor_phi_ak_post = c_cros_ante + ratio_phi_post * (c_self_crnt + c_cros_post);
    factor_phi_al_ante = -c_self_crnt;
    factor_phi_al_crnt = c_cros_ante + ratio_phi_post * c_cros_post;
    factor_phi_al_post = -ratio_phi_post * c_self_crnt;

    /* Computation of forces due to the derivatives of dihedral angle phi*/
    for (d = 0; d < DIM; d++)
    {
        f_phi_ai[d] = prefactor_phi * (factor_phi_ai_ante * delta_ante[d] + factor_phi_ai_crnt * delta_crnt[d] + factor_phi_ai_post * delta_post[d]);
        f_phi_aj[d] = prefactor_phi * (factor_phi_aj_ante * delta_ante[d] + factor_phi_aj_crnt * delta_crnt[d] + factor_phi_aj_post * delta_post[d]);
        f_phi_ak[d] = prefactor_phi * (factor_phi_ak_ante * delta_ante[d] + factor_phi_ak_crnt * delta_crnt[d] + factor_phi_ak_post * delta_post[d]);
        f_phi_al[d] = prefactor_phi * (factor_phi_al_ante * delta_ante[d] + factor_phi_al_crnt * delta_crnt[d] + factor_phi_al_post * delta_post[d]);
    }

    /* PART 3 - COMPUTES THE FORCE COMPONENTS DUE TO THE DERIVATIVES OF BENDING ANGLE THETHA_ANTHE */
    /*      Computation of ratios */
    ratio_theta_ante_ante = c_cros_ante / c_self_ante;
    ratio_theta_ante_crnt = c_cros_ante / c_self_crnt;

    /*      Computation of the prefactor */
    /*      Computing 2nd power */
    r1 = cosine_phi;
    /*      Computing 3rd power */
    r2 = cosine_phi;

    prefactor_theta_ante = -torsion_coef[0] * norm_theta_ante * ( torsion_coef[1] + torsion_coef[2] * cosine_phi + torsion_coef[3] * (r1 * r1) +
                                                                  torsion_coef[4] * (r2 * (r2 * r2))+ torsion_coef[5] * (r2 * (r2 * (r2 * r2)))) * (-3.0) * cosine_theta_ante * sine_theta_ante * sine_theta_post_sq * sine_theta_post;


    /*      Computation of forces due to the derivatives of bending angle theta_ante */
    for (d = 0; d < DIM; d++)
    {
        f_theta_ante_ai[d] = prefactor_theta_ante * (ratio_theta_ante_ante * delta_ante[d] - delta_crnt[d]);
        f_theta_ante_aj[d] = prefactor_theta_ante * ((ratio_theta_ante_crnt + 1.0) * delta_crnt[d] - (ratio_theta_ante_ante + 1.0) * delta_ante[d]);
        f_theta_ante_ak[d] = prefactor_theta_ante * (delta_ante[d] - ratio_theta_ante_crnt * delta_crnt[d]);
    }

    /* PART 4 - COMPUTES THE FORCE COMPONENTS DUE TO THE DERIVATIVES OF THE BENDING ANGLE THETA_POST */

    /*      Computation of ratios */
    ratio_theta_post_crnt = c_cros_post / c_self_crnt;
    ratio_theta_post_post = c_cros_post / c_self_post;

    /*     Computation of the prefactor */
    /*      Computing 2nd power */
    r1 = cosine_phi;
    /*      Computing 3rd power */
    r2 = cosine_phi;

    prefactor_theta_post = -torsion_coef[0] * norm_theta_post * (torsion_coef[1] + torsion_coef[2] * cosine_phi + torsion_coef[3] * (r1 * r1) +
                                                                 torsion_coef[4] * (r2 * (r2 * r2)) + torsion_coef[5] * (r2 * (r2 * (r2 * r2)))) * sine_theta_ante_sq * sine_theta_ante * (-3.0) * cosine_theta_post * sine_theta_post;


    /*      Computation of forces due to the derivatives of bending angle Theta_Post */
    for (d = 0; d < DIM; d++)
    {
        f_theta_post_aj[d] = prefactor_theta_post * (ratio_theta_post_crnt * delta_crnt[d] - delta_post[d]);
        f_theta_post_ak[d] = prefactor_theta_post * ((ratio_theta_post_post + 1.0) * delta_post[d] - (ratio_theta_post_crnt + 1.0) * delta_crnt[d]);
        f_theta_post_al[d] = prefactor_theta_post * (delta_crnt[d] - ratio_theta_post_post * delta_post[d]);
    }
    r1 = cosine_phi;
    r2 = cosine_phi;

    /* Contribution to energy - for formula see file "restcbt.h" */
    *v = torsion_coef[0] * (torsion_coef[1] + torsion_coef[2] * cosine_phi + torsion_coef[3] * (r1 * r1) +
                            torsion_coef[4] * (r2 * (r2 * r2)) + torsion_coef[5] * (r2 * (r2 * (r2 * r2)))) * sine_theta_ante_sq *
        sine_theta_ante * sine_theta_post_sq * sine_theta_post;


}
