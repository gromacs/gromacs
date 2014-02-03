/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#include <math.h>
#include <assert.h>
#include "physics.h"
#include "vec.h"
#include "gromacs/math/utilities.h"
#include "txtdump.h"
#include "bondf.h"
#include "smalloc.h"
#include "pbc.h"
#include "ns.h"
#include "macros.h"
#include "names.h"
#include "gmx_fatal.h"
#include "mshift.h"
#include "main.h"
#include "disre.h"
#include "orires.h"
#include "force.h"
#include "nonbonded.h"

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

    norm          = gmx_invsqrt(c_ante * c_post);
    cosine_theta  = c_cros * norm;
    sine_theta_sq = 1 - cosine_theta * cosine_theta;

    *ratio_ante = c_cros / c_ante;
    *ratio_post = c_cros / c_post;

    delta_cosine           = cosine_theta - cosine_theta_equil;
    term_theta_theta_equil = 1 - cosine_theta * cosine_theta_equil;
    *prefactor             = -(k_bending) * delta_cosine * norm * term_theta_theta_equil / (sine_theta_sq * sine_theta_sq);

    *v = k_bending * 0.5 * delta_cosine * delta_cosine / sine_theta_sq;

}



void compute_factors_restrdihs(int type,  const t_iparams forceparams[],
                               rvec delta_ante, rvec delta_crnt, rvec delta_post,
                               real *factor_phi_middle_ante_ante, real *factor_phi_middle_ante_crnt,
                               real *factor_phi_middle_ante_post, real *factor_phi_middle_post_ante,
                               real *factor_phi_extrem_ante_ante, real *factor_phi_middle_post_crnt,
                               real *factor_phi_extrem_ante_crnt, real *factor_phi_middle_post_post,
                               real *factor_phi_extrem_ante_post, real *factor_phi_extrem_post_ante,
                               real *factor_phi_extrem_post_crnt,  real *factor_phi_extrem_post_post,
                               real *prefactor_phi, real *v)
{
    real  phi, cos_phi;
    real  k_torsion;
    real  norm_phi, sine_phi0;
    rvec  vec_temp;
    real  c_self_ante, c_cros_ante, c_cros_acrs, c_self_crnt;
    real  c_self_post, sine_phi_sq, phi0, c_cros_post;
    real  d_ante, c_prod, d_post;
    real  cosine_phi, cosine_phi0;
    real  delta_cosine, term_phi_phi0;
    real  ratio_phi_ante, ratio_phi_post;

    phi0        = forceparams[type].pdihs.phiA * DEG2RAD;
    cosine_phi0 = cos(phi0);
    sine_phi0   = sin(phi0);
    k_torsion   = forceparams[type].pdihs.cpA;

    c_self_ante = iprod(delta_ante, delta_ante);
    c_self_crnt = iprod(delta_crnt, delta_crnt);
    c_self_post = iprod(delta_post, delta_post);
    c_cros_ante = iprod(delta_ante, delta_crnt);
    c_cros_acrs = iprod(delta_ante, delta_post);
    c_cros_post = iprod(delta_crnt, delta_post);
    c_prod      = c_cros_ante * c_cros_post - c_self_crnt * c_cros_acrs;
    d_ante      = c_self_ante * c_self_crnt - c_cros_ante * c_cros_ante;
    d_post      = c_self_post * c_self_crnt - c_cros_post * c_cros_post;

    /*      Take good care of the singularity !!! */
    if (d_ante < GMX_REAL_EPS)
    {
        d_ante = GMX_REAL_EPS;
    }
    if (d_post < GMX_REAL_EPS)
    {
        d_post = GMX_REAL_EPS;
    }

    norm_phi    = gmx_invsqrt(d_ante * d_post);
    cosine_phi  = c_prod * norm_phi;
    sine_phi_sq = 1.0 - cosine_phi * cosine_phi;

    /*      Take good care of trig functions !!! */
    if (sine_phi_sq < 0.0)
    {
        sine_phi_sq = 0.0;
    }

    delta_cosine  = cosine_phi - cosine_phi0;
    term_phi_phi0 = 1 - cosine_phi * cosine_phi0;


    /*      Computation of ratios */
    ratio_phi_ante = c_prod / d_ante;
    ratio_phi_post = c_prod / d_post;

    /*      Computation of the prefactor */
    *prefactor_phi = -(k_torsion) * delta_cosine * norm_phi * term_phi_phi0 / (sine_phi_sq * sine_phi_sq);

    /*      Computation of factors (important for gaining speed) */
    *factor_phi_extrem_ante_ante = ratio_phi_ante * c_self_crnt;
    *factor_phi_extrem_ante_crnt = -c_cros_post - ratio_phi_ante * c_cros_ante;
    *factor_phi_extrem_ante_post = c_self_crnt;
    *factor_phi_middle_ante_ante = -c_cros_post - ratio_phi_ante * (c_self_crnt + c_cros_ante);
    *factor_phi_middle_ante_crnt = c_cros_post + c_cros_acrs * 2.0 + ratio_phi_ante * (c_self_ante + c_cros_ante) + ratio_phi_post * c_self_post;
    *factor_phi_middle_ante_post = -(c_cros_ante + c_self_crnt) - ratio_phi_post * c_cros_post;
    *factor_phi_middle_post_ante = c_cros_post + c_self_crnt + ratio_phi_ante * c_cros_ante;
    *factor_phi_middle_post_crnt = -(c_cros_ante + c_cros_acrs * 2.0)- ratio_phi_ante * c_self_ante - ratio_phi_post * (c_self_post + c_cros_post);
    *factor_phi_middle_post_post = c_cros_ante + ratio_phi_post * (c_self_crnt + c_cros_post);
    *factor_phi_extrem_post_ante = -c_self_crnt;
    *factor_phi_extrem_post_crnt = c_cros_ante + ratio_phi_post * c_cros_post;
    *factor_phi_extrem_post_post = -ratio_phi_post * c_self_crnt;

    *v = k_torsion * 0.5 * delta_cosine * delta_cosine / sine_phi_sq;
}




void compute_factors_cbtdihs(int type,  const t_iparams forceparams[],
                             rvec delta_ante, rvec delta_crnt, rvec delta_post,
                             rvec f_theta_ante_middle_ante, rvec f_theta_ante_middle_post,
                             rvec f_theta_post_middle_ante, rvec f_theta_ante_extrem_ante,
                             rvec f_theta_post_middle_post, rvec f_phi_middle_ante,
                             rvec f_phi_middle_post, rvec f_phi_extrem_ante,
                             rvec f_phi_extrem_post, rvec f_theta_post_extrem_post, real * v)

{
    real c_cros_ante, c_cros_acrs, c_self_crnt;
    real sine_theta_post, norm_theta_post, norm_phi;
    real cosine_theta_ante, cosine_theta_post, cosine_phi;
    real sine_theta_ante_sq, sine_theta_post_sq, c_self_ante;
    real factor_phi_middle_ante_ante, factor_phi_middle_ante_crnt;
    real factor_phi_middle_ante_post;
    real factor_phi_middle_post_ante, factor_phi_extrem_ante_ante;
    real factor_phi_middle_post_crnt;
    real factor_phi_extrem_ante_crnt;
    real factor_phi_middle_post_post;
    real factor_phi_extrem_ante_post;
    real factor_phi_extrem_post_ante, c_self_post, c_cros_post;
    real factor_phi_extrem_post_crnt,  factor_phi_extrem_post_post, prefactor_theta_ante;
    real prefactor_theta_post, prefactor_phi, ratio_theta_ante_ante;
    real d_ante, ratio_theta_ante_crnt, c_prod, d_post;
    real ratio_theta_post_crnt, ratio_theta_post_post;
    real ratio_phi_ante, ratio_phi_post, sine_theta_ante;
    real norm_theta_ante;
    real torsion_coef[NR_CBTDIHS];
    real r1, r2;
    int  j, d;

    for (j = 0; (j < NR_CBTDIHS); j++)
    {
        torsion_coef[j]  = forceparams[type].cbtdihs.cbtcA[j];
    }

    c_self_ante = iprod(delta_ante, delta_ante);
    c_self_crnt = iprod(delta_crnt, delta_crnt);
    c_self_post = iprod(delta_post, delta_post);
    c_cros_ante = iprod(delta_ante, delta_crnt);
    c_cros_acrs = iprod(delta_ante, delta_post);
    c_cros_post = iprod(delta_crnt, delta_post);
    c_prod      = c_cros_ante * c_cros_post - c_self_crnt * c_cros_acrs;
    d_ante      = c_self_ante * c_self_crnt - c_cros_ante * c_cros_ante;
    d_post      = c_self_post * c_self_crnt - c_cros_post * c_cros_post;

    /*      Take good care of the singularity !!! */
    if (d_ante < GMX_REAL_EPS)
    {
        d_ante = GMX_REAL_EPS;
    }
    if (d_post < GMX_REAL_EPS)
    {
        d_post = GMX_REAL_EPS;
    }

    norm_phi           = gmx_invsqrt(d_ante * d_post);
    norm_theta_ante    = gmx_invsqrt(c_self_ante * c_self_crnt);
    norm_theta_post    = gmx_invsqrt(c_self_crnt * c_self_post);
    cosine_phi         = c_prod * norm_phi;
    cosine_theta_ante  = c_cros_ante * norm_theta_ante;
    cosine_theta_post  = c_cros_post * norm_theta_post;
    sine_theta_ante_sq = 1 - cosine_theta_ante * cosine_theta_ante;
    sine_theta_post_sq = 1 - cosine_theta_post * cosine_theta_post;

    /*      Take good care of trig functions !!! */
    if (sine_theta_ante_sq < 0.0)
    {
        sine_theta_ante_sq = 0.0;
    }
    if (sine_theta_post_sq < 0.0)
    {
        sine_theta_post_sq = 0.0;
    }

    sine_theta_ante = sqrt(sine_theta_ante_sq);
    sine_theta_post = sqrt(sine_theta_post_sq);

    /*      Computing forces due to the derivative with repect to the dehidral angle 'PHI': 4 beads */
    /*      Computation of ratios */
    ratio_phi_ante = c_prod / d_ante;
    ratio_phi_post = c_prod / d_post;

    /*       Computation of the prefactor */
    /*      Computing 2nd power */
    r1 = cosine_phi;

    prefactor_phi = -torsion_coef[0] * norm_phi * (torsion_coef[2] + torsion_coef[3] * 2.0 * cosine_phi + torsion_coef[4] * 3.0 * (r1 * r1) + 4*torsion_coef[5]*r1*r1*r1) *
        sine_theta_ante_sq * sine_theta_ante * sine_theta_post_sq * sine_theta_post;

    /*      Computation of factors (important for gaining speed) */
    factor_phi_extrem_ante_ante = ratio_phi_ante * c_self_crnt;
    factor_phi_extrem_ante_crnt = -c_cros_post - ratio_phi_ante * c_cros_ante;
    factor_phi_extrem_ante_post = c_self_crnt;
    factor_phi_middle_ante_ante = -c_cros_post - ratio_phi_ante * (c_self_crnt + c_cros_ante);
    factor_phi_middle_ante_crnt = c_cros_post + c_cros_acrs * 2.0 + ratio_phi_ante * (c_self_ante + c_cros_ante) +  ratio_phi_post * c_self_post;
    factor_phi_middle_ante_post = -(c_cros_ante + c_self_crnt) - ratio_phi_post * c_cros_post;
    factor_phi_middle_post_ante = c_cros_post + c_self_crnt + ratio_phi_ante * c_cros_ante;
    factor_phi_middle_post_crnt = -(c_cros_ante + c_cros_acrs * 2.0) - ratio_phi_ante * c_self_ante - ratio_phi_post * (c_self_post + c_cros_post);
    factor_phi_middle_post_post = c_cros_ante + ratio_phi_post * (c_self_crnt + c_cros_post);
    factor_phi_extrem_post_ante = -c_self_crnt;
    factor_phi_extrem_post_crnt = c_cros_ante + ratio_phi_post * c_cros_post;
    factor_phi_extrem_post_post = -ratio_phi_post * c_self_crnt;

    /* Computation of forces */
    for (d = 0; d < DIM; d++)
    {
        f_phi_extrem_ante[d] = prefactor_phi * (factor_phi_extrem_ante_ante * delta_ante[d] + factor_phi_extrem_ante_crnt * delta_crnt[d] + factor_phi_extrem_ante_post * delta_post[d]);
        f_phi_middle_ante[d] = prefactor_phi * (factor_phi_middle_ante_ante * delta_ante[d] + factor_phi_middle_ante_crnt * delta_crnt[d] + factor_phi_middle_ante_post * delta_post[d]);
        f_phi_middle_post[d] = prefactor_phi * (factor_phi_middle_post_ante * delta_ante[d] + factor_phi_middle_post_crnt * delta_crnt[d] + factor_phi_middle_post_post * delta_post[d]);
        f_phi_extrem_post[d] = prefactor_phi * (factor_phi_extrem_post_ante * delta_ante[d] + factor_phi_extrem_post_crnt * delta_crnt[d] + factor_phi_extrem_post_post * delta_post[d]);
    }

    /*      Computing forces due to the derivative with repect to the bending angle 'THETA_ANTE': 3 beads */
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


    /*      Computation of forces */
    for (d = 0; d < DIM; d++)
    {
        f_theta_ante_extrem_ante[d] = prefactor_theta_ante * (ratio_theta_ante_ante * delta_ante[d] - delta_crnt[d]);
        f_theta_ante_middle_ante[d] = prefactor_theta_ante * ((ratio_theta_ante_crnt + 1.0) * delta_crnt[d] - (ratio_theta_ante_ante + 1.0) * delta_ante[d]);
        f_theta_ante_middle_post[d] = prefactor_theta_ante * (delta_ante[d] - ratio_theta_ante_crnt * delta_crnt[d]);
    }

    /*      Computing forces due to the derivative with repect to the bending angle 'THETA_POST': 3 beads */
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


    /*      Computation of forces */
    for (d = 0; d < DIM; d++)
    {
        f_theta_post_middle_ante[d] = prefactor_theta_post * (ratio_theta_post_crnt * delta_crnt[d] - delta_post[d]);
        f_theta_post_middle_post[d] = prefactor_theta_post * ((ratio_theta_post_post + 1.0) * delta_post[d] - (ratio_theta_post_crnt + 1.0) * delta_crnt[d]);
        f_theta_post_extrem_post[d] = prefactor_theta_post * (delta_crnt[d] - ratio_theta_post_post * delta_post[d]);
    }
    r1 = cosine_phi;
    r2 = cosine_phi;
    *v = torsion_coef[0] * (torsion_coef[1] + torsion_coef[2] * cosine_phi + torsion_coef[3] * (r1 * r1) +
                            torsion_coef[4] * (r2 * (r2 * r2)) + torsion_coef[5] * (r2 * (r2 * (r2 * r2)))) * sine_theta_ante_sq *
        sine_theta_ante * sine_theta_post_sq * sine_theta_post;


}
