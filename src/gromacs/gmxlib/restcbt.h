/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _restcbt_h
#define _restcbt_h


#ifdef __cplusplus
extern "C" {
#endif


/* \brief This function computes factors needed for restricted angle potential.
*
*  \param[in]  type           type of force paramters
*  \param[in]  forceparams    array of parameters for force computations
*  \param[in]  lambda         parameter used for computation of bending constant 
*  \param[in]  delta_ante     position diference between the first two particles
*  \param[in]  delta_post     position difference between the last two particles
*  \param[out] prefactor      common term that comes in front of each force
*  \param[out] ratio_ante     ratio betwen scallar products of delta_ante with delta_post 
                              and delta_ante with delta_ante
*  \param[out] ration_post    ratio between scalar products of delta_ante  with delta_post 
                              and delta_post with delta_ante
*  \param[out] v              contribution to energy   
*/

void compute_factors_restangles(int type, const t_iparams forceparams[], real lambda,
                                rvec delta_ante,  rvec delta_post,
                                real *prefactor, real *ratio_ante, real *ratio_post, real *v);


/* Compute factors for restricted dihedral potential */
void compute_factors_restrdihs(int type,  const t_iparams forceparams[], real lambda,
                               rvec delta_ante, rvec delta_crnt, rvec delta_post,
                               real *factor_phi_middle_ante_ante, real *factor_phi_middle_ante_crnt,
                               real *factor_phi_middle_ante_post, real *factor_phi_middle_post_ante,
                               real *factor_phi_extrem_ante_ante, real *factor_phi_middle_post_crnt,
                               real *factor_phi_extrem_ante_crnt, real *factor_phi_middle_post_post,
                               real *factor_phi_extrem_ante_post, real *factor_phi_extrem_post_ante,
                               real *factor_phi_extrem_post_crnt,  real *factor_phi_extrem_post_post,
                               real *prefactor_phi, real *v);


/* Compute factors for CBT potential */
void compute_factors_cbtdihs(int type,  const t_iparams forceparams[],
                             rvec delta_ante, rvec delta_crnt, rvec delta_post,
                             rvec f_theta_ante_middle_ante, rvec f_theta_ante_middle_post,
                             rvec f_theta_post_middle_ante, rvec f_theta_ante_extrem_ante,
                             rvec f_theta_post_middle_post, rvec f_phi_middle_ante,
                             rvec f_phi_middle_post, rvec f_phi_extrem_ante,
                             rvec f_phi_extrem_post, rvec f_theta_post_extrem_post, real * v);


#ifdef __cplusplus
}
#endif

#endif  /* _restcbt_h */

