/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 *
 *
 * \brief
 * This file contains function declarations necessary
   for computations of forces due to restricted angle, restricted dihedral and
   combined bending-torsion potentials.
 *
 * \author Nicolae Goga
 *
 * \ingroup module_listed_forces
 */

#ifndef GMX_LISTED_FORCES_RESTCBT_H
#define GMX_LISTED_FORCES_RESTCBT_H

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/real.h"

/*! \brief This function computes factors needed for restricted angle potentials.
 *
 * The restricted angle potential is used in coarse-grained simulations to avoid singularities
 * when three particles align and the dihedral angle and dihedral potential cannot be calculated.
 * This potential is calculated using the formula:
 * \f[V_{\rm ReB}(\theta_i) = \frac{1}{2} k_{\theta} \frac{(\cos\theta_i - \cos\theta_0)^2}{\sin^2\theta_i}\f]
 * (see section "Restricted Bending Potential" from the manual).
 * The derivative of the restricted angle potential is calculated as:
 * \f[\frac{\partial V_{\rm ReB}(\theta_i)} {\partial \vec r_{k}} = \frac{dV_{\rm ReB}(\theta_i)}{dcos\theta_i} \frac{\partial cos\theta_i}{\partial \vec r_{k}}\f]
 * where all the derivatives of the bending angle with respect to Cartesian coordinates are calculated as in Allen & Tildesley (pp. 330-332)
 *
 *  \param[in]  type           type of force parameters
 *  \param[in]  forceparams    array of parameters for force computations
 *  \param[in]  delta_ante     distance between the first two particles
 *  \param[in]  delta_post     distance between the last two particles
 *  \param[out] prefactor      common term that comes in front of each force
 *  \param[out] ratio_ante     ratio of scalar products of delta_ante with delta_post
                              and delta_ante with delta_ante
 *  \param[out] ratio_post    ratio of scalar products of delta_ante with delta_post
                              and delta_post with delta_ante
 *  \param[out] v              contribution to energy   (see formula above)
 */


void compute_factors_restangles(int             type,
                                const t_iparams forceparams[],
                                rvec            delta_ante,
                                rvec            delta_post,
                                double*         prefactor,
                                double*         ratio_ante,
                                double*         ratio_post,
                                real*           v);


/*! \brief Compute factors for restricted dihedral potentials.
 *
 * The restricted dihedral potential is the equivalent of the restricted bending potential
 * for the dihedral angle. It imposes the dihedral angle to have only one equilibrium value.
 * This potential is calculated using the formula:
 * \f[V_{\rm ReT}(\phi_i) = \frac{1}{2} k_{\phi} \frac{(\cos\phi_i - \cos\phi_0)^2}{\sin^2\phi_i}\f]
 * (see section "Proper dihedrals: Restricted torsion potential" from the manual).
 * The derivative of the restricted dihedral potential is calculated as:
 * \f[\frac{\partial V_{\rm ReT}(\phi_i)} {\partial \vec r_{k}} = \frac{dV_{\rm ReT}(\phi_i)}{dcos\phi_i} \frac{\partial cos\phi_i}{\partial \vec r_{k}}\f]
 * where all the derivatives of the dihedral angle with respect to Cartesian coordinates
 * are calculated as in Allen & Tildesley (pp. 330-332). Factors factor_phi_*  are coming from the
 * derivatives of the torsion angle (phi) with respect to the beads ai, aj, ak, al, (four) coordinates
 * and they are multiplied in the force computations with the particle distance
 * stored in parameters delta_ante, delta_crnt, delta_post.
 *
 *  \param[in]  type                             type of force parameters
 *  \param[in]  forceparams                      array of parameters for force computations
 *  \param[in]  delta_ante                       distance between the first two particles
 *  \param[in]  delta_crnt                       distance between the middle pair of particles
 *  \param[in]  delta_post                       distance between the last two particles
 *  \param[out] factor_phi_ai_ante               force factor for particle ai multiplied with delta_ante
 *  \param[out] factor_phi_ai_crnt               force factor for particle ai multiplied with delta_crnt
 *  \param[out] factor_phi_ai_post               force factor for particle ai multiplied with delta_post
 *  \param[out] factor_phi_aj_ante               force factor for particle aj multiplied with delta_ante
 *  \param[out] factor_phi_aj_crnt               force factor for particle aj multiplied with delta_crnt
 *  \param[out] factor_phi_aj_post               force factor for particle aj multiplied with delta_post
 *  \param[out] factor_phi_ak_ante               force factor for particle ak multiplied with delta_ante
 *  \param[out] factor_phi_ak_crnt               force factor for particle ak multiplied with delta_crnt
 *  \param[out] factor_phi_ak_post               force factor for particle ak multiplied with delta_post
 *  \param[out] factor_phi_al_ante               force factor for particle al multiplied with delta_ante
 *  \param[out] factor_phi_al_crnt               force factor for particle al multiplied with delta_crnt
 *  \param[out] factor_phi_al_post               force factor for particle al multiplied with delta_post
 *  \param[out] prefactor_phi                    multiplication constant of the torsion force
 *  \param[out] v                                contribution to energy  (see formula above)
 */

void compute_factors_restrdihs(int             type,
                               const t_iparams forceparams[],
                               rvec            delta_ante,
                               rvec            delta_crnt,
                               rvec            delta_post,
                               real*           factor_phi_ai_ante,
                               real*           factor_phi_ai_crnt,
                               real*           factor_phi_ai_post,
                               real*           factor_phi_aj_ante,
                               real*           factor_phi_aj_crnt,
                               real*           factor_phi_aj_post,
                               real*           factor_phi_ak_ante,
                               real*           factor_phi_ak_crnt,
                               real*           factor_phi_ak_post,
                               real*           factor_phi_al_ante,
                               real*           factor_phi_al_crnt,
                               real*           factor_phi_al_post,
                               real*           prefactor_phi,
                               real*           v);

/*! \brief Compute factors for combined bending-torsion (CBT) potentials.
 *
 * The combined bending-torsion potential goes to zero in a very smooth manner, eliminating the numerical
 * instabilities, when three coarse-grained particles align and the dihedral angle and standard
 * dihedral potentials cannot be calculated. The CBT potential is calculated using the formula:
 * \f[V_{\rm CBT}(\theta_{i-1}, \theta_i, \phi_i) = k_{\phi} \sin^3\theta_{i-1} \sin^3\theta_{i}
 * \sum_{n=0}^4 { a_n \cos^n\phi_i}\f] (see section "Proper dihedrals: Combined bending-torsion potential" from the manual).
 * It contains in its expression not only the dihedral angle \f$\phi\f$
 * but also \f$\theta_{i-1}\f$ (denoted as theta_ante below) and \f$\theta_{i}\f$ (denoted as theta_post below)
 * --- the adjacent bending angles. The derivative of the CBT potential is calculated as:
 * \f[\frac{\partial V_{\rm CBT}(\theta_{i-1},\theta_i,\phi_i)} {\partial \vec r_{l}} = \frac{\partial V_
 * {\rm CBT}}{\partial \theta_{i-1}} \frac{\partial \theta_{i-1}}{\partial \vec r_{l}} +
 * \frac{\partial V_{\rm CBT}}{\partial \phi_{i    }} \frac{\partial \phi_{i    }}{\partial \vec r_{l}}\f]
 * where all the derivatives of the angles with respect to Cartesian coordinates are calculated as
 * in Allen & Tildesley (pp. 330-332). Factors f_phi_* come from  the derivatives of the torsion angle
 * with respect to the beads ai, aj, ak, al (four) coordinates; f_theta_ante_* come from the derivatives of
 * the bending angle theta_ante (theta_{i-1} in formula above) with respect to the beads ai, aj, ak (three
 * particles) coordinates and f_theta_post_* come from the derivatives of  the bending angle theta_post
 * (theta_{i} in formula above) with respect to the beads aj, ak, al (three particles) coordinates.
 *
 *  \param[in]  type                             type of force parameters
 *  \param[in]  forceparams                      array of parameters for force computations
 *  \param[in]  delta_ante                       distance between the first two particles
 *  \param[in]  delta_crnt                       distance between the middle pair of particles
 *  \param[in]  delta_post                       distance between the last two particles
 *  \param[out]  f_phi_ai                        force for particle ai due to derivative in phi angle
 *  \param[out]  f_phi_aj                        force for particle aj due to derivative in phi angle
 *  \param[out]  f_phi_ak                        force for particle ak due to derivative in phi angle
 *  \param[out]  f_phi_al                        force for particle al due to derivative in phi angle
 *  \param[out]  f_theta_ante_ai                 force for particle ai due to derivative in theta_ante angle
 *  \param[out]  f_theta_ante_aj                 force for particle aj due to derivative in theta_ante angle
 *  \param[out]  f_theta_ante_ak                 force for particle ak due to derivative in theta_ante angle
 *  \param[out]  f_theta_post_aj                 force for particle aj due to derivative in theta_post angle
 *  \param[out]  f_theta_post_ak                 force for particle ak due to derivative in theta_post angle
 *  \param[out]  f_theta_post_al                 force for particle al due to derivative in theta_psot angle
 *  \param[out] v                                contribution to energy (see formula above)
 */

void compute_factors_cbtdihs(int             type,
                             const t_iparams forceparams[],
                             rvec            delta_ante,
                             rvec            delta_crnt,
                             rvec            delta_post,
                             rvec            f_phi_ai,
                             rvec            f_phi_aj,
                             rvec            f_phi_ak,
                             rvec            f_phi_al,
                             rvec            f_theta_ante_ai,
                             rvec            f_theta_ante_aj,
                             rvec            f_theta_ante_ak,
                             rvec            f_theta_post_aj,
                             rvec            f_theta_post_ak,
                             rvec            f_theta_post_al,
                             real*           v);

#endif
