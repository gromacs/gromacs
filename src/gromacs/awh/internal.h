/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Functions and declarations for internal use in the AWH module.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_INTERNAL_H
#define GMX_AWH_INTERNAL_H

#include "types.h" /* This currently needed for awh_dvec */

//! Linewidth used for warning output
static const int linewidth = 78;

//! Indent used for warning output
static const int indent    = 0;

struct AwhBias;
struct awh_dim_params_t;
struct gmx_multisim_t;
struct pull_params_t;
struct t_inputrec;
struct t_commrec;

/*! \brief
 * Calculates the convolved bias for a given coordinate value.
 *
 * The convolved bias is the effective bias acting on the coordinate.
 * Since the bias here has arbitrary normalization, this only makes
 * sense as a relative, to other coordinate values, measure of the bias.
 *
 * \note If it turns out to be costly to calculate this pointwise
 * the convolved bias for the whole grid could be returned instead.
 *
 * \param[in] awh_bias       AWH bias.
 * \param[in] coord_value    Coordinate value.
 * \returns the convolved bias >= -GMX_DOUBLE_MAX.
 */
double calc_convolved_bias(const AwhBias *awh_bias, const awh_dvec coord_value);

/*! \brief
 * Sets the given array with PMF values.
 *
 * Points outside of the biasing target region will get PMF = GMX_DOUBLE_MAX.
 * In the simplest case the PMF is simply the negative of the PMF histogram.
 * If there are sharing replicas however, histograms need to be summed
 * across multiple simulations. The output PMF is not normalized.
 *
 * \param[in] awh_bias    The AWH bias.
 * \param[in] ms          Struct for multi-simulation communication, needed for bias sharing replicas.
 * \param[in,out] pmf     Array returned will be of the same length as the AWH grid to store the PMF in.
 */
void getPmf(const AwhBias &awh_bias, const gmx_multisim_t *ms,
            std::vector<float> *pmf);

/*! \brief Convolves the given PMF using the given AWH bias.
 *
 * \param[in] awh_bias          The AWH bias.
 * \param[in] ms                Struct for multi-simulation communication, needed for bias sharing replicas.
 * \param[in,out] convolvedPmf  Array returned will be of the same length as the AWH grid to store the convolved PMF in.
 */
void getConvolvedPmf(const AwhBias &awh_bias, const gmx_multisim_t *ms,
                     std::vector<float> *convolvedPmf);

/*! \brief Allocate, initialize and return an AWH struct.
 *
 * This is used to get a temporary AWH working struct at preprocessing time.
 *
 * \param[in,out] fplog           General output log file (or NULL).
 * \param[in]     ir              General input parameters.
 * \param[in]     cr              Struct for communication (or NULL).
 * \param[in]     awh_params      AWH input parameters.
 * \returns the initialized AWH struct.
 */
AwhBiasCollection *init_awh(FILE                    *fplog,
                            const t_inputrec        *ir,
                            const t_commrec         *cr,
                            const awh_params_t      *awh_params);

/*! \brief Register the AWH biased coordinates with pull.
 *
 * \param[in]     awh_params      The AWH parameters.
 * \param[in,out] pull_work       Pull struct which AWH will register the bias into.
 */
void register_bias_with_pull(const awh_params_t      *awh_params,
                             struct pull_t           *pull_work);


/*! \brief Initialize the target and bias values for the given AWH bias.
 *
 * \param[in,out]     awh_bias             The AWH bias.
 */
void initTargetAndBias(AwhBias *awh_bias);


#endif  /* GMX_AWH_INTERNAL_H */
