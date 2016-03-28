/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief
 * This file contains datatypes and function declarations necessary
   for the mdrun replica exchange code to interface with the AWH replica exchange code.
 *
 * \author Viveca Lindahl
 *
 * \inlibraryapi
 * \ingroup module_awh
 */

#ifndef GMX_AWH_REPLEX_H
#define GMX_AWH_REPLEX_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct awh_params_t;
struct awh_t;
struct gmx_multisim_t;
struct awh_replex_t;

/*! \brief Allocate, initialize and return an AWH replica exchange struct.
 *
 * \param[in] awh_params    AWH input parameters.
 * \param[in] awh           AWH working struct.
 * \param[in] ms                Struct for multisimulation communication.
 * \param[in] exchange_indices  Replica indices in the order of increasing exchange parameter.
 * \param[in] bMulti_exchange   If true, doing multiple exchanges per exchange interval, otherwise nearest neighbor exchange.
 * \returns the initialized AWH replica exchange struct.
 */
awh_replex_t *init_awh_replica_exchange(const awh_params_t *awh_params, const awh_t *awh, const gmx_multisim_t *ms,
                                        const int *exchange_indices, int bMulti_exchange);

/*! \brief Update the replica exchange struct with the current AWH bias.
 *
 * \param[in] awh           AWH working struct.
 * \param[in] replex        AWH replica exchange struct.
 * \param[in] ms            Struct for multisimulation communication.
 * \returns true if the update was successful.
 */
gmx_bool update_awh_replica_exchange_state(const awh_t *awh, awh_replex_t *replex, const gmx_multisim_t *ms);

/*! \brief Query if two replicas are up for exchanging.
 *
 * For replicas that are equivalent (i.e. have the same target distribution) it may be preferable to not attempt exchanges.
 *
 * \param[in] replex            AWH replica exchange struct.
 * \param[in] replica_id1       Id of the first replica.
 * \param[in] replica_id2       Id of the second replica.
 * \returns true if the replicas can exchange.
 */
gmx_bool awh_replicas_can_exchange(const awh_replex_t *replex, int replica_id1, int replica_id2);

/*! \brief Calculate the AWH bias exchange probability delta.
 *
 * The exchange probability for swapping the configurations of replica A and replica B
 * depends on the states (biases) and the configuration (coordinate values) of A and B.
 * The exchange probability is expressed as exp(-delta).
 *
 * Note that when calling this function the index given for the configuration of A and
 * the state of A is not generally the same. When multiple swaps are attempted per
 * exchange interval the indices may already have been permuted.
 *
 * \param[in] replex            AWH replica exchange struct.
 * \param[in] config_A          Replica id for configuration A.
 * \param[in] config_B          Replica id for configuration B.
 * \param[in] state_A           Replica id for state A.
 * \param[in] state_B           Replica id for state B.
 * \returns delta.
 */
real calc_awh_replica_exchange_delta(const awh_replex_t *replex, int config_A, int config_B, int state_A, int state_B);

#endif  /* GMX_AWH_REPLEX_H */
