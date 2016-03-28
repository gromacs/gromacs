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

/*! \libinternal \file
 *
 * \brief
 * Declares functions necessary for applying an AWH bias.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 */

#ifndef GMX_AWH_H
#define GMX_AWH_H

#include <cstdio>

#include "gromacs/fileio/enxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct t_commrec;
struct gmx_multisim_t;
struct t_mdatoms;
struct t_inputrec;
struct AwhBiasCollection;
struct awh_params_t;
struct awh_history_t;

/*! \brief Allocate, initialize and return an AWH working struct for mdrun.
 *
 * This function also handles initializing of the AWH history and restoring
 * a checkpointed history to the current state if needed.
 *
 * \param[in,out] fplog               General output file, normally md.log.
 * \param[in] ir                      General input parameters.
 * \param[in] cr                      Struct for communication.
 * \param[in] awh_params              AWH input parameters.
 * \param[in,out] awh_history_ptr     AWH bias history to initialize.
 * \param[in,out] pull_work           Pull struct which AWH will register the bias into.
 * \param[in] startingFromCheckpoint  True if this this is a continuation run.
 * \returns the initialized AWH struct.
 */
AwhBiasCollection *init_awh_md(FILE                    *fplog,
                               const t_inputrec        *ir,
                               const t_commrec         *cr,
                               const awh_params_t      *awh_params,
                               awh_history_t          **awh_history_ptr,
                               struct pull_t           *pull_work,
                               bool                     startingFromCheckpoint);

/*! \brief Do an AWH biasing update.
 *
 * An AWH update includes:
 * - sample observables related to the AWH coordinates
 * - feed the forces and potential resulting from the AWH bias to the pull struct.
 * - update the AWH bias
 * - prepare AWH output data (this is later written to file using write_awh_to_energyframe).
 *
 * Note: requires that pull_potential from pull.h has been called first
 * since AWH need the updated center-of-masses of the pull coordinates.
 *
 * \param[in,out] awh      AWH working struct.
 * \param[in] awh_params   AWH input parameters.
 * \param[in,out] pull_work    Pull working struct.
 * \param[in]   mdatoms        Atom properties.
 * \param[in] ePBC             Type of periodic boundary conditions.
 * \param[in] box              Box vectors.
 * \param[in,out] force        Forces.
 * \param[in,out] virial       The virial
 * \param[in] ms               Struct for multi-simulation communication.
 * \param[in]     t            Time.
 * \param[in]     step         Time step.
 * \param[in,out] wallcycle    Wallcycle counter.
 * \param[in,out] fplog        General output file, normally md.log.
 * \returns the potential energy for the bias.
 */
real update_awh(AwhBiasCollection      *awh,
                const awh_params_t     *awh_params,
                struct pull_t          *pull_work,
                int                     ePBC,
                const t_mdatoms        *mdatoms,
                const matrix            box,
                rvec                   *force,
                tensor                  virial,
                const gmx_multisim_t   *ms,
                double                  t,
                gmx_int64_t             step,
                struct gmx_wallcycle   *wallcycle,
                FILE                   *fplog);

/*! \brief
 * Update the AWH bias history for checkpointing.
 *
 * \param[in,out] awh_history  AWH bias history to set.
 * \param[in] awh          AWH state to read.
 */
void update_awh_history(awh_history_t            *awh_history,
                        const AwhBiasCollection  *awh);

/*! \brief Fills the AWH data block of an energy frame with data at certain steps.
 *
 * \param[in,out] fr      Energy data frame.
 * \param[in,out] awh     AWH working struct.
 */
void write_awh_to_energyframe(t_enxframe *fr, AwhBiasCollection *awh);

#endif /* GMX_AWH_H */
