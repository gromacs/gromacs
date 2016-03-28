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

/*! \internal \file
 *
 * \brief
 * This file contains function declarations needed internally by AWH
 * for preparing and writing output data to an energy frame.
 *
 * \author Viveca Lindahl
 * \ingroup module_awh
 */

#ifndef GMX_AWH_ENERGY_WRITER_H
#define GMX_AWH_ENERGY_WRITER_H

#include "gromacs/fileio/enxio.h"
#include "gromacs/utility/basedefinitions.h"

struct awh_energywriter_t;
struct awh_t;
struct awh_params_t;
struct gmx_multisim_t;

/*! \brief Allocate, initialize and return an AWH writer.
 *
 * \param[in] nstout            Number of steps per output writing.
 * \param[in] awh               AWH working struct.
 * \param[in] awh_params        AWH parameters.
 * \returns the initialized writer.
 */
awh_energywriter_t *init_awh_energywriter(int nstout,
                                          const awh_t *awh, const awh_params_t *awh_params);


/*! \brief Query if output should be written at the given step.
 *
 * \param[in,out] energywriter   The writer struct.
 * \param[in]     awh_params AWH parameters.
 * \param[in]     awh        AWH working struct.
 * \param[in]     ms             Struct for multi-simulation communication, needed for bias sharing replicas.
 */
void prep_awh_output(awh_energywriter_t *energywriter, const awh_params_t *awh_params,
                     const awh_t *awh, const gmx_multisim_t *ms);


/*! \brief Query if output should be written at the given step.
 *
 * \param[in] step          Time step.
 * \param[in] energywriter  Writer struct.
 * \returns true if the step is a writing step.
 */
bool time_to_write(gmx_int64_t step, const awh_energywriter_t *energywriter);


/*! \brief Fill the AWH data block of an energy frame with data (if there is any).
 *
 * \param[in,out] frame     Energy data frame.
 * \param[in,out] energywriter  Writer struct.
 */
void write_awh_to_frame(t_enxframe *frame, awh_energywriter_t *energywriter);

#endif  /* GMX_AWH_ENERGY_WRITER_H */
