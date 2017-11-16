/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Declares functions needed for reading, initializing and setting the AWH parameter data types.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 * \ingroup module_awh
 */

#ifndef GMX_AWH_READPARAMS_H
#define GMX_AWH_READPARAMS_H

#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vectypes.h"

struct t_grpopts;
struct t_inputrec;
struct pull_params_t;
struct pull_t;

namespace gmx
{
struct AwhParams;

/*! \brief Allocate, initialize and check the AWH parameters with values from the input file.
 *
 * \param[in,out] ninp_p       Number of read input file entries.
 * \param[in,out] inp_p        Input file entries.
 * \param[in]     inputrec     Input parameter struct.
 * \param[in,out] wi           Struct for bookeeping warnings.
 * \returns AWH parameters.
 */
AwhParams *readAndCheckAwhParams(int               *ninp_p,
                                 t_inpfile        **inp_p,
                                 const t_inputrec  *inputrec,
                                 warninp_t          wi);


/*! \brief
 * Sets AWH parameters that need state parameters such as the box vectors.
 *
 * \param[in,out] awhParams             AWH parameters.
 * \param[in]     pull_params           Pull parameters.
 * \param[in,out] pull_work             Pull working struct to register AWH bias in.
 * \param[in]     box                   Box vectors.
 * \param[in]     ePBC                  Periodic boundary conditions enum.
 * \param[in]     inputrecGroupOptions  Parameters for atom groups.
 * \param[in,out] wi                    Struct for bookeeping warnings.
 *
 * \note This function currently relies on the function set_pull_init to have been called.
 */
void setStateDependentAwhParams(AwhParams           *awhParams,
                                const pull_params_t *pull_params,
                                pull_t              *pull_work,
                                const matrix         box,
                                int                  ePBC,
                                const t_grpopts     *inputrecGroupOptions,
                                warninp_t            wi);

}      // namespace gmx

#endif /* GMX_AWH_READPARAMS_H */
