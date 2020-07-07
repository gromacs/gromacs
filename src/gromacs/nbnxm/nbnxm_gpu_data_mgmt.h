/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 *  \brief Declare common functions for NBNXM GPU data management.
 *
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 *  \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_NBNXM_GPU_DATA_MGMT_H
#define GMX_NBNXM_NBNXM_GPU_DATA_MGMT_H

struct interaction_const_t;
struct NBParamGpu;
struct PairlistParams;

namespace Nbnxm
{

struct gpu_plist;

/*! \brief Tabulates the Ewald Coulomb force and initializes the size/scale and the table GPU array.
 *
 * If called with an already allocated table, it just re-uploads the
 * table.
 */
void init_ewald_coulomb_force_table(const EwaldCorrectionTables& tables,
                                    NBParamGpu*                  nbp,
                                    const DeviceContext&         deviceContext);

/*! \brief Selects the Ewald kernel type, analytical or tabulated, single or twin cut-off. */
int nbnxn_gpu_pick_ewald_kernel_type(const interaction_const_t gmx_unused& ic);

/*! \brief Copies all parameters related to the cut-off from ic to nbp
 */
void set_cutoff_parameters(NBParamGpu* nbp, const interaction_const_t* ic, const PairlistParams& listParams);

/*! \brief Initializes the pair list data structure.
 */
void init_plist(gpu_plist* pl);

/*! \brief Initializes the timings data structure. */
void init_timings(gmx_wallclock_gpu_nbnxn_t* t);

} // namespace Nbnxm

#endif // GMX_NBNXM_NBNXM_GPU_DATA_MGMT_H
