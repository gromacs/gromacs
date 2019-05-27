/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief SETTLE tests header.
 *
 * Declares the functions that do the buffer management and apply
 * SETTLE constraints ("test runners").
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_SETTLE_RUNNERS_H
#define GMX_MDLIB_TESTS_SETTLE_RUNNERS_H

#include "config.h"

#include "gromacs/math/vectypes.h"

struct t_pbc;
struct gmx_mtop_t;
struct t_idef;
struct t_mdatoms;

namespace gmx
{
namespace test
{

#if GMX_GPU == GMX_GPU_CUDA

/*! \brief
 * Initialize and apply SETTLE constraints on CUDA-enabled GPU.
 *
 * \param[in]     numAtoms          Number of atoms.
 * \param[in]     h_x               Coordinates before timestep (in CPU memory).
 * \param[in,out] h_xp              Coordinates after timestep (in CPU memory). The
 *                                  resulting constrained coordinates will be saved here.
 * \param[in]     updateVelocities  If the velocities should be updated.
 * \param[in,out] h_v               Velocities to update (in CPU memory, can be nullptr
 *                                  if not updated).
 * \param[in]     invdt             Reciprocal timestep (to scale Lagrange
 *                                  multipliers when velocities are updated)
 * \param[in]     computeVirial     If virial should be updated.
 * \param[in,out] virialScaled      Scaled virial tensor to be updated.
 * \param[in]     pbc               Periodic boundary data.
 * \param[in]     mtop              Topology of the system to get the masses for O and
 *                                  H atoms target O-H and H-H distances.
 * \param[in]     idef              System topology.
 * \param[in]     mdatoms           Atoms data.
 */
void applySettleCuda(int                numAtoms,
                     const rvec        *h_x,
                     rvec              *h_xp,
                     bool               updateVelocities,
                     rvec              *h_v,
                     real               invdt,
                     bool               computeVirial,
                     tensor             virialScaled,
                     const t_pbc       *pbc,
                     const gmx_mtop_t  &mtop,
                     const t_idef      &idef,
                     const t_mdatoms   &mdatoms);

#endif // GMX_GPU == GMX_GPU_CUDA

}      // namespace test
}      // namespace gmx

#endif // GMX_MDLIB_TESTS_SETTLE_RUNNERS_H
