/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Declares functionality for deciding whether tasks will run on GPUs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 * \inlibraryapi
 */

#ifndef GMX_TASKASSIGNMENT_DECIDEGPUUSAGE_H
#define GMX_TASKASSIGNMENT_DECIDEGPUUSAGE_H

#include <utility>
#include <vector>

struct gmx_hw_info_t;

enum class EmulateGpuNonbonded : bool;

namespace gmx
{

/*! \brief Decide whether this thread-MPI simulation will run
 * nonbonded tasks on GPUs.
 *
 * The number of GPU tasks and devices influences both the choice of
 * the number of ranks, and checks upon any such choice made by the
 * user. So we need to consider this before any automated choice of
 * the number of thread-MPI ranks.
 *
 * \param[in]  nbOptionString            The user's choice for mdrun -nb for where to assign short-ranged nonbonded interaction tasks.
 * \param[in]  gpuIds                    The user-supplied GPU IDs.
 * \param[in]  emulateGpuNonbonded       Whether we will emulate GPU calculation of nonbonded interactions.
 * \param[in]  hardwareInfo              The detected hardware
 * \param[in]  usingVerletScheme         Whether the nonbondeds are using the Verlet scheme.
 * \param[in]  nonbondedOnGpuIsUseful    Whether computing nonbonded interactions on a GPU is useful for this calculation.
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 *
 * \returns    Whether the simulation will run nonbonded tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForNonbondedWithThreadMpi(const char               *nbOptionString,
                                                     const std::vector<int>   &gpuIds,
                                                     const EmulateGpuNonbonded emulateGpuNonbonded,
                                                     const gmx_hw_info_t      &hardwareInfo,
                                                     const bool                usingVerletScheme,
                                                     const bool                nonbondedOnGpuIsUseful,
                                                     const int                 numRanksPerSimulation);

/*! \brief Decide whether this thread-MPI simulation will run
 * PME tasks on GPUs.
 *
 * The number of GPU tasks and devices influences both the choice of
 * the number of ranks, and checks upon any such choice made by the
 * user. So we need to consider this before any automated choice of
 * the number of thread-MPI ranks.
 *
 * \param[in]  useGpuForNonbonded        Whether GPUs will be used for nonbonded interactions.
 * \param[in]  pmeOptionString           The user's choice for mdrun -pme for where to assign long-ranged PME nonbonded interaction tasks.
 * \param[in]  pmeFftOptionString        The user's choice for mdrun -pmefft for where to assign PME solve tasks.
 * \param[in]  gpuIds                    The user-supplied GPU IDs.
 * \param[in]  hardwareInfo              The detected hardware
 * \param[in]  usingAnyPme               Whether any form of PME is in use.
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 *
 * \returns    Whether the simulation will run PME tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForPmeWithThreadMpi(bool                    useGpuForNonbonded,
                                               const char             *pmeOptionString,
                                               const char             *pmeFftOptionString,
                                               const std::vector<int> &gpuIds,
                                               const gmx_hw_info_t    &hardwareInfo,
                                               // TODO should this be a call into the PME module instead? Should this arg appear in the function below?
                                               const bool              usingAnyPme,
                                               const int               numRanksPerSimulation);

/*! \brief Decide whether the simulation will run tasks of different types on GPUs.
 *
 * With thread-MPI, calls have been made to
 * decideWhetherToUseGpusForNonbondedWithThreadMpi() and
 * decideWhetherToUseGpusForPmeWithThreadMpi() to help determine
 * the number of ranks and run some checks, but the final
 * decision is made in this routine, along with many more
 * consistency checks.
 *
 * \param[in]  nbOptionString            The user's choice for mdrun -nb for where to assign short-ranged nonbonded interaction tasks.
 * \param[in]  pmeOptionString           The user's choice for mdrun -pme for where to assign long-ranged PME nonbonded interaction tasks.
 * \param[in]  gpuIds                    The user-supplied GPU IDs.
 * \param[in]  emulateGpuNonbonded       Whether we will emulate GPU calculation of nonbonded interactions.
 * \param[in]  hardwareInfo              The detected hardware
 * \param[in]  usingVerletScheme         Whether the nonbondeds are using the Verlet scheme.
 * \param[in]  nonbondedOnGpuIsUseful    Whether computing nonbonded interactions on a GPU is useful for this calculation.
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 *
 * \returns    Whether the simulation will run nonbonded and PME tasks, respectively, on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
std::pair<bool, bool> decideWhetherToUseGpus(const char                *nbOptionString,
                                             const char                *pmeOptionString,
                                             const std::vector<int>    &gpuIds,
                                             const EmulateGpuNonbonded  emulateGpuNonbonded,
                                             const gmx_hw_info_t       &hardwareInfo,
                                             const bool                 usingVerletScheme,
                                             const bool                 nonbondedOnGpuIsUseful,
                                             const int                  numRanksPerSimulation);

}

#endif
