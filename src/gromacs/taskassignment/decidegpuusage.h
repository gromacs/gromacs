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

#include <vector>

struct gmx_hw_info_t;

enum class EmulateGpuNonbonded : bool;

namespace gmx
{

//! Record where a compute task is targetted.
enum class TaskTarget : int
{
    Auto,
    Cpu,
    Gpu
};

/*! \brief Decide whether this thread-MPI simulation will run
 * nonbonded tasks on GPUs.
 *
 * The number of GPU tasks and devices influences both the choice of
 * the number of ranks, and checks upon any such choice made by the
 * user. So we need to consider this before any automated choice of
 * the number of thread-MPI ranks.
 *
 * \param[in]  nonbondedTarget           The user's choice for mdrun -nb for where to assign short-ranged nonbonded interaction tasks.
 * \param[in]  gpuIdsToUse               The compatible GPUs that the user permitted us to use.
 * \param[in]  userGpuTaskAssignment     The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  emulateGpuNonbonded       Whether we will emulate GPU calculation of nonbonded interactions.
 * \param[in]  usingVerletScheme         Whether the nonbondeds are using the Verlet scheme.
 * \param[in]  nonbondedOnGpuIsUseful    Whether computing nonbonded interactions on a GPU is useful for this calculation.
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 *
 * \returns    Whether the simulation will run nonbonded tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForNonbondedWithThreadMpi(const TaskTarget          nonbondedTarget,
                                                     const std::vector<int>   &gpuIdsToUse,
                                                     const std::vector<int>   &userGpuTaskAssignment,
                                                     const EmulateGpuNonbonded emulateGpuNonbonded,
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
 * \param[in]  pmeTarget                 The user's choice for mdrun -pme for where to assign long-ranged PME nonbonded interaction tasks.
 * \param[in]  gpuIdsToUse               The compatible GPUs that the user permitted us to use.
 * \param[in]  userGpuTaskAssignment     The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  canUseGpuForPme           Whether the form of PME chosen can run on a GPU
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 * \param[in]  numPmeRanksPerSimulation  The number of PME ranks in each simulation.
 *
 * \returns    Whether the simulation will run PME tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForPmeWithThreadMpi(const bool              useGpuForNonbonded,
                                               const TaskTarget        pmeTarget,
                                               const std::vector<int> &gpuIdsToUse,
                                               const std::vector<int> &userGpuTaskAssignment,
                                               const bool              canUseGpuForPme,
                                               const int               numRanksPerSimulation,
                                               const int               numPmeRanksPerSimulation);

/*! \brief Decide whether the simulation will try to run nonbonded
 * tasks on GPUs.
 *
 * The final decision cannot be made until after the duty of the rank
 * is known. But we need to know if nonbonded will run on GPUs for
 * setting up DD (particularly rlist) and determining duty. If the
 * user requires GPUs for the tasks of that duty, then it will be an
 * error when none are found.
 *
 * With thread-MPI, calls have been made to
 * decideWhetherToUseGpusForNonbondedWithThreadMpi() and
 * decideWhetherToUseGpusForPmeWithThreadMpi() to help determine
 * the number of ranks and run some checks, but the final
 * decision is made in this routine, along with many more
 * consistency checks.
 *
 * \param[in]  nonbondedTarget           The user's choice for mdrun -nb for where to assign short-ranged nonbonded interaction tasks.
 * \param[in]  userGpuTaskAssignment     The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  emulateGpuNonbonded       Whether we will emulate GPU calculation of nonbonded interactions.
 * \param[in]  usingVerletScheme         Whether the nonbondeds are using the Verlet scheme.
 * \param[in]  nonbondedOnGpuIsUseful    Whether computing nonbonded interactions on a GPU is useful for this calculation.
 * \param[in]  gpusWereDetected          Whether compatible GPUs were detected on any node.
 *
 * \returns    Whether the simulation will run nonbonded and PME tasks, respectively, on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForNonbonded(const TaskTarget           nonbondedTarget,
                                        const std::vector<int>    &userGpuTaskAssignment,
                                        const EmulateGpuNonbonded  emulateGpuNonbonded,
                                        const bool                 usingVerletScheme,
                                        const bool                 nonbondedOnGpuIsUseful,
                                        const bool                 gpusWereDetected);

/*! \brief Decide whether the simulation will try to run tasks of
 * different types on GPUs.
 *
 * The final decision cannot be made until after the duty of the rank
 * is known. But we need to know if nonbonded will run on GPUs for
 * setting up DD (particularly rlist) and determining duty. If the
 * user requires GPUs for the tasks of that duty, then it will be an
 * error when none are found.
 *
 * With thread-MPI, calls have been made to
 * decideWhetherToUseGpusForNonbondedWithThreadMpi() and
 * decideWhetherToUseGpusForPmeWithThreadMpi() to help determine
 * the number of ranks and run some checks, but the final
 * decision is made in this routine, along with many more
 * consistency checks.
 *
 * \param[in]  useGpuForNonbonded        Whether GPUs will be used for nonbonded interactions.
 * \param[in]  pmeTarget                 The user's choice for mdrun -pme for where to assign long-ranged PME nonbonded interaction tasks.
 * \param[in]  userGpuTaskAssignment     The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  canUseGpuForPme           Whether the form of PME chosen can run on a GPU
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 * \param[in]  numPmeRanksPerSimulation  The number of PME ranks in each simulation.
 * \param[in]  gpusWereDetected          Whether compatible GPUs were detected on any node.
 *
 * \returns    Whether the simulation will run nonbonded and PME tasks, respectively, on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForPme(const bool              useGpuForNonbonded,
                                  const TaskTarget        pmeTarget,
                                  const std::vector<int> &userGpuTaskAssignment,
                                  const bool              canUseGpuForPme,
                                  const int               numRanksPerSimulation,
                                  const int               numPmeRanksPerSimulation,
                                  const bool              gpusWereDetected);

}

#endif
