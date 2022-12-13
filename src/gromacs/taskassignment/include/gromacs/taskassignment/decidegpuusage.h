/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
struct gmx_mtop_t;
struct t_inputrec;
enum class PmeRunMode;

namespace gmx
{

class MDLogger;

//! Record where a compute task is targetted.
enum class TaskTarget : int
{
    Auto,
    Cpu,
    Gpu
};

//! Help pass GPU-emulation parameters with type safety.
enum class EmulateGpuNonbonded : bool
{
    //! Do not emulate GPUs.
    No,
    //! Do emulate GPUs.
    Yes
};

/*! \libinternal
 *  \brief Structure that holds boolean flags corresponding to the development
 *        features present enabled through environment variables.
 *
 */
struct DevelopmentFeatureFlags
{
    //! True if the Buffer ops development feature is enabled
    // TODO: when the trigger of the buffer ops offload is fully automated this should go away
    bool enableGpuBufferOps = false;
    //! True if the GPU-aware MPI can be used for GPU direct communication feature
    bool canUseGpuAwareMpi = false;
    //! True if GPU PME-decomposition is enabled
    bool enableGpuPmeDecomposition = false;
    //! True if CUDA Graphs are enabled
    bool enableCudaGraphs = false;
};


class MDAtoms;

/*! \brief Decide whether this thread-MPI simulation will run
 * nonbonded tasks on GPUs.
 *
 * The number of GPU tasks and devices influences both the choice of
 * the number of ranks, and checks upon any such choice made by the
 * user. So we need to consider this before any automated choice of
 * the number of thread-MPI ranks.
 *
 * \param[in] nonbondedTarget              The user's choice for mdrun -nb for where to assign
 *                                         short-ranged nonbonded interaction tasks.
 * \param[in] haveAvailableDevices         Whether there are available devices.
 * \param[in] userGpuTaskAssignment        The user-specified assignment of GPU tasks to device IDs.
 * \param[in] emulateGpuNonbonded          Whether we will emulate GPU calculation of nonbonded
 *                                         interactions.
 * \param[in] buildSupportsNonbondedOnGpu  Whether GROMACS was built with GPU support.
 * \param[in] nonbondedOnGpuIsUseful       Whether computing nonbonded interactions on a GPU is
 *                                         useful for this calculation.
 * \param[in] numRanksPerSimulation        The number of ranks in each simulation.
 *
 * \returns    Whether the simulation will run nonbonded tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForNonbondedWithThreadMpi(TaskTarget              nonbondedTarget,
                                                     bool                    haveAvailableDevices,
                                                     const std::vector<int>& userGpuTaskAssignment,
                                                     EmulateGpuNonbonded     emulateGpuNonbonded,
                                                     bool buildSupportsNonbondedOnGpu,
                                                     bool nonbondedOnGpuIsUseful,
                                                     int  numRanksPerSimulation);

/*! \brief Decide whether this thread-MPI simulation will run
 * PME tasks on GPUs.
 *
 * The number of GPU tasks and devices influences both the choice of
 * the number of ranks, and checks upon any such choice made by the
 * user. So we need to consider this before any automated choice of
 * the number of thread-MPI ranks.
 *
 * \param[in]  useGpuForNonbonded        Whether GPUs will be used for nonbonded interactions.
 * \param[in]  pmeTarget                 The user's choice for mdrun -pme for where to assign
 *                                       long-ranged PME nonbonded interaction tasks.
 * \param[in]  pmeFftTarget              The user's choice for mdrun -pmefft for where to run FFT.
 * \param[in]  numDevicesToUse           The number of compatible GPUs that the user permitted us to use.
 * \param[in]  userGpuTaskAssignment     The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  inputrec                  The user input
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 * \param[in]  numPmeRanksPerSimulation  The number of PME ranks in each simulation.
 *
 * \returns    Whether the simulation will run PME tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForPmeWithThreadMpi(bool                    useGpuForNonbonded,
                                               TaskTarget              pmeTarget,
                                               TaskTarget              pmeFftTarget,
                                               int                     numDevicesToUse,
                                               const std::vector<int>& userGpuTaskAssignment,
                                               const t_inputrec&       inputrec,
                                               int                     numRanksPerSimulation,
                                               int                     numPmeRanksPerSimulation);

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
 * \param[in]  nonbondedTarget             The user's choice for mdrun -nb for where to assign short-ranged nonbonded interaction tasks.
 * \param[in]  userGpuTaskAssignment       The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  emulateGpuNonbonded         Whether we will emulate GPU calculation of nonbonded interactions.
 * \param[in]  buildSupportsNonbondedOnGpu Whether GROMACS was build with GPU support.
 * \param[in]  nonbondedOnGpuIsUseful      Whether computing nonbonded interactions on a GPU is useful for this calculation.
 * \param[in]  gpusWereDetected            Whether compatible GPUs were detected on any node.
 *
 * \returns    Whether the simulation will run nonbonded and PME tasks, respectively, on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForNonbonded(TaskTarget              nonbondedTarget,
                                        const std::vector<int>& userGpuTaskAssignment,
                                        EmulateGpuNonbonded     emulateGpuNonbonded,
                                        bool                    buildSupportsNonbondedOnGpu,
                                        bool                    nonbondedOnGpuIsUseful,
                                        bool                    gpusWereDetected);

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
 * \param[in]  pmeFftTarget              The user's choice for mdrun -pmefft for where to do FFT for PME.
 * \param[in]  userGpuTaskAssignment     The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  inputrec                  The user input
 * \param[in]  numRanksPerSimulation     The number of ranks in each simulation.
 * \param[in]  numPmeRanksPerSimulation  The number of PME ranks in each simulation.
 * \param[in]  gpusWereDetected          Whether compatible GPUs were detected on any node.
 *
 * \returns    Whether the simulation will run nonbonded and PME tasks, respectively, on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForPme(bool                    useGpuForNonbonded,
                                  TaskTarget              pmeTarget,
                                  TaskTarget              pmeFftTarget,
                                  const std::vector<int>& userGpuTaskAssignment,
                                  const t_inputrec&       inputrec,
                                  int                     numRanksPerSimulation,
                                  int                     numPmeRanksPerSimulation,
                                  bool                    gpusWereDetected);

/*! \brief Determine PME run mode.
 *
 * Given the PME task assignment in \p useGpuForPme and the user-provided
 * FFT task target in \p pmeFftTarget, returns a PME run mode for the
 * current run. It also checks the compatibility of the two.
 *
 * \note Aborts the run upon incompatible values of \p useGpuForPme and \p pmeFftTarget.
 *
 * \param[in]  useGpuForPme              PME task assignment, true if PME task is mapped to the GPU.
 * \param[in]  pmeFftTarget              The user's choice for -pmefft for where to assign the FFT
 *                                       work of the PME task.
 * \param[in]  inputrec                  The user input record
 * */
PmeRunMode determinePmeRunMode(bool useGpuForPme, const TaskTarget& pmeFftTarget, const t_inputrec& inputrec);

/*! \brief Decide whether the simulation will try to run bonded tasks on GPUs.
 *
 * \param[in]  useGpuForNonbonded        Whether GPUs will be used for nonbonded interactions.
 * \param[in]  useGpuForPme              Whether GPUs will be used for PME interactions.
 * \param[in]  bondedTarget              The user's choice for mdrun -bonded for where to assign tasks.
 * \param[in]  inputrec                  The user input.
 * \param[in]  mtop                      The global topology.
 * \param[in]  numPmeRanksPerSimulation  The number of PME ranks in each simulation, can be -1 for auto.
 * \param[in]  gpusWereDetected          Whether compatible GPUs were detected on any node.
 *
 * \returns    Whether the simulation will run bondeded tasks on GPUs.
 *
 * \throws     std::bad_alloc          If out of memory
 *             InconsistentInputError  If the user requirements are inconsistent. */
bool decideWhetherToUseGpusForBonded(bool              useGpuForNonbonded,
                                     bool              useGpuForPme,
                                     TaskTarget        bondedTarget,
                                     const t_inputrec& inputrec,
                                     const gmx_mtop_t& mtop,
                                     int               numPmeRanksPerSimulation,
                                     bool              gpusWereDetected);

/*! \brief Decide whether to use GPU for update.
 *
 * \param[in]  isDomainDecomposition        Whether there more than one domain.
 * \param[in]  useUpdateGroups              If the constraints can be split across domains.
 * \param[in]  pmeRunMode                   PME running mode: CPU, GPU or mixed.
 * \param[in]  havePmeOnlyRank              If there is a PME-only rank in the simulation.
 * \param[in]  useGpuForNonbonded           Whether GPUs will be used for nonbonded interactions.
 * \param[in]  updateTarget                 User choice for running simulation on GPU.
 * \param[in]  gpusWereDetected             Whether compatible GPUs were detected on any node.
 * \param[in]  inputrec                     The user input.
 * \param[in]  mtop                         The global topology.
 * \param[in]  useEssentialDynamics         If essential dynamics is active.
 * \param[in]  doOrientationRestraints      If orientation restraints are enabled.
 * \param[in]  haveFrozenAtoms              If this simulation has frozen atoms (see Issue #3920).
 * \param[in]  useModularSimulator          Whether the modular simulator is used
 * \param[in]  doRerun                      It this is a rerun.
 * \param[in]  mdlog                        MD logger.
 *
 * \returns    Whether complete simulation can be run on GPU.
 * \throws     std::bad_alloc            If out of memory
 *             InconsistentInputError    If the user requirements are inconsistent.
 */
bool decideWhetherToUseGpuForUpdate(bool                 isDomainDecomposition,
                                    bool                 useUpdateGroups,
                                    PmeRunMode           pmeRunMode,
                                    bool                 havePmeOnlyRank,
                                    bool                 useGpuForNonbonded,
                                    TaskTarget           updateTarget,
                                    bool                 gpusWereDetected,
                                    const t_inputrec&    inputrec,
                                    const gmx_mtop_t&    mtop,
                                    bool                 useEssentialDynamics,
                                    bool                 doOrientationRestraints,
                                    bool                 haveFrozenAtoms,
                                    bool                 useModularSimulator,
                                    bool                 doRerun,
                                    const gmx::MDLogger& mdlog);

/*! \brief Decide whether direct GPU communication can be used.
 *
 * Takes into account the build type which determines feature support as well as GPU
 * development feature flags, determines whether this run can use direct GPU communication.
 * The final decision whether the run will use direct communication for either of the features
 * which rely on it is made during task assignment / simulationWorkload initialization.
 *
 * \param[in]  devFlags                     GPU development / experimental feature flags.
 * \param[in]  haveMts                      Whether the simulation uses multiple time stepping
 * \param[in]  haveSwapCoords               Whether the swap-coords functionality is active
 * \param[in]  mdlog                        MD logger.
 *
 * \returns    Whether the MPI-parallel runs can use direct GPU communication.
 */
bool decideWhetherDirectGpuCommunicationCanBeUsed(const DevelopmentFeatureFlags& devFlags,
                                                  bool                           haveMts,
                                                  bool                           haveSwapCoords,
                                                  const gmx::MDLogger&           mdlog);

/*! \brief Decide whether to use GPU for halo exchange.
 *
 * \param[in]  havePPDomainDecomposition    Whether PP domain decomposition is in use.
 * \param[in]  useGpuForNonbonded           Whether GPUs will be used for nonbonded interactions.
 * \param[in]  canUseDirectGpuComm          Whether direct GPU communication can be used.
 * \param[in]  useModularSimulator          Whether modularsimulator is in use.
 * \param[in]  doRerun                      Whether this is a rerun.
 * \param[in]  haveEnergyMinimization       Whether energy minimization is in use.
 * \param[in]  mdlog                        MD logger.
 *
 * \returns    Whether halo exchange can be run on GPU.
 */
bool decideWhetherToUseGpuForHalo(bool                 havePPDomainDecomposition,
                                  bool                 useGpuForNonbonded,
                                  bool                 canUseDirectGpuComm,
                                  bool                 useModularSimulator,
                                  bool                 doRerun,
                                  bool                 haveEnergyMinimization,
                                  const gmx::MDLogger& mdlog);

} // namespace gmx

#endif
