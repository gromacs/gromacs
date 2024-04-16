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
/*! \defgroup module_taskassignment Assigning simulation tasks to hardware (taskassignment)
 * \ingroup group_mdrun
 * \brief Provides code that manages assignment of simulation tasks to hardware.
 */
/*! \libinternal
 * \file
 * \brief Declares high-level functionality for managing assigning
 * tasks on ranks of a node to hardware on that node, and the factory
 * function to build the correct flavours of gmx::INodeTaskAssigner
 * required to implement the user's requirements.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 * \inlibraryapi
 */
#ifndef GMX_TASKASSIGNMENT_TASKASSIGNMENT_H
#define GMX_TASKASSIGNMENT_TASKASSIGNMENT_H

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"

struct DeviceInformation;
struct gmx_hw_info_t;
struct t_commrec;

enum class PmeRunMode;

namespace gmx
{

enum class TaskTarget;
class MDLogger;
class PhysicalNodeCommunicator;
class SimulationWorkload;

/*! \brief Types of compute tasks that can be run on a GPU.
 *
 * These names refer to existing practice in GROMACS, which is not
 * strictly accurate. */
enum class GpuTask : int
{
    //! Short-ranged interactions.
    Nonbonded,
    //! Long-ranged interactions.
    Pme,
    //! Number of possible tasks.
    Count
};

/*! \libinternal
 * \brief Specifies the GPU deviceID_ available for task_ to use. */
struct GpuTaskMapping
{
    //! The type of this GPU task.
    GpuTask task_;
    //! Device ID on this node to which this GPU task is mapped.
    int deviceId_;
};

//! Container of GPU tasks on a rank, specifying the task type and GPU device ID, e.g. potentially ready for consumption by the modules on that rank.
using GpuTaskAssignment = std::vector<GpuTaskMapping>;

class GpuTaskAssignments;

/*! \libinternal
 * \brief Builder for the GpuTaskAssignments for all ranks on this
 * node.
 *
 * This will coordinate the final stages of task assignment and
 * reporting, and build the GpuTaskAssignments object used to
 * configure the modules that might run tasks on GPUs.
 *
 * Communicates between ranks on a node to coordinate task assignment
 * between them onto available hardware, e.g. accelerators.
 *
 * \todo Later, this might become a loop over all registered modules
 * relevant to the mdp inputs, to find those that have such tasks.
 *
 * \todo Later we might need the concept of computeTasksOnThisRank,
 * from which we construct gpuTasksOnThisRank.
 *
 * Currently the DD code assigns duty to ranks that can
 * include PP work that currently can be executed on a single
 * GPU, if present and compatible.  This has to be coordinated
 * across PP ranks on a node, with possible multiple devices
 * or sharing devices on a node, either from the user
 * selection, or automatically. */
class GpuTaskAssignmentsBuilder
{
public:
    //! Constructor
    GpuTaskAssignmentsBuilder();

    /*! \brief Builds a GpuTaskAssignments
     *
     * This method reconciles
     *
     *   - user mdrun command-line options,
     *   - the results of hardware detection
     *   - the duty assigned by the DD setup,
     *   - the requested simulation modules, and
     *   - the possible existence of multi-simulations
     *
     * to assign the GPUs on each physical node to the tasks on
     * the ranks of that node. It throws InconsistentInputError
     * when a/the useful GPU task assignment is not possible.
     *
     * \param[in]  availableDevices       The compatible devices that the user permitted us to use.
     * \param[in]  userGpuTaskAssignment  The user-specified assignment of GPU tasks to device IDs.
     * \param[in]  hardwareInfo           The detected hardware
     * \param[in]  gromacsWorldComm       MPI communicator for all ranks in the current GROMACS run
     * \param[in]  physicalNodeComm       Communication object for this physical node.
     * \param[in]  nonbondedTarget        The user's choice for mdrun -nb for where to assign
     *                                    short-ranged nonbonded interaction tasks.
     * \param[in]  pmeTarget              The user's choice for mdrun -pme for where to assign
     *                                    long-ranged PME nonbonded interaction tasks.
     * \param[in]  bondedTarget           The user's choice for mdrun -bonded for where to assign tasks.
     * \param[in]  updateTarget           The user's choice for mdrun -update for where to assign tasks.
     * \param[in]  useGpuForNonbonded     Whether GPUs will be used for nonbonded interactions.
     * \param[in]  useGpuForPme           Whether GPUs will be used for PME interactions.
     * \param[in]  rankHasPpTask          Whether this rank has a PP task
     * \param[in]  rankHasPmeTask         Whether this rank has a PME task
     *
     * \throws   std::bad_alloc          If out of memory.
     *           InconsistentInputError  If user and/or detected inputs are inconsistent.
     */
    static GpuTaskAssignments build(ArrayRef<const int>             availableDevices,
                                    ArrayRef<const int>             userGpuTaskAssignment,
                                    const gmx_hw_info_t&            hardwareInfo,
                                    MPI_Comm                        gromacsWorldComm,
                                    const PhysicalNodeCommunicator& physicalNodeComm,
                                    TaskTarget                      nonbondedTarget,
                                    TaskTarget                      pmeTarget,
                                    TaskTarget                      bondedTarget,
                                    TaskTarget                      updateTarget,
                                    bool                            useGpuForNonbonded,
                                    bool                            useGpuForPme,
                                    bool                            rankHasPpTask,
                                    bool                            rankHasPmeTask);
};

/*! \libinternal
 * \brief Contains the GPU task assignment for all ranks on this
 * physical node.
 *
 * This can be used to configure the modules on this rank that might
 * run tasks on GPUs.
 *
 * This assignment is made by a GpuTaskAssignmentsBuilder object. */
class GpuTaskAssignments
{
public:
    //! Public move constructor to use with the builder
    GpuTaskAssignments(GpuTaskAssignments&& source) noexcept = default;

private:
    // Let the builder handle construction
    friend class GpuTaskAssignmentsBuilder;
    //! Private constructor so only the builder can construct
    GpuTaskAssignments(const gmx_hw_info_t& hardwareInfo);
    /*! \brief Information about hardware on this physical node
     *
     * The lifetime of the object referred to must exceed that
     * of this object. */
    const gmx_hw_info_t& hardwareInfo_;
    //! The GPU task assignment for all ranks on this node
    std::vector<GpuTaskAssignment> assignmentForAllRanksOnThisNode_;
    /*! \brief The index of this rank within those on this node.
     *
     * This is useful for indexing into \c
     * assignmentForAllRanksOnThisNode_. */
    Index indexOfThisRank_ = -1;
    //! Number of GPU tasks on this node.
    size_t numGpuTasksOnThisNode_ = 0;
    //! Number of ranks on this physical node.
    size_t numRanksOnThisNode_ = 0;

    //! Vector of device IDs assigned to this node
    std::vector<int> deviceIdsAssigned_;

public:
    /*! \brief Log a report on how GPUs are being used on
     * the ranks of the physical node of rank 0 of the simulation.
     *
     * \todo It could be useful to report also whether any nodes differed,
     * and in what way.
     *
     * \param[in]  mdlog           Logging object.
     * \param[in]  printHostName   Print the hostname in the usage information.
     * \param[in]  pmeRunMode      Describes the execution of PME tasks.
     * \param[in]  simulationWork  Simulation workload descriptor
     *
     * \throws     std::bad_alloc if out of memory
     */
    void reportGpuUsage(const MDLogger&           mdlog,
                        bool                      printHostName,
                        PmeRunMode                pmeRunMode,
                        const SimulationWorkload& simulationWork);

    /*! \brief Logs to \c mdlog information that may help a user
     * learn how to let mdrun make a task assignment that runs
     * faster.
     *
     * \param[in]  mdlog                         Logging object.
     * \param[in]  numAvailableDevicesOnThisNode The number of compatible devices on this node
     *                                           that the user permitted us to use.
     * */
    void logPerformanceHints(const MDLogger& mdlog, size_t numAvailableDevicesOnThisNode);
    /*! \brief Return handle to the initialized GPU to use in this rank.
     *
     * \param[out] deviceId Index of the assigned device.
     *
     * \returns Device information on the selected device. Returns nullptr if no GPU task
     *          is assigned to this rank.
     */
    DeviceInformation* initDevice(int* deviceId) const;
    //! Return whether this rank has a PME task running on a GPU
    bool thisRankHasPmeGpuTask() const;
    //! Return whether this rank has any task running on a GPU
    bool thisRankHasAnyGpuTask() const;
    //! Get the list of unique devices that have been assigned tasks on this physical node
    std::vector<int> deviceIdsAssigned() { return deviceIdsAssigned_; }
};

} // namespace gmx

#endif
