/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * \brief
 * Declares gmx::HardwareTopology
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inlibraryapi
 * \ingroup module_hardware
 */
#ifndef GMX_HARDWARE_HARDWARETOPOLOGY_H
#define GMX_HARDWARE_HARDWARETOPOLOGY_H

#include <cmath>
#include <cstddef>
#include <cstdint>

#include <array>
#include <map>
#include <string>
#include <vector>

namespace gmx
{

/*! \libinternal \brief Information about packages, cores, processing units, numa, caches
 *
 * This class is the main GROMACS interface to provide information about the
 * hardware of the system we are running on. Internally, it uses either
 * hwloc for full or almost-full information, or a fallback implementation
 * that relies on Linux sysinfo.
 *
 * You should always use this class to query the hardware layout in user code.
 * Note that you cannot rely on any information being present, but you must check
 * with the supportLevel() method before trying to access any information.
 */
class HardwareTopology
{
public:
    /*! \brief Amount of topology information present (incremental)
     *
     * For the LogicalProcessorCount alternative, the value of maxThreads()
     * will primarily reflect the allowed cpuLimit on the machine so we don't
     * overload it. If we could not find that, we will have tried to make it
     * reflect only the cpus on which we are allowed to run. Failing that too,
     * it will correspond to the total number of logical cpus.
     */
    enum class SupportLevel
    {
        None,                  //!< No hardware information whatsoever. Sorry.
        LogicalProcessorCount, //!< Only total processor count
        Basic,                 //!< Package, core and processing unit for allowed processors
        Full,                  //!< Cache, memory and numa node info
        FullWithDevices        //!< Information about devices on the PCI bus
    };

    /*! \libinternal \brief Information about a single cache level */
    struct Cache
    {
        int         level;         //!< Level relative to core (starts at 1)
        std::size_t size;          //!< size in bytes, 0 if unknown
        int         linesize;      //!< size of each cache line in bytes, 0 if unknown
        int         associativity; //!< associativity, -1 means fully associative
        int         shared;        //!< Number of processing units sharing this cache
    };

    /*! \libinternal \brief Information about a single processing unit (in a core)
     *
     * Gromacs now follows hwloc's nomenclature of calling the lowest-level
     * entity in the system "processing unit", which is typically a hardware thread,
     * while we call the (physical) grouping of many cores a package rather than
     * socket.
     * This avoids the confusion between "logical processor", "core", "thread",
     * and the indices assigned by the operating system, and it resolves the
     * situation where a single socket might host multiple chip modules.
     *
     * The id of the processing unit is specific only to this hardware topology.
     * The processing unit id typically increases continuously as you walk
     * through packages and cores in order of their ids. In general it
     * will be different from the processor id assigned by the operating
     * system. To achieve better load balancing when starting fewer threads than
     * hardware threads on SMT systems, Linux e.g. typically assigns processor ids
     * in a round-robin fashion over all cores.
     * In addition, the hardware topology now only contains the packages,
     * cores, and processing units on which we are permitted to run,
     * which in general means it will not be uniform.
     * This also means all the resources present in the entire (physical) node
     * will in general not be visible in the hardware topology structure, so you
     * cannot e.g. use the number of cores in the first package to check the
     * number of cores per package.
     *
     * In practice, even hwloc appears to be a bit inconsistent, and when testing
     * with a subset of cores the topology included empty packages, while
     * core Ids were changed so only permitted ones were included.
     */
    struct ProcessingUnit
    {
        int id;   //!<  index of this processing unit in hardware topology
        int osId; //!<  index assigned by the operating system
    };

    /*! \libinternal \brief Information about a single core in a package */
    struct Core
    {
        int                         id;              //!<  index of this core in hardware topology
        int                         numaNodeId;      //!< id of the numa node of this core
        std::vector<ProcessingUnit> processingUnits; //!< PUs in this core on which we can run
    };

    /*! \libinternal \brief Information about a single package in the system.
     *
     * \note The list of cores in the package might be empty if none of them contains
     * processing units on which we are allowed to run.
     */
    struct Package
    {
        int               id;    //!<  id of this package in hardware topology/system
        std::vector<Core> cores; //!< Cores in this package with processing units on which we can run.
    };

    /*! \libinternal \brief Information about each numa node in system
     *
     * \note The list of processing units might be empty if we are not allowed to
     * run on any in this numa node.
     */
    struct NumaNode
    {
        int              id;              //!<  id of numa node in hardware topology
        std::size_t      memory;          //!< Total detected memory in bytes
        std::vector<int> processingUnits; //!< PUs in this node on which we can run.
    };

    /*! \libinternal \brief Information about all numa nodes.
     *
     * \note We are usually allowed to use memory on all numa nodes even if we are restricted
     * to processing units on a single numa node. In this case one or more of the numa nodes might
     * have an empty list of processing unit ids on which we are allowed to run.
     */
    struct Numa
    {
        std::vector<NumaNode>           nodes;       //!< Information about each numa node
        float                           baseLatency; //!< Scale factor for relative latencies
        std::vector<std::vector<float>> relativeLatency; //!< 2D matrix of relative latencies between nodes
        float                           maxRelativeLatency; //!< Largest relative latency
    };

    /*! \libinternal \brief Information about a single PCI device.
     *
     *  \note On many systems the PCI bus is not directly connected to any numa node.
     *        For these systems the numaNodeId will be -1, so you cannot rely on this
     *        number reflecting a specific numa node.
     */
    struct Device
    {
        std::uint16_t vendorId;   //!< Vendor identification
        std::uint16_t deviceId;   //!< Vendor-specific device identification
        std::uint16_t classId;    //!< class (high 8 bits) and subclass (low 8 bits)
        std::uint16_t domain;     //!< Domain, usually 0 for PCI bus
        std::uint8_t  bus;        //!< Bus number in domain
        std::uint8_t  dev;        //!< Device on bus
        std::uint8_t  func;       //!< Function id for multi-function devices
        int           numaNodeId; //!< Numa node, -1 if the bus is not located inside a node
    };

    /*! \libinternal \brief Information about package, core and processing unit for a logical processor
     *
     * \note In general, not all logical processors in the system will be visible, so the osId
     * field can be larger than the total number of processing units in the topology.
     *
     * The indices of packages, cores and PUs correspond to the  indices in this hardware
     * topology, rather than the global/physical rank in the system, meaning they are consistent
     * with the contents of the vectors in the machine.packags[].cores[].processingUnits[] hierarchy.
     */
    struct LogicalProcessor
    {
        int puId;                     //!< Index of PU in hardware topology
        int osId;                     //!<  index assigned by the operating system
        int packageRankInTopology;    //!<  index of package in machine
        int coreRankInPackage;        //!<  index of core in package
        int processingUnitRankInCore; //!<  index of processing unit in core
        int numaNodeId;               //!< Index of numa node
    };

    /*! \libinternal \brief Hardware topology information about the entire machine
     *
     * The machine structure is a tree with top-down information about all
     * packages, cores, and processing units in the system. For example, an
     * operating system logical processor index can be found as
     * machine.packages[0].cores[2].processingUnits[1].osId.
     * In some cases you might need the opposite lookup, i.e. the processing
     * unit Id on the topology from the OS-assigned Id. This is present
     * in the map osIdToPuId, provided we have sufficient support level.
     *
     * It can happen that our process is only allowed
     * to run on a subset of processors, e.g. because we are running in a
     * container or because some cores are reserved for management,
     * which complicates core detection a bit. In this case, the hardware
     * topology will only list the PUs on which we are allowed to run,
     * meaning it might look quite unbalanced.
     *
     * First, remember that you must consider the supportLevel before using data in this class.
     * - If supportLevel is SupportLevel::None, we simply don't have any
     *  valid information - not even the number of logical processors.
     * - If we have at least SupportLevel::LogicalProcessorCount, we have
     *   detected the number of logical CPUs  in the system, but we might not have any
     *   information whether we are allowed to run on all of them or not.
     *   To decide how many threads to start, always check the maxThreads()
     *   method, since we might have been able to detect that we are only allowed to
     *   run on a subset, e.g. when executing in a container environment where
     *  the machine has many cores, but the cpu limit for our container is low.
     * - If supportLevel is at least SupportLevel::Basic, we have
     *  valid topology information about packages, cores and processing units both
     *  in the tree-like structure starting with machine().packages, in
     *  the linear logicalProcessors vector, and the map osIdToPuId is also valid.
     * - If supportLevel is at least SupportLevel::Full, the cache and
     *  ccNUMA information is also valid.
     * - Finally, for SupportLevel::FullWithDevices, the device structure is
     *  valid in addition to everything else.
     *
     * \note Never use the size of the logicalProcessors vector to decide
     *      how many threads to start, but consult maxThreads().
     * \note If you use cpu sets or other ways to decide what processors to run on,
     *       such OS-provided functions always refer to the osId of the processing unit.
     */
    struct Machine
    {
        Machine();

        std::vector<LogicalProcessor> logicalProcessors; //!< Logical processing unit info, indexed by PU Ids in vector.
        std::map<int, int>   osIdToPuId; //!< Map from OS to topology processing unit ids
        std::vector<Package> packages;   //!< All the packages in the system
        std::vector<Cache>   caches;     //!< Caches in increasing level order
        Numa                 numa;       //!< Structure with all numa information
        std::vector<Device>  devices;    //!< Devices on PCI bus
    };

    /*! \brief Detects the hardware topology. */
    static HardwareTopology detect();

    /*! \brief Creates mock topology with given number of logical cores.
     *
     * The support level will be None if the argument is 0 or smaller, otherwise LogicalProcessorCount.
     *
     * Intended for testing of code that uses the hardware topology.
     */
    explicit HardwareTopology(int logicalProcessorCount);

    /*! \brief Creates mock topology based on APIC (or similar) CPU indices
     *
     * This routine assembles a fake hardware topology based on a vector
     * containing indices describing each logical processor, a second
     * vector describing what processors we should be allowed to run on,
     * and a path to a (fake) filesystem optionally containing cgroup
     * information to detect the recommended maximum load.
     *
     * Intended for testing of code that uses the hardware topology.
     *
     * \param logicalProcessorIdMap Each key in this map is the (OS-provided) index
     *                              of a logical processor on which we are allowed
     *                              to run. The value is an array with three integers
     *                              describing (1) the packageId in the machine,
     *                              (2) the coreId in the package, and (3) the hardware
     *                              thread or "processing unit" id in the core.
     *                              Note that these indices do NOT have to be ranks, only
     *                              unique. In fact, it is common e.g. on x86 that the
     *                              low-level core indices are based on connectivity, so
     *                              the core ids in a 12-core CPU might be 0-5 and 8-13.
     *                              For x86 systems, you can extract an initializer list
     *                              for this parameter by running the standalone version of our
     *                              cpuinfo tool with the hidden/debug "-topology" option.
     * \param filesystemRoot        Path to (fake) filesystem where we attempt to parse
     *                              the allowed cpu load - all file paths mentioned should
     *                              be relative to the provided root. We will first try
     *                              to find if cgroups2 or cgroups1 is mounted by checking
     *                              /etc/mtab. The cgroup of the process will be parsed
     *                              from /proc/self/cgroup, and we also always add the
     *                              blank top-level cgroup ("/") if the specific cgroup
     *                              is not found in the next step. For this cgroup,
     *                              we parse the mount path specified in /etc/mtab. For
     *                              cgroups2, the cpu limit is specified in cpu.max, but
     *                              note that this file can be absent if no limit is set.
     *                              This file contains two numbers where the first is the
     *                              quota of this process, and the second the period, both
     *                              usually specified as microseconds. Note that the quota can
     *                              be larger than the period if we allow loads above 1.0.
     *                              For cgroups1, we first locate the cgroups subgroup,
     *                              but then instead look in the subdirectory "cpu,cpuacct"
     *                              and get the quota from cpu.cfs_quota_us and period
     *                              from cpu.cfs_period_us. The tests directory in the
     *                              hardware module contains a simple script that can
     *                              capture these files from a Linux system.
     */
    explicit HardwareTopology(const std::map<int, std::array<int, 3>>& logicalProcessorIdMap,
                              const std::string&                       filesystemRoot);

    /*! \brief Creates mock topology by parsing mock Linux sys/fs path
     *
     * Create mock hardware topology by attempting to parse processor
     * information from mock Linux sys/fs path.
     *
     * Intended for testing of code that uses the hardware topology.
     *
     * \param filesystemRoot      Path to (fake) filesystem where we will first find all
     *                            logical cpus from /sys/devices/system/cpu/possible,
     *                            after which the topology indices for processor XX are
     *                            read from the directory /sys/devices/system/cpu/cpuXX/topology .
     *                            The package id is read from the file physical_package_id,
     *                            the core id in the package from core_id, and then we
     *                            assume the hardware thread/processing unit ids are assigned
     *                            in the enumeration order of the logical processors.
     *                            After this, we also look for cpu load limits specified
     *                            with cgroups, as described in the other constructor above.
     *                            The tests directory in the hardware module contains a simple
     *                            script that can capture these files from a Linux system.
     * \param allowedProcessors   Vector containing the logical (OS) processor indices
     *                            that should be retained in the topology, mocking the
     *                            logical processors that are enabled in our cpu mask.
     */
    explicit HardwareTopology(const std::string& filesystemRoot, const std::vector<int>& allowedProcessors);

    /*! \brief Check what topology information is available and valid
     *
     *  The amount of hardware topology information that can be detected depends
     *  on both the hardware and whether GROMACS was linked with the external
     *  hwloc library.
     *
     * - If supportLevel is SupportLevel::None, we simply don't have any
     *  valid information - not even the number of logical processors.
     * - If we have at least SupportLevel::LogicalProcessorCount, we have
     *   detected the number of logical CPUs  in the system, but we might not have any
     *   information whether we are allowed to run on all of them or not.
     *   To decide how many threads to start, always use the maxThreads()
     *   method, since we might have been able to detect that we are only allowed to
     *   run on a subset, e.g. when executing in a container environment where
     *  the machine has many cores, but the cpu limit for our container is low.
     * - If supportLevel is at least SupportLevel::Basic, we have
     *  valid topology information about packages, cores and processing units both
     *  in the tree-like structure starting with machine().packages, in
     *  the linear logicalProcessors vector, and the map osIdToPuId is also valid.
     * - If supportLevel is at least SupportLevel::Full, the cache and
     *  ccNUMA information is also valid.
     * - Finally, for SupportLevel::FullWithDevices, the device structure is
     *  valid in addition to everything else.
     */
    SupportLevel supportLevel() const { return supportLevel_; }

    /*! \brief Return true if we actually detected hardware.
     *
     *  \return This method will normally return true, when we actually ran
     *          the hardware detection as part of this process to construct
     *          the object. It will be false when the object was constructed
     *          by reading a cached XML file, or possibly generated from
     *          synthetic data.
     */
    bool isThisSystem() const { return isThisSystem_; }

    /*! \brief Return the machine topology tree
     *
     *  You can always call this routine, but be aware that some or all contents
     *  will not be valid unless supportLevel() returns a sufficient level.
     *
     *  - With SupportLevel::LogicalProcessorCount, only the field
     *    machine.logicalProcessorCount is valid.
     *  - With SupportLevel::Basic, you can access the vectors of packages,
     *    cores, and hardware threads, and query what logical processorId
     *    each hardware thread corresponds to.
     *  - SupportLevel::Full adds cache, memory and ccNUMA information.
     *  - SupportLevel::FullWithDevices also adds the PCI express bus.
     *
     *  While data that is not valid has been initialized to special values,
     *  you should not rely on those but query the supportLevel() method before
     *  accessing it.
     */
    const Machine& machine() const { return machine_; }

    /*! \brief Practical max cpu load, as limited by the OS
     *
     * In some cases, in particular when running in containers, the total
     * number of logical processors on which we re allowed to execute can be
     * quite large, while there is a relatively low limit on the amount of CPU time
     * the process can consume. In this case it will be much better to limit ourselves
     * based on the amount of CPU we can use to improve scaling and avoid extra I/O.
     *
     * You can always call this routine, but if sufficient support is not
     * available, it may just return 0.0.
     */
    float cpuLimit() const { return cpuLimit_; }

    /*! \brief Recommended max number of active threads
     *
     * This method provides a recommendation for how many active cpu-consuming
     * threads to start while taking cpuLimit into account.
     * It is not a hard limit on threads, and since it is based on limits on the cpu load
     * there could be cases where it is wise to start additional threads e.g. for handling I/O.
     *
     *  \return Recommended number of threads to start, or 0 if we could not detect.
     */
    int maxThreads() const { return maxThreads_; }

private:
    HardwareTopology();

    SupportLevel supportLevel_; //!< Available topology information
    Machine      machine_;      //!< The machine map
    bool         isThisSystem_; //!< Machine map is real (vs. cached/synthetic)
    float        cpuLimit_;     //!< Max practical load as limited by OS
    int          maxThreads_;   //!< Recommended max # threads
};

} // namespace gmx

#endif // GMX_HARDWARE_HARDWARETOPOLOGY_H
