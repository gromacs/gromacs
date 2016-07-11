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

#include <cstdint>

#include <vector>

namespace gmx
{

/*! \libinternal \brief Information about sockets, cores, threads, numa, caches
 *
 * This class is the main GROMACS interface to provide information about the
 * hardware of the system we are running on. Internally, it uses either
 * hwloc for full or almost-full information, or a fallback implementation
 * that relies on CpuInfo on x86.
 *
 * You should always use this class (rather than CpuInfo directly) to query
 * the hardware layout in user code. Note that you cannot rely on any
 * information being present, but you must check with the supportLevel()
 * method before trying to access any information.
 */
class HardwareTopology
{
    public:

        /*! \brief Amount of topology information present (incremental) */
        enum class SupportLevel
        {
            None,                  //!< No hardware information whatsoever. Sorry.
            LogicalProcessorCount, //!< Only machine().logicalProcessorCount is valid
            Basic,                 //!< Socket, core and hardware thread info
            Full,                  //!< Cache, memory and numa node info
            FullWithDevices        //!< Information about devices on the PCI bus
        };

        /*! \libinternal \brief Information about a single cache level */
        struct Cache
        {
            int                    level;          //!< Level relative to core (starts at 1)
            std::size_t            size;           //!< size in bytes, 0 if unknown
            int                    linesize;       //!< size of each cache line in bytes, 0 if unknown
            int                    associativity;  //!< associativity, -1 means fully associative
            int                    shared;         //!< Number of logical processors sharing this cache
        };

        /*! \libinternal \brief Information about a single hardware thread in a core
         *
         * The id of the thread typically increases continuously as you walk
         * through sockets and cores in order of their ids. In general, this can
         * be different from the logical processor id provided by the operating
         * system. To achieve better load balancing when using SMT, Linux
         * typically assigns logical processors in a round-robin fashion
         * over all cores.
         */
        struct HWThread
        {
            int                    id;                 //!< Absolute id of this thread in hardware topology
            int                    logicalProcessorId; //!< Id of the operating system logical processor
        };

        /*! \libinternal \brief Information about a single core in a socket */
        struct Core
        {
            int                    id;              //!< Absolute id of this core in hardware topology
            int                    numaNodeId;      //!< id of the numa node of this core
            std::vector<HWThread>  hwThreads;       //!< All the hardware threads in this core
        };

        /*! \libinternal \brief Information about a single socket in the system */
        struct Socket
        {
            int                    id;             //!< Absolute id of this socket in hardware topology
            std::vector<Core>      cores;          //!< All the cores in this socket
        };

        /*! \libinternal \brief Information about each numa node in system */
        struct NumaNode
        {
            int                    id;                 //!< Absolute id of numa node in hardware topology
            std::size_t            memory;             //!< Total detected memory in bytes
            std::vector<int>       logicalProcessorId; //!< Vector of all the logical processors in this node
        };

        /*! \libinternal \brief Information about a single numa node */
        struct Numa
        {
            std::vector<NumaNode>               nodes;              //!< Information about each numa node
            float                               baseLatency;        //!< Scale factor for relative latencies
            std::vector< std::vector<float> >   relativeLatency;    //!< 2D matrix of relative latencies between nodes
            float                               maxRelativeLatency; //!< Largest relative latency
        };

        /*! \libinternal \brief Information about a single PCI device.
         *
         *  \note On many systems the PCI bus is not directly connected to any numa node.
         *        For these systems the numaNodeId will be -1, so you cannot rely on this
         *        number reflecting a specific numa node.
         */
        struct Device
        {
            std::uint16_t          vendorId;    //!< Vendor identification
            std::uint16_t          deviceId;    //!< Vendor-specific device identification
            std::uint16_t          classId;     //!< class (high 8 bits) and subclass (low 8 bits)
            std::uint16_t          domain;      //!< Domain, usually 0 for PCI bus
            std::uint8_t           bus;         //!< Bus number in domain
            std::uint8_t           dev;         //!< Device on bus
            std::uint8_t           func;        //!< Function id for multi-function devices
            int                    numaNodeId;  //!< Numa node, -1 if the bus is not located inside a node
        };

        /*! \libinternal \brief Information about socket, core and hwthread for a logical processor */
        struct LogicalProcessor
        {
            int                    socketRankInMachine; //!< Index of socket in machine
            int                    coreRankInSocket;    //!< Index of core in socket
            int                    hwThreadRankInCore;  //!< Index of hardware thread in core
            int                    numaNodeId;          //!< Index of numa node
        };

        /*! \libinternal \brief Hardware topology information about the entire machine
         *
         * The machine structure is a tree with top-down information about all
         * sockets, cores, and hardware threads in the system. For example, an
         * operating system logical processor index can be found as
         * machine.socket[0].core[1].hwthread[2].logicalProcessorId.
         * In some cases you might need the opposite lookup, i.e. the physical
         * hardware data for a specific logical processor. This is present in the
         * logicalProcessor vector for convenience.
         *
         * \note The logicalProcessor vector will only have non-zero length if the
         *       support level is SupportLevel::Basic or higher. You cannot use the
         *       size of this vector to query the number of logical processors on
         *       lower support levels.
         */
        struct Machine
        {
            Machine();

            int                            logicalProcessorCount; //!< Number of logical processors in system
            std::vector<LogicalProcessor>  logicalProcessors;     //!< Map logical processors to socket/core
            std::vector<Socket>            sockets;               //!< All the sockets in the system
            std::vector<Cache>             caches;                //!< Caches in increasing level order
            Numa                           numa;                  //!< Structure with all numa information
            std::vector<Device>            devices;               //!< Devices on PCI bus
        };

    public:

        /*! \brief Detects the hardware topology. */
        static HardwareTopology detect();

        /*! \brief Creates a topology with given number of logical cores.
         *
         * The support level will be either None or LogicalProcessorCount.
         *
         * Intended for testing of code that uses the hardware topology.
         */
        explicit HardwareTopology(int logicalProcessorCount);

        /*! \brief Check what topology information that is available and valid
         *
         *  The amount of hardware topology information that can be detected depends
         *  on both the hardware and whether GROMACS was linked with the external
         *  hwloc library. You cannot assume that any information is present,
         *  although we can almost always provide the number of logical processors.
         *  On x86 we can usually get basic information about how sockets, cores
         *  and hardware threads are ordered even without hwloc.
         *  With the hwloc library we can usually also get information about cache,
         *  memory and concepts such as core groups and ccNUMA nodes.
         *  Finally, if hwloc was built with support for libpci we can also
         *  detect how the PCI devices are connected.
         */
        SupportLevel
        supportLevel() const { return supportLevel_; }

        /*! \brief Return true if we actually detected hardware.
         *
         *  \return This method will normally return true, when we actually ran
         *          the hardware detection as part of this process to construct
         *          the object. It will be false when the object was constructed
         *          by reading a cached XML file, or possibly generated from
         *          synthetic data.
         */
        bool
        isThisSystem() const { return isThisSystem_; }

        /*! \brief Return the machine topology tree
         *
         *  You can always call this routine, but be aware that some or all contents
         *  will not be valid unless supportLevel() returns a sufficient level.
         *
         *  - With SupportLevel::LogicalProcessorCount, only the field
         *    machine.logicalProcessorCount is valid.
         *  - With SupportLevel::Basic, you can access the vectors of sockets,
         *    cores, and hardware threads, and query what logical processorId
         *    each hardware thread corresponds to.
         *  - SupportLevel::Full adds cache, memory and ccNUMA information.
         *  - SupportLevel::FullWithDevices also adds the PCI express bus.
         *
         *  While data that is not valid has been initialized to special values,
         *  you should not rely on those but query the supportLevel() method before
         *  accessing it.
         */
        const Machine &
        machine() const { return machine_; }

        /*! \brief Returns the number of cores.
         *
         * You can always call this routine, but if sufficient support is not
         * available, it may return the logical processor count or zero instead
         * of the physical core count.
         */
        int numberOfCores() const;

    private:

        HardwareTopology();

        SupportLevel        supportLevel_; //!< Available topology information
        Machine             machine_;      //!< The machine map
        bool                isThisSystem_; //!< Machine map is real (vs. cached/synthetic)
};

}

#endif // GMX_HARDWARE_HARDWARETOPOLOGY_H
