/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief
 * Implements gmx::HardwareTopology.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_hardware
 */

#include "gmxpre.h"

#include "hardwaretopology.h"

#include "config.h"

#include <cstdio>

#include <algorithm>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#if GMX_USE_HWLOC
#    include <hwloc.h>
#endif

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/utility/gmxassert.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h> // sysconf()
#endif
#if GMX_NATIVE_WINDOWS
#    include <windows.h> // GetSystemInfo()
#endif

//! Convenience macro to help us avoid ifdefs each time we use sysconf
#if !defined(_SC_NPROCESSORS_ONLN) && defined(_SC_NPROC_ONLN)
#    define _SC_NPROCESSORS_ONLN _SC_NPROC_ONLN
#endif

namespace gmx
{

namespace
{

/*****************************************************************************
 *                                                                           *
 *   Utility functions for extracting hardware topology from CpuInfo object  *
 *                                                                           *
 *****************************************************************************/

/*! \brief Initialize machine data from basic information in cpuinfo
 *
 *  \param  machine      Machine tree structure where information will be assigned
 *                       if the cpuinfo object contains topology information.
 *  \param  supportLevel If topology information is available in CpuInfo,
 *                       this will be updated to reflect the amount of
 *                       information written to the machine structure.
 */
void parseCpuInfo(HardwareTopology::Machine* machine, HardwareTopology::SupportLevel* supportLevel)
{
    CpuInfo cpuInfo(CpuInfo::detect());

    if (!cpuInfo.logicalProcessors().empty())
    {
        int nSockets   = 0;
        int nCores     = 0;
        int nHwThreads = 0;

        // Copy the logical processor information from cpuinfo
        for (auto& l : cpuInfo.logicalProcessors())
        {
            machine->logicalProcessors.push_back(
                    { l.socketRankInMachine, l.coreRankInSocket, l.hwThreadRankInCore, -1 });
            nSockets   = std::max(nSockets, l.socketRankInMachine);
            nCores     = std::max(nCores, l.coreRankInSocket);
            nHwThreads = std::max(nHwThreads, l.hwThreadRankInCore);
        }

        // Fill info form sockets/cores/hwthreads
        int socketId   = 0;
        int coreId     = 0;
        int hwThreadId = 0;

        machine->sockets.resize(nSockets + 1);
        for (auto& s : machine->sockets)
        {
            s.id = socketId++;
            s.cores.resize(nCores + 1);
            for (auto& c : s.cores)
            {
                c.id         = coreId++;
                c.numaNodeId = -1; // No numa information
                c.hwThreads.resize(nHwThreads + 1);
                for (auto& t : c.hwThreads)
                {
                    t.id                 = hwThreadId++;
                    t.logicalProcessorId = -1; // set as unassigned for now
                }
            }
        }

        // Fill the logical processor id in the right place
        for (std::size_t i = 0; i < machine->logicalProcessors.size(); i++)
        {
            const HardwareTopology::LogicalProcessor& l = machine->logicalProcessors[i];
            machine->sockets[l.socketRankInMachine]
                    .cores[l.coreRankInSocket]
                    .hwThreads[l.hwThreadRankInCore]
                    .logicalProcessorId = static_cast<int>(i);
        }
        machine->logicalProcessorCount = machine->logicalProcessors.size();
        *supportLevel                  = HardwareTopology::SupportLevel::Basic;
    }
    else
    {
        *supportLevel = HardwareTopology::SupportLevel::None;
    }
}

#if GMX_USE_HWLOC

#    if HWLOC_API_VERSION < 0x00010b00
#        define HWLOC_OBJ_PACKAGE HWLOC_OBJ_SOCKET
#        define HWLOC_OBJ_NUMANODE HWLOC_OBJ_NODE
#    endif

// Preprocessor variable for if hwloc api is version 1.x.x or 2.x.x
#    if HWLOC_API_VERSION >= 0x00020000
#        define GMX_HWLOC_API_VERSION_IS_2XX 1
#        if GMX_HWLOC_API_VERSION < 0x00020000
#            error "HWLOC library major version set during configuration is 1, but currently using version 2 headers"
#        endif
#    else
#        define GMX_HWLOC_API_VERSION_IS_2XX 0
#        if GMX_HWLOC_API_VERSION >= 0x00020000
#            error "HWLOC library major version set during configuration is 2, but currently using version 1 headers"
#        endif
#    endif

/*****************************************************************************
 *                                                                           *
 *   Utility functions for extracting hardware topology from hwloc library   *
 *                                                                           *
 *****************************************************************************/

// Compatibility function for accessing hwloc_obj_t object memory with different API versions of hwloc
std::size_t getHwLocObjectMemory(const hwloc_obj* obj)
{
#    if GMX_HWLOC_API_VERSION_IS_2XX
    return obj->total_memory;
#    else
    return obj->memory.total_memory;
#    endif
}

/*! \brief Return vector of all descendants of a given type in hwloc topology
 *
 *  \param topo  hwloc topology handle that has been initialized and loaded
 *  \param obj   Non-null hwloc object.
 *  \param type  hwloc object type to find. The routine will only search
 *               on levels below obj.
 *
 *  \return vector containing all the objects of given type that are
 *          descendants of the provided object. If no objects of this type
 *          were found, the vector will be empty.
 */
std::vector<const hwloc_obj*> getHwLocDescendantsByType(const hwloc_topology*  topo,
                                                        const hwloc_obj*       obj,
                                                        const hwloc_obj_type_t type)
{
    GMX_RELEASE_ASSERT(obj, "NULL hwloc object provided to getHwLocDescendantsByType()");

    std::vector<const hwloc_obj*> v;

    if (obj->type == type)
    {
        v.push_back(obj);
    }
    // Go through children; if this object has no children obj->arity is 0,
    // and we'll return an empty vector.
    hwloc_obj_t tempNode = nullptr;
    while ((tempNode = hwloc_get_next_child(const_cast<hwloc_topology_t>(topo),
                                            const_cast<hwloc_obj_t>(obj), tempNode))
           != nullptr)
    {
        std::vector<const hwloc_obj*> v2 = getHwLocDescendantsByType(topo, tempNode, type);
        v.insert(v.end(), v2.begin(), v2.end());
    }
    return v;
}

/*! \brief Read information about sockets, cores and threads from hwloc topology
 *
 *  \param topo    hwloc topology handle that has been initialized and loaded
 *  \param machine Pointer to the machine structure in the HardwareTopology
 *                 class, where the tree of sockets/cores/threads will be written.
 *
 *  \return If all the data is found
 */
bool parseHwLocSocketsCoresThreads(hwloc_topology_t topo, HardwareTopology::Machine* machine)
{
    const hwloc_obj*              root         = hwloc_get_root_obj(topo);
    std::vector<const hwloc_obj*> hwlocSockets = getHwLocDescendantsByType(topo, root, HWLOC_OBJ_PACKAGE);

    machine->logicalProcessorCount = hwloc_get_nbobjs_by_type(topo, HWLOC_OBJ_PU);
    machine->logicalProcessors.resize(machine->logicalProcessorCount);
    machine->sockets.resize(hwlocSockets.size());

    bool topologyOk = !hwlocSockets.empty(); // Fail if we have no sockets in machine

    for (std::size_t i = 0; i < hwlocSockets.size() && topologyOk; i++)
    {
        // Assign information about this socket
        machine->sockets[i].id = hwlocSockets[i]->logical_index;

        // Get children (cores)
        std::vector<const hwloc_obj*> hwlocCores =
                getHwLocDescendantsByType(topo, hwlocSockets[i], HWLOC_OBJ_CORE);
        machine->sockets[i].cores.resize(hwlocCores.size());

        topologyOk = topologyOk && !hwlocCores.empty(); // Fail if we have no cores in socket

        // Loop over child cores
        for (std::size_t j = 0; j < hwlocCores.size() && topologyOk; j++)
        {
            // Assign information about this core
            machine->sockets[i].cores[j].id         = hwlocCores[j]->logical_index;
            machine->sockets[i].cores[j].numaNodeId = -1;

            // Get children (hwthreads)
            std::vector<const hwloc_obj*> hwlocPUs =
                    getHwLocDescendantsByType(topo, hwlocCores[j], HWLOC_OBJ_PU);
            machine->sockets[i].cores[j].hwThreads.resize(hwlocPUs.size());

            topologyOk = topologyOk && !hwlocPUs.empty(); // Fail if we have no hwthreads in core

            // Loop over child hwthreads
            for (std::size_t k = 0; k < hwlocPUs.size() && topologyOk; k++)
            {
                // Assign information about this hwthread
                std::size_t logicalProcessorId               = hwlocPUs[k]->os_index;
                machine->sockets[i].cores[j].hwThreads[k].id = hwlocPUs[k]->logical_index;
                machine->sockets[i].cores[j].hwThreads[k].logicalProcessorId = logicalProcessorId;

                if (logicalProcessorId < machine->logicalProcessors.size())
                {
                    // Cross-assign data for this hwthread to the logicalprocess vector
                    machine->logicalProcessors[logicalProcessorId].socketRankInMachine =
                            static_cast<int>(i);
                    machine->logicalProcessors[logicalProcessorId].coreRankInSocket =
                            static_cast<int>(j);
                    machine->logicalProcessors[logicalProcessorId].hwThreadRankInCore =
                            static_cast<int>(k);
                    machine->logicalProcessors[logicalProcessorId].numaNodeId = -1;
                }
                else
                {
                    topologyOk = false;
                }
            }
        }
    }

    if (!topologyOk)
    {
        machine->logicalProcessors.clear();
        machine->sockets.clear();
    }
    return topologyOk;
}

/*! \brief Read cache information from hwloc topology
 *
 *  \param topo    hwloc topology handle that has been initialized and loaded
 *  \param machine Pointer to the machine structure in the HardwareTopology
 *                 class, where cache data will be filled.
 *
 *  \return If any cache data is found
 */
bool parseHwLocCache(hwloc_topology_t topo, HardwareTopology::Machine* machine)
{
    // Parse caches up to L5
    for (int cachelevel : { 1, 2, 3, 4, 5 })
    {
        int depth = hwloc_get_cache_type_depth(topo, cachelevel, HWLOC_OBJ_CACHE_DATA);

        if (depth >= 0)
        {
            hwloc_obj_t cache = hwloc_get_next_obj_by_depth(topo, depth, nullptr);
            if (cache != nullptr)
            {
                std::vector<const hwloc_obj*> hwThreads =
                        getHwLocDescendantsByType(topo, cache, HWLOC_OBJ_PU);

                machine->caches.push_back({ static_cast<int>(cache->attr->cache.depth),
                                            static_cast<std::size_t>(cache->attr->cache.size),
                                            static_cast<int>(cache->attr->cache.linesize),
                                            static_cast<int>(cache->attr->cache.associativity),
                                            std::max<int>(hwThreads.size(), 1) });
            }
        }
    }
    return !machine->caches.empty();
}


/*! \brief Read numa information from hwloc topology
 *
 *  \param topo    hwloc topology handle that has been initialized and loaded
 *  \param machine Pointer to the machine structure in the HardwareTopology
 *                 class, where numa information will be filled.
 *
 *  Hwloc should virtually always be able to detect numa information, but if
 *  there is only a single numa node in the system it is not reported at all.
 *  In this case we create a single numa node covering all cores.
 *
 *  This function uses the basic socket/core/thread information detected by
 *  parseHwLocSocketsCoresThreads(), which means that routine must have
 *  completed successfully before calling this one. If this is not the case,
 *  you will get an error return code.
 *
 *  \return If the data found makes sense (either in the numa node or the
 *          entire machine)
 */
bool parseHwLocNuma(hwloc_topology_t topo, HardwareTopology::Machine* machine)
{
    const hwloc_obj*              root = hwloc_get_root_obj(topo);
    std::vector<const hwloc_obj*> hwlocNumaNodes =
            getHwLocDescendantsByType(topo, root, HWLOC_OBJ_NUMANODE);
    bool topologyOk = true;

    if (!hwlocNumaNodes.empty())
    {
        machine->numa.nodes.resize(hwlocNumaNodes.size());

        for (std::size_t i = 0; i < hwlocNumaNodes.size(); i++)
        {
            machine->numa.nodes[i].id     = hwlocNumaNodes[i]->logical_index;
            machine->numa.nodes[i].memory = getHwLocObjectMemory(hwlocNumaNodes[i]);

            machine->numa.nodes[i].logicalProcessorId.clear();

            // Get list of PUs in this numa node. Get from numa node if v1.x.x, get from numa node's parent if 2.x.x
#    if GMX_HWLOC_API_VERSION_IS_2XX
            std::vector<const hwloc_obj*> hwlocPUs =
                    getHwLocDescendantsByType(topo, hwlocNumaNodes[i]->parent, HWLOC_OBJ_PU);
#    else
            std::vector<const hwloc_obj*> hwlocPUs =
                    getHwLocDescendantsByType(topo, hwlocNumaNodes[i], HWLOC_OBJ_PU);
#    endif
            for (auto& p : hwlocPUs)
            {
                machine->numa.nodes[i].logicalProcessorId.push_back(p->os_index);

                GMX_RELEASE_ASSERT(p->os_index < machine->logicalProcessors.size(),
                                   "OS index of PU in hwloc larger than processor count");

                machine->logicalProcessors[p->os_index].numaNodeId = static_cast<int>(i);
                std::size_t s = machine->logicalProcessors[p->os_index].socketRankInMachine;
                std::size_t c = machine->logicalProcessors[p->os_index].coreRankInSocket;

                GMX_RELEASE_ASSERT(s < machine->sockets.size(),
                                   "Socket index in logicalProcessors larger than socket count");
                GMX_RELEASE_ASSERT(c < machine->sockets[s].cores.size(),
                                   "Core index in logicalProcessors larger than core count");
                // Set numaNodeId in core too
                machine->sockets[s].cores[c].numaNodeId = i;
            }
        }
        // Getting the distance matrix
#    if GMX_HWLOC_API_VERSION_IS_2XX
        // with hwloc api v. 2.x.x, distances are no longer directly accessible. Need to retrieve and release hwloc_distances_s object
        // In addition, there can now be multiple types of distances, ie latency, bandwidth. We look only for latency, but have to check
        // if multiple distance matrices are returned.

        // If only 1 numa node exists, the v2.x.x hwloc api won't have a distances matrix, set manually
        if (hwlocNumaNodes.size() == 1)
        {
            machine->numa.relativeLatency = { { 1.0 } };
        }
        else
        {
            hwloc_distances_s* dist;
            // Set the number of distance matrices to return (1 in our case, but hwloc 2.x.x allows
            // for multiple distances types and therefore multiple distance matrices)
            unsigned nr = 1;
            hwloc_distances_get(topo, &nr, &dist, HWLOC_DISTANCES_KIND_MEANS_LATENCY, 0);
            // If no distances were found, nr will be 0, otherwise distances will be populated with
            // 1 hwloc_distances_s object
            if (nr > 0 && dist->nbobjs == hwlocNumaNodes.size())
            {

                machine->numa.relativeLatency.resize(dist->nbobjs);
                for (std::size_t i = 0; i < dist->nbobjs; i++)
                {
                    machine->numa.relativeLatency[i].resize(dist->nbobjs);
                    for (std::size_t j = 0; j < dist->nbobjs; j++)
                    {
                        machine->numa.relativeLatency[i][j] = dist->values[i * dist->nbobjs + j];
                    }
                }
            }
            else
            {
                topologyOk = false;
            }
            hwloc_distances_release(topo, dist);
        }

        // hwloc-2.x provides latencies as integers, but to make things more similar to the case of
        // a single numa node as well as hwloc-1.x, we rescale to relative floating-point values and
        // also set the largest relative latency value.

        // find smallest value in matrix
        float minLatency = std::numeric_limits<float>::max(); // large number
        float maxLatency = std::numeric_limits<float>::min(); // 0.0
        for (const auto& v : machine->numa.relativeLatency)
        {
            auto result = std::minmax_element(v.begin(), v.end());
            minLatency  = std::min(minLatency, *result.first);
            maxLatency  = std::max(maxLatency, *result.second);
        }

        // assign stuff
        for (auto& v : machine->numa.relativeLatency)
        {
            std::transform(v.begin(), v.end(), v.begin(),
                           std::bind(std::multiplies<float>(), std::placeholders::_1, 1.0 / minLatency));
        }
        machine->numa.baseLatency = 1.0; // latencies still do not have any units in hwloc-2.x
        machine->numa.maxRelativeLatency = maxLatency / minLatency;

#    else  // GMX_HWLOC_API_VERSION_IS_2XX == false, hwloc api is 1.x.x
        int                             depth = hwloc_get_type_depth(topo, HWLOC_OBJ_NUMANODE);
        const struct hwloc_distances_s* dist = hwloc_get_whole_distance_matrix_by_depth(topo, depth);
        if (dist != nullptr && dist->nbobjs == hwlocNumaNodes.size())
        {
            machine->numa.baseLatency        = dist->latency_base;
            machine->numa.maxRelativeLatency = dist->latency_max;
            machine->numa.relativeLatency.resize(dist->nbobjs);
            for (std::size_t i = 0; i < dist->nbobjs; i++)
            {
                machine->numa.relativeLatency[i].resize(dist->nbobjs);
                for (std::size_t j = 0; j < dist->nbobjs; j++)
                {
                    machine->numa.relativeLatency[i][j] = dist->latency[i * dist->nbobjs + j];
                }
            }
        }
        else
        {
            topologyOk = false;
        }
#    endif // end GMX_HWLOC_API_VERSION_IS_2XX == false
    }
    else
    // Deals with the case of no numa nodes found.
#    if GMX_HWLOC_API_VERSION_IS_2XX
    // If the hwloc version is 2.x.x, and there is no numa node, something went wrong
    {
        topologyOk = false;
    }
#    else
    {
        // No numa nodes found. Use the entire machine as a numa node.
        // Note that this should only be the case with hwloc api v 1.x.x,
        // a numa node is assigned to the machine by default in v 2.x.x
        const hwloc_obj* const hwlocMachine = hwloc_get_next_obj_by_type(topo, HWLOC_OBJ_MACHINE, nullptr);

        if (hwlocMachine != nullptr)
        {
            machine->numa.nodes.resize(1);
            machine->numa.nodes[0].id        = 0;
            machine->numa.nodes[0].memory    = hwlocMachine->memory.total_memory;
            machine->numa.baseLatency        = 10;
            machine->numa.maxRelativeLatency = 1;
            machine->numa.relativeLatency    = { { 1.0 } };

            for (int i = 0; i < machine->logicalProcessorCount; i++)
            {
                machine->numa.nodes[0].logicalProcessorId.push_back(i);
            }
            for (auto& l : machine->logicalProcessors)
            {
                l.numaNodeId = 0;
            }
            for (auto& s : machine->sockets)
            {
                for (auto& c : s.cores)
                {
                    c.numaNodeId = 0;
                }
            }
        }
        else
        {
            topologyOk = false;
        }
    }
#    endif // end if not GMX_HWLOC_API_VERSION_IS_2XX
    if (!topologyOk)
    {
        machine->numa.nodes.clear();
    }
    return topologyOk;
}

/*! \brief Read PCI device information from hwloc topology
 *
 *  \param topo    hwloc topology handle that has been initialized and loaded
 *  \param machine Pointer to the machine structure in the HardwareTopology
 *                 class, where PCI device information will be filled.
 * *
 *  \return If any devices were found
 */
bool parseHwLocDevices(hwloc_topology_t topo, HardwareTopology::Machine* machine)
{
    const hwloc_obj*              root = hwloc_get_root_obj(topo);
    std::vector<const hwloc_obj*> pcidevs = getHwLocDescendantsByType(topo, root, HWLOC_OBJ_PCI_DEVICE);

    for (auto& p : pcidevs)
    {
#    if GMX_HWLOC_API_VERSION_IS_2XX
        const hwloc_obj* ancestor = nullptr;
        // Numa nodes not directly part of tree. Walk up the tree until we find an ancestor with a numa node
        hwloc_obj_t parent = p->parent;
        while (parent && !parent->memory_arity)
        {
            parent = parent->parent;
        }
        if (parent)
        {
            ancestor = parent->memory_first_child;
        }
#    else  // GMX_HWLOC_API_VERSION_IS_2XX = false, api v 1.x.x
        // numa nodes are normal part of tree, can use hwloc ancestor function
        const hwloc_obj* const ancestor =
                hwloc_get_ancestor_obj_by_type(topo, HWLOC_OBJ_NUMANODE, const_cast<hwloc_obj_t>(p));
#    endif // end if GMX_HWLOC_API_VERSION_IS_2XX
        int numaId;
        if (ancestor != nullptr)
        {
            numaId = ancestor->logical_index;
        }
        else
        {
            // If we only have a single numa node we belong to it, otherwise set it to -1 (unknown)
            numaId = (machine->numa.nodes.size() == 1) ? 0 : -1;
        }

        GMX_RELEASE_ASSERT(p->attr, "Attributes should not be NULL for hwloc PCI object");

        machine->devices.push_back({ p->attr->pcidev.vendor_id, p->attr->pcidev.device_id,
                                     p->attr->pcidev.class_id, p->attr->pcidev.domain,
                                     p->attr->pcidev.bus, p->attr->pcidev.dev, p->attr->pcidev.func,
                                     numaId });
    }
    return !pcidevs.empty();
}

void parseHwLoc(HardwareTopology::Machine* machine, HardwareTopology::SupportLevel* supportLevel, bool* isThisSystem)
{
    hwloc_topology_t topo;

    // Initialize a hwloc object, set flags to request IO device information too,
    // try to load the topology, and get the root object. If either step fails,
    // return that we do not have any support at all from hwloc.
    if (hwloc_topology_init(&topo) != 0)
    {
        hwloc_topology_destroy(topo);
        return; // SupportLevel::None.
    }

    // Flags to look for io devices
#    if GMX_HWLOC_API_VERSION_IS_2XX
    GMX_RELEASE_ASSERT(
            (hwloc_get_api_version() >= 0x20000),
            "Mismatch between hwloc headers and library, using v2 headers with v1 library");
    hwloc_topology_set_io_types_filter(topo, HWLOC_TYPE_FILTER_KEEP_IMPORTANT);
#    else
    GMX_RELEASE_ASSERT(
            (hwloc_get_api_version() < 0x20000),
            "Mismatch between hwloc headers and library, using v1 headers with v2 library");
    hwloc_topology_set_flags(topo, HWLOC_TOPOLOGY_FLAG_IO_DEVICES);
#    endif

    if (hwloc_topology_load(topo) != 0 || hwloc_get_root_obj(topo) == nullptr)
    {
        hwloc_topology_destroy(topo);
        return; // SupportLevel::None.
    }

    // If we get here, we can get a valid root object for the topology
    *isThisSystem = hwloc_topology_is_thissystem(topo) != 0;

    // Parse basic information about sockets, cores, and hardware threads
    if (parseHwLocSocketsCoresThreads(topo, machine))
    {
        *supportLevel = HardwareTopology::SupportLevel::Basic;
    }
    else
    {
        hwloc_topology_destroy(topo);
        return; // SupportLevel::None.
    }

    // Get information about cache and numa nodes
    if (parseHwLocCache(topo, machine) && parseHwLocNuma(topo, machine))
    {
        *supportLevel = HardwareTopology::SupportLevel::Full;
    }
    else
    {
        hwloc_topology_destroy(topo);
        return; // SupportLevel::Basic.
    }

    // PCI devices
    if (parseHwLocDevices(topo, machine))
    {
        *supportLevel = HardwareTopology::SupportLevel::FullWithDevices;
    }

    hwloc_topology_destroy(topo);
    // SupportLevel::Full or SupportLevel::FullWithDevices.
}

#endif

/*! \brief Try to detect the number of logical processors.
 *
 *  \return The number of hardware processing units, or 0 if it fails.
 */
int detectLogicalProcessorCount()
{
    int count = 0;

    {
#if GMX_NATIVE_WINDOWS
        // Windows
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        count = sysinfo.dwNumberOfProcessors;
#elif defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_ONLN)
        // We are probably on Unix. Check if we have the argument to use before executing any calls
        count = sysconf(_SC_NPROCESSORS_ONLN);
#else
        count = 0; // Neither windows nor Unix.
#endif
    }

    return count;
}

} // namespace

// static
HardwareTopology HardwareTopology::detect()
{
    HardwareTopology result;

#if GMX_USE_HWLOC
    parseHwLoc(&result.machine_, &result.supportLevel_, &result.isThisSystem_);
#endif

    // If something went wrong in hwloc (or if it was not present) we might
    // have more information in cpuInfo
    if (result.supportLevel_ < SupportLevel::Basic)
    {
        // There might be topology information in cpuInfo
        parseCpuInfo(&result.machine_, &result.supportLevel_);
    }
    // If we did not manage to get anything from either hwloc or cpuInfo, find the cpu count at least
    if (result.supportLevel_ == SupportLevel::None)
    {
        // No topology information; try to detect the number of logical processors at least
        result.machine_.logicalProcessorCount = detectLogicalProcessorCount();
        if (result.machine_.logicalProcessorCount > 0)
        {
            result.supportLevel_ = SupportLevel::LogicalProcessorCount;
        }
    }
    return result;
}

HardwareTopology::Machine::Machine()
{
    logicalProcessorCount   = 0;
    numa.baseLatency        = 0.0;
    numa.maxRelativeLatency = 0.0;
}


HardwareTopology::HardwareTopology() :
    supportLevel_(SupportLevel::None),
    machine_(),
    isThisSystem_(true)
{
}

HardwareTopology::HardwareTopology(int logicalProcessorCount) :
    supportLevel_(SupportLevel::None),
    machine_(),
    isThisSystem_(true)
{
    if (logicalProcessorCount > 0)
    {
        machine_.logicalProcessorCount = logicalProcessorCount;
        supportLevel_                  = SupportLevel::LogicalProcessorCount;
    }
}

int HardwareTopology::numberOfCores() const
{
    if (supportLevel() >= SupportLevel::Basic)
    {
        // We assume all sockets have the same number of cores as socket 0.
        // Since topology information is present, we can assume there is at least one socket.
        return machine().sockets.size() * machine().sockets[0].cores.size();
    }
    else if (supportLevel() >= SupportLevel::LogicalProcessorCount)
    {
        return machine().logicalProcessorCount;
    }
    else
    {
        return 0;
    }
}

} // namespace gmx
