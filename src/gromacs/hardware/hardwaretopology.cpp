/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#include <cstdlib>

#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#if GMX_USE_HWLOC
#    include <hwloc.h>
#endif

#include <sys/types.h>

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/utility/gmxassert.h"

#ifdef HAVE_SCHED_H
#    include <sched.h>
#endif
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

/*! \brief Utlility function to renumber and translate low-level APIC info to topology
 *
 * \param logicalProcessors Logical processor information according to the CpuInfo
 *                          structure. Note that the indices refer
 *                          to low-level hardware, so they will
 *                          be renumbered to only consider
 *                          logical cores/processing units we use.
 * \param machine           Hardware topology machine structure where result is written.
 */
void translateCpuInfoLogicalProcessorsToMachine(const std::vector<CpuInfo::LogicalProcessor>& logicalProcessors,
                                                HardwareTopology::Machine* machine)
{
    // We will keep and report the os-provided indices for packages, since we need to be
    // able to identify e.g. that the GPU is connected to package 0 while we run on package 1,
    // even if we don't have any cores on package 0. Still, to avoid having empty entries in our
    // vector of packages, we need map to renumber from the seen packages to the logical indices.
    std::unordered_map<int, int> renumPkg;

    // When it comes to cores, they must be renumbered, including the index we list for them,
    // since at least Intel/AMD sometimes use non-consecutive core enumeration, and we want our
    // own core ids to be sequential even if we only run on a subset of them. This is a vector over
    // packages where each index is a map from hardware core id to our logical core id in that package.
    std::vector<std::unordered_map<int, int>> renumCorePerPkg;

    // For the processing units we will be A-OK without renumbering, since we don't use the actual
    // indices they have (we create new ones), and each PU will only appear once.
    for (const auto& p : logicalProcessors)
    {
        int oldPkg = p.packageIdInMachine;
        if (renumPkg.find(oldPkg) == renumPkg.end())
        {
            renumPkg[oldPkg] = renumPkg.size();
        }
        int newPkg = renumPkg[oldPkg];

        // Extend translation map if we have not seen this package before
        if (renumCorePerPkg.size() <= static_cast<std::size_t>(newPkg))
        {
            renumCorePerPkg.resize(newPkg + 1);
        }
        // Add new package in tree structure if we have not seen it before
        if (machine->packages.size() <= static_cast<std::size_t>(newPkg))
        {
            machine->packages.resize(newPkg + 1);
            machine->packages[newPkg].id = oldPkg;
        }

        // Check if this core has been seen before, add to translation map if not
        if (renumCorePerPkg[newPkg].find(p.coreIdInPackage) == renumCorePerPkg[newPkg].end())
        {
            renumCorePerPkg[newPkg][p.coreIdInPackage] = renumCorePerPkg[newPkg].size();
        }
        // Get translated core id
        int newCore = renumCorePerPkg[newPkg][p.coreIdInPackage];

        // Add new core in package (for the translated core id) if we have not already seen it
        if (machine->packages[newPkg].cores.size() <= static_cast<std::size_t>(newCore))
        {
            machine->packages[newPkg].cores.resize(newCore + 1);
            machine->packages[newPkg].cores[newCore].numaNodeId = -1;
            // We cannot assign global core id yet
        }
        // Note: We cannot assign global PU id yet, so use -1 for now
        machine->packages[newPkg].cores[newCore].processingUnits.push_back({ -1, p.osId });
    }

    // Fill linear structure of logical processors and assign global core/PU id in tree
    for (int pkg = 0, coreId = 0, puId = 0; static_cast<std::size_t>(pkg) < machine->packages.size(); pkg++)
    {
        for (int core = 0; static_cast<std::size_t>(core) < machine->packages[pkg].cores.size(); core++)
        {
            machine->packages[pkg].cores[core].id = coreId++; // Global core id
            for (int pu = 0; static_cast<std::size_t>(pu)
                             < machine->packages[pkg].cores[core].processingUnits.size();
                 pu++)
            {
                int osId = machine->packages[pkg].cores[core].processingUnits[pu].osId;
                // No numa info, set it to -1.
                machine->logicalProcessors.push_back({ puId, osId, pkg, core, pu, -1 });
                machine->osIdToPuId.insert({ osId, puId });
                machine->packages[pkg].cores[core].processingUnits[pu].id = puId++; // global PU id
            }
        }
    }
}

/*! \brief Initialize machine data from basic information in cpuinfo
 *
 *  \param  machine      Machine tree structure where information will be assigned
 *                       if the cpuinfo object contains topology information.
 *  \return   SupportLevel::Basic if topology information was found.
 */
HardwareTopology::SupportLevel parseCpuInfo(HardwareTopology::Machine* machine)
{
    CpuInfo cpuInfo(CpuInfo::detect());

    if (cpuInfo.supportLevel() >= CpuInfo::SupportLevel::LogicalProcessorInfo)
    {
        translateCpuInfoLogicalProcessorsToMachine(cpuInfo.logicalProcessors(), machine);
        return HardwareTopology::SupportLevel::Basic;
    }
    else
    {
        return HardwareTopology::SupportLevel::None;
    }
}

/*****************************************************************************
 *                                                                           *
 *   Utility functions for extracting hardware topology from hwloc library   *
 *                                                                           *
 *****************************************************************************/

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
 *          were found, or the input was invalid, the vector will be empty.
 */
std::vector<const hwloc_obj*> getHwLocDescendantsByType(const hwloc_topology*  topo,
                                                        const hwloc_obj*       obj,
                                                        const hwloc_obj_type_t type)
{
    std::vector<const hwloc_obj*> v;

    if (topo == nullptr || obj == nullptr)
    {
        return v;
    }

    if (obj->type == type)
    {
        v.push_back(obj);
    }
    // Go through children; if this object has no children obj->arity is 0,
    // and we'll return an empty vector.
    hwloc_obj_t tempNode = nullptr;
    while ((tempNode = hwloc_get_next_child(
                    const_cast<hwloc_topology_t>(topo), const_cast<hwloc_obj_t>(obj), tempNode))
           != nullptr)
    {
        std::vector<const hwloc_obj*> v2 = getHwLocDescendantsByType(topo, tempNode, type);
        v.insert(v.end(), v2.begin(), v2.end());
    }
    return v;
}

/*! \brief Read information about packages, cores and processing units from hwloc topology
 *
 *  \param topo    hwloc topology handle that has been initialized and loaded
 *  \param machine Pointer to the machine structure in the HardwareTopology
 *                 class, where the tree of packages/cores/PUs will be written.
 *
 *  \return true If all the data is found
 */
bool parseHwLocPackagesCoresProcessingUnits(hwloc_topology_t topo, HardwareTopology::Machine* machine)
{
    const hwloc_obj*              root          = hwloc_get_root_obj(topo);
    std::vector<const hwloc_obj*> hwlocPackages = getHwLocDescendantsByType(topo, root, HWLOC_OBJ_PACKAGE);

    std::size_t puCount = hwloc_get_nbobjs_by_type(topo, HWLOC_OBJ_PU);
    machine->logicalProcessors.resize(puCount);
    machine->packages.resize(hwlocPackages.size());

    for (std::size_t i = 0; i < hwlocPackages.size(); i++)
    {
        // Assign information about this package
        machine->packages[i].id = hwlocPackages[i]->logical_index;

        // Get children (cores)
        std::vector<const hwloc_obj*> hwlocCores =
                getHwLocDescendantsByType(topo, hwlocPackages[i], HWLOC_OBJ_CORE);
        machine->packages[i].cores.resize(hwlocCores.size());

        // Loop over child cores
        for (std::size_t j = 0; j < hwlocCores.size(); j++)
        {
            // Assign information about this core
            machine->packages[i].cores[j].id         = hwlocCores[j]->logical_index;
            machine->packages[i].cores[j].numaNodeId = -1;

            // Get children (PUs)
            std::vector<const hwloc_obj*> hwlocPUs =
                    getHwLocDescendantsByType(topo, hwlocCores[j], HWLOC_OBJ_PU);
            machine->packages[i].cores[j].processingUnits.resize(hwlocPUs.size());

            // Loop over child PUs
            for (std::size_t k = 0; k < hwlocPUs.size(); k++)
            {
                // Assign information about this processing unit
                const std::size_t puId                                = hwlocPUs[k]->logical_index;
                const std::size_t osId                                = hwlocPUs[k]->os_index;
                machine->packages[i].cores[j].processingUnits[k].id   = static_cast<int>(puId);
                machine->packages[i].cores[j].processingUnits[k].osId = static_cast<int>(osId);

                // Cross-assign data for this processing unit to the logicalProcessors vector
                machine->logicalProcessors[puId].puId                     = static_cast<int>(puId);
                machine->logicalProcessors[puId].osId                     = static_cast<int>(osId);
                machine->logicalProcessors[puId].packageRankInTopology    = static_cast<int>(i);
                machine->logicalProcessors[puId].coreRankInPackage        = static_cast<int>(j);
                machine->logicalProcessors[puId].processingUnitRankInCore = static_cast<int>(k);
                machine->logicalProcessors[puId].numaNodeId               = -1;

                // Use a convenience map so we can easily find the internal logical
                // processing unit id based on the OS-provided one if necessary.
                machine->osIdToPuId.insert({ osId, puId });
            }
        }
    }
    return true; // for now we can't really fail cleanly, but keep the option to signal failed detection
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
        const int depth = hwloc_get_cache_type_depth(topo, cachelevel, HWLOC_OBJ_CACHE_DATA);

        if (depth >= 0)
        {
            hwloc_obj_t cache = hwloc_get_next_obj_by_depth(topo, depth, nullptr);
            if (cache != nullptr)
            {
                std::vector<const hwloc_obj*> hwlocPUs =
                        getHwLocDescendantsByType(topo, cache, HWLOC_OBJ_PU);

                machine->caches.push_back({ static_cast<int>(cache->attr->cache.depth),
                                            static_cast<std::size_t>(cache->attr->cache.size),
                                            static_cast<int>(cache->attr->cache.linesize),
                                            static_cast<int>(cache->attr->cache.associativity),
                                            std::max<int>(hwlocPUs.size(), 1) });
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
 *  This function uses the basic package/core/processing unit information detected by
 *  parseHwLocPackagesCoresProcessingUnits(), which means that routine must have
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

            machine->numa.nodes[i].processingUnits.clear();

            // Get list of PUs in this numa node. Get from numa node if v1.x.x, get from numa node's parent if 2.x.x
#    if GMX_HWLOC_API_VERSION_IS_2XX
            std::vector<const hwloc_obj*> hwlocPUs =
                    getHwLocDescendantsByType(topo, hwlocNumaNodes[i]->parent, HWLOC_OBJ_PU);
#    else
            std::vector<const hwloc_obj*> hwlocPUs =
                    getHwLocDescendantsByType(topo, hwlocNumaNodes[i], HWLOC_OBJ_PU);
#    endif
            for (const auto& pu : hwlocPUs)
            {
                machine->numa.nodes[i].processingUnits.push_back(pu->logical_index);

                machine->logicalProcessors[pu->logical_index].numaNodeId = static_cast<int>(i);
                std::size_t pkgRank = machine->logicalProcessors[pu->logical_index].packageRankInTopology;
                std::size_t coreRank = machine->logicalProcessors[pu->logical_index].coreRankInPackage;

                if (pkgRank >= machine->packages.size()
                    || coreRank >= machine->packages[pkgRank].cores.size())
                {
                    // Inconsistent indices. Reset the entire numa structure
                    // and return false to indiciate that data was not valid.
                    machine->numa.nodes.clear();
                    return false;
                }

                // Set numaNodeId in core too
                machine->packages[pkgRank].cores[coreRank].numaNodeId = i;
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
            // Rescale the latencies to a relative float-point value
            for (float& value : v)
            {
                value /= minLatency;
            }
        }
        machine->numa.baseLatency = 1.0; // latencies still do not have any units in hwloc-2.x
        machine->numa.maxRelativeLatency = maxLatency / minLatency;

#    else  // GMX_HWLOC_API_VERSION_IS_2XX == false, hwloc api is 1.x.x
        const int                       depth = hwloc_get_type_depth(topo, HWLOC_OBJ_NUMANODE);
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

            for (auto& l : machine->logicalProcessors)
            {
                machine->numa.nodes[0].processingUnits.push_back(l.puId);
                l.numaNodeId = 0;
            }
            for (auto& s : machine->packages)
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

    for (const auto& p : pcidevs)
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

        if (p->attr == nullptr)
        {
            // Attributes should not be NULL for hwloc PCI object. Skip this device.
            continue;
        }

        machine->devices.push_back({ p->attr->pcidev.vendor_id,
                                     p->attr->pcidev.device_id,
                                     p->attr->pcidev.class_id,
                                     p->attr->pcidev.domain,
                                     p->attr->pcidev.bus,
                                     p->attr->pcidev.dev,
                                     p->attr->pcidev.func,
                                     numaId });
    }
    return !pcidevs.empty();
}

/*! \brief Populate hardware topology from hwloc information
 *
 * \param machine Pointer to the machine structure in the  topology  to be populated.
 * \param isThisSystem Pointer that will be set to true if detection is for this system.
 *
 * \return A supportLevel we managed to get from hwloc if it worked.
 */
HardwareTopology::SupportLevel parseHwLoc(HardwareTopology::Machine* machine, bool* isThisSystem)
{
    hwloc_topology_t topo;

    // Initialize a hwloc object, set flags to request IO device information too,
    // try to load the topology, and get the root object. If either step fails,
    // return that we do not have any support at all from hwloc.
    if (hwloc_topology_init(&topo) != 0)
    {
        hwloc_topology_destroy(topo);
        return HardwareTopology::SupportLevel::None;
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
        return HardwareTopology::SupportLevel::None;
    }

    // If we get here, we can get a valid root object for the topology
    *isThisSystem = hwloc_topology_is_thissystem(topo) != 0;

    HardwareTopology::SupportLevel supportLevel;
    // Parse basic information about packages, cores, and processing units
    if (parseHwLocPackagesCoresProcessingUnits(topo, machine))
    {
        supportLevel = HardwareTopology::SupportLevel::Basic;
    }
    else
    {
        hwloc_topology_destroy(topo);
        return HardwareTopology::SupportLevel::None;
    }

    // Get information about cache and numa nodes
    if (parseHwLocCache(topo, machine) && parseHwLocNuma(topo, machine))
    {
        supportLevel = HardwareTopology::SupportLevel::Full;
    }
    else
    {
        hwloc_topology_destroy(topo);
        // keep previously set SupportLevel::Basic.
        return supportLevel;
    }

    // PCI devices
    if (parseHwLocDevices(topo, machine))
    {
        supportLevel = HardwareTopology::SupportLevel::FullWithDevices;
    }

    hwloc_topology_destroy(topo);
    // SupportLevel::Full or SupportLevel::FullWithDevices.
    return supportLevel;
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

/*! \brief Parse cpu list (integers) from string separated by commas and ranges
 *
 * \param cpuString String with integers in the standard linux cpulist format, i.e.
 *                  integers and integer ranges (specified with dashes) separated
 *                  by commas like "0,3-5,7,8,10-15". If the string is ill-formed,
 *                  or If any integer is negative or the indices are
 *                  not strictly increasing, we return an empty vector.
 *
 * \return vector containing the integers specified, with ranges expanded
 *         to all entries.
 */
std::vector<int> parseCpuString(const std::string& cpuString)
{
    std::vector<int> cpus;

    // We first split the string based on commas into tokens
    // that are either single indices ('4') or intervals ('4-6')
    std::istringstream       cpuStringStream(cpuString);
    std::vector<std::string> cpuIntervals;
    std::string              s1;
    while (std::getline(cpuStringStream, s1, ','))
    {
        cpuIntervals.emplace_back(s1);
    }

    for (const std::string& singleInterval : cpuIntervals)
    {
        // Split each intervals into limits based on '-'.
        std::istringstream       singleIntervalStream(singleInterval);
        std::vector<std::string> singleIntervalParts;
        std::string              s2;
        while (std::getline(singleIntervalStream, s2, '-'))
        {
            singleIntervalParts.emplace_back(s2);
        }
        if (singleIntervalParts.size() == 1)
        {
            cpus.emplace_back(std::stoi(singleIntervalParts[0]));
        }
        else if (singleIntervalParts.size() == 2)
        {
            const int low  = std::stoi(singleIntervalParts[0]);
            const int high = std::stoi(singleIntervalParts[1]);
            for (int i = low; i <= high; i++)
            {
                if (i < 0 || (!cpus.empty() && cpus[cpus.size() - 1] >= i))
                {
                    return {}; // data in string is bad
                }
                cpus.push_back(i);
            }
        }
        else
        {
            return {}; // string not well-formatted; interval has 0 or >2 parts
        }
    }
    return cpus;
}

/*! \brief Attempt to read basic topology from Linux sysfs interface
 *
 * \param machine Pointer to machine structure in topology to be populated
 * \param root    Optional path to (mock) file system root. If this is
 *                provided, all file system access will be relative to this
 *                path instead. This allows us to create mock topologies
 *                based on saved data.
 * \param allowedCpus When creating mock topologies by setting
 *                    the root parameter, we will not attempt to
 *                    check if we can run on processors based on
 *                    cpuset affinity masks, but you should use this
 *                    parameter to provide a vector of logical (OS)
 *                    processor indices on which we are allowed to run.
 *
 * This is a poor man's version of the much fancier topology detection
 * available from hwloc, but given the compilcations of modern hardware
 * we are critically dependent on being able to understand at least
 * how many sockets/cores/threads we have, even when Gromacs is compiled
 * without hwloc support.
 *
 *  \return SupportLevel::Basic if topology information was found.
 */
HardwareTopology::SupportLevel parseSysFsCpuTopology(HardwareTopology::Machine* machine,
                                                     const std::string&         root        = "",
                                                     const std::vector<int>&    allowedCpus = {})
{
    std::string possibleCpuString;
    std::getline(std::ifstream(root + "/sys/devices/system/cpu/possible"), possibleCpuString);
    std::vector<int> rawCpuList =
            parseCpuString(possibleCpuString); // Vector will be empty if file did not exist
    std::vector<int> okCpuList;

    if (root.empty())
    {
#if HAVE_SCHED_AFFINITY
        // Presently, Linux lists all CPUs as possible, rather than the cpu set for the process.
        // Instead of going through the difficulties of trying to parse cgroups1/2, we just
        // check what CPUs we are allowed to run on in practice.
        cpu_set_t saveCpuSet;
        sched_getaffinity(0, sizeof(cpu_set_t), &saveCpuSet);
        // If we did not get any contents of the file above, we won't loop over anything here
        for (int cpu : rawCpuList)
        {
            cpu_set_t cpuSet;
            CPU_ZERO(&cpuSet);
            CPU_SET(cpu, &cpuSet);
            if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuSet) == 0)
            {
                okCpuList.push_back(cpu);
            }
        }
        sched_setaffinity(0, sizeof(cpu_set_t), &saveCpuSet);
#else
        // Bad luck - we couldn't test if we are actually allowed to run on the CPUs, so include all of them.
        okCpuList = rawCpuList;
#endif
    }
    else
    {
        // We are using a mock root path, so filter Cpus based on allowedCpus vector instead
        for (int cpu : rawCpuList)
        {
            if (std::find(allowedCpus.begin(), allowedCpus.end(), cpu) != allowedCpus.end())
            {
                okCpuList.push_back(cpu);
            }
        }
    }

    std::vector<CpuInfo::LogicalProcessor> logicalProcessors;

    for (int c : okCpuList)
    {
        const std::string logicalCpuTopologyPath(root + "/sys/devices/system/cpu/cpu"
                                                 + std::to_string(c) + "/topology/");
        std::string       packageString, coreString;
        std::getline(std::ifstream(logicalCpuTopologyPath + "physical_package_id"), packageString);
        std::getline(std::ifstream(logicalCpuTopologyPath + "core_id"), coreString);

        if (packageString.empty() || coreString.empty())
        {
            continue; // Failsafe so we don't try to parse integers for incorrectly configured cpus
        }
        const int packageId = std::stoi(packageString);
        const int coreId    = std::stoi(coreString);

        // We don't have the rank of the PU in the core, but since that won't be used when we
        // make the translation later, we just set it to -1 here.
        logicalProcessors.push_back({ packageId, coreId, -1, c });
    }

    if (!logicalProcessors.empty())
    {
        translateCpuInfoLogicalProcessorsToMachine(logicalProcessors, machine);
        return HardwareTopology::SupportLevel::Basic;
    }
    else
    {
        return HardwareTopology::SupportLevel::None;
    }
}

/*! \brief Find the active cgroup path for the current process
 *
 * \param mountPoints List of potential mount points
 * \param subGroups   List of potential subgroup paths
 * \param root        Optional mock root for creating topologies from
 *                    saved data. If provided, all paths will be relative
 *                    to this mock root.
 *
 * This routine will attempt to detect the path to the cgroup that's
 * active for the present process by searching for our own process id in
 * the cgroups.procs file in each combination of mountPoint and subGroup.
 *
 * \return Path to the active cgroup directory, or empty string if none found.
 */
std::string findCgroupPath(const std::vector<std::string>& mountPoints,
                           const std::vector<std::string>& subGroups,
                           const std::string&              root)
{
    // We read the pid from /proc/self/stat instead of calling getpid(),
    // so we can mock this properly in testing.
    std::ifstream procSelfStatStream(root + "/proc/self/stat");
    std::string   statString;
    std::getline(procSelfStatStream, statString);
    const int myPid = std::atoi(statString.c_str()); // will be 0 if file was not found
    for (const std::string_view mountPoint : mountPoints)
    {
        for (const std::string_view subGroup : subGroups)
        {
            std::string path = root;
            path.append(mountPoint).append(subGroup);
            std::ifstream cgroupProcsStream(path + "/cgroup.procs");
            std::string   s;
            while (std::getline(cgroupProcsStream, s))
            {
                const int thisPid = std::atoi(s.c_str());
                if (thisPid == myPid)
                {
                    return path;
                }
                else if (thisPid > myPid)
                {
                    break; // no point reading further in this file
                }
            }
        }
    }
    return "";
}

/*! \brief Parse cpu limits from cgroup version 1 file system on Linux
 *
 * \param mountPoints Top-level mount point for cgroups1
 * \param root        Optional mock root for creating topologies from
 *                    saved data. If provided, all paths will be relative
 *                    to this mock root.
 *
 * \return Allowed CPU limit. Note that this is often larger than 1,
 *                    meaning the limit is larger than 1 thread.
 */
float parseCgroup1CpuLimit(const std::vector<std::string>& mountPoints, const std::string& root = "")
{
    std::vector<std::string> subGroups;

    std::string line;
    bool        found = false;

    std::ifstream cgroupStream(root + "/proc/self/cgroup");
    while (!found && std::getline(cgroupStream, line))
    {
        std::istringstream       lineStream(line);
        std::vector<std::string> columns;
        std::string              s;
        while (std::getline(lineStream, s, ':'))
        {
            columns.emplace_back(s);
        }

        if (columns.size() < 3)
        {
            continue;
        }
        if (columns[1] == "cpu,cpuacct" || columns[1] == "cpu")
        {
            subGroups.emplace_back(columns[2]);
            found = true;
        }
    }
    // If the file above could not be read, the subGroups vector will simply be empty

    // Include the top-level group to handle bugs (?) for containers in Ubuntu-18.04
    if (subGroups.empty() || subGroups[subGroups.size() - 1] != "/")
    {
        subGroups.emplace_back("/");
    }

    const std::string cgroupPath = findCgroupPath(mountPoints, subGroups, root);
    if (!cgroupPath.empty())
    {
        // Reach one value from each of two files.
        std::string quotaString, periodString;
        std::getline(std::ifstream(cgroupPath + "/cpu.cfs_quota_us"), quotaString);
        std::getline(std::ifstream(cgroupPath + "/cpu.cfs_period_us"), periodString);
        const int quota  = std::atoi(quotaString.c_str()); // will be 0 if files did not exist
        const int period = std::atoi(periodString.c_str());
        // If conversions could not be done, quota & period will be 0. Quota -1 means no limit.
        if (quota > 0 && period > 0)
        {
            return static_cast<float>(quota) / static_cast<float>(period);
        }
    }
    return -1.0; // no limit found
}

/*! \brief Parse cpu limits from cgroup version 2 file system on Linux
 *
 * \param mountPoints Top-level mount point for cgroups2
 * \param root        Optional mock root for creating topologies from
 *                    saved data. If provided, all paths will be relative
 *                    to this mock root.
 *
 * \return Allowed CPU limit. Note that this is often larger than 1,
 *                    meaning the limit is larger than 1 thread.
 */
float parseCgroup2CpuLimit(const std::vector<std::string>& mountPoints, const std::string& root = "")
{
    std::vector<std::string> subGroups;

    std::string line;
    bool        found = false;

    std::ifstream cgroupStream(root + "/proc/self/cgroup");
    while (!found && std::getline(cgroupStream, line))
    {
        std::istringstream       lineStream(line);
        std::vector<std::string> columns;
        std::string              s;
        while (std::getline(lineStream, s, ':'))
        {
            columns.emplace_back(s);
        }

        if (columns.size() < 3)
        {
            continue;
        }
        // For cgroups2, we are looking for the blank top-level unified subgroup
        if (columns[1].empty())
        {
            subGroups.push_back(columns[2]);
            found = true;
        }
    }

    const std::string cgroupPath = findCgroupPath(mountPoints, subGroups, root);
    if (!cgroupPath.empty())
    {
        // Read two values from one file
        std::ifstream cpuMaxStream(cgroupPath + "/cpu.max");
        std::string   quotaString, periodString;

        if (std::getline(cpuMaxStream, quotaString, ' ') && std::getline(cpuMaxStream, periodString, ' '))
        {
            if (quotaString == "max")
            {
                return -1.0; // no limit
            }
            else
            {
                const int quota = std::atoi(quotaString.c_str()); // will not throw, but might return 0
                const int period = std::atoi(periodString.c_str());
                // If conversions could not be done, quota & period will be 0.
                if (quota != 0 && period > 0)
                {
                    return static_cast<float>(quota) / static_cast<float>(period);
                }
            }
        }
    }
    return -1.0; // no limit found
}

/*\! brief Find the practical cpu limit for the process
 *
 * \param root  Optional mock root for creating topologies from
 *              saved data. If provided, all paths will be relative
 *              to this mock root.
 *
 * \return Allowed CPU limit. Note that this is often larger than 1,
 *                    meaning the limit is larger than 1 thread.
 *
 * This function will attempt to detect the actual allowed CPU limit
 * from Linux kernel filesystems.
 */
float detectCpuLimit(const std::string& root = "")
{
    float cpuLimit = -1;

    std::string              line;
    bool                     found = false;
    std::vector<std::string> cgroups1Mounts;
    std::vector<std::string> cgroups2Mounts;

    // if /etc/mtab isn't present, std::getline will return 0.
    std::ifstream mtabStream(root + "/etc/mtab");
    while (!found && std::getline(mtabStream, line))
    {
        std::istringstream       lineStream(line);
        std::vector<std::string> columns{ std::istream_iterator<std::string>(lineStream),
                                          std::istream_iterator<std::string>() };

        if (columns.size() < 3)
        {
            continue;
        }
        // Column 1 is the fs type (we look for cgroup), column 1 the mount path
        if (columns[2] == "cgroup2")
        {
            // cgroup2 uses a unified mount, so there will only be one entry
            cgroups2Mounts.emplace_back(columns[1]);
        }
        // If token is cgroup1, look for cpu or cpu,cpuacct at end of second token
        else if (columns[2] == "cgroup")
        {
            const std::size_t pos  = columns[1].find_last_of('/');
            const std::string tail = columns[1].substr(pos + 1);
            if (tail == "cpu,cpuacct" || tail == "cpu")
            {
                cgroups1Mounts.emplace_back(columns[1]);
            }
        }
    }

    // Try cgroups2 first
    if (!cgroups2Mounts.empty())
    {
        cpuLimit = parseCgroup2CpuLimit(cgroups2Mounts, root);
    }

    // If we did not find any limits in cgroups2, try cgroups1 interface
    if (cpuLimit < 0 && !cgroups1Mounts.empty())
    {
        cpuLimit = parseCgroup1CpuLimit(cgroups1Mounts, root);
    }
    return cpuLimit;
}


/*\! brief Set the max number of threads in hardware topology
 *
 * The recommended number of threads to start depends on a number
 * of factors.
 * - First, if an explicit CPU limit has been set, we
 *   only round that upwards and use that as the thread limit.
 * - Without cpu limit set, if the support level means
 *   we have information about the packages/cores/processingunits
 *   we are allowed to run on, the total number of logical
 *   processing units is a natural choice for the thread limit.
 * - If we neither have any cpu limit nor topology information,
 *   we try to set the thread limit from the logical processor
 *   count on the hardware. The drawback with this is that it
 *   includes processors we might not be allowed to use, but we
 *   don't really have any alternative.
 *   Note that we provide this as an integer parameter rather
 *   than calling the detection from this routine, so that we
 *   can also use this routine when creating mock topologies.
 *
 * \param cpuLimit      Detected cpu limit (-1 means no limit)
 * \param topologyCpus  Number of allowed processors in topology
 * \param systemCpus    Total number of processors in system
 *
 * \return Recommended max number of threads to start.
 */
int setMaxThreads(float cpuLimit, int topologyCpus, int systemCpus)
{
    if (cpuLimit > 0)
    {
        return std::ceil(cpuLimit);
    }
    else if (topologyCpus > 0)
    {
        return topologyCpus;
    }
    else
    {
        return systemCpus;
    }
}

} // namespace


// static
HardwareTopology HardwareTopology::detect()
{
    HardwareTopology result;

#if GMX_USE_HWLOC
    result.supportLevel_ = parseHwLoc(&result.machine_, &result.isThisSystem_);
#endif
    // If something went wrong in hwloc (or if it was not present) we
    // first try to use CpuInfo, since this avoids file system access
    if (result.supportLevel_ < SupportLevel::Basic)
    {
        result.supportLevel_ = parseCpuInfo(&result.machine_);
    }
#ifdef __linux__
    // If neither hwloc nor CpuInfo parsing worked, try to parse topology from
    // the file system on Linux.
    if (result.supportLevel_ < SupportLevel::Basic)
    {
        result.supportLevel_ = parseSysFsCpuTopology(&result.machine_);
    }
#endif

    result.cpuLimit_   = detectCpuLimit();
    result.maxThreads_ = setMaxThreads(
            result.cpuLimit_, result.machine_.logicalProcessors.size(), detectLogicalProcessorCount());

    // If we have some hunch that more than 1 thread should be started, set the supportlevel to reflect this
    if (result.supportLevel_ == SupportLevel::None && result.maxThreads_ > 1)
    {
        result.supportLevel_ = SupportLevel::LogicalProcessorCount;
    }
    return result;
}

HardwareTopology::Machine::Machine()
{
    numa.baseLatency        = 0.0;
    numa.maxRelativeLatency = 0.0;
}

HardwareTopology::HardwareTopology() :
    supportLevel_(SupportLevel::None), machine_(), isThisSystem_(true), cpuLimit_(1.0)
{
}

HardwareTopology::HardwareTopology(int logicalProcessorCount) :
    supportLevel_(SupportLevel::LogicalProcessorCount),
    machine_(),
    isThisSystem_(false),
    cpuLimit_(logicalProcessorCount),
    maxThreads_(logicalProcessorCount)
{
    if (logicalProcessorCount <= 0)
    {
        cpuLimit_     = -1.0;
        supportLevel_ = SupportLevel::None;
    }
}

HardwareTopology::HardwareTopology(const std::map<int, std::array<int, 3>>& logicalProcessorIdMap,
                                   const std::string&                       filesystemRoot) :
    supportLevel_(SupportLevel::None), machine_(), isThisSystem_(false)
{
    // Create mock topology from saved APIC/core id data and (optionally) saved cgroups files.
    // Start by translating the map into the internal CpuInfo data structure
    // that we do not want to expose in the external class interface/definition.

    std::vector<CpuInfo::LogicalProcessor> logicalProcessors;

    for (const auto& [osId, apic] : logicalProcessorIdMap)
    {
        if (osId >= 0 && apic[0] >= 0 && apic[1] >= 0 && apic[2] >= 0)
        {
            logicalProcessors.push_back({ apic[0], apic[1], apic[2], osId });
        }
        else
        {
            // If any index is negative we return a topology with SupportLevel::None
            return;
        }
    }

    // Now we follow same logic as default detection
    if (!logicalProcessors.empty())
    {
        translateCpuInfoLogicalProcessorsToMachine(logicalProcessors, &machine_);
        supportLevel_ = HardwareTopology::SupportLevel::Basic;
    }

    cpuLimit_ = detectCpuLimit(filesystemRoot);
    // Since we are mocking a topology, use 1 for the system cpu count rather than detect it
    maxThreads_ = setMaxThreads(cpuLimit_, machine_.logicalProcessors.size(), 1);

    // If we have some hunch that more than 1 thread should be started, set the supportlevel to reflect this
    if (supportLevel_ == SupportLevel::None && maxThreads_ > 1)
    {
        supportLevel_ = SupportLevel::LogicalProcessorCount;
    }
}

HardwareTopology::HardwareTopology(const std::string&      filesystemRoot,
                                   const std::vector<int>& allowedProcessors) :
    supportLevel_(SupportLevel::None), machine_(), isThisSystem_(false)
{
    // Create mock topology by parsing saved sys/fs and (optionally) cgroups files
    supportLevel_ = parseSysFsCpuTopology(&machine_, filesystemRoot, allowedProcessors);

    cpuLimit_ = detectCpuLimit(filesystemRoot);
    // Since we are mocking a topology, use 1 for the system cpu count rather than detect it
    maxThreads_ = setMaxThreads(cpuLimit_, machine_.logicalProcessors.size(), 1);

    // If we have some hunch that more than 1 thread should be started, set the supportlevel to reflect this
    if (supportLevel_ == SupportLevel::None && maxThreads_ > 1)
    {
        supportLevel_ = SupportLevel::LogicalProcessorCount;
    }
}

} // namespace gmx
