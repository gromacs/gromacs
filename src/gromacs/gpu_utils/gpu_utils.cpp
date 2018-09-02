/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017,2018, by the GROMACS development team, led by
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
 *  \brief Function definitions for non-GPU builds
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gpu_utils.h"

#include "config.h"

#include <cassert>

#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/smalloc.h"

#if !GMX_GPU
/*! \brief Set allocation functions used by the GPU host
 *
 * Since GPU support is not configured, there is no host memory to
 * allocate. */
void gpu_set_host_malloc_and_free(bool /*unused*/,
                                  gmx_host_alloc_t **nb_alloc,
                                  gmx_host_free_t  **nb_free)
{
    *nb_alloc = nullptr;
    *nb_free  = nullptr;
}

int gpu_info_get_stat(const gmx_gpu_info_t & /*unused*/, int /*unused*/)
{
    return egpuNonexistent;
}
#endif

void free_gpu_info(const gmx_gpu_info_t *gpu_info)
{
    sfree(static_cast<void*>(gpu_info->gpu_dev)); //circumvent is_pod check in sfree
}

std::vector<int> getCompatibleGpus(const gmx_gpu_info_t &gpu_info)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<int> compatibleGpus;
    compatibleGpus.reserve(gpu_info.n_dev);
    for (int i = 0; i < gpu_info.n_dev; i++)
    {
        assert(gpu_info.gpu_dev);
        if (gpu_info_get_stat(gpu_info, i) == egpuCompatible)
        {
            compatibleGpus.push_back(i);
        }
    }
    return compatibleGpus;
}

const char *getGpuCompatibilityDescription(const gmx_gpu_info_t &gpu_info,
                                           int                   index)
{
    return (index >= gpu_info.n_dev ?
            gpu_detect_res_str[egpuNonexistent] :
            gpu_detect_res_str[gpu_info_get_stat(gpu_info, index)]);
}
