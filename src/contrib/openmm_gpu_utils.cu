/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010,2012 The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include <stdio.h>
#include <stdlib.h>

#include "smalloc.h"
#include "string2.h"
#include "types/hw_info.h"

#include "openmm_gpu_utils.h"
#include "../../src/gmxlib/cuda_tools/cudautils.cuh"

// TODO put this list into an external file and include it so that the list is easily accessible
/*! List of supported GPUs. */
static const char * const SupportedGPUs[] = {
    /* GT400 */
    "Geforce GTX 480",
    "Geforce GTX 470",
    "Geforce GTX 465",
    "Geforce GTX 460",

    "Tesla C2070",
    "Tesla C2050",
    "Tesla S2070",
    "Tesla S2050",
    "Tesla M2070",
    "Tesla M2050",

    "Quadro 5000",
    "Quadro 6000",

    /* GT200 */
    "Geforce GTX 295",
    "Geforce GTX 285",
    "Geforce GTX 280",
    "Geforce GTX 275",
    "Geforce GTX 260",
    "GeForce GTS 250",
    "GeForce GTS 150",

    "GeForce GTX 285M",
    "GeForce GTX 280M",

    "Tesla S1070",
    "Tesla C1060",
    "Tesla M1060",

    "Quadro FX 5800",
    "Quadro FX 4800",
    "Quadro CX",
    "Quadro Plex 2200 D2",
    "Quadro Plex 2200 S4",

    /* G90 */
    "GeForce 9800 G", /* GX2, GTX, GTX+, GT */
    "GeForce 9800M GTX",

    "Quadro FX 4700",
    "Quadro Plex 2100 D4"
};

/*! Number of supported GPUs */
#define NB_GPUS (sizeof(SupportedGPUs)/sizeof(SupportedGPUs[0]))

FUNC_QUALIFIER
gmx_bool is_gmx_openmm_supported_gpu(int dev_id, char *gpu_name) FUNC_TERM_INT

/*! 
 * \brief Checks whether the GPU with the given name is supported in Gromacs-OpenMM.
 * 
 * \param[in] gpu_name  the name of the CUDA device
 * \returns             TRUE if the device is supported, otherwise FALSE
 */
static bool is_gmx_openmm_supported_gpu_name(char *gpuName)
{
    size_t i;
    for (i = 0; i < NB_GPUS; i++)
    {
        trim(gpuName);
        if (gmx_strncasecmp(gpuName, SupportedGPUs[i], strlen(SupportedGPUs[i])) == 0)
            return 1;
    }
    return 0;
}

/*! \brief Checks whether the GPU with the given device id is supported in Gromacs-OpenMM.
 *
 * \param[in] dev_id    the device id of the GPU or -1 if the device has already been selected
 * \param[out] gpu_name Set to contain the name of the CUDA device, if NULL passed, no device name is set. 
 * \returns             TRUE if the device is supported, otherwise FALSE
 * 
 */
gmx_bool is_gmx_openmm_supported_gpu(int dev_id, char *gpu_name)
{
    cudaDeviceProp dev_prop;

    if (debug) fprintf(debug, "Checking compatibility with device #%d, %s\n", dev_id, gpu_name);

    if (do_sanity_checks(dev_id, &dev_prop) != 0)
        return -1;

    if (gpu_name != NULL)
    { 
        strcpy(gpu_name, dev_prop.name);
    }
    return is_gmx_openmm_supported_gpu_name(dev_prop.name);
}


