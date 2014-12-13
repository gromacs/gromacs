/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 *  \brief Defines functions that support JIT compilation (e.g. for OpenCL)
 *
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_mdlib
 */
#include "gmxpre.h"

#include <stdlib.h>

#include <cassert>

#include "gromacs/gmxlib/gpu_utils/gpu_utils.h"
#include "gromacs/gmxlib/gpu_utils/ocl_compiler.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/nbnxn_gpu_jit_support.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

#include "nbnxn_ocl_types.h"

//! This function is documented in the header file
void
nbnxn_ocl_map_interaction_types_to_gpu_kernel_flavors(const interaction_const_t *ic,
                                                      int                       *gpu_eeltype,
                                                      int                       *gpu_vdwtype)
{
    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodNONE:
            case eintmodPOTSHIFT:
                *gpu_vdwtype = evdwOclCUT;
                break;
            case eintmodFORCESWITCH:
                *gpu_vdwtype = evdwOclFSWITCH;
                break;
            case eintmodPOTSWITCH:
                *gpu_vdwtype = evdwOclPSWITCH;
                break;
            default:
                gmx_incons("The requested VdW interaction modifier is not implemented in the GPU accelerated kernels!");
                break;
        }
    }
    else if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == ljcrGEOM)
        {
            *gpu_vdwtype = evdwOclEWALDGEOM;
        }
        else
        {
            *gpu_vdwtype = evdwOclEWALDLB;
        }
    }
    else
    {
        gmx_incons("The requested VdW type is not implemented in the GPU accelerated kernels!");
    }

    if (ic->eeltype == eelCUT)
    {
        *gpu_eeltype = eelOclCUT;
    }
    else if (EEL_RF(ic->eeltype))
    {
        *gpu_eeltype = eelOclRF;
    }
    else if ((EEL_PME(ic->eeltype) || ic->eeltype == eelEWALD))
    {
        /* Initially rcoulomb == rvdw, so it's surely not twin cut-off. */
        *gpu_eeltype = nbnxn_gpu_pick_ewald_kernel_type(false);
    }
    else
    {
        /* Shouldn't happen, as this is checked when choosing Verlet-scheme */
        gmx_incons("The requested electrostatics type is not implemented in the GPU accelerated kernels!");
    }
}

/*! \brief Compiles nbnxn kernels for OpenCL GPU given by \p mygpu
 *
 * \param[in]  mygpu          Index of the GPU to initialize
 * \param[out] result_str     The message related to the error that occurred
 *                            during the initialization (if there was any).
 * \param[in] gpu_info        GPU info of all detected devices in the system.
 * \param[in] gpu_opt         Options for using the GPUs in gpu_info
 * \param[in] ic              Interaction constants used to determine which kernels will be used
 * \param[in] bOclDoFastGen   If true, only the requested kernels are compiled, significantly reducing
 * the total compilatin time. If false, all OpenCL kernels are compiled.
 * \returns                   true if no error occurs during initialization.
 */
static gmx_bool
nbnxn_ocl_compile_kernels_inner(int                        mygpu,
                                char                      *result_str,
                                const gmx_gpu_info_t      *gpu_info,
                                const gmx_gpu_opt_t       *gpu_opt,
                                const interaction_const_t *ic,
                                const gmx_bool             bOclDoFastGen)
{
    cl_context_properties     context_properties[3];
    struct gmx_device_info_t *selected_ocl_gpu;
    cl_platform_id            platform_id;
    cl_device_id              device_id;
    cl_context                context;
    cl_program                program;
    cl_int                    cl_error;

    gmx_algo_family_t         gmx_algo_family;

    nbnxn_ocl_map_interaction_types_to_gpu_kernel_flavors(ic,
                                                          &gmx_algo_family.eeltype,
                                                          &gmx_algo_family.vdwtype);
    assert(gmx_algo_family.eeltype < eelOclNR);
    assert(gmx_algo_family.vdwtype < evdwOclNR);

    assert(gpu_info);
    assert(result_str);

    result_str[0] = 0;

    selected_ocl_gpu = gpu_info->gpu_dev + gpu_opt->dev_use[mygpu];
    platform_id      = selected_ocl_gpu->ocl_gpu_id.ocl_platform_id;
    device_id        = selected_ocl_gpu->ocl_gpu_id.ocl_device_id;

    context_properties[0] = CL_CONTEXT_PLATFORM;
    context_properties[1] = (cl_context_properties)platform_id;
    context_properties[2] = 0;

    context = clCreateContext(context_properties, 1, &device_id, NULL, NULL, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        sprintf(result_str, "OpenCL error %d", cl_error);
        return FALSE;
    }

    cl_error =
        ocl_compile_program(_default_source_,
                            _auto_vendor_kernels_,
                            &gmx_algo_family,
                            bOclDoFastGen,
                            result_str,
                            context,
                            device_id,
                            selected_ocl_gpu->vendor_e,
                            &program
                            );
    if (cl_error != CL_SUCCESS)
    {
        sprintf(result_str, "OpenCL error %d", cl_error);
        return FALSE;
    }

    selected_ocl_gpu->context = context;
    selected_ocl_gpu->program = program;

    return TRUE;
}

/*! \brief Compiles nbnxn kernels for OpenCL GPU given by \p mygpu
 *
 * With OpenCL, a call to this function must precede nbnxn_gpu_init().
 *
 * Doing bOclDoFastGen means only the requested kernels are compiled,
 * significantly reducing the total compilation time. If false, all
 * OpenCL kernels are compiled.
 *
 * A fatal error results if compilation fails.
 *
 * \param[in]  mygpu          The intra-node PP rank, which is the index into the list of GPUs to use
 * \param[in]  rank           MPI rank (for error reporting)
 * \param[in] gpu_info        GPU info of all detected devices in the system.
 * \param[in] gpu_opt         Options for using the GPUs in gpu_info
 * \param[in] ic              Interaction constants used to determine which kernels will be used
 */
void
nbnxn_gpu_compile_kernels(int                        mygpu,
                          int                        rank,
                          const gmx_gpu_info_t      *gpu_info,
                          const gmx_gpu_opt_t       *gpu_opt,
                          const interaction_const_t *ic)
{
    char     gpu_err_str[STRLEN];
    gmx_bool bOclDoFastGen;

    bOclDoFastGen = (getenv("GMX_OCL_NOFASTGEN") == NULL);
    if (bOclDoFastGen)
    {
        bOclDoFastGen = FALSE;
#ifndef NDEBUG
        printf("\nFastGen temporary disabled. All kernel flavours will be generated.");
#endif
    }

    if (!nbnxn_ocl_compile_kernels_inner(mygpu, gpu_err_str,
                                         gpu_info, gpu_opt, ic, bOclDoFastGen))
    {
        gmx_fatal(FARGS, "On rank %d failed to compile NBNXN kernels for GPU #%s: %s",
                  rank,
                  get_ocl_gpu_device_name(gpu_info, gpu_opt, mygpu),
                  gpu_err_str);
    }
}
