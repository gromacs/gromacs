/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/nbnxn_gpu_jit_support.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

#include "nbnxn_ocl_types.h"

/*! \brief Stringifies the input argument
 */
#define STRINGIFY_PARAM(c) #c

/*! \brief Stringifies the result of expansion of a macro argument
 */
#define STRINGIFY_MACRO(c) STRINGIFY_PARAM(c)

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
 * With OpenCL, a call to this function must precede nbnxn_gpu_init().
 *
 * Doing bOclDoFastGen means only the requested kernels are compiled,
 * significantly reducing the total compilation time. If false, all
 * OpenCL kernels are compiled.
 *
 * A fatal error results if compilation fails.
 *
 * \param[inout] nb  Manages OpenCL non-bonded calculations; compiled kernels returned in dev_info members
 * \param[in]    ic  Interaction constants used to determine which kernels will be used
 */
void
nbnxn_gpu_compile_kernels(gmx_nbnxn_ocl_t           *nb,
                          const interaction_const_t *ic)
{
    char                      gpu_err_str[STRLEN];
    gmx_bool                  bOclDoFastGen = TRUE;
    cl_device_id              device_id;
    cl_context                context;
    cl_program                program;
    cl_int                    cl_error;
    char                      runtime_consts[256];
    gmx_algo_family_t         gmx_algo_family;

    if (getenv("GMX_OCL_NOFASTGEN") != NULL)
    {
        bOclDoFastGen = FALSE;
    }

    nbnxn_ocl_map_interaction_types_to_gpu_kernel_flavors(ic,
                                                          &gmx_algo_family.eeltype,
                                                          &gmx_algo_family.vdwtype);
    assert(gmx_algo_family.eeltype < eelOclNR);
    assert(gmx_algo_family.vdwtype < evdwOclNR);

    device_id        = nb->dev_info->ocl_gpu_id.ocl_device_id;
    context          = nb->dev_info->context;

    sprintf(runtime_consts,
            "-DCENTRAL=%d -DNBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER=%d -DNBNXN_GPU_CLUSTER_SIZE=%d -DNBNXN_GPU_JGROUP_SIZE=%d -DNBNXN_AVOID_SING_R2_INC=%s",
            CENTRAL,                                    /* Defined in ishift.h */
            NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER,        /* Defined in nbnxn_consts.h */
            NBNXN_GPU_CLUSTER_SIZE,                     /* Defined in nbnxn_consts.h */
            NBNXN_GPU_JGROUP_SIZE,                      /* Defined in nbnxn_consts.h */
            STRINGIFY_MACRO(NBNXN_AVOID_SING_R2_INC)    /* Defined in nbnxn_consts.h */
                                                        /* NBNXN_AVOID_SING_R2_INC passed as string to avoid
                                                           floating point representation problems with sprintf */
            );

    cl_error =
        ocl_compile_program(default_source,
                            auto_vendor_kernels,
                            &gmx_algo_family,
                            bOclDoFastGen,
                            gpu_err_str,
                            context,
                            device_id,
                            nb->dev_info->vendor_e,
                            &program,
                            runtime_consts
                            );
    if (cl_error != CL_SUCCESS)
    {
        gmx_fatal(FARGS, "Failed to compile NBNXN kernels for GPU #%s: %s",
                  nb->dev_info->device_name,
                  gpu_err_str);
    }

    nb->dev_info->program = program;
}
