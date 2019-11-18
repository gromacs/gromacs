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
 *  \brief Define OpenCL implementation of nbnxm_gpu_data_mgmt.h
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

// TODO We would like to move this down, but the way gmx_nbnxn_gpu_t
//      is currently declared means this has to be before gpu_types.h
#include "nbnxm_ocl_types.h"

// TODO Remove this comment when the above order issue is resolved
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/gpu_jit_support.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "nbnxm_ocl_internal.h"

namespace Nbnxm
{

/*! \brief Copies of values from cl_driver_diagnostics_intel.h,
 * which isn't guaranteed to be available. */
/**@{*/
#define CL_CONTEXT_SHOW_DIAGNOSTICS_INTEL 0x4106
#define CL_CONTEXT_DIAGNOSTICS_LEVEL_GOOD_INTEL 0x1
#define CL_CONTEXT_DIAGNOSTICS_LEVEL_BAD_INTEL 0x2
#define CL_CONTEXT_DIAGNOSTICS_LEVEL_NEUTRAL_INTEL 0x4
/**@}*/

/*! \brief This parameter should be determined heuristically from the
 * kernel execution times
 *
 * This value is best for small systems on a single AMD Radeon R9 290X
 * (and about 5% faster than 40, which is the default for CUDA
 * devices). Larger simulation systems were quite insensitive to the
 * value of this parameter.
 */
static unsigned int gpu_min_ci_balanced_factor = 50;


/*! \brief Returns true if LJ combination rules are used in the non-bonded kernels.
 *
 * Full doc in nbnxn_ocl_internal.h */
bool useLjCombRule(int vdwType)
{
    return (vdwType == evdwOclCUTCOMBGEOM || vdwType == evdwOclCUTCOMBLB);
}

/*! \brief Tabulates the Ewald Coulomb force and initializes the size/scale
 * and the table GPU array.
 *
 * If called with an already allocated table, it just re-uploads the
 * table.
 */
static void init_ewald_coulomb_force_table(const EwaldCorrectionTables&     tables,
                                           cl_nbparam_t*                    nbp,
                                           const gmx_device_runtime_data_t* runData)
{
    cl_mem coul_tab;

    cl_int cl_error;

    if (nbp->coulomb_tab_climg2d != nullptr)
    {
        freeDeviceBuffer(&(nbp->coulomb_tab_climg2d));
    }

    /* Switched from using textures to using buffers */
    // TODO: decide which alternative is most efficient - textures or buffers.
    /*
       cl_image_format array_format;

       array_format.image_channel_data_type = CL_FLOAT;
       array_format.image_channel_order     = CL_R;

       coul_tab = clCreateImage2D(runData->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
       &array_format, tabsize, 1, 0, ftmp, &cl_error);
     */

    coul_tab = clCreateBuffer(
            runData->context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
            tables.tableF.size() * sizeof(cl_float), const_cast<real*>(tables.tableF.data()), &cl_error);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

    nbp->coulomb_tab_climg2d = coul_tab;
    nbp->coulomb_tab_scale   = tables.scale;
}


/*! \brief Initializes the atomdata structure first time, it only gets filled at
    pair-search.
 */
static void init_atomdata_first(cl_atomdata_t* ad, int ntypes, gmx_device_runtime_data_t* runData)
{
    cl_int cl_error;

    ad->ntypes = ntypes;

    /* An element of the shift_vec device buffer has the same size as one element
       of the host side shift_vec buffer. */
    ad->shift_vec_elem_size = sizeof(*nbnxn_atomdata_t::shift_vec.data());

    ad->shift_vec = clCreateBuffer(runData->context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                   SHIFTS * ad->shift_vec_elem_size, nullptr, &cl_error);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
    ad->bShiftVecUploaded = CL_FALSE;

    /* An element of the fshift device buffer has the same size as one element
       of the host side fshift buffer. */
    ad->fshift_elem_size = sizeof(*cl_nb_staging_t::fshift);

    ad->fshift = clCreateBuffer(runData->context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
                                SHIFTS * ad->fshift_elem_size, nullptr, &cl_error);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

    ad->e_lj = clCreateBuffer(runData->context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
                              sizeof(float), nullptr, &cl_error);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

    ad->e_el = clCreateBuffer(runData->context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
                              sizeof(float), nullptr, &cl_error);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

    /* initialize to nullptr pointers to data that is not allocated here and will
       need reallocation in nbnxn_gpu_init_atomdata */
    ad->xq = nullptr;
    ad->f  = nullptr;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    ad->natoms = -1;
    ad->nalloc = -1;
}

/*! \brief Copies all parameters related to the cut-off from ic to nbp
 */
static void set_cutoff_parameters(cl_nbparam_t* nbp, const interaction_const_t* ic, const PairlistParams& listParams)
{
    nbp->ewald_beta        = ic->ewaldcoeff_q;
    nbp->sh_ewald          = ic->sh_ewald;
    nbp->epsfac            = ic->epsfac;
    nbp->two_k_rf          = 2.0 * ic->k_rf;
    nbp->c_rf              = ic->c_rf;
    nbp->rvdw_sq           = ic->rvdw * ic->rvdw;
    nbp->rcoulomb_sq       = ic->rcoulomb * ic->rcoulomb;
    nbp->rlistOuter_sq     = listParams.rlistOuter * listParams.rlistOuter;
    nbp->rlistInner_sq     = listParams.rlistInner * listParams.rlistInner;
    nbp->useDynamicPruning = listParams.useDynamicPruning;

    nbp->sh_lj_ewald   = ic->sh_lj_ewald;
    nbp->ewaldcoeff_lj = ic->ewaldcoeff_lj;

    nbp->rvdw_switch      = ic->rvdw_switch;
    nbp->dispersion_shift = ic->dispersion_shift;
    nbp->repulsion_shift  = ic->repulsion_shift;
    nbp->vdw_switch       = ic->vdw_switch;
}

/*! \brief Returns the kinds of electrostatics and Vdw OpenCL
 *  kernels that will be used.
 *
 * Respectively, these values are from enum eelOcl and enum
 * evdwOcl. */
static void map_interaction_types_to_gpu_kernel_flavors(const interaction_const_t* ic,
                                                        int                        combRule,
                                                        int*                       gpu_eeltype,
                                                        int*                       gpu_vdwtype)
{
    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodNONE:
            case eintmodPOTSHIFT:
                switch (combRule)
                {
                    case ljcrNONE: *gpu_vdwtype = evdwOclCUT; break;
                    case ljcrGEOM: *gpu_vdwtype = evdwOclCUTCOMBGEOM; break;
                    case ljcrLB: *gpu_vdwtype = evdwOclCUTCOMBLB; break;
                    default:
                        gmx_incons(
                                "The requested LJ combination rule is not implemented in the "
                                "OpenCL GPU accelerated kernels!");
                }
                break;
            case eintmodFORCESWITCH: *gpu_vdwtype = evdwOclFSWITCH; break;
            case eintmodPOTSWITCH: *gpu_vdwtype = evdwOclPSWITCH; break;
            default:
                gmx_incons(
                        "The requested VdW interaction modifier is not implemented in the GPU "
                        "accelerated kernels!");
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
        *gpu_eeltype = nbnxn_gpu_pick_ewald_kernel_type(*ic);
    }
    else
    {
        /* Shouldn't happen, as this is checked when choosing Verlet-scheme */
        gmx_incons(
                "The requested electrostatics type is not implemented in the GPU accelerated "
                "kernels!");
    }
}

/*! \brief Initializes the nonbonded parameter data structure.
 */
static void init_nbparam(cl_nbparam_t*                    nbp,
                         const interaction_const_t*       ic,
                         const PairlistParams&            listParams,
                         const nbnxn_atomdata_t::Params&  nbatParams,
                         const gmx_device_runtime_data_t* runData)
{
    cl_int cl_error;

    set_cutoff_parameters(nbp, ic, listParams);

    map_interaction_types_to_gpu_kernel_flavors(ic, nbatParams.comb_rule, &(nbp->eeltype), &(nbp->vdwtype));

    if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == ljcrGEOM)
        {
            GMX_ASSERT(nbatParams.comb_rule == ljcrGEOM, "Combination rule mismatch!");
        }
        else
        {
            GMX_ASSERT(nbatParams.comb_rule == ljcrLB, "Combination rule mismatch!");
        }
    }
    /* generate table for PME */
    nbp->coulomb_tab_climg2d = nullptr;
    if (nbp->eeltype == eelOclEWALD_TAB || nbp->eeltype == eelOclEWALD_TAB_TWIN)
    {
        GMX_RELEASE_ASSERT(ic->coulombEwaldTables, "Need valid Coulomb Ewald correction tables");
        init_ewald_coulomb_force_table(*ic->coulombEwaldTables, nbp, runData);
    }
    else
    // TODO: improvement needed.
    // The image2d is created here even if eeltype is not eelCuEWALD_TAB or eelCuEWALD_TAB_TWIN
    // because the OpenCL kernels don't accept nullptr values for image2D parameters.
    {
        /* Switched from using textures to using buffers */
        // TODO: decide which alternative is most efficient - textures or buffers.
        /*
           cl_image_format array_format;

           array_format.image_channel_data_type = CL_FLOAT;
           array_format.image_channel_order     = CL_R;

           nbp->coulomb_tab_climg2d = clCreateImage2D(runData->context, CL_MEM_READ_WRITE,
            &array_format, 1, 1, 0, nullptr, &cl_error);
         */

        nbp->coulomb_tab_climg2d = clCreateBuffer(runData->context, CL_MEM_READ_ONLY,
                                                  sizeof(cl_float), nullptr, &cl_error);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
    }

    const int nnbfp      = 2 * nbatParams.numTypes * nbatParams.numTypes;
    const int nnbfp_comb = 2 * nbatParams.numTypes;

    {
        /* Switched from using textures to using buffers */
        // TODO: decide which alternative is most efficient - textures or buffers.
        /*
           cl_image_format array_format;

           array_format.image_channel_data_type = CL_FLOAT;
           array_format.image_channel_order     = CL_R;

           nbp->nbfp_climg2d = clCreateImage2D(runData->context, CL_MEM_READ_ONLY |
           CL_MEM_COPY_HOST_PTR, &array_format, nnbfp, 1, 0, nbat->nbfp, &cl_error);
         */

        nbp->nbfp_climg2d = clCreateBuffer(
                runData->context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
                nnbfp * sizeof(cl_float), const_cast<float*>(nbatParams.nbfp.data()), &cl_error);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

        if (ic->vdwtype == evdwPME)
        {
            /* Switched from using textures to using buffers */
            // TODO: decide which alternative is most efficient - textures or buffers.
            /*  nbp->nbfp_comb_climg2d = clCreateImage2D(runData->context, CL_MEM_READ_WRITE |
               CL_MEM_COPY_HOST_PTR, &array_format, nnbfp_comb, 1, 0, nbat->nbfp_comb, &cl_error);*/
            nbp->nbfp_comb_climg2d = clCreateBuffer(
                    runData->context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
                    nnbfp_comb * sizeof(cl_float), const_cast<float*>(nbatParams.nbfp_comb.data()),
                    &cl_error);
            GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                               ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
        }
        else
        {
            // TODO: improvement needed.
            // The image2d is created here even if vdwtype is not evdwPME because the OpenCL kernels
            // don't accept nullptr values for image2D parameters.
            /* Switched from using textures to using buffers */
            // TODO: decide which alternative is most efficient - textures or buffers.
            /* nbp->nbfp_comb_climg2d = clCreateImage2D(runData->context, CL_MEM_READ_WRITE,
                &array_format, 1, 1, 0, nullptr, &cl_error);*/
            nbp->nbfp_comb_climg2d = clCreateBuffer(runData->context, CL_MEM_READ_ONLY,
                                                    sizeof(cl_float), nullptr, &cl_error);
            GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                               ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
        }
    }
}

//! This function is documented in the header file
void gpu_pme_loadbal_update_param(const nonbonded_verlet_t* nbv, const interaction_const_t* ic)
{
    if (!nbv || !nbv->useGpu())
    {
        return;
    }
    gmx_nbnxn_ocl_t* nb  = nbv->gpu_nbv;
    cl_nbparam_t*    nbp = nb->nbparam;

    set_cutoff_parameters(nbp, ic, nbv->pairlistSets().params());

    nbp->eeltype = nbnxn_gpu_pick_ewald_kernel_type(*ic);

    GMX_RELEASE_ASSERT(ic->coulombEwaldTables, "Need valid Coulomb Ewald correction tables");
    init_ewald_coulomb_force_table(*ic->coulombEwaldTables, nbp, nb->dev_rundata);
}

/*! \brief Initializes the pair list data structure.
 */
static void init_plist(cl_plist_t* pl)
{
    /* initialize to nullptr pointers to data that is not allocated here and will
       need reallocation in nbnxn_gpu_init_pairlist */
    pl->sci   = nullptr;
    pl->cj4   = nullptr;
    pl->imask = nullptr;
    pl->excl  = nullptr;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    pl->na_c          = -1;
    pl->nsci          = -1;
    pl->sci_nalloc    = -1;
    pl->ncj4          = -1;
    pl->cj4_nalloc    = -1;
    pl->nimask        = -1;
    pl->imask_nalloc  = -1;
    pl->nexcl         = -1;
    pl->excl_nalloc   = -1;
    pl->haveFreshList = false;
}

/*! \brief Initializes the timings data structure.
 */
static void init_timings(gmx_wallclock_gpu_nbnxn_t* t)
{
    int i, j;

    t->nb_h2d_t = 0.0;
    t->nb_d2h_t = 0.0;
    t->nb_c     = 0;
    t->pl_h2d_t = 0.0;
    t->pl_h2d_c = 0;
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            t->ktime[i][j].t = 0.0;
            t->ktime[i][j].c = 0;
        }
    }

    t->pruneTime.c        = 0;
    t->pruneTime.t        = 0.0;
    t->dynamicPruneTime.c = 0;
    t->dynamicPruneTime.t = 0.0;
}


//! OpenCL notification callback function
static void CL_CALLBACK ocl_notify_fn(const char* pErrInfo,
                                      const void gmx_unused* private_info,
                                      size_t gmx_unused cb,
                                      void gmx_unused* user_data)
{
    if (pErrInfo != nullptr)
    {
        printf("%s\n", pErrInfo); // Print error/hint
    }
}

/*! \brief Creates context for OpenCL GPU given by \p mygpu
 *
 * A fatal error results if creation fails.
 *
 * \param[inout] runtimeData runtime data including program and context
 * \param[in]    devInfo     device info struct
 * \param[in]    rank        MPI rank (for error reporting)
 */
static void nbnxn_gpu_create_context(gmx_device_runtime_data_t* runtimeData,
                                     const gmx_device_info_t*   devInfo,
                                     int                        rank)
{
    cl_context_properties context_properties[5];
    cl_platform_id        platform_id;
    cl_device_id          device_id;
    cl_context            context;
    cl_int                cl_error;

    assert(runtimeData != nullptr);
    assert(devInfo != nullptr);

    platform_id = devInfo->ocl_gpu_id.ocl_platform_id;
    device_id   = devInfo->ocl_gpu_id.ocl_device_id;

    int i                   = 0;
    context_properties[i++] = CL_CONTEXT_PLATFORM;
    context_properties[i++] = reinterpret_cast<cl_context_properties>(platform_id);
    if (getenv("GMX_OCL_SHOW_DIAGNOSTICS"))
    {
        context_properties[i++] = CL_CONTEXT_SHOW_DIAGNOSTICS_INTEL;
        context_properties[i++] =
                CL_CONTEXT_DIAGNOSTICS_LEVEL_BAD_INTEL | CL_CONTEXT_DIAGNOSTICS_LEVEL_NEUTRAL_INTEL;
    }
    context_properties[i++] = 0; /* Terminates the list of properties */

    context = clCreateContext(context_properties, 1, &device_id, ocl_notify_fn, nullptr, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS, "On rank %d failed to create context for GPU #%s:\n OpenCL error %d: %s",
                  rank, devInfo->device_name, cl_error, ocl_get_error_string(cl_error).c_str());
    }

    runtimeData->context = context;
}

/*! \brief Initializes the OpenCL kernel pointers of the nbnxn_ocl_ptr_t input data structure. */
static cl_kernel nbnxn_gpu_create_kernel(gmx_nbnxn_ocl_t* nb, const char* kernel_name)
{
    cl_kernel kernel;
    cl_int    cl_error;

    kernel = clCreateKernel(nb->dev_rundata->program, kernel_name, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS, "Failed to create kernel '%s' for GPU #%s: OpenCL error %d", kernel_name,
                  nb->dev_info->device_name, cl_error);
    }

    return kernel;
}

/*! \brief Clears nonbonded shift force output array and energy outputs on the GPU.
 */
static void nbnxn_ocl_clear_e_fshift(gmx_nbnxn_ocl_t* nb)
{

    cl_int           cl_error;
    cl_atomdata_t*   adat = nb->atdat;
    cl_command_queue ls   = nb->stream[InteractionLocality::Local];

    size_t local_work_size[3]  = { 1, 1, 1 };
    size_t global_work_size[3] = { 1, 1, 1 };

    cl_int shifts = SHIFTS * 3;

    cl_int arg_no;

    cl_kernel zero_e_fshift = nb->kernel_zero_e_fshift;

    local_work_size[0] = 64;
    // Round the total number of threads up from the array size
    global_work_size[0] = ((shifts + local_work_size[0] - 1) / local_work_size[0]) * local_work_size[0];

    arg_no   = 0;
    cl_error = clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_mem), &(adat->fshift));
    cl_error |= clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_mem), &(adat->e_lj));
    cl_error |= clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_mem), &(adat->e_el));
    cl_error |= clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_uint), &shifts);
    GMX_ASSERT(cl_error == CL_SUCCESS, ocl_get_error_string(cl_error).c_str());

    cl_error = clEnqueueNDRangeKernel(ls, zero_e_fshift, 3, nullptr, global_work_size,
                                      local_work_size, 0, nullptr, nullptr);
    GMX_ASSERT(cl_error == CL_SUCCESS, ocl_get_error_string(cl_error).c_str());
}

/*! \brief Initializes the OpenCL kernel pointers of the nbnxn_ocl_ptr_t input data structure. */
static void nbnxn_gpu_init_kernels(gmx_nbnxn_ocl_t* nb)
{
    /* Init to 0 main kernel arrays */
    /* They will be later on initialized in select_nbnxn_kernel */
    // TODO: consider always creating all variants of the kernels here so that there is no
    // need for late call to clCreateKernel -- if that gives any advantage?
    memset(nb->kernel_ener_noprune_ptr, 0, sizeof(nb->kernel_ener_noprune_ptr));
    memset(nb->kernel_ener_prune_ptr, 0, sizeof(nb->kernel_ener_prune_ptr));
    memset(nb->kernel_noener_noprune_ptr, 0, sizeof(nb->kernel_noener_noprune_ptr));
    memset(nb->kernel_noener_prune_ptr, 0, sizeof(nb->kernel_noener_prune_ptr));

    /* Init pruning kernels
     *
     * TODO: we could avoid creating kernels if dynamic pruning is turned off,
     * but ATM that depends on force flags not passed into the initialization.
     */
    nb->kernel_pruneonly[epruneFirst] = nbnxn_gpu_create_kernel(nb, "nbnxn_kernel_prune_opencl");
    nb->kernel_pruneonly[epruneRolling] =
            nbnxn_gpu_create_kernel(nb, "nbnxn_kernel_prune_rolling_opencl");

    /* Init auxiliary kernels */
    nb->kernel_zero_e_fshift = nbnxn_gpu_create_kernel(nb, "zero_e_fshift");
}

/*! \brief Initializes simulation constant data.
 *
 *  Initializes members of the atomdata and nbparam structs and
 *  clears e/fshift output buffers.
 */
static void nbnxn_ocl_init_const(gmx_nbnxn_ocl_t*                nb,
                                 const interaction_const_t*      ic,
                                 const PairlistParams&           listParams,
                                 const nbnxn_atomdata_t::Params& nbatParams)
{
    init_atomdata_first(nb->atdat, nbatParams.numTypes, nb->dev_rundata);
    init_nbparam(nb->nbparam, ic, listParams, nbatParams, nb->dev_rundata);
}


//! This function is documented in the header file
gmx_nbnxn_ocl_t* gpu_init(const gmx_device_info_t*   deviceInfo,
                          const interaction_const_t* ic,
                          const PairlistParams&      listParams,
                          const nbnxn_atomdata_t*    nbat,
                          const int                  rank,
                          const gmx_bool             bLocalAndNonlocal)
{
    gmx_nbnxn_ocl_t*            nb;
    cl_int                      cl_error;
    cl_command_queue_properties queue_properties;

    assert(ic);

    snew(nb, 1);
    snew(nb->atdat, 1);
    snew(nb->nbparam, 1);
    snew(nb->plist[InteractionLocality::Local], 1);
    if (bLocalAndNonlocal)
    {
        snew(nb->plist[InteractionLocality::NonLocal], 1);
    }

    nb->bUseTwoStreams = static_cast<cl_bool>(bLocalAndNonlocal);

    nb->timers = new cl_timers_t();
    snew(nb->timings, 1);

    /* set device info, just point it to the right GPU among the detected ones */
    nb->dev_info = deviceInfo;
    snew(nb->dev_rundata, 1);

    /* init nbst */
    pmalloc(reinterpret_cast<void**>(&nb->nbst.e_lj), sizeof(*nb->nbst.e_lj));
    pmalloc(reinterpret_cast<void**>(&nb->nbst.e_el), sizeof(*nb->nbst.e_el));
    pmalloc(reinterpret_cast<void**>(&nb->nbst.fshift), SHIFTS * sizeof(*nb->nbst.fshift));

    init_plist(nb->plist[InteractionLocality::Local]);

    /* OpenCL timing disabled if GMX_DISABLE_GPU_TIMING is defined. */
    nb->bDoTime = static_cast<cl_bool>(getenv("GMX_DISABLE_GPU_TIMING") == nullptr);

    /* Create queues only after bDoTime has been initialized */
    if (nb->bDoTime)
    {
        queue_properties = CL_QUEUE_PROFILING_ENABLE;
    }
    else
    {
        queue_properties = 0;
    }

    nbnxn_gpu_create_context(nb->dev_rundata, nb->dev_info, rank);

    /* local/non-local GPU streams */
    nb->stream[InteractionLocality::Local] = clCreateCommandQueue(
            nb->dev_rundata->context, nb->dev_info->ocl_gpu_id.ocl_device_id, queue_properties, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS, "On rank %d failed to create context for GPU #%s: OpenCL error %d", rank,
                  nb->dev_info->device_name, cl_error);
    }

    if (nb->bUseTwoStreams)
    {
        init_plist(nb->plist[InteractionLocality::NonLocal]);

        nb->stream[InteractionLocality::NonLocal] =
                clCreateCommandQueue(nb->dev_rundata->context, nb->dev_info->ocl_gpu_id.ocl_device_id,
                                     queue_properties, &cl_error);
        if (CL_SUCCESS != cl_error)
        {
            gmx_fatal(FARGS, "On rank %d failed to create context for GPU #%s: OpenCL error %d",
                      rank, nb->dev_info->device_name, cl_error);
        }
    }

    if (nb->bDoTime)
    {
        init_timings(nb->timings);
    }

    nbnxn_ocl_init_const(nb, ic, listParams, nbat->params());

    /* Enable LJ param manual prefetch for AMD or Intel or if we request through env. var.
     * TODO: decide about NVIDIA
     */
    nb->bPrefetchLjParam = (getenv("GMX_OCL_DISABLE_I_PREFETCH") == nullptr)
                           && ((nb->dev_info->vendor_e == OCL_VENDOR_AMD)
                               || (nb->dev_info->vendor_e == OCL_VENDOR_INTEL)
                               || (getenv("GMX_OCL_ENABLE_I_PREFETCH") != nullptr));

    /* NOTE: in CUDA we pick L1 cache configuration for the nbnxn kernels here,
     * but sadly this is not supported in OpenCL (yet?). Consider adding it if
     * it becomes supported.
     */
    nbnxn_gpu_compile_kernels(nb);
    nbnxn_gpu_init_kernels(nb);

    /* clear energy and shift force outputs */
    nbnxn_ocl_clear_e_fshift(nb);

    if (debug)
    {
        fprintf(debug, "Initialized OpenCL data structures.\n");
    }

    return nb;
}

/*! \brief Clears the first natoms_clear elements of the GPU nonbonded force output array.
 */
static void nbnxn_ocl_clear_f(gmx_nbnxn_ocl_t* nb, int natoms_clear)
{
    if (natoms_clear == 0)
    {
        return;
    }

    cl_int gmx_used_in_debug cl_error;

    cl_atomdata_t*   atomData = nb->atdat;
    cl_command_queue ls       = nb->stream[InteractionLocality::Local];
    cl_float         value    = 0.0F;

    cl_error = clEnqueueFillBuffer(ls, atomData->f, &value, sizeof(cl_float), 0,
                                   natoms_clear * sizeof(rvec), 0, nullptr, nullptr);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clEnqueueFillBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
}

//! This function is documented in the header file
void gpu_clear_outputs(gmx_nbnxn_ocl_t* nb, bool computeVirial)
{
    nbnxn_ocl_clear_f(nb, nb->atdat->natoms);
    /* clear shift force array and energies if the outputs were
       used in the current step */
    if (computeVirial)
    {
        nbnxn_ocl_clear_e_fshift(nb);
    }

    /* kick off buffer clearing kernel to ensure concurrency with constraints/update */
    cl_int gmx_unused cl_error;
    cl_error = clFlush(nb->stream[InteractionLocality::Local]);
    assert(CL_SUCCESS == cl_error);
}

//! This function is documented in the header file
void gpu_init_pairlist(gmx_nbnxn_ocl_t* nb, const NbnxnPairlistGpu* h_plist, const InteractionLocality iloc)
{
    char sbuf[STRLEN];
    // Timing accumulation should happen only if there was work to do
    // because getLastRangeTime() gets skipped with empty lists later
    // which leads to the counter not being reset.
    bool             bDoTime = ((nb->bDoTime == CL_TRUE) && !h_plist->sci.empty());
    cl_command_queue stream  = nb->stream[iloc];
    cl_plist_t*      d_plist = nb->plist[iloc];

    if (d_plist->na_c < 0)
    {
        d_plist->na_c = h_plist->na_ci;
    }
    else
    {
        if (d_plist->na_c != h_plist->na_ci)
        {
            sprintf(sbuf, "In cu_init_plist: the #atoms per cell has changed (from %d to %d)",
                    d_plist->na_c, h_plist->na_ci);
            gmx_incons(sbuf);
        }
    }

    gpu_timers_t::Interaction& iTimers = nb->timers->interaction[iloc];

    if (bDoTime)
    {
        iTimers.pl_h2d.openTimingRegion(stream);
        iTimers.didPairlistH2D = true;
    }

    // TODO most of this function is same in CUDA and OpenCL, move into the header
    DeviceContext context = nb->dev_rundata->context;

    reallocateDeviceBuffer(&d_plist->sci, h_plist->sci.size(), &d_plist->nsci, &d_plist->sci_nalloc, context);
    copyToDeviceBuffer(&d_plist->sci, h_plist->sci.data(), 0, h_plist->sci.size(), stream,
                       GpuApiCallBehavior::Async, bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);

    reallocateDeviceBuffer(&d_plist->cj4, h_plist->cj4.size(), &d_plist->ncj4, &d_plist->cj4_nalloc, context);
    copyToDeviceBuffer(&d_plist->cj4, h_plist->cj4.data(), 0, h_plist->cj4.size(), stream,
                       GpuApiCallBehavior::Async, bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);

    reallocateDeviceBuffer(&d_plist->imask, h_plist->cj4.size() * c_nbnxnGpuClusterpairSplit,
                           &d_plist->nimask, &d_plist->imask_nalloc, context);

    reallocateDeviceBuffer(&d_plist->excl, h_plist->excl.size(), &d_plist->nexcl,
                           &d_plist->excl_nalloc, context);
    copyToDeviceBuffer(&d_plist->excl, h_plist->excl.data(), 0, h_plist->excl.size(), stream,
                       GpuApiCallBehavior::Async, bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);

    if (bDoTime)
    {
        iTimers.pl_h2d.closeTimingRegion(stream);
    }

    /* need to prune the pair list during the next step */
    d_plist->haveFreshList = true;
}

//! This function is documented in the header file
void gpu_upload_shiftvec(gmx_nbnxn_ocl_t* nb, const nbnxn_atomdata_t* nbatom)
{
    cl_atomdata_t*   adat = nb->atdat;
    cl_command_queue ls   = nb->stream[InteractionLocality::Local];

    /* only if we have a dynamic box */
    if (nbatom->bDynamicBox || !adat->bShiftVecUploaded)
    {
        ocl_copy_H2D_async(adat->shift_vec, nbatom->shift_vec.data(), 0,
                           SHIFTS * adat->shift_vec_elem_size, ls, nullptr);
        adat->bShiftVecUploaded = CL_TRUE;
    }
}

//! This function is documented in the header file
void gpu_init_atomdata(gmx_nbnxn_ocl_t* nb, const nbnxn_atomdata_t* nbat)
{
    cl_int           cl_error;
    int              nalloc, natoms;
    bool             realloced;
    bool             bDoTime = nb->bDoTime == CL_TRUE;
    cl_timers_t*     timers  = nb->timers;
    cl_atomdata_t*   d_atdat = nb->atdat;
    cl_command_queue ls      = nb->stream[InteractionLocality::Local];

    natoms    = nbat->numAtoms();
    realloced = false;

    if (bDoTime)
    {
        /* time async copy */
        timers->atdat.openTimingRegion(ls);
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initialized yet, i.e d_atdat->natoms == -1 */
    if (natoms > d_atdat->nalloc)
    {
        nalloc = over_alloc_small(natoms);

        /* free up first if the arrays have already been initialized */
        if (d_atdat->nalloc != -1)
        {
            freeDeviceBuffer(&d_atdat->f);
            freeDeviceBuffer(&d_atdat->xq);
            freeDeviceBuffer(&d_atdat->lj_comb);
            freeDeviceBuffer(&d_atdat->atom_types);
        }

        d_atdat->f_elem_size = sizeof(rvec);

        d_atdat->f = clCreateBuffer(nb->dev_rundata->context, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
                                    nalloc * d_atdat->f_elem_size, nullptr, &cl_error);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

        d_atdat->xq = clCreateBuffer(nb->dev_rundata->context, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                     nalloc * sizeof(cl_float4), nullptr, &cl_error);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());

        if (useLjCombRule(nb->nbparam->vdwtype))
        {
            d_atdat->lj_comb = clCreateBuffer(nb->dev_rundata->context,
                                              CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                              nalloc * sizeof(cl_float2), nullptr, &cl_error);
            GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                               ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
        }
        else
        {
            d_atdat->atom_types = clCreateBuffer(nb->dev_rundata->context,
                                                 CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                                 nalloc * sizeof(int), nullptr, &cl_error);
            GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                               ("clCreateBuffer failed: " + ocl_get_error_string(cl_error)).c_str());
        }

        d_atdat->nalloc = nalloc;
        realloced       = true;
    }

    d_atdat->natoms       = natoms;
    d_atdat->natoms_local = nbat->natoms_local;

    /* need to clear GPU f output if realloc happened */
    if (realloced)
    {
        nbnxn_ocl_clear_f(nb, nalloc);
    }

    if (useLjCombRule(nb->nbparam->vdwtype))
    {
        ocl_copy_H2D_async(d_atdat->lj_comb, nbat->params().lj_comb.data(), 0, natoms * sizeof(cl_float2),
                           ls, bDoTime ? timers->atdat.fetchNextEvent() : nullptr);
    }
    else
    {
        ocl_copy_H2D_async(d_atdat->atom_types, nbat->params().type.data(), 0, natoms * sizeof(int),
                           ls, bDoTime ? timers->atdat.fetchNextEvent() : nullptr);
    }

    if (bDoTime)
    {
        timers->atdat.closeTimingRegion(ls);
    }

    /* kick off the tasks enqueued above to ensure concurrency with the search */
    cl_error = clFlush(ls);
    GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                       ("clFlush failed: " + ocl_get_error_string(cl_error)).c_str());
}

/*! \brief Releases an OpenCL kernel pointer */
static void free_kernel(cl_kernel* kernel_ptr)
{
    cl_int gmx_unused cl_error;

    assert(nullptr != kernel_ptr);

    if (*kernel_ptr)
    {
        cl_error = clReleaseKernel(*kernel_ptr);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clReleaseKernel failed: " + ocl_get_error_string(cl_error)).c_str());

        *kernel_ptr = nullptr;
    }
}

/*! \brief Releases a list of OpenCL kernel pointers */
static void free_kernels(cl_kernel* kernels, int count)
{
    int i;

    for (i = 0; i < count; i++)
    {
        free_kernel(kernels + i);
    }
}

/*! \brief Free the OpenCL runtime data (context and program).
 *
 *  The function releases the OpenCL context and program assuciated with the
 *  device that the calling PP rank is running on.
 *
 *  \param runData [in]  porinter to the structure with runtime data.
 */
static void free_gpu_device_runtime_data(gmx_device_runtime_data_t* runData)
{
    if (runData == nullptr)
    {
        return;
    }

    cl_int gmx_unused cl_error;

    if (runData->context)
    {
        cl_error = clReleaseContext(runData->context);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clReleaseContext failed: " + ocl_get_error_string(cl_error)).c_str());
        runData->context = nullptr;
    }

    if (runData->program)
    {
        cl_error = clReleaseProgram(runData->program);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clReleaseProgram failed: " + ocl_get_error_string(cl_error)).c_str());
        runData->program = nullptr;
    }
}

//! This function is documented in the header file
void gpu_free(gmx_nbnxn_ocl_t* nb)
{
    if (nb == nullptr)
    {
        return;
    }

    /* Free kernels */
    int kernel_count = sizeof(nb->kernel_ener_noprune_ptr) / sizeof(nb->kernel_ener_noprune_ptr[0][0]);
    free_kernels(nb->kernel_ener_noprune_ptr[0], kernel_count);

    kernel_count = sizeof(nb->kernel_ener_prune_ptr) / sizeof(nb->kernel_ener_prune_ptr[0][0]);
    free_kernels(nb->kernel_ener_prune_ptr[0], kernel_count);

    kernel_count = sizeof(nb->kernel_noener_noprune_ptr) / sizeof(nb->kernel_noener_noprune_ptr[0][0]);
    free_kernels(nb->kernel_noener_noprune_ptr[0], kernel_count);

    kernel_count = sizeof(nb->kernel_noener_prune_ptr) / sizeof(nb->kernel_noener_prune_ptr[0][0]);
    free_kernels(nb->kernel_noener_prune_ptr[0], kernel_count);

    free_kernel(&(nb->kernel_zero_e_fshift));

    /* Free atdat */
    freeDeviceBuffer(&(nb->atdat->xq));
    freeDeviceBuffer(&(nb->atdat->f));
    freeDeviceBuffer(&(nb->atdat->e_lj));
    freeDeviceBuffer(&(nb->atdat->e_el));
    freeDeviceBuffer(&(nb->atdat->fshift));
    freeDeviceBuffer(&(nb->atdat->lj_comb));
    freeDeviceBuffer(&(nb->atdat->atom_types));
    freeDeviceBuffer(&(nb->atdat->shift_vec));
    sfree(nb->atdat);

    /* Free nbparam */
    freeDeviceBuffer(&(nb->nbparam->nbfp_climg2d));
    freeDeviceBuffer(&(nb->nbparam->nbfp_comb_climg2d));
    freeDeviceBuffer(&(nb->nbparam->coulomb_tab_climg2d));
    sfree(nb->nbparam);

    /* Free plist */
    auto* plist = nb->plist[InteractionLocality::Local];
    freeDeviceBuffer(&plist->sci);
    freeDeviceBuffer(&plist->cj4);
    freeDeviceBuffer(&plist->imask);
    freeDeviceBuffer(&plist->excl);
    sfree(plist);
    if (nb->bUseTwoStreams)
    {
        auto* plist_nl = nb->plist[InteractionLocality::NonLocal];
        freeDeviceBuffer(&plist_nl->sci);
        freeDeviceBuffer(&plist_nl->cj4);
        freeDeviceBuffer(&plist_nl->imask);
        freeDeviceBuffer(&plist_nl->excl);
        sfree(plist_nl);
    }

    /* Free nbst */
    pfree(nb->nbst.e_lj);
    nb->nbst.e_lj = nullptr;

    pfree(nb->nbst.e_el);
    nb->nbst.e_el = nullptr;

    pfree(nb->nbst.fshift);
    nb->nbst.fshift = nullptr;

    /* Free command queues */
    clReleaseCommandQueue(nb->stream[InteractionLocality::Local]);
    nb->stream[InteractionLocality::Local] = nullptr;
    if (nb->bUseTwoStreams)
    {
        clReleaseCommandQueue(nb->stream[InteractionLocality::NonLocal]);
        nb->stream[InteractionLocality::NonLocal] = nullptr;
    }
    /* Free other events */
    if (nb->nonlocal_done)
    {
        clReleaseEvent(nb->nonlocal_done);
        nb->nonlocal_done = nullptr;
    }
    if (nb->misc_ops_and_local_H2D_done)
    {
        clReleaseEvent(nb->misc_ops_and_local_H2D_done);
        nb->misc_ops_and_local_H2D_done = nullptr;
    }

    free_gpu_device_runtime_data(nb->dev_rundata);
    sfree(nb->dev_rundata);

    /* Free timers and timings */
    delete nb->timers;
    sfree(nb->timings);
    sfree(nb);

    if (debug)
    {
        fprintf(debug, "Cleaned up OpenCL data structures.\n");
    }
}

//! This function is documented in the header file
gmx_wallclock_gpu_nbnxn_t* gpu_get_timings(gmx_nbnxn_ocl_t* nb)
{
    return (nb != nullptr && nb->bDoTime) ? nb->timings : nullptr;
}

//! This function is documented in the header file
void gpu_reset_timings(nonbonded_verlet_t* nbv)
{
    if (nbv->gpu_nbv && nbv->gpu_nbv->bDoTime)
    {
        init_timings(nbv->gpu_nbv->timings);
    }
}

//! This function is documented in the header file
int gpu_min_ci_balanced(gmx_nbnxn_ocl_t* nb)
{
    return nb != nullptr ? gpu_min_ci_balanced_factor * nb->dev_info->compute_units : 0;
}

//! This function is documented in the header file
gmx_bool gpu_is_kernel_ewald_analytical(const gmx_nbnxn_ocl_t* nb)
{
    return ((nb->nbparam->eeltype == eelOclEWALD_ANA) || (nb->nbparam->eeltype == eelOclEWALD_ANA_TWIN));
}

} // namespace Nbnxm
