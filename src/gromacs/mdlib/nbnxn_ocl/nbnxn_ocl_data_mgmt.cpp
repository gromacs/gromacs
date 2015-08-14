/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 *  \brief Define OpenCL implementation of nbnxn_gpu_data_mgmt.h
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 */
#include "gmxpre.h"

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/gmxlib/gpu_utils/gpu_utils.h"
#include "gromacs/gmxlib/ocl_tools/oclutils.h"
#include "gromacs/legacyheaders/gmx_detect_hardware.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/force_flags.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_gpu_jit_support.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "nbnxn_ocl_types.h"


/*! \brief This parameter should be determined heuristically from the
 * kernel execution times
 *
 * This value is best for small systems on a single AMD Radeon R9 290X
 * (and about 5% faster than 40, which is the default for CUDA
 * devices). Larger simulation systems were quite insensitive to the
 * value of this parameter.
 */
static unsigned int gpu_min_ci_balanced_factor = 50;

/*! \brief Helper function for warning output
 *
 * We should actually be using md_print_warn in md_logging.c,
 * but we can't include mpi.h in OpenCL code.
 */
static void md_print_warn(FILE       *fplog,
                          const char *fmt, ...)
{
    va_list ap;

    if (fplog != NULL)
    {
        /* We should only print to stderr on the master node,
         * in most cases fplog is only set on the master node, so this works.
         */
        va_start(ap, fmt);
        fprintf(stderr, "\n");
        vfprintf(stderr, fmt, ap);
        fprintf(stderr, "\n");
        va_end(ap);

        va_start(ap, fmt);
        fprintf(fplog, "\n");
        vfprintf(fplog, fmt, ap);
        fprintf(fplog, "\n");
        va_end(ap);
    }
}

/*! \brief Free device buffers
 *
 * If the pointers to the size variables are NULL no resetting happens.
 */
void ocl_free_buffered(cl_mem d_ptr, int *n, int *nalloc)
{
    cl_int gmx_unused cl_error;

    if (d_ptr)
    {
        cl_error = clReleaseMemObject(d_ptr);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors
    }

    if (n)
    {
        *n = -1;
    }

    if (nalloc)
    {
        *nalloc = -1;
    }
}

/*! \brief Reallocation device buffers
 *
 *  Reallocation of the memory pointed by d_ptr and copying of the data from
 *  the location pointed by h_src host-side pointer is done. Allocation is
 *  buffered and therefore freeing is only needed if the previously allocated
 *  space is not enough.
 *  The H2D copy is launched in command queue s and can be done synchronously or
 *  asynchronously (the default is the latter).
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying the H2D copy. The event can further be used to queue a wait
 *  for this operation or to query profiling information.
 *  OpenCL equivalent of cu_realloc_buffered.
 */
void ocl_realloc_buffered(cl_mem *d_dest, void *h_src,
                          size_t type_size,
                          int *curr_size, int *curr_alloc_size,
                          int req_size,
                          cl_context context,
                          cl_command_queue s,
                          bool bAsync = true,
                          cl_event *copy_event = NULL)
{
    cl_int cl_error;

    if (d_dest == NULL || req_size < 0)
    {
        return;
    }

    /* reallocate only if the data does not fit = allocation size is smaller
       than the current requested size */
    if (req_size > *curr_alloc_size)
    {
        /* only free if the array has already been initialized */
        if (*curr_alloc_size >= 0)
        {
            ocl_free_buffered(*d_dest, curr_size, curr_alloc_size);
        }

        *curr_alloc_size = over_alloc_large(req_size);

        *d_dest = clCreateBuffer(context, CL_MEM_READ_WRITE, *curr_alloc_size * type_size, NULL, &cl_error);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors, check clCreateBuffer flags
    }

    /* size could have changed without actual reallocation */
    *curr_size = req_size;

    /* upload to device */
    if (h_src)
    {
        if (bAsync)
        {
            ocl_copy_H2D_async(*d_dest, h_src, 0, *curr_size * type_size, s, copy_event);
        }
        else
        {
            ocl_copy_H2D(*d_dest, h_src,  0, *curr_size * type_size, s);
        }
    }
}

/*! \brief Releases the input OpenCL buffer */
static void free_ocl_buffer(cl_mem *buffer)
{
    cl_int gmx_unused cl_error;

    assert(NULL != buffer);

    if (*buffer)
    {
        cl_error = clReleaseMemObject(*buffer);
        assert(CL_SUCCESS == cl_error);
        *buffer = NULL;
    }
}

/*! \brief Tabulates the Ewald Coulomb force and initializes the size/scale
 * and the table GPU array.
 *
 * If called with an already allocated table, it just re-uploads the
 * table.
 */
static void init_ewald_coulomb_force_table(const interaction_const_t *ic,
                                           cl_nbparam_t              *nbp,
                                           const gmx_device_info_t   *dev_info)
{
    cl_mem       coul_tab;

    cl_int       cl_error;

    if (nbp->coulomb_tab_climg2d != NULL)
    {
        free_ocl_buffer(&(nbp->coulomb_tab_climg2d));
    }

    /* Switched from using textures to using buffers */
    // TODO: decide which alternative is most efficient - textures or buffers.
    /*
       cl_image_format array_format;

       array_format.image_channel_data_type = CL_FLOAT;
       array_format.image_channel_order     = CL_R;

       coul_tab = clCreateImage2D(dev_info->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
       &array_format, tabsize, 1, 0, ftmp, &cl_error);
     */

    coul_tab = clCreateBuffer(dev_info->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, ic->tabq_size*sizeof(cl_float), ic->tabq_coul_F, &cl_error);
    assert(cl_error == CL_SUCCESS);
    // TODO: handle errors, check clCreateBuffer flags

    nbp->coulomb_tab_climg2d  = coul_tab;
    nbp->coulomb_tab_size     = ic->tabq_size;
    nbp->coulomb_tab_scale    = ic->tabq_scale;
}


/*! \brief Initializes the atomdata structure first time, it only gets filled at
    pair-search.
 */
static void init_atomdata_first(cl_atomdata_t *ad, int ntypes, gmx_device_info_t *dev_info)
{
    cl_int cl_error;

    ad->ntypes  = ntypes;

    /* An element of the shift_vec device buffer has the same size as one element
       of the host side shift_vec buffer. */
    ad->shift_vec_elem_size = sizeof(*(((nbnxn_atomdata_t*)0)->shift_vec));

    // TODO: handle errors, check clCreateBuffer flags
    ad->shift_vec = clCreateBuffer(dev_info->context, CL_MEM_READ_WRITE, SHIFTS * ad->shift_vec_elem_size, NULL, &cl_error);
    assert(cl_error == CL_SUCCESS);
    ad->bShiftVecUploaded = false;

    /* An element of the fshift device buffer has the same size as one element
       of the host side fshift buffer. */
    ad->fshift_elem_size = sizeof(*(((cl_nb_staging_t*)0)->fshift));

    ad->fshift = clCreateBuffer(dev_info->context, CL_MEM_READ_WRITE, SHIFTS * ad->fshift_elem_size, NULL, &cl_error);
    assert(cl_error == CL_SUCCESS);
    // TODO: handle errors, check clCreateBuffer flags

    ad->e_lj = clCreateBuffer(dev_info->context, CL_MEM_READ_WRITE, sizeof(float), NULL, &cl_error);
    assert(cl_error == CL_SUCCESS);
    // TODO: handle errors, check clCreateBuffer flags

    ad->e_el = clCreateBuffer(dev_info->context, CL_MEM_READ_WRITE, sizeof(float), NULL, &cl_error);
    assert(cl_error == CL_SUCCESS);
    // TODO: handle errors, check clCreateBuffer flags

    /* initialize to NULL pointers to data that is not allocated here and will
       need reallocation in nbnxn_gpu_init_atomdata */
    ad->xq = NULL;
    ad->f  = NULL;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    ad->natoms = -1;
    ad->nalloc = -1;
}

/*! \brief Copies all parameters related to the cut-off from ic to nbp
 */
static void set_cutoff_parameters(cl_nbparam_t              *nbp,
                                  const interaction_const_t *ic)
{
    nbp->ewald_beta       = ic->ewaldcoeff_q;
    nbp->sh_ewald         = ic->sh_ewald;
    nbp->epsfac           = ic->epsfac;
    nbp->two_k_rf         = 2.0 * ic->k_rf;
    nbp->c_rf             = ic->c_rf;
    nbp->rvdw_sq          = ic->rvdw * ic->rvdw;
    nbp->rcoulomb_sq      = ic->rcoulomb * ic->rcoulomb;
    nbp->rlist_sq         = ic->rlist * ic->rlist;

    nbp->sh_lj_ewald      = ic->sh_lj_ewald;
    nbp->ewaldcoeff_lj    = ic->ewaldcoeff_lj;

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
static void
map_interaction_types_to_gpu_kernel_flavors(const interaction_const_t *ic,
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

/*! \brief Initializes the nonbonded parameter data structure.
 */
static void init_nbparam(cl_nbparam_t              *nbp,
                         const interaction_const_t *ic,
                         const nbnxn_atomdata_t    *nbat,
                         const gmx_device_info_t   *dev_info)
{
    int         ntypes, nnbfp, nnbfp_comb;
    cl_int      cl_error;


    ntypes = nbat->ntype;

    set_cutoff_parameters(nbp, ic);

    map_interaction_types_to_gpu_kernel_flavors(ic,
                                                &(nbp->eeltype),
                                                &(nbp->vdwtype));

    if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == ljcrGEOM)
        {
            assert(nbat->comb_rule == ljcrGEOM);
        }
        else
        {
            assert(nbat->comb_rule == ljcrLB);
        }
    }
    /* generate table for PME */
    nbp->coulomb_tab_climg2d = NULL;
    if (nbp->eeltype == eelOclEWALD_TAB || nbp->eeltype == eelOclEWALD_TAB_TWIN)
    {
        init_ewald_coulomb_force_table(ic, nbp, dev_info);
    }
    else
    // TODO: improvement needed.
    // The image2d is created here even if eeltype is not eelCuEWALD_TAB or eelCuEWALD_TAB_TWIN because the OpenCL kernels
    // don't accept NULL values for image2D parameters.
    {
        /* Switched from using textures to using buffers */
        // TODO: decide which alternative is most efficient - textures or buffers.
        /*
           cl_image_format array_format;

           array_format.image_channel_data_type = CL_FLOAT;
           array_format.image_channel_order     = CL_R;

           nbp->coulomb_tab_climg2d = clCreateImage2D(dev_info->context, CL_MEM_READ_WRITE,
            &array_format, 1, 1, 0, NULL, &cl_error);
         */

        nbp->coulomb_tab_climg2d = clCreateBuffer(dev_info->context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &cl_error);
        // TODO: handle errors
    }

    nnbfp      = 2*ntypes*ntypes;
    nnbfp_comb = 2*ntypes;

    {
        /* Switched from using textures to using buffers */
        // TODO: decide which alternative is most efficient - textures or buffers.
        /*
           cl_image_format array_format;

           array_format.image_channel_data_type = CL_FLOAT;
           array_format.image_channel_order     = CL_R;

           nbp->nbfp_climg2d = clCreateImage2D(dev_info->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            &array_format, nnbfp, 1, 0, nbat->nbfp, &cl_error);
         */

        nbp->nbfp_climg2d = clCreateBuffer(dev_info->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, nnbfp*sizeof(cl_float), nbat->nbfp, &cl_error);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors

        if (ic->vdwtype == evdwPME)
        {
            /* Switched from using textures to using buffers */
            // TODO: decide which alternative is most efficient - textures or buffers.
            /*  nbp->nbfp_comb_climg2d = clCreateImage2D(dev_info->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                &array_format, nnbfp_comb, 1, 0, nbat->nbfp_comb, &cl_error);*/
            nbp->nbfp_comb_climg2d = clCreateBuffer(dev_info->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, nnbfp_comb*sizeof(cl_float), nbat->nbfp_comb, &cl_error);


            assert(cl_error == CL_SUCCESS);
            // TODO: handle errors
        }
        else
        {
            // TODO: improvement needed.
            // The image2d is created here even if vdwtype is not evdwPME because the OpenCL kernels
            // don't accept NULL values for image2D parameters.
            /* Switched from using textures to using buffers */
            // TODO: decide which alternative is most efficient - textures or buffers.
            /* nbp->nbfp_comb_climg2d = clCreateImage2D(dev_info->context, CL_MEM_READ_WRITE,
                &array_format, 1, 1, 0, NULL, &cl_error);*/
            nbp->nbfp_comb_climg2d = clCreateBuffer(dev_info->context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &cl_error);


            assert(cl_error == CL_SUCCESS);
            // TODO: handle errors
        }
    }
}

//! This function is documented in the header file
void nbnxn_gpu_pme_loadbal_update_param(const nonbonded_verlet_t    *nbv,
                                        const interaction_const_t   *ic)
{
    if (!nbv || nbv->grp[0].kernel_type != nbnxnk8x8x8_GPU)
    {
        return;
    }
    gmx_nbnxn_ocl_t    *nb  = nbv->gpu_nbv;
    cl_nbparam_t       *nbp = nb->nbparam;

    set_cutoff_parameters(nbp, ic);

    nbp->eeltype = nbnxn_gpu_pick_ewald_kernel_type(ic->rcoulomb != ic->rvdw);

    init_ewald_coulomb_force_table(ic, nb->nbparam, nb->dev_info);
}

/*! \brief Initializes the pair list data structure.
 */
static void init_plist(cl_plist_t *pl)
{
    /* initialize to NULL pointers to data that is not allocated here and will
       need reallocation in nbnxn_gpu_init_pairlist */
    pl->sci     = NULL;
    pl->cj4     = NULL;
    pl->excl    = NULL;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    pl->na_c        = -1;
    pl->nsci        = -1;
    pl->sci_nalloc  = -1;
    pl->ncj4        = -1;
    pl->cj4_nalloc  = -1;
    pl->nexcl       = -1;
    pl->excl_nalloc = -1;
    pl->bDoPrune    = false;
}

/*! \brief Initializes the timer data structure.
 */
static void init_timers(cl_timers_t gmx_unused *t, bool gmx_unused bUseTwoStreams)
{
    /* Nothing to initialize for OpenCL */
}

/*! \brief Initializes the timings data structure.
 */
static void init_timings(gmx_wallclock_gpu_t *t)
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
}

/*! \brief Creates context for OpenCL GPU given by \p mygpu
 *
 * A fatal error results if creation fails.
 *
 * \param[inout] nb        Manages OpenCL non-bonded calculations;
 *                         contexts returned in dev_info members
 * \param[in]    rank      MPI rank (for error reporting)
 */
static void
nbnxn_gpu_create_context(gmx_nbnxn_ocl_t           *nb,
                         int                        rank)
{
    cl_context_properties     context_properties[3];
    cl_platform_id            platform_id;
    cl_device_id              device_id;
    cl_context                context;
    cl_int                    cl_error;

    platform_id      = nb->dev_info->ocl_gpu_id.ocl_platform_id;
    device_id        = nb->dev_info->ocl_gpu_id.ocl_device_id;

    context_properties[0] = CL_CONTEXT_PLATFORM;
    context_properties[1] = (cl_context_properties) platform_id;
    context_properties[2] = 0; /* Terminates the list of properties */

    context = clCreateContext(context_properties, 1, &device_id, NULL, NULL, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS, "On rank %d failed to create context for GPU #%s: OpenCL error %d",
                  rank,
                  nb->dev_info->device_name,
                  cl_error);
        return;
    }

    nb->dev_info->context = context;
}

/*! \brief Initializes the OpenCL kernel pointers of the nbnxn_ocl_ptr_t input data structure. */
static cl_kernel nbnxn_gpu_create_kernel(gmx_nbnxn_ocl_t *nb,
                                         const char      *kernel_name)
{
    cl_kernel kernel;
    cl_int    cl_error;

    kernel = clCreateKernel(nb->dev_info->program, kernel_name, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS, "Failed to create kernel '%s' for GPU #%s: OpenCL error %d",
                  kernel_name,
                  nb->dev_info->device_name,
                  cl_error);
    }

    return kernel;
}

/*! \brief Clears nonbonded shift force output array and energy outputs on the GPU.
 */
static void
nbnxn_ocl_clear_e_fshift(gmx_nbnxn_ocl_t *nb)
{

    cl_int               cl_error;
    cl_atomdata_t *      adat     = nb->atdat;
    cl_command_queue     ls       = nb->stream[eintLocal];

    size_t               local_work_size[3]   = {1, 1, 1};
    size_t               global_work_size[3]  = {1, 1, 1};

    cl_int               shifts   = SHIFTS*3;

    cl_int               arg_no;

    cl_kernel            zero_e_fshift = nb->kernel_zero_e_fshift;

    local_work_size[0]   = 64;
    global_work_size[0]  = ((shifts/64)*64) + ((shifts%64) ? 64 : 0);

    arg_no    = 0;
    cl_error  = clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_mem), &(adat->fshift));
    cl_error |= clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_mem), &(adat->e_lj));
    cl_error |= clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_mem), &(adat->e_el));
    cl_error |= clSetKernelArg(zero_e_fshift, arg_no++, sizeof(cl_uint), &shifts);
    assert(cl_error == CL_SUCCESS);

    cl_error = clEnqueueNDRangeKernel(ls, zero_e_fshift, 3, NULL, global_work_size, local_work_size, 0, NULL, NULL);
    assert(cl_error == CL_SUCCESS);

}

/*! \brief Initializes the OpenCL kernel pointers of the nbnxn_ocl_ptr_t input data structure. */
static void nbnxn_gpu_init_kernels(gmx_nbnxn_ocl_t *nb)
{
    /* Init to 0 main kernel arrays */
    /* They will be later on initialized in select_nbnxn_kernel */
    memset(nb->kernel_ener_noprune_ptr, 0, sizeof(nb->kernel_ener_noprune_ptr));
    memset(nb->kernel_ener_prune_ptr, 0, sizeof(nb->kernel_ener_prune_ptr));
    memset(nb->kernel_noener_noprune_ptr, 0, sizeof(nb->kernel_noener_noprune_ptr));
    memset(nb->kernel_noener_prune_ptr, 0, sizeof(nb->kernel_noener_prune_ptr));

    /* Init auxiliary kernels */
    nb->kernel_memset_f      = nbnxn_gpu_create_kernel(nb, "memset_f");
    nb->kernel_memset_f2     = nbnxn_gpu_create_kernel(nb, "memset_f2");
    nb->kernel_memset_f3     = nbnxn_gpu_create_kernel(nb, "memset_f3");
    nb->kernel_zero_e_fshift = nbnxn_gpu_create_kernel(nb, "zero_e_fshift");
}

//! This function is documented in the header file
void nbnxn_gpu_init(FILE gmx_unused           *fplog,
                    gmx_nbnxn_ocl_t          **p_nb,
                    const gmx_gpu_info_t      *gpu_info,
                    const gmx_gpu_opt_t       *gpu_opt,
                    const interaction_const_t *ic,
                    nonbonded_verlet_group_t  *nbv_grp,
                    int                        my_gpu_index,
                    int                        rank,
                    gmx_bool                   bLocalAndNonlocal)
{
    gmx_nbnxn_ocl_t            *nb;
    cl_int                      cl_error;
    /*
       bool gmx_unused             bStreamSync;
       bool gmx_unused             bNoStreamSync;
       bool gmx_unused             bTMPIAtomics;
       bool gmx_unused             bX86;
       bool gmx_unused             bOldDriver;
     */
    cl_command_queue_properties queue_properties;

    assert(gpu_info);
    assert(gpu_opt);
    assert(ic);

    if (p_nb == NULL)
    {
        return;
    }

    snew(nb, 1);
    snew(nb->atdat, 1);
    snew(nb->nbparam, 1);
    snew(nb->plist[eintLocal], 1);
    if (bLocalAndNonlocal)
    {
        snew(nb->plist[eintNonlocal], 1);
    }

    nb->bUseTwoStreams = bLocalAndNonlocal;

    snew(nb->timers, 1);
    snew(nb->timings, 1);

    /* set device info, just point it to the right GPU among the detected ones */
    nb->dev_info = gpu_info->gpu_dev + gpu_opt->dev_use[my_gpu_index];

    /* init to NULL the debug buffer */
    nb->debug_buffer = NULL;

    /* init nbst */
    ocl_pmalloc((void**)&nb->nbst.e_lj, sizeof(*nb->nbst.e_lj));
    ocl_pmalloc((void**)&nb->nbst.e_el, sizeof(*nb->nbst.e_el));
    ocl_pmalloc((void**)&nb->nbst.fshift, SHIFTS * sizeof(*nb->nbst.fshift));

    init_plist(nb->plist[eintLocal]);

    /* OpenCL timing disabled if GMX_DISABLE_OCL_TIMING is defined. */
    nb->bDoTime = (getenv("GMX_DISABLE_OCL_TIMING") == NULL);

    /* Create queues only after bDoTime has been initialized */
    if (nb->bDoTime)
    {
        queue_properties = CL_QUEUE_PROFILING_ENABLE;
    }
    else
    {
        queue_properties = 0;
    }

    nbnxn_gpu_create_context(nb, rank);

    /* local/non-local GPU streams */
    nb->stream[eintLocal] = clCreateCommandQueue(nb->dev_info->context, nb->dev_info->ocl_gpu_id.ocl_device_id, queue_properties, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS, "On rank %d failed to create context for GPU #%s: OpenCL error %d",
                  rank,
                  nb->dev_info->device_name,
                  cl_error);
        return;
    }

    if (nb->bUseTwoStreams)
    {
        init_plist(nb->plist[eintNonlocal]);

        nb->stream[eintNonlocal] = clCreateCommandQueue(nb->dev_info->context, nb->dev_info->ocl_gpu_id.ocl_device_id, queue_properties, &cl_error);
        if (CL_SUCCESS != cl_error)
        {
            gmx_fatal(FARGS, "On rank %d failed to create context for GPU #%s: OpenCL error %d",
                      rank,
                      nb->dev_info->device_name,
                      cl_error);
            return;
        }
    }

    if (nb->bDoTime)
    {
        init_timers(nb->timers, nb->bUseTwoStreams);
        init_timings(nb->timings);
    }

    // TODO: check if it's worth implementing for NVIDIA GPUs
    ///////////* set the kernel type for the current GPU */
    ///////////* pick L1 cache configuration */
    //////////nbnxn_gpu_set_cacheconfig(nb->dev_info);

    init_atomdata_first(nb->atdat, nbv_grp[0].nbat->ntype, nb->dev_info);
    init_nbparam(nb->nbparam, ic, nbv_grp[0].nbat, nb->dev_info);
    nbnxn_gpu_compile_kernels(nb);
    nbnxn_gpu_init_kernels(nb);
    // TODO put this elsewhere? also mirror it in cuda
    nbnxn_ocl_clear_e_fshift(nb);

    *p_nb = nb;

    if (debug)
    {
        fprintf(debug, "Initialized OpenCL data structures.\n");
    }
}

/*! \brief Clears the first natoms_clear elements of the GPU nonbonded force output array.
 */
static void nbnxn_ocl_clear_f(gmx_nbnxn_ocl_t *nb, int natoms_clear)
{

    cl_int               cl_error;
    cl_atomdata_t *      adat     = nb->atdat;
    cl_command_queue     ls       = nb->stream[eintLocal];
    cl_float             value    = 0.0f;

    size_t               local_work_size[3]  = {1, 1, 1};
    size_t               global_work_size[3] = {1, 1, 1};

    cl_int               arg_no;

    cl_kernel            memset_f = nb->kernel_memset_f;

    cl_uint              natoms_flat = natoms_clear * (sizeof(rvec)/sizeof(real));

    local_work_size[0]  = 64;
    global_work_size[0] = ((natoms_flat/local_work_size[0])*local_work_size[0]) + ((natoms_flat%local_work_size[0]) ? local_work_size[0] : 0);

    arg_no    = 0;
    cl_error  = clSetKernelArg(memset_f, arg_no++, sizeof(cl_mem), &(adat->f));
    cl_error |= clSetKernelArg(memset_f, arg_no++, sizeof(cl_float), &value);
    cl_error |= clSetKernelArg(memset_f, arg_no++, sizeof(cl_uint), &natoms_flat);
    assert(cl_error == CL_SUCCESS);

    cl_error = clEnqueueNDRangeKernel(ls, memset_f, 3, NULL, global_work_size, local_work_size, 0, NULL, NULL);
    assert(cl_error == CL_SUCCESS);
}

//! This function is documented in the header file
void
nbnxn_gpu_clear_outputs(gmx_nbnxn_ocl_t   *nb,
                        int                flags)
{
    nbnxn_ocl_clear_f(nb, nb->atdat->natoms);
    /* clear shift force array and energies if the outputs were
       used in the current step */
    if (flags & GMX_FORCE_VIRIAL)
    {
        nbnxn_ocl_clear_e_fshift(nb);
    }

    /* kick off buffer clearing kernel to ensure concurrency with constraints/update */
    cl_int gmx_unused cl_error;
    cl_error = clFlush(nb->stream[eintLocal]);
    assert(CL_SUCCESS == cl_error);
}

//! This function is documented in the header file
void nbnxn_gpu_init_pairlist(gmx_nbnxn_ocl_t        *nb,
                             const nbnxn_pairlist_t *h_plist,
                             int                     iloc)
{
    char             sbuf[STRLEN];
    cl_command_queue stream     = nb->stream[iloc];
    cl_plist_t      *d_plist    = nb->plist[iloc];

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

    ocl_realloc_buffered(&d_plist->sci, h_plist->sci, sizeof(nbnxn_sci_t),
                         &d_plist->nsci, &d_plist->sci_nalloc,
                         h_plist->nsci,
                         nb->dev_info->context,
                         stream, true, &(nb->timers->pl_h2d_sci[iloc]));

    ocl_realloc_buffered(&d_plist->cj4, h_plist->cj4, sizeof(nbnxn_cj4_t),
                         &d_plist->ncj4, &d_plist->cj4_nalloc,
                         h_plist->ncj4,
                         nb->dev_info->context,
                         stream, true, &(nb->timers->pl_h2d_cj4[iloc]));

    ocl_realloc_buffered(&d_plist->excl, h_plist->excl, sizeof(nbnxn_excl_t),
                         &d_plist->nexcl, &d_plist->excl_nalloc,
                         h_plist->nexcl,
                         nb->dev_info->context,
                         stream, true, &(nb->timers->pl_h2d_excl[iloc]));

    /* need to prune the pair list during the next step */
    d_plist->bDoPrune = true;
}

//! This function is documented in the header file
void nbnxn_gpu_upload_shiftvec(gmx_nbnxn_ocl_t        *nb,
                               const nbnxn_atomdata_t *nbatom)
{
    cl_atomdata_t   *adat  = nb->atdat;
    cl_command_queue ls    = nb->stream[eintLocal];

    /* only if we have a dynamic box */
    if (nbatom->bDynamicBox || !adat->bShiftVecUploaded)
    {
        ocl_copy_H2D_async(adat->shift_vec, nbatom->shift_vec, 0,
                           SHIFTS * adat->shift_vec_elem_size, ls, NULL);
        adat->bShiftVecUploaded = true;
    }
}

//! This function is documented in the header file
void nbnxn_gpu_init_atomdata(gmx_nbnxn_ocl_t               *nb,
                             const struct nbnxn_atomdata_t *nbat)
{
    cl_int           cl_error;
    int              nalloc, natoms;
    bool             realloced;
    bool             bDoTime = nb->bDoTime;
    cl_timers_t     *timers  = nb->timers;
    cl_atomdata_t   *d_atdat = nb->atdat;
    cl_command_queue ls      = nb->stream[eintLocal];

    natoms    = nbat->natoms;
    realloced = false;

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initialized yet, i.e d_atdat->natoms == -1 */
    if (natoms > d_atdat->nalloc)
    {
        nalloc = over_alloc_small(natoms);

        /* free up first if the arrays have already been initialized */
        if (d_atdat->nalloc != -1)
        {
            ocl_free_buffered(d_atdat->f, &d_atdat->natoms, &d_atdat->nalloc);
            ocl_free_buffered(d_atdat->xq, NULL, NULL);
            ocl_free_buffered(d_atdat->atom_types, NULL, NULL);
        }

        d_atdat->f_elem_size = sizeof(rvec);

        // TODO: handle errors, check clCreateBuffer flags
        d_atdat->f = clCreateBuffer(nb->dev_info->context, CL_MEM_READ_WRITE, nalloc * d_atdat->f_elem_size, NULL, &cl_error);
        assert(CL_SUCCESS == cl_error);

        d_atdat->xq = clCreateBuffer(nb->dev_info->context, CL_MEM_READ_WRITE, nalloc * sizeof(cl_float4), NULL, &cl_error);
        assert(CL_SUCCESS == cl_error);
        // TODO: handle errors, check clCreateBuffer flags

        d_atdat->atom_types = clCreateBuffer(nb->dev_info->context, CL_MEM_READ_WRITE, nalloc * sizeof(int), NULL, &cl_error);
        assert(CL_SUCCESS == cl_error);
        // TODO: handle errors, check clCreateBuffer flags

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

    ocl_copy_H2D_async(d_atdat->atom_types, nbat->type, 0,
                       natoms*sizeof(int), ls, bDoTime ? &(timers->atdat) : NULL);

    /* kick off the tasks enqueued above to ensure concurrency with the search */
    cl_error = clFlush(ls);
    assert(CL_SUCCESS == cl_error);
}

/*! \brief Releases an OpenCL kernel pointer */
void free_kernel(cl_kernel *kernel_ptr)
{
    cl_int gmx_unused cl_error;

    assert(NULL != kernel_ptr);

    if (*kernel_ptr)
    {
        cl_error = clReleaseKernel(*kernel_ptr);
        assert(cl_error == CL_SUCCESS);

        *kernel_ptr = NULL;
    }
}

/*! \brief Releases a list of OpenCL kernel pointers */
void free_kernels(cl_kernel *kernels, int count)
{
    int i;

    for (i = 0; i < count; i++)
    {
        free_kernel(kernels + i);
    }
}

//! This function is documented in the header file
void nbnxn_gpu_free(gmx_nbnxn_ocl_t *nb)
{
    int    kernel_count;

    /* Free kernels */
    kernel_count = sizeof(nb->kernel_ener_noprune_ptr) / sizeof(nb->kernel_ener_noprune_ptr[0][0]);
    free_kernels((cl_kernel*)nb->kernel_ener_noprune_ptr, kernel_count);

    kernel_count = sizeof(nb->kernel_ener_prune_ptr) / sizeof(nb->kernel_ener_prune_ptr[0][0]);
    free_kernels((cl_kernel*)nb->kernel_ener_prune_ptr, kernel_count);

    kernel_count = sizeof(nb->kernel_noener_noprune_ptr) / sizeof(nb->kernel_noener_noprune_ptr[0][0]);
    free_kernels((cl_kernel*)nb->kernel_noener_noprune_ptr, kernel_count);

    kernel_count = sizeof(nb->kernel_noener_prune_ptr) / sizeof(nb->kernel_noener_prune_ptr[0][0]);
    free_kernels((cl_kernel*)nb->kernel_noener_prune_ptr, kernel_count);

    free_kernel(&(nb->kernel_memset_f));
    free_kernel(&(nb->kernel_memset_f2));
    free_kernel(&(nb->kernel_memset_f3));
    free_kernel(&(nb->kernel_zero_e_fshift));

    /* Free atdat */
    free_ocl_buffer(&(nb->atdat->xq));
    free_ocl_buffer(&(nb->atdat->f));
    free_ocl_buffer(&(nb->atdat->e_lj));
    free_ocl_buffer(&(nb->atdat->e_el));
    free_ocl_buffer(&(nb->atdat->fshift));
    free_ocl_buffer(&(nb->atdat->atom_types));
    free_ocl_buffer(&(nb->atdat->shift_vec));
    sfree(nb->atdat);

    /* Free nbparam */
    free_ocl_buffer(&(nb->nbparam->nbfp_climg2d));
    free_ocl_buffer(&(nb->nbparam->nbfp_comb_climg2d));
    free_ocl_buffer(&(nb->nbparam->coulomb_tab_climg2d));
    sfree(nb->nbparam);

    /* Free plist */
    free_ocl_buffer(&(nb->plist[eintLocal]->sci));
    free_ocl_buffer(&(nb->plist[eintLocal]->cj4));
    free_ocl_buffer(&(nb->plist[eintLocal]->excl));
    sfree(nb->plist[eintLocal]);
    if (nb->bUseTwoStreams)
    {
        free_ocl_buffer(&(nb->plist[eintNonlocal]->sci));
        free_ocl_buffer(&(nb->plist[eintNonlocal]->cj4));
        free_ocl_buffer(&(nb->plist[eintNonlocal]->excl));
        sfree(nb->plist[eintNonlocal]);
    }

    /* Free nbst */
    ocl_pfree(nb->nbst.e_lj);
    nb->nbst.e_lj = NULL;

    ocl_pfree(nb->nbst.e_el);
    nb->nbst.e_el = NULL;

    ocl_pfree(nb->nbst.fshift);
    nb->nbst.fshift = NULL;

    /* Free debug buffer */
    free_ocl_buffer(&nb->debug_buffer);

    /* Free command queues */
    clReleaseCommandQueue(nb->stream[eintLocal]);
    nb->stream[eintLocal] = NULL;
    if (nb->bUseTwoStreams)
    {
        clReleaseCommandQueue(nb->stream[eintNonlocal]);
        nb->stream[eintNonlocal] = NULL;
    }
    /* Free other events */
    if (nb->nonlocal_done)
    {
        clReleaseEvent(nb->nonlocal_done);
        nb->nonlocal_done = NULL;
    }
    if (nb->misc_ops_and_local_H2D_done)
    {
        clReleaseEvent(nb->misc_ops_and_local_H2D_done);
        nb->misc_ops_and_local_H2D_done = NULL;
    }

    /* Free timers and timings */
    sfree(nb->timers);
    sfree(nb->timings);
    sfree(nb);

    if (debug)
    {
        fprintf(debug, "Cleaned up OpenCL data structures.\n");
    }
}

//! This function is documented in the header file
gmx_wallclock_gpu_t * nbnxn_gpu_get_timings(gmx_nbnxn_ocl_t *nb)
{
    return (nb != NULL && nb->bDoTime) ? nb->timings : NULL;
}

//! This function is documented in the header file
void nbnxn_gpu_reset_timings(nonbonded_verlet_t* nbv)
{
    if (nbv->gpu_nbv && nbv->gpu_nbv->bDoTime)
    {
        init_timings(nbv->gpu_nbv->timings);
    }
}

//! This function is documented in the header file
int nbnxn_gpu_min_ci_balanced(gmx_nbnxn_ocl_t *nb)
{
    return nb != NULL ?
           gpu_min_ci_balanced_factor * nb->dev_info->compute_units : 0;
}

//! This function is documented in the header file
gmx_bool nbnxn_gpu_is_kernel_ewald_analytical(const gmx_nbnxn_ocl_t *nb)
{
    return ((nb->nbparam->eeltype == eelOclEWALD_ANA) ||
            (nb->nbparam->eeltype == eelOclEWALD_ANA_TWIN));
}
