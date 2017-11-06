/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 *  \brief
 *  Data types used internally in the nbnxn_ocl module.
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Szilárd Páll <pszilard@kth.se>
 *  \ingroup module_mdlib
 */

#ifndef NBNXN_OPENCL_TYPES_H
#define NBNXN_OPENCL_TYPES_H

#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/mdlib/nbnxn_gpu_types_common.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/utility/real.h"

/* kernel does #include "gromacs/math/utilities.h" */
/* Move the actual useful stuff here: */

//! Define 1/sqrt(pi)
#define M_FLOAT_1_SQRTPI 0.564189583547756f

/*! \brief Macros defining platform-dependent defaults for the prune kernel's j4 processing concurrency.
 *
 *  The GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY macro allows compile-time override.
 */
/*! @{ */
#ifndef GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY
#define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_AMD       4
#define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_NVIDIA    4
#define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_DEFAULT   4
#else
#define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_AMD       GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY
#define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_NVIDIA    GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY
#define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_DEFAULT   GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY
#endif
/*! @} */
/*! \brief Constants for platform-dependent defaults for the prune kernel's j4 processing concurrency.
 *
 *  Initialized using macros that can be overridden at compile-time (using #GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY).
 */
/*! @{ */
const int c_oclPruneKernelJ4ConcurrencyAMD     = GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_AMD;
const int c_oclPruneKernelJ4ConcurrencyNVIDIA  = GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_NVIDIA;
const int c_oclPruneKernelJ4ConcurrencyDefault = GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY_DEFAULT;
/*! @} */

/*! \brief Returns the j4 processing concurrency parameter for the vendor \p vendorId
 *  \param vendorId takes values from #ocl_vendor_id_t.
 */
static inline int getOclPruneKernelJ4Concurrency(int vendorId)
{
    assert(vendorId < OCL_VENDOR_UNKNOWN);
    switch (vendorId)
    {
        case OCL_VENDOR_AMD:    return c_oclPruneKernelJ4ConcurrencyAMD;     break;
        case OCL_VENDOR_NVIDIA: return c_oclPruneKernelJ4ConcurrencyNVIDIA;  break;
        default:                return c_oclPruneKernelJ4ConcurrencyDefault; break;
    }
}


#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Electrostatic OpenCL kernel flavors.
 *
 *  Types of electrostatics implementations available in the OpenCL non-bonded
 *  force kernels. These represent both the electrostatics types implemented
 *  by the kernels (cut-off, RF, and Ewald - a subset of what's defined in
 *  enums.h) as well as encode implementation details analytical/tabulated
 *  and single or twin cut-off (for Ewald kernels).
 *  Note that the cut-off and RF kernels have only analytical flavor and unlike
 *  in the CPU kernels, the tabulated kernels are ATM Ewald-only.
 *
 *  The row-order of pointers to different electrostatic kernels defined in
 *  nbnxn_cuda.cu by the nb_*_kfunc_ptr function pointer table
 *  should match the order of enumerated types below.
 */
enum eelOcl {
    eelOclCUT, eelOclRF, eelOclEWALD_TAB, eelOclEWALD_TAB_TWIN, eelOclEWALD_ANA, eelOclEWALD_ANA_TWIN, eelOclNR
};

/*! \brief VdW OpenCL kernel flavors.
 *
 * The enumerates values correspond to the LJ implementations in the OpenCL non-bonded
 * kernels.
 *
 * The column-order of pointers to different electrostatic kernels defined in
 * nbnxn_cuda.cu by the nb_*_kfunc_ptr function pointer table
 * should match the order of enumerated types below.
 */
enum evdwOcl {
    evdwOclCUT, evdwOclCUTCOMBGEOM, evdwOclCUTCOMBLB, evdwOclFSWITCH, evdwOclPSWITCH, evdwOclEWALDGEOM, evdwOclEWALDLB, evdwOclNR
};

/*! \brief Pruning kernel flavors.
 *
 * The values correspond to the first call of the pruning post-list generation
 * and the rolling pruning, respectively.
 */
enum ePruneKind {
    epruneFirst, epruneRolling, ePruneNR
};

/*! \internal
 * \brief Staging area for temporary data downloaded from the GPU.
 *
 *  The energies/shift forces get downloaded here first, before getting added
 *  to the CPU-side aggregate values.
 */
typedef struct cl_nb_staging
{
    float    *e_lj;           /**< LJ energy                       */
    float    *e_el;           /**< electrostatic energy            */
    float   (*fshift)[3];     /**< float3 buffer with shift forces */
} cl_nb_staging_t;

/*! \internal
 * \brief Nonbonded atom data - both inputs and outputs.
 */
typedef struct cl_atomdata
{
    int         natoms;              /**< number of atoms                              */
    int         natoms_local;        /**< number of local atoms                        */
    int         nalloc;              /**< allocation size for the atom data (xq, f)    */

    cl_mem      xq;                  /**< float4 buffer with atom coordinates + charges, size natoms */

    cl_mem      f;                   /**< float3 buffer with force output array, size natoms         */
    size_t      f_elem_size;         /**< Size in bytes for one element of f buffer      */

    cl_mem      e_lj;                /**< LJ energy output, size 1                       */
    cl_mem      e_el;                /**< Electrostatics energy input, size 1            */

    cl_mem      fshift;              /**< float3 buffer with shift forces                */
    size_t      fshift_elem_size;    /**< Size in bytes for one element of fshift buffer */

    int         ntypes;              /**< number of atom types                           */
    cl_mem      atom_types;          /**< int buffer with atom type indices, size natoms */
    cl_mem      lj_comb;             /**< float2 buffer with sqrt(c6),sqrt(c12), size natoms */

    cl_mem      shift_vec;           /**< float3 buffer with shifts values               */
    size_t      shift_vec_elem_size; /**< Size in bytes for one element of shift_vec buffer */

    cl_bool     bShiftVecUploaded;   /**< true if the shift vector has been uploaded  */
} cl_atomdata_t;

/*! \internal
 * \brief Parameters required for the OpenCL nonbonded calculations.
 */
typedef struct cl_nbparam
{

    int             eeltype;           /**< type of electrostatics, takes values from #eelOcl */
    int             vdwtype;           /**< type of VdW impl., takes values from #evdwOcl     */

    float           epsfac;            /**< charge multiplication factor                      */
    float           c_rf;              /**< Reaction-field/plain cutoff electrostatics const. */
    float           two_k_rf;          /**< Reaction-field electrostatics constant            */
    float           ewald_beta;        /**< Ewald/PME parameter                               */
    float           sh_ewald;          /**< Ewald/PME correction term substracted from the direct-space potential */
    float           sh_lj_ewald;       /**< LJ-Ewald/PME correction term added to the correction potential        */
    float           ewaldcoeff_lj;     /**< LJ-Ewald/PME coefficient                          */

    float           rcoulomb_sq;       /**< Coulomb cut-off squared                           */

    float           rvdw_sq;           /**< VdW cut-off squared                               */
    float           rvdw_switch;       /**< VdW switched cut-off                              */
    float           rlistOuter_sq;     /**< Full, outer pair-list cut-off squared             */
    float           rlistInner_sq;     /**< Inner, dynamic pruned pair-list cut-off squared   */
    bool            useDynamicPruning; /**< True if we use dynamic pair-list pruning          */

    shift_consts_t  dispersion_shift;  /**< VdW shift dispersion constants           */
    shift_consts_t  repulsion_shift;   /**< VdW shift repulsion constants            */
    switch_consts_t vdw_switch;        /**< VdW switch constants                     */

    /* LJ non-bonded parameters - accessed through texture memory */
    cl_mem                  nbfp_climg2d;      /**< nonbonded parameter table with C6/C12 pairs per atom type-pair, 2*ntype^2 elements */
    cl_mem                  nbfp_comb_climg2d; /**< nonbonded parameter table per atom type, 2*ntype elements                          */

    /* Ewald Coulomb force table data - accessed through texture memory */
    float                  coulomb_tab_scale;   /**< table scale/spacing                        */
    cl_mem                 coulomb_tab_climg2d; /**< pointer to the table in the device memory  */
} cl_nbparam_t;

/*! \internal
 * \brief Data structure shared between the OpenCL device code and OpenCL host code
 *
 * Must not contain OpenCL objects (buffers)
 * TODO: review, improve */
typedef struct cl_nbparam_params
{

    int             eeltype;          /**< type of electrostatics, takes values from #eelCu */
    int             vdwtype;          /**< type of VdW impl., takes values from #evdwCu     */

    float           epsfac;           /**< charge multiplication factor                      */
    float           c_rf;             /**< Reaction-field/plain cutoff electrostatics const. */
    float           two_k_rf;         /**< Reaction-field electrostatics constant            */
    float           ewald_beta;       /**< Ewald/PME parameter                               */
    float           sh_ewald;         /**< Ewald/PME correction term substracted from the direct-space potential */
    float           sh_lj_ewald;      /**< LJ-Ewald/PME correction term added to the correction potential        */
    float           ewaldcoeff_lj;    /**< LJ-Ewald/PME coefficient                          */

    float           rcoulomb_sq;      /**< Coulomb cut-off squared                           */

    float           rvdw_sq;          /**< VdW cut-off squared                               */
    float           rvdw_switch;      /**< VdW switched cut-off                              */
    float           rlistOuter_sq;    /**< Full, outer pair-list cut-off squared             */
    float           rlistInner_sq;    /**< Inner, dynamic pruned pair-list cut-off squared   */

    shift_consts_t  dispersion_shift; /**< VdW shift dispersion constants           */
    shift_consts_t  repulsion_shift;  /**< VdW shift repulsion constants            */
    switch_consts_t vdw_switch;       /**< VdW switch constants                     */

    /* Ewald Coulomb force table data - accessed through texture memory */
    float                  coulomb_tab_scale;  /**< table scale/spacing                        */
} cl_nbparam_params_t;


/*! \internal
 * \brief Pair list data.
 */
typedef struct cl_plist
{
    int              na_c;         /**< number of atoms per cluster                  */

    int              nsci;         /**< size of sci, # of i clusters in the list     */
    int              sci_nalloc;   /**< allocation size of sci                       */
    cl_mem           sci;          /**< list of i-cluster ("super-clusters").
                                        It contains elements of type nbnxn_sci_t     */

    int              ncj4;         /**< total # of 4*j clusters                      */
    int              cj4_nalloc;   /**< allocation size of cj4                       */
    cl_mem           cj4;          /**< 4*j cluster list, contains j cluster number and
                                        index into the i cluster list.
                                        It contains elements of type nbnxn_cj4_t     */
    int              nimask;       /**< # of 4*j clusters * # of warps               */
    int              imask_nalloc; /**< allocation size of imask                     */
    cl_mem           imask;        /**< imask for 2 warps for each 4*j cluster group */
    cl_mem           excl;         /**< atom interaction bits
                                        It contains elements of type nbnxn_excl_t    */
    int              nexcl;        /**< count for excl                               */
    int              excl_nalloc;  /**< allocation size of excl                      */

    /* parameter+variables for normal and rolling pruning */
    bool             haveFreshList;          /**< true after search, indictes that initial pruning with outer prunning is needed */
    int              rollingPruningNumParts; /**< the number of parts/steps over which one cyle of roling pruning takes places */
    int              rollingPruningPart;     /**< the next part to which the roling pruning needs to be applied */
}cl_plist_t;


/** \internal
 * \brief Typedef of actual timer type.
 */
typedef struct nbnxn_gpu_timers_t cl_timers_t;

/*! \internal
 * \brief Main data structure for OpenCL nonbonded force calculations.
 */
struct gmx_nbnxn_ocl_t
{
    const gmx_device_info_t          *dev_info;    /**< OpenCL device information                              */
    struct gmx_device_runtime_data_t *dev_rundata; /**< OpenCL runtime data (context, kernels)                 */

    /**< Pointers to non-bonded kernel functions
     * organized similar with nb_kfunc_xxx arrays in nbnxn_ocl.cpp */
    ///@{
    cl_kernel           kernel_noener_noprune_ptr[eelOclNR][evdwOclNR];
    cl_kernel           kernel_ener_noprune_ptr[eelOclNR][evdwOclNR];
    cl_kernel           kernel_noener_prune_ptr[eelOclNR][evdwOclNR];
    cl_kernel           kernel_ener_prune_ptr[eelOclNR][evdwOclNR];
    ///@}
    cl_kernel           kernel_pruneonly[ePruneNR]; /**< prune kernels, ePruneKind defined the kernel kinds */

    bool                bPrefetchLjParam;           /**< true if prefetching fg i-atom LJ parameters should be used in the kernels */

    /**< auxiliary kernels implementing memset-like functions */
    ///@{
    cl_kernel           kernel_memset_f;
    cl_kernel           kernel_memset_f2;
    cl_kernel           kernel_memset_f3;
    cl_kernel           kernel_zero_e_fshift;
    ///@}

    cl_bool             bUseTwoStreams;        /**< true if doing both local/non-local NB work on GPU          */
    cl_bool             bNonLocalStreamActive; /**< true indicates that the nonlocal_done event was enqueued   */

    cl_atomdata_t      *atdat;                 /**< atom data                                                  */
    cl_nbparam_t       *nbparam;               /**< parameters required for the non-bonded calc.               */
    cl_plist_t         *plist[2];              /**< pair-list data structures (local and non-local)            */
    cl_nb_staging_t     nbst;                  /**< staging area where fshift/energies get downloaded          */

    cl_mem              debug_buffer;          /**< debug buffer */

    cl_command_queue    stream[2];             /**< local and non-local GPU queues                             */

    /** events used for synchronization */
    cl_event nonlocal_done;                     /**< event triggered when the non-local non-bonded kernel
                                                   is done (and the local transfer can proceed) */
    cl_event misc_ops_and_local_H2D_done;       /**< event triggered when the tasks issued in
                                                   the local stream that need to precede the
                                                   non-local force calculations are done
                                                   (e.g. f buffer 0-ing, local x/q H2D) */

    cl_bool                           bDoTime;  /**< True if event-based timing is enabled.                     */
    cl_timers_t                      *timers;   /**< OpenCL event-based timers.                                 */
    struct gmx_wallclock_gpu_nbnxn_t *timings;  /**< Timing data. TODO: deprecate this and query timers for accumulated data instead */
};

#ifdef __cplusplus
}
#endif

#endif  /* NBNXN_OPENCL_TYPES_H */
