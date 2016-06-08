/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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

#include <string>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/ocl_compiler.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/nbnxn_gpu_jit_support.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "nbnxn_ocl_types.h"

/*! \brief Stringifies the input argument
 */
#define STRINGIFY_PARAM(c) #c

/*! \brief Stringifies the result of expansion of a macro argument
 */
#define STRINGIFY_MACRO(c) STRINGIFY_PARAM(c)

/*! \brief Array of the defines needed to generate a specific eel flavour
 *
 * The twin-cutoff entries are not normally used, because those setups are
 * not available to the user. FastGen takes care of generating both
 * single- and twin-cutoff versions because PME tuning might need both.
 */
static const char * kernel_electrostatic_family_definitions[] =
{
    " -DEL_CUTOFF -DEELNAME=_ElecCut",
    " -DEL_RF -DEELNAME=_ElecRF",
    " -DEL_EWALD_TAB -DEELNAME=_ElecEwQSTab",
    " -DEL_EWALD_TAB -DVDW_CUTOFF_CHECK -DEELNAME=_ElecEwQSTabTwinCut",
    " -DEL_EWALD_ANA -DEELNAME=_ElecEw",
    " -DEL_EWALD_ANA -DVDW_CUTOFF_CHECK -DEELNAME=_ElecEwTwinCut"
};

/*! \brief Array of the defines needed to generate a specific vdw flavour
 */
static const char * kernel_VdW_family_definitions[] =
{
    " -DVDWNAME=_VdwLJ",
    " -DLJ_COMB_GEOM -DVDWNAME=_VdwLJCombGeom",
    " -DLJ_COMB_LB  -DVDWNAME=_VdwLJCombLB",
    " -DLJ_FORCE_SWITCH -DVDWNAME=_VdwLJFsw",
    " -DLJ_POT_SWITCH -DVDWNAME=_VdwLJPsw",
    " -DLJ_EWALD_COMB_GEOM -DVDWNAME=_VdwLJEwCombGeom",
    " -DLJ_EWALD_COMB_LB -DVDWNAME=_VdwLJEwCombLB"
};

/*! \brief Returns a string with the compiler defines required to avoid all flavour generation
 *
 * For example if flavour eelOclRF with evdwOclFSWITCH, the output will be such that the corresponding
 * kernel flavour is generated:
 * -DGMX_OCL_FASTGEN          (will replace flavour generator nbnxn_ocl_kernels.clh with nbnxn_ocl_kernels_fastgen.clh)
 * -DEL_RF                    (The eelOclRF flavour)
 * -DEELNAME=_ElecRF          (The first part of the generated kernel name )
 * -DLJ_EWALD_COMB_GEOM       (The evdwOclFSWITCH flavour)
 * -DVDWNAME=_VdwLJEwCombGeom (The second part of the generated kernel name )
 *
 * prune/energy are still generated as originally. It is only the the flavour-level that has changed, so that
 * only the required flavour for the simulation is compiled.
 *
 * If eeltype is single-range Ewald, then we need to add the
 * twin-cutoff flavour kernels to the JIT, because PME tuning might
 * need it. This path sets -DGMX_OCL_FASTGEN_ADD_TWINCUT, which
 * triggers the use of nbnxn_ocl_kernels_fastgen_add_twincut.clh. This
 * hard-codes the generation of extra kernels that have the same base
 * flavour, and add the required -DVDW_CUTOFF_CHECK and "TwinCut" to
 * the kernel name.
 *
 * If FastGen is not active, then nothing needs to be returned. The
 * JIT defaults to compiling all kernel flavours.
 *
 * \param[in]  bFastGen    Whether FastGen should be used
 * \param[in]  eeltype     Electrostatics kernel flavour for FastGen
 * \param[in]  vdwtype     VDW kernel flavour for FastGen
 * \return                 String with the defines if FastGen is active
 *
 * \throws std::bad_alloc if out of memory
 */
static std::string
makeDefinesForKernelTypes(bool bFastGen,
                          int  eeltype,
                          int  vdwtype)
{
    std::string defines_for_kernel_types;

    if (bFastGen)
    {
        bool bIsEwaldSingleCutoff = (eeltype == eelOclEWALD_TAB ||
                                     eeltype == eelOclEWALD_ANA);

        if (bIsEwaldSingleCutoff)
        {
            defines_for_kernel_types += "-DGMX_OCL_FASTGEN_ADD_TWINCUT";
        }
        else
        {
            /* This triggers the use of
               nbnxn_ocl_kernels_fastgen.clh. */
            defines_for_kernel_types += "-DGMX_OCL_FASTGEN";
        }
        defines_for_kernel_types += kernel_electrostatic_family_definitions[eeltype];
        defines_for_kernel_types += kernel_VdW_family_definitions[vdwtype];
    }

    return defines_for_kernel_types;
}

/*! \brief Compiles nbnxn kernels for OpenCL GPU given by \p mygpu
 *
 * With OpenCL, a call to this function must precede nbnxn_gpu_init().
 *
 * Doing bFastGen means only the requested kernels are compiled,
 * significantly reducing the total compilation time. If false, all
 * OpenCL kernels are compiled.
 *
 * A fatal error results if compilation fails.
 *
 * \param[inout] nb  Manages OpenCL non-bonded calculations; compiled kernels returned in dev_info members
 *
 * Does not throw
 */
void
nbnxn_gpu_compile_kernels(gmx_nbnxn_ocl_t *nb)
{
    gmx_bool                  bFastGen = TRUE;
    cl_program                program  = nullptr;

    if (getenv("GMX_OCL_NOFASTGEN") != NULL)
    {
        bFastGen = FALSE;
    }

    /* Need to catch std::bad_alloc here and during compilation string
       handling. */
    try
    {
        std::string extraDefines = makeDefinesForKernelTypes(bFastGen,
                                                             nb->nbparam->eeltype,
                                                             nb->nbparam->vdwtype);

        /* Here we pass macros and static const int variables defined in include
         * files outside the nbnxn_ocl as macros, to avoid including those files
         * in the JIT compilation that happens at runtime.
         */
        extraDefines += gmx::formatString(
                    " -DCENTRAL=%d -DNBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER=%d -DNBNXN_GPU_CLUSTER_SIZE=%d -DNBNXN_GPU_JGROUP_SIZE=%d -DNBNXN_MIN_RSQ=%s %s",
                    CENTRAL,                                 /* Defined in ishift.h */
                    c_nbnxnGpuNumClusterPerSupercluster,     /* Defined in nbnxn_pairlist.h */
                    c_nbnxnGpuClusterSize,                   /* Defined in nbnxn_pairlist.h */
                    c_nbnxnGpuJgroupSize,                    /* Defined in nbnxn_pairlist.h */
                    STRINGIFY_MACRO(NBNXN_MIN_RSQ)           /* Defined in nbnxn_consts.h */
                                                             /* NBNXN_MIN_RSQ passed as string to avoid
                                                                floating point representation problems with sprintf */
                    , (nb->bPrefetchLjParam) ? "-DIATYPE_SHMEM" : ""
                    );


        try
        {
            /* TODO when we have a proper MPI-aware logging module,
               the log output here should be written there */
            program = gmx::ocl::compileProgram(stderr,
                                               "nbnxn_ocl_kernels.cl",
                                               extraDefines,
                                               nb->dev_rundata->context,
                                               nb->dev_info->ocl_gpu_id.ocl_device_id,
                                               nb->dev_info->vendor_e);
        }
        catch (gmx::GromacsException &e)
        {
            e.prependContext(gmx::formatString("Failed to compile NBNXN kernels for GPU #%s\n",
                                               nb->dev_info->device_name));
            throw;
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    nb->dev_rundata->program = program;
}
