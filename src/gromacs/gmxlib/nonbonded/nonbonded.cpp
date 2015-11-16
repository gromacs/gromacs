/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "nonbonded.h"

#include "config.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "thread_mpi/threads.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"
#include "gromacs/gmxlib/nonbonded/nb_generic.h"
#include "gromacs/gmxlib/nonbonded/nb_generic_cg.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* Different default (c) and SIMD instructions interaction-specific kernels */
#include "gromacs/gmxlib/nonbonded/nb_kernel_c/nb_kernel_c.h"

#if GMX_SIMD_X86_SSE2 && !GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_sse2_single/nb_kernel_sse2_single.h"
#endif
#if GMX_SIMD_X86_SSE4_1 && !GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_sse4_1_single/nb_kernel_sse4_1_single.h"
#endif
#if GMX_SIMD_X86_AVX_128_FMA && !GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_avx_128_fma_single/nb_kernel_avx_128_fma_single.h"
#endif
#if (GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256) && !GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_avx_256_single/nb_kernel_avx_256_single.h"
#endif
#if GMX_SIMD_X86_SSE2 && GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_sse2_double/nb_kernel_sse2_double.h"
#endif
#if GMX_SIMD_X86_SSE4_1 && GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_sse4_1_double/nb_kernel_sse4_1_double.h"
#endif
#if GMX_SIMD_X86_AVX_128_FMA && GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_avx_128_fma_double/nb_kernel_avx_128_fma_double.h"
#endif
#if (GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256) && GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_avx_256_double/nb_kernel_avx_256_double.h"
#endif
#if GMX_SIMD_SPARC64_HPC_ACE && GMX_DOUBLE
#    include "gromacs/gmxlib/nonbonded/nb_kernel_sparc64_hpc_ace_double/nb_kernel_sparc64_hpc_ace_double.h"
#endif


static tMPI_Thread_mutex_t nonbonded_setup_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
static gmx_bool            nonbonded_setup_done  = FALSE;


void
gmx_nonbonded_setup(t_forcerec *   fr,
                    gmx_bool       bGenericKernelOnly)
{
    tMPI_Thread_mutex_lock(&nonbonded_setup_mutex);
    /* Here we are guaranteed only one thread made it. */
    if (nonbonded_setup_done == FALSE)
    {
        if (bGenericKernelOnly == FALSE)
        {
            /* Add the generic kernels to the structure stored statically in nb_kernel.c */
            nb_kernel_list_add_kernels(kernellist_c, kernellist_c_size);

            if (!(fr != NULL && fr->use_simd_kernels == FALSE))
            {
                /* Add interaction-specific kernels for different architectures */
                /* Single precision */
#if GMX_SIMD_X86_SSE2 && !GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_sse2_single, kernellist_sse2_single_size);
#endif
#if GMX_SIMD_X86_SSE4_1 && !GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_sse4_1_single, kernellist_sse4_1_single_size);
#endif
#if GMX_SIMD_X86_AVX_128_FMA && !GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_avx_128_fma_single, kernellist_avx_128_fma_single_size);
#endif
#if (GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256) && !GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_avx_256_single, kernellist_avx_256_single_size);
#endif
                /* Double precision */
#if GMX_SIMD_X86_SSE2 && GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_sse2_double, kernellist_sse2_double_size);
#endif
#if GMX_SIMD_X86_SSE4_1 && GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_sse4_1_double, kernellist_sse4_1_double_size);
#endif
#if GMX_SIMD_X86_AVX_128_FMA && GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_avx_128_fma_double, kernellist_avx_128_fma_double_size);
#endif
#if (GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256) && GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_avx_256_double, kernellist_avx_256_double_size);
#endif
#if GMX_SIMD_SPARC64_HPC_ACE && GMX_DOUBLE
                nb_kernel_list_add_kernels(kernellist_sparc64_hpc_ace_double, kernellist_sparc64_hpc_ace_double_size);
#endif
                ; /* empty statement to avoid a completely empty block */
            }
        }
        /* Create a hash for faster lookups */
        nb_kernel_list_hash_init();

        nonbonded_setup_done = TRUE;
    }
    tMPI_Thread_mutex_unlock(&nonbonded_setup_mutex);
}



void
gmx_nonbonded_set_kernel_pointers(FILE *log, t_nblist *nl, gmx_bool bElecAndVdwSwitchDiffers)
{
    const char *     elec;
    const char *     elec_mod;
    const char *     vdw;
    const char *     vdw_mod;
    const char *     geom;
    const char *     other;

    struct
    {
        const char *  arch;
        int           simd_padding_width;
    }
    arch_and_padding[] =
    {
        /* Single precision */
#if (GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256) && !GMX_DOUBLE
        { "avx_256_single", 8 },
#endif
#if GMX_SIMD_X86_AVX_128_FMA && !GMX_DOUBLE
        { "avx_128_fma_single", 4 },
#endif
#if GMX_SIMD_X86_SSE4_1 && !GMX_DOUBLE
        { "sse4_1_single", 4 },
#endif
#if GMX_SIMD_X86_SSE2 && !GMX_DOUBLE
        { "sse2_single", 4 },
#endif
        /* Double precision */
#if (GMX_SIMD_X86_AVX_256 || GMX_SIMD_X86_AVX2_256) && GMX_DOUBLE
        { "avx_256_double", 4 },
#endif
#if GMX_SIMD_X86_AVX_128_FMA && GMX_DOUBLE
        /* Sic. Double precision 2-way SIMD does not require neighbor list padding,
         * since the kernels execute a loop unrolled a factor 2, followed by
         * a possible single odd-element epilogue.
         */
        { "avx_128_fma_double", 1 },
#endif
#if GMX_SIMD_X86_SSE2 && GMX_DOUBLE
        /* No padding - see comment above */
        { "sse2_double", 1 },
#endif
#if GMX_SIMD_X86_SSE4_1 && GMX_DOUBLE
        /* No padding - see comment above */
        { "sse4_1_double", 1 },
#endif
#if GMX_SIMD_SPARC64_HPC_ACE && GMX_DOUBLE
        /* No padding - see comment above */
        { "sparc64_hpc_ace_double", 1 },
#endif
        { "c", 1 },
    };
    int              narch = asize(arch_and_padding);
    int              i;

    if (nonbonded_setup_done == FALSE)
    {
        /* We typically call this setup routine before starting timers,
         * but if that has not been done for whatever reason we do it now.
         */
        gmx_nonbonded_setup(NULL, FALSE);
    }

    /* Not used yet */
    other = "";

    nl->kernelptr_vf = NULL;
    nl->kernelptr_v  = NULL;
    nl->kernelptr_f  = NULL;

    elec     = gmx_nbkernel_elec_names[nl->ielec];
    elec_mod = eintmod_names[nl->ielecmod];
    vdw      = gmx_nbkernel_vdw_names[nl->ivdw];
    vdw_mod  = eintmod_names[nl->ivdwmod];
    geom     = gmx_nblist_geometry_names[nl->igeometry];

    if (nl->type == GMX_NBLIST_INTERACTION_FREE_ENERGY)
    {
        nl->kernelptr_vf       = (void *) gmx_nb_free_energy_kernel;
        nl->kernelptr_f        = (void *) gmx_nb_free_energy_kernel;
        nl->simd_padding_width = 1;
    }
    else if (!gmx_strcasecmp_min(geom, "CG-CG"))
    {
        nl->kernelptr_vf       = (void *) gmx_nb_generic_cg_kernel;
        nl->kernelptr_f        = (void *) gmx_nb_generic_cg_kernel;
        nl->simd_padding_width = 1;
    }
    else
    {
        /* Try to find a specific kernel first */

        for (i = 0; i < narch && nl->kernelptr_vf == NULL; i++)
        {
            nl->kernelptr_vf       = (void *) nb_kernel_list_findkernel(log, arch_and_padding[i].arch, elec, elec_mod, vdw, vdw_mod, geom, other, "PotentialAndForce");
            nl->simd_padding_width = arch_and_padding[i].simd_padding_width;
        }
        for (i = 0; i < narch && nl->kernelptr_f == NULL; i++)
        {
            nl->kernelptr_f        = (void *) nb_kernel_list_findkernel(log, arch_and_padding[i].arch, elec, elec_mod, vdw, vdw_mod, geom, other, "Force");
            nl->simd_padding_width = arch_and_padding[i].simd_padding_width;

            /* If there is not force-only optimized kernel, is there a potential & force one? */
            if (nl->kernelptr_f == NULL)
            {
                nl->kernelptr_f        = (void *) nb_kernel_list_findkernel(NULL, arch_and_padding[i].arch, elec, elec_mod, vdw, vdw_mod, geom, other, "PotentialAndForce");
                nl->simd_padding_width = arch_and_padding[i].simd_padding_width;
            }
        }

        /* For now, the accelerated kernels cannot handle the combination of switch functions for both
         * electrostatics and VdW that use different switch radius or switch cutoff distances
         * (both of them enter in the switch function calculation). This would require
         * us to evaluate two completely separate switch functions for every interaction.
         * Instead, we disable such kernels by setting the pointer to NULL.
         * This will cause the generic kernel (which can handle it) to be called instead.
         *
         * Note that we typically already enable tabulated coulomb interactions for this case,
         * so this is mostly a safe-guard to make sure we call the generic kernel if the
         * tables are disabled.
         */
        if ((nl->ielec != GMX_NBKERNEL_ELEC_NONE) && (nl->ielecmod == eintmodPOTSWITCH) &&
            (nl->ivdw  != GMX_NBKERNEL_VDW_NONE)  && (nl->ivdwmod  == eintmodPOTSWITCH) &&
            bElecAndVdwSwitchDiffers)
        {
            nl->kernelptr_vf = NULL;
            nl->kernelptr_f  = NULL;
        }

        /* Give up, pick a generic one instead.
         * We only do this for particle-particle kernels; by leaving the water-optimized kernel
         * pointers to NULL, the water optimization will automatically be disabled for this interaction.
         */
        if (nl->kernelptr_vf == NULL && !gmx_strcasecmp_min(geom, "Particle-Particle"))
        {
            nl->kernelptr_vf       = (void *) gmx_nb_generic_kernel;
            nl->kernelptr_f        = (void *) gmx_nb_generic_kernel;
            nl->simd_padding_width = 1;
            if (debug)
            {
                fprintf(debug,
                        "WARNING - Slow generic NB kernel used for neighborlist with\n"
                        "    Elec: '%s', Modifier: '%s'\n"
                        "    Vdw:  '%s', Modifier: '%s'\n",
                        elec, elec_mod, vdw, vdw_mod);
            }
        }
    }
    return;
}

void do_nonbonded(t_forcerec *fr,
                  rvec x[], rvec f_shortrange[], t_mdatoms *mdatoms, t_blocka *excl,
                  gmx_grppairener_t *grppener,
                  t_nrnb *nrnb, real *lambda, real *dvdl,
                  int nls, int eNL, int flags)
{
    t_nblist *        nlist;
    int               n, n0, n1, i, i0, i1;
    t_nblists *       nblists;
    nb_kernel_data_t  kernel_data;
    nb_kernel_t *     kernelptr = NULL;
    rvec *            f;

    kernel_data.flags                   = flags;
    kernel_data.exclusions              = excl;
    kernel_data.lambda                  = lambda;
    kernel_data.dvdl                    = dvdl;

    if (fr->bAllvsAll)
    {
        gmx_incons("All-vs-all kernels have not been implemented in version 4.6");
        return;
    }

    if (eNL >= 0)
    {
        i0 = eNL;
        i1 = i0+1;
    }
    else
    {
        i0 = 0;
        i1 = eNL_NR;
    }

    if (nls >= 0)
    {
        n0 = nls;
        n1 = nls+1;
    }
    else
    {
        n0 = 0;
        n1 = fr->nnblists;
    }

    for (n = n0; (n < n1); n++)
    {
        nblists = &fr->nblists[n];

        /* Tabulated kernels hard-code a lot of assumptions about the
         * structure of these tables, but that's not worth fixing with
         * the group scheme due for removal soon. As a token
         * improvement, this assertion will stop code segfaulting if
         * someone assumes that extending the group-scheme table-type
         * enumeration is something that GROMACS supports. */
        /* cppcheck-suppress duplicateExpression */
        assert(etiNR == 3);

        kernel_data.table_elec              = nblists->table_elec;
        kernel_data.table_vdw               = nblists->table_vdw;
        kernel_data.table_elec_vdw          = nblists->table_elec_vdw;

        {
            {
                /* Short-range */
                if (!(flags & GMX_NONBONDED_DO_SR))
                {
                    continue;
                }
                kernel_data.energygrp_elec          = grppener->ener[egCOULSR];
                kernel_data.energygrp_vdw           = grppener->ener[fr->bBHAM ? egBHAMSR : egLJSR];
                kernel_data.energygrp_polarization  = grppener->ener[egGB];
                nlist = nblists->nlist_sr;
                f                                   = f_shortrange;
            }

            for (i = i0; (i < i1); i++)
            {
                if (nlist[i].nri > 0)
                {
                    if (flags & GMX_NONBONDED_DO_POTENTIAL)
                    {
                        /* Potential and force */
                        kernelptr = (nb_kernel_t *)nlist[i].kernelptr_vf;
                    }
                    else
                    {
                        /* Force only, no potential */
                        kernelptr = (nb_kernel_t *)nlist[i].kernelptr_f;
                    }

                    if (nlist[i].type != GMX_NBLIST_INTERACTION_FREE_ENERGY && (flags & GMX_NONBONDED_DO_FOREIGNLAMBDA))
                    {
                        /* We don't need the non-perturbed interactions */
                        continue;
                    }
                    /* Neighborlists whose kernelptr==NULL will always be empty */
                    if (kernelptr != NULL)
                    {
                        (*kernelptr)(&(nlist[i]), x, f, fr, mdatoms, &kernel_data, nrnb);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Non-empty neighborlist does not have any kernel pointer assigned.");
                    }
                }
            }
        }
    }
}
