/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_THREAD_MPI
#include <thread_mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "txtdump.h"
#include "smalloc.h"
#include "ns.h"
#include "vec.h"
#include "maths.h"
#include "macros.h"
#include "string2.h"
#include "force.h"
#include "names.h"
#include "main.h"
#include "xvgr.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "force.h"
#include "bondf.h"
#include "nrnb.h"
#include "smalloc.h"
#include "nonbonded.h"

#include "nb_kernel.h"
#include "nb_free_energy.h"
#include "nb_generic.h"
#include "nb_generic_cg.h"
#include "nb_generic_adress.h"

/* Different default (c) and accelerated interaction-specific kernels */
#include "nb_kernel_c/nb_kernel_c.h"

#if (defined GMX_CPU_ACCELERATION_X86_SSE2) && !(defined GMX_DOUBLE)
#    include "nb_kernel_sse2_single/nb_kernel_sse2_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1) && !(defined GMX_DOUBLE)
#    include "nb_kernel_sse4_1_single/nb_kernel_sse4_1_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA) && !(defined GMX_DOUBLE)
#    include "nb_kernel_avx_128_fma_single/nb_kernel_avx_128_fma_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256) && !(defined GMX_DOUBLE)
#    include "nb_kernel_avx_256_single/nb_kernel_avx_256_single.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE2 && defined GMX_DOUBLE)
#    include "nb_kernel_sse2_double/nb_kernel_sse2_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1 && defined GMX_DOUBLE)
#    include "nb_kernel_sse4_1_double/nb_kernel_sse4_1_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA && defined GMX_DOUBLE)
#    include "nb_kernel_avx_128_fma_double/nb_kernel_avx_128_fma_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256 && defined GMX_DOUBLE)
#    include "nb_kernel_avx_256_double/nb_kernel_avx_256_double.h"
#endif
#if (defined GMX_CPU_ACCELERATION_SPARC64_HPC_ACE && defined GMX_DOUBLE)
#    include "nb_kernel_sparc64_hpc_ace_double/nb_kernel_sparc64_hpc_ace_double.h"
#endif


#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t nonbonded_setup_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif
static gmx_bool            nonbonded_setup_done  = FALSE;


void
gmx_nonbonded_setup(FILE *         fplog,
                    t_forcerec *   fr,
                    gmx_bool       bGenericKernelOnly)
{
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&nonbonded_setup_mutex);
#endif
    /* Here we are guaranteed only one thread made it. */
    if (nonbonded_setup_done == FALSE)
    {
        if (bGenericKernelOnly == FALSE)
        {
            /* Add the generic kernels to the structure stored statically in nb_kernel.c */
            nb_kernel_list_add_kernels(kernellist_c, kernellist_c_size);

            if (!(fr != NULL && fr->use_cpu_acceleration == FALSE))
            {
                /* Add interaction-specific kernels for different architectures */
                /* Single precision */
#if (defined GMX_CPU_ACCELERATION_X86_SSE2) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse2_single, kernellist_sse2_single_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse4_1_single, kernellist_sse4_1_single_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_128_fma_single, kernellist_avx_128_fma_single_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256) && !(defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_256_single, kernellist_avx_256_single_size);
#endif
                /* Double precision */
#if (defined GMX_CPU_ACCELERATION_X86_SSE2 && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse2_double, kernellist_sse2_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1 && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sse4_1_double, kernellist_sse4_1_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_128_fma_double, kernellist_avx_128_fma_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256 && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_avx_256_double, kernellist_avx_256_double_size);
#endif
#if (defined GMX_CPU_ACCELERATION_SPARC64_HPC_ACE && defined GMX_DOUBLE)
                nb_kernel_list_add_kernels(kernellist_sparc64_hpc_ace_double,kernellist_sparc64_hpc_ace_double_size);
#endif
                ; /* empty statement to avoid a completely empty block */
            }
        }
        /* Create a hash for faster lookups */
        nb_kernel_list_hash_init();

        nonbonded_setup_done = TRUE;
    }
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&nonbonded_setup_mutex);
#endif
}



void
gmx_nonbonded_set_kernel_pointers(FILE *log, t_nblist *nl)
{
    const char *     elec;
    const char *     elec_mod;
    const char *     vdw;
    const char *     vdw_mod;
    const char *     geom;
    const char *     other;
    const char *     vf;

    struct
    {
        const char *  arch;
        int           simd_padding_width;
    }
    arch_and_padding[] =
    {
        /* Single precision */
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256) && !(defined GMX_DOUBLE)
        { "avx_256_single", 8 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA) && !(defined GMX_DOUBLE)
        { "avx_128_fma_single", 4 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1) && !(defined GMX_DOUBLE)
        { "sse4_1_single", 4 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE2) && !(defined GMX_DOUBLE)
        { "sse2_single", 4 },
#endif
        /* Double precision */
#if (defined GMX_CPU_ACCELERATION_X86_AVX_256 && defined GMX_DOUBLE)
        { "avx_256_double", 4 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_AVX_128_FMA && defined GMX_DOUBLE)
        /* Sic. Double precision 2-way SIMD does not require neighbor list padding,
         * since the kernels execute a loop unrolled a factor 2, followed by
         * a possible single odd-element epilogue.
         */
        { "avx_128_fma_double", 1 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE2 && defined GMX_DOUBLE)
        /* No padding - see comment above */
        { "sse2_double", 1 },
#endif
#if (defined GMX_CPU_ACCELERATION_X86_SSE4_1 && defined GMX_DOUBLE)
        /* No padding - see comment above */
        { "sse4_1_double", 1 },
#endif
#if (defined GMX_CPU_ACCELERATION_SPARC64_HPC_ACE && defined GMX_DOUBLE)
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
        gmx_nonbonded_setup(NULL, NULL, FALSE);
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

    if (nl->type == GMX_NBLIST_INTERACTION_ADRESS)
    {
        nl->kernelptr_vf       = (void *) gmx_nb_generic_adress_kernel;
        nl->kernelptr_f        = (void *) gmx_nb_generic_adress_kernel;
        nl->simd_padding_width = 1;
        return;
    }

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

        /* Give up. If this was a water kernel, leave the pointer as NULL, which
         * will disable water optimization in NS. If it is a particle kernel, set
         * the pointer to the generic NB kernel.
         */
        if (nl->kernelptr_vf == NULL && !gmx_strcasecmp_min(geom,"Particle-Particle"))
        {
            nl->kernelptr_vf       = (void *) gmx_nb_generic_kernel;
            nl->kernelptr_f        = (void *) gmx_nb_generic_kernel;
            nl->simd_padding_width = 1;
            if (debug)
            {
                fprintf(debug,
                        "WARNING - Slow generic NB kernel used for neighborlist with\n"
                        "    Elec: '%s', Modifier: '%s'\n"
                        "    Vdw:  '%s', Modifier: '%s'\n"
                        "    Geom: '%s', Other: '%s'\n\n",
                        elec, elec_mod, vdw, vdw_mod, geom, other);
            }
        }
    }

    return;
}

void do_nonbonded(t_commrec *cr, t_forcerec *fr,
                  rvec x[], rvec f_shortrange[], rvec f_longrange[], t_mdatoms *mdatoms, t_blocka *excl,
                  gmx_grppairener_t *grppener, rvec box_size,
                  t_nrnb *nrnb, real *lambda, real *dvdl,
                  int nls, int eNL, int flags)
{
    t_nblist *        nlist;
    int               n, n0, n1, i, i0, i1, sz, range;
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

        kernel_data.table_elec              = &nblists->table_elec;
        kernel_data.table_vdw               = &nblists->table_vdw;
        kernel_data.table_elec_vdw          = &nblists->table_elec_vdw;

        for (range = 0; range < 2; range++)
        {
            /* Are we doing short/long-range? */
            if (range == 0)
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
            else if (range == 1)
            {
                /* Long-range */
                if (!(flags & GMX_NONBONDED_DO_LR))
                {
                    continue;
                }
                kernel_data.energygrp_elec          = grppener->ener[egCOULLR];
                kernel_data.energygrp_vdw           = grppener->ener[fr->bBHAM ? egBHAMLR : egLJLR];
                kernel_data.energygrp_polarization  = grppener->ener[egGB];
                nlist = nblists->nlist_lr;
                f                                   = f_longrange;
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
                    if(kernelptr != NULL)
                    {
                        (*kernelptr)(&(nlist[i]), x, f, fr, mdatoms, &kernel_data, nrnb);
                    }
                }
            }
        }
    }
}

static void
nb_listed_warning_rlimit(const rvec *x, int ai, int aj, int * global_atom_index, real r, real rlimit)
{
    gmx_warning("Listed nonbonded interaction between particles %d and %d\n"
                "at distance %.3f which is larger than the table limit %.3f nm.\n\n"
                "This is likely either a 1,4 interaction, or a listed interaction inside\n"
                "a smaller molecule you are decoupling during a free energy calculation.\n"
                "Since interactions at distances beyond the table cannot be computed,\n"
                "they are skipped until they are inside the table limit again. You will\n"
                "only see this message once, even if it occurs for several interactions.\n\n"
                "IMPORTANT: This should not happen in a stable simulation, so there is\n"
                "probably something wrong with your system. Only change the table-extension\n"
                "distance in the mdp file if you are really sure that is the reason.\n",
                glatnr(global_atom_index, ai), glatnr(global_atom_index, aj), r, rlimit);

    if (debug)
    {
        fprintf(debug,
                "%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. Ignored\n",
                x[ai][XX], x[ai][YY], x[ai][ZZ], x[aj][XX], x[aj][YY], x[aj][ZZ],
                glatnr(global_atom_index, ai), glatnr(global_atom_index, aj), r);
    }
}



/* This might logically belong better in the nb_generic.c module, but it is only
 * used in do_nonbonded_listed(), and we want it to be inlined there to avoid an
 * extra functional call for every single pair listed in the topology.
 */
static real
nb_evaluate_single(real r2, real tabscale, real *vftab,
                   real qq, real c6, real c12, real *velec, real *vvdw)
{
    real       rinv, r, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VVe, FFe, VVd, FFd, VVr, FFr, fscal;
    int        ntab;

    /* Do the tabulated interactions - first table lookup */
    rinv             = gmx_invsqrt(r2);
    r                = r2*rinv;
    rtab             = r*tabscale;
    ntab             = rtab;
    eps              = rtab-ntab;
    eps2             = eps*eps;
    ntab             = 12*ntab;
    /* Electrostatics */
    Y                = vftab[ntab];
    F                = vftab[ntab+1];
    Geps             = eps*vftab[ntab+2];
    Heps2            = eps2*vftab[ntab+3];
    Fp               = F+Geps+Heps2;
    VVe              = Y+eps*Fp;
    FFe              = Fp+Geps+2.0*Heps2;
    /* Dispersion */
    Y                = vftab[ntab+4];
    F                = vftab[ntab+5];
    Geps             = eps*vftab[ntab+6];
    Heps2            = eps2*vftab[ntab+7];
    Fp               = F+Geps+Heps2;
    VVd              = Y+eps*Fp;
    FFd              = Fp+Geps+2.0*Heps2;
    /* Repulsion */
    Y                = vftab[ntab+8];
    F                = vftab[ntab+9];
    Geps             = eps*vftab[ntab+10];
    Heps2            = eps2*vftab[ntab+11];
    Fp               = F+Geps+Heps2;
    VVr              = Y+eps*Fp;
    FFr              = Fp+Geps+2.0*Heps2;

    *velec           = qq*VVe;
    *vvdw            = c6*VVd+c12*VVr;

    fscal            = -(qq*FFe+c6*FFd+c12*FFr)*tabscale*rinv;

    return fscal;
}


real
do_nonbonded_listed(int ftype, int nbonds,
                    const t_iatom iatoms[], const t_iparams iparams[],
                    const rvec x[], rvec f[], rvec fshift[],
                    const t_pbc *pbc, const t_graph *g,
                    real *lambda, real *dvdl,
                    const t_mdatoms *md,
                    const t_forcerec *fr, gmx_grppairener_t *grppener,
                    int *global_atom_index)
{
    int              ielec, ivdw;
    real             qq, c6, c12;
    rvec             dx;
    ivec             dt;
    int              i, j, itype, ai, aj, gid;
    int              fshift_index;
    real             r2, rinv;
    real             fscal, velec, vvdw;
    real *           energygrp_elec;
    real *           energygrp_vdw;
    static gmx_bool  warned_rlimit = FALSE;
    /* Free energy stuff */
    gmx_bool         bFreeEnergy;
    real             LFC[2], LFV[2], DLF[2], lfac_coul[2], lfac_vdw[2], dlfac_coul[2], dlfac_vdw[2];
    real             qqB, c6B, c12B, sigma2_def, sigma2_min;


    switch (ftype)
    {
        case F_LJ14:
        case F_LJC14_Q:
            energygrp_elec = grppener->ener[egCOUL14];
            energygrp_vdw  = grppener->ener[egLJ14];
            break;
        case F_LJC_PAIRS_NB:
            energygrp_elec = grppener->ener[egCOULSR];
            energygrp_vdw  = grppener->ener[egLJSR];
            break;
        default:
            energygrp_elec = NULL; /* Keep compiler happy */
            energygrp_vdw  = NULL; /* Keep compiler happy */
            gmx_fatal(FARGS, "Unknown function type %d in do_nonbonded14", ftype);
            break;
    }

    if (fr->efep != efepNO)
    {
        /* Lambda factor for state A=1-lambda and B=lambda */
        LFC[0] = 1.0 - lambda[efptCOUL];
        LFV[0] = 1.0 - lambda[efptVDW];
        LFC[1] = lambda[efptCOUL];
        LFV[1] = lambda[efptVDW];

        /*derivative of the lambda factor for state A and B */
        DLF[0] = -1;
        DLF[1] = 1;

        /* precalculate */
        sigma2_def = pow(fr->sc_sigma6_def, 1.0/3.0);
        sigma2_min = pow(fr->sc_sigma6_min, 1.0/3.0);

        for (i = 0; i < 2; i++)
        {
            lfac_coul[i]  = (fr->sc_power == 2 ? (1-LFC[i])*(1-LFC[i]) : (1-LFC[i]));
            dlfac_coul[i] = DLF[i]*fr->sc_power/fr->sc_r_power*(fr->sc_power == 2 ? (1-LFC[i]) : 1);
            lfac_vdw[i]   = (fr->sc_power == 2 ? (1-LFV[i])*(1-LFV[i]) : (1-LFV[i]));
            dlfac_vdw[i]  = DLF[i]*fr->sc_power/fr->sc_r_power*(fr->sc_power == 2 ? (1-LFV[i]) : 1);
        }
    }
    else
    {
        sigma2_min = sigma2_def = 0;
    }

    bFreeEnergy = FALSE;
    for (i = 0; (i < nbonds); )
    {
        itype = iatoms[i++];
        ai    = iatoms[i++];
        aj    = iatoms[i++];
        gid   = GID(md->cENER[ai], md->cENER[aj], md->nenergrp);

        /* Get parameters */
        switch (ftype)
        {
            case F_LJ14:
                bFreeEnergy =
                    (fr->efep != efepNO &&
                     ((md->nPerturbed && (md->bPerturbed[ai] || md->bPerturbed[aj])) ||
                      iparams[itype].lj14.c6A != iparams[itype].lj14.c6B ||
                      iparams[itype].lj14.c12A != iparams[itype].lj14.c12B));
                qq               = md->chargeA[ai]*md->chargeA[aj]*fr->epsfac*fr->fudgeQQ;
                c6               = iparams[itype].lj14.c6A;
                c12              = iparams[itype].lj14.c12A;
                break;
            case F_LJC14_Q:
                qq               = iparams[itype].ljc14.qi*iparams[itype].ljc14.qj*fr->epsfac*iparams[itype].ljc14.fqq;
                c6               = iparams[itype].ljc14.c6;
                c12              = iparams[itype].ljc14.c12;
                break;
            case F_LJC_PAIRS_NB:
                qq               = iparams[itype].ljcnb.qi*iparams[itype].ljcnb.qj*fr->epsfac;
                c6               = iparams[itype].ljcnb.c6;
                c12              = iparams[itype].ljcnb.c12;
                break;
            default:
                /* Cannot happen since we called gmx_fatal() above in this case */
                qq = c6 = c12 = 0; /* Keep compiler happy */
                break;
        }

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0;
        c12 *= 12.0;

        /* Do we need to apply full periodic boundary conditions? */
        if (fr->bMolPBC == TRUE)
        {
            fshift_index = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            fshift_index = CENTRAL;
            rvec_sub(x[ai], x[aj], dx);
        }
        r2           = norm2(dx);

        if (r2 >= fr->tab14.r*fr->tab14.r)
        {
            if (warned_rlimit == FALSE)
            {
                nb_listed_warning_rlimit(x, ai, aj, global_atom_index, sqrt(r2), fr->tab14.r);
                warned_rlimit = TRUE;
            }
            continue;
        }

        if (bFreeEnergy)
        {
            /* Currently free energy is only supported for F_LJ14, so no need to check for that if we got here */
            qqB              = md->chargeB[ai]*md->chargeB[aj]*fr->epsfac*fr->fudgeQQ;
            c6B              = iparams[itype].lj14.c6B*6.0;
            c12B             = iparams[itype].lj14.c12B*12.0;

            fscal            = nb_free_energy_evaluate_single(r2, fr->sc_r_power, fr->sc_alphacoul, fr->sc_alphavdw,
                                                              fr->tab14.scale, fr->tab14.data, qq, c6, c12, qqB, c6B, c12B,
                                                              LFC, LFV, DLF, lfac_coul, lfac_vdw, dlfac_coul, dlfac_vdw,
                                                              fr->sc_sigma6_def, fr->sc_sigma6_min, sigma2_def, sigma2_min, &velec, &vvdw, dvdl);
        }
        else
        {
            /* Evaluate tabulated interaction without free energy */
            fscal            = nb_evaluate_single(r2, fr->tab14.scale, fr->tab14.data, qq, c6, c12, &velec, &vvdw);
        }

        energygrp_elec[gid]  += velec;
        energygrp_vdw[gid]   += vvdw;
        svmul(fscal, dx, dx);

        /* Add the forces */
        rvec_inc(f[ai], dx);
        rvec_dec(f[aj], dx);

        if (g)
        {
            /* Correct the shift forces using the graph */
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            fshift_index = IVEC2IS(dt);
        }
        if (fshift_index != CENTRAL)
        {
            rvec_inc(fshift[fshift_index], dx);
            rvec_dec(fshift[CENTRAL], dx);
        }
    }
    return 0.0;
}
