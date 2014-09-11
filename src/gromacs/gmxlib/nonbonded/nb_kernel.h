/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#ifndef _nb_kernel_h_
#define _nb_kernel_h_

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif


#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/simple.h"

/* Structure to collect kernel data not available in forcerec or mdatoms structures.
 * This is only used inside the nonbonded module.
 */
typedef struct
{
    int                flags;
    t_blocka *         exclusions;
    real *             lambda;
    real *             dvdl;

    /* pointers to tables */
    t_forcetable *     table_elec;
    t_forcetable *     table_vdw;
    t_forcetable *     table_elec_vdw;

    /* potentials */
    real *             energygrp_elec;
    real *             energygrp_vdw;
    real *             energygrp_polarization;
}
nb_kernel_data_t;


typedef void
    nb_kernel_t (t_nblist *                nlist,
                 rvec *                    x,
                 rvec *                    f,
                 t_forcerec *              fr,
                 t_mdatoms *               mdatoms,
                 nb_kernel_data_t *        kernel_data,
                 t_nrnb *                  nrnb);


/* Structure with a kernel pointer and settings. This cannot be abstract
 * since we define the kernel list statically for each architecture in a header,
 * and use it to set up the kernel hash functions to find kernels.
 *
 * The electrostatics/vdw names should be obvious and correspond to the
 * forms of the interactions calculated in this function, and the interaction
 * modifiers describe switch/shift and similar alterations. Geometry refers
 * to whether this kernel calculates interactions between single particles or
 * waters (groups of 3/4 atoms) for better performance. Finally, the VF string
 * selects whether the kernel calculates only potential, only force, or both
 *
 * The allowed values for kernel interactions are described by the
 * enumerated types gmx_nbkernel_elec and gmx_nbkernel_vdw (see types/enums.h).
 * Note that these are deliberately NOT identical to the interactions the
 * user can set, since some user-specified interactions will be tabulated, and
 * Lennard-Jones and Buckingham use different kernels while their setting in
 * the input is decided by nonbonded parameter formats rather than mdp options.
 *
 * The interaction modifiers are described by the eintmod enum type, while
 * the kernel geometry is decided from the neighborlist geometry, which is
 * described by the enum gmx_nblist_kernel_geometry (again, see types/enums.h).
 * The
 *
 * Note that any particular implementation of kernels might not support all of
 * these strings. In fact, some might not be supported by any architecture yet.
 * The whole point of using strings and hashes is that we do not have to define a
 * unique set of strings in a single place. Thus, as long as you implement a
 * corresponding kernel, you could in theory provide any string you want.
 */
typedef struct nb_kernel_info
{
    nb_kernel_t *   kernelptr;
    const char *    kernelname;
    const char *    architecture;     /* e.g. "C", "SSE", "BlueGene", etc. */

    const char *    electrostatics;
    const char *    electrostatics_modifier;
    const char *    vdw;
    const char *    vdw_modifier;
    const char *    geometry;
    const char *    other;  /* Any extra info you want/need to select a kernel */
    const char *    vf;     /* "PotentialAndForce", "Potential", or "Force" */
}
nb_kernel_info_t;


void
nb_kernel_list_add_kernels(nb_kernel_info_t *   new_kernelinfo,
                           int                  new_size);

int
nb_kernel_list_hash_init(void);

/* Return a function pointer to the nonbonded kernel pointer with
 * settings according to the text strings provided. GROMACS does not guarantee
 * the existence of accelerated kernels for any combination, so the return value
 * can be NULL.
 * In that case, you can try a different/lower-level acceleration, and
 * eventually you need to prepare to fall back to generic kernels or change
 * your settings and try again.
 *
 * The names of the text strings are obviously meant to reflect settings in
 * GROMACS, but inside this routine they are merely used as a set of text
 * strings not defined here. The routine will simply compare the arguments with
 * the contents of the corresponding strings in the nb_kernel_list_t structure.
 *
 * This function does not check whether the kernel in question can run on the
 * present architecture since that would require a slow cpuid call for every
 * single invocation.
 */
nb_kernel_t *
nb_kernel_list_findkernel(FILE *              log,
                          const char *        architecture,
                          const char *        electrostatics,
                          const char *        electrostatics_modifier,
                          const char *        vdw,
                          const char *        vdw_modifier,
                          const char *        geometry,
                          const char *        other,
                          const char *        vf);



#ifdef __cplusplus
}
#endif

#endif /* _nb_kernel_h_ */
