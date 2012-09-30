/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * This file is part of GROMACS.
 * Copyright (c) 2012-
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _nb_kernel_h_
#define _nb_kernel_h_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif


#include "types/simple.h"

/* Basic definition of the nonbonded kernel type and a datatype for lists
 * of kernel pointers together with text strings describing the interactions
 * computed by the kernel.
 */



/* Temporary structure to be able to pass 2 GB data through the
 * work data pointer argument.
 */
typedef struct
{
    real      gb_epsilon_solvent;  /* Epsilon for solvent */
    real      epsilon_r;           /* Epsilon for inner dielectric */
    real *    gpol;                /* Polarization energy group */
} gmx_gbdata_t;



typedef void
nb_kernel_t(int *             nri,
int *             iinr,
int *             jindex,
int *             jjnr,
int *             shift,
real *            shiftvec,
real *            fshift,
int *             gid,
real *            pos,
real *            faction,
real *            charge,
real *            facel,
real *            krf,
real *            crf,
real *            vc,
int *             type,
int *             ntype,
real *            vdwparam,
real *            vvdw,
real *            tabscale,
real *            vftab,
real *            invsqrta,
real *            dvda,
real *            gbtabscale,
real *            gbtab,
int *             nthreads,
int *             count,
void *            mtx,
int *             outeriter,
int *             inneriter,
real *            work);


/* Structure with a kernel pointer and settings. This cannot be abstract
 * since we define the kernel list statically for each architecture in a header,
 * and use it to set up the kernel hash functions to find kernels.
 *
 * The electrostatics/vdw names should be obvious and correspond to the
 * forms of the interactions calculated in this function. Geometry refers
 * to whether this kernel calculates interactions between single particles or
 * waters (groups of 3/4 atoms) for better performance. Finally, the VF string
 * selects whether the kernel calculates only potential, only force, or both.
 *
 * As of late summer 2012, the strings below have been defined/reserved.
 * This does not mean all architectures will support all options, or that
 * no other options exist in user-contributed kernels. In fact, some might
 * not be supported by any architecture yet. The whole point of
 * using strings and hashes is that we do not have to define a unique set of
 * strings in a single place. Finally, remember that these strings correspond
 * to how things are calculated in the *kernels*. Many user-level settings
 * might e.g. be translated to tabulated kernels, so there is no 1-to-1
 * correspondence with the mdp options.
 * Since we will move to CamelCase for Gromacs 5.0 (when we use C++), I have
 * used that for the strings already now to facilitate portability.
 *
 * Strings used for Electrostatics:
 *
 * "None"
 * "Coulomb"
 * "CoulombSwitch"
 * "ReactionField"
 * "ReactionFieldCutoff" (Sets interaction to zero outside cutoff)
 * "GeneralizedBorn"
 * "CubicSplineTable"
 * "LinearTable"
 * "Ewald" (Efficient real-space part of PME/Ewald - either analytical or table)
 *
 *
 * Strings used for Vdw:
 *
 * "None"
 * "LennardJones"
 * "LennardJonesSwitch"
 * "Buckingham"
 * "BuckinghamSwitch"
 * "CubicSplineTable"
 * "LinearTable"
 *
 *
 * Strings used for Geometry:
 *
 * "ParticleParticle"  (atom to atom)
 * "Water3Particle"    (spc or tip3p to atoms)
 * "Water3Water3"      (pairs of spc/tip3p)
 * "Water4Particle"    (tip4p to atoms)
 * "Water4Water4"      (pairs of tip4p)
 *
 *
 * Strings used for VF:
 *
 * "PotentialAndForce"
 * "Potential"
 * "Force"
 */
typedef struct nb_kernel_info
{
    nb_kernel_t *   kernelptr;
    const char *    kernelname;
    const char *    architecture;     /* e.g. "C", "SSE", "BlueGene", etc. */
    
    const char *    electrostatics;
    const char *    vdw;
    const char *    geometry;
    const char *    other;  /* Any extra info you want/need to select a kernel */
    const char *    vf;
}
nb_kernel_info_t;


void
nb_kernel_list_add_kernels(nb_kernel_info_t *   new_kernelinfo,
                           int                  new_size);

int
nb_kernel_list_hash_init(void);

/* Return a function pointer to the nonbonded kernel with settings according
 * to the text strings provided. GROMACS does not guarantee the existence of
 * accelerated kernels for any combination, so the return value can be NULL.
 * In that case, you can try a different/lower-level acceleration, and
 * eventually you need to prepare to fall back to generic kernels or change
 * your settings and try again.
 *
 * The names of the text strings are obviously meant to reflect settings in
 * GROMACS, but inside this routine they are merely used as four text
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
                          const char *        vdw,
                          const char *        geometry,
                          const char *        other,
                          const char *        vf);



#ifdef __cplusplus
}
#endif

#endif /* _nb_kernel_h_ */
