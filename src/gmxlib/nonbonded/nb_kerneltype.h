/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _NB_KERNEL_H_
#define _NB_KERNEL_H_

#include <types/simple.h>

/** \file
 *  \brief Function interfaces and nonbonded kernel meta-data.
 *
 *  \internal
 *
 *  This file contains the definitions of the call sequences and 
 *  parameter documetation for the nonbonded kernels, both the 
 *  optimized level1/level2 kernels in Fortran or C, and the more 
 *  general level0 normal and free energy kernels. 
 *  It also defines an information structure used to record
 *  data about nonbonded kernels, their type of call sequence, 
 *  and the number of flops used in the outer and inner loops.
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif




/* Temporary structure to be able to pass 2 GB data through the
 * work data pointer argument.
 */
typedef struct 
{
	real      gb_epsilon_solvent;  /* Epsilon for solvent */
	real      epsilon_r;           /* Epsilon for inner dielectric */
	real *    gpol;                /* Polarization energy group */
} gmx_gbdata_t;
	
	
/** Interface to level1 optimized nonbonded kernels.
 * 
 *  \internal
 *
 *  To simplify the jungle of nonbonded function calls in earlier
 *  versions of Gromacs, we have changed most nonbonded kernels
 *  to use the same call sequence, even if some of the parameters
 *  are ununsed in a particular kernel.
 *
 *  This function call sequence is used by all level1 kernels. When this is
 *  written it is actually used by the level2 kernels too, but that might
 *  change in the future - don't assume they are identical.
 *  All arguments are passed as pointers to make to interface identical across 
 *  the Fortran and C implementations. For the same reason (Fortran & assembly) 
 *  we cannot use structures.
 *
 *  The nonbonded kernels should never be called directly - use
 *  the setup and wrapper routines in gmx_nonbonded.h instead.
 *
 *  \param nri        Number of neighborlists/i-particles/outer-particles.
 *                    The outer loop of the nonbonded kernel will iterate
 *                    over these indices. This number must always
 *                    be provided,
 *  \param iinr       Atom number (starting on 0) for each index 0..(nri-1)
 *                    in the outer loop. These are the 'owner'/'parent' 
 *                    particles of each neighborlist. It is quite common
 *                    that a particle occurs several times with different
 *                    shift indices. This is the field 'ilist' in 
 *                    the neighborlist. This array must always be provided.
 *  \param jindex     This array specifies the index  of the first j neighbor
 *                    in the array jjnr[] which belongs to the corresponding
 *                    i atom. The j neighbors for all lists are merged into
 *                    one continous array, so we also get the first index for
 *                    the next list in the next element. The size of this list
 *                    is nri+1, and the last position contains the number of
 *                    entries in the j list. Confused? It is actually pretty
 *                    straightforward; for index <tt>i</tt> the i atom is 
 *                    <tt>iinr[i]</tt>, and the atom indices of the neighbors 
 *                    to interact with are in positions <tt>jindex[i]</tt> 
 *                    through <tt>(jindex[i+1]-1)</tt> in <tt>jjnr[]</tt>.
 *                    This corresponds to ilist_jindex in the neighborlist.
 *                    This array must always be provided.
 *  \param jjnr       Array with the j particle neighbors of list i in positions
 *                    <tt>jindex[i]</tt> through <tt>(jindex[i+1]-1)</tt>. This
 *                    field is called jlist in the neighborlist structure.
 *                    This array must always be provided.
 *  \param shift      Array with shift vector index corresponding to each 
 *                    iinr neighborlist. In most codes, periodic boundary 
 *                    conditions are applied by checking pairwise distances 
 *                    of particles, and correcting differences of more than 
 *                    half a box. This is costly (latencies) and unnecessary 
 *                    since we already decided which periodic copy to use during
 *                    neighborsearching. We solve this by introducing 
 *                    \a shiftvectors and a \a shiftindex. The shiftvectors is 
 *                    normally an array of 3*3*3 vectors, corresponding to 
 *                    displacements +- one box unit in one or more directions.
 *                    The central box (no displacement) is always in the middle,
 *                    index 13 in the case of 3*3*3 shift vectors.
 *                    Imagine a particle interacting with two particles in the
 *                    current box, and three particles in the box to the left. 
 *                    We represent this with two ilist entries, one with shift 
 *                    index=13 (central box) and one with shift index=12. 
 *                    all neighbors in each sublist will have the same shift,  
 *                    so in the nonbonded interaction, we only have to add the 
 *                    shift vector to the outer (i) particle coordinates before
 *                    starting the loop over neighbors. This extracts periodic 
 *                    boundary condition calculations from the inner loops, 
 *                    increasing performance significantly. This array must
 *                    always be provided, but if you do not use any periodic
 *                    boundary conditions you could set all elements to zero
 *                    and provide a shiftvec array of length three, with all
 *                    elements there too being zero.
 *  \param shiftvec   The shift vectors for each index. Since the code needs to
 *                    interface with Fortran or assembly, we have to use 
 *                    pointers to float/double instead of the gmx_rvec_t type.
 *                    The x/y/z shift distances for shift index idx are in
 *                    <tt>shiftvec[3*idx]</tt>, <tt>shiftvec[3*idx+1]</tt>, 
 *                    and <tt>shiftvec[3*idx+2]</tt>. This is fully binary 
 *                    compatible with an array of gmx_rvec_t, so you can simply
 *                    use a typecast when calling. The length of the array is
 *                    three times the number of shift indices. This array
 *                    must be provided, but see the comment for shift.
 *  \param fshift     Forces from 'other periodic boxes' for each shift index.
 *                    The shift concept does not affect the forces on individual
 *                    atoms, but it will have an effect on the virial. To 
 *                    account for this, we add the total force on the i particle 
 *                    in each outer loop to the current shift index in this 
 *                    array, which has exactly the same dimensions as the 
 *                    shiftvec array above. The Gromacs manual describes how 
 *                    this is used to calculate the effect in the virial.
 *                    Dimension is identical to shiftvec.
 *                    For kernels that do not calculate forces this can be NULL.
 *  \param gid        Energy group index for each i neighborlist index; this
 *                    corresponds to ilist_groupindex in the neighborlist and is
 *                    used to decompose energies into groupwise contributions.
 *                    If an i particle belongs to group A, and has interactions 
 *                    both with groups A,B, and C, you will get three different 
 *                    i list entries with energy group indices corresponding to 
 *                    A-A, A-B, and A-C neighbors, in addition to possible
 *                    copies due to different shift indices. The energies will 
 *                    be added to this index (i.e. index=gid[i]) in the vc and 
 *                    vnb arrays. This array must always be provided.
 *                    If you don't want to decompose energies all groups should 
 *                    be set to index 0 (length or the array nri).
 *  \param pos        Coordinates of atoms, starting on 0 as all arrays. To
 *                    interface with fortran and assembly this is a simple list
 *                    of floats, with x/y/z of atom n in positions
 *                    <tt>pos[3*n]</tt>, <tt>pos[3*n+1]</tt>, and
 *                    <tt>pos[3*n+2]</tt>. This is binary compatible with,
 *                    and can be typecast from, an array of gmx_rvec_t.
 *                    This array must always be provided.
 *  \param faction    Forces of atoms. To interface with fortran and assembly 
 *                    this is a simple list of floats, with x/y/z of atom n in 
 *                    positions <tt>pos[3*n]</tt>, <tt>pos[3*n+1]</tt>, and
 *                    <tt>pos[3*n+2]</tt>. This is binary compatible with,
 *                    and can be typecast from, an array of gmx_rvec_t.
 *                    For kernels that do not calculate forces this can be NULL.
 *  \param charge     Array with the charge of atom n in position n.
 *                    If you are calling a kernel without coulomb interaction
 *                    this parameter can safely be set to NULL.
 *  \param facel      Factor to multiply all electrostatic interactions with;
 *                    this is essentially 1/(4*pi*epsilon). For kernels without
 *                    coulomb interaction the value will not be used. In the
 *                    fortran version (pass-by-reference) it can be NULL in
 *                    that case.
 *  \param krf        The 'k' constant for reaction-field electrostatic
 *                    interactions. Can be zero/NULL for non-RF interactions.
 *  \param crf        The 'c' constant for reaction-field electrostatic
 *                    interactions. Can be zero/NULL for non-RF interactions.
 *  \param vc         List of Coulomb energies, the length is equal to the
 *                    number of energy group combinations. As described for
 *                    gid above, the coulomb energy of each outer/i-particle
 *                    list will be added to vc[gid[i]]. Can be NULL for
 *                    kernels without electrostatic interaction.
 *  \param type       Array with the index of the Van der Waals type for each
 *                    atom, starting on 0 . This is used to look up the VdW 
 *                    parameters. When no VdW interactions are calculated 
 *                    (i.e. only Coulomb), this can be set to NULL. 
 *  \param ntype      Number of VdW types, used to look up the VdW parameters. 
 *                    This parameter can be zero (NULL for pass-by-reference)
 *                    if the kernel does not include VdW interactions.
 *  \param vdwparam   Array with Van der Waals parameters. The contents and 
 *                    size depends on the number of parameters required for
 *                    the selected VdW interaction (3 for B-ham, 2 for LJ, etc).
 *                    There are a total of ntype*ntype pair combinations, and 
 *                    the size of this array is the number of such pairs times
 *                    the number of parameters.
 *                    With 'nparam' parameters, and 'ntype' types in total,
 *                    the first parameter for interactions between atoms of VdW 
 *                    type 'ti' and VdW type 'tj' is located in
 *                    <tt>vdwparam[nparam*(ntype*ti+tj)]</tt>, and the 
 *                    additional parameters follow directly in the next 
 *                    position(s). Note that the array corresponds to a 
 *                    symmetric matrix - you should get the same parameters if
 *                    you swap ti and tj. This parameter can be NULL 
 *                    if the kernel does not include VdW interactions.
 *  \param vvdw       List of Van der Waals energies, the length is equal to 
 *                    the number of energy group combinations. As described 
 *                    for gid above, the coulomb energy of each outer list
 *                    will be added to vc[gid[i]]. Can be NULL for
 *                    kernels without electrostatic interaction.
 *  \param tabscale   Distance between table points in the vftab table.
 *                    This is always the same for Coulomb and VdW when both
 *                    are tabulated. If the kernel does not include any
 *                    tabulated interactions it can be zero, or NULL when
 *                    passed-by-reference.
 *  \param vftab      Table data for Coulomb and/or VdW interactions. 
 *                    Each 'point' in the table consists of the four 
 *                    floating-point numbers Y,F,G,H as described in the
 *                    Gromacs manual. If the kernel only tabulates the 
 *                    Coulomb interaction this is followed by the next
 *                    point. If both Coulomb and Van der Waals interactions 
 *                    are tabulated it is instead followed by another
 *                    four numbers for the dispersion table at the same
 *                    point, and then finally the repulsion table, i.e.
 *                    a total of 12 numbers per table point. If only 
 *                    Van der Waals interactions are tabulated the 
 *                    Dispersion and Repulsion table (in that order) use
 *                    a total of 8 floating-point numbers per point.
 *                    This array only needs to be provided for kernels with
 *                    tabulated interactions - otherwise it can be NULL.
 *  \param invsqrta   Array with the inverse square root of the born radius
 *                    of each atom. This can safely be NULL when the kernel
 *                    does not calculate Generalized Born interactions.
 *  \param dvda       Array where the derivative of the potential with respect
 *                    to the born radius for each atom will be written. This
 *                    is necessary to calculate the effect of born radius
 *                    derivatives on forces. However, just as for the force
 *                    arrays it must be provided even when not calculating
 *                    forces, since an implementation might not provide a
 *                    non-force version of the routine. When the kernel does
 *                    not calculate Generalized Born interactions it can 
 *                    safely be set to NULL, though.
 *  \param gbtabscale Distance between (scaled) table points for tabulated 
 *                    Generalized Born interactions. If the kernel does not 
 *                    calculate Generalized Born interactions it can be zero, 
 *                    or NULL when passed-by-reference. Note that this is 
 *                    usually different from the standard (non-GB) table scale.
 *  \param gbtab      Table data for Generalized Born Coulomb interactions.
 *                    Since these interactions contain parameters that 
 *                    enter in a non-trivial way I have introduced a trick
 *                    where the tabulated interaction is a function of
 *                    the distance scaled with the square roots of the two
 *                    born radii. Apart from this it contains similar 
 *                    data (Y,F,G,H) as the normal Coulomb/VdW tables.
 *                    See comments in the table code for details. This
 *                    can safely be set to NULL if you are not using
 *                    Generalized Born interactions.
 *  \param nthreads   Number of threads calling the kernel concurrently.
 *                    This value is only used to optimize the chunk size
 *                    of the neighborlists processed by each thread, so
 *                    it won't cause incorrect results if you get it wrong.
 *  \param count      Pointer to a common counter to synchronize 
 *                    multithreaded execution. This is normally the counter
 *                    in the neighborlist. It must be set to zero before
 *                    calling the kernel to get correct results. Note that
 *                    this is necessary even for a single thread if Gromacs
 *                    is compiled with thread support (meaning it is
 *                    a mandatory parameter).
 *  \param mtx        Pointer to a mutex protecting the counter, masked as 
 *                    a pointer-to-void since thread support is optional. 
 *                    However, if Gromacs is compiled with threads you \a must 
 *                    provide a mutex - it would hurt performance to check
 *                    it for each iteration in the nonbonded kernel. The
 *                    mutex is normally created automatically by the 
 *                    neighborlist initialization routine gmx_nlist_init().
 *                    (Yet another reason not to call kernels directly).
 *  \param outeriter  The number of iterations performed in the outer loop
 *                    of the kernel will be written to this integer. If
 *                    multiple threads are calling the loop this will be
 *                    the part of nri handled by this thread. 
 *                    This value is only used to get accurate flop accounting.
 *                    Since we cannot check for NULL pointers in Fortran it must
 *                    always be provided, but you can discard
 *                    the result if you don't need it.
 *  \param inneriter  The number of iterations performed in the inner loop
 *                    of the kernel will be written to this integer; compare
 *                    with nouter. Since we cannot check for NULL pointers in
 *                    Fortran it must always be provided, but you can discard
 *                    the result if you don't need it.
 *  \param work       Double precision work array - to be used for optimization.
 *
 *  \warning          There is \a very little error control inside
 *                    the nonbonded kernels in order to maximize performance.
 *                    It cannot be overemphasized that you should never call
 *                    them directly, but rely on the higher-level routines.
 */
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



#ifdef __cplusplus
}
#endif


#endif /* _NB_KERNEL_H_ */
