/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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

/*! \libinternal \file
 * \brief Assemble atomic positions of a (small) subset of atoms and distribute to all nodes.
 *
 *  This file contains functions to assemble the positions of a subset of the
 *  atoms and to do operations on it like determining the center of mass, or
 *  doing translations and rotations. These functions are useful when
 *  a subset of the positions needs to be compared to some set of reference
 *  positions, as e.g. done for essential dynamics.
 *
 * \inlibraryapi
 */

#include <stdio.h>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Select local atoms of a group.
 *
 * Selects the indices of local atoms of a group and stores them in anrs_loc[0..nr_loc].
 * If you need the positions of the group's atoms on all nodes, provide a coll_ind[0..nr]
 * array and pass it on to communicate_group_positions. Thus the collective array
 * will always have the same atom order (ascending indices).
 *
 *  \param[in]     ga2la      Global to local atom index conversion data.
 *  \param[in]     nr         The total number of atoms that the group contains.
 *  \param[in]     anrs       The global atom number of the group's atoms.
 *  \param[out]    nr_loc     The number of group atoms present on the local node.
 *  \param[out]    anrs_loc   The local atom numbers of the group.
 *  \param[in,out] nalloc_loc Local allocation size of anrs_loc array.
 *  \param[out]    coll_ind   If not NULL this array must be of size nr. It stores
 *                            for each local atom where it belongs in the global
 *                            (collective) array such that it can be gmx_summed
 *                            in the communicate_group_positions routine.
 */

extern void dd_make_local_group_indices(gmx_ga2la_t ga2la,
                                        const int nr, int anrs[], int *nr_loc,
                                        int *anrs_loc[], int *nalloc_loc,
                                        int coll_ind[]);


/*! \brief Assemble local positions into a collective array present on all nodes.
 *
 * Communicate the positions of the group's atoms such that every node has all of
 * them. Unless running on huge number of cores, this is not a big performance impact
 * as long as the collective subset [0..nr] is kept small. The atom indices are
 * retrieved from anrs_loc[0..nr_loc]. If you call the routine for the serial case,
 * provide an array coll_ind[i] = i for i in 1..nr.
 *
 * If shifts != NULL, the PBC representation of each atom is chosen such that a
 * continuous trajectory results. Therefore, if the group is whole at the start
 * of the simulation, it will always stay whole.
 * If shifts = NULL, the group positions are not made whole again, but assembled
 * and distributed to all nodes. The variables marked "optional" are not used in
 * that case.
 *
 * \param[in]     cr           Pointer to MPI communication data.
 * \param[out]    xcoll        Collective array of positions, identical on all nodes
 *                             after this routine has been called.
 * \param[in,out] shifts       Collective array of shifts for xcoll, needed to make
 *                             the group whole. This array remembers the shifts
 *                             since the start of the simulation (where the group
 *                             is whole) and must therefore not be changed outside
 *                             of this routine! If NULL, the group will not be made
 *                             whole and the optional variables are ignored.
 * \param[out]    extra_shifts Extra shifts since last time step, only needed as
 *                             buffer variable [0..nr] (optional).
 * \param[in]     bNS          Neighbor searching / domain re-decomposition has been
 *                             performed at the begin of this time step such that
 *                             the shifts have changed and need to be updated
 *                             (optional).
 * \param[in]     x_loc        Pointer to the local atom positions this node has.
 * \param[in]     nr           Total number of atoms in the group.
 * \param[in]     nr_loc       Number of group atoms on the local node.
 * \param[in]     anrs_loc     Array of the local atom indices.
 * \param[in]     coll_ind     This array of size nr stores for each local atom where
 *                             it belongs in the collective array so that the local
 *                             contributions can be gmx_summed. It is provided by
 *                             dd_make_local_group_indices.
 * \param[in,out] xcoll_old    Positions from the last time step, used to make the
 *                             group whole (optional).
 * \param[in]     box          Simulation box matrix, needed to shift xcoll such that
 *                             the group becomes whole (optional).
 */
extern void communicate_group_positions(t_commrec *cr, rvec *xcoll, ivec *shifts,
                                        ivec *extra_shifts, const gmx_bool bNS,
                                        rvec *x_loc, const int nr, const int nr_loc,
                                        int *anrs_loc, int *coll_ind, rvec *xcoll_old,
                                        matrix box);

/*! \brief Calculates the center of the positions x locally.
 *
 * Calculates the center of mass (if masses are given in the weight array) or
 * the geometrical center (if NULL is passed as weight).
 *
 * \param[in]   x            Positions.
 * \param[in]   weight       Can be NULL or an array of weights. If masses are
 *                           given as weights, the COM is calculated.
 * \param[in]   nr           Number of positions and weights if present.
 * \param[out]  center       The (weighted) center of the positions.
 *
 */
extern void get_center(rvec x[], real weight[], const int nr, rvec center);


/*! \brief Calculates the sum of the positions x locally.
 *
 * Calculates the (weighted) sum of position vectors and returns the sum of
 * weights, which is needed when local contributions shall be summed to a
 * global weighted center.
 *
 * \param[in]   x            Array of positions.
 * \param[in]   weight       Can be NULL or an array of weights.
 * \param[in]   nr           Number of positions and weights if present.
 * \param[out]  dsumvec      The (weighted) sum of the positions.
 * \return Sum of weights.
 *
 */
extern double get_sum_of_positions(rvec x[], real weight[], const int nr, dvec dsumvec);


/*! \brief Calculates the global center of all local arrays x.
 *
 * Get the center from local positions [0..nr_loc], this involves communication.
 * Not that the positions must already have the correct PBC representation. Use
 * this routine if no collective coordinates are assembled from which the center
 * could be calculated without communication.
 *
 * \param[in]   cr           Pointer to MPI communication data.
 * \param[in]   x_loc        Array of local positions [0..nr_loc].
 * \param[in]   weight_loc   Array of local weights, these are the masses if the
 *                           center of mass is to be calculated.
 * \param[in]   nr_loc       The number of positions on the local node.
 * \param[in]   nr_group     The number of positions in the whole group. Since
 *                           this is known anyway, we do not need to communicate
 *                           and sum nr_loc if we pass it over.
 * \param[out]  center       The (weighted) center of all x_loc from all the
 *                           nodes.
 */
extern void get_center_comm(t_commrec *cr, rvec x_loc[], real weight_loc[],
                            int nr_loc, int nr_group, rvec center);


/*! \brief Translate positions.
 *
 * Add a translation vector to the positions x.
 *
 * \param[in,out] x          Array of positions.
 * \param[in]     nr         Number of entries in the position array.
 * \param[in]     transvec   Translation vector to be added to all positions.
 *
 */
extern void translate_x(rvec x[], const int nr, const rvec transvec);


/*! \brief Rotate positions.
 *
 * Rotate the positions with the rotation matrix.
 *
 * \param[in,out] x          Array of positions.
 * \param[in]     nr         Number of entries in the position array.
 * \param[in]     rmat       Rotation matrix to operate on all positions.
 *
 */
extern void rotate_x(rvec x[], const int nr, matrix rmat);

#ifdef __cplusplus
}
#endif
