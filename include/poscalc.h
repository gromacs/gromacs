/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
/*! \file
 * \brief API for structured and optimized calculation of positions.
 *
 * The functions in this header are used internally by the analysis library
 * to calculate positions.
 * They can also be used in user code, but in most cases there should be no
 * need. Instead, one should write an analysis tool such that it gets all
 * positions through selections.
 *
 * \internal
 *
 * The API is documented in more detail on a separate page:
 * \ref poscalcengine.
 */
#ifndef POSCALC_H
#define POSCALC_H

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \name Flags for position calculation.
 * \anchor poscalc_flags
 */
/*@{*/
/*! \brief
 * Use mass weighting.
 *
 * If this flag is set, the positions will be calculated using mass weighting,
 * i.e., one gets center-of-mass positions.
 * Without the flag, center-of-geometry positions are calculated.
 * Does not have any effect if the calculation type is \ref POS_ATOM.
 */
#define POS_MASS        1
/*! \brief
 * Calculate positions for the same atoms in residues/molecules.
 *
 * If this flag is set, the positions are always calculated using the same
 * atoms for each residue/molecule, even if the evaluation group contains only
 * some of the atoms for some frames.
 * The group passed to gmx_ana_poscalc_set_maxindex() is used to determine
 * the atoms to use for the calculation.
 *
 * Has no effect unless \ref POS_DYNAMIC is set or if the calculation type
 * is not \ref POS_RES of \ref POS_MOL.
 */
#define POS_COMPLMAX    2
/*! \brief
 * Calculate positions for whole residues/molecules.
 *
 * If this flag is set, the positions will be calculated for whole
 * residues/molecules, even if the group contains only some of the atoms in
 * the residue/molecule.
 *
 * Has no effect unless the calculation type is \ref POS_RES or \ref POS_MOL.
 */
#define POS_COMPLWHOLE  4
/*! \brief
 * Enable handling of changing calculation groups.
 *
 * Can be used for static calculations as well, but implies a small
 * performance penalty.
 */
#define POS_DYNAMIC     16
/*! \brief
 * Update \c gmx_ana_pos_t::m dynamically for an otherwise static
 * calculation.
 *
 * Has effect only if \ref POS_DYNAMIC is not set.
 */
#define POS_MASKONLY    32
/*! \brief
 * Calculate velocities of the positions.
 */
#define POS_VELOCITIES  64
/*! \brief
 * Calculate forces on the positions.
 */
#define POS_FORCES      128
/*@}*/

/** Specifies the type of positions to be calculated. */
typedef enum
{
    POS_ATOM,    /**< Copy atomic coordinates. */
    POS_RES,     /**< Calculate center for each residue. */
    POS_MOL,     /**< Calculate center for each molecule. */
    POS_ALL,     /**< Calculate center for the whole group. */
    POS_ALL_PBC  /**< Calculate center for the whole group with PBC. */
} e_poscalc_t;

/** Collection of \c gmx_ana_poscalc_t structures for the same topology. */
typedef struct gmx_ana_poscalc_coll_t gmx_ana_poscalc_coll_t;
/** Data structure for position calculation. */
typedef struct gmx_ana_poscalc_t gmx_ana_poscalc_t;

struct gmx_ana_index_t;
struct gmx_ana_pos_t;

/** Converts a string to parameters for gmx_ana_poscalc_create(). */
int
gmx_ana_poscalc_type_from_enum(const char *post, e_poscalc_t *type, int *flags);
/** Creates a list of strings for position enum parameter handling. */
const char **
gmx_ana_poscalc_create_type_enum(gmx_bool bAtom);

/** Creates a new position calculation collection object. */
int
gmx_ana_poscalc_coll_create(gmx_ana_poscalc_coll_t **pccp);
/** Sets the topology for a position calculation collection. */
void
gmx_ana_poscalc_coll_set_topology(gmx_ana_poscalc_coll_t *pcc, t_topology *top);
/** Frees memory allocated for a position calculation collection. */
void
gmx_ana_poscalc_coll_free(gmx_ana_poscalc_coll_t *pcc);
/** Prints information about calculations in a position calculation collection. */
void
gmx_ana_poscalc_coll_print_tree(FILE *fp, gmx_ana_poscalc_coll_t *pcc);

/** Creates a new position calculation. */
int
gmx_ana_poscalc_create(gmx_ana_poscalc_t **pcp, gmx_ana_poscalc_coll_t *pcc,
                       e_poscalc_t type, int flags);
/** Creates a new position calculation based on an enum value. */
int
gmx_ana_poscalc_create_enum(gmx_ana_poscalc_t **pcp, gmx_ana_poscalc_coll_t *pcc,
                            const char *post, int flags);
/** Sets the flags for position calculation. */
void
gmx_ana_poscalc_set_flags(gmx_ana_poscalc_t *pc, int flags);
/** Sets the maximum possible input index group for position calculation. */
void
gmx_ana_poscalc_set_maxindex(gmx_ana_poscalc_t *pc, struct gmx_ana_index_t *g);
/** Initializes positions for position calculation output. */
void
gmx_ana_poscalc_init_pos(gmx_ana_poscalc_t *pc, struct gmx_ana_pos_t *p);
/** Frees the memory allocated for position calculation. */
void
gmx_ana_poscalc_free(gmx_ana_poscalc_t *pc);
/** Returns TRUE if the position calculation requires topology information. */
gmx_bool
gmx_ana_poscalc_requires_top(gmx_ana_poscalc_t *pc);

/** Initializes evaluation for a position calculation collection. */
void
gmx_ana_poscalc_init_eval(gmx_ana_poscalc_coll_t *pcc);
/** Initializes a position calculation collection for a new frame. */
void
gmx_ana_poscalc_init_frame(gmx_ana_poscalc_coll_t *pcc);
/** Updates a single COM/COG structure for a frame. */
void
gmx_ana_poscalc_update(gmx_ana_poscalc_t *pc,
                       struct gmx_ana_pos_t *p, struct gmx_ana_index_t *g,
                       t_trxframe *fr, t_pbc *pbc);

#ifdef __cplusplus
}
#endif

#endif
