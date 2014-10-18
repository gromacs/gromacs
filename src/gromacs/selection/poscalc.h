/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * \brief
 * API for structured and optimized calculation of positions.
 *
 * This header declares an API for calculating positions in an automated way,
 * for internal use by the selection engine.  This is useful in particular with
 * dynamic selections, because the same COM/COG positions may be needed in
 * several contexts.  The API makes it possible to optimize the evaluation such
 * that any heavy calculation is only done once, and the results just copied if
 * needed more than once.  The functions also provide a convenient interface
 * for keeping the whole \c gmx_ana_pos_t structure up-to-date.
 *
 * The API is documented in more detail in gmx::PositionCalculationCollection.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_POSCALC_H
#define GMX_SELECTION_POSCALC_H

#include <cstdio>

#include "gromacs/utility/classhelpers.h"

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

/** Data structure for position calculation. */
struct gmx_ana_poscalc_t;

struct gmx_ana_index_t;
struct gmx_ana_pos_t;
struct t_pbc;
struct t_topology;
struct t_trxframe;

namespace gmx
{

/*! \internal
 * \brief
 * Collection of \c gmx_ana_poscalc_t structures for the same topology.
 *
 * Calculations within one collection share the same topology, and they are
 * optimized.  Calculations in different collections do not interact.
 * The topology for a collection can be set with setTopology().
 * This needs to be done before calling gmx_ana_poscalc_set_maxindex() for
 * any calculation in the collection, unless that calculation does not
 * require topology information.
 *
 * A new calculation is created with createCalculation().
 * If flags need to be adjusted later, gmx_ana_poscalc_set_flags() can be
 * used.
 * After the flags are final, the largest possible index group for which the
 * positions are needed has to be set with gmx_ana_poscalc_set_maxindex().
 * setTopology() should have been called before this function is called.
 * After the above calls, gmx_ana_poscalc_init_pos() can be used to initialize
 * output to a \c gmx_ana_pos_t structure.  Several different structures can be
 * initialized for the same calculation; the only requirement is that the
 * structure passed later to gmx_ana_poscalc_update() has been initialized
 * properly.
 * The memory allocated for a calculation can be freed with
 * gmx_ana_poscalc_free().
 *
 * The position evaluation is simple: initFrame() should be
 * called once for each frame, and gmx_ana_poscalc_update() can then be called
 * for each calculation that is needed for that frame.
 *
 * It is also possible to initialize the calculations based on a type provided
 * as a string.
 * The possible strings are listed in \ref typeEnumValues, and the string can
 * be converted to the parameters for createCalculation() using typeFromEnum().
 * createCalculationFromEnum() is also provided for convenience.
 *
 * \ingroup module_selection
 */
class PositionCalculationCollection
{
    public:
        /*! \brief
         * Array of strings acceptable for position calculation type enum.
         *
         * This array contains the acceptable values for typeFromEnum() and
         * createCalculationFromEnum().
         * The array contains a NULL pointer after the last item to indicate
         * the end of the list.
         */
        static const char * const typeEnumValues[];

        /*! \brief
         * Converts a string to parameters for createCalculationFromEnum().
         *
         * \param[in]     post  String (typically an enum argument).
         *     Allowed values: 'atom', 'res_com', 'res_cog', 'mol_com', 'mol_cog',
         *     or one of the last four prepended by 'whole_', 'part_', or 'dyn_'.
         * \param[out]    type  \c e_poscalc_t corresponding to \p post.
         * \param[in,out] flags Flags corresponding to \p post.
         *     On input, the flags should contain the default flags.
         *     On exit, the flags \ref POS_MASS, \ref POS_COMPLMAX and
         *     \ref POS_COMPLWHOLE have been set according to \p post
         *     (the completion flags are left at the default values if no
         *     completion prefix is given).
         * \throws  InternalError  if post is not recognized.
         *
         * \attention
         * Checking is not complete, and other values than those listed above
         * may be accepted for \p post, but the results are undefined.
         *
         * \see typeEnumValues
         */
        static void typeFromEnum(const char *post, e_poscalc_t *type, int *flags);

        /*! \brief
         * Creates a new position calculation collection object.
         *
         * \throws  std::bad_alloc if out of memory.
         */
        PositionCalculationCollection();
        /*! \brief
         * Destroys a position calculation collection and its calculations.
         *
         * Any calculations in the collection are also freed, even if
         * references to them are left.
         */
        ~PositionCalculationCollection();

        /*! \brief
         * Sets the topology used for the calculations.
         *
         * \param[in]     top   Topology data structure.
         *
         * This function should be called to set the topology before using
         * gmx_ana_poscalc_set_maxindex() for any calculation that requires
         * topology information.
         *
         * Does not throw.
         */
        void setTopology(t_topology *top);
        /*! \brief
         * Prints information about calculations.
         *
         * \param[in] fp    File handle to receive the output.
         *
         * The output is very technical, making this function mainly useful for
         * debugging purposes.
         *
         * Does not throw.
         */
        void printTree(FILE *fp) const;

        /*! \brief
         * Creates a new position calculation.
         *
         * \param[in]  type  Type of calculation.
         * \param[in]  flags Flags for setting calculation options
         *   (see \ref poscalc_flags "documentation of the flags").
         *
         * Does not throw currently, but may throw std::bad_alloc in the
         * future.
         */
        gmx_ana_poscalc_t *createCalculation(e_poscalc_t type, int flags);
        /*! \brief
         * Creates a new position calculation based on an enum value.
         *
         * \param[in]  post  One of the strings acceptable for
         *      typeFromEnum().
         * \param[in]  flags Flags for setting calculation options
         *      (see \ref poscalc_flags "documentation of the flags").
         * \throws     InternalError  if post is not recognized.
         *
         * This is a convenience wrapper for createCalculation().
         * \p flags sets the default calculation options if not overridden by
         * \p post; see typeFromEnum().
         *
         * May also throw std::bad_alloc in the future.
         *
         * \see createCalculation(), typeFromEnum()
         */
        gmx_ana_poscalc_t *createCalculationFromEnum(const char *post, int flags);

        /*! \brief
         * Computes the highest atom index required to evaluate this collection.
         *
         * Does not throw.
         */
        int getHighestRequiredAtomIndex() const;

        /*! \brief
         * Initializes evaluation for a position calculation collection.
         *
         * This function does some final initialization of the data structures
         * in the collection to prepare them for evaluation.
         * After this function has been called, it is no longer possible to add
         * new calculations to the collection.
         *
         * Multiple calls to the function are ignored.
         *
         * Does not throw currently, but may throw std::bad_alloc in the
         * future.
         */
        void initEvaluation();
        /*! \brief
         * Initializes a position calculation collection for a new frame.
         *
         * Clears the evaluation flag for all calculations.
         * Should be called for each frame before calling
         * gmx_ana_poscalc_update().
         *
         * This function calls initEvaluation() automatically if it has not
         * been called earlier.
         *
         * Does not throw, but may throw if initEvaluation() is changed to
         * throw.
         */
        void initFrame();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        /*! \brief
         * Needed to access the implementation class from the C code.
         */
        friend struct ::gmx_ana_poscalc_t;
};

} // namespace gmx

/** Sets the flags for position calculation. */
void
gmx_ana_poscalc_set_flags(gmx_ana_poscalc_t *pc, int flags);
/** Sets the maximum possible input index group for position calculation. */
void
gmx_ana_poscalc_set_maxindex(gmx_ana_poscalc_t *pc, gmx_ana_index_t *g);
/** Initializes positions for position calculation output. */
void
gmx_ana_poscalc_init_pos(gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p);
/** Frees the memory allocated for position calculation. */
void
gmx_ana_poscalc_free(gmx_ana_poscalc_t *pc);
/** Returns true if the position calculation requires topology information. */
bool
gmx_ana_poscalc_requires_top(gmx_ana_poscalc_t *pc);

/** Updates a single COM/COG structure for a frame. */
void
gmx_ana_poscalc_update(gmx_ana_poscalc_t *pc,
                       gmx_ana_pos_t *p, gmx_ana_index_t *g,
                       t_trxframe *fr, t_pbc *pbc);

#endif
