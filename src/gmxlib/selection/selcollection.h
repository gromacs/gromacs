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
/*! \internal \file
 * \brief Definition of \c gmx_ana_selcollection_t.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef SELECTION_COLLECTION_H
#define SELECTION_COLLECTION_H

#include <typedefs.h>

#include <indexutil.h>

/*! \internal
 * \brief
 * Information for a collection of selections.
 *
 * The functions to deal with the structure are defined in selection.h.
 * The structure is allocated with gmx_ana_selcollection_create() and
 * freed with gmx_ana_selcollection_free().
 * Some default values must then be set with
 * gmx_ana_selcollection_set_refpostype() and
 * gmx_ana_selcollection_set_outpostype().
 *
 * After setting the default values, one or more selections can be parsed
 * with gmx_ana_selcollection_parse_*().
 * At latest at this point, the topology must be set with
 * gmx_ana_selcollection_set_topology() unless
 * gmx_ana_selcollection_requires_top() returns FALSE.
 * Once all selections are parsed, they must be compiled all at once using
 * gmx_ana_selcollection_compile().
 * After these calls, gmx_ana_selcollection_get_count() and
 * gmx_ana_selcollection_get_selections() can be used
 * to get the compiled selections.
 * gmx_ana_selcollection_evaluate() can be used to update the selections for a
 * new frame.
 * gmx_ana_selcollection_evaluate_fin() can be called after all the frames have
 * been processed to restore the selection values back to the ones they were
 * after gmx_ana_selcollection_compile(), i.e., dynamic selections have the
 * maximal index group as their value.
 *
 * At any point, gmx_ana_selcollection_requires_top() can be called to see
 * whether the information provided so far requires loading the topology.
 * gmx_ana_selcollection_print_tree() can be used to print the internal
 * representation of the selections (mostly useful for debugging).
 */
struct gmx_ana_selcollection_t
{
    /** Default reference position type for selections. */
    const char                     *rpost;
    /** Default output position type for selections. */
    const char                     *spost;
    /** TRUE if \ref POS_MASKONLY should be used for output position evaluation. */
    gmx_bool                        bMaskOnly;
    /** TRUE if velocities should be evaluated for output positions. */
    gmx_bool                        bVelocities;
    /** TRUE if forces should be evaluated for output positions. */
    gmx_bool                        bForces;
    /** TRUE if debugging output should be printed during compilation. */
    gmx_bool                        bDebugCompile;

    /** Root of the selection element tree. */
    struct t_selelem              *root;
    /** Number of selections in \p sel. */
    int                            nr;
    /** Array of compiled selections. */
    struct gmx_ana_selection_t   **sel;
    /** Number of variables defined. */
    int                            nvars;
    /** Selection strings for variables. */
    char                         **varstrs;

    /** Topology for the collection. */
    t_topology                    *top;
    /** Index group that contains all the atoms. */
    struct gmx_ana_index_t         gall;
    /** Position calculation collection used for selection position evaluation. */
    struct gmx_ana_poscalc_coll_t *pcc;
    /** Memory pool used for selection evaluation. */
    struct gmx_sel_mempool_t      *mempool;
    /** Parser symbol table. */
    struct gmx_sel_symtab_t       *symtab;
};

/** Clears the symbol table in the collection */
void
_gmx_selcollection_clear_symtab(struct gmx_ana_selcollection_t *sc);

#endif
