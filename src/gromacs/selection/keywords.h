/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2012,2014,2015, by the GROMACS development team, led by
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
 * \brief Definitions of generic keyword evaluation structures.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_KEYWORDS_H
#define GMX_SELECTION_KEYWORDS_H

#include "parsetree.h"
#include "selelem.h"

struct gmx_ana_selmethod_t;

/** Selection method data for comparison expression evaluation. */
extern struct gmx_ana_selmethod_t sm_compare;

/** Selection method data for integer keyword evaluation. */
extern struct gmx_ana_selmethod_t sm_keyword_int;
/** Selection method data for real keyword evaluation. */
extern struct gmx_ana_selmethod_t sm_keyword_real;
/** Selection method data for string keyword evaluation. */
extern struct gmx_ana_selmethod_t sm_keyword_str;
/** Selection method data for position keyword evaluation. */
extern struct gmx_ana_selmethod_t sm_keyword_pos;

/** Prints information about a comparison expression. */
void
_gmx_selelem_print_compare_info(FILE *fp, void *data);

/*! \brief
 * Returns whether the selection element is a default position keyword.
 *
 * \param[in] sel   Selection element to query.
 * \returns   ``true`` if ``sel`` represents a position keyword evaluation that
 *     uses the default (implicit) position keyword.
 *
 * This method only works before the selection has been compiled.
 */
bool
_gmx_selelem_is_default_kwpos(const gmx::SelectionTreeElement &sel);
/** Sets the position type for position keyword evaluation. */
void
_gmx_selelem_set_kwpos_type(gmx::SelectionTreeElement *sel, const char *type);
/** Sets the flags for position keyword evaluation. */
void
_gmx_selelem_set_kwpos_flags(gmx::SelectionTreeElement *sel, int flags);

/** Sets the string match type for string keyword evaluation. */
void
_gmx_selelem_set_kwstr_match_type(const gmx::SelectionTreeElementPointer &sel,
                                  gmx::SelectionStringMatchType           matchType);

/** Does custom processing for parameters of the \c same selection method. */
void
_gmx_selelem_custom_init_same(struct gmx_ana_selmethod_t                    **method,
                              const gmx::SelectionParserParameterListPointer &params,
                              void                                           *scanner);

/*! \brief
 * Initializes a selection element for evaluating a keyword in a given group.
 *
 * \param[in]   method  Keyword selection method to evaluate.
 * \param[in]   child   The group/positions to evaluate \p method in.
 * \param[in]   scanner Scanner data structure.
 * \returns     Pointer to the created selection element.
 *
 * Creates a \ref SEL_EXPRESSION selection element that evaluates the keyword
 * method given by \p method in the group/positions given by \p child.
 *
 * \p child should be a selection tree that evaluates to \ref GROUP_VALUE or
 * \ref POS_VALUE.
 */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword_evaluator(struct gmx_ana_selmethod_t              *method,
                                const gmx::SelectionTreeElementPointer  &child,
                                void                                    *scanner);

#endif
