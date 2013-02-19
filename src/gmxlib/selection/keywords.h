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
 * \brief Definitions of generic keyword evaluation structures.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef SELECTION_KEYWORDS_H
#define SELECTION_KEYWORDS_H

struct gmx_ana_selmethod_t;
struct t_selelem;
struct t_selexpr_param;

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

/** Sets the position type for position keyword evaluation. */
void
_gmx_selelem_set_kwpos_type(struct t_selelem *sel, const char *type);
/** Sets the flags for position keyword evaluation. */
void
_gmx_selelem_set_kwpos_flags(struct t_selelem *sel, int flags);

/** Does custom processing for parameters of the \c same selection method. */
int
_gmx_selelem_custom_init_same(struct gmx_ana_selmethod_t **method,
                              struct t_selexpr_param *params, void *scanner);

/** Initializes a selection element for evaluating a keyword in a given group. */
int
_gmx_sel_init_keyword_evaluator(struct t_selelem **sel,
                                struct gmx_ana_selmethod_t *method,
                                struct t_selexpr_param *param, void *scanner);

#endif
