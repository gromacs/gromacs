/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief Definitions of generic keyword evaluation structures.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_KEYWORDS_H
#define GMX_SELECTION_KEYWORDS_H

#include "parsetree.h"
#include "selelem.h"

struct gmx_ana_selmethod_t;
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
_gmx_selelem_set_kwpos_type(gmx::SelectionTreeElement *sel, const char *type);
/** Sets the flags for position keyword evaluation. */
void
_gmx_selelem_set_kwpos_flags(gmx::SelectionTreeElement *sel, int flags);

/** Does custom processing for parameters of the \c same selection method. */
int
_gmx_selelem_custom_init_same(struct gmx_ana_selmethod_t **method,
                              const gmx::SelectionParserParameterListPointer &params,
                              void *scanner);

/** Initializes a selection element for evaluating a keyword in a given group. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword_evaluator(struct gmx_ana_selmethod_t *method,
                                const gmx::SelectionParserParameterList &params,
                                void *scanner);

#endif
