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
 * \brief
 * Parser/scanner interaction functions.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef SELECTION_SCANNER_H
#define SELECTION_SCANNER_H

#include "parser.h"

struct gmx_ana_indexgrps_t;
struct gmx_ana_selcollection_t;

#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void *yyscan_t;
#endif

/** Initializes the selection scanner. */
int
_gmx_sel_init_lexer(yyscan_t *scannerp, struct gmx_ana_selcollection_t *sc,
                    gmx_bool bInteractive, int maxnr,
                    struct gmx_ana_indexgrps_t *grps);
/** Frees memory allocated for the selection scanner. */
void
_gmx_sel_free_lexer(yyscan_t scanner);

/** Returns TRUE if the scanner is interactive. */
gmx_bool
_gmx_sel_is_lexer_interactive(yyscan_t scanner);
/** Returns the selection collection for the scanner. */
struct gmx_ana_selcollection_t *
_gmx_sel_lexer_selcollection(yyscan_t scanner);
/** Returns the external index groups for the scanner. */
struct gmx_ana_indexgrps_t *
_gmx_sel_lexer_indexgrps(yyscan_t scanner);
/** Returns the number of selections after which the parser should stop. */
int
_gmx_sel_lexer_exp_selcount(yyscan_t scanner);

/** Returns a pretty string of the current selection.  */
const char *
_gmx_sel_lexer_pselstr(yyscan_t scanner);
/** Clears the current selection string.  */
void
_gmx_sel_lexer_clear_pselstr(yyscan_t scanner);
/** Clears the method stack in the scanner in error situations. */
void
_gmx_sel_lexer_clear_method_stack(yyscan_t scanner);
/** Notifies the scanner that a complete method expression has been parsed. */
void
_gmx_sel_finish_method(yyscan_t scanner);
/** Initializes the scanner to scan a file. */
void
_gmx_sel_set_lex_input_file(yyscan_t scanner, FILE *fp);
/** Initializes the scanner to scan a string. */
void
_gmx_sel_set_lex_input_str(yyscan_t scanner, const char *str);

/** A wrapper for the actual scanner, used by the Bison parser. */
int
_gmx_sel_yyblex(YYSTYPE *yylval, yyscan_t yyscanner);

#endif
