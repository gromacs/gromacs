/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2014,2015, by the GROMACS development team, led by
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
 * Parser/scanner interaction functions.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef SELECTION_SCANNER_H
#define SELECTION_SCANNER_H

#include <string>

#include <boost/exception_ptr.hpp>

#include "parser.h"

struct gmx_ana_indexgrps_t;
struct gmx_ana_selcollection_t;

#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void *yyscan_t;
#endif

/** Initializes the selection scanner. */
void
_gmx_sel_init_lexer(yyscan_t *scannerp, struct gmx_ana_selcollection_t *sc,
                    bool bInteractive, int maxnr, bool bGroups,
                    struct gmx_ana_indexgrps_t *grps);
/** Frees memory allocated for the selection scanner. */
void
_gmx_sel_free_lexer(yyscan_t scanner);
/** Stores an exception that is caught during parsing. */
void
_gmx_sel_lexer_set_exception(yyscan_t                    scanner,
                             const boost::exception_ptr &ex);
/** Rethrows and clears the stored exception if one is present. */
// TODO: The semantics is a bit confusing, need to be thought more,
// but easier to do as part of larger refactoring of the parsing.
void
_gmx_sel_lexer_rethrow_exception_if_occurred(yyscan_t scanner);

/** Returns true if the scanner is interactive. */
bool
_gmx_sel_is_lexer_interactive(yyscan_t scanner);
/** Returns the selection collection for the scanner. */
struct gmx_ana_selcollection_t *
_gmx_sel_lexer_selcollection(yyscan_t scanner);
/** Returns true if the external index groups for the scanner are set. */
bool
_gmx_sel_lexer_has_groups_set(yyscan_t scanner);
/** Returns the external index groups for the scanner. */
struct gmx_ana_indexgrps_t *
_gmx_sel_lexer_indexgrps(yyscan_t scanner);
/** Returns the number of selections after which the parser should stop. */
int
_gmx_sel_lexer_exp_selcount(yyscan_t scanner);

/** Returns a pretty string of the current selection.  */
const char *
_gmx_sel_lexer_pselstr(yyscan_t scanner);
/*! \brief
 * Sets the current parser context location.
 *
 * This location is set while Bison reductions are being processed, and
 * identifies the location of the current rule/reduction.
 */
void
_gmx_sel_lexer_set_current_location(yyscan_t                      scanner,
                                    const gmx::SelectionLocation &location);
/*! \brief
 * Returns the current parser context location.
 *
 * This returns the location last set with
 * _gmx_sel_lexer_set_current_location().
 */
const gmx::SelectionLocation &
_gmx_sel_lexer_get_current_location(yyscan_t scanner);
/*! \brief
 * Returns the selection text for the current parser context.
 *
 * This returns the selection text that corresponds to the position set last
 * with _gmx_sel_lexer_set_current_location().
 */
std::string
_gmx_sel_lexer_get_current_text(yyscan_t scanner);
/*! \brief
 * Returns the selection text at the given location.
 */
std::string
_gmx_sel_lexer_get_text(yyscan_t                      scanner,
                        const gmx::SelectionLocation &location);
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

/** The scanning function generated by Flex. */
#define YY_DECL int _gmx_sel_yylex(YYSTYPE *yylval, YYLTYPE *yylloc, yyscan_t yyscanner)
YY_DECL;

#endif
