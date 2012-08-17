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
 * \brief Helper functions for the selection parser.
 *
 * This header is includes only from parser.cpp (generated from parser.y), and
 * it includes functions and macros used internally by the parser.
 * They are in a separate file to make then easier to edit (no need to
 * regenerate the parser), and to keep parser.y as simple as possible.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_PARSER_INTERNAL_H
#define GMX_SELECTION_PARSER_INTERNAL_H

#include <exception>

#include "parsetree.h"
#include "selelem.h"

#include "scanner.h"

//! Helper method to reorder a list of parameter values.
static t_selexpr_value *
process_value_list(t_selexpr_value *values)
{
    t_selexpr_value *val, *pval, *nval;

    /* Count values (if needed) and reverse list */
    pval = NULL;
    val  = values;
    while (val)
    {
        nval = val->next;
        val->next = pval;
        pval = val;
        val = nval;
    }
    values = pval;

    return values;
}

//! Error handler needed by Bison.
static void
yyerror(yyscan_t scanner, char const *s)
{
    _gmx_selparser_error(scanner, "%s", s);
}

/*! \name Exception handling macros for actions
 *
 * These macros should be used at the beginning and end of each semantic action
 * that may throw an exception. For robustness, it's best to wrap all actions
 * that call functions declared outside parser.y should be wrapped.
 * These macros take care to catch any exceptions, store the exception (or
 * handle it and allow the parser to continue), and terminate the parser
 * cleanly if necessary.
 * The code calling the parser should use
 * _gmx_sel_lexer_rethrow_exception_if_occurred() to rethrow any exceptions.
 * \{
 */
//! Starts an action that may throw exceptions.
#define BEGIN_ACTION \
    try {
//! Finishes an action that may throw exceptions.
#define END_ACTION \
    } \
    catch(const std::exception &ex) \
    { \
        if (_gmx_selparser_handle_exception(scanner, ex)) \
            YYERROR; \
        else \
            YYABORT; \
    }
//!\}

/*! \brief
 * Retrieves a smart pointer from a semantic value.
 *
 * \param[in] src  Semantic value to get the pointer from.
 * \returns   Retrieved smart pointer.
 *
 * There should be no statements that may throw exceptions in actions before
 * this function has been called for all semantic values that have a smart
 * pointer stored.  Together with set, this function abstracts away exception
 * safety issues that arise from the use of a plain pointer for storing the
 * semantic values.
 *
 * Does not throw.
 */
template <typename PointerType> static
PointerType get(PointerType *src)
{
    PointerType result;
    if (src != NULL)
    {
        result.swap(*src);
        delete src;
    }
    return move(result);
}
/*! \brief
 * Sets a smart pointer to a semantic value.
 *
 * \param[out] dest  Semantic value to set (typically $$).
 * \param[in]  value Pointer to put into the semantic value.
 * \throws     std::bad_alloc if out of memory.
 *
 * This should be the last statement before ::END_ACTION, except for a
 * possible ::CHECK_SEL.
 */
template <typename PointerType> static
void set(PointerType *&dest, PointerType value)
{
    dest = new PointerType(move(value));
}
/*! \brief
 * Checks that a valid tree was set.
 *
 * Should be called after set_sel() if it was used to set a value where NULL
 * pointer indicates an error.
 *
 * \todo
 * Get rid of this macro.  It should now be possible to handle all errors using
 * exceptions.
 */
#define CHECK_SEL(sel) \
    if (!*(sel)) { \
        delete sel; \
        YYERROR; \
    }

#endif
