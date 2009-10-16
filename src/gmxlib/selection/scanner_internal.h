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
 * \brief Internal header file used by the selection tokenizer.
 */
#ifndef SELECTION_SCANNER_INTERNAL_H
#define SELECTION_SCANNER_INTERNAL_H

#include "parser.h"

/* These need to be defined before including scanner_flex.h, because it
 * uses YY_EXTRA_TYPE. But we also need to include it before defining
 * gmx_sel_lexer_t; hence the forward declaration. */
struct gmx_sel_lexer_t;
#define YY_EXTRA_TYPE struct gmx_sel_lexer_t *

/* We cannot include scanner_flex.h from the scanner itself, because it
 * seems to break everything. */
/* And we need to define YY_NO_UNISTD_H here as well, otherwise unistd.h
 * gets included in other files than scanner.c... */
#ifndef FLEX_SCANNER
#define YY_NO_UNISTD_H
#include "scanner_flex.h"
#endif

/*! \brief
 * Internal data structure for the selection tokenizer state.
 */
typedef struct gmx_sel_lexer_t
{
    struct gmx_ana_selcollection_t  *sc;
    bool                             bPrompt;
    const char                      *prompt;

    char                            *pselstr;
    int                              pslen;
    int                              nalloc_psel;

    struct gmx_ana_selmethod_t     **mstack;
    int                              msp;
    int                              mstack_alloc;
    int                              neom;
    struct gmx_ana_selparam_t       *nextparam;
    bool                             bBoolNo;

    bool                             bMatchOf;
    bool                             bMatchBool;
    bool                             bCmdStart;

    bool                             bBuffer;
    YY_BUFFER_STATE                  buffer;
} gmx_sel_lexer_t;

/* Because Flex defines yylval, yytext, and yyleng as macros,
 * and this file is included from scanner.l,
 * we cannot have them here as parameter names... */
/** Internal function for cases where several tokens need to be returned for a
 * single parameter. */
int
_gmx_sel_lexer_process_next_param(YYSTYPE *, gmx_sel_lexer_t *state);
/** Internal function that processes identifier tokens. */
int
_gmx_sel_lexer_process_identifier(YYSTYPE *, char *, int,
                                  gmx_sel_lexer_t *state);
/** Internal helper function that prints a prompt if appropriate. */
void
_gmx_sel_lexer_prompt_print(gmx_sel_lexer_t *state);
/** Internal helper function that updates the prompt after a newline. */
void
_gmx_sel_lexer_prompt_newline(bool bContinue, gmx_sel_lexer_t *state);
/** Internal function to add a token to the pretty-printed selection text. */
void
_gmx_sel_lexer_add_token(const char *str, int len, gmx_sel_lexer_t *state);

#endif
