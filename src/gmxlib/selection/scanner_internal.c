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
 * \brief Helper functions for the selection tokenizer.
 *
 * This file implements the functions in the headers scanner.h and
 * scanner_internal.h.
 */
/*! \internal file scanner_flex.h
 * \brief Generated (from scanner.l) header file by Flex.
 *
 * This file contains definitions of functions that are needed in
 * scanner_internal.c.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <typedefs.h>
#include <smalloc.h>
#include <string.h>

#include "string2.h"
#include "gmx_fatal.h"

#include <selmethod.h>

#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"
#include "symrec.h"

#include "parser.h"
#include "scanner.h"
#include "scanner_internal.h"

#define STRSTORE_ALLOCSTEP 1000

/* These are defined as macros in the generated scanner_flex.h.
 * We undefine them here to have them as variable names in the subroutines.
 * There are other ways of doing this, but this is probably the easiest. */
#undef yylval
#undef yytext
#undef yyleng

static gmx_bool
read_stdin_line(gmx_sel_lexer_t *state)
{
    char *ptr     = state->inputstr;
    int   max_len = state->nalloc_input;
    int   totlen  = 0;

    if (feof(stdin))
    {
        return FALSE;
    }
    if (state->bInteractive)
    {
        fprintf(stderr, "> ");
    }
    /* For some reason (at least on my Linux), fgets() doesn't return until
     * the user presses Ctrl-D _twice_ at the end of a non-empty line.
     * This can be a bit confusing for users, but there's not much we can
     * do, and the chances of a normal user noticing this are not very big. */
    while (fgets(ptr, max_len, stdin))
    {
        int len = strlen(ptr);

        totlen += len;
        if (len >= 2 && ptr[len - 1] == '\n' && ptr[len - 2] == '\\')
        {
            if (state->bInteractive)
            {
                fprintf(stderr, "... ");
            }
        }
        else if (len >= 1 && ptr[len - 1] == '\n')
        {
            break;
        }
        else if (len < max_len - 1)
        {
            if (state->bInteractive)
            {
                fprintf(stderr, "\n");
            }
            break;
        }
        ptr     += len;
        max_len -= len;
        if (max_len <= 2)
        {
            max_len             += state->nalloc_input;
            state->nalloc_input *= 2;
            len                  = ptr - state->inputstr;
            srenew(state->inputstr, state->nalloc_input);
            ptr = state->inputstr + len;
        }
    }
    if (ferror(stdin))
    {
        gmx_input("selection reading failed");
    }
    return totlen > 0;
}

int
_gmx_sel_yyblex(YYSTYPE *yylval, yyscan_t yyscanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(yyscanner);
    gmx_bool         bCmdStart;
    int              token;

    if (!state->bBuffer && !state->inputstr)
    {
        state->nalloc_input = 1024;
        snew(state->inputstr, state->nalloc_input);
        read_stdin_line(state);
        _gmx_sel_set_lex_input_str(yyscanner, state->inputstr);
    }
    bCmdStart = state->bCmdStart;
    token     = _gmx_sel_yylex(yylval, yyscanner);
    while (state->inputstr && token == 0 && read_stdin_line(state))
    {
        _gmx_sel_set_lex_input_str(yyscanner, state->inputstr);
        token = _gmx_sel_yylex(yylval, yyscanner);
    }
    if (token == 0 && !bCmdStart)
    {
        token = CMD_SEP;
        rtrim(state->pselstr);
    }
    state->bCmdStart = (token == CMD_SEP);
    return token;
}

static int
init_param_token(YYSTYPE *yylval, gmx_ana_selparam_t *param, gmx_bool bBoolNo)
{
    if (bBoolNo)
    {
        snew(yylval->str, strlen(param->name) + 3);
        yylval->str[0] = 'n';
        yylval->str[1] = 'o';
        strcpy(yylval->str+2, param->name);
    }
    else
    {
        yylval->str = param->name ? strdup(param->name) : NULL;
    }
    return PARAM;
}

static int
init_method_token(YYSTYPE *yylval, gmx_ana_selmethod_t *method, gmx_bool bPosMod,
                  gmx_sel_lexer_t *state)
{
    /* If the previous token was not KEYWORD_POS, return EMPTY_POSMOD
     * before the actual method to work around a limitation in Bison. */
    if (!bPosMod && method->type != POS_VALUE)
    {
        state->nextmethod = method;
        return EMPTY_POSMOD;
    }
    yylval->meth = method;
    if (!(method->flags & SMETH_MODIFIER) && method->nparams == 0)
    {
        /* Keyword */
        switch (method->type)
        {
            case INT_VALUE:   return KEYWORD_NUMERIC;
            case REAL_VALUE:  return KEYWORD_NUMERIC;
            case STR_VALUE:   return KEYWORD_STR;
            case GROUP_VALUE: return KEYWORD_GROUP;
            default:          return INVALID;
        }
    }
    else
    {
        /* Method with parameters or a modifier */
        if (method->flags & SMETH_MODIFIER)
        {
            /* Remove all methods from the stack */
            state->msp = -1;
            if (method->param[1].name == NULL)
            {
                state->nextparam = &method->param[1];
            }
        }
        else
        {
            if (method->param[0].name == NULL)
            {
                state->nextparam = &method->param[0];
            }
        }
        ++state->msp;
        if (state->msp >= state->mstack_alloc)
        {
            state->mstack_alloc += 10;
            srenew(state->mstack, state->mstack_alloc);
        }
        state->mstack[state->msp] = method;
        if (method->flags & SMETH_MODIFIER)
        {
            return MODIFIER;
        }
        switch (method->type)
        {
            case INT_VALUE:   return METHOD_NUMERIC;
            case REAL_VALUE:  return METHOD_NUMERIC;
            case POS_VALUE:   return METHOD_POS;
            case GROUP_VALUE: return METHOD_GROUP;
            default:
                --state->msp;
                return INVALID;
        }
    }
    return INVALID; /* Should not be reached */
}

int
_gmx_sel_lexer_process_pending(YYSTYPE *yylval, gmx_sel_lexer_t *state)
{
    if (state->nextparam)
    {
        gmx_ana_selparam_t     *param   = state->nextparam;
        gmx_bool                bBoolNo = state->bBoolNo;

        if (state->neom > 0)
        {
            --state->neom;
            return END_OF_METHOD;
        }
        state->nextparam = NULL;
        state->bBoolNo   = FALSE;
        _gmx_sel_lexer_add_token(param->name, -1, state);
        return init_param_token(yylval, param, bBoolNo);
    }
    if (state->prev_pos_kw > 0)
    {
        --state->prev_pos_kw;
    }
    if (state->nextmethod)
    {
        gmx_ana_selmethod_t *method = state->nextmethod;

        state->nextmethod = NULL;
        return init_method_token(yylval, method, TRUE, state);
    }
    return 0;
}

int
_gmx_sel_lexer_process_identifier(YYSTYPE *yylval, char *yytext, size_t yyleng,
                                  gmx_sel_lexer_t *state)
{
    gmx_sel_symrec_t *symbol;
    e_symbol_t        symtype;

    /* Check if the identifier matches with a parameter name */
    if (state->msp >= 0)
    {
        gmx_ana_selparam_t     *param   = NULL;
        gmx_bool                bBoolNo = FALSE;
        int                     sp      = state->msp;
        while (!param && sp >= 0)
        {
            int             i;
            for (i = 0; i < state->mstack[sp]->nparams; ++i)
            {
                /* Skip NULL parameters and too long parameters */
                if (state->mstack[sp]->param[i].name == NULL
                    || strlen(state->mstack[sp]->param[i].name) > yyleng)
                {
                    continue;
                }
                if (!strncmp(state->mstack[sp]->param[i].name, yytext, yyleng))
                {
                    param = &state->mstack[sp]->param[i];
                    break;
                }
                /* Check separately for a 'no' prefix on gmx_boolean parameters */
                if (state->mstack[sp]->param[i].val.type == NO_VALUE
                    && yyleng > 2 && yytext[0] == 'n' && yytext[1] == 'o'
                    && !strncmp(state->mstack[sp]->param[i].name, yytext+2, yyleng-2))
                {
                    param   = &state->mstack[sp]->param[i];
                    bBoolNo = TRUE;
                    break;
                }
            }
            if (!param)
            {
                --sp;
            }
        }
        if (param)
        {
            if (param->val.type == NO_VALUE && !bBoolNo)
            {
                state->bMatchBool = TRUE;
            }
            if (sp < state->msp)
            {
                state->neom      = state->msp - sp - 1;
                state->nextparam = param;
                state->bBoolNo   = bBoolNo;
                return END_OF_METHOD;
            }
            _gmx_sel_lexer_add_token(param->name, -1, state);
            return init_param_token(yylval, param, bBoolNo);
        }
    }

    /* Check if the identifier matches with a symbol */
    symbol = _gmx_sel_find_symbol_len(state->sc->symtab, yytext, yyleng, FALSE);
    /* If there is no match, return the token as a string */
    if (!symbol)
    {
        yylval->str = gmx_strndup(yytext, yyleng);
        _gmx_sel_lexer_add_token(yytext, yyleng, state);
        return IDENTIFIER;
    }
    _gmx_sel_lexer_add_token(_gmx_sel_sym_name(symbol), -1, state);
    symtype = _gmx_sel_sym_type(symbol);
    /* Reserved symbols should have been caught earlier */
    if (symtype == SYMBOL_RESERVED)
    {
        return INVALID;
    }
    /* For variable symbols, return the type of the variable value */
    if (symtype == SYMBOL_VARIABLE)
    {
        t_selelem *var;

        var = _gmx_sel_sym_value_var(symbol);
        /* Return simple tokens for constant variables */
        if (var->type == SEL_CONST)
        {
            switch (var->v.type)
            {
                case INT_VALUE:
                    yylval->i = var->v.u.i[0];
                    return TOK_INT;
                case REAL_VALUE:
                    yylval->r = var->v.u.r[0];
                    return TOK_REAL;
                case POS_VALUE:
                    break;
                default:
                    return INVALID;
            }
        }
        yylval->sel = var;
        switch (var->v.type)
        {
            case INT_VALUE:   return VARIABLE_NUMERIC;
            case REAL_VALUE:  return VARIABLE_NUMERIC;
            case POS_VALUE:   return VARIABLE_POS;
            case GROUP_VALUE: return VARIABLE_GROUP;
            default:          return INVALID;
        }
        return INVALID;
    }
    /* For method symbols, return the correct type */
    if (symtype == SYMBOL_METHOD)
    {
        gmx_ana_selmethod_t *method;

        method = _gmx_sel_sym_value_method(symbol);
        return init_method_token(yylval, method, state->prev_pos_kw > 0, state);
    }
    /* For position symbols, we need to return KEYWORD_POS, but we also need
     * some additional handling. */
    if (symtype == SYMBOL_POS)
    {
        state->bMatchOf    = TRUE;
        yylval->str        = _gmx_sel_sym_name(symbol);
        state->prev_pos_kw = 2;
        return KEYWORD_POS;
    }
    /* Should not be reached */
    return INVALID;
}

void
_gmx_sel_lexer_add_token(const char *str, int len, gmx_sel_lexer_t *state)
{
    /* Do nothing if the string is empty, or if it is a space and there is
     * no other text yet, or if there already is a space. */
    if (!str || len == 0 || strlen(str) == 0
        || (str[0] == ' ' && str[1] == 0
            && (state->pslen == 0 || state->pselstr[state->pslen - 1] == ' ')))
    {
        return;
    }
    if (len < 0)
    {
        len = strlen(str);
    }
    /* Allocate more memory if necessary */
    if (state->nalloc_psel - state->pslen < len)
    {
        int incr = STRSTORE_ALLOCSTEP < len ? len : STRSTORE_ALLOCSTEP;
        state->nalloc_psel += incr;
        srenew(state->pselstr, state->nalloc_psel);
    }
    /* Append the token to the stored string */
    strncpy(state->pselstr + state->pslen, str, len);
    state->pslen                += len;
    state->pselstr[state->pslen] = 0;
}

int
_gmx_sel_init_lexer(yyscan_t *scannerp, struct gmx_ana_selcollection_t *sc,
                    gmx_bool bInteractive, int maxnr,
                    struct gmx_ana_indexgrps_t *grps)
{
    gmx_sel_lexer_t *state;
    int              rc;

    rc = _gmx_sel_yylex_init(scannerp);
    if (rc != 0)
    {
        return rc;
    }

    snew(state, 1);
    state->sc        = sc;
    state->grps      = grps;
    state->nexpsel   = (maxnr > 0 ? sc->nr + maxnr : -1);

    state->bInteractive = bInteractive;
    state->nalloc_input = 0;
    state->inputstr     = NULL;

    snew(state->pselstr, STRSTORE_ALLOCSTEP);
    state->pselstr[0]   = 0;
    state->pslen        = 0;
    state->nalloc_psel  = STRSTORE_ALLOCSTEP;

    snew(state->mstack, 20);
    state->mstack_alloc = 20;
    state->msp          = -1;
    state->neom         = 0;
    state->nextparam    = NULL;
    state->nextmethod   = NULL;
    state->prev_pos_kw  = 0;
    state->bBoolNo      = FALSE;
    state->bMatchOf     = FALSE;
    state->bMatchBool   = FALSE;
    state->bCmdStart    = TRUE;
    state->bBuffer      = FALSE;

    _gmx_sel_yyset_extra(state, *scannerp);
    return 0;
}

void
_gmx_sel_free_lexer(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);

    sfree(state->inputstr);
    sfree(state->pselstr);
    sfree(state->mstack);
    if (state->bBuffer)
    {
        _gmx_sel_yy_delete_buffer(state->buffer, scanner);
    }
    sfree(state);
    _gmx_sel_yylex_destroy(scanner);
}

gmx_bool
_gmx_sel_is_lexer_interactive(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);
    return state->bInteractive;
}

struct gmx_ana_selcollection_t *
_gmx_sel_lexer_selcollection(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);
    return state->sc;
}

struct gmx_ana_indexgrps_t *
_gmx_sel_lexer_indexgrps(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);
    return state->grps;
}

int
_gmx_sel_lexer_exp_selcount(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);
    return state->nexpsel;
}

const char *
_gmx_sel_lexer_pselstr(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);
    return state->pselstr;
}

void
_gmx_sel_lexer_clear_pselstr(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);
    state->pselstr[0] = 0;
    state->pslen      = 0;
}

void
_gmx_sel_lexer_clear_method_stack(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);

    state->msp = -1;
}

void
_gmx_sel_finish_method(yyscan_t scanner)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);

    if (state->msp >= 0)
    {
        --state->msp;
    }
}

void
_gmx_sel_set_lex_input_file(yyscan_t scanner, FILE *fp)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);

    state->bBuffer = TRUE;
    state->buffer  = _gmx_sel_yy_create_buffer(fp, YY_BUF_SIZE, scanner);
    _gmx_sel_yy_switch_to_buffer(state->buffer, scanner);
}

void
_gmx_sel_set_lex_input_str(yyscan_t scanner, const char *str)
{
    gmx_sel_lexer_t *state = _gmx_sel_yyget_extra(scanner);

    if (state->bBuffer)
    {
        _gmx_sel_yy_delete_buffer(state->buffer, scanner);
    }
    state->bBuffer = TRUE;
    state->buffer  = _gmx_sel_yy_scan_string((yyconst char*) str, scanner);
}
