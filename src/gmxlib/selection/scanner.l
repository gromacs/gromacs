/*
 * $Id$
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
/*! \internal \file scanner.l
 * \brief
 * Tokenizer for the selection language.
 *
 * \todo
 * With a recent version of Flex, it should be straightforward to make the
 * scanner reentrant (it is the only component of the selection parser that
 * still uses global variables).
 */
/*! \internal \file scanner.c
 * \brief
 * Generated (from scanner.l by Flex) tokenizer for the selection language.
 */
%{
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include <smalloc.h>
#include <string2.h>

#include <selmethod.h>

#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"
#include "symrec.h"

#include "parser.h"
#include "scanner.h"

#define DEFAULT_PROMPT     ">"
#define CONTINUE_PROMPT    "..."
#define STRSTORE_ALLOCSTEP 1000

struct gmx_sel_lexer_t
{
    struct gmx_ana_selcollection_t *sc;
    bool                  bPrompt;
    char                 *prompt;

    char                **strstore;
    int                   slen;
    int                   nalloc_str;

    gmx_ana_selmethod_t **mstack;
    int                   msp;
    int                   mstack_alloc;
    int                   neom;
    gmx_ana_selparam_t   *nextparam;

    bool                  bMatchOf;

    bool                  bBuffer;
    YY_BUFFER_STATE       buffer;
};

static gmx_sel_lexer_t *gstate;

static int
process_next_param(YYSTYPE *yylval, gmx_sel_lexer_t *state);
static int
process_identifier(YYSTYPE *yylval, char *yytext, int yyleng,
                   gmx_sel_lexer_t *state);

#define YY_INPUT(buf,result,max_size) \
    { \
        gmx_sel_lexer_t *state = gstate; \
        int c = '*', n; \
        if (state->bPrompt) \
        { \
            fprintf(stderr, "%s ", state->prompt); \
            state->bPrompt = FALSE; \
        } \
        for (n = 0; n < max_size && \
                    (c = getc(yyin)) != EOF && c != '\n'; ++n) \
        { \
            buf[n] = (char)c; \
        } \
        if (c == '\n') \
        { \
            buf[n++] = (char)c; \
            if (state->prompt) \
            { \
                state->prompt  = DEFAULT_PROMPT; \
                state->bPrompt = TRUE; \
            } \
        } \
        if (c == EOF && ferror(yyin)) \
        { \
            YY_FATAL_ERROR("input in flex scanner failed"); \
        } \
        result = n; \
        if (state->strstore) \
        { \
            while (n > state->nalloc_str - state->slen) \
            { \
                state->nalloc_str += STRSTORE_ALLOCSTEP; \
                srenew(*state->strstore, state->nalloc_str); \
            } \
            strncpy((*state->strstore)+state->slen, buf, n); \
            state->slen += n; \
            if (state->nalloc_str > 0) \
            { \
                (*state->strstore)[state->slen] = 0; \
            } \
        } \
    }
/*! \cond */
%}

INTEGER    [[:digit:]]+
DSEQ       ([[:digit:]]+)
FRAC       (([[:digit:]]*"."{DSEQ})|{DSEQ}".")
EXP        ([eE][+-]?{DSEQ})
REAL       ("-"?(({FRAC}{EXP}?)|({DSEQ}{EXP})))
STRING     (\"([^\"\\\n]|(\\\"))*\")
IDENTIFIER ([[:alpha:]][_[:alnum:]]*)
CMPOP      (([<>]=?)|([!=]=))
COMMENT    (#.*)

%option nodefault
%option noyywrap

%s matchof

%%

%{
    /* Return END_OF_METHOD/PARAM_* immediately if necessary */
    if (state->nextparam)
    {
        return process_next_param(yylval, state);
    }
    /* Handle the start conditions for 'of' matching */
    if (state->bMatchOf)
    {
        BEGIN(matchof);
        state->bMatchOf = FALSE;
    }
    else
    {
        BEGIN(0);
    }
%}

{COMMENT}
{INTEGER}       { yylval->i   = atoi(yytext);                return INT; }
{REAL}          { yylval->r   = atof(yytext);                return REAL; }
{STRING}        { yylval->str = strndup(yytext+1, yyleng-2); return STR; }

\\\n            { if (state->prompt) state->prompt = CONTINUE_PROMPT; }
";"|\n          { if (yytext[0] == ';' || state->prompt) return CMD_SEP; }
group           { return GROUP; }
to              { return TO; }
<matchof>of     { BEGIN(0); return OF; }
and|"&&"        { return AND; }
or|"||"         { return OR; }
xor             { return XOR; }
not|"!"         { return NOT; }
{CMPOP}         { yylval->str = strndup(yytext, yyleng); return CMP_OP; }

{IDENTIFIER}    { return process_identifier(yylval, yytext, yyleng, state); }

[[:blank:]]+
[_[:alnum:]]+   { yylval->str = strndup(yytext, yyleng);     return STR; }
.               { return yytext[0]; }

%%
/*! \endcond */

static int
process_next_param(YYSTYPE *yylval, gmx_sel_lexer_t *state)
{
    gmx_ana_selparam_t *param = state->nextparam;

    if (state->neom > 0)
    {
        --state->neom;
        return END_OF_METHOD;
    }
    state->nextparam = NULL;
    yylval->str = param->name;
    switch (param->val.type)
    {
        case NO_VALUE:    return PARAM_BOOL;
        case INT_VALUE:   return PARAM_INT;
        case REAL_VALUE:  return PARAM_REAL;
        case STR_VALUE:   return PARAM_STR;
        case POS_VALUE:   return PARAM_POS;
        case GROUP_VALUE: return PARAM_GROUP;
    }
    return INVALID; /* Should not be reached */
}

static int
process_identifier(YYSTYPE *yylval, char *yytext, int yyleng,
                   gmx_sel_lexer_t *state)
{
    gmx_sel_symrec_t *symbol;
    e_symbol_t        symtype;

    /* Check if the identifier matches with a parameter name */
    if (state->msp >= 0)
    {
        gmx_ana_selparam_t *param = NULL;
        int                 sp = state->msp;
        while (!param && sp >= 0)
        {
            int             i;
            for (i = 0; i < state->mstack[sp]->nparams; ++i)
            {
                if (state->mstack[sp]->param[i].name == NULL)
                {
                    continue;
                }
                if (!strncmp(state->mstack[sp]->param[i].name, yytext, yyleng))
                {
                    param = &state->mstack[sp]->param[i];
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
            if (sp < state->msp)
            {
                state->neom = state->msp - sp - 1;
                state->nextparam = param;
                return END_OF_METHOD;
            }
            yylval->str = param->name;
            switch (param->val.type)
            {
                case NO_VALUE:    return PARAM_BOOL;
                case INT_VALUE:   return PARAM_INT;
                case REAL_VALUE:  return PARAM_REAL;
                case STR_VALUE:   return PARAM_STR;
                case POS_VALUE:   return PARAM_POS;
                case GROUP_VALUE: return PARAM_GROUP;
            }
            return INVALID; /* Should not be reached */
        }
    }

    /* Check if the identifier matches with a symbol */
    symbol = _gmx_sel_find_symbol_len(state->sc->symtab, yytext, yyleng, FALSE);
    /* If there is no match, return the token as a string */
    if (!symbol)
    {
        yylval->str = strndup(yytext, yyleng);
        return IDENTIFIER;
    }
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
                    return INT;
                case REAL_VALUE:
                    yylval->r = var->v.u.r[0];
                    return REAL;
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
        yylval->meth = method;
        if (!(method->flags & SMETH_MODIFIER) && method->nparams == 0)
        {
            /* Keyword */
            switch (method->type)
            {
                case INT_VALUE:   return KEYWORD_INT;
                case REAL_VALUE:  return KEYWORD_REAL;
                case STR_VALUE:   return KEYWORD_STR;
                case GROUP_VALUE: return KEYWORD_GROUP;
                default:          return INVALID;
            }
        } else {
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
    }
    /* For position symbols, we need to return KEYWORD_POS, but we also need
     * some additional handling. */
    if (symtype == SYMBOL_POS)
    {
        state->bMatchOf = TRUE;
        yylval->str = _gmx_sel_sym_name(symbol);
        return KEYWORD_POS;
    }
    /* Should not be reached */
    return INVALID;
}

void
_gmx_sel_init_lexer(gmx_sel_lexer_t **statep, struct gmx_ana_selcollection_t *sc,
                    bool bInteractive)
{
    gmx_sel_lexer_t *state;

    snew(state, 1);
    state->sc        = sc;
    state->bPrompt   = bInteractive;
    state->prompt    = bInteractive ? DEFAULT_PROMPT : NULL;

    state->strstore  = &sc->selstr;
    if (sc->selstr)
    {
        state->slen       = strlen(sc->selstr);
        state->nalloc_str = state->slen + 1;
    }
    else
    {
        state->slen       = 0;
        state->nalloc_str = 0;
    }

    snew(state->mstack, 20);
    state->mstack_alloc = 20;
    state->msp          = -1;
    state->neom         = 0;
    state->nextparam    = NULL;
    state->bMatchOf     = FALSE;
    state->bBuffer      = FALSE;

    *statep = state;
    gstate  = state;
}

void
_gmx_sel_free_lexer(gmx_sel_lexer_t *state)
{
    sfree(state->mstack);
    if (state->bBuffer)
    {
        yy_delete_buffer(state->buffer);
    }
    sfree(state);
    gstate  = NULL;
}

bool
_gmx_sel_is_lexer_interactive(gmx_sel_lexer_t *state)
{
    return state->bPrompt;
}

struct gmx_ana_selcollection_t *
_gmx_sel_lexer_selcollection(gmx_sel_lexer_t *state)
{
    return state->sc;
}

void
_gmx_sel_finish_method(gmx_sel_lexer_t *state)
{
    if (state->msp >= 0)
    {
        --state->msp;
    }
}

void
_gmx_sel_set_lex_input_file(gmx_sel_lexer_t *state, FILE *fp)
{
    state->bBuffer = TRUE;
    state->buffer  = yy_create_buffer(fp, YY_BUF_SIZE);
    yy_switch_to_buffer(state->buffer);
}

void
_gmx_sel_set_lex_input_str(gmx_sel_lexer_t *state, const char *str)
{

    state->bBuffer = TRUE;
    state->buffer  = yy_scan_string(str);
    /* Append the buffer to the string store as YY_INPUT is not called */
    if (state->strstore)
    {
        int len, slen;

        if (*state->strstore)
        {
            slen = strlen(*state->strstore);
        }
        else
        {
            slen = 0;
        }
        len = strlen(str);
        snew(*state->strstore, len + slen + 1);
        strcpy((*state->strstore)+slen, str);
    }
}
