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
/*! \cond \internal \file scanner.l
 * \brief
 * Tokenizer for the selection language.
 * \endcond
 */
/*! \internal \file scanner.c
 * \brief
 * Generated (from scanner.l by Flex) tokenizer for the selection language.
 */
%{
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string2.h>

#include "parser.h"
#include "scanner.h"
#include "scanner_internal.h"

/* This macro is here to make the actions a bit shorter, since nearly every
 * action needs this call. */
#define ADD_TOKEN _gmx_sel_lexer_add_token(yytext, yyleng, state)

%}

INTEGER    [[:digit:]]+
DSEQ       ([[:digit:]]+)
FRAC       (([[:digit:]]*"."{DSEQ})|{DSEQ}".")
EXP        ([eE][+-]?{DSEQ})
REAL       (({FRAC}{EXP}?)|({DSEQ}{EXP}))
STRING     (\"([^\"\\\n]|(\\\"))*\")
IDENTIFIER ([[:alpha:]][_[:alnum:]]*)
CMPOP      (([<>]=?)|([!=]=))
COMMENT    (#.*)

%option nodefault
%option noyywrap
%option reentrant
%option prefix="_gmx_sel_yy"
%option header-file="scanner_flex.h"
%option nounistd
%option never-interactive

%s matchof
%s matchbool
%s cmdstart
%x help

%%

%{
    gmx_sel_lexer_t *state = yyget_extra(yyscanner);
    int              retval;
    /* Return a token if one is pending */
    retval = _gmx_sel_lexer_process_pending(yylval, state);
    if (retval != 0)
    {
        return retval;
    }
    /* Handle the start conditions for 'of' matching */
    if (state->bMatchOf)
    {
        BEGIN(matchof);
        state->bMatchOf = FALSE;
    }
    else if (state->bMatchBool)
    {
        BEGIN(matchbool);
        state->bMatchBool = FALSE;
    }
    else if (state->bCmdStart)
    {
        BEGIN(cmdstart);
    }
    else if (YYSTATE != help)
    {
        BEGIN(0);
    }
%}

{COMMENT}
{INTEGER}       { yylval->i   = strtol(yytext, NULL, 10);    ADD_TOKEN; return TOK_INT; }
{REAL}          { yylval->r   = strtod(yytext, NULL);        ADD_TOKEN; return TOK_REAL; }
{STRING}        { yylval->str = gmx_strndup(yytext+1, yyleng-2); ADD_TOKEN; return STR;  }

\\\n            { _gmx_sel_lexer_add_token(" ", 1, state); }
";"|\n          {
                    if (yytext[0] == ';' || state->bInteractive)
                    {
                        rtrim(state->pselstr);
                        return CMD_SEP;
                    }
                    else
                    {
                        _gmx_sel_lexer_add_token(" ", 1, state);
                    }
                }

<cmdstart>help  { BEGIN(help); return HELP; }
<help>{
[[:blank:]]+
{IDENTIFIER}    { yylval->str = gmx_strndup(yytext, yyleng); return HELP_TOPIC; }
";"|\n          { return CMD_SEP; }
.               { return INVALID; }
}

<matchbool>{
yes|on          { ADD_TOKEN; yylval->i = 1; return TOK_INT; }
no|off          { ADD_TOKEN; yylval->i = 0; return TOK_INT; }
}
group           { ADD_TOKEN; return GROUP; }
to              { ADD_TOKEN; return TO; }
<matchof>of     { ADD_TOKEN; BEGIN(0); return OF; }
and|"&&"        { ADD_TOKEN; return AND; }
or|"||"         { ADD_TOKEN; return OR; }
xor             { ADD_TOKEN; return XOR; }
not|"!"         { ADD_TOKEN; return NOT; }
{CMPOP}         { yylval->str = gmx_strndup(yytext, yyleng); ADD_TOKEN; return CMP_OP; }

{IDENTIFIER}    { return _gmx_sel_lexer_process_identifier(yylval, yytext, yyleng, state); }

[[:blank:]]+    { _gmx_sel_lexer_add_token(" ", 1, state); }
[_[:alnum:]]+   { yylval->str = gmx_strndup(yytext, yyleng); ADD_TOKEN; return STR; }
.               { ADD_TOKEN; return yytext[0]; }
