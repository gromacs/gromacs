/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \cond \internal \file scanner.l
 * \brief
 * Tokenizer for the selection language.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 * \endcond
 */
/*! \internal \file scanner.cpp
 * \brief
 * Generated (from scanner.l by Flex) tokenizer for the selection language.
 *
 * \ingroup module_selection
 */
%top{
#if !defined _gmx_sel_yyIN_HEADER
#include "gmxpre.h"
#endif

// Required before flex definitions, since it includes <stdint.h>.
// Otherwise, compilers not strictly C99 get macro redefinition errors,
// since flex defines INT32_MAX etc. in such cases.
#include "gromacs/utility/basedefinitions.h"
}
%{
#include "gromacs/utility/cstringutil.h"

#include "parser.h"
#include "scanner.h"
#include "scanner_internal.h"

// This macro makes the actions a bit shorter, since nearly every action needs
// this call.
#define ADD_TOKEN _gmx_sel_lexer_add_token(yylloc, yytext, yyleng, state)

// Set YY_BREAK to an empty value to avoid warnings (for the PGI compiler)
// when we have return statements followed by break. Instead, we add breaks
// manually.
#define YY_BREAK
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

%%

%{
    gmx_sel_lexer_t *state = yyget_extra(yyscanner);
    int              retval;
    /* Return a token if one is pending */
    retval = _gmx_sel_lexer_process_pending(yylval, yylloc, state);
    if (retval != 0)
    {
        return retval;
    }
    /* Handle the start conditions for 'of' matching */
    if (state->bMatchOf)
    {
        BEGIN(matchof);
        state->bMatchOf = false;
    }
    else if (state->bMatchBool)
    {
        BEGIN(matchbool);
        state->bMatchBool = false;
    }
    else if (state->bCmdStart)
    {
        BEGIN(cmdstart);
        state->bCmdStart = false;
    }
    else
    {
        BEGIN(0);
    }
%}

{COMMENT}       break;
{INTEGER}       { yylval->i   = strtol(yytext, NULL, 10);    ADD_TOKEN; return TOK_INT; }
{REAL}          { yylval->r   = strtod(yytext, NULL);        ADD_TOKEN; return TOK_REAL; }
{STRING}        { yylval->str = gmx_strndup(yytext+1, yyleng-2); ADD_TOKEN; return STR;  }

\\\n            { _gmx_sel_lexer_add_token(yylloc, " ", 1, state); break; }
";"|\n          {
                    if (yytext[0] == ';' || state->bInteractive)
                    {
                        rtrim(state->pselstr);
                        state->bCmdStart = true;
                        return CMD_SEP;
                    }
                    else
                    {
                        _gmx_sel_lexer_add_token(yylloc, " ", 1, state);
                    }
                    break;
                }

<cmdstart><<EOF>> { state->bCmdStart = true; yyterminate(); }
<<EOF>>         { state->bCmdStart = true; return CMD_SEP; }

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

{IDENTIFIER}    { return _gmx_sel_lexer_process_identifier(yylval, yylloc, yytext, yyleng, state); }

[[:space:]]+    { _gmx_sel_lexer_add_token(yylloc, " ", 1, state); break; }
[_[:alnum:]]+   { yylval->str = gmx_strndup(yytext, yyleng); ADD_TOKEN; return STR; }
.               { ADD_TOKEN; return yytext[0]; }
