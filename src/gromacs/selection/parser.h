/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY__GMX_SEL_YY_PARSER_H_INCLUDED
#define YY__GMX_SEL_YY_PARSER_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
#    define YYDEBUG 1
#endif
#if YYDEBUG
extern int _gmx_sel_yydebug;
#endif
/* "%code requires" blocks.  */
#line 1 "parser.y" /* yacc.c:1909  */

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#line 76 "parser.y" /* yacc.c:1909  */

#include "parsetree.h"
#include "selelem.h"
#include "gromacs/utility/real.h"

#define YYLTYPE ::gmx::SelectionLocation

#line 87 "parser.h" /* yacc.c:1909  */

/* Token type.  */
#ifndef YYTOKENTYPE
#    define YYTOKENTYPE
enum yytokentype
{
    INVALID          = 258,
    TOK_INT          = 259,
    TOK_REAL         = 260,
    STR              = 261,
    IDENTIFIER       = 262,
    CMD_SEP          = 263,
    GROUP            = 264,
    TO               = 265,
    VARIABLE_NUMERIC = 266,
    VARIABLE_GROUP   = 267,
    VARIABLE_POS     = 268,
    KEYWORD_NUMERIC  = 269,
    KEYWORD_STR      = 270,
    KEYWORD_POS      = 271,
    KEYWORD_GROUP    = 272,
    METHOD_NUMERIC   = 273,
    METHOD_GROUP     = 274,
    METHOD_POS       = 275,
    MODIFIER         = 276,
    EMPTY_POSMOD     = 277,
    PARAM            = 278,
    END_OF_METHOD    = 279,
    OF               = 280,
    CMP_OP           = 281,
    PARAM_REDUCT     = 282,
    OR               = 283,
    XOR              = 284,
    AND              = 285,
    NOT              = 286,
    UNARY_NEG        = 287,
    NUM_REDUCT       = 288
};
#endif

/* Value type.  */
#if !defined YYSTYPE && !defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#    line 83 "parser.y" /* yacc.c:1909  */

    int                         i;
    real                        r;
    char*                       str;
    struct gmx_ana_selmethod_t* meth;

    gmx::SelectionStringMatchType smt;

    gmx::SelectionTreeElementPointer*         sel;
    gmx::SelectionParserValue*                val;
    gmx::SelectionParserValueListPointer*     vlist;
    gmx::SelectionParserParameter*            param;
    gmx::SelectionParserParameterListPointer* plist;

#    line 148 "parser.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
#    define YYSTYPE_IS_TRIVIAL 1
#    define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if !defined YYLTYPE && !defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE           YYLTYPE;
struct YYLTYPE
{
    int first_line;
    int first_column;
    int last_line;
    int last_column;
};
#    define YYLTYPE_IS_DECLARED 1
#    define YYLTYPE_IS_TRIVIAL 1
#endif


#ifndef YYPUSH_MORE_DEFINED
#    define YYPUSH_MORE_DEFINED
enum
{
    YYPUSH_MORE = 4
};
#endif

typedef struct _gmx_sel_yypstate _gmx_sel_yypstate;

int _gmx_sel_yypush_parse(_gmx_sel_yypstate* ps,
                          int                pushed_char,
                          YYSTYPE const*     pushed_val,
                          YYLTYPE*           pushed_loc,
                          void*              scanner);

_gmx_sel_yypstate* _gmx_sel_yypstate_new();
void               _gmx_sel_yypstate_delete(_gmx_sel_yypstate* ps);

#endif /* !YY__GMX_SEL_YY_PARSER_H_INCLUDED  */
