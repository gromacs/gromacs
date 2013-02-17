/* A Bison parser, made by GNU Bison 2.6.2.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
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

#ifndef _GMX_SEL_YY_PARSER_H
# define _GMX_SEL_YY_PARSER_H
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int _gmx_sel_yydebug;
#endif
/* "%code requires" blocks.  */
/* Line 2049 of yacc.c  */
#line 1 "parser.y"

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012, by the GROMACS development team, led by
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

/* Line 2049 of yacc.c  */
#line 73 "parser.y"

#include "parsetree.h"
#include "selelem.h"


/* Line 2049 of yacc.c  */
#line 90 "parser.h"

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INVALID = 258,
     HELP = 259,
     HELP_TOPIC = 260,
     TOK_INT = 261,
     TOK_REAL = 262,
     STR = 263,
     IDENTIFIER = 264,
     CMD_SEP = 265,
     GROUP = 266,
     TO = 267,
     VARIABLE_NUMERIC = 268,
     VARIABLE_GROUP = 269,
     VARIABLE_POS = 270,
     KEYWORD_NUMERIC = 271,
     KEYWORD_STR = 272,
     KEYWORD_POS = 273,
     KEYWORD_GROUP = 274,
     METHOD_NUMERIC = 275,
     METHOD_GROUP = 276,
     METHOD_POS = 277,
     MODIFIER = 278,
     EMPTY_POSMOD = 279,
     PARAM = 280,
     END_OF_METHOD = 281,
     OF = 282,
     CMP_OP = 283,
     PARAM_REDUCT = 284,
     XOR = 285,
     OR = 286,
     AND = 287,
     NOT = 288,
     UNARY_NEG = 289,
     NUM_REDUCT = 290
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 2049 of yacc.c  */
#line 78 "parser.y"

    int                         i;
    real                        r;
    char                       *str;
    struct gmx_ana_selmethod_t *meth;

    gmx::SelectionStringMatchType                smt;

    gmx::SelectionTreeElementPointer            *sel;
    gmx::SelectionParserValue                   *val;
    gmx::SelectionParserValueListPointer        *vlist;
    gmx::SelectionParserParameter               *param;
    gmx::SelectionParserParameterListPointer    *plist;


/* Line 2049 of yacc.c  */
#line 156 "parser.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


#ifndef YYPUSH_MORE_DEFINED
# define YYPUSH_MORE_DEFINED
enum { YYPUSH_MORE = 4 };
#endif

typedef struct _gmx_sel_yypstate _gmx_sel_yypstate;

#if defined __STDC__ || defined __cplusplus
int _gmx_sel_yypush_parse (_gmx_sel_yypstate *ps, int pushed_char, YYSTYPE const *pushed_val, void *scanner);
#else
int _gmx_sel_yypush_parse ();
#endif

#if defined __STDC__ || defined __cplusplus
_gmx_sel_yypstate * _gmx_sel_yypstate_new (void);
#else
_gmx_sel_yypstate * _gmx_sel_yypstate_new ();
#endif
#if defined __STDC__ || defined __cplusplus
void _gmx_sel_yypstate_delete (_gmx_sel_yypstate *ps);
#else
void _gmx_sel_yypstate_delete ();
#endif

#endif /* !_GMX_SEL_YY_PARSER_H  */
