/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INVALID = 258,
     INT = 259,
     REAL = 260,
     STR = 261,
     IDENTIFIER = 262,
     CMD_SEP = 263,
     GROUP = 264,
     TO = 265,
     OF = 266,
     VARIABLE_NUMERIC = 267,
     VARIABLE_GROUP = 268,
     VARIABLE_POS = 269,
     KEYWORD_INT = 270,
     KEYWORD_REAL = 271,
     KEYWORD_STR = 272,
     KEYWORD_POS = 273,
     KEYWORD_GROUP = 274,
     METHOD_NUMERIC = 275,
     METHOD_GROUP = 276,
     METHOD_POS = 277,
     MODIFIER = 278,
     PARAM_BOOL = 279,
     PARAM_INT = 280,
     PARAM_REAL = 281,
     PARAM_STR = 282,
     PARAM_POS = 283,
     PARAM_GROUP = 284,
     END_OF_METHOD = 285,
     XOR = 286,
     OR = 287,
     AND = 288,
     NOT = 289,
     CMP_OP = 290
   };
#endif
/* Tokens.  */
#define INVALID 258
#define INT 259
#define REAL 260
#define STR 261
#define IDENTIFIER 262
#define CMD_SEP 263
#define GROUP 264
#define TO 265
#define OF 266
#define VARIABLE_NUMERIC 267
#define VARIABLE_GROUP 268
#define VARIABLE_POS 269
#define KEYWORD_INT 270
#define KEYWORD_REAL 271
#define KEYWORD_STR 272
#define KEYWORD_POS 273
#define KEYWORD_GROUP 274
#define METHOD_NUMERIC 275
#define METHOD_GROUP 276
#define METHOD_POS 277
#define MODIFIER 278
#define PARAM_BOOL 279
#define PARAM_INT 280
#define PARAM_REAL 281
#define PARAM_STR 282
#define PARAM_POS 283
#define PARAM_GROUP 284
#define END_OF_METHOD 285
#define XOR 286
#define OR 287
#define AND 288
#define NOT 289
#define CMP_OP 290




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 92 "parser.y"
{
    int                  i;
    real                 r;
    char                *str;
    gmx_ana_selmethod_t *meth;

    t_selelem        *sel;

    t_selexpr_value  *val;
    t_selexpr_param  *param;
}
/* Line 1529 of yacc.c.  */
#line 131 "parser.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



