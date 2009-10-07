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
     HELP = 259,
     HELP_TOPIC = 260,
     INTEGER = 261,
     REAL = 262,
     STR = 263,
     IDENTIFIER = 264,
     CMD_SEP = 265,
     GROUP = 266,
     TO = 267,
     OF = 268,
     VARIABLE_NUMERIC = 269,
     VARIABLE_GROUP = 270,
     VARIABLE_POS = 271,
     KEYWORD_INT = 272,
     KEYWORD_REAL = 273,
     KEYWORD_STR = 274,
     KEYWORD_POS = 275,
     KEYWORD_GROUP = 276,
     METHOD_NUMERIC = 277,
     METHOD_GROUP = 278,
     METHOD_POS = 279,
     MODIFIER = 280,
     PARAM_BOOL = 281,
     BOOL_VALUE = 282,
     PARAM_INT = 283,
     PARAM_REAL = 284,
     PARAM_STR = 285,
     PARAM_POS = 286,
     PARAM_GROUP = 287,
     END_OF_METHOD = 288,
     XOR = 289,
     OR = 290,
     AND = 291,
     NOT = 292,
     CMP_OP = 293
   };
#endif
/* Tokens.  */
#define INVALID 258
#define HELP 259
#define HELP_TOPIC 260
#define INTEGER 261
#define REAL 262
#define STR 263
#define IDENTIFIER 264
#define CMD_SEP 265
#define GROUP 266
#define TO 267
#define OF 268
#define VARIABLE_NUMERIC 269
#define VARIABLE_GROUP 270
#define VARIABLE_POS 271
#define KEYWORD_INT 272
#define KEYWORD_REAL 273
#define KEYWORD_STR 274
#define KEYWORD_POS 275
#define KEYWORD_GROUP 276
#define METHOD_NUMERIC 277
#define METHOD_GROUP 278
#define METHOD_POS 279
#define MODIFIER 280
#define PARAM_BOOL 281
#define BOOL_VALUE 282
#define PARAM_INT 283
#define PARAM_REAL 284
#define PARAM_STR 285
#define PARAM_POS 286
#define PARAM_GROUP 287
#define END_OF_METHOD 288
#define XOR 289
#define OR 290
#define AND 291
#define NOT 292
#define CMP_OP 293




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 71 "parser.y"
{
    int                  i;
    real                 r;
    char                *str;
    struct gmx_ana_selmethod_t *meth;

    struct t_selelem           *sel;

    struct t_selexpr_value     *val;
    struct t_selexpr_param     *param;
}
/* Line 1489 of yacc.c.  */
#line 137 "parser.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



