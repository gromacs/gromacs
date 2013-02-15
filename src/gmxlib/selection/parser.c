/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse _gmx_sel_yybparse
#define yylex   _gmx_sel_yyblex
#define yyerror _gmx_sel_yyberror
#define yylval  _gmx_sel_yyblval
#define yychar  _gmx_sel_yybchar
#define yydebug _gmx_sel_yybdebug
#define yynerrs _gmx_sel_yybnerrs


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
/* Tokens.  */
#define INVALID 258
#define HELP 259
#define HELP_TOPIC 260
#define TOK_INT 261
#define TOK_REAL 262
#define STR 263
#define IDENTIFIER 264
#define CMD_SEP 265
#define GROUP 266
#define TO 267
#define VARIABLE_NUMERIC 268
#define VARIABLE_GROUP 269
#define VARIABLE_POS 270
#define KEYWORD_NUMERIC 271
#define KEYWORD_STR 272
#define KEYWORD_POS 273
#define KEYWORD_GROUP 274
#define METHOD_NUMERIC 275
#define METHOD_GROUP 276
#define METHOD_POS 277
#define MODIFIER 278
#define EMPTY_POSMOD 279
#define PARAM 280
#define END_OF_METHOD 281
#define OF 282
#define CMP_OP 283
#define PARAM_REDUCT 284
#define XOR 285
#define OR 286
#define AND 287
#define NOT 288
#define UNARY_NEG 289
#define NUM_REDUCT 290




/* Copy the first part of user declarations.  */
#line 34 "parser.y"

/*! \internal \file parser.c
 * \brief Generated (from parser.y by Bison) parser for the selection language.
 */
/*! \internal \file parser.h
 * \brief Generated (from parser.y by Bison) parser include file.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string2.h>

#include "parsetree.h"
#include "selelem.h"

#include "scanner.h"

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
static t_selexpr_param *
process_param_list(t_selexpr_param *params);

static void
yyerror(yyscan_t, char const *s);

// Work around compiler warnings that result from bison not correctly
// dealing with stdlib.h with ICC on Windows.
#if (defined __INTEL_COMPILER && defined _WIN32)
#define YYMALLOC malloc
#define YYFREE free
#endif


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 68 "parser.y"
{
    int                         i;
    real                        r;
    char                       *str;
    struct gmx_ana_selmethod_t *meth;

    struct t_selelem           *sel;

    struct t_selexpr_value     *val;
    struct t_selexpr_param     *param;
}
/* Line 193 of yacc.c.  */
#line 220 "parser.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 233 "parser.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   417

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  49
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  26
/* YYNRULES -- Number of rules.  */
#define YYNRULES  91
/* YYNRULES -- Number of states.  */
#define YYNSTATES  150

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   290

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      42,    43,    36,    34,    45,    35,     2,    37,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    41,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    44,     2,    46,    39,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    47,     2,    48,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    38,
      40
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    14,    16,    18,
      20,    22,    25,    29,    33,    37,    39,    41,    44,    47,
      49,    51,    55,    59,    61,    64,    66,    69,    71,    73,
      75,    77,    80,    84,    88,    92,    96,    99,   102,   104,
     106,   109,   113,   117,   121,   123,   125,   128,   132,   136,
     140,   144,   148,   151,   155,   159,   161,   164,   172,   176,
     179,   183,   185,   187,   189,   191,   194,   195,   198,   201,
     202,   204,   208,   210,   213,   217,   219,   223,   225,   228,
     232,   234,   236,   238,   240,   242,   244,   246,   248,   250,
     254,   258
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      50,     0,    -1,    -1,    50,    51,    -1,    52,    10,    -1,
       1,    10,    -1,    -1,    53,    -1,     6,    -1,    59,    -1,
      55,    -1,    59,    55,    -1,     9,    41,    60,    -1,     9,
      41,    62,    -1,     9,    41,    64,    -1,     4,    -1,    54,
      -1,     4,     5,    -1,    54,     5,    -1,    64,    -1,    60,
      -1,    42,    55,    43,    -1,    55,    23,    65,    -1,     6,
      -1,    35,     6,    -1,     7,    -1,    35,     7,    -1,    56,
      -1,    57,    -1,     8,    -1,     9,    -1,    33,    60,    -1,
      60,    32,    60,    -1,    60,    31,    60,    -1,    42,    60,
      43,    -1,    62,    28,    62,    -1,    11,    59,    -1,    11,
       6,    -1,    24,    -1,    18,    -1,    61,    19,    -1,    61,
      17,    70,    -1,    61,    16,    70,    -1,    61,    21,    65,
      -1,     6,    -1,     7,    -1,    61,    16,    -1,    61,    20,
      65,    -1,    62,    34,    62,    -1,    62,    35,    62,    -1,
      62,    36,    62,    -1,    62,    37,    62,    -1,    35,    62,
      -1,    62,    39,    62,    -1,    42,    62,    43,    -1,    59,
      -1,    61,    17,    -1,    44,    58,    45,    58,    45,    58,
      46,    -1,    42,    64,    43,    -1,    22,    65,    -1,    18,
      27,    60,    -1,    14,    -1,    13,    -1,    15,    -1,    66,
      -1,    66,    26,    -1,    -1,    66,    67,    -1,    25,    68,
      -1,    -1,    69,    -1,    47,    69,    48,    -1,    72,    -1,
      69,    72,    -1,    69,    45,    72,    -1,    71,    -1,    47,
      71,    48,    -1,    73,    -1,    71,    73,    -1,    71,    45,
      73,    -1,    60,    -1,    64,    -1,    62,    -1,    63,    -1,
      74,    -1,    56,    -1,    57,    -1,    59,    -1,    74,    -1,
      56,    12,    56,    -1,    56,    12,    57,    -1,    57,    12,
      58,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   188,   188,   189,   198,   199,   219,   223,   224,   233,
     243,   245,   247,   249,   251,   257,   258,   261,   262,   266,
     267,   272,   273,   285,   286,   290,   291,   294,   295,   298,
     299,   307,   313,   319,   331,   335,   343,   349,   357,   358,
     362,   367,   372,   380,   392,   399,   409,   414,   422,   424,
     426,   428,   430,   432,   434,   441,   448,   460,   465,   469,
     477,   488,   492,   496,   505,   507,   512,   513,   518,   525,
     526,   527,   531,   532,   534,   539,   540,   544,   545,   547,
     551,   553,   555,   557,   559,   563,   568,   573,   578,   582,
     587,   592
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INVALID", "HELP", "HELP_TOPIC",
  "TOK_INT", "TOK_REAL", "STR", "IDENTIFIER", "CMD_SEP", "GROUP", "TO",
  "VARIABLE_NUMERIC", "VARIABLE_GROUP", "VARIABLE_POS", "KEYWORD_NUMERIC",
  "KEYWORD_STR", "KEYWORD_POS", "KEYWORD_GROUP", "METHOD_NUMERIC",
  "METHOD_GROUP", "METHOD_POS", "MODIFIER", "EMPTY_POSMOD", "PARAM",
  "END_OF_METHOD", "OF", "CMP_OP", "PARAM_REDUCT", "XOR", "OR", "AND",
  "NOT", "'+'", "'-'", "'*'", "'/'", "UNARY_NEG", "'^'", "NUM_REDUCT",
  "'='", "'('", "')'", "'['", "','", "']'", "'{'", "'}'", "$accept",
  "commands", "command", "cmd_plain", "help_request", "help_topic",
  "selection", "integer_number", "real_number", "number", "string",
  "sel_expr", "pos_mod", "num_expr", "str_expr", "pos_expr",
  "method_params", "method_param_list", "method_param", "value_list",
  "value_list_contents", "basic_value_list", "basic_value_list_contents",
  "value_item", "basic_value_item", "value_item_range", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,    43,    45,    42,    47,   289,    94,
     290,    61,    40,    41,    91,    44,    93,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    49,    50,    50,    51,    51,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    53,    53,    54,    54,    55,
      55,    55,    55,    56,    56,    57,    57,    58,    58,    59,
      59,    60,    60,    60,    60,    60,    60,    60,    61,    61,
      60,    60,    60,    60,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    62,    63,    63,    64,    64,    64,
      64,    60,    62,    64,    65,    65,    66,    66,    67,    68,
      68,    68,    69,    69,    69,    70,    70,    71,    71,    71,
      72,    72,    72,    72,    72,    73,    73,    73,    73,    74,
      74,    74
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     0,     1,     1,     1,
       1,     2,     3,     3,     3,     1,     1,     2,     2,     1,
       1,     3,     3,     1,     2,     1,     2,     1,     1,     1,
       1,     2,     3,     3,     3,     3,     2,     2,     1,     1,
       2,     3,     3,     3,     1,     1,     2,     3,     3,     3,
       3,     3,     2,     3,     3,     1,     2,     7,     3,     2,
       3,     1,     1,     1,     1,     2,     0,     2,     2,     0,
       1,     3,     1,     2,     3,     1,     3,     1,     2,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,    15,    44,    45,    29,    30,     0,
      62,    61,    63,    39,    66,    38,     0,     0,     0,     0,
       3,     0,     7,    16,    10,     9,    20,     0,     0,    19,
       5,    17,     0,    37,    30,    36,     0,    59,    64,    44,
      39,     0,    31,     0,     0,    52,     0,    20,     0,    19,
      23,    25,     0,    27,    28,     0,     4,    18,    66,    11,
       0,     0,    46,     0,    40,    66,    66,     0,     0,     0,
       0,     0,     0,     0,    12,    13,    14,    60,    69,    65,
      67,     0,     0,    46,    21,    34,    54,    58,    24,    26,
       0,    22,    33,    32,     0,    85,    86,    87,    42,    75,
      77,    88,    41,    47,    43,    35,    48,    49,    50,    51,
      53,     0,    44,    45,     0,     0,     0,     0,    55,    80,
       0,    82,    83,    81,    68,    70,    72,    84,     0,     0,
       0,     0,     0,    78,    44,    45,     0,    56,     0,    73,
       0,    76,    89,    90,    91,    79,    71,    74,     0,    57
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    20,    21,    22,    23,    24,    95,    96,    55,
      97,   119,    27,    28,   122,   123,    37,    38,    80,   124,
     125,   102,    99,   126,   100,   101
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -93
static const yytype_int16 yypact[] =
{
     -93,   155,   -93,    10,    19,    26,   -93,   -93,    -3,    73,
     -93,   -93,   -93,    22,   -93,   -93,   356,   372,   317,    11,
     -93,    79,   -93,    86,    70,   317,    29,   384,   180,   -93,
     -93,   -93,   342,   -93,   -93,   -93,   356,   -93,     6,   -93,
     -93,   356,   -93,   372,   -10,    57,   -20,   -17,   256,    54,
     -93,   -93,    88,   -93,   -93,    55,   -93,   -93,   -93,    70,
     356,   356,   197,   174,   -93,   -93,   -93,   372,   372,   372,
     372,   372,   372,   342,    29,   180,   -93,    29,   221,   -93,
     -93,   -17,   223,   -93,   -93,   -93,   -93,   -93,   -93,   -93,
      11,   -93,    69,   -93,    78,    90,    94,   -93,   -93,   244,
     -93,   -93,   -93,   -93,   -93,   267,    36,    36,    57,    57,
      57,    54,    95,    96,   375,   303,    90,    94,   -93,    29,
     392,   267,   -93,   -93,   -93,   263,   -93,   -93,    71,    35,
      11,    11,    78,   -93,   105,   106,   178,   174,   303,   -93,
      11,   -93,   -93,   -93,   -93,   -93,   -93,   -93,    80,   -93
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -93,   -93,   -93,   -93,   -93,   -93,   -13,     0,    14,   -81,
      -1,    87,    -4,   -16,   -93,     3,   -36,   -93,   -93,   -93,
      12,    81,    39,   -91,   -92,   -67
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -27
static const yytype_int16 yytable[] =
{
      25,    45,    48,    58,    29,    46,    83,   133,    35,   128,
      65,   127,    59,    44,    60,    61,    75,    50,    51,    53,
      30,    49,    91,    84,    31,    48,    85,    82,    29,   103,
     104,    78,    79,    54,   139,    76,    -8,   133,    32,    44,
     145,    50,    51,     7,    34,   139,    52,   147,   127,    36,
     144,   105,   106,   107,   108,   109,   110,    48,   127,   148,
      60,    61,   121,    44,    44,    44,    44,    44,    44,   127,
      52,   127,    70,    71,   120,    72,   111,   118,   116,    33,
     132,     7,    34,   141,    50,    51,     7,    34,    26,    56,
      53,    57,   117,    58,    88,    89,    72,    87,    45,   121,
      90,    61,   130,    42,    54,    47,   131,   -23,   -25,   121,
      44,   120,    26,    52,   118,   116,   140,   -24,   -26,    74,
     121,   120,   121,    77,   118,   116,   149,   136,    81,   117,
     142,    53,   120,   129,   120,   118,   116,   118,   116,   117,
      53,     0,     0,    98,   143,    54,     0,    92,    93,     0,
     117,     0,   117,     0,    54,     2,     3,     0,     0,     4,
      81,     5,     6,     7,     8,    -6,     9,     0,    10,    11,
      12,     0,     0,    13,     0,     0,     0,    14,     0,    15,
      50,    51,     7,    34,   112,   113,     7,    34,    16,     9,
      17,    10,    11,    12,     0,     0,    13,    18,     0,    19,
      14,     0,    15,    50,    51,     7,    34,     0,    67,    52,
       0,    16,     0,   114,    68,    69,    70,    71,     0,    72,
      73,    94,    19,   138,     0,     0,   146,   112,   113,     7,
      34,     0,     9,     0,    10,    11,    12,     0,     0,    13,
       0,     0,     0,    14,    94,    15,     0,     0,     0,     0,
      50,    51,     7,    34,    16,     0,   114,    68,    69,    70,
      71,     0,    72,    73,     0,    19,    86,     0,   115,   112,
     113,     7,    34,     0,     9,     0,    10,    11,    12,    52,
       0,    13,     0,     0,    67,    14,     0,    15,     0,   132,
      68,    69,    70,    71,     0,    72,    16,     0,   114,    86,
       0,    68,    69,    70,    71,    73,    72,    19,   138,   112,
     113,     7,    34,     0,     9,     0,    10,    11,    12,     0,
       0,    13,     0,    39,     6,    14,     0,    15,     9,     0,
      10,    11,    12,     0,     0,    13,    16,     0,   114,    14,
       0,    15,     0,     0,     0,    73,     0,    19,    39,     6,
      16,     0,    17,     9,     0,    10,    11,    12,     0,    18,
      13,    19,    39,     6,    14,     0,    15,     9,     0,    10,
      11,     0,     0,     0,    40,    16,     0,    17,    39,     6,
      15,   134,   135,     0,    73,    10,    19,     0,    10,    16,
      40,    17,     0,    40,     0,     0,    15,     0,    41,    15,
      62,    63,     0,    64,    65,    66,     0,    17,    62,   137,
      17,    64,    65,    66,    43,     0,     0,    43
};

static const yytype_int16 yycheck[] =
{
       1,    17,    18,    23,     1,    18,    16,    99,     9,    90,
      20,    78,    25,    17,    31,    32,    32,     6,     7,    19,
      10,    18,    58,    43,     5,    41,    43,    43,    25,    65,
      66,    25,    26,    19,   125,    32,    10,   129,    41,    43,
     132,     6,     7,     8,     9,   136,    35,   138,   115,    27,
     131,    67,    68,    69,    70,    71,    72,    73,   125,   140,
      31,    32,    78,    67,    68,    69,    70,    71,    72,   136,
      35,   138,    36,    37,    78,    39,    73,    78,    78,     6,
      45,     8,     9,    48,     6,     7,     8,     9,     1,    10,
      90,     5,    78,    23,     6,     7,    39,    43,   114,   115,
      45,    32,    12,    16,    90,    18,    12,    12,    12,   125,
     114,   115,    25,    35,   115,   115,    45,    12,    12,    32,
     136,   125,   138,    36,   125,   125,    46,   115,    41,   115,
     130,   131,   136,    94,   138,   136,   136,   138,   138,   125,
     140,    -1,    -1,    62,   130,   131,    -1,    60,    61,    -1,
     136,    -1,   138,    -1,   140,     0,     1,    -1,    -1,     4,
      73,     6,     7,     8,     9,    10,    11,    -1,    13,    14,
      15,    -1,    -1,    18,    -1,    -1,    -1,    22,    -1,    24,
       6,     7,     8,     9,     6,     7,     8,     9,    33,    11,
      35,    13,    14,    15,    -1,    -1,    18,    42,    -1,    44,
      22,    -1,    24,     6,     7,     8,     9,    -1,    28,    35,
      -1,    33,    -1,    35,    34,    35,    36,    37,    -1,    39,
      42,    47,    44,    45,    -1,    -1,    48,     6,     7,     8,
       9,    -1,    11,    -1,    13,    14,    15,    -1,    -1,    18,
      -1,    -1,    -1,    22,    47,    24,    -1,    -1,    -1,    -1,
       6,     7,     8,     9,    33,    -1,    35,    34,    35,    36,
      37,    -1,    39,    42,    -1,    44,    43,    -1,    47,     6,
       7,     8,     9,    -1,    11,    -1,    13,    14,    15,    35,
      -1,    18,    -1,    -1,    28,    22,    -1,    24,    -1,    45,
      34,    35,    36,    37,    -1,    39,    33,    -1,    35,    43,
      -1,    34,    35,    36,    37,    42,    39,    44,    45,     6,
       7,     8,     9,    -1,    11,    -1,    13,    14,    15,    -1,
      -1,    18,    -1,     6,     7,    22,    -1,    24,    11,    -1,
      13,    14,    15,    -1,    -1,    18,    33,    -1,    35,    22,
      -1,    24,    -1,    -1,    -1,    42,    -1,    44,     6,     7,
      33,    -1,    35,    11,    -1,    13,    14,    15,    -1,    42,
      18,    44,     6,     7,    22,    -1,    24,    11,    -1,    13,
      14,    -1,    -1,    -1,    18,    33,    -1,    35,     6,     7,
      24,     6,     7,    -1,    42,    13,    44,    -1,    13,    33,
      18,    35,    -1,    18,    -1,    -1,    24,    -1,    42,    24,
      16,    17,    -1,    19,    20,    21,    -1,    35,    16,    17,
      35,    19,    20,    21,    42,    -1,    -1,    42
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    50,     0,     1,     4,     6,     7,     8,     9,    11,
      13,    14,    15,    18,    22,    24,    33,    35,    42,    44,
      51,    52,    53,    54,    55,    59,    60,    61,    62,    64,
      10,     5,    41,     6,     9,    59,    27,    65,    66,     6,
      18,    42,    60,    42,    61,    62,    55,    60,    62,    64,
       6,     7,    35,    56,    57,    58,    10,     5,    23,    55,
      31,    32,    16,    17,    19,    20,    21,    28,    34,    35,
      36,    37,    39,    42,    60,    62,    64,    60,    25,    26,
      67,    60,    62,    16,    43,    43,    43,    43,     6,     7,
      45,    65,    60,    60,    47,    56,    57,    59,    70,    71,
      73,    74,    70,    65,    65,    62,    62,    62,    62,    62,
      62,    64,     6,     7,    35,    47,    56,    57,    59,    60,
      61,    62,    63,    64,    68,    69,    72,    74,    58,    71,
      12,    12,    45,    73,     6,     7,    69,    17,    45,    72,
      45,    48,    56,    57,    58,    73,    48,    72,    58,    46
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (scanner, YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval, scanner)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value, scanner); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, yyscan_t                 scanner)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    yyscan_t                 scanner;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (scanner);
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, yyscan_t                 scanner)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    yyscan_t                 scanner;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, yyscan_t                 scanner)
#else
static void
yy_reduce_print (yyvsp, yyrule, scanner)
    YYSTYPE *yyvsp;
    int yyrule;
    yyscan_t                 scanner;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       , scanner);
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, scanner); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, yyscan_t                 scanner)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, scanner)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    yyscan_t                 scanner;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (scanner);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {
      case 5: /* "HELP_TOPIC" */
#line 167 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1323 "parser.c"
	break;
      case 8: /* "STR" */
#line 167 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1328 "parser.c"
	break;
      case 9: /* "IDENTIFIER" */
#line 167 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1333 "parser.c"
	break;
      case 25: /* "PARAM" */
#line 168 "parser.y"
	{ if((yyvaluep->str)) free((yyvaluep->str));              };
#line 1338 "parser.c"
	break;
      case 28: /* "CMP_OP" */
#line 167 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1343 "parser.c"
	break;
      case 51: /* "command" */
#line 169 "parser.y"
	{ if((yyvaluep->sel)) _gmx_selelem_free((yyvaluep->sel)); };
#line 1348 "parser.c"
	break;
      case 52: /* "cmd_plain" */
#line 169 "parser.y"
	{ if((yyvaluep->sel)) _gmx_selelem_free((yyvaluep->sel)); };
#line 1353 "parser.c"
	break;
      case 55: /* "selection" */
#line 170 "parser.y"
	{ _gmx_selelem_free_chain((yyvaluep->sel));  };
#line 1358 "parser.c"
	break;
      case 59: /* "string" */
#line 167 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1363 "parser.c"
	break;
      case 60: /* "sel_expr" */
#line 171 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1368 "parser.c"
	break;
      case 62: /* "num_expr" */
#line 171 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1373 "parser.c"
	break;
      case 63: /* "str_expr" */
#line 171 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1378 "parser.c"
	break;
      case 64: /* "pos_expr" */
#line 171 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1383 "parser.c"
	break;
      case 65: /* "method_params" */
#line 172 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1388 "parser.c"
	break;
      case 66: /* "method_param_list" */
#line 172 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1393 "parser.c"
	break;
      case 67: /* "method_param" */
#line 172 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1398 "parser.c"
	break;
      case 68: /* "value_list" */
#line 173 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1403 "parser.c"
	break;
      case 69: /* "value_list_contents" */
#line 173 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1408 "parser.c"
	break;
      case 70: /* "basic_value_list" */
#line 174 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1413 "parser.c"
	break;
      case 71: /* "basic_value_list_contents" */
#line 174 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1418 "parser.c"
	break;
      case 72: /* "value_item" */
#line 173 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1423 "parser.c"
	break;
      case 73: /* "basic_value_item" */
#line 174 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1428 "parser.c"
	break;
      case 74: /* "value_item_range" */
#line 173 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1433 "parser.c"
	break;

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (yyscan_t                 scanner);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */






/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (yyscan_t                 scanner)
#else
int
yyparse (scanner)
    yyscan_t                 scanner;
#endif
#endif
{
  /* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 188 "parser.y"
    { (yyval.sel) = NULL ;}
    break;

  case 3:
#line 190 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_append_selection((yyvsp[(2) - (2)].sel), (yyvsp[(1) - (2)].sel), scanner);
                 if (_gmx_sel_parser_should_finish(scanner))
                     YYACCEPT;
             ;}
    break;

  case 4:
#line 198 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (2)].sel); ;}
    break;

  case 5:
#line 200 "parser.y"
    {
                 (yyval.sel) = NULL;
                 _gmx_selparser_error("invalid selection '%s'",
                                      _gmx_sel_lexer_pselstr(scanner));
                 _gmx_sel_lexer_clear_method_stack(scanner);
                 if (_gmx_sel_is_lexer_interactive(scanner))
                 {
                     _gmx_sel_lexer_clear_pselstr(scanner);
                     yyerrok;
                 }
                 else
                 {
                     YYABORT;
                 }
             ;}
    break;

  case 6:
#line 219 "parser.y"
    {
                 (yyval.sel) = NULL;
                 _gmx_sel_handle_empty_cmd(scanner);
             ;}
    break;

  case 7:
#line 223 "parser.y"
    { (yyval.sel) = NULL; ;}
    break;

  case 8:
#line 225 "parser.y"
    {
                 t_selelem *s, *p;
                 s = _gmx_sel_init_group_by_id((yyvsp[(1) - (1)].i), scanner);
                 if (s == NULL) YYERROR;
                 p = _gmx_sel_init_position(s, NULL, scanner);
                 if (p == NULL) YYERROR;
                 (yyval.sel) = _gmx_sel_init_selection(strdup(s->name), p, scanner);
             ;}
    break;

  case 9:
#line 234 "parser.y"
    {
                 t_selelem *s, *p;
                 s = _gmx_sel_init_group_by_name((yyvsp[(1) - (1)].str), scanner);
                 free((yyvsp[(1) - (1)].str));
                 if (s == NULL) YYERROR;
                 p = _gmx_sel_init_position(s, NULL, scanner);
                 if (p == NULL) YYERROR;
                 (yyval.sel) = _gmx_sel_init_selection(strdup(s->name), p, scanner);
             ;}
    break;

  case 10:
#line 244 "parser.y"
    { (yyval.sel) = _gmx_sel_init_selection(NULL, (yyvsp[(1) - (1)].sel), scanner); ;}
    break;

  case 11:
#line 246 "parser.y"
    { (yyval.sel) = _gmx_sel_init_selection((yyvsp[(1) - (2)].str), (yyvsp[(2) - (2)].sel), scanner);   ;}
    break;

  case 12:
#line 248 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel), scanner);  ;}
    break;

  case 13:
#line 250 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel), scanner);  ;}
    break;

  case 14:
#line 252 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel), scanner);  ;}
    break;

  case 15:
#line 257 "parser.y"
    { _gmx_sel_handle_help_cmd(NULL, scanner); ;}
    break;

  case 17:
#line 261 "parser.y"
    { _gmx_sel_handle_help_cmd((yyvsp[(2) - (2)].str), scanner); ;}
    break;

  case 18:
#line 262 "parser.y"
    { _gmx_sel_handle_help_cmd((yyvsp[(2) - (2)].str), scanner); ;}
    break;

  case 19:
#line 266 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 20:
#line 268 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_position((yyvsp[(1) - (1)].sel), NULL, scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 21:
#line 272 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); ;}
    break;

  case 22:
#line 274 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_modifier((yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].param), (yyvsp[(1) - (3)].sel), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 23:
#line 285 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].i); ;}
    break;

  case 24:
#line 286 "parser.y"
    { (yyval.r) = -(yyvsp[(2) - (2)].i); ;}
    break;

  case 25:
#line 290 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); ;}
    break;

  case 26:
#line 291 "parser.y"
    { (yyval.r) = -(yyvsp[(2) - (2)].r); ;}
    break;

  case 27:
#line 294 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); ;}
    break;

  case 28:
#line 295 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); ;}
    break;

  case 29:
#line 298 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); ;}
    break;

  case 30:
#line 299 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); ;}
    break;

  case 31:
#line 308 "parser.y"
    {
                 (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                 (yyval.sel)->u.boolt = BOOL_NOT;
                 (yyval.sel)->child = (yyvsp[(2) - (2)].sel);
             ;}
    break;

  case 32:
#line 314 "parser.y"
    {
                 (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                 (yyval.sel)->u.boolt = BOOL_AND;
                 (yyval.sel)->child = (yyvsp[(1) - (3)].sel); (yyval.sel)->child->next = (yyvsp[(3) - (3)].sel);
             ;}
    break;

  case 33:
#line 320 "parser.y"
    {
                 (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                 (yyval.sel)->u.boolt = BOOL_OR;
                 (yyval.sel)->child = (yyvsp[(1) - (3)].sel); (yyval.sel)->child->next = (yyvsp[(3) - (3)].sel);
             ;}
    break;

  case 34:
#line 331 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); ;}
    break;

  case 35:
#line 336 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_comparison((yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), (yyvsp[(2) - (3)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 36:
#line 344 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_group_by_name((yyvsp[(2) - (2)].str), scanner);
                 free((yyvsp[(2) - (2)].str));
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 37:
#line 350 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_group_by_id((yyvsp[(2) - (2)].i), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 38:
#line 357 "parser.y"
    { (yyval.str) = NULL; ;}
    break;

  case 39:
#line 358 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str);   ;}
    break;

  case 40:
#line 363 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 41:
#line 368 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_keyword((yyvsp[(2) - (3)].meth), process_value_list((yyvsp[(3) - (3)].val), NULL), (yyvsp[(1) - (3)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 42:
#line 373 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_keyword((yyvsp[(2) - (3)].meth), process_value_list((yyvsp[(3) - (3)].val), NULL), (yyvsp[(1) - (3)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 43:
#line 381 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_method((yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].param), (yyvsp[(1) - (3)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 44:
#line 393 "parser.y"
    {
                 (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                 _gmx_selelem_set_vtype((yyval.sel), INT_VALUE);
                 _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                 (yyval.sel)->v.u.i[0] = (yyvsp[(1) - (1)].i);
             ;}
    break;

  case 45:
#line 400 "parser.y"
    {
                 (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                 _gmx_selelem_set_vtype((yyval.sel), REAL_VALUE);
                 _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                 (yyval.sel)->v.u.r[0] = (yyvsp[(1) - (1)].r);
             ;}
    break;

  case 46:
#line 410 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 47:
#line 415 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_method((yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].param), (yyvsp[(1) - (3)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 48:
#line 423 "parser.y"
    { (yyval.sel) = _gmx_sel_init_arithmetic((yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), '+', scanner); ;}
    break;

  case 49:
#line 425 "parser.y"
    { (yyval.sel) = _gmx_sel_init_arithmetic((yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), '-', scanner); ;}
    break;

  case 50:
#line 427 "parser.y"
    { (yyval.sel) = _gmx_sel_init_arithmetic((yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), '*', scanner); ;}
    break;

  case 51:
#line 429 "parser.y"
    { (yyval.sel) = _gmx_sel_init_arithmetic((yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), '/', scanner); ;}
    break;

  case 52:
#line 431 "parser.y"
    { (yyval.sel) = _gmx_sel_init_arithmetic((yyvsp[(2) - (2)].sel), NULL, '-', scanner); ;}
    break;

  case 53:
#line 433 "parser.y"
    { (yyval.sel) = _gmx_sel_init_arithmetic((yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), '^', scanner); ;}
    break;

  case 54:
#line 434 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); ;}
    break;

  case 55:
#line 442 "parser.y"
    {
                 (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                 _gmx_selelem_set_vtype((yyval.sel), STR_VALUE);
                 _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                 (yyval.sel)->v.u.s[0] = (yyvsp[(1) - (1)].str);
             ;}
    break;

  case 56:
#line 449 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 57:
#line 461 "parser.y"
    { (yyval.sel) = _gmx_sel_init_const_position((yyvsp[(2) - (7)].r), (yyvsp[(4) - (7)].r), (yyvsp[(6) - (7)].r)); ;}
    break;

  case 58:
#line 465 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); ;}
    break;

  case 59:
#line 470 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_method((yyvsp[(1) - (2)].meth), (yyvsp[(2) - (2)].param), NULL, scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 60:
#line 478 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_position((yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].str), scanner);
                 if ((yyval.sel) == NULL) YYERROR;
             ;}
    break;

  case 61:
#line 489 "parser.y"
    { (yyval.sel) = _gmx_sel_init_variable_ref((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 62:
#line 493 "parser.y"
    { (yyval.sel) = _gmx_sel_init_variable_ref((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 63:
#line 497 "parser.y"
    { (yyval.sel) = _gmx_sel_init_variable_ref((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 64:
#line 506 "parser.y"
    { (yyval.param) = process_param_list((yyvsp[(1) - (1)].param)); ;}
    break;

  case 65:
#line 508 "parser.y"
    { (yyval.param) = process_param_list((yyvsp[(1) - (2)].param)); ;}
    break;

  case 66:
#line 512 "parser.y"
    { (yyval.param) = NULL;              ;}
    break;

  case 67:
#line 514 "parser.y"
    { (yyvsp[(2) - (2)].param)->next = (yyvsp[(1) - (2)].param); (yyval.param) = (yyvsp[(2) - (2)].param); ;}
    break;

  case 68:
#line 519 "parser.y"
    {
                 (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                 (yyval.param)->value = process_value_list((yyvsp[(2) - (2)].val), &(yyval.param)->nval);
             ;}
    break;

  case 69:
#line 525 "parser.y"
    { (yyval.val) = NULL; ;}
    break;

  case 70:
#line 526 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val);   ;}
    break;

  case 71:
#line 527 "parser.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val);   ;}
    break;

  case 72:
#line 531 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); ;}
    break;

  case 73:
#line 533 "parser.y"
    { (yyvsp[(2) - (2)].val)->next = (yyvsp[(1) - (2)].val); (yyval.val) = (yyvsp[(2) - (2)].val); ;}
    break;

  case 74:
#line 535 "parser.y"
    { (yyvsp[(3) - (3)].val)->next = (yyvsp[(1) - (3)].val); (yyval.val) = (yyvsp[(3) - (3)].val); ;}
    break;

  case 75:
#line 539 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); ;}
    break;

  case 76:
#line 540 "parser.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val); ;}
    break;

  case 77:
#line 544 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); ;}
    break;

  case 78:
#line 546 "parser.y"
    { (yyvsp[(2) - (2)].val)->next = (yyvsp[(1) - (2)].val); (yyval.val) = (yyvsp[(2) - (2)].val); ;}
    break;

  case 79:
#line 548 "parser.y"
    { (yyvsp[(3) - (3)].val)->next = (yyvsp[(1) - (3)].val); (yyval.val) = (yyvsp[(3) - (3)].val); ;}
    break;

  case 80:
#line 552 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value_expr((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 81:
#line 554 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value_expr((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 82:
#line 556 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value_expr((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 83:
#line 558 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value_expr((yyvsp[(1) - (1)].sel)); ;}
    break;

  case 84:
#line 559 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); ;}
    break;

  case 85:
#line 564 "parser.y"
    {
                 (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                 (yyval.val)->u.i.i1 = (yyval.val)->u.i.i2 = (yyvsp[(1) - (1)].r);
             ;}
    break;

  case 86:
#line 569 "parser.y"
    {
                 (yyval.val) = _gmx_selexpr_create_value(REAL_VALUE);
                 (yyval.val)->u.r.r1 = (yyval.val)->u.r.r2 = (yyvsp[(1) - (1)].r);
             ;}
    break;

  case 87:
#line 574 "parser.y"
    {
                 (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                 (yyval.val)->u.s = (yyvsp[(1) - (1)].str);
             ;}
    break;

  case 88:
#line 578 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); ;}
    break;

  case 89:
#line 583 "parser.y"
    {
                 (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                 (yyval.val)->u.i.i1 = (yyvsp[(1) - (3)].r); (yyval.val)->u.i.i2 = (yyvsp[(3) - (3)].r);
             ;}
    break;

  case 90:
#line 588 "parser.y"
    {
                 (yyval.val) = _gmx_selexpr_create_value(REAL_VALUE);
                 (yyval.val)->u.r.r1 = (yyvsp[(1) - (3)].r); (yyval.val)->u.r.r2 = (yyvsp[(3) - (3)].r);
             ;}
    break;

  case 91:
#line 593 "parser.y"
    {
                 (yyval.val) = _gmx_selexpr_create_value(REAL_VALUE);
                 (yyval.val)->u.r.r1 = (yyvsp[(1) - (3)].r); (yyval.val)->u.r.r2 = (yyvsp[(3) - (3)].r);
             ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2315 "parser.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (scanner, YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (scanner, yymsg);
	  }
	else
	  {
	    yyerror (scanner, YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval, scanner);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp, scanner);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, scanner);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, scanner);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 599 "parser.y"


static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr)
{
    t_selexpr_value *val, *pval, *nval;

    /* Count values (if needed) and reverse list */
    if (nr)
    {
        *nr  = 0;
    }
    pval = NULL;
    val  = values;
    while (val)
    {
        if (nr)
        {
            ++*nr;
        }
        nval = val->next;
        val->next = pval;
        pval = val;
        val = nval;
    }
    values = pval;

    return values;
}

static t_selexpr_param *
process_param_list(t_selexpr_param *params)
{
    t_selexpr_param *par, *ppar, *npar;

    /* Reverse list */
    ppar = NULL;
    par  = params;
    while (par)
    {
        npar = par->next;
        par->next = ppar;
        ppar = par;
        par = npar;
    }
    params = ppar;

    return params;
}

static void
yyerror(yyscan_t scanner, char const *s)
{
    _gmx_selparser_error("%s", s);
}



