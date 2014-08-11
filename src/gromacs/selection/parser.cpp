/* A Bison parser, made by GNU Bison 2.6.2.  */

/* Bison implementation for Yacc-like parsers in C
   
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
#define YYBISON_VERSION "2.6.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 1

/* Pull parsers.  */
#define YYPULL 0


/* Substitute the variable and function names.  */
#define yypush_parse    _gmx_sel_yypush_parse
#define yypstate_new    _gmx_sel_yypstate_new
#define yypstate_delete _gmx_sel_yypstate_delete
#define yypstate        _gmx_sel_yypstate
#define yylex           _gmx_sel_yylex
#define yyerror         _gmx_sel_yyerror
#define yylval          _gmx_sel_yylval
#define yychar          _gmx_sel_yychar
#define yydebug         _gmx_sel_yydebug
#define yynerrs         _gmx_sel_yynerrs

/* Copy the first part of user declarations.  */
/* Line 336 of yacc.c  */
#line 43 "parser.y"

/*! \internal \file parser.cpp
 * \brief Generated (from parser.y by Bison) parser for the selection language.
 *
 * \ingroup module_selection
 */
/*! \internal \file parser.h
 * \brief Generated (from parser.y by Bison) parser include file.
 *
 * \ingroup module_selection
 */
#include "gromacs/utility/scoped_ptr_sfree.h"

#include "parser_internal.h"

using gmx::scoped_ptr_sfree;
using gmx::SelectionParserValue;
using gmx::SelectionParserValueList;
using gmx::SelectionParserValueListPointer;
using gmx::SelectionParserParameter;
using gmx::SelectionParserParameterList;
using gmx::SelectionParserParameterListPointer;
using gmx::SelectionTreeElement;
using gmx::SelectionTreeElementPointer;

#ifdef _MSC_VER
#pragma warning(disable: 4065)
#endif

/* Line 336 of yacc.c  */
#line 108 "parser.cpp"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser.h".  */
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
/* Line 350 of yacc.c  */
#line 1 "parser.y"

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

/* Line 350 of yacc.c  */
#line 73 "parser.y"

#include "parsetree.h"
#include "selelem.h"


/* Line 350 of yacc.c  */
#line 184 "parser.cpp"

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INVALID = 258,
     TOK_INT = 259,
     TOK_REAL = 260,
     STR = 261,
     IDENTIFIER = 262,
     CMD_SEP = 263,
     GROUP = 264,
     TO = 265,
     VARIABLE_NUMERIC = 266,
     VARIABLE_GROUP = 267,
     VARIABLE_POS = 268,
     KEYWORD_NUMERIC = 269,
     KEYWORD_STR = 270,
     KEYWORD_POS = 271,
     KEYWORD_GROUP = 272,
     METHOD_NUMERIC = 273,
     METHOD_GROUP = 274,
     METHOD_POS = 275,
     MODIFIER = 276,
     EMPTY_POSMOD = 277,
     PARAM = 278,
     END_OF_METHOD = 279,
     OF = 280,
     CMP_OP = 281,
     PARAM_REDUCT = 282,
     XOR = 283,
     OR = 284,
     AND = 285,
     NOT = 286,
     UNARY_NEG = 287,
     NUM_REDUCT = 288
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 350 of yacc.c  */
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


/* Line 350 of yacc.c  */
#line 248 "parser.cpp"
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

/* Copy the second part of user declarations.  */

/* Line 353 of yacc.c  */
#line 285 "parser.cpp"

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
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
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
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   378

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  49
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  25
/* YYNRULES -- Number of rules.  */
#define YYNRULES  89
/* YYNRULES -- Number of states.  */
#define YYNSTATES  149

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   288

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      40,    41,    34,    32,    45,    33,     2,    35,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    39,     2,    43,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    44,     2,    46,    37,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    47,     2,    48,    42,     2,     2,     2,
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
      25,    26,    27,    28,    29,    30,    31,    36,    38
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    14,    16,    18,
      20,    23,    27,    31,    35,    37,    39,    43,    47,    49,
      52,    54,    57,    59,    61,    63,    65,    68,    72,    76,
      80,    84,    87,    90,    92,    94,    96,    98,   100,   103,
     107,   112,   116,   120,   122,   124,   127,   131,   135,   139,
     143,   147,   150,   154,   158,   160,   163,   171,   175,   178,
     182,   184,   186,   188,   190,   193,   194,   197,   200,   202,
     206,   207,   210,   214,   216,   220,   222,   225,   229,   231,
     233,   235,   237,   239,   241,   243,   245,   247,   251,   255
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      50,     0,    -1,    -1,    50,    51,    -1,    52,     8,    -1,
       1,     8,    -1,    -1,     4,    -1,    57,    -1,    53,    -1,
      57,    53,    -1,     7,    39,    58,    -1,     7,    39,    61,
      -1,     7,    39,    63,    -1,    63,    -1,    58,    -1,    40,
      53,    41,    -1,    53,    21,    64,    -1,     4,    -1,    33,
       4,    -1,     5,    -1,    33,     5,    -1,    54,    -1,    55,
      -1,     6,    -1,     7,    -1,    31,    58,    -1,    58,    30,
      58,    -1,    58,    29,    58,    -1,    40,    58,    41,    -1,
      61,    26,    61,    -1,     9,    57,    -1,     9,     4,    -1,
      22,    -1,    16,    -1,    42,    -1,    43,    -1,    39,    -1,
      59,    17,    -1,    59,    15,    69,    -1,    59,    15,    60,
      69,    -1,    59,    14,    69,    -1,    59,    19,    64,    -1,
       4,    -1,     5,    -1,    59,    14,    -1,    59,    18,    64,
      -1,    61,    32,    61,    -1,    61,    33,    61,    -1,    61,
      34,    61,    -1,    61,    35,    61,    -1,    33,    61,    -1,
      61,    37,    61,    -1,    40,    61,    41,    -1,    57,    -1,
      59,    15,    -1,    44,    56,    45,    56,    45,    56,    46,
      -1,    40,    63,    41,    -1,    20,    64,    -1,    16,    25,
      58,    -1,    12,    -1,    11,    -1,    13,    -1,    65,    -1,
      65,    24,    -1,    -1,    65,    66,    -1,    23,    67,    -1,
      68,    -1,    47,    68,    48,    -1,    -1,    68,    71,    -1,
      68,    45,    71,    -1,    70,    -1,    47,    70,    48,    -1,
      72,    -1,    70,    72,    -1,    70,    45,    72,    -1,    58,
      -1,    63,    -1,    61,    -1,    62,    -1,    73,    -1,    54,
      -1,    55,    -1,    57,    -1,    73,    -1,    54,    10,    54,
      -1,    54,    10,    55,    -1,    55,    10,    56,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   193,   193,   198,   209,   210,   232,   237,   248,   260,
     266,   273,   280,   287,   297,   298,   305,   306,   320,   321,
     325,   326,   329,   330,   333,   334,   342,   353,   364,   375,
     379,   390,   397,   406,   407,   412,   413,   414,   418,   426,
     434,   442,   453,   468,   479,   493,   501,   512,   518,   524,
     530,   536,   542,   548,   555,   566,   581,   590,   594,   604,
     618,   626,   634,   647,   649,   655,   660,   671,   680,   681,
     686,   691,   699,   710,   711,   715,   721,   729,   739,   745,
     751,   757,   763,   767,   773,   779,   786,   790,   796,   802
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INVALID", "TOK_INT", "TOK_REAL", "STR",
  "IDENTIFIER", "CMD_SEP", "GROUP", "TO", "VARIABLE_NUMERIC",
  "VARIABLE_GROUP", "VARIABLE_POS", "KEYWORD_NUMERIC", "KEYWORD_STR",
  "KEYWORD_POS", "KEYWORD_GROUP", "METHOD_NUMERIC", "METHOD_GROUP",
  "METHOD_POS", "MODIFIER", "EMPTY_POSMOD", "PARAM", "END_OF_METHOD", "OF",
  "CMP_OP", "PARAM_REDUCT", "XOR", "OR", "AND", "NOT", "'+'", "'-'", "'*'",
  "'/'", "UNARY_NEG", "'^'", "NUM_REDUCT", "'='", "'('", "')'", "'~'",
  "'?'", "'['", "','", "']'", "'{'", "'}'", "$accept", "commands",
  "command", "cmd_plain", "selection", "integer_number", "real_number",
  "number", "string", "sel_expr", "pos_mod", "str_match_type", "num_expr",
  "str_expr", "pos_expr", "method_params", "method_param_list",
  "method_param", "value_list", "value_list_contents", "basic_value_list",
  "basic_value_list_contents", "value_item", "basic_value_item",
  "value_item_range", YY_NULL
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
     285,   286,    43,    45,    42,    47,   287,    94,   288,    61,
      40,    41,   126,    63,    91,    44,    93,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    49,    50,    50,    51,    51,    52,    52,    52,    52,
      52,    52,    52,    52,    53,    53,    53,    53,    54,    54,
      55,    55,    56,    56,    57,    57,    58,    58,    58,    58,
      58,    58,    58,    59,    59,    60,    60,    60,    58,    58,
      58,    58,    58,    61,    61,    61,    61,    61,    61,    61,
      61,    61,    61,    61,    62,    62,    63,    63,    63,    63,
      58,    61,    63,    64,    64,    65,    65,    66,    67,    67,
      68,    68,    68,    69,    69,    70,    70,    70,    71,    71,
      71,    71,    71,    72,    72,    72,    72,    73,    73,    73
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     0,     1,     1,     1,
       2,     3,     3,     3,     1,     1,     3,     3,     1,     2,
       1,     2,     1,     1,     1,     1,     2,     3,     3,     3,
       3,     2,     2,     1,     1,     1,     1,     1,     2,     3,
       4,     3,     3,     1,     1,     2,     3,     3,     3,     3,
       3,     2,     3,     3,     1,     2,     7,     3,     2,     3,
       1,     1,     1,     1,     2,     0,     2,     2,     1,     3,
       0,     2,     3,     1,     3,     1,     2,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     3,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,    43,    44,    24,    25,     0,    61,
      60,    62,    34,    65,    33,     0,     0,     0,     0,     3,
       0,     9,     8,    15,     0,     0,    14,     5,     0,    32,
      25,    31,     0,    58,    63,    43,    34,     0,    26,     0,
       0,    51,     0,    15,     0,    14,    18,    20,     0,    22,
      23,     0,     4,    65,    10,     0,     0,    45,     0,    38,
      65,    65,     0,     0,     0,     0,     0,     0,     0,    11,
      12,    13,    59,    70,    64,    66,     0,     0,    45,    16,
      29,    53,    57,    19,    21,     0,    17,    28,    27,     0,
      83,    84,    85,    41,    73,    75,    86,    37,    35,    36,
       0,    39,    46,    42,    30,    47,    48,    49,    50,    52,
       0,    70,    67,    68,     0,     0,     0,     0,     0,    76,
      40,     0,    43,    44,     0,     0,     0,     0,    54,    78,
       0,    80,    81,    79,    71,    82,     0,    74,    87,    88,
      89,    77,    69,    43,    44,    72,    55,     0,    56
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,    19,    20,    21,    90,    91,    51,    92,   129,
      24,   100,    25,   132,   133,    33,    34,    75,   112,   113,
     101,    94,   134,    95,    96
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -84
static const yytype_int16 yypact[] =
{
     -84,   148,   -84,     4,     6,   -84,   -84,   -20,   189,   -84,
     -84,   -84,     1,   -84,   -84,   322,   216,   285,     0,   -84,
      27,    16,   285,    22,   122,   173,   -84,   -84,   308,   -84,
     -84,   -84,   322,   -84,    71,   -84,   -84,   322,   -84,   216,
      49,    24,   -11,   -13,   232,    23,   -84,   -84,   140,   -84,
     -84,    21,   -84,   -84,    16,   322,   322,    40,   197,   -84,
     -84,   -84,   216,   216,   216,   216,   216,   216,   308,    22,
     173,   -84,    22,    30,   -84,   -84,   -13,    47,   -84,   -84,
     -84,   -84,   -84,   -84,   -84,     0,   -84,    56,   -84,    85,
      93,    99,   -84,   -84,   210,   -84,   -84,   -84,   -84,   -84,
      36,   -84,   -84,   -84,   139,    94,    94,    24,    24,    24,
      23,   -84,   -84,   241,    78,    17,     0,     0,    85,   -84,
     -84,   178,   116,   125,   338,   271,    93,    99,   -84,    22,
     349,   139,   -84,   -84,   -84,   -84,     0,   -84,   -84,   -84,
     -84,   -84,   -84,   133,   136,   -84,   197,   104,   -84
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -84,   -84,   -84,   -84,    -4,   -17,   -15,   -83,    -1,   110,
       9,   -84,    -8,   -84,    10,    54,   -84,   -84,   -84,    51,
     -51,    62,    33,   -79,   -28
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -22
static const yytype_int16 yytable[] =
{
      22,    49,   114,    50,    46,    47,    93,    31,    41,    44,
      53,    26,    27,    42,    -7,   119,    55,    56,    54,    28,
      70,    46,    47,     6,    30,    40,    32,    45,    80,    44,
      79,    77,    26,    48,   140,    52,   119,    53,    71,   141,
      46,    47,     6,    30,    46,    47,     6,    30,    40,   120,
      48,    55,    56,   147,   104,   105,   106,   107,   108,   109,
      44,    67,   118,    78,    82,   137,    85,    60,    49,    48,
      50,    40,    40,    40,    40,    40,    40,   111,   110,    63,
      64,    65,    66,    89,    67,   135,    56,    89,    81,    46,
      47,     6,    30,   135,    73,    74,   126,   135,   127,   138,
      49,   139,    50,   116,   126,   131,   127,    86,   126,   117,
     127,    23,   128,   131,   102,   103,    41,   131,    48,    49,
     128,    50,   130,   136,   128,    38,   -18,    43,    65,    66,
     130,    67,    23,    40,   130,   -20,    57,    58,    69,    59,
      60,    61,    72,   -19,    83,    84,   -21,    76,     2,     3,
     148,   115,     4,     5,     6,     7,    -6,     8,   145,     9,
      10,    11,   121,     0,    12,    87,    88,     0,    13,     0,
      14,    63,    64,    65,    66,     0,    67,     0,    76,    15,
       0,    16,   122,   123,     6,    30,     0,     8,    17,     9,
      10,    11,    18,    29,    12,     6,    30,     0,    13,    62,
      14,    46,    47,     6,    30,    63,    64,    65,    66,    15,
      67,   124,     0,     0,    46,    47,     6,    30,    68,     0,
      35,     5,    18,   125,     0,     0,   142,     9,     0,     0,
      48,     0,    36,     0,     0,     0,    97,     0,    14,    98,
      99,     0,     0,    48,    89,   122,   123,     6,    30,    16,
       8,     0,     9,    10,    11,   118,    39,    12,    62,     0,
       0,    13,     0,    14,    63,    64,    65,    66,     0,    67,
       0,     0,    15,    81,   124,   122,   123,     6,    30,     0,
       8,    68,     9,    10,    11,    18,   125,    12,     0,    35,
       5,    13,     0,    14,     8,     0,     9,    10,    11,     0,
       0,    12,    15,     0,   124,    13,     0,    14,     0,     0,
       0,    68,    35,     5,     0,    18,    15,     8,    16,     9,
      10,    11,     0,     0,    12,    17,    35,     5,    13,    18,
      14,     8,     0,     9,    10,     0,     0,     0,    36,    15,
       0,    16,   143,   144,    14,     0,     0,     0,    68,     9,
       0,     0,    18,    15,    36,    16,     0,     0,     0,     0,
      14,     0,    37,    57,   146,     0,    59,    60,    61,     0,
       0,    16,     0,     0,     0,     0,     0,     0,    39
};

#define yypact_value_is_default(yystate) \
  ((yystate) == (-84))

#define yytable_value_is_error(yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
       1,    18,    85,    18,     4,     5,    57,     8,    16,    17,
      21,     1,     8,    17,     8,    94,    29,    30,    22,    39,
      28,     4,     5,     6,     7,    16,    25,    17,    41,    37,
      41,    39,    22,    33,   117,     8,   115,    21,    28,   118,
       4,     5,     6,     7,     4,     5,     6,     7,    39,   100,
      33,    29,    30,   136,    62,    63,    64,    65,    66,    67,
      68,    37,    45,    14,    41,    48,    45,    18,    85,    33,
      85,    62,    63,    64,    65,    66,    67,    47,    68,    32,
      33,    34,    35,    47,    37,   113,    30,    47,    41,     4,
       5,     6,     7,   121,    23,    24,   113,   125,   113,   116,
     117,   116,   117,    10,   121,   113,   121,    53,   125,    10,
     125,     1,   113,   121,    60,    61,   124,   125,    33,   136,
     121,   136,   113,    45,   125,    15,    10,    17,    34,    35,
     121,    37,    22,   124,   125,    10,    14,    15,    28,    17,
      18,    19,    32,    10,     4,     5,    10,    37,     0,     1,
      46,    89,     4,     5,     6,     7,     8,     9,   125,    11,
      12,    13,   111,    -1,    16,    55,    56,    -1,    20,    -1,
      22,    32,    33,    34,    35,    -1,    37,    -1,    68,    31,
      -1,    33,     4,     5,     6,     7,    -1,     9,    40,    11,
      12,    13,    44,     4,    16,     6,     7,    -1,    20,    26,
      22,     4,     5,     6,     7,    32,    33,    34,    35,    31,
      37,    33,    -1,    -1,     4,     5,     6,     7,    40,    -1,
       4,     5,    44,    45,    -1,    -1,    48,    11,    -1,    -1,
      33,    -1,    16,    -1,    -1,    -1,    39,    -1,    22,    42,
      43,    -1,    -1,    33,    47,     4,     5,     6,     7,    33,
       9,    -1,    11,    12,    13,    45,    40,    16,    26,    -1,
      -1,    20,    -1,    22,    32,    33,    34,    35,    -1,    37,
      -1,    -1,    31,    41,    33,     4,     5,     6,     7,    -1,
       9,    40,    11,    12,    13,    44,    45,    16,    -1,     4,
       5,    20,    -1,    22,     9,    -1,    11,    12,    13,    -1,
      -1,    16,    31,    -1,    33,    20,    -1,    22,    -1,    -1,
      -1,    40,     4,     5,    -1,    44,    31,     9,    33,    11,
      12,    13,    -1,    -1,    16,    40,     4,     5,    20,    44,
      22,     9,    -1,    11,    12,    -1,    -1,    -1,    16,    31,
      -1,    33,     4,     5,    22,    -1,    -1,    -1,    40,    11,
      -1,    -1,    44,    31,    16,    33,    -1,    -1,    -1,    -1,
      22,    -1,    40,    14,    15,    -1,    17,    18,    19,    -1,
      -1,    33,    -1,    -1,    -1,    -1,    -1,    -1,    40
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    50,     0,     1,     4,     5,     6,     7,     9,    11,
      12,    13,    16,    20,    22,    31,    33,    40,    44,    51,
      52,    53,    57,    58,    59,    61,    63,     8,    39,     4,
       7,    57,    25,    64,    65,     4,    16,    40,    58,    40,
      59,    61,    53,    58,    61,    63,     4,     5,    33,    54,
      55,    56,     8,    21,    53,    29,    30,    14,    15,    17,
      18,    19,    26,    32,    33,    34,    35,    37,    40,    58,
      61,    63,    58,    23,    24,    66,    58,    61,    14,    41,
      41,    41,    41,     4,     5,    45,    64,    58,    58,    47,
      54,    55,    57,    69,    70,    72,    73,    39,    42,    43,
      60,    69,    64,    64,    61,    61,    61,    61,    61,    61,
      63,    47,    67,    68,    56,    70,    10,    10,    45,    72,
      69,    68,     4,     5,    33,    45,    54,    55,    57,    58,
      59,    61,    62,    63,    71,    73,    45,    48,    54,    55,
      56,    72,    48,     4,     5,    71,    15,    56,    46
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
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (scanner, YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (YYID (N))                                                     \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (YYID (0))
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])



/* This macro is provided for backward compatibility. */

#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval)
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
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void *scanner)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    void *scanner;
#endif
{
  FILE *yyo gmx_unused = yyoutput;
  YYUSE (yyo);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, void *scanner)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    void *scanner;
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
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
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
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, void *scanner)
#else
static void
yy_reduce_print (yyvsp, yyrule, scanner)
    YYSTYPE *yyvsp;
    int yyrule;
    void *scanner;
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
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       , scanner);
      YYFPRINTF (stderr, "\n");
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

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  YYSIZE_T yysize1;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                if (! (yysize <= yysize1
                       && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                  return 2;
                yysize = yysize1;
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  yysize1 = yysize + yystrlen (yyformat);
  if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
    return 2;
  yysize = yysize1;

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, void *scanner)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, scanner)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    void *scanner;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (scanner);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {
      case 6: /* STR */
/* Line 1381 of yacc.c  */
#line 173 "parser.y"
	{ free(((*yyvaluep).str));        };
/* Line 1381 of yacc.c  */
#line 1380 "parser.cpp"
	break;
      case 7: /* IDENTIFIER */
/* Line 1381 of yacc.c  */
#line 173 "parser.y"
	{ free(((*yyvaluep).str));        };
/* Line 1381 of yacc.c  */
#line 1387 "parser.cpp"
	break;
      case 16: /* KEYWORD_POS */
/* Line 1381 of yacc.c  */
#line 173 "parser.y"
	{ free(((*yyvaluep).str));        };
/* Line 1381 of yacc.c  */
#line 1394 "parser.cpp"
	break;
      case 23: /* PARAM */
/* Line 1381 of yacc.c  */
#line 174 "parser.y"
	{ if(((*yyvaluep).str)) free(((*yyvaluep).str)); };
/* Line 1381 of yacc.c  */
#line 1401 "parser.cpp"
	break;
      case 26: /* CMP_OP */
/* Line 1381 of yacc.c  */
#line 173 "parser.y"
	{ free(((*yyvaluep).str));        };
/* Line 1381 of yacc.c  */
#line 1408 "parser.cpp"
	break;
      case 50: /* commands */
/* Line 1381 of yacc.c  */
#line 175 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1415 "parser.cpp"
	break;
      case 51: /* command */
/* Line 1381 of yacc.c  */
#line 175 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1422 "parser.cpp"
	break;
      case 52: /* cmd_plain */
/* Line 1381 of yacc.c  */
#line 175 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1429 "parser.cpp"
	break;
      case 53: /* selection */
/* Line 1381 of yacc.c  */
#line 175 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1436 "parser.cpp"
	break;
      case 57: /* string */
/* Line 1381 of yacc.c  */
#line 173 "parser.y"
	{ free(((*yyvaluep).str));        };
/* Line 1381 of yacc.c  */
#line 1443 "parser.cpp"
	break;
      case 58: /* sel_expr */
/* Line 1381 of yacc.c  */
#line 176 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1450 "parser.cpp"
	break;
      case 59: /* pos_mod */
/* Line 1381 of yacc.c  */
#line 174 "parser.y"
	{ if(((*yyvaluep).str)) free(((*yyvaluep).str)); };
/* Line 1381 of yacc.c  */
#line 1457 "parser.cpp"
	break;
      case 61: /* num_expr */
/* Line 1381 of yacc.c  */
#line 176 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1464 "parser.cpp"
	break;
      case 62: /* str_expr */
/* Line 1381 of yacc.c  */
#line 176 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1471 "parser.cpp"
	break;
      case 63: /* pos_expr */
/* Line 1381 of yacc.c  */
#line 176 "parser.y"
	{ delete ((*yyvaluep).sel);       };
/* Line 1381 of yacc.c  */
#line 1478 "parser.cpp"
	break;
      case 64: /* method_params */
/* Line 1381 of yacc.c  */
#line 177 "parser.y"
	{ delete ((*yyvaluep).plist);       };
/* Line 1381 of yacc.c  */
#line 1485 "parser.cpp"
	break;
      case 65: /* method_param_list */
/* Line 1381 of yacc.c  */
#line 177 "parser.y"
	{ delete ((*yyvaluep).plist);       };
/* Line 1381 of yacc.c  */
#line 1492 "parser.cpp"
	break;
      case 66: /* method_param */
/* Line 1381 of yacc.c  */
#line 177 "parser.y"
	{ delete ((*yyvaluep).param);       };
/* Line 1381 of yacc.c  */
#line 1499 "parser.cpp"
	break;
      case 67: /* value_list */
/* Line 1381 of yacc.c  */
#line 178 "parser.y"
	{ delete ((*yyvaluep).vlist);       };
/* Line 1381 of yacc.c  */
#line 1506 "parser.cpp"
	break;
      case 68: /* value_list_contents */
/* Line 1381 of yacc.c  */
#line 178 "parser.y"
	{ delete ((*yyvaluep).vlist);       };
/* Line 1381 of yacc.c  */
#line 1513 "parser.cpp"
	break;
      case 69: /* basic_value_list */
/* Line 1381 of yacc.c  */
#line 178 "parser.y"
	{ delete ((*yyvaluep).vlist);       };
/* Line 1381 of yacc.c  */
#line 1520 "parser.cpp"
	break;
      case 70: /* basic_value_list_contents */
/* Line 1381 of yacc.c  */
#line 178 "parser.y"
	{ delete ((*yyvaluep).vlist);       };
/* Line 1381 of yacc.c  */
#line 1527 "parser.cpp"
	break;
      case 71: /* value_item */
/* Line 1381 of yacc.c  */
#line 179 "parser.y"
	{ delete ((*yyvaluep).val);       };
/* Line 1381 of yacc.c  */
#line 1534 "parser.cpp"
	break;
      case 72: /* basic_value_item */
/* Line 1381 of yacc.c  */
#line 179 "parser.y"
	{ delete ((*yyvaluep).val);       };
/* Line 1381 of yacc.c  */
#line 1541 "parser.cpp"
	break;
      case 73: /* value_item_range */
/* Line 1381 of yacc.c  */
#line 179 "parser.y"
	{ delete ((*yyvaluep).val);       };
/* Line 1381 of yacc.c  */
#line 1548 "parser.cpp"
	break;

      default:
	break;
    }
}



struct yypstate
  {
    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;
    /* Used to determine if this is the first time this instance has
       been used.  */
    int yynew;
  };

/* Initialize the parser data structure.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
yypstate *
yypstate_new (void)
#else
yypstate *
yypstate_new ()

#endif
{
  yypstate *yyps;
  yyps = (yypstate *) malloc (sizeof *yyps);
  if (!yyps)
    return YY_NULL;
  yyps->yynew = 1;
  return yyps;
}

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void
yypstate_delete (yypstate *yyps)
#else
void
yypstate_delete (yyps)
    yypstate *yyps;
#endif
{
#ifndef yyoverflow
  /* If the stack was reallocated but the parse did not complete, then the
     stack still needs to be freed.  */
  if (!yyps->yynew && yyps->yyss != yyps->yyssa)
    YYSTACK_FREE (yyps->yyss);
#endif
  free (yyps);
}

#define _gmx_sel_yynerrs yyps->_gmx_sel_yynerrs
#define yystate yyps->yystate
#define yyerrstatus yyps->yyerrstatus
#define yyssa yyps->yyssa
#define yyss yyps->yyss
#define yyssp yyps->yyssp
#define yyvsa yyps->yyvsa
#define yyvs yyps->yyvs
#define yyvsp yyps->yyvsp
#define yystacksize yyps->yystacksize


/*---------------.
| yypush_parse.  |
`---------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yypush_parse (yypstate *yyps, int yypushed_char, YYSTYPE const *yypushed_val, void *scanner)
#else
int
yypush_parse (yyps, yypushed_char, yypushed_val, scanner)
    yypstate *yyps;
    int yypushed_char;
    YYSTYPE const *yypushed_val;
    void *scanner;
#endif
{
/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  if (!yyps->yynew)
    {
      yyn = yypact[yystate];
      goto yyread_pushed_token;
    }

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

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
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
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

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      if (!yyps->yynew)
        {
          YYDPRINTF ((stderr, "Return for a new token:\n"));
          yyresult = YYPUSH_MORE;
          goto yypushreturn;
        }
      yyps->yynew = 0;
yyread_pushed_token:
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yypushed_char;
      if (yypushed_val)
        yylval = *yypushed_val;
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
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
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
/* Line 1787 of yacc.c  */
#line 193 "parser.y"
    {
                 BEGIN_ACTION;
                 set_empty((yyval.sel));
                 END_ACTION;
             }
    break;

  case 3:
/* Line 1787 of yacc.c  */
#line 199 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_append_selection(get((yyvsp[(2) - (2)].sel)), get((yyvsp[(1) - (2)].sel)), scanner));
                 if (_gmx_sel_parser_should_finish(scanner))
                     YYACCEPT;
                 END_ACTION;
             }
    break;

  case 4:
/* Line 1787 of yacc.c  */
#line 209 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (2)].sel); }
    break;

  case 5:
/* Line 1787 of yacc.c  */
#line 211 "parser.y"
    {
                 BEGIN_ACTION;
                 _gmx_selparser_error(scanner, "invalid selection '%s'",
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
                 set_empty((yyval.sel));
                 END_ACTION;
             }
    break;

  case 6:
/* Line 1787 of yacc.c  */
#line 232 "parser.y"
    {
                 BEGIN_ACTION;
                 set_empty((yyval.sel));
                 END_ACTION;
             }
    break;

  case 7:
/* Line 1787 of yacc.c  */
#line 238 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer s
                        = _gmx_sel_init_group_by_id((yyvsp[(1) - (1)].i), scanner);
                 SelectionTreeElementPointer p
                        = _gmx_sel_init_position(s, NULL, scanner);
                 if (!p) YYERROR;
                 set((yyval.sel), _gmx_sel_init_selection(NULL, p, scanner));
                 END_ACTION;
             }
    break;

  case 8:
/* Line 1787 of yacc.c  */
#line 249 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(1) - (1)].str));
                 SelectionTreeElementPointer s
                        = _gmx_sel_init_group_by_name((yyvsp[(1) - (1)].str), scanner);
                 SelectionTreeElementPointer p
                        = _gmx_sel_init_position(s, NULL, scanner);
                 if (!p) YYERROR;
                 set((yyval.sel), _gmx_sel_init_selection(NULL, p, scanner));
                 END_ACTION;
             }
    break;

  case 9:
/* Line 1787 of yacc.c  */
#line 261 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_selection(NULL, get((yyvsp[(1) - (1)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 10:
/* Line 1787 of yacc.c  */
#line 267 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(1) - (2)].str));
                 set((yyval.sel), _gmx_sel_init_selection((yyvsp[(1) - (2)].str), get((yyvsp[(2) - (2)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 11:
/* Line 1787 of yacc.c  */
#line 274 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), get((yyvsp[(3) - (3)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 12:
/* Line 1787 of yacc.c  */
#line 281 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), get((yyvsp[(3) - (3)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 13:
/* Line 1787 of yacc.c  */
#line 288 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), get((yyvsp[(3) - (3)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 14:
/* Line 1787 of yacc.c  */
#line 297 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); }
    break;

  case 15:
/* Line 1787 of yacc.c  */
#line 299 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_position(get((yyvsp[(1) - (1)].sel)), NULL, scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 16:
/* Line 1787 of yacc.c  */
#line 305 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 17:
/* Line 1787 of yacc.c  */
#line 307 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_modifier((yyvsp[(2) - (3)].meth), get((yyvsp[(3) - (3)].plist)), get((yyvsp[(1) - (3)].sel)), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 18:
/* Line 1787 of yacc.c  */
#line 320 "parser.y"
    { (yyval.i) = (yyvsp[(1) - (1)].i); }
    break;

  case 19:
/* Line 1787 of yacc.c  */
#line 321 "parser.y"
    { (yyval.i) = -(yyvsp[(2) - (2)].i); }
    break;

  case 20:
/* Line 1787 of yacc.c  */
#line 325 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); }
    break;

  case 21:
/* Line 1787 of yacc.c  */
#line 326 "parser.y"
    { (yyval.r) = -(yyvsp[(2) - (2)].r); }
    break;

  case 22:
/* Line 1787 of yacc.c  */
#line 329 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].i); }
    break;

  case 23:
/* Line 1787 of yacc.c  */
#line 330 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); }
    break;

  case 24:
/* Line 1787 of yacc.c  */
#line 333 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); }
    break;

  case 25:
/* Line 1787 of yacc.c  */
#line 334 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); }
    break;

  case 26:
/* Line 1787 of yacc.c  */
#line 343 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg(get((yyvsp[(2) - (2)].sel)));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN));
                 sel->u.boolt = BOOL_NOT;
                 sel->child = arg;
                 set((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 27:
/* Line 1787 of yacc.c  */
#line 354 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg1(get((yyvsp[(1) - (3)].sel))), arg2(get((yyvsp[(3) - (3)].sel)));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN));
                 sel->u.boolt = BOOL_AND;
                 sel->child = arg1; sel->child->next = arg2;
                 set((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 28:
/* Line 1787 of yacc.c  */
#line 365 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg1(get((yyvsp[(1) - (3)].sel))), arg2(get((yyvsp[(3) - (3)].sel)));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN));
                 sel->u.boolt = BOOL_OR;
                 sel->child = arg1; sel->child->next = arg2;
                 set((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 29:
/* Line 1787 of yacc.c  */
#line 375 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 30:
/* Line 1787 of yacc.c  */
#line 380 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree opGuard((yyvsp[(2) - (3)].str));
                 set((yyval.sel), _gmx_sel_init_comparison(get((yyvsp[(1) - (3)].sel)), get((yyvsp[(3) - (3)].sel)), (yyvsp[(2) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 31:
/* Line 1787 of yacc.c  */
#line 391 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(2) - (2)].str));
                 set((yyval.sel), _gmx_sel_init_group_by_name((yyvsp[(2) - (2)].str), scanner));
                 END_ACTION;
             }
    break;

  case 32:
/* Line 1787 of yacc.c  */
#line 398 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_group_by_id((yyvsp[(2) - (2)].i), scanner));
                 END_ACTION;
             }
    break;

  case 33:
/* Line 1787 of yacc.c  */
#line 406 "parser.y"
    { (yyval.str) = NULL; }
    break;

  case 34:
/* Line 1787 of yacc.c  */
#line 407 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str);   }
    break;

  case 35:
/* Line 1787 of yacc.c  */
#line 412 "parser.y"
    { (yyval.smt) = gmx::eStringMatchType_RegularExpression; }
    break;

  case 36:
/* Line 1787 of yacc.c  */
#line 413 "parser.y"
    { (yyval.smt) = gmx::eStringMatchType_Wildcard; }
    break;

  case 37:
/* Line 1787 of yacc.c  */
#line 414 "parser.y"
    { (yyval.smt) = gmx::eStringMatchType_Exact; }
    break;

  case 38:
/* Line 1787 of yacc.c  */
#line 419 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (2)].str));
                 set((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), SelectionParserValueListPointer(), (yyvsp[(1) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 39:
/* Line 1787 of yacc.c  */
#line 427 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_init_keyword_strmatch((yyvsp[(2) - (3)].meth), gmx::eStringMatchType_Auto, get((yyvsp[(3) - (3)].vlist)), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 40:
/* Line 1787 of yacc.c  */
#line 435 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (4)].str));
                 set((yyval.sel), _gmx_sel_init_keyword_strmatch((yyvsp[(2) - (4)].meth), (yyvsp[(3) - (4)].smt), get((yyvsp[(4) - (4)].vlist)), (yyvsp[(1) - (4)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 41:
/* Line 1787 of yacc.c  */
#line 443 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (3)].meth), get((yyvsp[(3) - (3)].vlist)), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 42:
/* Line 1787 of yacc.c  */
#line 454 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_init_method((yyvsp[(2) - (3)].meth), get((yyvsp[(3) - (3)].plist)), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 43:
/* Line 1787 of yacc.c  */
#line 469 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST));
                 _gmx_selelem_set_vtype(sel, INT_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.i[0] = (yyvsp[(1) - (1)].i);
                 set((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 44:
/* Line 1787 of yacc.c  */
#line 480 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST));
                 _gmx_selelem_set_vtype(sel, REAL_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.r[0] = (yyvsp[(1) - (1)].r);
                 set((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 45:
/* Line 1787 of yacc.c  */
#line 494 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (2)].str));
                 set((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), SelectionParserValueListPointer(), (yyvsp[(1) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 46:
/* Line 1787 of yacc.c  */
#line 502 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_init_method((yyvsp[(2) - (3)].meth), get((yyvsp[(3) - (3)].plist)), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 47:
/* Line 1787 of yacc.c  */
#line 513 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_arithmetic(get((yyvsp[(1) - (3)].sel)), get((yyvsp[(3) - (3)].sel)), '+', scanner));
                 END_ACTION;
             }
    break;

  case 48:
/* Line 1787 of yacc.c  */
#line 519 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_arithmetic(get((yyvsp[(1) - (3)].sel)), get((yyvsp[(3) - (3)].sel)), '-', scanner));
                 END_ACTION;
             }
    break;

  case 49:
/* Line 1787 of yacc.c  */
#line 525 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_arithmetic(get((yyvsp[(1) - (3)].sel)), get((yyvsp[(3) - (3)].sel)), '*', scanner));
                 END_ACTION;
             }
    break;

  case 50:
/* Line 1787 of yacc.c  */
#line 531 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_arithmetic(get((yyvsp[(1) - (3)].sel)), get((yyvsp[(3) - (3)].sel)), '/', scanner));
                 END_ACTION;
             }
    break;

  case 51:
/* Line 1787 of yacc.c  */
#line 537 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_arithmetic(get((yyvsp[(2) - (2)].sel)), SelectionTreeElementPointer(), '-', scanner));
                 END_ACTION;
             }
    break;

  case 52:
/* Line 1787 of yacc.c  */
#line 543 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_arithmetic(get((yyvsp[(1) - (3)].sel)), get((yyvsp[(3) - (3)].sel)), '^', scanner));
                 END_ACTION;
             }
    break;

  case 53:
/* Line 1787 of yacc.c  */
#line 548 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 54:
/* Line 1787 of yacc.c  */
#line 556 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST));
                 _gmx_selelem_set_vtype(sel, STR_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.s[0] = (yyvsp[(1) - (1)].str);
                 set((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 55:
/* Line 1787 of yacc.c  */
#line 567 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree posmodGuard((yyvsp[(1) - (2)].str));
                 set((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), SelectionParserValueListPointer(), (yyvsp[(1) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 56:
/* Line 1787 of yacc.c  */
#line 582 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_const_position((yyvsp[(2) - (7)].r), (yyvsp[(4) - (7)].r), (yyvsp[(6) - (7)].r)));
                 END_ACTION;
             }
    break;

  case 57:
/* Line 1787 of yacc.c  */
#line 590 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 58:
/* Line 1787 of yacc.c  */
#line 595 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_method((yyvsp[(1) - (2)].meth), get((yyvsp[(2) - (2)].plist)), NULL, scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 59:
/* Line 1787 of yacc.c  */
#line 605 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree keywordGuard((yyvsp[(1) - (3)].str));
                 set((yyval.sel), _gmx_sel_init_position(get((yyvsp[(3) - (3)].sel)), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 60:
/* Line 1787 of yacc.c  */
#line 619 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_variable_ref(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 61:
/* Line 1787 of yacc.c  */
#line 627 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_variable_ref(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 62:
/* Line 1787 of yacc.c  */
#line 635 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.sel), _gmx_sel_init_variable_ref(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 63:
/* Line 1787 of yacc.c  */
#line 648 "parser.y"
    { (yyval.plist) = (yyvsp[(1) - (1)].plist); }
    break;

  case 64:
/* Line 1787 of yacc.c  */
#line 650 "parser.y"
    { (yyval.plist) = (yyvsp[(1) - (2)].plist); }
    break;

  case 65:
/* Line 1787 of yacc.c  */
#line 655 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.plist), SelectionParserParameter::createList());
                 END_ACTION;
             }
    break;

  case 66:
/* Line 1787 of yacc.c  */
#line 661 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionParserParameterListPointer list(get((yyvsp[(1) - (2)].plist)));
                 list->push_back(get((yyvsp[(2) - (2)].param)));
                 set((yyval.plist), move(list));
                 END_ACTION;
             }
    break;

  case 67:
/* Line 1787 of yacc.c  */
#line 672 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree nameGuard((yyvsp[(1) - (2)].str));
                 set((yyval.param), SelectionParserParameter::create((yyvsp[(1) - (2)].str), get((yyvsp[(2) - (2)].vlist))));
                 END_ACTION;
             }
    break;

  case 68:
/* Line 1787 of yacc.c  */
#line 680 "parser.y"
    { (yyval.vlist) = (yyvsp[(1) - (1)].vlist);   }
    break;

  case 69:
/* Line 1787 of yacc.c  */
#line 681 "parser.y"
    { (yyval.vlist) = (yyvsp[(2) - (3)].vlist);   }
    break;

  case 70:
/* Line 1787 of yacc.c  */
#line 686 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.vlist), SelectionParserValue::createList());
                 END_ACTION;
             }
    break;

  case 71:
/* Line 1787 of yacc.c  */
#line 692 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get((yyvsp[(1) - (2)].vlist)));
                 list->push_back(get((yyvsp[(2) - (2)].val)));
                 set((yyval.vlist), move(list));
                 END_ACTION;
             }
    break;

  case 72:
/* Line 1787 of yacc.c  */
#line 700 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get((yyvsp[(1) - (3)].vlist)));
                 list->push_back(get((yyvsp[(3) - (3)].val)));
                 set((yyval.vlist), move(list));
                 END_ACTION;
             }
    break;

  case 73:
/* Line 1787 of yacc.c  */
#line 710 "parser.y"
    { (yyval.vlist) = (yyvsp[(1) - (1)].vlist); }
    break;

  case 74:
/* Line 1787 of yacc.c  */
#line 711 "parser.y"
    { (yyval.vlist) = (yyvsp[(2) - (3)].vlist); }
    break;

  case 75:
/* Line 1787 of yacc.c  */
#line 716 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.vlist), SelectionParserValue::createList(get((yyvsp[(1) - (1)].val))));
                 END_ACTION;
             }
    break;

  case 76:
/* Line 1787 of yacc.c  */
#line 722 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get((yyvsp[(1) - (2)].vlist)));
                 list->push_back(get((yyvsp[(2) - (2)].val)));
                 set((yyval.vlist), move(list));
                 END_ACTION;
             }
    break;

  case 77:
/* Line 1787 of yacc.c  */
#line 730 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get((yyvsp[(1) - (3)].vlist)));
                 list->push_back(get((yyvsp[(3) - (3)].val)));
                 set((yyval.vlist), move(list));
                 END_ACTION;
             }
    break;

  case 78:
/* Line 1787 of yacc.c  */
#line 740 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createExpr(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 79:
/* Line 1787 of yacc.c  */
#line 746 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createExpr(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 80:
/* Line 1787 of yacc.c  */
#line 752 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createExpr(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 81:
/* Line 1787 of yacc.c  */
#line 758 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createExpr(get((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 82:
/* Line 1787 of yacc.c  */
#line 763 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 83:
/* Line 1787 of yacc.c  */
#line 768 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createInteger((yyvsp[(1) - (1)].i)));
                 END_ACTION;
             }
    break;

  case 84:
/* Line 1787 of yacc.c  */
#line 774 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createReal((yyvsp[(1) - (1)].r)));
                 END_ACTION;
             }
    break;

  case 85:
/* Line 1787 of yacc.c  */
#line 780 "parser.y"
    {
                 BEGIN_ACTION;
                 scoped_ptr_sfree stringGuard((yyvsp[(1) - (1)].str));
                 set((yyval.val), SelectionParserValue::createString((yyvsp[(1) - (1)].str)));
                 END_ACTION;
             }
    break;

  case 86:
/* Line 1787 of yacc.c  */
#line 786 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 87:
/* Line 1787 of yacc.c  */
#line 791 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createIntegerRange((yyvsp[(1) - (3)].i), (yyvsp[(3) - (3)].i)));
                 END_ACTION;
             }
    break;

  case 88:
/* Line 1787 of yacc.c  */
#line 797 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createRealRange((yyvsp[(1) - (3)].i), (yyvsp[(3) - (3)].r)));
                 END_ACTION;
             }
    break;

  case 89:
/* Line 1787 of yacc.c  */
#line 803 "parser.y"
    {
                 BEGIN_ACTION;
                 set((yyval.val), SelectionParserValue::createRealRange((yyvsp[(1) - (3)].r), (yyvsp[(3) - (3)].r)));
                 END_ACTION;
             }
    break;


/* Line 1787 of yacc.c  */
#line 2763 "parser.cpp"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
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
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (scanner, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (scanner, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
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

  /* Else will try to reuse lookahead token after shifting the error
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
      if (!yypact_value_is_default (yyn))
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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, scanner);
    }
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
  yyps->yynew = 1;

yypushreturn:
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


