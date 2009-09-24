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
#define yyparse _gmx_sel_yyparse
#define yylex   _gmx_sel_yylex
#define yyerror _gmx_sel_yyerror
#define yylval  _gmx_sel_yylval
#define yychar  _gmx_sel_yychar
#define yydebug _gmx_sel_yydebug
#define yynerrs _gmx_sel_yynerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INVALID = 258,
     HELP = 259,
     HELP_TOPIC = 260,
     INT = 261,
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
     PARAM_INT = 282,
     PARAM_REAL = 283,
     PARAM_STR = 284,
     PARAM_POS = 285,
     PARAM_GROUP = 286,
     END_OF_METHOD = 287,
     XOR = 288,
     OR = 289,
     AND = 290,
     NOT = 291,
     CMP_OP = 292
   };
#endif
/* Tokens.  */
#define INVALID 258
#define HELP 259
#define HELP_TOPIC 260
#define INT 261
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
#define PARAM_INT 282
#define PARAM_REAL 283
#define PARAM_STR 284
#define PARAM_POS 285
#define PARAM_GROUP 286
#define END_OF_METHOD 287
#define XOR 288
#define OR 289
#define AND 290
#define NOT 291
#define CMP_OP 292




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

#include <stdio.h>
#include <string.h>

#include <smalloc.h>
#include <vec.h>

#include <position.h>
#include <selection.h>
#include <selmethod.h>

#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"
#include "selhelp.h"

#include "scanner.h"

static void
yyerror(yyscan_t, int, gmx_ana_indexgrps_t *, char const *s);

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
static t_selexpr_param *
process_param_list(t_selexpr_param *params);


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
#line 71 "parser.y"
{
    int                  i;
    real                 r;
    char                *str;
    gmx_ana_selmethod_t *meth;

    t_selelem        *sel;

    t_selexpr_value  *val;
    t_selexpr_param  *param;
}
/* Line 187 of yacc.c.  */
#line 227 "parser.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 240 "parser.c"

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
# if YYENABLE_NLS
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
#define YYFINAL  44
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   165

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  42
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  21
/* YYNRULES -- Number of rules.  */
#define YYNRULES  70
/* YYNRULES -- Number of states.  */
#define YYNSTATES  115

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   292

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      39,    40,     2,     2,    41,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    38,     2,     2,     2,     2,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     9,    11,    13,    15,    17,    18,
      20,    21,    23,    25,    28,    32,    36,    40,    42,    44,
      46,    49,    52,    54,    58,    62,    65,    69,    73,    77,
      80,    83,    85,    87,    89,    92,    96,   100,   104,   108,
     110,   112,   115,   118,   122,   126,   129,   131,   135,   137,
     139,   143,   145,   147,   155,   159,   161,   164,   166,   169,
     171,   175,   177,   180,   181,   184,   186,   189,   192,   195,
     198
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      43,     0,    -1,    47,    -1,    43,    10,    47,    -1,     6,
      -1,     7,    -1,     8,    -1,     9,    -1,    -1,    20,    -1,
      -1,     6,    -1,    50,    -1,    45,    50,    -1,     9,    38,
      51,    -1,     9,    38,    52,    -1,     9,    38,    55,    -1,
      48,    -1,     4,    -1,    49,    -1,     4,     5,    -1,    49,
       5,    -1,    54,    -1,    50,    25,    60,    -1,    39,    50,
      40,    -1,    36,    51,    -1,    51,    35,    51,    -1,    51,
      34,    51,    -1,    39,    51,    40,    -1,    11,    45,    -1,
      11,     6,    -1,    15,    -1,    14,    -1,    16,    -1,    46,
      21,    -1,    46,    19,    57,    -1,    46,    17,    58,    -1,
      46,    23,    60,    -1,    52,    37,    52,    -1,     6,    -1,
       7,    -1,    46,    17,    -1,    46,    18,    -1,    46,    22,
      60,    -1,    39,    52,    40,    -1,    24,    60,    -1,    53,
      -1,    20,    13,    51,    -1,    51,    -1,    53,    -1,    20,
      13,    51,    -1,    55,    -1,    51,    -1,    39,    44,    41,
      44,    41,    44,    40,    -1,    39,    53,    40,    -1,    45,
      -1,    57,    45,    -1,    59,    -1,    58,    59,    -1,     6,
      -1,     6,    12,     6,    -1,    61,    -1,    61,    32,    -1,
      -1,    61,    62,    -1,    26,    -1,    27,    58,    -1,    28,
      44,    -1,    29,    45,    -1,    30,    56,    -1,    31,    51,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   167,   167,   176,   188,   189,   192,   193,   196,   197,
     201,   202,   213,   214,   215,   217,   219,   221,   226,   232,
     234,   241,   250,   251,   258,   262,   265,   269,   277,   281,
     287,   295,   303,   310,   322,   329,   336,   346,   357,   368,
     372,   376,   383,   390,   402,   407,   419,   420,   427,   438,
     439,   450,   451,   461,   471,   475,   477,   482,   483,   488,
     490,   495,   496,   501,   502,   507,   508,   511,   514,   517,
     520
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INVALID", "HELP", "HELP_TOPIC", "INT",
  "REAL", "STR", "IDENTIFIER", "CMD_SEP", "GROUP", "TO", "OF",
  "VARIABLE_NUMERIC", "VARIABLE_GROUP", "VARIABLE_POS", "KEYWORD_INT",
  "KEYWORD_REAL", "KEYWORD_STR", "KEYWORD_POS", "KEYWORD_GROUP",
  "METHOD_NUMERIC", "METHOD_GROUP", "METHOD_POS", "MODIFIER", "PARAM_BOOL",
  "PARAM_INT", "PARAM_REAL", "PARAM_STR", "PARAM_POS", "PARAM_GROUP",
  "END_OF_METHOD", "XOR", "OR", "AND", "NOT", "CMP_OP", "'='", "'('",
  "')'", "','", "$accept", "commands", "number", "string", "pos_mod",
  "command", "help_request", "help_topic", "selection", "sel_expr",
  "numeric_expr", "pos_expr", "pos_expr_sel", "pos_expr_nosel",
  "pos_expr_nosel_impl", "string_list", "int_list", "int_list_item",
  "method_params", "method_param_list", "method_param", 0
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
     285,   286,   287,   288,   289,   290,   291,   292,    61,    40,
      41,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    42,    43,    43,    44,    44,    45,    45,    46,    46,
      47,    47,    47,    47,    47,    47,    47,    47,    48,    48,
      49,    49,    50,    50,    50,    51,    51,    51,    51,    51,
      51,    51,    52,    53,    51,    51,    51,    51,    51,    52,
      52,    52,    52,    52,    52,    53,    54,    54,    54,    55,
      55,    56,    56,    53,    53,    57,    57,    58,    58,    59,
      59,    60,    60,    61,    61,    62,    62,    62,    62,    62,
      62
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     1,     1,     1,     1,     0,     1,
       0,     1,     1,     2,     3,     3,     3,     1,     1,     1,
       2,     2,     1,     3,     3,     2,     3,     3,     3,     2,
       2,     1,     1,     1,     2,     3,     3,     3,     3,     1,
       1,     2,     2,     3,     3,     2,     1,     3,     1,     1,
       3,     1,     1,     7,     3,     1,     2,     1,     2,     1,
       3,     1,     2,     0,     2,     1,     2,     2,     2,     2,
       2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       8,    18,    11,    40,     6,     7,     0,    32,    31,    33,
       9,    63,     8,     8,     0,     8,     0,     2,    17,    19,
      12,    48,     0,    46,    22,    20,     8,    30,     7,    29,
       8,    45,    61,    39,     9,     8,    25,    39,    40,     0,
       0,    48,     0,    46,     1,     8,    13,    41,    42,     0,
      34,    63,    63,    21,    63,     8,     8,     8,     9,     8,
      14,    15,    49,    16,    47,    65,     0,     0,     0,     8,
       8,    62,    64,     0,     0,    24,    28,    44,    54,     3,
      59,    36,    57,    55,    35,    43,    37,    23,    27,    26,
       8,     0,    38,     8,     0,    66,     4,     5,    67,    68,
      52,    51,    69,    70,     0,     0,    58,    56,     0,    41,
      50,     0,    60,     0,    53
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,    14,    39,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    63,   102,    84,    81,    82,    31,    32,
      72
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -69
static const yytype_int16 yypact[] =
{
      62,    20,   -17,   -69,   -69,     3,   138,   -69,   -69,   -69,
      21,   -69,   114,    15,     7,    88,   137,   -69,   -69,    48,
      36,    40,    30,   -69,   -69,   -69,    99,   -69,   -69,   -69,
     114,   -69,    61,   -69,   -69,   114,   -69,    22,    59,    66,
     -16,     2,   -25,    39,   -69,    62,    36,   103,   -69,    76,
     -69,   -69,   -69,   -69,   -69,   114,   114,    26,    98,   125,
      40,    30,   -69,   -69,    40,   -69,   103,    90,    76,    99,
     114,   -69,   -69,     2,    90,   -69,   -69,   -69,   -69,   -69,
     104,   103,   -69,   -69,    76,   -69,   -69,   -69,   -69,   -69,
      26,    38,   -69,   114,    39,   103,   -69,   -69,   -69,   -69,
      40,   -69,   -69,    40,    77,   116,   -69,   -69,    86,   -69,
      40,    90,   -69,    93,   -69
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -69,   -69,   -63,    -4,   -52,    72,   -69,   -69,    -5,   -12,
      -7,   -10,   -69,    68,   -69,   -69,    64,   -68,   111,   -69,
     -69
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -40
static const yytype_int8 yytable[] =
{
      36,    41,    29,    43,    98,    91,    42,    44,    40,    54,
      46,   104,    57,   106,    60,    77,    62,    45,    64,    61,
     -39,    37,    38,    73,    75,    25,     6,   106,    42,     7,
       8,     9,    33,     3,    30,    10,    55,    56,    91,    11,
       7,    26,    76,    88,    89,    83,    34,    73,   113,    94,
      92,    12,    42,    53,    13,   109,    48,   100,   103,    62,
      51,    54,   -10,    -4,    99,    90,     1,    57,     2,     3,
       4,     5,   -10,     6,    55,    56,     7,     8,     9,    78,
     107,   110,    10,   108,     4,    28,    11,    65,    66,    67,
      68,    69,    70,    71,    33,     3,    96,    97,    12,     6,
      -5,    13,     7,     8,     9,    33,     3,    74,    10,    80,
       6,    93,    11,     7,     8,     9,   105,    79,   111,    58,
      33,     3,   112,    11,    12,     6,    77,    13,     7,     8,
      95,    37,    38,   114,    34,    12,     6,   101,    59,     7,
       8,     9,     0,     0,    27,    34,     4,    28,     0,    11,
      12,     0,     0,    35,    47,    48,    49,     0,    50,    51,
      52,    12,    85,    86,    59,    87
};

static const yytype_int8 yycheck[] =
{
      12,    13,     6,    13,    67,    57,    13,     0,    13,    25,
      15,    74,    37,    81,    26,    40,    26,    10,    30,    26,
      37,     6,     7,    35,    40,     5,    11,    95,    35,    14,
      15,    16,     6,     7,    13,    20,    34,    35,    90,    24,
      14,    38,    40,    55,    56,    49,    20,    59,   111,    59,
      57,    36,    59,     5,    39,    17,    18,    69,    70,    69,
      22,    25,     0,    41,    68,    39,     4,    37,     6,     7,
       8,     9,    10,    11,    34,    35,    14,    15,    16,    40,
      84,    93,    20,    90,     8,     9,    24,    26,    27,    28,
      29,    30,    31,    32,     6,     7,     6,     7,    36,    11,
      41,    39,    14,    15,    16,     6,     7,    41,    20,     6,
      11,    13,    24,    14,    15,    16,    12,    45,    41,    20,
       6,     7,     6,    24,    36,    11,    40,    39,    14,    15,
      66,     6,     7,    40,    20,    36,    11,    69,    39,    14,
      15,    16,    -1,    -1,     6,    20,     8,     9,    -1,    24,
      36,    -1,    -1,    39,    17,    18,    19,    -1,    21,    22,
      23,    36,    51,    52,    39,    54
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,     6,     7,     8,     9,    11,    14,    15,    16,
      20,    24,    36,    39,    43,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,     5,    38,     6,     9,    45,
      13,    60,    61,     6,    20,    39,    51,     6,     7,    44,
      50,    51,    52,    53,     0,    10,    50,    17,    18,    19,
      21,    22,    23,     5,    25,    34,    35,    37,    20,    39,
      51,    52,    53,    55,    51,    26,    27,    28,    29,    30,
      31,    32,    62,    51,    41,    40,    40,    40,    40,    47,
       6,    58,    59,    45,    57,    60,    60,    60,    51,    51,
      39,    46,    52,    13,    53,    58,     6,     7,    44,    45,
      51,    55,    56,    51,    44,    12,    59,    45,    52,    17,
      51,    41,     6,    44,    40
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
      yyerror (scanner, nexp, grps, YY_("syntax error: cannot back up")); \
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
# if YYLTYPE_IS_TRIVIAL
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
		  Type, Value, scanner, nexp, grps); \
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
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, yyscan_t                 scanner, int                      nexp, gmx_ana_indexgrps_t     *grps)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner, nexp, grps)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    yyscan_t                 scanner;
    int                      nexp;
    gmx_ana_indexgrps_t     *grps;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (scanner);
  YYUSE (nexp);
  YYUSE (grps);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, yyscan_t                 scanner, int                      nexp, gmx_ana_indexgrps_t     *grps)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, scanner, nexp, grps)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    yyscan_t                 scanner;
    int                      nexp;
    gmx_ana_indexgrps_t     *grps;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, scanner, nexp, grps);
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
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, yyscan_t                 scanner, int                      nexp, gmx_ana_indexgrps_t     *grps)
#else
static void
yy_reduce_print (yyvsp, yyrule, scanner, nexp, grps)
    YYSTYPE *yyvsp;
    int yyrule;
    yyscan_t                 scanner;
    int                      nexp;
    gmx_ana_indexgrps_t     *grps;
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
		       		       , scanner, nexp, grps);
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, scanner, nexp, grps); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, yyscan_t                 scanner, int                      nexp, gmx_ana_indexgrps_t     *grps)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, scanner, nexp, grps)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    yyscan_t                 scanner;
    int                      nexp;
    gmx_ana_indexgrps_t     *grps;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (scanner);
  YYUSE (nexp);
  YYUSE (grps);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {
      case 5: /* "HELP_TOPIC" */
#line 145 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1267 "parser.c"
	break;
      case 8: /* "STR" */
#line 145 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1272 "parser.c"
	break;
      case 9: /* "IDENTIFIER" */
#line 145 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1277 "parser.c"
	break;
      case 45: /* "string" */
#line 145 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1282 "parser.c"
	break;
      case 47: /* "command" */
#line 146 "parser.y"
	{ if((yyvaluep->sel)) _gmx_selelem_free((yyvaluep->sel)); };
#line 1287 "parser.c"
	break;
      case 50: /* "selection" */
#line 147 "parser.y"
	{ _gmx_selelem_free_chain((yyvaluep->sel));  };
#line 1292 "parser.c"
	break;
      case 51: /* "sel_expr" */
#line 148 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1297 "parser.c"
	break;
      case 52: /* "numeric_expr" */
#line 148 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1302 "parser.c"
	break;
      case 53: /* "pos_expr" */
#line 149 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1307 "parser.c"
	break;
      case 54: /* "pos_expr_sel" */
#line 149 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1312 "parser.c"
	break;
      case 55: /* "pos_expr_nosel" */
#line 149 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1317 "parser.c"
	break;
      case 56: /* "pos_expr_nosel_impl" */
#line 149 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1322 "parser.c"
	break;
      case 57: /* "string_list" */
#line 150 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1327 "parser.c"
	break;
      case 58: /* "int_list" */
#line 150 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1332 "parser.c"
	break;
      case 59: /* "int_list_item" */
#line 150 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1337 "parser.c"
	break;
      case 60: /* "method_params" */
#line 151 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1342 "parser.c"
	break;
      case 61: /* "method_param_list" */
#line 151 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1347 "parser.c"
	break;
      case 62: /* "method_param" */
#line 151 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1352 "parser.c"
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
int yyparse (yyscan_t                 scanner, int                      nexp, gmx_ana_indexgrps_t     *grps);
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
yyparse (yyscan_t                 scanner, int                      nexp, gmx_ana_indexgrps_t     *grps)
#else
int
yyparse (scanner, nexp, grps)
    yyscan_t                 scanner;
    int                      nexp;
    gmx_ana_indexgrps_t     *grps;
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
#line 168 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_append_selection(sc, (yyvsp[(1) - (1)].sel), NULL);
                 if (sc->nr == nexp)
                     YYACCEPT;
                 _gmx_sel_lexer_clear_pselstr(scanner);
             ;}
    break;

  case 3:
#line 177 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_append_selection(sc, (yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].sel));
                 if (sc->nr == nexp)
                     YYACCEPT;
                 _gmx_sel_lexer_clear_pselstr(scanner);
             ;}
    break;

  case 4:
#line 188 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].i); ;}
    break;

  case 5:
#line 189 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); ;}
    break;

  case 6:
#line 192 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); ;}
    break;

  case 7:
#line 193 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); ;}
    break;

  case 8:
#line 196 "parser.y"
    { (yyval.str) = NULL; ;}
    break;

  case 9:
#line 197 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str);   ;}
    break;

  case 10:
#line 201 "parser.y"
    { (yyval.sel) = NULL;                            ;}
    break;

  case 11:
#line 203 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 t_selelem               *s, *p;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 s = _gmx_sel_init_group_by_id(grps, (yyvsp[(1) - (1)].i));
                 if (s == NULL) YYABORT;
                 p = _gmx_sel_init_position(sc, s, sc->spost, TRUE);
                 if (p == NULL) YYABORT;
                 (yyval.sel) = _gmx_sel_init_selection(strdup(s->name), p, scanner);
             ;}
    break;

  case 12:
#line 213 "parser.y"
    { (yyval.sel) = _gmx_sel_init_selection(NULL, (yyvsp[(1) - (1)].sel), scanner); ;}
    break;

  case 13:
#line 214 "parser.y"
    { (yyval.sel) = _gmx_sel_init_selection((yyvsp[(1) - (2)].str), (yyvsp[(2) - (2)].sel), scanner);  ;}
    break;

  case 14:
#line 216 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel), scanner); ;}
    break;

  case 15:
#line 218 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel), scanner); ;}
    break;

  case 16:
#line 220 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel), scanner); ;}
    break;

  case 17:
#line 221 "parser.y"
    { (yyval.sel) = NULL; ;}
    break;

  case 18:
#line 227 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 _gmx_sel_print_help(sc, NULL);
             ;}
    break;

  case 20:
#line 235 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 _gmx_sel_print_help(sc, (yyvsp[(2) - (2)].str));
                 sfree((yyvsp[(2) - (2)].str));
             ;}
    break;

  case 21:
#line 242 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 _gmx_sel_print_help(sc, (yyvsp[(2) - (2)].str));
                 sfree((yyvsp[(2) - (2)].str));
             ;}
    break;

  case 22:
#line 250 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 23:
#line 252 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_modifier(sc, (yyvsp[(2) - (3)].meth), process_param_list((yyvsp[(3) - (3)].param)), (yyvsp[(1) - (3)].sel));
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 24:
#line 258 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                              ;}
    break;

  case 25:
#line 262 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                                  (yyval.sel)->u.boolt = BOOL_NOT;
                                  (yyval.sel)->child = (yyvsp[(2) - (2)].sel);                       ;}
    break;

  case 26:
#line 266 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                                  (yyval.sel)->u.boolt = BOOL_AND;
                                  (yyval.sel)->child = (yyvsp[(1) - (3)].sel); (yyval.sel)->child->next = (yyvsp[(3) - (3)].sel); ;}
    break;

  case 27:
#line 270 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                                  (yyval.sel)->u.boolt = BOOL_OR;
                                  (yyval.sel)->child = (yyvsp[(1) - (3)].sel); (yyval.sel)->child->next = (yyvsp[(3) - (3)].sel); ;}
    break;

  case 28:
#line 277 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                               ;}
    break;

  case 29:
#line 282 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_group_by_name(grps, (yyvsp[(2) - (2)].str));
                 sfree((yyvsp[(2) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 30:
#line 288 "parser.y"
    {
                 (yyval.sel) = _gmx_sel_init_group_by_id(grps, (yyvsp[(2) - (2)].i));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 31:
#line 295 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  _gmx_selelem_set_vtype((yyval.sel), (yyvsp[(1) - (1)].sel)->v.type);
                                  (yyval.sel)->name   = (yyvsp[(1) - (1)].sel)->name;
                                  (yyval.sel)->child  = (yyvsp[(1) - (1)].sel);
                                  (yyvsp[(1) - (1)].sel)->refcount++;                        ;}
    break;

  case 32:
#line 303 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  _gmx_selelem_set_vtype((yyval.sel), (yyvsp[(1) - (1)].sel)->v.type);
                                  (yyval.sel)->name   = (yyvsp[(1) - (1)].sel)->name;
                                  (yyval.sel)->child  = (yyvsp[(1) - (1)].sel);
                                  (yyvsp[(1) - (1)].sel)->refcount++;                        ;}
    break;

  case 33:
#line 310 "parser.y"
    { if ((yyvsp[(1) - (1)].sel)->type == SEL_CONST) {
                                      (yyval.sel) = (yyvsp[(1) - (1)].sel);
                                  } else {
                                      (yyval.sel) = _gmx_selelem_create(SEL_SUBEXPRREF);
                                      _gmx_selelem_set_vtype((yyval.sel), (yyvsp[(1) - (1)].sel)->v.type);
                                      (yyval.sel)->name   = (yyvsp[(1) - (1)].sel)->name;
                                      (yyval.sel)->child  = (yyvsp[(1) - (1)].sel);
                                  }
                                  (yyvsp[(1) - (1)].sel)->refcount++;                        ;}
    break;

  case 34:
#line 323 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_keyword(sc, (yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 35:
#line 330 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_keyword(sc, (yyvsp[(2) - (3)].meth), process_value_list((yyvsp[(3) - (3)].val), NULL), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 36:
#line 337 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_keyword(sc, (yyvsp[(2) - (3)].meth), process_value_list((yyvsp[(3) - (3)].val), NULL), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 37:
#line 347 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_method(sc, (yyvsp[(2) - (3)].meth), process_param_list((yyvsp[(3) - (3)].param)), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 38:
#line 358 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_comparison(sc, (yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), (yyvsp[(2) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 39:
#line 368 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype((yyval.sel), INT_VALUE);
                                  _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                                  (yyval.sel)->v.u.i[0] = (yyvsp[(1) - (1)].i);                  ;}
    break;

  case 40:
#line 372 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype((yyval.sel), REAL_VALUE);
                                  _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                                  (yyval.sel)->v.u.r[0] = (yyvsp[(1) - (1)].r);                  ;}
    break;

  case 41:
#line 377 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_keyword(sc, (yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 42:
#line 384 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_keyword(sc, (yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 43:
#line 391 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_method(sc, (yyvsp[(2) - (3)].meth), process_param_list((yyvsp[(3) - (3)].param)), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 44:
#line 403 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                               ;}
    break;

  case 45:
#line 408 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_method(sc, (yyvsp[(1) - (2)].meth), process_param_list((yyvsp[(2) - (2)].param)), NULL);
                 if ((yyval.sel) == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 46:
#line 419 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 47:
#line 421 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].str), TRUE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 48:
#line 428 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(1) - (1)].sel), sc->spost, TRUE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 49:
#line 438 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 50:
#line 440 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].str), FALSE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 51:
#line 450 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 52:
#line 452 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(1) - (1)].sel), sc->rpost, FALSE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 53:
#line 462 "parser.y"
    { rvec x;
                                  (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype((yyval.sel), POS_VALUE);
                                  _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                                  x[XX] = (yyvsp[(2) - (7)].r); x[YY] = (yyvsp[(4) - (7)].r); x[ZZ] = (yyvsp[(6) - (7)].r);
                                  gmx_ana_pos_init_const((yyval.sel)->v.u.p, x);  ;}
    break;

  case 54:
#line 471 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                               ;}
    break;

  case 55:
#line 475 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                                  (yyval.val)->u.s = (yyvsp[(1) - (1)].str);                          ;}
    break;

  case 56:
#line 477 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                                  (yyval.val)->u.s = (yyvsp[(2) - (2)].str); (yyval.val)->next = (yyvsp[(1) - (2)].val);           ;}
    break;

  case 57:
#line 482 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val);                               ;}
    break;

  case 58:
#line 484 "parser.y"
    { (yyvsp[(2) - (2)].val)->next = (yyvsp[(1) - (2)].val); (yyval.val) = (yyvsp[(2) - (2)].val);                ;}
    break;

  case 59:
#line 488 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                                  (yyval.val)->u.i.i1 = (yyval.val)->u.i.i2 = (yyvsp[(1) - (1)].i);          ;}
    break;

  case 60:
#line 490 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                                  (yyval.val)->u.i.i1 = (yyvsp[(1) - (3)].i); (yyval.val)->u.i.i2 = (yyvsp[(3) - (3)].i);      ;}
    break;

  case 61:
#line 495 "parser.y"
    { (yyval.param) = (yyvsp[(1) - (1)].param); ;}
    break;

  case 62:
#line 497 "parser.y"
    { (yyval.param) = (yyvsp[(1) - (2)].param); ;}
    break;

  case 63:
#line 501 "parser.y"
    { (yyval.param) = NULL;              ;}
    break;

  case 64:
#line 503 "parser.y"
    { (yyvsp[(2) - (2)].param)->next = (yyvsp[(1) - (2)].param); (yyval.param) = (yyvsp[(2) - (2)].param); ;}
    break;

  case 65:
#line 507 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (1)].str));    ;}
    break;

  case 66:
#line 509 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = (yyvsp[(2) - (2)].val);                        ;}
    break;

  case 67:
#line 511 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value(REAL_VALUE);
                                  (yyval.param)->value->u.r = (yyvsp[(2) - (2)].r);                   ;}
    break;

  case 68:
#line 514 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value(STR_VALUE);
                                  (yyval.param)->value->u.s = (yyvsp[(2) - (2)].str);                   ;}
    break;

  case 69:
#line 518 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value_expr((yyvsp[(2) - (2)].sel)); ;}
    break;

  case 70:
#line 521 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value_expr((yyvsp[(2) - (2)].sel)); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2165 "parser.c"
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
      yyerror (scanner, nexp, grps, YY_("syntax error"));
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
	    yyerror (scanner, nexp, grps, yymsg);
	  }
	else
	  {
	    yyerror (scanner, nexp, grps, YY_("syntax error"));
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
		      yytoken, &yylval, scanner, nexp, grps);
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
		  yystos[yystate], yyvsp, scanner, nexp, grps);
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
  yyerror (scanner, nexp, grps, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, scanner, nexp, grps);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, scanner, nexp, grps);
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


#line 525 "parser.y"


/*!
 * \param[in,out] scanner Scanner data structure.
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     grps  External index groups (can be NULL).
 * \param[in]     maxnr Maximum number of selections to parse
 *   (if -1, parse as many as provided by the user).
 * \returns       0 on success, -1 on error.
 */
int
_gmx_sel_run_parser(yyscan_t scanner, gmx_ana_selcollection_t *sc,
                    gmx_ana_indexgrps_t *grps, int maxnr)
{
    bool bOk;
    int  nr;
    int  nexp;

    nr        = sc->nr;
    nexp      = (maxnr > 0) ? (sc->nr + maxnr) : -1;
    bOk = !_gmx_sel_yyparse(scanner, nexp, grps);
    _gmx_sel_free_lexer(scanner);
    if (sc->selstr)
    {
        srenew(sc->selstr, strlen(sc->selstr) + 1);
    }
    nr = sc->nr - nr;
    if (maxnr > 0 && nr != maxnr)
    {
        return -1;
    }
    return bOk ? 0 : -1;
}

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

    /* Reverse list and process values */
    ppar = NULL;
    par  = params;
    while (par)
    {
        par->value = process_value_list(par->value, &par->nval);

        npar = par->next;
        par->next = ppar;
        ppar = par;
        par = npar;
    }
    params = ppar;

    return params;
}

static void
yyerror(yyscan_t scanner, int nexp, gmx_ana_indexgrps_t *grps,
        char const *s)
{
    _gmx_selparser_error("%s", s);
}

