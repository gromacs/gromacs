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

#include <indexutil.h>
#include <position.h>
#include <selection.h>
#include <selmethod.h>

#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"

#include "scanner.h"

static void
yyerror(yyscan_t, int, gmx_ana_indexgrps_t *, char const *s);

static t_selelem *
get_group_by_name(gmx_ana_indexgrps_t *grps, char *name);
static t_selelem *
get_group_by_id(gmx_ana_indexgrps_t *grps, int id);

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
static t_selexpr_param *
process_param_list(t_selexpr_param *params);
static t_selelem *
init_keyword_expr(gmx_ana_selcollection_t *sc, gmx_ana_selmethod_t *method,
                  t_selexpr_value *values, char *rpost);


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
#line 79 "parser.y"
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
#line 231 "parser.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 244 "parser.c"

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
#define YYFINAL  40
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   189

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  40
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  65
/* YYNRULES -- Number of states.  */
#define YYNSTATES  110

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
      37,    38,     2,     2,    39,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    36,     2,     2,     2,     2,     2,     2,     2,     2,
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
      35
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     9,    11,    13,    15,    17,    18,
      20,    21,    23,    25,    28,    32,    36,    40,    42,    46,
      50,    53,    57,    61,    65,    68,    71,    73,    75,    77,
      80,    84,    88,    92,    96,    98,   100,   103,   106,   110,
     114,   117,   119,   121,   125,   127,   131,   133,   135,   143,
     147,   149,   152,   154,   157,   159,   163,   165,   168,   169,
     172,   174,   177,   180,   183,   186
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      41,     0,    -1,    45,    -1,    41,     8,    45,    -1,     4,
      -1,     5,    -1,     6,    -1,     7,    -1,    -1,    18,    -1,
      -1,     4,    -1,    46,    -1,    43,    46,    -1,     7,    36,
      47,    -1,     7,    36,    48,    -1,     7,    36,    51,    -1,
      50,    -1,    46,    23,    56,    -1,    37,    46,    38,    -1,
      34,    47,    -1,    47,    33,    47,    -1,    47,    32,    47,
      -1,    37,    47,    38,    -1,     9,    43,    -1,     9,     4,
      -1,    13,    -1,    12,    -1,    14,    -1,    44,    19,    -1,
      44,    17,    53,    -1,    44,    15,    54,    -1,    44,    21,
      56,    -1,    48,    35,    48,    -1,     4,    -1,     5,    -1,
      44,    15,    -1,    44,    16,    -1,    44,    20,    56,    -1,
      37,    48,    38,    -1,    22,    56,    -1,    49,    -1,    47,
      -1,    18,    11,    47,    -1,    49,    -1,    18,    11,    47,
      -1,    51,    -1,    47,    -1,    37,    42,    39,    42,    39,
      42,    38,    -1,    37,    49,    38,    -1,    43,    -1,    53,
      43,    -1,    55,    -1,    54,    55,    -1,     4,    -1,     4,
      10,     4,    -1,    57,    -1,    57,    30,    -1,    -1,    57,
      58,    -1,    24,    -1,    25,    54,    -1,    26,    42,    -1,
      27,    43,    -1,    28,    52,    -1,    29,    47,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   171,   171,   179,   190,   191,   194,   195,   198,   199,
     203,   204,   215,   216,   218,   220,   222,   227,   228,   235,
     239,   242,   246,   254,   258,   261,   266,   274,   281,   293,
     300,   307,   317,   328,   339,   343,   347,   354,   361,   373,
     378,   390,   391,   398,   409,   410,   421,   422,   432,   442,
     446,   448,   453,   454,   459,   461,   466,   467,   472,   473,
     478,   479,   482,   485,   488,   491
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INVALID", "INT", "REAL", "STR",
  "IDENTIFIER", "CMD_SEP", "GROUP", "TO", "OF", "VARIABLE_NUMERIC",
  "VARIABLE_GROUP", "VARIABLE_POS", "KEYWORD_INT", "KEYWORD_REAL",
  "KEYWORD_STR", "KEYWORD_POS", "KEYWORD_GROUP", "METHOD_NUMERIC",
  "METHOD_GROUP", "METHOD_POS", "MODIFIER", "PARAM_BOOL", "PARAM_INT",
  "PARAM_REAL", "PARAM_STR", "PARAM_POS", "PARAM_GROUP", "END_OF_METHOD",
  "XOR", "OR", "AND", "NOT", "CMP_OP", "'='", "'('", "')'", "','",
  "$accept", "commands", "number", "string", "pos_mod", "command",
  "selection", "sel_expr", "numeric_expr", "pos_expr", "pos_expr_sel",
  "pos_expr_nosel", "pos_expr_nosel_impl", "string_list", "int_list",
  "int_list_item", "method_params", "method_param_list", "method_param", 0
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
     285,   286,   287,   288,   289,   290,    61,    40,    41,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    40,    41,    41,    42,    42,    43,    43,    44,    44,
      45,    45,    45,    45,    45,    45,    45,    46,    46,    46,
      47,    47,    47,    47,    47,    47,    47,    48,    49,    47,
      47,    47,    47,    47,    48,    48,    48,    48,    48,    48,
      49,    50,    50,    50,    51,    51,    52,    52,    49,    49,
      53,    53,    54,    54,    55,    55,    56,    56,    57,    57,
      58,    58,    58,    58,    58,    58
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     1,     1,     1,     1,     0,     1,
       0,     1,     1,     2,     3,     3,     3,     1,     3,     3,
       2,     3,     3,     3,     2,     2,     1,     1,     1,     2,
       3,     3,     3,     3,     1,     1,     2,     2,     3,     3,
       2,     1,     1,     3,     1,     3,     1,     1,     7,     3,
       1,     2,     1,     2,     1,     3,     1,     2,     0,     2,
       1,     2,     2,     2,     2,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       8,    11,    35,     6,     7,     0,    27,    26,    28,     9,
      58,     8,     8,     0,     8,     0,     2,    12,    42,     0,
      41,    17,     8,    25,     7,    24,     8,    40,    56,    34,
       9,     8,    20,    34,    35,     0,     0,    42,     0,    41,
       1,     8,    13,    36,    37,     0,    29,    58,    58,    58,
       8,     8,     8,     9,     8,    14,    15,    44,    16,    43,
      60,     0,     0,     0,     8,     8,    57,    59,     0,     0,
      19,    23,    39,    49,     3,    54,    31,    52,    50,    30,
      38,    32,    18,    22,    21,     8,     0,    33,     8,     0,
      61,     4,     5,    62,    63,    47,    46,    64,    65,     0,
       0,    53,    51,     0,    36,    45,     0,    55,     0,    48
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,    13,    35,    14,    15,    16,    17,    18,    19,    20,
      21,    58,    97,    79,    76,    77,    27,    28,    67
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -68
static const yytype_int16 yypact[] =
{
      78,    -1,   -68,   -68,   -12,    52,   -68,   -68,   -68,    25,
     -68,   141,    89,    21,   115,    53,   -68,    18,    48,     3,
     -68,   -68,   126,   -68,   -68,   -68,   141,   -68,    37,   -68,
     -68,   141,   -68,     5,     8,    11,   -16,    -5,   -21,    19,
     -68,    78,    18,    67,   -68,    82,   -68,   -68,   -68,   -68,
     141,   141,    14,    41,   152,    48,     3,   -68,   -68,    48,
     -68,    67,   104,    82,   126,   141,   -68,   -68,    -5,   104,
     -68,   -68,   -68,   -68,   -68,    65,    67,   -68,   -68,    82,
     -68,   -68,   -68,   -68,   -68,    14,    15,   -68,   141,    19,
      67,   -68,   -68,   -68,   -68,    48,   -68,   -68,    48,    56,
      93,   -68,   -68,    61,   -68,    48,   104,   -68,    72,   -68
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -68,   -68,   -57,    -3,   -48,    73,    -4,   -11,    -6,    -9,
     -68,    49,   -68,   -68,    55,   -67,    57,   -68,   -68
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -35
static const yytype_int8 yytable[] =
{
      32,    37,    25,    39,    86,    93,    38,    49,    36,   101,
      42,    55,    99,    57,    52,    59,    56,    72,    29,     2,
      68,    40,    70,   101,    22,    38,     6,    50,    51,    41,
     104,    44,    30,    71,   -34,    47,    26,    86,    52,    83,
      84,    49,    78,    68,    -4,    89,    87,    -5,    38,   108,
      69,    85,    88,    95,    98,    57,    23,    73,     3,    24,
      94,    60,    61,    62,    63,    64,    65,    66,    43,    44,
      45,    75,    46,    47,    48,   100,   102,   105,   -10,   103,
      50,    51,     1,     2,     3,     4,   -10,     5,     3,    24,
       6,     7,     8,    33,    34,   106,     9,   107,     5,    72,
      10,     6,     7,     8,    80,    81,    82,     9,    91,    92,
     109,    10,    11,    96,    74,    12,    90,     0,     0,    29,
       2,     0,     0,    11,     5,     0,    12,     6,     7,     8,
      29,     2,     0,     9,     0,     5,     0,    10,     6,     7,
       8,     0,     0,     0,    53,    29,     2,     0,    10,    11,
       5,     0,    12,     6,     7,     0,    33,    34,     0,    30,
      11,     5,     0,    54,     6,     7,     8,     0,     0,     0,
      30,     0,     0,     0,    10,    11,     0,     0,    31,     0,
       0,     0,     0,     0,     0,     0,    11,     0,     0,    54
};

static const yytype_int8 yycheck[] =
{
      11,    12,     5,    12,    52,    62,    12,    23,    12,    76,
      14,    22,    69,    22,    35,    26,    22,    38,     4,     5,
      31,     0,    38,    90,    36,    31,    12,    32,    33,     8,
      15,    16,    18,    38,    35,    20,    11,    85,    35,    50,
      51,    23,    45,    54,    39,    54,    52,    39,    54,   106,
      39,    37,    11,    64,    65,    64,     4,    38,     6,     7,
      63,    24,    25,    26,    27,    28,    29,    30,    15,    16,
      17,     4,    19,    20,    21,    10,    79,    88,     0,    85,
      32,    33,     4,     5,     6,     7,     8,     9,     6,     7,
      12,    13,    14,     4,     5,    39,    18,     4,     9,    38,
      22,    12,    13,    14,    47,    48,    49,    18,     4,     5,
      38,    22,    34,    64,    41,    37,    61,    -1,    -1,     4,
       5,    -1,    -1,    34,     9,    -1,    37,    12,    13,    14,
       4,     5,    -1,    18,    -1,     9,    -1,    22,    12,    13,
      14,    -1,    -1,    -1,    18,     4,     5,    -1,    22,    34,
       9,    -1,    37,    12,    13,    -1,     4,     5,    -1,    18,
      34,     9,    -1,    37,    12,    13,    14,    -1,    -1,    -1,
      18,    -1,    -1,    -1,    22,    34,    -1,    -1,    37,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    34,    -1,    -1,    37
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,     5,     6,     7,     9,    12,    13,    14,    18,
      22,    34,    37,    41,    43,    44,    45,    46,    47,    48,
      49,    50,    36,     4,     7,    43,    11,    56,    57,     4,
      18,    37,    47,     4,     5,    42,    46,    47,    48,    49,
       0,     8,    46,    15,    16,    17,    19,    20,    21,    23,
      32,    33,    35,    18,    37,    47,    48,    49,    51,    47,
      24,    25,    26,    27,    28,    29,    30,    58,    47,    39,
      38,    38,    38,    38,    45,     4,    54,    55,    43,    53,
      56,    56,    56,    47,    47,    37,    44,    48,    11,    49,
      54,     4,     5,    42,    43,    47,    51,    52,    47,    42,
      10,    55,    43,    48,    15,    47,    39,     4,    42,    38
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
      case 6: /* "STR" */
#line 149 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1262 "parser.c"
	break;
      case 7: /* "IDENTIFIER" */
#line 149 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1267 "parser.c"
	break;
      case 43: /* "string" */
#line 149 "parser.y"
	{ free((yyvaluep->str));                     };
#line 1272 "parser.c"
	break;
      case 45: /* "command" */
#line 150 "parser.y"
	{ if((yyvaluep->sel)) _gmx_selelem_free((yyvaluep->sel)); };
#line 1277 "parser.c"
	break;
      case 46: /* "selection" */
#line 151 "parser.y"
	{ _gmx_selelem_free_chain((yyvaluep->sel));  };
#line 1282 "parser.c"
	break;
      case 47: /* "sel_expr" */
#line 152 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1287 "parser.c"
	break;
      case 48: /* "numeric_expr" */
#line 152 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1292 "parser.c"
	break;
      case 49: /* "pos_expr" */
#line 153 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1297 "parser.c"
	break;
      case 50: /* "pos_expr_sel" */
#line 153 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1302 "parser.c"
	break;
      case 51: /* "pos_expr_nosel" */
#line 153 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1307 "parser.c"
	break;
      case 52: /* "pos_expr_nosel_impl" */
#line 153 "parser.y"
	{ _gmx_selelem_free((yyvaluep->sel));        };
#line 1312 "parser.c"
	break;
      case 53: /* "string_list" */
#line 154 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1317 "parser.c"
	break;
      case 54: /* "int_list" */
#line 154 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1322 "parser.c"
	break;
      case 55: /* "int_list_item" */
#line 154 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };
#line 1327 "parser.c"
	break;
      case 56: /* "method_params" */
#line 155 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1332 "parser.c"
	break;
      case 57: /* "method_param_list" */
#line 155 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1337 "parser.c"
	break;
      case 58: /* "method_param" */
#line 155 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };
#line 1342 "parser.c"
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
#line 172 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_append_selection(sc, (yyvsp[(1) - (1)].sel), NULL);
                 if (sc->nr == nexp)
                     YYACCEPT;
             ;}
    break;

  case 3:
#line 180 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_append_selection(sc, (yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].sel));
                 if (sc->nr == nexp)
                     YYACCEPT;
             ;}
    break;

  case 4:
#line 190 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].i); ;}
    break;

  case 5:
#line 191 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); ;}
    break;

  case 6:
#line 194 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); ;}
    break;

  case 7:
#line 195 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); ;}
    break;

  case 8:
#line 198 "parser.y"
    { (yyval.str) = NULL; ;}
    break;

  case 9:
#line 199 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str);   ;}
    break;

  case 10:
#line 203 "parser.y"
    { (yyval.sel) = NULL;                            ;}
    break;

  case 11:
#line 205 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 t_selelem               *s, *p;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 s = get_group_by_id(grps, (yyvsp[(1) - (1)].i));
                 if (s == NULL) YYABORT;
                 p = _gmx_sel_init_position(sc, s, sc->spost, TRUE);
                 if (p == NULL) YYABORT;
                 (yyval.sel) = _gmx_sel_init_selection(scanner, p);
             ;}
    break;

  case 12:
#line 215 "parser.y"
    { (yyval.sel) = _gmx_sel_init_selection(scanner, (yyvsp[(1) - (1)].sel)); ;}
    break;

  case 13:
#line 216 "parser.y"
    { (yyval.sel) = _gmx_sel_init_selection(scanner, (yyvsp[(2) - (2)].sel));
                                  (yyval.sel)->name = (yyvsp[(1) - (2)].str); (yyval.sel)->u.cgrp.name = (yyvsp[(1) - (2)].str);   ;}
    break;

  case 14:
#line 219 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable(scanner, (yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel)); ;}
    break;

  case 15:
#line 221 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable(scanner, (yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel)); ;}
    break;

  case 16:
#line 223 "parser.y"
    { (yyval.sel) = _gmx_sel_assign_variable(scanner, (yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].sel)); ;}
    break;

  case 17:
#line 227 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 18:
#line 229 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_modifier(sc, (yyvsp[(2) - (3)].meth), process_param_list((yyvsp[(3) - (3)].param)), (yyvsp[(1) - (3)].sel));
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 19:
#line 235 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                              ;}
    break;

  case 20:
#line 239 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                                  (yyval.sel)->u.boolt = BOOL_NOT;
                                  (yyval.sel)->child = (yyvsp[(2) - (2)].sel);                       ;}
    break;

  case 21:
#line 243 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                                  (yyval.sel)->u.boolt = BOOL_AND;
                                  (yyval.sel)->child = (yyvsp[(1) - (3)].sel); (yyval.sel)->child->next = (yyvsp[(3) - (3)].sel); ;}
    break;

  case 22:
#line 247 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_BOOLEAN);
                                  (yyval.sel)->u.boolt = BOOL_OR;
                                  (yyval.sel)->child = (yyvsp[(1) - (3)].sel); (yyval.sel)->child->next = (yyvsp[(3) - (3)].sel); ;}
    break;

  case 23:
#line 254 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                               ;}
    break;

  case 24:
#line 258 "parser.y"
    { (yyval.sel) = get_group_by_name(grps, (yyvsp[(2) - (2)].str));
                                  sfree((yyvsp[(2) - (2)].str));
                                  if ((yyval.sel) == NULL) YYABORT;               ;}
    break;

  case 25:
#line 261 "parser.y"
    { (yyval.sel) = get_group_by_id(grps, (yyvsp[(2) - (2)].i));
                                  if ((yyval.sel) == NULL) YYABORT;               ;}
    break;

  case 26:
#line 266 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  _gmx_selelem_set_vtype((yyval.sel), (yyvsp[(1) - (1)].sel)->v.type);
                                  (yyval.sel)->name   = (yyvsp[(1) - (1)].sel)->name;
                                  (yyval.sel)->child  = (yyvsp[(1) - (1)].sel);
                                  (yyvsp[(1) - (1)].sel)->refcount++;                        ;}
    break;

  case 27:
#line 274 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  _gmx_selelem_set_vtype((yyval.sel), (yyvsp[(1) - (1)].sel)->v.type);
                                  (yyval.sel)->name   = (yyvsp[(1) - (1)].sel)->name;
                                  (yyval.sel)->child  = (yyvsp[(1) - (1)].sel);
                                  (yyvsp[(1) - (1)].sel)->refcount++;                        ;}
    break;

  case 28:
#line 281 "parser.y"
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

  case 29:
#line 294 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = init_keyword_expr(sc, (yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 30:
#line 301 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = init_keyword_expr(sc, (yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].val), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 31:
#line 308 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = init_keyword_expr(sc, (yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].val), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 32:
#line 318 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_method(sc, (yyvsp[(2) - (3)].meth), process_param_list((yyvsp[(3) - (3)].param)), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 33:
#line 329 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_comparison(sc, (yyvsp[(1) - (3)].sel), (yyvsp[(3) - (3)].sel), (yyvsp[(2) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 34:
#line 339 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype((yyval.sel), INT_VALUE);
                                  _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                                  (yyval.sel)->v.u.i[0] = (yyvsp[(1) - (1)].i);                  ;}
    break;

  case 35:
#line 343 "parser.y"
    { (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype((yyval.sel), REAL_VALUE);
                                  _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                                  (yyval.sel)->v.u.r[0] = (yyvsp[(1) - (1)].r);                  ;}
    break;

  case 36:
#line 348 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = init_keyword_expr(sc, (yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 37:
#line 355 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = init_keyword_expr(sc, (yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str));
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 38:
#line 362 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_method(sc, (yyvsp[(2) - (3)].meth), process_param_list((yyvsp[(3) - (3)].param)), (yyvsp[(1) - (3)].str));
                 if ((yyval.sel) == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 39:
#line 374 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                               ;}
    break;

  case 40:
#line 379 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_method(sc, (yyvsp[(1) - (2)].meth), process_param_list((yyvsp[(2) - (2)].param)), NULL);
                 if ((yyval.sel) == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             ;}
    break;

  case 41:
#line 390 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 42:
#line 392 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(1) - (1)].sel), sc->spost, TRUE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 43:
#line 399 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].str), TRUE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 44:
#line 409 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 45:
#line 411 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(3) - (3)].sel), (yyvsp[(1) - (3)].str), FALSE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 46:
#line 421 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); ;}
    break;

  case 47:
#line 423 "parser.y"
    {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 (yyval.sel) = _gmx_sel_init_position(sc, (yyvsp[(1) - (1)].sel), sc->rpost, FALSE);
                 if ((yyval.sel) == NULL) YYABORT;
             ;}
    break;

  case 48:
#line 433 "parser.y"
    { rvec x;
                                  (yyval.sel) = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype((yyval.sel), POS_VALUE);
                                  _gmx_selvalue_reserve(&(yyval.sel)->v, 1);
                                  x[XX] = (yyvsp[(2) - (7)].r); x[YY] = (yyvsp[(4) - (7)].r); x[ZZ] = (yyvsp[(6) - (7)].r);
                                  gmx_ana_pos_init_const((yyval.sel)->v.u.p, x);  ;}
    break;

  case 49:
#line 442 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel);                               ;}
    break;

  case 50:
#line 446 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                                  (yyval.val)->u.s = (yyvsp[(1) - (1)].str);                          ;}
    break;

  case 51:
#line 448 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                                  (yyval.val)->u.s = (yyvsp[(2) - (2)].str); (yyval.val)->next = (yyvsp[(1) - (2)].val);           ;}
    break;

  case 52:
#line 453 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val);                               ;}
    break;

  case 53:
#line 455 "parser.y"
    { (yyvsp[(2) - (2)].val)->next = (yyvsp[(1) - (2)].val); (yyval.val) = (yyvsp[(2) - (2)].val);                ;}
    break;

  case 54:
#line 459 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                                  (yyval.val)->u.i.i1 = (yyval.val)->u.i.i2 = (yyvsp[(1) - (1)].i);          ;}
    break;

  case 55:
#line 461 "parser.y"
    { (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                                  (yyval.val)->u.i.i1 = (yyvsp[(1) - (3)].i); (yyval.val)->u.i.i2 = (yyvsp[(3) - (3)].i);      ;}
    break;

  case 56:
#line 466 "parser.y"
    { (yyval.param) = (yyvsp[(1) - (1)].param); ;}
    break;

  case 57:
#line 468 "parser.y"
    { (yyval.param) = (yyvsp[(1) - (2)].param); ;}
    break;

  case 58:
#line 472 "parser.y"
    { (yyval.param) = NULL;              ;}
    break;

  case 59:
#line 474 "parser.y"
    { (yyvsp[(2) - (2)].param)->next = (yyvsp[(1) - (2)].param); (yyval.param) = (yyvsp[(2) - (2)].param); ;}
    break;

  case 60:
#line 478 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (1)].str));    ;}
    break;

  case 61:
#line 480 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = (yyvsp[(2) - (2)].val);                        ;}
    break;

  case 62:
#line 482 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value(REAL_VALUE);
                                  (yyval.param)->value->u.r = (yyvsp[(2) - (2)].r);                   ;}
    break;

  case 63:
#line 485 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value(STR_VALUE);
                                  (yyval.param)->value->u.s = (yyvsp[(2) - (2)].str);                   ;}
    break;

  case 64:
#line 489 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value_expr((yyvsp[(2) - (2)].sel)); ;}
    break;

  case 65:
#line 492 "parser.y"
    { (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                                  (yyval.param)->value = _gmx_selexpr_create_value_expr((yyvsp[(2) - (2)].sel)); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2116 "parser.c"
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


#line 496 "parser.y"


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

static t_selelem *
get_group_by_name(gmx_ana_indexgrps_t *grps, char *name)
{
    t_selelem *sel;

    if (!grps)
    {
        return NULL;
    }
    sel = _gmx_selelem_create(SEL_CONST);
    _gmx_selelem_set_vtype(sel, GROUP_VALUE);
    if (!gmx_ana_indexgrps_find(&sel->u.cgrp, grps, name))
    {
        _gmx_selelem_free(sel);
        return NULL;
    }
    sel->name = sel->u.cgrp.name;
    return sel;
}

static t_selelem *
get_group_by_id(gmx_ana_indexgrps_t *grps, int id)
{
    t_selelem *sel;

    if (!grps)
    {
        return NULL;
    }
    sel = _gmx_selelem_create(SEL_CONST);
    _gmx_selelem_set_vtype(sel, GROUP_VALUE);
    if (!gmx_ana_indexgrps_extract(&sel->u.cgrp, grps, id))
    {
        _gmx_selelem_free(sel);
        return NULL;
    }
    sel->name = sel->u.cgrp.name;
    return sel;
}

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr)
{
    t_selexpr_value *val, *pval, *nval;

    /* Count values and reverse list */
    *nr  = 0;
    pval = NULL;
    val  = values;
    while (val)
    {
        ++*nr;
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

static t_selelem *
init_keyword_expr(gmx_ana_selcollection_t *sc, gmx_ana_selmethod_t *method,
                  t_selexpr_value *values, char *rpost)
{
    t_selelem     *sel;
    int            nargs;

    values = process_value_list(values, &nargs);
    sel = _gmx_sel_init_keyword(sc, method, nargs, values, rpost);
    return sel;
}

static void
yyerror(yyscan_t scanner, int nexp, gmx_ana_indexgrps_t *grps,
        char const *s)
{
    _gmx_selparser_error("%s", s);
}

