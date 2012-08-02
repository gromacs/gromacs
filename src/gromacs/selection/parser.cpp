/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.5"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 1

/* Pull parsers.  */
#define YYPULL 0

/* Using locations.  */
#define YYLSP_NEEDED 0

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

/* Line 268 of yacc.c  */
#line 37 "parser.y"

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <exception>

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/string2.h"

#include "parsetree.h"
#include "selelem.h"

#include "scanner.h"

using gmx::sfree_guard;
using gmx::SelectionTreeElement;
using gmx::SelectionTreeElementPointer;

//! Helper method to reorder a list of parameter values and to count the values.
static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
//! Helper method to reorder a list of parameters.
static t_selexpr_param *
process_param_list(t_selexpr_param *params);

/*! \brief
 * Retrieves a selection tree pointer from a semantic value.
 *
 * \param[in] src  Semantic value to get the tree from.
 * \returns   Pointer to the selection tree.
 *
 * There should be no statements that may throw exceptions in actions before
 * this function has been called for all semantic values that have a tree
 * argument.  Together with set_sel(), this function abstracts away exception
 * safety issues that arise from the use of a plain pointer for storing the
 * selection tree semantic values.
 *
 * Does not throw.
 */
static SelectionTreeElementPointer
get_sel(SelectionTreeElementPointer *src)
{
    SelectionTreeElementPointer result;
    if (src != NULL)
    {
        result.swap(*src);
        delete src;
    }
    return result;
}
/*! \brief
 * Sets a selection tree pointer to a semantic value.
 *
 * \param[out] dest  Semantic value to set the tree to.
 * \param[in]  value Pointer to the selection tree to set.
 * \throws     std::bad_alloc if out of memory.
 *
 * This should be the last statement before ::END_ACTION, except for a
 * possible ::CHECK_SEL.
 */
static void
set_sel(SelectionTreeElementPointer *&dest,
        const SelectionTreeElementPointer &value)
{
    dest = new SelectionTreeElementPointer(value);
}
/*! \brief
 * Checks that a valid tree was set.
 *
 * Should be called after set_sel() if it was used to set a value where NULL
 * pointer indicates an error.
 *
 * \todo
 * Get rid of this macro.  It should now be possible to handle all errors using
 * exceptions.
 */
#define CHECK_SEL(sel) \
    if (!*(sel)) { \
        delete sel; \
        YYERROR; \
    }

//! Error handler needed by Bison.
static void
yyerror(yyscan_t, char const *s);

#ifdef _MSC_VER
#pragma warning(disable: 4065)
#endif

/*! \name Exception handling macros for actions
 *
 * These macros should be used at the beginning and end of each semantic action
 * that may throw an exception. For robustness, it's best to wrap all actions
 * that call functions declared outside parser.y should be wrapped.
 * These macros take care to catch any exceptions, store the exception (or
 * handle it and allow the parser to continue), and terminate the parser
 * cleanly if necessary.
 * The code calling the parser should use
 * _gmx_sel_lexer_rethrow_exception_if_occurred() to rethrow any exceptions.
 * \{
 */
//! Starts an action that may throw exceptions.
#define BEGIN_ACTION \
    try {
//! Finishes an action that may throw exceptions.
#define END_ACTION \
    } \
    catch(const std::exception &ex) \
    { \
        if (_gmx_selparser_handle_exception(scanner, ex)) \
            YYERROR; \
        else \
            YYABORT; \
    }
//!\}


/* Line 268 of yacc.c  */
#line 212 "parser.cpp"

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

/* "%code requires" blocks.  */

/* Line 288 of yacc.c  */
#line 166 "parser.y"

#include "selelem.h"



/* Line 288 of yacc.c  */
#line 242 "parser.cpp"

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

/* Line 293 of yacc.c  */
#line 170 "parser.y"

    int                         i;
    real                        r;
    char                       *str;
    struct gmx_ana_selmethod_t *meth;

    gmx::SelectionTreeElementPointer *sel;

    struct t_selexpr_value     *val;
    struct t_selexpr_param     *param;



/* Line 293 of yacc.c  */
#line 308 "parser.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

#ifndef YYPUSH_DECLS
#  define YYPUSH_DECLS
struct yypstate;
typedef struct yypstate yypstate;
enum { YYPUSH_MORE = 4 };

#if defined __STDC__ || defined __cplusplus
int yypush_parse (yypstate *yyps, int yypushed_char, YYSTYPE const *yypushed_val, void *scanner);
#else
int yypush_parse ();
#endif

#if defined __STDC__ || defined __cplusplus
yypstate * yypstate_new (void);
#else
yypstate * yypstate_new ();
#endif
#if defined __STDC__ || defined __cplusplus
void yypstate_delete (yypstate *yyps);
#else
void yypstate_delete ();
#endif
#endif


/* Copy the second part of user declarations.  */


/* Line 343 of yacc.c  */
#line 344 "parser.cpp"

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
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   417

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  49
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  26
/* YYNRULES -- Number of rules.  */
#define YYNRULES  90
/* YYNRULES -- Number of states.  */
#define YYNSTATES  149

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
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    14,    16,    18,
      20,    22,    25,    29,    33,    37,    40,    41,    44,    46,
      48,    52,    56,    58,    61,    63,    66,    68,    70,    72,
      74,    77,    81,    85,    89,    93,    96,    99,   101,   103,
     106,   110,   114,   118,   120,   122,   125,   129,   133,   137,
     141,   145,   148,   152,   156,   158,   161,   169,   173,   176,
     180,   182,   184,   186,   188,   191,   192,   195,   198,   199,
     201,   205,   207,   210,   214,   216,   220,   222,   225,   229,
     231,   233,   235,   237,   239,   241,   243,   245,   247,   251,
     255
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      50,     0,    -1,    -1,    50,    51,    -1,    52,    10,    -1,
       1,    10,    -1,    -1,    53,    -1,     6,    -1,    59,    -1,
      55,    -1,    59,    55,    -1,     9,    41,    60,    -1,     9,
      41,    62,    -1,     9,    41,    64,    -1,     4,    54,    -1,
      -1,    54,     5,    -1,    64,    -1,    60,    -1,    42,    55,
      43,    -1,    55,    23,    65,    -1,     6,    -1,    35,     6,
      -1,     7,    -1,    35,     7,    -1,    56,    -1,    57,    -1,
       8,    -1,     9,    -1,    33,    60,    -1,    60,    32,    60,
      -1,    60,    31,    60,    -1,    42,    60,    43,    -1,    62,
      28,    62,    -1,    11,    59,    -1,    11,     6,    -1,    24,
      -1,    18,    -1,    61,    19,    -1,    61,    17,    70,    -1,
      61,    16,    70,    -1,    61,    21,    65,    -1,     6,    -1,
       7,    -1,    61,    16,    -1,    61,    20,    65,    -1,    62,
      34,    62,    -1,    62,    35,    62,    -1,    62,    36,    62,
      -1,    62,    37,    62,    -1,    35,    62,    -1,    62,    39,
      62,    -1,    42,    62,    43,    -1,    59,    -1,    61,    17,
      -1,    44,    58,    45,    58,    45,    58,    46,    -1,    42,
      64,    43,    -1,    22,    65,    -1,    18,    27,    60,    -1,
      14,    -1,    13,    -1,    15,    -1,    66,    -1,    66,    26,
      -1,    -1,    66,    67,    -1,    25,    68,    -1,    -1,    69,
      -1,    47,    69,    48,    -1,    72,    -1,    69,    72,    -1,
      69,    45,    72,    -1,    71,    -1,    47,    71,    48,    -1,
      73,    -1,    71,    73,    -1,    71,    45,    73,    -1,    60,
      -1,    64,    -1,    62,    -1,    63,    -1,    74,    -1,    56,
      -1,    57,    -1,    59,    -1,    74,    -1,    56,    12,    56,
      -1,    56,    12,    57,    -1,    57,    12,    58,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   285,   285,   286,   297,   298,   320,   326,   327,   339,
     352,   358,   365,   372,   379,   390,   398,   399,   409,   410,
     417,   418,   432,   433,   437,   438,   441,   442,   445,   446,
     454,   465,   476,   487,   491,   501,   509,   519,   520,   524,
     531,   538,   548,   562,   573,   587,   594,   604,   610,   616,
     622,   628,   634,   640,   647,   658,   672,   681,   685,   695,
     708,   716,   724,   737,   739,   744,   745,   750,   759,   760,
     761,   765,   766,   768,   773,   774,   778,   779,   781,   785,
     791,   797,   803,   809,   813,   820,   827,   834,   838,   845,
     852
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
      52,    52,    52,    52,    52,    53,    54,    54,    55,    55,
      55,    55,    56,    56,    57,    57,    58,    58,    59,    59,
      60,    60,    60,    60,    60,    60,    60,    61,    61,    60,
      60,    60,    60,    62,    62,    62,    62,    62,    62,    62,
      62,    62,    62,    62,    63,    63,    64,    64,    64,    64,
      60,    62,    64,    65,    65,    66,    66,    67,    68,    68,
      68,    69,    69,    69,    70,    70,    71,    71,    71,    72,
      72,    72,    72,    72,    73,    73,    73,    73,    74,    74,
      74
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     0,     1,     1,     1,
       1,     2,     3,     3,     3,     2,     0,     2,     1,     1,
       3,     3,     1,     2,     1,     2,     1,     1,     1,     1,
       2,     3,     3,     3,     3,     2,     2,     1,     1,     2,
       3,     3,     3,     1,     1,     2,     3,     3,     3,     3,
       3,     2,     3,     3,     1,     2,     7,     3,     2,     3,
       1,     1,     1,     1,     2,     0,     2,     2,     0,     1,
       3,     1,     2,     3,     1,     3,     1,     2,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     3,
       3
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,    16,    43,    44,    28,    29,     0,
      61,    60,    62,    38,    65,    37,     0,     0,     0,     0,
       3,     0,     7,    10,     9,    19,     0,     0,    18,     5,
      15,     0,    36,    29,    35,     0,    58,    63,    43,    38,
       0,    30,     0,     0,    51,     0,    19,     0,    18,    22,
      24,     0,    26,    27,     0,     4,    65,    11,     0,     0,
      45,     0,    39,    65,    65,     0,     0,     0,     0,     0,
       0,    17,     0,    12,    13,    14,    59,    68,    64,    66,
       0,     0,    45,    20,    33,    53,    57,    23,    25,     0,
      21,    32,    31,     0,    84,    85,    86,    41,    74,    76,
      87,    40,    46,    42,    34,    47,    48,    49,    50,    52,
       0,    43,    44,     0,     0,     0,     0,    54,    79,     0,
      81,    82,    80,    67,    69,    71,    83,     0,     0,     0,
       0,     0,    77,    43,    44,     0,    55,     0,    72,     0,
      75,    88,    89,    90,    78,    70,    73,     0,    56
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    20,    21,    22,    30,    23,    94,    95,    54,
      96,   118,    26,    27,   121,   122,    36,    37,    79,   123,
     124,   101,    98,   125,    99,   100
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -72
static const yytype_int16 yypact[] =
{
     -72,   157,   -72,     4,   -72,     8,   -72,   -72,   -21,   116,
     -72,   -72,   -72,    20,   -72,   -72,   359,   243,   320,     6,
     -72,    25,   -72,    14,   320,    31,   286,   258,   -72,   -72,
      76,   345,   -72,   -72,   -72,   359,   -72,    52,   -72,   -72,
     359,   -72,   243,    18,    46,   -20,   -22,    75,    56,   -72,
     -72,    84,   -72,   -72,    50,   -72,   -72,    14,   359,   359,
      -2,   176,   -72,   -72,   -72,   243,   243,   243,   243,   243,
     243,   -72,   345,    31,   258,   -72,    31,   224,   -72,   -72,
     -22,   361,   -72,   -72,   -72,   -72,   -72,   -72,   -72,     6,
     -72,    69,   -72,    23,   108,   117,   -72,   -72,   210,   -72,
     -72,   -72,   -72,   -72,   114,    68,    68,    46,    46,    46,
      56,   118,   123,   375,   306,   108,   117,   -72,    31,   386,
     114,   -72,   -72,   -72,   266,   -72,   -72,    83,   199,     6,
       6,    23,   -72,   131,   132,   180,   176,   306,   -72,     6,
     -72,   -72,   -72,   -72,   -72,   -72,   -72,    99,   -72
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -72,   -72,   -72,   -72,   -72,   -72,    -7,     3,    17,   -46,
      -1,    24,     2,   -16,   -72,    15,    10,   -72,   -72,   -72,
      41,   100,    66,   -35,   -71,   -49
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -26
static const yytype_int16 yytable[] =
{
      24,    44,    47,    56,    49,    50,     7,    33,    34,    58,
      59,    45,    49,    50,    29,    74,    28,    57,    -8,    43,
      31,    84,    52,    83,    47,    25,    81,   132,   126,    49,
      50,     7,    33,    48,    82,    55,    53,    56,    63,    28,
      41,    51,    46,   127,    43,    93,    75,    35,    25,   104,
     105,   106,   107,   108,   109,    73,    47,   132,    51,    76,
     144,   120,    58,    59,    80,   126,    90,    43,    43,    43,
      43,    43,    43,   102,   103,   126,   117,    77,    78,   119,
     115,    71,    91,    92,   143,    70,   126,   110,   126,   138,
      87,    88,    52,   147,   116,    89,    80,    44,   120,    86,
     138,    59,   146,    65,    68,    69,    53,    70,   120,    66,
      67,    68,    69,   117,    70,    43,   119,   115,    85,   120,
     129,   120,    32,   117,     7,    33,   119,   115,   139,   130,
     -22,   116,   141,    52,   117,   -24,   117,   119,   115,   119,
     115,   116,    52,   -23,   -25,   148,   142,    53,    66,    67,
      68,    69,   116,    70,   116,   135,    53,     2,     3,   128,
      97,     4,     0,     5,     6,     7,     8,    -6,     9,     0,
      10,    11,    12,     0,     0,    13,     0,     0,     0,    14,
       0,    15,    49,    50,     7,    33,   111,   112,     7,    33,
      16,     9,    17,    10,    11,    12,     0,     0,    13,    18,
       0,    19,    14,     0,    15,    49,    50,     7,    33,     0,
       0,    51,     0,    16,     0,   113,    49,    50,     7,    33,
       0,     0,    72,    93,    19,   137,     0,     0,   145,     0,
     111,   112,     7,    33,    51,     9,     0,    10,    11,    12,
       0,     0,    13,     0,   131,    51,    14,   140,    15,    38,
       6,     0,     0,     0,     0,   131,    10,    16,     0,   113,
       0,    39,     0,     0,     0,     0,    72,    15,    19,     0,
       0,   114,   111,   112,     7,    33,     0,     9,    17,    10,
      11,    12,     0,     0,    13,    42,    65,     0,    14,     0,
      15,     0,    66,    67,    68,    69,     0,    70,     0,    16,
       0,   113,    60,    61,     0,    62,    63,    64,    72,     0,
      19,   137,   111,   112,     7,    33,     0,     9,     0,    10,
      11,    12,     0,     0,    13,     0,    38,     6,    14,     0,
      15,     9,     0,    10,    11,    12,     0,     0,    13,    16,
       0,   113,    14,     0,    15,     0,     0,     0,    72,     0,
      19,    38,     6,    16,     0,    17,     9,     0,    10,    11,
      12,     0,    18,    13,    19,    38,     6,    14,     0,    15,
       9,     0,    10,    11,     0,     0,     0,    39,    16,     0,
      17,   133,   134,    15,     0,     0,     0,    72,    10,    19,
       0,     0,    16,    39,    17,    66,    67,    68,    69,    15,
      70,    40,    60,   136,    85,    62,    63,    64,     0,     0,
      17,     0,     0,     0,     0,     0,     0,    42
};

#define yypact_value_is_default(yystate) \
  ((yystate) == (-72))

#define yytable_value_is_error(yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
       1,    17,    18,    23,     6,     7,     8,     9,     9,    31,
      32,    18,     6,     7,    10,    31,     1,    24,    10,    17,
      41,    43,    19,    43,    40,     1,    42,    98,    77,     6,
       7,     8,     9,    18,    16,    10,    19,    23,    20,    24,
      16,    35,    18,    89,    42,    47,    31,    27,    24,    65,
      66,    67,    68,    69,    70,    31,    72,   128,    35,    35,
     131,    77,    31,    32,    40,   114,    56,    65,    66,    67,
      68,    69,    70,    63,    64,   124,    77,    25,    26,    77,
      77,     5,    58,    59,   130,    39,   135,    72,   137,   124,
       6,     7,    89,   139,    77,    45,    72,   113,   114,    43,
     135,    32,   137,    28,    36,    37,    89,    39,   124,    34,
      35,    36,    37,   114,    39,   113,   114,   114,    43,   135,
      12,   137,     6,   124,     8,     9,   124,   124,    45,    12,
      12,   114,   129,   130,   135,    12,   137,   135,   135,   137,
     137,   124,   139,    12,    12,    46,   129,   130,    34,    35,
      36,    37,   135,    39,   137,   114,   139,     0,     1,    93,
      60,     4,    -1,     6,     7,     8,     9,    10,    11,    -1,
      13,    14,    15,    -1,    -1,    18,    -1,    -1,    -1,    22,
      -1,    24,     6,     7,     8,     9,     6,     7,     8,     9,
      33,    11,    35,    13,    14,    15,    -1,    -1,    18,    42,
      -1,    44,    22,    -1,    24,     6,     7,     8,     9,    -1,
      -1,    35,    -1,    33,    -1,    35,     6,     7,     8,     9,
      -1,    -1,    42,    47,    44,    45,    -1,    -1,    48,    -1,
       6,     7,     8,     9,    35,    11,    -1,    13,    14,    15,
      -1,    -1,    18,    -1,    45,    35,    22,    48,    24,     6,
       7,    -1,    -1,    -1,    -1,    45,    13,    33,    -1,    35,
      -1,    18,    -1,    -1,    -1,    -1,    42,    24,    44,    -1,
      -1,    47,     6,     7,     8,     9,    -1,    11,    35,    13,
      14,    15,    -1,    -1,    18,    42,    28,    -1,    22,    -1,
      24,    -1,    34,    35,    36,    37,    -1,    39,    -1,    33,
      -1,    35,    16,    17,    -1,    19,    20,    21,    42,    -1,
      44,    45,     6,     7,     8,     9,    -1,    11,    -1,    13,
      14,    15,    -1,    -1,    18,    -1,     6,     7,    22,    -1,
      24,    11,    -1,    13,    14,    15,    -1,    -1,    18,    33,
      -1,    35,    22,    -1,    24,    -1,    -1,    -1,    42,    -1,
      44,     6,     7,    33,    -1,    35,    11,    -1,    13,    14,
      15,    -1,    42,    18,    44,     6,     7,    22,    -1,    24,
      11,    -1,    13,    14,    -1,    -1,    -1,    18,    33,    -1,
      35,     6,     7,    24,    -1,    -1,    -1,    42,    13,    44,
      -1,    -1,    33,    18,    35,    34,    35,    36,    37,    24,
      39,    42,    16,    17,    43,    19,    20,    21,    -1,    -1,
      35,    -1,    -1,    -1,    -1,    -1,    -1,    42
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    50,     0,     1,     4,     6,     7,     8,     9,    11,
      13,    14,    15,    18,    22,    24,    33,    35,    42,    44,
      51,    52,    53,    55,    59,    60,    61,    62,    64,    10,
      54,    41,     6,     9,    59,    27,    65,    66,     6,    18,
      42,    60,    42,    61,    62,    55,    60,    62,    64,     6,
       7,    35,    56,    57,    58,    10,    23,    55,    31,    32,
      16,    17,    19,    20,    21,    28,    34,    35,    36,    37,
      39,     5,    42,    60,    62,    64,    60,    25,    26,    67,
      60,    62,    16,    43,    43,    43,    43,     6,     7,    45,
      65,    60,    60,    47,    56,    57,    59,    70,    71,    73,
      74,    70,    65,    65,    62,    62,    62,    62,    62,    62,
      64,     6,     7,    35,    47,    56,    57,    59,    60,    61,
      62,    63,    64,    68,    69,    72,    74,    58,    71,    12,
      12,    45,    73,     6,     7,    69,    17,    45,    72,    45,
      48,    56,    57,    58,    73,    48,    72,    58,    46
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

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
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
  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  YYSIZE_T yysize1;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = 0;
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
                yysize1 = yysize + yytnamerr (0, yytname[yyx]);
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
      case 5: /* "HELP_TOPIC" */

/* Line 1391 of yacc.c  */
#line 265 "parser.y"
	{ free((yyvaluep->str));                     };

/* Line 1391 of yacc.c  */
#line 1450 "parser.cpp"
	break;
      case 8: /* "STR" */

/* Line 1391 of yacc.c  */
#line 265 "parser.y"
	{ free((yyvaluep->str));                     };

/* Line 1391 of yacc.c  */
#line 1459 "parser.cpp"
	break;
      case 9: /* "IDENTIFIER" */

/* Line 1391 of yacc.c  */
#line 265 "parser.y"
	{ free((yyvaluep->str));                     };

/* Line 1391 of yacc.c  */
#line 1468 "parser.cpp"
	break;
      case 25: /* "PARAM" */

/* Line 1391 of yacc.c  */
#line 266 "parser.y"
	{ if((yyvaluep->str)) free((yyvaluep->str));              };

/* Line 1391 of yacc.c  */
#line 1477 "parser.cpp"
	break;
      case 28: /* "CMP_OP" */

/* Line 1391 of yacc.c  */
#line 265 "parser.y"
	{ free((yyvaluep->str));                     };

/* Line 1391 of yacc.c  */
#line 1486 "parser.cpp"
	break;
      case 50: /* "commands" */

/* Line 1391 of yacc.c  */
#line 267 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1495 "parser.cpp"
	break;
      case 51: /* "command" */

/* Line 1391 of yacc.c  */
#line 267 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1504 "parser.cpp"
	break;
      case 52: /* "cmd_plain" */

/* Line 1391 of yacc.c  */
#line 267 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1513 "parser.cpp"
	break;
      case 54: /* "help_topic" */

/* Line 1391 of yacc.c  */
#line 272 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1522 "parser.cpp"
	break;
      case 55: /* "selection" */

/* Line 1391 of yacc.c  */
#line 267 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1531 "parser.cpp"
	break;
      case 59: /* "string" */

/* Line 1391 of yacc.c  */
#line 265 "parser.y"
	{ free((yyvaluep->str));                     };

/* Line 1391 of yacc.c  */
#line 1540 "parser.cpp"
	break;
      case 60: /* "sel_expr" */

/* Line 1391 of yacc.c  */
#line 268 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1549 "parser.cpp"
	break;
      case 62: /* "num_expr" */

/* Line 1391 of yacc.c  */
#line 268 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1558 "parser.cpp"
	break;
      case 63: /* "str_expr" */

/* Line 1391 of yacc.c  */
#line 268 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1567 "parser.cpp"
	break;
      case 64: /* "pos_expr" */

/* Line 1391 of yacc.c  */
#line 268 "parser.y"
	{ delete (yyvaluep->sel);                    };

/* Line 1391 of yacc.c  */
#line 1576 "parser.cpp"
	break;
      case 65: /* "method_params" */

/* Line 1391 of yacc.c  */
#line 269 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };

/* Line 1391 of yacc.c  */
#line 1585 "parser.cpp"
	break;
      case 66: /* "method_param_list" */

/* Line 1391 of yacc.c  */
#line 269 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };

/* Line 1391 of yacc.c  */
#line 1594 "parser.cpp"
	break;
      case 67: /* "method_param" */

/* Line 1391 of yacc.c  */
#line 269 "parser.y"
	{ _gmx_selexpr_free_params((yyvaluep->param)); };

/* Line 1391 of yacc.c  */
#line 1603 "parser.cpp"
	break;
      case 68: /* "value_list" */

/* Line 1391 of yacc.c  */
#line 270 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1612 "parser.cpp"
	break;
      case 69: /* "value_list_contents" */

/* Line 1391 of yacc.c  */
#line 270 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1621 "parser.cpp"
	break;
      case 70: /* "basic_value_list" */

/* Line 1391 of yacc.c  */
#line 271 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1630 "parser.cpp"
	break;
      case 71: /* "basic_value_list_contents" */

/* Line 1391 of yacc.c  */
#line 271 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1639 "parser.cpp"
	break;
      case 72: /* "value_item" */

/* Line 1391 of yacc.c  */
#line 270 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1648 "parser.cpp"
	break;
      case 73: /* "basic_value_item" */

/* Line 1391 of yacc.c  */
#line 271 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1657 "parser.cpp"
	break;
      case 74: /* "value_item_range" */

/* Line 1391 of yacc.c  */
#line 270 "parser.y"
	{ _gmx_selexpr_free_values((yyvaluep->val)); };

/* Line 1391 of yacc.c  */
#line 1666 "parser.cpp"
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

       Refer to the stacks thru separate pointers, to allow yyoverflow
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
    return 0;
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

/* Line 1806 of yacc.c  */
#line 285 "parser.y"
    { (yyval.sel) = NULL; }
    break;

  case 3:

/* Line 1806 of yacc.c  */
#line 287 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_append_selection(get_sel((yyvsp[(2) - (2)].sel)), get_sel((yyvsp[(1) - (2)].sel)), scanner));
                 if (_gmx_sel_parser_should_finish(scanner))
                     YYACCEPT;
                 END_ACTION;
             }
    break;

  case 4:

/* Line 1806 of yacc.c  */
#line 297 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (2)].sel); }
    break;

  case 5:

/* Line 1806 of yacc.c  */
#line 299 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.sel) = NULL;
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
                 END_ACTION;
             }
    break;

  case 6:

/* Line 1806 of yacc.c  */
#line 320 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.sel) = NULL;
                 _gmx_sel_handle_empty_cmd(scanner);
                 END_ACTION;
             }
    break;

  case 7:

/* Line 1806 of yacc.c  */
#line 326 "parser.y"
    { (yyval.sel) = NULL; }
    break;

  case 8:

/* Line 1806 of yacc.c  */
#line 328 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer s
                        = _gmx_sel_init_group_by_id((yyvsp[(1) - (1)].i), scanner);
                 if (!s) YYERROR;
                 SelectionTreeElementPointer p
                        = _gmx_sel_init_position(s, NULL, scanner);
                 if (!p) YYERROR;
                 set_sel((yyval.sel), _gmx_sel_init_selection(s->name, p, scanner));
                 END_ACTION;
             }
    break;

  case 9:

/* Line 1806 of yacc.c  */
#line 340 "parser.y"
    {
                 BEGIN_ACTION;
                 sfree_guard nameGuard((yyvsp[(1) - (1)].str));
                 SelectionTreeElementPointer s
                        = _gmx_sel_init_group_by_name((yyvsp[(1) - (1)].str), scanner);
                 if (!s) YYERROR;
                 SelectionTreeElementPointer p
                        = _gmx_sel_init_position(s, NULL, scanner);
                 if (!p) YYERROR;
                 set_sel((yyval.sel), _gmx_sel_init_selection(s->name, p, scanner));
                 END_ACTION;
             }
    break;

  case 10:

/* Line 1806 of yacc.c  */
#line 353 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_selection(NULL, get_sel((yyvsp[(1) - (1)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 11:

/* Line 1806 of yacc.c  */
#line 359 "parser.y"
    {
                 BEGIN_ACTION;
                 sfree_guard nameGuard((yyvsp[(1) - (2)].str));
                 set_sel((yyval.sel), _gmx_sel_init_selection((yyvsp[(1) - (2)].str), get_sel((yyvsp[(2) - (2)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 12:

/* Line 1806 of yacc.c  */
#line 366 "parser.y"
    {
                 BEGIN_ACTION;
                 sfree_guard nameGuard((yyvsp[(1) - (3)].str));
                 set_sel((yyval.sel), _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), get_sel((yyvsp[(3) - (3)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 13:

/* Line 1806 of yacc.c  */
#line 373 "parser.y"
    {
                 BEGIN_ACTION;
                 sfree_guard nameGuard((yyvsp[(1) - (3)].str));
                 set_sel((yyval.sel), _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), get_sel((yyvsp[(3) - (3)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 14:

/* Line 1806 of yacc.c  */
#line 380 "parser.y"
    {
                 BEGIN_ACTION;
                 sfree_guard nameGuard((yyvsp[(1) - (3)].str));
                 set_sel((yyval.sel), _gmx_sel_assign_variable((yyvsp[(1) - (3)].str), get_sel((yyvsp[(3) - (3)].sel)), scanner));
                 END_ACTION;
             }
    break;

  case 15:

/* Line 1806 of yacc.c  */
#line 391 "parser.y"
    {
                 BEGIN_ACTION;
                 _gmx_sel_handle_help_cmd(process_value_list((yyvsp[(2) - (2)].val), NULL), scanner);
                 END_ACTION;
             }
    break;

  case 16:

/* Line 1806 of yacc.c  */
#line 398 "parser.y"
    { (yyval.val) = NULL; }
    break;

  case 17:

/* Line 1806 of yacc.c  */
#line 400 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                 (yyval.val)->u.s = (yyvsp[(2) - (2)].str); (yyval.val)->next = (yyvsp[(1) - (2)].val);
                 END_ACTION;
             }
    break;

  case 18:

/* Line 1806 of yacc.c  */
#line 409 "parser.y"
    { (yyval.sel) = (yyvsp[(1) - (1)].sel); }
    break;

  case 19:

/* Line 1806 of yacc.c  */
#line 411 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_position(get_sel((yyvsp[(1) - (1)].sel)), NULL, scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 20:

/* Line 1806 of yacc.c  */
#line 417 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 21:

/* Line 1806 of yacc.c  */
#line 419 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_modifier((yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].param), get_sel((yyvsp[(1) - (3)].sel)), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 22:

/* Line 1806 of yacc.c  */
#line 432 "parser.y"
    { (yyval.i) = (yyvsp[(1) - (1)].i); }
    break;

  case 23:

/* Line 1806 of yacc.c  */
#line 433 "parser.y"
    { (yyval.i) = -(yyvsp[(2) - (2)].i); }
    break;

  case 24:

/* Line 1806 of yacc.c  */
#line 437 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); }
    break;

  case 25:

/* Line 1806 of yacc.c  */
#line 438 "parser.y"
    { (yyval.r) = -(yyvsp[(2) - (2)].r); }
    break;

  case 26:

/* Line 1806 of yacc.c  */
#line 441 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].i); }
    break;

  case 27:

/* Line 1806 of yacc.c  */
#line 442 "parser.y"
    { (yyval.r) = (yyvsp[(1) - (1)].r); }
    break;

  case 28:

/* Line 1806 of yacc.c  */
#line 445 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); }
    break;

  case 29:

/* Line 1806 of yacc.c  */
#line 446 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str); }
    break;

  case 30:

/* Line 1806 of yacc.c  */
#line 455 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg(get_sel((yyvsp[(2) - (2)].sel)));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN));
                 sel->u.boolt = BOOL_NOT;
                 sel->child = arg;
                 set_sel((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 31:

/* Line 1806 of yacc.c  */
#line 466 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg1(get_sel((yyvsp[(1) - (3)].sel))), arg2(get_sel((yyvsp[(3) - (3)].sel)));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN));
                 sel->u.boolt = BOOL_AND;
                 sel->child = arg1; sel->child->next = arg2;
                 set_sel((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 32:

/* Line 1806 of yacc.c  */
#line 477 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg1(get_sel((yyvsp[(1) - (3)].sel))), arg2(get_sel((yyvsp[(3) - (3)].sel)));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN));
                 sel->u.boolt = BOOL_OR;
                 sel->child = arg1; sel->child->next = arg2;
                 set_sel((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 33:

/* Line 1806 of yacc.c  */
#line 487 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 34:

/* Line 1806 of yacc.c  */
#line 492 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_comparison(get_sel((yyvsp[(1) - (3)].sel)), get_sel((yyvsp[(3) - (3)].sel)), (yyvsp[(2) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 35:

/* Line 1806 of yacc.c  */
#line 502 "parser.y"
    {
                 BEGIN_ACTION;
                 sfree_guard nameGuard((yyvsp[(2) - (2)].str));
                 set_sel((yyval.sel), _gmx_sel_init_group_by_name((yyvsp[(2) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 36:

/* Line 1806 of yacc.c  */
#line 510 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_group_by_id((yyvsp[(2) - (2)].i), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 37:

/* Line 1806 of yacc.c  */
#line 519 "parser.y"
    { (yyval.str) = NULL; }
    break;

  case 38:

/* Line 1806 of yacc.c  */
#line 520 "parser.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str);   }
    break;

  case 39:

/* Line 1806 of yacc.c  */
#line 525 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 40:

/* Line 1806 of yacc.c  */
#line 532 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (3)].meth), process_value_list((yyvsp[(3) - (3)].val), NULL), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 41:

/* Line 1806 of yacc.c  */
#line 539 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (3)].meth), process_value_list((yyvsp[(3) - (3)].val), NULL), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 42:

/* Line 1806 of yacc.c  */
#line 549 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_method((yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].param), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 43:

/* Line 1806 of yacc.c  */
#line 563 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST));
                 _gmx_selelem_set_vtype(sel, INT_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.i[0] = (yyvsp[(1) - (1)].i);
                 set_sel((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 44:

/* Line 1806 of yacc.c  */
#line 574 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST));
                 _gmx_selelem_set_vtype(sel, REAL_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.r[0] = (yyvsp[(1) - (1)].r);
                 set_sel((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 45:

/* Line 1806 of yacc.c  */
#line 588 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 46:

/* Line 1806 of yacc.c  */
#line 595 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_method((yyvsp[(2) - (3)].meth), (yyvsp[(3) - (3)].param), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 47:

/* Line 1806 of yacc.c  */
#line 605 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_arithmetic(get_sel((yyvsp[(1) - (3)].sel)), get_sel((yyvsp[(3) - (3)].sel)), '+', scanner));
                 END_ACTION;
             }
    break;

  case 48:

/* Line 1806 of yacc.c  */
#line 611 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_arithmetic(get_sel((yyvsp[(1) - (3)].sel)), get_sel((yyvsp[(3) - (3)].sel)), '-', scanner));
                 END_ACTION;
             }
    break;

  case 49:

/* Line 1806 of yacc.c  */
#line 617 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_arithmetic(get_sel((yyvsp[(1) - (3)].sel)), get_sel((yyvsp[(3) - (3)].sel)), '*', scanner));
                 END_ACTION;
             }
    break;

  case 50:

/* Line 1806 of yacc.c  */
#line 623 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_arithmetic(get_sel((yyvsp[(1) - (3)].sel)), get_sel((yyvsp[(3) - (3)].sel)), '/', scanner));
                 END_ACTION;
             }
    break;

  case 51:

/* Line 1806 of yacc.c  */
#line 629 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_arithmetic(get_sel((yyvsp[(2) - (2)].sel)), SelectionTreeElementPointer(), '-', scanner));
                 END_ACTION;
             }
    break;

  case 52:

/* Line 1806 of yacc.c  */
#line 635 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_arithmetic(get_sel((yyvsp[(1) - (3)].sel)), get_sel((yyvsp[(3) - (3)].sel)), '^', scanner));
                 END_ACTION;
             }
    break;

  case 53:

/* Line 1806 of yacc.c  */
#line 640 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 54:

/* Line 1806 of yacc.c  */
#line 648 "parser.y"
    {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST));
                 _gmx_selelem_set_vtype(sel, STR_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.s[0] = (yyvsp[(1) - (1)].str);
                 set_sel((yyval.sel), sel);
                 END_ACTION;
             }
    break;

  case 55:

/* Line 1806 of yacc.c  */
#line 659 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_keyword((yyvsp[(2) - (2)].meth), NULL, (yyvsp[(1) - (2)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 56:

/* Line 1806 of yacc.c  */
#line 673 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_const_position((yyvsp[(2) - (7)].r), (yyvsp[(4) - (7)].r), (yyvsp[(6) - (7)].r)));
                 END_ACTION;
             }
    break;

  case 57:

/* Line 1806 of yacc.c  */
#line 681 "parser.y"
    { (yyval.sel) = (yyvsp[(2) - (3)].sel); }
    break;

  case 58:

/* Line 1806 of yacc.c  */
#line 686 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_method((yyvsp[(1) - (2)].meth), (yyvsp[(2) - (2)].param), NULL, scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 59:

/* Line 1806 of yacc.c  */
#line 696 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_position(get_sel((yyvsp[(3) - (3)].sel)), (yyvsp[(1) - (3)].str), scanner));
                 CHECK_SEL((yyval.sel));
                 END_ACTION;
             }
    break;

  case 60:

/* Line 1806 of yacc.c  */
#line 709 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_variable_ref(get_sel((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 61:

/* Line 1806 of yacc.c  */
#line 717 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_variable_ref(get_sel((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 62:

/* Line 1806 of yacc.c  */
#line 725 "parser.y"
    {
                 BEGIN_ACTION;
                 set_sel((yyval.sel), _gmx_sel_init_variable_ref(get_sel((yyvsp[(1) - (1)].sel))));
                 END_ACTION;
             }
    break;

  case 63:

/* Line 1806 of yacc.c  */
#line 738 "parser.y"
    { (yyval.param) = process_param_list((yyvsp[(1) - (1)].param)); }
    break;

  case 64:

/* Line 1806 of yacc.c  */
#line 740 "parser.y"
    { (yyval.param) = process_param_list((yyvsp[(1) - (2)].param)); }
    break;

  case 65:

/* Line 1806 of yacc.c  */
#line 744 "parser.y"
    { (yyval.param) = NULL;              }
    break;

  case 66:

/* Line 1806 of yacc.c  */
#line 746 "parser.y"
    { (yyvsp[(2) - (2)].param)->next = (yyvsp[(1) - (2)].param); (yyval.param) = (yyvsp[(2) - (2)].param); }
    break;

  case 67:

/* Line 1806 of yacc.c  */
#line 751 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.param) = _gmx_selexpr_create_param((yyvsp[(1) - (2)].str));
                 (yyval.param)->value = process_value_list((yyvsp[(2) - (2)].val), &(yyval.param)->nval);
                 END_ACTION;
             }
    break;

  case 68:

/* Line 1806 of yacc.c  */
#line 759 "parser.y"
    { (yyval.val) = NULL; }
    break;

  case 69:

/* Line 1806 of yacc.c  */
#line 760 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val);   }
    break;

  case 70:

/* Line 1806 of yacc.c  */
#line 761 "parser.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val);   }
    break;

  case 71:

/* Line 1806 of yacc.c  */
#line 765 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 72:

/* Line 1806 of yacc.c  */
#line 767 "parser.y"
    { (yyvsp[(2) - (2)].val)->next = (yyvsp[(1) - (2)].val); (yyval.val) = (yyvsp[(2) - (2)].val); }
    break;

  case 73:

/* Line 1806 of yacc.c  */
#line 769 "parser.y"
    { (yyvsp[(3) - (3)].val)->next = (yyvsp[(1) - (3)].val); (yyval.val) = (yyvsp[(3) - (3)].val); }
    break;

  case 74:

/* Line 1806 of yacc.c  */
#line 773 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 75:

/* Line 1806 of yacc.c  */
#line 774 "parser.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val); }
    break;

  case 76:

/* Line 1806 of yacc.c  */
#line 778 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 77:

/* Line 1806 of yacc.c  */
#line 780 "parser.y"
    { (yyvsp[(2) - (2)].val)->next = (yyvsp[(1) - (2)].val); (yyval.val) = (yyvsp[(2) - (2)].val); }
    break;

  case 78:

/* Line 1806 of yacc.c  */
#line 782 "parser.y"
    { (yyvsp[(3) - (3)].val)->next = (yyvsp[(1) - (3)].val); (yyval.val) = (yyvsp[(3) - (3)].val); }
    break;

  case 79:

/* Line 1806 of yacc.c  */
#line 786 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value_expr(get_sel((yyvsp[(1) - (1)].sel)));
                 END_ACTION;
             }
    break;

  case 80:

/* Line 1806 of yacc.c  */
#line 792 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value_expr(get_sel((yyvsp[(1) - (1)].sel)));
                 END_ACTION;
             }
    break;

  case 81:

/* Line 1806 of yacc.c  */
#line 798 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value_expr(get_sel((yyvsp[(1) - (1)].sel)));
                 END_ACTION;
             }
    break;

  case 82:

/* Line 1806 of yacc.c  */
#line 804 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value_expr(get_sel((yyvsp[(1) - (1)].sel)));
                 END_ACTION;
             }
    break;

  case 83:

/* Line 1806 of yacc.c  */
#line 809 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 84:

/* Line 1806 of yacc.c  */
#line 814 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                 (yyval.val)->u.i.i1 = (yyval.val)->u.i.i2 = (yyvsp[(1) - (1)].i);
                 END_ACTION;
             }
    break;

  case 85:

/* Line 1806 of yacc.c  */
#line 821 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(REAL_VALUE);
                 (yyval.val)->u.r.r1 = (yyval.val)->u.r.r2 = (yyvsp[(1) - (1)].r);
                 END_ACTION;
             }
    break;

  case 86:

/* Line 1806 of yacc.c  */
#line 828 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(STR_VALUE);
                 (yyval.val)->u.s = (yyvsp[(1) - (1)].str);
                 END_ACTION;
             }
    break;

  case 87:

/* Line 1806 of yacc.c  */
#line 834 "parser.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); }
    break;

  case 88:

/* Line 1806 of yacc.c  */
#line 839 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(INT_VALUE);
                 (yyval.val)->u.i.i1 = (yyvsp[(1) - (3)].i); (yyval.val)->u.i.i2 = (yyvsp[(3) - (3)].i);
                 END_ACTION;
             }
    break;

  case 89:

/* Line 1806 of yacc.c  */
#line 846 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(REAL_VALUE);
                 (yyval.val)->u.r.r1 = (yyvsp[(1) - (3)].i); (yyval.val)->u.r.r2 = (yyvsp[(3) - (3)].r);
                 END_ACTION;
             }
    break;

  case 90:

/* Line 1806 of yacc.c  */
#line 853 "parser.y"
    {
                 BEGIN_ACTION;
                 (yyval.val) = _gmx_selexpr_create_value(REAL_VALUE);
                 (yyval.val)->u.r.r1 = (yyvsp[(1) - (3)].r); (yyval.val)->u.r.r2 = (yyvsp[(3) - (3)].r);
                 END_ACTION;
             }
    break;



/* Line 1806 of yacc.c  */
#line 2935 "parser.cpp"
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

#if !defined(yyoverflow) || YYERROR_VERBOSE
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



/* Line 2067 of yacc.c  */
#line 861 "parser.y"


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
    _gmx_selparser_error(scanner, "%s", s);
}

