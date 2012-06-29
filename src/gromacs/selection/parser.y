/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief Grammar description and parser for the selection language.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
%{
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

#include "gromacs/legacyheaders/string2.h"

#include "parsetree.h"
#include "selelem.h"

#include "scanner.h"

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
static t_selexpr_param *
process_param_list(t_selexpr_param *params);

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
%}

%union{
    int                         i;
    real                        r;
    char                       *str;
    struct gmx_ana_selmethod_t *meth;

    struct t_selelem           *sel;

    struct t_selexpr_value     *val;
    struct t_selexpr_param     *param;
};

/* Invalid token to report lexer errors */
%token INVALID

/* Tokens for help requests */
%token         HELP
%token <str>   HELP_TOPIC

/* Simple input tokens */
%token <i>     TOK_INT
%token <r>     TOK_REAL
%token <str>   STR
%token <str>   IDENTIFIER
%token         CMD_SEP

/* Simple keyword tokens */
%token         GROUP
%token         TO

/* Variable tokens */
%token <sel>   VARIABLE_NUMERIC
%token <sel>   VARIABLE_GROUP
%token <sel>   VARIABLE_POS

/* Selection method tokens */
%token <meth>  KEYWORD_NUMERIC
%token <meth>  KEYWORD_STR
%token <str>   KEYWORD_POS
%token <meth>  KEYWORD_GROUP
%token <meth>  METHOD_NUMERIC
%token <meth>  METHOD_GROUP
%token <meth>  METHOD_POS
%token <meth>  MODIFIER
/* Empty token that should precede any non-position KEYWORD/METHOD token that
 * is not preceded by KEYWORD_POS. This is used to work around reduce/reduce
 * conflicts that appear when a lookahead token would require a reduction of
 * a rule with empty RHS before shifting, and there is an alternative reduction
 * available. Replacing the empty RHS with a dummy token makes these conflicts
 * only shift/reduce conflicts. Another alternative would be to remove the
 * pos_mod non-terminal completely and split each rule that uses it into two,
 * but this would require duplicating six rules in the grammar. */
%token         EMPTY_POSMOD

%token <str>   PARAM
%token         END_OF_METHOD

%token          OF
/* Comparison operators have lower precedence than parameter reduction
 * to make it possible to parse, e.g., "mindist from resnr 1 < 2" without
 * parenthesis. */
%nonassoc <str> CMP_OP
/* A dummy token that determines the precedence of parameter reduction */
%nonassoc       PARAM_REDUCT
/* Boolean operator tokens */
%left           OR XOR
%left           AND
%left           NOT
/* Arithmetic operator tokens */
%left           '+' '-'
%left           '*' '/'
%right          UNARY_NEG   /* Dummy token for unary negation precedence */
%right          '^'
%nonassoc       NUM_REDUCT  /* Dummy token for numerical keyword reduction precedence */

/* Simple non-terminals */
%type <i>     integer_number
%type <r>     real_number number
%type <str>   string
%type <str>   pos_mod

/* Expression non-terminals */
%type <sel>   commands command cmd_plain
%type <sel>   selection
%type <sel>   sel_expr
%type <sel>   num_expr
%type <sel>   str_expr
%type <sel>   pos_expr

/* Parameter/value non-terminals */
%type <param> method_params method_param_list method_param
%type <val>   value_list value_list_contents value_item value_item_range
%type <val>   basic_value_list basic_value_list_contents basic_value_item
%type <val>   help_topic

%destructor { free($$);                     } HELP_TOPIC STR IDENTIFIER CMP_OP string
%destructor { if($$) free($$);              } PARAM
%destructor { if($$) _gmx_selelem_free($$); } command cmd_plain
%destructor { _gmx_selelem_free_chain($$);  } selection
%destructor { _gmx_selelem_free($$);        } sel_expr num_expr str_expr pos_expr
%destructor { _gmx_selexpr_free_params($$); } method_params method_param_list method_param
%destructor { _gmx_selexpr_free_values($$); } value_list value_list_contents value_item value_item_range
%destructor { _gmx_selexpr_free_values($$); } basic_value_list basic_value_list_contents basic_value_item
%destructor { _gmx_selexpr_free_values($$); } help_topic

%expect 50
%debug
%pure-parser
%define api.push-pull push

%name-prefix="_gmx_sel_yy"
%parse-param { void *scanner }

%%

/* The start rule: allow one or more commands */
commands:    /* empty */        { $$ = NULL; }
           | commands command
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_append_selection($2, $1, scanner);
                 if (_gmx_sel_parser_should_finish(scanner))
                     YYACCEPT;
                 END_ACTION;
             }
;

/* A command is formed from an actual command and a separator */
command:     cmd_plain CMD_SEP  { $$ = $1; }
           | error CMD_SEP
             {
                 BEGIN_ACTION;
                 $$ = NULL;
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
;

/* Commands can be selections or variable assignments */
cmd_plain:   /* empty */
             {
                 BEGIN_ACTION;
                 $$ = NULL;
                 _gmx_sel_handle_empty_cmd(scanner);
                 END_ACTION;
             }
           | help_request       { $$ = NULL; }
           | TOK_INT
             {
                 BEGIN_ACTION;
                 t_selelem *s, *p;
                 s = _gmx_sel_init_group_by_id($1, scanner);
                 if (s == NULL) YYERROR;
                 p = _gmx_sel_init_position(s, NULL, scanner);
                 if (p == NULL) YYERROR;
                 $$ = _gmx_sel_init_selection(strdup(s->name), p, scanner);
                 END_ACTION;
             }
           | string
             {
                 BEGIN_ACTION;
                 t_selelem *s, *p;
                 s = _gmx_sel_init_group_by_name($1, scanner);
                 free($1);
                 if (s == NULL) YYERROR;
                 p = _gmx_sel_init_position(s, NULL, scanner);
                 if (p == NULL) YYERROR;
                 $$ = _gmx_sel_init_selection(strdup(s->name), p, scanner);
                 END_ACTION;
             }
           | selection
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_selection(NULL, $1, scanner);
                 END_ACTION;
             }
           | string selection
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_selection($1, $2, scanner);
                 END_ACTION;
             }
           | IDENTIFIER '=' sel_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_assign_variable($1, $3, scanner);
                 END_ACTION;
             }
           | IDENTIFIER '=' num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_assign_variable($1, $3, scanner);
                 END_ACTION;
             }
           | IDENTIFIER '=' pos_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_assign_variable($1, $3, scanner);
                 END_ACTION;
             }
;

/* Help requests */
help_request:
             HELP help_topic
             {
                 BEGIN_ACTION;
                 _gmx_sel_handle_help_cmd(process_value_list($2, NULL), scanner);
                 END_ACTION;
             }
;

help_topic:  /* empty */            { $$ = NULL; }
           | help_topic HELP_TOPIC
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(STR_VALUE);
                 $$->u.s = $2; $$->next = $1;
                 END_ACTION;
             }
;

/* Selection is made of an expression and zero or more modifiers */
selection:   pos_expr           { $$ = $1; }
           | sel_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_position($1, NULL, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
           | '(' selection ')'  { $$ = $2; }
           | selection MODIFIER method_params
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_modifier($2, $3, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/********************************************************************
 * BASIC NON-TERMINAL SYMBOLS
 ********************************************************************/

integer_number:
             TOK_INT            { $$ = $1; }
           | '-' TOK_INT        { $$ = -$2; }
;

real_number:
             TOK_REAL           { $$ = $1; }
           | '-' TOK_REAL       { $$ = -$2; }
;

number:      integer_number     { $$ = $1; }
           | real_number        { $$ = $1; }
;

string:      STR                { $$ = $1; }
           | IDENTIFIER         { $$ = $1; }
;

/********************************************************************
 * ATOM SELECTION EXPRESSIONS
 ********************************************************************/

/* Boolean expressions and grouping */
sel_expr:    NOT sel_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selelem_create(SEL_BOOLEAN);
                 $$->u.boolt = BOOL_NOT;
                 $$->child = $2;
                 END_ACTION;
             }
           | sel_expr AND sel_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selelem_create(SEL_BOOLEAN);
                 $$->u.boolt = BOOL_AND;
                 $$->child = $1; $$->child->next = $3;
                 END_ACTION;
             }
           | sel_expr OR  sel_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selelem_create(SEL_BOOLEAN);
                 $$->u.boolt = BOOL_OR;
                 $$->child = $1; $$->child->next = $3;
                 END_ACTION;
             }
           | '(' sel_expr ')'   { $$ = $2; }
;

/* Numeric comparisons */
sel_expr:    num_expr CMP_OP num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_comparison($1, $3, $2, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/* External groups */
sel_expr:    GROUP string
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_group_by_name($2, scanner);
                 free($2);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
           | GROUP TOK_INT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_group_by_id($2, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/* Position modifiers for selection methods */
pos_mod:     EMPTY_POSMOD       { $$ = NULL; }
           | KEYWORD_POS        { $$ = $1;   }
;

/* Keyword selections */
sel_expr:    pos_mod KEYWORD_GROUP
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_keyword($2, NULL, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
           | pos_mod KEYWORD_STR basic_value_list
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_keyword($2, process_value_list($3, NULL), $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
           | pos_mod KEYWORD_NUMERIC basic_value_list
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_keyword($2, process_value_list($3, NULL), $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/* Custom selection methods */
sel_expr:    pos_mod METHOD_GROUP method_params
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_method($2, $3, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/********************************************************************
 * NUMERICAL EXPRESSIONS
 ********************************************************************/

/* Basic numerical values */
num_expr:    TOK_INT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selelem_create(SEL_CONST);
                 _gmx_selelem_set_vtype($$, INT_VALUE);
                 _gmx_selvalue_reserve(&$$->v, 1);
                 $$->v.u.i[0] = $1;
                 END_ACTION;
             }
           | TOK_REAL
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selelem_create(SEL_CONST);
                 _gmx_selelem_set_vtype($$, REAL_VALUE);
                 _gmx_selvalue_reserve(&$$->v, 1);
                 $$->v.u.r[0] = $1;
                 END_ACTION;
             }
;

/* Numeric selection methods */
num_expr:    pos_mod KEYWORD_NUMERIC    %prec NUM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_keyword($2, NULL, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
           | pos_mod METHOD_NUMERIC method_params
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_method($2, $3, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/* Arithmetic evaluation and grouping */
num_expr:    num_expr '+' num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_arithmetic($1, $3, '+', scanner);
                 END_ACTION;
             }
           | num_expr '-' num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_arithmetic($1, $3, '-', scanner);
                 END_ACTION;
             }
           | num_expr '*' num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_arithmetic($1, $3, '*', scanner);
                 END_ACTION;
             }
           | num_expr '/' num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_arithmetic($1, $3, '/', scanner);
                 END_ACTION;
             }
           | '-' num_expr %prec UNARY_NEG
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_arithmetic($2, NULL, '-', scanner);
                 END_ACTION;
             }
           | num_expr '^' num_expr
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_arithmetic($1, $3, '^', scanner);
                 END_ACTION;
             }
           | '(' num_expr ')'   { $$ = $2; }
;

/********************************************************************
 * STRING EXPRESSIONS
 ********************************************************************/

str_expr:    string
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selelem_create(SEL_CONST);
                 _gmx_selelem_set_vtype($$, STR_VALUE);
                 _gmx_selvalue_reserve(&$$->v, 1);
                 $$->v.u.s[0] = $1;
                 END_ACTION;
             }
           | pos_mod KEYWORD_STR
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_keyword($2, NULL, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/********************************************************************
 * POSITION EXPRESSIONS
 ********************************************************************/

/* Constant position expressions */
pos_expr:    '[' number ',' number ',' number ']'
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_const_position($2, $4, $6);
                 END_ACTION;
             }
;

/* Grouping of position expressions */
pos_expr:    '(' pos_expr ')'   { $$ = $2; }
;

/* Expressions with a position value */
pos_expr:    METHOD_POS method_params
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_method($1, $2, NULL, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/* Evaluation of positions using a keyword */
pos_expr:    KEYWORD_POS OF sel_expr    %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_position($3, $1, scanner);
                 if ($$ == NULL) YYERROR;
                 END_ACTION;
             }
;

/********************************************************************
 * VARIABLES
 ********************************************************************/

sel_expr:    VARIABLE_GROUP
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_variable_ref($1);
                 END_ACTION;
             }
;

num_expr:    VARIABLE_NUMERIC
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_variable_ref($1);
                 END_ACTION;
             }
;

pos_expr:    VARIABLE_POS
             {
                 BEGIN_ACTION;
                 $$ = _gmx_sel_init_variable_ref($1);
                 END_ACTION;
             }
;

/********************************************************************
 * METHOD PARAMETERS
 ********************************************************************/

method_params:
             method_param_list
             { $$ = process_param_list($1); }
           | method_param_list END_OF_METHOD
             { $$ = process_param_list($1); }
;

method_param_list:
             /* empty */        { $$ = NULL;              }
           | method_param_list method_param
                                { $2->next = $1; $$ = $2; }
;

method_param:
             PARAM value_list
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_param($1);
                 $$->value = process_value_list($2, &$$->nval);
                 END_ACTION;
             }
;

value_list:  /* empty */                         { $$ = NULL; }
           | value_list_contents                 { $$ = $1;   }
           | '{' value_list_contents '}'         { $$ = $2;   }
;

value_list_contents:
             value_item          { $$ = $1; }
           | value_list_contents value_item
                                 { $2->next = $1; $$ = $2; }
           | value_list_contents ',' value_item
                                 { $3->next = $1; $$ = $3; }
;

basic_value_list:
             basic_value_list_contents           { $$ = $1; }
           | '{' basic_value_list_contents '}'   { $$ = $2; }
;

basic_value_list_contents:
             basic_value_item    { $$ = $1; }
           | basic_value_list_contents basic_value_item
                                 { $2->next = $1; $$ = $2; }
           | basic_value_list_contents ',' basic_value_item
                                 { $3->next = $1; $$ = $3; }
;

value_item:  sel_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value_expr($1);
                 END_ACTION;
             }
           | pos_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value_expr($1);
                 END_ACTION;
             }
           | num_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value_expr($1);
                 END_ACTION;
             }
           | str_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value_expr($1);
                 END_ACTION;
             }
           | value_item_range    { $$ = $1; }
;

basic_value_item:
             integer_number      %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(INT_VALUE);
                 $$->u.i.i1 = $$->u.i.i2 = $1;
                 END_ACTION;
             }
           | real_number         %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(REAL_VALUE);
                 $$->u.r.r1 = $$->u.r.r2 = $1;
                 END_ACTION;
             }
           | string              %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(STR_VALUE);
                 $$->u.s = $1;
                 END_ACTION;
             }
           | value_item_range    { $$ = $1; }
;

value_item_range:
             integer_number TO integer_number
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(INT_VALUE);
                 $$->u.i.i1 = $1; $$->u.i.i2 = $3;
                 END_ACTION;
             }
           | integer_number TO real_number
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(REAL_VALUE);
                 $$->u.r.r1 = $1; $$->u.r.r2 = $3;
                 END_ACTION;
             }
           | real_number TO number
             {
                 BEGIN_ACTION;
                 $$ = _gmx_selexpr_create_value(REAL_VALUE);
                 $$->u.r.r1 = $1; $$->u.r.r2 = $3;
                 END_ACTION;
             }
;

%%

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


