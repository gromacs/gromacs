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
 */
%{
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
yyerror(int, gmx_ana_indexgrps_t *, yyscan_t, char const *s);

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
static t_selexpr_param *
process_param_list(t_selexpr_param *params);
%}

%union{
    int                  i;
    real                 r;
    char                *str;
    struct gmx_ana_selmethod_t *meth;

    struct t_selelem           *sel;

    struct t_selexpr_value     *val;
    struct t_selexpr_param     *param;
}

/* Invalid token to report lexer errors */
%token INVALID

/* Tokens for help requests */
%token         HELP
%token <str>   HELP_TOPIC

/* Simple input tokens */
%token <i>     INT
%token <r>     REAL
%token <str>   STR
%token <str>   IDENTIFIER
%token         CMD_SEP

/* Simple keyword tokens */
%token GROUP
%token TO
%token OF

%token <sel>   VARIABLE_NUMERIC
%token <sel>   VARIABLE_GROUP
%token <sel>   VARIABLE_POS
%token <meth>  KEYWORD_INT
%token <meth>  KEYWORD_REAL
%token <meth>  KEYWORD_STR
%token <str>   KEYWORD_POS
%token <meth>  KEYWORD_GROUP
%token <meth>  METHOD_NUMERIC
%token <meth>  METHOD_GROUP
%token <meth>  METHOD_POS
%token <meth>  MODIFIER

%token <str>   PARAM_BOOL
%token <i>     BOOL_VALUE
%token <str>   PARAM_INT
%token <str>   PARAM_REAL
%token <str>   PARAM_STR
%token <str>   PARAM_POS
%token <str>   PARAM_GROUP
%token         END_OF_METHOD

/* Operator tokens */
%left           AND OR XOR
%left           NOT
%nonassoc <str> CMP_OP

/* Simple non-terminals */
%type <r>     number
%type <str>   string
%type <str>   pos_mod

/* Expression non-terminals */
%type <sel>   commands command
%type <sel>   selection
%type <sel>   sel_expr
%type <sel>   numeric_expr
%type <sel>   pos_expr pos_expr_sel pos_expr_nosel pos_expr_nosel_impl

/* Parameter/value non-terminals */
%type <val>   string_list
%type <val>   int_list     int_list_item
%type <param> method_params method_param_list method_param

%destructor { free($$);                     } HELP_TOPIC STR IDENTIFIER string
%destructor { if($$) _gmx_selelem_free($$); } command
%destructor { _gmx_selelem_free_chain($$);  } selection
%destructor { _gmx_selelem_free($$);        } sel_expr numeric_expr
%destructor { _gmx_selelem_free($$);        } pos_expr pos_expr_sel pos_expr_nosel pos_expr_nosel_impl
%destructor { _gmx_selexpr_free_values($$); } string_list int_list int_list_item
%destructor { _gmx_selexpr_free_params($$); } method_params method_param_list method_param

%expect 15
%debug
%pure-parser

/* If you change these, you also need to update the prototype in parsetree.c. */
%name-prefix="_gmx_sel_yy"
%parse-param { int                      nexp    }
%parse-param { gmx_ana_indexgrps_t     *grps    }
%parse-param { yyscan_t                 scanner }
%lex-param   { yyscan_t                 scanner }

%%

/* The start rule: allow one or more commands separated by semicolons */
commands:    command
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_append_selection($1, NULL, scanner);
                 if (sc->nr == nexp)
                     YYACCEPT;
                 _gmx_sel_lexer_clear_pselstr(scanner);
             }
           | commands CMD_SEP command
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_append_selection($3, $1, scanner);
                 if (sc->nr == nexp)
                     YYACCEPT;
                 _gmx_sel_lexer_clear_pselstr(scanner);
             }
;

/* Basic expressions */
number:      INT                { $$ = $1; }
           | REAL               { $$ = $1; }
;

string:      STR                { $$ = $1; }
           | IDENTIFIER         { $$ = $1; }
;

pos_mod:     /* empty */        { $$ = NULL; }
           | KEYWORD_POS        { $$ = $1;   }
;

/* Commands can be selections or variable assignments */
command:     /* empty */        { $$ = NULL;                            }
           | INT
             {
                 t_selelem               *s, *p;
                 s = _gmx_sel_init_group_by_id(grps, $1);
                 if (s == NULL) YYABORT;
                 p = _gmx_sel_init_position(s, NULL, TRUE, scanner);
                 if (p == NULL) YYABORT;
                 $$ = _gmx_sel_init_selection(strdup(s->name), p, scanner);
             }
           | selection          { $$ = _gmx_sel_init_selection(NULL, $1, scanner); }
           | string selection   { $$ = _gmx_sel_init_selection($1, $2, scanner);  }
           | IDENTIFIER '=' sel_expr
                                { $$ = _gmx_sel_assign_variable($1, $3, scanner); }
           | IDENTIFIER '=' numeric_expr
                                { $$ = _gmx_sel_assign_variable($1, $3, scanner); }
           | IDENTIFIER '=' pos_expr_nosel
                                { $$ = _gmx_sel_assign_variable($1, $3, scanner); }
           | help_request       { $$ = NULL; }
;

/* Help requests */
help_request:
             HELP
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 _gmx_sel_print_help(sc, NULL);
             }
           | help_topic

help_topic:  HELP HELP_TOPIC
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 _gmx_sel_print_help(sc, $2);
                 sfree($2);
             }
           | help_topic HELP_TOPIC
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 _gmx_sel_print_help(sc, $2);
                 sfree($2);
             }

/* Selection is made of an expression and zero or more modifiers */
selection:   pos_expr_sel       { $$ = $1; }
           | selection MODIFIER method_params
             {
                 $$ = _gmx_sel_init_modifier($2, process_param_list($3), $1, scanner);
                 _gmx_sel_finish_method(scanner);
             }
           | '(' selection ')'  { $$ = $2;                              }
;

/* Boolean expressions and grouping */
sel_expr:    NOT sel_expr       { $$ = _gmx_selelem_create(SEL_BOOLEAN);
                                  $$->u.boolt = BOOL_NOT;
                                  $$->child = $2;                       }
           | sel_expr AND sel_expr
                                { $$ = _gmx_selelem_create(SEL_BOOLEAN);
                                  $$->u.boolt = BOOL_AND;
                                  $$->child = $1; $$->child->next = $3; }
           | sel_expr OR  sel_expr
                                { $$ = _gmx_selelem_create(SEL_BOOLEAN);
                                  $$->u.boolt = BOOL_OR;
                                  $$->child = $1; $$->child->next = $3; }
/*           | sel_expr XOR sel_expr
                                { $$ = _gmx_selelem_create(SEL_BOOLEAN);
                                  $$->u.boolt = BOOL_XOR;
                                  $$->child = $1; $$->child->next = $3; }*/
           | '(' sel_expr ')'   { $$ = $2;                               }
;

/* External groups */
sel_expr:    GROUP string
             {
                 $$ = _gmx_sel_init_group_by_name(grps, $2);
                 sfree($2);
                 if ($$ == NULL) YYABORT;
             }
           | GROUP INT
             {
                 $$ = _gmx_sel_init_group_by_id(grps, $2);
                 if ($$ == NULL) YYABORT;
             }
;

/* External variables */
sel_expr:    VARIABLE_GROUP     { $$ = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  _gmx_selelem_set_vtype($$, $1->v.type);
                                  $$->name   = $1->name;
                                  $$->child  = $1;
                                  $1->refcount++;                        }
;

numeric_expr:
             VARIABLE_NUMERIC   { $$ = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  _gmx_selelem_set_vtype($$, $1->v.type);
                                  $$->name   = $1->name;
                                  $$->child  = $1;
                                  $1->refcount++;                        }
;

pos_expr:    VARIABLE_POS       { if ($1->type == SEL_CONST) {
                                      $$ = $1;
                                  } else {
                                      $$ = _gmx_selelem_create(SEL_SUBEXPRREF);
                                      _gmx_selelem_set_vtype($$, $1->v.type);
                                      $$->name   = $1->name;
                                      $$->child  = $1;
                                  }
                                  $1->refcount++;                        }
;

/* Keyword selections */
sel_expr:    pos_mod KEYWORD_GROUP
             {
                 $$ = _gmx_sel_init_keyword($2, NULL, $1, scanner);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod KEYWORD_STR string_list
             {
                 $$ = _gmx_sel_init_keyword($2, process_value_list($3, NULL), $1, scanner);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod KEYWORD_INT int_list
             {
                 $$ = _gmx_sel_init_keyword($2, process_value_list($3, NULL), $1, scanner);
                 if ($$ == NULL) YYABORT;
             }
;

/* Custom selection methods */
sel_expr:    pos_mod METHOD_GROUP method_params
             {
                 $$ = _gmx_sel_init_method($2, process_param_list($3), $1, scanner);
                 if ($$ == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             }
;

/* Numeric selections */
sel_expr:    numeric_expr CMP_OP numeric_expr
             {
                 $$ = _gmx_sel_init_comparison($1, $3, $2, scanner);
                 if ($$ == NULL) YYABORT;
             }
;

/* Expressions that can (and should) be compared numerically */
numeric_expr:
             INT                { $$ = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype($$, INT_VALUE);
                                  _gmx_selvalue_reserve(&$$->v, 1);
                                  $$->v.u.i[0] = $1;                  }
           | REAL               { $$ = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype($$, REAL_VALUE);
                                  _gmx_selvalue_reserve(&$$->v, 1);
                                  $$->v.u.r[0] = $1;                  }
           | pos_mod KEYWORD_INT
             {
                 $$ = _gmx_sel_init_keyword($2, NULL, $1, scanner);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod KEYWORD_REAL
             {
                 $$ = _gmx_sel_init_keyword($2, NULL, $1, scanner);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod METHOD_NUMERIC method_params
             {
                 $$ = _gmx_sel_init_method($2, process_param_list($3), $1, scanner);
                 if ($$ == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             }
;

/* Grouping of numeric expressions */
numeric_expr:
             '(' numeric_expr ')'
                                { $$ = $2;                               }
;

/* Expressions with a position value */
pos_expr:    METHOD_POS method_params
             {
                 $$ = _gmx_sel_init_method($1, process_param_list($2), NULL, scanner);
                 if ($$ == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             }
;

/* Evaluation of selection output positions */
pos_expr_sel:
             pos_expr           { $$ = $1; }
           | KEYWORD_POS OF sel_expr
             {
                 $$ = _gmx_sel_init_position($3, $1, TRUE, scanner);
                 if ($$ == NULL) YYABORT;
             }
           | sel_expr
             {
                 $$ = _gmx_sel_init_position($1, NULL, TRUE, scanner);
                 if ($$ == NULL) YYABORT;
             }
;

/* Evaluation of positions somewhere else */
pos_expr_nosel:
             pos_expr           { $$ = $1; }
           | KEYWORD_POS OF sel_expr
             {
                 $$ = _gmx_sel_init_position($3, $1, FALSE, scanner);
                 if ($$ == NULL) YYABORT;
             }
;

/* Evaluation of reference positions for a selection */
pos_expr_nosel_impl:
             pos_expr_nosel     { $$ = $1; }
           | sel_expr
             {
                 $$ = _gmx_sel_init_position($1, NULL, FALSE, scanner);
                 if ($$ == NULL) YYABORT;
             }
;

/* Constant position expressions */
pos_expr:    '(' number ',' number ',' number ')'
                                { rvec x;
                                  $$ = _gmx_selelem_create(SEL_CONST);
                                  _gmx_selelem_set_vtype($$, POS_VALUE);
                                  _gmx_selvalue_reserve(&$$->v, 1);
                                  x[XX] = $2; x[YY] = $4; x[ZZ] = $6;
                                  gmx_ana_pos_init_const($$->v.u.p, x);  }
;

/* Grouping of position expressions */
pos_expr:    '(' pos_expr ')'   { $$ = $2;                               }
;

string_list:
             string             { $$ = _gmx_selexpr_create_value(STR_VALUE);
                                  $$->u.s = $1;                          }
           | string_list string { $$ = _gmx_selexpr_create_value(STR_VALUE);
                                  $$->u.s = $2; $$->next = $1;           }
;

int_list:
             int_list_item      { $$ = $1;                               }
           | int_list int_list_item
                                { $2->next = $1; $$ = $2;                }
;

int_list_item:
             INT                { $$ = _gmx_selexpr_create_value(INT_VALUE);
                                  $$->u.i.i1 = $$->u.i.i2 = $1;          }
           | INT TO INT         { $$ = _gmx_selexpr_create_value(INT_VALUE);
                                  $$->u.i.i1 = $1; $$->u.i.i2 = $3;      }
;

method_params:
             method_param_list  { $$ = $1; }
           | method_param_list END_OF_METHOD
                                { $$ = $1; }
;

method_param_list:
             /* empty */        { $$ = NULL;              }
           | method_param_list method_param
                                { $2->next = $1; $$ = $2; }
;

method_param:
             PARAM_BOOL         { $$ = _gmx_selexpr_create_param($1);    }
           | PARAM_BOOL BOOL_VALUE
             {
                 $$ = _gmx_selexpr_create_param($1);
                 $$->value = _gmx_selexpr_create_value(INT_VALUE);
                 $$->value->u.i.i1 = $$->value->u.i.i2 = $2;
             }
           | PARAM_INT  int_list
                                { $$ = _gmx_selexpr_create_param($1);
                                  $$->value = $2;                        }
           | PARAM_REAL number  { $$ = _gmx_selexpr_create_param($1);
                                  $$->value = _gmx_selexpr_create_value(REAL_VALUE);
                                  $$->value->u.r = $2;                   }
           | PARAM_STR  string  { $$ = _gmx_selexpr_create_param($1);
                                  $$->value = _gmx_selexpr_create_value(STR_VALUE);
                                  $$->value->u.s = $2;                   }
           | PARAM_POS  pos_expr_nosel_impl
                                { $$ = _gmx_selexpr_create_param($1);
                                  $$->value = _gmx_selexpr_create_value_expr($2); }
           | PARAM_GROUP sel_expr
                                { $$ = _gmx_selexpr_create_param($1);
                                  $$->value = _gmx_selexpr_create_value_expr($2); }
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
yyerror(int nexp, gmx_ana_indexgrps_t *grps, yyscan_t scanner,
        char const *s)
{
    _gmx_selparser_error("%s", s);
}
