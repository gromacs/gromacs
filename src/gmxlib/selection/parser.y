/*
 * $Id$
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
/*! \internal \file parser.c
 * \brief Generated (from parser.y by Bison) parser for the selection language.
 */
/*! \internal \file parser.h
 * \brief Generated (from parser.y by Bison) parser include file.
 */
%{
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#include <futil.h>
#include <smalloc.h>
#include <vec.h>

#include <indexutil.h>
#include <position.h>
#include <selection.h>
#include <selmethod.h>

#include "keywords.h"
#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"
#include "symrec.h"

#include "scanner.h"

static void
yyerror(gmx_sel_lexer_t *, int, gmx_ana_indexgrps_t *, char const *s);

static t_selelem *
get_group_by_name(gmx_ana_indexgrps_t *grps, char *name);
static t_selelem *
get_group_by_id(gmx_ana_indexgrps_t *grps, int id);

static t_selelem *
append_selection(gmx_sel_lexer_t *scanner, t_selelem *sel, t_selelem *last);

static t_selelem *
assign_variable(gmx_sel_lexer_t *scanner, char *name, t_selelem *expr);
static t_selelem *
init_selection(gmx_sel_lexer_t *scanner, t_selelem *sel);

static t_selexpr_value *
process_value_list(t_selexpr_value *values, int *nr);
static t_selexpr_param *
process_param_list(t_selexpr_param *params);
static t_selelem *
init_keyword_expr(gmx_ana_selcollection_t *sc, gmx_ana_selmethod_t *method,
                  t_selexpr_value *values, char *rpost);
/*! \cond */
%}

%union{
    int                  i;
    real                 r;
    char                *str;
    gmx_ana_selmethod_t *meth;

    t_selelem        *sel;

    t_selexpr_value  *val;
    t_selexpr_param  *param;
}

/* Invalid token to report lexer errors */
%token INVALID

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

%destructor { free($$);                     } STR IDENTIFIER string
%destructor { if($$) _gmx_selelem_free($$); } command
%destructor { _gmx_selelem_free_chain($$);  } selection
%destructor { _gmx_selelem_free($$);        } sel_expr numeric_expr
%destructor { _gmx_selelem_free($$);        } pos_expr pos_expr_sel pos_expr_nosel pos_expr_nosel_impl
%destructor { _gmx_selexpr_free_values($$); } string_list int_list int_list_item
%destructor { _gmx_selexpr_free_params($$); } method_params method_param_list method_param

%expect 15
%debug
%pure-parser

%parse-param { gmx_sel_lexer_t         *scanner }
%lex-param   { gmx_sel_lexer_t         *scanner }
%parse-param { int                      nexp    }
%parse-param { gmx_ana_indexgrps_t     *grps    }

%%

/* The start rule: allow one or more commands separated by semicolons */
commands:    command
             {
                 $$ = append_selection(scanner, $1, NULL);
                 if (_gmx_sel_lexer_selcollection(scanner)->nr == nexp)
                     YYACCEPT;
             }
           | commands CMD_SEP command
             {
                 $$ = append_selection(scanner, $3, $1);
                 if (_gmx_sel_lexer_selcollection(scanner)->nr == nexp)
                     YYACCEPT;
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
                 gmx_ana_selcollection_t *sc;
                 t_selelem               *s, *p;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 s = get_group_by_id(grps, $1);
                 if (s == NULL) YYABORT;
                 p = _gmx_sel_init_position(sc, s, sc->spost, TRUE);
                 if (p == NULL) YYABORT;
                 $$ = init_selection(scanner, p);
             }
           | selection          { $$ = init_selection(scanner, $1);      }
           | string selection   { $$ = init_selection(scanner, $2);
                                  $$->name = $1; $$->u.cgrp.name = $1;   }
           | IDENTIFIER '=' sel_expr
                                { $$ = assign_variable(scanner, $1, $3); }
           | IDENTIFIER '=' numeric_expr
                                { $$ = assign_variable(scanner, $1, $3); }
           | IDENTIFIER '=' pos_expr_nosel
                                { $$ = assign_variable(scanner, $1, $3); }
;

/* Selection is made of an expression and zero or more modifiers */
selection:   pos_expr_sel       { $$ = $1; }
           | selection MODIFIER method_params
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_modifier(sc, $2, process_param_list($3), $1);
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
sel_expr:    GROUP string       { $$ = get_group_by_name(grps, $2);
                                  if ($$ == NULL) YYABORT;               }
           | GROUP INT          { $$ = get_group_by_id(grps, $2);
                                  if ($$ == NULL) YYABORT;               }
;

/* External variables */
sel_expr:    VARIABLE_GROUP     { $$ = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  $$->name   = $1->name;
                                  $$->v.type = $1->v.type;
                                  $$->child  = $1;
                                  $1->refcount++;                        }
;

numeric_expr:
             VARIABLE_NUMERIC   { $$ = _gmx_selelem_create(SEL_SUBEXPRREF);
                                  $$->name   = $1->name;
                                  $$->v.type = $1->v.type;
                                  $$->child  = $1;
                                  $1->refcount++;                        }
;

pos_expr:    VARIABLE_POS       { if ($1->type == SEL_CONST) {
                                      $$ = $1;
                                  } else {
                                      $$ = _gmx_selelem_create(SEL_SUBEXPRREF);
                                      $$->name   = $1->name;
                                      $$->v.type = $1->v.type;
                                      $$->child  = $1;
                                  }
                                  $1->refcount++;                        }
;

/* Keyword selections */
sel_expr:    pos_mod KEYWORD_GROUP
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = init_keyword_expr(sc, $2, NULL, $1);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod KEYWORD_STR string_list
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = init_keyword_expr(sc, $2, $3, $1);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod KEYWORD_INT int_list
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = init_keyword_expr(sc, $2, $3, $1);
                 if ($$ == NULL) YYABORT;
             }
;

/* Custom selection methods */
sel_expr:    pos_mod METHOD_GROUP method_params
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_method(sc, $2, process_param_list($3), $1);
                 if ($$ == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             }
;

/* Numeric selections */
sel_expr:    numeric_expr CMP_OP numeric_expr
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_comparison(sc, $1, $3, $2);
                 if ($$ == NULL) YYABORT;
             }
;

/* Expressions that can (and should) be compared numerically */
numeric_expr:
             INT                { $$ = _gmx_selelem_create(SEL_CONST);
                                  $$->v.type = INT_VALUE;
                                  snew($$->v.u.i, 1);
                                  $$->v.u.i[0] = $1;                  }
           | REAL               { $$ = _gmx_selelem_create(SEL_CONST);
                                  $$->v.type = REAL_VALUE;
                                  snew($$->v.u.r, 1);
                                  $$->v.u.r[0] = $1;                  }
           | pos_mod KEYWORD_INT
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = init_keyword_expr(sc, $2, NULL, $1);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod KEYWORD_REAL
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = init_keyword_expr(sc, $2, NULL, $1);
                 if ($$ == NULL) YYABORT;
             }
           | pos_mod METHOD_NUMERIC method_params
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_method(sc, $2, process_param_list($3), $1);
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
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_method(sc, $1, process_param_list($2), NULL);
                 if ($$ == NULL) YYABORT;
                 _gmx_sel_finish_method(scanner);
             }
;

/* Evaluation of selection output positions */
pos_expr_sel:
             pos_expr           { $$ = $1; }
           | sel_expr
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_position(sc, $1, sc->spost, TRUE);
                 if ($$ == NULL) YYABORT;
             }
           | KEYWORD_POS OF sel_expr
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_position(sc, $3, $1, TRUE);
                 if ($$ == NULL) YYABORT;
             }
;

/* Evaluation of positions somewhere else */
pos_expr_nosel:
             pos_expr           { $$ = $1; }
           | KEYWORD_POS OF sel_expr
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_position(sc, $3, $1, FALSE);
                 if ($$ == NULL) YYABORT;
             }
;

/* Evaluation of reference positions for a selection */
pos_expr_nosel_impl:
             pos_expr_nosel     { $$ = $1; }
           | sel_expr
             {
                 gmx_ana_selcollection_t *sc;
                 sc = _gmx_sel_lexer_selcollection(scanner);
                 $$ = _gmx_sel_init_position(sc, $1, sc->rpost, FALSE);
                 if ($$ == NULL) YYABORT;
             }
;

/* Constant position expressions */
pos_expr:    '(' number ',' number ',' number ')'
                                { rvec x;
                                  $$ = _gmx_selelem_create(SEL_CONST);
                                  $$->v.type = POS_VALUE;
                                  snew($$->v.u.p, 1);
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
/*! \endcond */

/*! \brief
 * Internal helper function used by parse_selection_*() to do the actual work.
 *
 * \param[in,out] scanner Scanner data structure.
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     grps  External index groups (can be NULL).
 * \param[in]     maxnr Maximum number of selections to parse
 *   (if -1, parse as many as provided by the user).
 * \returns       0 on success, -1 on error.
 */
static int
run_parser(gmx_sel_lexer_t *scanner, gmx_ana_selcollection_t *sc,
           gmx_ana_indexgrps_t *grps, int maxnr)
{
    bool bOk;
    int  nr;
    int  nexp;

    nr        = sc->nr;
    nexp      = (maxnr > 0) ? (sc->nr + maxnr) : -1;
    bOk = !yyparse(scanner, nexp, grps);
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

/*!
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     nr    Number of selections to parse
 *   (if -1, parse as many as provided by the user).
 * \param[in]     grps  External index groups (can be NULL).
 * \returns       0 on success, -1 on error.
 *
 * The number of selections parsed can be accessed with
 * gmx_ana_selcollection_get_count() (note that if you call the parser
 * multiple times, this function returns the total count).
 */
int
gmx_ana_selcollection_parse_stdin(gmx_ana_selcollection_t *sc, int nr,
                                  gmx_ana_indexgrps_t *grps)
{
    gmx_sel_lexer_t *scanner;

    _gmx_sel_init_lexer(&scanner, sc, isatty(fileno(stdin)));
    _gmx_sel_set_lex_input_file(scanner, stdin);
    return run_parser(scanner, sc, grps, nr);
}

/*!
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     fnm   Name of the file to parse selections from.
 * \param[in]     grps  External index groups (can be NULL).
 * \returns       0 on success, -1 on error.
 *
 * The number of selections parsed can be accessed with
 * gmx_ana_selcollection_get_count() (note that if you call the parser
 * multiple times, this function returns the total count).
 */
int
gmx_ana_selcollection_parse_file(gmx_ana_selcollection_t *sc, char *fnm,
                                 gmx_ana_indexgrps_t *grps)
{
    gmx_sel_lexer_t *scanner;
    FILE *fp;
    int   rc;

    _gmx_sel_init_lexer(&scanner, sc, FALSE);
    fp = ffopen(fnm, "r");
    _gmx_sel_set_lex_input_file(scanner, fp);
    rc = run_parser(scanner, sc, grps, -1);
    fclose(fp);
    return rc;
}

/*!
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     str   String to parse selections from.
 * \param[in]     grps  External index groups (can be NULL).
 * \returns       0 on success, -1 on error.
 *
 * The number of selections parsed can be accessed with
 * gmx_ana_selcollection_get_count() (note that if you call the parser
 * multiple times, this function returns the total count).
 */
int
gmx_ana_selcollection_parse_str(gmx_ana_selcollection_t *sc, char *str,
                                gmx_ana_indexgrps_t *grps)
{
    gmx_sel_lexer_t *scanner;

    _gmx_sel_init_lexer(&scanner, sc, FALSE);
    _gmx_sel_set_lex_input_str(scanner, str);
    return run_parser(scanner, sc, grps, -1);
}

static t_selelem *
get_group_by_name(gmx_ana_indexgrps_t *grps, char *name)
{
    t_selelem *sel;

    if (!grps)
    {
        sfree(name);
        return NULL;
    }
    sel = _gmx_selelem_create(SEL_CONST);
    sel->v.type = GROUP_VALUE;
    if (!gmx_ana_indexgrps_find(&sel->u.cgrp, grps, name))
    {
        sfree(name);
        _gmx_selelem_free(sel);
        return NULL;
    }
    sfree(name);
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
    sel->v.type = GROUP_VALUE;
    if (!gmx_ana_indexgrps_extract(&sel->u.cgrp, grps, id))
    {
        _gmx_selelem_free(sel);
        return NULL;
    }
    sel->name = sel->u.cgrp.name;
    return sel;
}

static t_selelem *
append_selection(gmx_sel_lexer_t *scanner, t_selelem *sel, t_selelem *last)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);

    if (last)
    {
        last->next = sel;
    }
    else
    {
        if (sc->root)
        {
            last = sc->root;
            while (last->next)
            {
                last = last->next;
            }
            last->next = sel;
        }
        else
        {
            sc->root = sel;
        }
    }
    if (sel)
    {
        last = sel;
        /* Add the new selection to the collection if it is not a variable. */
        if (sel->child->type != SEL_SUBEXPR)
        {
            int        i;

            sc->nr++;
            srenew(sc->sel, sc->nr);
            i = sc->nr - 1;
            snew(sc->sel[i], 1);

            if (sel->child->type == SEL_CONST)
            {
                gmx_ana_pos_copy(&sc->sel[i]->p, sel->child->v.u.p, TRUE);
                sc->sel[i]->bDynamic = FALSE;
            }
            else
            {
                t_selelem *child;

                child = sel->child;
                child->flags     &= ~SEL_ALLOCVAL;
                child->v.u.p      = &sc->sel[i]->p;
                /* We should also skip any modifiers to determine the dynamic
                 * status. */
                while (child->type == SEL_MODIFIER)
                {
                    child = child->child;
                }
                /* For variable references, we should skip the
                 * SEL_SUBEXPRREF and SEL_SUBEXPR elements. */
                if (child->type == SEL_SUBEXPRREF)
                {
                    child = child->child->child;
                }
                sc->sel[i]->bDynamic = (child->child->flags & SEL_DYNAMIC);
            }
            /* The group will be set after compilation */
            sc->sel[i]->g        = NULL;
            sc->sel[i]->selelem  = sel;
            gmx_ana_selection_init_coverfrac(sc->sel[i], CFRAC_NONE);
        }
    }
    return last;
}

static t_selelem *
assign_variable(gmx_sel_lexer_t *scanner, char *name, t_selelem *expr)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem   *root;
    int          rc;

    rc = _gmx_selelem_update_flags(expr);
    if (rc != 0)
    {
        sfree(name);
        _gmx_selelem_free(expr);
        return NULL;
    }
    /* Check if this is a constant non-group value */
    if (expr->type == SEL_CONST && expr->v.type != GROUP_VALUE)
    {
        /* If so, just assign the constant value to the variable */
        if (!_gmx_sel_add_var_symbol(sc->symtab, name, expr))
        {
            _gmx_selelem_free(expr);
            sfree(name);
            return NULL;
        }
        _gmx_selelem_free(expr);
        if (_gmx_sel_is_lexer_interactive(scanner))
        {
            fprintf(stderr, "Variable '%s' parsed\n", name);
        }
        sfree(name);
        return NULL;
    }
    /* Check if we are assigning a variable to another variable */
    if (expr->type == SEL_SUBEXPRREF)
    {
        /* If so, make a simple alias */
        if (!_gmx_sel_add_var_symbol(sc->symtab, name, expr->child))
        {
            _gmx_selelem_free(expr);
            sfree(name);
            return NULL;
        }
        _gmx_selelem_free(expr);
        if (_gmx_sel_is_lexer_interactive(scanner))
        {
            fprintf(stderr, "Variable '%s' parsed\n", name);
        }
        sfree(name);
        return NULL;
    }
    root = _gmx_selelem_create(SEL_ROOT);
    root->name          = name;
    root->u.cgrp.name   = name;
    root->child = _gmx_selelem_create(SEL_SUBEXPR);
    root->child->name   = name;
    root->child->v.type = expr->v.type;
    root->child->child  = expr;
    rc = _gmx_selelem_update_flags(root);
    if (rc != 0)
    {
        _gmx_selelem_free(root);
        return NULL;
    }
    if (!_gmx_sel_add_var_symbol(sc->symtab, name, root->child))
    {
        _gmx_selelem_free(root);
        return NULL;
    }
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        fprintf(stderr, "Variable '%s' parsed\n", name);
    }
    return root;
}

static t_selelem *
init_selection(gmx_sel_lexer_t *scanner, t_selelem *sel)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem               *root, *child;
    int                      rc;

    if (sel->v.type != POS_VALUE)
    {
        gmx_bug("each selection must evaluate to a position");
        /* FIXME: Better handling of this error */
        return NULL;
    }

    root = _gmx_selelem_create(SEL_ROOT);
    root->child = sel;
    /* Update the flags */
    rc = _gmx_selelem_update_flags(root);
    if (rc != 0)
    {
        _gmx_selelem_free(root);
        return NULL;
    }
    /* Set the reference position type if applicable */
    child = sel;
    while (child->type == SEL_MODIFIER)
    {
        /* We skip the modifier and the SEL_SUBEXPRREF after it */
        child = child->child->child;
    }

    /* Print out some information if the parser is interactive */
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        /* TODO: It would be nice to print the whole selection here */
        fprintf(stderr, "Selection parsed\n");
    }

    return root;
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
yyerror(gmx_sel_lexer_t *scanner, int nexp, gmx_ana_indexgrps_t *grps,
        char const *s)
{
    _gmx_selparser_error("%s", s);
}
