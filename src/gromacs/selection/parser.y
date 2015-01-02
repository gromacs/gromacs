%code requires {
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
}
/*! \internal \file
 * \brief Grammar description and parser for the selection language.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
%code top {
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
#include "gmxpre.h"
}
%{
#include "gromacs/utility/scoped_cptr.h"

#include "parser_internal.h"

using gmx::scoped_guard_sfree;
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
%}

%code requires{
#include "parsetree.h"
#include "selelem.h"

#define YYLTYPE ::gmx::SelectionLocation
}

%union{
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
};

/* Invalid token to report lexer errors */
%token INVALID

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
%type <smt>   str_match_type

/* Expression non-terminals */
%type <sel>   commands command cmd_plain
%type <sel>   selection
%type <sel>   sel_expr
%type <sel>   num_expr
%type <sel>   str_expr
%type <sel>   pos_expr

/* Parameter/value non-terminals */
%type <plist> method_params method_param_list
%type <param> method_param
%type <vlist> value_list value_list_contents basic_value_list basic_value_list_contents
%type <val>   value_item value_item_range basic_value_item

%destructor { free($$);        } STR IDENTIFIER KEYWORD_POS CMP_OP string
%destructor { if($$) free($$); } PARAM pos_mod
%destructor { delete $$;       } commands command cmd_plain selection
%destructor { delete $$;       } sel_expr num_expr str_expr pos_expr
%destructor { delete $$;       } method_params method_param_list method_param
%destructor { delete $$;       } value_list value_list_contents basic_value_list basic_value_list_contents
%destructor { delete $$;       } value_item value_item_range basic_value_item

%expect 35
%debug
%pure-parser
%define api.push-pull push
%locations

%name-prefix="_gmx_sel_yy"
%parse-param { void *scanner }

%%

/* The start rule: allow one or more commands */
commands:    /* empty */
             {
                 BEGIN_ACTION;
                 set_empty($$);
                 END_ACTION_TOPLEVEL;
             }
           | commands command
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_append_selection(get($2), get($1), scanner));
                 if (_gmx_sel_parser_should_finish(scanner))
                     YYACCEPT;
                 END_ACTION_TOPLEVEL;
             }
;

/* A command is formed from an actual command and a separator */
command:     cmd_plain CMD_SEP  { $$ = $1; }
           | error CMD_SEP
             {
                 BEGIN_ACTION;
                 _gmx_sel_lexer_clear_method_stack(scanner);
                 if (_gmx_selparser_handle_error(scanner))
                 {
                     yyerrok;
                 }
                 else
                 {
                     YYABORT;
                 }
                 _gmx_sel_lexer_clear_pselstr(scanner);
                 set_empty($$);
                 END_ACTION_TOPLEVEL;
             }
;

/* Commands can be selections or variable assignments */
cmd_plain:   /* empty */
             {
                 BEGIN_ACTION;
                 set_empty($$);
                 END_ACTION;
             }
           | TOK_INT
             {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer s
                        = _gmx_sel_init_group_by_id($1, scanner);
                 SelectionTreeElementPointer p
                        = _gmx_sel_init_position(s, NULL, scanner);
                 if (!p) YYERROR;
                 set($$, _gmx_sel_init_selection(NULL, p, scanner));
                 END_ACTION;
             }
           | string
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($1);
                 SelectionTreeElementPointer s
                        = _gmx_sel_init_group_by_name($1, scanner);
                 SelectionTreeElementPointer p
                        = _gmx_sel_init_position(s, NULL, scanner);
                 if (!p) YYERROR;
                 set($$, _gmx_sel_init_selection(NULL, p, scanner));
                 END_ACTION;
             }
           | selection
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_selection(NULL, get($1), scanner));
                 END_ACTION;
             }
           | STR selection
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($1);
                 set($$, _gmx_sel_init_selection($1, get($2), scanner));
                 END_ACTION;
             }
           | IDENTIFIER '=' sel_expr
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($1);
                 set($$, _gmx_sel_assign_variable($1, get($3), scanner));
                 END_ACTION;
             }
           | IDENTIFIER '=' num_expr
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($1);
                 set($$, _gmx_sel_assign_variable($1, get($3), scanner));
                 END_ACTION;
             }
           | IDENTIFIER '=' pos_expr
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($1);
                 set($$, _gmx_sel_assign_variable($1, get($3), scanner));
                 END_ACTION;
             }
;

/* Selection is made of an expression and zero or more modifiers */
selection:   pos_expr           { $$ = $1; }
           | sel_expr
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_position(get($1), NULL, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
           | '(' selection ')'  { $$ = $2; }
           | selection MODIFIER method_params
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_modifier($2, get($3), get($1), scanner));
                 CHECK_SEL($$);
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
                 SelectionTreeElementPointer arg(get($2));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN, @$));
                 sel->u.boolt = BOOL_NOT;
                 sel->child = arg;
                 set($$, sel);
                 END_ACTION;
             }
           | sel_expr AND sel_expr
             {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg1(get($1)), arg2(get($3));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN, @$));
                 sel->u.boolt = BOOL_AND;
                 sel->child = arg1; sel->child->next = arg2;
                 set($$, sel);
                 END_ACTION;
             }
           | sel_expr OR  sel_expr
             {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer arg1(get($1)), arg2(get($3));
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_BOOLEAN, @$));
                 sel->u.boolt = BOOL_OR;
                 sel->child = arg1; sel->child->next = arg2;
                 set($$, sel);
                 END_ACTION;
             }
           | '(' sel_expr ')'   { $$ = $2; }
;

/* Numeric comparisons */
sel_expr:    num_expr CMP_OP num_expr
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree opGuard($2);
                 set($$, _gmx_sel_init_comparison(get($1), get($3), $2, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
;

/* External groups */
sel_expr:    GROUP string
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($2);
                 set($$, _gmx_sel_init_group_by_name($2, scanner));
                 END_ACTION;
             }
           | GROUP TOK_INT
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_group_by_id($2, scanner));
                 END_ACTION;
             }
;

/* Position modifiers for selection methods */
pos_mod:     EMPTY_POSMOD       { $$ = NULL; }
           | KEYWORD_POS        { $$ = $1;   }
;

/* Matching mode forcing for keyword matching */
str_match_type:
             '~'                { $$ = gmx::eStringMatchType_RegularExpression; }
           | '?'                { $$ = gmx::eStringMatchType_Wildcard; }
           | '='                { $$ = gmx::eStringMatchType_Exact; }
;

/* Keyword selections */
sel_expr:    pos_mod KEYWORD_GROUP
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword($2, SelectionParserValueListPointer(), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
           | pos_mod KEYWORD_STR basic_value_list
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword_strmatch($2, gmx::eStringMatchType_Auto, get($3), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
           | pos_mod KEYWORD_STR str_match_type basic_value_list
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword_strmatch($2, $3, get($4), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
           | pos_mod KEYWORD_NUMERIC basic_value_list
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword($2, get($3), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
;

/* Custom selection methods */
sel_expr:    pos_mod METHOD_GROUP method_params
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_method($2, get($3), $1, scanner));
                 CHECK_SEL($$);
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
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST, @$));
                 _gmx_selelem_set_vtype(sel, INT_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.i[0] = $1;
                 set($$, sel);
                 END_ACTION;
             }
           | TOK_REAL
             {
                 BEGIN_ACTION;
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST, @$));
                 _gmx_selelem_set_vtype(sel, REAL_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.r[0] = $1;
                 set($$, sel);
                 END_ACTION;
             }
;

/* Numeric selection methods */
num_expr:    pos_mod KEYWORD_NUMERIC    %prec NUM_REDUCT
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword($2, SelectionParserValueListPointer(), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
           | pos_mod KEYWORD_NUMERIC OF pos_expr
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword_of($2, get($4), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
           | pos_mod METHOD_NUMERIC method_params
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_method($2, get($3), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
;

/* Arithmetic evaluation and grouping */
num_expr:    num_expr '+' num_expr
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_arithmetic(get($1), get($3), '+', scanner));
                 END_ACTION;
             }
           | num_expr '-' num_expr
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_arithmetic(get($1), get($3), '-', scanner));
                 END_ACTION;
             }
           | num_expr '*' num_expr
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_arithmetic(get($1), get($3), '*', scanner));
                 END_ACTION;
             }
           | num_expr '/' num_expr
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_arithmetic(get($1), get($3), '/', scanner));
                 END_ACTION;
             }
           | '-' num_expr %prec UNARY_NEG
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_arithmetic(get($2), SelectionTreeElementPointer(), '-', scanner));
                 END_ACTION;
             }
           | num_expr '^' num_expr
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_arithmetic(get($1), get($3), '^', scanner));
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
                 SelectionTreeElementPointer sel(
                        new SelectionTreeElement(SEL_CONST, @$));
                 _gmx_selelem_set_vtype(sel, STR_VALUE);
                 _gmx_selvalue_reserve(&sel->v, 1);
                 sel->v.u.s[0] = $1;
                 set($$, sel);
                 END_ACTION;
             }
           | pos_mod KEYWORD_STR
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree posmodGuard($1);
                 set($$, _gmx_sel_init_keyword($2, SelectionParserValueListPointer(), $1, scanner));
                 CHECK_SEL($$);
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
                 set($$, _gmx_sel_init_const_position($2, $4, $6, scanner));
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
                 set($$, _gmx_sel_init_method($1, get($2), NULL, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
;

/* Evaluation of positions using a keyword */
pos_expr:    KEYWORD_POS OF sel_expr    %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree keywordGuard($1);
                 set($$, _gmx_sel_init_position(get($3), $1, scanner));
                 CHECK_SEL($$);
                 END_ACTION;
             }
;

/********************************************************************
 * VARIABLES
 ********************************************************************/

sel_expr:    VARIABLE_GROUP
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_variable_ref(get($1), scanner));
                 END_ACTION;
             }
;

num_expr:    VARIABLE_NUMERIC
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_variable_ref(get($1), scanner));
                 END_ACTION;
             }
;

pos_expr:    VARIABLE_POS
             {
                 BEGIN_ACTION;
                 set($$, _gmx_sel_init_variable_ref(get($1), scanner));
                 END_ACTION;
             }
;

/********************************************************************
 * METHOD PARAMETERS
 ********************************************************************/

method_params:
             method_param_list
             { $$ = $1; }
           | method_param_list END_OF_METHOD
             { $$ = $1; }
;

method_param_list:
             /* empty */
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserParameter::createList());
                 END_ACTION;
             }
           | method_param_list method_param
             {
                 BEGIN_ACTION;
                 SelectionParserParameterListPointer list(get($1));
                 list->push_back(get($2));
                 set($$, move(list));
                 END_ACTION;
             }
;

method_param:
             PARAM value_list
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree nameGuard($1);
                 set($$, SelectionParserParameter::create($1, get($2), @$));
                 END_ACTION;
             }
;

value_list:  value_list_contents                 { $$ = $1;   }
           | '{' value_list_contents '}'         { $$ = $2;   }
;

value_list_contents:
             /* empty */
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createList());
                 END_ACTION;
             }
           | value_list_contents value_item
             {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get($1));
                 list->push_back(get($2));
                 set($$, move(list));
                 END_ACTION;
             }
           | value_list_contents ',' value_item
             {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get($1));
                 list->push_back(get($3));
                 set($$, move(list));
                 END_ACTION;
             }
;

basic_value_list:
             basic_value_list_contents           { $$ = $1; }
           | '{' basic_value_list_contents '}'   { $$ = $2; }
;

basic_value_list_contents:
             basic_value_item
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createList(get($1)));
                 END_ACTION;
             }
           | basic_value_list_contents basic_value_item
             {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get($1));
                 list->push_back(get($2));
                 set($$, move(list));
                 END_ACTION;
             }
           | basic_value_list_contents ',' basic_value_item
             {
                 BEGIN_ACTION;
                 SelectionParserValueListPointer list(get($1));
                 list->push_back(get($3));
                 set($$, move(list));
                 END_ACTION;
             }
;

value_item:  sel_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createExpr(get($1)));
                 END_ACTION;
             }
           | pos_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createExpr(get($1)));
                 END_ACTION;
             }
           | num_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createExpr(get($1)));
                 END_ACTION;
             }
           | str_expr            %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createExpr(get($1)));
                 END_ACTION;
             }
           | value_item_range    { $$ = $1; }
;

basic_value_item:
             integer_number      %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createInteger($1, @$));
                 END_ACTION;
             }
           | real_number         %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createReal($1, @$));
                 END_ACTION;
             }
           | string              %prec PARAM_REDUCT
             {
                 BEGIN_ACTION;
                 scoped_guard_sfree stringGuard($1);
                 set($$, SelectionParserValue::createString($1, @$));
                 END_ACTION;
             }
           | value_item_range    { $$ = $1; }
;

value_item_range:
             integer_number TO integer_number
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createIntegerRange($1, $3, @$));
                 END_ACTION;
             }
           | integer_number TO real_number
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createRealRange($1, $3, @$));
                 END_ACTION;
             }
           | real_number TO number
             {
                 BEGIN_ACTION;
                 set($$, SelectionParserValue::createRealRange($1, $3, @$));
                 END_ACTION;
             }
;
