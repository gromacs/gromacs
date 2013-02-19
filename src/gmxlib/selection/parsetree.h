/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
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
/*! \internal \file
 * \brief Handling of intermediate selection parser data.
 *
 * The data types declared in this header are used by the parser to store
 * intermediate data when constructing method expressions.
 * In particular, the parameters for the method are stored.
 * The intermediate data is freed once a \c t_selelem object can be
 * constructed.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef SELECTION_PARSETREE_H
#define SELECTION_PARSETREE_H

/*#include <typedefs.h>*/
#include <types/simple.h>


#include <selvalue.h>

struct t_selelem;
struct gmx_ana_indexgrps_t;
struct gmx_ana_selmethod_t;
struct gmx_ana_selparam_t;

/*! \internal \brief
 * Describes a parsed value, possibly resulting from expression evaluation.
 */
typedef struct t_selexpr_value
{
    /** Type of the value. */
    e_selvalue_t                type;
    /** TRUE if the value is the result of an expression. */
    gmx_bool                    bExpr;
    union {
        /** The integer value/range (\p type INT_VALUE); */
        struct {
            /** Beginning of the range. */
            int             i1;
            /** End of the range; equals \p i1 for a single integer. */
            int             i2;
        }                   i;
        /** The real value/range (\p type REAL_VALUE); */
        struct {
            /** Beginning of the range. */
            real            r1;
            /** End of the range; equals \p r1 for a single number. */
            real            r2;
        }                   r;
        /** The string value (\p type STR_VALUE); */
        char               *s;
        /** The position value (\p type POS_VALUE); */
        rvec                x;
        /** The expression if \p bExpr is TRUE. */
        struct t_selelem   *expr;
    }                       u;
    /** Pointer to the next value. */
    struct t_selexpr_value *next;
} t_selexpr_value;

/*! \internal \brief
 * Describes a parsed method parameter.
 */
typedef struct t_selexpr_param
{
    /** Name of the parameter. */
    char                   *name;
    /** Number of values given for this parameter. */
    int                     nval;
    /** Pointer to the first value. */
    struct t_selexpr_value *value;
    /** Pointer to the next parameter. */
    struct t_selexpr_param *next;
} t_selexpr_param;

/** Error reporting function for the selection parser. */
void
_gmx_selparser_error(const char *fmt, ...);

/** Allocates and initializes a constant \c t_selexpr_value. */
t_selexpr_value *
_gmx_selexpr_create_value(e_selvalue_t type);
/** Allocates and initializes an expression \c t_selexpr_value. */
t_selexpr_value *
_gmx_selexpr_create_value_expr(struct t_selelem *expr);
/** Allocates and initializes a \c t_selexpr_param. */
t_selexpr_param *
_gmx_selexpr_create_param(char *name);

/** Frees the memory allocated for a chain of values. */
void
_gmx_selexpr_free_values(t_selexpr_value *value);
/** Frees the memory allocated for a chain of parameters. */
void
_gmx_selexpr_free_params(t_selexpr_param *param);

/** Propagates the flags for selection elements. */
int
_gmx_selelem_update_flags(struct t_selelem *sel);

/** Initializes the method parameter data of \ref SEL_EXPRESSION and
 * \ref SEL_MODIFIER elements. */
void
_gmx_selelem_init_method_params(struct t_selelem *sel, void *scanner);
/** Initializes the method for a \ref SEL_EXPRESSION selection element. */
void
_gmx_selelem_set_method(struct t_selelem *sel,
                        struct gmx_ana_selmethod_t *method, void *scanner);

/** Creates a \c t_selelem for arithmetic expression evaluation. */
struct t_selelem *
_gmx_sel_init_arithmetic(struct t_selelem *left, struct t_selelem *right,
                         char op, void *scanner);
/** Creates a \c t_selelem for comparsion expression evaluation. */
struct t_selelem *
_gmx_sel_init_comparison(struct t_selelem *left, struct t_selelem *right,
                         char *cmpop, void *scanner);
/** Creates a \c t_selelem for a keyword expression from the parsed data. */
struct t_selelem *
_gmx_sel_init_keyword(struct gmx_ana_selmethod_t *method,
                      t_selexpr_value *args, const char *rpost, void *scanner);
/** Creates a \c t_selelem for a method expression from the parsed data. */
struct t_selelem *
_gmx_sel_init_method(struct gmx_ana_selmethod_t *method,
                     t_selexpr_param *params, const char *rpost,
                     void *scanner);
/** Creates a \c t_selelem for a modifier expression from the parsed data. */
struct t_selelem *
_gmx_sel_init_modifier(struct gmx_ana_selmethod_t *mod, t_selexpr_param *params,
                       struct t_selelem *sel, void *scanner);
/** Creates a \c t_selelem for evaluation of reference positions. */
struct t_selelem *
_gmx_sel_init_position(struct t_selelem *expr, const char *type, void *scanner);

/** Creates a \c t_selelem for a constant position. */
struct t_selelem *
_gmx_sel_init_const_position(real x, real y, real z);
/** Creates a \c t_selelem for a index group expression using group name. */
struct t_selelem *
_gmx_sel_init_group_by_name(const char *name, void *scanner);
/** Creates a \c t_selelem for a index group expression using group index. */
struct t_selelem *
_gmx_sel_init_group_by_id(int id, void *scanner);
/** Creates a \c t_selelem for a variable reference */
struct t_selelem *
_gmx_sel_init_variable_ref(struct t_selelem *sel);

/** Creates a root \c t_selelem for a selection. */
struct t_selelem *
_gmx_sel_init_selection(char *name, struct t_selelem *sel, void *scanner);
/** Creates a root \c t_selelem elements for a variable assignment. */
struct t_selelem *
_gmx_sel_assign_variable(char *name, struct t_selelem *expr, void *scanner);
/** Appends a root \c t_selelem to a selection collection. */
struct t_selelem *
_gmx_sel_append_selection(struct t_selelem *sel, struct t_selelem *last,
                          void *scanner);
/** Check whether the parser should finish. */
gmx_bool
_gmx_sel_parser_should_finish(void *scanner);

/** Handle empty commands. */
void
_gmx_sel_handle_empty_cmd(void *scanner);
/** Process help commands. */
void
_gmx_sel_handle_help_cmd(char *topic, void *scanner);

/* In params.c */
/** Initializes an array of parameters based on input from the selection parser. */
gmx_bool
_gmx_sel_parse_params(t_selexpr_param *pparams, int nparam,
                      struct gmx_ana_selparam_t *param, struct t_selelem *root,
                      void *scanner);

#endif
