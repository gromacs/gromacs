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
 * \brief
 * Handling of intermediate selection parser data.
 *
 * The data types declared in this header are used by the parser to store
 * intermediate data when constructing method expressions.
 * In particular, the parameters for the method are stored.
 * The intermediate data is freed once a gmx::SelectionTreeElement object can
 * be constructed.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_PARSETREE_H
#define GMX_SELECTION_PARSETREE_H

#include <exception>

#include "gromacs/legacyheaders/types/simple.h"

#include "selelem.h"
#include "selvalue.h"

struct gmx_ana_indexgrps_t;
struct gmx_ana_selmethod_t;
struct gmx_ana_selparam_t;

/*! \internal \brief
 * Describes a parsed value, possibly resulting from expression evaluation.
 *
 * \todo
 * Make this a proper class.
 */
typedef struct t_selexpr_value
{
    //! Returns true if the value comes from expression evaluation.
    bool hasExpressionValue() const { return expr; }

    //! Type of the value.
    e_selvalue_t            type;
    //! Expression pointer if the value is the result of an expression.
    gmx::SelectionTreeElementPointer expr;
    //! The actual value if \p expr is NULL.
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
_gmx_selparser_error(void *scanner, const char *fmt, ...);
/** Handle exceptions caught within the Bison code. */
bool
_gmx_selparser_handle_exception(void *scanner, const std::exception &ex);

/** Allocates and initializes a constant \c t_selexpr_value. */
t_selexpr_value *
_gmx_selexpr_create_value(e_selvalue_t type);
/** Allocates and initializes an expression \c t_selexpr_value. */
t_selexpr_value *
_gmx_selexpr_create_value_expr(const gmx::SelectionTreeElementPointer &expr);
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
void
_gmx_selelem_update_flags(const gmx::SelectionTreeElementPointer &sel,
                          void *scanner);

/** Initializes the method parameter data of \ref SEL_EXPRESSION and
 * \ref SEL_MODIFIER elements. */
void
_gmx_selelem_init_method_params(const gmx::SelectionTreeElementPointer &sel,
                                void *scanner);
/** Initializes the method for a \ref SEL_EXPRESSION selection element. */
void
_gmx_selelem_set_method(const gmx::SelectionTreeElementPointer &sel,
                        struct gmx_ana_selmethod_t *method, void *scanner);

/** Creates a gmx::SelectionTreeElement for arithmetic expression evaluation. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_arithmetic(const gmx::SelectionTreeElementPointer &left,
                         const gmx::SelectionTreeElementPointer &right,
                         char op, void *scanner);
/** Creates a gmx::SelectionTreeElement for comparsion expression evaluation. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_comparison(const gmx::SelectionTreeElementPointer &left,
                         const gmx::SelectionTreeElementPointer &right,
                         char *cmpop, void *scanner);
/** Creates a gmx::SelectionTreeElement for a keyword expression from the parsed data. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword(struct gmx_ana_selmethod_t *method,
                      t_selexpr_value *args, const char *rpost, void *scanner);
/** Creates a gmx::SelectionTreeElement for a method expression from the parsed data. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_method(struct gmx_ana_selmethod_t *method,
                     t_selexpr_param *params, const char *rpost,
                     void *scanner);
/** Creates a gmx::SelectionTreeElement for a modifier expression from the parsed data. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_modifier(struct gmx_ana_selmethod_t *mod, t_selexpr_param *params,
                       const gmx::SelectionTreeElementPointer &sel,
                       void *scanner);
/** Creates a gmx::SelectionTreeElement for evaluation of reference positions. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_position(const gmx::SelectionTreeElementPointer &expr,
                       const char *type, void *scanner);

/** Creates a gmx::SelectionTreeElement for a constant position. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_const_position(real x, real y, real z);
/** Creates a gmx::SelectionTreeElement for a index group expression using group name. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_group_by_name(const char *name, void *scanner);
/** Creates a gmx::SelectionTreeElement for a index group expression using group index. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_group_by_id(int id, void *scanner);
/** Creates a gmx::SelectionTreeElement for a variable reference */
gmx::SelectionTreeElementPointer
_gmx_sel_init_variable_ref(const gmx::SelectionTreeElementPointer &sel);

/** Creates a root gmx::SelectionTreeElement for a selection. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_selection(const char *name,
                        const gmx::SelectionTreeElementPointer &sel,
                        void *scanner);
/** Creates a root gmx::SelectionTreeElement elements for a variable assignment. */
gmx::SelectionTreeElementPointer
_gmx_sel_assign_variable(const char *name,
                         const gmx::SelectionTreeElementPointer &expr,
                         void *scanner);
/** Appends a root gmx::SelectionTreeElement to a selection collection. */
gmx::SelectionTreeElementPointer
_gmx_sel_append_selection(const gmx::SelectionTreeElementPointer &sel,
                          gmx::SelectionTreeElementPointer last,
                          void *scanner);
/** Check whether the parser should finish. */
bool
_gmx_sel_parser_should_finish(void *scanner);

/** Handle empty commands. */
void
_gmx_sel_handle_empty_cmd(void *scanner);
/** Process help commands. */
void
_gmx_sel_handle_help_cmd(t_selexpr_value *topic, void *scanner);

/* In params.c */
/** Initializes an array of parameters based on input from the selection parser. */
bool
_gmx_sel_parse_params(t_selexpr_param *pparams, int nparam,
                      struct gmx_ana_selparam_t *param,
                      const gmx::SelectionTreeElementPointer &root,
                      void *scanner);

#endif
