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
 * \brief Definition of \c t_selelem and related things.
 *
 * The selection element trees constructed by the parser and the compiler
 * are described on the respective pages:
 * \ref selparser and \ref selcompiler.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef SELECTION_ELEMENT_H
#define SELECTION_ELEMENT_H

#include <types/simple.h>

#include <indexutil.h>
#include <selvalue.h>

struct gmx_ana_poscalc_t;
struct gmx_ana_selparam_t;
struct gmx_ana_selmethod_t;

struct gmx_sel_evaluate_t;
struct gmx_sel_mempool_t;
struct t_selelem;

/********************************************************************/
/*! \name Enumerations for expression types
 ********************************************************************/
/*@{*/

/** Defines the type of a \c t_selelem object. */
typedef enum
{
    /** Constant-valued expression. */
    SEL_CONST,
    /** Method expression that requires evaluation. */
    SEL_EXPRESSION,
    /** Boolean expression. */
    SEL_BOOLEAN,
    /** Arithmetic expression. */
    SEL_ARITHMETIC,
    /** Root node of the evaluation tree. */
    SEL_ROOT,
    /** Subexpression that may be referenced several times. */
    SEL_SUBEXPR,
    /** Reference to a subexpression. */
    SEL_SUBEXPRREF,
    /** Post-processing of selection value. */
    SEL_MODIFIER
} e_selelem_t;

/** Defines the gmx_boolean operation of \c t_selelem objects with type \ref SEL_BOOLEAN. */
typedef enum
{
    BOOL_NOT,           /**< Not */
    BOOL_AND,           /**< And */
    BOOL_OR,            /**< Or */
    BOOL_XOR            /**< Xor (not implemented). */
} e_boolean_t;

/** Defines the arithmetic operation of \c t_selelem objects with type \ref SEL_ARITHMETIC. */
typedef enum
{
    ARITH_PLUS,         /**< + */
    ARITH_MINUS,        /**< - */
    ARITH_NEG,          /**< Unary - */
    ARITH_MULT,         /**< * */
    ARITH_DIV,          /**< / */
    ARITH_EXP           /**< ^ (to power) */
} e_arithmetic_t;

/** Returns a string representation of the type of a \c t_selelem. */
extern const char *
_gmx_selelem_type_str(struct t_selelem *sel);
/** Returns a string representation of the gmx_boolean type of a \ref SEL_BOOLEAN \c t_selelem. */
extern const char *
_gmx_selelem_gmx_boolean_type_str(struct t_selelem *sel);
/** Returns a string representation of the type of a \c gmx_ana_selvalue_t. */
extern const char *
_gmx_sel_value_type_str(gmx_ana_selvalue_t *val);

/*@}*/


/********************************************************************/
/*! \name Selection expression flags
 * \anchor selelem_flags
 ********************************************************************/
/*@{*/
/*! \brief
 * Selection value flags are set.
 *
 * If this flag is set, the flags covered by \ref SEL_VALFLAGMASK
 * have been set properly for the element.
 */
#define SEL_FLAGSSET    1
/*! \brief
 * The element evaluates to a single value.
 *
 * This flag is always set for \ref GROUP_VALUE elements.
 */
#define SEL_SINGLEVAL   2
/*! \brief
 * The element evaluates to one value for each input atom.
 */
#define SEL_ATOMVAL     4
/*! \brief
 * The element evaluates to an arbitrary number of values.
 */
#define SEL_VARNUMVAL   8
/*! \brief
 * The element (or one of its children) is dynamic.
 */
#define SEL_DYNAMIC     16
/*! \brief
 * Mask that covers the flags that describe the number of values.
 */
#define SEL_VALTYPEMASK (SEL_SINGLEVAL | SEL_ATOMVAL | SEL_VARNUMVAL)
/*! \brief
 * Mask that covers the flags that describe the value type.
 */
#define SEL_VALFLAGMASK (SEL_FLAGSSET | SEL_VALTYPEMASK | SEL_DYNAMIC)
/*! \brief
 * Data has been allocated for the \p v.u union.
 *
 * If not set, the \p v.u.ptr points to data allocated externally.
 * This is the case if the value of the element is used as a parameter
 * for a selection method or if the element evaluates the final value of
 * a selection.
 *
 * Even if the flag is set, \p v.u.ptr can be NULL during initialization.
 *
 * \todo
 * This flag overlaps with the function of \p v.nalloc field, and could
 * probably be removed, making memory management simpler. Currently, the
 * \p v.nalloc field is not kept up-to-date in all cases when this flag
 * is changed and is used in places where this flag is not, so this would
 * require a careful investigation of the selection code.
 */
#define SEL_ALLOCVAL    (1<<8)
/*! \brief
 * Data has been allocated for the group/position structure.
 *
 * If not set, the memory allocated for fields in \p v.u.g or \p v.u.p is
 * managed externally.
 *
 * This field has no effect if the value type is not \ref GROUP_VALUE or
 * \ref POS_VALUE, but should not be set.
 */
#define SEL_ALLOCDATA   (1<<9)
/*! \brief
 * \p method->init_frame should be called for the frame.
 */
#define SEL_INITFRAME   (1<<10)
/*! \brief
 * Parameter has been evaluated for the current frame.
 *
 * This flag is set for children of \ref SEL_EXPRESSION elements (which
 * describe method parameters) after the element has been evaluated for the
 * current frame.
 * It is not set for \ref SEL_ATOMVAL elements, because they may need to
 * be evaluated multiple times.
 */
#define SEL_EVALFRAME   (1<<11)
/*! \brief
 * \p method->init has been called.
 */
#define SEL_METHODINIT  (1<<12)
/*! \brief
 * \p method->outinit has been called.
 *
 * This flag is also used for \ref SEL_SUBEXPRREF elements.
 */
#define SEL_OUTINIT     (1<<13)
/*@}*/


/********************************************************************/
/*! \name Selection expression data structures and functions
 ********************************************************************/
/*@{*/

struct t_selelem;

/*! \brief
 * Function pointer for evaluating a \c t_selelem.
 */
typedef int (*sel_evalfunc)(struct gmx_sel_evaluate_t *data,
                            struct t_selelem *sel, gmx_ana_index_t *g);

/*! \internal \brief
 * Represents an element of a selection expression.
 */
typedef struct t_selelem
{
    /*! \brief Name of the element.
     *
     * This field is only used for informative purposes.
     * It is always either NULL or a pointer to a string.
     * Memory is never allocated for it directly.
     */
    const char                         *name;
    /** Type of the element. */
    e_selelem_t                         type;
    /*! \brief
     * Value storage of the element.
     *
     * This field contains the evaluated value of the element, as well as
     * the output value type.
     */
    gmx_ana_selvalue_t                  v;
    /*! \brief
     * Evaluation function for the element.
     *
     * Can be either NULL (if the expression is a constant and does not require
     * evaluation) or point to one of the functions defined in evaluate.h.
     */
    sel_evalfunc                        evaluate;
    /*! \brief
     * Information flags about the element.
     *
     * Allowed flags are listed here:
     * \ref selelem_flags "flags for \c t_selelem".
     */
    int                                 flags;
    /** Data required by the evaluation function. */
    union {
        /*! \brief Index group data for several element types.
         *
         *  - \ref SEL_CONST : if the value type is \ref GROUP_VALUE,
         *    this field holds the unprocessed group value.
         *  - \ref SEL_ROOT : holds the group value for which the
         *    selection subtree should be evaluated.
         *  - \ref SEL_SUBEXPR : holds the group for which the subexpression
         *    has been evaluated.
         */
        gmx_ana_index_t                 cgrp;
        /** Data for \ref SEL_EXPRESSION and \ref SEL_MODIFIER elements. */
        struct {
            /** Pointer the the method used in this expression. */
            struct gmx_ana_selmethod_t *method;
            /** Pointer to the data allocated by the method's \p init_data (see sel_datafunc()). */
            void                       *mdata;
            /** Pointer to the position data passed to the method. */
            struct gmx_ana_pos_t       *pos;
            /** Pointer to the evaluation data for \p pos. */
            struct gmx_ana_poscalc_t   *pc;
        }                               expr;
        /** Operation type for \ref SEL_BOOLEAN elements. */
        e_boolean_t                     boolt;
        /** Operation type for \ref SEL_ARITHMETIC elements. */
        struct {
            /** Operation type. */
            e_arithmetic_t              type;
            /** String representation. */
            char                       *opstr;
        }                               arith;
        /** Associated selection parameter for \ref SEL_SUBEXPRREF elements. */
        struct gmx_ana_selparam_t      *param;
    }                                   u;
    /** Memory pool to use for values, or NULL if standard memory handling. */
    struct gmx_sel_mempool_t           *mempool;
    /** Internal data for the selection compiler. */
    struct t_compiler_data             *cdata;

    /*! \brief The first child element.
     *
     * Other children can be accessed through the \p next field of \p child.
     */
    struct t_selelem                    *child;
    /** The next sibling element. */
    struct t_selelem                    *next;
    /*! \brief Number of references to this element.
     *
     * Should be larger than one only for \ref SEL_SUBEXPR elements.
     */
    int                                  refcount;
} t_selelem;

/* In evaluate.c */
/** Writes out a human-readable name for an evaluation function. */
extern void
_gmx_sel_print_evalfunc_name(FILE *fp, sel_evalfunc evalfunc);

/** Allocates memory and performs some common initialization for a \c t_selelem. */
extern t_selelem *
_gmx_selelem_create(e_selelem_t type);
/** Sets the value type of a \c t_selelem. */
extern int
_gmx_selelem_set_vtype(t_selelem *sel, e_selvalue_t vtype);
/** Reserves memory for value of a \c t_selelem from a memory pool. */
extern int
_gmx_selelem_mempool_reserve(t_selelem *sel, int count);
/** Releases memory pool used for value of a \c t_selelem. */
extern void
_gmx_selelem_mempool_release(t_selelem *sel);
/** Frees the memory allocated for a \c t_selelem structure and all its children. */
extern void
_gmx_selelem_free(t_selelem *sel);
/** Frees the memory allocated for a \c t_selelem structure, all its children, and also all structures referenced through t_selelem::next fields. */
extern void
_gmx_selelem_free_chain(t_selelem *first);

/** Frees the memory allocated for the \c t_selelem::d union. */
extern void
_gmx_selelem_free_values(t_selelem *sel);
/** Frees the memory allocated for a selection method. */
extern void
_gmx_selelem_free_method(struct gmx_ana_selmethod_t *method, void *mdata);
/** Frees the memory allocated for the \c t_selelem::u field. */
extern void
_gmx_selelem_free_exprdata(t_selelem *sel);
/* In compiler.c */
/** Frees the memory allocated for the selection compiler. */
extern void
_gmx_selelem_free_compiler_data(t_selelem *sel);

/** Prints a human-readable version of a selection element subtree. */
extern void
_gmx_selelem_print_tree(FILE *fp, t_selelem *root, gmx_bool bValues, int level);
/* In compile.c */
/** Prints a human-readable version of the internal compiler data structure. */
extern void
_gmx_selelem_print_compiler_info(FILE *fp, t_selelem *sel, int level);

/** Returns TRUE if the selection element subtree requires topology information for evaluation. */
extern gmx_bool
_gmx_selelem_requires_top(t_selelem *root);

/* In sm_insolidangle.c */
/** Returns TRUE if the covered fraction of the selection can be calculated. */
extern gmx_bool
_gmx_selelem_can_estimate_cover(t_selelem *sel);
/** Returns the covered fraction of the selection for the current frame. */
extern real
_gmx_selelem_estimate_coverfrac(t_selelem *sel);

/*@}*/

#endif
