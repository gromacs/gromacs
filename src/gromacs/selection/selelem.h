/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares gmx::SelectionTreeElement and related things.
 *
 * The selection element trees constructed by the parser and the compiler
 * are described on the respective pages:
 * \ref page_module_selection_parser and \ref page_module_selection_compiler.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELELEM_H
#define GMX_SELECTION_SELELEM_H

#include <string>

#include <boost/shared_ptr.hpp>

#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

#include "selvalue.h"

struct gmx_ana_poscalc_t;
struct gmx_ana_selparam_t;
struct gmx_ana_selmethod_t;

struct gmx_sel_evaluate_t;
struct gmx_sel_mempool_t;

struct t_compiler_data;

namespace gmx
{
class SelectionTreeElement;

//! Smart pointer type for selection tree element pointers.
typedef boost::shared_ptr<SelectionTreeElement> SelectionTreeElementPointer;
} // namespace gmx

/********************************************************************/
/*! \name Enumerations for expression types
 ********************************************************************/
//!\{

/** Defines the type of a gmx::SelectionTreeElement object. */
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
    /** Unresolved reference to an external group. */
    SEL_GROUPREF,
    /** Post-processing of selection value. */
    SEL_MODIFIER
} e_selelem_t;

/** Defines the boolean operation of gmx::SelectionTreeElement objects with type \ref SEL_BOOLEAN. */
typedef enum
{
    BOOL_NOT,           /**< Not */
    BOOL_AND,           /**< And */
    BOOL_OR,            /**< Or */
    BOOL_XOR            /**< Xor (not implemented). */
} e_boolean_t;

/** Defines the arithmetic operation of gmx::SelectionTreeElement objects with type \ref SEL_ARITHMETIC. */
typedef enum
{
    ARITH_PLUS,         /**< Addition (`+`) */
    ARITH_MINUS,        /**< Subtraction (`-`) */
    ARITH_NEG,          /**< Unary `-` */
    ARITH_MULT,         /**< Multiplication (`*`) */
    ARITH_DIV,          /**< Division (`/`) */
    ARITH_EXP           /**< Power (`^`) */
} e_arithmetic_t;

/** Returns a string representation of the type of a gmx::SelectionTreeElement. */
extern const char *
_gmx_selelem_type_str(const gmx::SelectionTreeElement &sel);
/** Returns a string representation of the boolean type of a \ref SEL_BOOLEAN gmx::SelectionTreeElement. */
extern const char *
_gmx_selelem_boolean_type_str(const gmx::SelectionTreeElement &sel);
/** Returns a string representation of the type of a \c gmx_ana_selvalue_t. */
extern const char *
_gmx_sel_value_type_str(const gmx_ana_selvalue_t *val);

//!\}


/********************************************************************/
/*! \name Selection expression flags
 * \anchor selelem_flags
 ********************************************************************/
//!\{
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
 * The element may contain atom indices in an unsorted order.
 */
#define SEL_UNSORTED    32
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
//!\}


namespace gmx
{

class ExceptionInitializer;

//! \cond internal
/*! \brief
 * Function pointer for evaluating a gmx::SelectionTreeElement.
 */
typedef void (*sel_evalfunc)(struct gmx_sel_evaluate_t         *data,
                             const SelectionTreeElementPointer &sel,
                             gmx_ana_index_t                   *g);
//! \endcond

/*! \internal
 * \brief
 * Stores the location of a selection element in the selection text.
 *
 * The location is stored as a range in the pretty-printed selection text
 * (where whitespace has been sanitized), and can be used to extract that text
 * for error messages and other diagnostic purposes.
 * During parsing, the extraction is done with _gmx_sel_lexer_get_text().
 *
 * This needs to be a plain C struct for Bison to properly deal with it.
 */
struct SelectionLocation
{
    //! Returns an empty location.
    static SelectionLocation createEmpty()
    {
        SelectionLocation empty = {0, 0};
        return empty;
    }

    //! Start index of the string where this element has been parsed from.
    int  startIndex;
    //! End index of the string where this element has been parsed from.
    int  endIndex;
};

/*! \internal \brief
 * Represents an element of a selection expression.
 */
class SelectionTreeElement
{
    public:
        /*! \brief
         * Allocates memory and performs common initialization.
         *
         * \param[in] type     Type of selection element to create.
         * \param[in] location Location of the element.
         *
         * \a type is set to \p type,
         * \a v::type is set to \ref GROUP_VALUE for boolean and comparison
         * expressions and \ref NO_VALUE for others, and
         * \ref SEL_ALLOCVAL is set for non-root elements (\ref SEL_ALLOCDATA
         * is also set for \ref SEL_BOOLEAN elements).
         * All the pointers are set to NULL.
         */
        SelectionTreeElement(e_selelem_t type, const SelectionLocation &location);
        ~SelectionTreeElement();

        //! Frees the memory allocated for the \a v union.
        void freeValues();
        //! Frees the memory allocated for the \a u union.
        void freeExpressionData();
        /* In compiler.cpp */
        /*! \brief
         * Frees the memory allocated for the selection compiler.
         *
         * This function only frees the data for the given selection, not its
         * children.  It is safe to call the function when compiler data has
         * not been allocated or has already been freed; in such a case,
         * nothing is done.
         */
        void freeCompilerData();

        /*! \brief
         * Reserves memory for value from a memory pool.
         *
         * \param[in]     count Number of values to reserve memory for.
         *
         * Reserves memory for the values of this element from the \a mempool
         * memory pool.
         * If no memory pool is set, nothing is done.
         */
        void mempoolReserve(int count);
        /*! \brief
         * Releases memory pool used for value.
         *
         * Releases the memory allocated for the values of this element from the
         * \a mempool memory pool.
         * If no memory pool is set, nothing is done.
         */
        void mempoolRelease();

        //! Returns the name of the element.
        const std::string &name() const { return name_; }
        //! Returns the location of the element.
        const SelectionLocation &location() const { return location_; }

        /*! \brief
         * Sets the name of the element.
         *
         * \param[in] name  Name to set (can be NULL).
         * \throws    std::bad_alloc if out of memory.
         */
        void setName(const char *name) { name_ = (name != NULL ? name : ""); }
        //! \copydoc setName(const char *)
        void setName(const std::string &name) { name_ = name; }
        /*! \brief
         * Sets the name of a root element if it is missing.
         *
         * \param[in] selectionText  Full selection text to use as a fallback.
         * \throws    std::bad_alloc if out of memory.
         *
         * If index groups have not yet been set and the selection is a result
         * of a group reference, the name may still be empty after this call.
         *
         * Strong exception safety guarantee.
         */
        void fillNameIfMissing(const char *selectionText);

        /*! \brief
         * Checks that this element and its children do not contain unsupported
         * elements with unsorted atoms.
         *
         * \param[in] bUnsortedAllowed Whether this element's parents allow it
         *     to have unsorted atoms.
         * \param     errors           Object for reporting any error messages.
         * \throws    std::bad_alloc if out of memory.
         *
         * Errors are reported as nested exceptions in \p errors.
         */
        void checkUnsortedAtoms(bool                  bUnsortedAllowed,
                                ExceptionInitializer *errors) const;
        /*! \brief
         * Resolves an unresolved reference to an index group.
         *
         * \param[in] grps   Index groups to use to resolve the reference.
         * \param[in] natoms Maximum number of atoms the selections can evaluate to
         *     (zero if the topology/atom count is not set yet).
         * \throws    std::bad_alloc if out of memory.
         * \throws    InconsistentInputError if the reference cannot be
         *     resolved.
         */
        void resolveIndexGroupReference(gmx_ana_indexgrps_t *grps, int natoms);
        /*! \brief
         * Checks that an index group has valid atom indices.
         *
         * \param[in] natoms Maximum number of atoms the selections can evaluate to.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InconsistentInputError if there are invalid atom indices.
         */
        void checkIndexGroup(int natoms);

        //! Type of the element.
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
         * Can be either NULL (if the expression is a constant and does not
         * require evaluation) or point to one of the functions defined in
         * evaluate.h.
         */
        sel_evalfunc                        evaluate;
        /*! \brief
         * Information flags about the element.
         *
         * Allowed flags are listed here:
         * \ref selelem_flags "flags for gmx::SelectionTreeElement".
         */
        int                                 flags;
        //! Data required by the evaluation function.
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
            //! Data for \ref SEL_EXPRESSION and \ref SEL_MODIFIER elements.
            struct {
                //! Pointer the the method used in this expression.
                struct gmx_ana_selmethod_t *method;
                //! Pointer to the data allocated by the method's \p init_data (see sel_datafunc()).
                void                       *mdata;
                //! Pointer to the position data passed to the method.
                struct gmx_ana_pos_t       *pos;
                //! Pointer to the evaluation data for \p pos.
                struct gmx_ana_poscalc_t   *pc;
            }                               expr;
            //! Operation type for \ref SEL_BOOLEAN elements.
            e_boolean_t                     boolt;
            //! Operation type for \ref SEL_ARITHMETIC elements.
            struct {
                //! Operation type.
                e_arithmetic_t              type;
                //! String representation.
                char                       *opstr;
            }                               arith;
            //! Associated selection parameter for \ref SEL_SUBEXPRREF elements.
            struct gmx_ana_selparam_t      *param;
            //! The string/number used to reference the group.
            struct {
                //! Name of the referenced external group.
                char                       *name;
                //! If \a name is NULL, the index number of the referenced group.
                int                         id;
            }                               gref;
        }                                   u;
        //! Memory pool to use for values, or NULL if standard memory handling.
        struct gmx_sel_mempool_t           *mempool;
        //! Internal data for the selection compiler.
        t_compiler_data                    *cdata;

        /*! \brief The first child element.
         *
         * Other children can be accessed through the \p next field of \p child.
         */
        SelectionTreeElementPointer         child;
        //! The next sibling element.
        SelectionTreeElementPointer         next;

    private:
        /*! \brief
         * Name of the element.
         *
         * This field is only used for diagnostic purposes.
         */
        std::string                         name_;
        /*! \brief
         * Location of the element in the selection text.
         *
         * This field is only used for diagnostic purposes (including error
         * messages).
         */
        SelectionLocation                   location_;

        GMX_DISALLOW_COPY_AND_ASSIGN(SelectionTreeElement);
};

} // namespace gmx

/********************************************************************/
/*! \name Selection expression functions
 */
//!\{

/* In evaluate.c */
/** Writes out a human-readable name for an evaluation function. */
void
_gmx_sel_print_evalfunc_name(FILE *fp, gmx::sel_evalfunc evalfunc);

/** Sets the value type of a gmx::SelectionTreeElement. */
void
_gmx_selelem_set_vtype(const gmx::SelectionTreeElementPointer &sel,
                       e_selvalue_t                            vtype);

/*! \brief
 * Frees the memory allocated for a selection method parameter.
 *
 * \param[in] param Parameter to free.
 */
void
_gmx_selelem_free_param(struct gmx_ana_selparam_t *param);
/*! \brief
 * Frees the memory allocated for a selection method.
 *
 * \param[in] method Method to free.
 * \param[in] mdata  Method data to free.
 */
void
_gmx_selelem_free_method(struct gmx_ana_selmethod_t *method, void *mdata);

/** Prints a human-readable version of a selection element subtree. */
void
_gmx_selelem_print_tree(FILE *fp, const gmx::SelectionTreeElement &sel,
                        bool bValues, int level);
/* In compiler.c */
/** Prints a human-readable version of the internal compiler data structure. */
void
_gmx_selelem_print_compiler_info(FILE *fp, const gmx::SelectionTreeElement &sel,
                                 int level);

/** Returns true if the selection element subtree requires topology information for evaluation. */
bool
_gmx_selelem_requires_top(const gmx::SelectionTreeElement &root);

/* In sm_insolidangle.c */
/** Returns true if the covered fraction of the selection can be calculated. */
bool
_gmx_selelem_can_estimate_cover(const gmx::SelectionTreeElement &sel);
/** Returns the covered fraction of the selection for the current frame. */
real
_gmx_selelem_estimate_coverfrac(const gmx::SelectionTreeElement &sel);

//!\}

#endif
