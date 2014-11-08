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
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_PARSETREE_H
#define GMX_SELECTION_PARSETREE_H

#include <exception>
#include <list>
#include <string>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/uniqueptr.h"

#include "selelem.h"
#include "selvalue.h"

struct gmx_ana_indexgrps_t;
struct gmx_ana_selmethod_t;
struct gmx_ana_selparam_t;

namespace gmx
{

//! \cond internal
/*! \internal \brief
 * String matching mode for string keyword expressions.
 *
 * \ingroup module_selection
 */
enum SelectionStringMatchType
{
    eStringMatchType_Auto,              //!< Deduce from the string.
    eStringMatchType_Exact,             //!< Match as a literal string.
    eStringMatchType_Wildcard,          //!< Match using ? and * as wildcards.
    eStringMatchType_RegularExpression  //!< Match using regular expressions.
};
/*! \endcond */

class SelectionParserValue;

//! Container for a list of SelectionParserValue objects.
typedef std::list<SelectionParserValue>
    SelectionParserValueList;
//! Smart pointer type for managing a SelectionParserValueList.
typedef gmx::gmx_unique_ptr<SelectionParserValueList>::type
    SelectionParserValueListPointer;

/*! \internal
 * \brief
 * Describes a parsed value, possibly resulting from expression evaluation.
 *
 * All factory methods and the constructors may throw an std::bad_alloc if
 * out of memory.
 *
 * \ingroup module_selection
 */
class SelectionParserValue
{
    public:
        //! Allocates and initializes an empty value list.
        static SelectionParserValueListPointer createList()
        {
            return SelectionParserValueListPointer(new SelectionParserValueList);
        }
        /*! \brief
         * Allocates and initializes a value list with a single value.
         *
         * \param[in] value  Initial value to put in the list.
         * \returns   Pointer to a new value list that contains \p value.
         */
        static SelectionParserValueListPointer
        createList(const SelectionParserValue &value)
        {
            SelectionParserValueListPointer list(new SelectionParserValueList);
            list->push_back(value);
            return move(list);
        }
        /*! \brief
         * Allocates and initializes an expression value.
         *
         * \param[in] expr  Root of the expression tree to assign to the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createExpr(const gmx::SelectionTreeElementPointer &expr)
        {
            return SelectionParserValue(expr);
        }
        /*! \brief
         * Allocates and initializes a constant integer value.
         *
         * \param[in] value    Integer value to assign to the value.
         * \param[in] location Location of the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createInteger(int value, const SelectionLocation &location)
        {
            SelectionParserValue result(INT_VALUE, location);
            result.u.i.i1 = result.u.i.i2 = value;
            return result;
        }
        /*! \brief
         * Allocates and initializes a constant integer range value.
         *
         * \param[in] from     Beginning of the range to assign to the value.
         * \param[in] to       End of the range to assign to the value.
         * \param[in] location Location of the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createIntegerRange(int from, int to, const SelectionLocation &location)
        {
            SelectionParserValue result(INT_VALUE, location);
            result.u.i.i1 = from;
            result.u.i.i2 = to;
            return result;
        }
        /*! \brief
         * Allocates and initializes a constant floating-point value.
         *
         * \param[in] value    Floating-point value to assign to the value.
         * \param[in] location Location of the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createReal(real value, const SelectionLocation &location)
        {
            SelectionParserValue result(REAL_VALUE, location);
            result.u.r.r1 = result.u.r.r2 = value;
            return result;
        }
        /*! \brief
         * Allocates and initializes a constant floating-point range value.
         *
         * \param[in] from     Beginning of the range to assign to the value.
         * \param[in] to       End of the range to assign to the value.
         * \param[in] location Location of the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createRealRange(real from, real to, const SelectionLocation &location)
        {
            SelectionParserValue result(REAL_VALUE, location);
            result.u.r.r1 = from;
            result.u.r.r2 = to;
            return result;
        }
        /*! \brief
         * Allocates and initializes a constant string value.
         *
         * \param[in] value    String to assign to the value.
         * \param[in] location Location of the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createString(const char *value, const SelectionLocation &location)
        {
            SelectionParserValue result(STR_VALUE, location);
            result.str = value;
            return result;
        }
        /*! \brief
         * Allocates and initializes a constant position value.
         *
         * \param[in] value    Position vector to assign to the value.
         * \param[in] location Location of the value.
         * \returns   The newly created value.
         */
        static SelectionParserValue
        createPosition(rvec value, const SelectionLocation &location)
        {
            SelectionParserValue result(POS_VALUE, location);
            copy_rvec(value, result.u.x);
            return result;
        }

        //! Returns the location of this value in the parsed selection text.
        const SelectionLocation &location() const { return location_; }
        //! Returns true if the value comes from expression evaluation.
        bool hasExpressionValue() const { return static_cast<bool>(expr); }

        //! Returns the string value (\a type must be ::STR_VALUE).
        const std::string &stringValue() const
        {
            GMX_ASSERT(type == STR_VALUE && !hasExpressionValue(),
                       "Attempted to retrieve string value from a non-string value");
            return str;
        }

        // TODO: boost::any or similar could be nicer for the implementation.
        //! Type of the value.
        e_selvalue_t                     type;
        //! Expression pointer if the value is the result of an expression.
        gmx::SelectionTreeElementPointer expr;
        //! String value for \a type ::STR_VALUE.
        std::string                      str;
        //! The actual value if \a expr is NULL and \a type is not ::STR_VALUE.
        union {
            //! The integer value/range (\a type ::INT_VALUE).
            struct {
                //! Beginning of the range.
                int             i1;
                //! End of the range; equals \a i1 for a single integer.
                int             i2;
            }                   i;
            //! The real value/range (\a type ::REAL_VALUE).
            struct {
                //! Beginning of the range.
                real            r1;
                //! End of the range; equals \a r1 for a single number.
                real            r2;
            }                   r;
            //! The position value (\a type ::POS_VALUE).
            rvec                x;
        }                       u;

    private:
        /*! \brief
         * Initializes a new value.
         *
         * \param[in] type     Type for the new value.
         * \param[in] location Location for the value.
         */
        SelectionParserValue(e_selvalue_t type, const SelectionLocation &location);
        /*! \brief
         * Initializes a new expression value.
         *
         * \param[in] expr  Expression for the value.
         */
        explicit SelectionParserValue(const gmx::SelectionTreeElementPointer &expr);

        //! Location of the value in the parsed text.
        SelectionLocation       location_;
};

class SelectionParserParameter;

//! Container for a list of SelectionParserParameter objects.
typedef std::list<SelectionParserParameter>
    SelectionParserParameterList;
//! Smart pointer type for managing a SelectionParserParameterList.
typedef gmx::gmx_unique_ptr<SelectionParserParameterList>::type
    SelectionParserParameterListPointer;

/*! \internal \brief
 * Describes a parsed method parameter.
 *
 * \ingroup module_selection
 */
class SelectionParserParameter
{
    public:
        //! Allocates and initializes an empty parameter list.
        static SelectionParserParameterListPointer createList()
        {
            return SelectionParserParameterListPointer(
                    new SelectionParserParameterList);
        }
        /*! \brief
         * Allocates and initializes a parsed method parameter.
         *
         * \param[in] name     Name for the new parameter (can be NULL).
         * \param[in] values   List of values for the parameter.
         * \param[in] location Location of the parameter.
         * \returns   Pointer to the newly allocated parameter.
         * \throws    std::bad_alloc if out of memory.
         */
        static SelectionParserParameter
        create(const char *name, SelectionParserValueListPointer values,
               const SelectionLocation &location)
        {
            return SelectionParserParameter(name, move(values), location);
        }
        //! \copydoc create(const char *, SelectionParserValueListPointer, const SelectionLocation &)
        static SelectionParserParameter
        create(const std::string &name, SelectionParserValueListPointer values,
               const SelectionLocation &location)
        {
            return SelectionParserParameter(name.c_str(), move(values), location);
        }
        /*! \brief
         * Allocates and initializes a parsed method parameter.
         *
         * \param[in] name     Name for the new parameter (can be NULL).
         * \param[in] value    Value for the parameter.
         * \param[in] location Location of the parameter.
         * \returns   Pointer to the newly allocated parameter.
         * \throws    std::bad_alloc if out of memory.
         *
         * This overload is a convenience wrapper for the case when creating
         * parameters outside the actual Bison parser and only a single value
         * is necessary.
         */
        static SelectionParserParameter
        create(const char *name, const SelectionParserValue &value,
               const SelectionLocation &location)
        {
            return create(name, SelectionParserValue::createList(value), location);
        }
        /*! \brief
         * Allocates and initializes a parsed method parameter.
         *
         * \param[in] name    Name for the new parameter (can be NULL).
         * \param[in] expr    Expression value for the parameter.
         * \returns   Pointer to the newly allocated parameter.
         * \throws    std::bad_alloc if out of memory.
         *
         * This overload is a convenience wrapper for the case when creating
         * parameters outside the actual Bison parser and only a single
         * expression value is necessary.
         */
        static SelectionParserParameter
        createFromExpression(const char                        *name,
                             const SelectionTreeElementPointer &expr)
        {
            return create(name, SelectionParserValue::createExpr(expr),
                          expr->location());
        }
        //! \copydoc createFromExpression(const char *, const SelectionTreeElementPointer &)
        static SelectionParserParameter
        createFromExpression(const std::string                 &name,
                             const SelectionTreeElementPointer &expr)
        {
            return create(name.c_str(), SelectionParserValue::createExpr(expr),
                          expr->location());
        }

        //! Returns the name of the parameter (may be empty).
        const std::string &name() const { return name_; }
        //! Returns the location of this parameter in the parsed selection text.
        const SelectionLocation        &location() const { return location_; }
        //! Returns the values for the parameter.
        const SelectionParserValueList &values() const { return *values_; }

    private:
        /*! \brief
         * Initializes a parsed method parameter.
         *
         * \param[in] name     Name for the new parameter (can be NULL).
         * \param[in] values   List of values for the parameter.
         * \param[in] location Location of the parameter.
         * \throws    std::bad_alloc if out of memory.
         */
        SelectionParserParameter(const char                      *name,
                                 SelectionParserValueListPointer  values,
                                 const SelectionLocation         &location);

        //! Name of the parameter.
        std::string                     name_;
        //! Location of the parameter in the parsed text.
        SelectionLocation               location_;

        // TODO: Make private, there is only one direct user.
    public:
        //! Values for this parameter.
        SelectionParserValueListPointer values_;
};

} // namespace gmx

/*! \brief
 * Handles exceptions caught within the Bison code.
 *
 * \retval `true`  if the parser should attempt error recovery.
 * \retval `false` if the parser should immediately abort.
 *
 * This function is called whenever an exception is caught within Bison
 * actions.  Since exceptions cannot propagate through Bison code, the
 * exception is saved (potentially with some extra context information) so that
 * the caller of the parser can rethrow the exception.
 *
 * If it is possible to recover from the exception, then the function returns
 * `true`, and Bison enters error recovery state.  At the end of the recovery,
 * _gmx_selparser_handle_error() is called.
 * If this function returns false, then Bison immediately aborts the parsing
 * so that the caller can rethrow the exception.
 */
bool
_gmx_selparser_handle_exception(void *scanner, std::exception *ex);
/*! \brief
 * Handles errors in the selection parser.
 *
 * \returns `true` if parsing can continue with the next selection.
 * \throws  std::bad_alloc if out of memory during the error processing.
 * \throws  unspecified    Can throw the stored exception if recovery from that
 *     exception is not possible.
 *
 * This function is called during error recovery, after Bison has discarded all
 * the symbols for the erroneous selection.
 * At this point, the full selection that caused the error is known, and can be
 * added to the error context.
 *
 * For an interactive parser, this function returns `true` to let the parsing
 * continue with the next selection, or to let the user enter the next
 * selection, if it was possible to recover from the exception.
 * For other cases, this will either rethrow the original exception with added
 * context, or return `false` after adding the context to the error reporter.
 * Any exceptions thrown from this method are again caught by Bison and result
 * in termination of the parsing; the caller can then rethrow them.
 */
bool
_gmx_selparser_handle_error(void *scanner);

/** Propagates the flags for selection elements. */
void
_gmx_selelem_update_flags(const gmx::SelectionTreeElementPointer &sel);

/** Initializes the method parameter data of \ref SEL_EXPRESSION and
 * \ref SEL_MODIFIER elements. */
void
_gmx_selelem_init_method_params(const gmx::SelectionTreeElementPointer &sel,
                                void                                   *scanner);
/** Initializes the method for a \ref SEL_EXPRESSION selection element. */
void
_gmx_selelem_set_method(const gmx::SelectionTreeElementPointer &sel,
                        struct gmx_ana_selmethod_t *method, void *scanner);

/* An opaque pointer. */
#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif
/** \brief Creates a gmx::SelectionTreeElement for arithmetic expression evaluation.
 *
 * \param[in]  left    Selection element for the left hand side.
 * \param[in]  right   Selection element for the right hand side.
 * \param[in]  op      String representation of the operator.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * arithmetic expressions.
 */
gmx::SelectionTreeElementPointer
_gmx_sel_init_arithmetic(const gmx::SelectionTreeElementPointer &left,
                         const gmx::SelectionTreeElementPointer &right,
                         char op, yyscan_t scanner);
/** Creates a gmx::SelectionTreeElement for comparsion expression evaluation. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_comparison(const gmx::SelectionTreeElementPointer &left,
                         const gmx::SelectionTreeElementPointer &right,
                         const char *cmpop, void *scanner);
/** Creates a gmx::SelectionTreeElement for a keyword expression from the parsed data. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword(struct gmx_ana_selmethod_t *method,
                      gmx::SelectionParserValueListPointer args,
                      const char *rpost, void *scanner);
/** Creates a gmx::SelectionTreeElement for string-matching keyword expression. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword_strmatch(struct gmx_ana_selmethod_t *method,
                               gmx::SelectionStringMatchType matchType,
                               gmx::SelectionParserValueListPointer args,
                               const char *rpost, void *scanner);
/** Creates a gmx::SelectionTreeElement for "keyword of" expression. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword_of(struct gmx_ana_selmethod_t *method,
                         const gmx::SelectionTreeElementPointer &group,
                         const char *rpost, void *scanner);
/** Creates a gmx::SelectionTreeElement for a method expression from the parsed data. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_method(struct gmx_ana_selmethod_t *method,
                     gmx::SelectionParserParameterListPointer params,
                     const char *rpost, void *scanner);
/** Creates a gmx::SelectionTreeElement for a modifier expression from the parsed data. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_modifier(struct gmx_ana_selmethod_t              *mod,
                       gmx::SelectionParserParameterListPointer params,
                       const gmx::SelectionTreeElementPointer  &sel,
                       void                                    *scanner);
/** Creates a gmx::SelectionTreeElement for evaluation of reference positions. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_position(const gmx::SelectionTreeElementPointer &expr,
                       const char *type, void *scanner);

/** Creates a gmx::SelectionTreeElement for a constant position. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_const_position(real x, real y, real z, void *scanner);
/** Creates a gmx::SelectionTreeElement for a index group expression using group name. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_group_by_name(const char *name, void *scanner);
/** Creates a gmx::SelectionTreeElement for a index group expression using group index. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_group_by_id(int id, void *scanner);
/** Creates a gmx::SelectionTreeElement for a variable reference */
gmx::SelectionTreeElementPointer
_gmx_sel_init_variable_ref(const gmx::SelectionTreeElementPointer &sel,
                           void                                   *scanner);

/** Creates a root gmx::SelectionTreeElement for a selection. */
gmx::SelectionTreeElementPointer
_gmx_sel_init_selection(const char                             *name,
                        const gmx::SelectionTreeElementPointer &sel,
                        void                                   *scanner);
/** Creates a root gmx::SelectionTreeElement elements for a variable assignment. */
gmx::SelectionTreeElementPointer
_gmx_sel_assign_variable(const char                             *name,
                         const gmx::SelectionTreeElementPointer &expr,
                         void                                   *scanner);
/** Appends a root gmx::SelectionTreeElement to a selection collection. */
gmx::SelectionTreeElementPointer
_gmx_sel_append_selection(const gmx::SelectionTreeElementPointer &sel,
                          gmx::SelectionTreeElementPointer        last,
                          void                                   *scanner);
/** Check whether the parser should finish. */
bool
_gmx_sel_parser_should_finish(void *scanner);

/* In params.c */
/** Initializes an array of parameters based on input from the selection parser. */
void
_gmx_sel_parse_params(const gmx::SelectionParserParameterList &params,
                      int nparam, struct gmx_ana_selparam_t *param,
                      const gmx::SelectionTreeElementPointer &root,
                      void *scanner);

#endif
