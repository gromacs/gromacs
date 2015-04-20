/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
 * \brief Helper functions for the selection parser.
 *
 * This header is includes only from parser.cpp (generated from parser.y), and
 * it includes functions and macros used internally by the parser.
 * They are in a separate file to make then easier to edit (no need to
 * regenerate the parser), and to keep parser.y as simple as possible.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_PARSER_INTERNAL_H
#define GMX_SELECTION_PARSER_INTERNAL_H

#include "config.h"

#include <exception>

#include <boost/exception_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "parsetree.h"
#include "scanner.h"
#include "selelem.h"

//! Error handler needed by Bison.
static void
yyerror(YYLTYPE *location, yyscan_t scanner, char const *s)
{
    try
    {
        std::string            context(_gmx_sel_lexer_get_text(scanner, *location));
        gmx::InvalidInputError ex(s);
        // TODO: Examine how to show better context information.
        if (!context.empty())
        {
            context = gmx::formatString("Near '%s'", context.c_str());
            ex.prependContext(context);
        }
        _gmx_sel_lexer_set_exception(scanner, boost::copy_exception(ex));
    }
    catch (const std::exception &)
    {
        _gmx_sel_lexer_set_exception(scanner, boost::current_exception());
    }
}

//! Logic for computing the location of the output of Bison reduction.
#define YYLLOC_DEFAULT(Current, Rhs, N)                          \
    do {                                                         \
        if (N != 0)                                              \
        {                                                        \
            (Current).startIndex = YYRHSLOC(Rhs, 1).startIndex;  \
            (Current).endIndex   = YYRHSLOC(Rhs, N).endIndex;    \
        }                                                        \
        else                                                     \
        {                                                        \
            (Current).startIndex = (Current).endIndex =          \
                    YYRHSLOC(Rhs, 0).endIndex;                   \
        }                                                        \
        _gmx_sel_lexer_set_current_location(scanner, (Current)); \
    } while (0)

/*! \brief
 * Custom macro to influence Bison behavior.
 *
 * This macro added to parser.cpp through our patch to force Bison to
 * use C-style logic for stack reallocation even though we have provided
 * YYLTYPE and are compiling the code in C++ (our YYLTYPE can safely be copied
 * this way).
 * An alternative would be to provide the whole reallocation logic through an
 * undocumented yyoverflow() macro, but that is probably also more trouble than
 * it is worth.
 */
#define GMX_YYFORCE_C_STACK_EXTENSION 1

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
#define END_ACTION                                              \
    }                                                           \
    catch (std::exception &ex)                                  \
    {                                                           \
        if (_gmx_selparser_handle_exception(scanner, &ex))      \
        {                                                       \
            YYERROR;                                            \
        }                                                       \
        else                                                    \
        {                                                       \
            YYABORT;                                            \
        }                                                       \
    }
//! Finishes an action that may throw exceptions and does not support resuming.
#define END_ACTION_TOPLEVEL                                     \
    }                                                           \
    catch (const std::exception &)                              \
    {                                                           \
        _gmx_sel_lexer_set_exception(scanner, boost::current_exception()); \
        YYABORT;                                                \
    }
//!\}

#ifndef GMX_CXX11
//! No-op to enable use of same get()/set() implementation as with C++11.
static gmx::SelectionParserValue &move(gmx::SelectionParserValue &src)
{
    return src;
}
//! No-op to enable use of same get()/set() implementation as with C++11.
static gmx::SelectionParserParameter &move(gmx::SelectionParserParameter &src)
{
    return src;
}
#endif

/*! \brief
 * Retrieves a semantic value.
 *
 * \param[in] src  Semantic value to get the value from.
 * \returns   Retrieved value.
 * \throws    unspecified  Any exception thrown by the move constructor of
 *      ValueType (copy constructor if GMX_CXX11 is not set).
 *
 * There should be no statements that may throw exceptions in actions before
 * this function has been called for all semantic values that have a C++ object
 * stored.  Together with set(), this function abstracts away exception
 * safety issues that arise from the use of a plain pointer for storing the
 * semantic values.
 *
 * Does not throw for smart pointer types.  If used with types that may throw,
 * the order of operations should be such that it is exception-safe.
 */
template <typename ValueType> static
ValueType get(ValueType *src)
{
    GMX_RELEASE_ASSERT(src != NULL, "Semantic value pointers should be non-NULL");
    boost::scoped_ptr<ValueType> srcGuard(src);
    return ValueType(move(*src));
}
/*! \brief
 * Sets a semantic value.
 *
 * \tparam     ValueType Type of value to set.
 * \param[out] dest  Semantic value to set (typically $$).
 * \param[in]  value Value to put into the semantic value.
 * \throws     std::bad_alloc if out of memory.
 * \throws     unspecified  Any exception thrown by the move constructor of
 *      ValueType (copy constructor if GMX_CXX11 is not set).
 *
 * This should be the last statement before ::END_ACTION, except for a
 * possible ::CHECK_SEL.
 */
template <typename ValueType> static
void set(ValueType * &dest, ValueType value)
{
    dest = new ValueType(move(value));
}
/*! \brief
 * Sets an empty semantic value.
 *
 * \tparam     ValueType Type of value to set (must be default constructible).
 * \param[out] dest  Semantic value to set (typically $$).
 * \throws     std::bad_alloc if out of memory.
 * \throws     unspecified  Any exception thrown by the default constructor of
 *      ValueType.
 *
 * This should be the last statement before ::END_ACTION, except for a
 * possible ::CHECK_SEL.
 */
template <typename ValueType> static
void set_empty(ValueType * &dest)
{
    dest = new ValueType;
}
/*! \brief
 * Checks that a valid tree was set.
 *
 * Should be called after set() if it was used to set a value where NULL
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

#endif
