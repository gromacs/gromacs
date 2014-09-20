/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2013,2014, by the GROMACS development team, led by
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
 * \brief API for handling parameters used in selections.
 *
 * There should be no need to use the data structures or call the
 * functions in this file directly unless implementing a custom selection
 * method.
 *
 * More details can be found on the page discussing
 * \ref page_module_selection_custom "custom selection methods".
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELPARAM_H
#define GMX_SELECTION_SELPARAM_H

#include "gromacs/selection/indexutil.h"

#include "selvalue.h"

/*! \name Parameter flags
 * \anchor selparam_flags
 */
/*@{*/
/*! \brief
 * This flag is set if the user has provided the parameter.
 *
 * This flag is set automatically, and should not be set by the user.
 */
#define SPAR_SET      1
/*! \brief
 * If not set, an error is reported if the parameter is not specified by the
 * user.
 */
#define SPAR_OPTIONAL 2
/*! \brief
 * If set, the parameter value can be dynamic, i.e., be different for
 * different frames.
 *
 * If set, the parameter value should only be accessed in the update function
 * of \c gmx_ana_selmethod_t.
 * The flag is cleared before sel_initfunc() if the value provided is actually
 * static.
 */
#define SPAR_DYNAMIC  4
/*! \brief
 * If set, the parameter value is parsed into sorted ranges.
 *
 * Can only be specified for integer parameters.
 * If specified, the value of the parameter (\c gmx_ana_selparam_t::val)
 * consists of sets of two integers, each specifying a range.
 * The values give the endpoints of the ranges (inclusive).
 * The ranges are sorted and overlapping/continuous ranges are merged into
 * a single range to minimize the number of ranges.
 *
 * If this flag is specified, \c gmx_ana_selparam_t::nval gives the number of
 * ranges. \p gmx_ana_selparam_t::nval should be 1 or \ref SPAR_VARNUM should be
 * specified; other values would lead to unpredictable behavior.
 */
#define SPAR_RANGES   8
/*! \brief
 * If set, the parameter can have any number of values.
 *
 * If specified, the data pointer in \c gmx_ana_selparam_t::val should be NULL;
 * the memory is allocated by the parameter parser.
 * The implementation of the method should ensure that the pointer to the
 * allocated memory is stored somewhere in sel_initfunc();
 * otherwise, the memory is lost.
 *
 * The initial value of \c gmx_ana_selparam_t::nval is not used with this flag.
 * Instead, it will give the number of actual values provided by the user
 * after the parameters have been parsed.
 * For consistency, it should be initialized to -1.
 *
 * Cannot be combined with \ref GROUP_VALUE parameters.
 */
#define SPAR_VARNUM   16
/*! \brief
 * If set, the parameter can have a separate value for each atom.
 *
 * The flag is cleared before sel_initfunc() if the value provided is actually
 * a single value.
 *
 * Cannot be combined with \ref POS_VALUE or \ref GROUP_VALUE parameters.
 */
#define SPAR_ATOMVAL  32
/*! \brief
 * If set, the parameter takes one of a set of predefined strings.
 *
 * Can only be specified for a \ref STR_VALUE parameter that takes a single
 * string.
 * The data pointer in \c gmx_ana_selparam_t::val should be initialized into an
 * array of strings such that the first and last elements are NULL, and the
 * rest give the possible values. For optional values, the second element in
 * the array should give the default value. The string given by the user is
 * matched against the beginnings of the given strings, and if a unique match
 * is found, the first pointer in the array will be initialized to point to
 * the matching string.
 * The data pointer can be initialized as a static array; duplication of the
 * array for multiple instances of the same method is automatically taken care
 * of.
 */
#define SPAR_ENUMVAL  128
/*@}*/

/*! \internal \brief
 * Describes a single parameter for a selection method.
 */
typedef struct gmx_ana_selparam_t
{
    /** Name of the parameter. */
    const char         *name;
    /*! \brief
     * The parameter value.
     *
     * Type \ref NO_VALUE can be used to define a boolean parameter.
     * The number of values should be 0 for boolean parameters.
     *
     * The value pointer be initialized to NULL in the definition of a
     * \c gmx_ana_selmethod_t and initialized in the
     * \c gmx_ana_selmethod_t::init_data call
     * (see sel_datafunc()).
     * However, if \ref SPAR_VARNUM is provided and the parameter is not
     * \ref POS_VALUE, this field should not be initialized. Instead,
     * sufficient memory is allocated automatically and the pointer should be
     * stored in \c gmx_ana_selmethod_t::init
     * (see sel_initfunc()).
     *
     * The values cannot be accessed outside these two functions: the compiler
     * makes a copy of the parameter structure for each instance of the
     * method, and the original parameter array is not changed.
     */
    gmx_ana_selvalue_t  val;
    /*! \brief
     * Pointer to store the number of values.
     *
     * If not NULL, the number of values for the parameter is stored in the
     * pointed value.
     * Should be specified if \ref SPAR_VARNUM and \ref SPAR_DYNAMIC are both
     * set.
     *
     * Should be initialized to NULL in the definition a \c gmx_ana_selmethod_t
     * and initialized in sel_datafunc().
     */
    int                *nvalptr;
    /*! \brief
     * Flags that alter the way the parameter is parsed/handled.
     *
     * See \ref selparam_flags for allowed values.
     */
    int                 flags;
} gmx_ana_selparam_t;

/** Finds a parameter from an array by name. */
gmx_ana_selparam_t *
gmx_ana_selparam_find(const char *name, int nparam, gmx_ana_selparam_t *param);

#endif
