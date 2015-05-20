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
/*! \internal
 * \page page_module_selection_custom Custom selection methods
 *
 * Custom selection methods are defined by creating a new instance of
 * \c gmx_ana_selmethod_t and filling it with the necessary data for handling
 * the selection.
 * The structure contains callback pointers that define the actual behavior
 * of the method.
 * The following sections discuss how the structure should be filled and how
 * to implement the callbacks.
 *
 *
 * \section selmethods_define \c gmx_ana_selmethod_t data structure
 *
 * An example \c gmx_ana_selmethod_t definition could look like this:
 *
 * \code
   gmx_ana_selmethod_t sm_example = {
       "example", GROUP_VALUE, 0,
       asize(sm_params_example), sm_params_example,
       &init_data_example,
        NULL,
       &init_example,
        NULL,
       &free_data_example,
       &init_frame_example,
       &evaluate_example,
        NULL,
       {"example from POS_EXPR [cutoff REAL]", NULL, 0, NULL},
   };
 * \endcode
 *
 * The first value defines the name of the method.
 * It is used mostly for informational purposes; the actual name(s) recognized
 * by the selection parser are defined by the call to
 * gmx_ana_selmethod_register() (see \ref selmethods_register).
 *
 * The second value defines the type of the value the method returns.
 * Possible values are
 *  - \ref NO_VALUE : This is allowed only for methods that have the flag
 *    \ref SMETH_MODIFIER set (see \ref selmethods_modifiers).
 *  - \ref INT_VALUE : The method returns one or more integer values.
 *  - \ref REAL_VALUE : The method returns one or more floating-point values.
 *  - \ref STR_VALUE : The method returns one or more strings.
 *  - \ref POS_VALUE : The method returns one or more 3D vectors.
 *  - \ref GROUP_VALUE : The method returns a single index group.
 *
 * The third value gives additional information about the method using
 * a combination of flags.
 * Possible flags are:
 *  - \ref SMETH_REQTOP : If set, the topology information is always loaded
 *    and the \p top pointer passed to the callbacks is guaranteed to be
 *    non-NULL. Should be set if the method requires topology information
 *    for evaluation.
 *  - \ref SMETH_DYNAMIC : If set, the method can only be evaluated dynamically,
 *    i.e., it requires data from the trajectory frame.
 *  - \ref SMETH_MODIFIER : If set, the method is a selection modifier and
 *    not an actual selection method.
 *    For more details, see \ref selmethods_modifiers.
 *  - \ref SMETH_ALLOW_UNSORTED : If set, the method supports unsorted atoms
 *    in its input parameters. \ref SMETH_MODIFIER methods are assumed to always
 *    support unsorted atoms, as their purpose is to affect the ordering.
 *
 * There are two additional flags that specify the number of values the
 * method returns. Only one of them can be set at a time.
 * If neither is set, the default behavior is to evaluate a value for each
 * input atom (except for \ref GROUP_VALUE methods, which always return a
 * single group).
 * Other behaviors can be specified with these flags:
 *  - \ref SMETH_SINGLEVAL : If set, the method evaluates to a single value.
 *    This is automatically set if the type is \ref GROUP_VALUE.
 *  - \ref SMETH_VARNUMVAL : If set, the method evaluates to an arbitrary
 *    number of values.
 *    The number of values is determined based on the values given by the user
 *    to the method parameters (see \ref selmethods_params).
 *  .
 * If either of these flags is specified (and the method type is not
 * \ref GROUP_VALUE), the group passed to the evaluation callback should not
 * be used as it can be NULL.
 * Currently, the above flags only work (have been tested) for \ref POS_VALUE
 * methods.
 *
 * There is one additional flag that can only be specified for \ref STR_VALUE
 * methods: \ref SMETH_CHARVAL . It is meant for to ease implementation of
 * methods that evaluate to strings consisting of single characters.
 *
 * The next two values determine the number of parameters and a pointer to
 * the parameter array. The contents of the parameter array are described in
 * \ref selmethods_params. If the method does not take parameters, the first
 * value should be 0 and the second can be NULL.
 * Currently, \ref STR_VALUE methods cannot take parameters, but this limitation
 * should be easy to lift if required.
 *
 * These are followed by function callbacks that determine the
 * actual behavior of the method. Any of these except the evaluation callback
 * can be NULL (the evaluation callback can also be NULL if \ref NO_VALUE is
 * specified for a selection modifier). However, the presence of parameters
 * can require some of the callbacks to be implemented.
 * The details are described in \ref selmethods_callbacks.
 *
 * Finally, there is a data structure that gives help texts for the method.
 *
 * The \c gmx_ana_selmethod_t variable should be declared as a global variable
 * or it should be otherwise ensured that the structure is not freed: only a
 * pointer to the structure is stored by the library.
 *
 *
 * \section selmethods_params Defining parameters
 *
 * Parameters to selection methods are defined in a separate array of
 * \c gmx_ana_selparam_t structures.
 * The order of the parameters does not matter (except possibly for callback
 * implementation), with one important exception:
 * If the method evaluates to a \ref POS_VALUE, the first parameter should
 * have \ref GROUP_VALUE and be the one that is used to calculate the
 * positions.
 *
 * An example parameter definition:
 * \code
   static gmx_ana_selparam_t sm_params_example[] = {
     {"cutoff", {REAL_VALUE, 1, {NULL}}, NULL, SPAR_OPTIONAL},
     {"from",   {POS_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
   };
 * \endcode
 *
 * The first value gives the name of the parameter.
 * The first parameter can have a NULL name, which means that the value should
 * immediately follow the method name. This can be used to specify methods
 * of the type 'within 5 of ...'.
 *
 * The second value specifies the type of the value that the parameter accepts.
 * \ref NO_VALUE can be used to specify a boolean parameter, other possibilities
 * are the same as for the selection method type.
 *
 * The third value gives the number of values that the parameter accepts.
 * For boolean parameters (\ref NO_VALUE), it should be 0.
 * For parameters with \ref SPAR_VARNUM of \ref SPAR_ATOMVAL, it should be set
 * to -1 for consistency (it is not used).
 * If \ref SPAR_RANGES is specified, it should be either 1 (to accept a single
 * continuous range) or -1 (if combined with \ref SPAR_VARNUM).
 * In all other cases, it should be a positive integer; in most cases, it
 * should be 1.
 *
 * The nest two pointers should always be NULL (they should be initialized in
 * the callbacks), except the first pointer in the case of \ref SPAR_ENUMVAL
 * (see below).
 *
 * The final value gives additional information about the acceptable values
 * for the parameter using a combination of flags.
 * The possible flags are:
 *  - \ref SPAR_OPTIONAL : If set, the user does not need to provide a value
 *    for the parameter. If not set, an error is reported if the parameter
 *    is not specified by the user.
 *  - \ref SPAR_DYNAMIC : If set, the method can handle dynamic values for
 *    the parameter, i.e., the value(s) can be given by an expression that
 *    evaluates to different values for different frames.
 *  - \ref SPAR_RANGES : Can be set only for \ref INT_VALUE and
 *    \ref REAL_VALUE parameters,
 *    and cannot be combined with \ref SPAR_DYNAMIC.
 *    If set, the parameter accepts ranges of values.
 *    The ranges are automatically sorted and compacted such that a minimum
 *    amount of non-overlapping ranges are given for the method.
 *  - \ref SPAR_VARNUM : If set, the parameter can have a variable number
 *    of values. These can be provided by the user as a list of values, or
 *    using a single \ref SMETH_VARNUMVAL (or a single \ref SMETH_SINGLEVAL)
 *    method.
 *  - \ref SPAR_ATOMVAL : If set, the parameter accepts either a single value
 *    or an expression that evaluates to a value for each input atom.
 *    The single input value is treated as if the same value was returned for
 *    each atom.
 *    Cannot be combined with \ref SPAR_RANGES or \ref SPAR_VARNUM.
 *  - \ref SPAR_ENUMVAL : Can only be set for \ref STR_VALUE parameters that
 *    take a single value, and cannot be combined with any other flag than
 *    \ref SPAR_OPTIONAL. If set, the parameter only accepts one of predefined
 *    string values. See \ref SPAR_ENUMVAL documentation for details on how
 *    to specify the acceptable values.
 *
 *
 * \section selmethods_callbacks Implementing callbacks
 *
 * There are eight differen callback functions that can be implemented for
 * selection methods: sel_datafunc(), sel_posfunc(), sel_initfunc(),
 * sel_outinitfunc(), sel_freefunc(), sel_framefunc(), and two update functions.
 * They are in this order in the \c gmx_ana_selmethod_t data structure.
 * In general, any of the callbacks can be NULL, but the presence of
 * parameters or other callbacks imposes some restrictions:
 *  - sel_datafunc() should be provided if the method takes parameters.
 *  - sel_initfunc() should be provided if the method takes
 *    any parameters with the \ref SPAR_VARNUM or \ref SPAR_ATOMVAL flags,
 *    except if those parameters have a \ref POS_VALUE.
 *  - sel_outinitfunc() should be provided for \ref POS_VALUE methods
 *    and \ref SMETH_VARNUMVAL methods.
 *  - sel_freefunc() should be provided if sel_datafunc() and/or
 *    sel_initfunc() allocate any dynamic memory in addition to the data
 *    structure itself (or allocates the data structure using some other means
 *    than malloc()).
 *  - sel_updatefunc_pos() only makes sense for methods with \ref SMETH_DYNAMIC
 *    set.
 *  - At least one update function should be provided unless the method type is
 *    \ref NO_VALUE.
 *
 * The documentations for the function pointer types provide more information
 * about how the callbacks should be implemented.
 *
 *
 * \section selmethods_modifiers Selection modifiers
 *
 * Selection modifiers are a special kind of selection methods that can be
 * appended to the end of a selection. They are specified by adding the
 * \ref SMETH_MODIFIER flag to the \c gmx_ana_selmethod_t.
 * They can have two different types:
 *  - \ref POS_VALUE : These modifiers are given the final positions
 *    as an input, and they can make modifications to the selection that are
 *    not possible otherwise (e.g., permute the atoms).
 *    The modifier should implement sel_updatefunc_pos() and also have
 *    one NULL parameter in the beginning of the parameter list that takes
 *    \ref POS_VALUE and is used to give the input positions.
 *  - \ref NO_VALUE : These modifiers do not modify the final selection, but
 *    can be used to implement per-selection options for analysis tools
 *    or to control the default behavior of the selection engine
 *    (currently, such a framework is not implemented, but should be easy to
 *    implement if required).
 *
 * In addition to restricting the type of the method, selection modifiers
 * do not allow the flags \ref SMETH_SINGLEVAL and \ref SMETH_VARNUMVAL
 * (they would not make sense).
 *
 * Parameters and callbacks should be implemented as with normal selection
 * method, but beware that very little of the functionality has been tested.
 *
 * \todo
 * The modifier handling could be made more flexible and more generic;
 * the current implementation does not allow many things which would be
 * possible with slight changes in the internals of the library.
 *
 *
 * \section selmethods_register Registering the method
 *
 * After defining the method with \c gmx_ana_selmethod_t, it should be
 * registered with the selection engine.
 * In analysis programs, this can be done by calling
 * gmx_ana_selmethod_register().
 * If adding the method to the library, you should add a pointer to the new
 * method structure into the \c smtable_def array (in selmethod.cpp), and it is
 * registered automatically.
 * In both cases, gmx_ana_selmethod_register() does several checks on the
 * structure and reports any errors or inconsistencies it finds.
 */
/*! \internal \file
 * \brief API for handling selection methods.
 *
 * There should be no need to use the data structures or call the
 * functions in this file directly unless implementing a custom selection
 * method.
 *
 * Instructions for implementing custom selection methods can be found
 * on a separate page: \ref page_module_selection_custom
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELMETHOD_H
#define GMX_SELECTION_SELMETHOD_H

#include "selparam.h"
#include "selvalue.h"

namespace gmx
{
class PositionCalculationCollection;
class SelectionParserSymbolTable;
} // namespace gmx

struct gmx_ana_index_t;
struct gmx_ana_pos_t;
struct gmx_ana_selcollection_t;
struct t_pbc;
struct t_topology;
struct t_trxframe;

/*! \name Selection method flags
 * \anchor selmethod_flags
 */
/*@{*/
/*! \brief
 * If set, the method requires topology information.
 */
#define SMETH_REQTOP     1
/*! \brief
 * If set, the method can only be evaluated dynamically.
 */
#define SMETH_DYNAMIC    2
/*! \brief
 * If set, the method evaluates to a single value.
 *
 * The default is that the method evaluates to a value for each input atom.
 * Cannot be combined with \ref SMETH_VARNUMVAL.
 */
#define SMETH_SINGLEVAL  4
/*! \brief
 * If set, the method evaluates to an arbitrary number of values.
 *
 * The default is that the method evaluates to a value for each input atom.
 * Cannot be combined with \ref SMETH_SINGLEVAL or with \ref GROUP_VALUE.
 */
#define SMETH_VARNUMVAL  8
/*! \brief
 * If set, the method evaluates to single-character strings.
 *
 * This flag can only be set for \ref STR_VALUE methods. If it is set, the
 * selection engine automatically allocates and frees the required strings.
 * The evaluation function should store the character values as the first
 * character in the strings in the output data structure and should not change
 * the string pointers.
 */
#define SMETH_CHARVAL    64
/*! \brief
 * If set, the method accepts unsorted atoms in its input parameters.
 *
 * Currently, the support for this functionality is fairly limited, and only
 * static index group references can actually contain unsorted atoms.
 * But to make this single case work, the position evaluation must support
 * unsorted atoms as well.
 */
#define SMETH_ALLOW_UNSORTED 128
/*! \brief
 * If set, the method is a selection modifier.
 *
 * The method type should be \ref GROUP_VALUE or \ref NO_VALUE .
 * Cannot be combined with \ref SMETH_SINGLEVAL or \ref SMETH_VARNUMVAL .
 */
#define SMETH_MODIFIER   256
/*@}*/

/*! \brief
 * Allocates and initializes internal data and parameter values.
 *
 * \param[in]     npar  Number of parameters in \p param.
 * \param[in,out] param Pointer to (a copy of) the method's
 *   \c gmx_ana_selmethod_t::param.
 * \returns       Pointer to method-specific data structure.
 *   This pointer will be passed as the last parameter of all other function
 *   calls.
 * \throws        unspecified Any errors should be indicated by throwing an
 *      exception.
 *
 * Should allocate and initialize any internal data required by the method.
 * Should also initialize the value pointers (\c gmx_ana_selparam_t::val) in
 * \p param to point to variables within the internal data structure,
 * with the exception of parameters that specify the \ref SPAR_VARNUM or
 * the \ref SPAR_ATOMVAL flag (these should be handled in sel_initfunc()).
 * However, parameters with a position value should be initialized.
 * It is also possible to initialize \ref SPAR_ENUMVAL statically outside
 * this function (see \ref SPAR_ENUMVAL).
 * The \c gmx_ana_selparam_t::nvalptr should also be initialized for
 * non-position-valued parameters that have both \ref SPAR_VARNUM and
 * \ref SPAR_DYNAMIC set (it can also be initialized for other parameters if
 * desired, but the same information will be available through other means).
 * For optional parameters, the default values can (and should) be initialized
 * here, as the parameter values are not changed if the parameter is not
 * provided.
 *
 * For boolean parameters (type equals \ref NO_VALUE), the default value
 * should be set here. The user can override the value by giving the parameter
 * either as 'NAME'/'noNAME', or as 'NAME on/off/yes/no'.
 *
 * If the method takes any parameters, this function must be provided.
 */
typedef void *(*sel_datafunc)(int npar, gmx_ana_selparam_t *param);
/*! \brief
 * Sets the position calculation collection for the method.
 *
 * \param[in]  pcc   Position calculation collection that the method should use
 *   for position calculations.
 * \param      data  Internal data structure from sel_datafunc().
 *
 * This function should be provided if the method uses the routines from
 * poscalc.h for calculating positions.
 * The pointer \p pcc should then be stored and used for initialization for
 * any position calculation structures.
 */
typedef void  (*sel_posfunc)(gmx::PositionCalculationCollection *pcc, void *data);
/*! \brief
 * Does initialization based on topology and/or parameter values.
 *
 * \param[in]  top   Topology structure
 *   (can be NULL if \ref SMETH_REQTOP is not set).
 * \param[in]  npar  Number of parameters in \p param.
 * \param[in]  param Pointer to (an initialized copy of) the method's
 *   \c gmx_ana_selmethod_t::param.
 * \param      data  Internal data structure from sel_datafunc().
 * \returns    0 on success, a non-zero error code on failure.
 *
 * This function is called after the parameters have been processed:
 * the values of the parameters are stored at the locations set in
 * sel_datafunc().
 * The flags \ref SPAR_DYNAMIC and \ref SPAR_ATOMVAL are cleared before
 * calling the function if the value is static or single-valued, respectively.
 * If a parameter had the \ref SPAR_VARNUM or \ref SPAR_ATOMVAL flag (and
 * is not \ref POS_VALUE), a pointer to the memory allocated for the values is
 * found in \c gmx_ana_selparam_t::val.
 * The pointer should be stored by this function, otherwise the values
 * cannot be accessed.
 * For \ref SPAR_VARNUM parameters, the number of values can be accessed
 * through \c gmx_ana_selparam_t::val. For parameters with \ref SPAR_DYNAMIC,
 * the number is the maximum number of values (the actual number can be
 * accessed in sel_framefunc() and in the update callback through the value
 * pointed by \c gmx_ana_selparam_t::nvalptr).
 * For \ref SPAR_ATOMVAL parameters, \c gmx_ana_selparam_t::val::nr is set to
 * 1 if a single value was provided, otherwise it is set to the maximum number
 * of values possibly passed to the method.
 * The value pointed by \c gmx_ana_selparam_t::nvalptr always contains the same
 * value as \c gmx_ana_selparam_t::val::nr.
 *
 * For dynamic \ref GROUP_VALUE parameters (\ref SPAR_DYNAMIC set), the value
 * will be the largest possible selection that may occur during the
 * evaluation. For other types of dynamic parameters, the values are
 * undefined.
 *
 * If the method takes any parameters with the \ref SPAR_VARNUM or
 * \ref SPAR_ATOMVAL flags, this function must be provided, except if these
 * parameters all have \ref POS_VALUE.
 *
 * This function may be called multiple times for the same method if the
 * method takes parameters with \ref SPAR_ATOMVAL set.
 */
typedef void  (*sel_initfunc)(t_topology *top, int npar,
                              gmx_ana_selparam_t *param, void *data);
/*! \brief
 * Initializes output data structure.
 *
 * \param[in]     top   Topology structure
 *   (can be NULL if \ref SMETH_REQTOP is not set).
 * \param[in,out] out   Output data structure.
 * \param[in]     data  Internal data structure from sel_datafunc().
 * \returns       0 on success, an error code on error.
 *
 * This function is called immediately after sel_initfunc().
 *
 * If the method evaluates to a position (\ref POS_VALUE), this function
 * should be provided, and it should initialize the \c gmx_ana_pos_t data
 * structure pointed by \p out.p (the pointer is guaranteed to be non-NULL).
 * The \p out.p->g pointer should be initialized to the group that is used
 * to evaluate positions in sel_updatefunc() or sel_updatefunc_pos().
 *
 * The function should also be provided for non-position-valued
 * \ref SMETH_VARNUMVAL methods. For these methods, it suffices to set the
 * \p out->nr field to reflect the maximum number of values returned by the
 * method.
 *
 * Currently, this function is not needed for other types of methods.
 *
 * This function may be called multiple times for the same method if the
 * method takes parameters with \ref SPAR_ATOMVAL set.
 */
typedef void  (*sel_outinitfunc)(t_topology *top, gmx_ana_selvalue_t *out,
                                 void *data);
/*! \brief
 * Frees the internal data.
 *
 * \param[in] data Internal data structure from sel_datafunc().
 *
 * This function should be provided if the internal data structure contains
 * dynamically allocated data, and should free any such data.
 * The data structure itself should also be freed.
 * For convenience, if there is no dynamically allocated data within the
 * structure and the structure is allocated using malloc()/snew(), this
 * function is not needed: the selection engine automatically frees the
 * structure using sfree().
 * Any memory pointers received as values of parameters are managed externally,
 * and should not be freed.
 * Pointers set as the value pointer of \ref SPAR_ENUMVAL parameters should not
 * be freed.
 */
typedef void  (*sel_freefunc)(void *data);

/*! \brief
 * Initializes the evaluation for a new frame.
 *
 * \param[in]  top  Topology structure
 *   (can be NULL if \ref SMETH_REQTOP is not set).
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  Initialized periodic boundary condition structure,
 *   or NULL if PBC should not be used.
 * \param      data Internal data structure from sel_datafunc().
 * \returns    0 on success, a non-zero error code on failure.
 *
 * This function should be implemented if the selection method needs to
 * do some preprocessing for each frame, and the preprocessing does not
 * depend on the evaluation group.
 * Because \p sel_updatefunc_* can be called more than once for a frame,
 * it is inefficient do the preprocessing there.
 * It is ensured that this function will be called before
 * \p sel_updatefunc_* for each frame, and that it will be called at most
 * once for each frame.
 * For static methods, it is called once, with \p fr and \p pbc set to
 * NULL.
 */
typedef void  (*sel_framefunc)(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                               void *data);
/*! \brief
 * Evaluates a selection method.
 *
 * \param[in]  top  Topology structure
 *   (can be NULL if \ref SMETH_REQTOP is not set).
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  Initialized periodic boundary condition structure,
 *   or NULL if PBC should not be used.
 * \param[in]  g    Index group for which the method should be evaluated.
 * \param[out] out  Output data structure.
 * \param      data Internal data structure from sel_datafunc().
 * \returns    0 on success, a non-zero error code on error.
 *
 * This function should evaluate the method for each atom included in \p g,
 * and write the output to \p out. The pointer in the union \p out->u that
 * corresponds to the type of the method should be used.
 * Enough memory has been allocated to store the output values.
 * The number of values in \p out should also be updated if necessary.
 * However, \ref POS_VALUE or \ref GROUP_VALUE methods should not touch
 * \p out->nr (it should be 1 anyways).
 *
 * For \ref STR_VALUE methods, the pointers stored in \p out->s are discarded
 * without freeing; it is the responsibility of this function to provide
 * pointers that can be discarded without memory leaks.
 *
 * If the method accesses \p fr outside the index group specified in \p g or
 * what it receives from its parameters, it must check that \p fr actually
 * contains such an atom in case the \p fr has been loaded from a trajectory
 * that only contains a subset of the system.
 */
typedef void  (*sel_updatefunc)(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                                gmx_ana_index_t *g, gmx_ana_selvalue_t *out,
                                void *data);
/*! \brief
 * Evaluates a selection method using positions.
 *
 * \param[in]  top  Topology structure
 *   (can be NULL if \ref SMETH_REQTOP is not set).
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  Initialized periodic boundary condition structure,
 *   or NULL if PBC should not be used.
 * \param[in]  pos  Positions for which the method should be evaluated.
 * \param[out] out  Output data structure.
 * \param      data Internal data structure from sel_datafunc().
 * \returns    0 on success, a non-zero error code on error.
 *
 * This function should evaluate the method for each position in \p pos,
 * and write the output values to \p out. The pointer in the union \p out->u
 * that corresponds to the type of the method should be used.
 * Enough memory has been allocated to store the output values.
 * The number of values in \p out should also be updated if necessary.
 * \ref POS_VALUE or \ref GROUP_VALUE methods should not touch
 * \p out->nr (it should be 1 anyways).  For other types of methods, the number
 * of output values should equal the number of positions in \p pos.
 *
 * For \ref STR_VALUE methods, the pointers stored in \p out->s are discarded
 * without freeing; it is the responsibility of this function to provide
 * pointers that can be discarded without memory leaks.
 *
 * If the method accesses \p fr outside the atoms referenced in \p pos or
 * what it receives from its parameters, it must check that \p fr actually
 * contains such an atom in case the \p fr has been loaded from a trajectory
 * that only contains a subset of the system.
 */
typedef void  (*sel_updatefunc_pos)(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                                    gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out,
                                    void *data);

/*! \internal
 * \brief
 * Help information for a selection method.
 *
 * If some information is not available, the corresponding field can be set to
 * 0/NULL.
 */
struct gmx_ana_selmethod_help_t
{
    /*! \brief
     * One-line description of the syntax of the method.
     *
     * If NULL, the name of the method is used.
     */
    const char         *syntax;
    /*! \brief
     * Title for the help text in \p help.
     *
     * If NULL, the name of the method is used.
     * Only used if `nlhelp > 0`.
     */
    const char         *helpTitle;
    /*! \brief
     * Number of strings in \p help.
     *
     * Set to 0 if \p help is NULL.
     */
    int                 nlhelp;
    /*! \brief
     * Detailed help for the method.
     *
     * If there is no help available in addition to \p syntax, this can be set
     * to NULL.
     */
    const char *const  *help;
};

/*! \internal
 * \brief
 * Describes a selection method.
 *
 * Any of the function pointers except the update call can be NULL if the
 * operation is not required or not supported. In this case,
 * corresponding function calls are skipped.
 *
 * See the function pointer type documentation for details of how the
 * functions should be implemented.
 * More details on implementing new selection methods can be found on a
 * separate page: \ref page_module_selection_custom.
 */
struct gmx_ana_selmethod_t
{
    /** Name of the method. */
    const char         *name;
    /** Type which the method returns. */
    e_selvalue_t        type;
    /*! \brief
     * Flags to specify how the method should be handled.
     *
     * See \ref selmethod_flags for allowed values.
     */
    int                 flags;
    /** Number of parameters the method takes. */
    int                 nparams;
    /** Pointer to the array of parameter descriptions. */
    gmx_ana_selparam_t *param;

    /** Function for allocating and initializing internal data and parameters. */
    sel_datafunc        init_data;
    /** Function to set the position calculation collection. */
    sel_posfunc         set_poscoll;
    /** Function to do initialization based on topology and/or parameter values. */
    sel_initfunc        init;
    /** Function to initialize output data structure. */
    sel_outinitfunc     outinit;
    /** Function to free the internal data. */
    sel_freefunc        free;

    /** Function to initialize the calculation for a new frame. */
    sel_framefunc       init_frame;
    /** Function to evaluate the value. */
    sel_updatefunc      update;
    /** Function to evaluate the value using positions. */
    sel_updatefunc_pos  pupdate;

    /** Help data for the method. */
    gmx_ana_selmethod_help_t help;
};

/** Registers a selection method. */
int
gmx_ana_selmethod_register(gmx::SelectionParserSymbolTable *symtab,
                           const char *name, gmx_ana_selmethod_t *method);
/** Registers all selection methods in the library. */
int
gmx_ana_selmethod_register_defaults(gmx::SelectionParserSymbolTable *symtab);

/** Finds a parameter from a selection method by name. */
gmx_ana_selparam_t *
gmx_ana_selmethod_find_param(const char *name, gmx_ana_selmethod_t *method);

#endif
