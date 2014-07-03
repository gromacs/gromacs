/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Defines flags used in option implementation.
 *
 * Symbols in this header are considered an implementation detail, and should
 * not be accessed outside the module.
 * Because of details in the implementation, it is still installed.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONFLAGS_H
#define GMX_OPTIONS_OPTIONFLAGS_H

#include "gromacs/utility/flags.h"

namespace gmx
{

/*! \cond libapi */
/*! \libinternal \brief
 * Flags for options.
 *
 * These flags are not part of the public interface, even though they are in an
 * installed header.  They are needed in a few template class implementations.
 *
 * \todo
 * The flags related to default values are confusing, consider reorganizing
 * them.
 */
enum OptionFlag
{
    //! %Option has been set.
    efOption_Set                        = 1<<0,
    //! The current value of the option is a programmatic default value.
    efOption_HasDefaultValue            = 1<<1,
    //! An explicit default value has been provided for the option.
    efOption_ExplicitDefaultValue       = 1<<2,
    /*! \brief
     * Next assignment to the option clears old values.
     *
     * This flag is set when a new option source starts, such that values
     * from the new source will overwrite old ones.
     */
    efOption_ClearOnNextSet             = 1<<3,
    //! %Option is required to be set.
    efOption_Required                   = 1<<4,
    //! %Option can be specified multiple times.
    efOption_MultipleTimes              = 1<<5,
    //! %Option is hidden from standard help.
    efOption_Hidden                     = 1<<6,
    /*! \brief
     * %Option value is a vector, but a single value is also accepted.
     *
     * \see AbstractOption::setVector()
     */
    efOption_Vector                     = 1<<8,
    //! %Option has a defaultValueIfSet() specified.
    efOption_DefaultValueIfSetExists    = 1<<11,
    //! %Option does not support default values.
    efOption_NoDefaultValue             = 1<<9,
    /*! \brief
     * Storage object does its custom checking for minimum value count.
     *
     * If this flag is set, the class derived from OptionStorageTemplate should
     * implement processSetValues(), processAll(), and possible other functions
     * it provides such that it always fails if not enough values are provided.
     * This is useful to override the default check, which is done in
     * OptionStorageTemplate::processSet().
     */
    efOption_DontCheckMinimumCount      = 1<<10
};

//! \libinternal Holds a combination of ::OptionFlag values.
typedef FlagsTemplate<OptionFlag> OptionFlags;
//! \endcond

} // namespace gmx

#endif
