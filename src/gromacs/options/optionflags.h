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
/*! \file
 * \brief
 * Defines flags used in option implementation.
 *
 * Symbols in this header are considered an implementation detail, and should
 * not be accessed outside the module.
 * Because of details in the implementation, it is still installed.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONFLAGS_H
#define GMX_OPTIONS_OPTIONFLAGS_H

#include "../utility/flags.h"

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
    efOption_Set = 1<<0,
    //! The current value of the option is a programmatic default value.
    efOption_HasDefaultValue       = 1<<1,
    //! An explicit default value has been provided for the option.
    efOption_ExplicitDefaultValue  = 1<<2,
    /*! \brief
     * Next assignment to the option clears old values.
     *
     * This flag is set when a new option source starts, such that values
     * from the new source will overwrite old ones.
     */
    efOption_ClearOnNextSet        = 1<<3,
    //! %Option is required to be set.
    efOption_Required              = 1<<4,
    //! %Option can be specified multiple times.
    efOption_MultipleTimes         = 1<<5,
    //! %Option is hidden from standard help.
    efOption_Hidden                = 1<<6,
    /*! \brief
     * %Option value is a vector, but a single value is also accepted.
     *
     * \see AbstractOption::setVector()
     */
    efOption_Vector                = 1<<8,
    //! %Option does not support default values.
    efOption_NoDefaultValue        = 1<<9,
    /*! \brief
     * Storage object does its custom checking for minimum value count.
     *
     * If this flag is set, the class derived from OptionStorageTemplate should
     * implement processSetValues(), processAll(), and possible other functions
     * it provides such that it always fails if not enough values are provided.
     * This is useful to override the default check, which is done in
     * OptionStorageTemplate::processSet().
     */
    efOption_DontCheckMinimumCount = 1<<10
};

//! \libinternal Holds a combination of ::OptionFlag values.
typedef FlagsTemplate<OptionFlag> OptionFlags;
//! \endcond

} // namespace gmx

#endif
