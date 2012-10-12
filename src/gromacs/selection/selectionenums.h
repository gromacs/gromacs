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
 * Declares common types used in selections.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONENUMS_H
#define GMX_SELECTION_SELECTIONENUMS_H

#include "../utility/flags.h"

/*! \brief
 * Defines the type of covered fraction.
 *
 * \inpublicapi
 */
typedef enum
{
    CFRAC_NONE,         /**< No covered fraction (everything covered). */
    CFRAC_SOLIDANGLE    /**< Fraction of a solid (3D) angle covered. */
} e_coverfrac_t;

namespace gmx
{

/*! \cond internal */
/*! \internal \brief
 * Flags for options.
 *
 * These flags are not part of the public interface, even though they are in an
 * installed header.  They are needed in the implementation of SelectionOption.
 */
enum SelectionFlag
{
    efSelection_OnlyStatic              = 1<<0,
    efSelection_OnlyAtoms               = 1<<1,
    //! Whether ::POS_MASKONLY should be used for output position evaluation.
    efSelection_DynamicMask             = 1<<2,
    efSelection_DynamicOnlyWhole        = 1<<3,
    //! Whether velocities of output positions should be evaluated.
    efSelection_EvaluateVelocities      = 1<<5,
    //! Whether forces on output positions should be evaluated.
    efSelection_EvaluateForces          = 1<<6,
};

//! \internal Holds a collection of ::SelectionFlag values.
typedef FlagsTemplate<SelectionFlag> SelectionFlags;
//! \endcond

}

#endif
