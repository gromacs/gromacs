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
 * \brief API for handling selection (the \c gmx_ana_selection_t structure and related functions).
 *
 * There should be no need to use the data structures or call the
 * functions in this file directly unless using the selection routines outside
 * the main trajectory analysis API.
 */
#ifndef GMX_SELECTION_SELECTIONENUMS_H
#define GMX_SELECTION_SELECTIONENUMS_H

/** Defines the type of covered fraction. */
typedef enum
{
    CFRAC_NONE,         /**< No covered fraction (everything covered). */
    CFRAC_SOLIDANGLE    /**< Fraction of a solid (3D) angle covered. */
} e_coverfrac_t;

namespace gmx
{

enum SelectionFlag
{
    efOnlyStatic        = 1<<0,
    efOnlyAtoms         = 1<<1,
    efDynamicMask       = 1<<2,
    efDynamicOnlyWhole  = 1<<3,
    efCollectRemaining  = 1<<4,
};

typedef unsigned long SelectionFlags;

}

#endif
