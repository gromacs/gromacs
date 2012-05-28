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
/*! \libinternal \file
 * \brief
 * Declares gmx::SelectionFileOption.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONFILEOPTION_H
#define GMX_SELECTION_SELECTIONFILEOPTION_H

#include "../options/abstractoption.h"

namespace gmx
{

/*! \libinternal \brief
 * Specifies a special option that provides selections from a file.
 *
 * This option is used internally by the command-line framework to implement
 * file input for selections.  The option takes a file name, and reads it in
 * using SelectionOptionManager::parseRequestedFromFile().  This means that
 * selections from the file are assigned to selection options that have been
 * explicitly provided without values earlier on the command line.
 *
 * Public methods in this class do not throw.
 *
 * \inlibraryapi
 * \ingroup module_selection
 */
class SelectionFileOption : public AbstractOption
{
    public:
        //! Initializes an option with the given name.
        explicit SelectionFileOption(const char *name);

    private:
        virtual AbstractOptionStoragePointer createStorage() const;
};

} // namespace gmx

#endif
