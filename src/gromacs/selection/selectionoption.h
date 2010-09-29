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
 * Declares gmx::SelectionOption.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONOPTION_H
#define GMX_SELECTION_SELECTIONOPTION_H

#include "../options/abstractoption.h"
#include "selectionenums.h"

namespace gmx
{

class Selection;
class SelectionOptionStorage;

/*! \brief
 * Specifies an option that provides selection(s).
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionOption : public OptionTemplate<Selection *, SelectionOption>
{
    public:
        //! Initializes an option with the given name.
        explicit SelectionOption(const char *name) : MyBase(name), _flags(0)
        {
            setFlag(efConversionMayNotAddValues);
        }

        /*! \brief
         * Only accept selections that evaluate to atom positions.
         */
        MyClass &onlyAtoms() { _flags |= efOnlyAtoms; return me(); }
        /*! \brief
         * Only accept static selections for this option.
         */
        MyClass &onlyStatic() { _flags |= efOnlyStatic; return me(); }
        /*! \brief
         * Handle dynamic selections for this option with position masks.
         */
        MyClass &dynamicMask() { _flags |= efDynamicMask; return me(); }
        /*! \brief
         * Disallow using atom coordinates as the reference positions.
         */
        MyClass &dynamicOnlyWhole() { _flags |= efDynamicOnlyWhole; return me(); }
        /*! \brief
         * Mark this option as receiving selections that can't be assigned to
         * other options.
         */
        MyClass &collectRemaining() { _flags |= efCollectRemaining; return me(); }

    private:
        virtual int createDefaultStorage(Options *options,
                                         AbstractOptionStorage **storage) const;

        SelectionFlags          _flags;

        /*! \brief
         * Needed to initialize SelectionOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class SelectionOptionStorage;
};

} // namespace gmx

#endif
