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
#include "selection.h"
#include "selectionenums.h"
#include "selectionoptioninfo.h"

namespace gmx
{

class SelectionOptionInfo;
class SelectionOptionStorage;

/*! \brief
 * Specifies an option that provides selection(s).
 *
 * Public methods in this class do not throw.
 *
 * \todo
 * Support for specifying that an option accepts, e.g., two to four selections.
 * Currently, only a fixed count or any number of selections is possible.
 * \if internal
 * In addition to allowing this in OptionTemplate, also SelectionOptionManager
 * needs to be updated.
 * \endif
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionOption : public OptionTemplate<Selection, SelectionOption>
{
    public:
        //! OptionInfo subclass corresponding to this option type.
        typedef SelectionOptionInfo InfoType;

        //! Initializes an option with the given name.
        explicit SelectionOption(const char *name) : MyBase(name) { }

        /*! \brief
         * Request velocity evaluation for output positions.
         *
         * Note that even with this flag set, velocities may not be available,
         * in which case SelectionPosition::hasVelocity() returns false.
         */
        MyClass &evaluateVelocities()
        { selectionFlags_.set(efSelection_EvaluateVelocities); return me(); }
        /*! \brief
         * Request force evaluation for output positions.
         *
         * Note that even with this flag set, forces may not be available,
         * in which case SelectionPosition::hasForce() returns false.
         */
        MyClass &evaluateForces()
        { selectionFlags_.set(efSelection_EvaluateForces); return me(); }
        /*! \brief
         * Only accept selections that evaluate to atom positions.
         *
         * TODO: This option is not yet implemented.
         */
        MyClass &onlyAtoms()
        { selectionFlags_.set(efSelection_OnlyAtoms); return me(); }
        /*! \brief
         * Only accept static selections for this option.
         */
        MyClass &onlyStatic()
        { selectionFlags_.set(efSelection_OnlyStatic); return me(); }
        /*! \brief
         * Handle dynamic selections for this option with position masks.
         *
         * \see Selection
         * \see SelectionPosition::selected()
         */
        MyClass &dynamicMask()
        { selectionFlags_.set(efSelection_DynamicMask); return me(); }
        /*! \brief
         * Disallow using atom coordinates as the reference positions.
         *
         * TODO: This option is not yet implemented.
         */
        MyClass &dynamicOnlyWhole()
        { selectionFlags_.set(efSelection_DynamicOnlyWhole); return me(); }

    private:
        // Disable possibility to allow multiple occurrences, since it isn't
        // implemented.
        using MyBase::allowMultiple;
        // Disable default value because it is impossible to provide a
        // Selection object.
        using MyBase::defaultValue;
        using MyBase::defaultValueIfSet;

        virtual AbstractOptionStoragePointer createStorage() const;

        SelectionFlags          selectionFlags_;

        /*! \brief
         * Needed to initialize SelectionOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class SelectionOptionStorage;
};

} // namespace gmx

#endif
