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
 * Declares gmx::SelectionOption and gmx::SelectionOptionInfo.
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

namespace gmx
{

class SelectionOptionInfo;
class SelectionOptionManager;
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

        SelectionFlags selectionFlags_;

        /*! \brief
         * Needed to initialize SelectionOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class SelectionOptionStorage;
};

/*! \brief
 * Wrapper class for accessing and modifying selection option information.
 *
 * Allows changes to a selection option after creation.
 *
 * This class provides the necessary interface for changing, e.g., the number
 * of allowed selections for a selection option after the option has been
 * created with Options::addOption().  This is needed if the number or other
 * flags are only known after other options have been parsed.  The main
 * advantage of this class over custom checks is that if used before
 * interactive selection prompt, the interactive prompt is updated accordingly.
 *
 * When using this class, the option should be initially created with the most
 * permissive flags, and this class should be used to place restrictions where
 * appropriate.  Otherwise, values that are provided before adjustments will
 * need to follow the more strict checks.  In most cases in trajectory analysis
 * (which is the main use case for selection options), the adjustments should
 * be done in TrajectoryAnalysisModule::initOptionsDone() for them to take
 * place before interactive selection prompts.
 *
 * An instance of this class for a selection option can be obtained with
 * SelectionOption::getAdjuster() when the option is created.
 *
 * Example use:
 * \code
   SelectionList sel;
   Options options("example", "Example options");
   SelectionOptionInfo *info;
   info = options.addOption(SelectionOption("sel").storeVector(&sel)
                             .multiValue());
   // < ... assign values to options ...>
   if ( condition )
   {
    // Put limitations on the selections based on the condition,
    // which can depend on other option values.
    // Throws if input given so far violates the limitations.
    info->setValueCount(2);
    info->setOnlyStatic(true);
   }
 * \endcode
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionOptionInfo : public OptionInfo
{
    public:
        /*! \brief
         * Creates option info object for given storage object.
         *
         * Does not throw.
         */
        explicit SelectionOptionInfo(SelectionOptionStorage *option);

        /*! \brief
         * Set manager for handling interaction with other options and the
         * selection collection.
         *
         * \param   manager  Selection manager to set.
         *
         * This must be called before the values are added.
         *
         * Typically it is called through setManagerForSelectionOptions(),
         * which recursively sets the manager for all selection options in
         * an Options object.
         *
         * Does not throw.
         */
        void setManager(SelectionOptionManager *manager);

        /*! \brief
         * Sets the number of selections allowed for the option.
         *
         * \param[in] count  Number of allowed selections.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if values have already been provided
         *      and their count does not match.
         */
        void setValueCount(int count);

        /*! \brief
         * Sets whether this option evaluates velocities for positions.
         *
         * \param[in] bEnabled  If true, velocities are evaluated.
         *
         * Does not throw.
         *
         * \see SelectionOption::evaluateVelocities()
         */
        void setEvaluateVelocities(bool bEnabled);
        /*! \brief
         * Sets whether this option evaluates forces for positions.
         *
         * \param[in] bEnabled  If true, forces are evaluated.
         *
         * Does not throw.
         *
         * \see SelectionOption::evaluateForces()
         */
        void setEvaluateForces(bool bEnabled);
        /*! \brief
         * Sets whether this option accepts positions that come from multiple
         * atoms.
         *
         * \param[in] bEnabled  If true, the option accepts only positions that
         *      evaluate to atom positions.
         *
         * TODO: This is not yet implemented.
         *
         * \see SelectionOption::onlyAtoms()
         */
        void setOnlyAtoms(bool bEnabled);
        /*! \brief
         * Sets whether this option accepts dynamic selections.
         *
         * \param[in] bEnabled  If true, the option accepts only static
         *      selections.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if dynamic selections have already been
         *      provided.
         *
         * Strong exception safety guarantee.
         *
         * \see SelectionOption::onlyStatic()
         */
        void setOnlyStatic(bool bEnabled);
        /*! \brief
         * Sets whether this option uses position masks for dynamic selections.
         *
         * \param[in] bEnabled  If true, the position masks are used.
         *
         * Does not throw.
         *
         * \see SelectionOption::dynamicMask()
         */
        void setDynamicMask(bool bEnabled);
        /*! \brief
         * Sets whether atom coordinates are allowed as reference positions.
         *
         * \param[in] bEnabled  If true, the option does not accept atom
         *      coordinates as reference positions.
         *
         * TODO: This is not yet implemented.
         *
         * \see SelectionOption::dynamicOnlyWhole()
         */
        void setDynamicOnlyWhole(bool bEnabled);

    private:
        SelectionOptionStorage &option();
        const SelectionOptionStorage &option() const;
};

} // namespace gmx

#endif
