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

#include <cassert>

#include "../options/abstractoption.h"
#include "selectionenums.h"

namespace gmx
{

class AbstractErrorReporter;
class Selection;
class SelectionOptionAdjuster;
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
        explicit SelectionOption(const char *name)
            : MyBase(name), _adjuster(NULL)
        { }

        /*! \brief
         * Request velocity evaluation for output positions.
         */
        MyClass &evaluateVelocities()
        { _selectionFlags.set(efEvaluateVelocities); return me(); }
        /*! \brief
         * Request force evaluation for output positions.
         */
        MyClass &evaluateForces()
        { _selectionFlags.set(efEvaluateForces); return me(); }
        /*! \brief
         * Only accept selections that evaluate to atom positions.
         *
         * TODO: This option is not yet implemented.
         */
        MyClass &onlyAtoms()
        { _selectionFlags.set(efOnlyAtoms); return me(); }
        /*! \brief
         * Only accept static selections for this option.
         */
        MyClass &onlyStatic()
        { _selectionFlags.set(efOnlyStatic); return me(); }
        /*! \brief
         * Handle dynamic selections for this option with position masks.
         *
         * Sets ::POS_MASKONLY on the positions for this selection.
         */
        MyClass &dynamicMask()
        { _selectionFlags.set(efDynamicMask); return me(); }
        /*! \brief
         * Disallow using atom coordinates as the reference positions.
         *
         * TODO: This option is not yet implemented.
         */
        MyClass &dynamicOnlyWhole()
        { _selectionFlags.set(efDynamicOnlyWhole); return me(); }

        /*! \brief
         * Get an adjuster that can be used to alter the option after creation.
         *
         * \see SelectionOptionAdjuster
         */
        MyClass &getAdjuster(SelectionOptionAdjuster **adjusterp)
        { _adjuster = adjusterp; return me(); }

    private:
        // Disable default value because it is impossible to provide a
        // Selection object.
        using MyBase::defaultValue;
        using MyBase::defaultValueIfSet;

        virtual int createDefaultStorage(Options *options,
                                         AbstractOptionStorage **storage) const;

        SelectionFlags          _selectionFlags;
        SelectionOptionAdjuster **_adjuster;

        /*! \brief
         * Needed to initialize SelectionOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class SelectionOptionStorage;
};


/*! \brief
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
std::vector<Selection *> sel;
Options options("example", "Example options");
SelectionOptionAdjuster *adjuster;
options.addOption(SelectionOption("sel").storeVector(&sel)
                      .multiValue().getAdjuster(&adjuster));
// < ... assign values to options ...>
if ( condition )
{
    OptionAdjusterErrorContext context(adjuster, errors);
    // Put limitations on the selections based on the condition,
    // which can depend on other option values.
    int rc = adjuster->setValueCount(2);
    if (rc == 0)
    {
        rc = adjuster->setOnlyStatic(true);
    }
    if (rc != 0)
    {
        return rc;
    }
}
 * \endcode
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionOptionAdjuster
{
    public:
        /*! \brief
         * Creates a new adjuster.
         *
         * Should only be called internally by the options module.
         */
        SelectionOptionAdjuster(SelectionOptionStorage *storage);

        /*! \brief
         * Sets an error reporter for all subsequent operations.
         *
         * \param[in] errors  Error reporter object.
         * \returns The previously set error reporter (may be NULL).
         *
         * Caller must ensure that the error reporter is valid as long as it
         * is assigned to the adjuster.
         *
         * There must be a valid error reporter associated with the adjuster
         * before any other method in the object is called.
         *
         * \see OptionAdjusterErrorContext
         */
        AbstractErrorReporter *setErrorReporter(AbstractErrorReporter *errors);

        /*! \brief
         * Sets the number of selections allowed for the option.
         *
         * \param[in] count  Number of allowed selections.
         * \retval 0 on success.
         */
        int setValueCount(int count);

        //! \copydoc SelectionOption::evaluateVelocities()
        int setEvaluateVelocities(bool bEnabled);
        //! \copydoc SelectionOption::evaluateForces()
        int setEvaluateForces(bool bEnabled);
        //! \copydoc SelectionOption::onlyAtoms()
        int setOnlyAtoms(bool bEnabled);
        //! \copydoc SelectionOption::onlyStatic()
        int setOnlyStatic(bool bEnabled);
        //! \copydoc SelectionOption::dynamicMask()
        int setDynamicMask(bool bEnabled);
        //! \copydoc SelectionOption::dynamicOnlyWhole()
        int setDynamicOnlyWhole(bool bEnabled);

    private:
        //! Returns the storage object associated with this adjuster.
        SelectionOptionStorage &storage() { return _storage; }
        //! Returns the storage object associated with this adjuster.
        const SelectionOptionStorage &storage() const { return _storage; }
        //! Returns the current error reporter object, asserts if there is none.
        AbstractErrorReporter *errors()
        { assert(_errors != NULL); return _errors; }

        SelectionOptionStorage &_storage;
        AbstractErrorReporter  *_errors;
};


/*! \brief
 * Convenience class for providing an error reporter for an option adjuster.
 *
 * This class implements a RAII-type interface over
 * SelectionOptionAdjuster::setErrorReporter(): the constructor sets a new
 * error reporter, and the destructor restores the old reporter.
 *
 * Example use:
 * \code
{
    OptionAdjusterErrorContext context(adjuster, errors);
    adjuster->setValueCount(2); // Errors are reported using 'errors'.
}  // Previous reporter is automatically restored on scope exit.
 * \endcode
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class OptionAdjusterErrorContext
{
    public:
        //! Sets error reporter for an adjuster and stores the old reporter.
        OptionAdjusterErrorContext(SelectionOptionAdjuster *adjuster,
                                   AbstractErrorReporter *errors)
            : _adjuster(adjuster)
        {
            _oldReporter = adjuster->setErrorReporter(errors);
        }
        //! Restores the old error reporter.
        ~OptionAdjusterErrorContext()
        {
            _adjuster->setErrorReporter(_oldReporter);
        }

    private:
        SelectionOptionAdjuster *_adjuster;
        AbstractErrorReporter   *_oldReporter;
};

} // namespace gmx

#endif
