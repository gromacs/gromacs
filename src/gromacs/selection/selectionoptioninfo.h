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
 * Declares gmx::SelectionOptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONOPTIONINFO_H
#define GMX_SELECTION_SELECTIONOPTIONINFO_H

#include "../options/optioninfo.h"

namespace gmx
{

class Options;
class SelectionCollection;
class SelectionOptionStorage;

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
options.addOption(SelectionOption("sel").storeVector(&sel)
                      .multiValue().getAdjuster(&info));
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
        //! Creates option info object for given storage object.
        explicit SelectionOptionInfo(SelectionOptionStorage *option);

        /*! \brief
         * Set selection collection into which this option adds selections.
         *
         * This must be called before the values are added.
         *
         * Typically it is called through setSelectionCollectionForOptions(),
         * which recursively sets the collection for all selection options in
         * an Options object.
         */
        void setSelectionCollection(SelectionCollection *selections);

        /*! \brief
         * Sets the number of selections allowed for the option.
         *
         * \param[in] count  Number of allowed selections.
         */
        void setValueCount(int count);

        //! \copydoc SelectionOption::evaluateVelocities()
        void setEvaluateVelocities(bool bEnabled);
        //! \copydoc SelectionOption::evaluateForces()
        void setEvaluateForces(bool bEnabled);
        //! \copydoc SelectionOption::onlyAtoms()
        void setOnlyAtoms(bool bEnabled);
        //! \copydoc SelectionOption::onlyStatic()
        void setOnlyStatic(bool bEnabled);
        //! \copydoc SelectionOption::dynamicMask()
        void setDynamicMask(bool bEnabled);
        //! \copydoc SelectionOption::dynamicOnlyWhole()
        void setDynamicOnlyWhole(bool bEnabled);

    private:
        SelectionOptionStorage &option();
        const SelectionOptionStorage &option() const;
};

/*! \cond libapi */
/*! \libinternal \brief
 * Set selection collection for all selection options.
 *
 * Recursively sets the selection collection to \p selections for all selection
 * options in \p options.
 * Must be called before value assignment starts for \p options.
 *
 * Does not throw.
 *
 * \inlibraryapi
 */
void setSelectionCollectionForOptions(Options *options,
                                      SelectionCollection *selections);
//! \endcond

} // namespace gmx

#endif
