/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Declares gmx::SelectionOption and gmx::SelectionOptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONOPTION_H
#define GMX_SELECTION_SELECTIONOPTION_H

#include "gromacs/options/abstractoption.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionenums.h"

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
 * To use options of this type, SelectionOptionManager must first be added to
 * the Options collection.  For trajectory analysis tools, the framework takes
 * care of this.
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
        explicit SelectionOption(const char *name)
            : MyBase(name), defaultText_(""),
              selectionFlags_(efSelection_DisallowEmpty)
        {
        }

        /*! \brief
         * Request velocity evaluation for output positions.
         *
         * \see Selection::setEvaluateVelocities()
         */
        MyClass &evaluateVelocities()
        { selectionFlags_.set(efSelection_EvaluateVelocities); return me(); }
        /*! \brief
         * Request force evaluation for output positions.
         *
         * \see Selection::setEvaluateForces()
         */
        MyClass &evaluateForces()
        { selectionFlags_.set(efSelection_EvaluateForces); return me(); }
        /*! \brief
         * Only accept selections that evaluate to atom positions.
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
         * Allow specifying an unconditionally empty selection for this option.
         *
         * If this option is not set, selections that are unconditionally empty
         * (i.e., can never match any atoms) result in errors.
         * Note that even without this option, it is still possible that a
         * dynamic selection evaluates to zero atoms for some frames.
         */
        MyClass &allowEmpty()
        { selectionFlags_.clear(efSelection_DisallowEmpty); return me(); }

        /*! \brief
         * Sets default selection text for the option.
         *
         * If the option is not set by the user, the provided text is parsed as
         * the value of the selection.
         */
        MyClass &defaultSelectionText(const char *text)
        { defaultText_ = text; return me(); }

    private:
        // Disable possibility to allow multiple occurrences, since it isn't
        // implemented.
        using MyBase::allowMultiple;
        // Disable default value because it is impossible to provide a
        // Selection object.
        using MyBase::defaultValue;
        using MyBase::defaultValueIfSet;

        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

        const char             *defaultText_;
        SelectionFlags          selectionFlags_;

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
         * \see Selection::setEvaluateVelocities()
         */
        void setEvaluateVelocities(bool bEnabled);
        /*! \brief
         * Sets whether this option evaluates forces for positions.
         *
         * \param[in] bEnabled  If true, forces are evaluated.
         *
         * Does not throw.
         *
         * \see Selection::setEvaluateForces()
         */
        void setEvaluateForces(bool bEnabled);
        /*! \brief
         * Sets whether this option accepts positions that come from multiple
         * atoms.
         *
         * \param[in] bEnabled  If true, the option accepts only positions that
         *      evaluate to atom positions.
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

    private:
        SelectionOptionStorage &option();
        const SelectionOptionStorage &option() const;
};

} // namespace gmx

#endif
