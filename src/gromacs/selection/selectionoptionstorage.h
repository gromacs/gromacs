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
/*! \internal \file
 * \brief
 * Declares gmx::SelectionOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONOPTIONSTORAGE_H
#define GMX_SELECTION_SELECTIONOPTIONSTORAGE_H

#include <string>

#include "gromacs/options/optionstoragetemplate.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selectionoption.h"

namespace gmx
{

class SelectionOption;
class SelectionOptionManager;

/*! \internal \brief
 * Converts, validates, and stores selection values.
 *
 * \see SelectionOptionManager
 *
 * \ingroup module_selection
 */
class SelectionOptionStorage : public OptionStorageTemplate<Selection>
{
    public:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings   Storage settings.
         * \param     manager    Manager for this object.
         */
        SelectionOptionStorage(const SelectionOption  &settings,
                               SelectionOptionManager *manager);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const { return "selection"; }
        virtual std::string formatSingleValue(const Selection &value) const;

        /*! \brief
         * Adds selections to the storage.
         *
         * \param[in] selections  List of selections to add.
         * \param[in] bFullValue  If true, the provided selections are the full
         *      value of the option, and additional checks are performed.
         * \throws  std::bad_alloc if out of memory.
         * \throws  InvalidInputError if
         *      - There is an incorrect number of selections in \p selections.
         *      - Any selection in \p selections is not allowed for this
         *        option.
         *
         * This function is used to add selections from SelectionOptionManager.
         * It is called with \p bFullValue set to false from
         * SelectionOptionManager::convertOptionValue(), and \p bFullValue set
         * to true when parsing requested selections.
         */
        void addSelections(const SelectionList &selections,
                           bool                 bFullValue);

        // Required to access the number of values in selection requests.
        // See SelectionCollection::Impl.
        using MyBase::maxValueCount;
        //! \copydoc SelectionOptionInfo::setValueCount()
        void setAllowedValueCount(int count);
        /*! \brief
         * Alters flags for the selections created by this option.
         *
         * \param[in] flag        Flag to change.
         * \param[in] bSet        Whether to set or clear the flag.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if selections have already been
         *      provided and conflict with the given flags.
         *
         * If selections have already been provided, it is checked that they
         * match the limitations enforced by the flags.  Pending requests are
         * also affected.
         *
         * Strong exception safety guarantee.
         */
        void setSelectionFlag(SelectionFlag flag, bool bSet);

    private:
        virtual void convertValue(const std::string &value);
        virtual void processSetValues(ValueList *values);
        virtual void processAll();

        SelectionOptionInfo     info_;
        SelectionOptionManager &manager_;
        std::string             defaultText_;
        SelectionFlags          selectionFlags_;
};

} // namespace gmx

#endif
