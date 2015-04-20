/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::SelectionFileOption and gmx::SelectionFileOptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONFILEOPTION_H
#define GMX_SELECTION_SELECTIONFILEOPTION_H

#include "gromacs/options/abstractoption.h"

namespace gmx
{

class SelectionFileOptionInfo;
class SelectionFileOptionStorage;
class SelectionOptionManager;

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
        //! OptionInfo subclass corresponding to this option type.
        typedef SelectionFileOptionInfo InfoType;

        //! Initializes an option with the given name.
        explicit SelectionFileOption(const char *name);

    private:
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;
};

/*! \libinternal \brief
 * Wrapper class for accessing and modifying selection file option information.
 *
 * \inlibraryapi
 * \ingroup module_selection
 */
class SelectionFileOptionInfo : public OptionInfo
{
    public:
        /*! \brief
         * Creates option info object for given storage object.
         *
         * Does not throw.
         */
        explicit SelectionFileOptionInfo(SelectionFileOptionStorage *option);
};

} // namespace gmx

#endif
