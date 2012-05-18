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
 * Declares gmx::SelectionFileOptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONFILEOPTIONINFO_H
#define GMX_SELECTION_SELECTIONFILEOPTIONINFO_H

#include "../options/optioninfo.h"

namespace gmx
{

class SelectionFileOptionStorage;
class SelectionOptionManager;

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

        //! \copydoc SelectionOptionInfo::setManager()
        void setManager(SelectionOptionManager *manager);

    private:
        SelectionFileOptionStorage &option();
        const SelectionFileOptionStorage &option() const;
};

} // namespace gmx

#endif
