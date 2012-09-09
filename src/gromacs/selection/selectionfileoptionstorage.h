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
/*! \internal \file
 * \brief
 * Declares gmx::SelectionFileOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONFILEOPTIONSTORAGE_H
#define GMX_SELECTION_SELECTIONFILEOPTIONSTORAGE_H

#include "../options/abstractoptionstorage.h"
#include "selectionfileoption.h"

namespace gmx
{

class SelectionFileOption;
class SelectionOptionManager;

/*! \internal \brief
 * Implementation for a special option for reading selections from files.
 *
 * \ingroup module_selection
 */
class SelectionFileOptionStorage : public AbstractOptionStorage
{
    public:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings   Storage settings.
         */
        SelectionFileOptionStorage(const SelectionFileOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual const char *typeString() const { return "file"; }
        virtual int valueCount() const { return 0; }
        virtual std::string formatValue(int /*i*/) const { return ""; }

        //! \copydoc SelectionFileOptionInfo::setManager()
        void setManager(SelectionOptionManager *manager)
        {
            manager_ = manager;
        }

    private:
        virtual void clearSet();
        virtual void convertValue(const std::string &value);
        virtual void processSet();
        virtual void processAll() {}

        SelectionFileOptionInfo info_;
        SelectionOptionManager *manager_;
        bool                    bValueParsed_;
};

} // namespace gmx

#endif
