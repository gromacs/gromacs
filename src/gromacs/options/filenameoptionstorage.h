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
 * Declares gmx::FileNameOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_FILENAMEOPTIONSTORAGE_H
#define GMX_OPTIONS_FILENAMEOPTIONSTORAGE_H

#include <string>

#include "filenameoption.h"
#include "optionfiletype.h"
#include "optionstoragetemplate.h"

namespace gmx
{

class FileNameOption;

/*! \internal \brief
 * Converts, validates, and stores file names.
 */
class FileNameOptionStorage : public OptionStorageTemplate<std::string>
{
    public:
        //! \copydoc StringOptionStorage::StringOptionStorage()
        explicit FileNameOptionStorage(const FileNameOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual const char *typeString() const { return "file"; }
        virtual std::string formatSingleValue(const std::string &value) const;

        //! \copydoc FileNameOptionInfo::isInputFile()
        bool isInputFile() const { return bRead_ && !bWrite_; }
        //! \copydoc FileNameOptionInfo::isOutputFile()
        bool isOutputFile() const { return !bRead_ && bWrite_; }
        //! \copydoc FileNameOptionInfo::isInputOutputFile()
        bool isInputOutputFile() const { return bRead_ && bWrite_; }
        //! \copydoc FileNameOptionInfo::isLibraryFile()
        bool isLibraryFile() const { return bLibrary_; }

    private:
        virtual void convertValue(const std::string &value);

        FileNameOptionInfo      info_;
        OptionFileType          filetype_;
        bool                    bRead_;
        bool                    bWrite_;
        bool                    bLibrary_;
};

} // namespace gmx

#endif
