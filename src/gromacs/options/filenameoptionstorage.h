/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Declares gmx::FileNameOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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
        virtual std::string typeString() const;
        virtual std::string formatExtraDescription() const;
        virtual std::string formatSingleValue(const std::string &value) const;

        //! \copydoc FileNameOptionInfo::isInputFile()
        bool isInputFile() const { return bRead_ && !bWrite_; }
        //! \copydoc FileNameOptionInfo::isOutputFile()
        bool isOutputFile() const { return !bRead_ && bWrite_; }
        //! \copydoc FileNameOptionInfo::isInputOutputFile()
        bool isInputOutputFile() const { return bRead_ && bWrite_; }
        //! \copydoc FileNameOptionInfo::isLibraryFile()
        bool isLibraryFile() const { return bLibrary_; }

        //! \copydoc FileNameOptionInfo::isDirectoryOption()
        bool isDirectoryOption() const;
        //! \copydoc FileNameOptionInfo::extensions()
        ConstArrayRef<const char *> extensions() const;

    private:
        virtual void convertValue(const std::string &value);

        FileNameOptionInfo      info_;
        OptionFileType          filetype_;
        int                     legacyType_;
        bool                    bRead_;
        bool                    bWrite_;
        bool                    bLibrary_;
};

} // namespace gmx

#endif
