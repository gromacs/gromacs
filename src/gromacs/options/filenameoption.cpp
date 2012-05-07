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
 * Implements classes in filenameoption.h, filenameoptioninfo.h and
 * filenameoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/filenameoption.h"

#include <string>

#include "gromacs/options/filenameoptioninfo.h"

#include "filenameoptionstorage.h"

namespace gmx
{

/********************************************************************
 * FileNameOptionStorage
 */

FileNameOptionStorage::FileNameOptionStorage(const FileNameOption &settings)
    : MyBase(settings), info_(this), filetype_(settings.filetype_),
      bRead_(settings.bRead_), bWrite_(settings.bWrite_),
      bLibrary_(settings.bLibrary_)
{
}

std::string FileNameOptionStorage::formatSingleValue(const std::string &value) const
{
    return value;
}

void FileNameOptionStorage::convertValue(const std::string &value)
{
    // TODO: Proper implementation.
    addValue(value);
}

/********************************************************************
 * FileNameOptionInfo
 */

FileNameOptionInfo::FileNameOptionInfo(FileNameOptionStorage *option)
    : OptionInfo(option)
{
}

const FileNameOptionStorage &FileNameOptionInfo::option() const
{
    return static_cast<const FileNameOptionStorage &>(OptionInfo::option());
}

bool FileNameOptionInfo::isInputFile() const
{
    return option().isInputFile();
}

bool FileNameOptionInfo::isOutputFile() const
{
    return option().isOutputFile();
}

bool FileNameOptionInfo::isInputOutputFile() const
{
    return option().isInputOutputFile();
}

bool FileNameOptionInfo::isLibraryFile() const
{
    return option().isLibraryFile();
}

/********************************************************************
 * FileNameOption
 */

AbstractOptionStoragePointer FileNameOption::createStorage() const
{
    return AbstractOptionStoragePointer(new FileNameOptionStorage(*this));
}

} // namespace gmx
