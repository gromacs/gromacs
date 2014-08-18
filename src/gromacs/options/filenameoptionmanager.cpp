/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Implements gmx::FileNameOptionManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "filenameoptionmanager.h"

#include <cstring>

#include <string>

#include "gromacs/fileio/filenm.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! Extensions that are recognized as compressed files.
const char *const c_compressedExtensions[] =
{ ".gz", ".Z" };

/********************************************************************
 * Helper functions
 */

/*! \brief
 * Adds an extension to \p prefix if it results in an existing file.
 *
 * Tries to add each extension for this file type to \p prefix and
 * checks whether this results in an existing file.
 * The first match is returned.
 * Returns an empty string if no existing file is found.
 */
std::string findExistingExtension(const std::string        &prefix,
                                  const FileNameOptionInfo &option)
{
    ConstArrayRef<int>                 types = option.fileTypes();
    ConstArrayRef<int>::const_iterator i;
    for (i = types.begin(); i != types.end(); ++i)
    {
        std::string testFilename(prefix + ftp2ext_with_dot(*i));
        if (File::exists(testFilename))
        {
            return testFilename;
        }
    }
    return std::string();
}

}   // namespace

/********************************************************************
 * FileNameOptionManager::Impl
 */

/*! \internal \brief
 * Private implemention class for FileNameOptionManager.
 *
 * \ingroup module_options
 */
class FileNameOptionManager::Impl
{
    public:
        Impl() : bInputCheckingDisabled_(false) {}

        //! Global default file name, if set.
        std::string     defaultFileName_;
        //! Whether input option processing has been disabled.
        bool            bInputCheckingDisabled_;
};

/********************************************************************
 * FileNameOptionManager
 */

FileNameOptionManager::FileNameOptionManager()
    : impl_(new Impl())
{
}

FileNameOptionManager::~FileNameOptionManager()
{
}

void FileNameOptionManager::disableInputOptionChecking(bool bDisable)
{
    impl_->bInputCheckingDisabled_ = bDisable;
}

void FileNameOptionManager::addDefaultFileNameOption(
        Options *options, const char *name)
{
    options->addOption(
            StringOption(name).store(&impl_->defaultFileName_)
                .description("Set the default filename for all file options"));
}

std::string FileNameOptionManager::completeFileName(
        const std::string &value, const FileNameOptionInfo &option)
{
    const bool bInput = option.isInputFile() || option.isInputOutputFile();
    // Currently, directory options are simple, and don't need any
    // special processing.
    // TODO: Consider splitting them into a separate DirectoryOption.
    if (option.isDirectoryOption())
    {
        if (!impl_->bInputCheckingDisabled_ && bInput && !Directory::exists(value))
        {
            std::string message
                = formatString("Directory '%s' does not exist or is not accessible.",
                               value.c_str());
            // TODO: Get actual errno value from the attempt to open the file
            // to provide better feedback to the user.
            GMX_THROW(InvalidInputError(message));
        }
        return value;
    }
    const int fileType = fn2ftp(value.c_str());
    if (bInput && !impl_->bInputCheckingDisabled_)
    {
        if (fileType == efNR && File::exists(value))
        {
            ConstArrayRef<const char *>                 compressedExtensions(c_compressedExtensions);
            ConstArrayRef<const char *>::const_iterator ext;
            for (ext = compressedExtensions.begin(); ext != compressedExtensions.end(); ++ext)
            {
                if (endsWith(value, *ext))
                {
                    std::string newValue = value.substr(0, value.length() - std::strlen(*ext));
                    if (option.isValidType(fn2ftp(newValue.c_str())))
                    {
                        return newValue;
                    }
                    else
                    {
                        return std::string();
                    }
                }
            }
            // VMD plugins may be able to read the file.
            if (option.isInputFile() && option.isTrajectoryOption())
            {
                return value;
            }
        }
        else if (fileType == efNR)
        {
            std::string processedValue = findExistingExtension(value, option);
            if (!processedValue.empty())
            {
                return processedValue;
            }
            if (option.isLibraryFile())
            {
                // TODO: Treat also library files here.
                return value + option.defaultExtension();
            }
            else
            {
                std::string message
                    = formatString("File '%s' does not exist or is not accessible.\n"
                                   "The following extensions were tried to complete the file name:\n  %s",
                                   value.c_str(), joinStrings(option.extensions(), ", ").c_str());
                GMX_THROW(InvalidInputError(message));
            }
        }
        else if (option.isValidType(fileType))
        {
            if (option.isLibraryFile())
            {
                // TODO: Treat also library files.
            }
            else if (!File::exists(value))
            {
                std::string message
                    = formatString("File '%s' does not exist or is not accessible.",
                                   value.c_str());
                // TODO: Get actual errno value from the attempt to open the file
                // to provide better feedback to the user.
                GMX_THROW(InvalidInputError(message));
            }
            return value;
        }
    }
    else // Not an input file
    {
        if (fileType == efNR)
        {
            return value + option.defaultExtension();
        }
        else if (option.isValidType(fileType))
        {
            return value;
        }
    }
    return std::string();
}

std::string FileNameOptionManager::completeDefaultFileName(
        const std::string &prefix, const FileNameOptionInfo &option)
{
    if (option.isDirectoryOption() || impl_->bInputCheckingDisabled_)
    {
        return std::string();
    }
    const bool        bInput = option.isInputFile() || option.isInputOutputFile();
    const std::string realPrefix
        = !impl_->defaultFileName_.empty() ? impl_->defaultFileName_ : prefix;
    if (bInput)
    {
        std::string completedName = findExistingExtension(realPrefix, option);
        if (!completedName.empty())
        {
            return completedName;
        }
        if (option.isLibraryFile())
        {
            // TODO: Treat also library files here.
            return realPrefix + option.defaultExtension();
        }
        else if (option.isSet())
        {
            std::string message
                = formatString("No file name was provided, and the default file "
                               "'%s' does not exist or is not accessible.\n"
                               "The following extensions were tried to complete the file name:\n  %s",
                               prefix.c_str(), joinStrings(option.extensions(), ", ").c_str());
            GMX_THROW(InvalidInputError(message));
        }
        else if (option.isRequired())
        {
            std::string message
                = formatString("Required option was not provided, and the default file "
                               "'%s' does not exist or is not accessible.\n"
                               "The following extensions were tried to complete the file name:\n  %s",
                               prefix.c_str(), joinStrings(option.extensions(), ", ").c_str());
            GMX_THROW(InvalidInputError(message));
        }
        // We get here with the legacy optional behavior.
    }
    return realPrefix + option.defaultExtension();
}

} // namespace gmx
