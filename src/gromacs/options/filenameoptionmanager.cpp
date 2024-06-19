/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::FileNameOptionManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/filenameoptionmanager.h"

#include <cstring>

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fileredirector.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! Extensions that are recognized as compressed files.
const char* const c_compressedExtensions[] = { ".gz", ".Z" };

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
std::string findExistingExtension(const std::string&          prefix,
                                  const FileNameOptionInfo&   option,
                                  const IFileInputRedirector* redirector)
{
    ArrayRef<const int>                 types = option.fileTypes();
    ArrayRef<const int>::const_iterator i;
    for (i = types.begin(); i != types.end(); ++i)
    {
        std::string testFilename(prefix + ftp2ext_with_dot(*i));
        if (redirector->fileExists(testFilename, File::throwOnError))
        {
            return testFilename;
        }
    }
    return std::string();
}

} // namespace

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
    Impl() : redirector_(&defaultFileInputRedirector()), bInputCheckingDisabled_(false) {}

    //! Redirector for file existence checks.
    const IFileInputRedirector* redirector_;
    //! Global default file name, if set.
    std::string defaultFileName_;
    //! Whether input option processing has been disabled.
    bool bInputCheckingDisabled_;
};

/********************************************************************
 * FileNameOptionManager
 */

FileNameOptionManager::FileNameOptionManager() : impl_(new Impl()) {}

FileNameOptionManager::~FileNameOptionManager() {}

void FileNameOptionManager::setInputRedirector(const IFileInputRedirector* redirector)
{
    impl_->redirector_ = redirector;
}

void FileNameOptionManager::disableInputOptionChecking(bool bDisable)
{
    impl_->bInputCheckingDisabled_ = bDisable;
}

void FileNameOptionManager::addDefaultFileNameOption(IOptionsContainer* options, const char* name)
{
    options->addOption(
            StringOption(name).store(&impl_->defaultFileName_).description("Set the default filename for all file options"));
}

std::string FileNameOptionManager::completeFileName(const std::string& value, const FileNameOptionInfo& option)
{
    const bool bAllowMissing = option.allowMissing();
    const bool bInput        = option.isInputFile() || option.isInputOutputFile();
    // Currently, directory options are simple, and don't need any
    // special processing.
    // TODO: Consider splitting them into a separate DirectoryOption.
    if (option.isDirectoryOption())
    {
        if (!impl_->bInputCheckingDisabled_ && bInput && !bAllowMissing
            && !std::filesystem::is_directory(value))
        {
            std::string message =
                    formatString("Directory '%s' does not exist or is not accessible.", value.c_str());
            // TODO: Get actual errno value from the attempt to open the file
            // to provide better feedback to the user.
            GMX_THROW(InvalidInputError(message));
        }
        return value;
    }
    const int fileType = fn2ftp(value.c_str());
    if (bInput && !impl_->bInputCheckingDisabled_)
    {
        if (fileType == efNR && impl_->redirector_->fileExists(value, File::throwOnError))
        {
            ArrayRef<const char* const> compressedExtensions(c_compressedExtensions);
            ArrayRef<const char* const>::const_iterator ext;
            for (ext = compressedExtensions.begin(); ext != compressedExtensions.end(); ++ext)
            {
                if (endsWith(value, *ext))
                {
                    std::string newValue = value.substr(0, value.length() - std::strlen(*ext));
                    if (option.isValidType(fn2ftp(newValue.c_str())))
                    {
                        // Note that here the "completed" filename no
                        // longer has the extension appropriate for a
                        // compressed file. When gmx_ffopen() sees
                        // that the uncompressed file does not exist,
                        // it will try to open the compressed
                        // versions, one of which exists, and thus
                        // will be able to fulfil the user's desire.
                        // This is not very robust; doing better might
                        // require a more elaborate abstraction around
                        // the concept of a filename returned by
                        // FileNameOptionManager.
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
            std::string processedValue = findExistingExtension(value, option, impl_->redirector_);
            if (!processedValue.empty())
            {
                return processedValue;
            }
            if (bAllowMissing)
            { // NOLINT bugprone-branch-clone
                return value + option.defaultExtension();
            }
            else if (option.isLibraryFile())
            {
                // TODO: Treat also library files here and remove the NOLINT.
                return value + option.defaultExtension();
            }
            else
            {
                std::string message = formatString(
                        "File '%s' does not exist or is not accessible.\n"
                        "The following extensions were tried to complete the file name:\n  %s",
                        value.c_str(),
                        joinStrings(option.extensions(), ", ").c_str());
                GMX_THROW(InvalidInputError(message));
            }
        }
        else if (option.isValidType(fileType))
        {
            if (option.isLibraryFile())
            {
                // TODO: Treat also library files.
            }
            else if (!bAllowMissing)
            {
                if (!impl_->redirector_->fileExists(value, File::throwOnNotFound))
                {
                    return std::string();
                }
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

std::string FileNameOptionManager::completeDefaultFileName(const std::string&        prefix,
                                                           const FileNameOptionInfo& option)
{
    if (option.isDirectoryOption())
    {
        return std::string();
    }
    const bool        bInput     = option.isInputFile() || option.isInputOutputFile();
    const std::string realPrefix = !impl_->defaultFileName_.empty() ? impl_->defaultFileName_ : prefix;
    if (bInput && !impl_->bInputCheckingDisabled_)
    {
        std::string completedName = findExistingExtension(realPrefix, option, impl_->redirector_);
        if (!completedName.empty())
        {
            return completedName;
        }
        if (option.allowMissing())
        { // NOLINT bugprone-branch-clone
            return realPrefix + option.defaultExtension();
        }
        else if (option.isLibraryFile()) // NOLINT bugprone-branch-clone
        {
            // TODO: Treat also library files here and remove the NOLINT
            return realPrefix + option.defaultExtension();
        }
        else if (option.isSet())
        {
            std::string message = formatString(
                    "No file name was provided, and the default file "
                    "'%s' does not exist or is not accessible.\n"
                    "The following extensions were tried to complete the file name:\n  %s",
                    prefix.c_str(),
                    joinStrings(option.extensions(), ", ").c_str());
            GMX_THROW(InvalidInputError(message));
        }
        else if (option.isRequired())
        {
            std::string message = formatString(
                    "Required option was not provided, and the default file "
                    "'%s' does not exist or is not accessible.\n"
                    "The following extensions were tried to complete the file name:\n  %s",
                    prefix.c_str(),
                    joinStrings(option.extensions(), ", ").c_str());
            GMX_THROW(InvalidInputError(message));
        }
        // We get here with the legacy optional behavior.
    }
    return realPrefix + option.defaultExtension();
}

} // namespace gmx
