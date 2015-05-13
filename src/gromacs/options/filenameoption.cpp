/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements classes in filenameoption.h and filenameoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "filenameoption.h"

#include <cstring>

#include <string>
#include <vector>

#include "gromacs/fileio/filenm.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "filenameoptionstorage.h"

namespace gmx
{

namespace
{

//! \addtogroup module_options
//! \{

/*! \brief
 * Mapping from OptionFileType to a file type in filenm.h.
 */
struct FileTypeMapping
{
    //! OptionFileType value to map.
    OptionFileType optionType;
    //! Corresponding file type from filenm.h.
    int            fileType;
};

//! Mappings from OptionFileType to file types in filenm.h.
const FileTypeMapping c_fileTypeMapping[] =
{
    { eftTopology,    efTPS },
    { eftTrajectory,  efTRX },
    { eftPDB,         efPDB },
    { eftIndex,       efNDX },
    { eftPlot,        efXVG },
    { eftGenericData, efDAT }
};

/********************************************************************
 * FileTypeHandler
 */

/*! \internal
 * \brief
 * Handles a single file type known to FileNameOptionStorage.
 *
 * Methods in this class do not throw, except for a possible std::bad_alloc
 * when constructing std::string return values.
 */
class FileTypeHandler
{
    public:
        /*! \brief
         * Returns a handler for a single file type.
         *
         * \param[in] fileType  File type (from filenm.h) to use.
         */
        explicit FileTypeHandler(int fileType);

        //! Returns the number of acceptable extensions for this file type.
        int extensionCount() const;
        //! Returns the extension with the given index.
        const char *extension(int i) const;

        //! Returns whether \p fileType (from filenm.h) is accepted for this type.
        bool isValidType(int fileType) const;

    private:
        /*! \brief
         * File type (from filenm.h) represented by this handler.
         *
         * -1 represents an unknown file type.
         */
        int        fileType_;
        //! Number of different extensions this type supports.
        int        extensionCount_;
        /*! \brief
         * List of simple file types that are included in this type.
         *
         * If `fileType_` represents a generic type in filenm.h, i.e., a type
         * that accepts multiple different types of files, then this is an
         * array of `extensionCount_` elements, each element specifying one
         * non-generic file type that this option accepts.
         * `NULL` for single-extension types.
         */
        const int *genericTypes_;
};

FileTypeHandler::FileTypeHandler(int fileType)
    : fileType_(fileType), extensionCount_(0), genericTypes_(NULL)
{
    if (fileType_ >= 0)
    {
        const int genericTypeCount = ftp2generic_count(fileType_);
        if (genericTypeCount > 0)
        {
            extensionCount_ = genericTypeCount;
            genericTypes_   = ftp2generic_list(fileType_);
        }
        else if (ftp2ext_with_dot(fileType_)[0] != '\0')
        {
            extensionCount_ = 1;
        }
    }
}

int FileTypeHandler::extensionCount() const
{
    return extensionCount_;
}

const char *FileTypeHandler::extension(int i) const
{
    GMX_ASSERT(i >= 0 && i < extensionCount_, "Invalid extension index");
    if (genericTypes_ != NULL)
    {
        return ftp2ext_with_dot(genericTypes_[i]);
    }
    return ftp2ext_with_dot(fileType_);
}

bool
FileTypeHandler::isValidType(int fileType) const
{
    if (genericTypes_ != NULL)
    {
        for (int i = 0; i < extensionCount(); ++i)
        {
            if (fileType == genericTypes_[i])
            {
                return true;
            }
        }
        return false;
    }
    else
    {
        return fileType == fileType_;
    }
}

//! \}

}   // namespace

/********************************************************************
 * FileNameOptionStorage
 */

FileNameOptionStorage::FileNameOptionStorage(const FileNameOption  &settings,
                                             FileNameOptionManager *manager)
    : MyBase(settings), info_(this), manager_(manager), fileType_(-1),
      defaultExtension_(""), bRead_(settings.bRead_), bWrite_(settings.bWrite_),
      bLibrary_(settings.bLibrary_), bAllowMissing_(settings.bAllowMissing_)
{
    GMX_RELEASE_ASSERT(!hasFlag(efOption_MultipleTimes),
                       "allowMultiple() is not supported for file name options");
    if (settings.optionType_ == eftUnknown && settings.legacyType_ >= 0)
    {
        fileType_ = settings.legacyType_;
    }
    else
    {
        ConstArrayRef<FileTypeMapping>                 map(c_fileTypeMapping);
        ConstArrayRef<FileTypeMapping>::const_iterator i;
        for (i = map.begin(); i != map.end(); ++i)
        {
            if (i->optionType == settings.optionType_)
            {
                fileType_ = i->fileType;
                break;
            }
        }
    }
    FileTypeHandler typeHandler(fileType_);
    if (settings.defaultType_ >= 0 && settings.defaultType_ < efNR)
    {
        // This also assures that the default type is not a generic type.
        GMX_RELEASE_ASSERT(typeHandler.isValidType(settings.defaultType_),
                           "Default type for a file option is not an accepted "
                           "type for the option");
        FileTypeHandler defaultHandler(settings.defaultType_);
        defaultExtension_ = defaultHandler.extension(0);
    }
    else if (typeHandler.extensionCount() > 0)
    {
        defaultExtension_ = typeHandler.extension(0);
    }
    if (settings.defaultBasename_ != NULL)
    {
        std::string defaultValue(settings.defaultBasename_);
        int         type = fn2ftp(settings.defaultBasename_);
        GMX_RELEASE_ASSERT(type == efNR || type == settings.defaultType_,
                           "Default basename has an extension that does not "
                           "match the default type");
        if (type == efNR)
        {
            defaultValue.append(defaultExtension());
        }
        setDefaultValueIfSet(defaultValue);
        if (isRequired() || settings.bLegacyOptionalBehavior_)
        {
            setDefaultValue(defaultValue);
        }
    }
}

std::string FileNameOptionStorage::typeString() const
{
    FileTypeHandler typeHandler(fileType_);
    std::string     result;
    int             count;
    for (count = 0; count < 2 && count < typeHandler.extensionCount(); ++count)
    {
        if (count > 0)
        {
            result.append("/");
        }
        result.append(typeHandler.extension(count));
    }
    if (count < typeHandler.extensionCount())
    {
        result.append("/...");
    }
    if (result.empty())
    {
        if (isDirectoryOption())
        {
            result = "dir";
        }
        else
        {
            result = "file";
        }
    }
    return result;
}

std::string FileNameOptionStorage::formatExtraDescription() const
{
    FileTypeHandler typeHandler(fileType_);
    std::string     result;
    if (typeHandler.extensionCount() > 2)
    {
        result.append(":");
        for (int i = 0; i < typeHandler.extensionCount(); ++i)
        {
            result.append(" [REF]");
            // Skip the dot.
            result.append(typeHandler.extension(i) + 1);
            result.append("[ref]");
        }
    }
    return result;
}

std::string FileNameOptionStorage::formatSingleValue(const std::string &value) const
{
    return value;
}

void FileNameOptionStorage::convertValue(const std::string &value)
{
    if (manager_ != NULL)
    {
        std::string processedValue = manager_->completeFileName(value, info_);
        if (!processedValue.empty())
        {
            // If the manager returns a value, use it without further checks,
            // except for sanity checking.
            if (!isDirectoryOption())
            {
                const int fileType = fn2ftp(processedValue.c_str());
                if (fileType == efNR)
                {
                    // If the manager returned an invalid file name, assume
                    // that it knows what it is doing.  But assert that it
                    // only does that for the only case that it is currently
                    // required for: VMD plugins.
                    GMX_ASSERT(isInputFile() && isTrajectoryOption(),
                               "Manager returned an invalid file name");
                }
                else
                {
                    GMX_ASSERT(isValidType(fileType),
                               "Manager returned an invalid file name");
                }
            }
            addValue(processedValue);
            return;
        }
    }
    // Currently, directory options are simple, and don't need any
    // special processing.
    // TODO: Consider splitting them into a separate DirectoryOption.
    if (isDirectoryOption())
    {
        addValue(value);
        return;
    }
    const int fileType = fn2ftp(value.c_str());
    if (fileType == efNR)
    {
        std::string message
            = formatString("File '%s' cannot be used by GROMACS because it "
                           "does not have a recognizable extension.\n"
                           "The following extensions are possible for this option:\n  %s",
                           value.c_str(), joinStrings(extensions(), ", ").c_str());
        GMX_THROW(InvalidInputError(message));
    }
    else if (!isValidType(fileType))
    {
        std::string message
            = formatString("File name '%s' cannot be used for this option.\n"
                           "Only the following extensions are possible:\n  %s",
                           value.c_str(), joinStrings(extensions(), ", ").c_str());
        GMX_THROW(InvalidInputError(message));
    }
    addValue(value);
}

void FileNameOptionStorage::processAll()
{
    if (manager_ != NULL && hasFlag(efOption_HasDefaultValue))
    {
        ValueList &valueList = values();
        GMX_RELEASE_ASSERT(valueList.size() == 1,
                           "There should be only one default value");
        if (!valueList[0].empty())
        {
            const std::string &oldValue = valueList[0];
            GMX_ASSERT(endsWith(oldValue, defaultExtension()),
                       "Default value does not have the expected extension");
            const std::string  prefix
                = stripSuffixIfPresent(oldValue, defaultExtension());
            const std::string  newValue
                = manager_->completeDefaultFileName(prefix, info_);
            if (!newValue.empty() && newValue != oldValue)
            {
                GMX_ASSERT(isValidType(fn2ftp(newValue.c_str())),
                           "Manager returned an invalid default value");
                valueList[0] = newValue;
                refreshValues();
            }
        }
    }
}

bool FileNameOptionStorage::isDirectoryOption() const
{
    return fileType_ == efRND;
}

bool FileNameOptionStorage::isTrajectoryOption() const
{
    return fileType_ == efTRX;
}

const char *FileNameOptionStorage::defaultExtension() const
{
    return defaultExtension_;
}

std::vector<const char *> FileNameOptionStorage::extensions() const
{
    FileTypeHandler           typeHandler(fileType_);
    std::vector<const char *> result;
    result.reserve(typeHandler.extensionCount());
    for (int i = 0; i < typeHandler.extensionCount(); ++i)
    {
        result.push_back(typeHandler.extension(i));
    }
    return result;
}

bool FileNameOptionStorage::isValidType(int fileType) const
{
    FileTypeHandler typeHandler(fileType_);
    return typeHandler.isValidType(fileType);
}

ConstArrayRef<int> FileNameOptionStorage::fileTypes() const
{
    if (fileType_ < 0)
    {
        return ConstArrayRef<int>();
    }
    const int genericTypeCount = ftp2generic_count(fileType_);
    if (genericTypeCount > 0)
    {
        return constArrayRefFromArray<int>(ftp2generic_list(fileType_), genericTypeCount);
    }
    return constArrayRefFromArray<int>(&fileType_, 1);
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

bool FileNameOptionInfo::allowMissing() const
{
    return option().allowMissing();
}

bool FileNameOptionInfo::isDirectoryOption() const
{
    return option().isDirectoryOption();
}

bool FileNameOptionInfo::isTrajectoryOption() const
{
    return option().isTrajectoryOption();
}

const char *FileNameOptionInfo::defaultExtension() const
{
    return option().defaultExtension();
}

FileNameOptionInfo::ExtensionList FileNameOptionInfo::extensions() const
{
    return option().extensions();
}

bool FileNameOptionInfo::isValidType(int fileType) const
{
    return option().isValidType(fileType);
}

ConstArrayRef<int> FileNameOptionInfo::fileTypes() const
{
    return option().fileTypes();
}

/********************************************************************
 * FileNameOption
 */

AbstractOptionStorage *
FileNameOption::createStorage(const OptionManagerContainer &managers) const
{
    return new FileNameOptionStorage(*this, managers.get<FileNameOptionManager>());
}

} // namespace gmx
