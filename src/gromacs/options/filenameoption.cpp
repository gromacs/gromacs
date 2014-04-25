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
 * Implements classes in filenameoption.h and filenameoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "filenameoption.h"
#include "filenameoptionstorage.h"

#include <cstring>

#include <string>
#include <vector>

#include "gromacs/fileio/filenm.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

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
    { eftEnergy,      efEDR },
    { eftPDB,         efPDB },
    { eftIndex,       efNDX },
    { eftPlot,        efXVG },
    { eftGenericData, efDAT }
};

//! Extensions that are recognized as compressed files.
const char *const c_compressedExtensions[] =
{ ".gz", ".Z" };

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

        //! Returns whether \p filename has a valid extension for this type.
        bool hasKnownExtension(const std::string &filename) const;
        //! Adds a default extension for this type to \p filename.
        std::string addExtension(const std::string &filename) const;
        /*! \brief
         * Adds an extension to \p filename if it results in an existing file.
         *
         * Tries to add each extension for this file type to \p filename and
         * checks whether this results in an existing file.
         * The first match is returned.
         * Returns an empty string if no existing file is found.
         */
        std::string findFileWithExtension(const std::string &filename) const;

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
FileTypeHandler::hasKnownExtension(const std::string &filename) const
{
    for (int i = 0; i < extensionCount(); ++i)
    {
        if (endsWith(filename, extension(i)))
        {
            return true;
        }
    }
    return false;
}

std::string
FileTypeHandler::addExtension(const std::string &filename) const
{
    if (extensionCount() == 0)
    {
        return filename;
    }
    return filename + extension(0);
}

std::string
FileTypeHandler::findFileWithExtension(const std::string &filename) const
{
    for (int i = 0; i < extensionCount(); ++i)
    {
        std::string testFilename(filename + extension(i));
        if (File::exists(testFilename))
        {
            return testFilename;
        }
    }
    return std::string();
}

/*! \brief
 * Helper method to complete a file name provided to a file name option.
 *
 * \param[in] value       Value provided to the file name option.
 * \param[in] typeHandler Handler for the file type.
 * \param[in] bCompleteToExisting
 *     Whether to check existing files when completing the extension.
 * \returns   \p value with possible extension added.
 */
std::string completeFileName(const std::string     &value,
                             const FileTypeHandler &typeHandler,
                             bool                   bCompleteToExisting)
{
    if (bCompleteToExisting && File::exists(value))
    {
        // TODO: This may not work as expected if the value is passed to a
        // function that uses fn2ftp() to determine the file type and the input
        // file has an unrecognized extension.
        ConstArrayRef<const char *>                 compressedExtensions(c_compressedExtensions);
        ConstArrayRef<const char *>::const_iterator ext;
        for (ext = compressedExtensions.begin(); ext != compressedExtensions.end(); ++ext)
        {
            if (endsWith(value, *ext))
            {
                return value.substr(0, value.length() - std::strlen(*ext));
            }
        }
        return value;
    }
    if (typeHandler.hasKnownExtension(value))
    {
        return value;
    }
    if (bCompleteToExisting)
    {
        std::string newValue = typeHandler.findFileWithExtension(value);
        if (!newValue.empty())
        {
            return newValue;
        }
    }
    return typeHandler.addExtension(value);
}

//! \}

}   // namespace

/********************************************************************
 * FileNameOptionStorage
 */

FileNameOptionStorage::FileNameOptionStorage(const FileNameOption  &settings,
                                             FileNameOptionManager *manager)
    : MyBase(settings), info_(this), manager_(manager), fileType_(-1),
      bRead_(settings.bRead_), bWrite_(settings.bWrite_),
      bLibrary_(settings.bLibrary_)
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
    if (settings.defaultBasename_ != NULL)
    {
        std::string defaultValue(settings.defaultBasename_);
        defaultValue.append(defaultExtension());
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
            result.append(" ");
            // Skip the dot.
            result.append(typeHandler.extension(i) + 1);
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
    const bool      bInput = isInputFile() || isInputOutputFile();
    FileTypeHandler typeHandler(fileType_);
    addValue(completeFileName(value, typeHandler, bInput));
}

void FileNameOptionStorage::processAll()
{
    FileTypeHandler typeHandler(fileType_);
    if (hasFlag(efOption_HasDefaultValue) && typeHandler.extensionCount() > 0)
    {
        const bool  bInput      = isInputFile() || isInputOutputFile();
        ValueList  &valueList   = values();
        GMX_RELEASE_ASSERT(valueList.size() == 1,
                           "There should be only one default value");
        const bool  bGlobalDefault =
            (manager_ != NULL && !manager_->defaultFileName().empty());
        if (!valueList[0].empty() && (typeHandler.extensionCount() > 1 || bGlobalDefault))
        {
            const std::string &oldValue = valueList[0];
            GMX_ASSERT(endsWith(oldValue, defaultExtension()),
                       "Default value does not have the expected extension");
            std::string prefix = stripSuffixIfPresent(oldValue, defaultExtension());
            if (bGlobalDefault)
            {
                prefix = manager_->defaultFileName();
            }
            std::string newValue = completeFileName(prefix, typeHandler, bInput);
            if (newValue != oldValue)
            {
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

const char *FileNameOptionStorage::defaultExtension() const
{
    FileTypeHandler typeHandler(fileType_);
    if (typeHandler.extensionCount() == 0)
    {
        return "";
    }
    return typeHandler.extension(0);
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

bool FileNameOptionInfo::isDirectoryOption() const
{
    return option().isDirectoryOption();
}

const char *FileNameOptionInfo::defaultExtension() const
{
    return option().defaultExtension();
}

FileNameOptionInfo::ExtensionList FileNameOptionInfo::extensions() const
{
    return option().extensions();
}

/********************************************************************
 * FileNameOption
 */

AbstractOptionStoragePointer
FileNameOption::createStorage(const OptionManagerContainer &managers) const
{
    return AbstractOptionStoragePointer(
            new FileNameOptionStorage(*this, managers.get<FileNameOptionManager>()));
}

} // namespace gmx
