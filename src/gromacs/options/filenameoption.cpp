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

#include <string>
#include <vector>

#include "gromacs/fileio/filenm.h"

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

class FileTypeRegistry;

//! \addtogroup module_options
//! \{

//! Shorthand for a list of file extensions.
typedef std::vector<const char *> ExtensionList;

/********************************************************************
 * FileTypeHandler
 */

/*! \internal \brief
 * Handles a single file type known to FileNameOptionStorage.
 */
class FileTypeHandler
{
    public:
        //! Returns the list of extensions for this file type.
        const ExtensionList &extensions() const { return extensions_; }

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
        //! Possible extensions for this file type.
        ExtensionList extensions_;

        /*! \brief
         * Needed for initialization; all initialization is handled by
         * FileTypeRegistry.
         */
        friend class FileTypeRegistry;
};

bool
FileTypeHandler::hasKnownExtension(const std::string &filename) const
{
    for (size_t i = 0; i < extensions_.size(); ++i)
    {
        if (endsWith(filename, extensions_[i]))
        {
            return true;
        }
    }
    return false;
}

std::string
FileTypeHandler::addExtension(const std::string &filename) const
{
    if (extensions_.empty())
    {
        return filename;
    }
    return filename + extensions_[0];
}

std::string
FileTypeHandler::findFileWithExtension(const std::string &filename) const
{
    for (size_t i = 0; i < extensions_.size(); ++i)
    {
        std::string testFilename(filename + extensions_[i]);
        if (File::exists(testFilename))
        {
            return testFilename;
        }
    }
    return std::string();
}

/********************************************************************
 * FileTypeRegistry
 */

/*! \internal \brief
 * Singleton for managing static file type info for FileNameOptionStorage.
 */
class FileTypeRegistry
{
    public:
        //! Returns a singleton instance of this class.
        static const FileTypeRegistry &instance();
        //! Returns a handler for a single file type.
        const FileTypeHandler &
        handlerForType(OptionFileType type, int legacyType) const;

    private:
        //! Initializes the file type registry.
        FileTypeRegistry();

        //! Registers a file type that corresponds to a ftp in filenm.h.
        void registerType(int type, int ftp);
        //! Registers a file type with a single extension.
        void registerType(int type, const char *extension);

        std::vector<FileTypeHandler> filetypes_;
};

// static
const FileTypeRegistry &
FileTypeRegistry::instance()
{
    static FileTypeRegistry singleton;
    return singleton;
}

const FileTypeHandler &
FileTypeRegistry::handlerForType(OptionFileType type, int legacyType) const
{
    int index = type;
    if (type == eftUnknown && legacyType >= 0)
    {
        index = eftOptionFileType_NR + legacyType;
    }
    GMX_RELEASE_ASSERT(index >= 0 && static_cast<size_t>(index) < filetypes_.size(),
                       "Invalid file type");
    return filetypes_[index];
}

FileTypeRegistry::FileTypeRegistry()
{
    filetypes_.resize(eftOptionFileType_NR + efNR);
    registerType(eftTopology,    efTPS);
    registerType(eftTrajectory,  efTRX);
    registerType(eftPDB,         efPDB);
    registerType(eftIndex,       efNDX);
    registerType(eftPlot,        efXVG);
    registerType(eftGenericData, efDAT);
    for (int i = 0; i < efNR; ++i)
    {
        registerType(eftOptionFileType_NR + i, i);
    }
}

void FileTypeRegistry::registerType(int type, int ftp)
{
    GMX_RELEASE_ASSERT(type >= 0 && static_cast<size_t>(type) < filetypes_.size(),
                       "Invalid file type");
    const int genericTypeCount = ftp2generic_count(ftp);
    if (genericTypeCount > 0)
    {
        const int *const genericTypes = ftp2generic_list(ftp);
        filetypes_[type].extensions_.clear();
        filetypes_[type].extensions_.reserve(genericTypeCount);
        for (int i = 0; i < genericTypeCount; ++i)
        {
            filetypes_[type].extensions_.push_back(ftp2ext_with_dot(genericTypes[i]));
        }
    }
    else
    {
        registerType(type, ftp2ext_with_dot(ftp));
    }
}

void FileTypeRegistry::registerType(int type, const char *extension)
{
    GMX_RELEASE_ASSERT(type >= 0 && static_cast<size_t>(type) < filetypes_.size(),
                       "Invalid file type");
    filetypes_[type].extensions_.assign(1, extension);
}

/*! \brief
 * Helper method to complete a file name provided to a file name option.
 *
 * \param[in] value      Value provided to the file name option.
 * \param[in] filetype   File type for the option.
 * \param[in] legacyType If \p filetype is eftUnknown, this gives the type as
 *     an enum value from filenm.h.
 * \param[in] bCompleteToExisting
 *     Whether to check existing files when completing the extension.
 * \returns   \p value with possible extension added.
 */
std::string completeFileName(const std::string &value, OptionFileType filetype,
                             int legacyType, bool bCompleteToExisting)
{
    if (bCompleteToExisting && File::exists(value))
    {
        // TODO: This may not work as expected if the value is passed to a
        // function that uses fn2ftp() to determine the file type and the input
        // file has an unrecognized extension.
        return value;
    }
    const FileTypeRegistry &registry    = FileTypeRegistry::instance();
    const FileTypeHandler  &typeHandler = registry.handlerForType(filetype, legacyType);
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

FileNameOptionStorage::FileNameOptionStorage(const FileNameOption &settings)
    : MyBase(settings), info_(this), filetype_(settings.filetype_),
      legacyType_(settings.legacyType_),
      bRead_(settings.bRead_), bWrite_(settings.bWrite_),
      bLibrary_(settings.bLibrary_)
{
    if (settings.defaultBasename_ != NULL)
    {
        std::string defaultValue =
            completeFileName(settings.defaultBasename_, filetype_,
                             legacyType_, false);
        setDefaultValueIfSet(defaultValue);
        if (isRequired())
        {
            setDefaultValue(defaultValue);
        }
    }
}

std::string FileNameOptionStorage::typeString() const
{
    const FileTypeRegistry       &registry    = FileTypeRegistry::instance();
    const FileTypeHandler        &typeHandler = registry.handlerForType(filetype_, legacyType_);
    const ExtensionList          &extensions  = typeHandler.extensions();
    std::string                   result;
    ExtensionList::const_iterator i;
    int                           count = 0;
    for (i = extensions.begin(); count < 2 && i != extensions.end(); ++i, ++count)
    {
        if (i != extensions.begin())
        {
            result.append("/");
        }
        result.append(*i);
    }
    if (i != extensions.end())
    {
        result.append("/...");
    }
    if (result.empty())
    {
        if (legacyType_ == efRND)
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
    const FileTypeRegistry       &registry    = FileTypeRegistry::instance();
    const FileTypeHandler        &typeHandler = registry.handlerForType(filetype_, legacyType_);
    const ExtensionList          &extensions  = typeHandler.extensions();
    std::string                   result;
    if (extensions.size() > 2)
    {
        result.append(":");
        ExtensionList::const_iterator i;
        for (i = extensions.begin(); i != extensions.end(); ++i)
        {
            result.append(" ");
            result.append((*i) + 1);
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
    bool bInput = isInputFile() || isInputOutputFile();
    addValue(completeFileName(value, filetype_, legacyType_, bInput));
}

bool FileNameOptionStorage::isDirectoryOption() const
{
    return legacyType_ == efRND;
}

ConstArrayRef<const char *> FileNameOptionStorage::extensions() const
{
    const FileTypeRegistry &registry    = FileTypeRegistry::instance();
    const FileTypeHandler  &typeHandler = registry.handlerForType(filetype_, legacyType_);
    const ExtensionList    &extensions  = typeHandler.extensions();
    return constArrayRefFromVector<const char *>(extensions.begin(), extensions.end());
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

FileNameOptionInfo::ExtensionList FileNameOptionInfo::extensions() const
{
    return option().extensions();
}

/********************************************************************
 * FileNameOption
 */

AbstractOptionStoragePointer FileNameOption::createStorage() const
{
    return AbstractOptionStoragePointer(new FileNameOptionStorage(*this));
}

} // namespace gmx
