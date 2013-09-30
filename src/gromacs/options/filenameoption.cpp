/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include "gromacs/legacyheaders/filenm.h"

#include "gromacs/utility/file.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

class FileTypeRegistry;

/********************************************************************
 * FileTypeHandler
 */

/*! \internal \brief
 * Handles a single file type known to FileNameOptionStorage.
 *
 * \ingroup module_options
 */
class FileTypeHandler
{
    public:
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
        std::vector<const char *> extensions_;

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
 *
 * \ingroup module_options
 */
class FileTypeRegistry
{
    public:
        //! Returns a singleton instance of this class.
        static const FileTypeRegistry &instance();
        //! Returns a handler for a single file type.
        const FileTypeHandler &handlerForType(OptionFileType type) const;

    private:
        //! Initializes the file type registry.
        FileTypeRegistry();

        //! Registers a file type that corresponds to a ftp in filenm.h.
        void registerType(OptionFileType type, int ftp);
        //! Registers a file type with a single extension.
        void registerType(OptionFileType type, const char *extension);

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
FileTypeRegistry::handlerForType(OptionFileType type) const
{
    GMX_RELEASE_ASSERT(type >= 0 && static_cast<size_t>(type) < filetypes_.size(),
                       "Invalid file type");
    return filetypes_[type];
}

FileTypeRegistry::FileTypeRegistry()
{
    filetypes_.resize(eftOptionFileType_NR);
    registerType(eftTopology,    efTPS);
    registerType(eftTrajectory,  efTRX);
    registerType(eftPDB,         efPDB);
    registerType(eftIndex,       efNDX);
    registerType(eftPlot,        efXVG);
    registerType(eftGenericData, efDAT);
}

void FileTypeRegistry::registerType(OptionFileType type, int ftp)
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

void FileTypeRegistry::registerType(OptionFileType type,
                                    const char    *extension)
{
    GMX_RELEASE_ASSERT(type >= 0 && static_cast<size_t>(type) < filetypes_.size(),
                       "Invalid file type");
    filetypes_[type].extensions_.assign(1, extension);
}

/*! \brief
 * Helper method to complete a file name provided to a file name option.
 *
 * \param[in] value     Value provided to the file name option.
 * \param[in] filetype  File type for the option.
 * \param[in] bCompleteToExisting
 *      Whether to check existing files when completing the extension.
 * \returns   \p value with possible extension added.
 */
std::string completeFileName(const std::string &value, OptionFileType filetype,
                             bool bCompleteToExisting)
{
    if (bCompleteToExisting && File::exists(value))
    {
        // TODO: This may not work as expected if the value is passed to a
        // function that uses fn2ftp() to determine the file type and the input
        // file has an unrecognized extension.
        return value;
    }
    const FileTypeRegistry &registry    = FileTypeRegistry::instance();
    const FileTypeHandler  &typeHandler = registry.handlerForType(filetype);
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

}   // namespace

/********************************************************************
 * FileNameOptionStorage
 */

FileNameOptionStorage::FileNameOptionStorage(const FileNameOption &settings)
    : MyBase(settings), info_(this), filetype_(settings.filetype_),
      bRead_(settings.bRead_), bWrite_(settings.bWrite_),
      bLibrary_(settings.bLibrary_)
{
    if (settings.defaultBasename_ != NULL)
    {
        if (isRequired())
        {
            setDefaultValue(completeFileName(settings.defaultBasename_,
                                             filetype_, false));
        }
        else
        {
            setDefaultValueIfSet(completeFileName(settings.defaultBasename_,
                                                  filetype_, false));
        }
    }
}

std::string FileNameOptionStorage::formatSingleValue(const std::string &value) const
{
    return value;
}

void FileNameOptionStorage::convertValue(const std::string &value)
{
    bool bInput = isInputFile() || isInputOutputFile();
    addValue(completeFileName(value, filetype_, bInput));
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
