/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::FileNameOption and gmx::FileNameOptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_FILENAMEOPTION_H
#define GMX_OPTIONS_FILENAMEOPTION_H

#include <string>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/optionfiletype.h"

namespace gmx
{

template <typename T> class ArrayRef;
class FileNameOptionInfo;
class FileNameOptionManager;
class FileNameOptionStorage;

/*! \brief
 * Specifies an option that provides file names.
 *
 * Public methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class FileNameOption : public OptionTemplate<std::string, FileNameOption>
{
    public:
        //! OptionInfo subclass corresponding to this option type.
        typedef FileNameOptionInfo InfoType;

        //! Initializes an option with the given name.
        explicit FileNameOption(const char *name)
            : MyBase(name), optionType_(eftUnknown), legacyType_(-1),
              defaultBasename_(nullptr), defaultType_(-1),
              bLegacyOptionalBehavior_(false),
              bRead_(false), bWrite_(false), bLibrary_(false),
              bAllowMissing_(false)
        {
        }

        /*! \brief
         * Sets the type of the file this option accepts.
         *
         * Either this attribute or legacyType() must be provided.
         */
        MyClass &filetype(OptionFileType type)
        { optionType_ = type; return me(); }
        /*! \brief
         * Sets the type of the file from an enum in filetypes.h.
         *
         * New code should prefer filetype(), extending the enumeration if
         * necessary.
         */
        MyClass &legacyType(int type)
        { legacyType_ = type; return me(); }
        /*! \brief
         * Changes the behavior of optional options to match old t_filenm.
         *
         * If this is not set, optional options return an empty string if not
         * set.  If this is set, a non-empty value is always returned.
         * In the latter case, whether the option is set only affects the
         * return value of OptionInfo::isSet() and Options::isSet().
         */
        MyClass &legacyOptionalBehavior()
        { bLegacyOptionalBehavior_ = true; return me(); }
        //! Tells that the file provided by this option is used for input only.
        MyClass &inputFile()
        { bRead_ = true; bWrite_ = false; return me(); }
        //! Tells that the file provided by this option is used for output only.
        MyClass &outputFile()
        { bRead_ = false; bWrite_ = true; return me(); }
        /*! \brief
         * Tells that the file provided by this option is used for input and
         * output both.
         */
        MyClass &inputOutputFile()
        { bRead_ = bWrite_ = true; return me(); }
        /*! \brief
         * Sets the read/write usage for this file from boolean flags.
         */
        MyClass &readWriteFlags(bool bRead, bool bWrite)
        { bRead_ = bRead; bWrite_ = bWrite; return me(); }
        /*! \brief
         * Tells that the file will be looked up in library directories in
         * addition to working directory.
         *
         * \todo
         * Currently, this flag only affects the help output.  Callers must
         * take care themselves to actually search the file in the library
         * directories.  It would be nicer to do this searching within the
         * file name option implementation.
         */
        MyClass &libraryFile(bool bLibrary = true)
        { bLibrary_ = bLibrary; return me(); }
        /*! \brief
         * Tells that missing file names explicitly provided by the user are
         * valid for this input option.
         *
         * If this method is not called, an error will be raised if the user
         * explicitly provides a file name that does not name an existing file,
         * or if the default value does not resolve to a valid file name for a
         * required option that the user has not set.
         *
         * This method only has effect with input files, and only if a
         * FileNameOptionManager is being used.
         */
        MyClass &allowMissing(bool bAllow = true)
        { bAllowMissing_ = bAllow; return me(); }
        /*! \brief
         * Sets a default basename for the file option.
         *
         * Use this method instead of defaultValue() or defaultValueIfSet() to
         * set a default value for a file name option.  No extension needs to
         * be provided; it is automatically added based on filetype() or
         * defaultType().
         * The behavior is also adjusted based on required(): if the option is
         * required, the value given to defaultBasename() is treated as for
         * both defaultValue() and defaultValueIfSet(), otherwise it is treated
         * as for defaultValueIfSet().
         *
         * For input files that accept multiple extensions, the extension is
         * completed to the default extension on creation of the option or at
         * time of parsing an option without a value.
         *
         * If FileNameOptionManager is used, the extension may change during
         * Options::finish(), as this is the time when the default names are
         * checked against the file system to provide an extension that matches
         * an existing file if that is possible.
         *
         * If FileNameOptionManager is used, and
         * FileNameOptionManager::addDefaultFileNameOption() is used, and the
         * user provides a global default file name using that option, then the
         * global default takes precedence over defaultBasename().
         */
        MyClass &defaultBasename(const char *basename)
        { defaultBasename_ = basename; return me(); }
        /*! \brief
         * Sets a default type/extension for the file option.
         *
         * For options that accept multiple types of files (e.g.,
         * eftTrajectory), this method sets the default extension used
         * for completing defaultBasename(), as well as the default extension
         * used by FileNameOptionManager to complete various file names.
         *
         * The value should be one of the enumerated `ef*` values from
         * filetypes.h, and be a valid type for the type specified with
         * filetype().
         */
        MyClass &defaultType(int filetype)
        { defaultType_ = filetype; return me(); }

    private:
        // Use defaultBasename() instead.
        using MyBase::defaultValue;
        using MyBase::defaultValueIfSet;

        //! Creates a FileNameOptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

        OptionFileType          optionType_;
        int                     legacyType_;
        const char             *defaultBasename_;
        int                     defaultType_;
        bool                    bLegacyOptionalBehavior_;
        bool                    bRead_;
        bool                    bWrite_;
        bool                    bLibrary_;
        bool                    bAllowMissing_;

        /*! \brief
         * Needed to initialize FileNameOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class FileNameOptionStorage;
};

/*! \brief
 * Wrapper class for accessing file name option information.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class FileNameOptionInfo : public OptionInfo
{
    public:
        //! Shorthand for a list of extensions.
        typedef std::vector<const char *> ExtensionList;

        //! Creates an option info object for the given option.
        explicit FileNameOptionInfo(FileNameOptionStorage *option);

        //! Whether the option specifies an input file.
        bool isInputFile() const;
        //! Whether the option specifies an output file.
        bool isOutputFile() const;
        //! Whether the option specifies a file used for both input and output.
        bool isInputOutputFile() const;
        /*! \brief
         * Whether the option specifies a library file.
         *
         * \see FileNameOption::libraryFile()
         */
        bool isLibraryFile() const;
        //! Whether the (input) option allows missing files to be provided.
        bool allowMissing() const;

        //! Whether the option specifies directories.
        bool isDirectoryOption() const;
        //! Whether the option specifies a generic trajectory file.
        bool isTrajectoryOption() const;
        //! Returns the default extension for this option.
        const char *defaultExtension() const;
        //! Returns the list of extensions this option accepts.
        ExtensionList extensions() const;
        //! Returns whether \p fileType (from filetypes.h) is accepted for this option.
        bool isValidType(int fileType) const;
        //! Returns the list of file types this option accepts.
        ArrayRef<const int> fileTypes() const;

    private:
        const FileNameOptionStorage &option() const;
};

} // namespace gmx

#endif
