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
/*! \file
 * \brief
 * Declares gmx::FileNameOption and gmx::FileNameOptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_FILENAMEOPTION_H
#define GMX_OPTIONS_FILENAMEOPTION_H

#include <string>

#include "abstractoption.h"
#include "optionfiletype.h"

namespace gmx
{

class FileNameOptionInfo;
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
            : MyBase(name), filetype_(eftUnknown), defaultBasename_(NULL),
            bRead_(false), bWrite_(false), bLibrary_(false)
        {
        }

        /*! \brief
         * Sets the type of the file this option accepts.
         *
         * This attribute must be provided.
         */
        MyClass &filetype(OptionFileType type)
        { filetype_ = type; return me(); }
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
         * Tells that the file will be looked up in library directories in
         * addition to working directory.
         *
         * \todo
         * Currently, this flag only affects the help output.  Callers must
         * take care themselves to actually search the file in the library
         * directories.  It would be nicer to do this searching within the
         * file name option implementation.
         */
        MyClass &libraryFile() { bLibrary_ = true; return me(); }
        /*! \brief
         * Sets a default basename for the file option.
         *
         * Use this method instead of defaultValue() or defaultValueIfSet() to
         * set a default value for a file name option.  No extension needs to
         * be provided; it is automatically added based on filetype().
         * The behavior is also adjusted based on required(): if the option is
         * required, the value given to defaultBasename() is treated as for
         * defaultValue(), otherwise it is treated as for defaultValueIfSet().
         */
        MyClass &defaultBasename(const char *basename)
        { defaultBasename_ = basename; return me(); }

    private:
        // Use defaultBasename() instead.
        using MyBase::defaultValue;
        using MyBase::defaultValueIfSet;

        //! Creates a FileNameOptionStorage object.
        virtual AbstractOptionStoragePointer createStorage() const;

        OptionFileType filetype_;
        const char             *defaultBasename_;
        bool bRead_;
        bool bWrite_;
        bool bLibrary_;

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

    private:
        const FileNameOptionStorage &option() const;
};

} // namespace gmx

#endif
