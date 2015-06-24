/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Declares gmx::FileNameOptionManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_FILENAMEOPTIONMANAGER_H
#define GMX_OPTIONS_FILENAMEOPTIONMANAGER_H

#include <string>

#include "gromacs/options/options.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class FileInputRedirectorInterface;
class FileNameOptionInfo;
class Options;

/*! \brief
 * Handles interaction of file name options with global options.
 *
 * This class contains all logic that completes file names based on user input
 * and file system contents.  Additionally, this class implements support for a
 * global default file name that overrides any option-specific default, as well
 * as additional control over how the completion is done.
 *
 * \todo
 * Most of the functionality in this class is specific to command line parsing,
 * so it would be cleaner to replace this with an interface, and have the
 * actual code in the `commandline` module.
 *
 * Adding a FileNameOptionManager for an Options object is optional, even if
 * the Options contains FileNameOption options.  Features from the manager are
 * not available if the manager is not created, but otherwise the options work:
 * the values provided to FileNameOption are used as they are, and exceptions
 * are thrown if they are no valid instead of attempting to complete them.
 *
 * \see Options::addManager()
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class FileNameOptionManager : public OptionManagerInterface
{
    public:
        FileNameOptionManager();
        virtual ~FileNameOptionManager();

        /*! \brief
         * Redirects file existence checks.
         *
         * \param[in] redirector  File redirector to use for existence checks.
         *
         * The manager checks for existence of various files on the file system
         * to complete file extensions.  This method can be used to redirect
         * those checks to an alternative implementation.
         *
         * This is used for unit tests to more easily control the result of the
         * checks and to keep the tests as fast as possible by avoiding real
         * file system access.  To keep implementation options open, behavior
         * with `redirector == NULL` is undefined and should not be relied on.
         * For tests, there should only be need to call this a single time,
         * right after creating the manager.
         */
        void setInputRedirector(const FileInputRedirectorInterface *redirector);

        /*! \brief
         * Disables special input file option handling.
         *
         * If disabled, this removes all file system calls from the file
         * name option parsing.
         * The values returned by FileNameOption for input and input/output
         * files are handled with the same simple rule as for output files:
         * the default extension is added if the file does not end in a
         * recognized extension, and no other checking is done.
         *
         * This changes the following behavior:
         *  - Providing non-existent files does not trigger errors.
         *  - Extensions for input files are not completed to an existing file.
         *  - Compressed input files do not work.
         */
        void disableInputOptionChecking(bool bDisable);

        /*! \brief
         * Adds an option for setting the default global file name.
         *
         * \param     options Options to add the option to.
         * \param[in] name    Name of the option to add.
         *
         * If the user sets the option, it affects all file name options that
         * would normally return a default value: the basename for the returned
         * value is taken from the value of the default file name option,
         * instead from an option-specific default
         * (FileNameOption::defaultBaseName()).
         */
        void addDefaultFileNameOption(Options *options, const char *name);

        /*! \brief
         * Completes file name option values.
         *
         * \param[in] value  Value provided by the user.
         * \param[in] option Option for which the value should be completed.
         * \returns   Value for the file name option.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if the value is not valid for this
         *     option.
         *
         * This method is called for each value that the user provides to
         * a FileNameOption.  The return value (if non-empty) is used as the
         * value of the option instead of the user-provided one.
         */
        std::string completeFileName(const std::string        &value,
                                     const FileNameOptionInfo &option);
        /*! \brief
         * Completes default values for file name options.
         *
         * \param[in] prefix Default prefix for the file name.
         * \param[in] option Option for which the value should be completed.
         * \returns   Value for the file name option.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if the value is not valid for this
         *     option.
         *
         * This method is called for each FileNameOption that has a default
         * value (either a standard default value, or if the user provided the
         * option without an explicit value).  \p prefix is the default value
         * without the default extension for the option.
         * If the return value is non-empty, it is used as the default value
         * for the option instead of \p prefix + default extension.
         */
        std::string completeDefaultFileName(const std::string        &prefix,
                                            const FileNameOptionInfo &option);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
