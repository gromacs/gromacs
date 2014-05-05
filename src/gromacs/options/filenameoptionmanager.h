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

#include "../utility/common.h"

namespace gmx
{

class Options;

/*! \brief
 * Handles interaction of file name options with global options.
 *
 * Currently, this class implements support for a global default file name
 * that overrides any option-specific default.
 *
 * \todo
 * Currently, this class has very little logic, and just provides the global
 * values to FileNameOptionStorage implementation.  A cleaner design would have
 * most of the non-trivial file name completion logic in this class, so that
 * the customizations would be centralized here.
 *
 * Creating a FileNameOptionManager for an Options object is optional, even if
 * the Options contains FileNameOption options.  Features from the manager are
 * not available if the manager is not created, but otherwise the options work.
 *
 * \see setManagerForFileNameOptions()
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class FileNameOptionManager
{
    public:
        FileNameOptionManager();
        ~FileNameOptionManager();

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

        //! Returns the currently set default file name.
        const std::string &defaultFileName() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \brief
 * Set manager for all file name options.
 *
 * Recursively sets the manager to \p manager for all file name options in
 * \p options.
 * Must be called before value assignment starts for \p options.
 *
 * Does not throw.
 *
 * \inpublicapi
 */
void setManagerForFileNameOptions(Options               *options,
                                  FileNameOptionManager *manager);

} // namespace gmx

#endif
