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
#include "gromacs/options/filenameoptionmanager.h"

#include <string>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"

namespace gmx
{

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
        //! Global default file name, if set.
        std::string     defaultFileName_;
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

void FileNameOptionManager::addDefaultFileNameOption(
        Options *options, const char *name)
{
    options->addOption(
            StringOption(name).store(&impl_->defaultFileName_)
                .description("Set the default filename for all file options"));
}

const std::string &FileNameOptionManager::defaultFileName() const
{
    return impl_->defaultFileName_;
}

/********************************************************************
 * Global functions
 */

namespace
{

/*! \internal \brief
 * Visitor that sets the manager for each file name option.
 *
 * \ingroup module_options
 */
class FileNameOptionManagerSetter : public OptionsModifyingVisitor
{
    public:
        //! Construct a visitor that sets given manager.
        explicit FileNameOptionManagerSetter(FileNameOptionManager *manager)
            : manager_(manager)
        {
        }

        void visitSubSection(Options *section)
        {
            OptionsModifyingIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }

        void visitOption(OptionInfo *option)
        {
            FileNameOptionInfo *fileOption
                = option->toType<FileNameOptionInfo>();
            if (fileOption != NULL)
            {
                fileOption->setManager(manager_);
            }
        }

    private:
        FileNameOptionManager *manager_;
};

}   // namespace

void setManagerForFileNameOptions(Options               *options,
                                  FileNameOptionManager *manager)
{
    FileNameOptionManagerSetter(manager).visitSubSection(options);
}

} // namespace gmx
