/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief This file declares a class for file handling, so that such
 * dependency may be injected into other code, and e.g. replaced when
 * testing.
 *
 * \todo Expand functionality here to absorb and replace old C-style
 * wrappers
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_FILEIOHANDLER_H
#define GMX_FILEIO_FILEIOHANDLER_H

#include "gromacs/fileio/fileiohandlerinterface.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

/*! \libinternal
 * \brief Wrapper used for dependency injection for objects that
 * need to get input from or about files.
 */
class FileIOHandler : public FileIOHandlerInterface
{
    public:
        //! Virtual destructor required
        virtual ~FileIOHandler();
        //! \copydoc gmx::FileIOHandlerInterface::gmx_fexist
        virtual gmx_bool gmx_fexist(const char *fname) const;
        //! \copydoc gmx::FileIOHandlerInterface::read_checkpoint_simulation_part_and_filenames
        virtual void
        read_checkpoint_simulation_part_and_filenames(const char           *filename,
                                                      int                  *simulation_part,
                                                      int                  *nfiles,
                                                      gmx_file_position_t **outputfiles) const;
};

//! Convenience typedef
typedef gmx_unique_ptr<FileIOHandler>::type FileIOHandlerPointer;

} // namespace

#endif
