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
 * \brief This file declares an interface for file handling, so that
 * such dependency may be injected into modules.
 *
 * It is intended that the gmx::FileIOHandlerInterface abstract base
 * class have multiple subclasses, e.g. gmx::FileIOHandler for use in
 * real GROMACS code, and other mock or fake subclasses used in
 * testing.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_FILEIOHANDLERINTERFACE_H
#define GMX_FILEIO_FILEIOHANDLERINTERFACE_H

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

/*! \libinternal
 * \brief Interface used for dependency injection of methods used to
 * get input from or about files. This permits mock subclasses to be
 * used during testing. */
class FileIOHandlerInterface
{
    public:
        //! Virtual destructor required
        virtual ~FileIOHandlerInterface() {};
        //! Returns whether the named file exists
        virtual gmx_bool gmx_fexist(const char *fname) const = 0;
        //! Reads and returns simulation metadata from a checkpoint file
        virtual void
        read_checkpoint_simulation_part_and_filenames(const char           *filename,
                                                      int                  *simulation_part,
                                                      int                  *nfiles,
                                                      gmx_file_position_t **outputfiles) const = 0;
};

} // namespace

#endif
