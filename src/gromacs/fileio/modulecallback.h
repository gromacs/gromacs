/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Interface for callback from output manager to decide if added modules
 * are feasable to be added to the output for a given filetype or not.
 *
 * \author
 * \inpublicapi
 * \ingroup fileio
 */
#ifndef GMX_FILEIO_MODULECALLBACK_H
#define GMX_FILEIO_MODULECALLBACK_H

#include <algorithm>
#include <memory>

#include "gromacs/fileio/trxio.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\brief
 * ModuleCallback class for checking that requested module registration in output manager
 * is valid.
 *
 * \inpublicapi
 * \ingroup fileio
 *
 */
class IModuleCallback
{
    public:
        //! Recognized flags.
        enum
        {
            /*! \brief
             * Base setting that says that the module has no requirements.
             *
             * Sets the flags to default setting to mkae sure all output methods
             * are supported.
             */
            efAnyOutputSupported = 1<<0,
            /*! \brief
             * Requires output method to support force output.
             *
             * If set, only output method supporting writing of forces to
             * the coordinate file will work, others will generate an invalid
             * input error.
             */
            efForceOutput     = 1<<1,
            /*! \brief
             * Requires output method to support velocity output.
             *
             * If set, only writing to files that support velocity output will succed.
             * Other writing methods will generate an error.
             *
             */
            efVelocityOutput        = 1<<2,
            /*! \brief
             * Requires output to support coordinate output.
             *
             * Default for most methods, will need to be able to write coordinates to
             * output file or generate an error.
             */
            efCoordinateOutput      = 1<<3,
            /*! \brief
             * Requires support for connection information in output format.
             *
             * If set, only file output that supports writing of connection information will succeed.
             * This means for now that only PDB and TNG files can be written. Other file writing
             * methods will fail.
             */
            efConnectionOutput    = 1<<4,
            /*! \brief
             * Requires that output format supports the writing of atom information to the file.
             *
             * If set, files will only be written if they can output the information from t_atoms
             * and otherwise create an error while writing.
             */
            efAtomOutput = 1<<5,
        };

        virtual ~IModuleCallback()
        {
        }
        //! Output format dependent flags to check if we can perform the module addition.
        unsigned long moduleFlags_ = 0;
        //! Return the flag status that decides if we can add the module or not.
        bool checkModuleFlag(unsigned long outputFlags) {return moduleFlags_ & outputFlags; }

    protected:
        /*! \brief
         * Default constructor for CoordinateOutput class.
         */
        IModuleCallback(unsigned long flag) : moduleFlags_(flag) {}
};

//! Smart pointer to manage the output manager object.
typedef std::shared_ptr<IModuleCallback>
    ModuleCallbackPointer;

} // namespace gmx

#endif
