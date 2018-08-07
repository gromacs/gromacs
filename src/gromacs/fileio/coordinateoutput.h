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
 * Interface for opening and modifying files for coordinate data output.
 *
 * \author
 * \inpublicapi
 * \ingroup fileio
 */
#ifndef GMX_FILEIO_COORDINATEOUTPUT_H
#define GMX_FILEIO_COORDINATEOUTPUT_H

#include <algorithm>

#include "gromacs/fileio/trxio.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\brief
 * CoordinateOutput class for handling trajectory file opening and initialization.
 *
 * \inpublicapi
 * \ingroup module_coordinatedata
 *
 */
class ICoordinateOutput
{
    public:
        /*! \brief
         * Default constructor for CoordinateOutput class.
         */
        ICoordinateOutput()
        {
            clear_trxframe(&coordinates_, true);
        }
        virtual ~ICoordinateOutput()
        {
        }
        //! Move constructor for old object.
        explicit ICoordinateOutput(ICoordinateOutput &&old)
            : coordinates_(std::move(old.coordinates_))
        {
        }

        /*! \brief
         * Change settings in t_trxframe according to user input.
         */
        virtual void processFrame(const int framenumber, const t_trxframe &input) = 0;

        //! Test if the atoms data is available for writing
        bool haveAtoms() const { return coordinates_.bAtoms; };
        //! Test if the topology data is available for writing
        bool haveMtop() const { return mtop_ ? true : false; };

        //! Return local coordinates in derived modules.
        t_trxframe getFrame() const { return coordinates_; };

        //! Private copy of coordinate frame data used in output modules.
        t_trxframe coordinates_;

        //! Pointer to topology information if available.
        gmx_mtop_t *mtop_;
};

//! Smart pointer to manage the output manager object.
typedef std::shared_ptr<ICoordinateOutput>
    CoordinateOutputPointer;

} // namespace gmx

#endif
