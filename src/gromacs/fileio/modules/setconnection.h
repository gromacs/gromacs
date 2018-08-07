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
 * Set frame connection information from user input.
 *
 * \author
 * \inpublicapi
 * \ingroup fileio
 */
#ifndef GMX_FILEIO_SETCONNECTION_H
#define GMX_FILEIO_SETCONNECTION_H

#include <algorithm>

#include "gromacs/fileio/coordinateoutput.h"
#include "gromacs/fileio/modulecallback.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\brief
 * SetConnection class allows changing writing of connection information to file.
 *
 * This class allows the user to define if connection information should be written
 * to the output coordinate file, and requests them from the
 * currently processed data.
 *
 * \inpublicapi
 * \ingroup fileio
 *
 */
class SetConnection : public ICoordinateOutput
{
    public:
        /*! \brief
         * Default constructor for SetConnection should not be used.
         *
         * Class should only be initialized with at least the base selection.
         */
        SetConnection() = delete;
        /*! \brief
         * Construct SetConnection object with choice for boolean value.
         *
         * Can be used to initialize SetConnection from outside of trajectoryanalysis
         * with the user specified option to write coordinate velocities or not.
         * framework.
         */
        explicit SetConnection(gmx_conect connections) : ICoordinateOutput(efConnectionOutput),
                                                         connections_(connections)
        {
        }
        /*! \brief
         * Copy constructor.
         */
        SetConnection(const SetConnection &old) = delete;
        /*! \brief
         * Assignment operator.
         */
        SetConnection &operator=(const SetConnection &old) = delete;
        /*! \brief
         * Move constructor for SetConnection.
         */
        SetConnection &operator=(SetConnection &&old)
        {
            connections_ = std::move(old.connections_);
            return *this;
        }
        /*! \brief
         *  Move constructor for SetConnection.
         */
        SetConnection(SetConnection &&old) : ICoordinateOutput(old.moduleFlags_),
                                             connections_(std::move(old.connections_))
        {
        }

        ~SetConnection() {}

        /*! \brief
         * Change coordinate frame information for output.
         */
        virtual void processFrame(const int /*framenumner*/, t_trxframe * /*input*/);

        gmx_conect getConnection() {return connections_; };

    private:
        //! Connection information that should be added to the output.
        gmx_conect connections_ = nullptr;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<SetConnection>
    SetConnectionPointer;

} // namespace gmx

#endif
