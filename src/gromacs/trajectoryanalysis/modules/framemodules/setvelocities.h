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
/*! \libinternal \file
 * \brief
 * Handler for writing coordinate data to files.
 * Takes care of all your coordinate file writing issues.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_SETVELOCITY_H
#define GMX_TRAJECTORYANALYSIS_MODULES_SETVELOCITY_H

#include <algorithm>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "framemanager.h"

namespace gmx
{

/*!\libinternal
 * \brief
 * SetVelocities class for handling trajectory file opening and data writing
 *
 * The filehandler keeps both track of the file being currently written too
 * and the correct number of coordinates that should be written, as well as
 * the output file type.
 * All the information is set by the user through the options mechanics
 * in the framework.
 *
 * \ingroup module_trajectoryanalysis
 *
 */
class SetVelocities : public IFrameManager
{
    public:
        /*! \brief
         * Default constructor for SetVelocities should not be used.
         *
         * Class should only be initialized with at least the base selection.
         */
        SetVelocities() = delete;
        /*! \brief
         * Construct SetVelocities object with choice for boolean value
         * for frame velocity writing.
         *
         * Can be used to initialize SetVelocities from outside of trajectoryanalysis
         * framework.
         */
        explicit SetVelocities(bool velocity) : velocity_(velocity)
        {
        }
        /*! \brief
         * Copy constructor.
         */
        SetVelocities(const SetVelocities &old) = delete;
        /*! \brief
         * Assignment operator.
         */
        SetVelocities &operator=(const SetVelocities &old) = delete;
        /*! \brief
         * Move constructor for SetVelocities.
         */
        SetVelocities &operator=(SetVelocities &&old)
        {
            velocity_ = std::move(old.velocity_);
            return *this;
        }
        /*! \brief
         *  Move constructor for SetVelocities.
         */
        SetVelocities(SetVelocities &&old) : velocity_(std::move(old.velocity_))
        {
        }

        ~SetVelocities() {}
        /*! \brief
         * Pass any user input options to the frame manager.
         *
         * Currently not used, will be useful to pass user input information to frame manager.
         */
        virtual void initFileOptions(IOptionsContainer * /*options*/);

        /*! \brief
         * Change coordinate frame information for output.
         *
         * Takes the previously internally stored coordinates and saves them again.
         * Applies correct number of atoms, as well as changing things such as
         * frame time or affect output of velocities or forces.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void modifyFrame(const t_trxframe &input);
        /*! \brief
         * Sanity check for user input options.
         *
         * This function performs the check of the user input for basic sanity issues
         * and should be called after option processing has been finished.
         */
        virtual void checkOptions();

    private:
        /*! \brief
         * Selection of atoms that will be written to disk.
         *
         * Internal selection of atoms chosen by the user that will be written
         * to disk during processing. All actions that the filehandler performs
         * will only be on those atoms, with the remaining ones being not affected.
         */
        bool                            velocity_;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<SetVelocities>
    SetVelocitiesPointer;

} // namespace gmx

#endif
