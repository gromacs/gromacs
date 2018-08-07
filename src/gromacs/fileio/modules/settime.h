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
 * Set flag for changing time setting in a coordinate frame.
 *
 * \author
 * \inpublicapi
 * \ingroup fileio
 */
#ifndef GMX_FILEIO_SETTIME_H
#define GMX_FILEIO_SETTIME_H

#include <algorithm>

#include "gromacs/fileio/coordinateoutput.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\brief
 * SetTime class allows changing writing of velocities to file.
 *
 * This class allows the user to define if velocities should be written
 * to the output coordinate file, and checks if they are available from the
 * currently processed data.
 *
 * \inpublicapi
 * \ingroup module_coordinatedata
 *
 */
class SetTime : public ICoordinateOutput
{
    public:
        /*! \brief
         * Default constructor for SetTime should not be used.
         *
         * Class should only be initialized with at least the base selection.
         */
        SetTime() = delete;
        /*! \brief
         * Construct SetTime object with choice for boolean value.
         *
         * Can be used to initialize SetTime from outside of trajectoryanalysis
         * with the user specified option to write coordinate velocities or not.
         * framework.
         */
        explicit SetTime(real startTime, real timeStep, bool bSetStartTime, bool bSetTimeStep) :
            startTime_(startTime),
            timeStep_(timeStep),
            bSetStartTime_(bSetStartTime),
            bSetTimeStep_(bSetTimeStep)
        {
        }
        /*! \brief
         * Copy constructor.
         */
        SetTime(const SetTime &old) = delete;
        /*! \brief
         * Assignment operator.
         */
        SetTime &operator=(const SetTime &old) = delete;
        /*! \brief
         * Move constructor for SetTime.
         */
        SetTime &operator=(SetTime &&old)
        {
            startTime_     = std::move(old.startTime_);
            timeStep_      = std::move(old.timeStep_);
            bSetStartTime_ = std::move(old.bSetStartTime_);
            bSetTimeStep_  = std::move(old.bSetTimeStep_);
            return *this;
        }
        /*! \brief
         *  Move constructor for SetTime.
         */
        SetTime(SetTime &&old) :
            startTime_(std::move(old.startTime_)),
            timeStep_(std::move(old.timeStep_)),
            bSetStartTime_(std::move(old.bSetStartTime_)),
            bSetTimeStep_(std::move(old.bSetTimeStep_))
        {
        }

        ~SetTime() {}

        /*! \brief
         * Change coordinate frame information for output.
         *
         * In this case, the correct flag for writing the velocities is applied
         * to the output frame, depending on user selection and availability
         * in the input data.
         *
         * \todo should this throw an error if velocity writing is desired
         * but not possible from input data?
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void processFrame(const int framenumber, const t_trxframe &input);

        /*! \brief
         * Initialize the timeshift variable.
         *
         * Set up the timeshift variable needed to compute the acutal current frame
         * time when it is supposed to get changed according to user input.
         * Only done once right at the first frame, never after.
         *
         * \param[in] time Input time from first frame.
         */
        void setTimeShift(real time);

    private:
        /*! \brief
         * Flag to specify if frame time will be changed.
         *
         * Internal storage for the user choice for new starting time.
         */
        real                            startTime_;
        //! New time difference between frames.
        real                            timeStep_;
        //! Has the start time been changed?
        bool                            bSetStartTime_;
        //! Has the time step been changed?
        bool                            bSetTimeStep_;
        //! Has the time shift been set from the first frame?
        bool                            bHaveSetTimeShift_ = false;
        //! What is the correct time shift?
        real                            timeShift_ = 0;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<SetTime>
    SetTimePointer;

} // namespace gmx

#endif
