/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Declares gmx::SetStartTime
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_SETSTARTTIME_H
#define GMX_FILEIO_SETSTARTTIME_H

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*!\brief
 * SetStartTime class allows changing trajectory time information.
 *
 * This class allows the user to set custom start time information for the
 * current frame in a trajectory.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class SetStartTime : public IOutputAdapter
{
public:
    /*! \brief
     * Construct object with choice for how to change initial time.
     *
     * \param[in] startTime User defined value for the initial time.
     */
    explicit SetStartTime(real startTime) :
        startTime_(startTime),
        haveProcessedFirstFrame_(false),
        differenceToInitialTime_(0)
    {
    }
    /*! \brief
     *  Move constructor for SetStartTime.
     */
    SetStartTime(SetStartTime&& old) noexcept = default;

    ~SetStartTime() override {}

    void processFrame(int /* framenumber */, t_trxframe* input) override;

    void checkAbilityDependencies(unsigned long /* abilities */) const override {}

private:
    /*! \brief
     * Set initial time from first processed frame.
     *
     * Calculates the time shift between the user set time and the time
     * in the coordinate frame being processed from the first processed coordinate
     * frame. This time shift is then used to calculate new frame times for each processed
     * coordinate frame.
     *
     * \param[in] initialTime Time value obtained from first frame.
     */
    void setInitialTime(real initialTime);

    /*! \brief
     * Stores the value of the initial time.
     *
     * In case users supply a new time step, the initial time of the
     * processed coordinate frame is stored here. In case the user also supplies
     * a new initial time, this variable is set to the user input instead.
     */
    real startTime_;
    //! Has the first frame been processed?
    bool haveProcessedFirstFrame_;
    /*! \brief
     * If the initial time is changed, we need to keep track of the initial
     * time difference to adjust the time of all following frames.
     */
    real differenceToInitialTime_;
};

//! Smart pointer to manage the object.
using SetStartTimePointer = std::unique_ptr<SetStartTime>;

} // namespace gmx

#endif
