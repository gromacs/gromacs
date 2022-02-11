/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal
 *
 * \brief Functionality for testing whether calls to mdrun produce the
 * same energy and force quantities when they should do so.
 */
#ifndef GMX_MDRUN_TESTS_MDRUNCOMPARISON_H
#define GMX_MDRUN_TESTS_MDRUNCOMPARISON_H

#include <functional>
#include <memory>

#include <gtest/gtest.h>

namespace gmx
{

namespace test
{

/*! \internal
 * \brief Manages returning a pair of frames from two
 * equivalent simulations that are meaningful to compare.
 *
 * \todo This is largely duplicated by ContinuationFrameManager. These
 * could be refactored into components that compare iterators to
 * frames.
 *
 * \tparam FrameReader  Has readNextFrame() and frame() methods
 *                      useful for returning successive Frame objects */
template<class FrameReader>
class FramePairManager
{
public:
    //! Convenience typedef
    typedef std::unique_ptr<FrameReader> FrameReaderPtr;
    //! Constructor
    FramePairManager(FrameReaderPtr first, FrameReaderPtr second) :
        first_(std::move(first)), second_(std::move(second))
    {
    }

private:
    /*! \brief Probe for a pair of valid frames, and return true if both are found.
     *
     * Give a test failure if exactly one frame is found, because
     * that file is longer than the other one, and this is not
     * expected behaviour. */
    bool shouldContinueComparing()
    {
        if (first_->readNextFrame())
        {
            if (second_->readNextFrame())
            {
                // Two valid next frames exist, so we should continue comparing.
                return true;
            }
            else
            {
                ADD_FAILURE() << "first file had at least one more frame than second file";
            }
        }
        else
        {
            if (second_->readNextFrame())
            {
                ADD_FAILURE() << "second file had at least one more frame than first file";
            }
            else
            {
                // Both files ran out of frames at the same time, which is the expected behaviour.
            }
        }
        // At least one file is out of frames, so should not continue comparing.
        return false;
    }

public:
    /*! \brief Compare all possible pairs of frames using \c compareTwoFrames.
     *
     * \tparam Frame  The type of frame used in the comparison (returned
     *                by FrameReader and used by compareTwoFrames). */
    template<class Frame>
    void compareAllFramePairs(std::function<void(const Frame&, const Frame&)> compareTwoFrames)
    {
        while (shouldContinueComparing())
        {
            Frame firstFrame  = first_->frame();
            Frame secondFrame = second_->frame();
            SCOPED_TRACE("Comparing frames from two runs '" + firstFrame.frameName() + "' and '"
                         + secondFrame.frameName() + "'");
            compareTwoFrames(firstFrame, secondFrame);
        }
    }

private:
    FrameReaderPtr first_;
    FrameReaderPtr second_;
};

} // namespace test
} // namespace gmx

#endif
