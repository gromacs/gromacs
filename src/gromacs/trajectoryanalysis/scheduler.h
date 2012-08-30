/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Declares gmx::TrajectoryAnalysisRunnerCommon.
 *
 * \author John Eblen <jeblen@acm.org>
 * \author Ryan Johnson <ryanphjohnson@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#pragma once

namespace gmx
{

/*! \internal \brief
 * Encapsulates a strategy for selecting the next frame to be analyzed.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisScheduler
{
public:
    virtual ~TrajectoryAnalysisScheduler() {}
    /*! \brief
     * Initialize scheduler - should be called after runner is initialized.
     */
    virtual void init() = 0;
    /*! \brief
     * Select a frame number to analyze. Every valid frame number should be
     * returned exactly once by one process. When there are no more valid
     * frames, the return value should be a number larger than the highest
     * possible frame.
     */
    virtual size_t selectNextFrameNumber() = 0;
    /*! \brief
     * Returns true if scheduler knows there are no more frames to process.
     * Some schedulers may always return false, since they may not have
     * enough information to know when frames have been exhausted.
     */
    virtual bool done() const = 0;
};

/*! \internal \brief
 * Block scheduler - must know the total number of frames but accurately reports
 * when frames are exhausted.
 */
class BlockScheduler : TrajectoryAnalysisScheduler
{
public:
    /*! \brief
     * Input runner for access to function getTotalNumberOfFrames()
     */
    // TODO: Make runner const once we fix the problems with counting number of frames.
    BlockScheduler(TrajectoryAnalysisRunnerCommon &runner);
    virtual void init();
    // Increments by one, even after frames exhausted.
    virtual size_t selectNextFrameNumber();
    // Returns true when frames exhausted.
    virtual bool done() const;
private:
    size_t nextFrame; // frame returned by selectNextFrameNumber()
    size_t endFrame; // last frame + 1
    bool isDone;
    // TODO: Make this const once we fix the problems with counting number of frames.
    TrajectoryAnalysisRunnerCommon &runner;
};

/*! \internal \brief
 * Round robin scheduler - does not need to know total number of frames but is
 * unable to report when frames are exhausted.
 */
class RoundRobinScheduler : TrajectoryAnalysisScheduler
{
public:
    virtual void init();
    // Increments by number of processes, even after frames exhausted.
    virtual size_t selectNextFrameNumber();
    // Always returns false.
    virtual bool done() const;
private:
    size_t nextFrame; // frame returned by selectNextFrameNumber()
    size_t inc; // increment amount
};

} // gmx
