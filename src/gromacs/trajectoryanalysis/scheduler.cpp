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
 * Implements gmx::TrajectoryAnalysisCommandLineRunner.
 *
 * \author John Eblen <jeblen@acm.org>
 * \author Ryan Johnson <ryanphjohnson@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "scheduler.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "external/mpp/gmxmpp.h"
#include "runnercommon.h"

namespace gmx
{
    BlockScheduler::BlockScheduler(TrajectoryAnalysisRunnerCommon &tarc)
        : nextFrame(0), endFrame(0), isDone(false), runner(tarc) {}

    void BlockScheduler::init()
    {
#ifdef GMX_LIB_MPI
        size_t rank        = mpi::comm::world.rank();
        size_t totalFrames = runner.getTotalNumberOfFrames();
        size_t groupSize   = totalFrames / mpi::comm::world.size();
        size_t mod         = totalFrames % mpi::comm::world.size();

        if (rank < mod)
        {
            groupSize++;
            nextFrame = rank * groupSize;
        }
        else
        {
            nextFrame = rank * groupSize + mod;
        }
#else
        nextFrame = 0;
        size_t groupSize = runner.getTotalNumberOfFrames();
#endif
        endFrame = nextFrame + groupSize;
    }

    size_t BlockScheduler::selectNextFrameNumber()
    {
        if (nextFrame >= endFrame-1)
        {
            isDone = true;
        }
        size_t frame = nextFrame;
        ++nextFrame;
        return frame;
    }

    bool BlockScheduler::done() const
    {
        return isDone;
    }

    RoundRobinScheduler::RoundRobinScheduler()
        : nextFrame(0), inc(0) {}

    void RoundRobinScheduler::init()
    {
#ifdef GMX_LIB_MPI
        nextFrame = mpi::comm::world.rank();
        inc       = mpi::comm::world.size();

#else
        nextFrame = 0;
        inc = 1;
#endif
    }

    size_t RoundRobinScheduler::selectNextFrameNumber()
    {
        size_t frame = nextFrame;
        nextFrame += inc;
        return frame;
    }

    bool RoundRobinScheduler::done() const
    {
        return false;
    }
}
