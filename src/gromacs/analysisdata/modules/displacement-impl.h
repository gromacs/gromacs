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
 * Declares private implementation class for
 * gmx::AnalysisDataDisplacementModule.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_DISPLACEMENT_IMPL_H
#define GMX_ANALYSISDATA_MODULES_DISPLACEMENT_IMPL_H

#include "displacement.h"

namespace gmx
{

class AbstractHistogramModule;

/*! \internal \brief
 * Private implementation class for AnalysisDataDisplacementModule.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataDisplacementModule::Impl
{
    public:
        Impl();
        ~Impl();

        //! Maximum number of particles for which the displacements are calculated.
        int                     nmax;
        //! Maximum time for which the displacements are needed.
        real                    tmax;
        //! Number of dimensions per data point.
        int                     ndim;

        //! TRUE if no frames have been read.
        bool                    bFirst;
        //! Stores the time of the first frame.
        real                    t0;
        //! Stores the time interval between frames.
        real                    dt;
        //! Stores the time of the current frame.
        real                    t;
        //! Stores the index in the store for the current positions.
        int                     ci;

        //! Maximum number of positions to store for a particle.
        int                     max_store;
        //! The total number of positions ever stored (can be larger than \p max_store).
        int                     nstored;
        //! Old values.
        real                   *oldval;
        //! The most recently calculated displacements.
        real                   *currd;

        //! Histogram module for calculating MSD histograms, or NULL if not set.
        AbstractHistogramModule *histm;
};

} // namespace gmx

#endif
