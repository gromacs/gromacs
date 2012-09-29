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
/*! \libinternal \file
 * \brief
 * Declares gmx::AnalysisDataParallelOptions.
 *
 * \if internal
 * Implementation of this class is currently in datastorage.cpp.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_PARALLELOPTIONS_H
#define GMX_ANALYSISDATA_PARALLELOPTIONS_H

namespace gmx
{

/*! \libinternal \brief
 * Parallelization options for analysis data objects.
 *
 * Methods in this class do not throw.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataParallelOptions
{
    public:
        //! Constructs options for serial execution.
        AnalysisDataParallelOptions();
        /*! \brief
         * Constructs options for parallel execution with given number of
         * concurrent frames.
         *
         * \param[in] parallelizationFactor
         *      Number of frames that may be constructed concurrently.
         *      Must be >= 1.
         */
        explicit AnalysisDataParallelOptions(int parallelizationFactor);

        //! Returns the number of frames that may be constructed concurrently.
        int parallelizationFactor() const { return parallelizationFactor_; }

        //! Returns whether MPI is used for parallelization
        bool useMPI() const;

    private:
        int                     parallelizationFactor_;
};

} // namespace gmx

#endif
