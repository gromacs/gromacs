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
 * Wrapper object to perform all basic convert functions for writing new coordinate files.
 *
 * \author
 * \inpublicapi
 * \ingroup fileio
 */
#ifndef GMX_FILEIO_OUTPUTFLAGS_H
#define GMX_FILEIO_OUTPUTFLAGS_H

#include <algorithm>

#include "gromacs/fileio/outputmanager.h"
#include "gromacs/options/ioptionscontainer.h"

namespace gmx
{

/*!\brief
 * OutputFlags class
 *
 * \inpublicapi
 * \ingroup fileio
 *
 */
class OutputFlags
{
    public:
        /*! \brief
         * Default constructor for OutputFlags.
         */
        OutputFlags() {}

        virtual ~OutputFlags() {}
        /*! \brief
         * Pass any user input options to the frame manager.
         *
         * Used to build the wrapper object used here to register all the modules and options.
         * \param[in] options Pointer to global options framework.
         */
        void initFileOptions(IOptionsContainer *options);

        /*! \brief
         * Sanity check for user input options.
         *
         * This function performs the check of the user input for basic sanity issues
         * and should be called after option processing has been finished.
         */
        void checkOptions();

        /*! \brief
         * Register modules used in the wrapper.
         *
         * Adds the convert related framemanager modules to the wrapper and primes them
         * for running in the analysis chain. Only registers modules that might actually
         * be useful.
         *
         * \param[in] output Pointer to method that registers the module.
         */
        void registerModules(OutputManagerPointer &output);

        //! Should velocities be written.
        ChangeSettingType                                  velocity_ = ChangeSettingType::efUnchanged;
        //! Should forces be written.
        ChangeSettingType                                  force_ = ChangeSettingType::efUnchanged;
        //! Should precision be changed.
        ChangeFrameUnchangedYesType                        precision_ = ChangeFrameUnchangedYesType::efUnchanged;
        //! Precision used in output file.
        int                                                prec_ = 3;
        //! Should frame start time be changed.
        ChangeFrameTimeType                                frameTime_ = ChangeFrameTimeType::efUnchanged;
        //! Time for first frame to start.
        double                                             startTimeValue_ = 0;
        //! Check if only start time will be changed.
        //! Time step to use between frames.
        double                                timeStepValue_  = 0;
        //! User supplied diagonal box vector.
        std::vector<double>                   newBoxVector_;
        //! Box vector converted to matrix format.
        matrix                                newBox_ = {{0}};
        //! Should frame box be changed.
        ChangeFrameUnchangedYesType                                  box_ = ChangeFrameUnchangedYesType::efUnchanged;
        //! Should frame atom setting be changed.
        ChangeSettingType atoms_ = ChangeSettingType::efUnchanged;
};

} // namespace gmx

#endif
