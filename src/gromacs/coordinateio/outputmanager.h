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
 * Output manager takes care of opening files and writing output to them.
 *
 * \author
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_OUTPUTMANAGER_H
#define GMX_COORDINATEIO_OUTPUTMANAGER_H

#include <algorithm>
#include <utility>

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/selection/selection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

/*!\brief
 * OutputManager class responsible for managing coordinate file output.
 *
 * The OutputManager has as its tasks the proper opening and closing of new
 * trajectory files written from provided trajectory data. It performs the actual
 * writing to disk and has the ability to recruit modules derived from
 * IOutputAdapter to perform changes to the meta information contained
 * within the t_trxframe datastructures.
 *
 * When writing files to disk, the manager decides if a deep copy of the
 * coordinate file information should be taken or not, depending on if changes
 * to the stored values are going to be done or not. If the values are changed,
 * the manager will clean up the information after each frame has been processed.
 *
 * \todo Does this make sense, or should the data be stored to avoid numerous
 * rounds of allocations?
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class OutputManager
{
    public:
        struct OutputManagerBuildHelper;

        ~OutputManager()
        {
            closeFile();
        }
        /*! \brief
         * Wrapper around processFrame to allow starting with a const
         * input of t_trxframe.
         *
         * \param[in] framenumber Number of frame being currently processed.
         * \param[in] input View of the constant t_trxframe object provided by the
         *                  method that calls the output manager.
         */
        void prepareFrame(int framenumber, const t_trxframe &input);

        //! Get registered filetype
        int getFiletype() const { return filetype_; }
    private:
        /*! \brief
         * Minimal used constructor with reference to the correct selection used for file output.
         */
        OutputManager(std::string               name,
                      int                       filetype,
                      const Selection          &sel,
                      const gmx_mtop_t         *mtop,
                      OutputAdapterContainer    adapters) :
            outputFileName_(std::move(name)),
            outputFile_(nullptr),
            filetype_(filetype),
            sel_(sel),
            mtop_(mtop),
            outputAdapters_(std::move(adapters))
        {
        }

        /*! \brief
         * Prepare file output.
         *
         * Opens the file on disk to allow writing of coordinate data.
         */
        void initOutput();

        /*! \brief
         * Closes new trajectory file after finishing the writing to it.
         *
         * Function to properly close opened output files after the file writing is finished.
         * Has to be called from parent or will be invoked when the filehandler itself is destroyed.
         */
        void closeFile();

        //! Name for the new coordinate file.
        std::string                          outputFileName_;

        //! File pointer to the coordinate file being written.
        t_trxstatus               *outputFile_;

        //! Storage of file type for determing what kind of file will be written to disk.
        int          filetype_;

        /*! \brief
         * Selection of atoms to write out.
         *
         * Currently, OutputManager expects that the lifetime of the selection is longer
         * than that of itself, and that the selection remains unchanged during this lifetime.
         * A better approach will be to pass the selection to it and expect it to
         * manage the lifetime instead.
         */
        const Selection   &sel_;

        //! Pointer to topology information if available.
        const gmx_mtop_t *mtop_;

        //! Storage for list of output adapters.
        OutputAdapterContainer outputAdapters_;

        //! Local storage for modified positions.
        std::vector<RVec> localX_;
        //! Local storage for modified velocities.
        std::vector<RVec> localV_;
        //! Local storage for modified forces.
        std::vector<RVec> localF_;
        //! Local storage for modified indices.
        std::vector<int>  localIndex_;
};

//! Smart pointer to manage the output manager object.
using OutputManagerPointer = std::unique_ptr<OutputManager>;

} // namespace gmx

#endif
