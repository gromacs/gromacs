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
 * Output manager takes care of opening files and writing output to them.
 * Takes care of file formatting and handling of file output sizes.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_OUTPUTMANAGER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_OUTPUTMANAGER_H

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
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "metadata.h"

namespace gmx
{

/*!\libinternal
 * \brief
 * OutputManager class for handling trajectory file opening and initialization.
 *
 * The fileopen class keeps track of the file being currently written too
 * and metadata needed to successfully finish writing to it. It also takes
 * care of setting the output file type.
 * All the information is set by the user through the options mechanics
 * in the framework.
 *
 * \ingroup module_trajectoryanalysis
 *
 */
class OutputManager : public TrajectoryOutputMetadata
{
    public:
        /*! \brief
         * Default constructor for OutputManager class should not be used and is thus deleted.
         */
        OutputManager() = delete;
        /*! \brief
         * Minimal used constructor with pointer to the correct selection used for file output.
         */
        OutputManager(Selection *sel) : trr_(nullptr), filetype_(efNR), sel_(sel)
        {
            // Set up the metadata object and initalize it to the default values.
            metadata_ = new TrajectoryOutputMetadata();
            strcpy(filemode_, "w");

        }

        ~OutputManager()
        {
            closeFile();
        }
        /*!\brief
         * Initialize settings through the options interface.
         *
         * The Filehandler is passed the information needed to properly set up the
         * new trajectory file directly from user information through the options
         * interface.
         *
         * \param[in,out] options Pointer to options interface passed from trajectoryanalysis framework.
         */
        void initFileOptions(IOptionsContainer *options);
        /*! \brief
         * Prepare file output.
         *
         * After processing user information, the new file can be opened, while also setting
         * legacy information and/or connections needed for plain text files.
         */
        void initOutput();

        /*! \brief
         * Closes new trajectory file after finishing the writing too it.
         *
         * Function to properly close opened output files after the file writing is finished.
         * Has to be called from parent or will be invoked when the filehandler itself is destroyed.
         */
        void closeFile();

        /*! \brief
         * Sanity check for user input options.
         *
         * This function performs the check of the user input for basic sanity issues
         * and should be called after option processing has been finished.
         *
         * \param[in] top global topology to check if information is avialable.
         */
        void checkOptions(const TopologyInformation *top);

        /*! \brief
         * Set the output filetype from the provided file name.
         *
         * Uses the filename provided by the user to set the correct output file type,
         * together with flags needed for checking the ability to provide the output later.
         */
        void setFiletype();

        /*! \brief
         * Write coordinates to disk.
         *
         * \param[in] write Coordinate frame object ready for writing.
         */
        void writeFrame(const t_trxframe *write) const;

        /*! \brief
         * Local storage for trajectory output metadata values.
         */
        TrajectoryOutputMetadata            *metadata_;
        /*! \brief
         * Get access to metadata object in outputmanager.
         */
        TrajectoryOutputMetadata *metadata()
        {
            return metadata_;
        }
        // Remaining functions are internal and called by functions above.
    private:
        /*! \brief
         * Individual function needed to open TNG files.
         *
         * Function used to open TNG formatted output files for writing. The function is
         * basically a wrapper about the lower level function that sets the necessary input
         * needed for opening the file.
         */
        t_trxstatus *trjOpenTng() const;
        /*! \brief
         * Function used to open all files execpt TNG files.
         *
         * This is just a wrapper around the open_trx function from the file handling routines.
         */
        t_trxstatus *trjOpenTrr() const;

        /*! \brief
         * Function used to request opening of a new PDB file.
         *
         * Needs to check for stuff such as
         * connections being needed or not, and has to determine how the writing takes place.
         */
        t_trxstatus *trjOpenPdb() const;

        /*! \brief
         * Function to request opening of nwe GRO output file.
         *
         * Has to perform some limited checks for how to open file and needs to validate that
         * t_atoms is available.
         */
        t_trxstatus *trjOpenGro() const;

        /*! \brief
         * Name for the new coordinate file.
         *
         * User supplied name for the new coordinate file. Also used to determine the
         * file type through checking of the file extension.
         */
        std::string                          name_;

        /*! \brief
         * File pointer to the coordinate file being written.
         */
        t_trxstatus               *trr_;

        /*! \brief
         * File mode for writing to coordinate files.
         *
         * TODO have this as a constant string.
         */
        char                       filemode_[5];

        /*! \brief
         * Storage of file type for determing what kind of file will be written to disk.
         */
        int          filetype_;
        /*! \brief
         * Selection of atoms to write out.
         */
        Selection   *sel_;

};

} // namespace gmx

#endif
