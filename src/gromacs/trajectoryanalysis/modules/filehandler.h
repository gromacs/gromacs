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
 * Handler for writing coordinate data to files.
 * Takes care of file formatting and handling of file correct output sizes.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_FILEHANDLER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_FILEHANDLER_H

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
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

/*! \libinternal
 * \brief
 * Container for double values for file output.
 *
 * Struct that holds user information for writing time and file output precision.
 * Can be extendet if needed to contain more values.
 */
typedef struct t_writeFileDoubles
{
    //! Precision for binary file formats.
    double prec;
    //! Non default starting time.
    double time;
} t_writeFileDoubles;

/*!\libinternal
 * \brief
 * Container for boolean values for file output.
 *
 * Struct that should hold all the boolean values for deciding what and how to write
 * the new coordinate files.
 */
typedef struct t_writeFileBool
{
    //! Do we need to write connection information.
    bool    bGC;
    //! Does the frame have non default precision.
    bool    bP;
    //! Does the frame have velocities.
    bool    bV;
    //! Does the frame have forces.
    bool    bF;
    //! Is the t_atoms struct included in the frame data.
    bool    bA;
    //! Are we using non default starting and frame times.
    bool    bT;
    //! Are we setting a new and different box.
    bool    bNewBox;
} t_writeFileBool;

/*!\libinternal
 * \brief
 * Filehandler class for handling trajectory file opening and data writing
 *
 * The filehandler keeps both track of the file being currently written too
 * and the correct number of coordinates that should be written, as well as
 * the output file type.
 * All the information is set by the user through the options mechanics
 * in the framework.
 *
 * \ingroup module_trajectoryanalysis
 *
 */
class Filehandler
{
    public:
        /*! \brief
         * Default constructor for Filehandler.
         */
        Filehandler() : sel_(nullptr), trr_(nullptr), filetype_(efNR), mtop_(nullptr)
        {
            initWriteBool(&writeBool_);
            initWriteDoubles(&writeDoubles_);
            connections_ = gmx_conect_init();
            init_atom(&atoms_);
            strcpy(filemode_, "w");
            clear_trxframe(&newFrame_, true);

        }
        /*! \brief
         * Construct Filehandler object with initial selection.
         *
         * Can be used to initialize Filehandler from outside of trajectoryanalysis
         * framework.
         * TODO Add initializers for the remaining fields.
         */
        Filehandler(Selection sel) : sel_(sel), trr_(nullptr), filetype_(efNR), mtop_(nullptr)
        {
            initWriteBool(&writeBool_);
            initWriteDoubles(&writeDoubles_);
            connections_ = gmx_conect_init();
            init_atom(&atoms_);
            strcpy(filemode_, "w");
            clear_trxframe(&newFrame_, true);
        }

        ~Filehandler()
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
         * Change coordinate frame information for output.
         *
         * Takes the `oldFrame` coordinates and changes them
         * internally for correct number of atoms, as well as changing things such as
         * frame time or affect output of velocities or forces.
         *
         * \param[in] oldFrame Coordinate frame to be modified for output.
         */
        void modifyFrame(const t_trxframe *oldFrame);
        /*! \brief
         * Writes new coordinate data from internal data to disk.
         */
        void writeFrame() const;
        /*! \brief
         * Closes new trajectory file after finishing the writing too it.
         *
         * Function to properly close opened output files after the file writing is finished.
         * Has to be called from parent or will be invoked when the filehandler itself is destroyed.
         */
        void closeFile();
        /*! \brief
         * Sets a new box for the frame output.
         *
         * In case the user wants to modify the box of the new frame, we set the new matrix for the box here.
         * \param[in] box New box in matrix format.
         */
        void setBox(const matrix box);

        /*! \brief
         * Get the currently processed coordinate frame for further processing.
         *
         * Function to access the internally stored trajectory frame for additional processing
         * after the initial information that will actually written to disk has been set.
         * This includes number of atoms, box, force and velocity output and frame times.
         *
         * \returns Trajectory frame pointer to internally stored frame.
         */
        t_trxframe *getFrameForModification();

        //! Returns if we need to set new box or not.
        bool getUseNewBox()
        {
            return writeBool_.bNewBox;
        }

        //! Returns the rectangular box dimensions as a vector.
        const std::vector<float> getNewBoxDimensions()
        {
            return newBox_;
        }

        //! Returns the internal Selection for the output.
        const Selection *getSel()
        {
            return &sel_;
        }

        //! Returns the pointer to the topology information obtained from the analysis framework.
        const gmx_mtop_t *getMtop() const
        {
            return mtop_;
        }

        /*! \brief
         * Set the pointer to the topology information.
         *
         * Makes the topology information available to the filehandler functions to set coordinates
         * and access molecule information for file writing.
         *
         * \param[in] mtop Pointer to topology from analysis framework.
         */
        void setMtop(const gmx_mtop_t *mtop)
        {
            mtop_ = mtop;
        }

        /*! \brief
         * Sanity check for user input options.
         *
         * This function performs the check of the user input for basic sanity issues
         * and should be called after option processing has been finished.
         */
        void checkOptions();

        // Remaining functions are internal and called by functions above.
    private:
        /*! \brief
         * Set the information needed for populating the atoms structure in coordinate files
         * through the `local` version of it.
         */
        void setLegacyInformation(t_atoms *local);
        /*! \brief
         * Individual function needed to open TNG files.
         *
         * Function used to open TNG formatted output files for writing. The function is
         * basically a wrapper about the lower level function that sets the necessary input
         * needed for opening the file.
         *
         * \param[in] filename Character string for the new filename.
         * \param[in] mode Character string indicating filemode, has to be 'w'.
         */
        t_trxstatus *trjOpenTng(const char *filename, const char *mode) const;
        /*! \brief
         * Function used to open all files execpt TNG files.
         *
         * This is just a wrapper around the open_trx function from the file handling routines.
         *
         * \param[in] filename Character string for the new filename.
         * \param[in] mode Character string indicating filemode, has to be 'w'.
         */
        t_trxstatus *trjOpenTrr(const char *filename, const char *mode) const;
        /*! \brief
         * Set connection information needed to write out plain text files.
         *
         * Plain text format files (pdb, gro, g96) need the information for which molecules are
         * connected for writing the new files. This function populates gmx_conect to make this
         * available for writing.
         */
        void setConnections();

        /*! \brief
         * Initialize boolean variables to false.
         *
         * \param[in,out] writeBool Boolean storage struct that will be initialized.
         */
        void initWriteBool(t_writeFileBool *writeBool);

        /*! \brief
         * Initialize double variables to 0.
         *
         * \param[in,out] writeDouble Double storage struct that will be initialized.
         */
        void initWriteDoubles(t_writeFileDoubles *writeDouble);


    private:
        /*! \brief
         * Selection of atoms that will be written to disk.
         *
         * Internal selection of atoms chosen by the user that will be written
         * to disk during processing. All actions that the filehandler performs
         * will only be on those atoms, with the remaining ones being not affected.
         */
        Selection                            sel_;
        /*! \brief
         * Name for the new coordinate file.
         *
         * User supplied name for the new coordinate file. Also used to determine the
         * file type through checking of the file extension.
         */
        std::string                          name_;

        /*! \brief
         * Vector of real values to set new box size if requested.
         *
         * Internal storage for new box formats if requested.
         */
        std::vector<float> newBox_;

        /*! \brief
         * Internal storage of modified coordinates.
         *
         * Local copy of the modified coordinates stored internally for writing and
         * possible use for additional modifications.
         */
        t_trxframe                 newFrame_;

        /*! \brief
         * Local storage for boolean options.
         */
        t_writeFileBool            writeBool_;
        /*! \brief
         * Local storage for double options.
         */
        t_writeFileDoubles         writeDoubles_;

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
         * Local storage of connection information for plain text files.
         *
         * The connections object needs to be populated for writing to plain text files
         * and is handled during the set-up of the file writing.
         */
        gmx_conect                 connections_;

        /*! \brief
         * Local copy of the atom information struct needed for writing to some coordinate file types.
         */
        t_atoms                    atoms_;

        /*! \brief
         * Storage of file type for determing what kind of file will be written to disk.
         */
        int filetype_;

        /*! \brief
         * Local copy of the topology so it is available for functions that rely on it.
         */
        const gmx_mtop_t          *mtop_;

};

} // namespace gmx

#endif
