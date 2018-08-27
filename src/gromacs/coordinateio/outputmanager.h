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
#ifndef GMX_FILEIO_OUTPUTMANAGER_H
#define GMX_FILEIO_OUTPUTMANAGER_H

#include <algorithm>

#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

/*!\brief
 * OutputManager class for handling trajectory file opening and initialization.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class OutputManager : public IModuleCallback
{
    public:
        /*! \brief
         * Minimal used constructor with pointer to the correct selection used for file output.
         */
        explicit OutputManager(std::string name, const Selection *sel, const gmx_mtop_t *mtop) :
            IModuleCallback(efAnyOutputSupported),
            name_(name), trr_(nullptr), filetype_(efNR), sel_(sel), mtop_(mtop)
        {
            setFiletype();
            flagModules_.clear();
        }
        /*! \brief
         * Move assignment constructor for OutputManager.
         */
        OutputManager &operator=(OutputManager &&old)
        {
            name_        = std::move(old.name_);
            trr_         = std::move(old.trr_);
            filetype_    = std::move(old.filetype_);
            sel_         = std::move(old.sel_);
            mtop_        = std::move(old.mtop_);
            flagModules_ = std::move(old.flagModules_);
            return *this;
        }
        /*! \brief
         * Move constructor for OutputManager.
         */
        OutputManager(OutputManager &&old) : IModuleCallback(old.moduleFlags_),
                                             name_(std::move(old.name_)), trr_(std::move(old.trr_)),
                                             filetype_(std::move(old.filetype_)), sel_(std::move(old.sel_)),
                                             flagModules_(std::move(old.flagModules_)),
                                             mtop_(std::move(old.mtop_))
        {
        }

        ~OutputManager()
        {
            closeFile();
        }
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
        void writeFrame(const t_trxframe write) const;

        /*! \brief
         * Wrapper around processFrame to allow starting with a const
         * input of t_trxframe.
         *
         * \param[in] framenumber Number of frame being currently processed.
         * \param[in] input View of the constant t_trxframe object provided by the
         *                  method that calls the output manager.
         */
        void prepareFrame(const int framenumber, const t_trxframe &input);

        // Remaining functions are internal and called by functions above.
        /*! \brief
         * Individual function needed to open TNG files.
         *
         * Function used to open TNG formatted output files for writing. The function is
         * basically a wrapper about the lower level function that sets the necessary input
         * needed for opening the file.
         */
        t_trxstatus *trjOpenTng();
        /*! \brief
         * Function used to open all files execpt TNG files.
         *
         * This is just a wrapper around the open_trx function from the file handling routines.
         */
        t_trxstatus *trjOpenTrr();

        /*! \brief
         * Function used to request opening of a new PDB file.
         *
         * Needs to check for stuff such as
         * connections being needed or not, and has to determine how the writing takes place.
         */
        t_trxstatus *trjOpenPdb();

        /*! \brief
         * Function to request opening of nwe GRO output file.
         *
         * Has to perform some limited checks for how to open file and needs to validate that
         * t_atoms is available.
         */
        t_trxstatus *trjOpenGro();

        /*! \brief
         * Costum method to clear frame coordinates allocated before.
         */
        void clearCoordinateFrame(t_trxframe *frame) const;

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
         */
        const std::string                       filemode_ = "w";

        /*! \brief
         * Storage of file type for determing what kind of file will be written to disk.
         */
        int          filetype_;
        /*! \brief
         * Selection of atoms to write out.
         */
        const Selection   *sel_;

        /*! \libinternal \brief
         * Storage for flag modules that modify the state of a t_trxframe object.
         */
        struct CoordinateFlagModule
        {
            //! Initialize module, same as for frame converters.
            explicit CoordinateFlagModule(CoordinateOutputPointer module)
                : module(module)
            {
            }
            //! Pointer to the module
            CoordinateOutputPointer module;
        };
        //! List of flag settings modules
        typedef std::vector<CoordinateFlagModule> CoordinateFlagModuleList;

        //! Adds new module to the process chain.
        void addFlagModule(CoordinateOutputPointer module);

        //! Storage for list of flag changing modules.
        CoordinateFlagModuleList flagModules_;

        //! Pointer to topology information if available.
        const gmx_mtop_t *mtop_;

        //! Test if the topology data is available for writing.
        bool haveMtop() const { return mtop_; }
        //! Return pointer to topology.
        const gmx_mtop_t *getMtop() const { return haveMtop() ? mtop_ : nullptr; }
};

//! Smart pointer to manage the output manager object.
typedef std::unique_ptr<OutputManager>
    OutputManagerPointer;

} // namespace gmx

#endif
