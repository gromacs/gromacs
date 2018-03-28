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
/*! \internal \file
 * \brief
 * Helper classes for file handling of trajectory files.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_FILEHANDLER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_FILEHANDLER_H

#include <algorithm>

#include "writesettings.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/unique_cptr.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"

namespace gmx
{

/*! \libinternal \brief
 * Filehandler classes for handling trajectory file opening and data writing
 *
 * The filehandler keeps both track of the file being currently written too
 * and the correct number of coordinates that should be written, as well as
 * the output file type.
 * It is passed the file name (used to derive file type) and a settings object containing
 * the appropriate fields to set the number of coordinates and internal state of
 * the final coordinate file.
 *
 * \ingroup module_trajectoryanalysis
 *
 */


class Filehandler
{
    public:
    Filehandler() : sel_(nullptr), trr_(nullptr), infile_(nullptr), filetype_(efNR), mtop_(nullptr)
        {
            strcpy(filemode_, "w");
            clear_trxframe(&newFrame_, true);
        }
    Filehandler(Selection sel) : trr_(nullptr), infile_(nullptr), filetype_(efNR), mtop_(nullptr)
        {
            sel_ = sel;
            strcpy(filemode_, "w");
        }

    ~Filehandler()
        {
            closeFile();
        }

        void initFileOptions(IOptionsContainer *options);

        /*! \brief
         * Initializes the file writer with the appropriate `settings`.
         */
        void initOutput(TrajectoryWriteSettings *settings);
        /*! \brief
         * Takes the `oldFrame` coordinates and changes them
         * internally for correct number of atoms
         */
        void modifyFrame(const t_trxframe *oldFrame);
        /*! \brief
         * Writes new coordinate data from internal data to disk.
         */
        void writeFrame() const;
        /*! \brief
         * Closes new trajectory file after finishing the writing too it.
         * Has to be called from parent before destroying the created filehandler object.
         */
	void closeFile();
        /*! \brief
         * Sets a new box for the frame output.
         */
        void setBox(const matrix box);
	/* \brief
	 * Get the currently processed coordinate frame for further processing
	 * after it has been passed through the base modifyFrame to set the attributes.
	 */
	t_trxframe *getFrameForModification();
        // Remaining functions are internal and called by functions above.
    private:
        /*! \brief
         * Set the information needed for populating the atoms structure in coordinate files
         * through the `local` version of it.
         */
        void setLegacyInformation(t_atoms *local);
        /*! \brief
         * Individual function needed to open TNG files. Needs `filename`, filemode `mode`
         */
        t_trxstatus *trjOpenTng(const char *filename, const char *mode) const;
        /*! \brief
         * Function used to open all files execpt TNG files. Needs `filename` and
         * filemode `mode`.
         */
        t_trxstatus *trjOpenTrr(const char *filename, const char *mode) const;
        /*! \brief
         * Set connection information needed to write out text based files.
         */
        void setConnections();

        Selection                            sel_;
        std::string                          name_;

        TrajectoryWriteSettings             *settings_;


        t_trxframe                 newFrame_;
        t_writeFileBool            writeBool_;
        t_writeFileDoubles         writeDoubles_;

        t_trxstatus               *trr_;

        char                       filemode_[5];
        t_trxstatus               *infile_;
        gmx_conect                 connections_;
        t_atoms                    atoms_;
        int filetype_;
        const gmx_mtop_t          *mtop_;

};

} // namespace gmx

#endif
