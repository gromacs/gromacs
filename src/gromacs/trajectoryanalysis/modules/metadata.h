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
 * Contains metadata struct declaration.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_METADATA_H
#define GMX_TRAJECTORYANALYSIS_MODULES_METADATA_H

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

/*!\libinternal
 * \brief
 * Container for input data values for file output.
 *
 * Struct that holds user information for writing time and file output precision,
 * as well as boolean values that decide on how file writing will be performed to new
 * coordinate files.
 * Can be extendet if needed to contain more values.
 *
 * Also contains flags for deciding which files will be written out.
 */
class TrajectoryOutputMetadata
{
    public:
        //! Default constructor creates object initialized to default values.
        TrajectoryOutputMetadata() : prec_(3), time_(0),
                                     bGC_(false), bP_(false), bV_(false), bF_(false), bA_(false), bT_(false),
                                     bNewBox_(false), bAtoms_(false), bMtop_(false), fwflags_(0)
        {
            init_atom(&atoms_);
        }

        ~TrajectoryOutputMetadata()
        {
            done_atom(&atoms_);
        }

        TrajectoryOutputMetadata(TrajectoryOutputMetadata &metadata)
            : prec_(std::move(metadata.prec_)), time_(std::move(metadata.time_)), bGC_(std::move(metadata.bGC_)),
              bP_(std::move(metadata.bP_)), bV_(std::move(metadata.bV_)), bF_(std::move(metadata.bF_)),
              bA_(std::move(metadata.bA_)), bT_(std::move(metadata.bT_)), bNewBox_(std::move(metadata.bNewBox_)),
              bAtoms_(std::move(metadata.bAtoms_)), bMtop_(std::move(metadata.bMtop_)),
              fwflags_(std::move(metadata.fwflags_))
        {
            init_atom(&atoms_);
        }

        //! Precision for binary file formats.
        double            prec_;
        //! Non default starting time.
        double            time_;
        //! Do we need to write connection information.
        bool              bGC_;
        //! Does the frame have non default precision.
        bool              bP_;
        //! Does the frame have velocities.
        bool              bV_;
        //! Does the frame have forces.
        bool              bF_;
        //! Is the t_atoms struct included in the frame data.
        bool              bA_;
        //! Are we using non default starting and frame times.
        bool              bT_;
        //! Are we setting a new and different box.
        bool              bNewBox_;
        //! Is the t_atoms struct available.
        bool              bAtoms_;
        //! Is the full topology available.
        bool              bMtop_;
        //! Flags containing information about the coordinate file writing.
        unsigned long     fwflags_;
        //! Local atoms struct for plain text files.
        t_atoms           atoms_;
        //! Local mtop topology information.
        gmx_mtop_t        mtop_;

        /*! \brief
         * Sets a coordinate frame writing flag.
         *
         * Sets or unsets a coordinate frame writing flag. Copy paste from
         * the routine in TrajectoryAnalysisSettings.
         */
        void setFlag(unsigned long flag, bool bSet)
        {
            if (bSet)
            {
                fwflags_ |= flag;
            }
            else
            {
                fwflags_ &= ~flag;
            }
        }

        /*! \brief
         * Returns true if the flag investigated has been set.
         */
        bool hasFlag(unsigned long flag)
        {
            return fwflags_ & flag;
        }

        /*! \brief
         * Set local topology information from Trajectoryanalysis framework.
         *
         * \param[in] atoms Pointer to atoms struct from topology.
         */
        void setAtoms(t_atoms atoms)
        {
            atoms_  = atoms;
            bAtoms_ = true;
        }

        /*! \brief
         * Set complete local topology. Also inits loval t_atoms during this.
         *
         * \param[in] mtop Pointer to global mtop topology.
         */
        void setMtop(const gmx_mtop_t *mtop)
        {
            mtop_ = *mtop;
            setAtoms(gmx_mtop_global_atoms(&mtop_));
            bMtop_ = true;
        }

        /*! \brief
         * Checks if t_atoms are available for file output.
         */
        bool haveAtoms()
        {
            return bAtoms_;
        }
        /*! \brief
         * Checks if full topology is available.
         */
        bool haveMtop()
        {
            return bMtop_;
        }

        /*! \brief
         * Get access to atom data for file writing.
         */
        t_atoms *getAtoms()
        {
            return &atoms_;
        }


        //! Recognized flags for coordinate file writing.
        enum
        {
            /*! \brief
             * Ensures topology information is present.
             *
             * Similar flag to \p efRequireTop in the global \p TrajectoryAnalysisSettings
             * class, if passed then the coordinate writer needs to have access to a molecular topology.
             */
            efRequireMtop = 1<<0,
            /*! \brief
             * Ensures atoms data is available to writer.
             *
             * When writing plain text files, the t_atoms data structure is needed for atom information.
             * The writer will use this flag to check if this is needed and inform the
             * user of the t_atoms data could not be accessed from either \p mtop or
             * the coordinate frame itself.
             */
            efRequireAtoms = 1<<1,
            /*! \brief
             * Request that availability of connection data is checked.
             *
             * Some plain text formats can make use of connection data for file write out.
             * If the user has requested that connections are written, this flag will be set by the
             * output manager. Can lead to other flags being set, as connection data has to
             * be derived from the topology.
             */
            efRequireConnectionCheck = 1<<2
        };



};

} // namespace gmx

#endif
