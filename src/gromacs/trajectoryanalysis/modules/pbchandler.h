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
 * Helper classes for PBC change
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_PBCHANDLER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_PBCHANDLER_H

#include <algorithm>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Enum value to store the selection for the PBC type from `-pbc`.
enum PBCType
{
    ePBCTypeMolecule,
    ePBCTypeResidue,
    ePBCTypeAtom,
    ePBCTypeNoJump,
    ePBCTypeCluster,
    ePBCTypeWhole,
    ePBCNotSet
};

//! Enum value to store the selection for the unit cell type from `-ur`.
enum UnitCellType
{
    eUnitCellTypeRectangular,
    eUnitCellTypeTriclinic,
    eUnitCellTypeCompact,
    eUnitCellNotSet
};

//! Enum value to store the selection for the centering from `-center`.
enum CenteringType
{
    eCenteringTypeTriclinic,
    eCenteringTypeRectangular,
    eCenteringTypeZero,
    eCenteringNotSet
};

//! Strings corresponding to PBCType.
const char *const    c_pbcTypes[] = {"mol", "res", "atom", "nojump", "cluster", "whole"};
//! Strings corresponding to UnitCellType.
const char *const    c_unitCellTypes[] = {"rect", "tric", "compact"};
//! Strings corresponding to CenteringType.
const char *const    c_centeringTypes[] = {"tric", "rect", "zero"};


/*! \libinternal
 * \brief
 * Allows changing the PBC representation of a system.
 *
 * This class makes it possible to change the PBC representation of a coordinate frame according
 * to user input specifications, and to hand the new coordiantes on for writing them to a new coordinate
 * file. The user has to provide the correct options to the class through the options interface in the
 * trajectoryanalysis framework to allow proper processing.
 */

class Pbchandler
{
    public:
        Pbchandler() :  selPBCType_(ePBCNotSet), selUnitCellType_(eUnitCellNotSet), selCenteringType_(eCenteringNotSet), selCOM_(nullptr), selCent_(nullptr), selClust_(nullptr), bCenter_(false), bHaveSelClust_(false),  bHaveSelCent_(false), bHaveSelCom_(false), bUseNewBox_(false)
        {
            clear_mat(box_);
        }
        ~Pbchandler()
        {}

        /*! \brief
         * Set user options through IOptions interface.
         *
         * Allows the setting of user facing options through the options interface
         * in the trajecotryanalysis framework.
         *
         * \param[in,out] options Pointer to options interface.
         */
        void initPBCOptions(IOptionsContainer *options);

        /*! \brief
         * Quick sanity check of user options.
         *
         * Makes sure that user options can be properly processed and all information needed to perform
         * the PBC change is available. Sets the selections that have been left unassigned to default
         * values if needed by the analysis for centering and COM based PBC actions.
         *
         * \param[in] sel Default output selection.
         * \param[in,out] settings Settings object to change pbc options.
         */
        void checkOptions(const Selection *sel, TrajectoryAnalysisSettings *settings);

        /*! \brief
         * Perform actual PBC changing operation.
         *
         * Changes coordinate representation of the output frame to obtain the representation
         * after the change to the PBC. The function relies on the proper setup of the class before
         * to be used with user supplied coordinates.
         *
         * \param[in] input Input coordinate frame used to obtain current coordinates.
         * \param[out] output Output coordinate frame that will be actually manipulated.
         * \param[in] mtop Topology object needed for a number of coordinate changes
         */
        void doPBC(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);

        /*! \brief
         * Set internal reference coordinates from topology.
         *
         * Makes the topology coordinates available to remove jumps of molecules across
         * PBC boundaries.
         *
         * \param[in] top Topology object from analysis framework.
         */
        void setReferenceCoordinates(const TopologyInformation &top);

        /*! \brief
         * Sets a new box for the frame output.
         *
         * In case the user wants to modify the box of the new frame, we set the new matrix for the box here.
         * \param[in] box New box in matrix format.
         */
        void setBox(const matrix box);

    private:
        /*! \brief
         * Remove molecule jumps across PBC boundaries
         *
         * \param[in,out] output Coordinate frame for manipulation.
         */
        void removeJump(t_trxframe *output);

        /*! \brief
         * Performs system centering according to user defined unit cell type.
         *
         * \param[in,out] output Coordinate frame for manipulation.
         * \param[in] centeringType How the system will be centered.
         */
        void centerSystem(const t_trxframe *input, t_trxframe *output, int centeringType);

        /*! \brief
         * Performs PBC operation accoridng to clustering of system molecules.
         *
         * \param[in,out] output Coordinate frame for manipulation.
         * \param[in] mtop Topology information.
         */
        void cluster(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);

        /*! \brief
         * Changes PBC representation according to Atom center of mass.
         *
         * \param[in,out] output Coordinate frame for manipulation.
         */
        void comAtom(t_trxframe *output);

        /*! \brief
         * Changes PBC representation according to residue center of mass.
         *
         * \param[in] input Unmodified coordinates from original coordinate frame.
         * \param[in,out] output Coordinate frame for manipulation.
         * \param[in] mtop Topology information.
         */
        void comResidue(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);

        /*! \brief
         * Changes PBC representation according to molecule center of mass.
         *
         * \param[in] input Unmodified coordiantes from original coordinate frame.
         * \param[in,out] output Coordinate frame for manipulation.
         * \param[in] mtop Topology information.
         */
        void comMolecule(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);

        /*! \brief
         * Performs the change of coordinates to change PBC representation.
         *
         * \param[in] atomStart Index for starting the shift of coordinates.
         * \param[in] atomEnd Index for ending the shift of coordinates.
         * \param[in] shift Change to coordinates.
         * \param[out] coord Coordinates to change.
         * \param[in] checkValid Test if the coordinates need to be checked before shifting.
         */
        void shiftCoord(const int atomStart, const int atomEnd, rvec shift, rvec *coord, const bool checkValid);

        /*! \brief
         * PBC type that the coordinates should be changed into.
         */
        PBCType              selPBCType_;
        /*! \brief
         * Unitcell type for changing representation.
         */
        UnitCellType         selUnitCellType_;
        /*! \brief
         * Centering type for changing represenation.
         */
        CenteringType        selCenteringType_;
        /*! \brief
         * Base selection for output.
         */
        Selection            sel_;
        /*! \brief
         * Selection for center of mass fitting.
         */
        Selection            selCOM_;
        /*! \brief
         * Selection for Centering.
         */
        Selection            selCent_;
        /*! \brief
         * Selection for clustering.
         */
        Selection            selClust_;
        //! Do we center something.
        bool                 bCenter_;
        //! Is clustering demanded.
        bool                 bHaveSelClust_;
        //! Is centering demanded.
        bool                 bHaveSelCent_;
        //! Is a center of mass action demanded.
        bool                 bHaveSelCom_;
        //! Reference coordinate array.
        std::vector<RVec>    refCoord_;
        //! Values for new diagonal box.
        std::vector<float>   newBox_;
        //! Are we using a new box.
        bool                 bUseNewBox_;
        //! New box matrix if set.
        matrix               box_;
};

} // namespace gmx

#endif
