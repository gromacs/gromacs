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
 * Helper classes for PBC change 
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_PBCHANDLER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_PBCHANDLER_H

#include <algorithm>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/unique_cptr.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"

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

class Pbchandler
{
    public:
        Pbchandler() :  selPBCType_(ePBCNotSet), selUnitCellType_(eUnitCellNotSet), selCenteringType_(eCenteringNotSet), selCOM_(nullptr), selCent_(nullptr), selClust_(nullptr)
    {}
        ~Pbchandler()
        {}

        void initPBCOptions(IOptionsContainer *options);
        void checkOptions(const Selection *sel, TrajectoryAnalysisSettings *settings);
        void doPBC(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);
        void setReferenceCoordinates(const TopologyInformation &top);
    private:
        void removeJump(t_trxframe *output);
        void centerSystem(t_trxframe *output, int centeringType);
        void cluster(t_trxframe *output, const gmx_mtop_t *mtop);
        void comAtom(t_trxframe *output);
        void comResidue(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);
        void comMolecule(const t_trxframe *input, t_trxframe *output, const gmx_mtop_t *mtop);
        void shiftCoord(const int atomStart, const int atomEnd, rvec shift, rvec *coord, const bool checkValid);

    PBCType    selPBCType_;
    UnitCellType    selUnitCellType_;
    CenteringType    selCenteringType_;
    Selection sel_;
    Selection selCOM_;
    Selection selCent_;
    Selection selClust_;
    bool    bCenter_;
    bool    bHaveSelClust_;
    bool    bHaveSelCent_;
    bool    bHaveSelCom_;
    std::vector<RVec>    refCoord_;
    std::vector<float>   newBox_;
    bool    bUseNewBox_;
    matrix  newMatrix_;


};

} // namespace gmx

#endif
