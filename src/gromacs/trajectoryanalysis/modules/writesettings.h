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
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_WRITESETTINGS_H
#define GMX_TRAJECTORYANALYSIS_MODULES_WRITESETTINGS_H

#include <algorithm>


#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/options/ioptionscontainer.h"
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

typedef struct t_writeFileDoubles
{
    double  prec;
    double  time;
} t_writeFileDoubles;
typedef struct t_writeFileBool
{
    bool    bGC;
    bool    bP;
    bool    bV;
    bool    bF;
    bool    bA;
    bool    bT;
    bool    bNewBox;
} t_writeFileBool;

class TrajectoryWriteSettings
{
    public:
        virtual ~TrajectoryWriteSettings();
        TrajectoryWriteSettings();

        t_writeFileBool   getWriteBool()
        {
            return writeBool_;

        }
        t_writeFileDoubles    getWriteDoubles()
        {
            return writeDoubles_;
        }
        const Selection *getSel()
        {
            return &sel_;
        }
        std::string getName()
        {
            return name_;
        }
        const gmx_mtop_t *getMtop() const
        {
            return mtop_;
        }

        bool getUseNewBox()
        {
            return writeBool_.bNewBox;
        }

        const std::vector<real> getNewBoxDimensions()
        {
            return newBox_;
        }
        void setName(const std::string name)
        {
            name_ = name;
        }
        void setMtop(const gmx_mtop_t *mtop)
        {
            mtop_ = mtop;
        }
        void setWriteBool(t_writeFileBool writeBool)
        {
            writeBool_ = writeBool;
        }
        void setWriteDoubles(t_writeFileDoubles writeDoubles)
        {
            writeDoubles_ = writeDoubles;
        }
        void setbGC(bool bGC)
        {
            writeBool_.bGC = bGC;
        }
        void setPrecision(double prec)
        {
            writeDoubles_.prec = prec;
        }
        void setbPrec(bool bPrec)
        {
            writeBool_.bP = bPrec;
        }
        void setTime(double time)
        {
            writeDoubles_.time = time;
        }
        void setbVel(bool bVel)
        {
            writeBool_.bV = bVel;
        }
        void setbForce(bool bForce)
        {
            writeBool_.bF = bForce;
        }
        void setbAtoms(bool bAtoms)
        {
            writeBool_.bA = bAtoms;
        }
        void setbTime(bool bTime)
        {
            writeBool_.bT = bTime;
        }

        void initWriteSettingsOptions(IOptionsContainer *options);
        void initWriteBool(t_writeFileBool *writeBool);
        void initWriteDouble(t_writeFileDoubles *writeDoubles);
        void checkOptions();
    private:
        t_writeFileDoubles    writeDoubles_;
        t_writeFileBool       writeBool_;
        std::string           name_;
        const gmx_mtop_t     *mtop_;
        Selection       sel_;
        std::vector<real> newBox_;

};

} // namespace gmx

#endif
