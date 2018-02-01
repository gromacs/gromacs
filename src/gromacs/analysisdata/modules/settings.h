/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Declares gmx::TrajectoryDataWriteModule for writing trajectory data (into a file).
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 * \author 
 */
#ifndef GMX_TRAJECTORYDATA_MODULES_SETTINGS_H
#define GMX_TRAJECTORYDATA_MODULES_SETTINGS_H

#include <memory>
#include <string>

#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class AnalysisDataValue;
class IOptionsContainer;
class SelectionCollection;
class Selection;

/*! \brief
 * Common settings for data plots.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class TrajectoryDataWriteSettings
{
    public:
        virtual ~TrajectoryDataWriteSettings();

        //! Constructs default analysis plot settings.
        TrajectoryDataWriteSettings();


        //! Returns the selection collection set with setSelectionCollection().
        const SelectionCollection *selectionCollection() const
        {
            return selections_;
        }
        /*! \brief
         * Set selection collection to print as comments into the output.
         *
         * Formatted selection text from all selections in \p selections is
         * printed as comments in the output file.
         * If this method is not called, no selection information is written
         * to the output.
         */
        void setSelectionCollection(const SelectionCollection *selections);
        /*! \brief
         * Adds common options for setting plot options.
         *
         * \param[in,out] options Options object to which options are added.
         */
        void initOptions(IOptionsContainer *options);

        void setName(const std::string name)
        {
            name_ = name;
            setFiletype();
        }
        void setInputSel(const Selection *sel)
        {
            sel_ = sel;
        }
        void setTopology(const gmx_mtop_t *mtop, const t_topology *top)
        {
            mtop_ = mtop;
            top_ = top;

        }
        void setFiletype()
        {
            filetype_ = fn2ftp(name_.c_str());
        }
        void setbGC(bool bgenCon)
        {
            bgenCon_ = bgenCon;
        }
        void setPrecision(double prec)
        {
            prec_ = prec;
        }
        void setbPrec(bool bPrec)
        {
            bPrec_ = bPrec;
        }
        void setTime(double time)
        {
            time_ = time;
        }
        void setbVel(bool bVel)
        {
            bVel_ = bVel;
        }
        void setbForce(bool bForce)
        {
            bForce_ = bForce;
        }
        void setbAtoms(bool bAtoms)
        {
            bAtoms_ = bAtoms;
        }
        void setAtoms(t_atoms *atoms)
        {
            atoms_ = *atoms;
        }
        void setConnections()
        {
            connections_ = gmx_conect_generate(top_);
        }

        const gmx_mtop_t *getMtop() const
        {
            return mtop_;
        }
        std::string getName() const
        {
            return name_;
        }
        const Selection *getInputSel() const
        {
            return sel_;
        }
        int getFiletype() const
        {
            return filetype_;
        }
        bool getbGC() const
        {
            return bgenCon_;
        }
        double getPrecision() const
        {
            return prec_;
        }
        bool getbPrec() const
        {
            return bPrec_;
        }
        double getTime() const
        {
            return time_;
        }
        bool getbTime() const
        {
            return bTime_;
        }
        bool getbVel() const
        {
            return bVel_;
        }
        bool getbForce() const
        {
            return bForce_;
        }
        bool getbAtoms() const
        {
            return bAtoms_;
        }

        t_atoms *getAtoms()
        {
            return &atoms_;
        }

        gmx_conect getConnections()
        {
            return connections_;
        }


    private:
        const SelectionCollection *selections_;
        std::string name_;
        const Selection *sel_;
        int filetype_;
        const gmx_mtop_t *mtop_;
        const t_topology *top_;
        bool    bgenCon_;
        double prec_;
        bool bPrec_;
        double time_;
        bool bTime_;
        bool bVel_;
        bool bForce_;
        bool bAtoms_;
        t_atoms atoms_;
        gmx_conect connections_;
};

} // namespace gmx

#endif
