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
 * Declares gmx::TrajectoryAnalysisSettings and gmx::TopologyInformation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_TOPOLOGYINFORMATION_H
#define GMX_TRAJECTORYANALYSIS_TOPOLOGYINFORMATION_H

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/classhelpers.h"

struct gmx_mtop_t;
struct t_topology;

namespace gmx
{

class TopologyInformation;
class TrajectoryAnalysisRunnerCommon;

/*! \brief Builder function to fill the contents of
 * TopologyInformation in \c topInfo from \c filename.
 *
 * \todo Once we use better data structures within TopologyInformation
 * then it should support copy and move, and we won't need to pass a
 * pointer here. */
void fillTopologyInformationFromTprFile(const char          *filename,
                                        TopologyInformation *topInfo);

/*! \brief
 * Topology information passed to a trajectory analysis module.
 *
 * This class is used to pass topology information to trajectory analysis
 * modules and to manage memory for them.  Having a single wrapper object
 * instead of passing each item separately makes TrajectoryAnalysisModule
 * interface simpler, and also reduces the need to change existing code if
 * additional information is added.
 *
 * Methods in this class do not throw if not explicitly stated.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class TopologyInformation
{
    public:
        //! Returns true if a topology file was loaded.
        bool hasTopology() const { return mtop_ != nullptr; }
        //! Returns true if a full topology file was loaded.
        bool hasFullTopology() const { return bTop_; }
        //! Returns the loaded topology, or NULL if not loaded.
        const gmx_mtop_t *mtop() const { return mtop_.get(); }
        //! Returns the loaded topology, or NULL if not loaded.
        t_topology *topology() const;
        //! Returns the ePBC field from the topology.
        int ePBC() const { return ePBC_; }
        /*! \brief
         * Gets the configuration from the topology.
         *
         * \param[out] x     Topology coordinate pointer to initialize.
         *      (can be NULL, in which case it is not used).
         * \param[out] box   Box size from the topology file
         *      (can be NULL, in which case it is not used).
         * \throws  APIError if topology coordinates are not available and
         *      \p x is not NULL.
         *
         * If TrajectoryAnalysisSettings::efUseTopX has not been specified,
         * \p x should be NULL.
         *
         * The pointer returned in \p *x should not be freed.
         */
        void getTopologyConf(rvec **x, matrix box) const;

    private:
        TopologyInformation();
        ~TopologyInformation();

        std::unique_ptr<gmx_mtop_t> mtop_;
        //! The topology structure, or NULL if no topology loaded.
        // TODO: Replace fully with mtop.
        mutable t_topology  *top_;
        //! true if full tpx file was loaded, false otherwise.
        bool                 bTop_;
        //! Coordinates from the topology (can be NULL).
        rvec                *xtop_;
        //! The box loaded from the topology file.
        matrix               boxtop_;
        //! The ePBC field loaded from the topology file.
        int                  ePBC_;

        GMX_DISALLOW_COPY_AND_ASSIGN(TopologyInformation);

        /*! \brief
         * Needed to initialize the data.
         */
        friend class TrajectoryAnalysisRunnerCommon;
        /*! \brief
         * Initialize \c topInfo from the tpr \c filename.
         */
        friend void fillTopologyInformationFromTprFile(const char          *filename,
                                                       TopologyInformation *topInfo);
};

} // namespace gmx

#endif
