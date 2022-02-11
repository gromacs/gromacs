/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \file
 * \brief
 * Declares gmx::TopologyInformation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_TOPOLOGYINFORMATION_H
#define GMX_TRAJECTORYANALYSIS_TOPOLOGYINFORMATION_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/classhelpers.h"

//! Forward declaration
typedef struct gmx_rmpbc* gmx_rmpbc_t;
enum class PbcType : int;

namespace gmx
{

template<typename T>
class ArrayRef;

class TopologyInformation;
class TrajectoryAnalysisRunnerCommon;

/*! \libinternal
 * \brief Topology information available to a trajectory analysis module.
 *
 * This class is used to provide topology information to trajectory
 * analysis modules and to manage memory for them. Having a single
 * wrapper object instead of passing each item separately makes
 * TrajectoryAnalysisModule interface simpler, and also reduces the
 * need to change existing code if additional information is added.
 *
 * It is intended that eventually most clients of this class will be
 * analysis tools ported to the new analysis framework, but we will
 * use this infrastructure also from the legacy analysis tools during
 * the transition period. That will make it easier to put those tools
 * under tests, and eventually port them.
 *
 * Methods in this class do not throw if not explicitly stated.
 *
 * The main data content is constant once loaded, but some content is
 * constructed only when required (e.g. atoms_ and
 * expandedTopology_). Their data members are mutable, so that the
 * lazy construction idiom works properly. Some clients wish to modify
 * the t_atoms, so there is an efficient mechanism for them to get a
 * copy they can modify without disturbing this class. (The
 * implementation releases the cached lazily constructed atoms_, but
 * from the point of view of the caller, that is a copy.) The getters
 * and copy functions are const for callers, which correctly expresses
 * that the topology information is not being changed, merely copied
 * and presented in a different form.
 *
 * \ingroup module_trajectoryanalysis
 */
class TopologyInformation
{
public:
    //! Returns true if a topology file was loaded.
    bool hasTopology() const { return hasLoadedMtop_; }
    //! Returns true if a full topology file was loaded.
    bool hasFullTopology() const { return bTop_; }
    /*! \brief Builder function to fill the contents of
     * TopologyInformation in \c topInfo from \c filename.
     *
     * Different tools require, might need, would benefit from, or
     * do not need topology information. This functions implements
     * the two-phase construction that is currently needed to
     * support that.
     *
     * Any coordinate or run input file format will work, but the
     * kind of data available from the getter methods afterwards
     * will vary. For example, the mtop() available after reading
     * a plain structure file will have a single molecule block and
     * molecule type, regardless of contents.
     *
     * After reading, this object can return many kinds of primary
     * and derived data structures to its caller.
     *
     * \todo This should throw upon error but currently does
     * not. */
    void fillFromInputFile(const std::string& filename);
    /*! \brief Returns the loaded topology, or nullptr if not loaded. */
    gmx_mtop_t* mtop() const { return mtop_.get(); }
    //! Returns the loaded topology fully expanded, or nullptr if no topology is available.
    const gmx_localtop_t* expandedTopology() const;
    /*! \brief Returns a read-only handle to the fully expanded
     * atom data arrays, which might be valid but empty if no
     * topology is available. */
    const t_atoms* atoms() const;
    /*! \brief Copies the fully expanded atom data arrays, which
     * might be valid but empty if no topology is available. */
    AtomsDataPtr copyAtoms() const;
    //! Returns the pbcType field from the topology.
    PbcType pbcType() const { return pbcType_; }
    /*! \brief
     * Gets the configuration positions from the topology file.
     *
     * If TrajectoryAnalysisSettings::efUseTopX has not been specified,
     * this method should not be called.
     *
     * \throws  APIError if topology position coordinates are not available
     */
    ArrayRef<const RVec> x() const;
    /*! \brief
     * Gets the configuration velocities from the topology file.
     *
     * If TrajectoryAnalysisSettings::efUseTopV has not been specified,
     * this method should not be called.
     *
     * \throws  APIError if topology velocity coordinates are not available
     */
    ArrayRef<const RVec> v() const;
    /*! \brief
     * Gets the configuration box from the topology file.
     *
     * \param[out] box   Box size from the topology file, must not be nullptr.
     */
    void getBox(matrix box) const;
    /*! \brief Returns a name for the topology.
     *
     * If a full topology was read from a a file, returns the name
     * it contained, otherwise the empty string. */
    const char* name() const;

    TopologyInformation();
    ~TopologyInformation();

private:
    //! The topology structure, or nullptr if no topology loaded.
    std::unique_ptr<gmx_mtop_t> mtop_;
    //! Whether a topology has been loaded.
    bool hasLoadedMtop_;
    //! The fully expanded topology structure, nullptr if not yet constructed.
    mutable ExpandedTopologyPtr expandedTopology_;
    //! The fully expanded atoms data structure, nullptr if not yet constructed.
    mutable AtomsDataPtr atoms_;
    //! true if full tpx file was loaded, false otherwise.
    bool bTop_;
    //! Position coordinates from the topology (can be nullptr).
    std::vector<RVec> xtop_;
    //! Velocity coordinates from the topology (can be nullptr).
    std::vector<RVec> vtop_;
    //! The box loaded from the topology file.
    matrix boxtop_{};
    //! The pbcType field loaded from the topology file.
    PbcType pbcType_;

    // TODO This type is probably movable if we need that.
    GMX_DISALLOW_COPY_AND_ASSIGN(TopologyInformation);

    /*! \brief
     * Needed to initialize the data.
     */
    friend class TrajectoryAnalysisRunnerCommon;
};

//! Convenience overload useful for implementing legacy tools.
gmx_rmpbc_t gmx_rmpbc_init(const gmx::TopologyInformation& topInfo);

} // namespace gmx

#endif
