/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "caching.h"

#include <memory>

#include "gromacs/fileio/trxio.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

namespace gmx
{
namespace trajectoryanalysis
{

const char CacheInfo::name[]             = "cache";
const char CacheInfo::shortDescription[] = "Cache a frame of trajectory data";

TrajectoryAnalysisModulePointer CacheInfo::create()
{
    return TrajectoryAnalysisModulePointer(new CachingTafModule());
}

// Implement required functions from base class
// Is there a Gromacs convention for tagging unused function parameters?
void CachingTafModule::initOptions(IOptionsContainer          gmx_unused *options,
                                   TrajectoryAnalysisSettings            *settings)
{
    // TODO: convert TRX_ flags in trxio.h to named enum
    // TODO: update gmx::TrajectoryAnalysisSettings::set{,Frame}Flags() call signature and enum to identifiable types.
    // TODO: get more helpful error from read_next_frame when _NEED_ flags aren't present?
    // Note that memory is allocated for v and f even if they are not available
    // for reading.
    settings->setFrameFlags(TRX_NEED_X | TRX_READ_V | TRX_READ_F);
}

// Is there a Gromacs convention for tagging unused function parameters?
void CachingTafModule::initAnalysis(const TrajectoryAnalysisSettings gmx_unused &settings,
                                    const TopologyInformation gmx_unused        &top)
{}

// Is there a Gromacs convention for tagging unused function parameters?
void CachingTafModule::analyzeFrame(int                          gmx_unused   frnr,
                                    const t_trxframe                         &fr,
                                    t_pbc                        gmx_unused  *pbc,
                                    TrajectoryAnalysisModuleData gmx_unused  *pdata)
{
    // Let's just grab a copy of the frame, using the trxframe interface.
    // We can stash it in an AnalysisData object, but don't need to for a first draft.
    // Copy assign the member shared_ptr with the input frame.
    last_frame_ = gmx::trajectory::trxframe_copy(fr);
    /*
       Note AnalysisData != TrajectoryAnalysisModuleData. The TrajectoryAnalysisModuleData object provided by the runner mediates access
       to the AnalysisData members of *this as configured with initAnalysis() and startFrames(). The Selection objects available via the TrajectoryAnalysisModuleData are probably more useful than direct trxframe access, anyway.
     */
    /*! TODO: use data module for analyzeFrame to retain a shared_ptr
     * to the trajectory frame data rather than copying fr.
     */
}

// Is there a Gromacs convention for tagging unused function parameters?
void CachingTafModule::finishAnalysis(int gmx_unused nframes)
{
    // If we're just caching trajectories, there is no post-processing.
}

/// Does not produce output unless requested.
void CachingTafModule::writeOutput() {}

/// Get a shared handle to the last frame processed.
std::shared_ptr<t_trxframe> CachingTafModule::frame() const
{
    return last_frame_;
}

} // end namespace trajectoryanalysis
} // end namespace gmx
