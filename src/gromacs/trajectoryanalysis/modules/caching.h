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

/*! \file
 * \brief Declares the CachingTafModule
 *
 * \internal
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_CACHING_H
#define GMX_TRAJECTORYANALYSIS_CACHING_H

#include "gromacs/trajectoryanalysis/analysismodule.h"

namespace gmx
{
namespace trajectoryanalysis
{
using gmx::TrajectoryAnalysisSettings;
using gmx::TopologyInformation;
using gmx::IOptionsContainer;

/*! \brief Provide a dummy module to grab copies of frames received.
 *
 * Objects of this class are useful for testing runners, pipelines,
 * and other proofs of concept. This rough draft should be replaced
 * soon with a class that uses the data modules to retrieve and
 * store trajectory info and to use selection processing.
 *
 * It may be helpful to reorganize a bit to allow modules to migrate from the
 * interface using unmanaged pointers to an alternative interface or API.
 * \internal
 * \ingroup module_trajectoryanalysis
 */
class CachingTafModule : public gmx::TrajectoryAnalysisModule
{
    public:
        // Implement required virtual functions from base class

        virtual void initOptions(IOptionsContainer         gmx_unused  *options,
                                 TrajectoryAnalysisSettings *settings);

        virtual void initAnalysis(const TrajectoryAnalysisSettings gmx_unused  &settings,
                                  const TopologyInformation      gmx_unused    &top);

        virtual void analyzeFrame(int                           gmx_unused frnr,
                                  const t_trxframe             &fr,
                                  t_pbc                        gmx_unused *pbc,
                                  TrajectoryAnalysisModuleData gmx_unused *pdata );

        virtual void finishAnalysis(int gmx_unused nframes);

        virtual void writeOutput();

        // Additional methods provided by this class.
        /*! \brief Get a shared pointer to the most recent frame.
         *
         * \return most recent frame seen by module.
         * \note This should not be a t_trxframe object, but an updated interface
         * to trajectory frame data. If managed data objects are not available, we
         * can use an AnalysisData object to keep shared pointers alive for selected
         * data.
         * \internal
         */
        std::shared_ptr<t_trxframe> frame() const;

    private:
        std::shared_ptr<t_trxframe> last_frame_; //!< cache the last frame read. \internal
};

/*! \brief Module info for Caching module
 *
 * Various code for registering modules requires a class
 * providing these three members.
 * TODO: this requirement may only be implicit in code that manages modules.
 * \internal
 * \ingroup module_trajectoryanalysis
 */
class CacheInfo
{
    public:
        /// Name to register for module \internal
        static const char name[];
        /// Description for registration \internal
        static const char shortDescription[];

        /*! \brief Get pointer for registering module
         *
         * \return pointer to a new CachingTafModule
         * \internal
         */
        static TrajectoryAnalysisModulePointer create();
};

}      // end namespace trajectoryanalysis
}      // end namespace gmx

#endif // GMX_TRAJECTORYANALYSIS_CACHING_H
