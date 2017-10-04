/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISSETTINGS_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISSETTINGS_H

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/utility/classhelpers.h"

struct gmx_mtop_t;
struct t_topology;

namespace gmx
{

template <typename T> class ArrayRef;

class AnalysisDataPlotSettings;
class ICommandLineOptionsModuleSettings;
class Options;
class TrajectoryAnalysisRunnerCommon;

/*! \brief
 * Trajectory analysis module configuration object.
 *
 * This class is used by trajectory analysis modules to inform the caller
 * about the requirements they have on the input (e.g., whether a topology is
 * required, or whether PBC removal makes sense).  It is also used to pass
 * similar information back to the analysis module after parsing user input.
 *
 * Having this functionality as a separate class makes the
 * TrajectoryAnalysisModule interface much cleaner, and also reduces the need to
 * change existing code when new options are added.
 *
 * Methods in this class do not throw, except for the constructor, which may
 * throw an std::bad_alloc.
 *
 * \todo
 * Remove plain flags from the public interface.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisSettings
{
    public:
        //! Recognized flags.
        enum
        {
            /*! \brief
             * Forces loading of a topology file.
             *
             * If this flag is not specified, the topology file is loaded only
             * if it is provided on the command line explicitly.
             */
            efRequireTop     = 1<<0,
            /*! \brief
             * Requests topology coordinates.
             *
             * If this flag is specified, the coordinates loaded from the
             * topology can be accessed, otherwise they are not loaded.
             *
             * \see TopologyInformation
             */
            efUseTopX        = 1<<1,
            /*! \brief
             * Disallows the user from changing PBC handling.
             *
             * If this option is not specified, the analysis module (see
             * TrajectoryAnalysisModule::analyzeFrame()) may be passed a NULL
             * PBC structure, and it should be able to handle such a situation.
             *
             * \see setPBC()
             */
            efNoUserPBC      = 1<<4,
            /*! \brief
             * Disallows the user from changing PBC removal.
             *
             * \see setRmPBC()
             */
            efNoUserRmPBC    = 1<<5,
        };

        //! Initializes default settings.
        TrajectoryAnalysisSettings();
        ~TrajectoryAnalysisSettings();

        //! Injects command line options module settings for some methods to use.
        void setOptionsModuleSettings(ICommandLineOptionsModuleSettings *settings);

        //! Returns the time unit the user has requested.
        TimeUnit timeUnit() const;
        //! Returns common settings for analysis data plot modules.
        const AnalysisDataPlotSettings &plotSettings() const;

        //! Returns the currently set flags.
        unsigned long flags() const;
        //! Tests whether a flag has been set.
        bool hasFlag(unsigned long flag) const;
        /*! \brief
         * Returns whether PBC should be used.
         *
         * Returns the value set with setPBC() and/or overridden by the user.
         * The user-provided value can be accessed in
         * TrajectoryAnalysisModule::optionsFinished(), and can be overridden
         * with a call to setPBC().
         */
        bool hasPBC() const;
        /*! \brief
         * Returns whether molecules should be made whole.
         *
         * See hasPBC() for information on accessing or overriding the
         * user-provided value.
         */
        bool hasRmPBC() const;
        //! Returns the currently set frame flags.
        int frflags() const;

        /*! \brief
         * Sets flags.
         *
         * Overrides any earlier set flags.
         * By default, no flags are set.
         */
        void setFlags(unsigned long flags);
        //! Sets or clears an individual flag.
        void setFlag(unsigned long flag, bool bSet = true);
        /*! \brief
         * Sets whether PBC are used.
         *
         * \param[in]  bPBC   true if PBC should be used.
         *
         * If called in TrajectoryAnalysisModule::initOptions(), this function
         * sets the default for whether PBC are used in the analysis.
         * If \ref efNoUserPBC is not set, a command-line option is provided
         * for the user to override the default value.
         * If called later, it overrides the setting provided by the user or an
         * earlier call.
         *
         * If this function is not called, the default is to use PBC.
         *
         * If PBC are not used, the \p pbc pointer passed to
         * TrajectoryAnalysisModule::analyzeFrame() is NULL.
         * The value of the flag can also be accessed with hasPBC().
         *
         * \see efNoUserPBC
         */
        void setPBC(bool bPBC);
        /*! \brief
         * Sets whether molecules are made whole.
         *
         * \param[in]     bRmPBC true if molecules should be made whole.
         *
         * If called in TrajectoryAnalysisModule::initOptions(), this function
         * sets the default for whether molecules are made whole.
         * If \ref efNoUserRmPBC is not set, a command-line option is provided
         * for the user to override the default value.
         * If called later, it overrides the setting provided by the user or an
         * earlier call.
         *
         * If this function is not called, the default is to make molecules
         * whole.
         *
         * The main use of this function is to call it with \c false if your
         * analysis program does not require whole molecules as this can
         * increase the performance.
         * In such a case, you can also specify \ref efNoUserRmPBC to not to
         * confuse the user with an option that would only slow the program
         * down.
         *
         * \see efNoUserRmPBC
         */
        void setRmPBC(bool bRmPBC);
        /*! \brief
         * Sets flags that determine what to read from the trajectory.
         *
         * \param[in]     frflags Flags for what to read from the trajectory file.
         *
         * If this function is not called, the flags default to TRX_NEED_X.
         * If the analysis module needs some other information (velocities,
         * forces), it can call this function to load additional information
         * from the trajectory.
         */
        void setFrameFlags(int frflags);

        //! \copydoc ICommandLineOptionsModuleSettings::setHelpText()
        void setHelpText(const ArrayRef<const char *const> &help);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        friend class TrajectoryAnalysisRunnerCommon;
};

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
        const gmx_mtop_t *mtop() const { return mtop_; }
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

        gmx_mtop_t          *mtop_;
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
};

} // namespace gmx

#endif
