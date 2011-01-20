/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \file
 * \brief
 * Declares gmx::TrajectoryAnalysisSettings and gmx::TopologyInformation.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISSETTINGS_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISSETTINGS_H

#include "../legacyheaders/typedefs.h"

namespace gmx
{

class Options;
class TrajectoryAnalysisRunnerCommon;

/*! \brief
 * Trajectory analysis module configuration object.
 *
 * This class is used by trajectory analysis modules to inform the caller
 * about the requirements they have on the input (e.g., whether a topology is
 * required, or whether PBC removal makes sense). It is also used to pass
 * similar information back to the analysis module after parsing user input.
 *
 * Having this functionality as a separate class makes the
 * TrajectoryAnalysisModule interface much cleaner, and also reduces the need to
 * change existing code when new options are added.
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
            /*! \brief
             * Requests dumps of parsed and compiled selection trees.
             *
             * This flag is used by internal debugging tools to request
             * the selection trees dumping to stderr.
             */
            efDebugSelection = 1<<16,
        };

        //! Initializes default settings.
        TrajectoryAnalysisSettings();
        ~TrajectoryAnalysisSettings();

        //! Returns the currently set flags.
        unsigned long flags() const;
        //! Tests whether a flag has been set.
        bool hasFlag(unsigned long flag) const;
        /*! \brief
         * Returns whether PBC should be used.
         *
         * Returns the value set with setPBC() and/or overridden by the user.
         * The user-provided value can be accessed in
         * TrajectoryAnalysisModule::initOptionsDone(), and can be overridden
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
        int setFlags(unsigned long flags);
        //! Sets or clears an individual flag.
        int setFlag(unsigned long flag, bool bSet = true);
        /*! \brief
         * Sets whether PBC are used.
         *
         * \param[in]  bPBC   TRUE if PBC should be used.
         * \returns    0 on success.
         *
         * If called in TrajectoryAnalysisModule::initOptions(), this function
         * sets the default for whether PBC are used in the analysis.
         * If ::efNoUserPBC is not set, a command-line option is provided
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
         * \see ::efNoUserPBC
         */
        int setPBC(bool bPBC);
        /*! \brief
         * Sets whether molecules are made whole.
         *
         * \param[in]     bRmPBC TRUE if molecules should be made whole.
         * \returns       0 on success.
         *
         * If called in TrajectoryAnalysisModule::initOptions(), this function
         * sets the default for whether molecules are made whole.
         * If ::efNoUserRmPBC is not set, a command-line option is provided
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
         * In such a case, you can also specify ::efNoUserRmPBC to not to
         * confuse the user with an option that would only slow the program
         * down.
         *
         * \see ::efNoUserRmPBC
         */
        int setRmPBC(bool bRmPBC);
        /*! \brief
         * Sets flags that determine what to read from the trajectory.
         *
         * \param[in]     frflags Flags for what to read from the trajectory file.
         * \returns       0 on success, an error code on error.
         *
         * If this function is not called, the flags default to TRX_NEED_X.
         * If the analysis module needs some other information (velocities,
         * forces), it can call this function to load additional information
         * from the trajectory.
         */
        int setFrameFlags(int frflags);

    private:
        class Impl;

        Impl                   *_impl;

        // Disallow copy and assign.
        TrajectoryAnalysisSettings(const TrajectoryAnalysisSettings &);
        void operator =(const TrajectoryAnalysisSettings &);

        friend class TrajectoryAnalysisRunnerCommon;
};

/*! \brief
 * Topology information passed to a trajectory analysis module.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class TopologyInformation
{
    public:
        //! Returns true if a topology file was loaded.
        bool hasTopology() const { return _top != NULL; }
        //! Returns true if a full topology file was loaded.
        bool hasFullTopology() const { return _bTop; }
        //! Returns the loaded topology, or NULL if not loaded.
        t_topology *topology() const { return _top; }
        //! Returns the ePBC field from the topology.
        int ePBC() const { return _ePBC; }
        /*! \brief
         * Gets the configuration from the topology.
         *
         * \param[out] x     Topology coordinate pointer to initialize.
         *      (can be NULL, in which case it is not used).
         * \param[out] box   Box size from the topology file
         *      (can be NULL, in which case it is not used).
         * \returns    0 on success, a non-zero error code on error.
         *
         * If TrajectoryAnalysisSettings::efUseTopX has not been specified,
         * \p x should be NULL.
         *
         * The pointer returned in \p *x should not be freed.
         */
        int getTopologyConf(rvec **x, matrix box) const;

    private:
        TopologyInformation();
        ~TopologyInformation();

        //! The topology structure, or NULL if no topology loaded.
        t_topology          *_top;
        //! true if full tpx file was loaded, false otherwise.
        bool                 _bTop;
        //! Coordinates from the topology (can be NULL).
        rvec                *_xtop;
        //! The box loaded from the topology file.
        matrix               _boxtop;
        //! The ePBC field loaded from the topology file.
        int                  _ePBC;

        // Disallow copy and assign.
        TopologyInformation(const TopologyInformation &);
        void operator =(const TopologyInformation &);

        friend class TrajectoryAnalysisRunnerCommon;
};

} // namespace gmx

#endif
