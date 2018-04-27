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
 * Interface class for frame handling, provides handles for all calls.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_FRAMEMANAGER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_FRAMEMANAGER_H

#include <algorithm>

#include "gromacs/fileio/trxio.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\libinternal
 * \brief
 * IFrameManager interface class taking care of all your frame handling needs.
 *
 * This interface is supposed to facilitate the usage of different frame handling
 * needs for writing and modifying trajectory coordinate frames. It is (supposed to be) written
 * in a way that makes it possible to chain different modules together and that all take
 * a trajectory frame as input and return a new trajectory frame after applying the given modification.
 *
 * \ingroup module_trajectoryanalysis
 *
 */
class IFrameManager
{
    public:
        /*! \brief
         * Construct IFrameManager object.
         *
         * Can be used to initialize IFrameManager from outside of trajectoryanalysis
         * framework.
         */
        IFrameManager() : prec_(3), time_(0),
                          bGC_(false), bP_(false), bV_(false), bF_(false), bA_(false), bT_(false),
                          bNewBox_(false), bAtoms_(false), fwflags_(0)
        {
            clear_trxframe(&newFrame_, true);
        }

        virtual ~IFrameManager()
        {
        }
        //! Move constructor for old object.
        explicit IFrameManager(IFrameManager &&old)
            : prec_(std::move(old.prec_)), time_(std::move(old.time_)), bGC_(std::move(old.bGC_)),
              bP_(std::move(old.bP_)), bV_(std::move(old.bV_)), bF_(std::move(old.bF_)),
              bA_(std::move(old.bA_)), bT_(std::move(old.bT_)), bNewBox_(std::move(old.bNewBox_)),
              bAtoms_(std::move(old.bAtoms_)), fwflags_(std::move(old.fwflags_)),
              newFrame_(std::move(old.newFrame_))
        {
        }

        //! Recognized flags for coordinate file writing.
        enum
        {
            /*! \brief
             * Ensures topology information is present.
             *
             * Similar flag to \p efRequireTop in the global \p TrajectoryAnalysisSettings
             * class, if passed then the coordinate writer needs to have access to a molecular topology.
             */
            efRequireMtop = 1<<0,
            /*! \brief
             * Ensures atoms data is available to writer.
             *
             * When writing plain text files, the t_atoms data structure is needed for atom information.
             * The writer will use this flag to check if this is needed and inform the
             * user of the t_atoms data could not be accessed from either \p mtop or
             * the coordinate frame itself.
             */
            efRequireAtoms = 1<<1,
            /*! \brief
             * Request that availability of connection data is checked.
             *
             * Some plain text formats can make use of connection data for file write out.
             * If the user has requested that connections are written, this flag will be set by the
             * output manager. Can lead to other flags being set, as connection data has to
             * be derived from the topology.
             */
            efRequireConnectionCheck = 1<<2
        };


        /*! \brief
         * Pass any user input options to the frame manager.
         *
         * Currently not used, will be useful to pass user input information to frame manager.
         */
        virtual void initFileOptions(IOptionsContainer *options) = 0;

        /*! \brief
         * Change coordinate frame information for output.
         *
         * Takes the previously internally stored coordinates and saves them again.
         * Applies correct number of atoms, as well as changing things such as
         * frame time or affect output of velocities or forces.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void modifyFrame(const t_trxframe &input) = 0;
        /*! \brief
         * Sanity check for user input options.
         *
         * This function performs the check of the user input for basic sanity issues
         * and should be called after option processing has been finished.
         */
        virtual void checkOptions() = 0;

        /*! \brief
         * Sets a coordinate frame writing flag.
         *
         * Sets or unsets a coordinate frame writing flag. Copy paste from
         * the routine in TrajectoryAnalysisSettings.
         */
        void setFlag(unsigned long flag, bool bSet)
        {
            if (bSet)
            {
                fwflags_ |= flag;
            }
            else
            {
                fwflags_ &= ~flag;
            }
        }

        /*! \brief
         * Returns true if the flag investigated has been set.
         */
        bool hasFlag(unsigned long flag)
        {
            return fwflags_ & flag;
        }

        /*! \brief
         * Checks if t_atoms are available for file output.
         */
        bool haveFrameAtoms() const
        {
            return newFrame_.bAtoms;
        }

        /*! \brief
         * Set the internal frame to the input value provided to the tool.
         */
        const t_trxframe getFrame() const
        {
            return newFrame_;
        }
        /*! \brief
         * Set the internal frame to the input value provided to the tool.
         */
        void setFrame(const t_trxframe &frame)
        {
            newFrame_ = frame;
        }
        //! Set the number of coordinates in the output.
        void setNatoms(int natoms)
        {
            newFrame_.natoms = natoms;
        }
        //! Set atoms data struct for current frame.
        void setFrameAtoms(const t_atoms *atoms, bool haveAtoms)
        {
            if (haveAtoms)
            {
                newFrame_.atoms     = const_cast<t_atoms*>(atoms);
            }
            else if (!haveFrameAtoms() && !haveAtoms)
            {
                newFrame_.atoms = nullptr;
            }
            newFrame_.bAtoms    = haveAtoms;
        }
        //! Get if we are having velocities in the frame.
        bool getbVel()
        {
            return newFrame_.bV;
        }
        /*! \brief
         * Change frame velocity settings.
         *
         * Corresponding function to \c getbVel to change frame settings for having
         * velocities saved.
         *
         * \param[in] velocity Boolean value to change frame setting.
         */
        void setbVel(bool velocity)
        {
            newFrame_.bV = velocity;
        }
        //! Get if we are having forces in the frame.
        bool getbForce()
        {
            return newFrame_.bF;
        }
        /*! \brief
         * Change frame force settings.
         *
         * Corresponding function to \c getbForce to change frame setting for having
         * forces saved.
         *
         * \param[in] force Boolean value to change frame setting.
         */
        void setbForce(bool force)
        {
            newFrame_.bF = force;
        }
        /*! \brief
         * Set frame coordinates from input.
         *
         * Sets the new frame coordinates to the user input given
         * by the tool.
         *
         * \param[in] x Pointer to coordinate array.
         */
        void setFrameCoordinates(rvec *x)
        {
            newFrame_.x = x;
        }
        /*! \brief
         * Set frame velocities from input.
         *
         * Sets new frame velocities (if possible) from input.
         * Asserts that the frame is meant to have velocities before doing things.
         *
         * \param[in] v Pointer to velocity array.
         */
        void setFrameVelocities(rvec *v)
        {
            GMX_RELEASE_ASSERT(newFrame_.bV, "Can not modify velocities in frame specified not to have them");
            newFrame_.v = v;
        }
        /*! \brief
         * Set frame forces from input.
         *
         * Sets new frame forces (if possible) from input.
         * Asserts that frame is meant to have forces.
         *
         * \param[in] f Pointer to force array.
         */
        void setFrameForces(rvec *f)
        {
            GMX_RELEASE_ASSERT(newFrame_.bF, "Can not modify forces in frame specified not to have them");
            newFrame_.f = f;
        }
        /*! \brief
         *  Resets the frame after it has been written out.
         *
         *  Because we (still) need to allocate the memory manually that is later
         *  written to the frame, we need to have a method that safely clears it after
         *  the frame has been written to disk.
         */
        void clearFrame()
        {
            sfree(newFrame_.x);
            if (getbVel())
            {
                sfree(newFrame_.v);
            }
            if (getbForce())
            {
                sfree(newFrame_.f);
            }
        }

        //! Precision for binary file formats.
        double                                prec_;
        //! Non default starting time.
        double                                time_;
        //! Do we need to write connection information.
        bool                                  bGC_;
        //! Does the frame have non default precision.
        bool                                  bP_;
        //! Does the frame have velocities.
        bool                                  bV_;
        //! Does the frame have forces.
        bool                                  bF_;
        //! Is the t_atoms struct included in the frame data.
        bool                                  bA_;
        //! Are we using non default starting and frame times.
        bool                                  bT_;
        //! Are we setting a new and different box.
        bool                                  bNewBox_;
        //! Is the t_atoms struct available.
        bool                                  bAtoms_;
        //! Flags containing information about the coordinate file writing.
        unsigned long                         fwflags_;
        //! Local storage of coordiante frame to work with.
        t_trxframe                            newFrame_;
};
/*! \brief
 * Typedef to have direct access to the individual FrameManager modules.
 */
typedef std::shared_ptr<IFrameManager>
    FrameManagerPointer;

} // namespace gmx

#endif
