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
/*! \inpublicapi \file
 * \brief
 * Interface class for frame handling, provides handles for all calls.
 *
 * \author
 * \ingroup module_coordinatedata
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_FRAMECONVERTER_H
#define GMX_TRAJECTORYANALYSIS_MODULES_FRAMECONVERTER_H

#include <algorithm>

#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\inpublicapi
 * \brief
 * IFrameConverter interface class taking care of all your frame handling needs.
 *
 * This interface is supposed to facilitate the usage of different frame handling
 * needs for writing and modifying trajectory coordinate frames. It is (supposed to be) written
 * in a way that makes it possible to chain different modules together and that all take
 * a trajectory frame as input and return a new trajectory frame after applying the given modification.
 *
 * \ingroup module_coordinatedata
 *
 */
class IFrameConverter
{
    public:
        /*! \brief
         * Construct IFrameConverter object.
         *
         * Can be used to initialize IFrameConverter from outside of trajectoryanalysis
         * framework.
         */
        IFrameConverter() : fwflags_(0)
        {
            clear_trxframe(&coordinateFrame_, true);
        }

        virtual ~IFrameConverter()
        {
        }
        //! Move constructor for old object.
        explicit IFrameConverter(IFrameConverter &&old)
            : fwflags_(std::move(old.fwflags_)),
              coordinateFrame_(std::move(old.coordinateFrame_))
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
        virtual void convertFrame(const t_trxframe &input) = 0;
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
         * Return the internally stored coordinates.
         */
        const t_trxframe getFrame() const
        {
            return coordinateFrame_;
        }
        /*! \brief
         * Set the internal frame to the input value provided to the tool.
         */
        void setFrame(const t_trxframe &frame)
        {
            coordinateFrame_ = frame;
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
            coordinateFrame_.x = x;
        }
        /*! \brief
         * Set new box for coordiante frame.
         *
         * Allows users to specifiy new box information for the output coordinate
         * frame dueing analysis.
         *
         * \param[in] box New box matrix to be copied over the original one.
         */
        void setFrameBox(matrix box)
        {
            copy_mat(box, coordinateFrame_.box);
        }

        //! Flags containing information about the coordinate file writing.
        unsigned long                         fwflags_;
        //! Local storage of coordiante frame to work with.
        t_trxframe                            coordinateFrame_;
};
/*! \brief
 * Typedef to have direct access to the individual FrameConverter modules.
 */
typedef std::shared_ptr<IFrameConverter>
    FrameConverterPointer;

} // namespace gmx

#endif
