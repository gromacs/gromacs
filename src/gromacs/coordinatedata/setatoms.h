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
 * Module to declare if atoms data is available for file writing.
 *
 * \author
 * \inpublicapi
 * \ingroup module_coordinatedata
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_SETATOMS_H
#define GMX_TRAJECTORYANALYSIS_MODULES_SETATOMS_H

#include <algorithm>

#include "gromacs/coordinatedata/framemanager.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\brief
 * SetAtoms class controls availability of atoms data.
 *
 * This modules allows the user to specify if a coordiante frame
 * should contain the t_atoms data structure or not, and sets it in the
 * new coordinate frame from either the current topology or from 
 * the data in the coordinate frame. The data is later used to identify
 * if certain output file types are legal or not.
 *
 * \inpublicapi
 * \ingroup module_coordinatedata
 *
 */
class SetAtoms : public IFrameManager
{
    public:
        /*! \brief
         * Default constructor for SetAtoms should not be used.
         *
         * Class should only be initialized with at least the base selection.
         */
        SetAtoms() = delete;
        /*! \brief
         * Construct SetAtoms object with choice for boolean value
         * for availability of the t_atoms struct.
         *
         * Can be used to initialize SetAtoms from outside of trajectoryanalysis
         * framework.
         */
        explicit SetAtoms(const t_atoms *atoms) : localAtoms_(atoms)
        {
        }
        /*! \brief
         * Copy constructor.
         */
        SetAtoms(const SetAtoms &old) = delete;
        /*! \brief
         * Assignment operator.
         */
        SetAtoms &operator=(const SetAtoms &old) = delete;
        /*! \brief
         * Move constructor for SetAtoms.
         */
        SetAtoms &operator=(SetAtoms &&old)
        {
            localAtoms_ = std::move(old.localAtoms_);
            return *this;
        }
        /*! \brief
         *  Move constructor for SetAtoms.
         */
        SetAtoms(SetAtoms &&old) : localAtoms_(std::move(old.localAtoms_))
        {
        }

        ~SetAtoms() {}
        /*! \brief
         * Pass any user input options to the frame manager.
         *
         * Currently not used, will be useful to pass user input information to frame manager.
         */
        virtual void initFileOptions(IOptionsContainer * /*options*/);

        /*! \brief
         * Change coordinate frame information for output.
         *
         * Changes the frame t_atoms struct according to user choice and
         * availability.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void modifyFrame(const t_trxframe &input);
        /*! \brief
         * Sanity check for user input options.
         *
         * This function performs the check of the user input for basic sanity issues
         * and should be called after option processing has been finished.
         */
        virtual void checkOptions();

    private:
        //! Local function to check that we have a proper t_atoms struct available.
        bool haveLocalAtoms()
        {
            return localAtoms_->nr > 0 ? true : false;
        }
        /*! \brief
         * Copy of t_atoms struct to add to frame.
         *
         * Local pointer to a t_atoms struct provided either from the current topology
         * or from the currently processed frame.
         */
        const t_atoms                           *localAtoms_;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<SetAtoms>
    SetAtomsPointer;

} // namespace gmx

#endif
