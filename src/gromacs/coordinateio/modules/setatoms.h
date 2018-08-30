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
 * Declares gmx::SetAtoms.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_SETATOMS_H
#define GMX_FILEIO_SETATOMS_H

#include <algorithm>

#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/topology/atoms.h"

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
 * \ingroup module_coordinateio
 *
 */
class SetAtoms : public ICoordinateOutput
{
    public:
        /*! \brief
         * Construct SetAtoms object with choice for boolean value
         * for availability of the t_atoms struct.
         *
         * Can be used to initialize SetAtoms from outside of trajectoryanalysis
         * framework.
         */
        explicit SetAtoms(ChangeSettingType atomFlag) :
            ICoordinateOutput(efAtomOutput),
            atomFlag_(atomFlag)
        {
            init_atom(&atoms_);
        }
        /*! \brief
         * Move constructor for SetAtoms.
         */
        SetAtoms &operator=(SetAtoms &&old) noexcept
        {
            atoms_    = old.atoms_;
            atomFlag_ = old.atomFlag_;
            return *this;
        }
        /*! \brief
         *  Move constructor for SetAtoms.
         */
        SetAtoms(SetAtoms &&old) noexcept :
            ICoordinateOutput(old.moduleFlags_),
            atomFlag_(old.atomFlag_),
            atoms_(old.atoms_)
        {
        }

        ~SetAtoms() {}

        /*! \brief
         * Change coordinate frame information for output.
         *
         * Changes the frame t_atoms struct according to user choice and
         * availability.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void processFrame(int /*framenumber*/, t_trxframe *input);

    private:
        //! Local function to check that we have a proper t_atoms struct available.
        bool haveLocalAtoms()
        {
            return haveAtoms();
        }
        /*! \brief
         *  Checking if t_trxframe has the atom information saved within.
         *
         *  \param[in] input t_trxframe before we start modifying it.
         */
        bool haveFrameAtoms(const t_trxframe &input)
        {
            return input.bAtoms;
        }
        //! Test if the atoms data is available for writing.
        bool haveAtoms() const { return atoms_.nr > 0; }
        //! Return pointer to t_atoms.
        t_atoms *atoms()
        {
            GMX_ASSERT(haveAtoms(), "No atoms information available");
            return &atoms_;
        }
        /*! \brief
         * Set atoms information from user input.
         *
         * \param[in] atoms Pointer to atoms information.
         */
        bool setAtoms(const t_atoms *atoms)
        {
            if (atoms != nullptr) {atoms_ = *atoms; }
            return atoms != nullptr;
        }
        //! Flag provided for setting atoms in coordinate frame from user.
        ChangeSettingType atomFlag_;
        //! Atoms information if available.
        t_atoms           atoms_;

};

//! Smart pointer to manage the object.
typedef std::unique_ptr<SetAtoms>
    SetAtomsPointer;

} // namespace gmx

#endif
