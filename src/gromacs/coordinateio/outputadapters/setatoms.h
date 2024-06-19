/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Declares gmx::SetAtoms.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_SETATOMS_H
#define GMX_COORDINATEIO_SETATOMS_H

#include <algorithm>
#include <memory>
#include <utility>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/unique_cptr.h"

struct t_trxframe;

namespace gmx
{

/*!\brief
 * SetAtoms class controls availability of atoms data.
 *
 * This modules allows the user to specify if a coordinate frame
 * should contain the t_atoms data structure or not, and sets it in the
 * new coordinate frame from either the current topology or from
 * the data in the coordinate frame. The data is later used to identify
 * if certain output file types are legal or not.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class SetAtoms : public IOutputAdapter
{
public:
    /*! \brief
     * Construct SetAtoms object with choice for boolean value
     * for availability of the t_atoms struct.
     *
     * Can be used to initialize SetAtoms from outside of trajectoryanalysis
     * framework.
     */
    explicit SetAtoms(ChangeAtomsType atomFlag, AtomsDataPtr inputAtoms) :
        atomFlag_(atomFlag), haveStructureFileAtoms_(false), atoms_(std::move(inputAtoms))
    {
        if (atoms_ != nullptr)
        {
            haveStructureFileAtoms_ = true;
        }
        if (atomFlag_ == ChangeAtomsType::Never)
        {
            moduleRequirements_ = CoordinateFileFlags::Base;
        }
        else
        {
            moduleRequirements_ = CoordinateFileFlags::RequireAtomInformation;
        }
    }
    /*! \brief
     *  Move constructor for SetAtoms.
     */
    SetAtoms(SetAtoms&& old) noexcept :
        atomFlag_(old.atomFlag_),
        haveStructureFileAtoms_(old.haveStructureFileAtoms_),
        atoms_(std::move(old.atoms_))
    {
    }

    ~SetAtoms() override {}

    /*! \brief
     * Change coordinate frame information for output.
     *
     * Changes the frame t_atoms struct according to user choice and
     * availability.
     *
     * \param[in] input Coordinate frame to be modified later.
     */
    void processFrame(int /*framenumber*/, t_trxframe* input) override;

    void checkAbilityDependencies(unsigned long abilities) const override;

private:
    //! Local function to check that we have a proper t_atoms struct available.
    bool haveStructureFileAtoms() const { return haveStructureFileAtoms_; }
    /*! \brief
     *  Checking if t_trxframe has the atom information saved within.
     *
     *  \param[in] input t_trxframe before we start modifying it.
     */
    static bool haveFrameAtoms(const t_trxframe& input);
    //! Test if the atoms data is available for writing.
    bool haveAtoms(const t_trxframe& input) const
    {
        return haveStructureFileAtoms() || haveFrameAtoms(input);
    }
    //! Return pointer to t_atoms.
    t_atoms* atoms()
    {
        GMX_RELEASE_ASSERT(haveStructureFileAtoms(), "No atoms information available");
        return atoms_.get();
    }
    //! Flag provided for setting atoms in coordinate frame from user.
    ChangeAtomsType atomFlag_;
    //! Flag set if input atoms have been valid from the beginning.
    bool haveStructureFileAtoms_;
    /*! \brief
     *  Atoms information if available.
     *
     *  Note, the module takes ownership of the information and
     *  will clean it up on exit.
     */
    AtomsDataPtr atoms_;
    //! Requirements obtained from user input.
    CoordinateFileFlags moduleRequirements_;
};

//! Smart pointer to manage the object.
using SetAtomsPointer = std::unique_ptr<SetAtoms>;

} // namespace gmx

#endif
