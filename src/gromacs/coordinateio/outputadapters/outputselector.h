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
 * Declares gmx::OutputSelector.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_OUTPUTSELECTOR_H
#define GMX_COORDINATEIO_OUTPUTSELECTOR_H

#include <algorithm>
#include <memory>
#include <vector>

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/gmxassert.h"

struct t_trxframe;

namespace gmx
{

/*!\brief
 * OutputSelector class controls setting which coordinates are actually written.
 *
 * This adapter selects a subset of the particles and their data for output.
 * The subset is specified via a selection object (see \ref module_selection).
 * The corresponding particle data (coordinates, velocities, forces) and
 * topology information is copied from as available from the input data.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class OutputSelector : public IOutputAdapter
{
public:
    /*! \brief
     * Construct OutputSelector object with initial selection.
     *
     * Can be used to initialize OutputSelector from outside of trajectoryanalysis
     * framework.
     */
    explicit OutputSelector(const Selection& sel) : sel_(sel), selectionAtoms_(nullptr)
    {
        GMX_RELEASE_ASSERT(sel.isValid() && sel.hasOnlyAtoms(),
                           "Need a valid selection out of simple atom indices");
    }
    /*! \brief
     *  Move assignment constructor for OutputSelector.
     */
    OutputSelector(OutputSelector&& old) noexcept = default;

    ~OutputSelector() override {}

    /*! \brief
     * Change coordinate frame information for output.
     *
     * Takes the previously internally stored coordinates and saves them again.
     * Applies correct number of atoms in this case.
     *
     * \param[in] input Coordinate frame to be modified later.
     */
    void processFrame(int /*framenumber*/, t_trxframe* input) override;

    void checkAbilityDependencies(unsigned long /* abilities */) const override {}

private:
    /*! \brief
     * Selection of atoms that will be written to disk.
     *
     * Internal selection of atoms chosen by the user that will be written
     * to disk during processing.
     */
    const Selection& sel_;
    /*! \brief
     * Local storage of modified atoms.
     *
     * When selection input information, we might will have to adjust the
     * atoms content to match the new output. To perform this, we keep track
     * of modified atoms information in this object that is not used by default.
     */
    AtomsDataPtr selectionAtoms_;
    //! Local storage for coordinates
    std::vector<RVec> localX_;
    //! Local storage for velocities
    std::vector<RVec> localV_;
    //! Local storage for forces
    std::vector<RVec> localF_;
    //! Local storage for atom indices
    std::vector<int> localIndex_;
};

//! Smart pointer to manage the object.
using OutputSelectorPointer = std::unique_ptr<OutputSelector>;

} // namespace gmx

#endif
