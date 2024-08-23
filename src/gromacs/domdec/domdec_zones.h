/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief Declares the DomdecZones class and helper structs
 *
 * A DomdecZones object is (or should be) the only data structure that
 * is shared from the domain decomposition module to the Nbnxm module.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef GMX_DOMDEC_ZONES_H
#define GMX_DOMDEC_ZONES_H

#include <array>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/range.h"

struct gmx_ddbox_t;
struct gmx_domdec_t;

namespace gmx
{

template<typename>
class ArrayRef;

//! The maximum possible number of zones, 2 along each dimension in the eighth shell method
static constexpr int sc_maxNumZones = 8;
//! The maximum possible number of i-zones, half of the zones are needed to cover all pairs
static constexpr int sc_maxNumIZones = sc_maxNumZones / 2;

/*! \internal
 * \brief Triclinic corners and the Cartesian bounding box
 *
 * Note that at decomposition steps atoms can be, slightly, ouf of these bounds
 * when update groups are used. At other steps somes atoms can (and usually will)
 * move/diffuse outside of these bounds.
 */
struct gmx_domdec_zone_size_t
{
    /* Zone lower corner in triclinic coordinates         */
    gmx::RVec x0 = { 0, 0, 0 };
    /* Zone upper corner in triclinic coordinates         */
    gmx::RVec x1 = { 0, 0, 0 };
    /* Zone bounding box lower corner in Cartesian coords */
    gmx::RVec bb_x0 = { 0, 0, 0 };
    /* Zone bounding box upper corner in Cartesian coords */
    gmx::RVec bb_x1 = { 0, 0, 0 };
};

/*! \internal
 * \brief Class for handling atom ranges and dimensions of domain decomposition zones
 *
 * In the eighth-shell domain decomposition there is a home zone and up to seven
 * zones in the halo. Along each dimension which is decomposed, there are two rows of zones,
 * resulting in a maximum of 2^3=8 zones.
 * For computing pair interactions, each pair only once, we take atom pairs between pairs
 * of zones where the two zones are labelled "i" and "j", like atom pairs in the pair list.
 * This class provides a list of zone pairs (i,j) from which to collect atom pairs for computing
 * all pair interactions within a given cut-off distance. These pairs are provided as
 * range of consecutive j-zones for each i-zone.
 */
class DomdecZones
{
public:
    DomdecZones(ArrayRef<const int> ddDims);

    //! Returns the number of zones
    int numZones() const { return num_; }

    //! Returns the number of "i-zones" for non-bonded pair interactions
    int numIZones() const { return numIZones_; }

    //! Returns the j-zone range for i-zone \p iZone
    const Range<int>& jZoneRange(int iZone) const { return jZoneRanges_[iZone]; }

    //! Returns the shift of zone \p zone with respect to the home zone
    const IVec& shift(int zone) const { return shifts_[zone]; }

    //! Returns the atom index range for zone \p zone
    const Range<int> atomRange(int zone) const
    {
        GMX_ASSERT(zone <= lastAtomRangeEndSet_, "Cannot access zones beyond the last one set");

        return { atomRanges_[zone], atomRanges_[zone + 1] };
    }

    //! Returns the atom index range end for atoms from directly neighboring domains for zone \p zone
    int directNeighborAtomRangeEnd(int zone) const
    {
        return atomRanges_[zone] + numAtomFromDirectNeighbors_[zone];
    }

    //! Returns the j-atom index range for i-zone \p iZone
    Range<int> jAtomRange(const int iZone) const
    {
        GMX_ASSERT(iZone < numIZones_, "iZone should be in range");

        const auto& jZoneRange = jZoneRanges_[iZone];

        return { atomRanges_[*jZoneRange.begin()], atomRanges_[*jZoneRange.end()] };
    }

    //! Returns the sizes for zone \p zone
    const gmx_domdec_zone_size_t& sizes(int zone) const { return sizes_[zone]; }

    /*! \brief Sets the end of the atom range for zone \p zone
     *
     * \param[in] zone                       The zone to set the range for
     * \param[in] atomRangeEnd               The end of the atom range
     * \param[in] rangeIsFromDirectNeighbor  Whether this range comes from directly neighboring domains only
     */
    void setAtomRangeEnd(const int zone, const int atomRangeEnd, const bool rangeIsFromDirectNeighbor)
    {
        GMX_ASSERT(zone <= lastAtomRangeEndSet_ + 1,
                   "The zone should be at most one beyond the last one that was set");

        GMX_ASSERT(atomRangeEnd >= atomRanges_[zone],
                   "atomRangeEnd should be >= than the previous end");

        atomRanges_[zone + 1] = atomRangeEnd;

        if (zone == 0 || rangeIsFromDirectNeighbor)
        {
            numAtomFromDirectNeighbors_[zone] = atomRangeEnd - atomRanges_[zone];
        }

        lastAtomRangeEndSet_ = zone;
    }

    /*! \brief Set zone dimensions for a range of zones
     *
     * \param[in,out] dd          The domain decomposition struct
     * \param[in]     box         The box
     * \param[in]     ddbox       The domain decomposition box struct
     * \param[in]     zoneRange   The range of zones to set sizes for
     */
    void setSizes(const gmx_domdec_t&   dd,
                  const matrix          box,
                  const gmx_ddbox_t*    ddbox,
                  const gmx::Range<int> zoneRange);

private:
    //! The number of zones including the home zone
    int num_ = 0;
    //! The number of i-zones for computing pair interactions
    int numIZones_ = 0;
    //! The j-zones ranges for each i-zone;
    std::array<Range<int>, sc_maxNumIZones> jZoneRanges_;
    //! The shift of the zones with respect to the home zone
    std::array<IVec, sc_maxNumZones> shifts_;
    //! The atom range boundaries for the zones
    std::array<int, sc_maxNumZones + 1> atomRanges_ = { 0 };
    //! The number of atoms in the zone that are communicated from direct neighbors (these come first)
    std::array<int, sc_maxNumZones> numAtomFromDirectNeighbors_ = { 0 };
    //! Boundaries of the zones
    std::array<gmx_domdec_zone_size_t, sc_maxNumZones> sizes_;

    //! The last zone for which the end was set, used for asserting consistency of atomRanges_
    int lastAtomRangeEndSet_ = -1;
};

} // namespace gmx

#endif
