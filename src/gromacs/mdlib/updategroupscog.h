/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares the UpdateGroupsCog class for managing centers of mass of update groups
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_UPDATEGROUPSCOG
#define GMX_MDLIB_UPDATEGROUPSCOG

#include <vector>

#include "gromacs/domdec/hashedmap.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/defaultinitializationallocator.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;

namespace gmx
{
template<typename>
class Range;
class RangePartitioning;

/*! \libinternal
 * \brief Class for managing and computing centers of geometry of update groups
 */
class UpdateGroupsCog
{
public:
    /*! \brief Constructor
     *
     * The temperature is used for computing the maximum update group
     * radius.
     *
     * \param[in] mtop                            The global topology
     * \param[in] updateGroupingsPerMoleculeType  List of update groups for each molecule type in \p mtop
     * \param[in] temperature                     The maximum reference temperature, pass -1 when unknown or not applicable
     */
    UpdateGroupsCog(const gmx_mtop_t&                           mtop,
                    gmx::ArrayRef<const gmx::RangePartitioning> updateGroupingsPerMoleculeType,
                    real                                        temperature);

    /*! \brief Compute centers of geometry for supplied coordinates
     *
     * Coordinates are processed starting after the last index
     * processed in the previous call to \p addCogs(), unless \p clear()
     * was called last, in which case processing starts at 0.
     * Coordinates are processed until \p globalAtomIndices.size().
     *
     * \param[in] globalAtomIndices  List of global atom indices for the atoms belonging to \p coordinates
     * \param[in] coordinates        List of coordinates to be processed, processing might not start at 0 (see above)
     * \param[in] numAtomsPerNbnxmGridColumn  A list of the number of atoms per NBNxM grid column,
     *                                        used for OpenMP parallelization, as update groups are
     *                                        never split over columns; can be empty
     */
    void addCogs(ArrayRef<const int>  globalAtomIndices,
                 ArrayRef<const RVec> coordinates,
                 ArrayRef<const int>  numAtomsPerNbnxmGridColumn);

    /*! \brief Returns the number of centers of geometry currently stored */
    int numCogs() const { return cogs_.size(); }

    /*! \brief Returns a reference to a center of geometry
     *
     * \param[in] cogIndex  The COG requested, should be; 0 <= \p cogIndex < cogs_.size()
     */
    RVec& cog(int cogIndex)
    {
        GMX_ASSERT(cogIndex >= 0 && static_cast<size_t>(cogIndex) < cogs_.size(),
                   "cogIndex should be in the range set in this object");

        return cogs_[cogIndex];
    }

    /*! \brief Returns whether a COG entry consists of a (single) filler particle
     *
     * \param[in] cogIndex  The COG requested, should be; 0 <= \p cogIndex < cogs_.size()
     */
    bool cogIsFillerParticle(int cogIndex) const { return numAtomsPerCog_[cogIndex] == 0; }

    /*! \brief Returns the COG index given an atom index
     *
     * \param[in] atomIndex  The local atom index that maps to the COG
     */
    int cogIndex(int atomIndex) const
    {
        GMX_ASSERT(atomIndex >= 0 && static_cast<size_t>(atomIndex) < cogIndices_.size(),
                   "atomIndex should be in the range set in this object");

        return cogIndices_[atomIndex];
    }

    /*! \brief Return a reference to a center of geometry given an atom index
     *
     * \param[in] atomIndex  The atom index that maps to the COG requested
     */
    const RVec& cogForAtom(int atomIndex) const { return cogs_[cogIndex(atomIndex)]; }

    /*! \brief Clear the lists of stored center of geometry coordinates */
    void clear();

    /*! \brief Returns the maximum radius over all update groups */
    real maxUpdateGroupRadius() const { return maxUpdateGroupRadius_; }

private:
    /*! \brief Like addCogs(), but performs the work for one OpenMP thread
     *
     * This method should be called from within an OpenMP region.
     *
     * \param[in] globalAtomIndices  The global atom indices
     * \param[in] coordinates        The coordinates
     * \param[in] cogBegin           The start of the COG range to assign for all threads
     * \param[in] numThreads         The number of threads used
     * \param[in] thread             The index of the calling thread
     * \param[in] threadAtomRange    The range range of atom to process on this thread
     *
     * \returns the end of the COG range assigned by this thread
     */
    int addCogsThread(ArrayRef<const int>  globalAtomIndices,
                      ArrayRef<const RVec> coordinates,
                      int                  cogBegin,
                      int                  numThreads,
                      int                  thread,
                      const Range<int>&    threadAtomRange);

    /*! \libinternal
     * \brief Helper struct for mapping atom indices to COG indices, used per molecule block in mtop
     */
    struct IndexToGroup
    {
        //! Starting index of groups for this molblock
        int groupStart_;
        //! The number of atoms in a molecule
        int numGroupsPerMolecule_;
        //! Map atom indices to group indices for one molecule
        std::vector<int> groupIndex_;
    };

    //! Thread-local working data for building the COGs
    struct ThreadData
    {
        //! The number of update groups on this thread
        int numUpdateGroups_ = 0;
        //! Maps global COG index to local COG index
        HashedMap<int> globalToLocalMap_;
    };

    //! Maps atom indices to COG indices
    FastVector<int> cogIndices_;
    //! List of COGs
    FastVector<RVec> cogs_;
    //! List of the number of atoms for each COG, for filler particles the value is zero
    FastVector<int> numAtomsPerCog_;
    //! Helper data for mapping atom indices to COG indices
    std::vector<IndexToGroup> indicesPerMoleculeblock_;
    //! The maximum radius over all update groups
    real maxUpdateGroupRadius_;
    //! Reference to mtop this object was constructed with
    const gmx_mtop_t& mtop_;
    //! Lists for each molecule type that tell whether an atom is the first in the update group
    std::vector<std::vector<bool>> isFirstAtomInUpdateGroup_;
    //! Thread local data for building COGs
    std::vector<ThreadData> threadData_;
};

} // namespace gmx

#endif // GMX_MDLIB_UPDATEGROUPSCOG
