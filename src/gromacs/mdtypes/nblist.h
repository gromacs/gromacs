/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * \brief
 * Declares atom pair list class t_nblist
 *
 * \author Berk Hess
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_NBLIST_H
#define GMX_MDTYPES_NBLIST_H

#include <vector>


/*! \brief A plain atom pair list, used only for perturbed non-bonded interactions
 *
 * A list of atom pairs stored as list of so-called j-atoms for each i-atom
 * in the list. The list is constructed by calling \p addIEntry() and specifying
 * a maximum number of j-entries that can be added. Then addJEntry() can be
 * called to add j-atoms. For performance it can be good to call
 * \p popIEntryWhenEmpty() after potentially adding j-atoms, to remove empty
 * i-entries.
 */
class t_nblist
{
public:
    //! An i-entry
    struct IEntry
    {
        //! The i-atom index
        int atom;
        //! The shift vector index
        int shiftIndex;
        //! The energy group pair index for the i-entry with all its j-entries
        int energyGroupPair;
    };

    //! A j-entry
    struct JEntry
    {
        //! The j-atom index
        int atom;
        //! Whether j interacts with i, i.e. this pair is not excluded
        bool interacts;
    };

    /*! \brief Adds a new i-entry with space for at most \p maxJEntries j-entries
     *
     * Note that this will replace the last entry instead when that is empty,
     * i.e. has an empty list of j-entries.
     */
    void addIEntry(const IEntry& iEntry, int maxJEntries)
    {
        // Is the last i-entry empty?
        if (iList_.empty() || jRanges_.back() > jRanges_[iList_.size() - 1])
        {
            iList_.push_back(iEntry);
        }
        else
        {
            // We have an empty i-entry at the end, overwrite it
            iList_.back() = iEntry;
        }
        const int numJEntries = jRanges_.back();
        jRanges_.push_back(numJEntries);
        if (numJEntries + maxJEntries > gmx::ssize(jList_))
        {
            jList_.resize(numJEntries + maxJEntries);
        }
    }

    //! Adds a j-entry for the last i-entry
    void addJEntry(const JEntry& jEntry)
    {
        GMX_ASSERT(jRanges_.back() < gmx::ssize(jList_),
                   "We should have reserved sufficient space for calling addJEntry()x times");

        jList_[jRanges_.back()] = jEntry;

        jRanges_.back()++;

        if (!jEntry.interacts)
        {
            numExclusionsWithinRlist_++;
        }
    }

    //! Checks whether the last i-entry has no j-entries and if so removes it
    void popIEntryWhenEmpty()
    {
        GMX_ASSERT(!iList_.empty(), "Can not prune with an empty list");

        if (jRanges_.back() == jRanges_[iList_.size() - 1])
        {
            iList_.pop_back();
            jRanges_.pop_back();
        }
    }

    //! Clears the pair list
    void clear()
    {
        iList_.clear();
        jRanges_.resize(1);

        numExclusionsWithinRlist_ = 0;
    }

    //! Returns the list of i-entries
    gmx::ArrayRef<const IEntry> iList() const { return iList_; }

    //! Returns the list of j-entries for i-entry \p iEntry
    gmx::ArrayRef<const JEntry> jList(gmx::Index iEntry) const
    {
        return gmx::constArrayRefFromArray(jList_.data() + jRanges_[iEntry],
                                           jRanges_[iEntry + 1] - jRanges_[iEntry]);
    }

    //! Returns a concatenated list of all j-entries for all i-entries
    gmx::ArrayRef<const JEntry> flatJList() const
    {
        return gmx::constArrayRefFromArray(jList_.data(), jRanges_.back());
    }

    //! Returns the number of excluded pairs within rlist at list creation
    int numExclusionsWithinRlist() const { return numExclusionsWithinRlist_; }

private:
    //! The list of i-entries
    std::vector<IEntry> iList_;
    //! The range of j-entries, i-entry with index i has \p jRanges_[i] to \p jRanges_[i+1]
    std::vector<int> jRanges_ = { 0 };
    //! Buffer for the j-entries; note that this buffer is only sized up, never down
    std::vector<JEntry> jList_;
    //! The number of exclusions at distance < rlist
    int numExclusionsWithinRlist_ = 0;
};

#endif /* GMX_MDTYPES_NBLIST_H */
