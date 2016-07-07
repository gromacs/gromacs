/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements classes in atomset.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include "atomset.h"

#include <algorithm>
#include <memory>

#include "gromacs/domdec/ga2la.h"
#include "gromacs/utility/classhelpers.h"


namespace gmx
{

/********************************************************************
 * AtomSet::Impl
 */

/*! \internal \brief
 * Private implementation class for AtomSet.
 */
class AtomSet::Impl
{
    public:
        Impl();
        std::vector<int> global_index_;     /**< Global indices of the atoms in this set, e.g., as in the index file.*/
        std::vector<int> collective_index_; /**< Maps indices on one node (0..num_atoms_local_) to global atom indicices. */
        std::vector<int> local_index_;      /**< Local indices of the external potential atoms, used to add to local force;
                                               set by set_indices. Access the coordinate of the i-th local atom of this set by x[local_index_[i]].
                                               Constructed and updated every domain-decomposition step by set_indices. */
};

AtomSet::Impl::Impl(){};

/********************************************************************
 * AtomSet
 */

std::unique_ptr<AtomSet> AtomSet::create()
{
    return std::unique_ptr<AtomSet>(new AtomSet());
};

void AtomSet::init(const int number_of_atoms, const int *index, bool bParallel)
{
    impl_->global_index_.assign(index, index+number_of_atoms);

    /* if not running in parallel, local atoms are global atoms and
     * the collective index runs from 0..number_of_atoms-1
     */
    if (!bParallel)
    {
        impl_->local_index_.assign(index, index+number_of_atoms);
        impl_->collective_index_.resize(number_of_atoms);
        std::iota(impl_->collective_index_.begin(), impl_->collective_index_.end(), 0);
    }
};

AtomSet::AtomSet() : impl_(PrivateImplPointer<Impl>(new Impl()))
{
}


AtomSet::~AtomSet()
{
}

void AtomSet::bparallel_set_local_and_collective_indices(const gmx_ga2la_t *ga2la)
{
    /* Loop over all the atom indices of the group to check
     * which ones are on the local node */

    int  i_local;
    int  nalloc_loc = 0;
    int  n_local    = 0;
    int  n_global   = impl_->global_index_.size();

    impl_->local_index_.clear();
    impl_->collective_index_.clear();

    for (int i_collective = 0; i_collective < n_global; i_collective++)
    {
        if (ga2la_get_home(ga2la, impl_->global_index_[i_collective], &i_local))
        {
            /* The atom with this index is a home atom ? */
            if (n_local >= nalloc_loc)  /* Check whether memory suffices */
            {
                nalloc_loc = over_alloc_dd(n_local+1);
                /* We never need more memory than the number of atoms in the group */
                nalloc_loc = std::min(nalloc_loc, n_local);

                impl_->local_index_.reserve(nalloc_loc);
                impl_->collective_index_.reserve(nalloc_loc);
            }
            /* Save the atoms index in the local atom numbers array */
            impl_->local_index_.push_back(i_local);

            /* Keep track of where this local atom belongs in the collective index array.
             * This is needed when reducing the local arrays to a collective/global array
             * in communicate_group_positions */
            impl_->collective_index_.push_back(i_collective);
            n_local++;
        }
    }
};


const std::vector<int> &AtomSet::global_index() const
{
    return impl_->global_index_;
}


const std::vector<int> &AtomSet::local_index() const
{
    return impl_->local_index_;
}

const std::vector<int> &AtomSet::collective_index() const
{
    return impl_->collective_index_;
}

int AtomSet::num_atoms_global() const
{
    return impl_->global_index_.size();
}

int AtomSet::num_atoms_local() const
{
    return impl_->local_index_.size();
};

std::vector<int>::const_iterator AtomSet::begin() const
{
    return local_index().begin();
};

std::vector<int>::const_iterator AtomSet::end() const
{
    return local_index().end();
};

} // namespace gmx
