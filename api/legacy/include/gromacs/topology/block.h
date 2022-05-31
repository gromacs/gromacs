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
#ifndef GMX_TOPOLOGY_BLOCK_H
#define GMX_TOPOLOGY_BLOCK_H

#include <stdio.h>

#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/range.h"

namespace gmx
{

template<typename>
class ListOfLists;

/*! \brief Division of a range of indices into consecutive blocks
 *
 * A range of consecutive indices 0 to full.range.end() is divided
 * into numBlocks() consecutive blocks of consecutive indices.
 * Block b contains indices i for which block(b).begin() <= i < block(b).end().
 */
class RangePartitioning
{
public:
    /*! \brief A block defined by a range of atom indices */
    using Block = Range<int>;

    /*! \brief Returns the number of blocks */
    int numBlocks() const { return static_cast<int>(index_.size()) - 1; }

    /*! \brief Returns the size of the block with index \p blockIndex */
    Block block(int blockIndex) const
    {
        return Block(index_[blockIndex], index_[blockIndex + 1LL]);
    }

    /*! \brief Returns the full range */
    Block fullRange() const { return Block(index_.front(), index_.back()); }

    /*! \brief Returns a range starting at \p blockIndexBegin and ending at \p blockIndexEnd */
    Block subRange(int blockIndexBegin, int blockIndexEnd) const
    {
        return Block(index_[blockIndexBegin], index_[blockIndexEnd]);
    }

    /*! \brief Returns true when all blocks have size 0 or numBlocks()=0 */
    bool allBlocksHaveSizeOne() const { return (index_.back() == numBlocks()); }

    /*! \brief Appends a block of size \p blockSize at the end of the range
     *
     * \note blocksize has to be >= 1
     */
    void appendBlock(int blockSize)
    {
        GMX_ASSERT(blockSize > 0, "block sizes should be >= 1");
        index_.push_back(index_.back() + blockSize);
    }

    /*! \brief Removes all blocks */
    void clear() { index_.resize(1); }

    /*! \brief Reduces the number of blocks to \p newNumBlocks
     *
     * \note \p newNumBlocks should be <= numBlocks().
     */
    void reduceNumBlocks(int newNumBlocks)
    {
        GMX_ASSERT(newNumBlocks <= numBlocks(), "Can only shrink to fewer blocks");
        index_.resize(newNumBlocks + 1LL);
    }

    /*! \brief Sets the partitioning to \p numBlocks blocks each of size 1 */
    void setAllBlocksSizeOne(int numBlocks);

    /*! \brief Returns the raw block index array, avoid using this */
    std::vector<int>& rawIndex() { return index_; }

private:
    std::vector<int> index_ = { 0 }; /**< The list of block begin/end indices */
};

} // namespace gmx

/* Deprecated, C-style version of RangePartitioning */
typedef struct t_block
{
    int blockSize(int blockIndex) const
    {
        GMX_ASSERT(blockIndex < nr, "blockIndex should be in range");
        return index[blockIndex + 1] - index[blockIndex];
    }

    int  nr;           /* The number of blocks          */
    int* index;        /* Array of indices (dim: nr+1)  */
    int  nalloc_index; /* The allocation size for index */
} t_block;

struct t_blocka
{
    int  nr;    /* The number of blocks              */
    int* index; /* Array of indices in a (dim: nr+1) */
    int  nra;   /* The number of atoms               */
    int* a;     /* Array of atom numbers in each group  */
    /* (dim: nra)                           */
    /* Block i (0<=i<nr) runs from          */
    /* index[i] to index[i+1]-1. There will */
    /* allways be an extra entry in index   */
    /* to terminate the table               */
    int nalloc_index; /* The allocation size for index        */
    int nalloc_a;     /* The allocation size for a            */
};

/*! \brief
 * Fully initialize t_block datastructure.
 *
 * Initializes a \p block and sets up the first index to zero.
 *
 * \param[in,out] block datastructure to initialize.
 */
void init_block(t_block* block);

/*! \brief
 * Fully initialize t_blocka datastructure.
 *
 * Initializes a \p block and sets up the first index to zero.
 * The atom number array is initialized to nullptr.
 *
 * \param[in,out] block datastructure to initialize.
 */
void init_blocka(t_blocka* block);

/* TODO
 * In general all t_block datastructures should be avoided
 * in favour of RangePartitioning. This here is a simple cludge
 * to use more modern initialization while we move to the use
 * of RangePartitioning.
 */

/*! \brief
 * Minimal initialization of t_block datastructure.
 *
 * Performs the equivalent to a snew on a t_block, setting all
 * values to zero or nullptr. Needed for some cases where the topology
 * handling expects a block to be valid initialized (e.g. during domain
 * decomposition) but without the first block set to zero.
 *
 * \param[in,out] block datastructure to initialize.
 */
void init_block_null(t_block* block);

/*! \brief
 * Minimal initialization of t_blocka datastructure.
 *
 * Performs the equivalent to a snew on a t_blocka, setting all
 * values to zero or nullptr. Needed for some cases where the topology
 * handling expects a block to be valid initialized (e.g. during domain
 * decomposition) but without the first block set to zero.
 *
 * \param[in,out] block datastructure to initialize.
 */
void init_blocka_null(t_blocka* block);

t_blocka* new_blocka();
/* allocate new block */

void done_block(t_block* block);
//! Deallocates memory within \c block
void done_blocka(t_blocka* block);

void copy_blocka(const t_blocka* src, t_blocka* dest);

void copy_block(const t_block* src, t_block* dst);

void stupid_fill_block(t_block* grp, int natom, gmx_bool bOneIndexGroup);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * If bOneIndexGroup, then all atoms are  lumped in one index group,
 * otherwise there is one atom per index entry
 */

void stupid_fill_blocka(t_blocka* grp, int natom);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * There is one atom per index entry
 */

void pr_block(FILE* fp, int indent, const char* title, const t_block* block, gmx_bool bShowNumbers);
void pr_blocka(FILE* fp, int indent, const char* title, const t_blocka* block, gmx_bool bShowNumbers);
void pr_listoflists(FILE* fp, int indent, const char* title, const gmx::ListOfLists<int>* block, gmx_bool bShowNumbers);

#endif
