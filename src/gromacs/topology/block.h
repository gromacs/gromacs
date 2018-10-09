/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2018, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_BLOCK_H
#define GMX_TOPOLOGY_BLOCK_H

#include <stdio.h>

#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Division of a range of indices into consecutive blocks
 *
 * A range of consecutive indices 0 to full.range.end() is divided
 * into numBlocks() consecutive blocks of consecutive indices.
 * Block b contains indices i for which block(b).begin() <= i < block(b).end().
 */
class RangePartitioning
{
    public:
        /*! \brief Struct for returning the range of a block.
         *
         * Can be used in a range loop.
         */
        struct Block
        {
            public:
                /*! \brief An iterator that loops over integers */
                struct iterator
                {
                    //! Constructor
                    iterator(int value) : value_(value) {}
                    //! Value
                    operator int () const { return value_; }
                    //! Reference
                    operator int &()      { return value_; }
                    //! Pointer
                    int operator* () const { return value_; }
                    //! Inequality comparison
                    bool operator!= (const iterator other) { return value_ != other; }
                    //! Increment operator
                    iterator &operator++() { ++value_; return *this; }
                    //! Increment operator
                    iterator operator++(int) { iterator tmp(*this); ++value_; return tmp; }
                    //! The actual value
                    int value_;
                };

                /*! \brief Constructor, constructs a range starting at 0 with 0 blocks */
                Block(int begin,
                      int end) :
                    begin_(begin),
                    end_(end)
                {
                }

                /*! \brief Begin iterator/value */
                const iterator begin() const { return begin_; }
                /*! \brief End iterator/value */
                const iterator end() const { return end_; }

                /*! \brief The number of items in the block */
                int size() const
                {
                    return end_ - begin_;
                }

                /*! \brief Returns whether \p index is within range of the block */
                bool inRange(int index) const
                {
                    return (begin_ <= index && index < end_);
                }

            private:
                const int begin_; /**< The start index of the block */
                const int end_;   /**< The end index of the block */
        };

        /*! \brief Returns the number of blocks */
        int numBlocks() const
        {
            return index_.size() - 1;
        }

        /*! \brief Returns the size of the block with index \p blockIndex */
        Block block(int blockIndex) const
        {
            return Block(index_[blockIndex], index_[blockIndex + 1]);
        }

        /*! \brief Returns the full range */
        Block fullRange() const
        {
            return Block(index_.front(), index_.back());
        }

        /*! \brief Returns a range starting at \p blockIndexBegin and ending at \p blockIndexEnd */
        Block subRange(int blockIndexBegin,
                       int blockIndexEnd) const
        {
            return Block(index_[blockIndexBegin], index_[blockIndexEnd]);
        }

        /*! \brief Returns true when all blocks have size 0 or numBlocks()=0 */
        bool allBlocksHaveSizeOne() const
        {
            return (index_.back() == numBlocks());
        }

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
        void clear()
        {
            index_.resize(1);
        }

        /*! \brief Reduces the number of blocks to \p newNumBlocks
         *
         * \note \p newNumBlocks should be <= numBlocks().
         */
        void reduceNumBlocks(int newNumBlocks)
        {
            GMX_ASSERT(newNumBlocks <= numBlocks(), "Can only shrink to fewer blocks");
            index_.resize(newNumBlocks + 1);
        }

        /*! \brief Sets the partitioning to \p numBlocks blocks each of size 1 */
        void setAllBlocksSizeOne(int numBlocks);

        /*! \brief Returns the raw block index array, avoid using this */
        std::vector<int> &rawIndex()
        {
            return index_;
        }

    private:
        std::vector<int> index_ = { 0 }; /**< The list of block begin/end indices */
};

}  // namespace gmx

/* Deprecated, C-style version of RangePartitioning */
typedef struct t_block
{
    int blockSize(int blockIndex) const
    {
        GMX_ASSERT(blockIndex < nr, "blockIndex should be in range");
        return index[blockIndex + 1] - index[blockIndex];
    }

    int      nr;           /* The number of blocks          */
    int     *index;        /* Array of indices (dim: nr+1)  */
    int      nalloc_index; /* The allocation size for index */
} t_block;

struct t_blocka
{
    int      nr;    /* The number of blocks              */
    int     *index; /* Array of indices in a (dim: nr+1) */
    int      nra;   /* The number of atoms               */
    int     *a;     /* Array of atom numbers in each group  */
    /* (dim: nra)                           */
    /* Block i (0<=i<nr) runs from          */
    /* index[i] to index[i+1]-1. There will */
    /* allways be an extra entry in index   */
    /* to terminate the table               */
    int nalloc_index;           /* The allocation size for index        */
    int nalloc_a;               /* The allocation size for a            */
};

void init_block(t_block *block);
void init_blocka(t_blocka *block);
t_blocka *new_blocka();
/* allocate new block */

void done_block(t_block *block);
void done_blocka(t_blocka *block);

void copy_blocka(const t_blocka *src, t_blocka *dest);

void copy_block(const t_block *src, t_block *dst);

void stupid_fill_block(t_block *grp, int natom, gmx_bool bOneIndexGroup);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * If bOneIndexGroup, then all atoms are  lumped in one index group,
 * otherwise there is one atom per index entry
 */

void stupid_fill_blocka(t_blocka *grp, int natom);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * There is one atom per index entry
 */

void pr_block(FILE *fp, int indent, const char *title, const t_block *block, gmx_bool bShowNumbers);
void pr_blocka(FILE *fp, int indent, const char *title, const t_blocka *block, gmx_bool bShowNumbers);

#endif
