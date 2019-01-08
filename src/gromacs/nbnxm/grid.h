/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

/*! \internal \file
 *
 * \brief
 * Declares the Grid class.
 *
 * This class provides functionality for setting up and accessing atoms
 * on a grid for one domain decomposition zone. This grid is use for
 * generating cluster pair lists for computing non-bonded pair interactions.
 * The grid consists of a regular array of columns along dimensions x and y.
 * Within each column the cells are irregular along dimension z.
 * Each cell can hold one or more clusters of atoms, depending on the grid
 * geometry, which in turn is set by the non-bonded kernel type.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GRID_H
#define GMX_NBNXM_GRID_H

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"


struct gmx_domdec_zones_t;
struct nbnxn_atomdata_t;
struct nbnxn_search;
enum class PairlistType;

namespace gmx
{
class UpdateGroupsCog;
}

namespace Nbnxn
{

#ifndef DOXYGEN

/* Pair search box lower and upper corner in x,y,z.
 * Store this in 4 iso 3 reals, which is useful with 4-wide SIMD.
 * To avoid complicating the code we also use 4 without 4-wide SIMD.
 */
#define NNBSBB_C         4
/* Pair search box lower and upper bound in z only. */
#define NNBSBB_D         2
/* Pair search box lower and upper corner x,y,z indices, entry 3 is unused */
#define BB_X  0
#define BB_Y  1
#define BB_Z  2

#endif // !DOXYGEN


/* Bounding box for a nbnxn atom cluster */
struct BoundingBox
{
    float lower[NNBSBB_C];
    float upper[NNBSBB_C];
};


#ifndef DOXYGEN

/* Bounding box calculations are (currently) always in single precision, so
 * we only need to check for single precision support here.
 * This uses less (cache-)memory and SIMD is faster, at least on x86.
 */
#if GMX_SIMD4_HAVE_FLOAT
#    define NBNXN_SEARCH_BB_SIMD4      1
/* Memory alignment in bytes as required by SIMD aligned loads/stores */
#    define NBNXN_SEARCH_BB_MEM_ALIGN  (GMX_SIMD4_WIDTH*sizeof(float))
#else
#    define NBNXN_SEARCH_BB_SIMD4      0
/* No alignment required, but set it so we can call the same routines */
#    define NBNXN_SEARCH_BB_MEM_ALIGN  32
#endif


#if NBNXN_SEARCH_BB_SIMD4
/* Always use 4-wide SIMD for bounding box calculations */

#    if !GMX_DOUBLE
/* Single precision BBs + coordinates, we can also load coordinates with SIMD */
#        define NBNXN_SEARCH_SIMD4_FLOAT_X_BB  1
#    else
#        define NBNXN_SEARCH_SIMD4_FLOAT_X_BB  0
#    endif

/* The packed bounding box coordinate stride is always set to 4.
 * With AVX we could use 8, but that turns out not to be faster.
 */
#    define STRIDE_PBB       GMX_SIMD4_WIDTH
#    define STRIDE_PBB_2LOG  2

/* Store bounding boxes corners as quadruplets: xxxxyyyyzzzz */
#    define NBNXN_BBXXXX  1
/* Size of bounding box corners quadruplet */
#    define NNBSBB_XXXX  (NNBSBB_D*DIM*STRIDE_PBB)

#else  /* NBNXN_SEARCH_BB_SIMD4 */

#    define NBNXN_SEARCH_SIMD4_FLOAT_X_BB  0
#    define NBNXN_BBXXXX                   0

#endif /* NBNXN_SEARCH_BB_SIMD4 */

#endif // !DOXYGEN


/*! \internal
 * \brief A pair-search grid struct for one domain decomposition zone
 *
 * Note that when atom groups, instead of individual atoms, are assigned
 * to grid cells, individual atoms can be geometrically outside the cell
 * and grid that they have been assigned to (as determined by the center
 * or geometry of the atom group they belong to).
 */
class Grid
{
    public:
        /*! \internal
         * \brief The cluster and cell geometry of a grid
         */
        struct Geometry
        {
            //! Constructs the cluster/cell geometry given the type of pairlist
            Geometry(PairlistType pairlistType);

            bool isSimple;             //!< Is this grid simple (CPU) or hierarchical (GPU)
            int  numAtomsICluster;     //!< Number of atoms per cluster
            int  numAtomsJCluster;     //!< Number of atoms for list j-clusters
            int  numAtomsPerCell;      //!< Number of atoms per super-cluster
            int  numAtomsICluster2Log; //!< 2log of na_c
        };

        // The physical dimensions of a grid
        struct Dimensions
        {
            //! The lower corner of the (local) grid
            rvec lowerCorner;
            //! The upper corner of the (local) grid
            rvec upperCorner;
            //! The physical grid size: upperCorner - lowerCorner
            rvec gridSize;
            //! The atom number density for the (local) grid
            real atomDensity;
            //! The maximum distance an atom can be outside of a cell and outside of the grid
            real maxAtomGroupRadius;
            //! Size of cell along dimension x and y
            real cellSize[DIM - 1];
            //! 1/size of a cell along dimensions x and y
            real invCellSize[DIM - 1];
            //! The number of grid cells along dimensions x and y
            int  numCells[DIM - 1];
        };

        
        //! Constructs a grid given the type of pairlist
        Grid(PairlistType pairlistType);

        //! Returns the geometry of the grid cells
        const Geometry &geometry() const
        {
            return geometry_;
        }

        //! Returns the dimensions of the grid
        const Dimensions &dimensions() const
        {
            return dimensions_;
        }

        //! Returns the total number of grid columns
        int numColumns() const
        {
            return dimensions_.numCells[XX]*dimensions_.numCells[YY];
        }

        //! Returns the total number of grid cells
        int numCells() const
        {
            return numCellsTotal_;
        }

        //! Returns the cell offset of (the first cell of) this grid in the list of cells combined over all grids
        int cellOffset() const
        {
            return cellOffset_;
        }

        //! Returns the first cell index in the grid, starting at 0 in this grid
        int firstCellInColumn(int columnIndex) const
        {
            return cxy_ind_[columnIndex];
        };

        //! Returns the number of cells in the column
        int numCellsInColumn(int columnIndex) const
        {
            return cxy_ind_[columnIndex + 1] - cxy_ind_[columnIndex];
        };

        //! Returns the index of the first atom in the column
        int firstAtomInColumn(int columnIndex) const
        {
            return (cellOffset_ + cxy_ind_[columnIndex])*geometry_.numAtomsPerCell;
        };

        //! Returns the number of real atoms in the column
        int numAtomsInColumn(int columnIndex) const
        {
            return cxy_na_[columnIndex];
        };

        //! Returns the number of atoms in the column including padding
        int paddedNumAtomsInColumn(int columnIndex) const
        {
            return numCellsInColumn(columnIndex)*geometry_.numAtomsPerCell;
        };

        //! Returns the end of the atom index range on the grid, including padding
        int atomIndexEnd() const
        {
            return (cellOffset_ + numCellsTotal_)*geometry_.numAtomsPerCell;
        }

        //! Returns whether any atom in the cluster is perturbed
        bool clusterIsPerturbed(int clusterIndex) const
        {
            return fep_[clusterIndex] != 0u;
        }

        //! Returns whether the given atom in the cluster is perturbed
        bool atomIsPerturbed(int clusterIndex,
                             int atomIndexInCluster) const
        {
            return (fep_[clusterIndex] & (1 << atomIndexInCluster)) != 0u;
        }

        //! Returns the free-energy perturbation bits for the cluster
        unsigned int fepBits(int clusterIndex) const
        {
            return fep_[clusterIndex];
        }

        //! Returns the i-bounding boxes for all clusters on the grid
        gmx::ArrayRef<const BoundingBox> boundingBoxes() const
        {
            return bb_;
        }

        //! Returns the j-bounding boxes for all clusters on the grid
        gmx::ArrayRef<const BoundingBox> jBoundingBoxes() const
        {
            return bbj_;
        }

        //! Returns the packed bounding boxes for all clusters on the grid, empty with a CPU list
        gmx::ArrayRef<const float> packedBoundingBoxes() const
        {
            return pbb_;
        }

        //! Returns the packed bounding boxes for all cells on the grid
        gmx::ArrayRef<const float> zBoundingBoxes() const
        {
            return bbcz_;
        }

        //! Returns the flags for all clusters on the grid
        gmx::ArrayRef<const int> clusterFlags() const
        {
            return flags_;
        }

        //! Returns the number of clusters for all cells on the grid, empty with CPU a list
        gmx::ArrayRef<const int> numClustersPerCell() const
        {
            return numClusters_;
        }

        //! Returns the cluster index for an atom
        int atomToCluster(int atomIndex) const
        {
            return (atomIndex >> geometry_.numAtomsICluster2Log);
        }

        //! Returns the total number of clusters on the grid
        int numClusters() const
        {
            if (geometry_.isSimple)
            {
                return numCellsTotal_;
            }
            else
            {
                return numClustersTotal_;
            }
        }

        //! Sets tge grid dimensions
        void setDimensions(const nbnxn_search   *nbs,
                           int                   ddZone,
                           int                   numAtoms,
                           const rvec            lowerCorner,
                           const rvec            upperCorner,
                           real                  atomDensity,
                           real                  maxAtomGroupRadius);

        //! Determine in which grid cells the atoms should go
        void calcCellIndices(nbnxn_search                   *nbs,
                             int                             ddZone,
                             int                             cellOffset,
                             const gmx::UpdateGroupsCog     *updateGroupsCog,
                             int                             atomStart,
                             int                             atomEnd,
                             const int                      *atinfo,
                             gmx::ArrayRef<const gmx::RVec>  x,
                             int                             numAtomsMoved,
                             const int                      *move,
                             nbnxn_atomdata_t               *nbat);

    private:
        /*! \brief Fill a pair search cell with atoms
         *
         * Potentially sorts atoms and sets the interaction flags.
         */
        void fillCell(nbnxn_search                   *nbs,
                      nbnxn_atomdata_t               *nbat,
                      int                             atomStart,
                      int                             atomEnd,
                      const int                      *atinfo,
                      gmx::ArrayRef<const gmx::RVec>  x,
                      BoundingBox gmx_unused         *bb_work_aligned);

        //! Spatially sort the atoms within one grid column
        void sortColumnsSimple(nbnxn_search *nbs,
                               int dd_zone,
                               int atomStart, int atomEnd,
                               const int *atinfo,
                               gmx::ArrayRef<const gmx::RVec> x,
                               nbnxn_atomdata_t *nbat,
                               int cxy_start, int cxy_end,
                               gmx::ArrayRef<int> sort_work);

        //! Spatially sort the atoms within one grid column
        void sortColumnsSuperSub(nbnxn_search *nbs,
                                 int dd_zone,
                                 int atomStart, int atomEnd,
                                 const int *atinfo,
                                 gmx::ArrayRef<const gmx::RVec> x,
                                 nbnxn_atomdata_t *nbat,
                                 int cxy_start, int cxy_end,
                                 gmx::ArrayRef<int> sort_work);

        /* Data members */
        //! The geometry of the grid clusters and cells
        Geometry   geometry_;
        //! The physical dimensions of the grid
        Dimensions dimensions_;

        //! The total number of cells in this grid
        int        numCellsTotal_;
        //! Index in nbs->cell corresponding to cell 0  */
        int        cellOffset_;

        /* Grid data */
        //! The number of, non-filler, atoms for each grid column
        std::vector<int> cxy_na_;
        //! The grid-local cell index for each grid column
        std::vector<int> cxy_ind_;

        //! The number of cluster for each cell
        std::vector<int> numClusters_;

        /* Bounding boxes */
        //! Bounding boxes in z for the cells
        std::vector<float>                                  bbcz_;
        //! 3D bounding boxes for the sub cells
        std::vector < BoundingBox, gmx::AlignedAllocator < BoundingBox>> bb_;
        //! 3D j-bounding boxes for the case where the i- and j-cluster sizes are different
        std::vector < BoundingBox, gmx::AlignedAllocator < BoundingBox>> bbjStorage_;
        //! 3D j-bounding boxes
        gmx::ArrayRef<BoundingBox>                           bbj_;
        //! 3D bounding boxes in packed xxxx format per cell
        std::vector < float, gmx::AlignedAllocator < float>> pbb_;

        /* Bit-flag information */
        //! Flags for properties of clusters in each cell
        std::vector<int>          flags_;
        //! Signal bits for atoms in each cell that tell whether an atom is perturbed
        std::vector<unsigned int> fep_;

        /* Statistics */
        //! Total number of clusters, used for printing
        int numClustersTotal_;
};

} // namespace Nbnxn

#endif
