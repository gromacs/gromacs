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

/*! \internal \file
 *
 * \brief
 * Declares the Grid class.
 *
 * This class provides functionality for setting up and accessing atoms
 * on a grid for one domain decomposition zone. This grid is used for
 * generating cluster pair lists for computing non-bonded pair interactions.
 * The grid consists of a regular array of columns along dimensions x and y.
 * Along z the number of cells and their boundaries vary between the columns.
 * Each cell can hold one or more clusters of atoms, depending on the grid
 * geometry, which is set by the pair-list type.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GRID_H
#define GMX_NBNXM_GRID_H

#include <cstdint>

#include <memory>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"

#include "boundingbox.h"

namespace gmx
{

struct nbnxn_atomdata_t;
enum class PairlistType;
class UpdateGroupsCog;
struct GridSetData;
struct GridWork;

/*! \internal
 * \brief Bounding box for one dimension of a grid cell
 */
struct BoundingBox1D
{
    //! lower bound
    float lower;
    //! upper bound
    float upper;
};

/*! \internal
 * \brief A pair-search grid object for one domain decomposition zone
 *
 * This is a rectangular 3D grid covering a potentially non-rectangular
 * volume which is either the whole unit cell or the local zone or part
 * of a non-local zone when using domain decomposition. The grid cells
 * are even spaced along x/y and irregular along z. Each cell is sub-divided
 * into atom clusters. With a CPU geometry, each cell contains 1 or 2 clusters.
 * With a GPU geometry, each cell contains up to 8 clusters. The geometry is
 * set by the pairlist type which is the only argument of the constructor.
 *
 * When multiple grids are used, i.e. with domain decomposition, we want
 * to avoid the overhead of multiple coordinate arrays or extra indexing.
 * Therefore each grid stores a cell offset, so a contiguous cell index
 * can be used to index atom arrays. All methods returning atom indices
 * return indices which index into a full atom array.
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

        //! Is this grid simple (CPU) or hierarchical (GPU)
        bool isSimple;
        //! Number of atoms per cluster
        int numAtomsICluster;
        //! Number of atoms for list j-clusters
        int numAtomsJCluster;
        //! Number of atoms per cell
        int numAtomsPerCell;
        //! 2log of na_c
        int numAtomsICluster2Log;
    };

    //! The physical dimensions of a grid \internal
    struct Dimensions
    {
        //! The lower corner of the (local) grid
        RVec lowerCorner;
        //! The upper corner of the (local) grid
        RVec upperCorner;
        //! The physical grid size: upperCorner - lowerCorner
        RVec gridSize;
        //! An estimate for the atom number density of the region targeted by the grid
        real atomDensity;
        //! The maximum distance an atom can be outside of a cell and outside of the grid
        real maxAtomGroupRadius;
        //! Size of cell along dimension x and y
        real cellSize[DIM - 1];
        //! 1/size of a cell along dimensions x and y
        real invCellSize[DIM - 1];
        //! The number of grid cells along dimensions x and y
        int numCells[DIM - 1];
    };

    //! Constructs a grid given the type of pairlist
    Grid(PairlistType pairlistType, const bool& haveFep, PinningPolicy pinningPolicy);

    //! Returns the geometry of the grid cells
    const Geometry& geometry() const { return geometry_; }

    //! Returns the dimensions of the grid
    const Dimensions& dimensions() const { return dimensions_; }

    //! Returns the total number of grid columns
    int numColumns() const { return dimensions_.numCells[XX] * dimensions_.numCells[YY]; }

    //! Returns the total number of grid cells
    int numCells() const { return numCellsTotal_; }

    //! Returns the cell offset of (the first cell of) this grid in the list of cells combined over all grids
    int cellOffset() const { return cellOffset_; }

    //! Returns the maximum number of grid cells in a column
    int numCellsColumnMax() const { return numCellsColumnMax_; }

    //! Returns the start of the source atom range mapped to this grid
    int srcAtomBegin() const { return srcAtomBegin_; }

    //! Returns the end of the source atom range mapped to this grid
    int srcAtomEnd() const { return srcAtomEnd_; }

    //! Returns the first cell index in the grid, starting at 0 in this grid
    int firstCellInColumn(int columnIndex) const { return cxy_ind_[columnIndex]; }

    //! Returns the number of cells in the column
    int numCellsInColumn(int columnIndex) const
    {
        return cxy_ind_[columnIndex + 1LL] - cxy_ind_[columnIndex];
    }

    //! Returns the index of the first atom in the column
    int firstAtomInColumn(int columnIndex) const
    {
        return (cellOffset_ + cxy_ind_[columnIndex]) * geometry_.numAtomsPerCell;
    }

    //! Returns the number of real atoms in the column
    int numAtomsInColumn(int columnIndex) const { return cxy_na_[columnIndex]; }

    /*! \brief Returns a view of the number of non-filler, atoms for each grid column
     *
     * \todo Needs a useful name. */
    ArrayRef<const int> cxy_na() const { return cxy_na_; }
    /*! \brief Returns a view of the grid-local cell index for each grid column
     *
     * \todo Needs a useful name. */
    ArrayRef<const int> cxy_ind() const { return cxy_ind_; }

    //! Returns the number of real atoms in the column
    int numAtomsPerCell() const { return geometry_.numAtomsPerCell; }

    //! Returns the number of atoms in the column including padding
    int paddedNumAtomsInColumn(int columnIndex) const
    {
        return numCellsInColumn(columnIndex) * geometry_.numAtomsPerCell;
    }

    //! Returns the end of the atom index range on the grid, including padding
    int atomIndexEnd() const { return (cellOffset_ + numCellsTotal_) * geometry_.numAtomsPerCell; }

    //! Returns whether any atom in the cluster is perturbed
    bool clusterIsPerturbed(int clusterIndex) const { return fep_[clusterIndex] != 0U; }

    //! Returns whether the given atom in the cluster is perturbed
    bool atomIsPerturbed(int clusterIndex, int atomIndexInCluster) const
    {
        return (fep_[clusterIndex] & (1 << atomIndexInCluster)) != 0U;
    }

    //! Returns the free-energy perturbation bits for the cluster
    unsigned int fepBits(int clusterIndex) const { return fep_[clusterIndex]; }

    //! Returns the i-bounding boxes for all clusters on the grid
    ArrayRef<const BoundingBox> iBoundingBoxes() const { return bb_; }

    //! Returns the j-bounding boxes for all clusters on the grid
    ArrayRef<const BoundingBox> jBoundingBoxes() const { return bbj_; }

    //! Returns the packed bounding boxes for all clusters on the grid, empty with a CPU list
    ArrayRef<const float> packedBoundingBoxes() const { return pbb_; }

    //! Returns the bounding boxes along z for all cells on the grid
    ArrayRef<const BoundingBox1D> zBoundingBoxes() const { return bbcz_; }

    //! Returns the flags for all clusters on the grid
    ArrayRef<const int> clusterFlags() const { return flags_; }

    //! Returns the number of clusters for all cells on the grid, empty with a CPU geometry
    ArrayRef<const int> numClustersPerCell() const { return numClusters_; }

    //! Returns the cluster index for an atom
    int atomToCluster(int atomIndex) const { return (atomIndex >> geometry_.numAtomsICluster2Log); }

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

    /*! \brief Sets the grid dimensions
     *
     * \param[in] ddZone           The domain decomposition zone index
     * \param[in] numAtoms         The number of atoms to put onto this grid
     * \param[in] lowerCorner      The minimum Cartesian coordinates of the grid
     * \param[in] upperCorner      The maximum Cartesian coordinates of the grid
     * \param[in,out] atomDensity  The atom density, will be computed when <= 0
     * \param[in] maxAtomGroupRadius  The maximum radius of atom groups
     */
    void setDimensions(int         ddZone,
                       int         numAtoms,
                       const RVec& lowerCorner,
                       const RVec& upperCorner,
                       real*       atomDensity,
                       real        maxAtomGroupRadius);

    //! Sets the cell indices using indices in \p gridSetData and \p gridWork
    void setCellIndices(int                     ddZone,
                        int                     cellOffset,
                        GridSetData*            gridSetData,
                        ArrayRef<GridWork>      gridWork,
                        Range<int>              atomRange,
                        int                     numGridAtomsWithoutFillers,
                        ArrayRef<const int32_t> atomInfo,
                        ArrayRef<const RVec>    x,
                        nbnxn_atomdata_t*       nbat);

    //! Determine in which grid columns atoms should go, store cells and atom counts in \p cell and \p cxy_na
    static void calcColumnIndices(const Grid::Dimensions& gridDims,
                                  const UpdateGroupsCog*  updateGroupsCog,
                                  Range<int>              atomRange,
                                  ArrayRef<const RVec>    x,
                                  int                     dd_zone,
                                  const int*              move,
                                  int                     thread,
                                  int                     nthread,
                                  ArrayRef<int>           cell,
                                  ArrayRef<int>           cxy_na);

private:
    /*! \brief Fill a pair search cell with atoms
     *
     * Potentially sorts atoms and sets the interaction flags.
     */
    void fillCell(GridSetData*            gridSetData,
                  nbnxn_atomdata_t*       nbat,
                  int                     atomStart,
                  int                     atomEnd,
                  ArrayRef<const int32_t> atomInfo,
                  ArrayRef<const RVec>    x);

    //! Spatially sort the atoms within the given column range, for CPU geometry
    void sortColumnsCpuGeometry(GridSetData*            gridSetData,
                                int                     dd_zone,
                                ArrayRef<const int32_t> atomInfo,
                                ArrayRef<const RVec>    x,
                                nbnxn_atomdata_t*       nbat,
                                Range<int>              columnRange,
                                ArrayRef<int>           sort_work);

    //! Spatially sort the atoms within the given column range, for GPU geometry
    void sortColumnsGpuGeometry(GridSetData*            gridSetData,
                                int                     dd_zone,
                                ArrayRef<const int32_t> atomInfo,
                                ArrayRef<const RVec>    x,
                                nbnxn_atomdata_t*       nbat,
                                Range<int>              columnRange,
                                ArrayRef<int>           sort_work);

    /* Data members */
    //! The geometry of the grid clusters and cells
    Geometry geometry_;
    //! The physical dimensions of the grid
    Dimensions dimensions_;

    //! The total number of cells in this grid
    int numCellsTotal_;
    //! Index in nbs->cell corresponding to cell 0
    int cellOffset_;
    //! The maximum number of cells in a column
    int numCellsColumnMax_;

    //! The start of the source atom range mapped to this grid
    int srcAtomBegin_;
    //! The end of the source atom range mapped to this grid
    int srcAtomEnd_;

    /* Grid data */
    /*! \brief The number of, non-filler, atoms for each grid column.
     *
     * \todo Needs a useful name. */
    HostVector<int> cxy_na_;
    /*! \brief The grid-local cell index for each grid column
     *
     * \todo Needs a useful name. */
    HostVector<int> cxy_ind_;

    //! The number of cluster for each cell
    std::vector<int> numClusters_;

    /* Bounding boxes */
    //! Bounding boxes in z for the cells
    std::vector<BoundingBox1D> bbcz_;
    //! 3D bounding boxes for the sub cells
    std::vector<BoundingBox, AlignedAllocator<BoundingBox>> bb_;
    //! 3D j-bounding boxes for the case where the i- and j-cluster sizes are different
    std::vector<BoundingBox, AlignedAllocator<BoundingBox>> bbjStorage_;
    //! 3D j-bounding boxes
    ArrayRef<BoundingBox> bbj_;
    //! 3D bounding boxes in packed xxxx format per cell
    std::vector<float, AlignedAllocator<float>> pbb_;

    //! Tells whether we have perturbed interactions, authorative source is in GridSet (never modified)
    const bool& haveFep_;

    /* Bit-flag information */
    //! Flags for properties of clusters in each cell
    std::vector<int> flags_;
    //! Signal bits for atoms in each cell that tell whether an atom is perturbed
    std::vector<unsigned int> fep_;

    /* Statistics */
    //! Total number of clusters, used for printing
    int numClustersTotal_;
};

/*! \brief Sets the 2D search grid dimensions puts the atoms on the 2D grid
 *
 * \param[in,out] grid      The pair search grid for one DD zone
 * \param[in,out] gridWork  Working data for each thread
 * \param[in,out] cells     The grid cell list
 * \param[in] lowerCorner   The minimum Cartesian coordinates of the grid
 * \param[in] upperCorner   The maximum Cartesian coordinates of the grid
 * \param[in] updateGroupsCog  The center of geometry of update groups, can be nullptr
 * \param[in] atomRange     The range of atoms to put on this grid, may include moved atoms
 * \param[in] numGridAtomsWithoutFillers  The number of non-filler atoms to put on this grid
 * \param[in,out] atomDensity  The atom density, will be computed when <= 0
 * \param[in] maxAtomGroupRadius  The maximum radius of atom groups
 * \param[in] x             The coordinates of the atoms
 * \param[in] ddZone        The domain decomposition zone
 * \param[in] move          Tells whether atoms have moved to another DD domain
 * \param[in] computeGridDensityRatio  When true, return the grid density ratio
 *
 * \returns When \p computeGridDensityRatio==true, the ratio of the effective 2D grid density and the uniform grid density
 */
real generateAndFill2DGrid(Grid*                  grid,
                           ArrayRef<GridWork>     gridWork,
                           HostVector<int>*       cells,
                           const rvec             lowerCorner,
                           const rvec             upperCorner,
                           const UpdateGroupsCog* updateGroupsCog,
                           Range<int>             atomRange,
                           int                    numGridAtomsWithoutFillers,
                           real*                  atomDensity,
                           real                   maxAtomGroupRadius,
                           ArrayRef<const RVec>   x,
                           int                    ddZone,
                           const int*             move,
                           bool                   computeGridDensityRatio);

} // namespace gmx

#endif
