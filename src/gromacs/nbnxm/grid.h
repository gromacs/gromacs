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
 * Along z the number of bins and their boundaries vary between the columns.
 * Each bin can hold one or more clusters of atoms, depending on the grid
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
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

#include "boundingbox.h"

namespace gmx
{

struct nbnxn_atomdata_t;
enum class PairlistType;
class UpdateGroupsCog;
struct GridSetData;
struct GridWork;

/*! \internal
 * \brief Bounding box for one dimension of a grid bin
 */
struct BoundingBox1D
{
    //! lower bound
    float lower;
    //! upper bound
    float upper;
};

//! The physical dimensions of a grid \internal
struct GridDimensions
{
    //! Returns the lower corner along dimension \p dim of the cell with index \p binIndex
    real cellLowerCorner(int dim, int cellIndex) const
    {
        return lowerCorner[dim] + cellIndex * cellSize[dim];
    }

    //! Return the index of the column on the grid given the x+y-indices
    int columnIndex(int columnIndexX, int columnIndexY) const
    {
        return columnIndexX * numCells[YY] + columnIndexY;
    }

    //! The lower corner of the (local) grid
    RVec lowerCorner;
    //! The upper corner of the (local) grid
    RVec upperCorner;
    //! The physical grid size: upperCorner - lowerCorner
    RVec gridSize;
    //! An estimate for the atom number density of the region targeted by the grid
    real atomDensity;
    //! The maximum distance an atom can be outside of a cell/bin and outside of the grid
    real maxAtomGroupRadius;
    //! Size of cell along dimension x and y
    real cellSize[DIM - 1];
    //! 1/size of a cell along dimensions x and y
    real invCellSize[DIM - 1];
    //! The number of grid cells along dimensions x and y
    int numCells[DIM - 1];
};

/*! \internal
 * \brief A pair-search grid object for one domain decomposition zone
 *
 * This is a rectangular 3D grid covering a potentially non-rectangular
 * volume which is either the whole unit cell or the local zone or part
 * of a non-local zone when using domain decomposition. An uniformly space
 * grid is used along dimensions x and y. Along dimension z atoms are
 * organized in bins with equal atom count. Each bin is sub-divided
 * into atom clusters. With a CPU geometry, each bin contains 1 or 2 clusters.
 * With a GPU geometry, each bin contains 8 clusters. The geometry is
 * set by the pairlist type which is the only argument of the constructor.
 *
 * When multiple grids are used, i.e. with domain decomposition, we want
 * to avoid the overhead of multiple coordinate arrays or extra indexing.
 * Therefore each grid stores a bin offset, so a contiguous bin index
 * can be used to index atom arrays. All methods returning atom indices
 * return indices which index into a full atom array.
 *
 * Note that when atom groups, instead of individual atoms, are assigned
 * to grid bins, individual atoms can be geometrically outside the bin
 * and grid that they have been assigned to (as determined by the center
 * or geometry of the atom group they belong to).
 */
class Grid
{
public:
    /*! \internal
     * \brief The cluster and bin geometry of a grid
     */
    struct Geometry
    {
        //! Constructs the cluster+bin geometry given the type of pairlist
        Geometry(PairlistType pairlistType);

        //! Is this grid simple (CPU) or hierarchical (GPU)
        bool isSimple_;
        //! Number of atoms per cluster
        int numAtomsICluster_;
        //! Number of atoms for list j-clusters
        int numAtomsJCluster_;
        //! Number of atoms per bin
        int numAtomsPerBin_;
        //! 2log of na_c
        int numAtomsICluster2Log_;
        //! What type of pairlist is in use.
        PairlistType pairlistType_;
    };

    //! Constructs a grid given the type of pairlist
    Grid(PairlistType pairlistType, int ddZone, const bool& haveFep, PinningPolicy pinningPolicy);

    //! Returns the geometry of the grid bins
    const Geometry& geometry() const { return geometry_; }

    //! Returns the zone this grid belong to
    int ddZone() const { return ddZone_; }

    //! Returns the dimensions of the grid
    const GridDimensions& dimensions() const { return dimensions_; }

    //! Returns the total number of grid columns
    int numColumns() const { return dimensions_.numCells[XX] * dimensions_.numCells[YY]; }

    //! Returns the total number of grid bins
    int numBins() const { return numBinsTotal_; }

    //! Returns the bin offset of (the first bin of) this grid in the list of bins combined over all grids
    int binOffset() const { return binOffset_; }

    //! Returns the maximum number of grid bins in a column
    int numBinsColumnMax() const { return numBinsColumnMax_; }

    //! Returns the first bin index in the grid, starting at 0 in this grid
    int firstBinInColumn(int columnIndex) const { return columnToBin_[columnIndex]; }

    //! Returns the number of bins in the column
    int numBinsInColumn(int columnIndex) const
    {
        return columnToBin_[columnIndex + 1LL] - columnToBin_[columnIndex];
    }

    //! Returns the index of the first atom in the column
    int firstAtomInColumn(int columnIndex) const
    {
        return (binOffset_ + columnToBin_[columnIndex]) * geometry_.numAtomsPerBin_;
    }

    //! Returns the number of real atoms in the column
    int numAtomsInColumn(int columnIndex) const { return numAtomsPerColumn_[columnIndex]; }

    //! Returns a view of the number of non-filler, atoms for each grid column
    ArrayRef<const int> numAtomsPerColumn() const { return numAtomsPerColumn_; }

    //! Returns a view of the grid-local bin index for each grid column
    ArrayRef<const int> columnToBin() const { return columnToBin_; }

    //! Returns the number of atoms in a bin
    int numAtomsPerBin() const { return geometry_.numAtomsPerBin_; }

    //! Returns the number of atoms in the column including padding
    int paddedNumAtomsInColumn(int columnIndex) const
    {
        return numBinsInColumn(columnIndex) * geometry_.numAtomsPerBin_;
    }

    //! Returns the end of the atom index range on the grid, including padding
    int atomIndexEnd() const { return (binOffset_ + numBinsTotal_) * geometry_.numAtomsPerBin_; }

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
    ArrayRef<const BoundingBox> iBoundingBoxes() const { return iClusterBoundingBoxes_; }

    //! Returns the j-bounding boxes for all clusters on the grid
    ArrayRef<const BoundingBox> jBoundingBoxes() const { return jClusterBoundingBoxes_; }

    //! Returns the packed bounding boxes for all clusters on the grid, empty with a CPU list
    ArrayRef<const float> packedBoundingBoxes() const { return packedClusterBoundingBoxes_; }

    //! Returns the bounding boxes along z for all i-bins on the grid
    ArrayRef<const BoundingBox1D> zBoundingBoxes() const { return binBoundingBoxesZ_; }

    //! Returns the flags for all clusters on the grid
    ArrayRef<const int> clusterFlags() const { return flags_; }

    //! Returns the number of clusters for all bins on the grid, empty with a CPU geometry
    ArrayRef<const int> numClustersPerBin() const { return numClusters_; }

    //! Returns the cluster index for an atom
    int atomToCluster(int atomIndex) const
    {
        return (atomIndex >> geometry_.numAtomsICluster2Log_);
    }

    //! Returns the total number of clusters on the grid
    int numClusters() const
    {
        if (geometry_.isSimple_)
        {
            return numBinsTotal_;
        }
        else
        {
            return numClustersTotal_;
        }
    }

    //! Returns the average spatial size of a grid bin
    RVec averageBinSize() const;

    //! Resizes the bouding box and FEP flag lists for at most \p maxNumBins
    void resizeBoundingBoxesAndFlags(int maxNumBins);

    /*! \brief Sets the grid dimensions
     *
     * \param[in] ddZone           The domain decomposition zone index
     * \param[in] numAtomsTotal    The total number of atoms to put onto this grid
     * \param[in] numAtomsWithoutFillers  The number of atoms that are not filler particles
     * \param[in] lowerCorner      The minimum Cartesian coordinates of the grid
     * \param[in] upperCorner      The maximum Cartesian coordinates of the grid
     * \param[in,out] atomDensity  The atom density, will be computed when <= 0
     * \param[in] maxAtomGroupRadius  The maximum radius of atom groups
     */
    void setDimensions(int         ddZone,
                       int         numAtomsTotal,
                       int         numAtomsWithoutFillers,
                       const RVec& lowerCorner,
                       const RVec& upperCorner,
                       real*       atomDensity,
                       real        maxAtomGroupRadius);

    //! Sets the bin indices using indices in \p gridSetData and \p gridWork
    void setBinIndices(int                     ddZone,
                       int                     binOffset,
                       GridSetData*            gridSetData,
                       ArrayRef<GridWork>      gridWork,
                       Range<int>              atomRange,
                       ArrayRef<const int32_t> atomInfo,
                       ArrayRef<const RVec>    x,
                       nbnxn_atomdata_t*       nbat);

    /*! \brief Sets a non-local grid using data communicated from a different domain
     *
     * \note The cluster indices passed in \p clusterRanges use a size of the maximum
     *       of the i- and j-cluster size in atoms, whereas \c Grid use the i-cluster size.
     *
     * \param[in] ddZone      The domain decomposition zone this grid belongs to
     * \param[in] dimensions  The dimensions of the grid
     * \param[in] clusterRanges  A list of column indices and number of clusters for a column,
     *                           the list should be ordered on column index, the same column
     *                           index can appear multiple times; note that the clusters
     *                           here have size of the maximum of the i- and j-sizes
     * \param[in] binOffset   The offset of this grid in the list of bins over all grids
     * \param[in] atomInfo    A list of information for all local and non-local atoms
     * \param[in] x           The coordinates for all local and non-local atoms
     * \param[in,out] gridSetData  The data shared over all grids
     * \param[in,out] nbat    The NBNxM atom data, used here for storing the atom coordinates
     */
    void setNonLocalGrid(int                                 ddZone,
                         const GridDimensions&               dimensions,
                         ArrayRef<const std::pair<int, int>> clusterRanges,
                         int                                 binOffset,
                         ArrayRef<const int32_t>             atomInfo,
                         ArrayRef<const RVec>                x,
                         GridSetData*                        gridSetData,
                         nbnxn_atomdata_t*                   nbat);

    //! Determine in which grid columns atoms should go, store bins and atom counts in \p bins and \p numAtomsPerColumn
    static void calcColumnIndices(const GridDimensions&  gridDims,
                                  const UpdateGroupsCog* updateGroupsCog,
                                  Range<int>             atomRange,
                                  ArrayRef<const RVec>   x,
                                  int                    dd_zone,
                                  const int*             move,
                                  int                    thread,
                                  int                    nthread,
                                  ArrayRef<int>          bins,
                                  ArrayRef<int>          numAtomsPerColumn);

private:
    /*! \brief Fill a pair search bin with atoms
     *
     * Optionally sorts atoms and sets the interaction flags.
     */
    void fillBin(GridSetData*            gridSetData,
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

    // Data members

    //! The geometry of the grid clusters and bins
    Geometry geometry_;

    // The DD zone the grid belong to
    int ddZone_;

    //! The physical dimensions of the grid
    GridDimensions dimensions_;

    //! The total number of bins in this grid
    int numBinsTotal_;
    //! Index in GridSetData corresponding to bin 0 of this grid
    int binOffset_;
    //! The maximum number of bins in a column
    int numBinsColumnMax_;

    /* Grid data */
    //! The number of, non-filler, atoms for each grid column.
    HostVector<int> numAtomsPerColumn_;
    //! The grid-local bin index for each grid column
    HostVector<int> columnToBin_;

    //! The number of clusters for each column
    std::vector<int> numClusters_;

    /* Bounding boxes */
    //! Bounding boxes in z for the i-bins
    std::vector<BoundingBox1D> binBoundingBoxesZ_;
    //! 3D bounding boxes for the clusters
    std::vector<BoundingBox, AlignedAllocator<BoundingBox>> iClusterBoundingBoxes_;
    //! 3D j-bounding boxes for the case where the i- and j-cluster sizes are different
    std::vector<BoundingBox, AlignedAllocator<BoundingBox>> jClusterBoundingBoxesStorage_;
    //! 3D j-bounding boxes
    ArrayRef<BoundingBox> jClusterBoundingBoxes_;
    //! 3D bounding boxes in packed xxxx format per cluster
    std::vector<float, AlignedAllocator<float>> packedClusterBoundingBoxes_;

    //! Tells whether we have perturbed interactions, authorative source is in GridSet (never modified)
    const bool& haveFep_;

    /* Bit-flag information */
    //! Flags for properties of clusters in each bin
    std::vector<int> flags_;
    //! Signal bits for atoms in each cluster that tell whether an atom is perturbed
    std::vector<unsigned int> fep_;

    /* Statistics */
    //! Total number of clusters, used for printing
    int numClustersTotal_;
};

/*! \brief Sets the 2D search grid dimensions puts the atoms on the 2D grid
 *
 * \param[in,out] grid      The pair search grid for one DD zone
 * \param[in,out] gridWork  Working data for each thread
 * \param[in,out] bins      The grid bin list
 * \param[in] lowerCorner   The minimum Cartesian coordinates of the grid
 * \param[in] upperCorner   The maximum Cartesian coordinates of the grid
 * \param[in] updateGroupsCog  The center of geometry of update groups, can be nullptr
 * \param[in] atomRange     The range of atoms to put on this grid, may include moved atoms
 * \param[in] numGridAtomsWithoutFillers  The number of atoms that are not filler particles
 *                                        and have not moved by to another domain by DD
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
                           HostVector<int>*       bins,
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
