Neighborhood search for analysis tools {#page_analysisnbsearch}
======================================

The header nbsearch.h declares a C++ interface to a relatively flexible and
efficient neighborhood search.  It is currently implemented within the
selection module where it originated, but it does not have any dependencies on
the other selection code and can be easily split out in the future.

The emphasis is on flexibility and ease of use; one main driver is to have
one common implementation of grid-based searching to avoid replicating this in
multiple tools (and to make more tools take advantage of the significant
performance improvement this allows).  The main features that it provides:

 - Grid-based searching with any triclinic box shape that \Gromacs supports
   (i.e., a triangular box matrix and not too skewed).
 - Grid-based searching with all PBC options except for screw boundary
   conditions.
 - With no PBC, grid-based searching where the grid is constructed based on the
   bounding box of the gridded atoms.
 - Efficient, rectangular grid cells whose size is determined by particle
   density and not limited by the cutoff.
 - Transparent fallback to a simple all-pairs search if the cutoff is too long
   for the algorithm or grid searching is not otherwise supported.
 - Support for computing all distances in the XY plane only (and still
   grid-based).
 - Convenience functions for finding the shortest distance or the nearest pair
   between two sets of positions.
 - Basic support for exclusions.
 - Thread-safe handling of multiple concurrent searches with the same cutoff
   with the same or different reference positions.

Usage
=====

The neighborhood search works conceptually with two different sets of
coordinates:

 - _reference positions_: When initiating the search, you provide one set of
   reference positions that get placed on the search grid and determine the
   size of the grid.
 - _test positions_: For each set of reference positions, you provide a set of
   test positions (or a single position).  The search is performed from each
   test position, finding the reference positions within the cutoff from this
   point.  It is possible to perform multiple searches against the same set of
   reference positions (and the same grid).

To start using the neighborhood search, you need to first create an instance of
gmx::AnalysisNeighborhood.  This class allows you to set some global properties
for the search (most notably, the cutoff distance).  Then you provide the
reference positions as a gmx::AnalysisNeighborhoodPositions and PBC information
to get a gmx::AnalysisNeighborhoodSearch instance.  You can then either use
methods directly in this class to find, e.g., the nearest reference point from
a test position, or you can do a full pair search that returns you all the
reference-test pairs within a cutoff.  The pair search is performed using an
instance of gmx::AnalysisNeighborhoodPairSearch that the search object returns.
Methods that return information about pairs return an instance of
gmx::AnalysisNeighborhoodPair, which can be used to access the indices of
the reference and test positions in the pair, as well as the computed distance.
See the class documentation for these classes for details.

For use together with selections, an instance of gmx::Selection or
gmx::SelectionPosition can be transparently passed as the positions for the
neighborhood search.

Implementation
==============

This section provides a high-level overview of the algorithm used.  It is not
necessary to understand all the details to use the API, but it can be useful to
get the best performance out of it.  The main audience is developers who may
need to extend the API to make it suitable for more cases.

The grid for the search is initialized based on the reference positions and the
PBC information:

 - The grid cells are always rectangular, even for fully triclinic boxes.
 - If there is no PBC, the grid edges are defined from the bounding box of the
   reference positions; with PBC, the grid covers the unit cell.
 - The grid cell size is determined such that on average, each cell contains
   ten particles.  Special considerations are in place for cases where the grid
   will only be one- or two-dimensional because of a flat box.
 - If the resulting grid has too few cells in some dimensions, the code
   falls back automatically to an all-pairs search.  For correct operation, the
   grid algorithm needs three cells in each dimension, but the code can fall
   back to a non-gridded search for each dimension separately.
 - If the resulting grid has so few cells that the search would anyways
   consider all (or nearly all) cell pairs, the search falls back to a
   simple search.
 - The initialization also pre-calculates the shifts required across the
   periodic boundaries for triclinic cells, i.e., the fractional number of
   cells that the grid origin is shifted when crossing the periodic boundary in
   Y or Z directions.
 - Finally, all the reference positions are mapped to the grid cells.

There are a few heuristic numbers in the above logic: the average number of
particles within a cell, and the cutover point from grid to an all-pairs
search.  These have not been particularly optimized for best performance.

When doing the search for test positions, each test position is considered
independently:

 - The coordinates of the test position are mapped to the grid coordinate
   system.  The coordinates here are fractional and may lay outside the grid
   for non-periodic dimensions.
 - The bounding box of the cutoff sphere centered at the mapped coordinates is
   determined, and each grid cell that intersects with this box is used for
   searching the reference positions.  So the searched grid cells may vary
   depending on the coordinates of the test position, even if the test position
   is within the same cell.
 - Possible triclinic shifts in the grid are considered when looping over the
   cells in the cutoff box if the coordinates wrap around a periodic dimension.
   This is done by shifting the search range in the other dimensions when the Z
   or Y dimension loop crosses the boundary.
