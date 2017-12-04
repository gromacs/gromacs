/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief
 * Implements gmx::analysismodules::tide.
 *
 * \author Camilo Aponte <ca.aponte.uniandes.edu.co>
 * \author Rodolfo Briones <fitobriones@gmail.com>
 * \author Carsten Kutzner <ckutzne@gwdg.de>

 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "tide.h"

#include <map>
#include <string.h>

#include "gromacs/trajectoryanalysis.h"
#include "gromacs/applied-forces/maputil.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


namespace gmx
{

namespace analysismodules
{

namespace
{

//! \brief Convert between XPLOR and GROMACS length units
const float NM2AA = 10.0;

//! \brief Accumulate block average density
void addToDensity( t_mapdata *targetMap, const t_mapdata &addedMap  )
{
    std::transform(std::begin(targetMap->vox), std::end(targetMap->vox), std::begin(addedMap.vox), std::begin(targetMap->vox), std::plus<float>());
}

//! \brief Set density of a map to zero
void setDensityToZero(t_mapdata *map)
{
    std::fill(std::begin(map->vox), std::end(map->vox), 0.f);
}

//! \brief Compute timed-average average density
void divideDensityByValue ( gmx_unused t_mapdata *map_av, float Nframes)
{
    std::transform(std::begin(map_av->vox), std::end(map_av->vox), std::begin(map_av->vox), [Nframes](float voxelValue){return voxelValue/Nframes; });
}

//! \brief Gets the index in an uni-dimensional array from the three coordinates of a three-dimensional array
int index_array(int i, int j, int k, int Nx, int Ny)
{
// The unidimensional array contains Nx*Ny*Nz elements:  0 ... Nx*Ny*Nz-1
    return k*(Nx*Ny) + j*Nx + i;
}

/*! \brief Compare variables relevant for tide (cell,grid, origin and map_dim) of map1 and map2.
 *
 *  \return true if their dimensions are equal and false otherwise===
 */
bool mapsHaveSameGrid(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2 )
{
    const float cell_tolerance = 0.0001;   // diff smaller than that will be treated as zero.

    for (int i = 0; i < 2*DIM; i++)
    {
        if ( (map1->cell[i] - map2->cell[i])*(map1->cell[i] - map2->cell[i]) > cell_tolerance*cell_tolerance)
        {
            return false;
        }
    }
    return (map1->grid == map2->grid) && (map1->origin == map2->origin) && (map1->map_dim == map2->map_dim);
}

//! \brief set mask
void set_mask(gmx_unused t_mapdata *mask, float mask_cutoff)
{
    std::transform(std::begin(mask->vox), std::end(mask->vox), std::begin(mask->vox), [mask_cutoff](float maskMapValue){return maskMapValue <= mask_cutoff ? 0 : 1; });
}

//! \brief Subtract maps: rho(map3) = rho(map2)-rho(map1)
void subtract_maps(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *map3, gmx_unused t_mapdata *mask)
{
    // do subtraction if maps have the same dimensions, else abort with fatal error
    if (!mapsHaveSameGrid(map1, map2) || !mapsHaveSameGrid(map1, map3) || !mapsHaveSameGrid(map1, mask))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    for (size_t l = 0; l < map1->vox.size(); l++)
    {
        if (mask->vox[l] > 0.5)
        {
            map3->vox[l] = map2->vox[l]-map1->vox[l];
        }
        else
        {
            map3->vox[l] = nanf("");         // masked point is not considered
        }
    }
}

//! \brief log of map: rho(map2) = - ln ( rho(map) )  PSEUDO FREE ENERGY
void minus_log_map(gmx_unused t_mapdata *map, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *mask)
{
    // do operation if maps have the same dimensions, else abort else abort with fatal error
    if (!mapsHaveSameGrid(map, map2) || !mapsHaveSameGrid(map, mask))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    for (size_t l = 0; l < map->vox.size(); l++)
    {
        if (mask->vox[l] > 0.5)                // discard masked points
        {
            if (map->vox[l] > 0)
            {
                map2->vox[l] = -log( map->vox[l] );
            }
            else
            {
                gmx_fatal(FARGS, "Density at %d is zero or less in the map. Discard this point by using -mask.\n", l);
            }
        }
        else
        {
            map2->vox[l] = nanf("");                 // masked point is not considered
        }
    }
}

//! \brief log of division of maps: rho(map3) = - ln ( rho(map2) / rho(map1) )  PSEUDO FREE ENERGY
void minus_log_division_maps(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *map3, gmx_unused t_mapdata *mask)
{
    // do operation if maps have the same dimensions, else abort else abort with fatal error
    if (!mapsHaveSameGrid(map1, map2) || !mapsHaveSameGrid(map1, map3) || !mapsHaveSameGrid(map1, mask))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    for (size_t l = 0; l < mask->vox.size(); l++)
    {
        if (mask->vox[l] > 0.5)            // discard masked points
        {
            if (map1->vox[l] != 0)
            {
                map3->vox[l] = -log( map2->vox[l] / map1->vox[l] );
            }
            else
            {
                gmx_fatal(FARGS, "Density at point %d is zero in the dividing map1: -ln[rho(map2)/rho(map1)]. You could discard this point from the calculation by using -mask.\n", l);
            }
        }
    }
}

//! \brief Shift map: rho(map) = rho(map) + S  , with S a shift scalar value
void shift_map(gmx_unused t_mapdata *map, float S )
{
    std::transform(std::begin(map->vox), std::end(map->vox), std::begin(map->vox), [S](float voxelValue){return voxelValue+S; });
}


//! \brief Compute spatial average over the map: <rho(map)> over all cell positions
void do_spatial_average_maps(gmx_unused t_mapdata *map, bool incl_zero, gmx_unused t_mapdata *mask  )
{
    if (!mapsHaveSameGrid(map, mask)) // map and mask must have same dimensions
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    map->mean = 0.0;
    map->rms  = 0.0;
    float nsamples      = 0;
    for (size_t l = 0; l < mask->vox.size(); l++)
    {
        if (mask->vox[l] > 0.5)                    // 0=masked point , 1=unmasked point to be considered
        {
            if (map->vox[l] != 0 || incl_zero)     // zero density point
            {
                map->mean += map->vox[l];
                map->rms  += map->vox[l]*map->vox[l];
                nsamples++;
            }
        }
    }
    if (nsamples > 0)
    {
        map->mean /= nsamples;
        map->rms   = sqrt(map->rms/(nsamples-1) - map->mean*map->mean * nsamples / (nsamples-1)   );
    }
}

//! \brief TODO
float calc_sum_rhoA_rhoB_mask (gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *mask)
{
    if (!mapsHaveSameGrid(map1, map2) || !mapsHaveSameGrid(map1, mask))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    double sum_ijk = 0.0;
    for (int l = 0; l < map1->map_dim[XX]*map1->map_dim[YY]*map1->map_dim[ZZ]; l++)
    {
        if (mask->vox[l] > 0.5)            // 0=masked point , 1=unmasked point to be considered
        {
            sum_ijk += map1->vox[l]*map2->vox[l];
        }
    }
    return sum_ijk;
}

//! \brief TODO
float global_correlation_coeff(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *mask)
{
    if (!mapsHaveSameGrid(map1, map2) || !mapsHaveSameGrid(map1, mask))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    float     numerator;
    float     denominator;
    float     cc;
    numerator   = calc_sum_rhoA_rhoB_mask(map1, map2, mask);
    denominator = sqrt(calc_sum_rhoA_rhoB_mask(map2, map2, mask) * calc_sum_rhoA_rhoB_mask(map1, map1, mask));
    cc          = numerator / denominator;
    return cc;
}

//! \brief Calculate the local correlation
void local_correlation_coeff(t_mapdata *map1, t_mapdata *map2, t_mapdata *mask, t_mapdata *loc_corr, int Nodd)
{
    int Nhalf = int(Nodd/2);    // center point of the sliding grid (e.g. if Nodd=5 ; then Nhalf=2, Note that it starts counting from 0)

    // This can only work for maps with equal dimensions
    if (!mapsHaveSameGrid(map1, map2) || !mapsHaveSameGrid(map1, mask) || !mapsHaveSameGrid(map1, loc_corr))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }

    // Initialize the output loc_corr map, setting the borders and masked points, with non-assigned NAN values
    std::fill(std::begin(loc_corr->vox), std::end(loc_corr->vox), nanf(""));

    // loop over the map1 (reduced by an amount Nhalf at both ends)
    for (int i = Nhalf; i < map1->map_dim[XX]-Nhalf; i++)
    {
        for (int j = Nhalf; j < map1->map_dim[YY]-Nhalf; j++)
        {
            for (int k = Nhalf; k < map1->map_dim[ZZ]-Nhalf; k++)
            {
                // reset correlation variables
                float sum12       = 0.0;
                float sum11       = 0.0;
                float sum22       = 0.0;
                //slide over Nodd points, centered at Nhalf to accumulate the density
                for (int i_sli = i - Nhalf; i_sli <=  i + Nhalf; i_sli++)
                {
                    for (int j_sli = j - Nhalf; j_sli <=  j + Nhalf; j_sli++)
                    {
                        for (int k_sli = k - Nhalf; k_sli <=  k + Nhalf; k_sli++)
                        {
                            int id_aux = index_array(i_sli, j_sli, k_sli, map1->map_dim[XX], map1->map_dim[YY]);
                            if (mask->vox[id_aux] > 0.5)        // skipped masked points
                            {
                                sum12 += map1->vox[id_aux]*map2->vox[id_aux];
                                sum11 += map1->vox[id_aux]*map1->vox[id_aux];
                                sum22 += map2->vox[id_aux]*map2->vox[id_aux];
                            }
                        }
                    }
                }

                // compute correlation
                float numerator   = sum12;
                float denominator = sqrt(sum22 * sum11);
                float cc          = numerator / denominator;

                // assign correlation to middle point of sliding grid (i,j,k)
                int id_aux            = index_array(i, j, k, map1->map_dim[XX], map1->map_dim[YY]);
                loc_corr->vox[id_aux] = cc;
            }
        }
    }
}

//! \brief Set the Gaussian coefficients A and B for the six atom types H, C, N, O, P, and S
void set_gaussian_coeff( std::string gauss_coeff_, std::map < std::string, std::vector < float>> *gaussianAmplitude, std::map < std::string, std::vector < float>> *gaussianPrecision, int Ngaussians)
{
    gaussianAmplitude->clear();
    gaussianPrecision->clear();

    // if there is an input Gaussian_coeff file, consider the values in there:
    if (!gauss_coeff_.empty())
    {

        fprintf(stdout, "Reading file with Gaussian coefficients\n");

        FILE *pfile;
        pfile = fopen(gauss_coeff_.c_str(), "r");

        // loop over the atoms
        char currentAtomType[1];
        while (fscanf(pfile, "%s", currentAtomType) == 1)
        {
            // gaussianAmplitude parameters
            for (int col = 0; col < Ngaussians; col++)
            {
                float currentgaussianAmplitudeFactor;
                fscanf(pfile, "%e", &currentgaussianAmplitudeFactor);
                (*gaussianAmplitude)[std::string(currentAtomType)].push_back(currentgaussianAmplitudeFactor);
            }

            // gaussianPrecision parameters
            for (int col = 0; col < Ngaussians; col++)
            {
                float currentgaussianPrecisionFactor;
                fscanf(pfile, "%e", &currentgaussianPrecisionFactor);
                (*gaussianPrecision)[std::string(currentAtomType)].push_back(currentgaussianPrecisionFactor);
            }
        }
        fclose (pfile);
    }
    else                                // take the defaults (2 Gaussians from cryoEM)
    {
        fprintf(stdout, "Taking default Gaussian coefficients from cryo EM for sums with 2 Gaussians \n");

        (*gaussianAmplitude)["H"] = {1.242258e-01, 1.133291e+00};
        (*gaussianAmplitude)["C"] = {5.948015e-01, 5.657122e+00};
        (*gaussianAmplitude)["N"] = {7.587542e-01, 7.705668e+00};
        (*gaussianAmplitude)["O"] = {9.061284e-01, 9.848492e+00};
        (*gaussianAmplitude)["P"] = {1.056651e+00, 1.083431e+01};
        (*gaussianAmplitude)["S"] = {1.335159e+00, 1.271311e+01};

        (*gaussianPrecision)["H"] = {1.718695e+00, 9.919200e+00};
        (*gaussianPrecision)["C"] = {1.700190e+00, 1.041647e+01};
        (*gaussianPrecision)["N"] = {2.236738e+00, 1.342803e+01};
        (*gaussianPrecision)["O"] = {2.786056e+00, 1.644934e+01};
        (*gaussianPrecision)["P"] = {1.390085e+00, 1.081600e+01};
        (*gaussianPrecision)["S"] = {1.711987e+00, 1.226038e+01};

    }
}

//! \brief  Output integer position in the array of atomtypes based on first character of atomname
std::vector<std::string> attypeStrings_from_top(const t_topology &top)
{
    std::vector<std::string> attypeString;
    for (int n = 0; n < top.atoms.nr; n++)
    {
        attypeString.push_back( std::string(*(top.atoms.atomname[n])));
    }
    return attypeString;
}


/*! \brief Outputs the density map to text file using XPLOR format
 *
 * For the format see: http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
void print_to3dmap(const char *inputFileName, gmx_unused t_mapdata *map, gmx_unused t_mapdata *mask )
{
    if (inputFileName == nullptr)
    {
        fprintf(stdout, "No filename given for X-PLOR output file, not outputting densities\n");
        return;
    }

    if (!mapsHaveSameGrid(map, mask))
    {
        gmx_fatal(FARGS, "Map and mask dimensions are different.");
    }

    FILE *fp = gmx_ffopen(inputFileName, "w");

    // === line 1 : Empty line ===
    fprintf(fp, "\n");

    // === Lines 2-4 : Title and remarks ===
    fprintf(fp, "       2\n");
    fprintf(fp, "REMARKS Lattice Nbins: %d %d %d\n", map->map_dim[XX], map->map_dim[YY], map->map_dim[ZZ]);
    float res[DIM]; for (int i = 0; i < DIM; i++)
    {
        res[i] = map->cell[i] / map->grid[i];
    }
    fprintf(fp, "REMARKS Lattice resolution (nm) %f %f %f\n", res[XX], res[YY], res[ZZ]);

    // === Line 5 : Grid points ===
    //  nbox = total number of grid points along the a,b, and c cell edges (of the simulation box)

    //  starting and stopping grid points along each cell edge in the portion of the map that is written
    int NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX;

    NA = map->grid[XX];
    NB = map->grid[YY];
    NC = map->grid[ZZ];;

    AMIN = map->origin[XX];
    BMIN = map->origin[YY];
    CMIN = map->origin[ZZ];

    AMAX = map->origin[XX]+map->map_dim[XX]-1;
    BMAX = map->origin[YY]+map->map_dim[YY]-1;
    CMAX = map->origin[ZZ]+map->map_dim[ZZ]-1;

    fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX  );

    // === Line 6 : box size (Angstrom) and angles ===
    fprintf(fp, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n", map->cell[XX]*NM2AA, map->cell[YY]*NM2AA, map->cell[ZZ]*NM2AA, map->cell[3], map->cell[4], map->cell[5]   );

    // === Line 7 : string ===
    fprintf(fp, "ZYX\n");

    for (int k = 0; k < map->map_dim[ZZ]; k++)     // loop over z
    {                                              // Print the k value
        fprintf(fp, "       %d\n", k);

        // variable to monitor the number of columns printed
        int ncols = 0;

        for (int j = 0; j < map->map_dim[YY]; j++)           // loop over y
        {
            for (int i = 0; i < map->map_dim[XX]; i++)       // loop over x
            {
                int index_aux = index_array(i, j, k, map->map_dim[XX], map->map_dim[YY]);

                if (mask->vox[index_aux] > 0.5)
                {
                    // Print the density values in 6 columns
                    fprintf(fp, "%12.5E", map->vox[index_aux] );
                }
                else                // print nan to indicate that this point is masked
                {
                    fprintf(fp, "%12.5E", nanf("") );
                }
                ncols++;
                // enter when the number of printed columns =6
                if (ncols == 6 || ( i == (map->map_dim[XX]-1) && j == (map->map_dim[YY]-1) ) )
                {
                    fprintf(fp, "\n");
                    ncols = 0;
                }
            }
        }
    }

    //  === Footer lines
    fprintf(fp, "   -9999\n");
    fprintf(fp, "%13.5E%13.5E\n", map->mean, map->rms);

    gmx_ffclose(fp);
}

void printMapInfo(const t_mapdata &map)
{
    printf( "\n");
    printf( "Header:\n");
    printf( "Title      : %s\n", map.title);
    printf( "Datamode   : %u\n", map.datamode);
    printf( "Cell       : %g A  %g A  %g A  %g %g %g\n", map.cell[0], map.cell[1], map.cell[2], map.cell[3], map.cell[4], map.cell[5]);
    printf( "Grid       : %d x %d x %d\n", map.grid[XX], map.grid[YY], map.grid[ZZ]);
    printf( "Origin     : %d A  %d A  %d A\n", map.origin[XX], map.origin[YY], map.origin[ZZ]);
    printf( "Axes order : %d %d %d\n", map.axes_order[0], map.axes_order[1], map.axes_order[2]);
    printf( "Map_dim    : %d %d %d\n", map.map_dim[0], map.map_dim[1], map.map_dim[2]);
    printf( "Spacegroup : %d\n", map.spacegroup);
    printf( "Min        : %g\n", map.min);
    printf( "Max        : %g\n", map.max);
    printf( "Mean       : %g\n", map.mean);
    printf( "RMS        : %g\n", map.rms);
    printf( "Skew_trans : %g %g %g\n", map.skew_trans[0], map.skew_trans[1], map.skew_trans[2]);
    printf( "Skew_mat   : ");
    for (int i = 0; i < 9; i++)
    {
        printf( "%g ", map.skew_mat[i]);
    }
    printf( "\n\n");
}

//! \brief Init map function
void init_map(gmx_unused t_mapdata *map, const std::string &title, int datamode, const std::array<float, 2*DIM> &cell, const std::array<int, DIM> &grid, const std::array<int, DIM> &origin, const std::array<int, DIM> &axes_order, const std::array<int, DIM> &map_dim, int spacegroup, double min, double max, float mean, float rms, const std::array<float, DIM*DIM> &skew_mat, const std::array<float, DIM> &skew_trans)
{
    map->title      = gmx_strdup(title.c_str());
    map->datamode   = datamode;   //stored as Reals
    map->cell       = cell;       // cell 0...2 are the box sizes // cell 3...5 are the angles
    map->grid       = grid;       // Number of grid points within box
    map->origin     = origin;     // Origin of the working lattice
    map->axes_order = axes_order; // x, y, z (guess)
    map->map_dim    = map_dim;    // the working lattice
    map->spacegroup = spacegroup; // this will be assigned arbitrarily. It won't be used in gtide anyways
    // minimum and maximum density values
    map->min = min;
    map->max = max;
    // average and stdev
    map->mean = mean;
    map->rms  = rms;
    // skew options will be ignored (0)
    map->skew_mat   = skew_mat;
    map->skew_trans = skew_trans;
}

/*!  \brief Read an input density map using XPLOR format
 *
 * For the format see: http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
void read_xplor_3dmap(const char *inputFileName, gmx_unused t_mapdata *ref_map )
{
    if (inputFileName == nullptr)
    {
        fprintf(stdout, "No filename given for X-PLOR input file\n");
        return;
    }

    FILE *fp = gmx_ffopen(inputFileName, "r");

    //  starting and stopping grid points along each cell edge in the portion of the map that is written
    rvec                 box, angle;
    std::array<int, DIM> N = {};
    std::array<int, DIM> MIN = {};
    std::array<int, DIM> MAX = {};
    std::array<int, DIM> bin = {};

    // ==== header ====

    // determine which line has the ZYX string
    char  tmp_char[80];    // auxiliary variable to copy tmp header information of the map
    int   row_ZYX = 0;
    do
    {
        fgets(tmp_char, 80, fp);
        row_ZYX++;
    }
    while (strncmp(tmp_char, "ZYX", 3 ) != 0);
    rewind ( fp );

    std::vector<float> lattice; // one-D array to save the lattice
    // loop over lines before string ZYX
    for (int row = 1; row <= row_ZYX; row++)
    {
        if (row < row_ZYX-2)
        {
            fgets(tmp_char, 80, fp);                      // === lines 1-4: empty space or header remarks
        }
        // === Line 5 : Grid points ===
        if (row == row_ZYX-2)
        {
            fscanf(fp, "%8d", &N[XX]);
            fscanf(fp, "%8d", &MIN[XX]);
            fscanf(fp, "%8d", &MAX[XX]);

            fscanf(fp, "%8d", &N[YY]);
            fscanf(fp, "%8d", &MIN[YY]);
            fscanf(fp, "%8d", &MAX[YY]);

            fscanf(fp, "%8d", &N[ZZ]);
            fscanf(fp, "%8d", &MIN[ZZ]);
            fscanf(fp, "%8d", &MAX[ZZ]);
        }
        // === Line 6 : box size (in Angstrom) and angles (degree)===
        if (row == row_ZYX-1)
        {
            fscanf(fp, "%12f", &box[XX]);
            fscanf(fp, "%12f", &box[YY]);
            fscanf(fp, "%12f", &box[ZZ]);

            fscanf(fp, "%12f", &angle[XX]);
            fscanf(fp, "%12f", &angle[YY]);
            fscanf(fp, "%12f", &angle[ZZ]);

            // having read N,MIN,MAX, box and angle, it is possible to compute bin
            for (int i = 0; i < DIM; i++)
            {
                // set box in nm
                box[i] /= NM2AA;

                //number of bins (working lattice)
                bin[i] = MAX[i]-MIN[i]+1;                    //Amax-Amin+1
            }
            // allocate space for the density one-dimensional array
            lattice.resize(bin[XX] * bin[YY] * bin[ZZ]);                        //lattice
        }
        // === Line 7 : ZYX string ===
        if (row == row_ZYX)
        {
            fscanf(fp, "%s", tmp_char);
        }
    }

    // ==== Line 8 onwards: the density map values =====
    for (int k = 0; k < bin[ZZ]; k++)              // loop over z
    {                                              // Print the k value
        int k_read;
        fscanf(fp, "%d", &k_read);                 //read value of k

        if (k != k_read)
        {
            printf("k!=k_read: %d %d\n", k, k_read);
            break;
        }
        else
        {
            for (int j = 0; j < bin[YY]; j++)                           // loop over y
            {
                for (int i = 0; i < bin[XX]; i++)                       // loop over x
                {
                    int index_aux = index_array(i, j, k, bin[XX], bin[YY]);
                    fscanf(fp, "%12f", &lattice[index_aux]);
                }
            }
        }
    }

    // === Second to last line : -9999 string===
    int cross_check_line;
    fscanf(fp, "%d", &cross_check_line);
    if (cross_check_line != -9999)
    {
        printf("There was an error: it reached second to last line and it is not -9999: %d\n", cross_check_line);
    }

    // === Last line : average and stdeviations ===
    float av;
    float stdev;
    fscanf(fp, "%f", &av);
    fscanf(fp, "%f", &stdev);

    fclose (fp);

//=== assign the information read in the xplor file to the map ===
    std::string              title = {"ref"};

    std::array<float, 2*DIM> cell;
    for (int i = 0; i < DIM; i++)
    {
        cell[i] = box[i];                        // cell 0...2 are the box sizes
    }
    for (int i = 3; i < 2*DIM; i++)
    {
        cell[i] = angle[i-3];                          // cell 3...5 are the angles
    }
    //  N Number of grid points within box
    //  MIN[i]  Origin of the working lattice
    std::array<int, DIM> axes_order = {{0, 1, 2}};              // x, y, z (guess) but  not used in tide anyways
    // bin = map_dim is the working lattice
    // this will be assigned arbitrarily. It won't be used in gtide anyways

    //skew options will be ignored (0)
    std::array<float, DIM*DIM> skew_mat = {};
    std::array<float, DIM>     skew_trans = {};

    auto minmax = std::minmax_element(std::begin(lattice), std::end(lattice));
    init_map(ref_map, title, 2, cell, N, MIN,  axes_order, bin, 1, *(minmax.first), *(minmax.second), av, stdev, skew_mat, skew_trans);

    /* ==== 1d continuous array for 3d voxel density  ====      */
    allocate_density_grid(ref_map->map_dim, &ref_map->vox);

    // copy lattice into vox
    ref_map->vox = lattice;
}

//! \brief Accummulate density every time step
void accumulate_denstiy( gmx_unused t_mapdata *map_av, gmx_unused t_mapdata *map_inst, gmx_unused t_mapdata *map_sd, Selection sel,  RVec res, RVec L0, IVec ncutoff, float cutoff, const std::map < std::string, std::vector < float>> &Ai, const std::map < std::string, std::vector < float>> &Bi,  const std::vector<std::string> &attypeString)
{
    // Coordinate vector with respect to the simulation box
    rvec R_trans;

    // vector coordinates of the (nx,ny,nz) grid point
    rvec grid_vec;

    // Distance between the atom and the grid_vector:  R_atom_grid= R_trans - grid_vec
    rvec  R_atom_grid;

    float dist_atom_grid;

    //=== INITIALIZE THE INSTANTANEOUS MAP LATTICE===
    setDensityToZero(map_inst);

    // === Count the visit of each one of the surrounding atoms (i.e. lipids) in the proper grid of the lattice

    // === Loop over the surrounding atoms (i.e. lipids )
    for (int atomIndex = 0; atomIndex < sel.atomCount(); atomIndex++)
    {
        // Select the first character of the atom name
        std::string atomName           = attypeString[sel.atomIndices()[atomIndex]].substr(0, 1);
        const auto &gaussianAmplitudes = Ai.at(atomName);
        const auto &gaussianPrecision  = Bi.at(atomName);

        // Change of system of coordinates to one where the origin is the origin of the lattice
        // L0 origin of the lattice coordinates
        // R_trans = coordinates of the atom with respect to the lattice
        // R_trans= R - L0
        rvec_sub( sel.position(atomIndex).x(), L0, R_trans );

        // === Discretize the Rtrans to locate it in the lattice
        const int i = floor( R_trans[XX] / res[XX]);
        const int j = floor( R_trans[YY] / res[YY]);
        const int k = floor( R_trans[ZZ] / res[ZZ]);

        // === only consider the atoms located in the lattice ===
        if (i >= 0 && i < map_av->map_dim[XX] && j >= 0 && j < map_av->map_dim[YY] && k >= 0 && k < map_av->map_dim[ZZ])
        {
            // Increment the density of the grid points to the atom, laying within a certain cutoff
            // loops over grid points around the (i,j,k) grid point
            for (int nx = i-ncutoff[XX]; nx <= i+ncutoff[XX]; nx++)
            {
                for (int ny = j-ncutoff[YY]; ny <= j+ncutoff[YY]; ny++)
                {
                    for (int nz = k-ncutoff[ZZ]; nz <= k+ncutoff[ZZ]; nz++)
                    {
                        // === only consider grid points lying in the lattice ===
                        if (nx >= 0 && nx < map_av->map_dim[XX] && ny >= 0 && ny < map_av->map_dim[YY] && nz >= 0 && nz < map_av->map_dim[ZZ])
                        {

                            // vector coordinates of the (nx,ny,nz) grid point
                            grid_vec[XX] = nx*res[XX];
                            grid_vec[YY] = ny*res[YY];
                            grid_vec[ZZ] = nz*res[ZZ];

                            // Distance between the atom and the grid_vector:
                            rvec_sub( R_trans, grid_vec, R_atom_grid ); // R_atom_grid  = R_trans - grid_vec

                            // Calculate the norm of R_atom_grid	(Angstrom) to coincide with the ai, bi and c factors
                            dist_atom_grid = NM2AA * norm(  R_atom_grid );

                            // === Only contribute to the density if the distance is smaller than the cutoff
                            if (dist_atom_grid < cutoff)
                            {

                                // === Increment the density at the position (nx,ny,nz)
                                real currentContribution = 0;
                                // === Gaussian contribution (two or four gaussians depending on the atom) to the density at this grid point
                                // contributions to the atomic potentials according to  equation (6) in J. Electron Microscopy, 56(4):131-140 (2007)
                                for (size_t g = 0; g < gaussianAmplitudes.size(); g++)
                                {
                                    currentContribution += gaussianAmplitudes[g]*exp(-gaussianPrecision[g]*dist_atom_grid * dist_atom_grid );
                                }
                                int index_aux = index_array(nx, ny, nz, map_av->map_dim[XX], map_av->map_dim[YY]);
                                map_inst->vox[index_aux] += currentContribution;

                            } // Finish if statement that only distances smaller than the cut off are considered for the density
                        }     // finish if statement that grid points are within the lattice
                    }         // finish over z grid points
                }             // finish over y grid points
            }                 // finish over x grid points
        }                     // finish if statement considering atoms lying in the lattice
    }                         // finish loop over surrounding atoms

    // ===  update average and stdev maps ===
    addToDensity(map_av, *map_inst);
    std::transform(std::begin(map_sd->vox), std::end(map_sd->vox), std::begin(map_inst->vox), std::begin(map_sd->vox), [](float oldSD, float voxelValue){return oldSD+voxelValue*voxelValue; });
}

//! \brief Compute timed-average standard deviation density
void sd_map( const  t_mapdata &map_av,  gmx_unused t_mapdata *map_sd, float Nframes  )
{
    if (Nframes > 1)
    {
        std::transform(std::begin(map_av.vox), std::end(map_av.vox), std::begin(map_sd->vox), std::begin(map_sd->vox), [Nframes](float av, float sd){
                           return sqrt(sd/(Nframes-1.0) - av*av*Nframes/(Nframes-1.0) );
                       }
                       );
    }
}

//! \brief Print the array in dat format: x, y, z, rho, sigma
void print_rho_sigma( const char *inputFileName,   rvec res, gmx_unused t_mapdata *map_av,  gmx_unused t_mapdata *map_sd  )
{
    // Nbins of the lattice: (bin[XX],bin[YY],bin[ZZ])
    // Resolution (size of each grid) (resx,resy,resz) (nm)
    FILE *pfile;
    pfile = fopen(inputFileName, "w");

    // === line 1 : Header lines ===
    fprintf(pfile, "# Rho ( and its standard dev) vs x y z.\n");
    fprintf(pfile, "# grid_res= %f %f %f (nm), dimensions: bin[XX]: %d bin[YY]: %d bin[ZZ]: %d\n", res[XX], res[YY], res[ZZ], map_av->map_dim[XX], map_av->map_dim[YY], map_av->map_dim[ZZ]);
    fprintf(pfile, "# x y z (Angstrom) rho sigma (Angstrom^-3)\n");

    // === Density array ===
    for (int i = 0; i < map_av->map_dim[XX]; i++)         // loop over x
    {
        for (int j = 0; j < map_av->map_dim[YY]; j++)     // loop over y
        {
            for (int k = 0; k < map_av->map_dim[ZZ]; k++) // loop over z

            {                                             //=== set the absolute coordinate system (the simulation box)by adding (minx,miny,minz) to each lattice value
                // x,y,z in Angstrom
                const float x = NM2AA*(i+map_av->origin[XX])*res[XX];
                const float y = NM2AA*(j+map_av->origin[YY])*res[YY];
                const float z = NM2AA*(k+map_av->origin[ZZ])*res[ZZ];

                int         index_aux = index_array(i, j, k, map_av->map_dim[XX], map_av->map_dim[YY]);
                fprintf(pfile, "%f %f %f %d %d %d %E %E\n", x, y, z, i, j, k, map_av->vox[index_aux], map_sd->vox[index_aux] );
            }
        }
    }
    fclose (pfile);
}


/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class TimeAveragedDensity final : public TrajectoryAnalysisModule
{
    public:
        TimeAveragedDensity();

        void initOptions(IOptionsContainer          *options,
                         TrajectoryAnalysisSettings *settings) override;
        void initAnalysis(const TrajectoryAnalysisSettings &settings,
                          const TopologyInformation        &top) override;

        void analyzeFrame(int                           frnr,
                          const t_trxframe             &fr,
                          t_pbc                        *pbc,
                          TrajectoryAnalysisModuleData *pdata ) override;

        void finishAnalysis(int nframes) override;
        void writeOutput() override;


    private:
        class ModuleData;

        // entire cell
        std::array<int, DIM> nGridPoints = {{250, 250, 250}}; // (NA,NB,NC)
        std::array<int, DIM> Nmin        = {{0, 0, 0}};       // (AMIN,BMIN,CMIN)
        std::array<int, DIM> Nmax        = {{249, 249, 249}}; // (AMAX,BMAX,CMAX)
        RVec                 box         = {10., 10., 10.};   //(boxA,boxB,boxC)

        // subcell where the density is actually computed
        std::array<int, DIM>             bin = {{0, 0, 0}};          // Grid points (0 0 0) (working lattice)

        RVec                             res         = {0., 0., 0.}; // grid resolution  (box/N)
        RVec                             L0          = {0., 0., 0.}; // origin of the lattice with respect to the origin of the map Nmin*res
        float                            cutoff      = 0.5;          // (nm) density will be incremented in the grid points within a cutoff distance of a certain atom
        float                            mask_cutoff = 1e-5;         //|rho|<mask_rheshold will be ignored from the calculation
        float                            time_window = 1.0;          // for time-resolved global correlation
        int                              Nslide      = 5;            // Size of the sliding grid for local-correlation calculation
        // ==== selection variables ===
        Selection                        refsel_;
        SelectionList                    sel_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;

        // === Definition of the parameters for the input/output files ===

        // strings to read the input and output file names
        std::string gauss_coef_; // input structure factor file
        std::string map_av_;     // average map (3d)
        std::string map_sd_;     // standard deviation (3d)
        std::string mi_;
        std::string diff_;
        std::string mask_;
        std::string mlog_;
        std::string loc_corr_;
        std::string rho_; // density (3d), same as the map but in a differnt format: x, y, z, rho, sigma
        std::string rho_diff_, rho_mlog_;
        std::string corr_;

        FILE       *pfile_corr = nullptr; // to print time-resolved correlation

        // === structure factor variables ===
        char               gauss_coef_atname[10] = {'\0'}; // atom names for crosscheck
        char               gauss_coef_attype[10] = {'\0'}; // atom types
        int                Ngaussians            = 2;      // Number of Gaussians that add up to the density
        bool               Norm                  = true;   // Normalization
        bool               incl_zero             = false;  // include zero-density values in the spatial average and stdev calculation (over the grid points)
        bool               shift_mlog            = true;   // remove mean value from mlog map

        std::map < std::string, std::vector < float>> A;   // Ai factors (Ngaussians x Natoms 1D array)
        std::map < std::string, std::vector < float>> B;   // Bi factors (Ngaussians x Natoms 1D array)
        std::vector<std::string> attypeString;             // array will store the index of the atom type for the A and B Gaussian-coefficient arrays

        // === define the lattice that will store the density for each point===
        t_mapdata          map_ref;             /* The reference=experimental map to compare to        */
        t_mapdata          map_av;              /* The computed timed-averaged map         */
        t_mapdata          map_inst;            /* The computed instantaneous         */
        t_mapdata          map_sd;              /* The computed timed-averaged standard deviation map */
        t_mapdata          map_diff;            /* diff map = computed_map - ref_map */
        t_mapdata          map_mask;            /* mask map*/
        t_mapdata          map_mlog;            /* minus log map*/
        t_mapdata          map_loc_corr;        /* local correlation map*/
        t_mapdata          map_block;           /* map for block averaging*/

        IVec               ncutoff = {0, 0, 0}; // === Definition of the cutoff and other gaussian-density related variables
        // Density will be increment in the grid points between:   (i-ncutoff, j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)

        // === Variables to be used in the trajectory loop ===
        float Nframes       = -1;   // Nframes
        float Nframes_block = -1;   // Nframes in the time_window for block averaging
};

TimeAveragedDensity::TimeAveragedDensity()
{
    registerAnalysisDataset(&data_, "density");
}

void
TimeAveragedDensity::initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes the time-averaged electron density of a group of atoms within a defined box."
        "The density is calculated as a linear combination of Gaussians:\n",
        "rho(R)=sum_{i=1}^{i=Ngaussians} Ai*exp(-Bi*R^2)\n",

        "=== Gaussian coefficients:\n",
        "By default, gtide considers 2-Gaussian coefficients (A and B), derived from Electron Microscopy (citation), for H,C,N,O,P,S.\n",
        "Non-default values can be specified via the -Gauss_coef input file.\n",
        "Format of Gauss_coef.dat:\n",
        "atomtype A_1 A_2 ...A_Ngaussian B1 B2 ... B_NGaussian \n",
        "Precomputed A and B parameters, for 2- and 4-Gaussian sums, from X-ray and EM, as well as for coarse-grained systems, are available at XXXX",
        "How many Gaussians can be changed with -Ngaussians\n",
        "How far each atom spreads in the space can be controlled with -cutoff.\n",

        "=== atom selection:\n",
        "The density map is computed for a specified selection of atoms. Both static and dynamic selections are supported\n",


        "=== Reference map:\n",
        "A reference map can be provided with the -mi option. Different comparison operations between the computed and the reference map are possible (see below)\n",

        "=== Masking:\n",
        "A masking map can be provided with the -mask option. Grid points in the mask, with density smaller than a user-specified"
        "threshold value (-mask_cutoff), will be ignored. In the output maps, the masked points will be labeled as NAN\n",

        "=== Output files:\n",
        "-map_av: time average of density in XPLOR format.\n",
        "-map_sd: time standard deviation of density in XPLOR format.\n",
        "These two maps can be normalized to sigma units and shifted to zero mean for easier visualization (-norm)",
        "and the spatial mean can (or not) include zero-value grid points (-incl_zero).\n",

        "-diff: Difference map:  rho(map_av) - rho(map_reference) in XPLOR format.\n",

        "-mlog: -log [ rho(map_av) ] or -log [ rho(map_av) /rho(map_reference]  if a reference map is provided.\n",
        "-shift_mlog will set the mean value of mlog to zero.\n"

        "-corr: Time-resolved global correlation between the computed and the reference map (xvg format).\n",
        "The correlation is computed in time windows defined by -time_window.\n",
        "-loc_corr: Local correlation between the computed and the reference map in XPLOR format. \n",
        "For a small cubic grid of dimensions (Nslide,Nslide,Nslide), the correlation between map_av ",
        "and the reference maps is calculated, and the resulting value is stored at the center point of this small grid.",
        "The process is repeated by sliding the grid over all the points of the lattice.\n",
        "The global correlation between the maps is written in the last line of the xplor map as the mean.",
        "The rms in the last line of the xplor map will be set to NAN.\n",
        "-rho_xxx: the output maps but in ASCII format. \n",

        "=== grid definition:\n"
        "Parameters are similar parameters as in an XPLOR map:\n",
        "If an input reference is provided, its grid parameters are considered. Otherwise, the following variables must be specified:\n",
        "    -size: Cell grid points (NA,NB,NC)\n",
        "    -box:  --(boxA, boxB,boxC) = crystal cell dimensions (nm,nm,nm)\n",
        "    -min: --(AMIN,BMIN,CMIN) = starting grid point in the considered portion of the map\n",
        "    -max: --(AMAX,BMAX,CMAX) = stopping grid point in the considered portion of the map\n",
        "    Note: The grid resolution is computed as boxi/Ni, independently for each axis i=a,b,c.\n",
        "    Note: The density will be ONLY computed on a sub-lattice going from (AMIN,BMIN,CMIN) to (AMAX,BMAX,CMAX).\n",
        "    Note: The density is absolute with respect to the grid defined by the integer vector Ni. There is no fit procedure.\n",
        "    Nore: If the density is intended to be calculated with respect to other reference group, a previous fitting step (e.g. with trjconv) should precede the gtide calculation.\n",

        "===Examples===:\n",
        "(1) Time-average density of the Protein in a box of size (10 nm,10 nm,10 nm), with its origin at (x,y,z)=(0,0,0) and with a grid resolution of 0.04 nm:\n",
        "gmx tide -f traj.xtc -s conf.gro -box 10 10 10 -size 250 250 250 -min 0 0 0 -max 249 249 249 -select \'group \"Protein\"\' -map_av\n",
        "(2) Same as example 1, but reading non-default Gaussian parameters, for 9 atom types, and 4 Gaussians contributing to the density, and expanding the Gaussians up to 0.8 nm:\n",
        "gmx tide -f traj.xtc -s conf.gro -Gauss_coef Gauss_coef.dat -box 10 10 10 -size 250 250 250 -min 0 0 0 -max 249 249 249 -Ngaussians 4 -cutoff 0.8 -select \'group \"Protein\"\' -map_av\n",
        "(3)Time-average density of the Protein, without normalization, including zeros in the spatial mean, printing as well the difference to a reference map",
        "(lattice parameters taken from the reference map):\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -nonorm -incl_zero -select \'group \"Protein\"\' -map_av -diff\n",
        "(4) Same as example 3, but considering a mask, to ignore grid points which have a density smaller than 1e-5 in the mask:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -mask mask.xplor -nonorm -incl_zero -mask_cutoff 1e-5 -select \'group \"Protein\"\' -map_av -diff\n",
        "(5) Same as example 4, but computing -log[rho(map_av)/rho(map_ref)], shifting the mlog map such that the background spatial mean is zero:",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -mask mask.xplor -mask_cutoff 1e-5 -shift_mlog -select \'group \"Protein\"\' -mlog \n",
        "(6) Time-resolved global correlation of the computed map to the reference map, computing the density in time blocks of 1ps:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -time_window 1.0 -select \'group \"Protein\"\' -corr \n",
        "(7) Local correlation of the computed map to the reference map, sliding a sub-lattice of size (5,5,5) bins through the main lattice:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -Nslide 5 -select \'group \"Protein\"\' -loc_corr \n",
        "(8) Print the density, difference and -log to reference map as a list in ASCII format instead of a XPLOR-format map:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -select \'group \"Protein\"\' -rho -rho_diff -rho_mlog\n",
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("Gauss_coef")
                           .filetype(eftGenericData).inputFile()
                           .store(&gauss_coef_).defaultBasename("Gauss_coef")
                           .description("File with Gaussian coefficients"));

    options->addOption(FileNameOption("mi")
                           .filetype(eftXPLOR).inputFile()
                           .store(&mi_).defaultBasename("mi")
                           .description("Input reference map (optional)"));

    options->addOption(FileNameOption("mask")
                           .filetype(eftXPLOR).inputFile()
                           .store(&mask_).defaultBasename("mask")
                           .description("Input mask map (optional)"));

    options->addOption(FileNameOption("map_av")
                           .filetype(eftXPLOR).outputFile()
                           .store(&map_av_).defaultBasename("map_av")
                           .description("Average density map in X-PLOR format"));

    options->addOption(FileNameOption("map_sd")
                           .filetype(eftXPLOR).outputFile()
                           .store(&map_sd_).defaultBasename("map_sd")
                           .description("St. dev. density map in X-PLOR format"));

    options->addOption(FileNameOption("diff")
                           .filetype(eftXPLOR).outputFile()
                           .store(&diff_).defaultBasename("diff")
                           .description("Difference densiy map in X-PLOR format"));

    options->addOption(FileNameOption("mlog")
                           .filetype(eftXPLOR).outputFile()
                           .store(&mlog_).defaultBasename("mlog")
                           .description("-ln(map_computed/map_reference) in X-PLOR format"));

    options->addOption(FileNameOption("rho")
                           .filetype(eftGenericData).outputFile()
                           .store(&rho_).defaultBasename("rho")
                           .description("av and sd of density in ASCII format"));

    options->addOption(FileNameOption("rho_diff")
                           .filetype(eftGenericData).outputFile()
                           .store(&rho_diff_).defaultBasename("rho_diff")
                           .description("difference density in ASCII format"));

    options->addOption(FileNameOption("rho_mlog")
                           .filetype(eftGenericData).outputFile()
                           .store(&rho_mlog_).defaultBasename("rho_mlog")
                           .description("-ln(map_computed/map_reference) density in ASCII format"));

    options->addOption(FileNameOption("corr").filetype(eftPlot).outputFile()
                           .store(&corr_).defaultBasename("corr")
                           .description("Time-resolved global correlation between computed and reference maps"));

    options->addOption(FileNameOption("loc_corr")
                           .filetype(eftXPLOR).outputFile()
                           .store(&loc_corr_).defaultBasename("loc_corr")
                           .description("Local correlation between computed and reference maps in X-PLOR format"));

    options->addOption(SelectionOption("select")
                           .storeVector(&sel_).required()
                           .description("Group to calculate the density"));

    options->addOption(IntegerOption("size").store(nGridPoints.data()).vector()
                           .description("(NA,NB,NC) = total number of grid points along each cell edge size"));

    options->addOption(IntegerOption("min").store(Nmin.data()).vector()
                           .description("(Amin,Bmin,Cmin) = starting grid points in the considered portion of the map"));

    options->addOption(IntegerOption("max").store(Nmax.data()).vector()
                           .description("(Amax,Bmax,Cmax) = stopping grid points in the considered portion of the map\nDensity will be calculated within min-max grid points"));

    options->addOption(RealOption("box").store(box).vector()
                           .description("crystal cell dimensions. Grid spacing computed as box/N for each axis"));

    options->addOption(FloatOption("cutoff").store(&cutoff)
                           .description("Cutoff distance to compute the density around a certain atom (nm)"));

    options->addOption(IntegerOption("Ngaussians").store(&Ngaussians)
                           .description("Number of Gaussians to use in the sum (e.g. 2,4, or 5)"));

    options->addOption(FloatOption("mask_cutoff").store(&mask_cutoff)
                           .description("Points with density smaller than mask_cutoff will be masked according to mask map (-mask)"));

    options->addOption(BooleanOption("norm").store(&Norm)
                           .description("Normalize the output maps. If yes, map in stdev units && shifted such that mean=0"));

    options->addOption(BooleanOption("incl_zero").store(&incl_zero)
                           .description("Include Zero-density points in the calculation of the spatial mean and stdev, over the grid"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

    options->addOption(BooleanOption("shift_mlog").store(&shift_mlog)
                           .description("remove mean value from mlog map"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);


    options->addOption(FloatOption("time_window").store(&time_window)
                           .description("Time window for block-averaging in time-resolved global correlation calculation(ps)"));

    options->addOption(IntegerOption("Nslide").store(&Nslide)
                           .description("Size of the sliding grid for local correlation calculation (Must be an odd number)"));
}

void
TimeAveragedDensity::initAnalysis(const gmx_unused TrajectoryAnalysisSettings  &settings,
                                  const TopologyInformation                    &top)
{

    data_.setColumnCount(0, sel_.size());
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    // === GAUSSIAN COEFFICIENTS BASED ON ATOM TYPES===
    fprintf(stdout, "Considering %d Gaussians \n", Ngaussians);

    set_gaussian_coeff(gauss_coef_, &A, &B, Ngaussians);
    attypeString = attypeStrings_from_top(*(top.topology()));

//=== INPUT REFERENCE MAP ====

    // === CREATE THE LATTICE ===
    if (!mi_.empty())
    {
        /* === Read input file (map) === */
        printf("=== Reading reference map\n");
        read_xplor_3dmap(mi_.c_str(), &map_ref );
        printf("Setting map parameters equal to those of the reference map. -size, -min, -max, and -box options ignored\n");
        nGridPoints = map_ref.grid; // cell grid points
        bin         = map_ref.map_dim;
        Nmin        = map_ref.origin;
        for (int i = 0; i < DIM; i++)
        {
            box[i] = map_ref.cell[i];                 // cell size
            res[i] = map_ref.cell[i]/map_ref.grid[i]; //resolution
            L0[i]  = map_ref.origin[i]*res[i];        // origin of the lattice
        }
    }
    else            // no reference map provided, read variables the -size, -min, -min , -box options
    {
        for (int i = 0; i < DIM; i++)
        {
            //number of bins (working lattice)
            bin[i] = Nmax[i]-Nmin[i]+1;         //Amax-Amin+1

            //resolution
            res[i] = box[i]/nGridPoints[i];

            // origin of the lattice
            L0[i] = Nmin[i]*res[i];
        }
    }

    // ====  initialize the maps =====
    std::string title {
        "map"
    };
    std::array<float, 2*DIM> cell    = {{static_cast<float>(box[XX]), static_cast<float>(box[YY]), static_cast<float>(box[ZZ]), 90., 90., 90.}};
    //N  Number of grid points within box
    //  MIN[i]; // Origin of the working lattice
    std::array<int, DIM> axes_order = {{0, 1, 2}};                  // x, y, z (guess) but  not used in tide anyways
    // bin = map_dim is the working lattice
    // spacegroup=1; // this will be assigned arbitrarily. It won't be used in gtide anyways

    // minimum and maximum density values
    double min  = 0;
    double max  = 0;
    float  mean = 0;
    float  rms  = 0;

    //skew options will be ignored (0)
    std::array<float, DIM*DIM> skew_mat;
    std::array<float, DIM>     skew_trans;
    std::fill(std::begin(skew_mat), std::end(skew_mat), 0.f);
    std::fill(std::begin(skew_trans), std::end(skew_trans), 0.f);

    printf("=== Initial MD map parameters:\n");
    init_map(&map_av, title, 2, cell, nGridPoints, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans);
    printMapInfo(map_av);
    //allocate memory for the lattices
    allocate_density_grid(map_av.map_dim, &(map_av.vox));
    map_inst     = map_av;
    map_sd       = map_av;
    map_diff     = map_av;
    map_mlog     = map_av;
    map_block    = map_av;
    map_loc_corr = map_av;

    /* === Mask map === */
    if (!mask_.empty())  // read map if provided
    {
        printf("=== Reading mask map\n");
        read_xplor_3dmap(mask_.c_str(), &map_mask );
        set_mask(&map_mask, mask_cutoff);
    }
    else
    {
        map_mask = map_av;
        std::fill(std::begin(map_mask.vox), std::end(map_mask.vox), 1.f);
    }

    // === Definition of the cutoff and other gaussian-density related variables
    ncutoff[XX] = floor(cutoff/res[XX]);
    ncutoff[YY] = floor(cutoff/res[YY]);
    ncutoff[ZZ] = floor(cutoff/res[ZZ]);
    // Density will be increment in the grid points between:   (i-ncutoff, j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)

    // === change the cutoff to Angstrom^2 ===
    cutoff *= NM2AA;

    // === NFRAMES = 0 ===
    Nframes       = 0;
    Nframes_block = 0;

    if (!corr_.empty())
    {
        pfile_corr = fopen(corr_.c_str(), "w");

        // header for xmgrace
        // To do: print properly the header with the addModule option
        fprintf(pfile_corr, "@    title \"Time-resolved correlation between computed and reference maps\"\n");
        fprintf(pfile_corr, "@    subtitle  \"Time-window for block-averaging: %f ps\"\n", time_window);
        fprintf(pfile_corr, "@    xaxis  label \"Time\"\n");
        fprintf(pfile_corr, "@    yaxis  label \"Correlation\"\n");

    }
}

void
TimeAveragedDensity::analyzeFrame(int frnr, const t_trxframe &fr, gmx_unused t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata)
{

    // === STUFF TO READ THE TRAJECTORY AND ATOM SELECTION DATA ===
    AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    dh.startFrame(frnr, fr.time);
    size_t                     group = 0; // selection-group index (only one considered thus it is 0)
    const Selection           &sel   = pdata->parallelSelection(sel_[group]);

    // === ACCUMMULATE DENSITY===
    accumulate_denstiy(&map_av, &map_inst, &map_sd, sel, res, L0, ncutoff, cutoff, A, B, attypeString);

    // === ACCUMULATE BLOCK AVERAGE ===
    addToDensity(&map_block, map_inst);

    // === Increment the number of frames variables ===
    Nframes++;
    Nframes_block++;

    // === calculate and print global correlation if requested ===
    if (!corr_.empty())
    {
        if (fr.time/time_window-floor(fr.time/time_window) == 0 && fr.time > 0)
        {

            divideDensityByValue(&map_block, Nframes_block );
            float cc = global_correlation_coeff(&map_ref, &map_block, &map_mask);
            fprintf(pfile_corr, "%f %f\n", fr.time, cc);
            setDensityToZero(&map_block);
            Nframes_block = 0;
        }
    }
    dh.finishFrame();
}


void
TimeAveragedDensity::finishAnalysis(int /*nframes*/)
{

    fprintf(stdout, "Nframes read: %f\n", Nframes);

    // calculate timed-average and timed-stdev of the cumulated map
    divideDensityByValue ( &map_av, Nframes);
    sd_map (map_av,  &map_sd,  Nframes);

    // calculate spatial average and spatial stdev
    do_spatial_average_maps( &map_av, incl_zero, &map_mask );
    do_spatial_average_maps( &map_sd, incl_zero, &map_mask );

    // calculate the difference to the reference map (if requested)
    if (!diff_.empty())
    {
        subtract_maps(&map_ref, &map_av, &map_diff, &map_mask ); // map_diff = map_av - map_ref
        do_spatial_average_maps( &map_diff, true, &map_mask);    // do include zero in the calculation of mean and rms
    }

    // calculate minus log (if requested)
    if (!mlog_.empty())
    {
        // if there is a ref map
        if (!mi_.empty())
        {
            minus_log_division_maps(&map_ref, &map_av, &map_mlog, &map_mask);               // map_mlog = - ln ( map_av /  map_ref )
        }
        // no ref map
        else
        {
            minus_log_map(&map_av, &map_mlog, &map_mask);      // map_mlog = -ln (map_av)
        }
        do_spatial_average_maps( &map_mlog, true, &map_mask);  // do include zero in the calculation of mean and rms

        // shift map if requested
        if (shift_mlog)
        {
            shift_map(&map_mlog, -map_mlog.mean);
            do_spatial_average_maps( &map_mlog, true, &map_mask);
        }
    }

// print global correlation if there is a reference map
    if (!mi_.empty())
    {
        float cc = global_correlation_coeff(&map_ref, &map_av, &map_mask);
        printf("The correlation coefficient between computed and reference maps is %f\n", cc);
        // calculate local correlation if requested
        if (!loc_corr_.empty())
        {
            local_correlation_coeff(&map_ref, &map_av, &map_mask, &map_loc_corr, Nslide);
            map_loc_corr.mean = cc;  // mean will be the global correlation coefficient
            map_loc_corr.rms  = nanf("");
        }
    }

    // shift the map by an amount equal to the average
    if (Norm)
    {
        shift_map(&map_av, -map_av.mean);
        shift_map(&map_sd, -map_sd.mean);

        // Normalize the map by an amount equal to the standard deviation
        divideDensityByValue(&map_av, map_av.rms);
        divideDensityByValue(&map_sd, map_sd.rms);
        map_av.mean = -1.0; // fixed values to indicate stdev units around the mean
        map_sd.mean = -1.0;
        map_av.rms  = -1.0;
        map_sd.rms  = -1.0;
    }
}

void
TimeAveragedDensity::writeOutput()
{

    // print the average to file 3dmap_average
    if (!map_av_.empty())
    {
        print_to3dmap(map_av_.c_str(), &map_av, &map_mask);
    }

    // print the standard deviation to file 3dmap_sigma
    if ( (!map_sd_.empty()) && (Nframes > 1) )
    {
        print_to3dmap(map_sd_.c_str(), &map_sd, &map_mask);
    }

    // print the difference map to file diff
    if (!diff_.empty())
    {
        print_to3dmap(diff_.c_str(), &map_diff, &map_mask);
    }

    // print the difference map to file diff
    if (!mlog_.empty())
    {
        print_to3dmap(mlog_.c_str(), &map_mlog, &map_mask);
    }

    //print to rho.dat file
    if (!rho_.empty())
    {
        print_rho_sigma(rho_.c_str(), res, &map_av, &map_sd );
    }

    //print to rho_diff.dat file
    if (!rho_diff_.empty())
    {
        print_rho_sigma(rho_diff_.c_str(), res, &map_diff, &map_diff );
    }


    //print to rho_mlog.dat file
    if (!rho_mlog_.empty())
    {
        print_rho_sigma(rho_mlog_.c_str(), res, &map_mlog, &map_mlog );
    }

    // close time-averaged correlation printing output file
    if (!corr_.empty())
    {
        fclose (pfile_corr);
    }

    // print the local correlation map
    if (!loc_corr_.empty())
    {
        print_to3dmap(loc_corr_.c_str(), &map_loc_corr, &map_mask);
    }

}

}       // namespace


const char TimeAveragedDensityInfo::name[]             = "tide";
const char TimeAveragedDensityInfo::shortDescription[] = "Calculate TIme averaged electron DEnsities";

TrajectoryAnalysisModulePointer TimeAveragedDensityInfo::create()
{
    return TrajectoryAnalysisModulePointer(new TimeAveragedDensity);
}

} // namespace analysismodules

} // namespace gmx
