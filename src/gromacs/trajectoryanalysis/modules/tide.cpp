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

//! \brief Gets the index in an uni-dimensional array from the three coordinates of a three-dimensional array
int index_array(int i, int j, int k, int Nx, int Ny, gmx_unused int Nz)
{
// The unidimensional array contains Nx*Ny*Nz elements:  0 ... Nx*Ny*Nz-1
// The element in the unidimensional array is given by
    return k*(Nx*Ny) + j*Nx + i;
}

/*! \brief Compare variables relevant for tide (cell,grid, origin and map_dim) of map1 and map2.
 *
 *  \return true if their dimensions are equal and false otherwise===
 */
bool compare_maps(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2 )
{
    float cell_tolerance = 0.0001;   // diff smaller than that will be treated as zero.

    bool  same = true;
    for (int i = 0; i < 2*DIM; i++)
    {
        if ( (map1->cell[i] - map2->cell[i])*(map1->cell[i] - map2->cell[i]) > cell_tolerance*cell_tolerance)
        {
            same = false; printf("cell %f %f\n", map1->cell[i], map2->cell[i]);
        }
    }
    for (int i = 0; i < DIM; i++)
    {
        if (map1->grid[i] != map2->grid[i])
        {
            same = false;  printf("grid %d %d\n", map1->grid[i], map2->grid[i]);
        }
    }
    for (int i = 0; i < DIM; i++)
    {
        if (map1->origin[i] != map2->origin[i])
        {
            same = false; printf("origin %d %d\n", map1->origin[i], map2->origin[i]);
        }
    }
    for (int i = 0; i < DIM; i++)
    {
        if (map1->map_dim[i] != map2->map_dim[i])
        {
            same = false; printf("dimensions %d %d\n", map1->map_dim[i], map2->map_dim[i]);
        }
    }

    return same;
}



//! \brief set mask
void set_mask( std::string map3d_mask_, gmx_unused t_mapdata *mask, float mask_cutoff)
{
    int l;


    for (l = 0; l < mask->map_dim[XX]*mask->map_dim[YY]*mask->map_dim[ZZ]; l++)
    {

        if (!map3d_mask_.empty())
        {

            // set to zero for the values for which mask is zero:
            if (mask->vox[l] <= mask_cutoff)
            {
                mask->vox[l] = 0;
            }
            else
            {
                mask->vox[l] = 1;
            }

        }

        else               // there is no file then set all mask values to 1
        {
            mask->vox[l] = 1;
        }

    }
}

//! \brief Subtract maps: rho(map3) = rho(map2)-rho(map1)
void subtract_maps(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *map3, gmx_unused t_mapdata *mask)
{
    //maps have to have the same dimensions

    bool same_dims = false;
    if (compare_maps(map1, map2) && compare_maps(map1, map3) && compare_maps(map1, mask) )
    {
        same_dims = true;
    }
    else
    {
        same_dims = false;
    }

    // do subtraction
    if (same_dims == true)
    {
        for (int l = 0; l < map1->map_dim[XX]*map1->map_dim[YY]*map1->map_dim[ZZ]; l++)
        {
            if (mask->vox[l] > 0.5)
            {
                map3->vox[l] = map2->vox[l]-map1->vox[l];
            }
            else
            {
                map3->vox[l] = nanf("");     // masked point is not considered
            }
        }
    }
    else            // map dimensions are different
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
}


//! \brief log of map: rho(map2) = - ln ( rho(map) )  PSEUDO FREE ENERGY
void minus_log_map(gmx_unused t_mapdata *map, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *mask)
{
    bool same_dims = false;
    if (compare_maps(map, map2) && compare_maps(map, mask) )
    {
        same_dims = true;
    }
    else
    {
        same_dims = false;
    }

    // do operation
    if (same_dims == true)
    {


        for (int l = 0; l < map->map_dim[XX]*map->map_dim[YY]*map->map_dim[ZZ]; l++)
        {
            if (mask->vox[l] > 0.5)            // discard masked points

            {
                if (map->vox[l] != 0)
                {
                    map2->vox[l] = -log( map->vox[l] );
                }
                else
                {
                    gmx_fatal(FARGS, "Density at point %d is zero in the map. You could discard this point from the calculation by using -mask.\n", l);
                }
            }
            else
            {
                map2->vox[l] = nanf("");             // masked point is not considered

            }
        }
    }
    else               // map dimensions are different
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
}

//! \brief log of division of maps: rho(map3) = - ln ( rho(map2) / rho(map1) )  PSEUDO FREE ENERGY
void minus_log_division_maps(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *map3, gmx_unused t_mapdata *mask)
{
    bool same_dims = false;
    if (compare_maps(map1, map2) && compare_maps(map1, map3) && compare_maps(map1, mask) )
    {
        same_dims = true;
    }
    else
    {
        same_dims = false;
    }


    // do operation
    if (same_dims == true)
    {


        for (int l = 0; l < map1->map_dim[XX]*map1->map_dim[YY]*map1->map_dim[ZZ]; l++)
        {
            if (mask->vox[l] > 0.5)        // discard masked points

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
    else            // map dimensions are different
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
}

//! \brief Normalize map: rho(map) = rho(map)/ N  , with N a normalization scalar
void normalize_map(gmx_unused t_mapdata *map, float N )
{
    int l;
    for (l = 0; l < map->map_dim[XX]*map->map_dim[YY]*map->map_dim[ZZ]; l++)
    {
        map->vox[l] = map->vox[l]/N;
    }
}

//! \brief Shift map: rho(map) = rho(map) + S  , with S a shift scalar value
void shift_map(gmx_unused t_mapdata *map, float S )
{
    int l;
    for (l = 0; l < map->map_dim[XX]*map->map_dim[YY]*map->map_dim[ZZ]; l++)
    {
        map->vox[l] = map->vox[l] + S;
    }
}


//! \brief Compute spatial average over the map: <rho(map)> over all cell positions
void do_spatial_average_maps(gmx_unused t_mapdata *map, bool incl_zero, gmx_unused t_mapdata *mask  )
{
    map->mean = 0.0;
    map->rms  = 0.0;
    float nsamples      = 0;

    if (compare_maps(map, mask))          // map and mask must have same dimensions
    {
        for (int l = 0; l < map->map_dim[XX]*map->map_dim[YY]*map->map_dim[ZZ]; l++)
        {
            if (mask->vox[l] > 0.5)                // 0=masked point , 1=unmasked point to be considered
            {
                if (map->vox[l] != 0 || incl_zero) // zero density point
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

}

//! \brief TODO
float calc_sum_rhoA_rhoB_mask (gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *mask)
{

    bool same_dims = false;
    if (compare_maps(map1, map2) && compare_maps(map1, mask) )
    {
        same_dims = true;
    }
    else
    {
        same_dims = false;
    }

    double sum_ijk = 0.0;
    int    l;

    if (same_dims)
    {
        for (l = 0; l < map1->map_dim[XX]*map1->map_dim[YY]*map1->map_dim[ZZ]; l++)
        {
            if (mask->vox[l] > 0.5)        // 0=masked point , 1=unmasked point to be considered
            {
                sum_ijk += map1->vox[l]*map2->vox[l];
            }
        }
        return sum_ijk;
    }
    else
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
}

//! \brief TODO
float global_correlation_coeff(gmx_unused t_mapdata *map1, gmx_unused t_mapdata *map2, gmx_unused t_mapdata *mask)
{
    bool same_dims = false;
    if (compare_maps(map1, map2) && compare_maps(map1, mask) )
    {
        same_dims = true;
    }
    else
    {
        same_dims = false;
    }

    if (same_dims)
    {
        float     numerator, denominator, cc;
        numerator   = calc_sum_rhoA_rhoB_mask(map1, map2, mask);
        denominator = sqrt(calc_sum_rhoA_rhoB_mask(map2, map2, mask) * calc_sum_rhoA_rhoB_mask(map1, map1, mask));
        cc          = numerator / denominator;
        return cc;
    }
    else
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
}



//! \brief Calculate the local correlation
void local_correlation_coeff(t_mapdata *map1, t_mapdata *map2, t_mapdata *mask, t_mapdata *loc_corr, int Nodd)
{
    int Nhalf = int(Nodd/2);    // center point of the sliding grid (e.g. if Nodd=5 ; then Nhalf=2, Note that it starts counting from 0)

    // This can only work for maps with equal dimensions
    if (!compare_maps(map1, map2) || !compare_maps(map1, mask) || !compare_maps(map1, loc_corr) )
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }

    // Initialize the output loc_corr map, setting the borders and masked points, with non-assigned NAN values
    for (int l = 0; l < map1->map_dim[XX]*map1->map_dim[YY]*map1->map_dim[ZZ]; l++)
    {
        loc_corr->vox[l] = nanf("");
    }

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
                            int id_aux = index_array(i_sli, j_sli, k_sli, map1->map_dim[XX], map1->map_dim[YY], map1->map_dim[ZZ]);
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
                int id_aux            = index_array(i, j, k, map1->map_dim[XX], map1->map_dim[YY], map1->map_dim[ZZ]);
                loc_corr->vox[id_aux] = cc;
            }
        }
    }
}

//! \brief Set the Gaussian coefficients A and B for the six atom types H, C, N, O, P, and S
void set_gaussian_coeff( std::string gauss_coeff_, std::vector<float> *A, std::vector<float> *B, std::vector<char> *attypes,  int Ngaussians, int Nattypes  )
{

    attypes->clear();
    A->clear();
    B->clear();

    // if there is an input Gaussian_coeff file, consider the values in there:
    if (!gauss_coeff_.empty())
    {

        fprintf(stdout, "Reading file with Gaussian coefficients\n");

        FILE *pfile;
        pfile = fopen(gauss_coeff_.c_str(), "r");

        // loop over the atoms
        for (int n = 0; n < Nattypes; n++)
        {

            // read the next parameter line
            char currentAtomType[1];
            fscanf(pfile, "%s", currentAtomType); // atomtype
            attypes->push_back(*currentAtomType);

            // A parameters
            for (int col = 0; col < Ngaussians; col++)
            {
                float currentAFactor;
                fscanf(pfile, "%e", &currentAFactor);
                A->push_back(currentAFactor);
            }

            // B parameters
            for (int col = 0; col < Ngaussians; col++)
            {
                float currentBFactor;
                fscanf(pfile, "%e", &currentBFactor);
                B->push_back(currentBFactor);
            }
        }
        fclose (pfile);
    }
    else                                // take the defaults (2 Gaussians from cryoEM)
    {
        A->resize(Nattypes*Ngaussians); // a factors
        B->resize(A->size());           // b factors
        fprintf(stdout, "Taking default Gaussian coefficients from cryo EM for sums with 2 Gaussians \n");

        *attypes = {'H', 'C', 'N', 'O', 'P', 'S'};

        //g=0               g=1                 g=0                 g=1
        (*A)[0]  = 1.242258e-01; (*A)[1] = 1.133291e+00; (*B)[0] = 1.718695e+00; (*B)[1] = 9.919200e+00;               //H  n=0
        (*A)[2]  = 5.948015e-01; (*A)[3] = 5.657122e+00; (*B)[2] = 1.700190e+00; (*B)[3] = 1.041647e+01;               //C  n=1
        (*A)[4]  = 7.587542e-01; (*A)[5] = 7.705668e+00; (*B)[4] = 2.236738e+00; (*B)[5] = 1.342803e+01;               //N  n=2
        (*A)[6]  = 9.061284e-01; (*A)[7] = 9.848492e+00; (*B)[6] = 2.786056e+00; (*B)[7] = 1.644934e+01;               //O  n=3
        (*A)[8]  = 1.056651e+00; (*A)[9] = 1.083431e+01; (*B)[8] = 1.390085e+00; (*B)[9] = 1.081600e+01;               //P  n=4
        (*A)[10] = 1.335159e+00; (*A)[11] = 1.271311e+01; (*B)[10] = 1.711987e+00; (*B)[11] = 1.226038e+01;            //S  n=5

    }
}

//! \brief Output integer position depending on atom type
int attype_index( char attype, const std::vector<char> &attypes, int Nattypes )
{
    int position = std::distance(std::begin(attypes), std::find(std::begin(attypes), std::end(attypes), attype));  // this will be the position number in the attypes array (e.g. H=0 , C=1, N=2 , O=3 ; P=4  ; S=5)

    if (position >= 0 && position < Nattypes)
    {
        return position;
    }
    else
    {
        gmx_fatal(FARGS, "atom type %.1s was not found in Gaussian-coefficient attype data base\nInclude this atom in the data base (-Gauss_coef)\nor exclude this atom from input structure and input trajectory", attype);
    }
}

//! \brief  Output integer position in the array of atomtypes based on first character of atomname
std::vector<int> attype_positions_from_top(const std::vector<char> &attypes, int Natoms, int Nattypes, const TopologyInformation &top)
{
    std::vector<int> attype_position;
    for (int n = 0; n < Natoms; n++)
    {
        char * atomname   = *(top.topology()->atoms.atomname[n]);
        attype_position.push_back( attype_index(atomname[0], attypes, Nattypes));
    }
    return attype_position;
}


/*! \brief Outputs the density map to text file using XPLOR format
 *
 * For the format see: http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
void print_to3dmap(const char *N_file, gmx_unused t_mapdata *map, gmx_unused t_mapdata *mask )
{
    if (N_file == nullptr)
    {
        fprintf(stdout, "No filename given for X-PLOR output file, not outputting densities\n");
        return;
    }


    if (compare_maps(map, mask))
    {

        FILE *fp = gmx_ffopen(N_file, "w");

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
        //  nBOX = total number of grid points along the a,b, and c cell edges (of the simulation box)

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

        // === Line 6 : BOX size (Angstrom) and angles ===
        fprintf(fp, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n", map->cell[XX]*NM2AA, map->cell[YY]*NM2AA, map->cell[ZZ]*NM2AA, map->cell[3], map->cell[4], map->cell[5]   );

        // === Line 7 : string ===
        fprintf(fp, "ZYX\n");

        for (int k = 0; k < map->map_dim[ZZ]; k++) // loop over z
        {                                          // Print the k value
            fprintf(fp, "       %d\n", k);

            // variable to monitor the number of columns printed
            int ncols = 0;

            for (int j = 0; j < map->map_dim[YY]; j++)       // loop over y

            {
                for (int i = 0; i < map->map_dim[XX]; i++)       // loop over x

                {
                    int index_aux = index_array(i, j, k, map->map_dim[XX], map->map_dim[YY], map->map_dim[ZZ]);

                    if (mask->vox[index_aux] > 0.5)
                    {
                        // Print the density values in 6 columns
                        fprintf(fp, "%12.5E", map->vox[index_aux] );
                    }
                    else            // print nan to indicate that this point is masked
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
    else
    {
        gmx_fatal(FARGS, "Map and mask dimensions are different.");

    }
}

//! \brief Init map function
void init_map(gmx_unused t_mapdata *map, char *title, int datamode, float cell[], int  grid[], int origin[], int axes_order[], int map_dim[], int spacegroup, double min, double max, float mean, float rms, float skew_mat[], float skew_trans[], bool printinfo)
{
    map->title    = title;
    map->datamode = datamode; //stored as Reals
    for (int i = 0; i < 2*DIM; i++)
    {
        map->cell[i] = cell[i];                       // cell 0...2 are the box sizes // cell 3...5 are the angles
    }
    for (int i = 0; i < DIM; i++)
    {
        map->grid[i] = grid[i];                     // Number of grid points within BOX
    }
    for (int i = 0; i < DIM; i++)
    {
        map->origin[i] = origin[i];                     // Origin of the working lattice
    }
    for (int i = 0; i < DIM; i++)
    {
        map->axes_order[i] = axes_order[i];                     // x, y, z (guess)
    }
    for (int i = 0; i < DIM; i++)
    {
        map->map_dim[i] = map_dim[i]; // the working lattice
    }
    map->spacegroup = spacegroup;     // this will be assigned arbitrarily. It won't be used in gtide anyways

    // minimum and maximum density values
    map->min = min;
    map->max = max;

    // average and stdev
    map->mean = mean;
    map->rms  = rms;


    // skew options will be ignored (0)
    for (int i = 0; i < 9; i++)
    {
        map->skew_mat[i] = skew_mat[i];
    }
    for (int i = 0; i < 3; i++)
    {
        map->skew_trans[i] = skew_trans[i];
    }


    // print crosscheck
    if (printinfo)
    {
        printf( "\n");
        printf( "Header:\n");
        printf( "Title      : %s\n", map->title);
        printf( "Datamode   : %u\n", map->datamode);
        printf( "Cell       : %g A  %g A  %g A  %g %g %g\n", map->cell[0], map->cell[1], map->cell[2], map->cell[3], map->cell[4], map->cell[5]);
        printf( "Grid       : %d x %d x %d\n", map->grid[XX], map->grid[YY], map->grid[ZZ]);
        printf( "Origin     : %d A  %d A  %d A\n", map->origin[XX], map->origin[YY], map->origin[ZZ]);
        printf( "Axes order : %d %d %d\n", map->axes_order[0], map->axes_order[1], map->axes_order[2]);
        printf( "Map_dim    : %d %d %d\n", map->map_dim[0], map->map_dim[1], map->map_dim[2]);
        printf( "Spacegroup : %d\n", map->spacegroup);
        printf( "Min        : %g\n", map->min);
        printf( "Max        : %g\n", map->max);
        printf( "Mean       : %g\n", map->mean);
        printf( "RMS        : %g\n", map->rms);
        printf( "Skew_trans : %g %g %g\n", map->skew_trans[0], map->skew_trans[1], map->skew_trans[2]);
        printf( "Skew_mat   : ");
        for (int i = 0; i < 9; i++)
        {
            printf( "%g ", map->skew_mat[i]);
        }
        printf( "\n\n");
    }
}

/*!  \brief Read an input density map using XPLOR format
 *
 * For the format see: http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
void read_xplor_3dmap(const char *N_file, gmx_unused t_mapdata *ref_map )
{
    if (N_file == nullptr)
    {
        fprintf(stdout, "No filename given for X-PLOR input file\n");
        return;
    }


    FILE *fp = gmx_ffopen(N_file, "r");

    //  starting and stopping grid points along each cell edge in the portion of the map that is written
    rvec   BOX, angle;
    ivec   N, MIN = {0, 0, 0}, MAX = {0, 0, 0}, bin = {0, 0, 0};
    float *lattice;  // one-D array to save the lattice

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
        // === Line 6 : BOX size (in Angstrom) and angles (degree)===
        if (row == row_ZYX-1)
        {
            fscanf(fp, "%12f", &BOX[XX]);
            fscanf(fp, "%12f", &BOX[YY]);
            fscanf(fp, "%12f", &BOX[ZZ]);

            fscanf(fp, "%12f", &angle[XX]);
            fscanf(fp, "%12f", &angle[YY]);
            fscanf(fp, "%12f", &angle[ZZ]);


            // having read N,MIN,MAX, BOX and angle, it is possible to compute bin
            for (int i = 0; i < DIM; i++)
            {

                // set box in nm
                BOX[i] /= NM2AA;

                //number of bins (working lattice)
                bin[i] = MAX[i]-MIN[i]+1;                    //Amax-Amin+1
            }
            // allocate space for the density one-dimensional array
            snew (lattice, bin[XX] * bin[YY] * bin[ZZ] );                        //lattice
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
            printf("k!=k_read: %d %d\n", k, k_read); break;
        }

        else
        {

            for (int j = 0; j < bin[YY]; j++)                           // loop over y
            {
                for (int i = 0; i < bin[XX]; i++)                       // loop over x
                {
                    int index_aux = index_array(i, j, k, bin[XX], bin[YY], bin[ZZ]);
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
    float   av, stdev;
    fscanf(fp, "%f", &av);
    fscanf(fp, "%f", &stdev);

    fclose (fp);



//=== assign the information read in the xplor file to the map ===
    char  title[] = {'r', 'e', 'f'};

    float cell[6];
    for (int i = 0; i < DIM; i++)
    {
        cell[i] = BOX[i];                        // cell 0...2 are the box sizes
    }
    for (int i = 3; i < 2*DIM; i++)
    {
        cell[i] = angle[i-3];                          // cell 3...5 are the angles
    }
    //N  Number of grid points within BOX
    //  MIN[i]; // Origin of the working lattice
    int axes_order[3] = {0, 1, 2};              // x, y, z (guess) but  not used in tide anyways
    // bin = map_dim is the working lattice
    // spacegroup=1; // this will be assigned arbitrarily. It won't be used in gtide anyways

    // minimum and maximum density values
    double min = 1e6, max = -1e6;
    for (int k = 0; k < bin[ZZ]; k++)
    {
        for (int j = 0; j < bin[YY]; j++)
        {
            for (int i = 0; i < bin[XX]; i++)
            {
                int index_aux = index_array(i, j, k, bin[XX], bin[YY], bin[ZZ]);
                if (lattice[index_aux] < min)
                {
                    min = lattice[index_aux];
                }
                if (lattice[index_aux] > max)
                {
                    max = lattice[index_aux];
                }
            }
        }
    }
    // average and stdev read them from last line in the input ref map

    //skew options will be ignored (0)
    float skew_mat[9];
    float skew_trans[3];
    for (int i = 0; i < 9; i++)
    {
        skew_mat[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        skew_trans[i] = 0;
    }

    init_map(ref_map, title, 2, cell, N, MIN,  axes_order, bin, 1, min, max, av, stdev, skew_mat, skew_trans, 0);


    /* ==== 1d continuous array for 3d voxel density  ====      */
    allocate_density_grid(ref_map->map_dim, &ref_map->vox);

    // copy lattice into vox
    for (int k = 0; k < bin[ZZ]; k++)
    {
        for (int j = 0; j < bin[YY]; j++)
        {
            for (int i = 0; i < bin[XX]; i++)
            {
                int index_aux               = index_array(i, j, k, bin[XX], bin[YY], bin[ZZ]);
                ref_map->vox[index_aux] = lattice[index_aux];
            }
        }
    }


}


//! \brief Accummulate density every time step
void accumulate_denstiy( gmx_unused t_mapdata *map_av, gmx_unused t_mapdata *map_inst, gmx_unused t_mapdata *map_sd, int natoms, Selection sel,  float res[], float L0[], int ncutoff[], float cutoff, int Ngaussians, const std::vector<float> &Ai, const std::vector<float> &Bi,  std::vector<int> attype_position)
{

    // Coordinate vector with respect to the simulation box
    rvec R_trans;

    // vector coordinates of the (nx,ny,nz) grid point
    rvec grid_vec;

    // Distance between the atom and the grid_vector:  R_atom_grid= R_trans - grid_vec
    rvec  R_atom_grid;

    float dist_atom_grid;

    //=== INITIALIZE THE INSTANTANEOUS MAP LATTICE===

    for (int l = 0; l < map_av->map_dim[XX]*map_av->map_dim[YY]*map_av->map_dim[ZZ]; l++)
    {
        map_inst->vox[l] = 0.0;
    }


    // === Count the visit of each one of the surrounding atoms (i.e. lipids) in the proper grid of the lattice

    // === Loop over the surrounding atoms (i.e. lipids )
    for (int n = 0; n < natoms; n++)
    {

        int atomindex = sel.atomIndices()[n];
        int position  = attype_position[atomindex]; // e.g. for default Gaussian-coeff array: H=0 , C=1 , N=2 , O=3, P=4, S=5
        //printf("%d %d %d\n",n,atomindex,position);


        // Coordinate vector with respect to the simulation box
        SelectionPosition p = sel.position(n);   // read position selection of n-th atom

        //	 printf("R at t=%8.3f : x= %8.5f y= %8.5f z= %8.5f\n",0.0,  R[XX] , R[YY] , R[ZZ] );

        // Change of system of coordinates to one where the origin is the origin of the lattice
        // R_trans = coordinates of the atom with respect to the lattice
        // R_trans= R - L0
        // L0 origin of the lattice coordinates
        rvec_sub( p.x(), L0, R_trans );      // R_trans = R - L0

        // === Discretize the Rtrans to locate it in the lattice
        const int i = floor( R_trans[XX] / res[XX]);
        const int j = floor( R_trans[YY] / res[YY]);
        const int k = floor( R_trans[ZZ] / res[ZZ]);

        // === only consider the atoms located in the lattice ===
        // 0 <= i <= bin[XX]
        // 0 <= j <= bin[YY]
        // 0 <= k <= bin[ZZ]

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
                        // 0 <= nx <= map_av->map_dim[XX]
                        // 0 <= ny <= map_av->map_dim[YY]
                        // 0 <= nz <= map_av->map_dim[ZZ]

                        if (nx >= 0 && nx < map_av->map_dim[XX] && ny >= 0 && ny < map_av->map_dim[YY] && nz >= 0 && nz < map_av->map_dim[ZZ])
                        {
                            //		printf("nx=%d ny=%d nz=%d\n",nx,ny,nz);

                            // vector coordinates of the (nx,ny,nz) grid point
                            grid_vec[XX] = nx*res[XX];
                            grid_vec[YY] = ny*res[YY];
                            grid_vec[ZZ] = nz*res[ZZ];

                            //printf("grid vec: (%f, %f, %f)\n", grid_vec[XX],grid_vec[YY] , grid_vec[ZZ] );


                            // Distance between the atom and the grid_vector:
                            rvec_sub( R_trans, grid_vec, R_atom_grid ); // R_atom_grid  = R_trans - grid_vec


                            // Calculate the norm of R_atom_grid	(Angstrom) to coincide with the ai, bi and c factors
                            dist_atom_grid = NM2AA * norm(  R_atom_grid );


                            // === Only contribute to the density if the distance is smaller than the cutoff
                            if (dist_atom_grid < cutoff)
                            {




                                // === Gaussian contribution (two or four gaussians depending on the atom) to the density at this grid point
                                real density_contrib = 0.0;


                                // contributions to the atomic potentials according to  equation (6) in J. Electron Microscopy, 56(4):131-140 (2007)
                                //gaussian contribution
                                for (int g = 0; g < Ngaussians; g++)
                                {
                                    int index_oneD_aux   = index_array(g, position, 0, Ngaussians, natoms, 0); // turn (position,g) into a one-dimensional array.
                                    density_contrib += Ai[index_oneD_aux]*exp(-Bi[index_oneD_aux]*dist_atom_grid * dist_atom_grid );
                                }



                                // === Increment the density at the position (nx,ny,nz)
                                int index_aux                 = index_array(nx, ny, nz, map_av->map_dim[XX], map_av->map_dim[YY], map_av->map_dim[ZZ]);
                                map_inst->vox[index_aux] += density_contrib;

                            } // Finish if statement that only distances smaller than the cut off are considered for the density
                        }     // finish if statement that grid points are within the lattice
                    }         // finish over z grid points
                }             // finish over y grid points
            }                 // finish over x grid points
        }                     // finish if statement considering atoms lying in the lattice
    }                         // finish loop over surrounding atoms



    // ===  update average and stdev maps ===
    for (int l = 0; l < map_av->map_dim[XX]*map_av->map_dim[YY]*map_av->map_dim[ZZ]; l++)
    {
        map_av->vox[l] += map_inst->vox[l];
        map_sd->vox[l] += map_inst->vox[l]*map_inst->vox[l];
    }
}


//! \brief Accumulate block average density
void accumulate_block_average( gmx_unused t_mapdata *map_block, gmx_unused t_mapdata *map_inst  )
{
    int l;
    for (l = 0; l < map_block->map_dim[XX]*map_block->map_dim[YY]*map_block->map_dim[ZZ]; l++)
    {
        map_block->vox[l] += map_inst->vox[l];
    }
}

//! \brief Set density of a map to zero
void reset_density( gmx_unused t_mapdata *map)
{
    int l;
    for (l = 0; l < map->map_dim[XX]*map->map_dim[YY]*map->map_dim[ZZ]; l++)
    {
        map->vox[l] = 0.0;
    }
}


//! \brief Compute timed-average average density
void av_map ( gmx_unused t_mapdata *map_av,  float Nframes  )
{

    int l;
    // === Calculate the density, normalizing the lattice values ===
    // density = Ncounts / (Nframes )
    for (l = 0; l < map_av->map_dim[XX]*map_av->map_dim[YY]*map_av->map_dim[ZZ]; l++)
    {
        // average
        map_av->vox[l] /= (Nframes);      // (Angstrom^-1)
    }

}

//! \brief Compute timed-average standard deviation density
void sd_map ( gmx_unused t_mapdata *map_av,  gmx_unused t_mapdata *map_sd, float Nframes  )
{

    int l;

    if (Nframes > 1)
    {
        for (l = 0; l < map_av->map_dim[XX]*map_av->map_dim[YY]*map_av->map_dim[ZZ]; l++)
        {
            //calculate standard deviation
            map_sd->vox[l] = sqrt(  map_sd->vox[l]/(Nframes-1.0) - map_av->vox[l]*map_av->vox[l]*Nframes/(Nframes-1.0) );

        }
    }
}

//! \brief Print the array in dat format: x, y, z, rho, sigma
void print_rho_sigma( const char *N_file,   rvec res, gmx_unused t_mapdata *map_av,  gmx_unused t_mapdata *map_sd  )
{
    // Nbins of the lattice: (bin[XX],bin[YY],bin[ZZ])
    // Resolution (size of each grid) (resx,resy,resz) (nm)
    // density array = density []
    // standard deviation array = density_sigma []
    FILE *pfile;
    pfile = fopen(N_file, "w");

    // === line 1 : Header lines ===
    fprintf(pfile, "# Rho ( and its standard dev) vs x y z.\n");
    fprintf(pfile, "# grid_res= %f %f %f (nm), dimensions: bin[XX]: %d bin[YY]: %d bin[ZZ]: %d\n", res[XX], res[YY], res[ZZ], map_av->map_dim[XX], map_av->map_dim[YY], map_av->map_dim[ZZ]);
    fprintf(pfile, "# x y z (Angstrom) rho sigma (Angstrom^-3)\n");


    // === Density array ===
    int   i, j, k;
    float x, y, z;


    for (i = 0; i < map_av->map_dim[XX]; i++)         // loop over x
    {
        for (j = 0; j < map_av->map_dim[YY]; j++)     // loop over y
        {
            for (k = 0; k < map_av->map_dim[ZZ]; k++) // loop over z

            {                                         //=== set the absolute coordinate system (the simulation box)by adding (minx,miny,minz) to each lattice value
                // x,y,z in Angstrom
                x = NM2AA*(i+map_av->origin[XX])*res[XX];
                y = NM2AA*(j+map_av->origin[YY])*res[YY];
                z = NM2AA*(k+map_av->origin[ZZ])*res[ZZ];

                int index_aux = index_array(i, j, k, map_av->map_dim[XX], map_av->map_dim[YY], map_av->map_dim[ZZ]);
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
        ivec N;    // (NA,NB,NC)
        ivec Nmin; // (AMIN,BMIN,CMIN)
        ivec Nmax; // (AMAX,BMAX,CMAX)
        rvec BOX;  //(BOXA,BOXB,BOXC)

        // subcell where the density is actually computed
        ivec    bin;         // Grid points (0 0 0) (working lattice)

        rvec    res;         // grid resolution  (box/N)
        rvec    L0;          // origin of the lattice with respect to the origin of the map Nmin*res
        float   cutoff;      // (nm) density will be incremented in the grid points within a cutoff distance of a certain atom
        float   mask_cutoff; //|rho|<mask_rheshold will be ignored from the calculation
        float   cc;          // global correlation coefficient
        float   time_window; // for time-resolved global correlation
        int     Nslide;      // Size of the sliding grid for local-correlation calculation
        // ==== selection variables ===
        Selection                        refsel_;
        SelectionList                    sel_;
        int     nr;       // number of atoms in the selection
        int     Ntotal;   // total number of atoms
        int     Nattypes; // Nummber of atomtypes (default = 6)

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;

        // === Definition of the parameters for the input/output files ===

        // strings to read the input and output file names
        std::string                      gauss_coef_; // input structure factor file
        std::string                      map_av_;     // average map (3d)
        std::string                      map_sd_;     // standard deviation (3d)
        std::string mi_;
        std::string diff_;
        std::string mask_;
        std::string mlog_;
        std::string loc_corr_;
        std::string                       rho_; // density (3d), same as the map but in a differnt format: x, y, z, rho, sigma
        std::string rho_diff_, rho_mlog_;
        std::string                       corr_;

        FILE *pfile_corr = nullptr; // to print time-resolved correlation


        // === structure factor variables ===
        char               gauss_coef_atname[10] = {'\0'}; // atom names for crosscheck
        char               gauss_coef_attype[10] = {'\0'}; // atom types
        int                Ngaussians;                     // Number of Gaussians that add up to the density
        bool               Norm;                           // Normalization
        bool               incl_zero;                      // include zero-density values in the spatial average and stdev calculation (over the grid points)
        bool               shift_mlog;                     // remove mean value from mlog map

        std::vector<float> A;                              // Ai factors (Ngaussians x Natoms 1D array)
        std::vector<float> B;                              // Bi factors (Ngaussians x Natoms 1D array)
        std::vector<char>  attypes;                        // atomtypes for A and B array
        std::vector<int>   attype_position;                // array will store the index of the atom type for the A and B Gaussian-coefficient arrays

        // === define the lattice that will store the density for each point===
        t_mapdata          map_ref;      /* The reference=experimental map to compare to        */
        t_mapdata          map_av;       /* The computed timed-averaged map         */
        t_mapdata          map_inst;     /* The computed instantaneous         */
        t_mapdata          map_sd;       /* The computed timed-averaged standard deviation map */
        t_mapdata          map_diff;     /* diff map = computed_map - ref_map */
        t_mapdata          map_mask;     /* mask map*/
        t_mapdata          map_mlog;     /* minus log map*/
        t_mapdata          map_loc_corr; /* local correlation map*/
        t_mapdata          map_block;    /* map for block averaging*/

        ivec               ncutoff;      // === Definition of the cutoff and other gaussian-density related variables
        // Density will be increment in the grid points between:   (i-ncutoff, j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)

        // === Variables to be used in the trajectory loop ===
        float Nframes;           // Nframes
        float Nframes_block;     // Nframes in the time_window for block averaging
};



TimeAveragedDensity::TimeAveragedDensity()
{
    registerAnalysisDataset(&data_, "density");

// Defaults
    N[XX]    = 250; N[YY] = 250; N[ZZ] = 250;
    Nmin[XX] = 0; Nmin[YY] = 0; Nmin[ZZ] = 0;
    Nmax[XX] = 249; Nmax[YY] = 249; Nmax[ZZ] = 249;
    BOX[XX]  = NM2AA; BOX[YY] = NM2AA; BOX[ZZ] = NM2AA;

    bin[XX]        = 0; bin[YY] = 0; bin[ZZ] = 0;
    res[XX]        = 0; res[YY] = 0; res[ZZ] = 0;
    L0[XX]         = 0; L0[YY] = 0; L0[ZZ] = 0;
    cutoff         = 0.5;
    mask_cutoff    = 1e-5;
    cc             = 0;
    time_window    = 1.0;
    Nslide         = 5;
    nr             = 0;
    Ntotal         = 0;
    Nattypes       = 6;
    Ngaussians     = 2;
    Norm           = true;
    incl_zero      = false;
    shift_mlog     = true;
    ncutoff[XX]    = 0; ncutoff[YY] = 0; ncutoff[ZZ] = 0;
    Nframes        = -1;
    Nframes_block  = -1;
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
        "How many atom types can be changed with -Nattypes\n",
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
        "    -box:  --(BOXA, BOXB,BOXC) = crystal cell dimensions (nm,nm,nm)\n",
        "    -min: --(AMIN,BMIN,CMIN) = starting grid point in the considered portion of the map\n",
        "    -max: --(AMAX,BMAX,CMAX) = stopping grid point in the considered portion of the map\n",
        "    Note: The grid resolution is computed as BOXi/Ni, independently for each axis i=a,b,c.\n",
        "    Note: The density will be ONLY computed on a sub-lattice going from (AMIN,BMIN,CMIN) to (AMAX,BMAX,CMAX).\n",
        "    Note: The density is absolute with respect to the grid defined by the integer vector Ni. There is no fit procedure.\n",
        "    Nore: If the density is intended to be calculated with respect to other reference group, a previous fitting step (e.g. with trjconv) should precede the gtide calculation.\n",

        "===Examples===:\n",
        "(1) Time-average density of the Protein in a box of size (10 nm,10 nm,10 nm), with its origin at (x,y,z)=(0,0,0) and with a grid resolution of 0.04 nm:\n",
        "gmx tide -f traj.xtc -s conf.gro -box 10 10 10 -size 250 250 250 -min 0 0 0 -max 249 249 249 -select \'group \"Protein\"\' -map_av\n",
        "(2) Same as example 1, but reading non-default Gaussian parameters, for 9 atom types, and 4 Gaussians contributing to the density, and expanding the Gaussians up to 0.8 nm:\n",
        "gmx tide -f traj.xtc -s conf.gro -Gauss_coef Gauss_coef.dat -box 10 10 10 -size 250 250 250 -min 0 0 0 -max 249 249 249 -Ngaussians 4 -Nattypes 9 -cutoff 0.8 -select \'group \"Protein\"\' -map_av\n",
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

    options->addOption(IntegerOption("size").store(N).vector()
                           .description("(NA,NB,NC) = total number of grid points along each cell edge size"));

    options->addOption(IntegerOption("min").store(Nmin).vector()
                           .description("(Amin,Bmin,Cmin) = starting grid points in the considered portion of the map"));

    options->addOption(IntegerOption("max").store(Nmax).vector()
                           .description("(Amax,Bmax,Cmax) = stopping grid points in the considered portion of the map\nDensity will be calculated within min-max grid points"));

    options->addOption(RealOption("box").store(BOX).vector()
                           .description("crystal cell dimensions. Grid spacing computed as BOX/N for each axis"));

    options->addOption(FloatOption("cutoff").store(&cutoff)
                           .description("Cutoff distance to compute the density around a certain atom (nm)"));

    options->addOption(IntegerOption("Ngaussians").store(&Ngaussians)
                           .description("Number of Gaussians to use in the sum (e.g. 2,4, or 5)"));

    options->addOption(IntegerOption("Nattypes").store(&Nattypes)
                           .description("Number of different atom types"));

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


    // === total Number of atoms===
    Ntotal    = top.topology()->atoms.nr;
    fprintf(stdout, "Total number of atoms: %d\n", Ntotal);

    // === GAUSSIAN COEFFICIENTS BASED ON ATOM TYPES===
    fprintf(stdout, "Considering %d Gaussians and %d different atom types\n", Ngaussians, Nattypes);

    set_gaussian_coeff(gauss_coef_, &A, &B, &attypes, Ngaussians, Nattypes);
    attype_position = attype_positions_from_top(attypes, Ntotal, Nattypes, top );


//=== INPUT REFERENCE MAP ====


    /* === Read input file (map) === */
    if (mi_.empty() == 0)
    {
        printf("=== Reading reference map\n");

        read_xplor_3dmap(mi_.c_str(), &map_ref );
    }


    /* === Mask map === */
    if (mask_.empty() == 0)  // read map if provided
    {
        printf("=== Reading mask map\n");
        read_xplor_3dmap(mask_.c_str(), &map_mask );
    }


    // === CREATE THE LATTICE ===
    if (mi_.empty() == 0)
    {
        printf("Setting map parameters equal to those of the reference map. -size, -min, -max, and -box options ignored\n");
        for (int i = 0; i < DIM; i++)
        {
            BOX[i] = map_ref.cell[i]; // cell size

            N[i] = map_ref.grid[i];   // cell grid points

            //number of bins (working lattice)
            bin[i] = map_ref.map_dim[i];

            //resolution
            res[i] = map_ref.cell[i]/map_ref.grid[i];

            // origin of the lattice
            L0[i] = map_ref.origin[i]*res[i];

            // origin in grid points
            Nmin[i] = map_ref.origin[i];
        }


    }
    else            // no reference map provided, read variables the -size, -min, -min , -box options

    {
        for (int i = 0; i < DIM; i++)
        {
            //number of bins (working lattice)
            bin[i] = Nmax[i]-Nmin[i]+1;         //Amax-Amin+1

            //resolution
            res[i] = BOX[i]/N[i];

            // origin of the lattice
            L0[i] = Nmin[i]*res[i];
        }
    }


    // ====  initialize the maps =====
    char  title[] = {'m', 'a', 'p'};
    float cell[6];
    for (int i = 0; i < DIM; i++)
    {
        cell[i] = BOX[i];                            // cell 0...2 are the box sizes
    }
    for (int i = 3; i < 2*DIM; i++)
    {
        cell[i] = 90.0;                              // cell 3...5 are the angles, always 90°
    }
    //N  Number of grid points within BOX
    //  MIN[i]; // Origin of the working lattice
    int axes_order[3] = {0, 1, 2};                  // x, y, z (guess) but  not used in tide anyways
    // bin = map_dim is the working lattice
    // spacegroup=1; // this will be assigned arbitrarily. It won't be used in gtide anyways

    // minimum and maximum density values
    double min  = 0, max = 0;
    float  mean = 0, rms = 0;

    //skew options will be ignored (0)
    float skew_mat[9];
    float skew_trans[3];
    for (int i = 0; i < 9; i++)
    {
        skew_mat[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        skew_trans[i] = 0;
    }


    printf("=== Initial MD map parameters:\n");
    init_map(&map_av, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 1);
    init_map(&map_inst, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0);
    init_map(&map_sd, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0);
    init_map(&map_diff, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0);
    init_map(&map_mlog, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0);
    init_map(&map_block, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0);
    init_map(&map_loc_corr, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0);


    // if mask was not provided:
    if (mask_.empty() == 1)
    {
        init_map(&map_mask, title, 2, cell, N, Nmin,  axes_order, bin, 1, min, max, mean, rms, skew_mat, skew_trans, 0 );
    }


    //allocate memory for the lattices
    allocate_density_grid(map_av.map_dim, &(map_av.vox));
    allocate_density_grid(map_inst.map_dim, &(map_inst.vox));
    allocate_density_grid(map_sd.map_dim, &(map_sd.vox));
    allocate_density_grid(map_diff.map_dim, &(map_diff.vox));
    allocate_density_grid(map_mlog.map_dim, &(map_mlog.vox));
    allocate_density_grid(map_block.map_dim, &(map_block.vox));
    allocate_density_grid(map_loc_corr.map_dim, &(map_loc_corr.vox));


    // if mask was not provided:
    if (mask_.empty() == 1)
    {
        allocate_density_grid(map_mask.map_dim, &(map_mask.vox));
    }

    // set mask
    //           set_mask(map3d_mask, map_mask, mask_cutoff);
    set_mask(mask_, &map_mask, mask_cutoff);



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
    nr    = sel.posCount();               // NUMBER OF ATOMS IN THE SELECTION (nr
//    fprintf(stdout, "Natoms in selection: %d\n", nr);


    // === ACCUMMULATE DENSITY===
    accumulate_denstiy(&map_av, &map_inst, &map_sd, nr, sel,  res, L0, ncutoff, cutoff, Ngaussians, A, B, attype_position);


    // === ACCUMULATE BLOCK AVERAGE ===
    accumulate_block_average(&map_block, &map_inst);


    // === Increment the number of frames variables ===
    Nframes++;
    Nframes_block++;


    // === calculate and print global correlation if requested ===
    if (!corr_.empty())
    {
        if (fr.time/time_window-floor(fr.time/time_window) == 0 && fr.time > 0)
        {

            av_map(&map_block, Nframes_block );
            cc = global_correlation_coeff(&map_ref, &map_block, &map_mask);
            fprintf(pfile_corr, "%f %f\n", fr.time, cc);
            reset_density(&map_block);
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
    av_map ( &map_av,    Nframes);
    sd_map ( &map_av,  &map_sd,  Nframes);


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
        cc = global_correlation_coeff(&map_ref, &map_av, &map_mask);
        printf("The correlation coefficient between computed and reference maps is %f\n", cc);
    }


    // calculate local correlation if requested
    if (!loc_corr_.empty())
    {
        local_correlation_coeff(&map_ref, &map_av, &map_mask, &map_loc_corr, Nslide);
        map_loc_corr.mean = cc;    // mean will be the global correlation coefficient
        map_loc_corr.rms  = nanf("");
    }


    // shift the map by an amount equal to the average
    if (Norm)
    {
        shift_map(&map_av, -map_av.mean);
        shift_map(&map_sd, -map_sd.mean);
    }

// Normalize the map by an amount equal to the standard deviation
    if (Norm)
    {
        normalize_map(&map_av, map_av.rms);
        normalize_map(&map_sd, map_sd.rms);
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
