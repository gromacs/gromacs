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

#include "gromacs/applied-forces/maputil.h"
#include "gromacs/math/container/mask.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/fileio/xplor.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/trajectoryanalysis.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

//! \brief Convert between XPLOR and GROMACS length units
const float NM2AA = 10.0;

struct GaussianParameters{
    float amplitude;
    float precision;
};

//! \brief Gets the index in an uni-dimensional array from the three coordinates
//! of a three-dimensional array
int index_array(int i, int j, int k, int Nx, int Ny)
{
    // The unidimensional array contains Nx*Ny*Nz elements:  0 ... Nx*Ny*Nz-1
    return k * (Nx * Ny) + j * Nx + i;
}

/*! \brief Compare variables relevant for tide (cell,grid, origin and map_dim)
 * of map1 and map2.
 *
 *  \return true if their dimensions are equal and false otherwise===
 */
bool mapsHaveSameGrid(const t_mapdata &map1, const t_mapdata &map2)
{
    const float cell_tolerance =
        0.0001; // diff smaller than that will be treated as zero.

    for (int i = 0; i < 2 * DIM; i++)
    {
        if ((map1.cell[i] - map2.cell[i]) * (map1.cell[i] - map2.cell[i]) >
            cell_tolerance * cell_tolerance)
        {
            return false;
        }
    }
    return (map1.grid == map2.grid) && (map1.origin == map2.origin) &&
           (map1.map_dim == map2.map_dim);
}

//! \brief Subtract maps: rho(map3) = rho(map2)-rho(map1)
std::vector<float> subtract_maps(const std::vector<float> &subtrahend, const std::vector<float> &minuend)
{
    std::vector<float> result;
    std::transform(std::begin(minuend), std::end(minuend), std::begin(subtrahend), std::back_inserter(result), [](float min, float sub){return min-sub; });
    return result;
}

//! \brief log of map: result = - ln ( rho(map) )  PSEUDO FREE ENERGY
std::vector<float> negLogVector(const std::vector<real> &input)
{
    std::vector<float> result;
    std::transform(std::begin(input), std::end(input), std::back_inserter(result), [](float value){return value > 0 ? -log(value) : nanf(""); });
    return result;
}

//! \brief log of division of maps: rho(map3) = - ln ( nominator / denominator )
std::vector<float> negLogRatioVector(const std::vector<float> &denominator,
                                     const std::vector<float> &nominator)
{
    std::vector<float> result;

    const auto        &negLogRatio =
        [] (float nom, float den) -> float {
            if ((den != 0) && (nom/den > 0))
            {
                return -log(nom/den);
            }
            return nanf("");
        };

    std::transform(std::begin(nominator), std::end(nominator), std::begin(denominator), std::back_inserter(result), negLogRatio);
    return result;
}

//! \brief Compute spatial average over the map: <rho(map)> over all cell
//! positions
void setMapMeanAndRms(t_mapdata *map, bool incl_zero, const Mask &mask)
{
    map->mean = 0.0;
    map->rms  = 0.0;
    float nsamples = 0;
    for (size_t l = 0; l < map->vox.size(); l++)
    {
        if (mask[l] && (map->vox[l] != 0 || incl_zero))
        {
            map->mean += map->vox[l];
            map->rms  += map->vox[l] * map->vox[l];
            nsamples++;
        }
    }
    if (nsamples > 0)
    {
        map->mean /= nsamples;
        map->rms   = sqrt(map->rms / (nsamples - 1) -
                          map->mean * map->mean * nsamples / (nsamples - 1));
    }
}

//! \brief TODO
float calc_sum_rhoA_rhoB_mask(const t_mapdata &map1,
                              const t_mapdata &map2,
                              const Mask      &mask)
{
    if (!mapsHaveSameGrid(map1, map2))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    double sum_ijk = 0.0;
    for (size_t l = 0; l < map1.vox.size(); l++)
    {
        if (mask[l])
        {
            sum_ijk += map1.vox[l] * map2.vox[l];
        }
    }
    return sum_ijk;
}

//! \brief TODO
float global_correlation_coeff(const t_mapdata &map1,
                               const t_mapdata &map2,
                               const Mask      &mask)
{
    if (!mapsHaveSameGrid(map1, map2))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }
    float numerator   = calc_sum_rhoA_rhoB_mask(map1, map2, mask);
    float denominator = sqrt(calc_sum_rhoA_rhoB_mask(map2, map2, mask) *
                             calc_sum_rhoA_rhoB_mask(map1, map1, mask));
    return numerator / denominator;
}

//! \brief Calculate the local correlation
void local_correlation_coeff(const t_mapdata &map1, const t_mapdata &map2, const Mask &mask,
                             t_mapdata *loc_corr, int Nodd)
{
    int Nhalf = int(Nodd / 2); // center point of the sliding grid (e.g. if Nodd=5
                               // ; then Nhalf=2, Note that it starts counting
                               // from 0)

    // This can only work for maps with equal dimensions
    if (!mapsHaveSameGrid(map1, map2) || !mapsHaveSameGrid(map1, *loc_corr))
    {
        gmx_fatal(FARGS, "Map dimensions must be identical");
    }

    // Initialize the output loc_corr map, setting the borders and masked points,
    // with non-assigned NAN values
    std::fill(std::begin(loc_corr->vox), std::end(loc_corr->vox), nanf(""));

    // loop over the map1 (reduced by an amount Nhalf at both ends)
    for (int i = Nhalf; i < map1.map_dim[XX] - Nhalf; i++)
    {
        for (int j = Nhalf; j < map1.map_dim[YY] - Nhalf; j++)
        {
            for (int k = Nhalf; k < map1.map_dim[ZZ] - Nhalf; k++)
            {
                // reset correlation variables
                float sum12 = 0.0;
                float sum11 = 0.0;
                float sum22 = 0.0;
                // slide over Nodd points, centered at Nhalf to accumulate the density
                for (int i_sli = i - Nhalf; i_sli <= i + Nhalf; i_sli++)
                {
                    for (int j_sli = j - Nhalf; j_sli <= j + Nhalf; j_sli++)
                    {
                        for (int k_sli = k - Nhalf; k_sli <= k + Nhalf; k_sli++)
                        {
                            int id_aux = index_array(i_sli, j_sli, k_sli, map1.map_dim[XX],
                                                     map1.map_dim[YY]);
                            if (mask[id_aux]) // skipped masked points
                            {
                                sum12 += map1.vox[id_aux] * map2.vox[id_aux];
                                sum11 += map1.vox[id_aux] * map1.vox[id_aux];
                                sum22 += map2.vox[id_aux] * map2.vox[id_aux];
                            }
                        }
                    }
                }

                // assign correlation to middle point of sliding grid (i,j,k)
                int id_aux = index_array(i, j, k, map1.map_dim[XX], map1.map_dim[YY]);
                loc_corr->vox[id_aux] = sum12 / sqrt(sum22 * sum11);
            }
        }
    }
}

//! \brief Set the Gaussian coefficients A and B for the six atom types H, C, N,
//! O, P, and S
std::map < std::string, std::vector < GaussianParameters>> set_gaussian_coeff( std::string parameterFilenName )
{
    // read gaussian parameters from file if provided, else use defaults
    FILE *parameterFile;
    if (!parameterFilenName.empty())
    {
        fprintf(stdout, "Reading file with Gaussian coefficients\n");
        parameterFile = fopen(parameterFilenName.c_str(), "r");
    }
    else
    {
        fprintf(stdout, "Reading default Gaussian coefficients.\n");
        parameterFile = libopen("electronscattering.dat");
    }
    std::map < std::string, std::vector < GaussianParameters>> gaussianParametersLookup;
    char line[1000];
    // all non-header lines
    while (get_a_line(parameterFile, line, 1000))
    {
        char  currentAtomType[100];
        float amplitude;
        float precision;
        if (sscanf(line, " %s %e %e", currentAtomType, &amplitude, &precision) == 3)
        {
            gaussianParametersLookup[currentAtomType].push_back({amplitude, precision});
        }
    }
    fclose(parameterFile);
    return gaussianParametersLookup;
}

//! \brief  Output integer position in the array of atomtypes based on first
//! character of atomname
std::vector<std::string> attypeStrings_from_top(const t_topology &top)
{
    std::vector<std::string> attypeString;
    for (int n = 0; n < top.atoms.nr; n++)
    {
        attypeString.push_back(std::string(*(top.atoms.atomname[n])));
    }
    return attypeString;
}

/*! \brief Outputs the density map to XPLOR format.
 */
void print_to3dmap(const std::string &inputFileName, const t_mapdata &map,
                   const Mask &mask)
{
    if (inputFileName.empty())
    {
        fprintf(
                stdout,
                "No filename given for X-PLOR output file, not outputting densities\n");
        return;
    }

    XplorData xplorOutput(map);
    mask.maskToNan(&(xplorOutput.data));
    FILE     *fp = gmx_ffopen(inputFileName.c_str(), "w");
    xplorOutput.write(fp);
    gmx_ffclose(fp);
}

void printMapInfo(const t_mapdata &map)
{
    printf("\n");
    printf("Header:\n");
    printf("Title      : %s\n", map.title);
    printf("Datamode   : %u\n", map.datamode);
    printf("Cell       : %g A  %g A  %g A  %g %g %g\n", map.cell[0], map.cell[1],
           map.cell[2], map.cell[3], map.cell[4], map.cell[5]);
    printf("Grid       : %d x %d x %d\n", map.grid[XX], map.grid[YY],
           map.grid[ZZ]);
    printf("Origin     : %d A  %d A  %d A\n", map.origin[XX], map.origin[YY],
           map.origin[ZZ]);
    printf("Axes order : %d %d %d\n", map.axes_order[0], map.axes_order[1],
           map.axes_order[2]);
    printf("Map_dim    : %d %d %d\n", map.map_dim[0], map.map_dim[1],
           map.map_dim[2]);
    printf("Spacegroup : %d\n", map.spacegroup);
    printf("Min        : %g\n", map.min);
    printf("Max        : %g\n", map.max);
    printf("Mean       : %g\n", map.mean);
    printf("RMS        : %g\n", map.rms);
    printf("Skew_trans : %g %g %g\n", map.skew_trans[0], map.skew_trans[1],
           map.skew_trans[2]);
    printf("Skew_mat   : ");
    for (int i = 0; i < 9; i++)
    {
        printf("%g ", map.skew_mat[i]);
    }
    printf("\n\n");
}

//! \brief Init map function
void init_map(gmx_unused t_mapdata *map, const std::string &title,
              const std::array<float, 2 * DIM> &cell,
              const std::array<int, DIM> &grid,
              const std::array<int, DIM> &origin,
              const std::array<int, DIM> &map_dim, double min,
              double max, float mean, float rms)
{
    map->title      = gmx_strdup(title.c_str());
    map->datamode   = 2;           // stored as Reals
    map->cell       = cell;        // cell 0...2 are the box sizes // cell 3...5 are the angles
    map->grid       = grid;        // Number of grid points within box
    map->origin     = origin;      // Origin of the working lattice
    map->axes_order = {{0, 1, 2}}; // x, y, z (guess)
    map->map_dim    = map_dim;     // the working lattice
    map->spacegroup = 1;           // this will be assigned arbitrarily. It won't
    // be used in gtide anyways
    // minimum and maximum density values
    map->min  = min;
    map->max  = max;
    map->mean = mean;
    map->rms  = rms;
    // skew options will be ignored (0)
    map->skew_mat = {};
    map->skew_trans = {};
}

/*!  \brief Read an input density map using XPLOR format
 *
 * For the format see:
 * http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
void read_xplor_3dmap(const char           *inputFileName,
                      gmx_unused t_mapdata *ref_map)
{
    if (inputFileName == nullptr)
    {
        fprintf(stdout, "No filename given for X-PLOR input file\n");
        return;
    }

    FILE     *fp = gmx_ffopen(inputFileName, "r");
    XplorData xplorData(fp);
    fclose(fp);

    //=== assign the information read in the xplor file to the map ===

    std::array<float, 2 * DIM> cell;
    for (int i = 0; i < DIM; i++)
    {
        cell[i] = xplorData.cell_length[i]; // cell 0...2 are the box sizes
    }
    for (int i = 3; i < 2 * DIM; i++)
    {
        cell[i] = xplorData.cell_angles[i - 3]; // cell 3...5 are the angles
    }
    std::array<int, DIM> bin = {{
                                    xplorData.crs_end[XX] - xplorData.crs_start[XX] + 1,
                                    xplorData.crs_end[YY] - xplorData.crs_start[YY] + 1,
                                    xplorData.crs_end[ZZ] - xplorData.crs_start[ZZ] + 1
                                }};

    auto                 minmax =
        std::minmax_element(std::begin(xplorData.data), std::end(xplorData.data));
    init_map(ref_map, "ref", cell, xplorData.extend, xplorData.crs_start,
             bin, *(minmax.first), *(minmax.second),
             xplorData.mean_value, xplorData.rms_value);
    // copy lattice into vox
    ref_map->vox = xplorData.data;
}

//! \brief Accummulate density every time step
void accumulate_denstiy(gmx_unused t_mapdata *map_av,
                        gmx_unused t_mapdata *map_inst,
                        gmx_unused t_mapdata *map_sd, Selection sel, RVec res,
                        RVec latticeOrigin_, IVec ncutoff, float cutoff,
                        const std::map < std::string, std::vector < GaussianParameters>> &gaussianParametersLookup,
                        const std::vector<std::string> &attypeString)
{

    //=== INITIALIZE THE INSTANTANEOUS MAP LATTICE===
    containeroperation::setZero(&(map_inst->vox));
    for (int atomIndex = 0; atomIndex < sel.atomCount(); atomIndex++)
    {
        // Select the first character of the atom name
        std::string atomName =
            attypeString[sel.atomIndices()[atomIndex]].substr(0, 1);
        const auto &gaussianParameter = gaussianParametersLookup.at(atomName);

        // Change of system of coordinates to one where the origin is the origin of
        // the lattice
        // latticeOrigin_ origin of the lattice coordinates
        // R_trans = coordinates of the atom with respect to the lattice
        // R_trans= R - latticeOrigin_
        rvec R_trans;
        rvec_sub(sel.position(atomIndex).x(), latticeOrigin_, R_trans);

        // === Discretize the Rtrans to locate it in the lattice
        const int i = floor(R_trans[XX] / res[XX]);
        const int j = floor(R_trans[YY] / res[YY]);
        const int k = floor(R_trans[ZZ] / res[ZZ]);

        // === only consider the atoms located in the lattice ===
        if (i >= 0 && i < map_av->map_dim[XX] && j >= 0 &&
            j < map_av->map_dim[YY] && k >= 0 && k < map_av->map_dim[ZZ])
        {
            for (auto g : gaussianParameter)
            {

                // Increment the density of the grid points to the atom, laying within a
                // certain cutoff
                // loops over grid points around the (i,j,k) grid point
                for (int nx = i - ncutoff[XX]; nx <= i + ncutoff[XX]; nx++)
                {
                    for (int ny = j - ncutoff[YY]; ny <= j + ncutoff[YY]; ny++)
                    {
                        for (int nz = k - ncutoff[ZZ]; nz <= k + ncutoff[ZZ]; nz++)
                        {
                            // === only consider grid points lying in the lattice ===
                            if (nx >= 0 && nx < map_av->map_dim[XX] && ny >= 0 &&
                                ny < map_av->map_dim[YY] && nz >= 0 &&
                                nz < map_av->map_dim[ZZ])
                            {
                                rvec grid_vec;
                                // vector coordinates of the (nx,ny,nz) grid point
                                grid_vec[XX] = nx * res[XX];
                                grid_vec[YY] = ny * res[YY];
                                grid_vec[ZZ] = nz * res[ZZ];

                                // Distance between the atom and the grid_vector:
                                rvec  R_atom_grid;
                                rvec_sub(R_trans, grid_vec,
                                         R_atom_grid); // R_atom_grid  = R_trans - grid_vec

                                // Calculate the norm of R_atom_grid	(Angstrom) to coincide
                                // with the ai, bi and c factors
                                float dist_atom_grid = norm(R_atom_grid);

                                // === Only contribute to the density if the distance is smaller
                                // than the cutoff
                                if (dist_atom_grid < cutoff)
                                {
                                    // === Increment the density at the position (nx,ny,nz)
                                    // === Gaussian contribution (two or four gaussians depending on
                                    // the atom) to the density at this grid point
                                    int index_aux = index_array(nx, ny, nz, map_av->map_dim[XX],
                                                                map_av->map_dim[YY]);
                                    map_inst->vox[index_aux] +=  g.amplitude * exp(-g.precision * dist_atom_grid * dist_atom_grid);

                                } // only distances smaller than the cut off are considered for the density
                            }     // grid points are within the lattice
                        }         // z grid points
                    }             // y grid points
                }                 // x grid points
            }
        }                         // if atom in lattice
    }                             // atoms

    // ===  update average and stdev maps ===
    containeroperation::add(&(map_av->vox), map_inst->vox);

    std::transform(std::begin(map_sd->vox), std::end(map_sd->vox),
                   std::begin(map_inst->vox), std::begin(map_sd->vox),
                   [](float oldSD, float voxelValue) {
                       return oldSD + voxelValue * voxelValue;
                   });
}

//! \brief Compute timed-average standard deviation density
void sd_map(const t_mapdata &map_av, gmx_unused t_mapdata *map_sd,
            float Nframes)
{
    if (Nframes > 1)
    {
        std::transform(std::begin(map_av.vox), std::end(map_av.vox),
                       std::begin(map_sd->vox), std::begin(map_sd->vox),
                       [Nframes](float av, float sd) {
                           return sqrt(sd / (Nframes - 1.0) -
                                       av * av * Nframes / (Nframes - 1.0));
                       });
    }
}

//! \brief Print the array in dat format: x, y, z, rho, sigma
void print_rho_sigma(const std::string &inputFileName, rvec res,
                     const t_mapdata &map_av,
                     const t_mapdata &map_sd)
{
    FILE *pfile;
    pfile = fopen(inputFileName.c_str(), "w");

    // === line 1 : Header lines ===
    fprintf(pfile, "# Rho ( and its standard dev) vs x y z.\n");
    fprintf(pfile, "# grid_res= %f %f %f (nm), dimensions: bin[XX]: %d bin[YY]: "
            "%d bin[ZZ]: %d\n",
            res[XX], res[YY], res[ZZ], map_av.map_dim[XX], map_av.map_dim[YY],
            map_av.map_dim[ZZ]);
    fprintf(pfile, "# x y z (Angstrom) rho sigma (Angstrom^-3)\n");

    // === Density array ===
    for (int i = 0; i < map_av.map_dim[XX]; i++)         // loop over x
    {
        for (int j = 0; j < map_av.map_dim[YY]; j++)     // loop over y
        {
            for (int k = 0; k < map_av.map_dim[ZZ]; k++) // loop over z

            {
                //=== set the absolute coordinate system (the simulation box)by adding
                //(minx,miny,minz) to each lattice value
                const float x = NM2AA * (i + map_av.origin[XX]) * res[XX];
                const float y = NM2AA * (j + map_av.origin[YY]) * res[YY];
                const float z = NM2AA * (k + map_av.origin[ZZ]) * res[ZZ];

                int         index_aux =
                    index_array(i, j, k, map_av.map_dim[XX], map_av.map_dim[YY]);
                fprintf(pfile, "%f %f %f %d %d %d %E %E\n", x, y, z, i, j, k,
                        map_av.vox[index_aux], map_sd.vox[index_aux]);
            }
        }
    }
    fclose(pfile);
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

        void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata) override;

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
        std::array<int, DIM>             bin = {
            {0, 0, 0}};                                                 // Grid points (0 0 0) (working lattice)

        RVec                             res            = {0., 0., 0.}; // grid resolution  (box/N)
        RVec                             latticeOrigin_ = {0., 0., 0.}; // origin of the lattice with respect to the origin of
        // the map Nmin*res
        float                            cutoff = 0.5;                  // (nm) density will be incremented in the grid points
                                                                        // within a cutoff distance of a certain atom
        float                            mask_cutoff =
            1e-5;                                                       //|rho|<mask_rheshold will be ignored from the calculation
        float                            time_window = 1.0;             // for time-resolved global correlation
        int                              Nslide      = 5;               // Size of the sliding grid for local-correlation calculation
        // ==== selection variables ===
        Selection                        refsel_;
        SelectionList                    sel_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;

        // === Definition of the parameters for the input/output files ===

        // strings to read the input and output file names
        std::string fnGaussCoefficients_;    // input structure factor file
        std::string fnMapAverage_;           // average map (3d)
        std::string fnMapStandardDeviation_; // standard deviation (3d)
        std::string fnInputMap_;
        std::string fnDifferenceMap_;
        std::string fnMaskMap_;
        std::string mlog_;
        std::string loc_corr_;
        std::string rho_; // density (3d), same as the map but in a differnt format:
                          // x, y, z, rho, sigma
        std::string rho_diff_;
        std::string rho_mlog_;
        std::string corr_;

        FILE       *pfile_corr = nullptr; // to print time-resolved correlation

        // === structure factor variables ===
        bool Norm                  = true;                                                        // Normalization
        bool incl_zero             = false;                                                       // include zero-density values in the spatial average
                                                                                                  // and stdev calculation (over the grid points)
        bool shift_mlog = true;                                                                   // remove mean value from mlog map

        std::map < std::string, std::vector < GaussianParameters>> gaussianParametersPerAtomType; // Ai factors (Ngaussians x Natoms 1D array)
        std::vector<std::string> attypeString;                                                    // array will store the index of the
                                                                                                  // atom type for the A and B
                                                                                                  // Gaussian-coefficient arrays

        // === define the lattice that will store the density for each point===
        t_mapdata map_ref;             /* The reference=experimental map to compare to        */
        t_mapdata map_av;              /* The computed timed-averaged map         */
        t_mapdata map_inst;            /* The computed instantaneous         */
        t_mapdata map_sd;              /* The computed timed-averaged standard deviation map */
        Mask      mask_;               //< Masking out density values
        t_mapdata map_block;           /* map for block averaging*/

        IVec      ncutoff = {0, 0, 0}; // === Definition of the cutoff and other
                                       // gaussian-density related variables
        // Density will be increment in the grid points between:   (i-ncutoff,
        // j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)

        // === Variables to be used in the trajectory loop ===
        float Nframes_block = -1; // Nframes in the time_window for block averaging
};

TimeAveragedDensity::TimeAveragedDensity()
{
    registerAnalysisDataset(&data_, "density");
}

void TimeAveragedDensity::initOptions(IOptionsContainer          *options,
                                      TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes the time-averaged electron density of a group of "
        "atoms within a defined box."
        "The density is calculated as a linear combination of Gaussians:\n",
        "rho(R)=sum_{i=1}^{i=Ngaussians} Ai*exp(-Bi*R^2)\n",

        "=== Gaussian coefficients:\n",
        "By default, gtide considers 2-Gaussian coefficients (A and B), derived "
        "from Electron Microscopy (citation), for H,C,N,O,P,S.\n",
        "Non-default values can be specified via the -Gauss_coef input file.\n",
        "Format of Gauss_coef.dat:\n",
        "atomtype A_1 A_2 ...A_Ngaussian B1 B2 ... B_NGaussian \n",
        "Precomputed A and B parameters, for 2- and 4-Gaussian sums, from X-ray "
        "and EM, as well as for coarse-grained systems, are available at XXXX",
        "How many Gaussians can be changed with -Ngaussians\n",
        "How far each atom spreads in the space can be controlled with "
        "-cutoff.\n",

        "=== atom selection:\n", "The density map is computed for a specified "
        "selection of atoms. Both static and dynamic "
        "selections are supported\n",

        "=== Reference map:\n", "A reference map can be provided with the -mi "
        "option. Different comparison operations between "
        "the computed and the reference map are possible "
        "(see below)\n",

        "=== Masking:\n",
        "A masking map can be provided with the -mask option. Grid points in the "
        "mask, with density smaller than a user-specified"
        "threshold value (-mask_cutoff), will be ignored. In the output maps, "
        "the masked points will be labeled as NAN\n",

        "=== Output files:\n",
        "-map_av: time average of density in XPLOR format.\n",
        "-map_sd: time standard deviation of density in XPLOR format.\n",
        "These two maps can be normalized to sigma units and shifted to zero "
        "mean for easier visualization (-norm)",
        "and the spatial mean can (or not) include zero-value grid points "
        "(-incl_zero).\n",

        "-diff: Difference map:  rho(map_av) - rho(map_reference) in XPLOR "
        "format.\n",

        "-mlog: -log [ rho(map_av) ] or -log [ rho(map_av) /rho(map_reference]  "
        "if a reference map is provided.\n",
        "-shift_mlog will set the mean value of mlog to zero.\n"

        "-corr: Time-resolved global correlation between the computed and the "
        "reference map (xvg format).\n",
        "The correlation is computed in time windows defined by -time_window.\n",
        "-loc_corr: Local correlation between the computed and the reference map "
        "in XPLOR format. \n",
        "For a small cubic grid of dimensions (Nslide,Nslide,Nslide), the "
        "correlation between map_av ",
        "and the reference maps is calculated, and the resulting value is stored "
        "at the center point of this small grid.",
        "The process is repeated by sliding the grid over all the points of the "
        "lattice.\n",
        "The global correlation between the maps is written in the last line of "
        "the xplor map as the mean.",
        "The rms in the last line of the xplor map will be set to NAN.\n",
        "-rho_xxx: the output maps but in ASCII format. \n",

        "=== grid definition:\n"
        "Parameters are similar parameters as in an XPLOR map:\n",
        "If an input reference is provided, its grid parameters are considered. "
        "Otherwise, the following variables must be specified:\n",
        "    -size: Cell grid points (NA,NB,NC)\n",
        "    -box:  --(boxA, boxB,boxC) = crystal cell dimensions (nm,nm,nm)\n",
        "    -min: --(AMIN,BMIN,CMIN) = starting grid point in the considered "
        "portion of the map\n",
        "    -max: --(AMAX,BMAX,CMAX) = stopping grid point in the considered "
        "portion of the map\n",
        "    Note: The grid resolution is computed as boxi/Ni, independently for "
        "each axis i=a,b,c.\n",
        "    Note: The density will be ONLY computed on a sub-lattice going from "
        "(AMIN,BMIN,CMIN) to (AMAX,BMAX,CMAX).\n",
        "    Note: The density is absolute with respect to the grid defined by "
        "the integer vector Ni. There is no fit procedure.\n",
        "    Nore: If the density is intended to be calculated with respect to "
        "other reference group, a previous fitting step (e.g. with trjconv) "
        "should precede the gtide calculation.\n",

        "===Examples===:\n", "(1) Time-average density of the Protein in a box "
        "of size (10 nm,10 nm,10 nm), with its origin at "
        "(x,y,z)=(0,0,0) and with a grid resolution of 0.04 "
        "nm:\n",
        "gmx tide -f traj.xtc -s conf.gro -box 10 10 10 -size 250 250 250 -min 0 "
        "0 0 -max 249 249 249 -select \'group \"Protein\"\' -map_av\n",
        "(2) Same as example 1, but reading non-default Gaussian parameters, for "
        "9 atom types, and 4 Gaussians contributing to the density, and "
        "expanding the Gaussians up to 0.8 nm:\n",
        "gmx tide -f traj.xtc -s conf.gro -Gauss_coef Gauss_coef.dat -box 10 10 "
        "10 -size 250 250 250 -min 0 0 0 -max 249 249 249 -Ngaussians 4 -cutoff "
        "0.8 -select \'group \"Protein\"\' -map_av\n",
        "(3)Time-average density of the Protein, without normalization, "
        "including zeros in the spatial mean, printing as well the difference to "
        "a reference map",
        "(lattice parameters taken from the reference map):\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -nonorm -incl_zero "
        "-select \'group \"Protein\"\' -map_av -diff\n",
        "(4) Same as example 3, but considering a mask, to ignore grid points "
        "which have a density smaller than 1e-5 in the mask:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -mask mask.xplor -nonorm "
        "-incl_zero -mask_cutoff 1e-5 -select \'group \"Protein\"\' -map_av "
        "-diff\n",
        "(5) Same as example 4, but computing -log[rho(map_av)/rho(map_ref)], "
        "shifting the mlog map such that the background spatial mean is zero:",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -mask mask.xplor "
        "-mask_cutoff 1e-5 -shift_mlog -select \'group \"Protein\"\' -mlog \n",
        "(6) Time-resolved global correlation of the computed map to the "
        "reference map, computing the density in time blocks of 1ps:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -time_window 1.0 -select "
        "\'group \"Protein\"\' -corr \n",
        "(7) Local correlation of the computed map to the reference map, sliding "
        "a sub-lattice of size (5,5,5) bins through the main lattice:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -Nslide 5 -select \'group "
        "\"Protein\"\' -loc_corr \n",
        "(8) Print the density, difference and -log to reference map as a list "
        "in ASCII format instead of a XPLOR-format map:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -select \'group "
        "\"Protein\"\' -rho -rho_diff -rho_mlog\n",
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("Gauss_coef")
                           .filetype(eftGenericData)
                           .inputFile()
                           .store(&fnGaussCoefficients_)
                           .defaultBasename("Gauss_coef")
                           .description("File with Gaussian coefficients"));

    options->addOption(FileNameOption("mi")
                           .filetype(eftXPLOR)
                           .inputFile()
                           .store(&fnInputMap_)
                           .defaultBasename("mi")
                           .description("Input reference map (optional)"));

    options->addOption(FileNameOption("mask")
                           .filetype(eftXPLOR)
                           .inputFile()
                           .store(&fnMaskMap_)
                           .defaultBasename("mask")
                           .description("Input mask map (optional)"));

    options->addOption(FileNameOption("map_av")
                           .filetype(eftXPLOR)
                           .outputFile()
                           .store(&fnMapAverage_)
                           .defaultBasename("map_av")
                           .description("Average density map in X-PLOR format"));

    options->addOption(FileNameOption("map_sd")
                           .filetype(eftXPLOR)
                           .outputFile()
                           .store(&fnMapStandardDeviation_)
                           .defaultBasename("map_sd")
                           .description("St. dev. density map in X-PLOR format"));

    options->addOption(
            FileNameOption("diff")
                .filetype(eftXPLOR)
                .outputFile()
                .store(&fnDifferenceMap_)
                .defaultBasename("diff")
                .description("Difference densiy map in X-PLOR format"));

    options->addOption(
            FileNameOption("mlog")
                .filetype(eftXPLOR)
                .outputFile()
                .store(&mlog_)
                .defaultBasename("mlog")
                .description("-ln(map_computed/map_reference) in X-PLOR format"));

    options->addOption(FileNameOption("rho")
                           .filetype(eftGenericData)
                           .outputFile()
                           .store(&rho_)
                           .defaultBasename("rho")
                           .description("av and sd of density in ASCII format"));

    options->addOption(FileNameOption("rho_diff")
                           .filetype(eftGenericData)
                           .outputFile()
                           .store(&rho_diff_)
                           .defaultBasename("rho_diff")
                           .description("difference density in ASCII format"));

    options->addOption(
            FileNameOption("rho_mlog")
                .filetype(eftGenericData)
                .outputFile()
                .store(&rho_mlog_)
                .defaultBasename("rho_mlog")
                .description(
                    "-ln(map_computed/map_reference) density in ASCII format"));

    options->addOption(FileNameOption("corr")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&corr_)
                           .defaultBasename("corr")
                           .description("Time-resolved global correlation "
                                        "between computed and reference maps"));

    options->addOption(FileNameOption("loc_corr")
                           .filetype(eftXPLOR)
                           .outputFile()
                           .store(&loc_corr_)
                           .defaultBasename("loc_corr")
                           .description("Local correlation between computed and "
                                        "reference maps in X-PLOR format"));

    options->addOption(
            SelectionOption("select").storeVector(&sel_).required().description(
                    "Group to calculate the density"));

    options->addOption(IntegerOption("size")
                           .store(nGridPoints.data())
                           .vector()
                           .description("(NA,NB,NC) = total number of grid "
                                        "points along each cell edge size"));

    options->addOption(IntegerOption("min")
                           .store(Nmin.data())
                           .vector()
                           .description("(Amin,Bmin,Cmin) = starting grid points "
                                        "in the considered portion of the map"));

    options->addOption(
            IntegerOption("max")
                .store(Nmax.data())
                .vector()
                .description("(Amax,Bmax,Cmax) = stopping grid points in the "
                             "considered portion of the map\nDensity will be "
                             "calculated within min-max grid points"));

    options->addOption(RealOption("box").store(box).vector().description(
                               "crystal cell dimensions. Grid spacing computed as box/N for each axis"));

    options->addOption(FloatOption("cutoff").store(&cutoff).description(
                               "Cutoff distance to compute the density around a certain atom (nm)"));

    options->addOption(
            FloatOption("mask_cutoff")
                .store(&mask_cutoff)
                .description("Points with density smaller than mask_cutoff will be "
                             "masked according to mask map (-mask)"));

    options->addOption(BooleanOption("norm").store(&Norm).description(
                               "Normalize the output maps. If yes, map in stdev units && shifted such "
                               "that mean=0"));

    options->addOption(
            BooleanOption("incl_zero")
                .store(&incl_zero)
                .description("Include Zero-density points in the calculation of the "
                             "spatial mean and stdev, over the grid"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

    options->addOption(BooleanOption("shift_mlog")
                           .store(&shift_mlog)
                           .description("remove mean value from mlog map"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

    options->addOption(
            FloatOption("time_window")
                .store(&time_window)
                .description("Time window for block-averaging in time-resolved "
                             "global correlation calculation(ps)"));

    options->addOption(IntegerOption("Nslide").store(&Nslide).description(
                               "Size of the sliding grid for local correlation calculation (Must be an "
                               "odd number)"));
}

void TimeAveragedDensity::initAnalysis(
        const gmx_unused TrajectoryAnalysisSettings &settings,
        const TopologyInformation                   &top)
{

    data_.setColumnCount(0, sel_.size());
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    // === GAUSSIAN COEFFICIENTS BASED ON ATOM TYPES===
    gaussianParametersPerAtomType = set_gaussian_coeff(fnGaussCoefficients_);
    attypeString                  = attypeStrings_from_top(*(top.topology()));

    //=== INPUT REFERENCE MAP ====

    // === CREATE THE LATTICE ===
    if (!fnInputMap_.empty())
    {
        /* === Read input file (map) === */
        printf("=== Reading reference map\n");
        read_xplor_3dmap(fnInputMap_.c_str(), &map_ref);
        printf("Setting map parameters equal to those of the reference map. -size, "
               "-min, -max, and -box options ignored\n");
        nGridPoints = map_ref.grid; // cell grid points
        bin         = map_ref.map_dim;
        Nmin        = map_ref.origin;
        for (int i = 0; i < DIM; i++)
        {
            box[i]             = map_ref.cell[i];                   // cell size
            res[i]             = map_ref.cell[i] / map_ref.grid[i]; // resolution
            latticeOrigin_[i]  = map_ref.origin[i] * res[i];        // origin of the lattice
        }
    }
    else // no reference map provided, read variables the -size, -min, -min ,
         // -box options
    {
        for (int i = 0; i < DIM; i++)
        {
            // number of bins (working lattice)
            bin[i] = Nmax[i] - Nmin[i] + 1;

            // resolution
            res[i] = box[i] / nGridPoints[i];

            // origin of the lattice
            latticeOrigin_[i] = Nmin[i] * res[i];
        }
    }

    // ====  initialize the maps =====
    std::array<float, 2 *DIM> cell = {
        {static_cast<float>(box[XX]), static_cast<float>(box[YY]),
         static_cast<float>(box[ZZ]), 90., 90., 90.}
    };
    // N  Number of grid points within box
    //  MIN[i]; // Origin of the working lattice
    // bin = map_dim is the working lattice
    // spacegroup=1; // this will be assigned arbitrarily. It won't be used in
    // gtide anyways

    printf("=== Initial MD map parameters:\n");
    init_map(&map_av, "map", cell, nGridPoints, Nmin, bin, 0, 0,
             0, 0);
    printMapInfo(map_av);
    // allocate memory for the lattices
    map_av.vox.resize(map_av.map_dim[XX]*map_av.map_dim[YY]*map_av.map_dim[ZZ]);
    map_inst     = map_av;
    map_sd       = map_av;
    map_block    = map_av;

    /* read mask map if provided, else set mask to always true */
    if (!fnMaskMap_.empty())
    {
        printf("=== Reading mask map\n");
        t_mapdata map_mask;
        read_xplor_3dmap(fnMaskMap_.c_str(), &map_mask);
        mask_ = Mask(map_mask.vox, mask_cutoff);
    }

    // === Definition of the cutoff and other gaussian-density related variables
    ncutoff[XX] = floor(cutoff / res[XX]);
    ncutoff[YY] = floor(cutoff / res[YY]);
    ncutoff[ZZ] = floor(cutoff / res[ZZ]);
    // Density will be increment in the grid points between:   (i-ncutoff,
    // j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)

    Nframes_block = 0;

    if (!corr_.empty())
    {
        pfile_corr = fopen(corr_.c_str(), "w");

        // header for xmgrace
        // To do: print properly the header with the addModule option
        fprintf(pfile_corr, "@    title \"Time-resolved correlation between "
                "computed and reference maps\"\n");
        fprintf(pfile_corr,
                "@    subtitle  \"Time-window for block-averaging: %f ps\"\n",
                time_window);
        fprintf(pfile_corr, "@    xaxis  label \"Time\"\n");
        fprintf(pfile_corr, "@    yaxis  label \"Correlation\"\n");
    }
}

void TimeAveragedDensity::analyzeFrame(int frnr, const t_trxframe &fr,
                                       gmx_unused t_pbc *pbc,
                                       TrajectoryAnalysisModuleData *pdata)
{

    // === STUFF TO READ THE TRAJECTORY AND ATOM SELECTION DATA ===
    AnalysisDataHandle dh = pdata->dataHandle(data_);
    dh.startFrame(frnr, fr.time);
    size_t             group = 0; // selection-group index (only one considered thus it is 0)
    const Selection   &sel   = pdata->parallelSelection(sel_[group]);

    // === ACCUMMULATE DENSITY===
    accumulate_denstiy(&map_av, &map_inst, &map_sd, sel, res, latticeOrigin_, ncutoff, cutoff,
                       gaussianParametersPerAtomType, attypeString);

    // === ACCUMULATE BLOCK AVERAGE ===
    containeroperation::add(&(map_block.vox), map_inst.vox);

    // === Increment the number of frames variables ===
    Nframes_block++;

    // === calculate and print global correlation if requested ===
    if (!corr_.empty())
    {
        if (fr.time / time_window - floor(fr.time / time_window) == 0 &&
            fr.time > 0)
        {
            containeroperation::divide(&(map_block.vox), Nframes_block);
            float cc = global_correlation_coeff(map_ref, map_block, mask_);
            fprintf(pfile_corr, "%f %f\n", fr.time, cc);
            containeroperation::setZero(&(map_block.vox));
            Nframes_block = 0;
        }
    }
    dh.finishFrame();
}

void TimeAveragedDensity::finishAnalysis(int nframes)
{

    fprintf(stdout, "Nframes read: %d\n", nframes);

    // calculate timed-average and timed-stdev of the cumulated map
    containeroperation::divide(&(map_av.vox), nframes);
    sd_map(map_av, &map_sd, nframes);

    // calculate spatial average and spatial stdev
    setMapMeanAndRms(&map_av, incl_zero, mask_);
    setMapMeanAndRms(&map_sd, incl_zero, mask_);

    // calculate the difference map_diff = map_av - map_ref to the reference map (if requested)
    if (!fnDifferenceMap_.empty() || !rho_diff_.empty())
    {
        t_mapdata map_diff = map_av;
        map_diff.vox = subtract_maps(map_ref.vox, map_av.vox);
        // include zero in the calculation of mean and rms
        setMapMeanAndRms( &map_diff, true, mask_);

        if (!fnDifferenceMap_.empty())
        {
            print_to3dmap(fnDifferenceMap_, map_diff, mask_);
        }

        if (!rho_diff_.empty())
        {
            print_rho_sigma(rho_diff_, res, map_diff, map_diff);
        }

    }

    // calculate minus log if requested
    if (!mlog_.empty() || !rho_mlog_.empty())
    {
        t_mapdata map_mlog = map_av;
        // with ref map, map_mlog = -ln( map_av / map_ref ), else -ln(map_av)
        if (!fnInputMap_.empty())
        {
            map_mlog.vox = negLogRatioVector(map_ref.vox, map_av.vox);
        }
        else
        {
            map_mlog.vox = negLogVector(map_av.vox);
        }
        // do include zero in the calculation of mean and rms
        setMapMeanAndRms( &map_mlog, true, mask_);

        // shift map if requested
        if (shift_mlog)
        {
            containeroperation::addScalar(&(map_mlog.vox), -map_mlog.mean);
            setMapMeanAndRms(&map_mlog, true, mask_);
        }

        if (!mlog_.empty())
        {
            print_to3dmap(mlog_, map_mlog, mask_);
        }

        if (!rho_mlog_.empty())
        {
            print_rho_sigma(rho_mlog_, res, map_mlog, map_mlog);
        }

    }

    // print global correlation if there is a reference map
    if (!fnInputMap_.empty())
    {
        float cc = global_correlation_coeff(map_ref, map_av, mask_);
        printf("The correlation coefficient between computed and reference maps is "
               "%f\n",
               cc);
        // calculate local correlation if requested
        if (!loc_corr_.empty())
        {
            auto map_loc_corr = map_av;
            local_correlation_coeff(map_ref, map_av, mask_, &map_loc_corr,
                                    Nslide);
            map_loc_corr.mean = cc; // mean will be the global correlation coefficient
            map_loc_corr.rms  = nanf("");
            print_to3dmap(loc_corr_, map_loc_corr, mask_);
        }
    }

    // shift the map by an amount equal to the average
    if (Norm)
    {
        containeroperation::addScalar(&(map_av.vox), -map_av.mean);
        containeroperation::addScalar(&(map_sd.vox), -map_sd.mean);

        // Normalize the map by an amount equal to the standard deviation
        containeroperation::divide(&(map_av.vox), map_av.rms);
        containeroperation::divide(&(map_sd.vox), map_sd.rms);
        map_av.mean = -1.0; // fixed values to indicate stdev units around the mean
        map_sd.mean = -1.0;
        map_av.rms  = -1.0;
        map_sd.rms  = -1.0;
    }
}

void TimeAveragedDensity::writeOutput()
{
    // print the average to file 3dmap_average
    if (!fnMapAverage_.empty())
    {
        print_to3dmap(fnMapAverage_, map_av, mask_);
    }

    // print the standard deviation to file 3dmap_sigma
    if (!fnMapStandardDeviation_.empty())
    {
        print_to3dmap(fnMapStandardDeviation_, map_sd, mask_);
    }

    // print to rho.dat file
    if (!rho_.empty())
    {
        print_rho_sigma(rho_, res, map_av, map_sd);
    }

    // close time-averaged correlation printing output file
    if (!corr_.empty())
    {
        fclose(pfile_corr);
    }

}

}       // namespace

const char TimeAveragedDensityInfo::name[]             = "tide";
const char TimeAveragedDensityInfo::shortDescription[] =
    "Calculate TIme averaged electron DEnsities";

TrajectoryAnalysisModulePointer TimeAveragedDensityInfo::create()
{
    return TrajectoryAnalysisModulePointer(new TimeAveragedDensity);
}

} // namespace analysismodules

} // namespace gmx
