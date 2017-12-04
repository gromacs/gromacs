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
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "tide.h"

#include <map>
#include <string.h>

#include "gromacs/applied-forces/densityfitting/densfit.h"
#include "gromacs/applied-forces/densityfitting/densityspreader.h"
#include "gromacs/math/container/mask.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/encompassinggrid.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/fileio/griddataio.h"
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
    // need comparsion operator for inclusion in map
    bool operator<(const GaussianParameters &other) const
    {
        return (other.amplitude < amplitude) || ((other.amplitude == amplitude) && (other.precision < precision));
    }
    float amplitude;
    float precision;
};

//! \brief TODO
float calc_sum_rhoA_rhoB_mask(const GridDataReal3D::container_type &vec1,
                              const GridDataReal3D::container_type &vec2,
                              const Mask                           &mask)
{
    double     sum_ijk    = 0.0;
    const auto minMapSize = std::min(vec1.size(), vec2.size());
    for (size_t l = 0; l < minMapSize; l++)
    {
        if (mask[l])
        {
            sum_ijk += vec1[l] * vec2[l];
        }
    }
    return sum_ijk;
}

//! \brief TODO
float global_correlation_coeff(const GridDataReal3D::container_type &map1,
                               const GridDataReal3D::container_type &map2,
                               const Mask                           &mask)
{
    float numerator   = calc_sum_rhoA_rhoB_mask(map1, map2, mask);
    float denominator = sqrt(calc_sum_rhoA_rhoB_mask(map2, map2, mask) *
                             calc_sum_rhoA_rhoB_mask(map1, map1, mask));
    return numerator / denominator;
}

//! \brief Calculate the local correlation
void local_correlation_coeff(const GridDataReal3D &map1, const GridDataReal3D &map2, const Mask &mask,
                             GridDataReal3D *loc_corr, int Nodd)
{
    int Nhalf = static_cast<int>(Nodd / 2); // center point of the sliding grid (e.g. if Nodd=5
    // ; then Nhalf=2, Note that it starts counting
    // from 0)

    // Initialize the output loc_corr map, setting the borders and masked points,
    // with non-assigned NAN values
    std::fill(std::begin(*loc_corr), std::end(*loc_corr), nanf(""));
    const auto &mapExtend = map1.getGrid().lattice().extend();
    const auto &lattice   = map1.getGrid().lattice();
    // loop over the map1 reduced by Nhalf at both ends
    for (int i = Nhalf; i < mapExtend[XX] - Nhalf; i++)
    {
        for (int j = Nhalf; j < mapExtend[YY] - Nhalf; j++)
        {
            for (int k = Nhalf; k < mapExtend[ZZ] - Nhalf; k++)
            {
                ColumnMajorLattice<DIM>::MultiIndex currentIndex({{i, j, k}});

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
                            ColumnMajorLattice<DIM>::MultiIndex localSurroundingIndex({{i_sli, j_sli, k_sli}});
                            if (mask[lattice.lineariseVectorIndex(localSurroundingIndex)]) // skip masked points
                            {
                                const auto v1 = *(map1.iteratorAtMultiIndex(localSurroundingIndex));
                                const auto v2 = *(map2.iteratorAtMultiIndex(localSurroundingIndex));
                                sum12 += v1 * v2;
                                sum11 += v1 * v1;
                                sum22 += v2 * v2;
                            }
                        }
                    }
                }

                *(loc_corr->iteratorAtMultiIndex(currentIndex)) = sum12 / sqrt(sum22 * sum11);
            }
        }
    }
}

//! \brief Set the Gaussian coefficients A and B for the six atom types H, C, N,
//! O, P, and S
std::map < std::string, std::vector < GaussianParameters>> set_gaussian_coeff( std::string parameterFileName )
{
    // read gaussian parameters from file if provided, else use defaults
    FILE *parameterFile;
    if (!parameterFileName.empty())
    {
        fprintf(stdout, "Reading file with Gaussian coefficients\n");
        parameterFile = fopen(parameterFileName.c_str(), "r");
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

//! \brief Compute timed-average standard deviation density
void sd_map(const GridDataReal3D &map_av, GridDataReal3D *map_sd,
            float Nframes)
{
    if (Nframes > 1)
    {
        std::transform(std::begin(map_av), std::end(map_av),
                       std::begin(*map_sd), std::begin(*map_sd),
                       [Nframes](float av, float sd) {
                           return sqrt(sd / (Nframes - 1.0) -
                                       av * av * Nframes / (Nframes - 1.0));
                       });
    }
}

//! \brief Print the array in dat format: x, y, z, rho, sigma
void print_rho_sigma(const std::string    &inputFileName,
                     const GridDataReal3D &map_av,
                     const GridDataReal3D &map_sd)
{
    FILE       *pfile;
    pfile = fopen(inputFileName.c_str(), "w");
    const auto &extend = map_av.getGrid().lattice().extend();
    const auto &res    = map_av.getGrid().unitCell().basisVectorLengths();
    // === line 1 : Header lines ===
    fprintf(pfile, "# Rho ( and its standard dev) vs x y z.\n");
    fprintf(pfile, "# grid_res= %f %f %f (nm), dimensions: bin[XX]: %d bin[YY]: "
            "%d bin[ZZ]: %d\n",
            res[XX], res[YY], res[ZZ], extend[XX], extend[YY],
            extend[ZZ]);
    fprintf(pfile, "# x y z (Angstrom) rho sigma (Angstrom^-3)\n");

    // === Density array ===
    ColumnMajorLattice<DIM>::MultiIndex multiIndex;
    for (multiIndex[XX] = 0; multiIndex[XX] < extend[XX]; multiIndex[XX]++)
    {
        for (multiIndex[YY] = 0; multiIndex[YY] < extend[YY]; multiIndex[YY]++)
        {
            for (multiIndex[ZZ] = 0; multiIndex[ZZ] < extend[ZZ]; multiIndex[ZZ]++)

            {
                const auto &r = map_av.getGrid().multiIndexToCoordinate(multiIndex);
                fprintf(pfile, "%f %f %f %d %d %d %E %E\n", NM2AA * r[XX], NM2AA * r[YY], NM2AA * r[ZZ], multiIndex[XX], multiIndex[YY], multiIndex[ZZ],
                        *(map_av.iteratorAtMultiIndex(multiIndex)), *(map_sd.iteratorAtMultiIndex(multiIndex)));
            }
        }
    }
    fclose(pfile);
}


struct MapDimensions
{
    std::array<int, DIM> nGridPoints = {{250, 250, 250}};
    std::array<int, DIM> Nmin        = {{0, 0, 0}};       // (AMIN,BMIN,CMIN)
    std::array<int, DIM> Nmax        = {{249, 249, 249}}; // (AMAX,BMAX,CMAX)
    RVec                 box         = {10., 10., 10.};   //(boxA,boxB,boxC)
};


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
        std::map < std::string, std::vector < GaussianParameters>> setGaussianCoeffFromSigma( );
        class ModuleData;
        Grid<DIM>     referenceCell; //< The cell read from the xplor file is lost in the conversion to GridDataReal3D
        MapDimensions mapDimensions;
        // entire cell


        // subcell where the density is actually computed
        float                            mask_cutoff = 1e-5; //|rho|<mask_rheshold will be ignored from the calculation
        float                            time_window = 1.0;  // for time-resolved global correlation
        float                            sigma_      = 0.3;  // spread with in nm
        bool                             hasSigma_   = false;
        int                              Nslide      = 5;    // Size of the sliding grid for local-correlation calculation
        float spacing_ = 0.2;                                // grid spacing in nm
        float margin_  = 0.5;                                // grid margin around atoms in nm
        Selection                        refsel_;
        SelectionList                    sel_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;

        // === Definition of the parameters for the input/output files ===

        std::string              fnGaussCoefficients_;    // input structure factor file
        std::string              fnMapAverage_;           // average map (3d)
        std::string              fnMapStandardDeviation_; // standard deviation (3d)
        std::string              fnInputMap_;
        std::string              fnDifferenceMap_;
        std::string              fnMaskMap_;
        std::string              mlog_;
        std::string              loc_corr_;
        std::string              rho_; // density (3d), same as the map but in a differnt format: // x, y, z, rho, sigma
        std::string              rho_diff_;
        std::string              rho_mlog_;
        std::string              corr_;

        bool                     Norm                  = true;                                     // Normalization
        bool                     incl_zero             = false;                                    // include zero-density values in the spatial average
                                                                                                   // and stdev calculation (over the grid points)
        bool                     shift_mlog        = true;                                         // remove mean value from mlog map
        bool                     integrateDensity_ = false;
        std::map < std::string, std::vector < GaussianParameters>> gaussianParametersPerAtomType_; // Ai factors (Ngaussians x Natoms 1D array)
        std::vector<std::string> attypeString;                                                     // array will store the index of the

        // === define the lattice that will store the density for each point===
        GridDataReal3D map_ref;             /* The reference=experimental map to compare to        */
        GridDataReal3D map_av;              /* The computed timed-averaged map         */
        GridDataReal3D map_sd;              /* The computed timed-averaged standard deviation map */
        GridDataReal3D map_block;           /* map for block averaging*/
        Mask           mask_;               //< Masking out density values

        void writeMaskedMap(const std::string &outputFileName, const GridDataReal3D &map,
                            const Mask &mask);

        // === Variables to be used in the trajectory loop ===
        float Nframes_block = -1; // Nframes in the time_window for block averaging
};

//! \brief Set the Gaussian coefficients for cryo-EM scattering,
std::map < std::string, std::vector < GaussianParameters>> TimeAveragedDensity::setGaussianCoeffFromSigma( )
{
    // read gaussian parameters from file if provided, else use defaults
    auto parameterFile = libopen("cryoemcrosssections.dat");
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
            gaussianParametersLookup[currentAtomType].push_back({amplitude, static_cast<float>(precision *(1/(2.*sigma_*sigma_)))});
        }
    }
    fclose(parameterFile);
    return gaussianParametersLookup;
}


/*! \brief Outputs the density map.
 */
void TimeAveragedDensity::writeMaskedMap(const std::string &outputFileName, const GridDataReal3D &map,
                                         const Mask &mask)
{
    if (outputFileName.empty())
    {
        InconsistentInputError("No filename given for writing masked file.");
    }

    auto mapCopy = map; // copy map before masking
    mask.maskToNan(&mapCopy);
    MrcFile().write(outputFileName, mapCopy);
}


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
        "and EM, as well as for coarse-grained systems, are available at XXXX\n",

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
        "10 -size 250 250 250 -min 0 0 0 -max 249 249 249 -Ngaussians 4 "
        " -select \'group \"Protein\"\' -map_av\n",
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

    options->addOption(FileNameOption("Gauss_coef").filetype(eftGenericData)
                           .inputFile().store(&fnGaussCoefficients_)
                           .defaultBasename("Gauss_coef")
                           .description("Gaussian coefficients"));

    options->addOption(FileNameOption("mi").filetype(eftCCP4).inputFile()
                           .store(&fnInputMap_).defaultBasename("mi")
                           .description("Reference map."));

    options->addOption(FileNameOption("mask").filetype(eftCCP4).inputFile()
                           .store(&fnMaskMap_).defaultBasename("mask")
                           .description("Input mask map."));

    options->addOption(FileNameOption("map_av").filetype(eftCCP4).outputFile()
                           .store(&fnMapAverage_).defaultBasename("map_av")
                           .description("Average density map."));

    options->addOption(FileNameOption("map_sd").filetype(eftCCP4).outputFile()
                           .store(&fnMapStandardDeviation_).defaultBasename("map_sd")
                           .description("St. dev. density map in X-PLOR format"));

    options->addOption(FileNameOption("diff").filetype(eftCCP4).outputFile()
                           .store(&fnDifferenceMap_).defaultBasename("diff")
                           .description("Difference densiy map format"));

    options->addOption(FileNameOption("mlog").filetype(eftCCP4).outputFile()
                           .store(&mlog_).defaultBasename("mlog")
                           .description("-ln(map_computed/map_reference) format"));

    options->addOption(FileNameOption("corr").filetype(eftPlot).outputFile()
                           .store(&corr_).defaultBasename("corr")
                           .description("Time-resolved global correlation "
                                        "between computed and reference maps"));

    options->addOption(FileNameOption("loc_corr").filetype(eftCCP4).outputFile()
                           .store(&loc_corr_).defaultBasename("loc_corr")
                           .description("Local correlation between computed and "
                                        "reference maps format"));

    options->addOption(SelectionOption("select").storeVector(&sel_).required()
                           .description("Group to calculate the density"));

    options->addOption(FloatOption("mask_cutoff").store(&mask_cutoff)
                           .description("Points with density smaller than mask_cutoff will be masked according to mask map (-mask)"));

    options->addOption(BooleanOption("norm").store(&Norm)
                           .description( "Normalize the output maps. If yes, map in stdev units and shifted such that mean=0"));

    options->addOption(BooleanOption("integrate").store(&integrateDensity_)
                           .description("Integrate density over voxel Volumes instead of default Gauss transform."));

    options->addOption( BooleanOption("incl_zero").store(&incl_zero)
                            .description("Include Zero-density points in the calculation of the "
                                         "spatial mean and stdev, over the grid"));

    options->addOption(BooleanOption("shift_mlog").store(&shift_mlog)
                           .description("remove mean value from mlog map"));
    options->addOption(FloatOption("spacing").store(&spacing_)
                           .description("Map spacing in nm."));

    options->addOption(FloatOption("time_window").store(&time_window)
                           .description("Time window for block-averaging in time-resolved "
                                        "global correlation calculation(ps)"));
    options->addOption(FloatOption("sigma").store(&sigma_).storeIsSet(&hasSigma_)
                           .description("Spreadwidth for atoms in nm."));

    options->addOption(IntegerOption("Nslide").store(&Nslide)
                           .description("Size of the sliding grid for local correlation calculation (Must be an odd number)"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}

void TimeAveragedDensity::initAnalysis(
        const gmx_unused TrajectoryAnalysisSettings &settings,
        const TopologyInformation                   &top)
{

    data_.setColumnCount(0, sel_.size());
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    // === GAUSSIAN COEFFICIENTS BASED ON ATOM TYPES===
    if (hasSigma_)
    {
        gaussianParametersPerAtomType_ = setGaussianCoeffFromSigma();
    }
    else
    {
        gaussianParametersPerAtomType_ = set_gaussian_coeff(fnGaussCoefficients_);
    }

    attypeString                  = attypeStrings_from_top(*(top.topology()));
    std::unique_ptr < IGrid < DIM>> referenceGrid;
    // create grid from reference map if provided, alternatively from structure
    if (!fnInputMap_.empty())
    {
        map_ref       = MrcFile().read(fnInputMap_);
        referenceGrid = map_ref.getGrid().duplicate();
    }
    else
    {
        rvec      * coordinates;
        matrix      box;
        top.getTopologyConf(&coordinates, box);
        t_topology *topol = top.topology();
        referenceGrid = encompassingGridFromCoordinates({coordinates, coordinates+topol->atoms.nr}, spacing_, margin_).duplicate();
    }

    map_av.setGrid(referenceGrid->duplicate());
    map_sd.setGrid(referenceGrid->duplicate());
    map_block.setGrid(referenceGrid->duplicate());

    /* read mask map if provided */
    if (!fnMaskMap_.empty())
    {
        mask_ = Mask(MrcFile().read(fnMaskMap_), mask_cutoff);
    }

    Nframes_block = 0;

    if (!corr_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(corr_);
        plotm->setTitle("Time-resolved correlation between "
                        "computed and reference maps");
        plotm->setSubtitle("Time-window for block-averaging:" + std::to_string(time_window)+" ps" );
        plotm->setXAxisIsTime();
        plotm->setYLabel("Correlation");
        data_.addModule(plotm);
    }
}

void TimeAveragedDensity::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /*pbc*/, TrajectoryAnalysisModuleData *pdata)
{

    AnalysisDataHandle dh = pdata->dataHandle(data_);
    dh.startFrame(frnr, fr.time);

    const size_t       group = 0; // selection-group index (only one considered thus it is 0)
    const Selection   &sel   = pdata->parallelSelection(sel_[group]);

    // Gather all atoms with the same Gaussian Parameters, also considering multiple Gaussian Parameters per atom
    std::map < GaussianParameters, std::vector < RVec>> atomIndexPerParameter;
    for (int atomIndex = 0; atomIndex < sel.atomCount(); atomIndex++)
    {
        // Select the first character of the atom name
        std::string atomName =
            attypeString[sel.atomIndices()[atomIndex]].substr(0, 1);
        const auto &gaussianParameter = gaussianParametersPerAtomType_.at(atomName);
        for (const auto &gp : gaussianParameter)
        {
            atomIndexPerParameter[gp].push_back(sel.position(atomIndex).x());
        }
    }

    DensitySpreader spreader(map_av.getGrid(), std::max(1, gmx_omp_nthreads_get(emntDefault)), 6, 1., integrateDensity_);

    for (const auto &para : atomIndexPerParameter)
    {
        if (para.first.amplitude > 0 && para.first.precision > 0)
        {
            const auto        sigma = sqrt(1/(2*para.first.precision));
            spreader.setSigma(sigma);
            std::vector<real> weights(para.second.size(), para.first.amplitude * sqrt(8*M_PI*M_PI*M_PI*sigma*sigma*sigma));
            ArrayRef<real>    weightsRef {
                weights.data(), weights.data()+weights.size()
            };
            PaddedArrayRef<const RVec> xRef {
                para.second.data(), para.second.data()+para.second.size()
            };
            spreader.spreadLocalAtoms(xRef, weightsRef);
        }
    }

    // ===  update average and stdev maps ===
    containeroperation::add(&map_av, spreader.getSpreadGrid());
    containeroperation::add(&map_block, spreader.getSpreadGrid());
    Nframes_block++;
    std::transform(std::begin(map_sd), std::end(map_sd), std::begin(spreader.getSpreadGrid()), std::begin(map_sd), [](float oldSD, float voxelValue) { return oldSD + voxelValue * voxelValue; });

    // === calculate and print global correlation if requested ===
    if (!corr_.empty())
    {
        if (fr.time / time_window == floor(fr.time / time_window) && fr.time > 0)
        {
            containeroperation::divide(&map_block, Nframes_block);
            dh.setPoint(0, global_correlation_coeff(map_ref, map_block, mask_));
            containeroperation::setZero(&map_block);
            Nframes_block = 0;
        }
    }

    dh.finishFrame();
}

void TimeAveragedDensity::finishAnalysis(int nframes)
{

    fprintf(stdout, "Nframes read: %d\n", nframes);

    // calculate timed-average and timed-stdev of the cumulated map
    containeroperation::divide(&map_av, nframes);
    sd_map(map_av, &map_sd, nframes);

    // calculate the difference map_diff = map_av - map_ref to the reference map (if requested)
    if (!fnDifferenceMap_.empty() || !rho_diff_.empty())
    {
        GridDataReal3D map_diff(map_av.getGrid());
        map_diff = containermeasure::subtract(map_ref, map_av);
        // include zero in the calculation of mean and rms

        if (!fnDifferenceMap_.empty())
        {
            writeMaskedMap(fnDifferenceMap_, map_diff, mask_);
        }

        if (!rho_diff_.empty())
        {
            print_rho_sigma(rho_diff_,  map_diff, map_diff);
        }

    }

    // calculate minus log if requested
    if (!mlog_.empty() || !rho_mlog_.empty())
    {
        GridDataReal3D map_mlog(map_av.getGrid());
        // with ref map, map_mlog = -ln( map_av / map_ref ), else -ln(map_av)
        if (!fnInputMap_.empty())
        {
            // TODO: mask division by zero first
            map_mlog = containermeasure::negLogRatio(map_ref, map_av);
        }
        else
        {
            // TODO: mask division by zero first
            map_mlog = containermeasure::negLog(map_av);
        }
        // do include zero in the calculation of mean and rms

        // shift map if requested
        if (shift_mlog)
        {
            containeroperation::addScalar(&map_mlog, -containermeasure::mean(map_mlog));
        }

        if (!mlog_.empty())
        {
            writeMaskedMap(mlog_, map_mlog, mask_);
        }

        if (!rho_mlog_.empty())
        {
            print_rho_sigma(rho_mlog_, map_mlog, map_mlog);
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
            writeMaskedMap(loc_corr_, map_loc_corr, mask_);
        }
    }

    // rescale map to zero average and unit standard deviation
    if (Norm)
    {
        containeroperation::addScalar(&map_av, -containermeasure::mean(map_av));
        containeroperation::divide(&map_av, containermeasure::rms(map_av));

        containeroperation::addScalar(&map_sd, -containermeasure::mean(map_sd));
        containeroperation::divide(&map_sd, containermeasure::rms(map_sd));
    }
}

void TimeAveragedDensity::writeOutput()
{
    // print the average to file 3dmap_average
    if (!fnMapAverage_.empty())
    {
        writeMaskedMap(fnMapAverage_, map_av, mask_);
    }

    // print the standard deviation to file 3dmap_sigma
    if (!fnMapStandardDeviation_.empty())
    {
        writeMaskedMap(fnMapStandardDeviation_, map_sd, mask_);
    }

    // print to rho.dat file
    if (!rho_.empty())
    {
        print_rho_sigma(rho_, map_av, map_sd);
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
