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

#include <string>
#include <vector>

#include "gromacs/trajectoryanalysis.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


namespace gmx
{

namespace analysismodules
{

namespace
{

//! \brief Gets the index in an uni-dimensional array from the three coordinates of a three-dimensional array
int index_array(int i, int j, int k, int Nx, int Ny, gmx_unused int Nz)
{
// The unidimensional array contains Nx*Ny*Nz elements:  0 ... Nx*Ny*Nz-1
// The element in the unidimensional array is given by
    return k*(Nx*Ny) + j*Nx + i;
}


/*!  \brief Outputs the density map to text file using XPLOR format
 *
 * For the format see: http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
void print_to3dmap(const char *N_file, rvec N, rvec Nmin, rvec Nmax, rvec BOX, ivec bin, rvec res, float density[])
{
    if (N_file == nullptr)
    {
        fprintf(stdout, "No filename given for X-PLOR output file, not outputting densities\n");
        return;
    }

    FILE *fp = gmx_ffopen(N_file, "w");

    // === line 1 : Empty line ===
    fprintf(fp, "\n");

    // === Lines 2-5 : Title and remarks ===
    fprintf(fp, "       2\n");
    fprintf(fp, "REMARKS Lattice Nbins: %d %d %d\n", bin[XX], bin[YY], bin[ZZ]);
    fprintf(fp, "REMARKS Lattice resolution (nm) %f %f %f\n", res[XX], res[YY], res[ZZ]);

    // === Line 6 : Grid points ===
    //  nBOX = total number of grid points along the a,b, and c cell edges (of the simulation box)

    //  starting and stopping grid points along each cell edge in the portion of the map that is written
    int NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX;

    NA = floor(N[XX]);
    NB = floor(N[YY]);
    NC = floor(N[ZZ]);

    AMIN = floor(Nmin[XX]);
    BMIN = floor(Nmin[YY]);
    CMIN = floor(Nmin[ZZ]);

    AMAX = floor(Nmax[XX]);
    BMAX = floor(Nmax[YY]);
    CMAX = floor(Nmax[ZZ]);


    fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX  );

    // === Line 7 : BOX size (Angstrom) and angles ===
    fprintf(fp, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n", BOX[XX]*10.0, BOX[YY]*10.0, BOX[ZZ]*10.0, 90.0, 90.0, 90.0   );

    // === Line 8 : string ===
    fprintf(fp, "ZYX\n");

    // === Density array ===
    int    i, j, k;

    double average = 0, sigma = 0, counts = 0; // average density and standard deviation

    for (k = 0; k < bin[ZZ]; k++)              // loop over z
    {                                          // Print the k value
        fprintf(fp, "       %d\n", k);

        // variable to monitor the number of columns printed
        int ncols = 0;

        for (j = 0; j < bin[YY]; j++)           // loop over y

        {
            for (i = 0; i < bin[XX]; i++)           // loop over x

            {
                int index_aux = index_array(i, j, k, bin[XX], bin[YY], bin[ZZ]);

                // Calculation of the average and the sigma density
                // if density is not zero
                if (density[index_aux] > 0)
                {
                    average += density[index_aux];
                    sigma   += density[index_aux]*density[index_aux];
                    counts++;
                }


                // Print the density values in 6 columns
                fprintf(fp, "%12.5E", density[index_aux] );
                ncols++;

                // enter when the number of printed columns =6
                if (ncols == 6 || ( i == (bin[XX]-1) && j == (bin[YY]-1) ) )
                {
                    fprintf(fp, "\n");
                    ncols = 0;
                }
            }
        }
    }

    //  === Footer lines
    fprintf(fp, "   -9999\n");

    // Calculate and print the average and the standard deviation
    average /= counts;
    sigma    = sqrt( sigma/(counts-1) - counts*average*average/(counts-1));

    fprintf(fp, "%13.5E%13.5E\n", average, sigma);

    gmx_ffclose(fp);
}

//! \brief Print the array in dat format: x, y, z, rho, sigma
void print_rho_sigma(const char *N_file, ivec bin, rvec res,  rvec min, float density[], float density_sigma[] )
{
    // Nbins of the lattice: (bin[XX],bin[YY],bin[ZZ])
    // Resolution (size of each grid) (resx,resy,resz) (nm)
    // density array = density []
    // standard deviation array = density_sigma []
    FILE *pfile;
    pfile = fopen(N_file, "w");

    // === line 1 : Header lines ===
    fprintf(pfile, "# Rho ( and its standard dev) vs x y z.\n");
    fprintf(pfile, "# grid_res= %f %f %f (nm), nbins: bin[XX]: %d bin[YY]: %d bin[ZZ]: %d\n", res[XX], res[YY], res[ZZ], bin[XX], bin[YY], bin[ZZ]);
    fprintf(pfile, "# x y z (Angstrom) rho sigma (Angstrom^-3)\n");


    // === Density array ===
    int   i, j, k;
    float x, y, z;


    for (i = 0; i < bin[XX]; i++)         // loop over x
    {
        for (j = 0; j < bin[YY]; j++)     // loop over y
        {
            for (k = 0; k < bin[ZZ]; k++) // loop over z

            {                             //=== set the absolute coordinate system (the simulation box)by adding (minx,miny,minz) to each lattice value
                // x,y,z in Angstrom
                x = 10*(i+min[XX])*res[XX];
                y = 10*(j+min[YY])*res[YY];
                z = 10*(k+min[ZZ])*res[ZZ];

                int index_aux = index_array(i, j, k, bin[XX], bin[YY], bin[ZZ]);

                // Print the density values only if they are > 0
                if (density[index_aux] > 0)
                {
                    fprintf(pfile, "%f %f %f %d %d %d %E %E\n", x, y, z, i, j, k, density[index_aux], density_sigma[index_aux] );
                }
            }
        }
    }
    fclose (pfile);
}


/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class TimeAveragedDensity : public TrajectoryAnalysisModule
{
    public:
        TimeAveragedDensity();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();


    private:
        class ModuleData;

        // entire cell
        rvec N;    // (NA,NB,NC)
        rvec Nmin; // (AMIN,BMIN,CMIN)
        rvec Nmax; // (AMAX,BMAX,CMAX)
        rvec BOX;  //(BOXA,BOXB,BOXC)

        // subcell where the density is actually computed
        ivec  bin;       // Grid points (0 0 0) (working lattice)
        rvec  res;       // grid resolution  (box/N)
        rvec  rmax;      // Lattice dimensions / 2
        rvec  L_lattice; // Lattice dimensions (L_lattice= 2*rmax)
        rvec  L0;        // origin of the lattice with respect to the origin of the map Nmin*res
        float cutoff;    // (nm) density will be incremented in the grid points within a cutoff distance of a certain atom
        float sigmaf;    // Assume that the density of an atom is gaussian distributed.


        // ==== selection variables ===
        Selection                        refsel_;
        SelectionList                    sel_;
        int     nr;     // number of atoms in the selection

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;

        // === Definition of the parameters for the input/output files ===

        // Input and output file names
        std::string                      gauss_coef_;
        std::string                      map_av_; // average map (3d)
        std::string                      map_sd_; // standard deviation (3d)
        std::string                      rho_;    // density (3d), same as the map but in a different format: x, y, z, rho, sigma
        char *gauss_coef;                         // input structure factor file

        // === structure factor variables ===
        char   gauss_coef_atname[10]; // atom names for crosscheck
        char   gauss_coef_attype[10]; // atom types
        int    Ngaussians;            // Number of Gaussians that add up to the density
        float *Ai;                    // Ai factors (Ngaussians x Natoms 1D array)
        float *Bi;                    // Bi factors (Ngaussians x Natoms 1D array)
        float *C;                     // c factor

        // === parse trajectory ===
        int parse_traj;

        // === define the lattice that will store the density for each point===
        float  *lattice;         //structure to locate the index of the protein
        float  *lattice_old;     //lattice old for standard deviation calculation
        float  *lattice_squared; // define the lattice that will store the density squared for each point (for the further standard deviation calculation)
        ivec    ncutoff;         // === Definition of the cutoff and other gaussian-density related variables
        // Density will be increment in the grid points between:   (i-ncutoff, j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)





        // === Variables to be used in the trajectory loop ===
        float Nframes;           // Nframes
        int   index_aux;


        float m;          // mass of the n-th atom
        float X;          // X coordinate of the n-th atom
        float Y;          // Y coordinate of the n-th atom
        float Z;          // Z coordinate of the n-th atom
        float onesMtotal; // 1/Mtotal

        // Coordinate vector with respect to the simulation box
        rvec R;
        rvec R_trans;
        int  nx, ny, nz;  // integer variables to parse the grid point around (i,j,k)

        // vector coordinates of the (nx,ny,nz) grid point
        rvec grid_vec;

        // Distance between the atom and the grid_vector:  R_atom_grid= R_trans - grid_vec
        rvec  R_atom_grid;

        float dist_atom_grid;

        // === Gaussian contribution (two or four gaussians depending on the atom) to the density at this grid point
        float density_contrib;

        // gaussians parsing integer
        int g;



        // === AUXILIARY PARSING VARIABLES ===
        int i, j, k, n, index_oneD_aux, l;


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
        "rho(R)=C+sum_{i=1}^{i=Ngaussians} Ai*exp(-Bi*R^2)\n",

        "=== Input files:\n",
        "-trajectory\n",
        "- reference structure \n",
        "- index file with the desired group to compute the density\n",
        "- gauss_coeff.dat file with Ai and Bi constants for each atom in the selected group\n",
        "Format of Gauss_coef.dat:\n",
        "atomname atomtype A1 A2 A3 A4...ANgaussian B1 B2 B3 B4... BNGaussian C\n",
        "Ai, Bi, and C constants can be can be computed from the structure factors using the"
        "scripts provided in the web page: XXXX\n",
        "WARNING: No dynamic group selection supported yet.\n",

        "=== Output files:\n",
        "-map_av time average of density in XPLOR format.\n",
        "-map_sd time standard deviation of density in XPLOR format.\n",
        "The average and standard deviation of the densities over the whole considered space\n",
        "are computed (discarding points with zero density) and are written at the end of\n",
        "the XPLOR files.\n\n",
        "-rho = writes both the time average and standard deviation of the density\n",
        "for each position in ASCII format\n",


        "Parameters (similar parameters as in an XPLOR map):\n",
        "    --(NA,NB,NC) = total number of grid points along each cell edge size (-size)\n",
        "    --(BOXA, BOXB,BOXC) = crystal cell dimensions (nm,nm,nm) (-box)\n",
        "    --(AMIN,BMIN,CMIN) = starting grid point in the considered portion of the map (-min)\n",
        "    --(AMAX,BMAX,CMAX) = stopping grid point in the considered portion of the map (-max)\n",
        "    NOTE: The grid resolution is computed as BOXi/Ni, independently for each axis i=a,b,c.\n",
        "    NOTE: The density will be ONLY computed on a sub-lattice going from (AMIN,BMIN,CMIN) to (AMAX,BMAX,CMAX).\n",
        "    NOTE: The density is absolute with respect to the grid defined by the integer vector Ni. There is no fit.\n",
        "    If the density is intended to be calculated with respect to other reference group, a previous fitting step (e.g. with trjconv) should precede the gtide calculation.\n",
        "    Example: NA NB NC = 100 100 200 ; AMIN BMIN CMIN = 10 10 40  ; AMAX BMAX CMAX = 20 30 70\n",
        "    It will compute the density on a lattice of size (11,21,31) grid points \n",
        "    with its origin at the grid point (10,10,40), embedded in a whole grid of size (100,100,200) (grid points)",
        "    --cutoff: Cutoff distance to compute the density around a certain atom (nm)\n",
        "    --Ngaussians: Number of Gaussians to add to rho (e.g. 2,4 or 5)\n"
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("Gauss_coef")
                           .filetype(eftGenericData).inputFile()
                           .store(&gauss_coef_).defaultBasename("Gauss_coef")
                           .description("File with Gaussian coefficients"));


    options->addOption(FileNameOption("map_av")
                           .filetype(eftXPLOR).outputFile()
                           .store(&map_av_).defaultBasename("map_av")
                           .description("Average density map in X-PLOR format"));

    options->addOption(FileNameOption("map_sd")
                           .filetype(eftXPLOR).outputFile()
                           .store(&map_sd_).defaultBasename("map_sd")
                           .description("St. dev. density map in X-PLOR format"));

    options->addOption(FileNameOption("rho")
                           .filetype(eftGenericData).outputFile()
                           .store(&rho_).defaultBasename("rho")
                           .description("av and sd of density in ASCII format"));

    options->addOption(SelectionOption("select")
                           .storeVector(&sel_).required()
                           .description("Group to calculate the density"));


    options->addOption(RealOption("size").store(N).vector()
                           .description("(NA,NB,NC) = total number of grid points along each cell edge size"));

    options->addOption(RealOption("min").store(Nmin).vector()
                           .description("(Amin,Bmin,Cmin) = starting grid points in the considered portion of the map"));

    options->addOption(RealOption("max").store(Nmax).vector()
                           .description("(Amax,Bmax,Cmax) = stopping grid points in the considered portion of the map\nDensity will be calculated within min-max grid points"));

    options->addOption(RealOption("box").store(BOX).vector()
                           .description("crystal cell dimensions. Grid spacing computed as BOX/N for each axis"));

    options->addOption(FloatOption("cutoff").store(&cutoff)
                           .description("Cutoff distance to compute the density around a certain atom (nm)"));

    options->addOption(IntegerOption("Ngaussians").store(&Ngaussians)
                           .description("Number of Gaussians to use in the sum (e.g. 2,4, or 5)"));



    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void
TimeAveragedDensity::initAnalysis(const gmx_unused TrajectoryAnalysisSettings  &settings,
                                  const TopologyInformation                    &top)
{



// set parse_traj to one so the analysis is done
    parse_traj = 1;


    data_.setColumnCount(0, sel_.size());
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);


    // === NUMBER OF ATOMS IN THE SELECTION (nr)===
    size_t           g     = 0;
    const Selection &sel   = sel_[g];
    nr    = sel.posCount();
    fprintf(stdout, "Natoms in selection: %d\n", nr);


    // === GAUSSIAN COEFFICIENTS===

    fprintf(stdout, "Considering %d Gaussians\n", Ngaussians);

    // allocation of arrays
    snew(Ai, nr*Ngaussians); // a factors
    snew(Bi, nr*Ngaussians); // bfactors
    snew(C, nr);             // c factor

    gauss_coef = new char [gauss_coef_.size()+1];
    strcpy (gauss_coef, gauss_coef_.c_str());

    //Create an array with the structure factors;
    //Electron density is approximated with two or four gaussians:
    // a1*exp(-b1*s*s) + a2*exp(-b2*s*s) + a3*exp(-b3*s*s) + a4*exp(-b4*s*s) +  c


    //scan the gauss_coeftor file
    FILE *pfile;
    pfile = fopen ( gauss_coef, "r");
    //printf("%s\n",gauss_coef);
    int col;

    // loop over the lipid atoms
    for (n = 0; n < nr; n++)
    {
        // read the following parameters
        fscanf (pfile, "%s", gauss_coef_atname);         // atomname
        fscanf (pfile, "%s", gauss_coef_attype);         // atomtype
        for (col = 0; col < Ngaussians; col++)
        {
            index_oneD_aux = index_array(col, n, 0, Ngaussians, nr, 0);
            fscanf (pfile, "%f", &Ai[index_oneD_aux]);         // Ai factors
        }
        for (col = 0; col < Ngaussians; col++)
        {
            index_oneD_aux = index_array(col, n, 0, Ngaussians, nr, 0);
            fscanf (pfile, "%f", &Bi[index_oneD_aux]);         // Bi factors

        }
        fscanf (pfile, "%f", &C[n]);         // c factors


        // crosscheck that the atomname is the same than in the topology (or pdb) input file
        int atomindex = sel.atomIndices()[n];

        //Declare atom name from top
        char *atomname;
        atomname = *(top.topology()->atoms.atomname[ atomindex  ]);

        //fprintf(stdout, "n=%d atid=%d atname=%s\n",n,atomindex,atomname);

        // compare with the atom name in the gauss_coeft file
        int comparison = strcmp(atomname, gauss_coef_atname );
        if (comparison != 0)
        {
            fprintf(stderr, "\n == ERROR ===\n : name(top): %s != name(gauss_coef_file): %s \n", atomname, gauss_coef_atname);
            parse_traj = 0;
        }
/*
                // print crosschecks
                fprintf(stdout,"CROSSCHECK %s %s ",gauss_coef_atname,gauss_coef_attype) ;
                for(col=0;col<Ngaussians;col++) { index_oneD_aux=index_array(col,n,0,Ngaussians,nr,0); fprintf(stdout,"%e ",Ai[index_oneD_aux]) ; }
                for(col=0;col<Ngaussians;col++) { index_oneD_aux=index_array(col,n,0,Ngaussians,nr,0); fprintf(stdout,"%e ",Bi[index_oneD_aux]) ; }
                fprintf(stdout,"%e ",C[n]) ;
                fprintf(stdout,"\n");
 */
    }
    fclose (pfile);


    // === CREATE THE LATTICE ===


    // === number of grid points;
    for (i = 0; i < DIM; i++)
    {

        //number of bins (working lattice)
        bin[i] = floor(Nmax[i]-Nmin[i]+1);         //Amax-Amin+1

        //resolution
        res[i] = BOX[i]/N[i];

        // working lattice size
        rmax[i]      = bin[i]*res[i]/2.0;     // half the size of the lattice
        L_lattice[i] = bin[i]*res[i];         // complete size

        // origin of the lattice
        L0[i] = Nmin[i]*res[i];
    }


    // === Print the parameters that it found ===
    fprintf(stdout, "#Atoms found in the group: %d\n", nr);
    fprintf(stdout, "#NA AMIN AMAX   NB BMIN BMAX   NC CMIN CMAX: %f %f %f   %f %f %f   %f %f %f\n", floor(N[XX]), floor(Nmin[XX]), floor(Nmax[XX]), floor(N[YY]), floor(Nmin[YY]), floor(Nmax[YY]), floor( N[ZZ]), floor(Nmin[ZZ]), floor(Nmax[ZZ]));
    fprintf(stdout, "#BOX: %f %f %f\n", BOX[XX], BOX[YY], BOX[ZZ]);

    fprintf(stdout, "#Lattice where density is computed. Dimensions (nm): %f %f %f, Nbins per coordinate: %d %d %d, resolution per coordinate: %f %f %f, origin: %f %f %f\n", L_lattice[XX], L_lattice[YY], L_lattice[ZZ], bin[XX], bin[YY], bin[ZZ], res[XX], res[YY], res[ZZ], L0[XX], L0[YY], L0[ZZ]);
    fprintf(stdout, "#Cutoff (nm): %f\n", cutoff);

    //allocate memory for the lattices
    snew (lattice, bin[XX] * bin[YY] * bin[ZZ] );             //lattice
    snew (lattice_old, bin[XX] * bin[YY] * bin[ZZ] );         //lattice old for stdev calculation
    snew (lattice_squared, bin[XX] * bin[YY] * bin[ZZ] );     //squared for stdev calculation

    // === Definition of the cutoff and other gaussian-density related variables
    ncutoff[XX] = floor(cutoff/res[XX]);
    ncutoff[YY] = floor(cutoff/res[YY]);
    ncutoff[ZZ] = floor(cutoff/res[ZZ]);
    // Density will be increment in the grid points between:   (i-ncutoff, j-ncutoff, k-ncutoff) and (i+ncutoff, j+cutoff, k+ncutoff)


    // === change the cutoff to Angstrom^2 ===
    cutoff *= 10.0;


    // === NFRAMES = 0 ===
    Nframes = 0;
}


void
TimeAveragedDensity::analyzeFrame(int frnr, const t_trxframe &fr, gmx_unused t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata)
{


    // === STUFF TO READ THE TRAJECTORY AND ATOM SELECTION DATA ===
    AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    dh.startFrame(frnr, fr.time);
    size_t                     group = 0;   // selection-group index (only one considered thus it is 0)
    const Selection           &sel   = pdata->parallelSelection(sel_[group]);

    //------------------------------------------------------------------------------

    //=== STANDARD DEVIATION CALCULATION : Save the lattice value in lattice_old ===
    for (l = 0; l < bin[XX]*bin[YY]*bin[ZZ]; l++)
    {
        lattice_old[l] = lattice[l];
    }


    // === Count the visit of each one of the surrounding atoms (i.e. lipids) in the proper grid of the lattice

    // === Loop over the surrounding atoms (i.e. lipids )
    for (n = 0; n < nr; n++)
    {

        // Coordinate vector with respect to the simulation box
        SelectionPosition p = sel.position(n);   // read position selection of n-th atom
        R[XX] = p.x()[0];
        R[YY] = p.x()[1];
        R[ZZ] = p.x()[2];


        //		 printf("R at t=%8.3f : x= %8.5f y= %8.5f z= %8.5f\n",fr.time,  R[XX] , R[YY] , R[ZZ] );

        // Change of system of coordinates to one where the origin is the origin of the lattice
        // R_trans = coordinates of the atom with respect to the lattice
        // R_trans= R - L0
        // L0 origin of the lattice coordinates


        rvec_sub( R, L0, R_trans );      // R_trans = R - L0
        //	         printf("Rtrans at t=%8.3f : x= %8.5f y= %8.5f z= %8.5f\n",fr.time,  R_trans[XX] , R_trans[YY] , R_trans[ZZ] );


        // === Discretize the Rtrans to locate it in the lattice
        i = floor( R_trans[XX] / res[XX]);
        j = floor( R_trans[YY] / res[YY]);
        k = floor( R_trans[ZZ] / res[ZZ]);

        // === only consider the atoms located in the lattice ===
        // 0 <= i <= bin[XX]
        // 0 <= j <= bin[YY]
        // 0 <= k <= bin[ZZ]

        if (i >= 0 && i < bin[XX] && j >= 0 && j < bin[YY] && k >= 0 && k < bin[ZZ])
        {
            //printf("Rtrans at t=%8.3f : x= %8.5f y= %8.5f z= %8.5f\n",fr.time,  R_trans[XX] , R_trans[YY] , R_trans[ZZ] );



            //printf("i=%d j=%d k=%d\n",i,j,k);



            // Increment the density of the grid points to the atom, laying within a certain cutoff
            // loops over grid points around the (i,j,k) grid point
            for (nx = i-ncutoff[XX]; nx <= i+ncutoff[XX]; nx++)
            {
                for (ny = j-ncutoff[YY]; ny <= j+ncutoff[YY]; ny++)
                {
                    for (nz = k-ncutoff[ZZ]; nz <= k+ncutoff[ZZ]; nz++)
                    {

                        // === only consider grid points lying in the lattice ===
                        // 0 <= nx <= bin[XX]
                        // 0 <= ny <= bin[YY]
                        // 0 <= nz <= bin[ZZ]

                        if (nx >= 0 && nx < bin[XX] && ny >= 0 && ny < bin[YY] && nz >= 0 && nz < bin[ZZ])
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
                            dist_atom_grid = 10.0 * norm(  R_atom_grid );


                            // === Only contribute to the density if the distance is smaller than the cutoff
                            if (dist_atom_grid < cutoff)
                            {


                                // === Gaussian contribution (two or four gaussians depending on the atom) to the density at this grid point
                                density_contrib = C[n];


                                // contributions to the atomic potentials according to  equation (6) in J. Electron Microscopy, 56(4):131-140 (2007)
                                //gaussian contribution
                                for (g = 0; g < Ngaussians; g++)
                                {
                                    index_oneD_aux   = index_array(g, n, 0, Ngaussians, nr, 0); // turn (n,g) into a one-dimensional array.
                                    density_contrib += Ai[index_oneD_aux]*exp(-Bi[index_oneD_aux]*dist_atom_grid * dist_atom_grid );
                                }

                                // print crosscheck
                                //					if(density_contrib>=0)	 printf("%d %f %f\n",n,dist_atom_grid,density_contrib);




                                // printf("dist_atom-grid: (%f, %f, %f) norm=%f density_contrib=%f\n", R_atom_grid[XX],R_atom_grid[YY] , R_atom_grid[ZZ], dist_atom_grid ,density_contrib);

                                // === Increment the density at the position (nx,ny,nz)
                                index_aux           = index_array(nx, ny, nz, bin[XX], bin[YY], bin[ZZ]);
                                lattice[index_aux] += density_contrib;




                            } // Finish if statement that only distances smaller than the cut off are considered for the density
                        }     // finish if statement that grid points are within the lattice
                    }         // finish over z grid points
                }             // finish over y grid points
            }                 // finish over x grid points
        }                     // finish if statement considering atoms lying in the lattice
    }                         // finish loop over surrounding atoms



    // ===  STANDARD DEVIATION CALCULATION: store the increment in the density ===

    for (l = 0; l < bin[XX]*bin[YY]*bin[ZZ]; l++)
    {
        float increment = lattice[l]-lattice_old[l];

        // increment latice_squared for the standard deviation calculation
        lattice_squared[l] += increment*increment;

        //print crosscheck
        //   int index_aux=index_array(38,38,75,bin[XX],bin[YY],bin[ZZ]) ;
        //if (i==index_aux){   printf("test %f %d %f %f %f \n",fr.time,i,increment,lattice[i],lattice_squared[i]) ;}

    }

    // Increment the number of frames variable
    Nframes++;
    //------------------------------------------------------------------------------


    dh.finishFrame();
}


void
TimeAveragedDensity::finishAnalysis(int /*nframes*/)
{

    fprintf(stdout, "Nframes read: %f\n", Nframes);

    // === Size of the box at the last frame ===


    // === Calculate the density, normalizing the lattice values ===
    // density = Ncounts / (Nframes )
    for (l = 0; l < bin[XX]*bin[YY]*bin[ZZ]; l++)
    {
        // average
        lattice[l] /= (Nframes);  // (Angstrom^-1)
    }

    if (Nframes > 1)
    {
        for (l = 0; l < bin[XX]*bin[YY]*bin[ZZ]; l++)
        {
            //calculate standard deviation
            lattice_squared[l] = sqrt(  lattice_squared[l]/(Nframes-1.0) - lattice[l]*lattice[l]*Nframes/(Nframes-1.0) );

            // print crosscheck
            //    if (i==1401321)   printf("%f %f nframes=%f\n",lattice[i],lattice_squared[i],Nframes) ;
        }
    }
}


void
TimeAveragedDensity::writeOutput()
{
    // output the average to file
    if (!map_av_.empty())
    {
        print_to3dmap(map_av_.c_str(), N, Nmin, Nmax, BOX, bin, res, lattice);
    }

    // output the standard deviation to file
    if (!map_sd_.empty() && (Nframes > 1) )
    {
        print_to3dmap(map_sd_.c_str(), N, Nmin, Nmax, BOX, bin, res, lattice_squared);
    }

    // output x, y, z, rho, and sigma to file
    if (!rho_.empty())
    {
        print_rho_sigma(rho_.c_str(), bin, res,  Nmin, lattice, lattice_squared );
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
