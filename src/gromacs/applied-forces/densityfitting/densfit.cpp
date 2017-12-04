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


#include "densfit.h"
#include "../maputil.h"

#include <cstdio>
#include <vector>
#include <cassert>

#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/topology/block.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/fileio/xvgr.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/topology/topology.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/pleasecite.h"

#include "gromacs/mdlib/sim_util.h"
#include "gromacs/commandline/filenm.h"

#include "gromacs/mdlib/groupcoord.h"

#include "gromacs/domdec/domdec_struct.h"

#include "gromacs/timing/cyclecounter.h"

#include "gromacs/domdec/domdec.h"

#include "gromacs/fileio/griddataio.h"

#include "gromacs/pbcutil/pbc.h"

static const char * FitStrMDRUN = "Density fitting: ";  //!< Used if density fitting output is written from mdrun
static const char * FitStr      = "";                   //!< Used in gmx map output
//! \brief Environment variable for setting nstfit
static const char*  NSTFIT_ENVVAR  =  "GMX_NSTFIT";

namespace gmx
{

namespace
{
//! \brief Sum grid contributions from all ranks
void sum_grid(t_mapdata *map, t_commrec *cr)
{
    int nx, ny, nz;


    nx = map->map_dim[XX];
    ny = map->map_dim[YY];
    nz = map->map_dim[ZZ];

    gmx_sumf(nx*ny*nz, map->vox.data(), cr);
}


//! \brief Adds 'buf' to 'str'
static void add_to_string(char **str, char *buf)
{
    int len;


    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}

//! \brief TODO
static void add_to_string_aligned(char **str, const char *buf)
{
    char buf_aligned[STRLEN];

    sprintf(buf_aligned, "%12s", buf);

    add_to_string(str, buf_aligned);
}


//! \brief Exponential function exp for a vector of n elements
void gmx_exp_vector(int n, real x[])
{
//    gmx_cycles_t c1 = gmx_cycles_read();

#if GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal sx = load<SimdReal>(x+i);
        sx = gmx::exp(sx);
        store(x+i, sx);
    }
#else       // GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i++)
    {
        x[i] = std::exp(x[i]);
    }
#endif      // GMX_SIMD_HAVE_REAL
//    gmx_cycles_t c2 = gmx_cycles_read();
//    fprintf(stderr, "exp vector %15d elements %15lld cycles  %15.3f per exp\n", n, c2-c1, (1.0*(c2-c1))/n);
}

/*! \brief Calculate one of the terms of the correlation coefficient
 *
 * Multiply each voxel of the map with the corresponding voxel of the other
 * map and sum up everything, i.e. compute
 *
 *  sum_ijk [ rhoA_ijk * rhoB_ijk ]
 *
 */
static double calc_sum_rhoA_rhoB(
        const std::vector<float> &vecA,
        const std::vector<float> &vecB)
{
    double sum_ijk = 0.;
    for (size_t i = 0; i < vecA.size(); i++)
    {
        sum_ijk += vecA[i] * vecB[i];
    }
    return sum_ijk;
}

//! \brief For a position vector x return all direct 27 image positions given the box
int get_images(rvec x, matrix box,
               rvec img[], ivec img_startindex[], ivec img_stop_index[],
               real grid_spacing, int voxrange, const std::array<int, DIM> &map_dim)
{
    int    i, j, k, m, count, icenter;
    rvec   image, shiftvecX, shiftvecY, shiftvecZ;
    ivec   istart, istop;

    double OOgrid_spacing = 1.0 / grid_spacing;

    /* Loop over all possible 3x3x3 neighboring boxes */
    /* Minimum image convention? TODO: loop over NTRICIMG = 14, not 27 */
    count = 0;
    for (i = -1; i <= 1; i++)
    {
        for (j = -1; j <= 1; j++)
        {
            for (k = -1; k <= 1; k++)
            {
                svmul(k, box[XX], shiftvecX);
                svmul(j, box[YY], shiftvecY);
                svmul(i, box[ZZ], shiftvecZ);

                copy_rvec(x, image);

                /* Move in x, y, and z directions */
                rvec_inc(image, shiftvecX);
                rvec_inc(image, shiftvecY);
                rvec_inc(image, shiftvecZ);

                /* Now check whether we will get any density from that image */
                icenter     = roundf(image[XX] * OOgrid_spacing);
                istart [XX] = icenter - voxrange;
                istop  [XX] = icenter + voxrange;

                if ( (istart[XX] <= map_dim[XX]) && (istop[XX] >= 0) )
                {
                    icenter     = roundf(image[YY] * OOgrid_spacing);
                    istart [YY] = icenter - voxrange;
                    istop  [YY] = icenter + voxrange;

                    if ( (istart[YY] <= map_dim[YY]) && (istop[YY] >= 0) )
                    {

                        icenter     = roundf(image[ZZ] * OOgrid_spacing);
                        istart [ZZ] = icenter - voxrange;
                        istop  [ZZ] = icenter + voxrange;

                        if ( (istart[ZZ] <= map_dim[ZZ]) && (istop[ZZ] >= 0) )
                        {
                            /* This is our new image */
                            copy_rvec(image, img[count]);

                            /* Cut away the volume where we do not have a grid */
                            for (m = 0; m < DIM; m++)
                            {
                                istart[m] = std::max(istart[m], 0         ); /* istart >= 0       */
                                istop [m] = std::min(istop [m], map_dim[m]); /* istop  <= map_dim */
                            }
                            copy_ivec(istart, img_startindex[count]);
                            copy_ivec(istop, img_stop_index[count]);

                            count++;

                        } /* ZZ */
                    }     /* YY */
                }         /* XX */
            }
        }
    }

    return count;
}

//! \brief The error function erf for a vector of n elements
void gmx_erf_vector(int n, real x[])
{
//    gmx_cycles_t c1 = gmx_cycles_read();

#if GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal sx = load<SimdReal>(x+i);
        sx = gmx::erf(sx);
        store(x+i, sx);
    }
#else       // GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i++)
    {
        x[i] = std::erf(x[i]);
    }
#endif      // GMX_SIMD_HAVE_REAL
//    gmx_cycles_t c2 = gmx_cycles_read();
//    fprintf(stderr, "erf vector %15d elements %15lld cycles  %15.3f per erf (SIMD width is %d)\n", n, c2-c1, (1.0*(c2-c1))/n, GMX_SIMD_REAL_WIDTH);
}
}

//! \brief maputil.cpp local density fitting data
struct gmxDensfit
{
    gmxDensfit(const Densfit &densfit,
               rvec          *x_densfit_whole, /* Can be NULL!                       */
               real           grid_spacing,
               gmx_bool       bVerbose,
               gmx_bool       bAppend,
               gmx_bool       bParallel);
    ~gmxDensfit() = default;
    real k;                     /**< Current value k of strength constant of map
                                     potential V_fit = k*(1 - cc) with cc being the
                                     correlation coefficient                        */
    real         sigma;         /**< Current value of Gaussian width used for
                                     spreading atomic positions on the map          */
    real         temp;          /**< Current value of temperature                   */
    FILE        *out_dens;      /**< Output file for density fitting data           */
    const char  *fn_map;        /**< Output filename for simulated density maps     */
    real         spacing;       /**< Map spacing in nm, derived from map data       */
    gmx_bool     bAppend;       /**< Are we appending to existing output files?     */
    gmx_bool     bVerbose;      /**< -v flag from command line                      */
    real         cc;            /**< Correlation coefficient of the two maps        */
    rvec        *x_assembled;   /**< Just the positions that are spread             */
    rvec        *f_loc;         /**< Forces on the local atoms of the
                                     density fitting group                          */
    real        *w_assembled;   /**< Weights for the assembled positions            */
    rvec        *x_old;         /**< x_assembled from last step                     */
    ivec        *x_shifts;      /**< To make the molecule whole                     */
    ivec        *extra_shifts;  /**< Shifts added since last NS                     */
    int         *c_ind;         /**< Where is this atom stored in the coll array?   */
    gmx_bool     bUpdateShifts; /**< Do we have to calculate new shift vectors
                                     from x_assembled and x_old?                    */
    int          voxrange;      /**< Max. number of voxels to be computed for a
                                     single atom in a single dimension x, y, or z   */
    int          vox_per_atom;  /**< Total number of voxels for one atom (x,y, + z) */
    double       sum_rho_exp2;  /**< Term needed for the calculation of the
                                     correlation coefficient and the forces         */

    /* The following two temporary vectors (one for each OpenMP thread) store
     * erf values around a single atoms, thus we can compute them all in one go,
     * and with SIMD acceleration */
    std::vector<real, gmx::AlignedAllocator<real> > *erfVector;  /**< pointer to vector of erf values */
    std::vector<real, gmx::AlignedAllocator<real> > *expVector;  /**< same for exp values             */
    private:
        /* TODO: crosscheck with code in new_t_gmx_densfit(), which depends on sigma!! */
        //! \brief TODO
        void initParametersThatCanBeTimeDependent(const Densfit &densfit);
};


//! \brief Allocate and initialize the arrays for assembled positions, forces and atomic weights
gmxDensfit::gmxDensfit(
        const Densfit &densfit,
        rvec          *x_densfit_whole, /* Can be NULL!                       */
        real           grid_spacing,
        gmx_bool       bVerbose,
        gmx_bool       bAppend,
        gmx_bool       bParallel)
{

    this->bVerbose = bVerbose;
    this->spacing  = grid_spacing;

    snew(this->x_assembled, densfit.nat);
    snew(this->f_loc, densfit.nat);
    snew(this->w_assembled, densfit.nat);
    snew(this->x_old, densfit.nat);
    snew(this->x_shifts, densfit.nat);
    snew(this->extra_shifts, densfit.nat);
    snew(this->c_ind, densfit.nat);

    if (this->bVerbose)
    {
        fprintf(stdout, "%sSetting all %d atomic weights to unity.\n", FitStr, densfit.nat);
    }
    for (int i = 0; i < densfit.nat; i++)
    {
        this->w_assembled[i] = 1.0;
    }

    this->bUpdateShifts = TRUE;
    this->bAppend       = bAppend;

    if (!bParallel)
    {
        for (int i = 0; i < densfit.nat; i++)
        {
            this->c_ind[i] = i;
        }
    }

    /* In mdruns we have to make the density fitting structure whole */
    if (NULL != x_densfit_whole)
    {
        for (int i = 0; i < densfit.nat; i++)
        {
            copy_rvec(x_densfit_whole[i], this->x_old[i]);
        }
    }

    /* Initialize the possibly time-dependent parameters */
    initParametersThatCanBeTimeDependent(densfit);
}

void gmxDensfit::initParametersThatCanBeTimeDependent(const Densfit &densfit)
{
    int           ompthreads;

    this->sigma = densfit.sigma_values[0];
    this->k     = densfit.k_values[0];
    this->temp  = densfit.temp_values[0];

/* For erf() and exp() performance, we want all values for one atom in one
 * linear array (not in three). Therefore, we compute the maximum number of
 * voxels to be computed for a single atom. How many voxels around each atom
 * we have to spread the density on is saved in voxrange (for each dimension
 * x, y, and z */
    this->voxrange     = ceil(0.5 * this->sigma * densfit.dist / this->spacing);
    this->vox_per_atom = 3*(2*this->voxrange + 1);

    ompthreads         = std::max(1, gmx_omp_nthreads_get(emntDefault));
    this->erfVector    = new std::vector<real, gmx::AlignedAllocator<real> >[ompthreads];
    this->expVector    = new std::vector<real, gmx::AlignedAllocator<real> >[ompthreads];
    for (int i = 0; i < ompthreads; i++)
    {
        int extra = 0;
#if GMX_SIMD_HAVE_REAL
        extra = GMX_SIMD_REAL_WIDTH;
#endif
        this->erfVector[i].reserve(this->vox_per_atom + extra);
        this->expVector[i].reserve(this->vox_per_atom + extra);
        // + GMX_SIMD_REAL_WIDTH adds some extra elements at the end of the array so that
        // we can always read GMX_SIMD_REAL_WIDTH elements, even near the end of the vector.
    }
}

Densfit::Densfit(){}

Densfit::Densfit(
        real      sigma,
        real      sigma_dist,
        real      k,        /* Spring constant */
        int       nat,
        int      *ind,      /* Which of the atoms should be used for spreading */
        real      grid_spacing,
        bool      bVerbose)
{
    int        i;
    this->npoints = 1;
    snew(this->time_values, this->npoints);
    snew(this->temp_values, this->npoints);
    snew(this->k_values, this->npoints);
    snew(this->sigma_values, this->npoints);
    this->sigma_values[0] = sigma;
    this->k_values[0]     = k;
    this->dist            = sigma_dist;
    this->nat             = nat;

    if (NULL == ind)
    {
        if (bVerbose)
        {
            fprintf(stdout, "%sYou did not provide an index group - will use all atoms for spreading.\n", FitStr);
        }

        snew(this->ind, this->nat);
        for (i = 0; i < this->nat; i++)
        {
            this->ind[i] = i;
        }
    }
    else
    {
        this->ind = ind;
    }
    setPrivateData(nullptr, grid_spacing, bVerbose, 0, 0);
};

Densfit::~Densfit(){};

void Densfit::makeGroups(char *groupname, t_blocka *grps, char **gnames)
{
    /* Group index from group name */
    int ig           = search_string(groupname, grps->nr, gnames);
    this->nat = grps->index[ig + 1] - grps->index[ig];

    fprintf(stderr, "Density fitting group '%s' has %d atoms\n", groupname,
            this->nat);
    snew(this->ind, this->nat);
    for (int i = 0; i < this->nat; i++)
    {
        this->ind[i] = grps->a[grps->index[ig] + i];
    }
}


//! \brief TODO
real Densfit::interpolateLinearly(real time, int npoints, real* time_values, real* y_values)
{
    real x, y;
    int  j;


    /* If we are left of any time point, return the zero'th y value */
    if (time <= time_values[0])
    {
        return y_values[0];
    }

    /* If we are to the right of any time point, return the last y value */
    if (time >= time_values[npoints-1])
    {
        return y_values[npoints-1];
    }

    /* We are somewhere in between, that means we need to interpolate y linearly */
    j = 0;
    while (time > time_values[j+1])
    {
        j++;
        assert(j < npoints);
    }

    /* Found our position between points j and j+1.
     * Interpolate: x is the amount from j+1, (1-x) from point j.
     */
    x = ((time - time_values[j]) /
         (time_values[j+1] - time_values[j]));
    y = x * y_values[j+1] + (1-x) * y_values[j];

    return y;
}
void Densfit::updateParametersThatAreTimeDependent(real time)
{
    int           ompthreads;

    df->temp  = interpolateLinearly(time, this->npoints, this->time_values, this->temp_values );
    assert(df->temp >= -1);
    df->sigma = interpolateLinearly(time, this->npoints, this->time_values, this->sigma_values);
    assert(df->sigma > 0.0);
    df->k     = interpolateLinearly(time, this->npoints, this->time_values, this->k_values    );

    /* For erf() and exp() performance, we want all values for one atom in one
     * linear array (not in three). Therefore, we compute the maximum number of
     * voxels to be computed for a single atom. How many voxels around each atom
     * we have to spread the density on is saved in voxrange (for each dimension
     * x, y, and z */
    df->voxrange     = ceil(0.5 * df->sigma * this->dist / df->spacing);
    df->vox_per_atom = 3*(2*df->voxrange + 1);

    ompthreads       = std::max(1, gmx_omp_nthreads_get(emntDefault));
    // TODO: initialize with the maximum vox_per_atom
    for (int i = 0; i < ompthreads; i++)
    {
        int extra = 0;
    #if GMX_SIMD_HAVE_REAL
        extra = GMX_SIMD_REAL_WIDTH;
    #endif
        // Give each OpenMP thread its own local storage; erfVector[i] is for OpenMP thread number i
        df->erfVector[i].reserve(df->vox_per_atom + extra);
        df->expVector[i].reserve(df->vox_per_atom + extra);
    }
}

void Densfit::dump_x(int nodeid, gmx_int64_t step)
{
    FILE          *fpout = nullptr;
    char           fn[STRLEN];
    char           buf[256];

    gmx_step_str(step, buf);
    sprintf(fn, "x_step%s_node%d.txt", buf, nodeid);
    fpout = fopen(fn, "w");

    for (int i = 0; i < this->nat; i++)
    {
        fprintf(fpout, "%4d %12.3e %12.3e %12.3e\n", i,
                df->x_assembled[i][XX],
                df->x_assembled[i][YY],
                df->x_assembled[i][ZZ]);
    }

    fclose(fpout);
}

void Densfit::dump_f(const char *fn, t_commrec * cr)
{
    int            i;
    rvec          *f;
    gmx_bool       bWrite;

    FILE          *fpout = NULL;



    /* Assemble the forces */
    snew(f, this->nat);
    /* Each node copies its local forces into the collective array */
    for (i = 0; i < this->nat_loc; i++)
    {
        copy_rvec(df->f_loc[i], f[df->c_ind[i]]);
    }

    /* Reduce */
    if (cr && PAR(cr))
    {
        gmx_sum(3*this->nat, f[0], cr);
    }

    /* Master outputs */
    bWrite = FALSE;
    if (cr)
    {
        if (MASTER(cr))
        {
            bWrite = TRUE;
        }
    }
    else
    {
        bWrite = TRUE;
    }

    if (bWrite)
    {
        if (fn)
        {
            fpout = gmx_fio_fopen(fn, "w");
        }
        else
        {
            fpout = df->out_dens;
        }

        fprintf(stderr, "%sDumping forces to file %s.\n", FitStr, fn ? fn : "");

        for (i = 0; i < this->nat; i++)
        {
            fprintf(fpout, "#     f[%5d]={%12.5e, %12.5e, %12.5e}\n", i,
                    f[i][XX],
                    f[i][YY],
                    f[i][ZZ]);
        }

        if (fn)
        {
            gmx_fio_fclose(fpout);
        }
    }
    sfree(f);
}


void Densfit::setPrivateData(rvec *x_densfit_whole, real grid_spacing, bool bVerbose, bool bAppend, bool bParallel)
{
    /* Initialize also the private data structure: */
    this->df = std::unique_ptr<gmxDensfit>( new gmxDensfit(*this, x_densfit_whole, grid_spacing, bVerbose, bAppend, bParallel));
}

void Densfit::couple_map_spacing_to_box(matrix box)
{
    int        i;
    rvec       spacing;

    for (i = 0; i < DIM; i++)
    {
        map_ref.cell[i] = NM2A * box[i][i];
        spacing[i]      = A2NM * map_ref.cell[i] / map_ref.map_dim[i];
    }

    /* Save the average grid spacing */
    df->spacing = (spacing[XX]+spacing[YY]+spacing[ZZ]) / 3.0;
}

void Densfit::spread_atoms_low( rvec x[], int        natoms, matrix     box)
{

    int                iz, iy, ix, i, j, n, th;
    int                nx, ny, nz, ngrid;
    double             V_vox;
    double             prefactor;
    double             left;
    real               ax, ay, az, ayaz;
    float             *ptr;             /* Pointer in that array                          */
    ivec               elements;
    rvec               pos;
    real              *aarr[3];
    int                image, num_images;
    rvec               img[27];        /* Stores PBC images of a position x as well as x */
    ivec               istart[27];     /* Store for each image the part of the grid ...  */
    ivec               istop[27];      /* ... where stuff has to be calculated           */

    double             sigma2     = df->sigma*df->sigma;
    double             c_fac      = sqrt(1.5/sigma2);
    int                nth        = std::max(1, gmx_omp_nthreads_get(emntDefault));

    //    gmx_cycles_t c1 = gmx_cycles_read();

    auto &sim_dens_1d = this->map_sim.vox;                  /* Calculated density map */
    V_vox       = pow(df->spacing, 3);
    prefactor   = pow(M_PI*sigma2/6.0, 1.5);
    prefactor  /= V_vox;
    nx          = this->map_sim.map_dim[XX];
    ny          = this->map_sim.map_dim[YY];
    nz          = this->map_sim.map_dim[ZZ];
    ngrid       = nx*ny*nz;



    /* Zero out the grid */
    #pragma omp parallel for num_threads(nth) schedule(static) \
    shared(sim_dens_1d, ngrid, nth) \
    default(none)
    for (i = 0; i < ngrid*nth; i++)
    {
        sim_dens_1d[i] = 0.0;
    }

    const auto dfp        = this->df.get();
    const auto mapsim_dim = map_sim.map_dim;
    /* Spread the local atoms on the grid, therefore loop over all atoms of the
     * density fitting group */
    #pragma omp parallel for num_threads(nth) schedule(static) \
    shared(nx, ny, nz, ngrid, natoms, x, box, dfp, mapsim_dim, c_fac, sim_dens_1d) \
    private(ptr, aarr, istart, istop, image, pos, num_images, img, i, ix, iy, iz, th, elements, ax, ay, az, ayaz, left) \
    default(none)
    for (n = 0; n < natoms; n++)
    {
        th = gmx_omp_get_thread_num();     /* Thread ID */

        /* Determine the grid points around the atom and its images for which we need to calculate stuff */
        num_images = get_images(x[n], box, img, istart, istop,
                                dfp->spacing, dfp->voxrange, mapsim_dim);

        /* Loop over all periodic images of the atom */
        for (image = 0; image < num_images; image++)
        {
            copy_rvec(img[image], pos);     // periodic image of the position

            elements[XX] = istop[image][XX] + 1 - istart[image][XX];
            elements[YY] = istop[image][YY] + 1 - istart[image][YY];
            elements[ZZ] = istop[image][ZZ] + 1 - istart[image][ZZ];

            //            int sum = elements[XX] + elements[YY] + elements[ZZ];

            /* For erf() performance, we want all values for one atom in one linear array (not
             * int three). However, we build pointers to access it as x, y, and z arrays.
             *
             * Each OpenMP thread th accesses its own temporary storage vector
             */
            aarr[XX] = &dfp->erfVector[th][0                              ];
            aarr[YY] = &dfp->erfVector[th][0 + elements[XX]               ];
            aarr[ZZ] = &dfp->erfVector[th][0 + elements[XX] + elements[YY]];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in a temporary array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            i = 0;
            for (ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                left          = (ix-0.5)*dfp->spacing - pos[XX];
                aarr[XX][i++] = c_fac*left;
            }

            i = 0;
            for (iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                left          = (iy-0.5)*dfp->spacing - pos[YY];
                aarr[YY][i++] = c_fac*left;
            }

            i = 0;
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                left          = (iz-0.5)*dfp->spacing - pos[ZZ];
                aarr[ZZ][i++] = c_fac*left;
            }
            gmx_erf_vector(elements[XX]+elements[YY]+elements[ZZ], aarr[XX]);

            /* Transform (in place) into a new array that directly stores the differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix+1] - aarr[XX][ix];
            }
            for (iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy+1] - aarr[YY][iy];
            }
            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz+1] - aarr[ZZ][iz];
            }

            /* Spread the density for this atom */
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ]; iz++)
            {
                az = aarr[ZZ][iz-istart[image][ZZ]];

                for (iy = istart[image][YY]; iy < istop[image][YY]; iy++)
                {
                    ay    = aarr[YY][iy-istart[image][YY]];
                    ayaz  = ay*az;
                    ptr   = &sim_dens_1d[th*ngrid +  nx*(ny*iz + iy) + istart[image][XX] ];
                    for (ix = istart[image][XX]; ix < istop[image][XX]; ix++)
                    {
                        ax = aarr[XX][ix-istart[image][XX]];

                        /* Assign value to voxel */
                        //                        sim_dens_1d[nz*ny*ix + nz*iy + iz] += ax*ay*az;
                        *ptr += ax*ayaz;
                        ptr++;
                    }
                }
            }
        }     /* End of loop over images of one atom */
    }         /* End of loop over atoms */

    #pragma omp parallel for num_threads(nth) schedule(static) \
    shared(nth, sim_dens_1d, ngrid, prefactor) \
    private(j) \
    default(none)
    for (i = 0; i < ngrid; i++)
    {
        for (j = 1; j < nth; j++)
        {
            sim_dens_1d[i] += sim_dens_1d[i + j*ngrid];
        }
        sim_dens_1d[i] *= prefactor;
    }

}

void Densfit::spread_atoms(matrix box)
{
//    double         t0 = 0;
//    long int       c1 = 0; /* make compiler happy */

    int           *map_dim = this->map_sim.map_dim.data();
    if (df->bVerbose)
    {
        fprintf(stderr, "%sSpreading %d atomic positions  on a %d x %d x %d grid ...\n",
                FitStr, this->nat, map_dim[XX], map_dim[YY], map_dim[ZZ]);
//        t0 = gettime();
//        c1 = gmx_cycles_read();
    }

    if (df->spacing <= 0)
    {
        gmx_fatal(FARGS, "%sSpacing must not be zero!", FitStr);
    }

    /* Make the density map */
    spread_atoms_low(df->x_assembled, this->nat, box);

//    if (df->bVerbose)
//    {
//        fprintf(stderr, "%sSpreading took %g seconds (%g M cycles).\n", FitStr,
//                gettime() - t0, 0.000001*(gmx_cycles_read()-c1));
//    }
}

void Densfit::do_forces(matrix box)
{
    int                 ix, iy, iz, i, l, ind;
    int                 nx, ny;
    double              sum_rho_exp2;
    double              sum_rho_sim2;
    double              sum_rho_exp_rho_sim;
    double              term1_prefactor, term2_prefactor;
    dvec                term1_sum, term2_sum;
    dvec                tmp1, tmp2;
    dvec                drho_dxl;
    dvec                force;
    rvec                pos;
    int                 th;    /* ID of the thread in this team */
    ivec                elements;
    double              V_vox; //TODO: with pressure coupling, the grid spacing needs to be adjusted!!!
    double              prefactor;
    double              left;
    double              OOsigma2_fac;     /* 3 / (2*sigma*sigma) */
    double              sqrtOOsigma2_fac; /* sqrt of the above */
    real                erfx, erfy, erfz; /* erf( sqrt( 3 / (2*sigma*sigma) * x)  for x */
    real                expx, expy, expz; /* exp( -3 x*x / (2*sigma*sigma) */
    real                erfy_erfz, expy_erfz, erfy_expz;

    real               *aarr[3], *garr[3];
    std::vector<float>  sim_dens_1d, ref_dens_1d;
    float              *ptr_ref, *ptr_sim;
    double              ref_dens, sim_dens;

    int                 image, num_images;
    rvec                img[27];         /* Stores PBC images of a position x as well as x */
    ivec                istart[27];      /* Store for each image the part of the grid ...  */
    ivec                istop[27];       /* ... where stuff has to be calculated           */

    sim_dens_1d     = this->map_sim.vox; /* Calculated density map */
    ref_dens_1d     = this->map_ref.vox; /* Reference density map */


    if (this->nat_loc <= 0)
    {
        return;
    }

    // cppcheck-suppress unreadVariable
    int gmx_unused nth = std::max(1, gmx_omp_nthreads_get(emntDefault)); // number of threads

    /* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (l = 0; l < this->nat_loc; l++)
    {
        clear_rvec(df->f_loc[l]);
    }

    /* Calculate various prefactors */
    sum_rho_exp2        = this->df->sum_rho_exp2;
    sum_rho_sim2        = calc_sum_rhoA_rhoB(this->map_sim.vox, this->map_sim.vox);
    sum_rho_exp_rho_sim = calc_sum_rhoA_rhoB(this->map_ref.vox, this->map_sim.vox);

    term1_prefactor     =   df->k / sqrt(sum_rho_exp2 * sum_rho_sim2);
    term2_prefactor     = -df->k * sum_rho_exp_rho_sim / (sqrt(sum_rho_exp2)*pow(sum_rho_sim2, 1.5));

    V_vox               = df->spacing*df->spacing*df->spacing;
    prefactor           = -M_PI*df->sigma*df->sigma/(6.0*V_vox);
    OOsigma2_fac        = 1.5 / (df->sigma*df->sigma);
    sqrtOOsigma2_fac    = sqrtf(OOsigma2_fac);

    nx  = this->map_sim.map_dim[XX];
    ny  = this->map_sim.map_dim[YY];

    auto     * dfp        = df.get();
    const auto mapsim_dim = this->map_sim.map_dim;
#pragma omp parallel for num_threads(nth) schedule(static) \
    private(term1_sum, term2_sum, pos, istart, istop, elements, aarr, garr, i, \
    ix, iy, iz, left, ind, erfx, erfy, erfz, expx, expy, expz, erfy_erfz, expy_erfz, erfy_expz, \
    ptr_ref, ptr_sim, ref_dens, sim_dens, drho_dxl, force, tmp1, tmp2, th, num_images, image, img) \
    shared(dfp, ref_dens_1d, sim_dens_1d, box, nx, ny, \
    OOsigma2_fac, sqrtOOsigma2_fac, prefactor, term1_prefactor, term2_prefactor) \
    default(none)
    /* Loop over all atoms of the density fitting group */
    for (l = 0; l < this->nat_loc; l++)
    {
        th = gmx_omp_get_thread_num();

        /* Loop over all voxels and calculate term1 and term2 for the force
         * simultaneously */
        clear_dvec(term1_sum);
        clear_dvec(term2_sum);

        /* Determine the grid points around the atom and its images which
         * have a contribution to the force on this atom.
         * Get the position of this atom of the density fitting group, we must
         * use the positions from the assembled array, because only these have
         * been put in the box in do_this() */
        num_images = get_images(dfp->x_assembled[dfp->c_ind[l]], box, img, istart, istop,
                                dfp->spacing, dfp->voxrange, mapsim_dim);

        /* Loop over all periodic images of the atom */
        for (image = 0; image < num_images; image++)
        {
            copy_rvec(img[image], pos); // periodic image of the position

            /* The calculation of drho_dxl is expensive. We only calculate it, where
             * there is actually simulated density, i.e. in the given sigma-range applied
             * during spreading. */
            elements[XX] = istop[image][XX] + 1 - istart[image][XX];
            elements[YY] = istop[image][YY] + 1 - istart[image][YY];
            elements[ZZ] = istop[image][ZZ] + 1 - istart[image][ZZ];

            /* For erf() performance, we want all values for one atom in one linear array (not
             * int three). However, we build pointers to access it as x, y, and z arrays.
             *
             * Each OpenMP thread th accesses its own temporary storage vector
             */
            aarr[XX] = &dfp->erfVector[th][0                              ];
            aarr[YY] = &dfp->erfVector[th][0 + elements[XX]               ];
            aarr[ZZ] = &dfp->erfVector[th][0 + elements[XX] + elements[YY]];

            /* Same for the temporary array that keeps the exp values */
            garr[XX] = &dfp->expVector[th][0                              ];
            garr[YY] = &dfp->expVector[th][0 + elements[XX]               ];
            garr[ZZ] = &dfp->expVector[th][0 + elements[XX] + elements[YY]];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in the tmpgrid[] array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            i = 0;
            for (ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                left          = (ix-0.5)*dfp->spacing - pos[XX];
                garr[XX][i  ] = -OOsigma2_fac*left*left;
                aarr[XX][i++] = sqrtOOsigma2_fac*left;
            }

            i = 0;
            for (iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                left          = (iy-0.5)*dfp->spacing - pos[YY];
                garr[YY][i  ] = -OOsigma2_fac*left*left;
                aarr[YY][i++] = sqrtOOsigma2_fac*left;
            }

            i = 0;
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                left          = (iz-0.5)*dfp->spacing - pos[ZZ];
                garr[ZZ][i  ] = -OOsigma2_fac*left*left;
                aarr[ZZ][i++] = sqrtOOsigma2_fac*left;
            }

            /* Call erf() and exp() for all input values for this atom in one go */
            gmx_erf_vector(elements[XX]+elements[YY]+elements[ZZ], aarr[XX]);
            gmx_exp_vector(elements[XX]+elements[YY]+elements[ZZ], garr[XX]);

            /* Transform (in place) into a new array that directly stores the differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix+1] - aarr[XX][ix];
                garr[XX][ix] = garr[XX][ix+1] - garr[XX][ix];
            }
            for (iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy+1] - aarr[YY][iy];
                garr[YY][iy] = garr[YY][iy+1] - garr[YY][iy];
            }
            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz+1] - aarr[ZZ][iz];
                garr[ZZ][iz] = garr[ZZ][iz+1] - garr[ZZ][iz];
            }

            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                erfz = aarr[ZZ][iz];
                expz = garr[ZZ][iz];

                for (iy = 0; iy < elements[YY] - 1; iy++)
                {
                    erfy = aarr[YY][iy];
                    expy = garr[YY][iy];

                    erfy_erfz = erfy * erfz;
                    expy_erfz = expy * erfz;
                    erfy_expz = erfy * expz;

                    ind     = ny*nx*(iz+istart[image][ZZ]) + nx*(iy+istart[image][YY]) + istart[image][XX];
                    ptr_ref = &ref_dens_1d[ind];
                    ptr_sim = &sim_dens_1d[ind];

                    /* Calculate d/dx_l rho_sim_ijk */
                    for (ix = 0; ix < elements[XX] - 1; ix++)
                    {
                        ref_dens = (double) *ptr_ref;
                        sim_dens = (double) *ptr_sim;

                        erfx = aarr[XX][ix];
                        expx = garr[XX][ix];

                        drho_dxl[XX] = expx * erfy_erfz;
                        drho_dxl[YY] = erfx * expy_erfz;
                        drho_dxl[ZZ] = erfx * erfy_expz;

                        term1_sum[XX] += ref_dens * drho_dxl[XX];
                        term1_sum[YY] += ref_dens * drho_dxl[YY];
                        term1_sum[ZZ] += ref_dens * drho_dxl[ZZ];

                        term2_sum[XX] += sim_dens * drho_dxl[XX];
                        term2_sum[YY] += sim_dens * drho_dxl[YY];
                        term2_sum[ZZ] += sim_dens * drho_dxl[ZZ];

                        ptr_ref++;
                        ptr_sim++;
                    }
                }
            } /* end of loop over voxels of this image */
        }     /* end of loop over images */

        dsvmul(prefactor*term1_prefactor, term1_sum, tmp1);
        dsvmul(prefactor*term2_prefactor, term2_sum, tmp2);

        dvec_add(tmp1, tmp2, force);

        dfp->f_loc[l][XX] = force[XX];
        dfp->f_loc[l][YY] = force[YY];
        dfp->f_loc[l][ZZ] = force[ZZ];

    }   /* end of loop over atoms */

}

real Densfit::calc_correlation_coeff(FILE *log)
{
    int              i;
    double           numerator, denominator;
    real             cc;

    if (log)
    {
        fprintf(log, "%sCalculating the correlation coefficient of the two maps.\n", FitStr);
    }

    for (i = 0; i < 3; i++)
    {
        if (this->map_ref.map_dim[i] != this->map_sim.map_dim[i])
        {
            gmx_fatal(FARGS, "Map dimensions must agree");
        }
    }

    numerator = calc_sum_rhoA_rhoB(this->map_ref.vox, this->map_sim.vox);

    denominator = sqrt(calc_sum_rhoA_rhoB(this->map_ref.vox, this->map_ref.vox) * calc_sum_rhoA_rhoB(this->map_sim.vox, this->map_sim.vox));

    cc = numerator / denominator;

    if (log)
    {
        fprintf(log, "%sThe correlation coefficient is %15.10f\n", FitStr, cc);
    }

    return cc;
}

/*! \brief Open output file and print some general information about density fitting
 *
 * Call on master only
 */
FILE * Densfit::open_out(const char *fn, const gmx_output_env_t *oenv)
{
    FILE        *fp;
    int          nsets;
    const char **setname;
    char         buf[50];
    char        *LegendStr = NULL;

    if (df->bAppend)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = xvgropen(fn, "Density fitting correlation coefficient", "Time (ps)", "correlation coefficient", oenv);
        fprintf(fp, "# Density fitting forces for %d atoms are updated every %d time step%s.\n", this->nat, this->nstfit, this->nstfit != 1 ? "s" : "");
        fprintf(fp, "#   Note that the force update frequency can be overwritten by the environment variable %s.\n", NSTFIT_ENVVAR);
        fprintf(fp, "# Output is written in intervals of %d time step%s.\n", this->nstout, this->nstout > 1 ? "s" : "");
        fprintf(fp, "# Number of refinement steps is %d.\n", this->npoints);
        fprintf(fp, "# Parameter sigma_dist = %g\n", this->dist);
        if (0 != this->nstmapout)
        {
            fprintf(fp, "# Calculated map will be saved every %d steps\n", this->nstmapout);
        }

        /* Print a nice legend */
        snew(LegendStr, 1);
        LegendStr[0] = '\0';
        sprintf(buf, "#     %6s", "time");
        add_to_string_aligned(&LegendStr, buf);

        nsets = 0;
        snew(setname, 5);
        setname[nsets++] = strdup("temperature (K)");
        add_to_string_aligned(&LegendStr, "T");
        setname[nsets++] = strdup("k (kJ/mol)");
        add_to_string_aligned(&LegendStr, "k");
        setname[nsets++] = strdup("sigma (nm)");
        add_to_string_aligned(&LegendStr, "sigma");
        setname[nsets++] = strdup("correlation coefficient");
        add_to_string_aligned(&LegendStr, "c.c.");
        setname[nsets++] = strdup("V_dens (kJ/mol)");
        add_to_string_aligned(&LegendStr, "V_dens");
        xvgr_legend(fp, nsets, setname, oenv);
        sfree(setname);

        fprintf(fp, "#\n# Legend for the following data columns:\n");
        fprintf(fp, "%s\n", LegendStr);
        sfree(LegendStr);

        fflush(fp);
    }

    return fp;

}


void Densfit::get_nstfit_from_env(t_commrec *cr, int nstlist)
{

    char *env;
    char  buf[STRLEN];
    int   env_nstfit;

    FILE *fp = df->out_dens;


    if (MASTER(cr))
    {
        env = getenv(NSTFIT_ENVVAR);

        if (env != NULL)
        {
            sprintf(buf, "Getting nstfit from environment variable %s=%s", NSTFIT_ENVVAR, env);
            fprintf(stderr, "%s%s\n", FitStr, buf);

            sscanf(env, "%d", &env_nstfit);

            if (env_nstfit >= 0)
            {
                this->nstfit = env_nstfit;
                if (fp)
                {
                    fprintf(fp, "# %s\n", buf);
                }
            }
            else
            {
                fprintf(stderr, "WARNING: Could not get a meaningful value for environment variable %s!\n", NSTFIT_ENVVAR);
            }
        }
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(this->nstfit), &this->nstfit, cr);
    }

    /* We check nstfit versus nstlist here so we can output a warning on stderr
     * AND a note in the density fitting output file */
    if (MASTER(cr))
    {
        if (this->nstfit > 1 && (this->nstfit % nstlist != 0) )
        {
            fprintf(stderr, "\n%s\n"
                    "WARNING: The frequency at which the fitting forces are calculated (nstfit=%d)\n"
                    "         is not a multiple of the neighbor searching frequency (nstlist=%d).\n"
                    "         Fitting forces will be calculated when the time step is a multiple\n"
                    "         of nstfit or nstlist, that could also be irregular intervals.\n\n", FitStr,
                    this->nstfit, nstlist);

            fprintf(fp, "# Note: will update fitting forces when the time step is a multiple of nstlist=%d or nstfit=%d.\n",
                    nstlist, this->nstfit);
        }

        if (this->nstfit < 1)
        {
            fprintf(stderr, "\n%s\n"
                    "WARNING: The frequency at which the density fitting potential is calculated (nstfit=%d) is < 1!\n"
                    "         No density fitting will be performed and no fitting output will be written.\n\n",
                    FitStr, this->nstfit);
            fprintf(fp, "# No output follows since nstfit < 1\n");
        }
    }
}

void Densfit::finish()
{
    if (df->out_dens)     // Close the density fitting output file
    {
        if (df->bAppend)
        {
            gmx_fio_fclose(df->out_dens);
        }
        else
        {
            xvgrclose(df->out_dens);
        }
    }
}

void Densfit::init(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
                   gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr, const gmx_output_env_t *oenv,
                   gmx_bool bAppend, gmx_bool bVerbose)
{
    int            i;
    real           grid_spacing;


    /* To be able to make the density fitting molecule whole: */
    rvec *x_pbc         = NULL;
    rvec *x_densfit_pbc = NULL;


    FitStr = FitStrMDRUN;

    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
    {
        gmx_fatal(FARGS, "%sOnly implemented for domain decomposition.\n", FitStr);
    }

    if (MASTER(cr) && bVerbose)
    {
        fprintf(stdout, "%sInitializing ...\n", FitStr);
    }

    /* In serial runs, the local atom indices are equal to the global ones */
    if (!PAR(cr))
    {
        this->nat_loc = this->nat;
        this->ind_loc = this->ind;
    }

    new_map_sim_from_ref(this, this->map_ref);

    /* Make space for the density fitting private data, allocate x, w, and f arrays, ... */
    snew(x_densfit_pbc, this->nat);                               /* There ... */
    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms);                                    /* There ... */
        //        m_rveccopy(mtop->natoms,x,x_pbc);
        copy_rvecn(x, x_pbc, 0, mtop->natoms);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        for (i = 0; i < this->nat; i++)
        {
            copy_rvec(x_pbc[this->ind[i]], x_densfit_pbc[i]);
        }
        sfree(x_pbc);                                       /* ... and back again */
    }
    if (PAR(cr))
    {
        gmx_bcast(this->nat*sizeof(rvec), x_densfit_pbc, cr);
    }

    grid_spacing = get_map_spacing(this->map_ref, MASTER(cr) ? stdout : NULL);


    setPrivateData(x_densfit_pbc, grid_spacing, MASTER(cr) && bVerbose, bAppend, PAR(cr));

    sfree(x_densfit_pbc);                                   /* ... and back again */

    if (MASTER(cr))
    {
        df->out_dens = open_out(opt2fn("-d", nfile, fnm), oenv);
        df->fn_map   = opt2fn("-mo", nfile, fnm);
    }
    else
    {
        df->out_dens = NULL;
        df->fn_map   = NULL;
    }

    /* Check whether densfit->nstfit is overwritten by environment variable */
    get_nstfit_from_env(cr, ir->nstlist);

    if (MASTER(cr) && !df->bAppend)
    {
        please_cite(fplog, "Tama2008");
    }

    /* This sum does not change, therefore we precalculate it */
    df->sum_rho_exp2 = calc_sum_rhoA_rhoB(this->map_ref.vox, this->map_ref.vox);

    //checkPerformance(ir);
}

void Densfit::assemble_atoms_for_spread( rvec       x[])
{
    int          i;
    for (i = 0; i < this->nat; i++)
    {
        copy_rvec(x[this->ind[i]], df->x_assembled[i]);
    }
}

void Densfit::dd_make_local_df_indices(gmx_domdec_t *dd)
{


    /* Local atoms of the density fitting group (these atoms will be spread) */
    dd_make_local_group_indices(dd->ga2la, this->nat, this->ind,
                                &this->nat_loc, &this->ind_loc, &this->nalloc_loc, this->df->c_ind);

    /* Indicate that the group's shift vectors for this structure need to be updated
     * at the next call to communicate_group_positions, since obviously we are in a NS step */
    df->bUpdateShifts = TRUE;
}

real Densfit::add_forces(rvec *f, gmx_unused t_commrec *cr, gmx_int64_t step, real time)
{
    int          ii, l;


    if (this->nstfit < 1)
    {
        return -1.0;
    }

/* Diagnostic output */
    if (step == -99)
    {
        char fn[STRLEN];
        make_filename("forces.txt", 0, step, fn);
        dump_x(cr->nodeid, step);
        dump_f(fn, cr);
    }

    if (df->out_dens && do_per_step(step, nstout))
    {
        fprintf(df->out_dens, "%12.5e%12g%12.5e%12.5e%12.7f%12.5e\n", time, df->temp, df->k, df->sigma, df->cc, df->k * (1.0 - df->cc));
    }

    for (l = 0; l < this->nat_loc; l++)
    {
        /* Get the right index of the local force, since typically not all local
         * atoms are subject to density fitting forces */
        ii = this->ind_loc[l];

        /* Add to local force */
        rvec_inc(f[ii], df->f_loc[l]);
    }

// return the correlation coefficient
    return df->cc;
}

void Densfit::do_densfit(real t, gmx_int64_t step, gmx_bool bOutputMap, t_inputrec *ir, t_commrec *cr, PaddedArrayRef< RVec> x, matrix box, gmx_wallcycle *wcycle)
{
    const char     *fn       = NULL;
    int             istart, nspread, nnodes;
    gmx_cycles_t    cycles_comp;                /* Cycles for the density fitting computations
                                                   only, does not count communication. This
                                                   counter is used for load-balancing              */

    /* If the user requested to have the density fitting forces calculated only
     * every N steps, we can skip the expensive do_densfit routine. However, we
     * must always calculate the fitting forces
     * - at the first step (obviously)
     * - after domain decompositions since the local atoms have changed
     */
    if (   ( !do_per_step(step, this->nstfit) && (FALSE == df->bUpdateShifts) )
           || (  this->nstfit < 1 ) )
    {
        return;
    }

    /* Update the time-dependent parameters like k, sigma, etc. */
    if (this->npoints > 1)
    {
        updateParametersThatAreTimeDependent(t);
    }

    /**************************************************************************/
    /* ASSEMBLE THE POSITIONS OF THE ATOMS USED FOR DENSITY FITTING           */

    /* Broadcast the DF positions such that every node has all of them
     * Every node contributes its local positions x and stores it in
     * the collective df->x_assembled array. Do all the communication first!  */
    wallcycle_start(wcycle, ewcDENSFIT_COMM);

    communicate_group_positions(cr, df->x_assembled, df->x_shifts, df->extra_shifts,
                                df->bUpdateShifts, as_rvec_array(x.data()), this->nat, this->nat_loc,
                                this->ind_loc, df->c_ind, df->x_old, box);

    /* If bUpdateShifts was TRUE then the shifts have just been updated in
     * communicate_group_positions. We do not need to update the shifts until
     * the next NS step */
    df->bUpdateShifts = FALSE;

    /* Put all atoms in the box (TODO: if we do that, we do not need to construct
     * a whole DF group before with communicate_group_positions!)
     */
    if (ir->ePBC != epbcNONE)
    {
        auto xAssembledArrayRef = arrayRefFromArray(reinterpret_cast<gmx::RVec *>(df->x_assembled), this->nat);
        put_atoms_in_box_omp(ir->ePBC, box, xAssembledArrayRef);
    }

    wallcycle_stop(wcycle, ewcDENSFIT_COMM);
    /* Now all nodes have all of the DF positions in df->x_assembled */

    /*                                                                        */
    /**************************************************************************/

    /**************************************************************************/
    /* SPREAD THE ATOMS ON THE DENSITY GRID                                   */
    /* Produces the density map. Each node spreads an equal part of all atoms,
     * this way we do not need load balancing here                            */
    wallcycle_start(wcycle, ewcDENSFIT_SPREAD);

    if (PAR(cr))
    {
        nnodes = cr->nnodes - cr->npmenodes;

        nspread = ceil( (real) this->nat / nnodes);
        istart  = cr->nodeid * nspread;

        if ( (nnodes - 1) == cr->nodeid)
        {
            nspread = this->nat - nspread*(nnodes-1);
        }
    }
    else
    {
        istart  = 0;
        nspread = this->nat;
    }

    assert(istart >= 0);
    assert(nspread >= 0);
    /* The reference map grows and shrinks with the box */
    //    get_reference_map_spacing(densfit, TRUE);
    //    df->spacing = couple_map_to_box(box, densfit);
    //    get_map_spacing(densfit, TRUE); TODO

    spread_atoms_low(&df->x_assembled[istart], nspread, box);

    wallcycle_stop(wcycle, ewcDENSFIT_SPREAD);

    if (PAR(cr))
    {
        wallcycle_start(wcycle, ewcDENSFIT_SUM_GRID);
        sum_grid(&(this->map_sim), cr);     /* Sum the density grid across all nodes */
        wallcycle_stop(wcycle, ewcDENSFIT_SUM_GRID);
    }
    /*                                                                        */
    /**************************************************************************/


    if (MASTER(cr))
    {
        if (bOutputMap)
        {
            char fn_with_step[STRLEN];

            if (this->bKeepAndNumberMaps)
            {
                make_filename(df->fn_map, 0, step, fn_with_step);
                fn = fn_with_step;
            }
            else
            {
                fn = df->fn_map;
            }
            MrcFile().write(fn, mapDataToGridData(this->map_sim));
            MrcFile().write("map_ref", mapDataToGridData(this->map_ref));
        }

        /* Calculate the correlation coefficient versus the reference map */
        df->cc = calc_correlation_coeff(nullptr);
        //        if (do_per_step(step, 10))
        //        {
        //            fprintf(stderr, "\rDensity fitting map correlation coefficient: %g\n", df->cc);
        //        }
    }


    /**************************************************************************/
    /* CALCULATE THE FORCES FROM THE DENSITY FITTING POTENTIAL                */
    /* Each rank calculates the forces on its home atoms                      */

    /* Start to count cycles for the dynamic load balancing now */
    cycles_comp = gmx_cycles_read();

    //    long int c1 = gmx_cycles_read();
    wallcycle_start(wcycle, ewcDENSFIT_FORCES);
    do_forces(box);
    wallcycle_stop(wcycle, ewcDENSFIT_FORCES);
    //    fprintf(stderr, "--- Forces took %g M cycles\n", 0.000001*(gmx_cycles_read()-c1));

    /* Stop the density fitting cycle counter and add the computation-only
     * cycles to the force cycles for load balancing */
    cycles_comp  = gmx_cycles_read() - cycles_comp;

    if (DOMAINDECOMP(cr) && wcycle)
    {
        dd_cycles_add(cr->dd, cycles_comp, ddCyclF);
    }
}
}  // namespace gmx
