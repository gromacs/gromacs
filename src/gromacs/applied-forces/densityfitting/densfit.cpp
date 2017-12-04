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

// #include "config.h"

#include "densfit.h"
#include "densfitdata.h"

#include <cassert>
#include <cstdio>
#include <vector>

#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdlib/sim_util.h"

#include "gromacs/mdlib/groupcoord.h"

#include "gromacs/domdec/domdec_struct.h"

#include "gromacs/timing/cyclecounter.h"

#include "gromacs/domdec/domdec.h"

#include "gromacs/fileio/griddataio.h"

#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/pbcutil/pbc.h"

#include "gromacs/fileio/mrcmetadata.h"

#include "gromacs/utility/txtdump.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/utility/strconvert.h"

#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/fileio/readinp.h"

static const char *FitStrMDRUN = "Density fitting: "; //!< Used if density
                                                      //!fitting output is
                                                      //!written from mdrun
static const char *FitStr = "";                       //!< Used in gmx map output
//! \brief Environment variable for setting nstfit
static const char *NSTFIT_ENVVAR = "GMX_NSTFIT";

#define MAXPTR 254



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

    gmx_sumf(nx * ny * nz, map->vox.data(), cr);
}

//! \brief Exponential function exp for a vector of n elements
void gmx_exp_vector(int n, real x[])
{
//    gmx_cycles_t c1 = gmx_cycles_read();

#if GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal sx = load<SimdReal>(x + i);
        sx = gmx::exp(sx);
        store(x + i, sx);
    }
#else       // GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i++)
    {
        x[i] = std::exp(x[i]);
    }
#endif      // GMX_SIMD_HAVE_REAL
    //    gmx_cycles_t c2 = gmx_cycles_read();
    //    fprintf(stderr, "exp vector %15d elements %15lld cycles  %15.3f per
    //    exp\n", n, c2-c1, (1.0*(c2-c1))/n);
}

/*! \brief Calculate one of the terms of the correlation coefficient
 *
 * Multiply each voxel of the map with the corresponding voxel of the other
 * map and sum up everything, i.e. compute
 *
 *  sum_ijk [ rhoA_ijk * rhoB_ijk ]
 *
 */
static double calc_sum_rhoA_rhoB(const std::vector<float> &vecA,
                                 const std::vector<float> &vecB)
{
    double sum_ijk = 0.;
    for (size_t i = 0; i < vecA.size(); i++)
    {
        sum_ijk += vecA[i] * vecB[i];
    }
    return sum_ijk;
}

//! \brief For a position vector x return all direct 27 image positions given
//! the box
int get_images(rvec x, matrix box, rvec img[], ivec img_startindex[],
               ivec img_stop_index[], real grid_spacing, int voxrange,
               const std::array<int, DIM> &map_dim)
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
                icenter    = roundf(image[XX] * OOgrid_spacing);
                istart[XX] = icenter - voxrange;
                istop[XX]  = icenter + voxrange;

                if ((istart[XX] <= map_dim[XX]) && (istop[XX] >= 0))
                {
                    icenter    = roundf(image[YY] * OOgrid_spacing);
                    istart[YY] = icenter - voxrange;
                    istop[YY]  = icenter + voxrange;

                    if ((istart[YY] <= map_dim[YY]) && (istop[YY] >= 0))
                    {

                        icenter    = roundf(image[ZZ] * OOgrid_spacing);
                        istart[ZZ] = icenter - voxrange;
                        istop[ZZ]  = icenter + voxrange;

                        if ((istart[ZZ] <= map_dim[ZZ]) && (istop[ZZ] >= 0))
                        {
                            /* This is our new image */
                            copy_rvec(image, img[count]);

                            /* Cut away the volume where we do not have a grid */
                            for (m = 0; m < DIM; m++)
                            {
                                istart[m] = std::max(istart[m], 0); /* istart >= 0       */
                                istop[m]  =
                                    std::min(istop[m], map_dim[m]); /* istop  <= map_dim */
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
        SimdReal sx = load<SimdReal>(x + i);
        sx = gmx::erf(sx);
        store(x + i, sx);
    }
#else       // GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i++)
    {
        x[i] = std::erf(x[i]);
    }
#endif      // GMX_SIMD_HAVE_REAL
    //    gmx_cycles_t c2 = gmx_cycles_read();
    //    fprintf(stderr, "erf vector %15d elements %15lld cycles  %15.3f per erf
    //    (SIMD width is %d)\n", n, c2-c1, (1.0*(c2-c1))/n, GMX_SIMD_REAL_WIDTH);
}
}

//! \brief local density fitting data
class LocalDensfitData
{
    public:
        LocalDensfitData(const DensfitData &parameters,
                         rvec *x_densfit_whole, /* Can be NULL!                       */
                         real grid_spacing, gmx_bool bVerbose, gmx_bool bAppend,
                         gmx_bool bParallel);
        ~LocalDensfitData() = default;
        real k;                    /**< Current value k of strength constant of map
                                        potential V_fit = k*(1 - cc) with cc being the
                                        correlation coefficient                        */
        real        sigma;         /**< Current value of Gaussian width used for
                                        spreading atomic positions on the map          */
        real        temp;          /**< Current value of temperature                   */
        FILE       *out_dens;      /**< Output file for density fitting data           */
        const char *fn_map;        /**< Output filename for simulated density maps     */
        real        spacing;       /**< Map spacing in nm, derived from map data       */
        gmx_bool    bAppend;       /**< Are we appending to existing output files?     */
        gmx_bool    bVerbose;      /**< -v flag from command line                      */
        real        cc;            /**< Correlation coefficient of the two maps        */
        rvec       *x_assembled;   /**< Just the positions that are spread             */
        rvec       *f_loc;         /**< Forces on the local atoms of the
                                        density fitting group                          */
        real       *w_assembled;   /**< Weights for the assembled positions            */
        rvec       *x_old;         /**< x_assembled from last step                     */
        ivec       *x_shifts;      /**< To make the molecule whole                     */
        ivec       *extra_shifts;  /**< Shifts added since last NS                     */
        int        *c_ind;         /**< Where is this atom stored in the coll array?   */
        gmx_bool    bUpdateShifts; /**< Do we have to calculate new shift vectors
                                        from x_assembled and x_old?                    */
        int         voxrange;      /**< Max. number of voxels to be computed for a
                                        single atom in a single dimension x, y, or z   */
        int         vox_per_atom;  /**< Total number of voxels for one atom (x,y, + z) */
        double      sum_rho_exp2;  /**< Term needed for the calculation of the
                                        correlation coefficient and the forces         */

        /* The following two temporary vectors (one for each OpenMP thread) store
         * erf values around a single atoms, thus we can compute them all in one go,
         * and with SIMD acceleration */
        std::vector < real, gmx::AlignedAllocator < real>>
        *erfVector; /**< pointer to vector of erf values */
        std::vector < real, gmx::AlignedAllocator < real>>
        *expVector; /**< same for exp values             */
    private:
        /* TODO: crosscheck with code in new_t_gmx_densfit(), which depends on sigma!!
         */
        //! \brief TODO
        void initParametersThatCanBeTimeDependent(
            const DensfitData &densfitParameters);
};

//! \brief Allocate and initialize the arrays for assembled positions, forces
//! and atomic weights
LocalDensfitData::LocalDensfitData(const DensfitData &parameters,
                                   rvec *x_densfit_whole, /* Can be NULL! */
                                   real grid_spacing, gmx_bool bVerbose, gmx_bool bAppend,
                                   gmx_bool bParallel)
{

    this->bVerbose = bVerbose;
    this->spacing  = grid_spacing;

    snew(this->x_assembled, parameters.nAtoms());
    snew(this->f_loc, parameters.nAtoms());
    snew(this->w_assembled, parameters.nAtoms());
    snew(this->x_old, parameters.nAtoms());
    snew(this->x_shifts, parameters.nAtoms());
    snew(this->extra_shifts, parameters.nAtoms());
    snew(this->c_ind, parameters.nAtoms());

    if (this->bVerbose)
    {
        fprintf(stdout, "%sSetting all %d atomic weights to unity.\n", FitStr,
                parameters.nAtoms());
    }
    for (int i = 0; i < parameters.nAtoms(); i++)
    {
        this->w_assembled[i] = 1.0;
    }

    this->bUpdateShifts = TRUE;
    this->bAppend       = bAppend;

    if (!bParallel)
    {
        for (int i = 0; i < parameters.nAtoms(); i++)
        {
            this->c_ind[i] = i;
        }
    }

    /* In mdruns we have to make the density fitting structure whole */
    if (NULL != x_densfit_whole)
    {
        for (int i = 0; i < parameters.nAtoms(); i++)
        {
            copy_rvec(x_densfit_whole[i], this->x_old[i]);
        }
    }

    /* Initialize the possibly time-dependent parameters */
    initParametersThatCanBeTimeDependent(parameters);
}

void LocalDensfitData::initParametersThatCanBeTimeDependent(
        const DensfitData &densfitParameters)
{
    int ompthreads;

    this->sigma = densfitParameters.currentSigma(0.);
    this->k     = densfitParameters.currentK(0.);
    this->temp  = densfitParameters.currentT(0.);

    /* For erf() and exp() performance, we want all values for one atom in one
     * linear array (not in three). Therefore, we compute the maximum number of
     * voxels to be computed for a single atom. How many voxels around each atom
     * we have to spread the density on is saved in voxrange (for each dimension
     * x, y, and z */
    this->voxrange =
        ceil(0.5 * this->sigma * densfitParameters.dist() / this->spacing);
    this->vox_per_atom = 3 * (2 * this->voxrange + 1);

    ompthreads      = std::max(1, gmx_omp_nthreads_get(emntDefault));
    this->erfVector =
        new std::vector < real, gmx::AlignedAllocator < real>>[ompthreads];
    this->expVector =
        new std::vector < real, gmx::AlignedAllocator < real>>[ompthreads];
    for (int i = 0; i < ompthreads; i++)
    {
        int extra = 0;
#if GMX_SIMD_HAVE_REAL
        extra = GMX_SIMD_REAL_WIDTH;
#endif
        this->erfVector[i].reserve(this->vox_per_atom + extra);
        this->expVector[i].reserve(this->vox_per_atom + extra);
        // + GMX_SIMD_REAL_WIDTH adds some extra elements at the end of the array so
        // that
        // we can always read GMX_SIMD_REAL_WIDTH elements, even near the end of the
        // vector.
    }
}

Densfit::Densfit() {}

Densfit::Densfit(real sigma, real sigma_dist, real k, /* Spring constant */
                 int nat,
                 int *ind,                            /* Which of the atoms should be used for spreading */
                 real grid_spacing, bool bVerbose)
{
    parameters_ = std::unique_ptr<DensfitData>(new DensfitData(1, sigma, k, sigma_dist, nat, ind, bVerbose));
    setPrivateData(nullptr, grid_spacing, bVerbose, 0, 0);
};

Densfit::Densfit(const DensfitData &parameters)
{
    parameters_ = std::unique_ptr<DensfitData>(new DensfitData(parameters));
}

Densfit::~Densfit(){};

void Densfit::makeGroups(char *groupname, t_blocka *grps, char **gnames)
{
    parameters_->makeGroups(groupname, grps, gnames);
}
void Densfit::updateParametersThatAreTimeDependent(real time)
{
    int ompthreads;

    df_->temp = parameters_->currentT(time);

    assert(df_->temp >= -1);
    df_->sigma = parameters_->currentSigma(time);

    assert(df_->sigma > 0.0);
    df_->k = parameters_->currentSigma(time);

    /* For erf() and exp() performance, we want all values for one atom in one
     * linear array (not in three). Therefore, we compute the maximum number of
     * voxels to be computed for a single atom. How many voxels around each atom
     * we have to spread the density on is saved in voxrange (for each dimension
     * x, y, and z */
    df_->voxrange     = ceil(0.5 * df_->sigma * parameters_->dist() / df_->spacing);
    df_->vox_per_atom = 3 * (2 * df_->voxrange + 1);

    ompthreads = std::max(1, gmx_omp_nthreads_get(emntDefault));
    // TODO: initialize with the maximum vox_per_atom
    for (int i = 0; i < ompthreads; i++)
    {
        int extra = 0;
#if GMX_SIMD_HAVE_REAL
        extra = GMX_SIMD_REAL_WIDTH;
#endif
        // Give each OpenMP thread its own local storage; erfVector[i] is for OpenMP
        // thread number i
        df_->erfVector[i].reserve(df_->vox_per_atom + extra);
        df_->expVector[i].reserve(df_->vox_per_atom + extra);
    }
}

void Densfit::dump_x(int nodeid, gmx_int64_t step)
{
    FILE *fpout = nullptr;
    char  fn[STRLEN];
    char  buf[256];

    gmx_step_str(step, buf);
    sprintf(fn, "x_step%s_node%d.txt", buf, nodeid);
    fpout = fopen(fn, "w");

    for (int i = 0; i < parameters_->nAtoms(); i++)
    {
        fprintf(fpout, "%4d %12.3e %12.3e %12.3e\n", i, df_->x_assembled[i][XX],
                df_->x_assembled[i][YY], df_->x_assembled[i][ZZ]);
    }

    fclose(fpout);
}

void Densfit::dump_f(const char *fn, t_commrec *cr)
{
    int      i;
    rvec    *f;
    gmx_bool bWrite;

    FILE    *fpout = NULL;

    /* Assemble the forces */
    snew(f, parameters_->nAtoms());
    /* Each node copies its local forces into the collective array */
    for (i = 0; i < nat_loc_; i++)
    {
        copy_rvec(df_->f_loc[i], f[df_->c_ind[i]]);
    }

    /* Reduce */
    if (cr && PAR(cr))
    {
        gmx_sum(3 * parameters_->nAtoms(), f[0], cr);
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
            fpout = df_->out_dens;
        }

        fprintf(stderr, "%sDumping forces to file %s.\n", FitStr, fn ? fn : "");

        for (i = 0; i < parameters_->nAtoms(); i++)
        {
            fprintf(fpout, "#     f[%5d]={%12.5e, %12.5e, %12.5e}\n", i, f[i][XX],
                    f[i][YY], f[i][ZZ]);
        }

        if (fn)
        {
            gmx_fio_fclose(fpout);
        }
    }
    sfree(f);
}

void Densfit::setPrivateData(rvec *x_densfit_whole, real grid_spacing,
                             bool bVerbose, bool bAppend, bool bParallel)
{
    /* Initialize also the private data structure: */
    this->df_ = std::unique_ptr<LocalDensfitData>(
                new LocalDensfitData(*parameters_, x_densfit_whole, grid_spacing, bVerbose,
                                     bAppend, bParallel));
}

void Densfit::spread_atoms_low(rvec x[], int natoms, matrix box)
{

    int    iz, iy, ix, i, j, n, th;
    int    nx, ny, nz, ngrid;
    double V_vox;
    double prefactor;
    double left;
    real   ax, ay, az, ayaz;
    float *ptr; /* Pointer in that array                          */
    ivec   elements;
    rvec   pos;
    real  *aarr[3];
    int    image, num_images;
    rvec   img[27];    /* Stores PBC images of a position x as well as x */
    ivec   istart[27]; /* Store for each image the part of the grid ...  */
    ivec   istop[27];  /* ... where stuff has to be calculated           */

    double sigma2 = df_->sigma * df_->sigma;
    double c_fac  = sqrt(1.5 / sigma2);
    int    nth    = std::max(1, gmx_omp_nthreads_get(emntDefault));

    //    gmx_cycles_t c1 = gmx_cycles_read();

    auto &sim_dens_1d = this->map_sim_.vox; /* Calculated density map */
    V_vox      = pow(df_->spacing, 3);
    prefactor  = pow(M_PI * sigma2 / 6.0, 1.5);
    prefactor /= V_vox;
    nx         = this->map_sim_.map_dim[XX];
    ny         = this->map_sim_.map_dim[YY];
    nz         = this->map_sim_.map_dim[ZZ];
    ngrid      = nx * ny * nz;

/* Zero out the grid */
#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(sim_dens_1d, ngrid, nth) default(none)
    for (i = 0; i < ngrid * nth; i++)
    {
        sim_dens_1d[i] = 0.0;
    }

    const auto dfp        = df_.get();
    const auto mapsim_dim = map_sim_.map_dim;
/* Spread the local atoms on the grid, therefore loop over all atoms of the
 * density fitting group */
#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(nx, ny, nz, ngrid, natoms, x, box, c_fac, \
    sim_dens_1d) private(ptr, aarr, istart, istop, image, pos, \
    num_images, img, i, ix, iy, iz, th, elements, \
    ax, ay, az, ayaz, left) default(none)
    for (n = 0; n < natoms; n++)
    {
        th = gmx_omp_get_thread_num(); /* Thread ID */

        /* Determine the grid points around the atom and its images for which we
         * need to calculate stuff */
        num_images = get_images(x[n], box, img, istart, istop, dfp->spacing,
                                dfp->voxrange, mapsim_dim);

        /* Loop over all periodic images of the atom */
        for (image = 0; image < num_images; image++)
        {
            copy_rvec(img[image], pos); // periodic image of the position

            elements[XX] = istop[image][XX] + 1 - istart[image][XX];
            elements[YY] = istop[image][YY] + 1 - istart[image][YY];
            elements[ZZ] = istop[image][ZZ] + 1 - istart[image][ZZ];

            //            int sum = elements[XX] + elements[YY] + elements[ZZ];

            /* For erf() performance, we want all values for one atom in one linear
             * array (not
             * int three). However, we build pointers to access it as x, y, and z
             * arrays.
             *
             * Each OpenMP thread th accesses its own temporary storage vector
             */
            aarr[XX] = &dfp->erfVector[th][0];
            aarr[YY] = &dfp->erfVector[th][0 + elements[XX]];
            aarr[ZZ] = &dfp->erfVector[th][0 + elements[XX] + elements[YY]];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in a temporary array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            i = 0;
            for (ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                left          = (ix - 0.5) * dfp->spacing - pos[XX];
                aarr[XX][i++] = c_fac * left;
            }

            i = 0;
            for (iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                left          = (iy - 0.5) * dfp->spacing - pos[YY];
                aarr[YY][i++] = c_fac * left;
            }

            i = 0;
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                left          = (iz - 0.5) * dfp->spacing - pos[ZZ];
                aarr[ZZ][i++] = c_fac * left;
            }
            gmx_erf_vector(elements[XX] + elements[YY] + elements[ZZ], aarr[XX]);

            /* Transform (in place) into a new array that directly stores the
             * differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix + 1] - aarr[XX][ix];
            }
            for (iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy + 1] - aarr[YY][iy];
            }
            for (iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz + 1] - aarr[ZZ][iz];
            }

            /* Spread the density for this atom */
            for (iz = istart[image][ZZ]; iz < istop[image][ZZ]; iz++)
            {
                az = aarr[ZZ][iz - istart[image][ZZ]];

                for (iy = istart[image][YY]; iy < istop[image][YY]; iy++)
                {
                    ay   = aarr[YY][iy - istart[image][YY]];
                    ayaz = ay * az;
                    ptr  = &sim_dens_1d[th * ngrid + nx * (ny * iz + iy) +
                                        istart[image][XX]];
                    for (ix = istart[image][XX]; ix < istop[image][XX]; ix++)
                    {
                        ax = aarr[XX][ix - istart[image][XX]];

                        /* Assign value to voxel */
                        //                        sim_dens_1d[nz*ny*ix + nz*iy + iz] +=
                        //                        ax*ay*az;
                        *ptr += ax * ayaz;
                        ptr++;
                    }
                }
            }
        } /* End of loop over images of one atom */
    }     /* End of loop over atoms */

#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(nth, sim_dens_1d, ngrid, prefactor) private(j) default(none)
    for (i = 0; i < ngrid; i++)
    {
        for (j = 1; j < nth; j++)
        {
            sim_dens_1d[i] += sim_dens_1d[i + j * ngrid];
        }
        sim_dens_1d[i] *= prefactor;
    }
}

void Densfit::setSpacingFromMap(const t_mapdata &map)
{
    df_->spacing = 0.;
    for (int i = 0; i < DIM; i++)
    {
        df_->spacing += A2NM * map.cell[i] / map.map_dim[i];
    }
    df_->spacing /= 3.;
}

const t_mapdata &Densfit::simulatedMap() { return map_sim_; }

void Densfit::spread_atoms(matrix box)
{
    //    double         t0 = 0;
    //    long int       c1 = 0; /* make compiler happy */

    int *map_dim = this->map_sim_.map_dim.data();
    if (df_->bVerbose)
    {
        fprintf(stderr,
                "%sSpreading %d atomic positions  on a %d x %d x %d grid ...\n",
                FitStr, parameters_->nAtoms(), map_dim[XX], map_dim[YY], map_dim[ZZ]);
        //        t0 = gettime();
        //        c1 = gmx_cycles_read();
    }

    if (df_->spacing <= 0)
    {
        gmx_fatal(FARGS, "%sSpacing must not be zero!", FitStr);
    }

    /* Make the density map */
    spread_atoms_low(df_->x_assembled, parameters_->nAtoms(), box);

    //    if (df->bVerbose)
    //    {
    //        fprintf(stderr, "%sSpreading took %g seconds (%g M cycles).\n",
    //        FitStr,
    //                gettime() - t0, 0.000001*(gmx_cycles_read()-c1));
    //    }
}

void Densfit::do_forces(matrix box)
{
    real                     *aarr[3], *garr[3];
    double                    ref_dens, sim_dens;

    rvec                      img[27];                                       /* Stores PBC images of a position x as well as x */
    ivec                      istart[27];                                    /* Store for each image the part of the grid ...  */
    ivec                      istop[27];                                     /* ... where stuff has to be calculated           */

    const std::vector<float> &sim_dens_1d = map_sim_.vox;                    /* Calculated density map */
    const std::vector<float> &ref_dens_1d = parameters_->referenceMap().vox; /* Reference density map */

    if (nat_loc_ <= 0)
    {
        return;
    }

    // cppcheck-suppress unreadVariable
    int gmx_unused nth =
        std::max(1, gmx_omp_nthreads_get(emntDefault)); // number of threads

/* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int l = 0; l < nat_loc_; l++)
    {
        clear_rvec(df_->f_loc[l]);
    }

    /* Calculate various prefactors */
    double sum_rho_exp2        = df_->sum_rho_exp2;
    double sum_rho_sim2        = calc_sum_rhoA_rhoB(map_sim_.vox, map_sim_.vox);
    double sum_rho_exp_rho_sim =
        calc_sum_rhoA_rhoB(parameters_->referenceMap().vox, map_sim_.vox);

    double term1_prefactor = df_->k / sqrt(sum_rho_exp2 * sum_rho_sim2);
    double term2_prefactor = -df_->k * sum_rho_exp_rho_sim /
        (sqrt(sum_rho_exp2) * pow(sum_rho_sim2, 1.5));

    double     V_vox            = df_->spacing * df_->spacing * df_->spacing;
    double     prefactor        = -M_PI * df_->sigma * df_->sigma / (6.0 * V_vox);
    double     OOsigma2_fac     = 1.5 / (df_->sigma * df_->sigma);
    double     sqrtOOsigma2_fac = sqrtf(OOsigma2_fac);

    int        nx = map_sim_.map_dim[XX];
    int        ny = map_sim_.map_dim[YY];

    auto      *dfp        = df_.get();
    const auto mapsim_dim = map_sim_.map_dim;
#pragma omp parallel for num_threads(nth) schedule(static) private( \
    istart, istop, aarr, garr, ref_dens, sim_dens, img) \
    shared(dfp, ref_dens_1d, sim_dens_1d, box, nx, ny, OOsigma2_fac, \
    sqrtOOsigma2_fac, prefactor, term1_prefactor, \
    term2_prefactor) default(none)
    /* Loop over all atoms of the density fitting group */
    for (int l = 0; l < nat_loc_; l++)
    {
        int  th = gmx_omp_get_thread_num();

        dvec term1_sum;
        dvec term2_sum;

        /* Loop over all voxels and calculate term1 and term2 for the force
         * simultaneously */
        clear_dvec(term1_sum);
        clear_dvec(term2_sum);

        /* Determine the grid points around the atom and its images which
         * have a contribution to the force on this atom.
         * Get the position of this atom of the density fitting group, we must
         * use the positions from the assembled array, because only these have
         * been put in the box in do_this() */
        int num_images = get_images(dfp->x_assembled[dfp->c_ind[l]], box, img, istart,
                                    istop, dfp->spacing, dfp->voxrange, mapsim_dim);

        /* Loop over all periodic images of the atom */
        for (int image = 0; image < num_images; image++)
        {
            rvec pos;
            copy_rvec(img[image], pos); // periodic image of the position

            ivec elements;
            /* The calculation of drho_dxl is expensive. We only calculate it, where
             * there is actually simulated density, i.e. in the given sigma-range
             * applied
             * during spreading. */
            elements[XX] = istop[image][XX] + 1 - istart[image][XX];
            elements[YY] = istop[image][YY] + 1 - istart[image][YY];
            elements[ZZ] = istop[image][ZZ] + 1 - istart[image][ZZ];

            /* For erf() performance, we want all values for one atom in one linear
             * array (not
             * int three). However, we build pointers to access it as x, y, and z
             * arrays.
             *
             * Each OpenMP thread th accesses its own temporary storage vector
             */
            aarr[XX] = &dfp->erfVector[th][0];
            aarr[YY] = &dfp->erfVector[th][0 + elements[XX]];
            aarr[ZZ] = &dfp->erfVector[th][0 + elements[XX] + elements[YY]];

            /* Same for the temporary array that keeps the exp values */
            garr[XX] = &dfp->expVector[th][0];
            garr[YY] = &dfp->expVector[th][0 + elements[XX]];
            garr[ZZ] = &dfp->expVector[th][0 + elements[XX] + elements[YY]];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in the tmpgrid[] array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            int i = 0;
            for (int ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                double left          = (ix - 0.5) * dfp->spacing - pos[XX];
                garr[XX][i]   = -OOsigma2_fac * left * left;
                aarr[XX][i++] = sqrtOOsigma2_fac * left;
            }

            i = 0;
            for (int iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                double left          = (iy - 0.5) * dfp->spacing - pos[YY];
                garr[YY][i]   = -OOsigma2_fac * left * left;
                aarr[YY][i++] = sqrtOOsigma2_fac * left;
            }

            i = 0;
            for (int iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                double left          = (iz - 0.5) * dfp->spacing - pos[ZZ];
                garr[ZZ][i]   = -OOsigma2_fac * left * left;
                aarr[ZZ][i++] = sqrtOOsigma2_fac * left;
            }

            /* Call erf() and exp() for all input values for this atom in one go */
            gmx_erf_vector(elements[XX] + elements[YY] + elements[ZZ], aarr[XX]);
            gmx_exp_vector(elements[XX] + elements[YY] + elements[ZZ], garr[XX]);

            /* Transform (in place) into a new array that directly stores the
             * differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (int ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix + 1] - aarr[XX][ix];
                garr[XX][ix] = garr[XX][ix + 1] - garr[XX][ix];
            }
            for (int iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy + 1] - aarr[YY][iy];
                garr[YY][iy] = garr[YY][iy + 1] - garr[YY][iy];
            }
            for (int iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz + 1] - aarr[ZZ][iz];
                garr[ZZ][iz] = garr[ZZ][iz + 1] - garr[ZZ][iz];
            }

            for (int iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                real erfz = aarr[ZZ][iz];
                real expz = garr[ZZ][iz];

                for (int iy = 0; iy < elements[YY] - 1; iy++)
                {
                    real erfy = aarr[YY][iy];
                    real expy = garr[YY][iy];

                    real erfy_erfz = erfy * erfz;
                    real expy_erfz = expy * erfz;
                    real erfy_expz = erfy * expz;

                    int  ind = ny * nx * (iz + istart[image][ZZ]) +
                        nx * (iy + istart[image][YY]) + istart[image][XX];

                    auto ptr_ref = ref_dens_1d.begin()+ind;
                    auto ptr_sim = sim_dens_1d.begin()+ind;

                    /* Calculate d/dx_l rho_sim_ijk */
                    for (int ix = 0; ix < elements[XX] - 1; ix++)
                    {
                        ref_dens = (double)*ptr_ref;
                        sim_dens = (double)*ptr_sim;

                        real erfx = aarr[XX][ix];
                        real expx = garr[XX][ix];

                        rvec drho_dxl;

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

        dvec tmp1;
        dvec tmp2;
        dsvmul(prefactor * term1_prefactor, term1_sum, tmp1);
        dsvmul(prefactor * term2_prefactor, term2_sum, tmp2);

        dvec force;
        dvec_add(tmp1, tmp2, force);

        dfp->f_loc[l][XX] = force[XX];
        dfp->f_loc[l][YY] = force[YY];
        dfp->f_loc[l][ZZ] = force[ZZ];

    }   /* end of loop over atoms */
}

real Densfit::calc_correlation_coeff(FILE *log)
{
    int    i;
    double numerator, denominator;
    real   cc;

    if (log)
    {
        fprintf(log, "%sCalculating the correlation coefficient of the two maps.\n",
                FitStr);
    }

    for (i = 0; i < 3; i++)
    {
        if (parameters_->referenceMap().map_dim[i] != map_sim_.map_dim[i])
        {
            gmx_fatal(FARGS, "Map dimensions must agree");
        }
    }

    numerator = calc_sum_rhoA_rhoB(parameters_->referenceMap().vox, map_sim_.vox);

    denominator = sqrt(calc_sum_rhoA_rhoB(parameters_->referenceMap().vox, parameters_->referenceMap().vox) *
                       calc_sum_rhoA_rhoB(map_sim_.vox, map_sim_.vox));

    cc = numerator / denominator;

    if (log)
    {
        fprintf(log, "%sThe correlation coefficient is %15.10f\n", FitStr, cc);
    }

    return cc;
}

/*! \brief Open output file and print some general information about density
 * fitting
 *
 * Call on master only
 */
FILE *Densfit::open_out(const char *fn, const gmx_output_env_t *oenv)
{
    FILE        *fp;

    if (df_->bAppend)
    {
        fp = gmx_fio_fopen(fn, "a");
    }
    else
    {
        fp = xvgropen(fn, "Density fitting correlation coefficient", "Time (ps)",
                      "correlation coefficient", oenv);
        parameters_->printToOut(fp, oenv);


        fflush(fp);
    }

    return fp;
}

void Densfit::get_nstfit_from_env(t_commrec *cr, int nstlist)
{
    parameters_->setNstFitFromEnv(cr, nstlist, df_->out_dens);
}

void Densfit::finish()
{
    if (df_->out_dens) // Close the density fitting output file
    {
        if (df_->bAppend)
        {
            gmx_fio_fclose(df_->out_dens);
        }
        else
        {
            xvgrclose(df_->out_dens);
        }
    }
}

void Densfit::init(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
                   gmx_mtop_t *mtop, rvec *x, matrix box, t_commrec *cr,
                   const gmx_output_env_t *oenv, gmx_bool bAppend,
                   gmx_bool bVerbose)
{
    int  i;
    real grid_spacing;

    /* To be able to make the density fitting molecule whole: */
    rvec *x_pbc         = NULL;
    rvec *x_densfit_pbc = NULL;

    FitStr = FitStrMDRUN;

    if ((PAR(cr)) && !DOMAINDECOMP(cr))
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
        nat_loc_ = parameters_->nAtoms();
        ind_loc_ = parameters_->indices();
    }

    new_map_sim_from_ref(parameters_->referenceMap());

    /* Make space for the density fitting private data, allocate x, w, and f
     * arrays, ... */
    snew(x_densfit_pbc, parameters_->nAtoms()); /* There ... */
    if (MASTER(cr))
    {
        /* Remove pbc and prepare a whole molecule for new_t_gmx_densfit */
        snew(x_pbc, mtop->natoms); /* There ... */
        //        m_rveccopy(mtop->natoms,x,x_pbc);
        copy_rvecn(x, x_pbc, 0, mtop->natoms);
        do_pbc_first_mtop(NULL, ir->ePBC, box, mtop, x_pbc);
        for (i = 0; i < parameters_->nAtoms(); i++)
        {
            copy_rvec(x_pbc[parameters_->indices()[i]], x_densfit_pbc[i]);
        }
        sfree(x_pbc); /* ... and back again */
    }
    if (PAR(cr))
    {
        gmx_bcast(parameters_->nAtoms() * sizeof(rvec), x_densfit_pbc, cr);
    }

    grid_spacing = get_map_spacing(parameters_->referenceMap(), MASTER(cr) ? stdout : NULL);

    setPrivateData(x_densfit_pbc, grid_spacing, MASTER(cr) && bVerbose, bAppend,
                   PAR(cr));

    sfree(x_densfit_pbc); /* ... and back again */

    if (MASTER(cr))
    {
        df_->out_dens = open_out(opt2fn("-d", nfile, fnm), oenv);
        df_->fn_map   = opt2fn("-mo", nfile, fnm);
    }
    else
    {
        df_->out_dens = NULL;
        df_->fn_map   = NULL;
    }

    /* Check whether densfit->nstfit is overwritten by environment variable */
    get_nstfit_from_env(cr, ir->nstlist);

    if (MASTER(cr) && !df_->bAppend)
    {
        please_cite(fplog, "Tama2008");
    }

    /* This sum does not change, therefore we precalculate it */
    df_->sum_rho_exp2 = calc_sum_rhoA_rhoB(parameters_->referenceMap().vox, parameters_->referenceMap().vox);

    // checkPerformance(ir);
}

void Densfit::assemble_atoms_for_spread(rvec x[])
{
    int i;
    for (i = 0; i < parameters_->nAtoms(); i++)
    {
        copy_rvec(x[parameters_->indices()[i]], df_->x_assembled[i]);
    }
}

void Densfit::dd_make_local_df_indices(gmx_domdec_t *dd)
{

    /* Local atoms of the density fitting group (these atoms will be spread) */
    dd_make_local_group_indices(dd->ga2la, parameters_->nAtoms(), parameters_->indices(),
                                &nat_loc_, &ind_loc_, &nalloc_loc_, df_->c_ind);

    /* Indicate that the group's shift vectors for this structure need to be
     * updated
     * at the next call to communicate_group_positions, since obviously we are in
     * a NS step */
    df_->bUpdateShifts = TRUE;
}

real Densfit::add_forces(rvec *f, gmx_unused t_commrec *cr, gmx_int64_t step,
                         real time)
{
    int ii, l;

    if (parameters_->nStepsFit() < 1)
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

    if (df_->out_dens && do_per_step(step, parameters_->nStepsMapOutput()))
    {
        fprintf(df_->out_dens, "%12.5e%12g%12.5e%12.5e%12.7f%12.5e\n", time,
                df_->temp, df_->k, df_->sigma, df_->cc, df_->k * (1.0 - df_->cc));
    }

    for (l = 0; l < nat_loc_; l++)
    {
        /* Get the right index of the local force, since typically not all local
         * atoms are subject to density fitting forces */
        ii = ind_loc_[l];

        /* Add to local force */
        rvec_inc(f[ii], df_->f_loc[l]);
    }

    // return the correlation coefficient
    return df_->cc;
}

void Densfit::do_densfit(real t, gmx_int64_t step, t_inputrec *ir,
                         t_commrec *cr, PaddedArrayRef<RVec> x, matrix box,
                         gmx_wallcycle *wcycle)
{
    const char  *fn = NULL;
    int          istart, nspread, nnodes;
    gmx_cycles_t cycles_comp; /* Cycles for the density fitting computations
                                 only, does not count communication. This
                                 counter is used for load-balancing */

    /* If the user requested to have the density fitting forces calculated only
     * every N steps, we can skip the expensive do_densfit routine. However, we
     * must always calculate the fitting forces
     * - at the first step (obviously)
     * - after domain decompositions since the local atoms have changed
     */
    if ((!do_per_step(step, parameters_->nStepsFit()) && (!df_->bUpdateShifts)) || (parameters_->nStepsFit() < 1))
    {
        return;
    }

    /* Update the time-dependent parameters like k, sigma, etc. */
    if (parameters_->timePoints() > 1)
    {
        updateParametersThatAreTimeDependent(t);
    }

    /**************************************************************************/
    /* ASSEMBLE THE POSITIONS OF THE ATOMS USED FOR DENSITY FITTING           */

    /* Broadcast the DF positions such that every node has all of them
     * Every node contributes its local positions x and stores it in
     * the collective df->x_assembled array. Do all the communication first!  */
    wallcycle_start(wcycle, ewcDENSFIT_COMM);

    communicate_group_positions(
            cr, df_->x_assembled, df_->x_shifts, df_->extra_shifts, df_->bUpdateShifts,
            as_rvec_array(x.data()), parameters_->nAtoms(), nat_loc_, ind_loc_,
            df_->c_ind, df_->x_old, box);

    /* If bUpdateShifts was TRUE then the shifts have just been updated in
     * communicate_group_positions. We do not need to update the shifts until
     * the next NS step */
    df_->bUpdateShifts = FALSE;

    /* Put all atoms in the box (TODO: if we do that, we do not need to construct
     * a whole DF group before with communicate_group_positions!)
     */
    if (ir->ePBC != epbcNONE)
    {
        auto xAssembledArrayRef = arrayRefFromArray(
                    reinterpret_cast<gmx::RVec *>(df_->x_assembled), parameters_->nAtoms());
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

        nspread = ceil((real)parameters_->nAtoms() / nnodes);
        istart  = cr->nodeid * nspread;

        if ((nnodes - 1) == cr->nodeid)
        {
            nspread = parameters_->nAtoms() - nspread * (nnodes - 1);
        }
    }
    else
    {
        istart  = 0;
        nspread = parameters_->nAtoms();
    }

    assert(istart >= 0);
    assert(nspread >= 0);
    /* The reference map grows and shrinks with the box */
    //    get_reference_map_spacing(densfit, TRUE);
    //    df->spacing = couple_map_to_box(box, densfit);
    //    get_map_spacing(densfit, TRUE); TODO

    spread_atoms_low(&df_->x_assembled[istart], nspread, box);

    wallcycle_stop(wcycle, ewcDENSFIT_SPREAD);

    if (PAR(cr))
    {
        wallcycle_start(wcycle, ewcDENSFIT_SUM_GRID);
        sum_grid(&(map_sim_), cr); /* Sum the density grid across all nodes */
        wallcycle_stop(wcycle, ewcDENSFIT_SUM_GRID);
    }
    /*                                                                        */
    /**************************************************************************/

    if (MASTER(cr))
    {
        if (do_per_step(step, parameters_->nStepsMapOutput()))
        {
            char fn_with_step[STRLEN];

            if (parameters_->keepAndNumberMaps())
            {
                make_filename(df_->fn_map, 0, step, fn_with_step);
                fn = fn_with_step;
            }
            else
            {
                fn = df_->fn_map;
            }
            MrcFile().write(fn, mapDataToGridData(map_sim_));
            MrcFile().write("map_ref", mapDataToGridData(parameters_->referenceMap()));
        }

        /* Calculate the correlation coefficient versus the reference map */
        df_->cc = calc_correlation_coeff(nullptr);
        //        if (do_per_step(step, 10))
        //        {
        //            fprintf(stderr, "\rDensity fitting map correlation
        //            coefficient: %g\n", df->cc);
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
    //    fprintf(stderr, "--- Forces took %g M cycles\n",
    //    0.000001*(gmx_cycles_read()-c1));

    /* Stop the density fitting cycle counter and add the computation-only
     * cycles to the force cycles for load balancing */
    cycles_comp = gmx_cycles_read() - cycles_comp;

    if (DOMAINDECOMP(cr) && wcycle)
    {
        dd_cycles_add(cr->dd, cycles_comp, ddCyclF);
    }
}

t_mapdata Densfit::referenceMapCopy() const { return parameters_->referenceMap(); };

const DensfitData &Densfit::parameters() const { return *parameters_; }

void Densfit::setParameters(const DensfitData &parameters)
{
    parameters_ = std::unique_ptr<DensfitData>(new DensfitData(parameters));
}

void Densfit::setReferenceMap(const gmx::t_mapdata &referenceMap)
{
    parameters_->setReferenceMap(referenceMap);
}


void Densfit::broadcast(const t_commrec *cr)
{
    block_bc(cr, *this);
    /* Now broadcast all remaining data. This is everything
     * not covered by the above statement, i.e. all data pointed to
     * by pointers.
     */
    parameters_->broadcast(cr);
    if (nweight_ > 0)
    {
        snew_bc(cr, weight_, parameters_->nAtoms());
        nblock_bc(cr, parameters_->nAtoms(), weight_);
    }
}

void Densfit::new_map_sim_from_ref(const t_mapdata &map_ref)
{
    map_sim_.datamode = static_cast<int>(MrcMetaData::MrcDataMode::float32);

    /* The simulated map copies the map_dim and grid settings from the reference
     * map. */
    for (int i = 0; i < 3; i++)
    {
        map_sim_.map_dim[i]    = map_ref.map_dim[i];
        map_sim_.grid[i]       = map_ref.grid[i];
        map_sim_.axes_order[i] = map_ref.axes_order[i];
    }
    for (int i = 0; i < 6; i++)
    {
        map_sim_.cell[i] = map_ref.cell[i];
    }
    map_sim_.title = strdup("Calculated map");

    allocate_density_grid(map_ref.map_dim, &map_sim_.vox);
}


void make_filename(const char *outf_name, int ndigit, int file_nr, char newname[STRLEN])
{
    char  fmt[128], nbuf[128];
    int   nd = 0, fnr;
    char *outf_base;
    char *extpos;


    /* Strip away the extension */
    outf_base = strdup(outf_name);
    extpos    = strrchr(outf_base, '.');
    if (extpos)
    {
        *extpos = '\0';
    }

    fnr = file_nr;
    do
    {
        fnr /= 10;
        nd++;
    }
    while (fnr > 0);

    if (ndigit > 0)
    {
        snprintf(fmt, 128, "%%0%dd", ndigit );
        snprintf(nbuf, 128, fmt, file_nr);
    }
    else
    {
        snprintf(nbuf, 128, "%d", file_nr);
    }

    if (extpos)
    {
        snprintf(newname, STRLEN, "%s%s.%s", outf_base, nbuf, extpos + 1);
    }
    else
    {
        snprintf(newname, STRLEN, "%s%s", outf_base, nbuf);
    }

    free(outf_base);
}


} // namespace gmx
