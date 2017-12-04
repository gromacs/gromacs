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


#include "localdensfitdata.h"

#include <cassert>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdlib/sim_util.h"


static const char *FitStr = "";                       //!< Used in gmx map output

namespace gmx
{

namespace
{
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

//! \brief For a position vector x return all direct 27 image positions given
//! the box
int get_images(rvec x, matrix box, rvec img[], ivec img_startindex[],
               ivec img_stop_index[], real grid_spacing, int voxrange,
               const std::array<int, DIM> &map_dim)
{

    const double OOgrid_spacing = 1.0 / grid_spacing;

    /* Loop over all possible 3x3x3 neighboring boxes */
    /* Minimum image convention? TODO: loop over NTRICIMG = 14, not 27 */
    int count = 0;
    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            for (int k = -1; k <= 1; k++)
            {
                rvec shiftvecX, shiftvecY, shiftvecZ;
                svmul(k, box[XX], shiftvecX);
                svmul(j, box[YY], shiftvecY);
                svmul(i, box[ZZ], shiftvecZ);

                rvec   image;
                copy_rvec(x, image);

                /* Move in x, y, and z directions */
                rvec_inc(image, shiftvecX);
                rvec_inc(image, shiftvecY);
                rvec_inc(image, shiftvecZ);

                /* Now check whether we will get any density from that image */
                int    icenter    = roundf(image[XX] * OOgrid_spacing);
                ivec   istart, istop;

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
                            for (int m = 0; m < DIM; m++)
                            {
                                istart[m] = std::max(istart[m], 0);  /* istart >= 0       */
                                istop[m]  =
                                    std::min(istop[m], map_dim[m]);  /* istop  <= map_dim */
                            }
                            copy_ivec(istart, img_startindex[count]);
                            copy_ivec(istop, img_stop_index[count]);

                            count++;

                        }  /* ZZ */
                    }      /* YY */
                }          /* XX */
            }
        }
    }

    return count;
}

}

/*! \brief Calculate one of the terms of the correlation coefficient
 *
 * Multiply each voxel of the map with the corresponding voxel of the other
 * map and sum up everything, i.e. compute
 *
 *  sum_ijk [ rhoA_ijk * rhoB_ijk ]
 *
 */
double calc_sum_rhoA_rhoB(const std::vector<float> &vecA,
                          const std::vector<float> &vecB)
{
    double sum_ijk = 0.;
    for (size_t i = 0; i < vecA.size(); i++)
    {
        sum_ijk += vecA[i] * vecB[i];
    }
    return sum_ijk;
}


//! \brief Allocate and initialize the arrays for assembled positions, forces
//! and atomic weights
LocalDensfitData::LocalDensfitData(LocalAtomSet localAtomSet, const DensfitData &parameters,
                                   rvec *x_densfit_whole) :
    localAtoms(localAtomSet)
{
    const auto &grid = parameters.referenceMap().getGrid();
    map_sim_.setGrid(grid.duplicate());

    /* Initialize the possibly time-dependent parameters */
    updateTimeDependentData(0., parameters);

    const auto nAtoms = parameters.fittingGroup().ind_.size();

    const RVec zeroRVecVector({0, 0, 0});
    x_assembled.resize(nAtoms, zeroRVecVector);
    f_loc.resize(nAtoms, zeroRVecVector);
    x_old.resize(nAtoms, zeroRVecVector);

    const IVec zeroIVecVector({0, 0, 0});
    x_shifts.resize(nAtoms, zeroIVecVector);
    extra_shifts.resize(nAtoms, zeroIVecVector);

    bUpdateShifts  = true;

    /* In mdruns we have to make the density fitting structure whole */
    if (x_densfit_whole != nullptr)
    {
        for (size_t i = 0; i < nAtoms; i++)
        {
            copy_rvec(x_densfit_whole[i], this->x_old[i]);
        }
    }
    /* For erf() and exp() performance, we want all values for one atom in one
     * linear array (not in three). Therefore, we compute the maximum number of
     * voxels to be computed for a single atom. How many voxels around each atom
     * we have to spread the density on is saved in voxrange (for each dimension
     * x, y, and z */
    voxrange =
        ceil(0.5 * sigma * parameters.dist() / grid.unitCell().basisVectorLength(0));
    vox_per_atom = 3 * (2 * voxrange + 1);
    const int ompthreads      = std::max(1, gmx_omp_nthreads_get(emntDefault));
    const int nVoxels         = map_sim_.getGrid().lattice().getNumLatticePoints();
    auto      voxBuf          = &openMPVoxelBuffer;
#pragma omp parallel for ordered num_threads(ompthreads) schedule(static) \
    shared(voxBuf) default(none)
    for (int i = 0; i < ompthreads; i++)
    {
#pragma omp ordered
        {
            voxBuf->emplace_back(nVoxels);
            int extra = 0;
#if GMX_SIMD_HAVE_REAL
            extra = GMX_SIMD_REAL_WIDTH;
#endif
            erfVector.emplace_back(AlignedRealVector(vox_per_atom+extra));
            expVector.emplace_back(AlignedRealVector(vox_per_atom+extra));
        }
        // + GMX_SIMD_REAL_WIDTH adds extra elements at the array end
        // so reading GMX_SIMD_REAL_WIDTH elements is always valid
        // even near the end of the vector.
    }

    sum_rho_exp2 = calc_sum_rhoA_rhoB(parameters.referenceMap(), parameters.referenceMap());
}

void LocalDensfitData::updateTimeDependentData(real time, const DensfitData &parameters)
{
    sigma = parameters.timeDependent().currentSigma(time);
    k     = parameters.timeDependent().currentForceConstant(time);
}

void LocalDensfitData::dump_x(int nodeid, gmx_int64_t step)
{
    FILE *fpout = nullptr;
    char  fn[STRLEN];
    char  buf[256];

    gmx_step_str(step, buf);
    sprintf(fn, "x_step%s_node%d.txt", buf, nodeid);
    fpout = fopen(fn, "w");

    for (size_t i = 0; i < localAtoms.numAtomsGlobal(); i++)
    {
        fprintf(fpout, "%4lu %12.3e %12.3e %12.3e\n", i, x_assembled[i][XX],
                x_assembled[i][YY], x_assembled[i][ZZ]);
    }

    fclose(fpout);
}

void LocalDensfitData::dump_f(const char *fn, const t_commrec *cr)
{
    rvec *f;
    FILE *fpout = NULL;

    /* Assemble the forces */
    snew(f, localAtoms.numAtomsGlobal());
    /* Each node copies its local forces into the collective array */
    for (size_t i = 0; i < localAtoms.numAtomsLocal(); i++)
    {
        copy_rvec(f_loc[i], f[localAtoms.collectiveIndex()[i]]);
    }

    /* Reduce */
    if (cr && PAR(cr))
    {
        gmx_sum(3 * localAtoms.numAtomsGlobal(), f[0], cr);
    }

    /* Master outputs */
    if (cr == nullptr || MASTER(cr))
    {
        if (fn)
        {
            fpout = gmx_fio_fopen(fn, "w");
        }

        fprintf(stderr, "%sDumping forces to file %s.\n", FitStr, fn ? fn : "");

        for (size_t i = 0; i < localAtoms.numAtomsGlobal(); i++)
        {
            fprintf(fpout, "#     f[%5lu]={%12.5e, %12.5e, %12.5e}\n", i, f[i][XX],
                    f[i][YY], f[i][ZZ]);
        }
        if (fn)
        {
            gmx_fio_fclose(fpout);
        }
    }
    sfree(f);
}

void LocalDensfitData::spread_atoms_low(int istart, int natoms, matrix box)
{
    //    gmx_cycles_t c1 = gmx_cycles_read();
    const int    nth        = std::max(1, gmx_omp_nthreads_get(emntDefault));
/* Zero out the grid */
#pragma omp parallel for num_threads(nth) schedule(static) \
    default(none)
    for (int thread = 0; thread < nth; thread++)
    {
        std::fill(std::begin(openMPVoxelBuffer[thread]), std::end(openMPVoxelBuffer[thread]), 0.);
    }
    const auto   x            = x_assembled.data() + istart;
    const int    nx           = map_sim_.getGrid().lattice().extend()[XX];
    const int    ny           = map_sim_.getGrid().lattice().extend()[YY];
    const int    nz           = map_sim_.getGrid().lattice().extend()[ZZ];
    const int    ngrid        = nx * ny * nz;
    const double sigma2       = sigma * sigma;
    const double c_fac        = sqrt(1.5 / sigma2);
    const auto   dfp          = this;
    const auto   mapsim_dim   = map_sim_.getGrid().lattice().extend();
    const auto   spacing      = (map_sim_.getGrid().unitCell().basisVectorLength(0)+map_sim_.getGrid().unitCell().basisVectorLength(1)+map_sim_.getGrid().unitCell().basisVectorLength(2))/3.0;
    const double V_vox        = pow(spacing, 3);
    const double prefactor    = pow(M_PI * sigma2 / 6.0, 1.5) / V_vox;

/* Spread the local atoms on the grid, therefore loop over all atoms of the
 * density fitting group */
#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(natoms, box) default(none)
    for (int n = 0; n < natoms; n++)
    {
        int th = gmx_omp_get_thread_num(); /* Thread ID */

        /* Determine the grid points around the atom and its images for which we
         * need to calculate stuff */
        rvec      img[27];      /* Stores PBC images of a position x as well as x */
        ivec      istart[27];   /* Store for each image the part of the grid ...  */
        ivec      istop[27];    /* ... where stuff has to be calculated           */

        const int num_images = get_images(x[n], box, img, istart, istop, spacing,
                                          dfp->voxrange, mapsim_dim);

        /* Loop over all periodic images of the atom */
        for (int image = 0; image < num_images; image++)
        {
            rvec pos;
            copy_rvec(img[image], pos); // periodic image of the position
            ivec elements;
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
            real  *aarr[3];

            aarr[XX] = &dfp->erfVector[th][0];
            aarr[YY] = &dfp->erfVector[th][0 + elements[XX]];
            aarr[ZZ] = &dfp->erfVector[th][0 + elements[XX] + elements[YY]];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in a temporary array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            int i = 0;
            for (int ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                const real left          = (ix - 0.5) * spacing - pos[XX];
                aarr[XX][i++] = c_fac * left;
            }

            i = 0;
            for (int iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                const real left          = (iy - 0.5) * spacing - pos[YY];
                aarr[YY][i++] = c_fac * left;
            }

            i = 0;
            for (int iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                const real left          = (iz - 0.5) * spacing - pos[ZZ];
                aarr[ZZ][i++] = c_fac * left;
            }
            gmx_erf_vector(elements[XX] + elements[YY] + elements[ZZ], aarr[XX]);

            /* Transform (in place) into a new array that directly stores the
             * differences of
             * the error functions, i.e. the area under the Gaussian curve */
            for (int ix = 0; ix < elements[XX] - 1; ix++)
            {
                aarr[XX][ix] = aarr[XX][ix + 1] - aarr[XX][ix];
            }
            for (int iy = 0; iy < elements[YY] - 1; iy++)
            {
                aarr[YY][iy] = aarr[YY][iy + 1] - aarr[YY][iy];
            }
            for (int iz = 0; iz < elements[ZZ] - 1; iz++)
            {
                aarr[ZZ][iz] = aarr[ZZ][iz + 1] - aarr[ZZ][iz];
            }

            /* Spread the density for this atom */
            for (int iz = istart[image][ZZ]; iz < istop[image][ZZ]; iz++)
            {
                const real az = aarr[ZZ][iz - istart[image][ZZ]];

                for (int iy = istart[image][YY]; iy < istop[image][YY]; iy++)
                {
                    const real ay   = aarr[YY][iy - istart[image][YY]];
                    const real ayaz = ay * az;
                    float    * ptr  = &(openMPVoxelBuffer[th][nx * (ny * iz + iy) +
                                                              istart[image][XX]]);
                    for (int ix = istart[image][XX]; ix < istop[image][XX]; ix++)
                    {
                        const real ax = aarr[XX][ix - istart[image][XX]];

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

    auto map_sim_p = &map_sim_;
#pragma omp parallel for num_threads(nth) schedule(static) shared(map_sim_p) default(none)
    for (int i = 0; i < ngrid; i++)
    {
        for (int j = 0; j < nth; j++)
        {
            (*map_sim_p)[i] += openMPVoxelBuffer[j][i];
        }
        (*map_sim_p)[i] *= prefactor;
    }
}

void LocalDensfitData::do_forces(matrix box, const GridDataReal3D &referenceMap)
{

    const std::vector<float> &sim_dens_1d = map_sim_;

    if (localAtoms.numAtomsLocal() <= 0)
    {
        return;
    }

    // cppcheck-suppress unreadVariable
    int gmx_unused nth =
        std::max(1, gmx_omp_nthreads_get(emntDefault)); // number of threads

/* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        clear_rvec(f_loc[l]);
    }

    /* Calculate various prefactors */
    const double sum_rho_sim2 =
        calc_sum_rhoA_rhoB(map_sim_, map_sim_);
    const double sum_rho_exp_rho_sim =
        calc_sum_rhoA_rhoB(referenceMap, map_sim_);

    const double term1_prefactor = k / sqrt(sum_rho_exp2 * sum_rho_sim2);
    const double term2_prefactor = -k * sum_rho_exp_rho_sim /
        (sqrt(sum_rho_exp2) * pow(sum_rho_sim2, 1.5));
    const auto   spacing          = (map_sim_.getGrid().unitCell().basisVectorLength(0)+map_sim_.getGrid().unitCell().basisVectorLength(1)+map_sim_.getGrid().unitCell().basisVectorLength(2))/3.0;
    const double V_vox            = spacing * spacing * spacing;
    const double prefactor        = -M_PI * sigma * sigma / (6.0 * V_vox);
    const double OOsigma2_fac     = 1.5 / (sigma * sigma);
    const double sqrtOOsigma2_fac = sqrtf(OOsigma2_fac);

    const int    nx = map_sim_.getGrid().lattice().extend()[XX];
    const int    ny = map_sim_.getGrid().lattice().extend()[YY];

    const auto   mapsim_dim = map_sim_.getGrid().lattice().extend();
#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(sim_dens_1d, referenceMap, box) default(none)
    /* Loop over all atoms of the density fitting group */
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        real                     *aarr[3], *garr[3];

        rvec                      img[27];    /* Stores PBC images of a position x as well as x */
        ivec                      istart[27]; /* Store for each image the part of the grid ...  */
        ivec                      istop[27];  /* ... where stuff has to be calculated           */
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
        int num_images =
            get_images(x_assembled[localAtoms.collectiveIndex()[l]], box, img, istart, istop,
                       spacing, voxrange, mapsim_dim);

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
            aarr[XX] = &erfVector[th][0];
            aarr[YY] = &erfVector[th][0 + elements[XX]];
            aarr[ZZ] = &erfVector[th][0 + elements[XX] + elements[YY]];

            /* Same for the temporary array that keeps the exp values */
            garr[XX] = &expVector[th][0];
            garr[YY] = &expVector[th][0 + elements[XX]];
            garr[ZZ] = &expVector[th][0 + elements[XX] + elements[YY]];

            /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
             * boundaries in the tmpgrid[] array. This way we can compute the erf's
             * for this atom all in one (vector) call. */
            int i = 0;
            for (int ix = istart[image][XX]; ix < istop[image][XX] + 1; ix++)
            {
                double left = (ix - 0.5) * spacing - pos[XX];
                garr[XX][i]   = -OOsigma2_fac * left * left;
                aarr[XX][i++] = sqrtOOsigma2_fac * left;
            }

            i = 0;
            for (int iy = istart[image][YY]; iy < istop[image][YY] + 1; iy++)
            {
                double left = (iy - 0.5) * spacing - pos[YY];
                garr[YY][i]   = -OOsigma2_fac * left * left;
                aarr[YY][i++] = sqrtOOsigma2_fac * left;
            }

            i = 0;
            for (int iz = istart[image][ZZ]; iz < istop[image][ZZ] + 1; iz++)
            {
                double left = (iz - 0.5) * spacing - pos[ZZ];
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

                    auto ptr_ref = referenceMap.begin() + ind;
                    auto ptr_sim = sim_dens_1d.begin() + ind;

                    /* Calculate d/dx_l rho_sim_ijk */
                    for (int ix = 0; ix < elements[XX] - 1; ix++)
                    {
                        const double ref_dens = (double)*ptr_ref;
                        const double sim_dens = (double)*ptr_sim;

                        const real   erfx = aarr[XX][ix];
                        const real   expx = garr[XX][ix];

                        rvec         drho_dxl;

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

        f_loc[l][XX] = force[XX];
        f_loc[l][YY] = force[YY];
        f_loc[l][ZZ] = force[ZZ];

    }   /* end of loop over atoms */
}

const GridDataReal3D &LocalDensfitData::simulatedMap() const
{
    return map_sim_;
}


std::string LocalDensfitData::infoString() const
{
    char buffer[STRLEN];
    sprintf(buffer, "%12g%12.5e%12.5e%12.7f%12.5e", 0., k, sigma, cc, k * (1.0 - cc));
    return buffer;
}


void LocalDensfitData::add_forces(rvec * f) const
{
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        /* Get the right index of the local force, since typically not all local
         * atoms are subject to density fitting forces */
        rvec_inc(f[localAtoms.localIndex()[l]], f_loc[l]);
    }
}

bool LocalDensfitData::afterDomainDecomposition() const
{
    return bUpdateShifts;
}

void LocalDensfitData::communicate( PaddedArrayRef<RVec> x, const t_commrec * cr, matrix box, int ePBC)
{
    communicate_group_positions(
            cr, as_vec_array(x_assembled.data()), as_vec_array(x_shifts.data()), as_vec_array(extra_shifts.data()),
            bUpdateShifts, as_rvec_array(x.data()), localAtoms.numAtomsGlobal(),
            localAtoms.numAtomsLocal(), localAtoms.localIndex().data(), localAtoms.collectiveIndex().data(), as_vec_array(x_old.data()), box);

/* If bUpdateShifts was TRUE then the shifts have just been updated in
 * communicate_group_positions. We do not need to update the shifts until
 * the next NS step */
    bUpdateShifts = false;

/* Put all atoms in the box (TODO: if we do that, we do not need to construct
 * a whole DF group before with communicate_group_positions!)
 */
    if (ePBC != epbcNONE)
    {
        auto xAssembledArrayRef = arrayRefFromArray(
                    x_assembled.data(), localAtoms.numAtomsGlobal());
        put_atoms_in_box_omp(ePBC, box, xAssembledArrayRef);
    }

}

void LocalDensfitData::triggerShiftUpdate()
{
    bUpdateShifts = true;
}


void LocalDensfitData::spreadAtoms(const t_commrec * cr, matrix box)
{

    if (cr == nullptr)
    {
        spread_atoms_low(0, localAtoms.numAtomsGlobal(), box);
        return;
    }

    int  istart;
    int  nspread;
    int  nnodes;

    if (PAR(cr))
    {
        nnodes = cr->nnodes - cr->npmenodes;

        nspread = ceil((real) localAtoms.numAtomsGlobal() / nnodes);
        istart  = cr->nodeid * nspread;

        if ((nnodes - 1) == cr->nodeid)
        {
            nspread = localAtoms.numAtomsGlobal()  - nspread * (nnodes - 1);
        }
    }
    else
    {
        istart  = 0;
        nspread = localAtoms.numAtomsGlobal();
    }

    assert(istart >= 0);
    assert(nspread >= 0);

    spread_atoms_low(istart, nspread, box);


}

void LocalDensfitData::sumSimulatedMap(const t_commrec * cr)
{
    gmx_sumf(map_sim_.size(), map_sim_.data(), cr);
}


void LocalDensfitData::calculateCorrelationCoefficient(const GridDataReal3D &referenceMap)
{
    cc = calc_sum_rhoA_rhoB(referenceMap, map_sim_) / sqrt(sum_rho_exp2 *
                                                           calc_sum_rhoA_rhoB(map_sim_, map_sim_));
}

}
