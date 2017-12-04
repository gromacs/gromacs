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
 * Declares data structure and utilities for density fitting
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "maputil.h"
#include "densityfitting/densfit.h"

#include <assert.h>
#include <array>
#include <stdio.h>

#include <cstring>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

#include "gromacs/fileio/mrcmetadata.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/math/griddata/grid.h"
#include "gromacs/mdtypes/inputrec.h"


#ifndef _maputil_h
#include "maputil.h"
#endif

static const char* FitStr      = "";                    //!< Used in gmx map output



namespace gmx
{

/*!\brief Converts mapdata structs to GridDataReal3D class and vice versa.
 *
 * This class is a temporary clutch until mapdata and GridDataReal3D functionality are merged .
 */
t_mapdata gridDataToMapData(const GridDataReal3D &inputGrid)
{
    t_mapdata   mapdata;
    MrcMetaData meta;
    meta.setEMDBDefaults();
    meta.fromGrid(inputGrid.getGrid());
    std::copy(std::begin(meta.crs_to_xyz), std::end(meta.crs_to_xyz), std::begin(mapdata.axes_order));
    std::copy(std::begin(meta.cell_length.as_vec()), std::end(meta.cell_length.as_vec()), std::begin(mapdata.cell));
    std::copy(std::begin(meta.cell_angles.as_vec()), std::end(meta.cell_angles.as_vec()), std::begin(mapdata.cell)+3);
    mapdata.datamode = meta.mrc_data_mode;
    std::copy(std::begin( meta.extend), std::end( meta.extend), std::begin(mapdata.grid));
    std::copy(std::begin( meta.extend), std::end( meta.extend), std::begin(mapdata.map_dim));
    mapdata.max  = meta.max_value;
    mapdata.mean = meta.mean_value;
    mapdata.min  = meta.min_value;
    std::copy(std::begin( meta.crs_start), std::end( meta.crs_start), std::begin(mapdata.origin));
    mapdata.rms = meta.rms_value;
    std::copy(std::begin( meta.skew_matrix), std::end( meta.skew_matrix), std::begin(mapdata.skew_mat));
    std::copy(std::begin(meta.skew_translation.as_vec()), std::end(meta.skew_translation.as_vec()), std::begin(mapdata.skew_trans));
    mapdata.spacegroup = meta.space_group;
    mapdata.title      = gmx_strdup(meta.labels[0].c_str());
    mapdata.vox        = inputGrid;
    return mapdata;
}

GridDataReal3D mapDataToGridData(const t_mapdata &mapdata)
{
    MrcMetaData meta;
    meta.setEMDBDefaults();
    std::copy(std::begin(mapdata.axes_order), std::end(mapdata.axes_order), std::begin(meta.crs_to_xyz));

    std::copy(std::begin(mapdata.cell), std::begin(mapdata.cell)+3, std::begin(meta.cell_length.as_vec()));
    std::copy(std::begin(mapdata.cell)+3, std::end(mapdata.cell), std::begin(meta.cell_angles.as_vec()));
    meta.mrc_data_mode = mapdata.datamode;
    std::copy(std::begin(mapdata.grid), std::end(mapdata.grid), std::begin(meta.extend));
    // std::copy(std::begin(mapdata.grid), std::end(meta.extend),std::begin(mapdata.map_dim));
    meta.max_value  = mapdata.max;
    meta.mean_value = mapdata.mean;
    meta.min_value  = mapdata.min;
    std::copy(std::begin(mapdata.origin), std::end(mapdata.origin), std::begin(meta.crs_start));
    meta.rms_value = mapdata.rms;
    std::copy(std::begin(mapdata.skew_mat), std::end(mapdata.skew_mat), std::begin(meta.skew_matrix));
    std::copy(std::begin(mapdata.skew_trans), std::end(mapdata.skew_trans), std::begin(meta.skew_translation.as_vec()));
    meta.space_group = mapdata.spacegroup;

    meta.labels[0] = std::string(mapdata.title);
    GridDataReal3D resultingGridData(meta.toGrid());
    std::copy(std::begin(mapdata.vox), std::end(mapdata.vox), std::begin(resultingGridData));
    return resultingGridData;
}

extern void allocate_density_grid(const std::array<int, DIM> &nValues, std::vector<float> *vox)
{
    int ompthreads = std::max(1, gmx_omp_nthreads_get(emntDefault));
    vox->resize(ompthreads * nValues[XX] * nValues[YY] * nValues[ZZ]);
}


extern void new_map_sim_from_ref(
        Densfit         *densfit,
        const t_mapdata &map_ref)
{
    int i;

    densfit->map_sim.datamode = static_cast<int>(MrcMetaData::MrcDataMode::float32);

    /* The simulated map copies the map_dim and grid settings from the reference map. */
    for (i = 0; i < 3; i++)
    {
        densfit->map_sim.map_dim   [i] = map_ref.map_dim   [i];
        densfit->map_sim.grid      [i] = map_ref.grid      [i];
        densfit->map_sim.axes_order[i] = map_ref.axes_order[i];
    }
    for (i = 0; i < 6; i++)
    {
        densfit->map_sim.cell[i] = map_ref.cell[i];
    }
    densfit->map_sim.title = strdup("Calculated map");

    allocate_density_grid(map_ref.map_dim, &densfit->map_sim.vox);
}


extern void translate_pdb(matrix box, rvec x[], int natoms, gmx_bool bTranslate)
{
    int  i;
    real x_max = -GMX_REAL_MAX, y_max = -GMX_REAL_MAX, z_max = -GMX_REAL_MAX;
    real x_min = GMX_REAL_MAX, y_min = GMX_REAL_MAX, z_min = GMX_REAL_MAX;


    if (   (box[XX][YY] != 0.0)
           || (box[XX][ZZ] != 0.0)
           || (box[YY][XX] != 0.0)
           || (box[YY][ZZ] != 0.0)
           || (box[ZZ][XX] != 0.0)
           || (box[ZZ][YY] != 0.0) )
    {
        gmx_fatal(FARGS, "All off-diagonal box elements must be zero.\n");
    }

    for (i = 0; i < natoms; i++)
    {
        x_min = std::min(x_min, x[i][XX]);
        y_min = std::min(y_min, x[i][YY]);
        z_min = std::min(z_min, x[i][ZZ]);

        x_max = std::max(x_max, x[i][XX]);
        y_max = std::max(y_max, x[i][YY]);
        z_max = std::max(z_max, x[i][ZZ]);
    }
    fprintf(stdout, "Positions were in range %g <= x <= %g\n", x_min, x_max);
    fprintf(stdout, "Positions were in range %g <= y <= %g\n", y_min, y_max);
    fprintf(stdout, "Positions were in range %g <= z <= %g\n", z_min, z_max);

    /* Optionally translate all positions to first quadrant */
    if (bTranslate)
    {
        for (i = 0; i < natoms; i++)
        {
            x[i][XX] -= x_min;
            x[i][YY] -= y_min;
            x[i][ZZ] -= z_min;
        }

        box[XX][XX] = x_max - x_min;
        box[YY][YY] = y_max - y_min;
        box[ZZ][ZZ] = z_max - z_min;
    }
    else
    {
        box[XX][XX] = x_max;
        box[YY][YY] = y_max;
        box[ZZ][ZZ] = z_max;
    }

    fprintf(stdout, "Assigning rectangular box xmax=%f ymax=%f zmax=%f nm\n",
            box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
}

extern real get_map_spacing(const t_mapdata &map, FILE *log)
{
    rvec spacing;
    real av_spacing;
    ivec points;
    int  i;

    /* Number of grid points in each dimension */
    for (i = 0; i < DIM; i++)
    {
        points [i] = map.map_dim[i];
        spacing[i] = A2NM * map.cell[i] / points[i];
    }

    if  (    ( fabs(spacing[XX] - spacing[YY]) > 0.01 )
             || ( fabs(spacing[XX] - spacing[ZZ]) > 0.01 ) )
    {
        fprintf(stderr, "WARING: Reference grid spacing (%g %g %g) should be identical in XYZ !",
                spacing[XX], spacing[YY], spacing[ZZ]);
    }

    av_spacing = (spacing[XX]+spacing[YY]+spacing[ZZ]) / 3.0;

    if (log)
    {
        fprintf(log, "%sGrid spacing (nm) is %g in x-, %g in y-, and %g in z-dimension (average %g).\n",
                FitStr, spacing[XX], spacing[YY], spacing[ZZ], av_spacing);
    }

    /* Return the grid spacing */
    return av_spacing;
}


//! \brief Get the length of the intersection of the two segments a and b
static inline double get_intersecting_length(double x11, double x12, double x21, double x22)
{
    return std::max(0.0, std::min(x12, x22) - std::max(x11, x21));
}

//! \brief Get the volume of the intersection of new_vox and ref_vox
static double get_intersecting_volume(
        dvec v1_start,
        dvec v1_stop,
        dvec v2_start,
        dvec v2_stop)
{
    int  d;
    dvec intersection;


    for (d = 0; d < DIM; d++)
    {
        intersection[d] = get_intersecting_length(v1_start[d], v1_stop[d], v2_start[d], v2_stop[d]);
    }

    return intersection[XX] * intersection[YY] * intersection[ZZ];
}


/*! \brief Calculates density for a rescaled map
 *
 * Return the average (integrated) density of the new, scaled voxel at ix, iy, iz which
 * is made of several smaller voxels of the original density map
 *
 * Variable naming convention: i for map_ref indices,
 *                             j for map_scaled indices
 */
static double get_rho_av(const t_mapdata &map_ref, const t_mapdata &map_new, ivec j)
{
    int    d, ix, iy, iz, ind;
    ivec   n_ref, n_new;
    dvec   mapsize;
    double rho;                               /* Average / integrated density in Angstrom^(-3) */
    double vol;
    dvec   ref_vox_size, new_vox_size;        /* Angstrom */
    dvec   new_vox_start, new_vox_stop;       /* Angstrom */
    dvec   ref_vox_start, ref_vox_stop;       /* Angstrom */
    ivec   ref_vox_startind, ref_vox_stopind; /* First and last ref voxel index per dimension
                                                 that contributes to new voxel density */
    float  refdensity;


    for (d = 0; d < DIM; d++)
    {
        /* Number of voxels per dimension of the reference and scaled map */
        n_ref[d]         = map_ref.map_dim[d];
        n_new[d]         = map_new.map_dim[d];

        /* Size (Angstrom) of the maps (both have the same size!) */
        mapsize[d]       = map_ref.cell[d];

        /* Voxel size per dimension of reference and scaled map */
        ref_vox_size[d]  = mapsize[d] / n_ref[d];
        new_vox_size[d]  = mapsize[d] / n_new[d];

        /* Start and stop position (Angstrom) of the new voxel */
        new_vox_start[d] = new_vox_size[d] *  j[d];
        new_vox_stop [d] = new_vox_size[d] * (j[d] + 1);

        /* Get the first and the last voxel index of the ref map which
         * contributes to the voxel density at this point of the scaled map */
        ref_vox_startind[d] = floor(new_vox_start[d] / ref_vox_size[d]);
        ref_vox_stopind [d] = floor(new_vox_stop [d] / ref_vox_size[d]);

        ref_vox_stopind[d] = std::min(ref_vox_stopind[d], n_ref[d] - 1);
    }

    /* Now loop over all small subvoxels, get the subvoxel volume plus density value,
     * and add to resulting voxel density */
    rho = 0.0;
    for (iz = ref_vox_startind[ZZ]; iz <= ref_vox_stopind[ZZ]; iz++)
    {
        for (iy = ref_vox_startind[YY]; iy <= ref_vox_stopind[YY]; iy++)
        {
            for (ix = ref_vox_startind[XX]; ix <= ref_vox_stopind[XX]; ix++)
            {
                ind        = n_ref[XX]*n_ref[YY]*iz + n_ref[XX]*iy + ix;
                refdensity = map_ref.vox[ind];

                if (refdensity != 0)
                {
                    /* Start and stop position (Angstrom) of the reference voxel with which to intersect */
                    ref_vox_start[XX] = ref_vox_size[XX] *  ix;
                    ref_vox_stop [XX] = ref_vox_size[XX] * (ix + 1);
                    ref_vox_start[YY] = ref_vox_size[YY] *  iy;
                    ref_vox_stop [YY] = ref_vox_size[YY] * (iy + 1);
                    ref_vox_start[ZZ] = ref_vox_size[ZZ] *  iz;
                    ref_vox_stop [ZZ] = ref_vox_size[ZZ] * (iz + 1);

                    /* Get the intersecting volume of the ref voxel at ix,iy,iz with
                     * the new voxel */
                    vol  =  get_intersecting_volume(new_vox_start, new_vox_stop, ref_vox_start, ref_vox_stop);
                    rho += refdensity * vol;
                }
            }
        }
    }

    return rho;
}


extern void make_positive(t_mapdata *map_ref)
{
    float     minimum = GMX_FLOAT_MAX;
    int       nx, ny, nz, i, count;


    nx  = map_ref->map_dim[XX];
    ny  = map_ref->map_dim[YY];
    nz  = map_ref->map_dim[ZZ];

    /* Determine the minimum value of the voxel density */
    count = 0;
    for (i = 0; i < nx*ny*nz; i++)
    {
        minimum = std::min(minimum, map_ref->vox[count] );
        count++;
    }

    /* Subtract the minimum */
    fprintf(stdout, "Subtracting the minimum vale of %g from all voxels.\n", minimum);
    count = 0;
    for (i = 0; i < nx*ny*nz; i++)
    {
        map_ref->vox[count] -= minimum;
        count++;
    }

    map_ref->min  = 0.0;
    map_ref->max  = 0.0;
    map_ref->mean = 0.0;
    map_ref->rms  = 0.0;
}


//TODO: scale sigma by the same factor?
//TODO: allow only values > 0?
extern t_mapdata rescale_map(const t_mapdata &map_ref, real scale)
{
    int                   i, ix, iy, iz;
    ivec                  ind;
    t_mapdata             map_scaled;
    int                   nx, ny, nz;
    int                   mx, my, mz;
    std::vector<float>    scaledDens;
    double                rho_av;

    map_scaled.datamode = static_cast<int>(MrcMetaData::MrcDataMode::float32);

    for (i = 0; i < 3; i++)
    {
        map_scaled.map_dim[i] = round(map_ref.map_dim[i] * scale);
        map_scaled.grid   [i] = round(map_ref.grid   [i] * scale);
    }
    for (i = 0; i < 6; i++)
    {
        map_scaled.cell[i] = map_ref.cell[i];
    }
    map_scaled.title = strdup("Scaled map");

    allocate_density_grid(map_scaled.map_dim, &map_scaled.vox);

    nx  = map_scaled.map_dim[XX];
    ny  = map_scaled.map_dim[YY];
    nz  = map_scaled.map_dim[ZZ];

    mx  = map_ref.map_dim[XX];
    my  = map_ref.map_dim[YY];
    mz  = map_ref.map_dim[ZZ];

    double volratio = (double)(nx*ny*nz) / (double)(mx*my*mz);
    fprintf(stderr, "volratio = %g\n", volratio);

    scaledDens = map_scaled.vox;
    /* Determine the average density for the new voxel per dimension */
    for (iz = 0; iz < nz; iz++)
    {
        ind[ZZ] = iz;
        for (iy = 0; iy < ny; iy++)
        {
            ind[YY] = iy;
            for (ix = 0; ix < nx; ix++)
            {
                /* Now determine the average density (in Angstrom^(-3)) of the original map at the ix,iy,iz voxel */
                ind[XX] = ix;
                rho_av  = get_rho_av(map_ref, map_scaled, ind);

                /* Assign average density of original map to new voxel */
                scaledDens[nx*ny*iz + nx*iy + ix] = rho_av * volratio;
            }
        }
    }

    return map_scaled;
}




extern void make_filename(const char *outf_name, int ndigit, int file_nr, char newname[STRLEN])
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


//static void checkPerformance(t_inputrec *ir)
//{
//    int i, j, n = 100000;
//    real *values = NULL;
//    real res;
//    gmx_rng_t rng;
//
//
//    fprintf(stderr, "\n=== Testing the performance of the erf and exp functions ===\n");
//    rng = gmx_rng_init(ir->ld_seed);
//
//
//    /* Fill the array with random values */
//    snew(values, n);
//    for (i = 0; i < n; i++)
//    {
//        values[i] = 10*(gmx_rng_uniform_real(rng) - 0.5);
//    }
//
//    real *x = NULL;
//    snew(x, 2);
//    gmx_cycles_t   c1, c2;
//
//    for (i = 1; i < 10; i++)
//    {
//        x[0] = values[0];
//        x[1] = values[1];
//
//        c1 = gmx_cycles_read();
//        res = erf(x[1]) - erf(x[0]);
//        c2 = gmx_cycles_read();
//        fprintf(stderr, "Result: %15.10f   %lld cycles  normal erf()\n", res, c2 - c1);
//
//        c1 = gmx_cycles_read();
//        res = std::erf(x[1]) - std::erf(x[0]);
//        c2 = gmx_cycles_read();
//        fprintf(stderr, "Result: %15.10f   %lld cycles  gmx_erf()\n", res, c2 - c1);
//
//        gmx_erf_vector(2, x);
//        fprintf(stderr, "Result: %15.10f  vector\n\n", x[1] - x[0]);
//    }
//
//
//    for (j = 0; j < 2; j++)
//    {
//        for (i = 1; i < n; i *= 2)
//        {
//            gmx_erf_vector(i, values);
//        }
//        fprintf(stderr, "\n");
//
//        for (i = 1; i < n; i *= 2)
//        {
//            gmx_exp_vector(i, values);
//        }
//        fprintf(stderr, "\n");
//    }
//    sfree(values);
//}

//static void dump_local_indices(Densfit *densfit, t_commrec *cr)
//{
//    FILE          *fpout = NULL;
//    char           fn[STRLEN];
//    int            i;
//    t_gmx_densfit *df;
//
//
//    df = densfit->df;
//    sprintf(fn, "ind_node%d.txt", cr->nodeid);
//    fpout = fopen(fn, "w");
//
//    for (i = 0; i < densfit->nat; i++)
//    {
//        fprintf(fpout, "%4d %4d %12.3e %12.3e %12.3e\n", i, df->c_ind[i],
//                df->x_assembled[i][XX],
//                df->x_assembled[i][YY],
//                df->x_assembled[i][ZZ]);
//    }
//
//    fclose(fpout);
//}

extern real add_densfit_forces(t_inputrec *ir, rvec *f, gmx_unused t_commrec *cr, gmx_int64_t step, real time)
{
    return ir->densfit->add_forces(f, cr, step, time);
}

extern void do_densfit(
        real                           t,
        gmx_int64_t                    step,
        gmx_bool                       bOutputMap,
        t_inputrec                    *ir,
        t_commrec                     *cr,
        gmx::PaddedArrayRef<gmx::RVec> x,
        matrix                         box,
        gmx_wallcycle                 *wcycle)
{
    ir->densfit->do_densfit(t, step, bOutputMap, ir, cr, x, box, wcycle);
}

} // namespace gmx
