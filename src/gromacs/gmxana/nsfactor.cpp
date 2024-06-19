/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "nsfactor.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdio>

#include <filesystem>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

void check_binwidth(real binwidth)
{
    real smallest_bin = 0.1;
    if (binwidth < smallest_bin)
    {
        gmx_fatal(FARGS,
                  "Binwidth shouldn't be smaller then smallest bond length (H-H bond ~0.1nm) in a "
                  "box");
    }
}

void check_mcover(real mcover)
{
    if (mcover > 1.0)
    {
        gmx_fatal(FARGS, "mcover should be -1 or (0,1]");
    }
    else if ((mcover < 0) && (mcover != -1))
    {
        gmx_fatal(FARGS, "mcover should be -1 or (0,1]");
    }
    else
    {
        return;
    }
}

void normalize_probability(int n, double* a)
{
    int    i;
    double norm = 0.0;
    for (i = 0; i < n; i++)
    {
        norm += a[i];
    }
    for (i = 0; i < n; i++)
    {
        a[i] /= norm;
    }
}

gmx_neutron_atomic_structurefactors_t* gmx_neutronstructurefactors_init(const char* datfn)
{
    /* read nsfactor.dat */
    char                                   line[STRLEN];
    int                                    nralloc = 10;
    int                                    n, p;
    int                                    i, line_no;
    char                                   atomnm[8];
    double                                 slength;
    gmx_neutron_atomic_structurefactors_t* gnsf;

    gmx::FilePtr fp = gmx::openLibraryFile(datfn);
    line_no         = 0;
    /* allocate memory for structure */
    snew(gnsf, nralloc);
    snew(gnsf->atomnm, nralloc);
    snew(gnsf->p, nralloc);
    snew(gnsf->n, nralloc);
    snew(gnsf->slength, nralloc);

    gnsf->nratoms = line_no;

    while (get_a_line(fp.get(), line, STRLEN))
    {
        i = line_no;
        if (sscanf(line, "%s %d %d %lf", atomnm, &p, &n, &slength) == 4)
        {
            gnsf->atomnm[i]  = gmx_strdup(atomnm);
            gnsf->n[i]       = n;
            gnsf->p[i]       = p;
            gnsf->slength[i] = slength;
            line_no++;
            gnsf->nratoms = line_no;
            if (line_no == nralloc)
            {
                nralloc++;
                srenew(gnsf->atomnm, nralloc);
                srenew(gnsf->p, nralloc);
                srenew(gnsf->n, nralloc);
                srenew(gnsf->slength, nralloc);
            }
        }
        else
        {
            fprintf(stderr, "WARNING: Error in file %s at line %d ignored\n", datfn, line_no);
        }
    }
    srenew(gnsf->atomnm, gnsf->nratoms);
    srenew(gnsf->p, gnsf->nratoms);
    srenew(gnsf->n, gnsf->nratoms);
    srenew(gnsf->slength, gnsf->nratoms);

    return gnsf;
}

gmx_sans_t* gmx_sans_init(const t_topology* top, gmx_neutron_atomic_structurefactors_t* gnsf)
{
    gmx_sans_t* gsans = nullptr;
    int         i, j;
    /* Try to assing scattering length from nsfactor.dat */
    snew(gsans, 1);
    snew(gsans->slength, top->atoms.nr);
    /* copy topology data */
    gsans->top = top;
    for (i = 0; i < top->atoms.nr; i++)
    {
        for (j = 0; j < gnsf->nratoms; j++)
        {
            if (top->atoms.atom[i].atomnumber == gnsf->p[j])
            {
                /* we need special case for H and D */
                if (top->atoms.atom[i].atomnumber == 1)
                {
                    if (gmx_within_tol(top->atoms.atom[i].m, 1.008000, 1e-3))
                    {
                        gsans->slength[i] = gnsf->slength[0];
                    }
                    else
                    {
                        gsans->slength[i] = gnsf->slength[1];
                    }
                }
                else
                {
                    gsans->slength[i] = gnsf->slength[j];
                }
            }
        }
    }

    return gsans;
}

gmx_radial_distribution_histogram_t* calc_radial_distribution_histogram(gmx_sans_t*  gsans,
                                                                        rvec*        x,
                                                                        matrix       box,
                                                                        const int*   index,
                                                                        int          isize,
                                                                        double       binwidth,
                                                                        gmx_bool     bMC,
                                                                        gmx_bool     bNORM,
                                                                        real         mcover,
                                                                        unsigned int seed)
{
    gmx_radial_distribution_histogram_t* pr = nullptr;
    rvec                                 dist;
    double                               rmax;
    int                                  i, j;
#if GMX_OPENMP
    double**                  tgr;
    int                       tid;
    int                       nthreads;
    gmx::DefaultRandomEngine* trng = nullptr;
#endif
    int64_t                  mc_max;
    gmx::DefaultRandomEngine rng(seed);

    /* allocate memory for pr */
    snew(pr, 1);
    /* set some fields */
    pr->binwidth = binwidth;

    /*
     * create max dist rvec
     * dist = box[xx] + box[yy] + box[zz]
     */
    rvec_add(box[XX], box[YY], dist);
    rvec_add(box[ZZ], dist, dist);

    rmax = norm(dist);

    pr->grn = static_cast<int>(std::floor(rmax / pr->binwidth) + 1);

    snew(pr->gr, pr->grn);

    if (bMC)
    {
        /* Special case for setting automaticaly number of mc iterations to 1% of total number of direct iterations */
        if (mcover == -1)
        {
            mc_max = static_cast<int64_t>(std::floor(0.5 * 0.01 * isize * (isize - 1)));
        }
        else
        {
            mc_max = static_cast<int64_t>(std::floor(0.5 * mcover * isize * (isize - 1)));
        }
#if GMX_OPENMP
        nthreads = gmx_omp_get_max_threads();
        snew(tgr, nthreads);
        trng = new gmx::DefaultRandomEngine[nthreads];
        for (i = 0; i < nthreads; i++)
        {
            snew(tgr[i], pr->grn);
            trng[i].seed(rng());
        }
#    pragma omp parallel shared(tgr, trng) private(tid, i, j)
        {
            gmx::UniformIntDistribution<int> tdist(0, isize - 1);
            tid = gmx_omp_get_thread_num();
            /* now starting parallel threads */
#    pragma omp for
            for (int64_t mc = 0; mc < mc_max; mc++)
            {
                try
                {
                    i = tdist(trng[tid]); // [0,isize-1]
                    j = tdist(trng[tid]); // [0,isize-1]
                    if (i != j)
                    {
                        tgr[tid][static_cast<int>(std::floor(std::sqrt(distance2(x[index[i]], x[index[j]])) / binwidth))] +=
                                gsans->slength[index[i]] * gsans->slength[index[j]];
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            }
        }
        /* collecting data from threads */
        for (i = 0; i < pr->grn; i++)
        {
            for (j = 0; j < nthreads; j++)
            {
                pr->gr[i] += tgr[j][i];
            }
        }
        /* freeing memory for tgr and destroying trng */
        for (i = 0; i < nthreads; i++)
        {
            sfree(tgr[i]);
        }
        sfree(tgr);
        delete[] trng;
#else
        gmx::UniformIntDistribution<int> dist(0, isize - 1);
        for (int64_t mc = 0; mc < mc_max; mc++)
        {
            i = dist(rng); // [0,isize-1]
            j = dist(rng); // [0,isize-1]
            if (i != j)
            {
                pr->gr[static_cast<int>(std::floor(std::sqrt(distance2(x[index[i]], x[index[j]])) / binwidth))] +=
                        gsans->slength[index[i]] * gsans->slength[index[j]];
            }
        }
#endif
    }
    else
    {
#if GMX_OPENMP
        nthreads = gmx_omp_get_max_threads();
        /* Allocating memory for tgr arrays */
        snew(tgr, nthreads);
        for (i = 0; i < nthreads; i++)
        {
            snew(tgr[i], pr->grn);
        }
#    pragma omp parallel shared(tgr) private(tid, i, j)
        {
            tid = gmx_omp_get_thread_num();
/* starting parallel threads */
#    pragma omp for
            for (i = 0; i < isize; i++)
            {
                try
                {
                    for (j = 0; j < i; j++)
                    {
                        tgr[tid][static_cast<int>(std::floor(std::sqrt(distance2(x[index[i]], x[index[j]])) / binwidth))] +=
                                gsans->slength[index[i]] * gsans->slength[index[j]];
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            }
        }
        /* collecating data for pr->gr */
        for (i = 0; i < pr->grn; i++)
        {
            for (j = 0; j < nthreads; j++)
            {
                pr->gr[i] += tgr[j][i];
            }
        }
        /* freeing memory for tgr */
        for (i = 0; i < nthreads; i++)
        {
            sfree(tgr[i]);
        }
        sfree(tgr);
#else
        for (i = 0; i < isize; i++)
        {
            for (j = 0; j < i; j++)
            {
                pr->gr[static_cast<int>(std::floor(std::sqrt(distance2(x[index[i]], x[index[j]])) / binwidth))] +=
                        gsans->slength[index[i]] * gsans->slength[index[j]];
            }
        }
#endif
    }

    /* normalize if needed */
    if (bNORM)
    {
        normalize_probability(pr->grn, pr->gr);
    }

    snew(pr->r, pr->grn);
    for (i = 0; i < pr->grn; i++)
    {
        pr->r[i] = (pr->binwidth * i + pr->binwidth * 0.5);
    }

    return pr;
}

gmx_static_structurefactor_t* convert_histogram_to_intensity_curve(gmx_radial_distribution_histogram_t* pr,
                                                                   double start_q,
                                                                   double end_q,
                                                                   double q_step)
{
    gmx_static_structurefactor_t* sq = nullptr;
    int                           i, j;
    /* init data */
    snew(sq, 1);
    sq->qn = static_cast<int>(std::floor((end_q - start_q) / q_step));
    snew(sq->q, sq->qn);
    snew(sq->s, sq->qn);
    for (i = 0; i < sq->qn; i++)
    {
        sq->q[i] = start_q + i * q_step;
    }

    if (start_q == 0.0)
    {
        sq->s[0] = 1.0;
        for (i = 1; i < sq->qn; i++)
        {
            for (j = 0; j < pr->grn; j++)
            {
                sq->s[i] += (pr->gr[j] / pr->r[j]) * std::sin(sq->q[i] * pr->r[j]);
            }
            sq->s[i] /= sq->q[i];
        }
    }
    else
    {
        for (i = 0; i < sq->qn; i++)
        {
            for (j = 0; j < pr->grn; j++)
            {
                sq->s[i] += (pr->gr[j] / pr->r[j]) * std::sin(sq->q[i] * pr->r[j]);
            }
            sq->s[i] /= sq->q[i];
        }
    }

    return sq;
}
