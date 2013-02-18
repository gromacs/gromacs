/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "maths.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "macros.h"
#include "index.h"
#include "strdb.h"
#include "copyrite.h"
#include "tpxio.h"
#include "typedefs.h"
#include "statutil.h"
#include "oenv.h"
#include "gmxfio.h"
#include "xvgr.h"
#include "matio.h"
#include "gmx_ana.h"
#include "names.h"
#include "sfactor.h"
#include "gmx_omp.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Commont procedures for SANS and SAXS
 */

gmx_structurefactors_t *gmx_structurefactors_init(const char *datfn)
{

    /* Read the database for the structure factor of the different atoms */

    FILE                   *fp;
    char                    line[STRLEN];
    gmx_structurefactors_t *gsf;
    double                  a1, a2, a3, a4, b1, b2, b3, b4, c, coh_b;
    int                     p;
    int                     i;
    int                     nralloc = 10;
    int                     line_no;
    char                    atomn[32];
    fp      = libopen(datfn);
    line_no = 0;
    snew(gsf, 1);

    snew(gsf->atomnm, nralloc);
    snew(gsf->a, nralloc);
    snew(gsf->b, nralloc);
    snew(gsf->c, nralloc);
    snew(gsf->p, nralloc);
    snew(gsf->coh_b, nralloc);
    gsf->n       = NULL;
    gsf->nratoms = line_no;
    while (get_a_line(fp, line, STRLEN))
    {
        i = line_no;
        if (sscanf(line, "%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   atomn, &p, &a1, &a2, &a3, &a4, &b1, &b2, &b3, &b4, &c, &coh_b) == 12)
        {
            gsf->atomnm[i] = strdup(atomn);
            gsf->p[i]      = p;
            snew(gsf->a[i], 4);
            snew(gsf->b[i], 4);
            gsf->a[i][0]  = a1;
            gsf->a[i][1]  = a2;
            gsf->a[i][2]  = a3;
            gsf->a[i][3]  = a4;
            gsf->b[i][0]  = b1;
            gsf->b[i][1]  = b2;
            gsf->b[i][2]  = b3;
            gsf->b[i][3]  = b4;
            gsf->c[i]     = c;
            gsf->coh_b[i] = coh_b;
            line_no++;
            gsf->nratoms = line_no;
            if (line_no == nralloc)
            {
                nralloc += 10;
                srenew(gsf->atomnm, nralloc);
                srenew(gsf->a, nralloc);
                srenew(gsf->b, nralloc);
                srenew(gsf->c, nralloc);
                srenew(gsf->p, nralloc);
                srenew(gsf->coh_b, nralloc);
            }
        }
        else
        {
            fprintf(stderr, "WARNING: Error in file %s at line %d ignored\n",
                    datfn, line_no);
        }
    }

    srenew(gsf->atomnm, gsf->nratoms);
    srenew(gsf->a, gsf->nratoms);
    srenew(gsf->b, gsf->nratoms);
    srenew(gsf->c, gsf->nratoms);
    srenew(gsf->p, gsf->nratoms);
    srenew(gsf->coh_b, gsf->nratoms);

    fclose(fp);

    return (gmx_structurefactors_t *) gsf;

}

int *create_indexed_atom_type (reduced_atom_t * atm, int size)
{
    /*
     * create an index of the atom types found in a  group
     * i.e.: for water index_atp[0]=type_number_of_O and
     *                 index_atp[1]=type_number_of_H
     *
     * the last element is set to 0
     */
    int            *index_atp, i, i_tmp, j;

    reduced_atom_t *att = (reduced_atom_t *)atm;

    snew (index_atp, 1);
    i_tmp        = 1;
    index_atp[0] = att[0].t;
    for (i = 1; i < size; i++)
    {
        for (j = 0; j < i_tmp; j++)
        {
            if (att[i].t == index_atp[j])
            {
                break;
            }
        }
        if (j == i_tmp)     /* i.e. no indexed atom type is  == to atm[i].t */
        {
            i_tmp++;
            srenew (index_atp, i_tmp * sizeof (int));
            index_atp[i_tmp - 1] = att[i].t;
        }
    }
    i_tmp++;
    srenew (index_atp, i_tmp * sizeof (int));
    index_atp[i_tmp - 1] = 0;
    return index_atp;
}

void rearrange_atoms (reduced_atom_t * positions, t_trxframe *fr, atom_id * index,
                      int isize, t_topology * top, gmx_bool flag, gmx_structurefactors_t *gsf)
/* given the group's index, return the (continuous) array of atoms */
{
    int             i;

    reduced_atom_t *pos = (reduced_atom_t *)positions;

    if (flag)
    {
        for (i = 0; i < isize; i++)
        {
            pos[i].t =
                return_atom_type (*(top->atoms.atomname[index[i]]), gsf);
        }
    }
    for (i = 0; i < isize; i++)
    {
        copy_rvec (fr->x[index[i]], pos[i].x);
    }

    positions = (reduced_atom_t *)pos;
}


int return_atom_type (const char *name, gmx_structurefactors_t *gsf)
{
    typedef struct {
        const char *name;
        int         nh;
    } t_united_h;
    t_united_h              uh[] = {
        { "CH1", 1 }, { "CH2", 2 }, { "CH3", 3 },
        { "CS1", 1 }, { "CS2", 2 }, { "CS3", 3 },
        { "CP1", 1 }, { "CP2", 2 }, { "CP3", 3 }
    };
    int                     i, cnt = 0;
    int                    *tndx;
    int                     nrc;
    int                     fndx = 0;
    int                     NCMT;

    gmx_structurefactors_t *gsft = (gmx_structurefactors_t *)gsf;

    NCMT = gsft->nratoms;

    snew(tndx, NCMT);

    for (i = 0; (i < asize(uh)); i++)
    {
        if (strcmp(name, uh[i].name) == 0)
        {
            return NCMT-1+uh[i].nh;
        }
    }

    for (i = 0; (i < NCMT); i++)
    {
        if (strncmp (name, gsft->atomnm[i], strlen(gsft->atomnm[i])) == 0)
        {
            tndx[cnt] = i;
            cnt++;
        }
    }

    if (cnt == 0)
    {
        gmx_fatal(FARGS, "\nError: atom (%s) not in list (%d types checked)!\n",
                  name, i);
    }
    else
    {
        nrc = 0;
        for (i = 0; i < cnt; i++)
        {
            if (strlen(gsft->atomnm[tndx[i]]) > (size_t)nrc)
            {
                nrc  = strlen(gsft->atomnm[tndx[i]]);
                fndx = tndx[i];
            }
        }

        return fndx;
    }

    return 0;
}

int gmx_structurefactors_get_sf(gmx_structurefactors_t *gsf, int elem, real a[4], real b[4], real *c)
{

    int                     success;
    int                     i;
    gmx_structurefactors_t *gsft = (gmx_structurefactors_t *)gsf;
    success = 0;

    for (i = 0; i < 4; i++)
    {
        a[i] = gsft->a[elem][i];
        b[i] = gsft->b[elem][i];
        *c   = gsft->c[elem];
    }

    success += 1;
    return success;
}

void done_gmx_structurefactors(gmx_structurefactors_t *gsf)
{

    int                     i;
    gmx_structurefactors_t *sf;
    sf = (gmx_structurefactors_t *) gsf;

    for (i = 0; i < sf->nratoms; i++)
    {
        sfree(sf->a[i]);
        sfree(sf->b[i]);
        sfree(sf->atomnm[i]);
    }

    sfree(sf->a);
    sfree(sf->b);
    sfree(sf->atomnm);
    sfree(sf->p);
    sfree(sf->c);
    sfree(sf->coh_b);

    sfree(sf);

}

void check_binwidth(real binwidth)
{
    real smallest_bin = 0.1;
    if (binwidth < smallest_bin)
    {
        gmx_fatal(FARGS, "Binwidth shouldnt be smaller then smallest bond length (H-H bond ~0.1nm) in a box");
    }
}

void check_mcover(real mcover)
{
    if (mcover > 1.0)
    {
        gmx_fatal(FARGS, "mcover should be -1 or (0,1]");
    }
    else if ((mcover < 0)&(mcover != -1))
    {
        gmx_fatal(FARGS, "mcover should be -1 or (0,1]");
    }
    else
    {
        return;
    }
}

void normalize_probability(int n, double *a)
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

real max_dx(matrix box)
{
    rvec dist;
    /*
     * create max dist rvec
     * dist = box[xx] + box[yy] + box[zz]
     */
    rvec_add(box[XX], box[YY], dist);
    rvec_add(box[ZZ], dist, dist);

    return norm(dist);
}

/*
 * SAXS funtions
 */

t_complex *** rc_tensor_allocation(int x, int y, int z)
{
    t_complex ***t;
    int          i, j;

    snew(t, x);
    t = (t_complex ***)calloc(x, sizeof(t_complex**));
    if (!t)
    {
        exit(fprintf(stderr, "\nallocation error"));
    }
    t[0] = (t_complex **)calloc(x*y, sizeof(t_complex*));
    if (!t[0])
    {
        exit(fprintf(stderr, "\nallocation error"));
    }
    t[0][0] = (t_complex *)calloc(x*y*z, sizeof(t_complex));
    if (!t[0][0])
    {
        exit(fprintf(stderr, "\nallocation error"));
    }

    for (j = 1; j < y; j++)
    {
        t[0][j] = t[0][j-1] + z;
    }
    for (i = 1; i < x; i++)
    {
        t[i]    = t[i-1] + y;
        t[i][0] = t[i-1][0] + y*z;
        for (j = 1; j < y; j++)
        {
            t[i][j] = t[i][j-1] + z;
        }
    }
    return t;
}

void compute_structure_factor (structure_factor_t * sft, matrix box,
                               reduced_atom_t * red, int isize, real start_q,
                               real end_q, int group, real **sf_table)
{
    structure_factor_t *sf   = (structure_factor_t *)sft;
    reduced_atom_t     *redt = (reduced_atom_t *)red;

    t_complex        ***tmpSF;
    rvec                k_factor;
    real                kdotx, asf, kx, ky, kz, krr;
    int                 kr, maxkx, maxky, maxkz, i, j, k, p, *counter;


    k_factor[XX] = 2 * M_PI / box[XX][XX];
    k_factor[YY] = 2 * M_PI / box[YY][YY];
    k_factor[ZZ] = 2 * M_PI / box[ZZ][ZZ];

    maxkx = (int) (end_q / k_factor[XX] + 0.5);
    maxky = (int) (end_q / k_factor[YY] + 0.5);
    maxkz = (int) (end_q / k_factor[ZZ] + 0.5);

    snew (counter, sf->n_angles);

    tmpSF = rc_tensor_allocation(maxkx, maxky, maxkz);
    /*
     * The big loop...
     * compute real and imaginary part of the structure factor for every
     * (kx,ky,kz))
     */
    fprintf(stderr, "\n");
    for (i = 0; i < maxkx; i++)
    {
        fprintf (stderr, "\rdone %3.1f%%     ", (double)(100.0*(i+1))/maxkx);
        kx = i * k_factor[XX];
        for (j = 0; j < maxky; j++)
        {
            ky = j * k_factor[YY];
            for (k = 0; k < maxkz; k++)
            {
                if (i != 0 || j != 0 || k != 0)
                {
                    kz  = k * k_factor[ZZ];
                    krr = sqrt (sqr (kx) + sqr (ky) + sqr (kz));
                    if (krr >= start_q && krr <= end_q)
                    {
                        kr = (int) (krr/sf->ref_k + 0.5);
                        if (kr < sf->n_angles)
                        {
                            counter[kr]++;  /* will be used for the copmutation
                                               of the average*/
                            for (p = 0; p < isize; p++)
                            {
                                asf = sf_table[redt[p].t][kr];

                                kdotx = kx * redt[p].x[XX] +
                                    ky * redt[p].x[YY] + kz * redt[p].x[ZZ];

                                tmpSF[i][j][k].re += cos (kdotx) * asf;
                                tmpSF[i][j][k].im += sin (kdotx) * asf;
                            }
                        }
                    }
                }
            }
        }
    }               /* end loop on i */
    /*
     *  compute the square modulus of the structure factor, averaging on the surface
     *  kx*kx + ky*ky + kz*kz = krr*krr
     *  note that this is correct only for a (on the macroscopic scale)
     *  isotropic system.
     */
    for (i = 0; i < maxkx; i++)
    {
        kx = i * k_factor[XX]; for (j = 0; j < maxky; j++)
        {
            ky = j * k_factor[YY]; for (k = 0; k < maxkz; k++)
            {
                kz = k * k_factor[ZZ]; krr = sqrt (sqr (kx) + sqr (ky)
                                                   + sqr (kz)); if (krr >= start_q && krr <= end_q)
                {
                    kr = (int) (krr / sf->ref_k + 0.5);
                    if (kr < sf->n_angles && counter[kr] != 0)
                    {
                        sf->F[group][kr] +=
                            (sqr (tmpSF[i][j][k].re) +
                             sqr (tmpSF[i][j][k].im))/ counter[kr];
                    }
                }
            }
        }
    }
    sfree (counter); free(tmpSF[0][0]); free(tmpSF[0]); free(tmpSF);
}

void save_data (structure_factor_t *sft, const char *file, int ngrps,
                real start_q, real end_q, const output_env_t oenv)
{

    FILE               *fp;
    int                 i, g = 0;
    double             *tmp, polarization_factor, A;

    structure_factor_t *sf = (structure_factor_t *)sft;

    fp = xvgropen (file, "Scattering Intensity", "q (1/nm)",
                   "Intensity (a.u.)", oenv);

    snew (tmp, ngrps);

    for (g = 0; g < ngrps; g++)
    {
        for (i = 0; i < sf->n_angles; i++)
        {

            /*
             *          theta is half the angle between incoming and scattered vectors.
             *
             *          polar. fact. = 0.5*(1+cos^2(2*theta)) = 1 - 0.5 * sin^2(2*theta)
             *
             *          sin(theta) = q/(2k) := A  ->  sin^2(theta) = 4*A^2 (1-A^2) ->
             *          -> 0.5*(1+cos^2(2*theta)) = 1 - 2 A^2 (1-A^2)
             */
            A                   = (double) (i * sf->ref_k) / (2.0 * sf->momentum);
            polarization_factor = 1 - 2.0 * sqr (A) * (1 - sqr (A));
            sf->F[g][i]        *= polarization_factor;
        }
    }
    for (i = 0; i < sf->n_angles; i++)
    {
        if (i * sf->ref_k >= start_q && i * sf->ref_k <= end_q)
        {
            fprintf (fp, "%10.5f  ", i * sf->ref_k);
            for (g = 0; g < ngrps; g++)
            {
                fprintf (fp, "  %10.5f ", (sf->F[g][i]) /( sf->total_n_atoms*
                                                           sf->nSteps));
            }
            fprintf (fp, "\n");
        }
    }

    ffclose (fp);
}


double CMSF (gmx_structurefactors_t *gsf, int type, int nh, double lambda, double sin_theta)
/*
 * return Cromer-Mann fit for the atomic scattering factor:
 * sin_theta is the sine of half the angle between incoming and scattered
 * vectors. See g_sq.h for a short description of CM fit.
 */
{
    int    i, success;
    double tmp = 0.0, k2;
    real  *a, *b;
    real   c;

    snew(a, 4);
    snew(b, 4);

    /*
     *
     * f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
     *             i=1,4
     */

    /*
     *  united atoms case
     *  CH2 / CH3 groups
     */
    if (nh > 0)
    {
        tmp = (CMSF (gsf, return_atom_type ("C", gsf), 0, lambda, sin_theta) +
               nh*CMSF (gsf, return_atom_type ("H", gsf), 0, lambda, sin_theta));
    }
    /* all atom case */
    else
    {
        k2      = (sqr (sin_theta) / sqr (10.0 * lambda));
        success = gmx_structurefactors_get_sf(gsf, type, a, b, &c);
        tmp     = c;
        for (i = 0; (i < 4); i++)
        {
            tmp += a[i] * exp (-b[i] * k2);
        }
    }
    return tmp;
}



real **gmx_structurefactors_table(gmx_structurefactors_t *gsf, real momentum, real ref_k, real lambda, int n_angles)
{

    int                     NCMT;
    int                     nsftable;
    int                     i, j;
    double                  q, sin_theta;
    real                  **sf_table;
    gmx_structurefactors_t *gsft = (gmx_structurefactors_t *)gsf;

    NCMT     = gsft->nratoms;
    nsftable = NCMT+3;

    snew (sf_table, nsftable);
    for (i = 0; (i < nsftable); i++)
    {
        snew (sf_table[i], n_angles);
        for (j = 0; j < n_angles; j++)
        {
            q = ((double) j * ref_k);
            /* theta is half the angle between incoming
               and scattered wavevectors. */
            sin_theta = q / (2.0 * momentum);
            if (i < NCMT)
            {
                sf_table[i][j] = CMSF (gsf, i, 0, lambda, sin_theta);
            }
            else
            {
                sf_table[i][j] = CMSF (gsf, i, i-NCMT+1, lambda, sin_theta);
            }
        }
    }
    return sf_table;
}


extern real **compute_scattering_factor_table (gmx_structurefactors_t *gsf, structure_factor_t *sft, int *nsftable)
{
    /*
     *  this function build up a table of scattering factors for every atom
     *  type and for every scattering angle.
     */

    double              hc = 1239.842;
    real             ** sf_table;

    structure_factor_t *sf = (structure_factor_t *)sft;


    /* \hbar \omega \lambda = hc = 1239.842 eV * nm */
    sf->momentum = ((double) (2. * 1000.0 * M_PI * sf->energy) / hc);
    sf->lambda   = hc / (1000.0 * sf->energy);
    fprintf (stderr, "\nwavelenght = %f nm\n", sf->lambda);

    sf_table = gmx_structurefactors_table(gsf, sf->momentum, sf->ref_k, sf->lambda, sf->n_angles);

    return sf_table;
}

/*
 * SANS
 */

gmx_sans_t *gmx_sans_init (t_topology *top, gmx_structurefactors_t *sf)
{
    gmx_sans_t    *gsans = NULL;
    int            i, j;
    /* Try to assing scattering length from nsfactor.dat */
    snew(gsans, 1);
    snew(gsans->slength, top->atoms.nr);
    /* copy topology data */
    gsans->top = top;
    for (i = 0; i < top->atoms.nr; i++)
    {
        for (j = 0; j < sf->nratoms; j++)
        {
            if (top->atoms.atom[i].atomnumber == sf->p[j])
            {
                /* we need special case for H and D */
                if (top->atoms.atom[i].atomnumber == 1)
                {
                    if (top->atoms.atom[i].m == 1.008000)
                    {
                        gsans->slength[i] = sf->coh_b[0];
                    }
                    else
                    {
                        gsans->slength[i] = sf->coh_b[1];
                    }
                }
                else
                {
                    gsans->slength[i] = sf->coh_b[j];
                }
            }
        }
    }

    return (gmx_sans_t *) gsans;
}

void done_sans(gmx_sans_t *gsans)
{
    sfree(gsans->slength);
    done_top(gsans->top);
    sfree(gsans);
}

gmx_radial_distribution_histogram_t *calc_radial_distribution_histogram (
        gmx_sans_t  *gsans,
        rvec        *x,
        rvec        *xf,
        matrix       box,
        matrix       boxf,
        atom_id     *index,
        int          isize,
        double       binwidth,
        gmx_bool     bMC,
        gmx_bool     bNSE,
        real         mcover,
        unsigned int seed)
{
    gmx_radial_distribution_histogram_t    *pr = NULL;
    rvec              dist;
    double            rmax;
    int               i, j;
#ifdef GMX_OPENMP
    double          **tgr;
    int               tid;
    int               nthreads;
    gmx_rng_t        *trng = NULL;
#endif
    gmx_large_int_t   mc  = 0, max;
    gmx_rng_t         rng = NULL;

    /* allocate memory for pr */
    snew(pr, 1);
    /* set some fields */
    pr->binwidth = binwidth;
    if (bNSE)
    {
        if (max_dx(box) - max_dx(boxf) > 0)
        {
            rmax = max_dx(box);
        }
        else
        {
            rmax = max_dx(boxf);
        }
    }
    else
    {
        rmax = max_dx(box);
    }

    pr->grn = (int)floor(rmax/pr->binwidth)+1;
    rmax    = pr->grn*pr->binwidth;

    snew(pr->gr, pr->grn);

    if (bMC)
    {
        /* Special case for setting automaticaly number of mc iterations to 1% of total number of direct iterations */
        if (mcover == -1)
        {
            if (bNSE)
            {
                max = (gmx_large_int_t)floor(0.01*isize*isize);
            }
            else
            {
                max = (gmx_large_int_t)floor(0.5*0.01*isize*(isize-1));
            }
        }
        else
        {
            if (bNSE)
            {
                max = (gmx_large_int_t)floor(mcover*isize*isize);
            }
            else
            {
                max = (gmx_large_int_t)floor(0.5*mcover*isize*(isize-1));
            }
        }
        rng = gmx_rng_init(seed);
#ifdef GMX_OPENMP
        nthreads = gmx_omp_get_max_threads();
        snew(tgr, nthreads);
        snew(trng, nthreads);
        for (i = 0; i < nthreads; i++)
        {
            snew(tgr[i], pr->grn);
            trng[i] = gmx_rng_init(gmx_rng_uniform_uint32(rng));
        }
#pragma omp parallel shared(tgr,trng,mc) private(tid,i,j)
        {
            tid = gmx_omp_get_thread_num();
            /* now starting parallel threads */
#pragma omp for
            for (mc = 0; mc < max; mc++)
            {
                i = (int)floor(gmx_rng_uniform_real(trng[tid])*isize);
                j = (int)floor(gmx_rng_uniform_real(trng[tid])*isize);
                if (bNSE)
                {
                    /* we already copyed x[index[i]] to gnse->x[frame] */
                    tgr[tid][(int)floor(sqrt(distance2(x[i], xf[j]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                }
                else
                {
                    if (i != j)
                    {
                        tgr[tid][(int)floor(sqrt(distance2(x[index[i]], x[index[j]]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                    }
                }
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
            gmx_rng_destroy(trng[i]);
        }
        sfree(tgr);
        sfree(trng);
#else
        for (mc = 0; mc < max; mc++)
        {
            i = (int)floor(gmx_rng_uniform_real(rng)*isize);
            j = (int)floor(gmx_rng_uniform_real(rng)*isize);
            if (bNSE)
            {
                /* we already copyed x[index[i]] to gnse->x[frame] */
                pr->gr[(int)floor(sqrt(distance2(x[i], xf[j]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
            }
            else
            {
                if (i != j)
                {
                    pr->gr[(int)floor(sqrt(distance2(x[index[i]], x[index[j]]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                }
            }
        }
#endif
        gmx_rng_destroy(rng);
    }
    else
    {
#ifdef GMX_OPENMP
        nthreads = gmx_omp_get_max_threads();
        /* Allocating memory for tgr arrays */
        snew(tgr, nthreads);
        for (i = 0; i < nthreads; i++)
        {
            snew(tgr[i], pr->grn);
        }
#pragma omp parallel shared(tgr) private(tid,i,j)
        {
            tid = gmx_omp_get_thread_num();
            /* starting parallel threads */
            if (bNSE)
            {
#pragma omp for
                /* we already copyed x[index[i]] to gnse->x[frame] */
                for (i = 0; i < isize; i++)
                {
                    for (j = 0; j < isize; j++)
                    {
                        tgr[tid][(int)floor(sqrt(distance2(x[i], xf[j]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                    }
                }
            }
            else
            {
#pragma omp for
                for (i = 0; i < isize; i++)
                {
                    for (j = 0; j < i; j++)
                    {
                        tgr[tid][(int)floor(sqrt(distance2(x[index[i]], x[index[j]]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                    }
                }
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
        if (bNSE)
        {
            /* we already copyed x[index[i]] to gnse->x[frame] */
            for (i = 0; i < isize; i++)
            {
                for (j = 0; j < isize; j++)
                {
                    pr->gr[(int)floor(sqrt(distance2(x[i], xf[j]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                }
            }
        }
        else
        {
            for (i = 0; i < isize; i++)
            {
                for (j = 0; j < i; j++)
                {
                    pr->gr[(int)floor(sqrt(distance2(x[index[i]], x[index[j]]))/binwidth)] += gsans->slength[index[i]]*gsans->slength[index[j]];
                }
            }
        }
#endif
    }

    snew(pr->r, pr->grn);
    for (i = 0; i < pr->grn; i++)
    {
        pr->r[i] = (pr->binwidth*i+pr->binwidth*0.5);
    }

    return (gmx_radial_distribution_histogram_t *) pr;
}

gmx_sans_structurefactor_t *convert_histogram_to_intensity_curve (gmx_radial_distribution_histogram_t *pr, double start_q, double end_q, double q_step)
{
    gmx_sans_structurefactor_t    *sq = NULL;
    int                            i, j;
    /* init data */
    snew(sq, 1);
    sq->qn = (int)floor((end_q-start_q)/q_step);
    snew(sq->q, sq->qn);
    snew(sq->s, sq->qn);
    for (i = 0; i < sq->qn; i++)
    {
        sq->q[i] = start_q+i*q_step;
    }

    if (start_q == 0.0)
    {
        sq->s[0] = 1.0;
        for (i = 1; i < sq->qn; i++)
        {
            for (j = 0; j < pr->grn; j++)
            {
                sq->s[i] += (pr->gr[j]/pr->r[j])*sin(sq->q[i]*pr->r[j]);
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
                sq->s[i] += (pr->gr[j]/pr->r[j])*sin(sq->q[i]*pr->r[j]);
            }
            sq->s[i] /= sq->q[i];
        }
    }

    return (gmx_sans_structurefactor_t *) sq;
}


#ifdef __cplusplus
}
#endif
