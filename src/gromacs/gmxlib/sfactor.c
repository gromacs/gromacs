/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/math/utilities.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "macros.h"
#include "index.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "typedefs.h"
#include "oenv.h"
#include "gromacs/fileio/gmxfio.h"
#include "xvgr.h"
#include "gromacs/fileio/matio.h"
#include "names.h"
#include "sfactor.h"


typedef struct gmx_structurefactors {
    int    nratoms;
    int   *p;      /* proton number */
    int   *n;      /* neutron number */
    /* Parameters for the Cromer Mann fit */
    real **a;      /* parameter a */
    real **b;      /* parameter b */
    real  *c;      /* parameter c */
    char **atomnm; /* atomname */

} gmx_structurefactors;

typedef struct reduced_atom{
    rvec x;
    int  t;
} reduced_atom;


typedef struct structure_factor
{
    int       n_angles;
    int       n_groups;
    double    lambda;
    double    energy;
    double    momentum;
    double    ref_k;
    double  **F;
    int       nSteps;
    int       total_n_atoms;
} structure_factor;


extern int * create_indexed_atom_type (reduced_atom_t * atm, int size)
{
/*
 * create an index of the atom types found in a  group
 * i.e.: for water index_atp[0]=type_number_of_O and
 *                 index_atp[1]=type_number_of_H
 *
 * the last element is set to 0
 */
    int          *index_atp, i, i_tmp, j;

    reduced_atom *att = (reduced_atom *)atm;

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
        if (j == i_tmp) /* i.e. no indexed atom type is  == to atm[i].t */
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



extern t_complex *** rc_tensor_allocation(int x, int y, int z)
{
    t_complex ***t;
    int          i, j;

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


extern void compute_structure_factor (structure_factor_t * sft, matrix box,
                                      reduced_atom_t * red, int isize, real start_q,
                                      real end_q, int group, real **sf_table)
{
    structure_factor *sf   = (structure_factor *)sft;
    reduced_atom     *redt = (reduced_atom *)red;

    t_complex      ***tmpSF;
    rvec              k_factor;
    real              kdotx, asf, kx, ky, kz, krr;
    int               kr, maxkx, maxky, maxkz, i, j, k, p, *counter;


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
                            counter[kr]++; /* will be used for the copmutation
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


extern gmx_structurefactors_t *gmx_structurefactors_init(const char *datfn)
{

    /* Read the database for the structure factor of the different atoms */

    FILE                 *fp;
    char                  line[STRLEN];
    gmx_structurefactors *gsf;
    double                a1, a2, a3, a4, b1, b2, b3, b4, c;
    int                   p;
    int                   i;
    int                   nralloc = 10;
    int                   line_no;
    char                  atomn[32];
    fp      = libopen(datfn);
    line_no = 0;
    snew(gsf, 1);

    snew(gsf->atomnm, nralloc);
    snew(gsf->a, nralloc);
    snew(gsf->b, nralloc);
    snew(gsf->c, nralloc);
    snew(gsf->p, nralloc);
    gsf->n       = NULL;
    gsf->nratoms = line_no;
    while (get_a_line(fp, line, STRLEN))
    {
        i = line_no;
        if (sscanf(line, "%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   atomn, &p, &a1, &a2, &a3, &a4, &b1, &b2, &b3, &b4, &c) == 11)
        {
            gsf->atomnm[i] = strdup(atomn);
            gsf->p[i]      = p;
            snew(gsf->a[i], 4);
            snew(gsf->b[i], 4);
            gsf->a[i][0] = a1;
            gsf->a[i][1] = a2;
            gsf->a[i][2] = a3;
            gsf->a[i][3] = a4;
            gsf->b[i][0] = b1;
            gsf->b[i][1] = b2;
            gsf->b[i][2] = b3;
            gsf->b[i][3] = b4;
            gsf->c[i]    = c;
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

    fclose(fp);

    return (gmx_structurefactors_t *) gsf;

}


extern void rearrange_atoms (reduced_atom_t * positions, t_trxframe *fr, atom_id * index,
                             int isize, t_topology * top, gmx_bool flag, gmx_structurefactors_t *gsf)
/* given the group's index, return the (continuous) array of atoms */
{
    int           i;

    reduced_atom *pos = (reduced_atom *)positions;

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


extern int return_atom_type (const char *name, gmx_structurefactors_t *gsf)
{
    typedef struct {
        const char *name;
        int         nh;
    } t_united_h;
    t_united_h            uh[] = {
        { "CH1", 1 }, { "CH2", 2 }, { "CH3", 3 },
        { "CS1", 1 }, { "CS2", 2 }, { "CS3", 3 },
        { "CP1", 1 }, { "CP2", 2 }, { "CP3", 3 }
    };
    int                   i, cnt = 0;
    int                  *tndx;
    int                   nrc;
    int                   fndx = 0;
    int                   NCMT;

    gmx_structurefactors *gsft = (gmx_structurefactors *)gsf;

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

extern int gmx_structurefactors_get_sf(gmx_structurefactors_t *gsf, int elem, real a[4], real b[4], real *c)
{

    int                   success;
    int                   i;
    gmx_structurefactors *gsft = (gmx_structurefactors *)gsf;
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

extern int do_scattering_intensity (const char* fnTPS, const char* fnNDX,
                                    const char* fnXVG, const char *fnTRX,
                                    const char* fnDAT,
                                    real start_q, real end_q,
                                    real energy, int ng, const output_env_t oenv)
{
    int                     i, *isize, flags = TRX_READ_X, **index_atp;
    t_trxstatus            *status;
    char                  **grpname, title[STRLEN];
    atom_id               **index;
    t_topology              top;
    int                     ePBC;
    t_trxframe              fr;
    reduced_atom_t        **red;
    structure_factor       *sf;
    rvec                   *xtop;
    real                  **sf_table;
    int                     nsftable;
    matrix                  box;
    double                  r_tmp;

    gmx_structurefactors_t *gmx_sf;
    real                   *a, *b, c;
    int                     success;

    snew(a, 4);
    snew(b, 4);


    gmx_sf = gmx_structurefactors_init(fnDAT);

    success = gmx_structurefactors_get_sf(gmx_sf, 0, a, b, &c);

    snew (sf, 1);
    sf->energy = energy;

    /* Read the topology informations */
    read_tps_conf (fnTPS, title, &top, &ePBC, &xtop, NULL, box, TRUE);
    sfree (xtop);

    /* groups stuff... */
    snew (isize, ng);
    snew (index, ng);
    snew (grpname, ng);

    fprintf (stderr, "\nSelect %d group%s\n", ng,
             ng == 1 ? "" : "s");
    if (fnTPS)
    {
        get_index (&top.atoms, fnNDX, ng, isize, index, grpname);
    }
    else
    {
        rd_index (fnNDX, ng, isize, index, grpname);
    }

    /* The first time we read data is a little special */
    read_first_frame (oenv, &status, fnTRX, &fr, flags);

    sf->total_n_atoms = fr.natoms;

    snew (red, ng);
    snew (index_atp, ng);

    r_tmp = max (box[XX][XX], box[YY][YY]);
    r_tmp = (double) max (box[ZZ][ZZ], r_tmp);

    sf->ref_k = (2.0 * M_PI) / (r_tmp);
    /* ref_k will be the reference momentum unit */
    sf->n_angles = (int) (end_q / sf->ref_k + 0.5);

    snew (sf->F, ng);
    for (i = 0; i < ng; i++)
    {
        snew (sf->F[i], sf->n_angles);
    }
    for (i = 0; i < ng; i++)
    {
        snew (red[i], isize[i]);
        rearrange_atoms (red[i], &fr, index[i], isize[i], &top, TRUE, gmx_sf);
        index_atp[i] = create_indexed_atom_type (red[i], isize[i]);
    }

    sf_table = compute_scattering_factor_table (gmx_sf, (structure_factor_t *)sf);


    /* This is the main loop over frames */

    do
    {
        sf->nSteps++;
        for (i = 0; i < ng; i++)
        {
            rearrange_atoms (red[i], &fr, index[i], isize[i], &top, FALSE, gmx_sf);

            compute_structure_factor ((structure_factor_t *)sf, box, red[i], isize[i],
                                      start_q, end_q, i, sf_table);
        }
    }

    while (read_next_frame (oenv, status, &fr));

    save_data ((structure_factor_t *)sf, fnXVG, ng, start_q, end_q, oenv);


    sfree(a);
    sfree(b);

    gmx_structurefactors_done(gmx_sf);

    return 0;
}


extern void save_data (structure_factor_t *sft, const char *file, int ngrps,
                       real start_q, real end_q, const output_env_t oenv)
{

    FILE             *fp;
    int               i, g = 0;
    double           *tmp, polarization_factor, A;

    structure_factor *sf = (structure_factor *)sft;

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

    gmx_ffclose (fp);
}


extern double CMSF (gmx_structurefactors_t *gsf, int type, int nh, double lambda, double sin_theta)
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



extern real **gmx_structurefactors_table(gmx_structurefactors_t *gsf, real momentum, real ref_k, real lambda, int n_angles)
{

    int                   NCMT;
    int                   nsftable;
    int                   i, j;
    double                q, sin_theta;
    real                **sf_table;
    gmx_structurefactors *gsft = (gmx_structurefactors *)gsf;

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

extern void gmx_structurefactors_done(gmx_structurefactors_t *gsf)
{

    int                   i;
    gmx_structurefactors *sf;
    sf = (gmx_structurefactors *) gsf;

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

    sfree(sf);

}

extern real **compute_scattering_factor_table (gmx_structurefactors_t *gsf, structure_factor_t *sft)
{
/*
 *  this function build up a table of scattering factors for every atom
 *  type and for every scattering angle.
 */

    double            hc = 1239.842;
    real           ** sf_table;

    structure_factor *sf = (structure_factor *)sft;


    /* \hbar \omega \lambda = hc = 1239.842 eV * nm */
    sf->momentum = ((double) (2. * 1000.0 * M_PI * sf->energy) / hc);
    sf->lambda   = hc / (1000.0 * sf->energy);
    fprintf (stderr, "\nwavelenght = %f nm\n", sf->lambda);

    sf_table = gmx_structurefactors_table(gsf, sf->momentum, sf->ref_k, sf->lambda, sf->n_angles);

    return sf_table;
}
