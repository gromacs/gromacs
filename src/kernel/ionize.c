/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include <string.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "random.h"
#include "gmx_random.h"
#include "physics.h"
#include "xvgr.h"
#include "vec.h"
#include "pbc.h"
#include "txtdump.h"
#include "ionize.h"
#include "names.h"
#include "futil.h"
#include "network.h"
#include "mtop_util.h"
#include "gmxfio.h"

typedef struct {
    real photo, coh, incoh, incoh_abs;
} t_cross;

/* THIS TABLE HAS ADDED A 12 keV COLUMN TO HYDROGEN, CARBON,  */
/* OXYGEN, NITROGEN AND SULPHUR BY FITTING A QUADRATIC TO THE */
/* POINTS 8keV, 10keV and 12keV - now contains 6, 8, 10, 12,  */
/* 15 and 20 keV                                              */
/* Units are barn. They are converted to nm^2 by multiplying  */
/* by 1e-10, which is done in Imax (ionize.c)                 */
/* Update: contains energy 2 KeV and B, Na, Li, Al, Mg        */
/* Update: data taken from NIST XCOM                          */

static const t_cross   cross_sec_h[] = {
    { 5.21E-01,    3.39E-01,     3.21E-01,        -1 },
    { 2.63e-2,     1.01e-1,      5.49e-1,         7.12e-3 },
    { 9.79e-3,     6.18e-2,      5.83e-1,         9.60e-3 },
    { 4.55e-3,     4.16e-2,      5.99e-1,         1.19e-2 },
    { 1.52e-3,     2.79e-2,      6.08e-1,         1.41e-2 },
    { 1.12e-3,     1.96e-2,      6.09e-1,         1.73e-2 },
    { 4.16e-4,     1.13e-2,      6.07e-1,         2.23e-2 }
};
static const t_cross   cross_sec_c[] = {
    { 3.10E+3,     1.42E+1,      1.03E+0,         -1 },
    { 1.99e+2,     5.88e+0,      2.29e+0,         3.06e-2 },
    { 8.01e+1,     4.22e+0,      2.56e+0,         4.38e-2 },
    { 3.92e+1,     3.26e+0,      2.74e+0,         5.72e-2 },
    { 2.19e+1,     2.55e+0,      2.88e+0,         7.07e-2 },
    { 1.06e+1,     1.97e+0,      3.04e+0,         9.15e-2 },
    { 4.15e+0,     1.30e+0,      3.20e+0,         1.24e-1 }
};
static const t_cross   cross_sec_n[] = {
    { 5.78E+3,     2.13E+1,      1.11E+0,         -1 },
    { 3.91e+2,     8.99e+0,      2.49e+0,         3.43e-2 },
    { 1.59e+2,     6.29e+0,      2.86e+0,         5.01e-2 },
    { 7.88e+1,     4.76e+0,      3.10e+0,         6.57e-2 },
    { 4.42e+1,     3.66e+0,      3.28e+0,         8.13e-2 },
    { 2.16e+1,     2.82e+0,      3.46e+0,         1.05e-1 },
    { 8.52e+0,     1.88e+0,      3.65e+0,         1.43e-1 }
};
static const t_cross   cross_sec_o[] = {
    { 9.74E+3,     3.00E+1,      1.06E+0,         -1 },
    { 6.90e+2,     1.33e+1,      2.66e+0,         3.75e-2 },
    { 2.84e+2,     9.21e+0,      3.14e+0,         5.62e-2 },
    { 1.42e+2,     6.85e+0,      3.44e+0,         7.43e-2 },
    { 8.01e+1,     5.18e+0,      3.66e+0,         9.20e-2 },
    { 3.95e+1,     3.97e+0,      3.87e+0,         1.18e-1 },
    { 1.57e+1,     2.64e+0,      4.10e+0,         1.61e-1 }
};
static const t_cross   cross_sec_s[] = {
    { 1.07E+5,      1.15E+2,     2.03E+0,         -1 },
    { 1.10e+4,      5.54e+1,     3.98e+0,         5.42e-2 },
    { 4.91e+3,      4.29e+1,     4.71e+0,         8.38e-2 },
    { 2.58e+3,      3.36e+1,     5.32e+0,         1.16e-1 },
    { 1.52e+3,      2.64e+1,     5.81e+0,         1.48e-1 },
    { 7.82e+2,      1.97e+1,     6.36e+0,         2.00e-1 },
    { 3.29e+2,      1.29e+1,     6.94e+0,         2.80e-1 }
};
static const t_cross   cross_sec_mg[] = {
    { 7.79E+4,      7.19E+1,     1.34E+0,         -1 },
    { 3.75E+3,      3.75E+1,     3.18E+0,         -1 },
    { 1.61E+3,      2.75E+1,     3.91E+0,         -1 },
    { 8.25E+2,      2.06E+1,     4.47E+0,         -1 },
    { 4.75E+2,      1.61E+1,     4.88E+0,         -1 },
    { 2.40E+2,      1.16E+1,     5.32E+0,         -1 },
    { 9.82E+1,      7.59E+0,     5.74E+0,         -1 }
};
static const t_cross   cross_sec_al[] = {
    { 1.01E+5,      8.24E+1,     1.51E+0,         -1 },
    { 5.12E+3,      4.32E+1,     3.45E+0,         -1 },
    { 2.22E+3,      3.24E+1,     4.16E+0,         -1 },
    { 1.14E+3,      2.47E+1,     4.74E+0,         -1 },
    { 6.63E+2,      1.93E+1,     5.19E+0,         -1 },
    { 3.37E+2,      1.41E+1,     5.67E+0,         -1 },
    { 1.39E+2,      9.17E+0,     6.14E+0,         -1 }
};
static const t_cross   cross_sec_b[] = {
    { 2.86E+3,      1.05E+1,     8.20E-1,         -1 },
    { 9.38E+1,      3.76E+0,     1.92E+0,         -1 },
    { 3.72E+1,      2.81E+0,     2.15E+0,         -1 },
    { 1.80E+1,      2.20E+0,     2.31E+0,         -1 },
    { 9.92E+0,      1.77E+0,     2.44E+0,         -1 },
    { 4.77E+0,      1.32E+0,     2.58E+0,         -1 },
    { 1.84E+0,      8.56E-1,     2.71E+0,         -1 }
};
static const t_cross   cross_sec_na[] = {
    { 5.80E+4,      6.27E+1,     1.01E+0,         -1 },
    { 2.65E+3,      3.17E+1,     2.95E+0,         -1 },
    { 1.13E+3,      2.26E+1,     3.67E+0,         -1 },
    { 5.74E+2,      1.68E+1,     4.20E+0,         -1 },
    { 3.28E+2,      1.30E+1,     4.58E+0,         -1 },
    { 1.65E+2,      9.43E+0,     4.96E+0,         -1 },
    { 6.71E+1,      6.16E+0,     5.34E+0,         -1 }
};
static const t_cross   cross_sec_li[] = {
    { 3.08E+2,      3.37E+0,     6.38E-1,         -1 },
    { 8.60E+0,      1.60E+0,     1.18E+0,         -1 },
    { 3.31E+0,      1.16E+0,     1.36E+0,         -1 },
    { 1.57E+0,      8.63E-1,     1.48E+0,         -1 },
    { 8.50E-1,      6.59E-1,     1.57E+0,         -1 },
    { 4.01E-1,      4.63E-1,     1.64E+0,         -1 },
    { 1.52E-1,      2.85E-1,     1.70E+0,         -1 }
};

typedef struct {
    const char    *name;
    int            nel;
    const t_cross *cross;
} t_element;

static const t_element element[] = {
    { "H",   1, cross_sec_h },
    { "C",   6, cross_sec_c },
    { "N",   7, cross_sec_n },
    { "O",   8, cross_sec_o },
    { "S",  16, cross_sec_s },
    { "LI",  3, cross_sec_li },
    { "B",   5, cross_sec_b },
    { "NA",  11, cross_sec_na },
    { "MG",  12, cross_sec_mg },
    { "AL",  13, cross_sec_al },
    { "AR", 20, cross_sec_s }, /* This is not correct! */
    { "EL",  1, cross_sec_h }  /* This is not correct! */
};
#define NELEM asize(element)

/*
 * In the first column the binding energy of the K-electrons;
 * THIS IS IN eV,  which matches the photon energies.
 * In the second column the binding energy of the outer shell electrons
 * The third column describes the photoelectric cross sections,
 * where this now gives the fraction of photoelectric events
 * which correspond to K-shell events, I called f_j in my
 * notes:
 * The final column (a new column) now gives the values for the lifetimes
 * in ps.
 */
typedef struct {
    real E_K, E_L, Prob_K, tau;
} t_recoil;

const t_recoil recoil[] = {
    { 0.0,    0.0,   0.0,   0},
    { 0.0136, 0.0,   0.0,   0},
    { 0.0246, 0.0,   0.0,   0},
    { 0.055,  0.005, 0.960, 0.012},
    { 0.117,  0.009, 0.956, 0.012},
    { 0.192,  0.008, 0.952, 0.012},
    { 0.284,  0.011, 0.950, 0.0113},
    { 0.402,  0.015, 0.950, 0.0083},
    { 0.532,  0.014, 0.936, 0.0066},
    { 0.687,  0.017, 0.928, 0.0045},
    { 0.874,  0.031, 0.922, 0.0033},
    { 1.072,  0.041, 0.933, 0.0028},
    { 1.305,  0.054, 0.927, 0.0022},
    { 1.560,  0.077, 0.922, 0.0019},
    { 1.839,  0.105, 0.918, 0.00165},
    { 2.146,  0.133, 0.921, 0.00145},
    { 2.472,  0.166, 0.908, 0.00130},
    { 2.822,  0.212, 0.902, 0.0012},
    { 3.203,  0.247, 0.902, 0.0010},
    { 3.607,  0.298, 0.894, 0.00095},
    { 4.038,  0.348, 0.890, 0.00085},
    { 4.490,  0.404, 0.886, 0.00078},
    { 4.966,  0.458, 0.882, 0.00073},
    { 5.465,  0.516, 0.885, 0.00062},
    { 5.989,  0.578, 0.883, 0.00055},
    { 6.539,  0.645, 0.880, 0.00049},
    { 7.112,  0.713, 0.877, 0.00044}
};

#define PREFIX "IONIZE: "

enum {
    eionCYL, eionSURF, eionGAUSS, eionNR
};

enum {
    ecollPHOTO, ecollINELASTIC, ecollNR
};

typedef struct {
    int  z, n, k;
    real fj, sigPh, sigIn, vAuger;
} t_cross_atom;

/* BEGIN GLOBAL VARIABLES */

/*
   UGLY HACK
   The 2 in this list doesn't really mean 2, but 2.5 keV as
   it's checked inside the code and added 0.5 when needed.
 */

static int   Energies[] = { 2, 6, 8, 10, 12, 15, 20 };
static int   ionize_seed = 1993;
#define NENER asize(Energies)
/* END GLOBAL VARIABLES */

void dump_ca(FILE *fp, t_cross_atom *ca, int i, const char *file, int line)
{
    fprintf(fp, PREFIX "(line %d) atom %d, z = %d, n = %d, k = %d\n",
            line, i, ca->z, ca->n, ca->k);
}

t_cross_atom *mk_cross_atom(FILE *log, t_mdatoms *md,
                            gmx_mtop_t *mtop, int Eindex)
{
    int           elem_index[] = { 0, 0, 0, 5, 0, 6, 1, 2, 3, 0, 0, 7, 8, 9, 0, 0, 4, 0, 0, 0, 10, 11 };
    t_cross_atom *ca;
    int          *elemcnt;
    char         *cc, *resname;
    int           i, j, resnr;

    fprintf(log, PREFIX "Filling data structure for ionization\n");
    fprintf(log, PREFIX "Warning: all fj values set to 0.95 for now\n");
    snew(ca, md->nr);
    snew(elemcnt, NELEM+1);
    for (i = 0; (i < md->nr); i++)
    {
        ca[i].n = 0;
        ca[i].k = 0;
        /* This code does not work for domain decomposition */
        gmx_mtop_atominfo_global(mtop, i, &cc, &resnr, &resname);
        for (j = 0; (j < NELEM); j++)
        {
            if (strncmp(cc, element[j].name, strlen(element[j].name)) == 0)
            {
                ca[i].z = element[j].nel;
                break;
            }
        }
        if (j == NELEM)
        {
            gmx_fatal(FARGS, PREFIX "Don't know number of electrons for %s", cc);
        }

        elemcnt[j]++;

        ca[i].sigPh = element[elem_index[ca[i].z]].cross[Eindex].photo;
        ca[i].sigIn = element[elem_index[ca[i].z]].cross[Eindex].incoh;
        ca[i].fj    = recoil[ca[i].z].Prob_K;
        switch (ca[i].z)
        {
            case 6:
                ca[i].vAuger  = 0.904;
                break;
            case 7:
                ca[i].vAuger  = 0.920;
                break;
            case 8:
                ca[i].vAuger  = 0.929;
                break;
            case 3:    /*  probably not correct! */
                ca[i].vAuger  = 0.9;
                break;
            case 5:    /*  probably not correct! */
                ca[i].vAuger  = 0.9;
                break;
            case 11:   /*  probably not correct! */
            case 12:   /*  probably not correct! */
            case 13:   /*  probably not correct! */
            case 16:
            case 20:
                ca[i].vAuger = 1.0;
                break;
            default:
                ca[i].vAuger = -1;
        }
    }

    fprintf(log, PREFIX "You have the following elements in your system (%d atoms):\n"PREFIX, md->nr);
    for (j = 0; (j < NELEM); j++)
    {
        if (elemcnt[j] > 0)
        {
            fprintf(log, "  %s: %d", element[j].name, elemcnt[j]);
        }
    }
    fprintf(log, " atoms\n");

    sfree(elemcnt);

    return ca;
}

int number_K(t_cross_atom *ca)
{
    if (ca->z <= 2)
    {
        return ca->z-ca->n;
    }
    else
    {
        return 2-ca->k;
    }
}

int number_L(t_cross_atom *ca)
{
    return ca->k-2+ca->z-ca->n;
}

real xray_cross_section(int eColl, t_cross_atom *ca)
{
    real c = 0;
    int  nK, nL;

    switch (eColl)
    {
        case ecollPHOTO:
            nK = number_K(ca);
            nL = number_L(ca);
            if (ca->z == 1)
            {
                c = ca->sigPh;
            }
            else if (ca->z == 2)
            {
                c = ca->sigPh*0.5;
            }
            else
            {
                c = (nK*0.5*ca->fj + nL/(ca->z-2)*(1-ca->fj))*ca->sigPh;
            }
            break;
        case ecollINELASTIC:
            c = (ca->z-ca->n)*ca->sigIn/ca->z;
            break;
        default:
            gmx_fatal(FARGS, "No such collision type %d\n", eColl);
    }
    return c;
}

real prob_K(int eColl, t_cross_atom *ca)
{
    real Pl, Pk, P = 0;

    if ((ca->z <= 2) || (ca->z == ca->n))
    {
        return 0;
    }

    switch (eColl)
    {
        case ecollPHOTO:
            Pl = (ca->k-2+ca->z-ca->n)*(1-ca->fj)/(ca->z-2);
            Pk = (2-ca->k)*ca->fj*0.5;
            P  = Pk/(Pl+Pk);
            break;
        case ecollINELASTIC:
            P = (2-ca->k)/(ca->z-ca->n);
            break;
        default:
            gmx_fatal(FARGS, "No such collision type %d\n", eColl);
    }
    return P;
}

double myexp(double x)
{
    if (x < -70)
    {
        return 0.0;
    }
    else
    {
        return exp(x);
    }
}

real ptheta_incoh(int Eindex, real theta)
/* theta should be in degrees */
{
    /* These numbers generated by fitting 5 gaussians to the real function
     * that describes the probability for theta.
     * We use symmetry in the gaussian (see 180-angle) therefore there
     * are fewer parameters (only 8 per energylevel).
     */
    static double ppp[NENER][8] = {
        { -0.00295169, 10.4847, 0.0341099, /*-*/ 43.1963,
          -0.0164054,  30.2452, 71.0311,    2.50282 },
        { -0.00370852, 9.02037, 0.100559,  /*-*/ 42.9962,
          -0.0537891,  35.5077, 71.4305,    1.05515 },
        { -0.00427039, 7.86831, 0.118042,  /*-*/ 45.9846,
          -0.0634505,  38.6134, 70.3857,    0.240082 },
        { -0.004514,   7.0728,  0.13464,  /*-*/ 48.213,
          -0.0723,     41.06,   69.38,     -0.02 },
        { -0.00488796, 5.87988, 0.159574,  /*-*/ 51.5556,
          -0.0855767,  44.7307, 69.0251,   -0.414604 },
        { -0.00504604, 4.56299, 0.201064,  /*-*/ 54.8599,
          -0.107153,   48.7016, 68.8212,   -0.487699 }
    };
    double        g1, g2, g3, g4, g5, ptheta;

    g1 = myexp(-0.5*sqr((theta-ppp[Eindex][7])/ppp[Eindex][1]));
    g2 = myexp(-0.5*sqr((theta-180+ppp[Eindex][7])/ppp[Eindex][1]));
    g3 = myexp(-0.5*sqr((theta-90)/ppp[Eindex][3]));
    g4 = myexp(-0.5*sqr((theta-ppp[Eindex][6])/ppp[Eindex][5]));
    g5 = myexp(-0.5*sqr((theta-180+ppp[Eindex][6])/ppp[Eindex][5]));

    ptheta = ppp[Eindex][0]*(g1+g2) + ppp[Eindex][2]*g3 + ppp[Eindex][4]*(g4+g5);

    return ptheta;
}

real rand_theta_incoh(int Eindex, int *seed)
{
#define NINTP 450
#define prev (1-cur)
    static gmx_bool bFirst = TRUE;
    static real   **intp;
    static int      i, j, cur = 1;
    real            rrr, dx;
    real            y[2];

    dx = 90.0/(real)NINTP;
    if (bFirst)
    {
        /* Compute cumulative integrals of all probability distributions */
        snew(intp, NENER);
        for (i = 0; (i < NENER); i++)
        {
            snew(intp[i], NINTP+1);
            y[prev]    = ptheta_incoh(i, 0.0);
            /*sum        = y[prev];*/
            for (j = 1; (j <= NINTP); j++)
            {
                y[cur]     = ptheta_incoh(i, j*dx);
                /*sum       += y[cur];*/
                intp[i][j] = intp[i][j-1] + (y[cur]+y[prev])*dx;
                cur        = prev;
            }
        }
        if (debug)
        {
            fprintf(debug, "Integrated probability functions for theta incoherent\n");
            for (j = 0; (j < NINTP); j++)
            {
                fprintf(debug, "%10f", dx*j);
                for (i = 0; (i < NENER); i++)
                {
                    fprintf(debug, "  %10f", intp[i][j]);
                }
                fprintf(debug, "\n");
            }
        }
        bFirst = FALSE;
    }

    rrr = rando(seed);
    for (j = 0; (j < NINTP) && (rrr > intp[Eindex][j]); j++)
    {
        ;
    }

    return (j-1+(rrr-intp[Eindex][j-1])/(intp[Eindex][j]-intp[Eindex][j-1]))*dx;
}

static void polar2cart(real phi, real theta, rvec v)
{
    v[XX] = cos(phi)*sin(theta);
    v[YY] = sin(phi)*sin(theta);
    v[ZZ] = cos(theta);
}

void rand_vector(rvec v, int *seed)
{
    real theta, phi;

    theta = 180.0*rando(seed);
    phi   = 360.0*rando(seed);
    polar2cart(phi, theta, v);
}

gmx_bool khole_decay(FILE *fp, t_cross_atom *ca, rvec x[], rvec v[], int ion,
                     int *seed, real dt)
{
    rvec dv;
    real factor;

    if ((ca->vAuger < 0) || (recoil[ca->z].tau == 0))
    {
        dump_ca(stderr, ca, ion, __FILE__, __LINE__);
        exit(1);
    }
    if ((rando(seed) < dt/recoil[ca->z].tau) && (number_L(ca) > 1))
    {
        if (debug)
        {
            fprintf(debug, "DECAY: Going to decay a k hole\n");
        }
        ca->n++;
        ca->k--;
        /* Generate random vector */
        rand_vector(dv, seed);

        factor = ca->vAuger;
        if (debug)
        {
            fprintf(debug, "DECAY: factor=%10g, dv = (%8.3f, %8.3f, %8.3f)\n",
                    factor, dv[XX], dv[YY], dv[ZZ]);
        }
        svmul(factor, dv, dv);
        rvec_inc(v[ion], dv);

        return TRUE;
    }
    else
    {
        return FALSE;
    }
}

void ionize(FILE *fp, const output_env_t oenv, t_mdatoms *md, gmx_mtop_t *mtop,
            real t, t_inputrec *ir, rvec x[], rvec v[], int start, int end,
            matrix box, t_commrec *cr)
{
    static FILE         *xvg, *ion;
    static const char   *const_leg[] = { "Probability", "Primary Ionization", "Integral over PI", "KHole-Decay", "Integral over KD" };
    static gmx_bool      bFirst      = TRUE;
    static real          t0, imax, width, rho, nphot;
    static real          interval;
    static int           dq_tot, nkd_tot, mode, ephot;
    static t_cross_atom *ca;
    static int           Eindex    = -1;
    static gmx_rng_t     gaussrand = NULL;

    real                 factor, E_lost = 0;
    real                 pt, ptot, pphot, pcoll[ecollNR], tmax;
    real                 hboxx, hboxy, rho2;
    rvec                 dv, ddv;
    gmx_bool             bIonize = FALSE, bKHole, bL, bDOIT;
    int                  i, k, kk, m, nK, nL, dq, nkh, nkdecay;
    int                 *nionize, *nkhole, *ndecay, nbuf[2];
    char               **leg;


    if (bFirst)
    {
        /* Get parameters for gaussian photon pulse from inputrec */
        t0          = ir->userreal1;      /* Peak of the gaussian pulse            */
        nphot       = ir->userreal2;      /* Intensity                             */
        width       = ir->userreal3;      /* Width of the peak (in time)           */
        rho         = ir->userreal4;      /* Diameter of the focal spot (nm)       */
        ionize_seed = ir->userint1;       /* Random seed for stochastic ionization */
        ephot       = ir->userint2;       /* Energy of the photons                 */
        mode        = ir->userint3;       /* Mode of ionizing                      */
        interval    = 0.001*ir->userint4; /* Interval between pulses (ps)    */
        gaussrand   = gmx_rng_init(ionize_seed);

        snew(leg, asize(const_leg));
        for (i = 0; i < asize(const_leg); i++)
        {
            leg[i] = strdup(const_leg[i]);
        }

        if ((width <= 0) || (nphot <= 0))
        {
            gmx_fatal(FARGS, "Your parameters for ionization are not set properly\n"
                      "width (userreal3) = %f,  nphot (userreal2) = %f",
                      width, nphot);
        }

        if ((mode < 0) || (mode >= eionNR))
        {
            gmx_fatal(FARGS, "Ionization mode (userint3)"
                      " should be in the range 0 .. %d", eionNR-1);
        }

        switch (mode)
        {
            case eionCYL:
                imax  = (nphot/(M_PI*sqr(rho/2)))*1e-10*1.0/(width*sqrt(2.0*M_PI));
                break;
            case eionSURF:
                imax  = (nphot/(M_PI*sqr(rho/2)))*1e-10*1.0/(width*sqrt(2.0*M_PI));
                break;
        }
        if (ionize_seed == 0)
        {
            ionize_seed = make_seed();
        }
        if (PAR(cr))
        {
            for (i = 0; (i < cr->nodeid); i++)
            {
                ionize_seed = INT_MAX*rando(&ionize_seed);
            }
            fprintf(fp, PREFIX "Modifying seed on parallel processor to %d\n",
                    ionize_seed);
        }

        for (Eindex = 0; (Eindex < NENER) && (Energies[Eindex] != ephot); Eindex++)
        {
            ;
        }
        if (Eindex == NENER)
        {
            gmx_fatal(FARGS, PREFIX "Energy level of %d keV not supported", ephot);
        }

        /* Initiate cross section data etc. */
        ca      = mk_cross_atom(fp, md, mtop, Eindex);

        dq_tot  = 0;
        nkd_tot = 0;

        xvg   = xvgropen("ionize.xvg", "Ionization Events", "Time (ps)", "()", oenv);
        xvgr_legend(xvg, asize(leg), (const char**)leg, oenv);
        ion   = gmx_fio_fopen("ionize.log", "w");

        fprintf(fp, PREFIX "Parameters for ionization events:\n");
        fprintf(fp, PREFIX "Imax = %g, t0 = %g, width = %g, seed = %d\n"
                PREFIX "# Photons = %g, rho = %g, ephot = %d (keV)\n",
                imax, t0, width, ionize_seed, nphot, rho, ephot);
        fprintf(fp, PREFIX "Electron_mass: %10.3e(keV) Atomic_mass: %10.3e(keV)\n"
                PREFIX "Speed_of_light: %10.3e(nm/ps)\n",
                ELECTRONMASS_keV, ATOMICMASS_keV, SPEED_OF_LIGHT);
        fprintf(fp, PREFIX "Interval between shots: %g ps\n", interval);
        fprintf(fp, PREFIX "Eindex = %d\n", Eindex);
        fprintf(fp, PREFIX "Doing ionizations for atoms %d - %d\n", start, end);

        fflush(fp);

        bFirst = FALSE;
    }

    /******************************************************
     *
     *    H E R E    S T A R T S   I O N I Z A T I O N
     *
     ******************************************************/

    /* Calculate probability */
    tmax        = t0;
    if (interval > 0)
    {
        while (t > (tmax+interval*0.5))
        {
            tmax += interval;
        }
    }
    /*  End when t <= t0 + (N+0.5) interval */

    pt          = imax*ir->delta_t*exp(-0.5*sqr((t-tmax)/width));
    dq          = 0;
    nkdecay     = 0;

    hboxx       = 0.5*box[XX][XX];
    hboxy       = 0.5*box[YY][YY];
    rho2        = sqr(rho);

    /* Arrays for ionization statistics */
    snew(nionize, md->nr);
    snew(nkhole, md->nr);
    snew(ndecay, md->nr);

    /* Loop over atoms */
    for (i = start; (i < end); i++)
    {
        /* Loop over collision types */
        bKHole  = FALSE;
        bIonize = FALSE;
        for (k = 0; (k < ecollNR); k++)
        {
            /* Determine cross section for this collision type */
            pcoll[k] = pt*xray_cross_section(k, &(ca[i]));
        }

        /* Total probability of ionisation */
        ptot = 1 - (1-pcoll[ecollPHOTO])*(1-pcoll[ecollINELASTIC]);
        if (debug && (i == 0))
        {
            fprintf(debug, PREFIX "Ptot = %g, t = %g\n", ptot, t);
        }

        /* Check whether to ionize this guy */
        bDOIT = FALSE;
        switch (mode)
        {
            case eionCYL:
                bDOIT = (((rando(&ionize_seed) < ptot) && (ca[i].n < ca[i].z)) &&
                         ((sqr(x[i][XX] - hboxx) + sqr(x[i][YY] - hboxy)) < rho2));
                break;
            case eionSURF:
                bDOIT = FALSE;
                break;
            default:
                gmx_fatal(FARGS, "Unknown ionization mode %d (%s, line %d)", mode,
                          __FILE__, __LINE__);
        }

        if (bDOIT)
        {
            clear_rvec(dv);

            /* The relative probability for a photoellastic event is given by: */
            pphot = pcoll[ecollPHOTO]/(pcoll[ecollPHOTO]+pcoll[ecollINELASTIC]);

            if (rando(&ionize_seed) < pphot)
            {
                k = ecollPHOTO;
            }
            else
            {
                k = ecollINELASTIC;
            }

            /* If a random number is smaller than the probability for
             * an L ionization than do that. Note that the probability
             * may be zero (H, He), but the < instead of <= covers that.
             */
            nK = number_K(&ca[i]);
            nL = number_L(&ca[i]);
            bL = (nK == 0) || ( (nL > 0) && (rando(&ionize_seed) > prob_K(k, &(ca[i]))));

            switch (k)
            {
                case ecollPHOTO:
                {
                    /* Select which one to take by yet another random numer */
                    real theta, phi;

                    /* Get parameters for photoelestic effect */
                    /* Note that in the article this is called 2 theta */
                    theta = DEG2RAD*gmx_rng_gaussian_table(gaussrand)*26.0+70.0;

                    phi   = 2*M_PI*rando( &ionize_seed);

                    if (bL)
                    {
                        E_lost = ephot-recoil[ca[i].z].E_L*(ca[i].n+1);
                    }
                    if (ephot == 2)
                    {
                        E_lost += 0.5; /* Real energy should be 2.5 KeV*/
                    }
                    else
                    {
                        E_lost = ephot-recoil[ca[i].z].E_K;
                        if (ephot == 2)
                        {
                            E_lost += 0.5; /* Real energy should be 2.5 KeV*/
                        }
                        if ((ca[i].z > 2) && (nL > 0))
                        {
                            bKHole = TRUE;
                        }
                    }
                    if (debug)
                    {
                        fprintf(debug, "i = %d, nK = %d, nL = %d, bL = %s, bKHole = %s\n",
                                i, nK, nL, EBOOL(bL), EBOOL(bKHole));
                    }
                    if (E_lost < 0)
                    {
                        E_lost  = 0.0;
                        bIonize = FALSE;
                        bKHole  = FALSE;
                    }
                    else
                    {
                        /* Compute the components of the velocity vector */
                        factor = ((ELECTRONMASS_keV/(ATOMICMASS_keV*md->massT[i]))*
                                  (SPEED_OF_LIGHT*sqrt(2*E_lost/ELECTRONMASS_keV)));

                        /* Subtract momentum of recoiling electron */
                        polar2cart(phi, theta, ddv);
                        for (m = 0; (m < DIM); m++)
                        {
                            dv[m] -= factor*ddv[m];
                        }

                        if (debug)
                        {
                            pr_rvec(debug, 0, "ELL", dv, DIM, TRUE);
                        }

                        bIonize = TRUE;
                    }
                    break;
                }
                case ecollINELASTIC:
                {
                    real theta, Ebind, Eelec;

                    if (bL)
                    {
                        Ebind = (ca[i].n+1)*recoil[ca[i].z].E_L;
                    }
                    else
                    {
                        Ebind  = recoil[ca[i].z].E_K;
                        if ((ca[i].z > 2) && (nL > 0))
                        {
                            bKHole = TRUE;
                        }
                    }
                    theta      = DEG2RAD*rand_theta_incoh(Eindex, &ionize_seed);
                    Eelec      = (sqr(ephot)/512)*(1-cos(2*theta));
                    if (ephot == 2)
                    {
                        Eelec = (sqr(ephot+.5)/512)*(1-cos(2*theta)); /* Real energy should be 2.5 KeV*/
                    }
                    bIonize    = (Ebind <= Eelec);
                    bKHole     = bKHole && bIonize;
                    if (debug)
                    {
                        fprintf(debug, PREFIX "Ebind: %g, Eelectron: %g\n", Ebind, Eelec);
                    }
                    if (!bIonize)
                    {
                        /* Subtract momentum of recoiling photon */
                        /*phi     = 2*M_PI*rando(&ionize_seed);
                           bKHole  = FALSE;
                           factor  = ephot*438;
                           dv[XX] -= factor*cos(phi)*sin(theta);
                           dv[YY] -= factor*sin(phi)*sin(theta);
                           dv[ZZ] -= factor*cos(theta);
                         */
                        if (debug)
                        {
                            pr_rvec(debug, 0, "INELL", dv, DIM, TRUE);
                        }
                    }
                    break;
                }
                default:
                    gmx_fatal(FARGS, "Ga direct naar de gevangenis. Ga niet langs start");
            }
            if (bIonize)
            {
                /* First increase the charge */
                if (ca[i].n < ca[i].z)
                {
                    md->chargeA[i] += 1.0;
                    md->chargeB[i] += 1.0;
                    ca[i].n++;
                    dq++;
                }
                if (debug)
                {
                    fprintf(debug, "Random-dv[%3d] = %10.3e,%10.3e,%10.3e,"
                            " ephot = %d, Elost=%10.3e\n",
                            i, dv[XX], dv[YY], dv[ZZ], ephot, E_lost);
                }
            }
            /* Now actually add the impulse to the velocities */
            for (m = 0; (m < DIM); m++)
            {
                v[i][m] += dv[m];
            }
            if (bKHole)
            {
                ca[i].k++;
                nkhole[i]++;
            }
            else if (bIonize)
            {
                nionize[i]++;
            }
        }

        /* Now check old event: Loop over k holes! */
        nkh = ca[i].k;
        for (kk = 0; (kk < nkh); kk++)
        {
            if (khole_decay(fp, &(ca[i]), x, v, i, &ionize_seed, ir->delta_t))
            {
                nkdecay++;
                ndecay[i]++;
                md->chargeA[i] += 1.0;
                md->chargeB[i] += 1.0;
            }
        }

        if (debug && (ca[i].n > 0))
        {
            dump_ca(debug, &(ca[i]), i, __FILE__, __LINE__);
        }
    }

    /* Sum events for statistics if necessary */
    if (PAR(cr))
    {
        gmx_sumi(md->nr, nionize, cr);
        gmx_sumi(md->nr, nkhole, cr);
        gmx_sumi(md->nr, ndecay, cr);
        nbuf[0] = dq; nbuf[1] = nkdecay;
        gmx_sumi(2, nbuf, cr);
        dq = nbuf[0]; nkdecay = nbuf[1];
    }
    /* Now sum global events on this timestep to cumulative numbers */
    dq_tot  += dq;
    nkd_tot += nkdecay;

    /* Printing time */
    if (MASTER(cr))
    {
        /* Print data to the file that holds ionization events per atom */
        fprintf(ion, "%12.8f", t);
        for (i = 0; (i < md->nr); i++)
        {
            if (nionize[i])
            {
                fprintf(ion, "  I:%d", i+1);
            }
            if (nkhole[i])
            {
                fprintf(ion, "  K:%d", i+1);
            }
            if (ndecay[i])
            {
                fprintf(ion, "  D:%d", i+1);
            }
        }
        fprintf(ion, "\n");
        if (debug)
        {
            fflush(ion);
        }

        /* Print statictics to file */
        fprintf(xvg, "%10.5f  %10.3e  %6d  %6d  %6d  %6d",
                t, pt, dq, dq_tot, nkdecay, nkd_tot);
        fprintf(xvg, "\n");
        if (debug)
        {
            fflush(xvg);
        }
    }
    sfree(nionize);
    sfree(nkhole);
    sfree(ndecay);
}
