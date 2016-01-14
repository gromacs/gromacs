/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/copyrite.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

// Alexandria stuff
#include "molprop_util.h"
#include "mymol.h"
#include "poldata_xml.h"
#include "stringutil.h"

typedef struct {
    char       *a1, *a2;
    int         order;
    int         nhisto;
    int        *histo;
    gmx_stats_t lsq;
} t_bond;

typedef struct {
    char       *a1, *a2, *a3;
    int         nhisto;
    int        *histo;
    gmx_stats_t lsq;
} t_angle;

typedef struct {
    char       *a1, *a2, *a3, *a4;
    int         nhisto;
    int        *histo;
    gmx_stats_t lsq;
} t_dih;

typedef struct {
    int      nbond;
    t_bond  *bond;
    int      nangle;
    t_angle *angle;
    int      ndih, nimp;
    t_dih   *dih, *imp;
} t_bonds;

int bondcmp(const void *a, const void *b)
{
    t_bond *ba = (t_bond *)a;
    t_bond *bb = (t_bond *)b;
    int     d;

    if ((d = strcmp(ba->a1, bb->a1)) == 0)
    {
        if ((d = strcmp(ba->a2, bb->a2)) == 0)
        {
            d = ba->order-bb->order;
        }
    }
    return d;
}

int anglecmp(const void *a, const void *b)
{
    t_angle *ba = (t_angle *)a;
    t_angle *bb = (t_angle *)b;
    int      d;

    if ((d = strcmp(ba->a1, bb->a1)) == 0)
    {
        if ((d = strcmp(ba->a2, bb->a2)) == 0)
        {
            return strcmp(ba->a3, bb->a3);
        }
    }

    return d;
}

int dihcmp(const void *a, const void *b)
{
    t_dih *ba = (t_dih *)a;
    t_dih *bb = (t_dih *)b;
    int    d;

    if ((d = strcmp(ba->a1, bb->a1)) == 0)
    {
        if ((d = strcmp(ba->a2, bb->a2)) == 0)
        {
            if ((d = strcmp(ba->a3, bb->a3)) == 0)
            {
                return strcmp(ba->a4, bb->a4);
            }
        }
    }

    return d;
}

void sort_bonds(t_bonds *b)
{
    qsort(b->bond, b->nbond, sizeof(b->bond[0]), bondcmp);
    qsort(b->angle, b->nangle, sizeof(b->angle[0]), anglecmp);
    qsort(b->dih, b->ndih, sizeof(b->dih[0]), dihcmp);
    qsort(b->imp, b->nimp, sizeof(b->imp[0]), dihcmp);
}

void add_bond(FILE *fplog, const char *molname, t_bonds *b, const char *a1, const char *a2,
              double blen, double spacing, int order)
{
    int i, j, index;

    index = std::lround(blen/spacing);
    for (i = 0; (i < b->nbond); i++)
    {
        if (((((strcmp(a1, b->bond[i].a1) == 0) && (strcmp(a2, b->bond[i].a2)) == 0)) ||
             (((strcmp(a2, b->bond[i].a1) == 0) && (strcmp(a1, b->bond[i].a2)) == 0))) &&
            (b->bond[i].order == order))
        {
            break;
        }
    }
    if (i == b->nbond)
    {
        b->nbond++;
        srenew(b->bond, b->nbond);
        if (strcmp(a1, a2) < 0)
        {
            b->bond[i].a1     = strdup(a1);
            b->bond[i].a2     = strdup(a2);
        }
        else
        {
            b->bond[i].a1     = strdup(a2);
            b->bond[i].a2     = strdup(a1);
        }
        b->bond[i].order  = order;
        b->bond[i].nhisto = 2*index+1;
        snew(b->bond[i].histo, b->bond[i].nhisto);
        b->bond[i].lsq = gmx_stats_init();
    }
    if (index >= b->bond[i].nhisto)
    {
        srenew(b->bond[i].histo, index+100);
        for (j = b->bond[i].nhisto; (j < index+100); j++)
        {
            b->bond[i].histo[j] = 0;
        }
        b->bond[i].nhisto = index+100;
    }
    gmx_stats_add_point(b->bond[i].lsq, 0, blen, 0, 0);
    b->bond[i].histo[index]++;
    if (NULL != fplog)
    {
        fprintf(fplog, "%s bond-%s-%s-%d %g\n", molname,
                b->bond[i].a1, b->bond[i].a2, order, blen);
    }
}

void add_angle(FILE *fplog, const char *molname, t_bonds *b,
               const char *a1, const char *a2, const char *a3, double angle, double spacing)
{
    int i, index;

    index = std::lround(angle/spacing);
    for (i = 0; (i < b->nangle); i++)
    {
        if ((strcmp(a2, b->angle[i].a2) == 0) &&
            (((strcmp(a1, b->angle[i].a1) == 0) && (strcmp(a3, b->angle[i].a3) == 0)) ||
             ((strcmp(a3, b->angle[i].a1) == 0) && (strcmp(a1, b->angle[i].a3) == 0))))
        {
            break;
        }
    }
    if (i == b->nangle)
    {
        b->nangle++;
        srenew(b->angle, b->nangle);
        b->angle[i].a2     = strdup(a2);
        if (strcmp(a1, a3) < 0)
        {
            b->angle[i].a1     = strdup(a1);
            b->angle[i].a3     = strdup(a3);
        }
        else
        {
            b->angle[i].a1     = strdup(a3);
            b->angle[i].a3     = strdup(a1);
        }
        b->angle[i].nhisto = (int) (180/spacing) + 1;
        snew(b->angle[i].histo, b->angle[i].nhisto);
        b->angle[i].lsq = gmx_stats_init();
    }
    gmx_stats_add_point(b->angle[i].lsq, 0, angle, 0, 0);
    b->angle[i].histo[index]++;
    if (NULL != fplog)
    {
        fprintf(fplog, "%s angle-%s-%s-%s %g\n", molname,
                b->angle[i].a1, b->angle[i].a2, b->angle[i].a3, angle);
    }
}

void add_dih(FILE *fplog, const char *molname, t_bonds *b,
             const char *a1, const char *a2, const char *a3, const char *a4,
             double angle, double spacing, int egd)
{
    int     i, *nd, index;
    t_dih **ddd;

    if (angle < 0)
    {
        angle += 360;
    }
    if (egd == egdPDIHS)
    {
        ddd = &b->dih;
        nd  = &b->ndih;
    }
    else
    {
        ddd = &b->imp;
        nd  = &b->nimp;
        while (angle > 176)
        {
            angle -= 180;
        }
    }
    index = std::lround(angle/spacing);
    if (index < 0)
    {
        index = 0;
    }
    for (i = 0; (i < (*nd)); i++)
    {
        if (((strcmp(a1, (*ddd)[i].a1) == 0) && (strcmp(a2, (*ddd)[i].a2) == 0) &&
             (strcmp(a3, (*ddd)[i].a3) == 0) && (strcmp(a4, (*ddd)[i].a4) == 0)) ||
            ((strcmp(a1, (*ddd)[i].a4) == 0) && (strcmp(a2, (*ddd)[i].a3) == 0) &&
             (strcmp(a3, (*ddd)[i].a2) == 0) && (strcmp(a4, (*ddd)[i].a1) == 0)))
        {
            break;
        }
    }
    if (i == (*nd))
    {
        (*nd)++;
        srenew((*ddd), (*nd));
        if (strcmp(a1, a4) < 0)
        {
            (*ddd)[i].a1     = strdup(a1);
            (*ddd)[i].a2     = strdup(a2);
            (*ddd)[i].a3     = strdup(a3);
            (*ddd)[i].a4     = strdup(a4);
        }
        else
        {
            (*ddd)[i].a4     = strdup(a1);
            (*ddd)[i].a3     = strdup(a2);
            (*ddd)[i].a2     = strdup(a3);
            (*ddd)[i].a1     = strdup(a4);
        }
        (*ddd)[i].nhisto = (int) (360/spacing) + 1;
        snew((*ddd)[i].histo, (*ddd)[i].nhisto);
        (*ddd)[i].lsq = gmx_stats_init();
    }
    gmx_stats_add_point((*ddd)[i].lsq, 0, angle, 0, 0);
    (*ddd)[i].histo[index]++;
    if (NULL != fplog)
    {
        fprintf(fplog, "%s %s-%s-%s-%s-%s %g\n", molname, (egd == egdPDIHS) ? "dih" : "imp",
                (*ddd)[i].a1, (*ddd)[i].a2, (*ddd)[i].a3, (*ddd)[i].a4, angle);
    }
}

static void lo_dump_histo(char *fn, char *xaxis, const gmx_output_env_t *oenv, int Nsample,
                          int n, int histo[], double spacing)
{
    FILE  *fp;
    int    j, j0, j1;
    double sum;

    for (j0 = 0; (j0 < n) && (histo[j0] == 0); j0++)
    {
        ;
    }
    j0 = std::max(j0-1, 0);
    for (j1 = n-1; (j1 > 0) && (histo[j1] == 0); j1--)
    {
        ;
    }
    j1  = std::min(j1+1, n-1);
    sum = 0;
    for (j = j0; (j <= j1); j++)
    {
        sum += histo[j];
    }
    if (sum > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s N = %d", fn, Nsample);
        fp = xvgropen(fn, buf, xaxis, "P (a.u.)", oenv);
        for (j = j0; (j <= j1); j++)
        {
            fprintf(fp, "%g  %g\n", spacing*j, histo[j]/sum);
        }
        fclose(fp);
    }
}

void dump_histo(t_bonds *b, double bspacing, double aspacing, const gmx_output_env_t *oenv)
{
    int  i, N;
    char buf[256];

    for (i = 0; (i < b->nbond); i++)
    {
        if ((gmx_stats_get_npoints(b->bond[i].lsq, &N) == 0) &&
            (b->bond[i].nhisto > 0))
        {
            sprintf(buf, "bond-%s-%s-%d.xvg", b->bond[i].a1, b->bond[i].a2, b->bond[i].order);
            lo_dump_histo(buf, (char *)"Distance (pm)", oenv, N,
                          b->bond[i].nhisto, b->bond[i].histo, bspacing);
        }
    }
    for (i = 0; (i < b->nangle); i++)
    {
        if ((gmx_stats_get_npoints(b->angle[i].lsq, &N) == 0) &&
            (b->angle[i].nhisto > 0))
        {
            sprintf(buf, "angle-%s-%s-%s.xvg",
                    b->angle[i].a1, b->angle[i].a2, b->angle[i].a3);
            lo_dump_histo(buf, (char *)"Angle (deg.)", oenv, N,
                          b->angle[i].nhisto, b->angle[i].histo, aspacing);
        }
    }
    for (i = 0; (i < b->ndih); i++)
    {
        if ((gmx_stats_get_npoints(b->dih[i].lsq, &N) == 0)  &&
            (b->dih[i].nhisto > 0))
        {
            sprintf(buf, "dih-%s-%s-%s-%s.xvg",
                    b->dih[i].a1, b->dih[i].a2, b->dih[i].a3, b->dih[i].a4);
            lo_dump_histo(buf, (char *)"Dihedral angle (deg.)", oenv, N,
                          b->dih[i].nhisto, b->dih[i].histo, aspacing);
        }
    }
    for (i = 0; (i < b->nimp); i++)
    {
        if ((gmx_stats_get_npoints(b->imp[i].lsq, &N) == 0)  &&
            (b->imp[i].nhisto > 0))
        {
            sprintf(buf, "imp-%s-%s-%s-%s.xvg",
                    b->imp[i].a1, b->imp[i].a2, b->imp[i].a3, b->imp[i].a4);
            lo_dump_histo(buf, (char *)"Improper angle (deg.)", oenv, N,
                          b->imp[i].nhisto, b->imp[i].histo, aspacing);
        }
    }
}

static void round_numbers(real *av, real *sig)
{
    *av  = ((int)(*av*100))/100.0;
    *sig = ((int)(*sig*100+50))/100.0;
}

void update_pd(FILE *fp, t_bonds *b, Poldata * pd,
               real Dm, real beta, real kt, real klin, real kp)
{
    int    i, N;
    real   av, sig;
    char   pbuf[256];
    double bondorder;

    pd->setLengthUnit( unit2string(eg2cPm));
    for (i = 0; (i < b->nbond); i++)
    {
        gmx_stats_get_average(b->bond[i].lsq, &av);
        gmx_stats_get_sigma(b->bond[i].lsq, &sig);
        gmx_stats_get_npoints(b->bond[i].lsq, &N);
        sprintf(pbuf, "%g  %g", Dm, beta);
        bondorder = b->bond[i].order;
        // Rounding the numbers to 1/10 pm and 1/10 degree
        round_numbers(&av, &sig);
        pd->addBond( b->bond[i].a1, b->bond[i].a2, av, sig, N, bondorder, pbuf);
        fprintf(fp, "bond-%s-%s len %g sigma %g (pm) N = %d%s\n",
                b->bond[i].a1, b->bond[i].a2, av, sig, N,
                (sig > 1.5) ? " WARNING" : "");

    }
    for (i = 0; (i < b->nangle); i++)
    {
        gmx_stats_get_average(b->angle[i].lsq, &av);
        gmx_stats_get_sigma(b->angle[i].lsq, &sig);
        gmx_stats_get_npoints(b->angle[i].lsq, &N);
        // Rounding the numbers to 1/10 pm and 1/10 degree
        round_numbers(&av, &sig);
        if ((av > 175) || (av < 5))
        {
            sprintf(pbuf, "%g", kt);
        }
        else
        {
            sprintf(pbuf, "%g", klin);
        }
        pd->addAngle(
                b->angle[i].a1, b->angle[i].a2,
                b->angle[i].a3, av, sig, N, pbuf);
        fprintf(fp, "angle-%s-%s-%s angle %g sigma %g (deg) N = %d%s\n",
                b->angle[i].a1, b->angle[i].a2, b->angle[i].a3, av, sig, N,
                (sig > 3) ? " WARNING" : "");
    }
    for (i = 0; (i < b->ndih); i++)
    {
        gmx_stats_get_average(b->dih[i].lsq, &av);
        gmx_stats_get_sigma(b->dih[i].lsq, &sig);
        gmx_stats_get_npoints(b->dih[i].lsq, &N);
        sprintf(pbuf, "%g  1", kp);
        // Rounding the numbers to 1/10 pm and 1/10 degree
        round_numbers(&av, &sig);
        pd->addDihedral( egdPDIHS,
                         b->dih[i].a1, b->dih[i].a2,
                         b->dih[i].a3, b->dih[i].a4, av, sig, N, pbuf);
        fprintf(fp, "dihedral-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                b->dih[i].a1, b->dih[i].a2, b->dih[i].a3, b->dih[i].a4, av, sig);
    }
    for (i = 0; (i < b->nimp); i++)
    {
        gmx_stats_get_average(b->imp[i].lsq, &av);
        gmx_stats_get_sigma(b->imp[i].lsq, &sig);
        gmx_stats_get_npoints(b->imp[i].lsq, &N);
        sprintf(pbuf, "%g  1", kp);
        // Rounding the numbers to 1/10 pm and 1/10 degree
        round_numbers(&av, &sig);
        pd->addDihedral( egdIDIHS,
                         b->imp[i].a1, b->imp[i].a2,
                         b->imp[i].a3, b->imp[i].a4, av, sig, N, pbuf);
        fprintf(fp, "improper-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                b->imp[i].a1, b->imp[i].a2, b->imp[i].a3, b->imp[i].a4, av, sig);
    }
}

int alex_bastat(int argc, char *argv[])
{
    static const char               *desc[] = {
        "bastat read a series of molecules and extracts average geometries from",
        "those. First atomtypes are determined and then bond-lengths, bond-angles",
        "and dihedral angles are extracted. The results are stored in a gentop.dat file."
    };

    t_filenm                         fnm[] = {
        { efDAT, "-f", "allmols",    ffRDMULT  },
        { efDAT, "-d", "gentop",     ffOPTRD },
        { efDAT, "-o", "bastat",    ffWRITE },
        { efDAT, "-sel", "molselect", ffREAD },
        { efLOG, "-g", "bastat",    ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])
    static int                       compress       = 0;
    static int                       maxwarn        = 0;
    static gmx_bool                  bHisto         = FALSE, bBondOrder = TRUE, bDih = FALSE;
    static real                      Dm             = 400, kt = 400, kp = 5, beta = 20, klin = 20;
    static char                     *lot            = (char *)"B3LYP/aug-cc-pVTZ";
    static const char               *cqdist[]       = {
        NULL, "AXp", "AXs", "AXg", "Yang", "Bultinck", "Rappe", NULL
    };
    t_pargs                          pa[]     = {
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-Dm",    FALSE, etREAL, {&Dm},
          "Dissociation energy (kJ/mol)" },
        { "-beta",    FALSE, etREAL, {&beta},
          "Steepness of the Morse potential (1/nm)" },
        { "-kt",    FALSE, etREAL, {&kt},
          "Angle force constant (kJ/mol/rad^2)" },
        { "-klin",  FALSE, etREAL, {&klin},
          "Linear angle force constant (kJ/mol/nm^2)" },
        { "-kp",    FALSE, etREAL, {&kp},
          "Dihedral angle force constant (kJ/mol/rad^2)" },
        { "-dih",   FALSE, etBOOL, {&bDih},
          "Generate proper dihedral terms" },
        { "-histo", FALSE, etBOOL, {&bHisto},
          "Print (hundreds of) xvg files containing histograms for bonds, angles and dihedrals" },
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Charge distribution used" },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" },
        { "-bondorder", FALSE, etBOOL, {&bBondOrder},
          "Make separate bonds for different bond orders" }
    };
    //alexandria::MolDip    md;
    FILE                            *fp;
    ChargeDistributionModel          iDistributionModel;
    gmx_molselect  *                 gms;
    time_t                           my_t;
    t_bonds                         *b;
    rvec                             dx, dx2, r_ij, r_kj, r_kl, mm, nn;
    t_pbc                            pbc;
    int                              t1, t2, t3;
    int                              ftb, fta, ftd, fti;
    matrix                           box;
    real                             sign;
    double                           bspacing = 1;   /* pm */
    double                           aspacing = 0.5; /* degree */
    double                           dspacing = 1;   /* degree */
    gmx_output_env_t                *oenv     = NULL;
    Poldata                         *pd;
    gmx_atomprop_t                   aps;
    int                              nfiles;
    char                           **fns;

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW,
                           NFILE, fnm,
                           sizeof(pa)/sizeof(pa[0]), pa,
                           sizeof(desc)/sizeof(desc[0]), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    iDistributionModel = Poldata::name2eemtype(cqdist[0]);

    fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# alexandria is part of G R O M A C S:\n#\n");
    fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());

    gms = gmx_molselect_init(opt2fn("-sel", NFILE, fnm));

    /* Read standard atom properties */
    aps = gmx_atomprop_init();

    /* Read polarization stuff */
    if ((pd = alexandria::PoldataXml::read(opt2fn_null("-d", NFILE, fnm), aps)) == NULL)
    {
        gmx_fatal(FARGS, "Can not read the force field information. File missing or incorrect.");
    }

    /* read Molprops */
    nfiles = opt2fns(&fns, "-f", NFILE, fnm);
    std::vector<alexandria::MolProp> mp;
    int nwarn = merge_xml(nfiles, fns, mp, NULL, NULL, NULL, aps, pd, TRUE);
    if (nwarn > maxwarn)
    {
        printf("Too many warnings (%d). Terminating.\n", nwarn);
        return 0;
    }
    ftb = F_BONDS;
    fta = F_ANGLES;
    ftd = F_PDIHS;
    fti = F_IDIHS;
    snew(b, 1);
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if (gmx_molselect_status(gms, mpi->getIupac().c_str()) == imsTrain)
        {
            alexandria::MyMol mmi;
            int               i;
            mmi.molProp()->Merge(mpi);
            if (mmi.molProp()->getMolname().size() == 0)
            {
                printf("Empty molname for molecule with formula %s\n", mmi.molProp()->formula().c_str());
                continue;
            }
            immStatus imm = mmi.GenerateTopology(aps, pd, lot, iDistributionModel, 2,
                                                 false, false, bDih);

            if (immOK != imm)
            {
                if (NULL != debug)
                {
                    fprintf(debug, "Could not make topology for %s\n", mmi.molProp()->getMolname().c_str());
                }
                continue;
            }
#define BTP(ii) pd->atypeToBtype( *mmi.topology_->atoms.atomtype[ii])
            for (i = 0; (i < mmi.topology_->atoms.nr); i++)
            {
                if (0 == BTP(i).size())
                {
                    if (NULL != debug)
                    {
                        fprintf(debug, "No bond-type support for atom %s in %s\n",
                                *mmi.topology_->atoms.atomtype[i], mmi.molProp()->getMolname().c_str());
                    }
                    break;
                }
            }
            if ((mmi.topology_->atoms.nr <= 0) || (i < mmi.topology_->atoms.nr))
            {
                continue;
            }
            for (int j = 0; (j < mmi.ltop_->idef.il[ftb].nr); j += interaction_function[ftb].nratoms+1)
            {
                int         ai = mmi.ltop_->idef.il[ftb].iatoms[j+1];
                int         aj = mmi.ltop_->idef.il[ftb].iatoms[j+2];
                rvec_sub(mmi.x_[ai], mmi.x_[aj], dx);
                std::string cai = BTP(ai);
                std::string caj = BTP(aj);
                if ((0 != cai.size()) && (0 != caj.size()))
                {
                    for (alexandria::BondIterator bi = mmi.molProp()->BeginBond(); (bi < mmi.molProp()->EndBond()); bi++)
                    {
                        int xi, xj, xb;
                        bi->get(&xi, &xj, &xb);
                        xi--;
                        xj--;
                        if (!bBondOrder)
                        {
                            xb = 1;
                        }
                        if (((xi == ai) && (xj == aj)) || ((xj == ai) && (xi == aj)))
                        {
                            add_bond(fp, mmi.molProp()->getMolname().c_str(), b, cai.c_str(), caj.c_str(), 1000*norm(dx),
                                     bspacing, xb);
                            break;
                        }
                    }
                }
                else
                {
                    fprintf(stderr, "No bond_atom type for either %s or %s\n",
                            BTP(ai).c_str(), BTP(aj).c_str());
                }
            }
            for (int j = 0; (j < mmi.ltop_->idef.il[fta].nr); j += interaction_function[fta].nratoms+1)
            {
                int ai = mmi.ltop_->idef.il[fta].iatoms[j+1];
                int aj = mmi.ltop_->idef.il[fta].iatoms[j+2];
                int ak = mmi.ltop_->idef.il[fta].iatoms[j+3];
                rvec_sub(mmi.x_[ai], mmi.x_[aj], dx);
                rvec_sub(mmi.x_[ak], mmi.x_[aj], dx2);
                double      ang = RAD2DEG*gmx_angle(dx, dx2);
                std::string cai = BTP(ai);
                std::string caj = BTP(aj);
                std::string cak = BTP(ak);
                if ((0 != cai.size()) && (0 != caj.size()) && (0 != cak.size()))
                {
                    add_angle(fp, mmi.molProp()->getMolname().c_str(), b, cai.c_str(), caj.c_str(), cak.c_str(), ang, aspacing);
                }
                else
                {
                    fprintf(stderr, "No bond_atom type for either %s, %s or %s\n",
                            BTP(ai).c_str(), BTP(aj).c_str(), BTP(ak).c_str());
                }
            }
            clear_mat(box);
            set_pbc(&pbc, epbcNONE, box);

            for (int j = 0; (j < mmi.ltop_->idef.il[ftd].nr); j += interaction_function[ftd].nratoms+1)
            {
                int    ai  = mmi.ltop_->idef.il[ftd].iatoms[j+1];
                int    aj  = mmi.ltop_->idef.il[ftd].iatoms[j+2];
                int    ak  = mmi.ltop_->idef.il[ftd].iatoms[j+3];
                int    al  = mmi.ltop_->idef.il[ftd].iatoms[j+4];
                double ang = RAD2DEG*dih_angle(mmi.x_[ai], mmi.x_[aj],
                                               mmi.x_[ak], mmi.x_[al],
                                               &pbc, r_ij, r_kj, r_kl, mm, nn, /* out */
                                               &sign, &t1, &t2, &t3);
                std::string cai = BTP(ai);
                std::string caj = BTP(aj);
                std::string cak = BTP(ak);
                std::string cal = BTP(al);
                if ((0 != cai.size()) && (0 != caj.size()) && (0 != cak.size()) && (0 != cal.size()))
                {
                    add_dih(fp, mmi.molProp()->getMolname().c_str(), b, cai.c_str(), caj.c_str(), cak.c_str(), cal.c_str(), ang, dspacing, egdPDIHS);
                }
                else
                {
                    fprintf(stderr, "No bond_atom type for either %s, %s, %s or %s\n",
                            BTP(ai).c_str(), BTP(aj).c_str(), BTP(ak).c_str(), BTP(al).c_str());
                }
            }

            for (int j = 0; (j < mmi.ltop_->idef.il[fti].nr); j += interaction_function[fti].nratoms+1)
            {
                int    ai  = mmi.ltop_->idef.il[fti].iatoms[j+1];
                int    aj  = mmi.ltop_->idef.il[fti].iatoms[j+2];
                int    ak  = mmi.ltop_->idef.il[fti].iatoms[j+3];
                int    al  = mmi.ltop_->idef.il[fti].iatoms[j+4];
                double ang = RAD2DEG*dih_angle(mmi.x_[ai], mmi.x_[aj],
                                               mmi.x_[ak], mmi.x_[al],
                                               &pbc, r_ij, r_kj, r_kl, mm, nn, /* out */
                                               &sign, &t1, &t2, &t3);
                std::string cai = BTP(ai);
                std::string caj = BTP(aj);
                std::string cak = BTP(ak);
                std::string cal = BTP(al);
                if ((0 != cai.size()) && (0 != caj.size()) && (0 != cak.size()) && (0 != cal.size()))
                {
                    add_dih(fp, mmi.molProp()->getMolname().c_str(), b, cai.c_str(), caj.c_str(), cak.c_str(), cal.c_str(), ang, dspacing, egdIDIHS);
                }
                else
                {
                    fprintf(stderr, "No bond_atom type for either %s, %s, %s or %s\n",
                            BTP(ai).c_str(), BTP(aj).c_str(), BTP(ak).c_str(), BTP(al).c_str());
                }
            }
        }
    }
    sort_bonds(b);
    if (bHisto)
    {
        dump_histo(b, bspacing, aspacing, oenv);
    }
    update_pd(fp, b, pd, Dm, beta, kt, klin, kp);

    alexandria::PoldataXml::write(opt2fn("-o", NFILE, fnm), pd, compress);

    printf("Extracted %d bondtypes, %d angletypes, %d dihedraltypes and %d impropertypes.\n",
           b->nbond, b->nangle, b->ndih, b->nimp);

    gmx_ffclose(fp);

    return 0;
}
