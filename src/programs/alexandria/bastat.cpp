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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
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
    std::string      a1, a2;
    int              order;
    std::vector<int> histo;
    gmx_stats_t      lsq;
} t_bond;

typedef struct {
    std::string      a1, a2, a3;
    std::vector<int> histo;
    gmx_stats_t      lsq;
} t_angle;

typedef struct {
    std::string      a1, a2, a3, a4;
    std::vector<int> histo;
    gmx_stats_t      lsq;
} t_dih;

typedef struct {
    std::vector<t_bond>  bond;
    std::vector<t_angle> angle;
    std::vector<t_dih>   dih;
    std::vector<t_dih>   imp;
} t_bonds;

static void sort_dihs(std::vector<t_dih> &dih)
{
    std::sort(dih.begin(), dih.end(),
              [](const t_dih &a, const t_dih &b)
              {
                  int d = a.a1.compare(b.a1);
                  if (0 == d)
                  {
                      if ((d = a.a2.compare(b.a2)) == 0)
                      {
                          if ((d = a.a3.compare(b.a3)) == 0)
                          {
                              return a.a4.compare(b.a4);
                          }
                      }
                  }
                  return d;
              });
}

static void sort_bonds(t_bonds *b)
{
    std::sort(b->bond.begin(), b->bond.end(),
              [](const t_bond &a, const t_bond &b)
              { 
                  int d = a.a1.compare(b.a1);
                  if (d == 0)
                  {
                      if ((d = a.a2.compare(b.a2)) == 0)
                      {
                          d = a.order-b.order;
                      }
                  }
                  return d;
              });
    std::sort(b->angle.begin(), b->angle.end(),
              [](const t_angle &a, const t_angle &b)
              {
                  int d = a.a1.compare(b.a1);
                  if (0 == d)
                  {
                      if ((d = a.a2.compare(b.a2)) == 0)
                      {
                          return a.a3.compare(b.a3);
                      }
                  }
                  return d;
              });
    sort_dihs(b->dih);
    sort_dihs(b->imp);
}

void add_bond(FILE *fplog, const char *molname, t_bonds *bonds, 
              const char *a1, const char *a2,
              double blen, double spacing, int order)
{
    GMX_RELEASE_ASSERT(strlen(a1) > 0, "atom name a1 is empty");
    GMX_RELEASE_ASSERT(strlen(a2) > 0, "atom name a2 is empty");
    size_t index = std::lround(blen/spacing);
    auto b       = std::find_if(bonds->bond.begin(), bonds->bond.end(),
                                [a1, a2, order](t_bond b)
                                {
                                    return ((((b.a1.compare(a1) == 0) && (b.a2.compare(a2) == 0)) ||
                                             ((b.a1.compare(a2) == 0) && (b.a2.compare(a1) == 0))) &&
                                            (b.order == order));
                                });
    if (b == bonds->bond.end())
    {
        t_bond bb;
        if (strcmp(a1, a2) < 0)
        {
            bb.a1 = a1;
            bb.a2 = a2;
        }
        else
        {
            bb.a1 = a2;
            bb.a2 = a1;
        }
        bb.order  = order;
        bb.histo.resize(2*index+1, 0);
        bb.lsq    = gmx_stats_init();
        bonds->bond.push_back(bb);
        b = bonds->bond.end()-1;
    }
    if (index >= b->histo.size())
    {
        b->histo.resize(index+100, 0);
    }
    gmx_stats_add_point(b->lsq, 0, blen, 0, 0);
    b->histo[index]++;
    if (NULL != fplog)
    {
        fprintf(fplog, "%s bond-%s-%s-%d %g\n", molname,
                b->a1.c_str(), b->a2.c_str(), order, blen);
    }
}

void add_angle(FILE *fplog, const char *molname, t_bonds *bonds,
               const char *a1, const char *a2, const char *a3, 
               double angle, double spacing)
{
    GMX_RELEASE_ASSERT(strlen(a1) > 0, "atom name a1 is empty");
    GMX_RELEASE_ASSERT(strlen(a2) > 0, "atom name a2 is empty");
    GMX_RELEASE_ASSERT(strlen(a3) > 0, "atom name a3 is empty");

    size_t index = std::lround(angle/spacing);
    auto   a     = std::find_if(bonds->angle.begin(), bonds->angle.end(),
                                [a1, a2, a3](const t_angle &a)
                                {
                                    int d = a.a2.compare(a2);
                                    if (0 == d)
                                    {
                                        return ((a.a1.compare(a1) == 0 && a.a3.compare(a3) == 0) ||
                                                (a.a1.compare(a3) == 0 && a.a3.compare(a1) == 0));
                                    }
                                    return false;
                                });

    if (a == bonds->angle.end())
    {
        t_angle aa;
        aa.a2 = a2;
        if (strcmp(a1, a3) < 0)
        {
            aa.a1 = a1;
            aa.a3 = a3;
        }
        else
        {
            aa.a1 = a3;
            aa.a3 = a1;
        }
        aa.histo.resize((int) (180/spacing) + 1, 0);
        aa.lsq = gmx_stats_init();
        bonds->angle.push_back(aa);
        a = bonds->angle.end()-1;
    }
    gmx_stats_add_point(a->lsq, 0, angle, 0, 0);
    a->histo[index]++;
    if (NULL != fplog)
    {
       fprintf(fplog, "%s angle-%s-%s-%s %g\n", molname,
	       a->a1.c_str(), a->a2.c_str(), a->a3.c_str(), angle);
    }
}

static void lo_add_dih(FILE *fplog, const char *molname, 
                       std::vector<t_dih> &dih,
                       const char *a1, const char *a2, const char *a3, const char *a4,
                       double angle, double spacing, unsigned int funcType)
                    
{
    if (angle < 0)
    {
        angle += 360;
    }
    if (funcType == F_IDIHS)
    {
        while (angle > 176)
        {
            angle -= 180;
        }
    }
    
    int index = std::lround(angle/spacing);
    if (index < 0)
    {
        index = 0;
    }
    auto d = std::find_if(dih.begin(), dih.end(),
                          [a1, a2, a3, a4](const t_dih &d)
                          {
                              return ((d.a1.compare(a1) == 0 && d.a2.compare(a2) == 0 &&
                                       d.a3.compare(a3) == 0 && d.a4.compare(a4) == 0) ||
                                      (d.a1.compare(a4) == 0 && d.a2.compare(a3) == 0 &&
                                       d.a3.compare(a2) == 0 && d.a4.compare(a1) == 0));
                          });
    
    if (dih.end() == d)
    {
        t_dih ddd;
        if (strcmp(a1, a4) < 0)
        {
            ddd.a1 = a1;
            ddd.a2 = a2;
            ddd.a3 = a3;
            ddd.a4 = a4;
        }
        else
        {
            ddd.a4 = a1;
            ddd.a3 = a2;
            ddd.a2 = a3;
            ddd.a1 = a4;
        }
        if (NULL != debug)
        {
            fprintf(debug, "NEWDIH  %5s  %5s  %5s  %5s\n",
                    a1, a2, a3, a4);
        }
        ddd.histo.resize((int) (360/spacing) + 1, 0);
        ddd.lsq = gmx_stats_init();
        dih.push_back(ddd);
        d = dih.end()-1;
    }
    gmx_stats_add_point(d->lsq, 0, angle, 0, 0);
    d->histo[index]++;
    if (NULL != fplog)
    {
        fprintf(fplog, "%s dihedral-%s-%s-%s-%s %g\n", molname,
                d->a1.c_str(), d->a2.c_str(), d->a3.c_str(), d->a4.c_str(), angle);
    }
}

static void add_dih(FILE *fplog, const char *molname, t_bonds *b,
                    const char *a1, const char *a2, const char *a3, const char *a4,
                    double angle, double spacing, unsigned int funcType)
{
    lo_add_dih(fplog, molname, 
               (F_PDIHS == funcType) ? b->dih : b->imp, 
               a1, a2, a3, a4, angle, spacing, funcType);
}

static void lo_dump_histo(char *fn, char *xaxis, const gmx_output_env_t *oenv, int Nsample,
                          int n, const int histo[], double spacing)
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
    int  N;
    char buf[256];

    for (const auto &i : b->bond)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0) && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "bond-%s-%s-%d.xvg", i.a1.c_str(), i.a2.c_str(), i.order);
            lo_dump_histo(buf, (char *)"Distance (pm)", oenv, N,
                          i.histo.size(), i.histo.data(), bspacing);
        }
    }
    for (const auto &i : b->angle)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0) && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "angle-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str());
            lo_dump_histo(buf, (char *)"Angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
    for (const auto &i : b->dih)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0)  && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "dih-%s-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str());
            lo_dump_histo(buf, (char *)"Dihedral angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
    for (const auto &i : b->imp)
    {
        if ((gmx_stats_get_npoints(i.lsq, &N) == 0)  && (i.histo.size() > 0))
        {
            snprintf(buf, sizeof(buf), "imp-%s-%s-%s-%s.xvg",
                     i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str());
            lo_dump_histo(buf, (char *)"Improper angle (deg.)", oenv, N,
                          i.histo.size(), i.histo.data(), aspacing);
        }
    }
}

static void round_numbers(real *av, real *sig)
{
    *av  = ((int)(*av*100))/100.0;
    *sig = ((int)(*sig*100+50))/100.0;
}

void update_pd(FILE *fp, t_bonds *b, Poldata &pd,
               real Dm, real beta, real kt, 
	       real klin, real kp, real kub)
{   
    int         N;
    real        av, sig;
    char        pbuf[256];
    std::vector<std::string> atoms;

    auto        morse             = pd.findForces(eitBONDS);
    auto        angle             = pd.findForces(eitANGLES);
    auto        linear_angle      = pd.findForces(eitLINEAR_ANGLES);
    auto        proper_dihedral   = pd.findForces(eitPROPER_DIHEDRALS);
    auto        improper_dihedral = pd.findForces(eitIMPROPER_DIHEDRALS);

    morse->eraseListedForce();
    angle->eraseListedForce();
    linear_angle->eraseListedForce();
    proper_dihedral->eraseListedForce();
    improper_dihedral->eraseListedForce();

    for (auto &i : b->bond)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        sprintf(pbuf, "%g  %g", Dm, beta);
        round_numbers(&av, &sig); // Rounding the numbers to 1/10 pm and 1/10 degree
	atoms = {i.a1, i.a2};
        morse->addForce(atoms, pbuf, av, sig, N);

        fprintf(fp, "bond-%s-%s len %g sigma %g (pm) N = %d%s\n",
                i.a1.c_str(),i.a2.c_str(), av, sig, N,
                (sig > 1.5) ? " WARNING" : "");
    }

    for (auto &i : b->angle)
    {

        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        round_numbers(&av, &sig);      
	atoms = {i.a1, i.a2, i.a3};
        if ((av > 175) || (av < 5))
        {
	    sprintf(pbuf, "%g",  klin);
	    linear_angle->addForce(atoms, pbuf, av, sig, N);
          
	    fprintf(fp, "linear_angle-%s-%s-%s angle %g sigma %g (deg) N = %d%s\n",
		    i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), av, sig, N,
		    (sig > 3) ? " WARNING" : "");
        }
        else
        {
	    sprintf(pbuf, "%g  %g", kt, kub);
            angle->addForce(atoms, pbuf, av, sig, N);

	    fprintf(fp, "harmonic_angle-%s-%s-%s angle %g sigma %g (deg) N = %d%s\n",
		    i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), av, sig, N,
		    (sig > 3) ? " WARNING" : "");
        }	
    }
    
    
    for (auto &i : b->dih)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        sprintf(pbuf, "%g  1", kp);   
        round_numbers(&av, &sig);
	atoms = {i.a1, i.a2, i.a3, i.a4};
        proper_dihedral->addForce(atoms, pbuf, av, sig, N);

        fprintf(fp, "dihedral-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str(), av, sig);
    }

    for (auto &i : b->imp)
    {
        gmx_stats_get_average(i.lsq, &av);
        gmx_stats_get_sigma(i.lsq, &sig);
        gmx_stats_get_npoints(i.lsq, &N);
        sprintf(pbuf, "%g  1", kp);  
        round_numbers(&av, &sig);
        atoms = {i.a1, i.a2, i.a3, i.a4};
        improper_dihedral->addForce(atoms, pbuf, av, sig, N);

        fprintf(fp, "improper-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                i.a1.c_str(), i.a2.c_str(), i.a3.c_str(), i.a4.c_str(), av, sig);
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
    static real                      Dm             = 0, kt = 0, kp = 0, beta = 0, klin = 0, kub = 0;
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
	{ "-kub",   FALSE, etREAL, {&kub},
	  "Urey_Bradley force constant" },
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

    FILE                            *fp;
    ChargeDistributionModel          iDistributionModel;
    time_t                           my_t;
    t_bonds                         *bonds = new(t_bonds);
    rvec                             dx, dx2, r_ij, r_kj, r_kl, mm, nn;
    t_pbc                            pbc;
    int                              t1, t2, t3;
    matrix                           box;
    real                             sign;
    double                           bspacing = 1;   /* pm */
    double                           aspacing = 0.5; /* degree */
    double                           dspacing = 1;   /* degree */
    gmx_output_env_t                *oenv     = NULL;
    Poldata                          pd;
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

    iDistributionModel = name2eemtype(cqdist[0]);

    fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# alexandria is part of G R O M A C S:\n#\n");
    fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());

    MolSelect gms;
    gms.read(opt2fn("-sel", NFILE, fnm));
    printf("There are %d molecules in the training group\n",
           gms.count(imsTrain));

    /* Read standard atom properties */
    aps = gmx_atomprop_init();

    /* Read force field stuff */
    try
    {
        readPoldata(opt2fn_null("-d", NFILE, fnm), pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    /* read Molprops */
    nfiles = opt2fns(&fns, "-f", NFILE, fnm);
    std::vector<alexandria::MolProp> mp;
    int nwarn = merge_xml(nfiles, fns, mp, NULL, NULL, NULL, aps, pd, TRUE);
    if (nwarn > maxwarn)
    {
        printf("Too many warnings (%d). Terminating.\n", nwarn);
        return 0;
    }

    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if (gms.status(mpi->getIupac()) == imsTrain)
        {
            alexandria::MyMol mmi;
            int               i;
            mmi.molProp()->Merge(mpi);
            if (mmi.molProp()->getMolname().size() == 0)
            {
                printf("Empty molname for molecule with formula %s\n", 
                       mmi.molProp()->formula().c_str());
                continue;
            }
            immStatus imm = mmi.GenerateTopology(aps, pd, lot, iDistributionModel,
                                                 false, false, bDih, false);

            if (immOK != imm)
            {
                if (NULL != debug)
                {
                    fprintf(debug, "Could not make topology for %s, reason %s\n", 
                            mmi.molProp()->getMolname().c_str(),
                            immsg(imm) );
                }
                continue;
            }
#define ATP(ii) (*mmi.topology_->atoms.atomtype[ii])
            for (i = 0; (i < mmi.topology_->atoms.nr); i++)
            {
                std::string btpi;
                if (!pd.atypeToBtype(*mmi.topology_->atoms.atomtype[i], btpi))
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
	        fprintf(debug, "You may need to check the number of atoms for %s\n",
                        mmi.molProp()->getMolname().c_str());
                continue;
            }

	    for (auto fs = pd.forcesBegin(); fs != pd.forcesEnd(); fs++)
	    {
	        if (eitBONDS == fs->iType())
		{
		    unsigned int funcType = fs->fType();

		    for (int j = 0; (j < mmi.ltop_->idef.il[funcType].nr); 
			 j += interaction_function[funcType].nratoms+1)
		    {
		        int         ai = mmi.ltop_->idef.il[funcType].iatoms[j+1];
			int         aj = mmi.ltop_->idef.il[funcType].iatoms[j+2];
			rvec_sub(mmi.x_[ai], mmi.x_[aj], dx);
			std::string cai, caj;
			if (pd.atypeToBtype(*mmi.topology_->atoms.atomtype[ai], cai) &&
			    pd.atypeToBtype(*mmi.topology_->atoms.atomtype[aj], caj))
			{
			    for (alexandria::BondIterator bi = mmi.molProp()->BeginBond(); 
				 (bi < mmi.molProp()->EndBond()); bi++)
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
				    add_bond(fp, mmi.molProp()->getMolname().c_str(), bonds, 
					     cai.c_str(), caj.c_str(), 1000*norm(dx),
					     bspacing, xb);
				    break;
				}
			    }
			}
			else
			{
			    fprintf(stderr, "No bond_atom type for either %s or %s\n",
				    ATP(ai), ATP(aj));
			}
		    }
		}
		else if (eitANGLES == fs->iType() || 
			 eitLINEAR_ANGLES == fs->iType())
		{
		    unsigned int funcType = fs->fType();

		    for (int j = 0; (j < mmi.ltop_->idef.il[funcType].nr); 
			 j+= interaction_function[funcType].nratoms+1)
		    {
		        int ai = mmi.ltop_->idef.il[funcType].iatoms[j+1];
			int aj = mmi.ltop_->idef.il[funcType].iatoms[j+2];
			int ak = mmi.ltop_->idef.il[funcType].iatoms[j+3];
			rvec_sub(mmi.x_[ai], mmi.x_[aj], dx);
			rvec_sub(mmi.x_[ak], mmi.x_[aj], dx2);
			double      ang = RAD2DEG*gmx_angle(dx, dx2);
			std::string cai, caj, cak;
			if (pd.atypeToBtype(*mmi.topology_->atoms.atomtype[ai], cai) &&
			    pd.atypeToBtype(*mmi.topology_->atoms.atomtype[aj], caj) &&
			    pd.atypeToBtype(*mmi.topology_->atoms.atomtype[ak], cak))
			{
			    add_angle(fp, mmi.molProp()->getMolname().c_str(), bonds, 
				      cai.c_str(), caj.c_str(), cak.c_str(), ang, aspacing);
			    if (NULL != debug)
			    {
			        fprintf(debug, "Molname: %s  btype1: %s  btype2: %s  btype3: %s  angle: %0.2f\n", 
					mmi.molProp()->getMolname().c_str(), cai.c_str(), caj.c_str(), cak.c_str(), ang);
			    }
			}
			else
			{
			    fprintf(stderr, "No bond_atom type for either %s, %s or %s in molecule %s\n",
				    ATP(ai), ATP(aj), ATP(ak), mmi.molProp()->getMolname().c_str());
			}
		    }
		}
		else if (eitPROPER_DIHEDRALS   == fs->iType()|| 
			 eitIMPROPER_DIHEDRALS == fs->iType())
		{
		    unsigned int funcType = fs->fType();
		
		    for (int j = 0; (j < mmi.ltop_->idef.il[funcType].nr); 
			 j += interaction_function[funcType].nratoms+1)
		    {
		        int    ai  = mmi.ltop_->idef.il[funcType].iatoms[j+1];
			int    aj  = mmi.ltop_->idef.il[funcType].iatoms[j+2];
			int    ak  = mmi.ltop_->idef.il[funcType].iatoms[j+3];
			int    al  = mmi.ltop_->idef.il[funcType].iatoms[j+4];
			double ang = RAD2DEG*dih_angle(mmi.x_[ai], mmi.x_[aj],
						       mmi.x_[ak], mmi.x_[al],
						       &pbc, r_ij, r_kj, r_kl, mm, nn, /* out */
						       &sign, &t1, &t2, &t3);
			std::string cai, caj, cak, cal;
			if (pd.atypeToBtype(*mmi.topology_->atoms.atomtype[ai], cai) &&
			    pd.atypeToBtype(*mmi.topology_->atoms.atomtype[aj], caj) &&
			    pd.atypeToBtype(*mmi.topology_->atoms.atomtype[ak], cak) &&
			    pd.atypeToBtype(*mmi.topology_->atoms.atomtype[al], cal))
			{
                            add_dih(fp, mmi.molProp()->getMolname().c_str(), bonds, 
				    cai.c_str(), caj.c_str(), cak.c_str(), cal.c_str(),
				    ang, dspacing, funcType);
			}
			else
			{
			    fprintf(stderr, "No bond_atom type for either %s, %s, %s or %s\n",
				    ATP(ai), ATP(aj), ATP(ak), ATP(al));
			}
		    }
		}
		else
	        {
		    fprintf(stderr, "Alexandria does not support the interaction type of %s\n",
			    iType2string(fs->iType()));
	        }
	    }

            clear_mat(box);
            set_pbc(&pbc, epbcNONE, box);
        }
    }
    sort_bonds(bonds);
    if (bHisto)
    {
        dump_histo(bonds, bspacing, aspacing, oenv);
    }
    update_pd(fp, bonds, pd, Dm, beta, kt, klin, kp, kub);

    writePoldata(opt2fn("-o", NFILE, fnm), pd, compress);

    printf("Extracted %d bondtypes, %d angletypes, %d dihedraltypes and %d impropertypes.\n",
           static_cast<int>(bonds->bond.size()), static_cast<int>(bonds->angle.size()),
           static_cast<int>(bonds->dih.size()), static_cast<int>(bonds->imp.size()));

    gmx_ffclose(fp);

    return 0;
}
