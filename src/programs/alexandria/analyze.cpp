/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"

#include "categories.h"
#include "composition.h"
#include "gmxpre.h"
#include "molprop.h"
#include "molprop_tables.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "poldata.h"
#include "poldata_xml.h"

class RefCount
{
    private:
        int         count_;
        std::string ref_;
    public:
        RefCount(const std::string ref) : count_(1), ref_(ref) {};

        int count() const { return count_; }

        void increment() { count_++; }

        const std::string ref() const { return ref_; }
};

static void add_refc(std::vector<RefCount> &rc, std::string ref)
{
    for (auto &r : rc)
    {
        if (r.ref().compare(ref) == 0)
        {
            r.increment();
            return;
        }
    }
    rc.push_back(RefCount(ref));
}

static void calc_frag_miller(alexandria::Poldata              &pd,
                             std::vector<alexandria::MolProp> &mp,
                             const alexandria::MolSelect      &gms)
{
    double                       bos0      = 0;
    double                       polar     = 0; 
    double                       sig_pol   = 0;
    char                        *null      = (char *)"0";
    char                        *empirical = (char *)"empirical";
    char                        *minimum   = (char *)"minimum";
    char                        *minus     = (char *)"-";
    char                        *nofile    = (char *)"none";
    const char                  *program   = "alexandria";
    const char                  *ang3;
    iMolSelect                   ims;
    alexandria::jobType          jobtype   = alexandria::JOB_UNKNOWN;
    alexandria::CompositionSpecs cs;

    ang3    = unit2string(eg2cAngstrom3);
    if (0 == pd.getBosquePol( null, &bos0))
    {
        gmx_fatal(FARGS, "Cannot find Bosque polarizability for %s", null);
    }

    for (auto &mpi : mp)
    {
        ims = gms.status(mpi.getIupac());
        if ((ims == imsTrain) || (ims == imsTest))
        {
            for (auto csi = cs.beginCS(); csi < cs.endCS(); ++csi)
            {
                alexandria::iComp ic = csi->iC();
                auto mci = mpi.SearchMolecularComposition(csi->name());
                if (mci != mpi.EndMolecularComposition())
                {
                    double p         = 0, sp = 0;
                    double ahc       = 0, ahp = 0;
                    bool   bSupport  = true;
                    int    natom_tot = 0, Nelec = 0;
                    for (auto ani = mci->BeginAtomNum(); bSupport && (ani < mci->EndAtomNum()); ani++)
                    {
                        const char *atomname = ani->getAtom().c_str();
                        int         natom    = ani->getNumber();
                        switch (ic)
                        {
                            case alexandria::iCalexandria:
                                bSupport = (pd.getAtypePol( (char *)atomname, &polar, &sig_pol) == 1);
                                break;
                            case alexandria::iCbosque:
                                bSupport = (pd.getBosquePol( (char *)atomname, &polar) == 1);
                                sig_pol  = 0;
                                break;
                            case alexandria::iCmiller:
                            {
                                double      tau_ahc, alpha_ahp;
                                int         atomnumber;
                                std::string aequiv;
                                bSupport = (pd.getMillerPol((char *)atomname,
                                                            &atomnumber,
                                                            &tau_ahc,
                                                            &alpha_ahp,
                                                            aequiv) == 1);

                                ahc   += tau_ahc*natom;
                                ahp   += alpha_ahp*natom;
                                Nelec += atomnumber*natom;
                            }
                            break;
                            default:
                                gmx_incons("Out of range");
                        }

                        if (bSupport && (ic != alexandria::iCmiller))
                        {
                            p         += polar*natom;
                            sp        += gmx::square(sig_pol)*natom;
                            natom_tot += natom;
                        }
                    }
                    sp = sqrt(sp/natom_tot);
                    if (bSupport)
                    {
                        const char *type = csi->abbreviation();
                        const char *ref  = csi->reference();
                        if (ic == alexandria::iCmiller)
                        {
                            alexandria::Experiment calc1(program, type, (char *)"ahc",
                                                         ref, minimum, nofile, jobtype);
                            ahc = 4*gmx::square(ahc)/Nelec;
                            alexandria::MolecularPolarizability md1(empirical, ang3, 0, 0, 0, 0, 0, 0, 0, ahc, 0);
                            calc1.AddPolar(md1);
                            mpi.AddExperiment(calc1);

                            alexandria::Experiment              calc2(program, type, (char *)"ahp",
                                                                      ref, minimum, nofile, jobtype);
                            alexandria::MolecularPolarizability md2(empirical, ang3, 0, 0, 0, 0, 0, 0, 0, ahp, 0);
                            calc2.AddPolar(md2);
                            mpi.AddExperiment(calc2);
                        }
                        else
                        {
                            alexandria::Experiment              calc(program, type, minus,
                                                                     ref, minimum, nofile, jobtype);
                            alexandria::MolecularPolarizability md(empirical, ang3, 0, 0, 0, 0, 0, 0, 0, p, sp);
                            calc.AddPolar(md);
                            mpi.AddExperiment(calc);
                            if (nullptr != debug)
                            {
                                fprintf(debug, "Added polarizability %g for %s\n", p, mpi.getIupac().c_str());
                            }
                        }
                    }
                }
            }
        }
    }
}

static void write_corr_xvg(FILE                             *fplog,
                           const char                       *fn,
                           std::vector<alexandria::MolProp> &mp,
                           MolPropObservable                 mpo,
                           const alexandria::QmCount        &qmc,
                           real                              rtoler,
                           real                              atoler,
                           const gmx_output_env_t           *oenv,
                           const alexandria::MolSelect      &gms,
                           char                             *exp_type)
{
    int          i         = 0;
    int          n         = 0;
    real         a         = 0;
    real         da        = 0;
    real         b         = 0;
    real         db        = 0;
    real         mse       = 0;
    real         mae       = 0;
    real         chi2      = 0;
    real         rmsd      = 0;
    real         Rfit      = 0;   
    double       exp_val   = 0;
    double       exp_error = 0;
    double       qm_val    = 0;
    double       qm_error  = 0;
    double       diff      = 0;
    double       Texp      = -1;
    double       Tqm       = -1;
    bool         bExp      = false;
    bool         bQM       = false;
    
    FILE        *fp;
    gmx_stats_t  lsq[qmc.nCalc()];    
    
    fp  = xvgropen(fn, "", "Exper.", "Calc. - Exper.", oenv);   
    for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q, ++i)
    {
        lsq[i] = gmx_stats_init();
        fprintf(fp, "@s%d legend \"%s/%s-%s\"\n", i,
                q->method().c_str(),
                q->basis().c_str(),
                q->type().c_str());
        fprintf(fp, "@s%d line linestyle 0\n", i);
        fprintf(fp, "@s%d symbol %d\n", i, i+1);
        fprintf(fp, "@s%d symbol size %g\n", i, 0.5);
        fprintf(fp, "@s%d symbol fill color %d\n", i, i+1);
        fprintf(fp, "@s%d symbol fill pattern 1\n", i);
    }
    i = 0;
    for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q, ++i)
    {
        char LevelOfTheory[256];
        snprintf(LevelOfTheory, sizeof(LevelOfTheory), "%s/%s", q->method().c_str(), q->basis().c_str());
        if (debug)
        {
            fprintf(debug, "QM: %s LoT: %s\n", q->lot().c_str(), LevelOfTheory);
        }
        int nout = 0;
        fprintf(fp, "@type xydy\n");
        for (auto &mpi : mp)
        {
            exp_val   = 0;
            exp_error = 0;
            bExp      = mpi.getProp(mpo, 
                                    iqmExp, 
                                    "", 
                                    "",
                                    exp_type,
                                    &exp_val, 
                                    &exp_error, 
                                    &Texp);
            auto ims = gms.status(mpi.getIupac());
            if ((ims == imsTrain) || (ims == imsTest))
            {
                
                qm_val   = 0;
                qm_error = 0;
                bQM      = mpi.getProp(mpo,                
                                       iqmQM, 
                                       LevelOfTheory, 
                                       "",
                                       q->type().c_str(),
                                       &qm_val,
                                       &qm_error, 
                                       &Tqm);
                if (bExp && bQM)
                {
                    gmx_stats_add_point(lsq[i], exp_val, qm_val, 0, 0);
                    fprintf(fp, "%8.3f  %8.3f  %8.3f\n", exp_val, qm_val-exp_val, qm_error);
                    diff = fabs(qm_val-exp_val);                  
                    if (((atoler > 0) && (diff >= atoler)) || (atoler == 0 && exp_val != 0 && (fabs(diff/exp_val) > rtoler)))
                    {
                        fprintf(fplog, "OUTLIER: %s Exp: %g, Calc: %g +/- %g Method:%s\n",
                                mpi.getIupac().c_str(), exp_val, qm_val, qm_error, q->method().c_str());
                        nout++;
                    }
                }
                else if (nullptr != debug)
                {
                    fprintf(debug, "%s bQM = %d bExp = %d\n", mpi.getMolname().c_str(), bQM ? 1 : 0, bExp ? 1 : 0);
                }
            }
        }
        fprintf(fp, "&\n");
        if (nout)
        {
            printf("There are %3d outliers for %s. Check the analyze log file for details.\n", nout, q->lot().c_str());
        }
    }
    
    fprintf(fplog, "Fitting %s data to y = ax + b\n", exp_type);
    fprintf(fplog, "%-12s %5s %13s %13s %8s %8s %8s %8s\n", "Method", "N", "a", "b", "R(%)", "RMSD", "MSE", "MAE");
    fprintf(fplog, "-----------------------------------------------------------------------------\n");              
    i = 0;
    for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q, ++i)
    {
        gmx_stats_get_npoints(lsq[i], &n);
        gmx_stats_get_ab(lsq[i], elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
        gmx_stats_get_rmsd(lsq[i],    &rmsd);
        gmx_stats_get_mse_mae(lsq[i], &mse, &mae);        
        fprintf(fplog, "%-12s %5d %6.3f(%5.4f) %6.3f(%5.4f) %7.2f %8.4f %8.4f %8.4f\n", 
                q->method().c_str(), n, a, da, b, db, Rfit*100, rmsd, mse, mae);
        gmx_stats_free(lsq[i]);
    }    
    fclose(fp);
    do_view(oenv, fn, nullptr);
}

static void alexandria_molprop_analyze(FILE                              *fplog,
                                std::vector<alexandria::MolProp>  &mp,
                                alexandria::Poldata               &pd,
                                gmx_bool                           bCalcPol,
                                MolPropObservable                  mpo, 
                                char                              *exp_type,
                                char                              *lot,
                                real                               rtoler,
                                real                               atoler, 
                                real                               outlier,
                                char                              *fc_str, 
                                gmx_bool                           bPrintAll,
                                gmx_bool                           bStatsTable,
                                int                                catmin,
                                const char                        *categoryfn,
                                gmx_bool                           bPropTable, 
                                gmx_bool                           bCompositionTable,
                                gmx_bool                           bPrintBasis, 
                                gmx_bool                           bPrintMultQ,
                                const char                        *texfn,
                                const char                        *xvgfn,
                                const gmx_output_env_t            *oenv,
                                const alexandria::MolSelect       &gms,
                                const char                        *selout)
{  
    int                       ntot;
    FILE                     *fp, *gp;
    double                    T, value, error, vec[3];
    tensor                    quadrupole;
    const char               *iupac;
    alexandria::QmCount       qmc;
    alexandria::CategoryList  cList;
    std::vector<RefCount>     rc;
    
    if (bCalcPol)
    {
        calc_frag_miller(pd, mp, gms);
    }
    find_calculations(mp, mpo, fc_str, &qmc);
    for (auto &mpi : mp)
    {
        for (auto ei = mpi.BeginExperiment(); ei < mpi.EndExperiment(); ++ei)
        {
            T = -1;
            if (ei->getVal(exp_type, mpo, &value, &error, &T, vec, quadrupole))
            {
                add_refc(rc, ei->getReference().c_str());
            }
        }
    }
    printf("--------------------------------------------------\n");
    printf("      Statistics for %s\n", mpo_name[mpo]);
    ntot = 0;
    for (const auto &r : rc)
    {
        printf("There are %d experiments with %s as reference\n",
               r.count(), r.ref().c_str());
        ntot += r.count();
    }
    printf("There are %d entries with experimental %s of type %s\n", ntot,
           mpo_name[mpo], exp_type);
    if (0 == ntot)
    {
        printf("   did you forget to pass the -exp_type flag?\n");
    }
    for (auto q = qmc.beginCalc(); q < qmc.endCalc(); ++q)
    {
        printf("There are %d calculation results using %s/%s type %s\n",
               q->count(), q->method().c_str(),
               q->basis().c_str(), q->type().c_str());
    }
    printf("--------------------------------------------------\n");
    makeCategoryList(cList, mp, gms, imsTrain);
    fp = gmx_ffopen(texfn, "w");
    if (bStatsTable)
    {
        alexandria_molprop_stats_table(fp, mpo, mp, qmc, exp_type,
                                outlier, cList, gms, imsTrain);
        if (0)
        {
            alexandria::CategoryList cListTest;
            makeCategoryList(cListTest, mp, gms, imsTest);
            alexandria_molprop_stats_table(fp, mpo, mp, qmc, exp_type,
                                           outlier, cListTest, gms, imsTest);
        }
    }
    if (bPropTable)
    {
        alexandria_molprop_prop_table(fp, mpo, rtoler, atoler, mp, qmc, exp_type, bPrintAll, bPrintBasis,
                                      bPrintMultQ, gms, imsTrain);
        alexandria_molprop_prop_table(fp, mpo, rtoler, atoler, mp, qmc, exp_type, bPrintAll, bPrintBasis,
                                      bPrintMultQ, gms, imsTest);
        if (nullptr != selout)
        {
            gp = fopen(selout, "w");
            for (auto mpi = mp.begin(); mpi < mp.end(); mpi++)
            {
                iupac = mpi->getIupac().c_str();
                if ((nullptr != iupac) && (strlen(iupac) > 0))
                {
                    std::string myref, mylot;
                    if (mpi->getPropRef(mpo, iqmBoth, lot, "", "",
                                        &value, &error, &T, myref, mylot,
                                        vec, quadrupole))
                    {
                        fprintf(gp, "%s|Train\n", iupac);
                    }
                }
            }
            fclose(gp);
        }
    }
    if (bCompositionTable)
    {
        alexandria_molprop_composition_table(fp, mp, gms, imsTrain);
        alexandria_molprop_composition_table(fp, mp, gms, imsTest);
    }
    fclose(fp);
    if (nullptr != categoryfn)
    {
        fp = fopen(categoryfn, "w");
        alexandria_molprop_category_table(fp, catmin, cList);
        fclose(fp);
    }
    if (nullptr != xvgfn)
    {
        write_corr_xvg(fplog, xvgfn, mp, mpo, qmc, rtoler, atoler, oenv, gms, exp_type);
    }
}


int alex_analyze(int argc, char *argv[])
{
    static const char               *desc[] = {
        "analyze reads a molecule database",
        "and produces tables and figures to describe the data.",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-m[tt] option). Missing molecules will be ignored. You can also write a",
        "selection file ([TT]-selout[tt]) that contains all the molecules in the",
        "output corresponding to the [TT]-prop[tt] flag. You can steer the content of",
        "this file e.g. by running the program with an empty [TT]-lot[tt] flag,",
        "yielding only those molecules for which experimental data is available."
    };
    t_filenm                         fnm[] = {
        { efDAT, "-d",      "gentop",    ffREAD   },
        { efDAT, "-m",      "allmols",   ffRDMULT },
        { efTEX, "-t",      "table",     ffWRITE  },
        { efTEX, "-cat",    "category",  ffOPTWR  },
        { efDAT, "-sel",    "molselect", ffREAD   },
        { efDAT, "-selout", "selout",    ffOPTWR  },
        { efXVG, "-c",      "correl",    ffWRITE  },
        { efLOG, "-g",      "analyze",   ffWRITE  }
    };
    int                 NFILE             = asize(fnm);
      
    static int          catmin            = 1;
    static int          maxwarn           = 0;
    static char        *fc_str            = (char *)"";
    static char        *exp_type          = (char *)"";
    static char        *lot               = (char *)"B3LYP/aug-cc-pVTZ";
    static real         rtoler            = 0.15;
    static real         atoler            = 0;
    static real         outlier           = 1;
    static real         th_toler          = 170;
    static real         ph_toler          = 5;
    static gmx_bool     bMerge            = true;
    static gmx_bool     bAll              = false;
    static gmx_bool     bCalcPol          = true;
    static gmx_bool     bPrintBasis       = true;
    static gmx_bool     bPrintMultQ       = false;
    static gmx_bool     bStatsTable       = true;
    static gmx_bool     bCompositionTable = false;
    static gmx_bool     bPropTable        = true;
    
    static char        *sort[]            = {nullptr, (char *)"molname", (char *)"formula", (char *)"composition", (char *)"selection", nullptr};
    static char        *prop[]            = {nullptr, (char *)"potential", (char *)"dipole", (char *)"quadrupole", (char *)"polarizability", (char *)"energy", (char *)"entropy", nullptr};

    t_pargs                          pa[]        = {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-rtol",    FALSE, etREAL, {&rtoler},
          "Relative tolerance for printing in bold in tables (see next option)" },
        { "-atol",    FALSE, etREAL, {&atoler},
          "Absolute tolerance for printing in bold in tables. If non-zero, an absolute tolerance in appropriate units, depending on property, will be used rather than a relative tolerance." },
        { "-outlier", FALSE, etREAL, {&outlier},
          "Outlier indicates a level (in units of sigma, one standard deviation). Calculations that deviate more than this level from the experiment are not taken into account when computing statistics. Moreover, the outliers are printed to the standard error. If outlier is 0, no action is taken. " },
        { "-merge",  FALSE, etBOOL, {&bMerge},
          "Merge molecule records in the input file and generate atomtype compositions based on the calculated geometry." },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Indicate the method and level of theory that were used together with experimental data in refining polarizabilities. If empty, is is assumed that only experimental data were used." },
        { "-prop",   FALSE, etENUM, {prop},
          "Property to print" },
        { "-catmin", FALSE, etINT, {&catmin},
          "Mininum number of molecules in a category for it to be printed" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-all",    FALSE, etBOOL, {&bAll},
          "Print calculated results for properties even if no experimental data is available to compare to" },
        { "-printbasis", FALSE, etBOOL, {&bPrintBasis},
          "Print the basis set in the property table" },
        { "-printmultq", FALSE, etBOOL, {&bPrintMultQ},
          "Print the multiplicity and charge in the property table" },
        { "-calcpol", FALSE, etBOOL, {&bCalcPol},
          "Calculate polarizabilities based on empirical methods" },
        { "-composition", FALSE, etBOOL, {&bCompositionTable},
          "Print a table of composition of the molecules" },
        { "-proptable", FALSE, etBOOL, {&bPropTable},
          "Print a table of properties (slect which with the [TT]-prop[tt] flag)." },
        { "-statstable", FALSE, etBOOL, {&bStatsTable},
          "Print a table of statistics per category" },
        { "-th_toler", FALSE, etREAL, {&th_toler},
          "If bond angles are larger than this value the group will be treated as a linear one and a virtual site will be created to keep the group linear" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "If dihedral angles are less than this (in absolute value) the atoms will be treated as a planar group with an improper dihedral being added to keep the group planar" },
        { "-fc_str", FALSE, etSTR, {&fc_str},
          "Selection of the stuff you want in the tables, given as a single string with spaces like: method1/basis1/type1:method2/basis2/type2 (you may have to put quotes around the whole thing in order to prevent the shell from interpreting it)." },
        { "-exp_type", FALSE, etSTR, {&exp_type},
          "The experimental property type that all the stuff above has to be compared to." }
    };
    
    FILE                            *fplog;   
    int                              npa;
    int                              i;
    alexandria::Poldata              pd;
    alexandria::MolSelect            gms;
    std::vector<alexandria::MolProp> mp;
    MolPropSortAlgorithm             mpsa;
    MolPropObservable                mpo;
    gmx_atomprop_t                   ap;
    gmx_output_env_t                *oenv;
    char                           **mpname = nullptr;
    int                              nmpfile;

    npa = asize(pa);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           npa, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    ap = gmx_atomprop_init();
    nmpfile = opt2fns(&mpname, "-m", NFILE, fnm);    
    gms.read(opt2fn("-sel", NFILE, fnm));
    mpsa = MPSA_NR;
    if (opt2parg_bSet("-sort", npa, pa))
    {
        for (i = 0; (i < MPSA_NR); i++)
        {
            if (strcasecmp(sort[0], sort[i+1]) == 0)
            {
                mpsa = (MolPropSortAlgorithm) i;
                break;
            }
        }
    }
    mpo = MPO_NR;
    if (opt2parg_bSet("-prop", npa, pa))
    {
        for (i = 0; (i < MPO_NR); i++)
        {
            if (strcasecmp(prop[0], prop[i+1]) == 0)
            {
                mpo = (MolPropObservable) i;
                break;
            }
        }
    }
    if (mpo == MPO_NR)
    {
        mpo = MPO_DIPOLE;
    }
    try
    {
        alexandria::readPoldata(opt2fn("-d", NFILE, fnm), pd, ap);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if (bMerge)
    {
        int nwarn = merge_xml(nmpfile, mpname, mp, nullptr, nullptr, nullptr, ap, pd, TRUE);
        if (nwarn > maxwarn)
        {
            printf("Too many warnings (%d). Terminating.\n", nwarn);
            return 0;
        }
    }
    else
    {
        MolPropRead((const char *)mpname[0], mp);
        generate_composition(mp, pd);
        generate_formula(mp, ap);
    }
    if (mpsa != MPSA_NR)
    {
        MolPropSort(mp, mpsa, ap, gms);
    }   
    fplog  = opt2FILE("-g", NFILE, fnm, "w");
    alexandria_molprop_analyze(fplog,
                               mp, 
                               pd, 
                               bCalcPol,
                               mpo, 
                               exp_type, 
                               lot, 
                               rtoler, 
                               atoler, 
                               outlier, 
                               fc_str, 
                               bAll,
                               bStatsTable,
                               catmin,
                               opt2fn_null("-cat", NFILE, fnm),
                               bPropTable, 
                               bCompositionTable,
                               bPrintBasis, 
                               bPrintMultQ,
                               opt2fn("-t", NFILE, fnm), 
                               opt2fn("-c", NFILE, fnm),
                               oenv, 
                               gms,
                               opt2fn_null("-selout", NFILE, fnm));
    gmx_ffclose(fplog);

    return 0;
}
