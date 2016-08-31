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
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

#include "categories.h"
#include "composition.h"
#include "molprop.h"
#include "molprop_tables.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "poldata.h"
#include "poldata_xml.h"

static void calc_frag_miller(alexandria::Poldata              &pd,
                             std::vector<alexandria::MolProp> &mp,
                             const alexandria::MolSelect      &gms)
{
    double                       bos0, polar, sig_pol;

    iMolSelect                   ims;
    char                        *null    = (char *)"0", *empirical = (char *)"empirical", *minimum = (char *)"minimum", *minus = (char *)"-", *nofile = (char *)"none";
    char                        *jobtype = (char *)"unknown";
    const char                  *program = "alexandria";
    const char                  *ang3;
    alexandria::CompositionSpecs cs;

    ang3    = unit2string(eg2cAngstrom3);
    if (0 == pd.getBosquePol( null, &bos0))
    {
        gmx_fatal(FARGS, "Can not find Bosque polarizability for %s", null);
    }

    for (auto &mpi : mp)
    {
        ims = gms.status(mpi.getIupac());
        if ((ims == imsTrain) || (ims == imsTest))
        {
            for (alexandria::CompositionSpecIterator csi = cs.beginCS(); (csi < cs.endCS()); ++csi)
            {
                alexandria::iComp ic = csi->iC();
                alexandria::MolecularCompositionIterator mci =
                    mpi.SearchMolecularComposition(csi->name());
                if (mci != mpi.EndMolecularComposition())
                {
                    double p         = 0, sp = 0;
                    double ahc       = 0, ahp = 0;
                    bool   bSupport  = true;
                    int    natom_tot = 0, Nelec = 0;
                    for (alexandria::AtomNumIterator ani = mci->BeginAtomNum(); bSupport && (ani < mci->EndAtomNum()); ani++)
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
                            if (NULL != debug)
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

static void write_corr_xvg(const char *fn,
                           std::vector<alexandria::MolProp> &mp,
                           MolPropObservable mpo, 
                           alexandria::t_qmcount *qmc,
                           real rtoler, real atoler,
                           const gmx_output_env_t *oenv, 
                           const alexandria::MolSelect &gms,
                           char *exp_type)
{
    FILE *fp;

    fp  = xvgropen(fn, "", "Exper.", "Calc. - Exper.", oenv);
    for (int i = 0; (i < qmc->n); i++)
    {
        fprintf(fp, "@s%d legend \"%s/%s-%s\"\n", i, qmc->method[i], qmc->basis[i],
                qmc->type[i]);
        fprintf(fp, "@s%d line linestyle 0\n", i);
        fprintf(fp, "@s%d symbol %d\n", i, i+1);
        fprintf(fp, "@s%d symbol size %g\n", i, 0.5);
        fprintf(fp, "@s%d symbol fill color %d\n", i, i+1);
        fprintf(fp, "@s%d symbol fill pattern 1\n", i);
    }
    for (int i = 0; (i < qmc->n); i++)
    {
        char LevelOfTheory[256];
        snprintf(LevelOfTheory, sizeof(LevelOfTheory), "%s/%s", qmc->method[i], qmc->basis[i]);
        if (debug)
        {
            fprintf(debug, "QM: %s LoT: %s\n", qmc->lot[i], LevelOfTheory);
        }
        int nout = 0;
        fprintf(fp, "@type xydy\n");
        for (auto &mpi : mp)
        {
            if (mpi.getMolname().compare("water") == 0)
            {
                printf("%s\n", mpi.getMolname().c_str());
            }
            double exp_val, exp_error;
            double Texp    = -1;
            bool   bExp    = mpi.getProp(mpo, iqmExp, NULL, NULL,
                                      exp_type,
                                      &exp_val, &exp_error, &Texp);
            iMolSelect ims = gms.status(mpi.getIupac());
            if ((ims == imsTrain) || (ims == imsTest))
            {
                double Tqm  = -1;
                double qm_val, qm_error;
                bool   bQM  = mpi.getProp(mpo, iqmQM, LevelOfTheory, NULL, //qmc->conf[k],
                                          qmc->type[i],
                                          &qm_val, &qm_error, &Tqm);
                if (bExp && bQM) // && (strcmp(exp_type, qmc->type[i]) == 0))
                {
                    fprintf(fp, "%8.3f  %8.3f  %8.3f\n", exp_val, qm_val-exp_val, qm_error);
                    double diff = fabs(qm_val-exp_val);
                    if (debug &&
                        (((atoler > 0) && (diff >= atoler)) ||
                         ((exp_val != 0) && (fabs(diff/exp_val) > rtoler))))
                    {
                        fprintf(debug, "OUTLIER: %s Exp: %g, Calc: %g +/- %g\n",
                                mpi.getIupac().c_str(), exp_val, qm_val, qm_error);
                        nout++;
                    }
                }
                else if (NULL != debug)
                {
                    fprintf(debug, "%s bQM = %d bExp = %d\n", mpi.getMolname().c_str(),
                            bQM ? 1 : 0, bExp ? 1 : 0);
                }
            }
        }
        fprintf(fp, "&\n");
        if (debug)
        {
            fprintf(debug, "There were %3d outliers for %s.\n",
                    nout, qmc->lot[i]);
        }
    }
    fclose(fp);
    do_view(oenv, fn, NULL);
}

typedef struct {
    int    nref;
    int   *rcount;
    char **refs;
} t_refcount;

static void add_refc(t_refcount *rc, const char *ref)
{
    int i;

    for (i = 0; (i < rc->nref); i++)
    {
        if (strcmp(ref, rc->refs[i]) == 0)
        {
            rc->rcount[i]++;
            break;
        }
    }
    if (i == rc->nref)
    {
        rc->nref++;
        srenew(rc->rcount, rc->nref);
        srenew(rc->refs, rc->nref);
        rc->rcount[i] = 1;
        rc->refs[i]   = strdup(ref);
    }
}

static void gmx_molprop_analyze(std::vector<alexandria::MolProp> &mp,
                                alexandria::Poldata &pd,
                                gmx_bool bCalcPol,
                                MolPropObservable mpo, char *exp_type,
                                char *lot,
                                real rtoler, real atoler, real outlier,
                                char *fc_str, gmx_bool bPrintAll,
                                gmx_bool bStatsTable,
                                int catmin,
                                const char *categoryfn,
                                gmx_bool bPropTable, gmx_bool bCompositionTable,
                                gmx_bool bPrintBasis, gmx_bool bPrintMultQ,
                                const char *texfn,
                                const char *xvgfn,
                                const gmx_output_env_t *oenv,
                                const alexandria::MolSelect &gms,
                                const char *selout)
{
    alexandria::CategoryList       cList;
    FILE                          *fp, *gp;
    int                            i, ntot;
    alexandria::t_qmcount         *qmc;
    t_refcount                    *rc;
    double                         T, value, error, vec[3];
    tensor                         quadrupole;
    const char                    *iupac;

    if (bCalcPol)
    {
        calc_frag_miller(pd, mp, gms);
    }

    qmc = find_calculations(mp, mpo, fc_str);

    snew(rc, 1);
    for (auto &mpi : mp)
    {
        for (alexandria::ExperimentIterator ei = mpi.BeginExperiment(); (ei < mpi.EndExperiment()); ++ei)
        {
            T = -1;
            if (ei->getVal(exp_type, mpo, &value, &error, &T, vec, quadrupole))
            {
                add_refc(rc, ei->getReference().c_str());
            }
        }
    }

    printf("--------------------------------------------------\n");
    printf("      Some statistics for %s\n", mpo_name[mpo]);
    ntot = 0;
    for (i = 0; (i < rc->nref); i++)
    {
        printf("There are %d experiments with %s as reference\n",
               rc->rcount[i], rc->refs[i]);
        ntot += rc->rcount[i];
    }
    printf("There are %d entries with experimental %s of type %s\n", ntot,
           mpo_name[mpo], exp_type);
    if (0 == ntot)
    {
        printf("   did you forget to pass the -exp_type flag?\n");
    }
    for (i = 0; (i < qmc->n); i++)
    {
        printf("There are %d calculation results using %s/%s type %s\n",
               qmc->count[i], qmc->method[i], qmc->basis[i], qmc->type[i]);
    }
    printf("--------------------------------------------------\n");

    makeCategoryList(cList, mp, gms, imsTrain);
    fp = gmx_ffopen(texfn, "w");
    if (bStatsTable)
    {
        gmx_molprop_stats_table(fp, mpo, mp, qmc, exp_type,
                                outlier, cList, gms, imsTrain);
        if (0)
        {
            alexandria::CategoryList cListTest;
            makeCategoryList(cListTest, mp, gms, imsTest);
            gmx_molprop_stats_table(fp, mpo, mp, qmc, exp_type,
                                    outlier, cListTest, gms, imsTest);
        }
    }
    if (bPropTable)
    {
        gmx_molprop_prop_table(fp, mpo, rtoler, atoler, mp, qmc, exp_type, bPrintAll, bPrintBasis,
                               bPrintMultQ, gms, imsTrain);
        gmx_molprop_prop_table(fp, mpo, rtoler, atoler, mp, qmc, exp_type, bPrintAll, bPrintBasis,
                               bPrintMultQ, gms, imsTest);
        if (NULL != selout)
        {
            gp = fopen(selout, "w");
            for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
            {
                iupac = mpi->getIupac().c_str();
                if ((NULL != iupac) && (strlen(iupac) > 0))
                {
                    std::string myref, mylot;
                    if (mpi->getPropRef(mpo, iqmBoth, lot, NULL, NULL,
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
        gmx_molprop_composition_table(fp, mp, gms, imsTrain);
        gmx_molprop_composition_table(fp, mp, gms, imsTest);
    }
    fclose(fp);
    if (NULL != categoryfn)
    {
        fp = fopen(categoryfn, "w");
        gmx_molprop_category_table(fp, catmin, cList);
        fclose(fp);
    }
    if (NULL != xvgfn)
    {
        write_corr_xvg(xvgfn, mp, mpo, qmc, rtoler, atoler, oenv, gms, exp_type);
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
        { efXVG, "-c",      "correl",    ffWRITE  }
    };
#define NFILE (sizeof(fnm)/sizeof(fnm[0]))
    static char                     *sort[]      = { NULL, (char *)"molname", (char *)"formula", (char *)"composition", (char *)"selection", NULL };
    static char                     *prop[]      = { NULL, (char *)"potential", (char *)"dipole", (char *)"quadrupole", (char *)"polarizability", (char *)"energy", (char *)"entropy", NULL };
    static char                     *fc_str      = (char *)"";
    static char                     *exp_type    = (char *)"";
    static char                     *lot         = (char *)"B3LYP/aug-cc-pVTZ";
    static real                      rtoler      = 0.15, atoler = 0, outlier = 1;
    static real                      th_toler    = 170, ph_toler = 5;
    static gmx_bool                  bMerge      = TRUE, bAll = FALSE, bCalcPol = TRUE, bPrintBasis = TRUE, bPrintMultQ = FALSE;
    static gmx_bool                  bStatsTable = TRUE, bCompositionTable = FALSE, bPropTable = TRUE;
    static int                       catmin      = 1;
    static int                       maxwarn     = 0;
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
    int                              npa;
    int                              i;

    std::vector<alexandria::MolProp> mp;
    MolPropSortAlgorithm             mpsa;
    MolPropObservable                mpo;
    gmx_atomprop_t                   ap;
    gmx_output_env_t                *oenv;
    char                           **mpname = NULL;
    int                              nmpfile;

    npa = sizeof(pa)/sizeof(pa[0]);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           npa, pa,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }

    ap = gmx_atomprop_init();

    nmpfile = opt2fns(&mpname, "-m", NFILE, fnm);
    alexandria::MolSelect gms;
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
    
    alexandria::Poldata pd;
    try 
    {
        alexandria::readPoldata(opt2fn("-d", NFILE, fnm), pd, ap);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if (bMerge)
    {
        int nwarn = merge_xml(nmpfile, mpname, mp, NULL, NULL, NULL, ap, pd, TRUE);
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
    gmx_molprop_analyze(mp, pd, bCalcPol,
                        mpo, exp_type, lot, rtoler, atoler, outlier, fc_str, bAll,
                        bStatsTable,
                        catmin,
                        opt2fn_null("-cat", NFILE, fnm),
                        bPropTable, bCompositionTable,
                        bPrintBasis, bPrintMultQ,
                        opt2fn("-t", NFILE, fnm), opt2fn("-c", NFILE, fnm),
                        oenv, gms,
                        opt2fn_null("-selout", NFILE, fnm));

    return 0;
}
