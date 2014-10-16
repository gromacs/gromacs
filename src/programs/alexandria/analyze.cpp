/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/math/vec.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/legacyheaders/viewit.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molselect.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molprop_tables.h"
#include "composition.h"

static void calc_frag_miller(gmx_poldata_t                     pd,
                             std::vector<alexandria::MolProp> &mp,
                             gmx_molselect_t                   gms)
{
    double                       bos0, polar, sig_pol;

    iMolSelect                   ims;
    char                        *null = (char *)"0", *empirical = (char *)"empirical", *minimum = (char *)"minimum", *minus = (char *)"-", *nofile = (char *)"none";
    const char                  *program;
    const char                  *ang3;
    alexandria::CompositionSpecs cs;

    ang3    = unit2string(eg2cAngstrom3);
    program = ShortProgram();
    if (0 == gmx_poldata_get_bosque_pol(pd, null, &bos0))
    {
        gmx_fatal(FARGS, "Can not find Bosque polarizability for %s", null);
    }

    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        const char *iupac = mpi->GetIupac().c_str();
        ims = gmx_molselect_status(gms, iupac);
        if ((ims == imsTrain) || (ims == imsTest))
        {
            for (alexandria::CompositionSpecIterator csi = cs.beginCS(); (csi < cs.endCS()); ++csi)
            {
                alexandria::iComp ic = csi->iC();
                alexandria::MolecularCompositionIterator mci =
                    mpi->SearchMolecularComposition(csi->name());
                if (mci != mpi->EndMolecularComposition())
                {
                    double p         = 0, sp = 0;
                    double ahc       = 0, ahp = 0;
                    bool   bSupport  = true;
                    int    natom_tot = 0, Nelec = 0;
                    for (alexandria::AtomNumIterator ani = mci->BeginAtomNum(); bSupport && (ani < mci->EndAtomNum()); ani++)
                    {
                        const char *atomname = ani->GetAtom().c_str();
                        int         natom    = ani->GetNumber();
                        switch (ic)
                        {
                            case alexandria::iCalexandria:
                                bSupport = (gmx_poldata_get_atype_pol(pd, (char *)atomname, &polar, &sig_pol) == 1);
                                break;
                            case alexandria::iCbosque:
                                bSupport = (gmx_poldata_get_bosque_pol(pd, (char *)atomname, &polar) == 1);
                                sig_pol  = 0;
                                break;
                            case alexandria::iCmiller:
                                double tau_ahc, alpha_ahp;
                                int    atomnumber;
                                bSupport = (gmx_poldata_get_miller_pol(pd, (char *)atomname, &atomnumber, &tau_ahc, &alpha_ahp) == 1);

                                ahc   += tau_ahc*natom;
                                ahp   += alpha_ahp*natom;
                                Nelec += atomnumber*natom;
                                break;
                            default:
                                gmx_incons("Out of range");
                        }

                        if (bSupport && (ic != alexandria::iCmiller))
                        {
                            p         += polar*natom;
                            sp        += sqr(sig_pol)*natom;
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
                            alexandria::Calculation calc1(program, type, (char *)"ahc",
                                                          ref, minimum, nofile);
                            ahc = 4*sqr(ahc)/Nelec;
                            alexandria::MolecularDipPolar md1(empirical, ang3, 0, 0, 0, ahc, 0);
                            calc1.AddPolar(md1);
                            mpi->AddCalculation(calc1);

                            alexandria::Calculation       calc2(program, type, (char *)"ahp",
                                                                ref, minimum, nofile);
                            alexandria::MolecularDipPolar md2(empirical, ang3, 0, 0, 0, ahp, 0);
                            calc2.AddPolar(md2);
                            mpi->AddCalculation(calc2);
                        }
                        else
                        {
                            alexandria::Calculation       calc(program, type, minus,
                                                               ref, minimum, nofile);
                            alexandria::MolecularDipPolar md(empirical, ang3, 0, 0, 0, p, sp);
                            calc.AddPolar(md);
                            mpi->AddCalculation(calc);
                            if (NULL != debug)
                            {
                                fprintf(debug, "Added polarizability %g for %s\n", p, iupac);
                            }
                        }
                    }
                }
            }
        }
    }
}

static void write_corr_xvg(const char *fn,
                           std::vector<alexandria::MolProp> mp,
                           MolPropObservable mpo, t_qmcount *qmc,
                           real rtoler, real atoler,
                           output_env_t oenv, gmx_molselect_t gms,
                           char *exp_type)
{
    alexandria::MolPropIterator mpi;
    FILE                       *fp;
    int    i, k, nout;
    iMolSelect                  ims;
    char   lbuf[256];
    double exp_val, qm_val, diff, qm_error;

    fp  = xvgropen(fn, "", "Exper.", "Calc. - Exper.", oenv);
    for (i = 0; (i < qmc->n); i++)
    {
        fprintf(fp, "@s%d legend \"%s/%s-%s\"\n", i, qmc->method[i], qmc->basis[i],
                qmc->type[i]);
        fprintf(fp, "@s%d line linestyle 0\n", i);
        fprintf(fp, "@s%d symbol %d\n", i, i+1);
        fprintf(fp, "@s%d symbol size %g\n", i, 0.5);
        fprintf(fp, "@s%d symbol fill color %d\n", i, i+1);
        fprintf(fp, "@s%d symbol fill pattern 1\n", i);
    }
    for (i = 0; (i < qmc->n); i++)
    {
        fprintf(fp, "@type xydy\n");
        sprintf(lbuf, "%s/%s", qmc->method[i], qmc->basis[i]);
        if (debug)
        {
            fprintf(debug, "QM: %s\n", qmc->lot[i]);
        }
        nout = 0;
        for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            ims = gmx_molselect_status(gms, mpi->GetIupac().c_str());
            if ((ims == imsTrain) || (ims == imsTest))
            {
                for (k = 0; (k < qmc->nconf); k++)
                {
                    bool bExp = mpi->GetProp(mpo, iqmExp, NULL, NULL, exp_type, &exp_val, NULL);
                    bool bQM  = mpi->GetProp(mpo, iqmQM, lbuf, qmc->conf[k],
                                             qmc->type[i], &qm_val, &qm_error);
                    if (bExp && bQM)
                    {
                        fprintf(fp, "%8.3f  %8.3f  %8.3f\n", exp_val, qm_val-exp_val, qm_error);
                        diff = fabs(qm_val-exp_val);
                        if (debug &&
                            (((atoler > 0) && (diff >= atoler)) ||
                             ((exp_val != 0) && (fabs(diff/exp_val) > rtoler))))
                        {
                            fprintf(debug, "OUTLIER: %s Exp: %g, Calc: %g +/- %g\n",
                                    mpi->GetIupac().c_str(), exp_val, qm_val, qm_error);
                            nout++;
                        }
                    }
                    else if (NULL != debug)
                    {
                        fprintf(debug, "%s bQM = %d bExp = %d\n", mpi->GetMolname().c_str(),
                                bQM ? 1 : 0, bExp ? 1 : 0);
                    }
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
                                gmx_poldata_t pd,
                                gmx_bool bCalcPol,
                                MolPropObservable prop, char *exp_type,
                                char *lot,
                                real rtoler, real atoler, real outlier,
                                char *fc_str, gmx_bool bPrintAll,
                                gmx_bool bStatsTable,
                                const char *categoryfn,
                                gmx_bool bPropTable, gmx_bool bCompositionTable,
                                gmx_bool bPrintBasis, gmx_bool bPrintMultQ,
                                const char *texfn,
                                const char *xvgfn,
                                output_env_t oenv,
                                gmx_molselect_t gms,
                                const char *selout)
{
    alexandria::MolPropIterator    mpi;
    alexandria::ExperimentIterator ei;

    FILE                          *fp, *gp;
    int                            i, ntot, cur = 0;
    t_qmcount                     *qmc;
    t_refcount                    *rc;
    const char                    *molname[2];
    double                         value, error, vec[3];
    tensor                         quadrupole;
#define prev (1-cur)
    const char                    *iupac;
    char                          *ref;

    if (bCalcPol)
    {
        calc_frag_miller(pd, mp, gms);
    }

    qmc = find_calculations(mp, prop, fc_str);

    snew(rc, 1);
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        molname[cur] = mpi->GetMolname().c_str();
        for (ei = mpi->BeginExperiment(); (ei < mpi->EndExperiment()); ei++)
        {
            if (ei->GetVal(exp_type, prop, &value, &error, vec, quadrupole))
            {
                add_refc(rc, ei->GetReference().c_str());
            }
        }
        if (debug && ((mpi > mp.begin()) && (strcasecmp(molname[cur], molname[prev]) == 0)))
        {
            fprintf(debug, "Double entry %s\n", molname[cur]);
        }
        cur = prev;
    }

    printf("--------------------------------------------------\n");
    printf("      Some statistics for %s\n", mpo_name[prop]);
    ntot = 0;
    for (i = 0; (i < rc->nref); i++)
    {
        printf("There are %d experiments with %s as reference\n",
               rc->rcount[i], rc->refs[i]);
        ntot += rc->rcount[i];
    }
    printf("There are %d entries with experimental data\n", ntot);
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

    fp = gmx_ffopen(texfn, "w");
    if (bStatsTable)
    {
        gmx_molprop_stats_table(fp, prop, mp, qmc, exp_type, outlier, gms, imsTrain);
        gmx_molprop_stats_table(fp, prop, mp, qmc, exp_type, outlier, gms, imsTest);
    }
    if (bPropTable)
    {
        gmx_molprop_prop_table(fp, prop, rtoler, atoler, mp, qmc, bPrintAll, bPrintBasis,
                               bPrintMultQ, gms, imsTrain);
        gmx_molprop_prop_table(fp, prop, rtoler, atoler, mp, qmc, bPrintAll, bPrintBasis,
                               bPrintMultQ, gms, imsTest);
        if (NULL != selout)
        {
            gp = fopen(selout, "w");
            for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
            {
                iupac = mpi->GetIupac().c_str();
                if ((NULL != iupac) && (strlen(iupac) > 0))
                {
                    if (mpi->GetPropRef(prop, iqmBoth, lot, NULL, NULL,
                                        &value, &error, &ref, NULL,
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
        gmx_molprop_category_table(fp, mp, gms, imsTrain);
        fclose(fp);
    }
    if (NULL != xvgfn)
    {
        write_corr_xvg(xvgfn, mp, prop, qmc, rtoler, atoler, oenv, gms, exp_type);
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
    static char                     *prop[]      = { NULL, (char *)"potential", (char *)"dipole", (char *)"quadrupole", (char *)"polarizability", (char *)"energy", (char *)"DHf(298.15K)", NULL };
    static char                     *fc_str      = (char *)"";
    static char                     *exp_type    = (char *)"";
    static char                     *lot         = (char *)"B3LYP/aug-cc-pVTZ";
    static real                      rtoler      = 0.15, atoler = 0, outlier = 1;
    static real                      th_toler    = 170, ph_toler = 5;
    static gmx_bool                  bMerge      = TRUE, bAll = FALSE, bCalcPol = TRUE, bPrintBasis = TRUE, bPrintMultQ = FALSE;
    static gmx_bool                  bStatsTable = TRUE, bCompositionTable = FALSE, bPropTable = TRUE;
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
    gmx_poldata_t                    pd;
    output_env_t                     oenv;
    gmx_molselect_t                  gms;
    char                           **mpname = NULL, **fns = NULL;
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
    gms     = gmx_molselect_init(opt2fn("-sel", NFILE, fnm));

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
    if (NULL == (pd = gmx_poldata_read(opt2fn("-d", NFILE, fnm), ap)))
    {
        gmx_fatal(FARGS, "Can not read the force field information. File %s missing or incorrect.", fns[i]);
    }

    if (bMerge)
    {
        merge_xml(nmpfile, mpname, mp, NULL, NULL, NULL, ap, pd, TRUE);
    }
    else
    {
        MolPropRead((const char *)mpname[0], mp);
    }
    generate_composition(mp, pd);
    generate_formula(mp, ap);

    if (mpsa != MPSA_NR)
    {
        MolPropSort(mp, mpsa, ap, gms);
    }
    gmx_molprop_analyze(mp, pd, bCalcPol,
                        mpo, exp_type, lot, rtoler, atoler, outlier, fc_str, bAll,
                        bStatsTable,
                        opt2fn_null("-cat", NFILE, fnm),
                        bPropTable, bCompositionTable,
                        bPrintBasis, bPrintMultQ,
                        opt2fn("-t", NFILE, fnm), opt2fn("-c", NFILE, fnm),
                        oenv, gms,
                        opt2fn_null("-selout", NFILE, fnm));

    return 0;
}
