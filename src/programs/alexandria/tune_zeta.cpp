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

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "gentop_core.h"
#include "getmdlogger.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"

namespace alexandria
{

class OptZeta : public MolGen
{
    using param_type = std::vector<double>;

    private:
        gmx_bool       bDipole_;
        gmx_bool       bQuadrupole_;
        gmx_bool       bFullTensor_;
        gmx_bool       bCharge_;
        gmx_bool       bFitAlpha_;
        
        Bayes <double> TuneZeta_;
        param_type     param_, lower_, upper_, best_;
        param_type     orig_, psigma_, pmean_;

    public:

        OptZeta() 
            : 
                bDipole_(false), 
                bQuadrupole_(false), 
                bFullTensor_(false), 
                bCharge_(false), 
                bFitAlpha_(false)
            {}

        ~OptZeta() {}

        gmx_bool dipole() const { return bDipole_; }

        gmx_bool quadrupole() const { return bQuadrupole_; }

        gmx_bool fullTensor() const { return bFullTensor_; }

        void add_pargs(std::vector<t_pargs> *pargs)
        {
            t_pargs pa[] =
            {
                { "-fullTensor", FALSE, etBOOL, {&bFullTensor_},
                  "consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization" },
                { "-dipole", FALSE, etBOOL, {&bDipole_},
                  "Calibrate parameters to reproduce dipole moment." },
                { "-quadrupole", FALSE, etBOOL, {&bQuadrupole_},
                  "Calibrate parameters to reproduce quadrupole tensor." },
                { "-fitalpha", FALSE, etBOOL, {&bFitAlpha_},
                  "Calibrate atomic polarizability." },
                { "-charge", FALSE, etBOOL, {&bCharge_},
                  "Calibrate parameters to keep reasonable charges (do not use with ESP)." },
            };
            for (size_t i = 0; i < asize(pa); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs);
            TuneZeta_.add_pargs(pargs);
        }

        double l2_regularizer (double x,
                               double min,
                               double max)
        {
            return (x < min) ? (0.5 * gmx::square(x-min)) : ((x > max) ? (0.5 * gmx::square(x-max)) : 0);
        }

        void polData2TuneZeta();

        void tuneZeta2PolData();

        void InitOpt(real factor);

        void calcDeviation();

        double objFunction(const double v[]);

        void optRun(FILE                   *fp,
                    FILE                   *fplog,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);
};

void OptZeta::polData2TuneZeta()
{
    param_.clear();
    auto          *ic = indexCount();
    const Poldata &pd = poldata();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto ei    = pd.findEem(iChargeDistributionModel(), ai->name());
            GMX_RELEASE_ASSERT(ei != pd.EndEemprops(), "Cannot find eemprops");
            auto nzeta = ei->getNzeta();
            auto zeta  = ei->getZeta(nzeta-1);
            if (0 != zeta)
            {
                param_.push_back(std::move(zeta));
            }
            else
            {
                gmx_fatal(FARGS, "Zeta is zero for atom %s in model %s\n",
                          ai->name().c_str(), getEemtypeName(iChargeDistributionModel()));
            }

            if (bFitAlpha_)
            {
                auto alpha = 0.0;
                auto sigma = 0.0;
                if (pd.getAtypePol(ai->name(), &alpha, &sigma))
                {
                    if (0 != alpha)
                    {
                        param_.push_back(std::move(alpha));
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Polarizability is zero for atom %s\n",
                                  ai->name().c_str());
                    }
                }
            }
        }
    }
}

void OptZeta::tuneZeta2PolData()
{
    int      n = 0;
    char     zstr[STRLEN];
    char     z_sig[STRLEN];
    char     buf[STRLEN];
    char     buf_sig[STRLEN];

    Poldata &pd    = poldata();
    auto    *ic    = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto        ei = pd.findEem(iChargeDistributionModel(), ai->name());
            GMX_RELEASE_ASSERT(ei != pd.EndEemprops(), "Cannot find eemprops");
            std::string qstr   = ei->getQstr();
            std::string rowstr = ei->getRowstr();
            zstr[0]  = '\0';
            z_sig[0] = '\0';            
            if (iChargeDistributionModel() == eqdAXps || iChargeDistributionModel() == eqdAXpg)
            {
                auto nzeta   = ei->getNzeta();
                double zeta  = ei->getZeta(0);
                double sigma = 0;
                for (auto i = 0; i < nzeta; i++)
                {
                    if (i > 0)
                    {
                        zeta   = param_[n];
                        sigma  = psigma_[n++];
                    }
                    sprintf(buf, "%g ", zeta);
                    sprintf(buf_sig, "%g ", sigma);
                    strcat(zstr, buf);
                    strcat(z_sig, buf_sig);
                }
                ei->setRowZetaQ(rowstr, zstr, qstr);
                ei->setZetastr(zstr);
                ei->setZeta_sigma(z_sig);
            }
            else
            {
                auto nzeta  = ei->getNzeta();
                auto zeta   = param_[n];
                auto sigma  = psigma_[n++];
                for (auto i = 0; i < nzeta; i++)
                {
                    sprintf(buf, "%g ", zeta);
                    sprintf(buf_sig, "%g ", sigma);
                    strcat(zstr, buf);
                    strcat(z_sig, buf_sig);
                }
                ei->setRowZetaQ(rowstr, zstr, qstr);
                ei->setZetastr(zstr);
                ei->setZeta_sigma(z_sig);
            }
            if (bFitAlpha_)
            {
                std::string ptype;
                if (pd.atypeToPtype(ai->name(), ptype))
                {
                    pd.setPtypePolarizability(ptype, param_[n], psigma_[n]);
                    n++;
                }
                else
                {
                    gmx_fatal(FARGS, "No Ptype for atom type %s\n", ai->name().c_str());
                }
            }
        }
    }
}

void OptZeta::InitOpt(real  factor)
{
    polData2TuneZeta();

    orig_.resize(param_.size(), 0);
    best_.resize(param_.size(), 0);
    lower_.resize(param_.size(), 0);
    upper_.resize(param_.size(), 0);
    psigma_.resize(param_.size(), 0);
    pmean_.resize(param_.size(), 0);
    if (factor < 1)
    {
        factor = 1/factor;
    }
    for (size_t i = 0; (i < param_.size()); i++)
    {
        best_[i]  = orig_[i] = param_[i];
        lower_[i] = orig_[i]/factor;
        upper_[i] = orig_[i]*factor;
    }
}

void OptZeta::calcDeviation()
{
    int    j;
    double qtot = 0;
    real   rrms = 0;
    real   wtot = 0;

    if (PAR(commrec()))
    {
        bool bFinal = final();
        gmx_bcast(sizeof(final()), &bFinal, commrec());
        if (bFinal)
        {
            setFinal();
        }
    }
    if (PAR(commrec()) && !final())
    {
        poldata().broadcast(commrec());
    }
    resetEnergies();
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (final() && (mymol.eSupp_ == eSupportRemote)))
        {
            mymol.Qgresp_.updateZeta(&mymol.topology_->atoms, poldata());
            mymol.Qgresp_.optimizeCharges();
            if (nullptr != mymol.shellfc_)
            {
                if (bFitAlpha_)
                {
                    mymol.UpdateIdef(poldata(), eitPOLARIZATION);
                }
                mymol.computeForces(nullptr, commrec());
                mymol.Qgresp_.updateAtomCoords(mymol.x());
            }
            if (bCharge_)
            {
                qtot = 0;
                for (j = 0; j < mymol.topology_->atoms.nr; j++)
                {
                    auto atomnr = mymol.topology_->atoms.atom[j].atomnumber;
                    auto qq     = mymol.Qgresp_.getAtomCharge(j);
                    qtot       += qq;
                    mymol.mtop_->moltype[0].atoms.atom[j].q              =
                        mymol.mtop_->moltype[0].atoms.atom[j].qB         =
                            mymol.topology_->atoms.atom[j].q             =
                                mymol.topology_->atoms.atom[j].qB        = qq;
                    if (mymol.topology_->atoms.atom[j].ptype == eptAtom ||
                        mymol.topology_->atoms.atom[j].ptype == eptNucleus)
                    {
                        auto q_H        = 0 ? (nullptr != mymol.shellfc_) : 1;
                        auto q_OFSClBrI = 0 ? (nullptr != mymol.shellfc_) : 2;
                        if (((qq < q_H) && (atomnr == 1)) ||
                            ((qq > q_OFSClBrI) && ((atomnr == 8)  || (atomnr == 9) ||
                                                   (atomnr == 16) || (atomnr == 17) ||
                                                   (atomnr == 35) || (atomnr == 53))))
                        {
                            increaseEnergy(ermsCHARGE, gmx::square(qq));
                        }
                    }

                }
                increaseEnergy(ermsCHARGE,
                         gmx::square(qtot - mymol.molProp()->getCharge()));
            }
            mymol.Qgresp_.calcPot();
            increaseEnergy(ermsESP,
                     convert2gmx(mymol.Qgresp_.getRms(&wtot, &rrms), eg2cHartree_e));
            if (bDipole_)
            {
                mymol.CalcDipole();
                if (bQM())
                {
                    rvec dmu;
                    rvec_sub(mymol.muQM(qtCalc), mymol.muQM(qtElec), dmu);
                    increaseEnergy(ermsMU, iprod(dmu, dmu));
                }
                else
                {
                    increaseEnergy(ermsMU, gmx::square(mymol.dipQM(qtCalc) - mymol.dipExper()));
                }
            }
            if (bQuadrupole_)
            {
                mymol.CalcQuadrupole();
                for (auto mm = 0; mm < DIM; mm++)
                {
                    for (auto nn = 0; nn < DIM; nn++)
                    {
                        if (fullTensor() || mm == nn)
                        {
                            increaseEnergy(ermsQUAD, 
                                     gmx::square(mymol.QQM(qtCalc)[mm][nn] - mymol.QQM(qtElec)[mm][nn]));
                        }
                    }
                }
            }
        }
    }
    sumEnergies();
    normalizeEnergies();
    printEnergies(debug);
}

double OptZeta::objFunction(const double v[])
{
    double bounds  = 0;
    int    n       = 0;

    auto   np = param_.size();
    for (size_t i = 0; i < np; i++)
    {
        param_[i] = v[i];
    }
    tuneZeta2PolData();
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto zeta = param_[n++];
            bounds   += l2_regularizer(zeta, zetaMin(), zetaMax());

        }
    }
    calcDeviation();
    increaseEnergy(ermsBOUNDS, bounds);
    increaseEnergy(ermsTOT, bounds);
    return energy(ermsTOT);
}

void OptZeta::optRun(FILE                   *fp,
                     FILE                   *fplog,
                     int                     nrun,
                     const gmx_output_env_t *oenv,
                     const char             *xvgconv,
                     const char             *xvgepot)
{
    std::vector<double> optb, opts, optm;
    double              chi2, chi2_min;
    gmx_bool            bMinimum = false;

    auto  func = [&] (const double v[]) {
        return objFunction(v);
    };
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, (nrun*TuneZeta_.maxIter()*param_.size()));
            }
        }
        chi2 = chi2_min = GMX_REAL_MAX;
        TuneZeta_.setFunc(func, param_, lower_, upper_, &chi2);
        TuneZeta_.Init(xvgconv, xvgepot, oenv);

        for (auto n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            TuneZeta_.simulate();
            TuneZeta_.getBestParam(optb);
            TuneZeta_.getPsigma(opts);
            TuneZeta_.getPmean(optm);
            if (chi2 < chi2_min)
            {
                bMinimum = true;
                for (size_t k = 0; k < param_.size(); k++)
                {
                    best_[k]   = optb[k];
                    pmean_[k]  = optm[k];
                    psigma_[k] = opts[k];
                }
                chi2_min = chi2;
            }
            TuneZeta_.setParam(best_);
        }
        if (bMinimum)
        {
            param_    = best_;
            auto emin = objFunction(best_.data());
            if (fplog)
            {
                fprintf(fplog, "\nMinimum rmsd value during optimization: %.3f.\n", sqrt(emin));
                fprintf(fplog, "Statistics of parameters after optimization\n");
                for (size_t k = 0; k < param_.size(); k++)
                {
                    fprintf(fplog, "Parameter %3zu  Best value:%10g  Mean value:%10g  Sigma:%10g\n",
                            k, best_[k], pmean_[k], psigma_[k]);
                }
            }
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        auto niter = gmx_recv_int(commrec(), 0);
        for (auto n = 0; n < niter + 2; n++)
        {
            calcDeviation();
        }
    }
    setFinal();
    if (MASTER(commrec()))
    {
        chi2 = objFunction(best_.data());
        printEnergies(fp);
        printEnergies(fplog);
    }
}
}

int alex_tune_zeta(int argc, char *argv[])
{
    const char                 *desc[] = {
        "tune_zeta reads a series of molecules and corresponding experimental",
        "dipole moments from a file, and tunes parameters in an algorithm",
        "until the experimental dipole moments are reproduced by the",
        "charge generating algorithm AX as implemented in the gentop program.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-random[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run. At reinit steps however, the parameters are only changed",
        "slightly, in order to speed-up local search but not global search."
        "In other words, complete random starts are done only at the beginning of each",
        "run, and only when explicitly requested.[PAR]",
        "The absolut dipole moment of a molecule remains unchanged if all the",
        "atoms swap the sign of the charge. To prevent this kind of mirror",
        "effects a penalty is added to the square deviation ",
        "if hydrogen atoms have a negative charge. Similarly a penalty is",
        "added if atoms from row VI or VII in the periodic table have a positive",
        "charge. The penalty is equal to the force constant given on the command line",
        "time the square of the charge.[PAR]",
        "One of the electronegativities (chi) is redundant in the optimization,",
        "only the relative values are meaningful.",
        "Therefore by default we fix the value for hydrogen to what is written",
        "in the eemprops.dat file (or whatever is given with the [tt]-d[TT] flag).",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm                    fnm[] = {
        { efDAT, "-f",         "allmols",           ffREAD  },
        { efDAT, "-d",         "gentop",            ffOPTRD },
        { efDAT, "-o",         "tunezeta",          ffWRITE },
        { efDAT, "-sel",       "molselect",         ffREAD  },
        { efXVG, "-table",     "table",             ffOPTRD },
        { efLOG, "-g",         "charges",           ffWRITE },
        { efXVG, "-qhisto",    "q_histo",           ffWRITE },
        { efXVG, "-dipcorr",   "dip_corr",          ffWRITE },
        { efXVG, "-mucorr",    "mu_corr",           ffWRITE },
        { efXVG, "-thetacorr", "theta_corr",        ffWRITE },
        { efXVG, "-espcorr",   "esp_corr",          ffWRITE },
        { efXVG, "-alphacorr", "alpha_corr",        ffWRITE },
        { efXVG, "-qcorr",     "q_corr",            ffWRITE },
        { efXVG, "-isopol",    "isopol_corr",       ffWRITE },
        { efXVG, "-anisopol",  "anisopol_corr",     ffWRITE },
        { efXVG, "-conv",      "param-conv",        ffWRITE },
        { efXVG, "-epot",      "param-epot",        ffWRITE },
        { efTEX, "-latex",     "orbital_exponents", ffWRITE }
    };

    const int                   NFILE         = asize(fnm);

    int                         nrun          = 1;
    int                         reinit        = 0;
    real                        th_toler      = 170;
    real                        ph_toler      = 5;
    real                        dip_toler     = 0.5;
    real                        quad_toler    = 5;
    real                        alpha_toler   = 3;
    real                        factor        = 0.8;
    real                        efield        = 1;
    char                       *opt_elem      = nullptr;
    char                       *const_elem    = nullptr;
    gmx_bool                    bcompress     = false;
    gmx_bool                    bZPE          = false;
    gmx_bool                    bZero         = true;
    gmx_bool                    bPrintTable   = false;
    gmx_bool                    bOptimize     = true;

    t_pargs                     pa[]          = {
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to optimize, e.g. \"H C Br\". The other available atom types in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of atom types to include but keep constant, e.g. \"O N\". These atom types from gentop.dat are left unmodified" },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule" },
        { "-dip_toler", FALSE, etREAL, {&dip_toler},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-alpha_toler", FALSE, etREAL, {&alpha_toler},
          "Tolerance (A^3) for marking polarizability as an outlier in the log file" },
        { "-th_toler", FALSE, etREAL, {&th_toler},
          "Minimum angle to be considered a linear A-B-C bond" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "Maximum angle to be considered a planar A-B-C/B-C-D torsion" },
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML file" },
        { "-btex", FALSE, etBOOL, {&bPrintTable},
          "[HIDDEN]Print the latex table for the Gaussian and Slater exponents" },
        { "-factor", FALSE, etREAL, {&factor},
          "Factor for generating random parameters. Parameters will be taken within the limit factor*x - x/factor" },
        { "-efield",  FALSE, etREAL, {&efield},
          "The magnitude of the external electric field to calculate polarizability tensor." },
        { "-optimize",     FALSE, etBOOL, {&bOptimize},
          "Optimize zeta values" },
    };

    FILE                       *fp;
    gmx_output_env_t           *oenv;
    time_t                      my_t;
    MolSelect                   gms;

    std::vector<t_pargs>        pargs;
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::OptZeta opt;
    opt.add_pargs(&pargs);

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           pargs.size(), pargs.data(),
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    opt.optionsFinished();

    if (MASTER(opt.commrec()))
    {
        fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fp, "# This file was created %s", ctime(&my_t));
        fprintf(fp, "# alexandria is part of GROMACS:\n#\n");
        fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
    }
    else
    {
        fp = nullptr;
    }
    if (MASTER(opt.commrec()))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
    }

    const char *tabfn = opt2fn_null("-table", NFILE, fnm);

    opt.Read(fp ? fp : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             opt_elem,
             const_elem,
             gms,
             true,
             false,
             false,
             bZPE,
             true,
             tabfn);

    if (nullptr != fp)
    {
        fprintf(fp, "In the total data set of %zu molecules we have:\n",
                opt.mymols().size());
    }
    
    if (bOptimize)
    {
        if (MASTER(opt.commrec()))
        {
            opt.InitOpt(factor);
        }

        opt.optRun(MASTER(opt.commrec()) ? stderr : nullptr,
                   fp,
                   nrun,
                   oenv,
                   opt2fn("-conv", NFILE, fnm),
                   opt2fn("-epot", NFILE, fnm));
    }
    
    if (MASTER(opt.commrec()))
    {
        gmx_bool bPolar = (opt.iChargeDistributionModel() == eqdAXpp  ||
                           opt.iChargeDistributionModel() == eqdAXpg  ||
                           opt.iChargeDistributionModel() == eqdAXps);
                           
        auto *ic = opt.indexCount();
        print_electric_props(fp,
                             opt.mymols(),
                             opt.poldata(),
                             opt2fn("-qhisto",    NFILE, fnm),
                             opt2fn("-dipcorr",   NFILE, fnm),
                             opt2fn("-mucorr",    NFILE, fnm),
                             opt2fn("-thetacorr", NFILE, fnm),
                             opt2fn("-espcorr",   NFILE, fnm),
                             opt2fn("-alphacorr", NFILE, fnm),
                             opt2fn("-isopol",    NFILE, fnm),
                             opt2fn("-anisopol",  NFILE, fnm),
                             opt2fn("-qcorr",     NFILE, fnm),
                             dip_toler,
                             quad_toler,
                             alpha_toler,
                             oenv,
                             bPolar,
                             opt.dipole(),
                             opt.quadrupole(),
                             opt.fullTensor(),
                             ic,
                             opt.commrec(),
                             efield);

        writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), bcompress);
        gmx_ffclose(fp);
        if (bPrintTable)
        {
            FILE        *tp;
            tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
            alexandria_poldata_eemprops_table(tp, true, false, opt.poldata());
            gmx_ffclose(tp);
        }
        done_filenms(NFILE, fnm);
    }
    return 0;
}
