/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
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
#include "gromacs/hardware/detecthardware.h"
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

#include "alex_modules.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "mymol_low.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"
#include "tuning_utility.h"

namespace alexandria
{

class OptACM : public MolGen, Bayes
{
    using param_type = std::vector<double>;

    private:
        bool       bFullTensor_;
        bool       bFitAlpha_;
        bool       bFitZeta_;
        bool       bSameZeta_;
        bool       bFitChi_;
        bool       bUseCM5_;
        real       penalty_;

    public:

        OptACM()
            :
              bFullTensor_(false),
              bFitAlpha_(false),
              bFitZeta_(false),
              bSameZeta_(false),
              bFitChi_(true),
              bUseCM5_(false),
              penalty_(0)
        {}

        ~OptACM() {}

        bool bESP() const { return weight(ermsESP); }

        bool dipole() const { return weight(ermsMU); }

        bool quadrupole() const { return weight(ermsQUAD); }

        bool fullTensor() const { return bFullTensor_; }

        bool fitZeta() const { return bFitZeta_; }
       
        bool sameZeta() const { return bSameZeta_; }
        
        bool fitChi() const { return bFitChi_; }

        bool useCM5() const {return bUseCM5_; }

        bool penalize() const {return penalty_ > 0; }
        
        void add_pargs(std::vector<t_pargs> *pargs)
        {
            t_pargs pa[] =
            {
                { "-fullTensor", FALSE, etBOOL, {&bFullTensor_},
                  "Consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization" },
                { "-fitalpha", FALSE, etBOOL, {&bFitAlpha_},
                  "Calibrate atomic polarizability." },
                { "-fitzeta", FALSE, etBOOL, {&bFitZeta_},
                  "Calibrate orbital exponent." },
                { "-samezeta", FALSE, etBOOL, {&bSameZeta_},
                  "Use the same zeta for both the core and the shell of the Drude model." },
                { "-fitchi", FALSE, etBOOL, {&bFitChi_},
                  "Calibrate electronegativity and hardness." },
                { "-penalty", FALSE, etREAL, {&penalty_},
                  "penalty to keep the Chi0 and J0 in order." },
                { "-cm5", FALSE, etBOOL, {&bUseCM5_},
                  "Reproduce CM5 charges in fitting." },
            };
            for (size_t i = 0; i < asize(pa); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs, etuneEEM);
            Bayes::add_pargs(pargs);
        }
        
        void optionsFinished()
        {
            MolGen::optionsFinished();
            setBounds(weight(ermsBOUNDS) > 0);
        }

        double l2_regularizer (double x,
                               double min,
                               double max)
        {
            if (x < min)
            {
                return (0.5 * gmx::square(x-min));
            }
            else if (x > max)
            { 
                return (0.5 * gmx::square(x-max));
            }
            else
            { 
                return 0;
            }
        }

        void initChargeGeneration();

        /*! \brief
         *
         * Fill parameter vector based on Poldata.
         * \param[in] factor Scaling factor for parameters
         */
        void polData2TuneACM(real factor);

        /*! \brief
         * Copy the optimization parameters to the poldata structure
         * \param[in] List over the parameters that have changed.
         */
        virtual void toPolData(const std::vector<bool> &changed);

        void InitOpt(real factor);

        virtual double calcDeviation();

        double calcPenalty(AtomIndexIterator ai);

        /*! \brief
         * Do the actual optimization.
         * \param[in] fp     FILE pointer for logging
         * \param[in] fplog  FILE pointer for logging
         * \param[in] oenv   Output environment for managing xvg files etc.
         * \param[in] xvgconv Output file monitoring parameters
         * \param[in] xvgepot Output file monitoring penalty function
         * \return true if better parameters were found.
         */
        bool optRun(FILE                   *fp,
                    FILE                   *fplog,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);
};

void OptACM::initChargeGeneration()
{
    for (auto &mymol : mymols())
    {
        if (mymol.eSupp_ != eSupportNo)
        {
            // If using ESP for fitting we need to be able to compute the
            // electrostatic potential, however we always want to report it
            // so have to initialize the data anyway.
            mymol.initQgresp(poldata(), 
                             lot(),
                             watoms(), 
                             maxPot());
            // ACM is needed always as well in this program
            mymol.Qgacm_.setInfo(poldata(), 
                                 mymol.atoms_,
                                 hfac(),
                                 mymol.molProp()->getCharge());
        }
    }
}

static void dumpQ(FILE *fp, const std::string &molname,
                  t_atoms *atoms)
{
    if (fp)
    {
        fprintf(fp, "%s q:", molname.c_str());
        for(int i = 0; i < atoms->nr; i++)
        {
            fprintf(fp, " %.3f", atoms->atom[i].q);
        }
        fprintf(fp, "\n");
    }
}

double OptACM::calcDeviation()
{
    int                  n         = 0;
    double               EemRms    = 0;
    double               bound     = 0;
    double               penalty   = 0;
    std::vector<double>  qq;
    const param_type    &param     = Bayes::getParam();

    if (MASTER(commrec()))
    {
        auto *ic = indexCount();
        for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
        {
            if (!ai->isConst())
            {
                auto name = ai->name();
                
                if (bFitChi_)
                {
                    auto J00  = param[n++];
                    bound    += l2_regularizer(J00, J0Min(), J0Max());
                    if (strcasecmp(name.c_str(), fixchi()) != 0)
                    {
                        auto Chi0 = param[n++];
                        bound    += l2_regularizer(Chi0, chi0Min(), chi0Max());
                    }
                    if (penalize())
                    {
                        penalty += calcPenalty(ai);
                    }
                }
                if (bFitZeta_)
                {
                    auto nzeta = poldata()->getNzeta(ai->name());
                    for (auto zz = 0; zz < (nzeta-1); zz++)
                    {
                        auto zeta = param[n++];
                        bound += l2_regularizer(zeta, zetaMin(), zetaMax());
                    }
                }
            }
        }
        if (optHfac())
        {
            setHfac(param[n++]);
            bound += 100*gmx::square(hfacDiff());
        }
    }

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
        poldata()->broadcast_eemprop(commrec());
    }
    resetEnergies();
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (final() && (mymol.eSupp_ == eSupportRemote)))
        {
            auto q     = mymol.Qgacm_.q();
            auto natom = mymol.Qgacm_.natom();

            qq.resize(natom + 1, 0);
            for (auto i = 0; i < natom + 1; i++)
            {
                qq[i] = q[i][0];
            }

            bool converged = false;
            int  iter      = 0;
            do
            {
                // Update charges in mtop before doing
                // shell optimization.
                for (int i = 0; i < mymol.mtop_->natoms; i++)
                {
                    mymol.mtop_->moltype[0].atoms.atom[i].q =
                        mymol.mtop_->moltype[0].atoms.atom[i].qB =
                        mymol.atoms_->atom[i].q;
                }

                if (nullptr != mymol.shellfc_)
                {
                    if (bFitAlpha_)
                    {
                        mymol.UpdateIdef(poldata(), eitPOLARIZATION);
                    }
                    mymol.computeForces(nullptr, commrec());
                }

                auto qgen =  mymol.Qgacm_.generateCharges(debug,
                                                          mymol.molProp()->getMolname().c_str(),
                                                          poldata(),
                                                          mymol.atoms_,
                                                          mymol.x());
                if (qgen != eQGEN_OK)
                {
                    gmx_fatal(FARGS, "Could not generate charges for %s: %s",
                              mymol.molProp()->getMolname().c_str(),
                              mymol.Qgacm_.message());
                }
                q       = mymol.Qgacm_.q();
                EemRms  = 0;
                for (int i = 0; i < natom + 1; i++)
                {
                    EemRms   += gmx::square(qq[i] - q[i][0]);
                    qq[i]     = q[i][0];
                }
                EemRms   /= natom;
                converged = (EemRms < qtol()) || (nullptr == mymol.shellfc_);
                iter++;
            }
            while ((!converged) && (iter < qcycle()));
            for (int i = 0; i < mymol.mtop_->natoms; i++)
            {
                mymol.mtop_->moltype[0].atoms.atom[i].q      =
                    mymol.mtop_->moltype[0].atoms.atom[i].qB =
                    mymol.atoms_->atom[i].q;
            }
            if (debug)
            {
                dumpQ(debug, mymol.molProp()->getMolname(), mymol.atoms_);
            }
            if (weight(ermsCHARGE))
            {
                int    nChargeResidual = 0; // number of charge residuals added per molecule
                double ChargeResidual  = 0;
                bool   isPolarizable   = (nullptr != mymol.shellfc_);
                double qtot = 0;
                int    i, j;
                for (j = i = 0; j < mymol.atoms_->nr; j++)
                {
                    auto atomnr = mymol.atoms_->atom[j].atomnumber;
                    auto qq     = mymol.atoms_->atom[j].q;
                    qtot       += qq;
                    if (mymol.atoms_->atom[j].ptype == eptAtom ||
                        mymol.atoms_->atom[j].ptype == eptNucleus)
                    {
                        double qref = (isPolarizable ? mymol.atoms_->atom[j+1].q : 0);
                        double dq   = 0;
                        if (atomnr == 1)
                        {
                            // Penalty if qH < 0
                            dq = qq + qref;
                        }
                        else if ((atomnr == 8)  || (atomnr == 9) ||
                                 (atomnr == 16) || (atomnr == 17) ||
                                 (atomnr == 35) || (atomnr == 53))
                        {
                            // Penalty if qO > 0, therefore we reverse the sign
                            dq = -(qq + qref);
                        }
                        if (dq < 0)
                        {
                            ChargeResidual += gmx::square(dq);
                            nChargeResidual++;
                        }
                        if (useCM5())
                        {
                            ChargeResidual += gmx::square(qq + qref - mymol.chargeQM(qtCM5)[i++]);
                            nChargeResidual++;
                        }
                    }
                }
                ChargeResidual += gmx::square(qtot - mymol.molProp()->getCharge());
                nChargeResidual++;
                increaseEnergy(ermsCHARGE, (ChargeResidual/nChargeResidual));
            }
            if (weight(ermsESP))
            {
                real rrms     = 0;
                real wtot     = 0;
                real cosangle = 0;
                if (nullptr != mymol.shellfc_)
                {
                    mymol.Qgresp_.updateAtomCoords(mymol.x());
                }
                if (bFitZeta_)
                {
                    mymol.Qgresp_.updateZeta(mymol.atoms_, poldata());
                }
                mymol.Qgresp_.updateAtomCharges(mymol.atoms_);
                mymol.Qgresp_.calcPot();
                auto myRms = 
                    convert2gmx(mymol.Qgresp_.getRms(&wtot, &rrms, &cosangle), 
                                eg2cHartree_e);
                increaseEnergy(ermsESP, gmx::square(myRms));
                if (debug)
                {
                    fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                            mymol.molProp()->getMolname().c_str(), 
                            myRms, cosangle);
                }
            }
            if (weight(ermsMU))
            {
                mymol.CalcDipole();
                mymol.rotateDipole(mymol.muQM(qtCalc), mymol.muQM(qtElec));
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
            if (weight(ermsQUAD))
            {
                mymol.CalcQuadrupole();
                for (int mm = 0; mm < DIM; mm++)
                {
                    for (int nn = 0; nn < DIM; nn++)
                    {
                        if (bFullTensor_ || mm == nn)
                        {
                            increaseEnergy(ermsQUAD, gmx::square(mymol.QQM(qtCalc)[mm][nn] - mymol.QQM(qtElec)[mm][nn]));
                        }
                    }
                }
            }
        }
    }
    increaseEnergy(ermsBOUNDS, bound);
    if (penalize())
    {
        increaseEnergy(ermsPENALTY, penalty);
    }
    sumEnergies();
    printEnergies(debug);
    return energy(ermsTOT);
}

void OptACM::polData2TuneACM(real factor)
{
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto ei   = poldata()->ztype2Eem(ai->name());
            GMX_RELEASE_ASSERT(ei != poldata()->EndEemprops(), 
                               gmx::formatString("Cannot find eemprops for %s", 
                               ai->name().c_str()).c_str());
            ai->setEemProps(ei);
            if (bFitChi_)
            {
                auto J00  = ei->getJ0();
                Bayes::addParam(J00, factor);
                
                if (ai->name().compare(fixchi()) != 0)
                {
                    auto Chi0 = ei->getChi0();
                    Bayes::addParam(Chi0, factor);
                }
            }
            if (bFitZeta_)
            {
                auto nzeta = ei->getNzeta();
                auto zeta  = ei->getZeta(nzeta-1); // We only optimize zeta for shell.
                if (0 != zeta)
                {
                    Bayes::addParam(zeta, factor);
                }
                else
                {
                    gmx_fatal(FARGS, "Zeta is zero for atom %s in model %s\n",
                              ai->name().c_str(), getEemtypeName(poldata()->getChargeModel()));
                }
            }
            if (bFitAlpha_)
            {
                auto alpha = 0.0;
                auto sigma = 0.0;
                if (poldata()->getAtypePol(ai->name(), &alpha, &sigma))
                {
                    if (0 != alpha)
                    {
                        Bayes::addParam(alpha, factor);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Polarizability is zero for atom %s\n", ai->name().c_str());
                    }
                }
            }
        }
    }
    if (optHfac())
    {
        Bayes::addParam(hfac(), factor);
    }
    
}

void OptACM::toPolData(const std::vector<bool> &changed)
{
    size_t   n           = 0;
    auto     pd          = poldata();
    bool     distributed = getEemtypeDistributed(pd->getChargeModel());
    auto    *ic          = indexCount();
    auto     param       = Bayes::getParam();
    auto     psigma      = Bayes::getPsigma();
    if (psigma.empty())
    {
        psigma.resize(param.size(), 0);
    }
    Bayes::dumpParam(debug);
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto ei = ai->eemProps();
            if (bFitChi_)
            {
                ei->setJ0(param[n]);
                ei->setJ0_sigma(psigma[n++]);
                if (ai->name().compare(fixchi()) != 0)
                {
                    ei->setChi0(param[n]);
                    ei->setChi0_sigma(psigma[n++]);
                }
            }
            if (bFitZeta_)
            {
                std::string zstr, z_sig;
                std::string qstr   = ei->getQstr();
                std::string rowstr = ei->getRowstr();
                auto   nZeta  = ei->getNzeta();
                double zeta   = 0;
                double sigma  = 0;
                if (distributed)
                {
                    if (bSameZeta_ && nZeta == 2)
                    {
                        // Same zeta will be used for both core and shell
                        zeta   = param[n];
                        sigma  = psigma[n++];
                        zstr.assign(gmx::formatString("%g %g ", zeta, zeta));
                        z_sig.assign(gmx::formatString("%g %g ", sigma, sigma));
                    }
                    else
                    {
                        for (auto i = 0; i < nZeta; i++)
                        {
                            zeta   = param[n];
                            sigma  = psigma[n++];
                            zstr.assign(gmx::formatString("%g ", zeta));
                            z_sig.assign(gmx::formatString("%g ", sigma));
                        }
                    }
                }
                ei->setRowZetaQ(rowstr, zstr, qstr);
                ei->setZetastr(zstr);
                ei->setZeta_sigma(z_sig);
            }
            if (bFitAlpha_)
            {
                std::string ptype;
                if (pd->atypeToPtype(ai->name(), ptype))
                {
                    pd->setPtypePolarizability(ptype, param[n], psigma[n]);
                    n++;
                }
                else
                {
                    gmx_fatal(FARGS, "No Ptype for atom type %s\n", ai->name().c_str());
                }
            }
        }
    }
    if (optHfac())
    {
        setHfac(param[n++]);
    }
    GMX_RELEASE_ASSERT(n == changed.size(),
                       gmx::formatString("n = %zu changed.size() = %zu",
                                         n, changed.size()).c_str());
}

void OptACM::InitOpt(real factor)
{
    polData2TuneACM(factor);
}

double OptACM::calcPenalty(AtomIndexIterator ai)
{
    double         penalty = 0;
    const auto     pd      = poldata();
    auto           ei      = ai->eemProps();
    auto           ai_elem = pd->ztype2elem(ei->getName());
    auto           ai_chi  = ei->getChi0();
    auto           ai_J0   = ei->getJ0();
    auto           ai_atn  = gmx_atomprop_atomnumber(atomprop(), ai_elem.c_str());

    if (strlen(fixchi()) != 0)
    {
        const auto ref_eem  = pd->atype2Eem(fixchi());
        if (ai_chi < ref_eem->getChi0())
        {
            penalty += penalty_;
        }
    }
    
    if (ai->name() == "z_c3" && (ai_chi < 5 or ai_chi > 8))
    {
        penalty += (ai_atn * penalty_);
    }

    if (ai->name() == "z_h1" && ai_chi > 2.5)
    {
       penalty += (6 * penalty_);
    }

    auto *ic = indexCount();
    for (auto aj = ic->beginIndex(); aj < ic->endIndex(); ++aj)
    {
        if (!aj->isConst())
        {
            const auto ej      = aj->eemProps();
            const auto aj_elem = pd->ztype2elem(ej->getName());
            auto       aj_atn  = gmx_atomprop_atomnumber(atomprop(), aj_elem.c_str());

            if (ai_atn != aj_atn)
            {
                //Penalize if HeavyAtoms_chi <= H_chi or HeavyAtoms_J0 <= H_J0                
                auto aj_chi = ej->getChi0();
                auto aj_J0  = ej->getJ0();
                if ((ai_atn == 1 && aj_atn > 1  && (aj_chi <= ai_chi || aj_J0 <= ai_J0)) ||
                    (ai_atn > 1  && aj_atn == 1 && (aj_chi <= ai_chi || aj_J0 <= ai_J0)))
                {
                    penalty += (std::abs((aj_atn - ai_atn)) * penalty_);
                }
            }
        }
    }
    return penalty;
}

bool OptACM::optRun(FILE                   *fp,
                    FILE                   *fplog,
                    int                     nrun,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot)
{
    bool bMinimum = false;

    auto func = [&] (const double v[]) {
                    return Bayes::objFunction(v);
        };

    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            int niter = 3 + nrun*Bayes::maxIter()*Bayes::nParam();
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, niter);
            }
        }
        double chi2;
        Bayes::setFunc(func, &chi2);
        Bayes::setOutputFiles(xvgconv, xvgepot, oenv);
        param_type param = Bayes::getParam();
        double chi2_min  = Bayes::objFunction(param.data());
        chi2             = chi2_min;

        for (auto n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            Bayes::MCMC();
            if (chi2 < chi2_min)
            {
                bMinimum = true;
                chi2_min = chi2;
            }
        }
        if (bMinimum)
        {
            auto best   = Bayes::getBestParam();
            auto pmean  = Bayes::getPmean();
            auto psigma = Bayes::getPsigma();
            auto emin   = Bayes::objFunction(best.data());
            if (fplog)
            {
                fprintf(fplog, "\nMinimum RMSD value during optimization: %.3f.\n", sqrt(emin));
                fprintf(fplog, "Statistics of parameters after optimization\n");
                for (size_t k = 0; k < Bayes::nParam(); k++)
                {
                    fprintf(fplog, "Parameter %3zu  Best value:%10g  Mean value:%10g  Sigma:%10g\n",
                            k, best[k], pmean[k], psigma[k]);
                }
            }
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        auto niter = gmx_recv_int(commrec(), 0);
        for (auto n = 0; n < niter; n++)
        {
            (void) calcDeviation();
        }
    }
    setFinal();
    if (MASTER(commrec()))
    {
        param_type best = Bayes::getBestParam();
        (void) Bayes::objFunction(best.data());
        printEnergies(fp);
        printEnergies(fplog);
    }
    return bMinimum;
}
}

int alex_tune_eem(int argc, char *argv[])
{
    static const char          *desc[] = {
        "tune_eem read a series of molecules and corresponding experimental",
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
        "A suitable value would be 2.3, the original, value due to Pauling,",
        "this can by overridden by setting the [tt]-fixchi[TT] flag to something else (e.g. a non-existing atom).[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm                    fnm[] = {
        { efDAT, "-f",         "allmols",       ffREAD  },
        { efDAT, "-d",         "gentop",        ffOPTRD },
        { efDAT, "-o",         "tunedip",       ffWRITE },
        { efDAT, "-sel",       "molselect",     ffREAD  },
        { efXVG, "-table",     "table",         ffOPTRD },
        { efLOG, "-g",         "charges",       ffWRITE },
        { efXVG, "-qhisto",    "q_histo",       ffWRITE },
        { efXVG, "-dipcorr",   "dip_corr",      ffWRITE },
        { efXVG, "-mucorr",    "mu_corr",       ffWRITE },
        { efXVG, "-thetacorr", "theta_corr",    ffWRITE },
        { efXVG, "-espcorr",   "esp_corr",      ffWRITE },
        { efXVG, "-alphacorr", "alpha_corr",    ffWRITE },
        { efXVG, "-qcorr",     "q_corr",        ffWRITE },
        { efXVG, "-isopol",    "isopol_corr",   ffWRITE },
        { efXVG, "-anisopol",  "anisopol_corr", ffWRITE },
        { efXVG, "-conv",      "param-conv",    ffWRITE },
        { efXVG, "-epot",      "param-epot",    ffWRITE },
        { efTEX, "-latex",     "eemprop",       ffWRITE }
    };

    const int NFILE         = asize(fnm);

    int       nrun          = 1;
    int       reinit        = 0;
    real      th_toler      = 170;
    real      ph_toler      = 5;
    real      dip_toler     = 0.5;
    real      quad_toler    = 5;
    real      alpha_toler   = 3;
    real      isopol_toler  = 2;
    real      factor        = 0.8;
    real      efield        = 10;
    char     *opt_elem      = nullptr;
    bool      bRandom       = false;
    bool      bcompress     = false;
    bool      bPrintTable   = false;
    bool      bZero         = true;
    bool      bOptimize     = true;
    bool      bForceOutput  = true;

    t_pargs                     pa[]         = {
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to optimize, e.g. \"H C Br\". The other available atom types in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-dip_toler", FALSE, etREAL, {&dip_toler},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-alpha_toler", FALSE, etREAL, {&alpha_toler},
          "Tolerance (A^3) for marking diagonal elements of the polarizability tensor as an outlier in the log file" },
        { "-isopol_toler", FALSE, etREAL, {&isopol_toler},
          "Tolerance (A^3) for marking isotropic polarizability as an outlier in the log file" },
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
          "The magnitude of the external electeric field to calculate polarizability tensor." },
        { "-optimize",     FALSE, etBOOL, {&bOptimize},
          "Do parameter optimization when true, or a single calculation otherwise." },
        { "-force_output", FALSE, etBOOL, {&bForceOutput},
          "Write output even if no new minimum is found" }
          
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
    alexandria::OptACM opt;
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
             gms,
             true,
             false,
             false,
             false,
             opt.fitZeta(),
             false,
             tabfn);

    opt.initChargeGeneration();

    bool bMinimum = false;
    if (bOptimize)
    {
        if (MASTER(opt.commrec()))
        {
            opt.InitOpt(factor);
        }
        
        bMinimum = opt.optRun(MASTER(opt.commrec()) ? stderr : nullptr,
                              fp,
                              nrun,
                              oenv,
                              opt2fn("-conv", NFILE, fnm),
                              opt2fn("-epot", NFILE, fnm));                   
    }

    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput)
        {
            auto iModel = opt.poldata()->getChargeModel();
            bool bPolar = getEemtypePolarizable(iModel);
            
            auto *ic = opt.indexCount();
            print_electric_props(fp,
                                 opt.mymols(),
                                 opt.poldata(),
                                 opt.mdlog(),
                                 opt.atomprop(),
                                 opt.watoms(),
                                 opt.hfac(),
                                 opt.lot(),
                                 tabfn,
                                 opt.hwinfo(),
                                 opt.qcycle(),
                                 opt.maxPot(),
                                 opt.qtol(),                          
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
                                 isopol_toler,
                                 oenv,
                                 bPolar,
                                 opt.dipole(),
                                 opt.quadrupole(),
                                 opt.fullTensor(),
                                 ic,
                                 opt.commrec(),
                                 efield);
            writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), bcompress);
            if (bPrintTable)
            {
                FILE        *tp;
                tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
                alexandria_poldata_eemprops_table(tp, opt.poldata());
                gmx_ffclose(tp);
            }
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
        gmx_ffclose(fp);
    }
    return 0;
}
