/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011-2016, by the GROMACS development team, led by
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
 * \author  Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author  David van der Spoel <david.vanderspoel@icm.uu.se>
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

#include "gentop_core.h"
#include "getmdlogger.h"
#include "gmx_simple_comm.h"
#include "moldip.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_xml.h"


static void print_stats(FILE *fp, const char *prop, gmx_stats_t lsq, gmx_bool bHeader,
                        char *xaxis, char *yaxis)
{
    real a, da, b, db, chi2, rmsd, Rfit;
    int  n;

    if (bHeader)
    {
        fprintf(fp, "Fitting data to y = ax+b, where x = %s and y = %s\n", xaxis, yaxis);
        fprintf(fp, "%-12s %5s %13s %13s %8s %8s\n",
                "Property", "N", "a", "b", "R", "RMSD");
        fprintf(fp, "---------------------------------------------------------------\n");
    }
    gmx_stats_get_ab(lsq, elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
    gmx_stats_get_rmsd(lsq, &rmsd);
    gmx_stats_get_npoints(lsq, &n);
    fprintf(fp, "%-12s %5d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f%% %8.4f\n",
            prop, n, a, da, b, db, Rfit*100, rmsd);
}

static void print_lsq_set(FILE *fp, gmx_stats_t lsq)
{
    real   x, y;

    fprintf(fp, "@type xy\n");
    while (gmx_stats_get_point(lsq, &x, &y, nullptr, nullptr, 0) == estatsOK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
}

static void xvgr_symbolize(FILE *xvgf, int nsym, const char *leg[],
                           const gmx_output_env_t * oenv)
{
    int i;

    xvgr_legend(xvgf, nsym, leg, oenv);
    for (i = 0; (i < nsym); i++)
    {
        xvgr_line_props(xvgf, i, elNone, ecBlack+i, oenv);
        fprintf(xvgf, "@ s%d symbol %d\n", i, i+1);
    }
}

namespace alexandria
{

class OPtimization : public MolDip
{
    using param_type = std::vector<double>;
    
    private:    
    public:
    
        OPtimization(bool  bfitESP, 
                     bool  bfitDipole, 
                     bool  bfitQuadrupole,
                     bool  bFitAlpha,
                     real  watoms, 
                     char *lot) 
           :
               bESP_(bfitESP),
               bDipole_(bfitDipole),
               bQuadrupole_(bfitQuadrupole),
               bFitAlpha_(bFitAlpha),
               watoms_(watoms),
               lot_(lot)
          {};
        
       ~OPtimization() {};
       
        param_type    param_, lower_, upper_, best_;
        param_type    orig_, psigma_, pmean_;
        bool          bESP_, bDipole_, bQuadrupole_; 
        bool          bFitAlpha_;
        real          watoms_;
        char         *lot_;
    
        double harmonic (double x, 
                         double min, 
                         double max)
       {
           return (x < min) ? (gmx::square(x-min)) : ((x > max) ? (gmx::square(x-max)) : 0);
       }
     
        void print_results(FILE *fp, 
                           const char             *qhisto,
                           const char             *dipcorr,
                           const char             *mucorr, 
                           const char             *Qcorr,
                           const char             *espcorr,
                           const char             *alphacorr, 
                           real                    dip_toler, 
                           real                    quad_toler, 
                           const gmx_output_env_t *oenv,
                           bool                    bPolar);
        
        void print_dipole(FILE  *fp, 
                          MyMol *mol, 
                          char  *calc_name, 
                          real   toler);
    
        void print_quadrapole(FILE  *fp, 
                              MyMol *mol, 
                              char  *calc_name, 
                              real   toler);
                              
        void addEspPoint();
        
        void setEEM();
        
        void polData2TuneDip();
        
        void tuneDip2PolData();
        
        void InitOpt(real factor);
                                 
        void calcDeviation();
    
        double objFunction(const double v[]);
        
        double calcPenalty(AtomIndexIterator ai);
        
        void optRun(FILE *fp, FILE *fplog, int maxiter,
                    int nrun, real stepsize, int seed,
                    const gmx_output_env_t *oenv,
                    int nprint, const char *xvgconv, 
                    const char *xvgepot, real temperature, 
                    bool bBound);
};

void OPtimization::addEspPoint()
{
    for (auto &mymol : mymol_)
    {
        if (mymol.eSupp_ == eSupportLocal)
        {
            mymol.Qgresp_.setChargeDistributionModel(iChargeDistributionModel_);
            mymol.Qgresp_.setAtomWeight(watoms_);
            mymol.Qgresp_.setAtomInfo(&mymol.topology_->atoms, pd_, mymol.state_->x, mymol.molProp()->getCharge());
            mymol.Qgresp_.setAtomSymmetry(mymol.symmetric_charges_);
            mymol.Qgresp_.setMolecularCharge(mymol.molProp()->getCharge());
            mymol.Qgresp_.summary(debug);
            
            auto ci = mymol.molProp()->getLotPropType(lot_, MPO_POTENTIAL, nullptr);
            if (ci != mymol.molProp()->EndExperiment())
            {
                size_t iesp = 0;
                for (auto epi = ci->BeginPotential(); epi < ci->EndPotential(); ++epi, ++iesp)
                {
                    if (mymol.Qgresp_.myWeight(iesp) == 0)
                    {
                        continue;
                    }
                    int xu = string2unit(epi->getXYZunit().c_str());
                    int vu = string2unit(epi->getVunit().c_str());
                    if (-1 == xu)
                    {
                        gmx_fatal(FARGS, "No such length unit '%s' for potential",
                                  epi->getXYZunit().c_str());
                    }
                    if (-1 == vu)
                    {
                        gmx_fatal(FARGS, "No such potential unit '%s' for potential",
                                  epi->getVunit().c_str());
                    }
                    mymol.Qgresp_.addEspPoint(convert2gmx(epi->getX(), xu),
                                              convert2gmx(epi->getY(), xu),
                                              convert2gmx(epi->getZ(), xu),
                                              convert2gmx(epi->getV(), vu));
                }
                if (debug)
                {
                    fprintf(debug, "Added %zu ESP points to the RESP structure.\n", mymol.Qgresp_.nEsp());
                }
            }
        }
    }
}

void OPtimization::setEEM()
{
    for (auto& mymol : mymol_)
    {
        if (mymol.eSupp_ != eSupportNo)
        {
            bool bHaveShells = false;
            if (nullptr != mymol.shellfc_)
            {
                bHaveShells = true;
            }         
            mymol.Qgeem_.setInfo(pd_, &(mymol.topology_->atoms), 
                                 iChargeDistributionModel_,
                                 hfac_,
                                 mymol.molProp()->getCharge(),
                                 bHaveShells);
        }
    }
}

void OPtimization::calcDeviation()
{
    int    j;
    double qtot = 0;

    if (PAR(cr_))
    {
        gmx_bcast(sizeof(bFinal_), &bFinal_, cr_);
    }
    if (PAR(cr_) && !bFinal_)
    {
        pd_.broadcast(cr_);
    }   
    for (j = 0; j < ermsNR; j++)
    {
        ener_[j] = 0;
    }    
    for (auto &mymol : mymol_)
    {        
        if ((mymol.eSupp_ == eSupportLocal) ||
            (bFinal_ && (mymol.eSupp_ == eSupportRemote)))
        {               
            mymol.Qgeem_.generateChargesSm(debug,
                                           pd_, 
                                           &(mymol.topology_->atoms),
                                           &mymol.chieq_,
                                           mymol.state_->x);           
            if (nullptr != mymol.shellfc_)
            {      
                if (bFitAlpha_)
                {
                    mymol.UpdateIdef(pd_, eitPOLARIZATION);  
                }             
                for (j = 0; j < mymol.topology_->atoms.nr; j++)
                {
                    mymol.mtop_->moltype[0].atoms.atom[j].q = 
                        mymol.mtop_->moltype[0].atoms.atom[j].qB = mymol.topology_->atoms.atom[j].q;     
                    
                }               
                mymol.computeForces(nullptr, cr_);
            }                                    
            qtot = 0;                       
            for (j = 0; j < mymol.topology_->atoms.nr; j++)
            {
                auto atomnr = mymol.topology_->atoms.atom[j].atomnumber;
                auto qq     = mymol.topology_->atoms.atom[j].q;
                qtot       += qq;
                if (mymol.topology_->atoms.atom[j].ptype == eptAtom)
                {
                    if (((qq < 0) && (atomnr == 1)) ||
                        ((qq > 0) && ((atomnr == 8)  || (atomnr == 9) ||
                                      (atomnr == 16) || (atomnr == 17) ||
                                      (atomnr == 35) || (atomnr == 53))))
                    {
                        ener_[ermsBOUNDS] += fabs(qq);
                    }                    
                }
            }
            if (fabs(qtot - mymol.molProp()->getCharge()) > 1e-2)
            {
                fprintf(stderr, "Warning qtot for %s is %g, should be %d\n",
                        mymol.molProp()->getMolname().c_str(),
                        qtot, mymol.molProp()->getCharge());
            }
            if (bESP_)
            {
                real  rrms   = 0;
                real  wtot   = 0;              
                if (nullptr != mymol.shellfc_)
                {
                    mymol.Qgresp_.updateAtomCoords(mymol.state_->x);
                }
                mymol.Qgresp_.updateAtomCharges(&mymol.topology_->atoms);
                mymol.Qgresp_.calcPot();
                ener_[ermsESP] += convert2gmx(mymol.Qgresp_.getRms(&wtot, &rrms), eg2cHartree_e);               
            }            
            if (bDipole_)
            {
                mymol.CalcDipole();
                if (bQM_)
                {
                    rvec dmu;                    
                    rvec_sub(mymol.mu_calc_, mymol.mu_elec_, dmu);
                    ener_[ermsMU]  += iprod(dmu, dmu);
                }
                else
                {
                    ener_[ermsMU]  += gmx::square(mymol.dip_calc_ - mymol.dip_exp_);
                }
            }                          
            if (bQuadrupole_)
            {
                mymol.CalcQuadrupole();
                for (auto mm = 0; mm < DIM; mm++)
                {
                    if (bfullTensor_)
                    {
                        for (auto nn = 0; nn < DIM; nn++)
                        {
                            ener_[ermsQUAD] += gmx::square(mymol.Q_calc_[mm][nn] - mymol.Q_elec_[mm][nn]);
                        }
                    }
                    else
                    {
                        ener_[ermsQUAD] += gmx::square(mymol.Q_calc_[mm][mm] - mymol.Q_elec_[mm][mm]);
                    }
                }
            }
        }
    }    
    if (PAR(cr_) && !bFinal_)
    {
        gmx_sum(ermsNR, ener_, cr_);
    }    
    if (MASTER(cr_))
    {
        for (j = 0; j < ermsTOT; j++)
        {
            ener_[ermsTOT] += ((fc_[j]*ener_[j])/nmol_support_);
        }
    }   
    if (nullptr != debug && MASTER(cr_))
    {
        fprintf(debug, "ENER:");
        for (j = 0; j < ermsNR; j++)
        {
            fprintf(debug, "  %8.3f", ener_[j]);
        }
        fprintf(debug, "\n");
    }    
}

void OPtimization::polData2TuneDip()
{
    param_.clear();    
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {   
        if (!ai->isConst())
        {
            auto ei   = pd_.findEem(iChargeDistributionModel_, ai->name());
            GMX_RELEASE_ASSERT(ei != pd_.EndEemprops(), "Cannot find eemprops");
            
            auto J00  = ei->getJ0();
            param_.push_back(std::move(J00));
            
            if (ai->name().compare(fixchi_) != 0)
            {
                auto Chi0 = ei->getChi0();
                param_.push_back(std::move(Chi0));
            }        
            if (bFitZeta_)
            {
                auto nzeta = ei->getNzeta();
                for (int k = 0; k < nzeta; k++)
                {
                    auto zeta = ei->getZeta(k);
                    if (0 != zeta)
                    {
                        param_.push_back(std::move(zeta));
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Zeta is zero for atom %s in %d model\n",
                                  ai->name().c_str(), iChargeDistributionModel_);
                    }
                }
            }
            if(bFitAlpha_)
            {
                auto alpha = 0.0;
                auto sigma = 0.0;
                if (pd_.getAtypePol(ai->name(), &alpha, &sigma))
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
    if (bOptHfac_)
    {
        param_.push_back(std::move(hfac_));
    }       
}

void OPtimization::tuneDip2PolData()
{
    int   n = 0;
    char  zstr[STRLEN];
    char  buf[STRLEN];
        
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto ei = pd_.findEem(iChargeDistributionModel_, ai->name());
            
            std::string qstr   = ei->getQstr();
            std::string rowstr = ei->getRowstr();
            
            if (qstr.size() == 0 || rowstr.size() == 0)
            {
                gmx_fatal(FARGS, "No qstr/rowstr for atom %s in %d model\n",
                          ai->name().c_str(), iChargeDistributionModel_);
            }
                       
            ei->setJ0(param_[n]);
            ei->setJ0_sigma(psigma_[n++]);
            
            if (ai->name().compare(fixchi_) != 0)
            {      
                ei->setChi0(param_[n]); 
                ei->setChi0_sigma(psigma_[n++]);          
            }       
            if (bFitZeta_)
            {
                zstr[0] = '\0';
                auto nzeta = ei->getNzeta();
                for (auto zz = 0; zz < nzeta; zz++)
                {
                    auto zeta = param_[n++];
                    sprintf(buf, "%g ", zeta);
                    strcat(zstr, buf);
                }                
                ei->setRowZetaQ(rowstr, zstr, qstr);
                ei->setZetastr(zstr);
            }
            if (bFitAlpha_)
            {
                std::string ptype;
                if (pd_.atypeToPtype(ai->name(), ptype))
                {
                    pd_.setPtypePolarizability(ptype, param_[n], psigma_[n]);
                    n++;
                }
                else
                {
                    gmx_fatal(FARGS, "No Ptype for atom type %s\n",
                              ai->name().c_str());
                }
            }
        }               
    }
    if (bOptHfac_)
    {
        hfac_ = param_[n++];
    }   
}

void OPtimization::InitOpt(real  factor)
{
    polData2TuneDip();

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

double OPtimization::calcPenalty(AtomIndexIterator ai)
{
    double p       = 0;
    double ref_chi = 0;
    
    ref_chi = pd_.getChi0(iChargeDistributionModel_, fixchi_); 

    auto ei      = pd_.findEem(iChargeDistributionModel_, ai->name());
    auto ai_elem = pd_.getElem(ai->name());
    auto ai_row  = ei->getRowstr();
    auto ai_chi  = ei->getChi0();
    auto ai_atn  = gmx_atomprop_atomnumber(atomprop_, ai_elem.c_str());
       
    if (ai_chi < ref_chi)
    {
        p += 1e5;
    }

    auto *ic = indexCount();
    for (auto aj = ic->beginIndex(); aj < ic->endIndex(); ++aj)
    {  
        if (!aj->isConst())
        {
            auto ej      = pd_.findEem(iChargeDistributionModel_, aj->name());
            auto aj_elem = pd_.getElem(aj->name());        
            auto aj_row  = ej->getRowstr();      
            auto aj_atn  = gmx_atomprop_atomnumber(atomprop_, aj_elem.c_str());
            
            if ((ai_row == aj_row) && (ai_atn != aj_atn))
            {           
                auto aj_chi = ej->getChi0();
                
                if (ai_atn > aj_atn)
                {
                    if (ai_chi < aj_chi)
                    {
                        p += 1e5;
                    }
                }
                else if (aj_atn > ai_atn) 
                {
                    if (aj_chi < ai_chi)
                    {
                        p += 1e5;
                    }
                }     
            }
        }
    }    
    return p;
}

double OPtimization::objFunction(const double v[])
{    
    double bounds  = 0;
    double penalty = 0;
    int    n       = 0;
    
    auto np = param_.size();
    for (size_t i = 0; i < np; i++)
    {
        param_[i] = v[i];
    }

    tuneDip2PolData();
    
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        if (!ai->isConst())
        {
            auto name = ai->name();
            auto J00  = param_[n++];
            bounds   += harmonic(J00, J0_min_, J0_max_);
            
            if (strcasecmp(name.c_str(), fixchi_) != 0)
            {
                auto Chi0 = param_[n++];
                bounds   += harmonic(Chi0, Chi0_min_, Chi0_max_);
            }
            
            if (bFitZeta_)
            {
                auto nzeta = pd_.getNzeta(iChargeDistributionModel_, ai->name());
                for (auto zz = 0; zz < nzeta; zz++)
                {
                    auto zeta = param_[n++];
                    bounds += harmonic(zeta, zeta_min_, zeta_max_);
                }
            }       
            penalty += calcPenalty(ai);
        }
    }   
    if (bOptHfac_)
    {
        hfac_ = param_[n++];
        if (hfac_ > hfac0_)
        {
            bounds += 100*gmx::square(hfac_ - hfac0_);
        }
        else if (hfac_ < -(hfac0_))
        {
            bounds += 100*gmx::square(hfac_ + hfac0_);
        }
    }      
    calcDeviation();

    ener_[ermsBOUNDS] += bounds;
    ener_[ermsTOT]    += bounds;
    ener_[ermsTOT]    += penalty;
    
    return ener_[ermsTOT];
}

void OPtimization::optRun(FILE *fp, FILE *fplog, int maxiter,
                          int nrun, real stepsize, int seed,
                          const gmx_output_env_t *oenv,
                          int nprint, const char *xvgconv, 
                          const char *xvgepot, real temperature, 
                          bool bBound)
{
    std::vector<double> optb, opts, optm;
    double              chi2, chi2_min;
    gmx_bool            bMinimum = false;
    
    auto func = [&] (const double v[]) {
        return objFunction(v);
    };
    
    if (MASTER(cr_))
    {    
        if (PAR(cr_))
        {
            for (int dest = 1; dest < cr_->nnodes; dest++)
            {
                gmx_send_int(cr_, dest, (nrun*maxiter*param_.size()));
            }
        }
        
        chi2 = chi2_min = GMX_REAL_MAX;
        Bayes <double> TuneDip(func, param_, lower_, upper_, &chi2);
        TuneDip.Init(xvgconv, xvgepot, oenv, seed, stepsize, 
                     maxiter, nprint,temperature, bBound);
                     
        for (auto n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            
            TuneDip.simulate();
            TuneDip.getBestParam(optb);
            TuneDip.getPsigma(opts);
            TuneDip.getPmean(optm);

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
            TuneDip.setParam(best_);
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
        auto niter = gmx_recv_int(cr_, 0);
        for (auto n = 0; n < niter + 2; n++)
        {
            calcDeviation();
        }
    }    
    bFinal_ = true;
    if(MASTER(cr_))
    {
        chi2 = objFunction(best_.data());;
        if (nullptr != fp)
        {
            fprintf(fp, "rmsd: %4.3f  ermsBOUNDS: %4.3f  after %d run(s)\n",
                    sqrt(chi2), ener_[ermsBOUNDS], nrun);
        }
        if (nullptr != fplog)
        {
            fprintf(fplog, "rmsd: %4.3f   ermsBOUNDS: %4.3f  after %d run(s)\n",
                    sqrt(chi2), ener_[ermsBOUNDS], nrun);
            fflush(fplog);
        }
    }
}

void OPtimization::print_quadrapole(FILE  *fp, 
                                    MyMol *mol,
                                    char  *calc_name,
                                    real   q_toler)
{
    tensor dQ;
    real   delta = 0;
    
    if (nullptr != calc_name)
    {
        if (strcmp(calc_name, (char *)"EEM") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_calc_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_calc_[XX][XX], mol->Q_calc_[XX][YY], mol->Q_calc_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_calc_[YY][YY], mol->Q_calc_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_calc_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"ESP") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_esp_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_esp_[XX][XX], mol->Q_esp_[XX][YY], mol->Q_esp_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_esp_[YY][YY], mol->Q_esp_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_esp_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"MPA") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_mulliken_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_mulliken_[XX][XX], mol->Q_mulliken_[XX][YY], mol->Q_mulliken_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_mulliken_[YY][YY], mol->Q_mulliken_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_mulliken_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"HPA") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_hirshfeld_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_hirshfeld_[XX][XX], mol->Q_hirshfeld_[XX][YY], mol->Q_hirshfeld_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_hirshfeld_[YY][YY], mol->Q_hirshfeld_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_hirshfeld_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"CM5") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_cm5_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_cm5_[XX][XX], mol->Q_cm5_[XX][YY], mol->Q_cm5_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_cm5_[YY][YY], mol->Q_cm5_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_cm5_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else
        {
            fprintf(fp, "Quadrupole analysis (6 independent components only)\n");
            fprintf(fp,
                    "QM   (%6.2f %6.2f %6.2f)\n"
                    "     (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)\n",
                    mol->Q_elec_[XX][XX], mol->Q_elec_[XX][YY], mol->Q_elec_[XX][ZZ],
                    "", mol->Q_elec_[YY][YY], mol->Q_elec_[YY][ZZ],
                    "", "", mol->Q_calc_[ZZ][ZZ]);
        }
    }
}

void OPtimization::print_dipole(FILE  *fp, 
                                MyMol *mol, 
                                char  *calc_name, 
                                real   toler)
{
    rvec dmu;
    real ndmu, cosa;
    char ebuf[32];

    if (nullptr != calc_name)
    {
        if (strcmp(calc_name, (char *)"EEM") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_calc_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_calc_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_calc_[XX], mol->mu_calc_[YY], mol->mu_calc_[ZZ], 
                    norm(mol->mu_calc_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"ESP") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_esp_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_esp_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_esp_[XX], mol->mu_esp_[YY], mol->mu_esp_[ZZ], 
                    norm(mol->mu_esp_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"MPA") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_mulliken_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_mulliken_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_mulliken_[XX], mol->mu_mulliken_[YY], mol->mu_mulliken_[ZZ], 
                    norm(mol->mu_mulliken_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"HPA") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_hirshfeld_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_hirshfeld_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_hirshfeld_[XX], mol->mu_hirshfeld_[YY], mol->mu_hirshfeld_[ZZ], 
                    norm(mol->mu_hirshfeld_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"CM5") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_cm5_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_cm5_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_cm5_[XX], mol->mu_cm5_[YY], mol->mu_cm5_[ZZ], 
                    norm(mol->mu_cm5_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else
        {
            fprintf(fp, "Dipole analysis\n");
            fprintf(fp, "QM   (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f\n",
                    mol->mu_elec_[XX], mol->mu_elec_[YY], mol->mu_elec_[ZZ], norm(mol->mu_elec_));
        }
    }
}

void OPtimization::print_results(FILE                   *fp, 
                                 const char             *qhisto,
                                 const char             *DipCorr,
                                 const char             *MuCorr, 
                                 const char             *Qcorr,
                                 const char             *EspCorr, 
                                 const char             *alphaCorr,
                                 real                    dip_toler, 
                                 real                    quad_toler, 
                                 const gmx_output_env_t *oenv,
                                 bool                    bPolar)
{
    int           i    = 0, j     = 0, n     = 0;
    int           nout = 0, mm    = 0, nn    = 0;
    real          sse  = 0, rms   = 0, sigma = 0;
    real          aver = 0, error = 0, qEEM  = 0;
    
    FILE          *dipc, *muc,  *Qc;
    FILE          *hh,   *espc, *alphac;   
        
    struct AtomTypeLsq {
        std::string atomtype;
        gmx_stats_t lsq;
    };    
    enum {
        eprEEM, eprESP, eprMPA, eprHPA, eprCM5, eprNR
    };
    
    gmx_stats_t               lsq_mu[eprNR], lsq_dip[eprNR], lsq_quad[eprNR], lsq_esp, lsq_alpha;
    const char               *eprnm[eprNR] = {"EEM", "ESP", "MPA", "HPA", "CM5"};
    std::vector<AtomTypeLsq>  lsqt;

    for (int i = 0; i < eprNR; i++)
    {
        lsq_quad[i] = gmx_stats_init();
        lsq_dip[i]  = gmx_stats_init();
        lsq_mu[i]   = gmx_stats_init();
    }
    lsq_esp     = gmx_stats_init();
    lsq_alpha   = gmx_stats_init();
    n           = 0;
    
    auto *ic = indexCount();
    for (auto ai = ic->beginIndex(); ai < ic->endIndex(); ++ai)
    {
        AtomTypeLsq k;
        k.atomtype.assign(ai->name());
        k.lsq = gmx_stats_init();
        lsqt.push_back(std::move(k));
    }
        
    for (auto &mol: mymol_)
    {
        if (mol.eSupp_ != eSupportNo)
        {
            fprintf(fp, "Molecule %d: %s Qtot: %d, Multiplicity %d\n", n+1,
                    mol.molProp()->getMolname().c_str(),
                    mol.molProp()->getCharge(),
                    mol.molProp()->getMultiplicity());
            
            mol.CalcDipole();
            print_dipole(fp, &mol, (char *)"QM",  dip_toler);
            print_dipole(fp, &mol, (char *)"EEM", dip_toler);
            print_dipole(fp, &mol, (char *)"ESP", dip_toler);
            print_dipole(fp, &mol, (char *)"MPA", dip_toler);
            print_dipole(fp, &mol, (char *)"HPA", dip_toler);
            print_dipole(fp, &mol, (char *)"CM5", dip_toler);

            sse += gmx::square(mol.dip_elec_ - mol.dip_calc_);

            mol.CalcQuadrupole();
            print_quadrapole(fp, &mol, (char *)"QM",  quad_toler);
            print_quadrapole(fp, &mol, (char *)"EEM", quad_toler);
            print_quadrapole(fp, &mol, (char *)"ESP", quad_toler);
            print_quadrapole(fp, &mol, (char *)"MPA", quad_toler);
            print_quadrapole(fp, &mol, (char *)"HPA", quad_toler);
            print_quadrapole(fp, &mol, (char *)"CM5", quad_toler);
            
            rms = mol.espRms();
            fprintf(fp,   "ESP rms: %g (Hartree/e)\n", rms);           
                        
            auto nEsp     = mol.Qgresp_.nEsp();
            auto EspPoint = mol.Qgresp_.espPoint();
            for (size_t i = 0; i < nEsp; i++)
            {
                gmx_stats_add_point(lsq_esp, gmx2convert(EspPoint[i].v(),eg2cHartree_e), gmx2convert(EspPoint[i].vCalc(), eg2cHartree_e), 0, 0);
            }
            
            gmx_stats_add_point(lsq_dip[eprEEM], mol.dip_elec_, mol.dip_calc_, 0, 0);
            gmx_stats_add_point(lsq_dip[eprESP], mol.dip_elec_, mol.dip_esp_, 0, 0);
            gmx_stats_add_point(lsq_dip[eprMPA], mol.dip_elec_, mol.dip_mulliken_, 0, 0);
            gmx_stats_add_point(lsq_dip[eprHPA], mol.dip_elec_, mol.dip_hirshfeld_, 0, 0);
            gmx_stats_add_point(lsq_dip[eprCM5], mol.dip_elec_, mol.dip_cm5_, 0, 0);
            
            for (mm = 0; mm < DIM; mm++)
            {
                gmx_stats_add_point(lsq_mu[eprEEM], mol.mu_elec_[mm], mol.mu_calc_[mm], 0, 0);
                gmx_stats_add_point(lsq_mu[eprESP], mol.mu_elec_[mm], mol.mu_esp_[mm], 0, 0);
                gmx_stats_add_point(lsq_mu[eprMPA], mol.mu_elec_[mm], mol.mu_mulliken_[mm], 0, 0);
                gmx_stats_add_point(lsq_mu[eprHPA], mol.mu_elec_[mm], mol.mu_hirshfeld_[mm], 0, 0);
                gmx_stats_add_point(lsq_mu[eprCM5], mol.mu_elec_[mm], mol.mu_cm5_[mm], 0, 0);
                
                if (bfullTensor_)
                {
                    for (nn = 0; nn < DIM; nn++)
                    {
                        gmx_stats_add_point(lsq_quad[eprEEM], mol.Q_elec_[mm][nn], mol.Q_calc_[mm][nn], 0, 0);
                        gmx_stats_add_point(lsq_quad[eprESP], mol.Q_elec_[mm][nn], mol.Q_esp_[mm][nn], 0, 0);
                        gmx_stats_add_point(lsq_quad[eprMPA], mol.Q_elec_[mm][nn], mol.Q_mulliken_[mm][nn], 0, 0);
                        gmx_stats_add_point(lsq_quad[eprHPA], mol.Q_elec_[mm][nn], mol.Q_hirshfeld_[mm][nn], 0, 0);
                        gmx_stats_add_point(lsq_quad[eprCM5], mol.Q_elec_[mm][nn], mol.Q_cm5_[mm][nn], 0, 0);
                    }
                }
                else
                {
                    gmx_stats_add_point(lsq_quad[eprEEM], mol.Q_elec_[mm][mm], mol.Q_calc_[mm][mm], 0, 0);
                    gmx_stats_add_point(lsq_quad[eprESP], mol.Q_elec_[mm][mm], mol.Q_esp_[mm][mm], 0, 0);
                    gmx_stats_add_point(lsq_quad[eprMPA], mol.Q_elec_[mm][mm], mol.Q_mulliken_[mm][nn], 0, 0);
                    gmx_stats_add_point(lsq_quad[eprHPA], mol.Q_elec_[mm][mm], mol.Q_hirshfeld_[mm][nn], 0, 0);
                    gmx_stats_add_point(lsq_quad[eprCM5], mol.Q_elec_[mm][mm], mol.Q_cm5_[mm][nn], 0, 0);

                }
            }
            if(bPolar)
            {
                mol.CalcPolarizability(10, cr_, nullptr);
                for (mm = 0; mm < DIM; mm++)
                {
                    gmx_stats_add_point(lsq_alpha, mol.alpha_elec_[mm][mm], mol.alpha_calc_[mm][mm], 0, 0);
                }
            }
            
            fprintf(fp, "Atom   Type      q_EEM     q_ESP     q_MPA     q_HPA     q_CM5       x       y       z\n");
            for (j = i = 0; j < mol.topology_->atoms.nr; j++)
            {
                if (mol.topology_->atoms.atom[j].ptype == eptAtom)
                {               
                    const char *at = *(mol.topology_->atoms.atomtype[j]);
                    if(indexCount_.isOptimized(at))
                    {
                        auto        k  = std::find_if(lsqt.begin(), lsqt.end(),
                                                      [at](const AtomTypeLsq &atlsq)
                                                      {
                                                          return atlsq.atomtype.compare(at) == 0;
                                                      });                                                 
                        if (k != lsqt.end())
                        {
                            qEEM = mol.topology_->atoms.atom[j].q;
                            if(nullptr != mol.shellfc_)
                            {
                                qEEM += mol.topology_->atoms.atom[j+1].q;
                            }
                            gmx_stats_add_point(k->lsq, mol.qESP_[i], qEEM, 0, 0);
                        } 
                        
                        fprintf(fp, "%-2d%3d  %-5s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f%8.3f%8.3f%8.3f\n",
                                mol.topology_->atoms.atom[j].atomnumber,
                                j+1,
                                *(mol.topology_->atoms.atomtype[j]),
                                qEEM, 
                                mol.qESP_[i],
                                mol.qMulliken_[i],
                                mol.qHirshfeld_[i],
                                mol.qCM5_[i],
                                mol.state_->x[j][XX], 
                                mol.state_->x[j][YY], 
                                mol.state_->x[j][ZZ]); 
                    }
                    i++;
                }
            }
            fprintf(fp, "\n");
            n++;
        }
    }

    fprintf(fp, "Dipoles are %s in EEM Parametrization.\n",     (bDipole_ ?     "used" : "not used"));
    fprintf(fp, "Quadrupoles are %s in EEM Parametrization.\n", (bQuadrupole_ ? "used" : "not used"));
    fprintf(fp, "ESP points are %s in EEM Parametrization.\n",  (bESP_ ?        "used" : "not used"));
    fprintf(fp, "\n");
    
    print_stats(fp, (char *)"Dipoles       (Debye)",       lsq_mu[eprEEM],   true,  (char *)"QM", (char *)"EEM");
    print_stats(fp, (char *)"Dipole Moment (Debye)",       lsq_dip[eprEEM],  false, (char *)"QM", (char *)"EEM");
    print_stats(fp, (char *)"Quadrupoles   (Buckingham)",  lsq_quad[eprEEM], false, (char *)"QM", (char *)"EEM");
    if (bESP_ || (!bESP_ && iChargeGenerationAlgorithm_ == eqgEEM))
    {
        print_stats(fp, (char *)"ESP (Hartree/e)", lsq_esp, false, (char *)"QM", (char *)"EEM");
    }        
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprESP],   true,  (char *)"QM", (char *)"ESP");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprESP],  false, (char *)"QM", (char *)"ESP");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprESP], false, (char *)"QM", (char *)"ESP");    
    if (!bESP_ && iChargeGenerationAlgorithm_ == eqgESP)
    {
        print_stats(fp, (char *)"ESP (Hartree/e)", lsq_esp, false, (char *)"QM", (char *)"ESP");
    }   
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprMPA],   true,  (char *)"QM", (char *)"MPA");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprMPA],  false, (char *)"QM", (char *)"MPA");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprMPA], false, (char *)"QM", (char *)"MPA");   
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprHPA],   true,  (char *)"QM", (char *)"HPA");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprHPA],  false, (char *)"QM", (char *)"HPA");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprHPA], false, (char *)"QM", (char *)"HPA");   
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprCM5],   true,  (char *)"QM", (char *)"CM5");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprCM5],  false, (char *)"QM", (char *)"CM5");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprCM5], false, (char *)"QM", (char *)"CM5");   
    fprintf(fp, "\n");

    
    std::vector<const char*> atypes;
    for (const auto &k : lsqt)
    {
        atypes.push_back(k.atomtype.c_str());
    }
    
    hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
    xvgr_legend(hh, atypes.size(), atypes.data(), oenv);
    
    fprintf(fp, "\nEEM parameters are optimized for %zu atom types:\n", atypes.size());    
    for (auto k = lsqt.begin(); k < lsqt.end(); ++k)
    {
        int   nbins;
        if (gmx_stats_get_npoints(k->lsq, &nbins) == estatsOK)
        {
            real *x, *y;
            fprintf(fp, "%-4d copies for %4s\n", nbins, k->atomtype.c_str());
            if (gmx_stats_make_histogram(k->lsq, 0, &nbins, ehistoY, 1, &x, &y) == estatsOK)
            {
                fprintf(hh, "@type xy\n");
                for (int i = 0; i < nbins; i++)
                {
                    fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                }
                fprintf(hh, "&\n");
                
                free(x);
                free(y);
            }
        }
        gmx_stats_free(k->lsq);
    }
    fclose(hh);
    fprintf(fp, "\n");
    
    dipc = xvgropen(DipCorr, "Dipole Moment (Debye)", "QM", "Empirical", oenv);
    xvgr_symbolize(dipc, 5, eprnm, oenv);
    print_lsq_set(dipc, lsq_dip[eprEEM]);
    print_lsq_set(dipc, lsq_dip[eprESP]);
    print_lsq_set(dipc, lsq_dip[eprMPA]);
    print_lsq_set(dipc, lsq_dip[eprHPA]);
    print_lsq_set(dipc, lsq_dip[eprCM5]);
    fclose(dipc);
    
    muc = xvgropen(MuCorr, "Dipoles (Debye)", "QM", "Empirical", oenv);
    xvgr_symbolize(muc, 5, eprnm, oenv);
    print_lsq_set(muc, lsq_mu[eprEEM]);
    print_lsq_set(muc, lsq_mu[eprESP]);
    print_lsq_set(muc, lsq_mu[eprMPA]);
    print_lsq_set(muc, lsq_mu[eprHPA]);
    print_lsq_set(muc, lsq_mu[eprCM5]);
    fclose(muc);

    Qc = xvgropen(Qcorr, "Quadrupoles (Buckingham)", "QM", "Empirical", oenv);
    xvgr_symbolize(Qc, 5, eprnm, oenv);
    print_lsq_set(Qc, lsq_quad[eprEEM]);
    print_lsq_set(Qc, lsq_quad[eprESP]);
    print_lsq_set(Qc, lsq_quad[eprMPA]);
    print_lsq_set(Qc, lsq_quad[eprHPA]);
    print_lsq_set(Qc, lsq_quad[eprCM5]);
    fclose(Qc);
    
    if (bESP_)
    {
        espc = xvgropen(EspCorr, "Electrostatic Potential (Hartree/e)", "QM", "EEM", oenv);
        xvgr_symbolize(espc, 1, eprnm, oenv);
        print_lsq_set(espc, lsq_esp);
        fclose(espc);
    }    
    if (bPolar)
    {
        alphac = xvgropen(alphaCorr, "Isotropic Polarizability (A^3)", "QM", "EEM", oenv);
        xvgr_symbolize(alphac, 1, eprnm, oenv);
        print_lsq_set(alphac, lsq_alpha);
        fclose(alphac);
    }

    fprintf(fp, "hfac = %g\n", hfac_);
    gmx_stats_get_ase(lsq_mu[0], &aver, &sigma, &error);
    sigma = sqrt(sse/n);
    nout  = 0;
    fprintf(fp, "Overview of dipole moment outliers (> %.3f off)\n", 2*sigma);
    fprintf(fp, "----------------------------------\n");
    fprintf(fp, "%-20s  %12s  %12s  %12s\n", "Name", "EEM", "QM", "Deviation (Debye)");
            
    for (auto &mol : mymol_)
    {
      auto deviation = std::abs(mol.dip_calc_ - mol.dip_elec_);
        if ((mol.eSupp_ != eSupportNo) &&
            (mol.dip_elec_ > sigma) &&
            (deviation > 2*sigma))
        {
            fprintf(fp, "%-20s  %12.3f  %12.3f  %12.3f\n",
                    mol.molProp()->getMolname().c_str(),
                    mol.dip_calc_, mol.dip_elec_, deviation);
            nout++;
        }
    }
    if (nout)
    {
        printf("There were %d outliers. See at the very bottom of the log file\n", nout);
    }
    else
    {
        printf("No outliers! Well done.\n");
    }

    for (int i = 0; i < eprNR; i++)
    {
        gmx_stats_free(lsq_quad[i]);
        gmx_stats_free(lsq_mu[i]);
        gmx_stats_free(lsq_dip[i]);    
    }
    gmx_stats_free(lsq_esp);
    gmx_stats_free(lsq_alpha);
    
}
}

int alex_tune_dip(int argc, char *argv[])
{
    static const char          *desc[] = {
        "tune_dip read a series of molecules and corresponding experimental",
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
        { efXVG, "-conv",      "param-conv",    ffWRITE },
        { efXVG, "-epot",      "param-epot",    ffWRITE }
    };
    
    const  int                  NFILE         = asize(fnm);

    static int                  nrun          = 1;
    static int                  nprint        = 10;
    static int                  maxiter       = 100;
    static int                  reinit        = 0;
    static int                  mindata       = 3;
    static int                  seed          = -1;
    static int                  qcycle        = 1000;
    static real                 qtol          = 1e-6;
    static real                 watoms        = 0;
    static real                 J0_min        = 5;
    static real                 Chi0_min      = 1;
    static real                 zeta_min      = 5;
    static real                 step          = 0.01;
    static real                 hfac          = 0;
    static real                 rDecrZeta     = -1;
    static real                 J0_max        = 30;
    static real                 Chi0_max      = 30;
    static real                 zeta_max      = 50;
    static real                 fc_mu         = 1;
    static real                 fc_bound      = 1;
    static real                 fc_quad       = 1;
    static real                 fc_charge     = 0;
    static real                 fc_esp        = 1;
    static real                 th_toler      = 170;
    static real                 ph_toler      = 5;
    static real                 dip_toler     = 0.5;
    static real                 quad_toler    = 5;
    static real                 factor        = 0.8;
    static real                 temperature   = 300;
    static char                *opt_elem      = nullptr;
    static char                *const_elem    = nullptr;
    static char                *fixchi        = (char *)"";
    static char                *lot           = (char *)"B3LYP/aug-cc-pVTZ";
    static gmx_bool             bRandom       = false;
    static gmx_bool             bOptHfac      = false;
    static gmx_bool             bcompress     = false;
    static gmx_bool             bQM           = false;
    static gmx_bool             bPolar        = false;
    static gmx_bool             bZPE          = false;
    static gmx_bool             bfullTensor   = false;
    static gmx_bool             bBound        = false;
    static gmx_bool             bQuadrupole   = false;
    static gmx_bool             bDipole       = false;
    static gmx_bool             bESP          = false;
    static gmx_bool             bFitZeta      = false; 
    static gmx_bool             bFitAlpha     = false;
    static gmx_bool             bZero         = true;  
    static gmx_bool             bGaussianBug  = true;     
    static const char          *cqdist[]      = {nullptr, "AXp", "AXg", "AXs", "AXpp", "AXpg", "AXps", nullptr};
    static const char          *cqgen[]       = {nullptr, "None", "EEM", "ESP", "RESP", nullptr};
    
    t_pargs                     pa[]         = {
        { "-maxiter", FALSE, etINT, {&maxiter},
          "Max number of iterations for optimization" },
        { "-mindata", FALSE, etINT, {&mindata},
          "Minimum number of data points to optimize a polarizability value" },
        { "-nprint",  FALSE, etINT, {&nprint},
          "How often to print the parameters during the simulation" },
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },          
        { "-fullTensor", FALSE, etBOOL, {&bfullTensor},
          "consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization" },        
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Model used for charge distribution" },
        { "-qgen",   FALSE, etENUM, {cqgen},
          "Algorithm used for charge generation" },
        { "-qtol",   FALSE, etREAL, {&qtol},
          "Tolerance for assigning charge generation algorithm." },
        { "-qcycle", FALSE, etINT, {&qcycle},
          "Max number of tries for optimizing the charges." },
        { "-fixchi", FALSE, etSTR,  {&fixchi},
          "Electronegativity for this atom type is fixed. Set to FALSE if you want this variable as well, but read the help text above." },
        { "-seed",   FALSE, etINT,  {&seed},
          "Random number seed. If zero, a seed will be generated." },
        { "-j0",    FALSE, etREAL, {&J0_min},
          "Minimum value that J0 (eV) can obtain in fitting" },
        { "-chi0",    FALSE, etREAL, {&Chi0_min},
          "Minimum value that Chi0 (eV) can obtain in fitting" },
        { "-z0",    FALSE, etREAL, {&zeta_min},
          "Minimum value that inverse radius (1/nm) can obtain in fitting" },
        { "-j1",    FALSE, etREAL, {&J0_max},
          "Maximum value that J0 (eV) can obtain in fitting" },
        { "-chi1",    FALSE, etREAL, {&Chi0_max},
          "Maximum value that Chi0 (eV) can obtain in fitting" },
        { "-z1",    FALSE, etREAL, {&zeta_max},
          "Maximum value that inverse radius (1/nm) can obtain in fitting" },
        { "-decrzeta", FALSE, etREAL, {&rDecrZeta},
          "Generate decreasing zeta with increasing row numbers for atoms that have multiple distributed charges. In this manner the 1S electrons are closer to the nucleus than 2S electrons and so on. If this number is < 0, nothing is done, otherwise a penalty is imposed in fitting if the Z2-Z1 < this number." },
        { "-fc_bound",    FALSE, etREAL, {&fc_bound},
          "Force constant in the penalty function for going outside the borders given with the above six options." },
        { "-fc_mu",    FALSE, etREAL, {&fc_mu},
          "Force constant in the penalty function for the magnitude of the dipole components." },
        { "-fc_quad",  FALSE, etREAL, {&fc_quad},
          "Force constant in the penalty function for the magnitude of the quadrupole components." },
        { "-fc_esp",   FALSE, etREAL, {&fc_esp},
          "Force constant in the penalty function for the magnitude of the electrostatic potential." },
        { "-fc_charge",  FALSE, etREAL, {&fc_charge},
          "Force constant in the penalty function for the magnitude of the charges with respect to the ESP charges." },
        { "-step",  FALSE, etREAL, {&step},
          "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to optimize, e.g. \"H C Br\". The other available atom types in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of atom types to include but keep constant, e.g. \"O N\". These atom types from gentop.dat are left unmodified" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-watoms", FALSE, etREAL, {&watoms},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use zero. For point+smeared charges 1 is recommended (the default)." },
        { "-fitzeta", FALSE, etBOOL, {&bFitZeta},
          "Controls whether or not the Gaussian/Slater widths are optimized." },
        { "-fitalpha", FALSE, etBOOL, {&bFitAlpha},
          "Controls whether or not the atomic polarizability is optimized." },
        { "-esp", FALSE, etBOOL, {&bESP},
          "Calibrate EEM paramters to reproduce QM electrostatic potential." },
        { "-dipole", FALSE, etBOOL, {&bDipole},
          "Calibrate EEM paramters to reproduce dipole moment." },
        { "-quadrupole", FALSE, etBOOL, {&bQuadrupole},
          "Calibrate EEM paramters to reproduce quadrupole tensor." },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule" },
        { "-hfac",  FALSE, etREAL, {&hfac},
          "Fudge factor to scale the J00 of hydrogen by (1 + hfac * qH). Default hfac is 0, means no fudging." },
        { "-opthfac",  FALSE, etBOOL, {&bOptHfac},
          "[HIDDEN]Optimize the fudge factor to scale the J00 of hydrogen (see above). If set, then [TT]-hfac[tt] set the absolute value of the largest hfac. Above this, a penalty is incurred." },
        { "-dip_toler", FALSE, etREAL, {&dip_toler},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-th_toler", FALSE, etREAL, {&th_toler},
          "Minimum angle to be considered a linear A-B-C bond" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "Maximum angle to be considered a planar A-B-C/B-C-D torsion" },
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML file" },
        { "-bgaussquad", FALSE, etBOOL, {&bGaussianBug},
          "[HIDDEN]Work around a bug in the off-diagonal quadrupole components in Gaussian" },
        { "-factor", FALSE, etREAL, {&factor},
          "Factor for generating random parameters. Parameters will be taken within the limit factor*x - x/factor" },
        { "-bound", FALSE, etBOOL, {&bBound},
          "Impose box-constrains for the optimization. Box constraints give lower and upper bounds for each parameter seperately." },
        { "-temp",    FALSE, etREAL, {&temperature},
          "'Temperature' for the Monte Carlo simulation" }
    };

    FILE                 *fp;
    gmx_output_env_t     *oenv;
    time_t                my_t;
    MolSelect             gms;
    
    t_commrec     *cr     = init_commrec(); 
    gmx::MDLogger  mdlog  = getMdLogger(cr, stdout);
    gmx_hw_info_t *hwinfo = gmx_detect_hardware(mdlog, cr, false);
    
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(cr);
        return 0;
    }
    if (MASTER(cr))
    {
        printf("There are %d threads/processes.\n", cr->nnodes);
    }    
    if (MASTER(cr))
    {
        fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fp, "# This file was created %s", ctime(&my_t));
        fprintf(fp, "# alexandria is part of G R O M A C S:\n#\n");
        fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
    }
    else
    {
        fp = nullptr;
    }    
    if (MASTER(cr))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
    }
    
    alexandria::OPtimization       opt(bESP, bDipole, bQuadrupole, bFitAlpha, watoms, lot);
    ChargeDistributionModel        iChargeDistributionModel   = name2eemtype(cqdist[0]);
    ChargeGenerationAlgorithm      iChargeGenerationAlgorithm = (ChargeGenerationAlgorithm) get_option(cqgen);
    const char                    *tabfn                      = opt2fn_null("-table", NFILE, fnm);

    if (iChargeDistributionModel == eqdAXps  &&  nullptr == tabfn)
    {
        gmx_fatal(FARGS, "Cannot generate charges with the %s charge model without a potential table. "
                  "Please supply a table file.", getEemtypeName(iChargeDistributionModel));
    }
    
    if (iChargeDistributionModel == eqdAXpp  || 
        iChargeDistributionModel == eqdAXpg  || 
        iChargeDistributionModel == eqdAXps)
    {
        bPolar = true;
    }
    
    opt.Init(cr, bQM, bGaussianBug,
             iChargeDistributionModel,
             iChargeGenerationAlgorithm,
             rDecrZeta,
             J0_min,
             Chi0_min,
             zeta_min, 
             J0_max,
             Chi0_max,
             zeta_max, fc_bound, 
             fc_mu,
             fc_quad,
             fc_charge,
             fc_esp,
             1,
             1,
             fixchi,
             bOptHfac,
             hfac, 
             bPolar,
             bFitZeta,
             hwinfo,
             bfullTensor,
             mindata);
            
    opt.Read(fp ? fp : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             opt_elem,
             const_elem,
             lot,
             gms,
             watoms,
             true,
             false,
             false,
             bPolar,
             bZPE,
             opt2fn_null("-table", NFILE, fnm),
             qcycle,
             qtol);
            
    if (nullptr != fp)
    {
        fprintf(fp, "In the total data set of %zu molecules we have:\n", opt.mymol_.size());
    }
    if (MASTER(cr))
    {
        opt.InitOpt(factor);
    }    
    if (bESP && iChargeGenerationAlgorithm != eqgESP)
    {
        opt.addEspPoint();
    }  
    if (iChargeGenerationAlgorithm != eqgEEM)
    {
        opt.setEEM();
    }
       
    opt.optRun(MASTER(cr) ? stderr : nullptr,
               fp,
               maxiter,
               nrun, step, seed,
               oenv,
               nprint,
               opt2fn("-conv", NFILE, fnm),
               opt2fn("-epot", NFILE, fnm),
               temperature,
               bBound);
               
    if (MASTER(cr))
    {
        opt.print_results(fp,  
                          opt2fn("-qhisto",    NFILE, fnm),
                          opt2fn("-dipcorr",   NFILE, fnm),
                          opt2fn("-mucorr",    NFILE, fnm),
                          opt2fn("-thetacorr", NFILE, fnm), 
                          opt2fn("-espcorr",   NFILE, fnm),
                          opt2fn("-alphacorr", NFILE, fnm),
                          dip_toler, 
                          quad_toler, 
                          oenv,
                          bPolar);
                            
        writePoldata(opt2fn("-o", NFILE, fnm), opt.pd_, bcompress);
        done_filenms(NFILE, fnm);
        gmx_ffclose(fp);
    }
    return 0;
}
