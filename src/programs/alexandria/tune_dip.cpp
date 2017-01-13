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
        fprintf(fp, "Fitting data to y = ax+b, where x = %s and y = %s\n",
                xaxis, yaxis);
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
    
        OPtimization() {};
        
       ~OPtimization() {};
       
        param_type    param_, lower_, upper_, best_;
        param_type    orig_, psigma_, pmean_;
          
        double harmonic (double x, 
                         double min, 
                         double max)
       {
           return (x < min) ? (gmx::square(x-min)) : ((x > max) ? (gmx::square(x-max)) : 0);
       }
     
        void print_molecules(FILE *fp, const char *xvgfn, const char *qhisto,
                             const char *cdiff, const char *mudiff, const char *Qdiff,
                             const char *espdiff, real dip_toler, real quad_toler, 
                             real q_toler, const gmx_output_env_t * oenv);
        
        void print_dipole(FILE  *fp, 
                          MyMol *mol, 
                          char  *calc_name, 
                          real   toler);
    
        void print_quadrapole(FILE  *fp, 
                              MyMol *mol, 
                              char  *calc_name, 
                              real   toler);
        
        void polData2TuneDip();
        
        void tuneDip2PolData();
        
        void InitOpt(real factor);
        
        void split_shell_charges(gmx_mtop_t *mtop,
                                 t_idef     *idef,
                                 t_topology *topology);
                                 
        void calcDeviation();
    
        double objFunction(const double v[]);
        
        void optRun(FILE *fp, FILE *fplog, int maxiter,
                    int nrun, real stepsize, int seed,
                    const gmx_output_env_t *oenv,
                    int nprint, const char *xvgconv, 
                    const char *xvgepot, real temperature, 
                    bool bBound);
};

void OPtimization::split_shell_charges(gmx_mtop_t *mtop,
                                       t_idef     *idef,
                                       t_topology *topology)
{
    int                     k, ai, aj;
    real                    q, Z;
    gmx_mtop_atomloop_all_t aloop;
    const t_atom           *atom;
    t_atom                 *atom_i, *atom_j;
    int                     at_global;

    for (k = 0; k < idef->il[F_POLARIZATION].nr; )
    {
        k++; // Skip over the type.
        ai = idef->il[F_POLARIZATION].iatoms[k++];
        aj = idef->il[F_POLARIZATION].iatoms[k++];

        atom_i = &topology->atoms.atom[ai];
        atom_j = &topology->atoms.atom[aj];

        if ((atom_i->ptype == eptAtom) &&
            (atom_j->ptype == eptShell))
        {
            q         = atom_i->q;
            Z         = atom_i->atomnumber;
            atom_i->q = Z;
            atom_j->q = q-Z;
        }
        else if ((atom_j->ptype == eptAtom) &&
                 (atom_i->ptype == eptShell))
        {
            q         = atom_j->q;
            Z         = atom_j->atomnumber;
            atom_j->q = Z;
            atom_i->q = q-Z;
        }
        else
        {
            gmx_incons("Polarization entry does not have one atom and one shell");
        }
    }
    q     = 0;
    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
    {
        q += atom->q;
    }
    Z = std::lround(q);
    if (fabs(q-Z) > 1e-3)
    {
        gmx_fatal(FARGS, "Total charge in molecule is not zero, but %f", q-Z);
    }
}

void OPtimization::calcDeviation()
{
    int  j;
    bool bHaveShells = false;

    if (PAR(_cr))
    {
        gmx_bcast(sizeof(_bFinal), &_bFinal, _cr);
    }

    if (PAR(_cr) && !_bFinal)
    {
        pd_.broadcast(_cr);
    }
    
    for (j = 0; j < ermsNR; j++)
    {
        _ener[j] = 0;
    }
    
    for (auto &mymol : _mymol)
    {
        GMX_RELEASE_ASSERT(mymol.mtop_->natoms == mymol.topology_->atoms.nr, "Inconsistency 3 in moldip.cpp");
        
        if ((mymol.eSupp == eSupportLocal) ||
            (_bFinal && (mymol.eSupp == eSupportRemote)))
        {
            if (nullptr != mymol.shellfc_)
            {
                bHaveShells = true;
                mymol.computeForces(nullptr, _cr);
            }
            
            QgenEem qgen(pd_, &(mymol.topology_->atoms), 
                         _iChargeDistributionModel,
                         _hfac,
                         mymol.molProp()->getCharge(),
                         bHaveShells);

            double chieq = 0;
            qgen.generateChargesSm(debug,
                                   pd_, &(mymol.topology_->atoms),
                                   &chieq,
                                   mymol.x_);
            mymol.chieq = chieq;
            
            double qtot = 0;            
            mymol.CalcMultipoles();
            for (j = 0; j < mymol.topology_->atoms.nr; j++)
            {
                auto atomnr = mymol.topology_->atoms.atom[j].atomnumber;
                auto qq     = mymol.topology_->atoms.atom[j].q;
                qtot       += qq;
                if (mymol.topology_->atoms.atom[j].ptype != eptShell)
                {
                    if (((qq < 0) && (atomnr == 1)) ||
                        ((qq > 0) && ((atomnr == 8)  || (atomnr == 9) ||
                                      (atomnr == 16) || (atomnr == 17) ||
                                      (atomnr == 35) || (atomnr == 53))))
                    {
                        _ener[ermsBOUNDS] += fabs(qq);
                    }                    
                }
            }
            if (fabs(qtot - mymol.molProp()->getCharge()) > 1e-2)
            {
                fprintf(stderr, "Warning qtot for %s is %g, should be %d\n",
                        mymol.molProp()->getMolname().c_str(),
                        qtot, mymol.molProp()->getCharge());
            }
            if (_bQM)
            {
                rvec dmu;
                rvec_sub(mymol.mu_calc, mymol.mu_exp, dmu);
                _ener[ermsMU]  += iprod(dmu, dmu);
                for (int mm = 0; mm < DIM; mm++)
                {
                    if (bfullTensor_)
                    {
                        for (int nn = 0; nn < DIM; nn++)
                        {
                            _ener[ermsQUAD] += gmx::square(mymol.Q_calc[mm][nn] - mymol.Q_exp[mm][nn]);
                        }
                    }
                    else
                    {
                        _ener[ermsQUAD] += gmx::square(mymol.Q_calc[mm][mm] - mymol.Q_exp[mm][mm]);
                    }
                }
            }
            else
            {
                _ener[ermsMU]  += gmx::square(mymol.dip_calc - mymol.dip_exp);
            }
        }
    }
    
    if (PAR(_cr) && !_bFinal)
    {
        gmx_sum(ermsNR, _ener, _cr);
    }    
    if (MASTER(_cr))
    {
        for (j = 0; j < ermsTOT; j++)
        {
            _ener[ermsTOT] += ((_fc[j]*_ener[j])/_nmol_support);
        }
    }   
    if (nullptr != debug && MASTER(_cr))
    {
        fprintf(debug, "ENER:");
        for (j = 0; j < ermsNR; j++)
        {
            fprintf(debug, "  %8.3f", _ener[j]);
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
        auto J00  = pd_.getJ00(_iChargeDistributionModel, ai->name());
        param_.push_back(std::move(J00));

        if (ai->name().compare(_fixchi) != 0)
        {
            auto Chi0 = pd_.getChi0(_iChargeDistributionModel, ai->name());
            param_.push_back(std::move(Chi0));
        }        
        if (_bFitZeta)
        {
            auto nzeta = pd_.getNzeta(_iChargeDistributionModel, ai->name());
            for (int zz = 0; zz < nzeta; zz++)
            {
                auto zeta = pd_.getZeta(_iChargeDistributionModel, ai->name(), zz);
                if (0 != zeta)
                {
                    param_.push_back(std::move(zeta));
                }
                else
                {
                    gmx_fatal(FARGS, "Zeta is zero for atom %s in %d model\n",
                              ai->name().c_str(), _iChargeDistributionModel);
                }
            }
        }       
    }
    if (_bOptHfac)
    {
        param_.push_back(std::move(_hfac));
    }       
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
        std::string qstr   = pd_.getQstr(_iChargeDistributionModel, ai->name());
        std::string rowstr = pd_.getRowstr(_iChargeDistributionModel, ai->name());
        
        if (qstr.size() == 0 || rowstr.size() == 0)
        {
            gmx_fatal(FARGS, "No qstr/rowstr for atom %s in %d model\n",
                      ai->name().c_str(), _iChargeDistributionModel);
        }
        
        auto ei = pd_.findEem(_iChargeDistributionModel, ai->name());
        GMX_RELEASE_ASSERT(ei != pd_.EndEemprops(), "Cannot find eemprops");
        
        ei->setJ0(param_[n++]);
        
        if (ai->name().compare(_fixchi) != 0)
        {      
            ei->setChi0(param_[n++]);           
        }       
        if (_bFitZeta)
        {
            zstr[0] = '\0';
            auto nzeta = pd_.getNzeta(_iChargeDistributionModel, ai->name());
            for (int zz = 0; zz < nzeta; zz++)
            {
                auto zeta = param_[n++];
                sprintf(buf, "  %10g", zeta);
                strcat(zstr, buf);
            }
            
            ei->setRowZetaQ(rowstr, zstr, qstr);
        }               
    }
    if (_bOptHfac)
    {
        _hfac = param_[n++];
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

double OPtimization::objFunction(const double v[])
{    
    double bounds = 0;
    int    n      = 0;
    
    auto *iCount = indexCount();
    for (auto ai = iCount->beginIndex(); ai < iCount->endIndex(); ++ai)
    {
        auto name = ai->name();
        auto J00  = v[n++];
        bounds   += harmonic(J00, _J0_0, _J0_1);
        
        if (strcasecmp(name.c_str(), _fixchi) != 0)
        {
            auto Chi0 = v[n++];
            bounds   += harmonic(Chi0, _Chi0_0, _Chi0_1);
        }
        
        if (_bFitZeta)
        {
            auto nzeta = pd_.getNzeta(_iChargeDistributionModel, ai->name());
            for (int zz = 0; zz < nzeta; zz++)
            {
                auto zeta = v[n++];
                bounds += harmonic(zeta, _w_0, _w_1);
            }
        }
    }
    
    if (_bOptHfac)
    {
        _hfac = v[n++];
        if (_hfac > _hfac0)
        {
            bounds += 100*gmx::square(_hfac - _hfac0);
        }
        else if (_hfac < -(_hfac0))
        {
            bounds += 100*gmx::square(_hfac + _hfac0);
        }
    }
    
    size_t np  = param_.size();
    for (size_t i = 0; i < np; i++)
    {
        param_[i] = v[i];
    }

    tuneDip2PolData();
    calcDeviation();

    _ener[ermsBOUNDS] += bounds;
    _ener[ermsTOT]    += bounds;

    return _ener[ermsTOT];
}

void OPtimization::optRun(FILE *fp, FILE *fplog, int maxiter,
                          int nrun, real stepsize, int seed,
                          const gmx_output_env_t *oenv,
                          int nprint, const char *xvgconv, 
                          const char *xvgepot, real temperature, 
                          bool bBound)
{
    std::vector<double> optx, opts, optm;
    double              chi2, chi2_min;
    gmx_bool            bMinimum = false;
    
    auto func = [&] (const double v[]) {
        return objFunction(v);
    };
    
    if (MASTER(_cr))
    {    
        if (PAR(_cr))
        {
            for (int dest = 1; dest < _cr->nnodes; dest++)
            {
                gmx_send_int(_cr, dest, (nrun*maxiter*param_.size()));
            }
        }
        
        chi2 = chi2_min  = GMX_REAL_MAX;
        Bayes <double> TuneDip(func, param_, lower_, upper_, &chi2);
        TuneDip.Init(xvgconv, xvgepot, oenv, seed, stepsize, 
                     maxiter, nprint,temperature, bBound);
                     
        for (int n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            
            TuneDip.simulate();
            TuneDip.getBestParam(optx);
            TuneDip.getPsigma(opts);
            TuneDip.getPmean(optm);

            if (chi2 < chi2_min)
            {
                bMinimum = true;
                for (size_t k = 0; k < param_.size(); k++)
                {
                    best_[k]   = optx[k];
                    psigma_[k] = opts[k];
                    pmean_[k]  = optm[k];
                }
                chi2_min = chi2;
            }
            TuneDip.setParam(best_);
        }
        if (bMinimum)
        {
            param_  = best_;
            double emin = objFunction(best_.data());
            if (fplog)
            {
                fprintf(fplog, "\nMinimum rmsd value during optimization: %.3f.\n", sqrt(emin));
                fprintf(fplog, "Average and standard deviation of parameters\n");
                for (size_t k = 0; k < param_.size(); k++)
                {
                    fprintf(fplog, "%5zu  %10g  %10g\n", k, pmean_[k], psigma_[k]);
                }
            }
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        int niter = gmx_recv_int(_cr, 0);
        for (int n = 0; n < niter + 2; n++)
        {
            calcDeviation();
        }
    }
    
    _bFinal = true;
    if(MASTER(_cr))
    {
        chi2 = objFunction(best_.data());;
        if (nullptr != fp)
        {
            fprintf(fp, "rmsd: %4.3f  ermsBOUNDS: %4.3f  after %d run(s)\n",
                    sqrt(chi2), _ener[ermsBOUNDS], nrun);
        }
        if (nullptr != fplog)
        {
            fprintf(fplog, "rmsd: %4.3f   ermsBOUNDS: %4.3f  after %d run(s)\n",
                    sqrt(chi2), _ener[ermsBOUNDS], nrun);
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
    real   delta;
    if (nullptr != calc_name)
    {
        m_sub(mol->Q_exp, mol->Q_calc, dQ);
        delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                     gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
        fprintf(fp,
                "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n",
                calc_name,
                mol->Q_calc[XX][XX], mol->Q_calc[XX][YY], mol->Q_calc[XX][ZZ],
                dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                "", mol->Q_calc[YY][YY], mol->Q_calc[YY][ZZ], "", dQ[YY][YY], dQ[YY][ZZ]);
    }
    else
    {
        fprintf(fp, "Quadrupole analysis (5 independent components only)\n");
        fprintf(fp,
                "Exp  (%6.2f %6.2f %6.2f)\n"
                "     (%6s %6.2f %6.2f)\n",
                mol->Q_exp[XX][XX], mol->Q_exp[XX][YY], mol->Q_exp[XX][ZZ],
                "", mol->Q_exp[YY][YY], mol->Q_exp[YY][ZZ]);
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
        rvec_sub(mol->mu_exp, mol->mu_calc, dmu);
        ndmu = norm(dmu);
        cosa = cos_angle(mol->mu_exp, mol->mu_calc);
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
                calc_name, mol->mu_calc[XX], mol->mu_calc[YY], mol->mu_calc[ZZ], 
                norm(mol->mu_calc), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
    }
    else
    {
        fprintf(fp, "Dipole analysis\n");
        fprintf(fp, "Exp  (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f\n",
                mol->mu_exp[XX], mol->mu_exp[YY], mol->mu_exp[ZZ], norm(mol->mu_exp));
    }
}

void OPtimization::print_molecules(FILE *fp, const char *xvgfn, const char *qhisto,
                                   const char *cdiff, const char *mudiff, const char *Qdiff,
                                   const char *espdiff, real dip_toler, real quad_toler, 
                                   real q_toler, const gmx_output_env_t * oenv)
{
    FILE         *xvgf, *qdiff, *mud, *tdiff, *hh, *espd;
    double        d2 = 0;
    real          rms, sigma, aver, error, qq, chi2;
    int           j, n, nout, mm, nn;
     
    struct AtomTypeLsq {
        std::string atomtype;
        gmx_stats_t lsq;
    };
    
    enum {
        eprEEM, eprESP, eprNR
    };
    
    gmx_stats_t               lsq_q, lsq_mu[eprNR], lsq_quad[eprNR], lsq_esp;
    std::vector<AtomTypeLsq>  lsqt;
    const char               *eprnm[eprNR] = { "EEM", "ESP" };
    

    xvgf  = xvgropen(xvgfn, "Correlation between dipoles",
                     "Experimental", "Predicted", oenv);
    xvgr_symbolize(xvgf, 2, eprnm, oenv);
    lsq_q       = gmx_stats_init();
    lsq_quad[0] = gmx_stats_init();
    lsq_quad[1] = gmx_stats_init();
    lsq_mu[0]   = gmx_stats_init();
    lsq_mu[1]   = gmx_stats_init();
    lsq_esp     = gmx_stats_init();
    n           = 0;
    
    for (auto &mol: _mymol)
    {
        if (mol.eSupp != eSupportNo)
        {
            fprintf(fp, "Molecule %d: %s. Qtot: %d, Multiplicity %d\n", n+1,
                    mol.molProp()->getMolname().c_str(),
                    mol.molProp()->getCharge(),
                    mol.molProp()->getMultiplicity());

            print_dipole(fp, &mol, nullptr, dip_toler);
            print_dipole(fp, &mol, (char *)"EEM", dip_toler);
            print_dipole(fp, &mol, (char *)"ESP", dip_toler);

            print_quadrapole(fp, &mol, nullptr, quad_toler);
            print_quadrapole(fp, &mol, (char *)"EEM", quad_toler);
            print_quadrapole(fp, &mol, (char *)"ESP", quad_toler);
            
            chi2 = mol.espRms();
            fprintf(fp, "ESP chi2 %g Hartree/e\n", chi2);           
            fprintf(xvgf, "%10g  %10g\n", mol.dip_exp, mol.dip_calc);
            
            for (mm = 0; mm < DIM; mm++)
            {
                gmx_stats_add_point(lsq_mu[0], mol.mu_exp[mm], mol.mu_calc[mm], 0, 0);
                if (0)
                {
                    for (nn = mm; nn < DIM; nn++)
                    {
                        if (mm < ZZ)
                        {
                            gmx_stats_add_point(lsq_quad[0], mol.Q_exp[mm][nn], mol.Q_calc[mm][nn], 0, 0);
                        }
                    }
                }
                else
                {
                    /* Ignore off-diagonal components */
                    gmx_stats_add_point(lsq_quad[0], mol.Q_exp[mm][mm], mol.Q_calc[mm][mm], 0, 0);

                }
            }

            d2 += gmx::square(mol.dip_exp - mol.dip_calc);
            fprintf(fp, "Atom   Type      q_EEM     q_ESP       x       y       z\n");
            for (j = 0; j < mol.topology_->atoms.nr; j++)
            {
                const char *at = *(mol.topology_->atoms.atomtype[j]);
                auto        k  = std::find_if(lsqt.begin(), lsqt.end(),
                                              [at](const AtomTypeLsq &atlsq)
                                              {
                                                  return atlsq.atomtype.compare(at) == 0;
                                              });
                if (k == lsqt.end())
                {
                    lsqt.resize(1+lsqt.size());
                    lsqt[lsqt.size()-1].atomtype.assign(at);
                    lsqt[lsqt.size()-1].lsq = gmx_stats_init();
                    k = lsqt.end() - 1;
                }
                
                qq = mol.topology_->atoms.atom[j].q;
                fprintf(fp, "%-2d%3d  %-5s  %8.4f  %8.4f%8.3f%8.3f%8.3f %s\n",
                        mol.topology_->atoms.atom[j].atomnumber,
                        j+1,
                        *(mol.topology_->atoms.atomtype[j]),
                        qq,mol.qESP[j],
                        mol.x_[j][XX], mol.x_[j][YY], mol.x_[j][ZZ],
                        fabs(qq - mol.qESP[j]) > q_toler ? "ZZZ" : "");
                gmx_stats_add_point(k->lsq, mol.qESP[j], qq, 0, 0);
                gmx_stats_add_point(lsq_q, mol.qESP[j], qq, 0, 0);
            }
            fprintf(fp, "\n");
            n++;
        }
    }
    fclose(xvgf);

    print_stats(fp, (char *)"dipoles", lsq_mu[0], true, (char *)"Elec", (char *)"EEM");
    print_stats(fp, (char *)"quadrupoles", lsq_quad[0], false, (char *)"Elec", (char *)"EEM");
    print_stats(fp, (char *)"charges", lsq_q, false, (char *)"ESP", (char *)"EEM");
    print_stats(fp, (char *)"esp", lsq_esp, false, (char *)"Elec", (char *)"EEM");
    fprintf(fp, "\n");

    print_stats(fp, (char *)"dipoles", lsq_mu[1], true, (char *)"Elec", (char *)"ESP");
    print_stats(fp, (char *)"quadrupoles", lsq_quad[1], false, (char *)"Elec", (char *)"ESP");

    mud = xvgropen(mudiff, "Correlation between Mu Elec and others",
                   "muElec", "mu", oenv);
    xvgr_symbolize(mud, 2, eprnm, oenv);
    print_lsq_set(mud, lsq_mu[0]);
    print_lsq_set(mud, lsq_mu[1]);
    fclose(mud);

    espd = xvgropen(espdiff, "Correlation between Esp Elec and others",
                    "ESP (Hartree/e)", "ESP (Hartree/e)", oenv);
    xvgr_symbolize(espd, 2, eprnm, oenv);
    print_lsq_set(espd, lsq_esp);
    fclose(espd);

    tdiff = xvgropen(Qdiff, "Correlation between Theta Elec and others",
                     "thetaElec", "theta", oenv);
    xvgr_symbolize(tdiff, 2, eprnm, oenv);
    print_lsq_set(tdiff, lsq_quad[0]);
    print_lsq_set(tdiff, lsq_quad[1]);
    fclose(tdiff);
    qdiff = xvgropen(cdiff, "Correlation between ESP and EEM", "qESP", "qEEM", oenv);
    
    std::vector<const char *> atypes;
    for (const auto &k : lsqt)
    {
        atypes.push_back(k.atomtype.c_str());
    }
    
    xvgr_legend(qdiff, atypes.size(), atypes.data(), oenv);
    xvgr_symbolize(qdiff, atypes.size(), atypes.data(), oenv);
    hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
    xvgr_legend(hh, atypes.size(), atypes.data(), oenv);
    fprintf(fp, "\nDeviations of the charges separated per atomtype:\n");
    
    for (auto k = lsqt.begin(); k < lsqt.end(); ++k)
    {
        int   N;
        real *x, *y;

        print_stats(fp, k->atomtype.c_str(), k->lsq, (k == lsqt.begin()), (char *)"ESP", (char *)"EEM");
        print_lsq_set(qdiff, k->lsq);
        if (gmx_stats_get_npoints(k->lsq, &N) == estatsOK)
        {
            N = N/4;
            if (gmx_stats_make_histogram(k->lsq, 0, &N, ehistoY, 0, &x, &y) == estatsOK)
            {
                fprintf(hh, "@type xy\n");
                for (int i = 0; i < N; i++)
                {
                    fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                }
                fprintf(hh, "&\n");
                free(x);
                free(y);
            }
        }
    }
    fclose(qdiff);
    fclose(hh);

    rms = sqrt(d2/n);
    fprintf(fp, "RMSD = %.3f D\n", rms);
    fprintf(fp, "hfac = %g\n", _hfac);
    gmx_stats_get_ase(lsq_mu[0], &aver, &sigma, &error);
    sigma = rms;
    nout  = 0;
    fprintf(fp, "Overview of outliers (> %.3f off)\n", 2*sigma);
    fprintf(fp, "----------------------------------\n");
    fprintf(fp, "%-20s  %12s  %12s  %12s\n",
            "Name", "Predicted", "Experimental", "Mu-Deviation");
            
    for (auto &mol : _mymol)
    {
        rvec dmu;
        rvec_sub(mol.mu_exp, mol.mu_calc, dmu);
        if ((mol.eSupp != eSupportNo) &&
            (mol.dip_exp > sigma) &&
            (norm(dmu) > 2*sigma))
        {
            fprintf(fp, "%-20s  %12.3f  %12.3f  %12.3f\n",
                    mol.molProp()->getMolname().c_str(),
                    mol.dip_calc, mol.dip_exp,
                    mol.dip_calc - mol.dip_exp);
            nout++;
        }
    }
    if (nout)
    {
        printf("There were %d outliers. See at the very bottom of the log file\n",
               nout);
    }
    else
    {
        printf("No outliers! Well done.\n");
    }
    
    do_view(oenv, xvgfn, nullptr);
    gmx_stats_free(lsq_q);
    gmx_stats_free(lsq_quad[0]);
    gmx_stats_free(lsq_quad[1]);
    gmx_stats_free(lsq_mu[0]);
    gmx_stats_free(lsq_mu[1]);
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
        { efDAT, "-f",         "allmols",     ffREAD  },
        { efDAT, "-d",         "gentop",      ffOPTRD },
        { efDAT, "-o",         "tunedip",     ffWRITE },
        { efDAT, "-sel",       "molselect",   ffREAD  },
        { efXVG, "-table",     "table",       ffOPTRD },
        { efLOG, "-g",         "charges",     ffWRITE },
        { efXVG, "-x",         "dipcorr",     ffWRITE },
        { efXVG, "-qhisto",    "q_histo",     ffWRITE },
        { efXVG, "-qdiff",     "q_diff",      ffWRITE },
        { efXVG, "-mudiff",    "mu_diff",     ffWRITE },
        { efXVG, "-thetadiff", "theta_diff",  ffWRITE },
        { efXVG, "-espdiff",   "esp_diff",    ffWRITE },
        { efXVG, "-conv",      "param-conv",  ffWRITE },
        { efXVG, "-epot",      "param-epot",  ffWRITE }
    };
    
    const  int                  NFILE         = asize(fnm);

    static int                  nrun          = 1;
    static int                  nprint        = 10;
    static int                  maxiter       = 100;
    static int                  reinit        = 0;
    static int                  seed          = 0;
    static int                  minimum_data  = 3;
    static real                 tol           = 1e-3;
    static real                 stol          = 1e-6;
    static real                 watoms        = 0;
    static real                 J0_0          = 5;
    static real                 Chi0_0        = 1;
    static real                 w_0           = 5;
    static real                 step          = 0.01;
    static real                 hfac          = 0;
    static real                 rDecrZeta     = -1;
    static real                 J0_1          = 30;
    static real                 Chi0_1        = 30;
    static real                 w_1           = 50;
    static real                 fc_mu         = 1;
    static real                 fc_bound      = 1;
    static real                 fc_quad       = 1;
    static real                 fc_charge     = 0;
    static real                 fc_esp        = 0;
    static real                 th_toler      = 170;
    static real                 ph_toler      = 5;
    static real                 dip_toler     = 0.5;
    static real                 quad_toler    = 5;
    static real                 q_toler       = 0.25;
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
    static gmx_bool             bZero         = true;
    static gmx_bool             bWeighted     = true;   
    static gmx_bool             bCharged      = true;
    static gmx_bool             bGaussianBug  = true;   
    static gmx_bool             bFitZeta      = true;   
    static const char          *cqdist[]      = {nullptr, "AXp", "AXg", "AXs", nullptr};
    static const char          *cqgen[]       = {nullptr, "None", "EEM", "ESP", "RESP", nullptr};
    
    t_pargs                     pa[]         = {
        { "-tol",   FALSE, etREAL, {&tol},
          "Tolerance for convergence in optimization" },
        { "-maxiter", FALSE, etINT, {&maxiter},
          "Max number of iterations for optimization" },
        { "-nprint",  FALSE, etINT, {&nprint},
          "How often to print the parameters during the simulation" },
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-stol",   FALSE, etREAL, {&stol},
          "If reinit is -1 then a reinit will be done as soon as the simplex size is below this treshold." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },
        { "-charged", FALSE, etBOOL, {&bCharged},
          "Use charged molecules in the parameter tuning as well" },          
        { "-fullTensor", FALSE, etBOOL, {&bfullTensor},
          "consider both diagonal and off-diagonal elements of the Q_Calc matrix for optimization" },        
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Model used for charge distribution" },
        { "-qgen",   FALSE, etENUM, {cqgen},
          "Algorithm used for charge generation" },
        { "-fixchi", FALSE, etSTR,  {&fixchi},
          "Electronegativity for this atom type is fixed. Set to FALSE if you want this variable as well, but read the help text above." },
        { "-seed",   FALSE, etINT,  {&seed},
          "Random number seed. If zero, a seed will be generated." },
        { "-j0",    FALSE, etREAL, {&J0_0},
          "Minimum value that J0 (eV) can obtain in fitting" },
        { "-chi0",    FALSE, etREAL, {&Chi0_0},
          "Minimum value that Chi0 (eV) can obtain in fitting" },
        { "-z0",    FALSE, etREAL, {&w_0},
          "Minimum value that inverse radius (1/nm) can obtain in fitting" },
        { "-j1",    FALSE, etREAL, {&J0_1},
          "Maximum value that J0 (eV) can obtain in fitting" },
        { "-chi1",    FALSE, etREAL, {&Chi0_1},
          "Maximum value that Chi0 (eV) can obtain in fitting" },
        { "-z1",    FALSE, etREAL, {&w_1},
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
        { "-min_data",  FALSE, etINT, {&minimum_data},
          "Minimum number of data points in order to be able to optimize the parameters for a given atomtype" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to optimize, e.g. \"H C Br\". The other available atom types in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of atom types to include but keep constant, e.g. \"O N\". These atom types from gentop.dat are left unmodified" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-watoms", FALSE, etREAL, {&watoms},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use zero. For point+smeared charges 1 is recommended (the default)." },
        { "-polar",    FALSE, etBOOL, {&bPolar},
          "Add polarizabilities to the topology and coordinate file" },
        { "-fitzeta", FALSE, etBOOL, {&bFitZeta},
          "Controls whether or not the Gaussian/Slater widths are optimized." },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule" },
        { "-weight", FALSE, etBOOL, {&bWeighted},
          "Perform a weighted fit, by using the errors in the dipoles presented in the input file. This may or may not improve convergence." },
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
    
    alexandria::OPtimization       opt;
    ChargeDistributionModel        iChargeDistributionModel   = name2eemtype(cqdist[0]);
    ChargeGenerationAlgorithm      iChargeGenerationAlgorithm = (ChargeGenerationAlgorithm) get_option(cqgen);
    const char                    *tabfn                      = opt2fn_null("-table", NFILE, fnm);

    if (iChargeDistributionModel == eqdAXs && nullptr == tabfn)
    {
        gmx_fatal(FARGS, "Cannot generate charges with the %s charge model without a potential table. "
                  "Please supply a table file.", getEemtypeName(iChargeDistributionModel));
    }
    
    opt.Init(cr, bQM, bGaussianBug,
             iChargeDistributionModel,
             iChargeGenerationAlgorithm,
             rDecrZeta, J0_0, Chi0_0, w_0, 
             J0_1, Chi0_1, w_1, fc_bound, 
             fc_mu, fc_quad, fc_charge,
             fc_esp, 1, 1, fixchi, bOptHfac, hfac, 
             bPolar, bFitZeta, hwinfo, bfullTensor);
            
    opt.Read(fp ? fp : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             minimum_data, bZero,
             opt_elem, const_elem,
             lot, gms, watoms, true,
             false, false, bPolar, bZPE,
             opt2fn_null("-table", NFILE, fnm));
            
    if (nullptr != fp)
    {
        fprintf(fp, "In the total data set of %zu molecules we have:\n", opt._mymol.size());
    }

    if (MASTER(cr))
    {
        opt.InitOpt(factor);
    }
    
    opt.optRun(MASTER(cr) ? stderr : nullptr, fp,
               maxiter, nrun, step, seed,
               oenv, nprint,
               opt2fn("-conv", NFILE, fnm),
               opt2fn("-epot", NFILE, fnm),
               temperature, bBound);
               
    if (MASTER(cr))
    {
        opt.print_molecules(fp, opt2fn("-x", NFILE, fnm), opt2fn("-qhisto", NFILE, fnm),
                            opt2fn("-qdiff", NFILE, fnm), opt2fn("-mudiff", NFILE, fnm),
                            opt2fn("-thetadiff", NFILE, fnm), opt2fn("-espdiff", NFILE, fnm),
                            dip_toler, quad_toler, q_toler, oenv);
                            
        writePoldata(opt2fn("-o", NFILE, fnm), opt.pd_, bcompress);
        done_filenms(NFILE, fnm);
        gmx_ffclose(fp);
    }
    return 0;
}
