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

#include "moldip.h"

#include <cmath>

#include <vector>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/smalloc.h"

#include "fill_inputrec.h"
#include "getmdlogger.h"
#include "gmx_simple_comm.h"
#include "molprop_xml.h"
#include "poldata_xml.h"


#define STRLEN 256

namespace alexandria
{

void IndexCount::addName(const std::string &name,
                         bool               bConst)
{
    auto ai = std::find_if(atomIndex_.begin(), atomIndex_.end(),
                           [name](const AtomIndex a)
                           {
                               return a.name().compare(name) == 0;
                           });
    if (atomIndex_.end() == ai)
    {
        AtomIndex aaa(name, bConst);
        atomIndex_.push_back(aaa);
    }
    else
    {
        if (ai->isConst() == bConst)
        {
            fprintf(stderr, "Trying to add %s twice\n", name.c_str());
            // ai.increment();
        }
        else
        {
            gmx_fatal(FARGS, "Trying to add atom %s as both constant and optimized",
                      name.c_str());
        }
    }
}

void IndexCount::sumCount(t_commrec *cr)
{
    totCount_.resize(atomIndex_.size(), 0);
    int i = 0;
    for (const auto &ai : atomIndex_)
    {
        totCount_[i++] = ai.count();
    }
    if (cr->nnodes > 1)
    {
        gmx_sumi(totCount_.size(), totCount_.data(), cr);
    }
}

void IndexCount::incrementName(const std::string &name)
{
    auto ai = findName(name);
    if (ai == atomIndex_.end())
    {
        gmx_fatal(FARGS, "No such atom %s", name.c_str());
    }
    ai->increment();
}

void IndexCount::decrementName(const std::string &name)
{
    auto ai = findName(name);
    if (ai == atomIndex_.end())
    {
        gmx_fatal(FARGS, "No such atom %s", name.c_str());
    }
    ai->decrement();
}

int IndexCount::count(const std::string &name)
{
    auto ai = findName(name);
    if (ai == atomIndex_.end())
    {
        return ai->count();
    }
    return 0;
}

static void dump_index_count(const IndexCount       *ic,
                             FILE                   *fp,
                             ChargeDistributionModel iDistributionModel,
                             const Poldata          &pd,
                             gmx_bool                bFitZeta)
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "Atom index for this optimization.\n");
    fprintf(fp, "Name  Number  Action   #Zeta\n");
    for (auto i = ic->beginIndex(); i < ic->endIndex(); ++i)
    {
        int nZeta = pd.getNzeta(iDistributionModel,
                                i->name());
        int nZopt = 0;
        for (int j = 0; (j < nZeta); j++)
        {
            if (pd.getZeta(iDistributionModel,
                           i->name(), j) > 0)
            {
                nZopt++;
            }
        }
        if (i->isConst())
        {
            fprintf(fp, "%-4s  %6d  Constant\n",
                    i->name().c_str(), i->count());
        }
        else
        {
            fprintf(fp, "%-4s  %6d  Optimized %4d%s\n",
                    i->name().c_str(),
                    i->count(), nZopt,
                    bFitZeta ? " optimized" : " constant");
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
}

int IndexCount::cleanIndex(int   minimum_data,
                           FILE *fp)
{
    int nremove = 0;

    for (auto i = atomIndex_.begin(); i < atomIndex_.end(); ++i)
    {
        if (!i->isConst() && (i->count() < minimum_data))
        {
            if (fp)
            {
                fprintf(fp, "Not enough support in data set for optimizing %s\n",
                        i->name().c_str());
            }
            i = atomIndex_.erase(i);
            nremove++;
        }
    }
    return nremove;
}

static void update_index_count_bool(IndexCount             *ic,
                                    const Poldata          &pd,
                                    const char             *string,
                                    gmx_bool                bSet,
                                    gmx_bool                bAllowZero,
                                    ChargeDistributionModel iDistributionModel)
{
    std::vector<std::string> ptr = gmx::splitString(string);
    for (auto &k : ptr)
    {
        if (pd.haveEemSupport(iDistributionModel, k, bAllowZero))
        {
            ic->addName(k, bSet);
        }
    }
}

static void make_index_count(IndexCount               *ic,
                             const Poldata             &pd,
                             char                      *opt_elem,
                             char                      *const_elem,
                             ChargeDistributionModel   iDistributionModel,
                             gmx_bool                  bFitZeta)
{
    if (nullptr != const_elem)
    {
        update_index_count_bool(ic, pd, const_elem, true, true, iDistributionModel);
    }
    if (nullptr != opt_elem)
    {
        update_index_count_bool(ic, pd, opt_elem, false, true, iDistributionModel);
    }
    else
    {
        for (auto eep = pd.BeginEemprops(); eep != pd.EndEemprops(); ++eep)
        {
            auto ai = ic->findName(eep->getName());
            if ((eep->getEqdModel() == iDistributionModel) &&
                (ai != ic->endIndex()) && !ai->isConst() &&
                pd.haveEemSupport(iDistributionModel, eep->getName(), false))
            {
                ic->addName(eep->getName(), false);
            }
        }
    }
    dump_index_count(ic, debug, iDistributionModel, pd, bFitZeta);
}

immStatus MolDip::check_data_sufficiency(alexandria::MyMol    mymol, 
                                         IndexCount          *ic)
{
    immStatus imm = immOK;

    for (int i = 0; i < mymol.topology_->atoms.nr; i++)
    {
        if ((mymol.topology_->atoms.atom[i].atomnumber > 0) && 
            (mymol.topology_->atoms.atom[i].ptype == eptAtom))
        {
            auto ai = ic->findName(*(mymol.topology_->atoms.atomtype[i]));
            if (ic->endIndex() == ai)
            {
                if (debug)
                {
                    fprintf(debug, "Removing %s because of lacking support for atom %s\n",
                            mymol.molProp()->getMolname().c_str(),
                            *(mymol.topology_->atoms.atomtype[i]));
                }
                imm = immInsufficientDATA;
            }
        }
    }    
    if (imm == immOK)
    {
        for (int i = 0; i < mymol.topology_->atoms.nr; i++)
        {
            if ((mymol.topology_->atoms.atom[i].atomnumber > 0) && 
                (mymol.topology_->atoms.atom[i].ptype == eptAtom))
            {
                ic->incrementName(*(mymol.topology_->atoms.atomtype[i]));
            }
        }
    }
    return imm;
}

MolDip::MolDip()
{
    _cr     = nullptr;
    _fixchi = nullptr;
    for (int i = 0; (i < ermsNR); i++)
    {
        _fc[i]   = 0;
        _ener[i] = 0;
    }

    inputrec_ = mdModules_.inputrec();
    fill_inputrec(inputrec_);
}

void MolDip::Init(t_commrec *cr, gmx_bool bQM, gmx_bool bGaussianBug,
                  ChargeDistributionModel iChargeDistributionModel,
                  ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                  real rDecrZeta,
                  real J0_0, real Chi0_0, real w_0,
                  real J0_1, real Chi0_1, real w_1,
                  real fc_bound, real fc_mu, real fc_quad, real fc_charge,
                  real fc_esp, real fc_epot, real fc_force, char *fixchi,
                  gmx_bool bOptHfac, real hfac,
                  gmx_bool bPol, gmx_bool bFitZeta, 
                  gmx_hw_info_t *hwinfo,
                  gmx_bool bfullTensor)
{
    _cr                         = cr;
    _bQM                        = bQM;
    _bDone                      = false;
    _bFinal                     = false;
    _bGaussianBug               = bGaussianBug;
    _bFitZeta                   = bFitZeta;
    _iChargeDistributionModel   = iChargeDistributionModel;
    _iChargeGenerationAlgorithm = iChargeGenerationAlgorithm;
    _decrzeta                   = rDecrZeta;
    _J0_0                       = J0_0;
    _Chi0_0                     = Chi0_0;
    _w_0                        = w_0;
    _J0_1                       = J0_1;
    _Chi0_1                     = Chi0_1;
    _w_1                        = w_1;
    _fc[ermsMU]                 = fc_mu;
    _fc[ermsBOUNDS]             = fc_bound;
    _fc[ermsQUAD]               = fc_quad;
    _fc[ermsCHARGE]             = fc_charge;
    _fc[ermsESP]                = fc_esp;
    _fc[ermsForce2]             = fc_force;
    _fc[ermsEPOT]               = fc_epot;
    _fixchi                     = fixchi;
    _hfac                       = hfac;
    _hfac0                      = hfac;
    _bOptHfac                   = bOptHfac;
    _bPol                       = bPol;
    hwinfo_                     = hwinfo;
    bfullTensor_                = bfullTensor;
}

void MolDip::Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  char            *opt_elem,
                  char            *const_elem,
                  char            *lot,
                  const MolSelect &gms,
                  real             watoms,
                  gmx_bool         bCheckSupport,
                  bool             bPairs,
                  bool             bDihedral,
                  bool             bPolar,
                  bool             bZPE,
                  const char      *tabfn)
{
    int                              nwarn = 0, nmol_cpu;
    int                              imm_count[immNR];
    immStatus                        imm;
    std::vector<alexandria::MolProp> mp;

    for (int i = 0; (i < immNR); i++)
    {
        imm_count[i] = 0;
    }

    /* Read the EEM parameters */
    _atomprop   = gmx_atomprop_init();

    /* Force field data */
    if (MASTER(_cr))
    {
        try
        {
            alexandria::readPoldata(pd_fn, pd_, _atomprop);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    
    if (PAR(_cr))
    {
        pd_.broadcast(_cr);
    }
    
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s:\n---\n",
                static_cast<int>(pd_.getNatypes()), pd_fn);
        fprintf(fp, "---\n\n");
    }

    /* Read other stuff */
    if (MASTER(_cr))
    {
        /* Now read the molecules */
        MolPropRead(fn, mp);
        for (auto mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            if (false == mpi->GenerateComposition(pd_))
            {
                mpi = mp.erase(mpi);
            }
            mpi->CheckConsistency();
        }
        nmol_cpu = mp.size()/_cr->nnodes + 1;
    }
    else
    {
        nmol_cpu = 0;
    }
    if (PAR(_cr))
    {
        gmx_sumi(1, &nmol_cpu, _cr);
    }    
    if (bCheckSupport && MASTER(_cr))
    {                     
        make_index_count(&indexCount_, pd_,
                         opt_elem, const_elem, 
                         _iChargeDistributionModel, 
                         _bFitZeta);
    }

    int ntopol = 0;
    if (MASTER(_cr))
    {
        for (auto mpi = mp.begin(); (mpi < mp.end()); ++mpi)
        {
            if (imsTrain == gms.status(mpi->getIupac()))
            {
                int               dest = (ntopol % _cr->nnodes);
                alexandria::MyMol mpnew;
                printf("%s\n", mpi->getMolname().c_str());
                mpnew.molProp()->Merge(mpi);

                mpnew.setInputrec(inputrec_);

                imm = mpnew.GenerateTopology(_atomprop, pd_, lot,
                                             _iChargeDistributionModel,
                                             false, bPairs, bDihedral, 
                                             bPolar, tabfn);
                                             
                if (bCheckSupport && immOK == imm)
                {
                    imm = check_data_sufficiency(mpnew, &indexCount_);
                }

                if (immOK == imm)
                {
                    gmx::MDLogger mdlog = getMdLogger(_cr, stdout);
                    imm = mpnew.GenerateCharges(pd_, mdlog, _atomprop,
                                                _iChargeDistributionModel,
                                                _iChargeGenerationAlgorithm,
                                                watoms, _hfac, lot, true,
                                                nullptr, _cr, tabfn, hwinfo_);
                    (void) mpnew.espRms();
                }
                if (immOK == imm)
                {
                    imm = mpnew.GenerateChargeGroups(ecgGroup, false);
                }
                if (immOK == imm)
                {
                    imm = mpnew.getExpProps(_bQM, bZero, bZPE, lot, pd_);
                }

                if (immOK == imm)
                {
                    if (dest > 0)
                    {
                        mpnew.eSupp_ = eSupportRemote;
                        /* Send another molecule */
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Going to send %s to cpu %d\n",
                                    mpi->getMolname().c_str(), dest);
                        }
                        gmx_send_int(_cr, dest, 1);
                        CommunicationStatus cs = mpi->Send(_cr, dest);
                        if (CS_OK != cs)
                        {
                            imm = immCommProblem;
                        }
                        else
                        {
                            imm = (immStatus)gmx_recv_int(_cr, dest);
                        }
                        if (imm != immOK)
                        {
                            fprintf(stderr, "Molecule %s was not accepted on node %d - error %s\n",
                                    mpnew.molProp()->getMolname().c_str(), dest, alexandria::immsg(imm));
                        }
                        else if (nullptr != debug)
                        {
                            fprintf(debug, "Succesfully beamed over %s\n", mpi->getMolname().c_str());
                        }

                    }
                    else
                    {
                        mpnew.eSupp_ = eSupportLocal;
                    }
                    if (immOK == imm)
                    {
                        _mymol.push_back(std::move(mpnew));
                        ntopol++;
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Added %s, ntopol = %d\n",
                                    mpnew.molProp()->getMolname().c_str(),
                                    ntopol);
                        }
                    }
                }
                if ((immOK != imm) && (nullptr != debug))
                {
                    fprintf(debug, "IMM: Dest: %d %s - %s\n",
                            dest, mpi->getMolname().c_str(), immsg(imm));
                }
            }
            else
            {
                imm = immTest;
            }
            imm_count[imm]++;
        }
        /* Send signal done with transferring molecules */
        for (int i = 1; (i < _cr->nnodes); i++)
        {
            gmx_send_int(_cr, i, 0);
        }
    }
    else
    {
        /***********************************************
         *
         *           S L A V E   N O D E S
         *
         ***********************************************/
        ntopol = 0;
        while (gmx_recv_int(_cr, 0) == 1)
        {
            /* Receive another molecule */
            alexandria::MyMol mpnew;

            if (nullptr != debug)
            {
                fprintf(debug, "Going to retrieve new molecule\n");
            }
            CommunicationStatus cs = mpnew.molProp()->Receive(_cr, 0);
            if (CS_OK != cs)
            {
                imm = immCommProblem;
            }
            else if (nullptr != debug)
            {
                fprintf(debug, "Succesfully retrieved %s\n", mpnew.molProp()->getMolname().c_str());
                fflush(debug);
            }

            mpnew.setInputrec(inputrec_);

            imm = mpnew.GenerateTopology(_atomprop, pd_, lot, 
                                         _iChargeDistributionModel,
                                         false, false, bDihedral, 
                                         bPolar, tabfn);

            if (immOK == imm)
            {
                gmx::MDLogger mdlog = getMdLogger(_cr, stdout);
                imm = mpnew.GenerateCharges(pd_, mdlog, _atomprop, 
                                            _iChargeDistributionModel,
                                            _iChargeGenerationAlgorithm, 
                                            watoms, _hfac, lot, true, nullptr, 
                                            _cr, tabfn, hwinfo_);
                (void) mpnew.espRms();
            }
            if (immOK == imm)
            {
                imm = mpnew.GenerateChargeGroups(ecgAtom, false);
            }
            if (immOK == imm)
            {
                imm = mpnew.getExpProps(_bQM, bZero, bZPE, lot, pd_);
            }

            mpnew.eSupp_ = eSupportLocal;
            imm_count[imm]++;
            if (immOK == imm)
            {
                _mymol.push_back(std::move(mpnew));
                if (nullptr != debug)
                {
                    fprintf(debug, "Added molecule %s\n", mpnew.molProp()->getMolname().c_str());
                }
            }
            gmx_send_int(_cr, 0, imm);
        }
    }
    int              nnn     = nmol_cpu;
    std::vector<int> nmolpar;

    if (PAR(_cr))
    {
        nmolpar.resize(_cr->nnodes, 0);
        nmolpar[_cr->nodeid] = nnn;
        gmx_sumi(_cr->nnodes, nmolpar.data(), _cr);
    }

    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        int nmoltot = 0;
        if (PAR(_cr))
        {
            for (int i = 0; (i < _cr->nnodes); i++)
            {
                fprintf(fp, "node %d has %d molecules\n", i, nmolpar[i]);
                nmoltot += nmolpar[i];
            }
        }
        else
        {
            nmoltot = mp.size();
        }
        fprintf(fp, "Made topologies for %d out of %d molecules.\n",
                ntopol, (MASTER(_cr)) ? nmoltot : nmol_cpu);

        for (int i = 0; (i < immNR); i++)
        {
            if (imm_count[i] > 0)
            {
                fprintf(fp, "%d molecules - %s.\n", imm_count[i], alexandria::immsg((immStatus)i));
            }
        }
        if (imm_count[immOK] != (int)mp.size())
        {
            fprintf(fp, "Check alexandria.debug for more information.\nYou may have to use the -debug 1 flag.\n\n");
        }
    }
    _nmol_support = _mymol.size();
    if (_nmol_support == 0)
    {
        gmx_fatal(FARGS, "No support for any molecule!");
    }
}
}
