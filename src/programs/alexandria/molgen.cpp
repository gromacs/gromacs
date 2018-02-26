/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/smalloc.h"

#include "fill_inputrec.h"
#include "getmdlogger.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

#define STRLEN 256

const char *rmsName(int e)
{
    static const char *rms[ermsNR] =
    {
        "BOUNDS", "MU", "QUAD", "CHARGE", "ESP",
        "EPOT", "Force2", "Polar", "TOT"
    };
    if (e >= 0 && e < ermsNR)
    {
        return rms[e];
    }
    else
    {
        return "Incorrect index in rmsName";
    }
}

namespace alexandria
{

static const char *cqdist[] = {nullptr, "AXp", "AXg", "AXs", "AXpp", "AXpg", "AXps", nullptr};
static const char *cqgen[]  = {nullptr, "None", "EEM", "ESP", "RESP", nullptr};

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
        int nZeta  = pd.getNzeta(iDistributionModel, i->name());
        int nZopt  = 0;
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

static void make_index_count(IndexCount                *ic,
                             const Poldata             &pd,
                             char                      *opt_elem,
                             char                      *const_elem,
                             ChargeDistributionModel    iDistributionModel,
                             gmx_bool                   bFitZeta)
{
    if (nullptr != const_elem)
    {
        update_index_count_bool(ic, pd, const_elem, true, true, iDistributionModel);
    }
    if (nullptr != opt_elem)
    {
        update_index_count_bool(ic, pd, opt_elem, false, false, iDistributionModel);
    }
    else
    {
        for (auto eep = pd.BeginEemprops(); eep != pd.EndEemprops(); ++eep)
        {
            if ((eep->getEqdModel() == iDistributionModel) &&
                pd.haveEemSupport(iDistributionModel, eep->getName(), false))
            {
                ic->addName(eep->getName(), false);
            }
        }
    }
    dump_index_count(ic, debug, iDistributionModel, pd, bFitZeta);
}

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
            gmx_fatal(FARGS, "Trying to add atom %s as both constant and optimized",
                      name.c_str());
        }
        else
        {
            fprintf(stderr, "Trying to add %s twice\n", name.c_str());
            // ai.increment();
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

bool IndexCount::isOptimized(const std::string &name)
{
    bool isObtimized = false;
    auto ai          = findName(name);
    if (ai != atomIndex_.end())
    {
        if (!ai->isConst())
        {
            isObtimized = true;
        }
    }
    return isObtimized;
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

int IndexCount::cleanIndex(int   minimum_data,
                           FILE *fp)
{
    int nremove = 0;

    for (auto ai = atomIndex_.begin(); ai < atomIndex_.end(); )
    {
        if (!ai->isConst() && (ai->count() < minimum_data))
        {
            if (fp)
            {
                fprintf(fp, "Not enough support in data set for optimizing %s\n",
                        ai->name().c_str());
            }
            ai = atomIndex_.erase(ai);
            nremove++;
        }
        else
        {
            ++ai;
        }
    }
    return nremove;
}

MolGen::MolGen()
{
    cr_        = nullptr;
    bFinal_    = false;
    bDone_     = false;
    bGenVsite_ = false;
    bOptHfac_  = false;
    qsymm_     = false;
    J0_min_    = 5;
    Chi0_min_  = 1;
    zeta_min_  = 2;
    J0_max_    = 30;
    Chi0_max_  = 30;
    zeta_max_  = 30;
    watoms_    = 0;
    qtol_      = 1e-6;
    qcycle_    = 1000;
    mindata_   = 3;
    nexcl_     = 2;
    hfac_      = 0;
    maxESP_    = 100;
    fixchi_    = (char *)"";
    lot_       = "B3LYP/aug-cc-pVTZ";
    inputrec_  = mdModules_.inputrec();
    fill_inputrec(inputrec_);
    for (int i = 0; i < ermsNR; i++)
    {
        fc_[i]   = 0;
        ener_[i] = 0;
    }
    fc_[ermsTOT] = 1;
}

MolGen::~MolGen()
{
    if (cr_)
    {
        done_commrec(cr_);
    }
}

void MolGen::addOptions(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] =
    {
        { "-mindata", FALSE, etINT, {&mindata_},
          "Minimum number of data points to optimize force field parameters" },
        { "-maxpot", FALSE, etINT, {&maxESP_},
          "Maximum percent of the electrostatic potential points that will be used to fit partial charges." },          
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Model used for charge distribution" },
        { "-qgen",   FALSE, etENUM, {cqgen},
          "Algorithm used for charge generation" },
        { "-lot",    FALSE, etSTR,  {&lot_},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },
        { "-watoms", FALSE, etREAL, {&watoms_},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use zero. For point+smeared charges 1 is recommended." },
        { "-fixchi", FALSE, etSTR,  {&fixchi_},
          "Electronegativity for this atom type is fixed. Set to FALSE if you want this variable as well, but read the help text above (or somewhere)." },
        { "-j0",    FALSE, etREAL, {&J0_min_},
          "Minimum value that J0 (eV) can obtain in fitting" },
        { "-chi0",    FALSE, etREAL, {&Chi0_min_},
          "Minimum value that Chi0 (eV) can obtain in fitting" },
        { "-z0",    FALSE, etREAL, {&zeta_min_},
          "Minimum value that inverse radius (1/nm) can obtain in fitting" },
        { "-j1",    FALSE, etREAL, {&J0_max_},
          "Maximum value that J0 (eV) can obtain in fitting" },
        { "-chi1",    FALSE, etREAL, {&Chi0_max_},
          "Maximum value that Chi0 (eV) can obtain in fitting" },
        { "-z1",    FALSE, etREAL, {&zeta_max_},
          "Maximum value that inverse radius (1/nm) can obtain in fitting" },
        { "-fc_bound",    FALSE, etREAL, {&fc_[ermsBOUNDS]},
          "Force constant in the penalty function for going outside the borders given with the above six options." },
        { "-fc_mu",    FALSE, etREAL, {&fc_[ermsMU]},
          "Force constant in the penalty function for the magnitude of the dipole components." },
        { "-fc_quad",  FALSE, etREAL, {&fc_[ermsQUAD]},
          "Force constant in the penalty function for the magnitude of the quadrupole components." },
        { "-fc_esp",   FALSE, etREAL, {&fc_[ermsESP]},
          "Force constant in the penalty function for the magnitude of the electrostatic potential." },
        { "-fc_charge",  FALSE, etREAL, {&fc_[ermsCHARGE]},
          "Force constant in the penalty function for the magnitude of the charges with respect to the ESP charges." },
        { "-fc_epot",  FALSE, etREAL, {&fc_[ermsEPOT]},
          "Force constant in the penalty function for the magnitude of the potential energy." },
        { "-fc_force",  FALSE, etREAL, {&fc_[ermsForce2]},
          "Force constant in the penalty function for the magnitude of the force." },
        { "-fc_polar",  FALSE, etREAL, {&fc_[ermsPolar]},
          "Force constant in the penalty function for polarizability." },
        { "-qtol",   FALSE, etREAL, {&qtol_},
          "Tolerance for assigning charge generation algorithm." },
        { "-qcycle", FALSE, etINT, {&qcycle_},
          "Max number of tries for optimizing the charges." },
        { "-hfac",  FALSE, etREAL, {&hfac_},
          "[HIDDEN]Fudge factor to scale the J00 of hydrogen by (1 + hfac * qH). Default hfac is 0, means no fudging." },
        { "-nexcl",  FALSE, etINT, {&nexcl_},
          "[HIDDEN]Exclusion number." },
        { "-opthfac",  FALSE, etBOOL, {&bOptHfac_},
          "[HIDDEN]Optimize the fudge factor to scale the J00 of hydrogen (see above). If set, then [TT]-hfac[tt] set the absolute value of the largest hfac. Above this, a penalty is incurred." },
        { "-qm",     FALSE, etBOOL, {&bQM_},
          "Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" },
        { "-qsymm",  FALSE, etBOOL, {&qsymm_},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2." },
        { "-genvsites", FALSE, etBOOL, {&bGenVsite_},
          "Generate virtual sites. Check and double check." }
    };
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void MolGen::optionsFinished()
{
    iChargeDistributionModel_   = name2eemtype(cqdist[0]);
    iChargeGenerationAlgorithm_ = (ChargeGenerationAlgorithm) get_option(cqgen);
    hfac0_                      = hfac_;
    cr_                         = init_commrec();
    mdlog_                      = getMdLogger(cr_, stdout);
    hwinfo_                     = gmx_detect_hardware(mdlog_, cr_, false);
    if (MASTER(cr_))
    {
        printf("There are %d threads/processes.\n", cr_->nnodes);
    }

}

immStatus MolGen::check_data_sufficiency(alexandria::MyMol mymol,
                                         IndexCount       *ic)
{
    immStatus imm = immOK;

    for (int i = 0; i < mymol.topology_->atoms.nr; i++)
    {
        if ((mymol.topology_->atoms.atom[i].atomnumber > 0) &&
            (mymol.topology_->atoms.atom[i].ptype == eptAtom))
        {
            auto fa = pd_.findAtype(*(mymol.topology_->atoms.atomtype[i]));
            if (pd_.getAtypeEnd() != fa)
            {
                const std::string &ztype = fa->getZtype();
                auto               ai    = ic->findName(ztype);
                if (ic->endIndex() == ai)
                {
                    if (debug)
                    {
                        fprintf(debug, "Removing %s because of lacking support for atom %s\n",
                                mymol.molProp()->getMolname().c_str(),
                                ztype.c_str());
                    }
                    imm = immInsufficientDATA;
                }
            }
            else
            {
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
                auto fa = pd_.findAtype(*(mymol.topology_->atoms.atomtype[i]));
                ic->incrementName(fa->getZtype());
            }
        }
    }
    return imm;
}

void MolGen::Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  char            *opt_elem,
                  char            *const_elem,
                  const MolSelect &gms,
                  gmx_bool         bCheckSupport,
                  bool             bPairs,
                  bool             bDihedral,
                  bool             bZPE,
                  bool             bFitZeta,
                  const char      *tabfn)
{
    int                              nwarn    = 0;
    int                              nmol_cpu = 0;
    int                              imm_count[immNR];
    immStatus                        imm      = immOK;
    std::vector<alexandria::MolProp> mp;

    atomprop_  = gmx_atomprop_init();
    for (int i = 0; i < immNR; i++)
    {
        imm_count[i] = 0;
    }
    /*Reading Force Field Data from gentop.dat*/
    if (MASTER(cr_))
    {
        try
        {
            alexandria::readPoldata(pd_fn, pd_, atomprop_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        if (pd_.getNexcl() != nexcl_)
        {
            fprintf(fp, "Exclusion number changed from %d to %d read from the command line.\n", 
                    pd_.getNexcl(), nexcl_);
            pd_.setNexcl(nexcl_);
        }
    }
    /*Broadcasting Force Field Data from Master to Slave nodes*/
    if (PAR(cr_))
    {
        pd_.broadcast(cr_);
    }
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s:\n---\n",
                static_cast<int>(pd_.getNatypes()), pd_fn);
        fprintf(fp, "---\n\n");
    }
    /*Reading Molecules from allmols.dat*/
    if (MASTER(cr_))
    {
        MolPropRead(fn, mp);
        for (auto mpi = mp.begin(); mpi < mp.end(); )
        {
            mpi->CheckConsistency();
            if (false == mpi->GenerateComposition(pd_) || imsTrain != gms.status(mpi->getIupac()))
            {
                mpi = mp.erase(mpi);
            }
            else
            {
                ++mpi;
            }
        }
        nmol_cpu = mp.size()/cr_->nnodes + 1;
    }
    else
    {
        nmol_cpu = 0;
    }
    /*Sort Molecules based on the number of atoms*/
    if (MASTER(cr_))
    {
        std::sort(mp.begin(), mp.end(),
                  [](alexandria::MolProp &mp1,
                     alexandria::MolProp &mp2)
                  {
                      return (mp1.NAtom() < mp2.NAtom());
                  });
    }
    if (PAR(cr_))
    {
        gmx_sumi(1, &nmol_cpu, cr_);
    }
    if (bCheckSupport && MASTER(cr_))
    {
        make_index_count(&indexCount_,
                         pd_,
                         opt_elem,
                         const_elem,
                         iChargeDistributionModel_,
                         bFitZeta);
    }
    /*Generate topology for Molecules and distribute them among the nodes*/
    bool bPolar = (iChargeDistributionModel() == eqdAXpp  ||
                   iChargeDistributionModel() == eqdAXpg  ||
                   iChargeDistributionModel() == eqdAXps);

    int ntopol = 0;
    if (MASTER(cr_))
    {
        for (auto mpi = mp.begin(); mpi < mp.end(); ++mpi)
        {
            if (imsTrain == gms.status(mpi->getIupac()))
            {
                int               dest = (ntopol % cr_->nnodes);
                alexandria::MyMol mymol;
                printf("%s\n", mpi->getMolname().c_str());
                mymol.molProp()->Merge(mpi);
                mymol.setInputrec(inputrec_);
                imm = mymol.GenerateTopology(atomprop_,
                                             pd_,
                                             lot_,
                                             iChargeDistributionModel_,
                                             bGenVsite_,
                                             bPairs,
                                             bDihedral,
                                             bPolar,
                                             false,
                                             tabfn);
                if (bCheckSupport && immOK == imm)
                {
                    imm = check_data_sufficiency(mymol, &indexCount_);
                }
                if (immOK == imm)
                {
                    gmx::MDLogger mdlog = getMdLogger(cr_, stdout);
                    imm = mymol.GenerateCharges(pd_,
                                                mdlog,
                                                atomprop_,
                                                iChargeDistributionModel_,
                                                iChargeGenerationAlgorithm_,
                                                watoms_,
                                                hfac_,
                                                lot_,
                                                qsymm_,
                                                nullptr,
                                                cr_,
                                                tabfn,
                                                hwinfo_,
                                                qcycle_,
                                                maxESP_,
                                                qtol_,
                                                nullptr);
                    (void) mymol.espRms();
                }
                if (immOK == imm)
                {
                    imm = mymol.GenerateChargeGroups(ecgGroup, false);
                }
                if (immOK == imm)
                {
                    imm = mymol.getExpProps(bQM_, bZero, bZPE, lot_, pd_);
                }
                if (immOK == imm)
                {
                    if (dest > 0)
                    {
                        mymol.eSupp_ = eSupportRemote;
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Going to send %s to cpu %d\n",
                                    mpi->getMolname().c_str(), dest);
                        }
                        gmx_send_int(cr_, dest, 1);
                        CommunicationStatus cs = mpi->Send(cr_, dest);
                        if (CS_OK != cs)
                        {
                            imm = immCommProblem;
                        }
                        else
                        {
                            imm = (immStatus)gmx_recv_int(cr_, dest);
                        }
                        if (imm != immOK)
                        {
                            fprintf(stderr, "Molecule %s was not accepted on node %d - error %s\n",
                                    mymol.molProp()->getMolname().c_str(), dest, alexandria::immsg(imm));
                        }
                        else if (nullptr != debug)
                        {
                            fprintf(debug, "Succesfully beamed over %s\n", mpi->getMolname().c_str());
                        }

                    }
                    else
                    {
                        mymol.eSupp_ = eSupportLocal;
                    }
                    if (immOK == imm)
                    {
                        mymol_.push_back(std::move(mymol));
                        ntopol++;
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Added %s, ntopol = %d\n", mymol.molProp()->getMolname().c_str(), ntopol);
                        }
                    }
                }
                if ((immOK != imm) && (nullptr != debug))
                {
                    fprintf(debug, "IMM: Dest: %d %s - %s\n", dest, mpi->getMolname().c_str(), immsg(imm));
                }
            }
            else
            {
                imm = immTest;
            }
            imm_count[imm]++;
        }
        /* Send signal done with transferring molecules */
        for (int i = 1; i < cr_->nnodes; i++)
        {
            gmx_send_int(cr_, i, 0);
        }
    }
    else
    {
        /***********************************************
         *                                             *
         *           S L A V E   N O D E S             *
         *                                             *
         ***********************************************/
        ntopol = 0;
        while (gmx_recv_int(cr_, 0) == 1)
        {
            alexandria::MyMol mymol;
            if (nullptr != debug)
            {
                fprintf(debug, "Going to retrieve new molecule\n");
            }
            CommunicationStatus cs = mymol.molProp()->Receive(cr_, 0);
            if (CS_OK != cs)
            {
                imm = immCommProblem;
            }
            else if (nullptr != debug)
            {
                fprintf(debug, "Succesfully retrieved %s\n", mymol.molProp()->getMolname().c_str());
                fflush(debug);
            }
            mymol.setInputrec(inputrec_);
            imm = mymol.GenerateTopology(atomprop_,
                                         pd_,
                                         lot_,
                                         iChargeDistributionModel_,
                                         bGenVsite_,
                                         bPairs,
                                         bDihedral,
                                         bPolar,
                                         false,
                                         tabfn);

            if (immOK == imm)
            {
                gmx::MDLogger mdlog = getMdLogger(cr_, stdout);
                imm = mymol.GenerateCharges(pd_,
                                            mdlog,
                                            atomprop_,
                                            iChargeDistributionModel_,
                                            iChargeGenerationAlgorithm_,
                                            watoms_,
                                            hfac_,
                                            lot_,
                                            qsymm_,
                                            nullptr,
                                            cr_,
                                            tabfn,
                                            hwinfo_,
                                            qcycle_,
                                            maxESP_,
                                            qtol_,
                                            nullptr);
                (void) mymol.espRms();
            }
            if (immOK == imm)
            {
                imm = mymol.GenerateChargeGroups(ecgAtom, false);
            }
            if (immOK == imm)
            {
                imm = mymol.getExpProps(bQM_, bZero, bZPE, lot_, pd_);
            }
            mymol.eSupp_ = eSupportLocal;
            imm_count[imm]++;
            if (immOK == imm)
            {
                mymol_.push_back(std::move(mymol));
                if (nullptr != debug)
                {
                    fprintf(debug, "Added molecule %s\n", mymol.molProp()->getMolname().c_str());
                }
            }
            gmx_send_int(cr_, 0, imm);
        }
    }
    int              nnn = nmol_cpu;
    std::vector<int> nmolpar;
    if (PAR(cr_))
    {
        nmolpar.resize(cr_->nnodes, 0);
        nmolpar[cr_->nodeid] = nnn;
        gmx_sumi(cr_->nnodes, nmolpar.data(), cr_);
    }
    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        int nmoltot = 0;
        if (PAR(cr_))
        {
            for (int i = 0; i < cr_->nnodes; i++)
            {
                fprintf(fp, "Node %d has %d molecules\n", i, nmolpar[i]);
                nmoltot += nmolpar[i];
            }
        }
        else
        {
            nmoltot = mp.size();
        }
        fprintf(fp, "Made topologies for %d out of %d molecules.\n",
                ntopol, (MASTER(cr_)) ? nmoltot : nmol_cpu);

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
    if (bCheckSupport && MASTER(cr_))
    {
        indexCount_.cleanIndex(mindata_, fp);
    }
    nmol_support_ = mymol_.size();
    if (nmol_support_ == 0)
    {
        gmx_fatal(FARGS, "No support for any molecule!");
    }
}
}
