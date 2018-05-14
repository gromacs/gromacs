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
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "communication.h"
#include "gentop_core.h"
#include "getmdlogger.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "mymol.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "regression.h"

/*! \brief Write a csv file containing molecule names and bond energy
 *
 * Writes the whole bond energy matrix.
 */
void dump_csv(const std::vector<std::string>        &ctest,
              const std::vector<alexandria::MyMol>  &mm,
              const std::vector<int>                &ntest,
              const std::vector<double>             &Edissoc,
              const MatrixWrapper                   &a,
              const double                           x[])
{
    FILE *csv = gmx_ffopen("tune_fc.csv", "w");
    fprintf(csv, ",");
    for (auto j : ctest)
    {
        fprintf(csv, "%s,", j.c_str());
    }
    fprintf(csv, "\n");
    int i = 0;
    for (auto &mymol : mm)
    {
        fprintf(csv, "%s,", mymol.molProp()->getMolname().c_str());
        for (size_t j = 0; (j < ctest.size()); j++)
        {
            fprintf(csv, "%g,", a.get(j, i));
        }
        fprintf(csv, "%.3f\n", x[i]);
        i++;
    }
    fprintf(csv, "Total,");
    for (auto j : ntest)
    {
        fprintf(csv, "%d,", j);
    }
    fprintf(csv, "\n");
    fprintf(csv, "Edissoc,");
    for (auto j : Edissoc)
    {
        fprintf(csv, "%.3f,", j);
    }
    fprintf(csv, "\n");
    fclose(csv);
}

namespace alexandria
{

class AtomTypes
{
    public:

        AtomTypes () {}

        AtomTypes(int                ncopies,
                  const std::string &name,
                  const std::string &vdwParams,
                  int                index)
            :
              ncopies_(ncopies),
              name_(name),
              vdwParams_(vdwParams),
              poldataIndex_(index)
        {
            extractParams();
        }

        void inc() { ncopies_++; }

        int nCopies() const { return ncopies_; }

        void setParamString(const std::string &params);

        int poldataIndex() const { return poldataIndex_; }

        const std::string &paramString() const { return vdwParams_; }

        const std::vector<double> &paramValues() const { return p_; }

        const std::string &name() const { return name_; }

        size_t nParams() const { return p_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        int                 ncopies_;
        std::string         name_;
        std::string         vdwParams_;
        std::vector<double> p_;
        int                 poldataIndex_;
        void extractParams();
};

using AtomTypesIterator = typename std::vector<AtomTypes>::iterator;

CommunicationStatus AtomTypes::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ncopies_);
        gmx_send_str(cr, dest, &name_);
        gmx_send_str(cr, dest, &vdwParams_);
        gmx_send_int(cr, dest, poldataIndex_);
        gmx_send_int(cr, dest, p_.size());
        for (auto &p : p_)
        {
            gmx_send_double(cr, dest, p);
        }
    }
    return cs;
}

CommunicationStatus AtomTypes::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        ncopies_      = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &name_);
        gmx_recv_str(cr, src, &vdwParams_);
        poldataIndex_ = gmx_recv_int(cr, src);
        int np        = gmx_recv_int(cr, src);
        for (int n = 0; n < np; n++)
        {
            auto p = gmx_recv_double(cr, src);
            p_.push_back(p);
        }
    }
    return cs;
}

void AtomTypes::setParamString(const std::string &params)
{
    vdwParams_ = params;
    extractParams();
}

void AtomTypes::extractParams()
{
    const auto p = gmx::splitString(vdwParams_);
    p_.clear();
    for (const auto &d : p)
    {
        p_.push_back(std::move(atof(d.c_str())));
    }
}

class NonBondParams
{
    public:

        NonBondParams () {}

        NonBondParams(bool bOpt, InteractionType  itype)

            :
              bOpt_(bOpt),
              itype_(itype)
        {}


        void addNonBonded(AtomTypes at) { at_.push_back(std::move(at)); }

        void analyzeIdef(std::vector<MyMol> &mm,
                         const Poldata      &pd);

        void makeReverseIndex();

        int reverseIndex(int poldataIndex)
        {
            GMX_RELEASE_ASSERT(poldataIndex >= 0 && poldataIndex < static_cast<int>(reverseIndex_.size()), "Incorrect poldataIndex");
            GMX_RELEASE_ASSERT(reverseIndex_[poldataIndex] != -1, "The reverseIndex is incorrect");

            return reverseIndex_[poldataIndex];
        }

        AtomTypesIterator beginAT() { return at_.begin(); }

        AtomTypesIterator endAT() { return at_.end(); }

        InteractionType interactionType() const { return itype_; }

        size_t nAT() const { return at_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        bool                   bOpt_;
        InteractionType        itype_;
        std::vector<AtomTypes> at_;
        std::vector<int>       reverseIndex_;
        std::vector<double>    params_;

};

CommunicationStatus NonBondParams::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    std::string         itype;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, bOpt_);
        itype.assign(iType2string(itype_));
        gmx_send_str(cr, dest, &itype);
        gmx_send_int(cr, dest, at_.size());
        gmx_send_int(cr, dest, reverseIndex_.size());
        gmx_send_int(cr, dest, params_.size());
        for (auto &at : at_)
        {
            at.Send(cr, dest);
        }
        for (auto &ri : reverseIndex_)
        {
            gmx_send_int(cr, dest, ri);
        }
        for (auto &param : params_)
        {
            gmx_send_double(cr, dest, param);
        }
    }
    return cs;
}

CommunicationStatus NonBondParams::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    std::string         itype;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        bOpt_      = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &itype);
        itype_     = string2iType(itype.c_str());
        int nat    = gmx_recv_int(cr, src);
        int nri    = gmx_recv_int(cr, src);
        int nparam = gmx_recv_int(cr, src);

        for (int n = 0; n < nat; n++)
        {
            AtomTypes at;
            cs = at.Receive(cr, src);
            if (CS_OK == cs)
            {
                at_.push_back(at);
            }
        }
        for (int n = 0; n < nri; n++)
        {
            auto ri = gmx_recv_int(cr, src);
            reverseIndex_.push_back(ri);
        }
        for (int n = 0; n < nparam; n++)
        {
            auto param = gmx_recv_double(cr, src);
            params_.push_back(param);
        }
    }
    return cs;
}

void NonBondParams::analyzeIdef(std::vector<MyMol> &mm,
                                const Poldata      &pd)
{
    if (!bOpt_)
    {
        return;
    }
    for (auto &mymol : mm)
    {
        for (int i = 0; i < mymol.topology_->atoms.nr; i++)
        {
            if (mymol.topology_->atoms.atom[i].ptype == eptAtom ||
                mymol.topology_->atoms.atom[i].ptype == eptNucleus)
            {
                std::string params;
                int         index = 0;
                char        buf[STRLEN];
                bool        found = false;
                auto        at    = *mymol.topology_->atoms.atomtype[i];
                auto        fat   = pd.findAtype(at);
                if (fat != pd.getAtypeEnd())
                {
                    sprintf(buf, "%s", at);
                    params    = fat->getVdwparams();
                    index     = fat - pd.getAtypeBegin();
                    found     = true;
                }
                if (found)
                {
                    auto c = std::find_if(at_.begin(), at_.end(),
                                          [buf](const AtomTypes &at)
                                          {
                                              return (at.name().compare(buf) == 0);
                                          });
                    if (c != at_.end())
                    {
                        c->inc();
                    }
                    else
                    {
                        AtomTypes at(1, buf, params, index);
                        addNonBonded(at);
                    }
                }
            }
        }
    }
}

void NonBondParams::makeReverseIndex()
{
    int N = 0;
    for (const auto &i : at_)
    {
        N = std::max(N, i.poldataIndex());
    }
    reverseIndex_.resize(N+1, -1);
    int j = 0;
    for (const auto &i : at_)
    {
        reverseIndex_[i.poldataIndex()] = j++;
    }
}

/*! \brief Helper class storing bond/angle/dihedral names
 *
 * For one bond/angle/dihedral here the name of the bondtypes
 * are stored as in e.g. c c h for an angle, along with the number
 * of occurrences in the force field.
 */
class BondNames
{
    public:

        BondNames () {}

        BondNames(int                ncopies,
                  const std::string &name,
                  const std::string &params,
                  int                index,
                  double             bondorder = 0)

            :
              ncopies_(ncopies),
              name_(name),
              params_(params),
              bondorder_(bondorder),
              poldataIndex_(index)
        {
            extractParams();
        }

        void inc() { ncopies_++; }

        int nCopies() const { return ncopies_; }

        void setParamString(const std::string &params);

        const std::string &name() const { return name_; }

        double bondorder() const { return bondorder_; }

        int poldataIndex() const { return poldataIndex_; }

        const std::string &paramString() const { return params_; }

        const std::vector<double> &paramValues() const { return p_; }

        size_t nParams() const { return p_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        //! Number of copies in the molecule data set
        int                 ncopies_;
        //! Name of this bond/angle/dihedral
        std::string         name_;
        //! String holding all the parameters
        std::string         params_;
        //! Vector containing all the parameters
        std::vector<double> p_;
        //! The bond order in case this is a bond
        double              bondorder_;
        //! Index in Poldata structure
        int                 poldataIndex_;
        //! Internal routine to extract the parameters
        void extractParams();
};

using BondNamesIterator = typename std::vector<BondNames>::iterator;

CommunicationStatus BondNames::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ncopies_);
        gmx_send_str(cr, dest, &name_);
        gmx_send_str(cr, dest, &params_);
        gmx_send_double(cr, dest, bondorder_);
        gmx_send_int(cr, dest, poldataIndex_);
        gmx_send_int(cr, dest, p_.size());
        for (auto &p : p_)
        {
            gmx_send_double(cr, dest, p);
        }
    }
    return cs;
}

CommunicationStatus BondNames::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        ncopies_      = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &name_);
        gmx_recv_str(cr, src, &params_);
        bondorder_    = gmx_recv_double(cr, src);
        poldataIndex_ = gmx_recv_int(cr, src);
        int np        = gmx_recv_int(cr, src);
        for (int n = 0; n < np; n++)
        {
            auto p = gmx_recv_double(cr, src);
            p_.push_back(p);
        }
    }
    return cs;
}

void BondNames::setParamString(const std::string &params)
{
    params_ = params;
    extractParams();
}

void BondNames::extractParams()
{
    const auto p = gmx::splitString(params_);
    p_.clear();
    for (const auto &d : p)
    {
        p_.push_back(std::move(atof(d.c_str())));
    }
}

/*! \brief Class holding for one type of interactions all names
 *
 * Class holding the OptNames for each interaction type.
 */
class ForceConstants
{

    public:

        ForceConstants () {}

        ForceConstants(int bt, int ftype, InteractionType itype, bool bOpt)
            :
              bt_(bt),
              ftype_(ftype),
              itype_(itype),
              bOpt_(bOpt)
        {}

        void addForceConstant(BondNames bn) { bn_.push_back(std::move(bn)); }

        void analyzeIdef(std::vector<MyMol> &mm,
                         const Poldata      &pd);

        /*! \brief Make reverse index from Poldata to BondNames
         *
         * The BondNames structure stores the Poldata index for
         * all interactions. This routine makes an index to convert
         * the Poldata index to the index in BondNames.
         */
        void makeReverseIndex();

        int reverseIndex(int poldataIndex)
        {
            GMX_RELEASE_ASSERT(poldataIndex >= 0 && poldataIndex < static_cast<int>(reverseIndex_.size()), "Incorrect poldataIndex");
            GMX_RELEASE_ASSERT(reverseIndex_[poldataIndex] != -1, "The reverseIndex is incorrect");

            return reverseIndex_[poldataIndex];
        }

        int bt2() const { return bt_; }

        int ftype() const { return ftype_; }

        InteractionType interactionType() const { return itype_; }

        void dump(FILE *fp) const;

        BondNamesIterator beginBN() { return bn_.begin(); }

        BondNamesIterator endBN() { return bn_.end(); }

        size_t nbad() const { return bn_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        int                    bt_;
        int                    ftype_;
        InteractionType        itype_;
        bool                   bOpt_;
        std::vector<BondNames> bn_;
        std::vector<int>       reverseIndex_;
        std::vector<double>    params_;
};

CommunicationStatus ForceConstants::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    std::string         itype;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, bt_);
        gmx_send_int(cr, dest, ftype_);
        itype.assign(iType2string(itype_));
        gmx_send_str(cr, dest, &itype);
        gmx_send_int(cr, dest, bOpt_);
        gmx_send_int(cr, dest, bn_.size());
        gmx_send_int(cr, dest, reverseIndex_.size());
        gmx_send_int(cr, dest, params_.size());

        for (auto &bn : bn_)
        {
            bn.Send(cr, dest);
        }
        for (auto &ri : reverseIndex_)
        {
            gmx_send_int(cr, dest, ri);
        }
        for (auto &param : params_)
        {
            gmx_send_double(cr, dest, param);
        }
    }
    return cs;
}

CommunicationStatus ForceConstants::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    std::string         itype;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        bt_           = gmx_recv_int(cr, src);
        ftype_        = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &itype);
        itype_        = string2iType(itype.c_str());
        bOpt_         = gmx_recv_int(cr, src);
        int nbn       = gmx_recv_int(cr, src);
        int nri       = gmx_recv_int(cr, src);
        int nparam    = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < nbn); n++)
        {
            BondNames bn;
            cs = bn.Receive(cr, src);
            if (CS_OK == cs)
            {
                bn_.push_back(bn);
            }
        }
        for (int n = 0; (CS_OK == cs) && (n < nri); n++)
        {
            auto ri = gmx_recv_int(cr, src);
            reverseIndex_.push_back(ri);
        }
        for (int n = 0; (CS_OK == cs) && (n < nparam); n++)
        {
            auto param = gmx_recv_double(cr, src);
            params_.push_back(param);
        }
    }
    return cs;
}

void ForceConstants::analyzeIdef(std::vector<MyMol> &mm,
                                 const Poldata      &pd)
{
    std::string  aai, aaj, aak, aal;

    if (!bOpt_)
    {
        return;
    }
    for (auto &mymol : mm)
    {
        for (int i = 0; (i < mymol.ltop_->idef.il[ftype_].nr);
             i += interaction_function[ftype_].nratoms+1)
        {
            std::vector<std::string> atoms;
            std::string              params;
            bool                     found     = false;
            int                      ai        = mymol.ltop_->idef.il[ftype_].iatoms[i+1];
            int                      aj        = mymol.ltop_->idef.il[ftype_].iatoms[i+2];
            if (pd.atypeToBtype( *mymol.topology_->atoms.atomtype[ai], aai) &&
                pd.atypeToBtype( *mymol.topology_->atoms.atomtype[aj], aaj))
            {
                int  index = 0;
                char buf[STRLEN];
                char buf_reverse[STRLEN];
                auto iType = static_cast<InteractionType>(bt_);
                switch (iType)
                {
                    case eitBONDS:
                    {
                        atoms   = {aai, aaj};
                        auto fs = pd.findForces(iType);
                        auto f  = fs->findForce(atoms);

                        if (fs->forceEnd() != f)
                        {
                            sprintf(buf, "%s %s", aai.c_str(), aaj.c_str());
                            sprintf(buf_reverse, "%s %s", aaj.c_str(), aai.c_str());
                            params    = f->params();
                            index     = f - fs->forceBegin();
                            found     = true;
                        }
                    }
                    break;
                    case eitANGLES:
                    case eitLINEAR_ANGLES:
                    {
                        int ak  = mymol.ltop_->idef.il[ftype_].iatoms[i+3];
                        if (pd.atypeToBtype( *mymol.topology_->atoms.atomtype[ak], aak))
                        {
                            atoms   = {aai, aaj, aak};
                            auto fs = pd.findForces(iType);
                            auto f  = fs->findForce(atoms);

                            if (fs->forceEnd() != f)
                            {
                                sprintf(buf, "%s %s %s", aai.c_str(),
                                        aaj.c_str(), aak.c_str());
                                sprintf(buf_reverse, "%s %s %s", aak.c_str(),
                                        aaj.c_str(), aai.c_str());
                                params = f->params();
                                index  = f - fs->forceBegin();
                                found  = true;
                            }
                        }
                    }
                    break;
                    case eitPROPER_DIHEDRALS:
                    case eitIMPROPER_DIHEDRALS:
                    {
                        int ak  = mymol.ltop_->idef.il[ftype_].iatoms[i+3];
                        int al  = mymol.ltop_->idef.il[ftype_].iatoms[i+4];
                        if (pd.atypeToBtype( *mymol.topology_->atoms.atomtype[ak], aak) &&
                            pd.atypeToBtype( *mymol.topology_->atoms.atomtype[al], aal))
                        {
                            atoms   = {aai, aaj, aak, aal};
                            auto fs = pd.findForces(iType);
                            auto f  = fs->findForce(atoms);

                            if (fs->forceEnd() != f)
                            {
                                sprintf(buf, "%s %s %s %s", aai.c_str(),
                                        aaj.c_str(), aak.c_str(), aal.c_str());
                                sprintf(buf_reverse, "%s %s %s %s", aal.c_str(),
                                        aaj.c_str(), aak.c_str(), aai.c_str());
                                params = f->params();
                                index  = f - fs->forceBegin();
                                found  = true;
                            }
                        }
                    }
                    break;
                    case eitPOLARIZATION:
                    case eitVDW:
                    case eitLJ14:
                    case eitVSITE2:
                    case eitVSITE3FAD:
                    case eitVSITE3OUT:
                    case eitCONSTR:
                    case eitNR:
                        break;
                }
                if (found)
                {
                    auto c = std::find_if(bn_.begin(), bn_.end(),
                                          [buf, buf_reverse](const BondNames &bn)
                                          {
                                              return (bn.name().compare(buf) == 0 ||
                                                      bn.name().compare(buf_reverse) == 0);
                                          });
                    if (c != bn_.end())
                    {
                        c->inc();
                    }
                    else
                    {
                        BondNames bn(1, buf, params, index, 0);
                        addForceConstant(bn);
                    }
                }
            }
        }
    }
}

void ForceConstants::makeReverseIndex()
{
    int N = 0;
    for (const auto &i : bn_)
    {
        N = std::max(N, i.poldataIndex());
    }
    reverseIndex_.resize(N+1, -1);
    int j = 0;
    for (const auto &i : bn_)
    {
        reverseIndex_[i.poldataIndex()] = j++;
    }
}

void ForceConstants::dump(FILE *fp) const
{
    const char  *btsnames[ebtsNR] =
    { "bond", "angle", "proper", "improper", nullptr, nullptr };

    if (bOpt_)
    {
        int ntot = 0;
        fprintf(fp, "Interaction  Bondtypes             Copies Poldata entry\n");
        for (const auto &i : bn_)
        {
            fprintf(fp, "%-10s  %-20s  %5d  %5d\n",
                    btsnames[bt_], i.name().c_str(), i.nCopies(), i.poldataIndex());
            ntot += i.nCopies();
        }
        fprintf(fp, "%-8s %d of %4d types\n", btsnames[bt_], ntot,
                static_cast<int>(bn_.size()));
    }
}

class Optimization : public MolGen
{
    using param_type = std::vector<double>;

    private:
        std::vector<ForceConstants> ForceConstants_;
        std::vector<NonBondParams>  NonBondParams_;
        param_type                  param_, lower_, upper_, best_;
        param_type                  De_, orig_, psigma_, pmean_;
        std::vector<gmx_bool>       bOpt_;
        real                        factor_;
        real                        beta0_, D0_, beta_min_, D0_min_;
        const char                 *lot_;
        Bayes <double>              TuneFc_;
    public:


        /*! \brief
         *
         * Constructor
         */
        Optimization()
        {
            bOpt_.resize(eitNR, false);
            bOpt_[eitBONDS] = true;
            factor_         = 1;
            beta0_          = 0;
            D0_             = 0;
            beta_min_       = 10;
            D0_min_         = 50;
            lot_            = nullptr;
        };

        /*! \brief
         *
         * Destructor
         */
        ~Optimization() {};

        bool bOpt(int i) const { return bOpt_[i]; }
        void add_pargs(std::vector<t_pargs> *pargs)
        {
            t_pargs pa[] =
            {
                { "-beta0", FALSE, etREAL, {&beta0_},
                  "Reset the initial beta for Morse potentials to this value, independent of gentop.dat. If value is <= 0 gentop.dat value is used." },
                { "-D0", FALSE, etREAL, {&D0_},
                  "Reset the initial D for Morse potentials to this value, independent of gentop.dat. If value is <= 0 gentop.dat value is used." },
                { "-beta_min", FALSE, etREAL, {&beta_min_},
                  "Minimum value for beta in Morse potential" },
                { "-DO_min", FALSE, etREAL, {&D0_min_},
                  "Minimum value for D0 in Morse potential" },
                { "-bonds",   FALSE, etBOOL, {&bOpt_[eitBONDS]},
                  "Optimize bond parameters" },
                { "-angles",  FALSE, etBOOL, {&bOpt_[eitANGLES]},
                  "Optimize angle parameters" },
                { "-langles", FALSE, etBOOL, {&bOpt_[eitLINEAR_ANGLES]},
                  "Optimize linear angle parameters" },
                { "-dihedrals", FALSE, etBOOL, {&bOpt_[eitPROPER_DIHEDRALS]},
                  "Optimize proper dihedral parameters" },
                { "-impropers", FALSE, etBOOL, {&bOpt_[eitIMPROPER_DIHEDRALS]},
                  "Optimize improper dihedral parameters" },
                { "-vdw", FALSE, etBOOL, {&bOpt_[eitVDW]},
                  "Optimize van der Waals parameters" },
                { "-pairs",  FALSE, etBOOL, {&bOpt_[eitLJ14]},
                  "Optimize 1-4 interaction parameters" },
                { "-factor", FALSE, etREAL, {&factor_},
                  "Factor for something I forgot. Sorry." }
            };
            for (size_t i = 0; i < sizeof(pa)/sizeof(pa[0]); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs);
            TuneFc_.add_pargs(pargs);
        }

        const char *lot() const { return lot_; }

        /*! \brief
         *
         * Check whether molecules are supported by the force field.
         * Check whether all the atomtypes in the molecules are present
         * in the force field file. If not the molecules are removed.
         *
         * Routine also divides molecules over processors.
         */
        void checkSupport(FILE *fp);

        /*! \brief
         *
         * Fill parameter vector using ForceConstants which
         * is built based on Poldata.
         */
        void polData2TuneFc();

        /*! \brief
         *
         * Copy the optimized parameters back to Poldata
         */
        void tuneFc2PolData();

        /*! \brief
         *
         * Compute the dissociation energies for all the bonds.
         * Given all the bonds and the enthalpies of formation of all
         * molecules, we can approximate the dissociation enthalpy (D0 in the
         * Morse potential by least squares fitting the D0 to reproduce the
         * molecular energy (Delta H formation of molecule - Delta H formation of
         * the atoms). This is a crude approximation since all other energy
         * terms in the force field are ignored, however the dissociation
         * energy is the largest contribution to the molecular energy.
         */
        void getDissociationEnergy(FILE *fplog);

        void InitOpt(FILE *fplog);

        void optRun(FILE *fp, FILE *fplog,
                    int nrun,
                    const gmx_output_env_t *oenv,
                    const char *xvgconv,
                    const char *xvgepot);

        void Print(FILE *fp);

        void printResults(FILE                   *fp, 
                          char                   *title,
                          const char             *xvg,
                          const char             *HF_xvg,
                          const gmx_output_env_t *oenv,
                          bool                    bCheckOutliers);

        double calcDeviation();

        double objFunction(const double v[]);

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

        void broadcast(t_commrec *cr);
};

CommunicationStatus Optimization::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ForceConstants_.size());
        gmx_send_int(cr, dest, NonBondParams_.size());
        for (auto &fc : ForceConstants_)
        {
            fc.Send(cr, dest);
        }
        for (auto &nb : NonBondParams_)
        {
            nb.Send(cr, dest);
        }
    }
    return cs;
}

CommunicationStatus Optimization::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        int nfc = gmx_recv_int(cr, src);
        int nnb = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < nfc); n++)
        {
            alexandria::ForceConstants fc;
            cs = fc.Receive(cr, src);
            if (CS_OK == cs)
            {
                ForceConstants_.push_back(fc);
            }
        }
        for (int n = 0; (CS_OK == cs) && (n < nnb); n++)
        {
            alexandria::NonBondParams nb;
            cs = nb.Receive(cr, src);
            if (CS_OK == cs)
            {
                NonBondParams_.push_back(nb);
            }
        }
    }
    return cs;
}

void Optimization::broadcast(t_commrec *cr)
{
    const int src = 0;
    if (MASTER(cr))
    {
        for (int dest = 1; dest < cr->nnodes; dest++)
        {
            Send(cr, dest);
        }
    }
    else
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Updating ForceConstants and NonBondParams on node %d\n",
                    cr->nodeid);
        }
        Receive(cr, src);
    }
}

void Optimization::checkSupport(FILE *fp)
{
    auto ntotal  = mymols().size();
    auto nlocal  = 0;

    for (auto mymol =  mymols().begin(); mymol <  mymols().end(); )
    {
        if (mymol->eSupp_ != eSupportLocal)
        {
            mymol++;
            continue;
        }

        bool bSupport = true;
        for (int bt = 0; bSupport && (bt < eitNR); bt++)
        {
            int  ft;
            if (bOpt_[bt])
            {
                auto iType = static_cast<InteractionType>(bt);
                ft         = poldata().findForces(iType)->fType();
                bSupport   = (mymol->ltop_ != nullptr);
                for (int i = 0; bSupport && (i < mymol->ltop_->idef.il[ft].nr);
                     i += interaction_function[ft].nratoms+1)
                {
                    int                      ai, aj, ak, al;
                    std::string              aai, aaj, aak, aal;
                    std::vector<std::string> atoms;

                    ai  = mymol->ltop_->idef.il[ft].iatoms[i+1];
                    aj  = mymol->ltop_->idef.il[ft].iatoms[i+2];
                    if (!(poldata().atypeToBtype(*mymol->topology_->atoms.atomtype[ai], aai) &&
                          poldata().atypeToBtype(*mymol->topology_->atoms.atomtype[aj], aaj)))
                    {
                        bSupport = false;
                        if (debug)
                        {
                            fprintf(debug, "Cannot find bond types %s and %s in %s.\n",
                                    aai.c_str(), aaj.c_str(),
                                    mymol->molProp()->getMolname().c_str());
                        }
                    }
                    switch (iType)
                    {
                        case eitBONDS:
                        {
                            atoms = {aai, aaj};
                            auto fs = poldata().findForces(iType);
                            auto f  = fs->findForce(atoms);
                            if (fs->forceEnd() == f)
                            {
                                bSupport = false;
                                if (debug)
                                {
                                    fprintf(debug, "Cannot find bond %s-%s in %s.\n",
                                            aai.c_str(), aaj.c_str(),
                                            mymol->molProp()->getMolname().c_str());
                                }
                            }
                            break;
                        }
                        case eitANGLES:
                        case eitLINEAR_ANGLES:
                        {
                            ak  = mymol->ltop_->idef.il[ft].iatoms[i+3];
                            if (!poldata().atypeToBtype( *mymol->topology_->atoms.atomtype[ak], aak))
                            {
                                bSupport = false;
                                if (debug)
                                {
                                    fprintf(debug, "Cannot find bond types %s, %s and %s in %s.\n",
                                            aai.c_str(), aaj.c_str(), aak.c_str(),
                                            mymol->molProp()->getMolname().c_str());
                                }
                            }
                            else
                            {
                                atoms   = {aai, aaj, aak};
                                auto fs = poldata().findForces(iType);
                                auto f  = fs->findForce(atoms);
                                if (fs->forceEnd() == f)
                                {
                                    bSupport = false;
                                    if (debug)
                                    {
                                        fprintf(debug, "Cannot find %s %s-%s-%s in %s.\n",
                                                eitANGLES ? "angles" : "linear_angles",
                                                aai.c_str(), aaj.c_str(), aak.c_str(),
                                                mymol->molProp()->getMolname().c_str());
                                    }
                                }
                            }
                        }
                        break;
                        case eitPROPER_DIHEDRALS:
                        case eitIMPROPER_DIHEDRALS:
                        {
                            ak  = mymol->ltop_->idef.il[ft].iatoms[i+3];
                            al  = mymol->ltop_->idef.il[ft].iatoms[i+4];
                            if (!(poldata().atypeToBtype( *mymol->topology_->atoms.atomtype[ak], aak) &&
                                  poldata().atypeToBtype( *mymol->topology_->atoms.atomtype[al], aal)))
                            {
                                bSupport = false;
                                if (debug)
                                {
                                    fprintf(debug, "Cannot find bond types %s, %s, %s, and %s in %s\n",
                                            aai.c_str(), aaj.c_str(), aak.c_str(), aal.c_str(),
                                            mymol->molProp()->getMolname().c_str());
                                }
                            }
                            else
                            {
                                atoms   = {aai, aaj, aak, aal};
                                auto fs = poldata().findForces(iType);
                                auto f  = fs->findForce(atoms);
                                if (fs->forceEnd() == f)
                                {
                                    bSupport = false;
                                    if (debug)
                                    {
                                        fprintf(debug, "Cannot find %s dihedral %s-%s-%s-%s in %s\n",
                                                eitPROPER_DIHEDRALS ? "proper" : "improper",
                                                aai.c_str(), aaj.c_str(), aak.c_str(), aal.c_str(),
                                                mymol->molProp()->getMolname().c_str());
                                    }
                                }
                            }
                            break;
                            case eitPOLARIZATION:
                            case eitVDW:
                            case eitLJ14:
                            case eitVSITE2:
                            case eitVSITE3FAD:
                            case eitVSITE3OUT:
                            case eitCONSTR:
                            case eitNR:
                                break;
                        }
                    }
                }
            }
        }
        if (!bSupport)
        {
            fprintf(stderr, "No force field support for %s\n",
                    mymol->molProp()->getMolname().c_str());
            mymol =  mymols().erase(mymol);
        }
        else
        {
            mymol++;
            nlocal++;
        }
    }
    if (PAR(commrec()))
    {
        gmx_sumi(1, &nlocal, commrec());
    }
    if (nullptr != fp)
    {
        fprintf(fp, "%d out of %zu molecules have support in the force field.\n",
                nlocal, ntotal);
    }
}

void Optimization::polData2TuneFc()
{
    param_.clear();
    for (auto &fc : ForceConstants_)
    {
        if (eitBONDS == fc.interactionType())
        {
            for (auto b = fc.beginBN(); b  < fc.endBN(); ++b)
            {
                int n = 0;
                for (const auto &p : b->paramValues())
                {
                    if (n == 0)
                    {
                        De_.push_back(std::move(p));
                    }
                    else
                    {
                        param_.push_back(std::move(p));
                    }
                    n++;
                }
            }
        }
        else
        {
            for (auto b = fc.beginBN(); b  < fc.endBN(); ++b)
            {
                for (const auto &p : b->paramValues())
                {
                    param_.push_back(std::move(p));
                }
            }
        }
    }
    for (auto nbp : NonBondParams_)
    {
        for (auto at = nbp.beginAT(); at < nbp.endAT(); ++at)
        {
            for (const auto &p : at->paramValues())
            {
                param_.push_back(std::move(p));
            }
        }
    }
}

void Optimization::tuneFc2PolData()
{
    int n = 0;
    int m = 0;
    std::vector<std::string> atoms;

    for (auto &fc : ForceConstants_)
    {
        const auto iType = fc.interactionType();
        for (auto b = fc.beginBN(); b  < fc.endBN(); ++b)
        {
            char buf[STRLEN];
            buf[0] = '\0';

            if (iType == eitBONDS)
            {
                int j = 0;
                for (size_t p = 0; p < b->nParams(); p++)
                {
                    strncat(buf, " ", sizeof(buf)-1);
                    if (j == 0)
                    {
                        strncat(buf, gmx_ftoa(De_[m++]).c_str(), sizeof(buf)-1);
                    }
                    else
                    {
                        strncat(buf, gmx_ftoa(param_[n++]).c_str(), sizeof(buf)-1);
                    }
                    j++;
                }
            }
            else
            {
                for (size_t p = 0; p < b->nParams(); p++)
                {
                    strncat(buf, " ", sizeof(buf)-1);
                    strncat(buf, gmx_ftoa(param_[n++]).c_str(), sizeof(buf)-1);
                }
            }

            b->setParamString(buf);
            const auto bondtypes = gmx::splitString(b->name());
            switch (iType)
            {
                case eitBONDS:
                {
                    atoms    = {bondtypes[0], bondtypes[1]};
                    auto  fs = poldata().findForces(iType);
                    auto  f  = fs->findForce(atoms);
                    if (fs->forceEnd() != f)
                    {
                        f->setParams(buf);
                    }
                }
                break;
                case eitANGLES:
                case eitLINEAR_ANGLES:
                {
                    atoms    = {bondtypes[0], bondtypes[1], bondtypes[2]};
                    auto  fs = poldata().findForces(iType);
                    auto  f  = fs->findForce(atoms);
                    if (fs->forceEnd() != f)
                    {
                        f->setParams(buf);
                    }
                }
                break;
                case eitPROPER_DIHEDRALS:
                case eitIMPROPER_DIHEDRALS:
                {
                    atoms    = {bondtypes[0], bondtypes[1], bondtypes[2], bondtypes[3]};
                    auto  fs = poldata().findForces(iType);
                    auto  f  = fs->findForce(atoms);
                    if (fs->forceEnd() != f)
                    {
                        f->setParams(buf);
                    }
                }
                break;
                default:
                    gmx_fatal(FARGS, "Unsupported InteractionType %d",
                              static_cast<int>(fc.interactionType()));
            }
        }
    }
    for (auto &nbp : NonBondParams_)
    {
        for (auto at = nbp.beginAT(); at < nbp.endAT(); ++at)
        {
            char buf[STRLEN];
            buf[0] = '\0';
            for (size_t p = 0; p < at->nParams(); p++)
            {
                strncat(buf, " ", sizeof(buf)-1);
                strncat(buf, gmx_ftoa(param_[n++]).c_str(), sizeof(buf)-1);
            }
            at->setParamString(buf);
            const auto atype       = at->name();
            auto       fat         = poldata().findAtype(atype);
            if (fat != poldata().getAtypeEnd())
            {
                fat->setVdwparams(buf);
            }
        }
    }
}

void Optimization::getDissociationEnergy(FILE *fplog)
{
    std::vector<double>         rhs;
    std::vector<int>            ntest;
    std::vector<std::string>    ctest;

    int nD   = ForceConstants_[eitBONDS].nbad();
    int nMol =  mymols().size();

    if ((0 == nD) || (0 == nMol))
    {
        gmx_fatal(FARGS, "Number of variables is %d and number of molecules is %d",
                  nD, nMol);
    }

    MatrixWrapper a(nD, nMol);
    MatrixWrapper a_copy(nD, nMol);
    ntest.resize(nD, 0);
    ctest.resize(nD);

    fprintf(fplog, "There are %d different bondtypes to optimize the heat of formation\n", nD);
    fprintf(fplog, "There are %d (experimental) reference heat of formation.\n", nMol);

    auto fs  = poldata().findForces(eitBONDS);
    auto ftb = fs->fType();
    auto j   = 0;

    for (auto mymol =  mymols().begin(); mymol <  mymols().end(); mymol++, j++)
    {
        for (auto i = 0; i < mymol->ltop_->idef.il[ftb].nr;
             i += interaction_function[ftb].nratoms+1)
        {
            auto                     ai = mymol->ltop_->idef.il[ftb].iatoms[i+1];
            auto                     aj = mymol->ltop_->idef.il[ftb].iatoms[i+2];
            std::string              aai, aaj;
            std::vector<std::string> atoms;
            if (poldata().atypeToBtype(*mymol->topology_->atoms.atomtype[ai], aai) &&
                poldata().atypeToBtype(*mymol->topology_->atoms.atomtype[aj], aaj))
            {
                atoms  = {aai, aaj};
                auto f = fs->findForce(atoms);
                if (fs->forceEnd() != f)
                {
                    auto  gt  = f - fs->forceBegin();
                    auto  gti = ForceConstants_[eitBONDS].reverseIndex(gt);
                    a.set(gti, j, a.get(gti, j) + 1);
                    a_copy.set(gti, j, a.get(gti, j));
                    ntest[gti]++;
                    if (ctest[gti].empty())
                    {
                        char buf[STRLEN];
                        snprintf(buf, sizeof(buf), "%s-%s", aai.c_str(), aaj.c_str());
                        ctest[gti].assign(buf);
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "No parameters for bond %s-%s in the force field, atoms %s-%s mol %s",
                          aai.c_str(), aaj.c_str(),
                          *mymol->topology_->atoms.atomtype[ai],
                          *mymol->topology_->atoms.atomtype[aj],
                          mymol->molProp()->getIupac().c_str());
            }
        }
        rhs.push_back(std::move(-mymol->Emol_));
    }

    char buf[STRLEN];
    snprintf(buf, sizeof(buf), "Inconsistency in number of energies nMol %d != #rhs %d", nMol, static_cast<int>(rhs.size()));
    GMX_RELEASE_ASSERT(static_cast<int>(rhs.size()) == nMol, buf);

    auto nzero = std::count_if(ntest.begin(), ntest.end(), [](const int n)
                               {
                                   return n == 0;
                               });

    GMX_RELEASE_ASSERT(nzero == 0, "Inconsistency in the number of bonds in poldata and ForceConstants_");

    std::vector<double> Edissoc(nD);
    a.solve(rhs, &Edissoc);
    dump_csv(ctest,  mymols(), ntest, Edissoc, a_copy, rhs.data());

    for (size_t i = 0; i < ctest.size(); i++)
    {
        if (fplog)
        {
            fprintf(fplog, "Optimized dissociation energy for %8s with %4d copies to %g\n",
                    ctest[i].c_str(), ntest[i], Edissoc[i]);
        }
    }

    auto i = 0;
    for (auto b = ForceConstants_[eitBONDS].beginBN();
         b < ForceConstants_[eitBONDS].endBN(); ++b)
    {
        const auto       atoms = gmx::splitString(b->name());
        auto             fs    = poldata().findForces(eitBONDS);
        auto             f     = fs->findForce(atoms);
        GMX_RELEASE_ASSERT(fs->forceEnd() != f, "Cannot find my bonds");
        const auto       pp    = gmx::splitString(b->paramString());
        char             buf[256];

        /*De is assumed to be the first parameter. Not good!*/
        snprintf(buf, sizeof(buf), "%.2f  %s", Edissoc[i], pp[1].c_str());
        f->setParams(buf);
        b->setParamString(buf);
        i++;
    }
}

void Optimization::InitOpt(FILE *fplog)
{
    std::vector<unsigned int> fts;
    
    for (auto fs = poldata().forcesBegin();
         fs != poldata().forcesEnd(); fs++)
    {
        fts.push_back(std::move(fs->fType()));
    }

    for (size_t bt = 0; bt < fts.size(); bt++)
    {
        ForceConstants fc(bt, fts[bt], static_cast<InteractionType>(bt), bOpt_[bt]);
        fc.analyzeIdef(mymols(), poldata());
        fc.makeReverseIndex();
        fc.dump(fplog);
        ForceConstants_.push_back(std::move(fc));
    }

    if (ForceConstants_[eitBONDS].nbad() <= mymols().size())
    {
        getDissociationEnergy(fplog);
    }
    else
    {
        printf("\n"
               "WARNING: %zu molecule(s) is (are) not enough to calculate dissociation\n"
               "         energy for %zu bond type(s) using linear regression. Defualt\n"
               "         values from gentop.dat will be used as the initial guess.\n"
               "         Recomendation is to add more molecules having the same bond types.\n\n",
               mymols().size(), ForceConstants_[eitBONDS].nbad());
    }

    NonBondParams nbp(bOpt_[eitVDW], eitVDW);
    nbp.analyzeIdef(mymols(), poldata());
    nbp.makeReverseIndex();
    NonBondParams_.push_back(std::move(nbp));

    polData2TuneFc();

    orig_.resize(param_.size(), 0);
    best_.resize(param_.size(), 0);
    lower_.resize(param_.size(), 0);
    upper_.resize(param_.size(), 0);
    psigma_.resize(param_.size(), 0);
    pmean_.resize(param_.size(), 0);

    if (factor_ < 1)
    {
        factor_ = 1/factor_;
    }
    for (size_t i = 0; (i < param_.size()); i++)
    {
        best_[i]  = orig_[i] = param_[i];
        lower_[i] = orig_[i]/factor_;
        upper_[i] = orig_[i]*factor_;
    }
}

void Optimization::Print(FILE *fp)
{
    fprintf(fp, "Param        Orig        Best\n");
    for (size_t k = 0; k < param_.size(); k++)
    {
        fprintf(fp, "%-5zu  %10g  %10g\n", k, orig_[k], best_[k]);
    }
}

double Optimization::calcDeviation()
{
    int     j;
    FILE   *dbcopy;
    double  ener    = 0;
    double  optHF   = 0;
    double  spHF    = 0;
    double  deltaEn = 0;
    double  Emol    = 0;

    if (PAR(commrec()))
    {
        gmx_bool bFinal = final();
        gmx_bcast(sizeof(final()), &bFinal, commrec());
        if (bFinal)
        {
            setFinal();
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Begin communicating force parameters\n");
        fflush(debug);
    }
    if (PAR(commrec()) && !final())
    {
        poldata().broadcast(commrec());
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Done communicating force parameters\n");
        fflush(debug);
    }
    resetEnergies();

    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (final() && (mymol.eSupp_ == eSupportRemote)))
        {
            int natoms      = mymol.topology_->atoms.nr;
            int nOptSP      = mymol.molProp()->NOptSP();
            gmx_bool bpolar = (mymol.shellfc_ != nullptr);
            
            if (mymol.molProp()->getOptHF(&optHF))
            {
                for (const auto fc : ForceConstants_)
                {
                    if (fc.nbad() > 0)
                    {
                        mymol.UpdateIdef(poldata(), fc.interactionType());
                    }
                }

                for (const auto nbp : NonBondParams_)
                {
                    if (nbp.nAT() > 0)
                    {
                        mymol.UpdateIdef(poldata(), nbp.interactionType());
                    }
                }

                for (auto ei = mymol.molProp()->BeginExperiment();
                     ei < mymol.molProp()->EndExperiment(); ++ei)
                {
                    auto jtype = ei->getJobtype();

                    if (jtype == JOB_OPT || jtype == JOB_SP)
                    {
                        ei->getHF(&spHF);

                        deltaEn = spHF - optHF;
                        Emol    = mymol.Emol_ + deltaEn;

                        mymol.f_.clear();
                        mymol.f_.resize(2*natoms);

                        dbcopy = debug;
                        debug  = nullptr;
                       
                        mymol.changeCoordinate(ei, bpolar);
                        mymol.computeForces(debug, commrec());

                        debug          = dbcopy;
                        mymol.Force2_  = 0.0;

                        mymol.Ecalc_  = mymol.enerd_->term[F_EPOT];
                        ener          = gmx::square(mymol.Ecalc_ - Emol);

                        for (j = 0; j < natoms; j++)
                        {
                            mymol.Force2_ += iprod(mymol.f_[j], mymol.f_[j]);
                        }

                        mymol.Force2_ /= natoms;

                        if (jtype == JOB_OPT)
                        {
                            mymol.OptForce2_ = 0.0;
                            for (j = 0; j < natoms; j++)
                            {
                                mymol.OptForce2_ += iprod(mymol.f_[j], mymol.f_[j]);
                                copy_rvec(mymol.f_[j], mymol.optf_[j]);
                            }
                            mymol.OptForce2_   /= natoms;
                            increaseEnergy(ermsForce2, mymol.OptForce2_);
                            mymol.OptEcalc_     = mymol.enerd_->term[F_EPOT];
                        }

                        increaseEnergy(ermsEPOT, ener);

                        if (nullptr != debug)
                        {
                            fprintf(debug, "spHF: %g  optHF: %g  DeltaEn: %g\n", spHF, optHF, deltaEn);

                            fprintf(debug, "%s Chi2 %g Hform %g Emol %g  Ecalc %g Morse %g  "
                                    "Hangle %g Langle %g PDIHS %g IDIHS %g Coul %g LJ %g BHAM %g POL %g  Force2 %g\n",
                                    mymol.molProp()->getMolname().c_str(), ener, mymol.Hform_, Emol, mymol.Ecalc_,
                                    mymol.enerd_->term[F_MORSE], mymol.enerd_->term[F_UREY_BRADLEY],
                                    mymol.enerd_->term[F_LINEAR_ANGLES], mymol.enerd_->term[F_PDIHS],
                                    mymol.enerd_->term[F_IDIHS], mymol.enerd_->term[F_COUL_SR],
                                    mymol.enerd_->term[F_LJ], mymol.enerd_->term[F_BHAM],
                                    mymol.enerd_->term[F_POLARIZATION], mymol.Force2_);
                        }
                    }
                }
                setEnergy(ermsEPOT, energy(ermsEPOT) / nOptSP);
            }
            else
            {
                gmx_fatal(FARGS, "There is no optimized structure for %s\n",
                          mymol.molProp()->getMolname().c_str());
            }
        }
    }
    /* Compute E-bounds */
    if (MASTER(commrec()))
    {
        for (size_t j = 0; j < param_.size(); j++)
        {
            if (param_[j] < lower_[j])
            {
                increaseEnergy(ermsBOUNDS, gmx::square(param_[j]-lower_[j]));
            }
            else if (param_[j] > upper_[j])
            {
                increaseEnergy(ermsBOUNDS, gmx::square(param_[j]-upper_[j]));
            }
        }
    }
    /* Global sum energies */
    sumEnergies();
    normalizeEnergies();
    printEnergies(debug);
    return energy(ermsTOT);
}

double Optimization::objFunction(const double v[])
{
    double rms = 0;
    size_t np  = param_.size();

    for (size_t i = 0; i < np; i++)
    {
        param_[i] = v[i];
    }

    tuneFc2PolData(); /* Copy parameters to topologies */
    rms = calcDeviation();

    return rms;
}

void Optimization::optRun(FILE *fp, FILE *fplog,
                          int nrun,
                          const gmx_output_env_t *oenv,
                          const char *xvgconv,
                          const char *xvgepot)
{

    std::vector<double> optx, opts, optm;
    double              chi2, chi2_min;
    gmx_bool            bMinimum = false;
    
    auto func = [&] (const double v[]) {
        return objFunction(v);
    };
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, (nrun*TuneFc_.maxIter()*param_.size()));
            }
            
        }
        chi2 = chi2_min  = GMX_REAL_MAX;
        TuneFc_.setFunc(func, param_, lower_, upper_, &chi2);
        TuneFc_.Init(xvgconv, xvgepot, oenv);

        for (int n = 0; n < nrun; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun);
            }
            TuneFc_.simulate();
            TuneFc_.getBestParam(optx);
            TuneFc_.getPsigma(opts);
            TuneFc_.getPmean(optm);

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
            TuneFc_.setParam(best_);
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
        int niter = gmx_recv_int(commrec(), 0);
        for (int n = 0; n < niter + 2; n++)
        {
            chi2 = calcDeviation();
        }
    }
    setFinal();
    if (MASTER(commrec()))
    {
        chi2 = calcDeviation();
        printEnergies(fp);
        printEnergies(fplog);
    }
}

void print_mols(FILE                          *fp, 
                std::vector<alexandria::MyMol> mol,
                gmx_bool                       bForce, 
                gmx_bool                       bMtop)
{
    int j, k;

    for (auto mi = mol.begin(); (mi < mol.end()); mi++)
    {
        fprintf(fp, "%-30s  %d\n", mi->molProp()->getMolname().c_str(), mi->topology_->atoms.nr);
        for (j = 0; j < mi->topology_->atoms.nr; j++)
        {
            fprintf(fp, "  %-5s  %-5s  q = %10g", 
                    *(mi->topology_->atoms.atomname[j]),
                    *(mi->topology_->atoms.atomtype[j]), 
                    mi->topology_->atoms.atom[j].q);
            if (bForce)
            {
                fprintf(fp, "   f = %8.3f  %8.3f  %8.3f",
                        mi->optf_[j][XX], mi->optf_[j][YY], mi->optf_[j][ZZ]);
            }
            fprintf(fp, "\n");
        }
        if (bForce)
        {
            for (k = 0; k < F_NRE; k++)
            {
                if ((mi->enerd_->term[k] != 0) ||
                    (mi->mtop_->moltype[0].ilist[k].nr > 0))
                {
                    fprintf(fp, "%s %d %g\n", interaction_function[k].name,
                            mi->mtop_->moltype[0].ilist[k].nr,
                            mi->enerd_->term[k]);
                }
            }
        }
        if (bMtop)
        {
            pr_mtop(fp, 0, mi->molProp()->getMolname().c_str(), mi->mtop_, true, false);
        }
    }
}


void Optimization::printResults(FILE                   *fp, 
                                char                   *title,
                                const char             *hform_xvg,
                                const char             *HF_xvg,
                                const gmx_output_env_t *oenv,
                                bool                    bCheckOutliers)
{
    int         N, i;
    real        a, b;
    real        da, db;
    real        chi2, Rfit;
    double      msd;
    FILE       *xfp, *hfp;
    gmx_stats_t gst;

    gst = gmx_stats_init();
    if (nullptr != hform_xvg)
    {
        xfp = xvgropen(hform_xvg, "Entalpy of Formation", "Experiment (kJ/mol)", "Calculated (kJ/mol)", oenv);
    }
    if (nullptr != HF_xvg)
    {
        hfp = xvgropen(HF_xvg, "Eenergy", "B3LYP/aug-cc-pVTZ (kJ/mol)", "Alexandria (kJ/mol)", oenv);
    }
    fprintf(fp, "%s\n", title);
    fprintf(fp, "Nr.   %-30s %10s %10s %10s %10s %10s\n", "Molecule", "DHf@298K", "Emol@0K", "Delta E", "rms F", "Outlier?");
    msd = 0;
    i   = 0;
    for (auto mi = mymols().begin(); mi < mymols().end(); mi++, i++)
    {
        real DeltaE = mi->OptEcalc_ - mi->Emol_;
        fprintf(fp, "%-5d %-30s %10g %10g %10g %10g %-10s\n",
                i, 
                mi->molProp()->getMolname().c_str(),
                mi->Hform_, 
                mi->Emol_, 
                DeltaE, 
                sqrt(mi->OptForce2_),
                (bCheckOutliers && (fabs(DeltaE) > 1000)) ? "XXX" : "");
        msd += gmx::square(DeltaE);
        gmx_stats_add_point(gst, mi->Hform_, mi->Hform_ + DeltaE, 0, 0);
        if (nullptr != hform_xvg)
        {
            fprintf(xfp, "%10g  %10g\n", mi->Hform_, mi->Hform_ + DeltaE);
        }
        if (nullptr != HF_xvg)
        {
            double   spHF     = 0;
            double   optHF    = 0;
            double   deltaEn  = 0;
            gmx_bool bpolar   = (mi->shellfc_ != nullptr);
            mi->molProp()->getOptHF(&optHF);
            for (auto ei = mi->molProp()->BeginExperiment(); ei < mi->molProp()->EndExperiment(); ++ei)
            {
                auto jtype = ei->getJobtype();
                if (jtype == JOB_SP)
                {
                    ei->getHF(&spHF);
                    deltaEn = spHF - optHF;
                    mi->f_.clear();
                    mi->f_.resize(2*mi->topology_->atoms.nr);                  
                    mi->changeCoordinate(ei, bpolar);
                    mi->computeForces(nullptr, commrec());                                   
                    fprintf(hfp, "%10g  %10g\n", mi->Emol_ + deltaEn, mi->enerd_->term[F_EPOT]);
                }
            }
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "RMSD is %g kJ/mol for %zu molecules.\n\n", sqrt(msd/mymols().size()), mymols().size());
    fflush(fp);
    if (nullptr != hform_xvg)
    {
        xvgrclose(xfp);
        do_view(oenv, hform_xvg, nullptr);
    }
    if (nullptr != HF_xvg)
    {
        xvgrclose(hfp);
        do_view(oenv, HF_xvg, nullptr);
    }

    gmx_stats_get_ab(gst, 1, &a, &b, &da, &db, &chi2, &Rfit);
    gmx_stats_get_npoints(gst, &N);
    fprintf(fp, "Regression analysis fit to y = ax + b:\n");
    fprintf(fp, "a = %.3f  b = %3f  R2 = %.1f%%  chi2 = %.1f N = %d\n", a, b, Rfit*100, chi2, N);
    gmx_stats_free(gst);
    fflush(fp);
}
}

int alex_tune_fc(int argc, char *argv[])
{
    const char           *desc[] = {
        "tune_fc read a series of molecules and corresponding experimental",
        "heats of formation from a file, and tunes parameters in an algorithm",
        "until the experimental energies are reproduced by the force field.[PAR]",
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

    t_filenm              fnm[] = {
        { efDAT, "-f",     "allmols",    ffREAD  },
        { efDAT, "-d",     "gentop",     ffOPTRD },
        { efDAT, "-o",     "tune_fc",    ffWRITE },
        { efDAT, "-sel",   "molselect",  ffREAD  },
        { efXVG, "-table", "table",      ffOPTRD },
        { efLOG, "-g",     "tune_fc",    ffWRITE },
        { efXVG, "-x",     "hform-corr", ffWRITE },
        { efXVG, "-hf",    "hf-corr",    ffWRITE },
        { efXVG, "-conv",  "param-conv", ffWRITE },
        { efXVG, "-epot",  "param-epot", ffWRITE }
    };

    const int             NFILE         = asize(fnm);

    int                   nrun          = 1;
    int                   reinit        = 0;
    int                   compress      = 0;
    int                   nmultisim     = 0;
    char                 *opt_elem      = nullptr;
    char                 *const_elem    = nullptr;
    gmx_bool              bWeighted     = false;
    gmx_bool              bZPE          = false;
    gmx_bool              bZero         = true;

    t_pargs               pa[]         = {
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do optimization in multiple simulation" },
        { "-reinit",  FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-nrun",    FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of elements to optimize, e.g. \"H C Br\". The other available elements in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of elements to include but keep constant, e.g. \"O N\". These elements from gentop.dat are left unmodified" },
        { "-weight", FALSE, etBOOL, {&bWeighted},
          "Perform a weighted fit, by using the errors in the dipoles presented in the input file. This may or may not improve convergence." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" }
    };

    FILE                 *fp;
    gmx_output_env_t     *oenv;
    time_t                my_t;
    MolSelect             gms;

    std::vector<t_pargs>  pargs;
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::Optimization opt;
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
        fprintf(fp, "# alexandria is part of G R O M A C S:\n#\n");
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
             false,
             opt.bOpt(eitLJ14),
             opt.bOpt(eitPROPER_DIHEDRALS),
             bZPE,
             false,
             tabfn);
    
    opt.checkSupport(fp);
    
    if (nullptr != fp)
    {
        fprintf(fp, "In the total data set of %d molecules we have:\n",
                static_cast<int>(opt.mymols().size()));
    }
   
    if (MASTER(opt.commrec()))
    {
        opt.InitOpt(fp);
        print_mols(fp, opt.mymols(), false, false);
    }

    if (PAR(opt.commrec()))
    {
        opt.broadcast(opt.commrec());
    }

    opt.calcDeviation();
    
    if (MASTER(opt.commrec()))
    {
        opt.printResults(fp, (char *)"Before optimization", nullptr, nullptr, oenv, false);
    }

    opt.optRun(MASTER(opt.commrec()) ? stderr : nullptr,
               fp,
               nrun,
               oenv,
               opt2fn("-conv", NFILE, fnm),
               opt2fn("-epot", NFILE, fnm));

    if (MASTER(opt.commrec()))
    {
        opt.tuneFc2PolData();
        print_mols(fp, opt.mymols(), true, false);
        opt.printResults(fp, (char *)"After optimization", 
                         opt2fn("-x", NFILE, fnm), 
                         opt2fn("-hf", NFILE, fnm), 
                         oenv, true);
        writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), compress);
        done_filenms(NFILE, fnm);
        gmx_ffclose(fp);
    }

    return 0;
}
