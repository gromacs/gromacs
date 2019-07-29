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
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

#include "alex_modules.h"
#include "communication.h"
#include "gentop_core.h"
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
static void dump_csv(const std::vector<std::string>        &ctest,
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
        p_.resize(np, 0.0);
        for (int n = 0; n < np; n++)
        {
            p_[n] = gmx_recv_double(cr, src);
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
    p_.resize(p.size(), 0.0);
    for (size_t d = 0; d < p.size(); d++)
    {
        p_[d] = gmx::doubleFromString(p[d].c_str());
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

        void analyzeIdef(const std::vector<MyMol> &mm,
                         const Poldata            &pd);

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

void NonBondParams::analyzeIdef(const std::vector<MyMol> &mm,
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
                auto at  = *mymol.topology_->atoms.atomtype[i];
                auto fat = pd.findAtype(at);
                if (fat != pd.getAtypeEnd())
                {
                    std::string typeName    = at;
                    std::string paramString = fat->getVdwparams();
                    int         index       = fat - pd.getAtypeBegin();
                    auto        myat        = std::find_if(at_.begin(), at_.end(),
                                                           [typeName](const AtomTypes &at)
                                                           {
                                                               return (at.name().compare(typeName) == 0);
                                                           });
                    if (myat != at_.end())
                    {
                        myat->inc();
                    }
                    else
                    {
                        AtomTypes at(1, typeName, paramString, index);
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
                  int                ftype,
                  const std::string &name,
                  const std::string &params,
                  int                index,
                  double             bondorder = 0)

            :
              ncopies_(ncopies),
              ftype_(ftype),
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
        //! Fucntion type for this particular bond
        int                 ftype_;
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
        p_.resize(np, 0.0);
        for (int n = 0; n < np; n++)
        {
            p_[n] = gmx_recv_double(cr, src);
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
    int np = 1+static_cast<int>(p.size());
    if ((ftype_ == F_UREY_BRADLEY && np != 3) ||
        (ftype_ != F_UREY_BRADLEY && np != NRFPA(ftype_)))
    {
        gmx_fatal(FARGS, "Number of parameters (%d) in gentop.dat does not match function type %s (%d)", np, 
                  interaction_function[ftype_].name, NRFPA(ftype_));
    }
    p_.resize(p.size(), 0.0);
    for (size_t d = 0; d < p.size(); d++)
    {
        p_[d] = gmx::doubleFromString(p[d].c_str());
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

        void analyzeIdef(const std::vector<MyMol> &mm,
                         const Poldata            &pd);

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

void ForceConstants::analyzeIdef(const std::vector<MyMol> &mm,
                                 const Poldata            &pd)
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
                        BondNames bn(1, ftype_, buf, params, index, 0);
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
    if (bOpt_)
    {
        int ntot = 0;
        fprintf(fp, "Interaction  Bondtypes             Copies Poldata entry\n");
        const char *name = interaction_function[ftype_].name;
        for (const auto &i : bn_)
        {
            fprintf(fp, "%-10s  %-20s  %5d  %5d\n",
                    name, i.name().c_str(), i.nCopies(), i.poldataIndex());
            ntot += i.nCopies();
        }
        fprintf(fp, "%d out of %d %s types will be optimized.\n",
                static_cast<int>(bn_.size()), ntot, name);
    }
}

class PoldataUpdate
{
public:
    PoldataUpdate() {}
    PoldataUpdate(InteractionType     iType,
                  int                 index,
                  std::vector<double> params,
                  std::string         paramString) : iType_(iType), index_(index), params_(params), paramString_(paramString)
    {}
    
    /*! \brief
     * Implement the changes in Poldata.
     *
     * \param[inout] pd The Poldata structure
     */
    void execute(Poldata &pd);
    /*! \brief
     * Dump the contents of the structure to a file.
     * \param[in] fp File pointer to dumpt to if not nullptr
     */
    void dump(FILE *fp) const;
    
    CommunicationStatus Send(t_commrec *cr, int dest);
    
    CommunicationStatus Receive(t_commrec *cr, int src);

private:
    InteractionType     iType_;
    int                 index_;
    std::vector<double> params_;
    std::string         paramString_;
};

void PoldataUpdate::execute(Poldata &pd)
{
    if (iType_ == eitVDW)
    {
        auto fat = pd.getAtypeBegin() + index_;
        fat->setVdwparams(paramString_);
        fat->setModified(true);
    }
    else
    {
        auto fs = pd.findForces(iType_);
        auto f  = fs->forceBegin() + index_;
        f->setParams(paramString_);
        f->setModified(true);
    }
}

void PoldataUpdate::dump(FILE *fp) const
{
    if (nullptr != fp)
    {
        fprintf(fp, "iType: %s index: %d paramString: %s params:",
                iType2string(iType_), index_, paramString_.c_str());
        for (auto &p : params_)
        {
            fprintf(fp, " %g", p);
        }
        fprintf(fp, "\n");        
    }
}

CommunicationStatus PoldataUpdate::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        dump(debug);
        gmx_send_int(cr, dest, static_cast<int>(iType_));
        gmx_send_int(cr, dest, index_);
        gmx_send_str(cr, dest, &paramString_);
        gmx_send_int(cr, dest, static_cast<int>(params_.size()));
        for(const auto &p : params_)
        {
            gmx_send_double(cr, dest, p);
        }
    }
    return cs;
}

CommunicationStatus PoldataUpdate::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        iType_ = static_cast<InteractionType>(gmx_recv_int(cr, src));
        index_ = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &paramString_);
        int n  = gmx_recv_int(cr, src);
        params_.clear();
        for(int i = 0; i < n; i++)
        {
            params_.push_back(gmx_recv_double(cr, src));
        }
        dump(debug);
    }
    return cs;
}

class Optimization : public MolGen
{
    using param_type = std::vector<double>;

    private:
        std::vector<ForceConstants> ForceConstants_;
        std::vector<NonBondParams>  NonBondParams_;
        param_type                  param_;
        param_type                  De_;
        std::vector<int>            iOpt_;
        std::vector<PoldataUpdate>  poldataUpdates;
        real                        factor_     = 1;
        bool                        OptimizeDe_ = true;
        bool                        bDissoc_    = true;
        real                        beta0_      = 20;
        real                        D0_         = 350;
        real                        beta_min_   = 10;
        real                        D0_min_     = 50;
        real                        w_dhf_      = 1;
        const char                 *lot_        = nullptr;
        int                         nrun_       = 1;
        Bayes <double>              TuneFc_;
    public:


        /*! \brief
         *
         * Constructor
         */
        Optimization()
        {
            iOpt_.resize(eitNR, false);
            iOpt_[eitBONDS] = true;
        };

        /*! \brief
         *
         * Destructor
         */
        ~Optimization() {};

        bool iOpt(int i) const { return iOpt_[i]; }
        void add_pargs(std::vector<t_pargs> *pargs)
        {
            t_pargs pa[] =
            {
                { "-beta0", FALSE, etREAL, {&beta0_},
                  "Reset the initial beta for Morse potentials to this value, independent of gentop.dat. If value is <= 0 gentop.dat value is used." },
                { "-beta_min", FALSE, etREAL, {&beta_min_},
                  "Minimum value for beta in Morse potential" },
                { "-D0", FALSE, etREAL, {&D0_},
                  "Reset the initial D for Morse potentials to this value, independent of gentop.dat. If value is <= 0 gentop.dat value is used." },
                { "-DO_min", FALSE, etREAL, {&D0_min_},
                  "Minimum value for D0 in Morse potential" },
                { "-dissoc",  FALSE, etBOOL, {&bDissoc_},
                  "Derive dissociation energy from the enthalpy of formation. If not chosen, the dissociation energy will be read from the gentop.dat file." },
                { "-optimizeDe", FALSE, etBOOL, {&OptimizeDe_},
                  "Optimize the dissociation energy or keep it constant." },
                { "-weight_dhf", FALSE, etREAL, {&w_dhf_},
                  "Fitting weight of the minimum energy structure, representing the enthalpy of formation relative to high energy structures." },
                { "-nrun",    FALSE, etINT,  {&nrun_},
                  "This many runs will be done, before each run a complete randomization will be done" },
                { "-bonds",   FALSE, etINT, {&(iOpt_[eitBONDS])},
                  "Optimize bond parameters" },
                { "-angles",  FALSE, etINT, {&(iOpt_[eitANGLES])},
                  "Optimize angle parameters" },
                { "-langles", FALSE, etINT, {&(iOpt_[eitLINEAR_ANGLES])},
                  "Optimize linear angle parameters" },
                { "-dihedrals", FALSE, etINT, {&(iOpt_[eitPROPER_DIHEDRALS])},
                  "Optimize proper dihedral parameters" },
                { "-impropers", FALSE, etINT, {&(iOpt_[eitIMPROPER_DIHEDRALS])},
                  "Optimize improper dihedral parameters" },
                { "-vdw", FALSE, etINT, {&(iOpt_[eitVDW])},
                  "Optimize van der Waals parameters" },
                { "-pairs",  FALSE, etINT, {&(iOpt_[eitLJ14])},
                  "Optimize 1-4 interaction parameters" },
                { "-factor", FALSE, etREAL, {&factor_},
                  "Parameters will be taken within the limit factor*x - x/factor." }
            };
            for (size_t i = 0; i < sizeof(pa)/sizeof(pa[0]); i++)
            {
                pargs->push_back(pa[i]);
            }
            addOptions(pargs, etuneFC);
            TuneFc_.add_pargs(pargs);
        }
        void optionsFinished()
        {
            MolGen::optionsFinished();
            TuneFc_.setBounds(weight(ermsBOUNDS) > 0);
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
         *
         * \param[in] changed Boolean list stating whether a parameter has changed
         */
        void tuneFc2PolData(const std::vector<bool> &changed);

        /*! \brief
         * Broadcast changes in Poldata to the 
         * slaves when in parallel.
         */
        void broadcastPoldataUpdate();
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

        /*! \brief
         * Initialiaze the optimization algorithm.
         * \param[in] fplog Log file for dumping information.
         */
        void InitOpt(FILE *fplog);

        /*! \brief
         * Do the actual optimization.
         * \param[in] fp     FILE pointer for urgent messages
         * \param[in] fplog  FILE pointer for logging
         * \param[in] oenv   Output environment for managing xvg files etc.
         * \param[in] xvgconv Output file monitoring parameters
         * \param[in] xvgepot Output file monitoring penalty function
         * \return true if better parameters were found.
         */
        bool optRun(FILE                   *fp,
                    FILE                   *fplog,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);

        void printResults(FILE                   *fp,
                          const char             *title,
                          const char             *xvg,
                          const char             *HF_xvg,
                          const gmx_output_env_t *oenv);

        void calcDeviation();

        double objFunction(const double v[]);
        
        size_t nParams() const { return param_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

        void broadcast();
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
                NonBondParams_.push_back(std::move(nb));
            }
        }
    }
    return cs;
}

void Optimization::broadcastPoldataUpdate()
{
    t_commrec *cr = commrec();
    if (!PAR(cr))
    {
        return;
    }
    
    if (MASTER(cr))
    {
        for(int i = 1; i < cr->nnodes; i++)
        {
            gmx_send_int(cr, i, static_cast<int>(poldataUpdates.size()));
            for (auto &p : poldataUpdates)
            {
                p.Send(cr, i);
            }
        }
    }
    else
    {
        int n = gmx_recv_int(cr, 0);
        poldataUpdates.resize(n);
        for (int i = 0; i < n; i++)
        {
            poldataUpdates[i].Receive(cr, 0);
        }
        if (debug)
        {
            fprintf(debug, "broadcastPoldataUpdate: received %d updates\n",
                    n);
        }
    }
}

void Optimization::broadcast()
{
    t_commrec *cr = commrec();
    if (!PAR(cr))
    {
        return;
    }
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
    poldata().broadcast(cr);
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
            if (iOpt_[bt])
            {
                auto iType = static_cast<InteractionType>(bt);
                auto fs    = poldata().findForces(iType);
                if (fs == poldata().forcesEnd())
                {
                    continue;
                }
                ft         = fs->fType();
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
                                                (iType == eitPROPER_DIHEDRALS) ? "proper" : "improper",
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
    De_.clear();
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
                    if (n == 0 && !OptimizeDe_)
                    {
                        De_.push_back(p);
                    }
                    else
                    {
                        param_.push_back(p);
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
                    param_.push_back(p);
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
                param_.push_back(p);
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "poldata2TuneFc: Copied %zu parameters.\n",
                param_.size());
    }
}

void Optimization::tuneFc2PolData(const std::vector<bool> &changed)
{
    size_t n   = 0;
    size_t nDe = 0;

    poldataUpdates.clear();
    for (auto &fc : ForceConstants_)
    {
        const auto iType = fc.interactionType();
        for (auto b = fc.beginBN(); b  < fc.endBN(); ++b)
        {
            std::vector<double> thisParam;
            std::string         paramString;
            bool                bondChanged = false;
            for (size_t p = 0; p < b->nParams(); p++)
            {
                if (iType == eitBONDS && p == 0 && !OptimizeDe_)
                {
                    paramString.append(gmx::formatString(" %g", De_[nDe++]));
                }
                else
                {
                    bondChanged = bondChanged || changed[n];
                    thisParam.push_back(param_[n]);
                    paramString.append(gmx::formatString(" %g", param_[n++]));
                }
            }
            if (true || bondChanged)
            {
                b->setParamString(paramString);
                poldataUpdates.push_back(PoldataUpdate(iType, b->poldataIndex(),
                                                       thisParam, paramString));
            }
        }
    }
    for (auto &nbp : NonBondParams_)
    {
        for (auto at = nbp.beginAT(); at < nbp.endAT(); ++at)
        {
            std::vector<double> thisParam;
            bool                bondChanged = false;
            std::string         paramString;
            for (size_t p = 0; p < at->nParams(); p++)
            {
                bondChanged = bondChanged || changed[n];
                thisParam.push_back(param_[n]);
                paramString.append(gmx::formatString(" %g", param_[n++]));
            }
            if (true || bondChanged)
            {
                at->setParamString(paramString);
                poldataUpdates.push_back(PoldataUpdate(eitVDW, 
                                                       at->poldataIndex(),
                                                       thisParam, paramString));
            }
        }
    }
    GMX_RELEASE_ASSERT(n == param_.size(), "Number of parameters set should be equal to the length of the parameter array");
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
        rhs.push_back(-mymol->Emol_);
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
    if (debug)
    {
        dump_csv(ctest,  mymols(), ntest, Edissoc, a_copy, rhs.data());
    }
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

        /* TODO: De is assumed to be the first parameter. Please fix. */
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
         fs != poldata().forcesEnd(); ++fs)
    {
        fts.push_back(fs->fType());
    }

    for (size_t bt = 0; bt < fts.size(); bt++)
    {
        ForceConstants fc(bt, fts[bt], static_cast<InteractionType>(bt), iOpt_[bt]);
        fc.analyzeIdef(mymols(), poldata());
        fc.makeReverseIndex();
        fc.dump(fplog);
        ForceConstants_.push_back(std::move(fc));
    }

    if (bDissoc_)
    {
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
    }
    NonBondParams nbp(iOpt_[eitVDW], eitVDW);
    nbp.analyzeIdef(mymols(), poldata());
    nbp.makeReverseIndex();
    NonBondParams_.push_back(std::move(nbp));

    polData2TuneFc();

    if (factor_ < 1)
    {
        factor_ = 1/factor_;
    }
}

void Optimization::calcDeviation()
{
    if (PAR(commrec()))
    {
        gmx_bool bFinal = final();
        gmx_bcast(sizeof(final()), &bFinal, commrec());
        if (bFinal)
        {
            setFinal();
        }
        broadcastPoldataUpdate();
    }
    for (auto &pUpd : poldataUpdates)
    {
        pUpd.execute(poldata());
    }
    if (debug)
    {
        static bool writtenPd = false;
        if (!writtenPd)
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "pd%d.dat", commrec()->nodeid);
            writePoldata(buf, poldata(), false);
            writtenPd = true;
        }
    }
    resetEnergies();
    double nCalc = 0;
    double ePot2 = 0;
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (final() && (mymol.eSupp_ == eSupportRemote)))
        {
            int      natoms = mymol.topology_->atoms.nr;
            gmx_bool bpolar = (mymol.shellfc_ != nullptr);
            double   optHF;
            if (mymol.molProp()->getOptHF(&optHF))
            {
                for (const auto fc : ForceConstants_)
                {
                    if (fc.nbad() > 0)
                    {
                        mymol.UpdateIdef(poldata(), fc.interactionType(), true);
                    }
                }

                for (const auto nbp : NonBondParams_)
                {
                    if (nbp.nAT() > 0)
                    {
                        mymol.UpdateIdef(poldata(), nbp.interactionType(), true);
                    }
                }
                mymol.f_.resizeWithPadding(natoms);
                mymol.optf_.resizeWithPadding(natoms);
                for (auto ei = mymol.molProp()->BeginExperiment();
                     ei < mymol.molProp()->EndExperiment(); ++ei)
                {
                    auto jtype = ei->getJobtype();
                    // Exclude other experimental data points!
                    if (jtype == JOB_OPT || jtype == JOB_SP)
                    {
                        double spHF;
                        if (!ei->getHF(&spHF))
                        {
                            continue;
                        }

                        double deltaEn = spHF - optHF;
                        double Emol    = mymol.Emol_ + deltaEn;

                        FILE *dbcopy = debug;
                        debug  = nullptr;
                        mymol.changeCoordinate(ei, bpolar);
                        mymol.computeForces(debug, commrec());
                        debug  = dbcopy;

                        real Ecalc     = mymol.enerd_->term[F_EPOT];
                        real EnerDiff2 = gmx::square(Ecalc - Emol);

                        // Force is added only for the opt geometry
                        if (jtype == JOB_OPT)
                        {
                            mymol.OptForce2_ = 0.0;
                            for (int j = 0; j < natoms; j++)
                            {
                                mymol.OptForce2_ += iprod(mymol.f_[j], mymol.f_[j]);
                                copy_rvec(mymol.f_[j], mymol.optf_[j]);
                            }
                            mymol.OptForce2_   /= natoms;
                            increaseEnergy(ermsForce2, mymol.OptForce2_);
                            mymol.OptEcalc_     = Ecalc;
                            // Weighted energy penalty for opt geometry
                            ePot2 += w_dhf_ * EnerDiff2; 
                            nCalc += w_dhf_;
                        }
                        else
                        {
                            // Energy penalty for sp geometry
                            ePot2 += EnerDiff2; 
                            nCalc += 1.0;
                        }
                        if (nullptr != debug)
                        {
                            int angleType = poldata().findForces(eitANGLES)->fType();
                            int vdwType   = poldata().getVdwFtype();
                            fprintf(debug, "spHF: %g  optHF: %g  DeltaEn: %g\n",
                                    spHF, optHF, deltaEn);
                            fprintf(debug, "%s Chi2 %g Hform %g Eqm %g  Ecalc %g Morse %g  "
                                    "Angle %g Langle %g PDIHS %g IDIHS %g Coul %g VdW %g POL %g  Force %g\n",
                                    mymol.molProp()->getMolname().c_str(),
                                    EnerDiff2,
                                    mymol.Hform_,
                                    Emol,
                                    Ecalc,
                                    mymol.enerd_->term[F_MORSE],
                                    mymol.enerd_->term[angleType],
                                    mymol.enerd_->term[F_LINEAR_ANGLES],
                                    mymol.enerd_->term[F_PDIHS],
                                    mymol.enerd_->term[F_IDIHS],
                                    mymol.enerd_->term[F_COUL_SR],
                                    mymol.enerd_->term[vdwType],
                                    mymol.enerd_->term[F_POLARIZATION],
                                    std::sqrt(mymol.OptForce2_));
                        }
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "There is no optimized structure for %s\n",
                          mymol.molProp()->getMolname().c_str());
            }
        }
    }
    if (PAR(commrec()))
    {
        gmx_sumd(1, &nCalc, commrec());
    }
    setEnergy(ermsEPOT, nMolSupport()*ePot2/nCalc);
    normalizeEnergies();
    sumEnergies();
    printEnergies(debug);
}

double Optimization::objFunction(const double v[])
{
    size_t np = param_.size();
    std::vector<bool> changed(np, false);
    for (size_t i = 0; i < np; i++)
    {
        if (param_[i] != v[i])
        {
            param_[i]  = v[i];
            changed[i] = true;
        }
    }
    tuneFc2PolData(changed);
    calcDeviation();

    return energy(ermsTOT);
}

bool Optimization::optRun(FILE                   *fp,
                          FILE                   *fplog,
                          const gmx_output_env_t *oenv,
                          const char             *xvgconv,
                          const char             *xvgepot)
{
    bool bMinimum = false;
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            // Tell the slave nodes how many times they have
            // to run calcDeviation.
            int niter = 3+nrun_*TuneFc_.maxIter()*param_.size();
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, niter);
            }
        }
        std::vector<double> pmean, psigma, best, lower, upper;
        pmean.resize(param_.size(), 0.0);
        psigma.resize(param_.size(), 0.0);
        best.resize(param_.size(), 0.0);
        lower.resize(param_.size(), 0);
        upper.resize(param_.size(), 0);
        for (size_t i = 0; (i < param_.size()); i++)
        {
            lower[i] = param_[i]/factor_;
            upper[i] = param_[i]*factor_;
        }
    
        double   chi2_min = objFunction(param_.data());
        double   chi2     = chi2_min;

        auto func = [&] (const double v[]) {
            return objFunction(v);
        };
        TuneFc_.setFunc(func, param_, lower, upper, &chi2);
        TuneFc_.Init(xvgconv, xvgepot, oenv);

        for (int n = 0; n < nrun_; n++)
        {
            if ((nullptr != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n, nrun_);
            }
            TuneFc_.MCMC();
            if (chi2 < chi2_min)
            {
                chi2_min = chi2;
                bMinimum = true;
                TuneFc_.getBestParam(best);
                TuneFc_.getPsigma(psigma);
                TuneFc_.getPmean(pmean);
                // For next run of optimization we start from the best.
                TuneFc_.setParam(best);
            }
        }
        if (bMinimum)
        {
            setFinal();
            // This call copies data to poldata as well.
            double chi2 = objFunction(best.data());
            if (fplog)
            {
                fprintf(fplog, "\nLowest RMSD value during optimization: %g.\n",
                        std::sqrt(chi2));
                fprintf(fplog, "Parameters after the optimization:\n");
                fprintf(fplog, "%-5s  %10s  %10s  %10s\n", "Index",
                        "Average", "Std. Dev.", "Optimum");
                for (size_t k = 0; k < param_.size(); k++)
                {
                    fprintf(fplog, "%5zu  %10g  %10g  %10g\n",
                            k, pmean[k], psigma[k], best[k]);
                }
            }
            printEnergies(fplog);
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        int niter = gmx_recv_int(commrec(), 0);
        for (int n = 0; n < niter; n++)
        {
            calcDeviation();
        }
    }
    return bMinimum;
}

static void print_mols(FILE                          *fp,
                       std::vector<alexandria::MyMol> mol,
                       gmx_bool                       bForce,
                       gmx_bool                       bMtop)
{
    int j, k;
    for (auto mi = mol.begin(); (mi < mol.end()); mi++)
    {
        int nSP = 0, nOpt = 0;
        for (auto ei = mi->molProp()->BeginExperiment(); ei < mi->molProp()->EndExperiment(); ++ei)
        {
            auto jtype = ei->getJobtype();
            if (jtype == JOB_SP)
            {
                nSP += 1;
            }
            else if (jtype == JOB_OPT)
            {
                nOpt += 1;
            }
        }
        if (nOpt != 1)
        {
            fprintf(stderr, "Number of optimized conformations is %d. Check your input.\n", nOpt);
        }
        fprintf(fp, "%s natoms: %d Opt conformations: %d SP conformations: %d\n",
                mi->molProp()->getMolname().c_str(),
                mi->topology_->atoms.nr,
                nOpt,
                nSP);
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
        
        if (bMtop)
        {
            pr_mtop(fp, 0, mi->molProp()->getMolname().c_str(), mi->mtop_, true, false);
        }
    }
    if (bForce)
    {
        for (auto mi = mol.begin(); (mi < mol.end()); mi++)
        {
            fprintf(fp, "%-20s", mi->molProp()->getMolname().c_str());
            for (k = 0; k < F_NRE; k++)
            {
                if ((mi->enerd_->term[k] != 0) ||
                    (mi->mtop_->moltype[0].ilist[k].size() > 0))
                {
                    fprintf(fp, " %s: %.2f", interaction_function[k].name,
                            mi->enerd_->term[k]);
                }
            }
            fprintf(fp, "\n");
        }
    }
}

void Optimization::printResults(FILE                   *fp,
                                const char             *title,
                                const char             *hform_xvg,
                                const char             *HF_xvg,
                                const gmx_output_env_t *oenv)
{
    FILE       *xfp, *hfp = nullptr;
    gmx_stats_t gst;

    gst = gmx_stats_init();
    {
        const char *title = "Enthalpy of Formation";
        const char *yaxis = "Alexandria (kJ/mol)";
        if (nullptr != hform_xvg)
        {
            xfp = xvgropen(hform_xvg, title, "Experiment (kJ/mol)", yaxis, oenv);
        }
        if (nullptr != HF_xvg)
        {
            hfp = xvgropen(HF_xvg, title,
                           "Experiment + \\f{12}DD\\f{4}E(B3LYP/aug-cc-pVTZ) (kJ/mol)",
                           yaxis, oenv);
            xvgr_view(hfp, 0.15, 0.15, 0.75, 0.85, oenv);
        }
    }
    fprintf(fp, "%s\n", title);
    fprintf(fp, "Fit of energy at different conformations to y = ax\n");
    fprintf(fp, "Nr.   %-30s %10s %10s %7s %7s %7s %7s %7s %4s\n",
            "Molecule", "DHf@298K", "Emol@0K", "rmsF@0K",
            "rms E", "MSE E", "a", "r2", "N");

    double msdTot        = 0;
    int    nconformation = 0;
    int    imol          = 0;
    for (auto mi = mymols().begin(); mi < mymols().end(); mi++, imol++)
    {
        real DeltaE = mi->OptEcalc_ - mi->Emol_;
        if (nullptr != hform_xvg)
        {
            fprintf(xfp, "%10g  %10g\n", mi->Hform_, mi->Hform_ + DeltaE);
        }
        gmx_stats_t gmol = gmx_stats_init();
        double   spHF    = 0;
        double   optHF   = 0;
        double   deltaEn = 0;
        gmx_bool bpolar  = (mi->shellfc_ != nullptr);
        GMX_RELEASE_ASSERT(mi->molProp()->getOptHF(&optHF), 
                           "No optimized energy");
        if (nullptr != hfp)
        {
            fprintf(hfp, "@ s%d legend \"%s\"\n", imol, 
                    mi->molProp()->getMolname().c_str());
            fprintf(hfp, "@type xy\n");
        }
        double rmsForce = 0;
        for (auto ei = mi->molProp()->BeginExperiment(); ei < mi->molProp()->EndExperiment(); ++ei)
        {
            auto jtype = ei->getJobtype();
            // Exclude other experimental data points!
            if (jtype == JOB_OPT || jtype == JOB_SP)
            {
                ei->getHF(&spHF);
                deltaEn = spHF - optHF;
                mi->changeCoordinate(ei, bpolar);
                mi->computeForces(nullptr, commrec());
                real Hexper = mi->Hform_ + deltaEn;
                real Halex  = mi->enerd_->term[F_EPOT] + mi->Hform_ - mi->Emol_;
                if (nullptr != hfp)
                {
                    fprintf(hfp, "%10g  %10g\n", Hexper, Halex);
                }
                gmx_stats_add_point(gmol, Hexper, Halex, 0, 0);
                if (jtype == JOB_OPT)
                {
                    double f2 = 0;
                    int    natoms = mi->mtop_->natoms;
                    for (int j = 0; j < natoms; j++)
                    {
                        f2 += iprod(mi->f_[j], mi->f_[j]);
                    }
                    rmsForce = std::sqrt(f2 / natoms);
                    gmx_stats_add_point(gst, Hexper, Halex, 0, 0);
                }
            }
        }
        if (nullptr != hfp)
        {
            fprintf(hfp, "&\n");
        }
        real a, chi2, da, rmsd, Rfit, Rdata, mse, mae;
        int  N;
        // Note we are ignoring the return value for these functions
        (void) gmx_stats_get_a(gmol, 0, &a, &da, &chi2, &Rfit);
        (void) gmx_stats_get_rmsd(gmol, &rmsd);
        (void) gmx_stats_get_mse_mae(gmol, &mse, &mae);
        (void) gmx_stats_get_npoints(gmol, &N);
        nconformation += N;
        msdTot += gmx::square(rmsd)*N;
        gmx_stats_get_corr_coeff(gmol, &Rdata);
        fprintf(fp, "%-5d %-30s %10g %10g %7.1f %7.3f %7.3f %7.3f %6.1f%% %4d\n",
                imol,
                mi->molProp()->getMolname().c_str(),
                mi->Hform_,
                mi->Emol_,
                rmsForce, rmsd, 
                mse, a, Rdata*100, N);
        gmx_stats_free(gmol);
    }
    fprintf(fp, "RMSD from target energies for %zu compounds and %d conformation is %g\n",
            mymols().size(), nconformation,
            std::sqrt(msdTot/nconformation)); 
    fprintf(fp, "\n");
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

    real a, da, chi2, rmsd, Rfit, Rdata;
    int  N;
    gmx_stats_get_a(gst, 1, &a, &da, &chi2, &Rfit);
    gmx_stats_get_npoints(gst, &N);
    gmx_stats_get_corr_coeff(gst, &Rdata);
    gmx_stats_get_rmsd(gst, &rmsd);
    fprintf(fp, "Regression analysis fit to y = ax of %d minimum energy structures:\n", N);
    fprintf(fp, "a = %.3f  r2 = %.1f%%  rmsd = %.2f kJ/mol\n",
            a, Rdata*100, rmsd);
    gmx_stats_free(gst);
    fflush(fp);
}

} // namespace alexandria

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

    int                   reinit        = 0;
    int                   compress      = 0;
    int                   nmultisim     = 0;
    gmx_bool              bZPE          = false;
    gmx_bool              bZero         = true;
    t_pargs               pa[]         = {
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do optimization in multiple simulation" },
        { "-reinit",  FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A value of 0 means this is never done at all." },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule. If set, the zero point energy will be subtracted from the target energy when optimizing the force field model." },
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
             nullptr,
             gms,
             false,
             opt.iOpt(eitLJ14),
             opt.iOpt(eitPROPER_DIHEDRALS),
             bZPE,
             false,
             tabfn);


    opt.checkSupport(fp);

    if (nullptr != fp)
    {
        fprintf(fp, "In the total data set of %zu molecules we have:\n", opt.mymols().size());
    }

    if (MASTER(opt.commrec()))
    {
        opt.InitOpt(fp);
        print_mols(fp, opt.mymols(), false, false);
    }
    
    opt.broadcast();
    opt.calcDeviation();

    if (MASTER(opt.commrec()))
    {
        opt.printResults(fp, (char *)"Before optimization", 
                         nullptr, nullptr, oenv);
    }

    // Optimize the parameters to minimize the penalty function.
    bool bMinimum = opt.optRun(MASTER(opt.commrec()) ? stderr : nullptr,
                               fp,
                               oenv,
                               opt2fn("-conv", NFILE, fnm),
                               opt2fn("-epot", NFILE, fnm));

    if (MASTER(opt.commrec()))
    {
        if (bMinimum)
        {
            // Copy the parameters from the tune_fc optimization
            // to the poldata array.
            // std::vector<bool> changed(opt.nParams(), true);
            // opt.tuneFc2PolData(changed);
            // Now print the output.
            print_mols(fp, opt.mymols(), true, false);
            opt.printResults(fp, (char *)"After optimization",
                             opt2fn("-x", NFILE, fnm),
                             opt2fn("-hf", NFILE, fnm),
                             oenv);
            writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), compress);
        }
        else
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
        gmx_ffclose(fp);
    }

    return 0;
}
