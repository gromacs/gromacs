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
 
#include "tune_fc_utils.h"

#include "gromacs/topology/topology.h"
#include "gromacs/utility/strconvert.h"

#include "gmx_simple_comm.h"

namespace alexandria
{

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
                if (fat != pd.getAtypeEnd() && !fat->fixed())
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
                int         index = 0;
                std::string buf;
                std::string buf_reverse;
                auto        iType = static_cast<InteractionType>(bt_);
                switch (iType)
                {
                    case eitBONDS:
                    {
                        atoms   = {aai, aaj};
                        auto fs = pd.findForces(iType);
                        auto f  = fs->findForce(atoms);

                        if (fs->forceEnd() != f && !f->fixed())
                        {
                            buf         = gmx::formatString("%s %s", aai.c_str(), aaj.c_str());
                            buf_reverse = gmx::formatString("%s %s", aaj.c_str(), aai.c_str());
                            params      = f->params();
                            index       = f - fs->forceBegin();
                            found       = true;
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

                            if (fs->forceEnd() != f && !f->fixed())
                            {
                                buf         = gmx::formatString("%s %s %s", aai.c_str(),
                                                                aaj.c_str(), aak.c_str());
                                buf_reverse = gmx::formatString("%s %s %s", aak.c_str(),
                                                                aaj.c_str(), aai.c_str());
                                params      = f->params();
                                index       = f - fs->forceBegin();
                                found       = true;
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

                            if (fs->forceEnd() != f && !f->fixed())
                            {
                                buf          = gmx::formatString("%s %s %s %s", aai.c_str(),
                                                                 aaj.c_str(), aak.c_str(), aal.c_str());
                                buf_reverse  = gmx::formatString("%s %s %s %s", aal.c_str(),
                                                                 aaj.c_str(), aak.c_str(), aai.c_str());
                                params       = f->params();
                                index        = f - fs->forceBegin();
                                found        = true;
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

} // namespace alexandria
