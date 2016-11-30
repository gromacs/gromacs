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
/*! \internal \briefassignstr
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"

#include "gmx_simple_comm.h"
#include "poldata.h"
#include "stringutil.h"

namespace alexandria
{

void Poldata::setFilename(const std::string &fn2)
{
    GMX_RELEASE_ASSERT((fn2.size() > 0),
                       "Trying to set empty Poldata filename");

    if (0 != filename_.size())
    {
        fprintf(stderr, "Changing Poldata filename_ from %s to %s\n",
                filename_.c_str(), fn2.c_str());

    }
    filename_ = fn2;
}

void Poldata::addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       std::string       &vdwparams,
                       const std::string &refEnthalpy)
{

    unsigned int        i;

    for (i = 0; (i < alexandria_.size()); i++)
    {
        if (alexandria_[i].getType().compare(atype) == 0)
        {
            break;
        }
    }
    if (i == alexandria_.size())
    {
        Ffatype sp(desc, atype, ptype, btype,
                   elem, vdwparams, refEnthalpy);

        alexandria_.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Atom type %s was already added to Poldata record\n", atype.c_str());
    }
}

gmx_bool Poldata::strcasestrStart(std::string needle, std::string haystack)
{
    std::string ptr;

    ptr = strcasestr(haystack.c_str(), needle.c_str());

    return (ptr == haystack);
}

bool Poldata::getAtypeRefEnthalpy(const std::string &atype,
                                  double            *Href) const
{
    auto fa = findAtype(atype);
    if (alexandria_.end() != fa)
    {
        *Href = atof(fa->getRefEnthalpy().c_str());
        return true;
    }
    return false;
}

std::string Poldata::getDesc( std::string atype)
{
    unsigned int i;

    if (atype.size() != 0)
    {
        for (i = 0; (i < alexandria_.size()); i++)
        {
            if (alexandria_[i].getType().compare(atype) == 0)
            {
                return alexandria_[i].getDesc();
            }
        }
    }

    return "";
}

/*
 *-+-+-+-+-+-+-+-+-+-+-+
 * Polarizability STUFF
 *-+-+-+-+-+-+-+-+-+-+-+
 */

bool Poldata::atypeToPtype(const std::string &atype,
                           std::string       &ptype) const
{
    if (atype.size() == 0)
    {
        return false;
    }
    auto ai = std::find_if(alexandria_.begin(), alexandria_.end(),
                           [atype](Ffatype const &fa)
        {
            return fa.getType().compare(atype) == 0;
        });
    if (ai != alexandria_.end() && ai->getPtype().size() > 0)
    {
        ptype = ai->getPtype();
        return true;
    }
    return false;
}

bool Poldata::getAtypePol(const std::string &atype,
                          double            *polar,
                          double            *sigPol) const
{
    auto fa = findAtype(atype);
    if (alexandria_.end() != fa)
    {
        return getPtypePol(fa->getPtype(), polar, sigPol);
    }
    return false;
}

bool Poldata::getPtypePol(const std::string &ptype,
                          double            *polar,
                          double            *sigPol) const
{
    unsigned int j;

    for (j = 0; (j < ptype_.size()); j++)
    {
        if (ptype.compare(ptype_[j].getType()) == 0)
        {
            *polar   = ptype_[j].getPolarizability();
            *sigPol  = ptype_[j].getSigPol();

            return true;
        }
    }
    return false;
}

bool Poldata::setPtypePolarizability(const std::string &ptype,
                                     double             polarizability,
                                     double             sigPol)
{
    PtypeIterator sp = findPtype(ptype);

    if (ptype_.end() != sp)
    {
        sp->setPolarizability(polarizability);
        sp->setSigPol(sigPol);

        return true;
    }
    return false;
}

void Poldata::addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double             polarizability,
                       double             sigPol)
{

    unsigned int      i;

    for (i = 0; (i < ptype_.size()); i++)
    {
        if (ptype_[i].getType().compare(ptype) == 0)
        {
            break;
        }
    }
    if (i == ptype_.size())
    {
        Ptype sp(ptype, miller, bosque, polarizability, sigPol);
        ptype_.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Polarizability type %s was already added to Poldata record\n", ptype.c_str());
    }
}

/*
 *-+-+-+-+-+-+-+-+
 * Bosque STUFF
 *-+-+-+-+-+-+-+-+
 */

bool Poldata::getBosquePol(const std::string &bosque,
                           double            *polarizability) const
{
    BosqueConstIterator bb  = bosque_.begin(), be = bosque_.end();
    BosqueConstIterator bos = std::find_if(bb, be, [bosque](Bosque const &b)
        {
            return (bosque.compare(b.getBosque()) == 0);
        });
    if (bosque_.end() != bos)
    {
        *polarizability = bos->getPolarizability();

        return true;
    }
    return false;
}

bool Poldata::ptypeToBosque(const std::string &ptype,
                            std::string       &bosque) const
{
    for (const auto &i : ptype_)
    {
        if (ptype.compare(i.getType()) == 0)
        {
            bosque = i.getBosque();
            return true;
        }
    }
    return false;
}

/*
 *-+-+-+-+-+-+-+-+
 * Miller STUFF
 *-+-+-+-+-+-+-+-+
 */

bool Poldata::ptypeToMiller(const std::string &ptype,
                            std::string       &miller) const
{
    for (const auto &i : ptype_)
    {
        if (ptype.compare(i.getType()) == 0)
        {
            miller = i.getMiller();
            return true;
        }
    }
    return false;
}

void Poldata::addMiller(const std::string &miller,
                        int                atomnumber,
                        double             tauAhc,
                        double             alphaAhp,
                        const std::string &alexandria_equiv)
{
    Miller mil(miller, atomnumber, tauAhc, alphaAhp, alexandria_equiv);

    miller_.push_back(mil);
}


bool Poldata::getMillerPol(const std::string &miller,
                           int               *atomnumber,
                           double            *tauAhc,
                           double            *alphaAhp,
                           std::string       &alexandria_equiv) const
{
    MillerConstIterator mb  = miller_.begin(), me = miller_.end();
    MillerConstIterator mil = std::find_if(mb, me, [miller](Miller const &m)
        {
            return (miller.compare(m.getMiller()) == 0);
        });
    if (miller_.end() != mil)
    {
        *atomnumber      = mil->getAtomnumber();
        *tauAhc          = mil->getTauAhc();
        *alphaAhp        = mil->getAlphaAhp();
        alexandria_equiv = mil->getAlexandriaEquiv();

        return true;
    }
    return false;
}

/*
 *-+-+-+-+-+-+-+-+-+-+
 * LISTED FORCES
 *-+-+-+-+-+-+-+-+-+-+
 */

bool Poldata::atypeToBtype(const std::string &atype,
                           std::string       &btype) const
{
    if (atype.size() == 0)
    {
        return false;
    }
    auto ai = std::find_if(alexandria_.begin(), alexandria_.end(),
                           [atype](Ffatype const &fa)
        {
            return fa.getType().compare(atype) == 0;
        });
    if (ai != alexandria_.end() && ai->getBtype().size() > 0)
    {
        btype = ai->getBtype();
        return true;
    }
    return false;
}

void Poldata::addBtype(const std::string &btype)
{
    size_t i;

    for (i = 0; (i < btype_.size()); i++)
    {
        if (btype.compare(btype_[i]) == 0)
        {
            break;
        }
    }
    if (i == btype_.size())
    {
        btype_.push_back(btype);
    }
}

bool Poldata::findForce(std::vector<std::string> &atoms,
                        ListedForceIterator      *force)
{
    for (auto &f : forces_)
    {
        auto tmp = f.findForce(atoms);
        if (f.forceEnd() != tmp)
        {
            *force = tmp;
            return true;
        }
    }
    return false;
}

bool Poldata::findForce(const std::vector<std::string> &atoms,
                        ListedForceConstIterator       *force) const
{
    for (const auto &f : forces_)
    {
        auto tmp = f.findForce(atoms);
        if (f.forceEnd() != tmp)
        {
            *force = tmp;
            return true;
        }
    }
    return false;
}

bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain) const
{
    for (auto &f : forces_)
    {
        if (f.searchForce(atoms, params, refValue,
                          sigma, ntrain))
        {
            return true;
        }
    }
    return false;
}

bool Poldata::searchForce(std::vector<std::string> &atoms,
                          std::string              &params,
                          double                   *refValue,
                          double                   *sigma,
                          size_t                   *ntrain,
                          InteractionType           iType) const
{
    auto f  = findForces(iType);

    if (forcesEnd() != f)
    {
        if (f->searchForce(atoms, params, refValue,
                           sigma, ntrain))
        {
            return true;
        }
    }
    return false;
}


/*
 *-+-+-+-+-+-+-+-+
 * VDW STUFF
 *-+-+-+-+-+-+-+-+
 */

void Poldata::setVdwFunction(const std::string &func)
{
    unsigned int i;

    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, func.c_str()) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Van der Waals function '%s' does not exist in gromacs", func.c_str());
    }
    gtVdwFtype_    = i;
    gtVdwFunction_ = func;
}

void Poldata::setCombinationRule(const std::string &func)
{
    unsigned int i;

    for (i = 0; (i < eCOMB_NR); i++)
    {
        if (strcasecmp(ecomb_names[i], func.c_str()) == 0)
        {
            break;
        }
    }
    if (i == eCOMB_NR)
    {
        gmx_fatal(FARGS, "Combination rule '%s' does not exist in gromacs", func.c_str());
    }
    gtCombRule_        = i;
    gtCombinationRule_ = func;
}

/*
 *-+-+-+-+-+-+-+-+
 * COULOMB STUFF
 *-+-+-+-+-+-+-+-+
 */

void Poldata::addSymcharges(const std::string &central,
                            const std::string &attached,
                            int                numattach)
{
    unsigned int           i;
    Symcharges           * sc;
    for (i = 0; (i < symcharges_.size()); i++)
    {
        sc = &(symcharges_[i]);
        if ((strcasecmp(sc->getCentral().c_str(), central.c_str()) == 0)   &&
            (strcasecmp(sc->getAttached().c_str(), attached.c_str()) == 0) &&
            (sc->getNumattach() == numattach))
        {
            break;
        }
    }
    if (i == symcharges_.size())
    {
        Symcharges symcharges(central, attached, numattach);

        symcharges_.push_back(symcharges);
    }
}

int Poldata::getNumprops(ChargeDistributionModel eqdModel) const
{
    unsigned int i, n = 0;

    for (i = 0; (i < eep_.size()); i++)
    {
        if (eep_[i].getEqdModel() == eqdModel)
        {
            n++;
        }
    }

    return n;
}

int Poldata::havePolSupport(const std::string &atype) const
{
    unsigned int i;

    for (i = 0; (i < alexandria_.size()); i++)
    {
        if (atype.compare(alexandria_[i].getType()) == 0)
        {
            return 1;
        }
    }
    return 0;
}

bool Poldata::haveEemSupport(ChargeDistributionModel  eqdModel,
                             const std::string       &name,
                             gmx_bool                 bAllowZeroParameters) const
{
    EempropsConstIterator eep = findEem(eqdModel, name);

    return (eep != EndEemprops() &&
            (bAllowZeroParameters || ((eep->getJ0() > 0) && (eep->getChi0() > 0))));
}

double Poldata::getJ00(ChargeDistributionModel  eqdModel,
                       const std::string       &name) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        return eer->getJ0();
    }
    return -1;
}

const char *Poldata::getQstr(ChargeDistributionModel  eqdModel,
                             const std::string       &name) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        return eer->getQstr();
    }
    return nullptr;
}

const char *Poldata::getRowstr(ChargeDistributionModel  eqdModel,
                               const std::string       &name) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        return eer->getRowstr();
    }
    return nullptr;
}

int Poldata::getRow(ChargeDistributionModel  eqdModel,
                    const std::string       &name,
                    int                      zz) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getRow(zz);
    }
    return -1;
}

double Poldata::getZeta(ChargeDistributionModel  eqdModel,
                        const std::string       &name,
                        int                      zz) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        if ((zz < 0) || (zz >= eer->getNzeta()))
        {
            printf("Bleh\n");
        }
        range_check(zz, 0, eer->getNzeta());
        return eer->getZeta(zz);
    }
    return -1;
}

int Poldata::getNzeta(ChargeDistributionModel  eqdModel,
                      const std::string       &name) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        return eer->getNzeta();
    }
    return 0;
}

double Poldata::getQ(ChargeDistributionModel  eqdModel,
                     const std::string       &name,
                     int                      zz) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getQ(zz);
    }
    return -1;
}

double Poldata::getChi0(ChargeDistributionModel   eqdModel,
                        const  std::string       &name) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        return eer->getChi0();
    }
    else
    {
        gmx_fatal(FARGS, "No chi0 data for eqdModel %d and name %s", eqdModel, name.c_str());
    }
    return -1;
}

void Poldata::setEpref(ChargeDistributionModel  eqdModel,
                       const std::string       &epref)
{
    for (auto &i : epr_)
    {
        if (i.getEqdModel() == eqdModel)
        {
            i.setEpref(epref);
            return;
        }
    }
    Epref epr(eqdModel, epref);
    epr_.push_back(epr);
}

const char *Poldata::getEpref(ChargeDistributionModel eqdModel) const
{
    for (auto &i : epr_)
    {
        if (i.getEqdModel() == eqdModel)
        {
            return i.getEpref();
        }
    }
    return nullptr;
}

CommunicationStatus Poldata::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, filename_.c_str());
        gmx_send_str(cr, dest, alexandriaPolarUnit_.c_str());
        gmx_send_str(cr, dest, alexandriaPolarRef_.c_str());
        gmx_send_str(cr, dest, alexandriaForcefield_.c_str());
        gmx_send_int(cr, dest, nexcl_);
        gmx_send_double(cr, dest, fudgeQQ_);
        gmx_send_double(cr, dest, fudgeLJ_);
        gmx_send_str(cr, dest, gtVdwFunction_.c_str());
        gmx_send_str(cr, dest, gtCombinationRule_.c_str());
        gmx_send_int(cr, dest, gtVdwFtype_);
        gmx_send_int(cr, dest, gtCombRule_);
        gmx_send_str(cr, dest, millerTauUnit_.c_str());
        gmx_send_str(cr, dest, millerAhpUnit_.c_str());
        gmx_send_str(cr, dest, millerRef_.c_str());
        gmx_send_str(cr, dest, bosquePolarUnit_.c_str());
        gmx_send_str(cr, dest, bosqueRef_.c_str());       
        gmx_send_int(cr, dest, ptype_.size());
        gmx_send_int(cr, dest, alexandria_.size());
        gmx_send_int(cr, dest, btype_.size());
        gmx_send_int(cr, dest, forces_.size());
        gmx_send_int(cr, dest, miller_.size());
        gmx_send_int(cr, dest, bosque_.size());
        gmx_send_int(cr, dest, symcharges_.size());
        gmx_send_int(cr, dest, eep_.size());
        gmx_send_int(cr, dest, epr_.size());
        
        /*Send ptype*/
        for (auto &ptype : ptype_)
        {
            cs = ptype.Send(cr, dest);
        }
        
        /*Send Ffatype*/
        for (auto &alexandria : alexandria_)
        {
            cs = alexandria.Send(cr, dest);           
        }
        
        /*Send btype*/
        for (auto &btype : btype_)
        {
            gmx_send_str(cr, dest, btype.c_str());
        }
        
        /*Send Listed Forces*/
        for (auto &force : forces_)
        {
            cs = force.Send(cr, dest);
        }
        
        /*Send Miller*/
        for (auto &miller : miller_)
        {
            cs = miller.Send(cr, dest);
        }
        
        /*Send Bosque*/
        for (auto &bosque : bosque_)
        {
            cs = bosque.Send(cr, dest);
        }
        
        /*Send Symcharge*/
        for (auto &symcharges : symcharges_)
        {
            cs = symcharges.Send(cr, dest);
        }
        
        /*Send Eemprops*/
        for (auto &eep : eep_)
        {
            cs = eep.Send(cr, dest);
        }
        
        /*Send Epref*/
        for (auto &epr : epr_)
        {
            cs = epr.Send(cr, dest);
        }
    }
    return cs;
}
        
CommunicationStatus Poldata::Receive(t_commrec *cr, int src)
{
    size_t nptype, nalexandria, nbtype, nforces;
    size_t nmiller, nbosque, nsymcharges, neep, nepr;
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        filename_             = gmx_recv_str(cr, src);
        alexandriaPolarUnit_  = gmx_recv_str(cr, src);
        alexandriaPolarRef_   = gmx_recv_str(cr, src);
        alexandriaForcefield_ = gmx_recv_str(cr, src);
        nexcl_                = gmx_recv_int(cr, src);
        fudgeQQ_              = gmx_recv_double(cr, src);
        fudgeLJ_              = gmx_recv_double(cr, src);
        gtVdwFunction_        = gmx_recv_str(cr, src);
        gtCombinationRule_    = gmx_recv_str(cr, src);
        gtVdwFtype_           = gmx_recv_int(cr, src);
        gtCombRule_           = gmx_recv_int(cr, src);
        millerTauUnit_        = gmx_recv_str(cr, src);
        millerAhpUnit_        = gmx_recv_str(cr, src);
        millerRef_            = gmx_recv_str(cr, src);
        bosquePolarUnit_      = gmx_recv_str(cr, src);
        bosqueRef_            = gmx_recv_str(cr, src);     
        nptype                = gmx_recv_int(cr, src);
        nalexandria           = gmx_recv_int(cr, src);
        nbtype                = gmx_recv_int(cr, src);
        nforces               = gmx_recv_int(cr, src);
        nmiller               = gmx_recv_int(cr, src);
        nbosque               = gmx_recv_int(cr, src);
        nsymcharges           = gmx_recv_int(cr, src);
        neep                  = gmx_recv_int(cr, src);
        nepr                  = gmx_recv_int(cr, src);
        
        
        /*Receive ptype*/
        ptype_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nptype); n++)
        {
            Ptype ptype;
            cs = ptype.Receive(cr, src);
            if (CS_OK == cs)
            {
                ptype_.push_back(ptype);
            }
        }
        
        /*Receive Ffatype*/
        alexandria_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nalexandria); n++)
        {
            Ffatype alexandria;
            cs = alexandria.Receive(cr, src);
            if (CS_OK == cs)
            {
                alexandria_.push_back(alexandria);
            }
        }
        
        /*Receive btype*/
        btype_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nbtype); n++)
        {
            char *btype = gmx_recv_str(cr, src);
            if (nullptr == btype)
            {
                btype_.push_back(btype);
            }
        }
        
        /*Receive Listed Forces*/
        forces_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nforces); n++)
        {
            ListedForces fs;
            cs = fs.Receive(cr, src);
            if (CS_OK == cs)
            {
                forces_.push_back(fs);
            }
        }
        
        /*Receive Miller*/
        miller_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nmiller); n++)
        {
            Miller miller;
            cs = miller.Receive(cr, src);
            if (CS_OK == cs)
            {
                miller_.push_back(miller);
            }
        }

        /*Receive Bosque*/
        bosque_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nbosque); n++)
        {
            Bosque bosque;
            cs = bosque.Receive(cr, src);
            if (CS_OK == cs)
            {
                bosque_.push_back(bosque);
            }
        }
        
        /*Receive Symcharges*/
        symcharges_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nsymcharges); n++)
        {
            Symcharges symcharges;
            cs = symcharges.Receive(cr, src);
            if (CS_OK == cs)
            {
                symcharges_.push_back(symcharges);
            }
        }
        
        /*Receive Eemprops*/
        eep_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < neep); n++)
        {
            Eemprops eep;
            cs = eep.Receive(cr, src);
            if (CS_OK == cs)
            {
                eep_.push_back(eep);
            }
        }
        
        /*Receive Epref*/
        epr_.clear();
        for (size_t n = 0; (CS_OK == cs) && (n < nepr); n++)
        {
            Epref epr;
            cs = epr.Receive(cr, src);
            if (CS_OK == cs)
            {
                epr_.push_back(epr);
            }
        }     
    }
    return cs;
}

void Poldata::broadcast(t_commrec *cr)
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
            fprintf(debug, "Going to update poldata on node %d\n", cr->nodeid);
        }
        Receive(cr, src);
    }
}

EempropsConstIterator Poldata::findEem(ChargeDistributionModel  eqdModel,
                                       const std::string       &name) const
{
    std::string nn = name;
    if (eqdModel == eqdRappe || eqdModel == eqdBultinck || eqdModel == eqdYang)
    {
        FfatypeConstIterator fa = findAtype(name);
        if (fa != getAtypeEnd())
        {
            nn = fa->getElem();
        }
        else
        {
            return eep_.end();
        }
    }
    return std::find_if(eep_.begin(), eep_.end(),
                        [eqdModel, nn](Eemprops const &eep)
        {
            return (strcasecmp(eep.getName(), nn.c_str()) == 0 &&
                    (eep.getEqdModel() == eqdModel));
        });
}

EempropsIterator Poldata::findEem(ChargeDistributionModel  eqdModel,
                                  const std::string       &name)
{
    std::string nn = name;
    if (eqdModel == eqdRappe || eqdModel == eqdBultinck || eqdModel == eqdYang)
    {
        FfatypeConstIterator fa = findAtype(name);
        if (fa != getAtypeEnd())
        {
            nn = fa->getElem();
        }
        else
        {
            return eep_.end();
        }
    }
    return std::find_if(eep_.begin(), eep_.end(),
                        [eqdModel, nn](Eemprops const &eep)
        {
            return (strcasecmp(eep.getName(), nn.c_str()) == 0 &&
                    (eep.getEqdModel() == eqdModel));
        });
}

typedef struct {
    ChargeDistributionModel eqd;
    const char             *name, *ref;
    gmx_bool                bWeight;
} t_eemtype_props;

t_eemtype_props eemtype_props[eqdNR] = {
    { eqdAXp,      "AXp",      "Maaren2016a",   FALSE },
    { eqdAXg,      "AXg",      "Maaren2016a",   TRUE },
    { eqdAXs,      "AXs",      "Maaren2016a",   TRUE },
    { eqdYang,     "Yang",     "Yang2006b",     TRUE },
    { eqdBultinck, "Bultinck", "Bultinck2002a", FALSE },
    { eqdRappe,    "Rappe",    "Rappe1991a",    TRUE }
};

ChargeDistributionModel name2eemtype(const std::string name)
{
    for (int i = 0; (i < eqdNR); i++)
    {
        if (strcasecmp(name.c_str(), eemtype_props[i].name) == 0)
        {
            return eemtype_props[i].eqd;
        }
    }
    return eqdNR;
}

const char *getEemtypeName(ChargeDistributionModel eem)
{
    for (int i = 0; (i < eqdNR); i++)
    {
        if (eem == eemtype_props[i].eqd)
        {
            return eemtype_props[i].name;
        }
    }
    return nullptr;
}
} // namespace alexandria
