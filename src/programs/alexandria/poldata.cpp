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
#include "gmxpre.h"

#include "poldata.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
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
                       const std::string &vdwparams,
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
                           { return fa.getType().compare(atype) == 0; });
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
                                           { return (bosque.compare(b.getBosque()) == 0); });
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
                                           { return (miller.compare(m.getMiller()) == 0); });
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
                           { return fa.getType().compare(atype) == 0; });
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

void Poldata::broadcast(t_commrec *cr)
{
    unsigned int          i, j, nep;
    std::vector<Eemprops> ep;

    if (NULL != debug)
    {
        fprintf(debug, "Going to update poldata on node %d\n", cr->nodeid);
    }
    fprintf(stderr, "Fixme: Poldata::broadcast is broken\n");
    if (MASTER(cr))
    {
        for (i = 1; ((int)i < cr->nnodes); i++)
        {
            gmx_send_int(cr, i, eep_.size());
            gmx_send(cr, i, eep_.data(), eep_.size()*sizeof(eep_[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != eep_.size())
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        eep_.clear();
        gmx_recv(cr, 0, ep.data(), eep_.size()*sizeof(ep[0]));
        for (i = 0; (i < nep); i++)
        {
            eep_.push_back(ep[i]);
        }
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    getEemtypeName(eep_[i].getEqdModel()),
                    eep_[i].getName(), eep_[i].getChi0(),
                    eep_[i].getJ0());
            for (j = 0; ((int)j < eep_[i].getNzeta()); j++)
            {
                fprintf(debug, " %8.3f", eep_[i].getZeta(j));
            }
            fprintf(debug, "\n");
        }
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
                        { return (strcasecmp(eep.getName(), nn.c_str()) == 0 &&
                                  (eep.getEqdModel() == eqdModel)); });
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
                        { return (strcasecmp(eep.getName(), nn.c_str()) == 0 &&
                                  (eep.getEqdModel() == eqdModel)); });
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
