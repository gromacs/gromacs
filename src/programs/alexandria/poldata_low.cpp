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

#include "poldata_low.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
#include "plistwrapper.h"

namespace alexandria
{

static const char * eit_names[eitNR] = {
    "BONDS",
    "ANGLES",
    "LINEAR_ANGLES",
    "PROPER_DIHEDRALS",
    "IMPROPER_DIHEDRALS",
    "VDW",
    "LJ14",
    "POLARIZATION",
    "CONSTR",
    "VSITE2",
    "VSITE3FAD",
    "VSITE3OUT"
};

static const char * evt_names[evtNR] = {
    "linear",
    "planar",
    "ring_planar",
    "in_plane",
    "out_of_plane",
    "all"
};

const char *iType2string(InteractionType iType)

{
    if (iType < eitNR)
    {
        return eit_names[iType];
    }
    return nullptr;
}

InteractionType string2iType(const char *string)
{
    int i;
    for (i = 0; i < eitNR; i++)
    {
        if (gmx_strcasecmp(string, eit_names[i]) == 0)
        {
            return static_cast<InteractionType>(i);
        }
    }
    return eitNR;
}

const char *vsiteType2string(VsiteType vType)
{
    if (vType < evtNR)
    {
        return evt_names[vType];
    }
    return nullptr;
}

VsiteType string2vsiteType(const char *string)
{
    int i;
    for (i = 0; i < evtNR; i++)
    {
        if (gmx_strcasecmp(string, evt_names[i]) == 0)
        {
            return static_cast<VsiteType>(i);
        }
    }
    return evtNR;
}

Ffatype::Ffatype(const std::string &desc,
                 const std::string &type,
                 const std::string &ptype,
                 const std::string &btype,
                 const std::string &ztype,
                 const std::string &elem,
                 std::string       &vdwparams,
                 const std::string &refEnthalpy)
    :
      desc_(desc),
      type_(type),
      ptype_(ptype),
      btype_(btype),
      ztype_(ztype),
      elem_(elem),
      vdwparams_(vdwparams),
      refEnthalpy_(refEnthalpy)
{}

CommunicationStatus Ffatype::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &desc_);
        gmx_send_str(cr, dest, &type_);
        gmx_send_str(cr, dest, &ptype_);
        gmx_send_str(cr, dest, &btype_);
        gmx_send_str(cr, dest, &ztype_);
        gmx_send_str(cr, dest, &elem_);
        gmx_send_str(cr, dest, &vdwparams_);
        gmx_send_str(cr, dest, &refEnthalpy_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Fftype %s %s %s %s %s %s %s, status %s\n",
                    desc_.c_str(), type_.c_str(), ptype_.c_str(),
                    btype_.c_str(), elem_.c_str(), vdwparams_.c_str(),
                    refEnthalpy_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;

}

CommunicationStatus Ffatype::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &desc_);
        gmx_recv_str(cr, src, &type_);
        gmx_recv_str(cr, src, &ptype_);
        gmx_recv_str(cr, src, &btype_);
        gmx_recv_str(cr, src, &ztype_);
        gmx_recv_str(cr, src, &elem_);
        gmx_recv_str(cr, src, &vdwparams_);
        gmx_recv_str(cr, src, &refEnthalpy_);

        if (nullptr != debug)
        {
            fprintf(debug, "Received Fftype %s %s %s %s %s %s %s, status %s\n",
                    desc_.c_str(), type_.c_str(), ptype_.c_str(),
                    btype_.c_str(), elem_.c_str(), vdwparams_.c_str(),
                    refEnthalpy_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Ptype::Ptype(const std::string &ptype,
             const std::string &miller,
             const std::string &bosque,
             double             polarizability,
             double             sigPol)
    :
      type_(ptype),
      miller_(miller),
      bosque_(bosque),
      polarizability_(polarizability),
      sigPol_(sigPol)
{}

CommunicationStatus Ptype::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &type_);
        gmx_send_str(cr, dest, &miller_);
        gmx_send_str(cr, dest, &bosque_);
        gmx_send_double(cr, dest, polarizability_);
        gmx_send_double(cr, dest, sigPol_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Ptype %s %s %s %g %g, status %s\n",
                    type_.c_str(), miller_.c_str(),
                    bosque_.c_str(), polarizability_,
                    sigPol_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Ptype::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &type_);
        gmx_recv_str(cr, src, &miller_);
        gmx_recv_str(cr, src, &bosque_);
        polarizability_ = gmx_recv_double(cr, src);
        sigPol_         = gmx_recv_double(cr, src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Ptype %s %s %s %g %g, status %s\n",
                    type_.c_str(), miller_.c_str(),
                    bosque_.c_str(), polarizability_,
                    sigPol_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Vsite::Vsite(const std::string &atype,
             const std::string &type,
             int                number,
             double             distance,
             double             angle,
             int                ncontrolatoms)
    :
      atype_(atype),
      type_(string2vsiteType(type.c_str())),
      number_(number),
      distance_(distance),
      angle_(angle),
      ncontrolatoms_(ncontrolatoms)
{}

CommunicationStatus Vsite::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        std::string vtype;
        vtype.assign(vsiteType2string(type_));
        gmx_send_str(cr, dest, &atype_);
        gmx_send_str(cr, dest, &vtype);
        gmx_send_int(cr, dest, number_);
        gmx_send_double(cr, dest, distance_);
        gmx_send_double(cr, dest, angle_);
        gmx_send_int(cr, dest, ncontrolatoms_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Vsite %s %s %d %g %g %d, status %s\n",
                    atype_.c_str(), vsiteType2string(type_), number_,
                    distance_, angle_, ncontrolatoms_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Vsite::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        std::string type;
        gmx_recv_str(cr, src, &atype_);
        gmx_recv_str(cr, src, &type);
        type_          = string2vsiteType(type.c_str());
        number_        = gmx_recv_int(cr, src);
        distance_      = gmx_recv_double(cr, src);
        angle_         = gmx_recv_double(cr, src);
        ncontrolatoms_ = gmx_recv_int(cr, src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Vsite %s %s %d %g %g %d, status %s\n",
                    atype_.c_str(), vsiteType2string(type_), number_,
                    distance_, angle_, ncontrolatoms_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

ListedForce::ListedForce(const std::vector<std::string> atoms,
                         std::string                    params,
                         double                         refValue,
                         double                         sigma,
                         size_t                         ntrain)
    :
      atoms_(atoms),
      params_(params),
      refValue_(refValue),
      sigma_(sigma),
      ntrain_(ntrain)

{}

CommunicationStatus ListedForce::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &params_);
        gmx_send_double(cr, dest, refValue_);
        gmx_send_double(cr, dest, sigma_);
        gmx_send_int(cr, dest, static_cast<int>(ntrain_));
        gmx_send_int(cr, dest, atoms_.size());

        for (auto &atom : atoms_)
        {
            gmx_send_str(cr, dest, &atom);
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Sent ListedForce %s %g %g %zu, status %s\n",
                    params_.c_str(), refValue_, sigma_,
                    ntrain_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus ListedForce::Receive(const t_commrec *cr, int src)
{
    int                 natom;
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &params_);
        refValue_ = gmx_recv_double(cr, src);
        sigma_    = gmx_recv_double(cr, src);
        ntrain_   = static_cast<size_t>(gmx_recv_double(cr, src));
        natom     = gmx_recv_int(cr, src);

        for (auto n = 0; n < natom; n++)
        {
            std::string atom;
            gmx_recv_str(cr, src, &atom);
            if (!atom.empty())
            {
                const_cast<std::vector<std::string> &>(atoms_).push_back(atom);
            }
            else
            {
                gmx_fatal(FARGS, "A category was promised but I got a nullptr pointer");
            }
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Received ListedForce %s %g %g %zu, status %s\n",
                    params_.c_str(), refValue_, sigma_,
                    ntrain_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

ListedForces::ListedForces(const std::string   iType,
                           const std::string  &function,
                           const std::string  &unit)
    :
      iType_(string2iType(iType.c_str())),
      function_(function),
      unit_(unit)
{
    size_t funcType;
    for (funcType = 0; funcType < F_NRE; funcType++)
    {
        if (strcasecmp(interaction_function[funcType].name, function_.c_str()) == 0)
        {
            break;
        }
    }
    if (funcType == F_NRE)
    {
        gmx_fatal(FARGS, "Force function '%s' does not exist in gromacs", function_.c_str());
    }

    fType_ = funcType;
}

CommunicationStatus ListedForces::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    std::string         itype;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        itype.assign(iType2string(iType_));
        gmx_send_str(cr, dest, &itype);
        gmx_send_str(cr, dest, &function_);
        gmx_send_str(cr, dest, &unit_);
        gmx_send_int(cr, dest, static_cast<int>(fType_));
        gmx_send_int(cr, dest, force_.size());

        for (auto &f : force_)
        {
            f.Send(cr, dest);
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Sent ListedForces %s %s %s %d, status %s\n",
                    iType2string(iType_), function_.c_str(), unit_.c_str(),
                    fType_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus ListedForces::Receive(const t_commrec *cr, int src)
{
    size_t              nforce;
    CommunicationStatus cs;
    std::string         iType, function, unit;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &iType);
        gmx_recv_str(cr, src, &function);
        gmx_recv_str(cr, src, &unit);
        iType_    = string2iType(iType.c_str());
        const_cast<std::string &>(function_) = function;
        const_cast<std::string &>(unit_)     = unit;
        fType_    = static_cast<unsigned int>(gmx_recv_int(cr, src));
        nforce    = gmx_recv_int(cr, src);

        force_.clear();
        for (size_t n = 0; n < nforce; n++)
        {
            ListedForce f;
            cs = f.Receive(cr, src);
            if (CS_OK == cs)
            {
                force_.push_back(f);
            }
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Received ListedForces %s %s %s %d, status %s\n",
                    iType2string(iType_), function_.c_str(), unit_.c_str(),
                    fType_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

ListedForceIterator ListedForces::findForce(const std::vector<std::string> &atoms)
{

    ListedForceIterator fb = forceBegin(), fe = forceEnd();
    return std::find_if(fb, fe, [atoms](const ListedForce &force)
                        {
                            std::vector<std::string> atoms_re(atoms.rbegin(), atoms.rend());
                            return (atoms == force.atoms() || atoms_re == force.atoms());
                        });
}

ListedForceConstIterator ListedForces::findForce(const std::vector<std::string> &atoms) const
{

    ListedForceConstIterator fb = forceBegin(), fe = forceEnd();
    return std::find_if(fb, fe, [atoms](const ListedForce &force)
                        {
                            std::vector<std::string> atoms_re(atoms.rbegin(), atoms.rend());
                            return (atoms == force.atoms() || atoms_re == force.atoms());
                        });
}

bool ListedForces::setForceParams(const std::vector<std::string> &atoms,
                                  const std::string              &params,
                                  double                          refValue,
                                  double                          sigma,
                                  size_t                          ntrain)
{
    auto force = findForce(atoms);
    if (forceEnd() != force)
    {
        force->setRefValue(refValue);
        force->setSigma(sigma);
        force->setNtrain(ntrain);
        force->setParams(params);

        return true;
    }
    return false;
}

void ListedForces::addForce(const std::vector<std::string> &atoms,
                            const std::string              &params,
                            double                          refValue,
                            double                          sigma,
                            size_t                          ntrain)
{
    if (setForceParams(atoms, params, refValue,
                       sigma, ntrain))
    {
        return;
    }
    ListedForce force(atoms, params, refValue,
                      sigma, ntrain);
    force_.push_back(force);
}

bool ListedForces::searchForce(std::vector<std::string> &atoms,
                               std::string              &params,
                               double                   *refValue,
                               double                   *sigma,
                               size_t                   *ntrain) const
{
    auto force = findForce(atoms);
    if (forceEnd() != force)
    {
        *refValue = force->refValue();
        *sigma    = force->sigma();
        *ntrain   = force->ntrain();
        params    = force->params();

        return true;
    }
    return false;
}

Bosque::Bosque(const std::string &bosque, double polarizability)
    :
      bosque_(bosque),
      polarizability_(polarizability)
{}

CommunicationStatus Bosque::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &bosque_);
        gmx_send_double(cr, dest, polarizability_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Bosque %s %g, status %s\n",
                    bosque_.c_str(), polarizability_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Bosque::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &bosque_);
        polarizability_ = gmx_recv_double(cr, src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Bosque %s %g, status %s\n",
                    bosque_.c_str(), polarizability_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Miller::Miller(const std::string &miller,
               int                atomnumber,
               double             tauAhc,
               double             alphaAhp,
               const std::string &alexandria_equiv)
    :
      miller_(miller),
      atomnumber_(atomnumber),
      tauAhc_(tauAhc),
      alphaAhp_(alphaAhp),
      alexandria_equiv_(alexandria_equiv)
{}

CommunicationStatus Miller::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &miller_);
        gmx_send_int(cr, dest, atomnumber_);
        gmx_send_double(cr, dest, tauAhc_);
        gmx_send_double(cr, dest, alphaAhp_);
        gmx_send_str(cr, dest, &alexandria_equiv_);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent Miller %s %d %g %g %s, status %s\n",
                    miller_.c_str(), atomnumber_, tauAhc_,
                    alphaAhp_, alexandria_equiv_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Miller::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &miller_);
        atomnumber_ = gmx_recv_int(cr, src);
        tauAhc_     = gmx_recv_double(cr, src);
        alphaAhp_   = gmx_recv_double(cr, src);
        gmx_recv_str(cr, src, &alexandria_equiv_);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Miller %s %d %g %g %s, status %s\n",
                    miller_.c_str(), atomnumber_, tauAhc_,
                    alphaAhp_, alexandria_equiv_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Symcharges::Symcharges(const std::string &central,
                       const std::string &attached,
                       int                numattach)
    :
      central_(central),
      attached_(attached),
      numattach_(numattach)
{}

CommunicationStatus Symcharges::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &central_);
        gmx_send_str(cr, dest, &attached_);
        gmx_send_int(cr, dest, numattach_);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent Symcharges %s %s %d, status %s\n",
                    central_.c_str(), attached_.c_str(), numattach_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Symcharges::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    std::string         central, attached;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &central);
        gmx_recv_str(cr, src, &attached);
        const_cast<std::string &>(central_)   = central;
        const_cast<std::string &>(attached_)  = attached;
        numattach_ = gmx_recv_int(cr, src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Symcharges %s %s %d, status %s\n",
                    central_.c_str(), attached_.c_str(), numattach_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}


Epref::Epref(ChargeDistributionModel  eqdModel,
             const std::string       &epref)
    :
      eqdModel_(eqdModel),
      epref_(epref)
{}

CommunicationStatus Epref::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    std::string         eqdModel;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        eqdModel.assign(getEemtypeName(eqdModel_));
        gmx_send_str(cr, dest, &eqdModel);
        gmx_send_str(cr, dest, &epref_);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent Epref %s %s, status %s\n",
                    getEemtypeName(eqdModel_), epref_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Epref::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    std::string         eqdModel;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &eqdModel);
        eqdModel_ = name2eemtype(eqdModel);
        gmx_recv_str(cr, src, &epref_);

        if (nullptr != debug)
        {
            fprintf(debug, "Received Epref %s %s, status %s\n",
                    getEemtypeName(eqdModel_), epref_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

RowZetaQ::RowZetaQ(int row, double zeta, double q)

    :

      row_(row),
      zeta_(zeta),
      q_(q),
      zetaRef_(zeta)


{
    zindex_ = -1;
    char buf[256];
    row_ = std::min(row_, SLATER_MAX);
    if (row_ < row && debug)
    {
        fprintf(debug, "Reducing row from %d to %d\n", row, row_);
    }
    snprintf(buf, sizeof(buf), "Row (%d) in the periodic table must be > 0 and <= %d",
             row_, SLATER_MAX);
    GMX_RELEASE_ASSERT(row_ > 0 && row_ <= SLATER_MAX, buf);
    fixedQ_ = (q != 0);
}

CommunicationStatus RowZetaQ::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, row_);
        gmx_send_double(cr, dest, zeta_);
        gmx_send_double(cr, dest, q_);
        gmx_send_double(cr, dest, zetaRef_);
        gmx_send_int(cr, dest, zindex_);
        gmx_send_int(cr, dest, fixedQ_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent RowZetaQ %d %g %g %g %d %d, status %s\n",
                    row_, zeta_, q_, zetaRef_, zindex_, fixedQ_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus RowZetaQ::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        row_     = gmx_recv_int(cr, src);
        zeta_    = gmx_recv_double(cr, src);
        q_       = gmx_recv_double(cr, src);
        zetaRef_ = gmx_recv_double(cr, src);
        zindex_  = gmx_recv_int(cr, src);
        fixedQ_  = gmx_recv_int(cr, src);

        if (nullptr != debug)
        {
            fprintf(debug, "Received RowZetaQ %d %g %g %g %d %d, status %s\n",
                    row_, zeta_, q_, zetaRef_, zindex_, fixedQ_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Eemprops::Eemprops(ChargeDistributionModel   eqdModel,
                   const std::string        &name,
                   const std::string        &rowstr,
                   const std::string        &zetastr,
                   const std::string        &zeta_sigma,
                   const std::string        &qstr,
                   double                    J0,
                   double                    J0_sigma,
                   double                    chi0,
                   double                    chi0_sigma)
    :
      eqdModel_(eqdModel),
      name_(name),
      rowstr_(rowstr),
      zetastr_(zetastr),
      zeta_sigma_(zeta_sigma),
      qstr_(qstr),
      J0_(J0),
      J0_sigma_(J0_sigma),
      chi0_(chi0),
      chi0_sigma_(chi0_sigma)
{
    setRowZetaQ(rowstr, zetastr, qstr);
}

CommunicationStatus Eemprops::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    std::string         eqdModel;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        eqdModel.assign(getEemtypeName(eqdModel_));
        gmx_send_str(cr, dest, &eqdModel);
        gmx_send_str(cr, dest, &name_);
        gmx_send_str(cr, dest, &rowstr_);
        gmx_send_str(cr, dest, &zetastr_);
        gmx_send_str(cr, dest, &zeta_sigma_);
        gmx_send_str(cr, dest, &qstr_);
        gmx_send_double(cr, dest, J0_);
        gmx_send_double(cr, dest, J0_sigma_);
        gmx_send_double(cr, dest, chi0_);
        gmx_send_double(cr, dest, chi0_sigma_);
        gmx_send_int(cr, dest, rzq_.size());

        for (auto &rzq : rzq_)
        {
            cs = rzq.Send(cr, dest);
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Eemprops %s %s %s %s %s %g %g, status %s\n",
                    getEemtypeName(eqdModel_), name_.c_str(), rowstr_.c_str(),
                    zetastr_.c_str(), qstr_.c_str(), J0_, chi0_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Eemprops::Receive(const t_commrec *cr, int src)
{
    size_t              nrzq;
    CommunicationStatus cs;
    std::string         eqdModel;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &eqdModel);
        eqdModel_   = name2eemtype(eqdModel);
        gmx_recv_str(cr, src, &name_);
        gmx_recv_str(cr, src, &rowstr_);
        gmx_recv_str(cr, src, &zetastr_);
        gmx_recv_str(cr, src, &zeta_sigma_);
        gmx_recv_str(cr, src, &qstr_);
        J0_         = gmx_recv_double(cr, src);
        J0_sigma_   = gmx_recv_double(cr, src);
        chi0_       = gmx_recv_double(cr, src);
        chi0_sigma_ = gmx_recv_double(cr, src);
        nrzq        = gmx_recv_int(cr, src);

        rzq_.clear();
        for (size_t n = 0; n < nrzq; n++)
        {
            RowZetaQ rzq;
            cs = rzq.Receive(cr, src);
            if (CS_OK == cs)
            {
                rzq_.push_back(rzq);
            }
        }
    }
    return cs;
}

void Eemprops::setRowZetaQ(const std::string &rowstr,
                           const std::string &zetastr,
                           const std::string &qstr)
{
    rzq_.clear();
    std::vector<std::string>  sz, sq, sr;
    sr = gmx::splitString(rowstr);
    sz = gmx::splitString(zetastr);
    sq = gmx::splitString(qstr);
    size_t nn = std::min(sz.size(), std::min(sq.size(), sr.size()));

    char   buf[256];
    snprintf(buf, sizeof(buf), "More zeta values than q/row values for %s n = %d\n",
             getName(), static_cast<int>(nn));
    GMX_RELEASE_ASSERT(sz.size() <= nn, buf);
    snprintf(buf, sizeof(buf), "More q values than zeta/row values for %s n = %d\n",
             getName(), static_cast<int>(nn));
    GMX_RELEASE_ASSERT(sq.size() <= nn, buf);
    snprintf(buf, sizeof(buf), "More row values than q/zeta values for %s n = %d\n",
             getName(), static_cast<int>(nn));
    GMX_RELEASE_ASSERT(sr.size() <= nn, buf);

    for (size_t n = 0; n < nn; n++)
    {
        RowZetaQ rzq(atoi(sr[n].c_str()), atof(sz[n].c_str()), atof(sq[n].c_str()));
        rzq_.push_back(rzq);
    }
}

typedef struct {
    ChargeDistributionModel eqd;
    const char             *name;
    const char             *ref;
    gmx_bool                bWeight;
} t_eemtype_props;

t_eemtype_props eemtype_props[eqdNR] = {
    { eqdAXp,      "AXp",      "Ghahremanpour2017a",   false },
    { eqdAXg,      "AXg",      "Ghahremanpour2017a",   true },
    { eqdAXs,      "AXs",      "Ghahremanpour2017a",   false },
    { eqdAXpp,     "AXpp",     "Ghahremanpour2017a",   false },
    { eqdAXpg,     "AXpg",     "Ghahremanpour2017a",   true },
    { eqdAXps,     "AXps",     "Ghahremanpour2017a",   false },
    { eqdYang,     "Yang",     "Yang2006b",            true },
    { eqdBultinck, "Bultinck", "Bultinck2002a",        false },
    { eqdRappe,    "Rappe",    "Rappe1991a",           true }
};

ChargeDistributionModel name2eemtype(const std::string name)
{
    for (auto i = 0; i < eqdNR; i++)
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
    for (auto i = 0; i < eqdNR; i++)
    {
        if (eem == eemtype_props[i].eqd)
        {
            return eemtype_props[i].name;
        }
    }
    return nullptr;
}

} // namespace alexandria
