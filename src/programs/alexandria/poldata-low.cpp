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
 
#include "poldata-low.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
#include "plistwrapper.h"


namespace alexandria
{

static const char * eit_names[eitNR] = {
    "BONDS", "ANGLES", "LINEAR_ANGLES",
    "PROPER_DIHEDRALS", "IMPROPER_DIHEDRALS", "VDW",
    "LJ14", "POLARIZATION", "CONSTR", "VSITE2"
};

static const char * eqd_names[eqdNR] = {
     "eqdAXp", "eqdAXg", "eqdAXs", "eqdAXpp", "eqdAXpg", "eqdAXps", 
     "eqdYang", "eqdBultinck", "eqdRappe"
};

const char *eqd2string(ChargeDistributionModel eqd)
{
    if (eqd < eqdNR)
    {
        return eqd_names[eqd];
    }
    return nullptr;
}

ChargeDistributionModel string2eqd(const char *string)
{
    int i;    
    for (i = 0; i < eqdNR; i++)
    {
        if (gmx_strcasecmp(string, eqd_names[i]) == 0)
        {
            return static_cast<ChargeDistributionModel>(i);
        }
    }
    return eqdNR;
}

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
    for (i = 0; (i < eitNR); i++)
    {
        if (gmx_strcasecmp(string, eit_names[i]) == 0)
        {
            return static_cast<InteractionType>(i);
        }
    }
    return eitNR;
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

CommunicationStatus Ptype::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, type_.c_str());
        gmx_send_str(cr, dest, miller_.c_str());
        gmx_send_str(cr, dest, bosque_.c_str());
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

CommunicationStatus Ptype::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        type_.assign(gmx_recv_str(cr, src));
        miller_.assign(gmx_recv_str(cr, src));
        bosque_.assign(gmx_recv_str(cr, src));
        polarizability_ = gmx_recv_double(cr, src);
        sigPol_ = gmx_recv_double(cr, src);
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

Ffatype::Ffatype(const std::string &desc,
                 const std::string &type,
                 const std::string &ptype,
                 const std::string &btype,
                 const std::string &elem,
                 std::string       &vdwparams,
                 const std::string &refEnthalpy)
    :
      desc_(desc),
      type_(type),
      ptype_(ptype),
      btype_(btype),
      elem_(elem),
      vdwparams_(vdwparams),
      refEnthalpy_(refEnthalpy)
{}

CommunicationStatus Ffatype::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, desc_.c_str());
        gmx_send_str(cr, dest, type_.c_str());
        gmx_send_str(cr, dest, ptype_.c_str());
        gmx_send_str(cr, dest, btype_.c_str());
        gmx_send_str(cr, dest, elem_.c_str());
        gmx_send_str(cr, dest, vdwparams_.c_str());
        gmx_send_str(cr, dest, refEnthalpy_.c_str());
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Fftype %s %s %s %s %s %s %s, status %s\n",
                    desc_.c_str(), type_.c_str(), ptype_.c_str(),
                    btype_.c_str(), elem_.c_str(), vdwparams_.c_str() ,
                    refEnthalpy_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;

}
        
CommunicationStatus Ffatype::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        desc_.assign(gmx_recv_str(cr, src));
        type_.assign(gmx_recv_str(cr, src));
        ptype_.assign(gmx_recv_str(cr, src));        
        btype_.assign(gmx_recv_str(cr, src));
        elem_.assign(gmx_recv_str(cr, src));
        vdwparams_.assign(gmx_recv_str(cr, src));
        refEnthalpy_.assign(gmx_recv_str(cr, src));
        if (nullptr != debug)
        {
            fprintf(debug, "Receive Fftype %s %s %s %s %s %s %s, status %s\n",
                    desc_.c_str(), type_.c_str(), ptype_.c_str(),
                    btype_.c_str(), elem_.c_str(), vdwparams_.c_str() ,
                    refEnthalpy_.c_str(), cs_name(cs));
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

CommunicationStatus ListedForce::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, params_.c_str());
        gmx_send_double(cr, dest, refValue_);
        gmx_send_double(cr, dest, sigma_);
        gmx_send_int(cr, dest, static_cast<int>(ntrain_));
        gmx_send_int(cr, dest, atoms_.size());
        
        for (auto &atom : atoms_)
        {
            gmx_send_str(cr, dest, atom.c_str());
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
        
CommunicationStatus ListedForce::Receive(t_commrec *cr, int src)
{
    int natom;
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        params_.assign(gmx_recv_str(cr, src));
        refValue_ = gmx_recv_double(cr, src);
        sigma_    = gmx_recv_double(cr, src);       
        ntrain_   = static_cast<size_t>(gmx_recv_double(cr, src));
        natom     = gmx_recv_int(cr, src);
        
        for(int n = 0; n < natom; n++)
        {
            char *atom = gmx_recv_str(cr, src);
            if (nullptr != atom)
            {
                const_cast<std::vector<std::string>&>(atoms_).push_back(atom);
                free(atom);
            }
            else
            {
                gmx_fatal(FARGS, "A category was promised but I got a NULL pointer");
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
    unsigned int funcType;

    for (funcType = 0; (funcType < F_NRE); funcType++)
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


CommunicationStatus ListedForces::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, iType2string(iType_));
        gmx_send_str(cr, dest, function_.c_str());
        gmx_send_str(cr, dest, unit_.c_str());
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
        
CommunicationStatus ListedForces::Receive(t_commrec *cr, int src)
{
    size_t nforce;
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        iType_    = string2iType(gmx_recv_str(cr, src));
        const_cast<std::string&>(function_) = gmx_recv_str(cr, src);
        const_cast<std::string&>(unit_)     = gmx_recv_str(cr, src);
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

CommunicationStatus Bosque::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, bosque_.c_str());
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
        
CommunicationStatus Bosque::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        bosque_ = gmx_recv_str(cr, src);
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
   
CommunicationStatus Miller::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, miller_.c_str());
        gmx_send_int(cr, dest, atomnumber_);
        gmx_send_double(cr, dest, tauAhc_);
        gmx_send_double(cr, dest, alphaAhp_);
        gmx_send_str(cr, dest, alexandria_equiv_.c_str());
       
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
        
CommunicationStatus Miller::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        miller_     = gmx_recv_str(cr, src);
        atomnumber_ = gmx_recv_int(cr, src);
        tauAhc_     = gmx_recv_double(cr, src);
        alphaAhp_   = gmx_recv_double(cr, src);
        alexandria_equiv_ = gmx_recv_str(cr, src);
        
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

CommunicationStatus Symcharges::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, central_.c_str());
        gmx_send_str(cr, dest, attached_.c_str());
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
        
CommunicationStatus Symcharges::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        const_cast<std::string&>(central_)   = gmx_recv_str(cr, src);
        const_cast<std::string&>(attached_)  = gmx_recv_str(cr, src);
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

CommunicationStatus Epref::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, eqd2string(eqdModel_));
        gmx_send_str(cr, dest, epref_.c_str());
        
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Epref %s %s, status %s\n",
                    eqd2string(eqdModel_), epref_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}
        
CommunicationStatus Epref::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        eqdModel_ = string2eqd(gmx_recv_str(cr, src));
        epref_    = gmx_recv_str(cr, src);
        
        if (nullptr != debug)
        {
            fprintf(debug, "Received Epref %s %s, status %s\n",
                    eqd2string(eqdModel_), epref_.c_str(), cs_name(cs));
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

CommunicationStatus RowZetaQ::Send(t_commrec *cr, int dest)
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
        
CommunicationStatus RowZetaQ::Receive(t_commrec *cr, int src)
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

Eemprops::Eemprops(ChargeDistributionModel  eqdModel,
                   const std::string        &name,
                   const std::string        &rowstr,
                   const std::string        &zetastr,
                   const std::string        &qstr,
                   double                    J0,
                   double                    chi0)
    :
      eqdModel_(eqdModel),
      name_(name),
      rowstr_(rowstr),
      zetastr_(zetastr),
      qstr_(qstr),
      J0_(J0),
      chi0_(chi0)
{
    setRowZetaQ(rowstr, zetastr, qstr);
}

CommunicationStatus Eemprops::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, eqd2string(eqdModel_));
        gmx_send_str(cr, dest, name_.c_str());
        gmx_send_str(cr, dest, rowstr_.c_str());
        gmx_send_str(cr, dest, zetastr_.c_str());
        gmx_send_str(cr, dest, qstr_.c_str());
        gmx_send_double(cr, dest, J0_);
        gmx_send_double(cr, dest, chi0_);
        gmx_send_int(cr, dest, rzq_.size());
        
        for (auto &rzq : rzq_)
        {
            cs = rzq.Send(cr, dest);
        } 
        
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Eemprops %s %s %s %s %s %g %g, status %s\n",
                    eqd2string(eqdModel_), name_.c_str(), rowstr_.c_str(), 
                    zetastr_.c_str(), qstr_.c_str(), J0_, chi0_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}
        
CommunicationStatus Eemprops::Receive(t_commrec *cr, int src)
{
    size_t nrzq;
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        eqdModel_ = string2eqd(gmx_recv_str(cr, src));
        name_     = gmx_recv_str(cr, src);
        rowstr_   = gmx_recv_str(cr, src);
        zetastr_  = gmx_recv_str(cr, src);
        qstr_     = gmx_recv_str(cr, src);
        J0_       = gmx_recv_double(cr, src);
        chi0_     = gmx_recv_double(cr, src);
        nrzq      = gmx_recv_int(cr, src);
        
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

} // namespace alexandria
