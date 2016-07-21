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

#include "poldata-low.h"

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
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
#include "plistwrapper.h"
#include "stringutil.h"

namespace alexandria
{

static const char * eit_names[eitNR] = {
    "BONDS", "ANGLES", "LINEAR_ANGLES",
    "PROPER_DIHEDRALS", "IMPROPER_DIHEDRALS",
    "LJ14", "POLARIZATION", "CONSTR", "VSITE2"
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

Ffatype::Ffatype(const std::string &desc,
                 const std::string &type,
                 const std::string &ptype,
                 const std::string &btype,
                 const std::string &elem,
                 const std::string &vdwparams,
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

Symcharges::Symcharges(const std::string &central,
                       const std::string &attached,
                       int                numattach)
    :
    central_(central),
    attached_(attached),
    numattach_(numattach)
{}

Epref::Epref(ChargeDistributionModel  eqdModel,
             const std::string       &epref)
    :
    eqdModel_(eqdModel),
    epref_(epref)
{}

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

Eemprops::Eemprops(ChargeDistributionModel   eqdModel,
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
