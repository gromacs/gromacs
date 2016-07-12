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
#include "stringutil.h"


namespace alexandria
{

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

GtBond::GtBond(const std::string btype1,
               const std::string btype2,
               const std::string params,
               double            length,
               double            sigma,
               double            bondorder,
               int               ntrain)
    :
      atom1_(btype1),
      atom2_(btype2),
      params_(params),
      length_(length),
      sigma_(sigma),
      bondorder_(bondorder),
      ntrain_(ntrain)
{}

GtBonds::GtBonds(const std::string function,
                 const std::string unit)
    :
      gtBondFunction_(function),
      gtLengthUnit_(unit)
{
    unsigned int Functype;

    for (Functype = 0; (Functype < F_NRE); Functype++)
    {
        if (strcasecmp(interaction_function[Functype].name, gtBondFunction_.c_str()) == 0)
        {
            break;
        }
    }
    if (Functype == F_NRE)
    {
        gmx_fatal(FARGS, "Bond function '%s' does not exist in gromacs", gtBondFunction_.c_str());
    }
    gtBondFtype_ = Functype;
}

GtBondIterator GtBonds::findBond(const std::string &btype1,
                                 const std::string &btype2,
                                 double             bondorder)
{
    GtBondIterator bb = getBondBegin(), be = getBondEnd();
    return std::find_if(bb, be, [btype1, btype2, bondorder](const GtBond &gtb)
                        { return ((((btype1.compare(gtb.getAtom1()) == 0)   &&
                                    (btype2.compare(gtb.getAtom2()) == 0))  ||
                                   ((btype1.compare(gtb.getAtom2()) == 0)   &&
                                    (btype2.compare(gtb.getAtom1()) == 0))) &&
                                  (0 == bondorder || bondorder == gtb.getBondorder())); });
}

GtBondConstIterator GtBonds::findBond(const std::string &btype1,
                                      const std::string &btype2,
                                      double             bondorder) const
{
    GtBondConstIterator bb = getBondBegin(), be = getBondEnd();
    return std::find_if(bb, be, [btype1, btype2, bondorder](const GtBond &gtb)
                        { return ((((btype1.compare(gtb.getAtom1()) == 0)   &&
                                    (btype2.compare(gtb.getAtom2()) == 0))  ||
                                   ((btype1.compare(gtb.getAtom2()) == 0)   &&
                                    (btype2.compare(gtb.getAtom1()) == 0))) &&
                                  (0 == bondorder || bondorder == gtb.getBondorder())); });
}

bool GtBonds::setBondParams(const std::string &btype1,
                            const std::string &btype2,
                            double             length,
                            double             sigma,
                            int                ntrain,
                            double             bondorder,
                            const std::string &params)
{
    auto gtb = findBond(btype1, btype2, bondorder);

    if (getBondEnd() != gtb)
    {
        gtb->setLength(length);
        gtb->setSigma(sigma);
        gtb->setNtrain(ntrain);
        gtb->setParams(params);

        return true;
    }
    return false;
}

void GtBonds::addBond(const std::string &btype1,
                      const std::string &btype2,
                      double             length,
                      double             sigma,
                      int                ntrain,
                      double             bondorder,
                      const std::string &params)
{

    if (setBondParams(btype1, btype2, length, sigma,
                      ntrain, bondorder, params))
    {
        return;
    }

    GtBond gtb(btype1, btype2, params, length,
               sigma, bondorder, ntrain);

    gtBond_.push_back(gtb);

}

bool GtBonds::searchBond(const std::string &btype1,
                         const std::string &btype2,
                         double            *length,
                         double            *sigma,
                         int               *ntrain,
                         double            *bondorder,
                         std::string       &params) const
{
    int minBondOrder = 1, maxBondOrder = 6;

    if (*bondorder > 0)
    {
        minBondOrder = maxBondOrder = *bondorder;
    }
    for (int i = minBondOrder; (i <= maxBondOrder); i++)
    {
        auto gtb = findBond(btype1, btype2, i);
        if (getBondEnd() != gtb)
        {
            *length    = gtb->getLength();
            *sigma     = gtb->getSigma();
            *ntrain    = gtb->getNtrain();
            *bondorder = gtb->getBondorder();
            params.assign(gtb->getParams());

            return true;
        }
    }
    return false;
}

GtAngle::GtAngle(const std::string &btype1,
                 const std::string &btype2,
                 const std::string &btype3,
                 const std::string &params,
                 double             angle,
                 double             sigma,
                 int                ntrain)
    :
      atom1_(btype1),
      atom2_(btype2),
      atom3_(btype3),
      params_(params),
      angle_(angle),
      sigma_(sigma),
      ntrain_(ntrain)
{}

GtAngles::GtAngles(const std::string function,
                   const std::string unit)
    :
      gtAngleFunction_(function),
      gtAngleUnit_(unit)
{
    unsigned int funcType;

    for (funcType = 0; (funcType < F_NRE); funcType++)
    {
        if (strcasecmp(interaction_function[funcType].name, gtAngleFunction_.c_str()) == 0)
        {
            break;
        }
    }
    if (funcType == F_NRE)
    {
        gmx_fatal(FARGS, "Angle function '%s' does not exist in gromacs", gtAngleFunction_.c_str());
    }
    gtAngleFtype_ = funcType;
}

GtAngleIterator GtAngles::findAngle(const std::string &btype1,
                                    const std::string &btype2,
                                    const std::string &btype3)
{
    GtAngleIterator ab = getAngleBegin(), ae = getAngleEnd();
    return std::find_if(ab, ae, [btype1, btype2, btype3](const GtAngle &gta)
                        { return ((((btype1.compare(gta.getAtom1()) == 0)   &&
                                    (btype3.compare(gta.getAtom3()) == 0))  ||
                                   ((btype1.compare(gta.getAtom3()) == 0)   &&
                                    (btype3.compare(gta.getAtom1()) == 0))) &&
                                  (btype2.compare(gta.getAtom2()) == 0)); });
}

GtAngleConstIterator GtAngles::findAngle(const std::string &btype1,
                                         const std::string &btype2,
                                         const std::string &btype3) const
{
    GtAngleConstIterator ab = getAngleBegin(), ae = getAngleEnd();
    return std::find_if(ab, ae, [btype1, btype2, btype3](const GtAngle &gta)
                        { return ((((btype1.compare(gta.getAtom1()) == 0)   &&
                                    (btype3.compare(gta.getAtom3()) == 0))  ||
                                   ((btype1.compare(gta.getAtom3()) == 0)   &&
                                    (btype3.compare(gta.getAtom1()) == 0))) &&
                                  (btype2.compare(gta.getAtom2()) == 0)); });
}

bool GtAngles::setAngleParams(const std::string &btype1,
                              const std::string &btype2,
                              const std::string &btype3,
                              double             angle,
                              double             sigma,
                              int                ntrain,
                              const std::string &params)
{
    auto gta = findAngle(btype1, btype2, btype3);

    if (getAngleEnd() != gta)
    {
        gta->setAngle(angle);
        gta->setSigma(sigma);
        gta->setNtrain(ntrain);
        gta->setParams(params);

        return true;
    }
    return false;
}

void GtAngles::addAngle(const std::string &btype1,
                        const std::string &btype2,
                        const std::string &btype3,
                        double             angle,
                        double             sigma,
                        int                ntrain,
                        const std::string &params)
{

    if (setAngleParams(btype1, btype2, btype3, angle,
                       sigma, ntrain, params))
    {
        return;
    }

    GtAngle gta(btype1, btype2, btype3, params,
                angle, sigma, ntrain);

    gtAngle_.push_back(gta);

}

bool GtAngles::searchAngle(const std::string &btype1,
                           const std::string &btype2,
                           const std::string &btype3,
                           double            *angle,
                           double            *sigma,
                           int               *ntrain,
                           std::string       &params) const
{
    auto gta = findAngle(btype1, btype2, btype3);

    if (getAngleEnd() != gta)
    {
        *angle  = gta->getAngle();
        *sigma  = gta->getSigma();
        *ntrain = gta->getNtrain();
        params.assign(gta->getParams());

        return true;
    }
    return false;
}

GtDihedral::GtDihedral(const std::string &atom1,
                       const std::string &atom2,
                       const std::string &atom3,
                       const std::string &atom4,
                       const std::string &params,
                       double             dihedral,
                       double             sigma,
                       int                ntrain)
    :
      atom1_(atom1),
      atom2_(atom2),
      atom3_(atom3),
      atom4_(atom4),
      params_(params),
      dihedral_(dihedral),
      sigma_(sigma),
      ntrain_(ntrain)
{}

GtDihedrals::GtDihedrals(const std::string function,
                         const std::string unit)
    :
      gtDihedralFunction_(function),
      gtDihedralUnit_(unit)
{
    unsigned int funcType;

    for (funcType = 0; (funcType < F_NRE); funcType++)
    {
        if (strcasecmp(interaction_function[funcType].name, gtDihedralFunction_.c_str()) == 0)
        {
            break;
        }
    }
    if (funcType == F_NRE)
    {
        gmx_fatal(FARGS, "Angle function '%s' does not exist in gromacs", gtDihedralFunction_.c_str());
    }
    gtDihedralFtype_ = funcType;
}

GtDihedralIterator GtDihedrals::findDihedral(const std::string &btype1,
                                             const std::string &btype2,
                                             const std::string &btype3,
                                             const std::string &btype4)
{
    GtDihedralIterator db = getDihedralBegin(), de = getDihedralEnd();
    return std::find_if(db, de, [btype1, btype2, btype3, btype4](const GtDihedral &gtd)
                        { return ((btype1.compare(gtd.getAtom1())  &&
                                   btype2.compare(gtd.getAtom2())  &&
                                   btype3.compare(gtd.getAtom3())  &&
                                   btype4.compare(gtd.getAtom4())) ||
                                  (btype1.compare(gtd.getAtom4())  &&
                                   btype2.compare(gtd.getAtom3())  &&
                                   btype3.compare(gtd.getAtom2())  &&
                                   btype4.compare(gtd.getAtom1()))); });
}

GtDihedralConstIterator GtDihedrals::findDihedral(const std::string &btype1,
                                                  const std::string &btype2,
                                                  const std::string &btype3,
                                                  const std::string &btype4) const
{
    GtDihedralConstIterator db = getDihedralBegin(), de = getDihedralEnd();
    return std::find_if(db, de, [btype1, btype2, btype3, btype4](const GtDihedral &gtd)
                        { return ((btype1.compare(gtd.getAtom1())  &&
                                   btype2.compare(gtd.getAtom2())  &&
                                   btype3.compare(gtd.getAtom3())  &&
                                   btype4.compare(gtd.getAtom4())) ||
                                  (btype1.compare(gtd.getAtom4())  &&
                                   btype2.compare(gtd.getAtom3())  &&
                                   btype3.compare(gtd.getAtom2())  &&
                                   btype4.compare(gtd.getAtom1()))); });
}

bool GtDihedrals::setDihedralParams(const std::string &btype1,
                                    const std::string &btype2,
                                    const std::string &btype3,
                                    const std::string &btype4,
                                    double             dihedral,
                                    double             sigma,
                                    int                ntrain,
                                    const std::string &params)
{
    auto gtd = findDihedral(btype1, btype2, btype3, btype4);
    if (getDihedralEnd() != gtd)
    {
        gtd->setDihedral(dihedral);
        gtd->setSigma(sigma);
        gtd->setNtrain(ntrain);
        gtd->setParams(params);

        return true;
    }
    return false;
}

void GtDihedrals::addDihedral(const std::string &btype1,
                              const std::string &btype2,
                              const std::string &btype3,
                              const std::string &btype4,
                              double             dihedral,
                              double             sigma,
                              int                ntrain,
                              const std::string &params)
{

    if (setDihedralParams(btype1, btype2, btype3, btype4,
                          dihedral, sigma, ntrain, params))
    {
        return;
    }

    GtDihedral gtd(btype1, btype2, btype3, btype4,
                   params, dihedral, sigma, ntrain);

    gtDihedral_.push_back(gtd);

}

bool GtDihedrals::searchDihedral(const std::string &btype1,
                                 const std::string &btype2,
                                 const std::string &btype3,
                                 const std::string &btype4,
                                 double            *dihedral,
                                 double            *sigma,
                                 int               *ntrain,
                                 std::string       &params) const
{
    auto gtd = findDihedral(btype1, btype2, btype3, btype4);
    if (getDihedralEnd() != gtd)
    {
        *dihedral = gtd->getDihedral();
        *sigma    = gtd->getSigma();
        *ntrain   = gtd->getNtrain();
        params    = gtd->getParams();

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

Eemprops::Eemprops(ChargeDistributionModel eqdModel,
                   const std::string      &name,
                   const std::string      &rowstr,
                   const std::string      &zetastr,
                   const std::string      &qstr,
                   double                  J0,
                   double                  chi0)
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
