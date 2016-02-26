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
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
#include "stringutil.h"

namespace alexandria
{

Poldata::Poldata()
    :
      _gtDihedralFunction(egdNR),
      _gtDihedralFtype(egdNR),
      _gtDihedral(egdNR)
{
    _nexcl        = 0;
    _gtBondFtype  = 0;
    _gtAngleFtype = 0;

    /* Initiate some crucial variables */
    _nexcl          = 0;
    _fudgeQQ        = 1;
    _fudgeLJ        = 1;
    _gtCombRule     = 0;
    _gtVdwFtype     = F_LJ;
    _gtBondFtype    = F_BONDS;
    _gtAngleFtype   = F_ANGLES;
    _gtDihedralFtype[egdPDIHS] = F_PDIHS;
    _gtDihedralFtype[egdIDIHS] = F_IDIHS;
}

void Poldata::setFilename(const std::string &fn2)
{
    GMX_RELEASE_ASSERT((fn2.size() > 0),
                       "Trying to set empty Poldata filename");
        
    if (0 != _filename.size())
    {
        fprintf(stderr, "Changing Poldata _filename from %s to %s\n",
                _filename.c_str(), fn2.c_str());

    }
    _filename = fn2;
}

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
    _gtVdwFtype    = i;
    _gtVdwFunction = func;
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
    _gtCombRule        = i;
    _gtCombinationRule = func;
}

void Poldata::addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double             polarizability,
                       double             sigPol)
{

    unsigned int      i;

    for (i = 0; (i < _ptype.size()); i++)
    {
        if (_ptype[i].getType().compare(ptype) == 0)
        {
            break;
        }
    }
    if (i == _ptype.size())
    {
        Ptype sp(ptype, miller, bosque, polarizability, sigPol);
        _ptype.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Polarizability type %s was already added to Poldata record\n", ptype.c_str());
    }
}

void Poldata::addBtype(const std::string &btype)
{
    size_t i;

    for (i = 0; (i < _btype.size()); i++)
    {
        if (btype.compare(_btype[i]) == 0)
        {
            break;
        }
    }
    if (i == _btype.size())
    {
        _btype.push_back(btype);
    }
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

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (_alexandria[i].getType().compare(atype) == 0)
        {
            break;
        }
    }
    if (i == _alexandria.size())
    {
        Ffatype sp(desc, atype, ptype, btype,
                   elem, vdwparams, refEnthalpy);

        _alexandria.push_back(sp);
    }
    else
    {
        fprintf(stderr, "Atom type %s was already added to Poldata record\n", atype.c_str());
    }
}

void Poldata::setBondFunction(const std::string &fn)
{
    unsigned int i;

    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, fn.c_str()) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Bond function '%s' does not exist in gromacs", fn.c_str());
    }
    _gtBondFtype    = i;
    _gtBondFunction = fn;
}


void Poldata::setAngleFunction(const std::string &fn)
{
    unsigned int i;

    for (i = 0; (i < F_NRE); i++)
    {
        if (strcasecmp(interaction_function[i].name, fn.c_str()) == 0)
        {
            break;
        }
    }
    if (i == F_NRE)
    {
        gmx_fatal(FARGS, "Angle function '%s' does not exist in gromacs", fn.c_str());
    }
    _gtAngleFtype    = i;
    _gtAngleFunction.assign(fn);
}


void Poldata::setDihedralFunction(int egd, const std::string &func)
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
        gmx_fatal(FARGS, "Dihedral function '%s' does not exist in gromacs", func.c_str());
    }
    _gtDihedralFtype[egd]    = i;
    _gtDihedralFunction[egd].assign(func);
}

bool Poldata::setPtypePolarizability(const std::string &ptype,
                                     double polarizability,
                                     double sigPol)
{
    PtypeIterator sp = findPtype(ptype);
    
    if (_ptype.end() != sp)
    {
        sp->setPolarizability(polarizability);
        sp->setSigPol(sigPol);
        
        return true;
    }
    return false;
}

std::string Poldata::getDesc( std::string atype)
{
    unsigned int i;

    if (atype.size() != 0)
    {
        for (i = 0; (i < _alexandria.size()); i++)
        {
            if (_alexandria[i].getType().compare(atype) == 0)
            {
                return _alexandria[i].getDesc();
            }
        }
    }

    return "";
}

gmx_bool Poldata::strcasestrStart(std::string needle, std::string haystack)
{
    std::string ptr;

    ptr = strcasestr(haystack.c_str(), needle.c_str());

    return (ptr == haystack);
}

bool Poldata::getPtypePol(const std::string &ptype,
                          double *polar, 
                          double *sigPol) const
{
    unsigned int j;

    for (j = 0; (j < _ptype.size()); j++)
    {
        if (ptype.compare(_ptype[j].getType()) == 0)
        {
            *polar   = _ptype[j].getPolarizability();
            *sigPol = _ptype[j].getSigPol();
            
            return true;
        }
    }
    return false;
}

bool Poldata::getAtypePol(const std::string &atype,
                          double *polar, 
                          double *sigPol) const
{
    auto fa = findAtype(atype);
    if (_alexandria.end() != fa)
    {
        return getPtypePol(fa->getPtype(), polar, sigPol);
    }
    return false;
}


bool Poldata::getAtypeRefEnthalpy(const std::string &atype,
                                  double            *Href) const
{
    auto fa = findAtype(atype);
    if (_alexandria.end() != fa)
    {
        *Href = atof(fa->getRefEnthalpy().c_str());
        return true;
    }
    return false;
}

bool Poldata::ptypeToMiller(const std::string &ptype,
                            std::string       &miller) const
{
    for (const auto &i : _ptype)
    {
        if (ptype.compare(i.getType()) == 0)
        {
            miller = i.getMiller();
            return true;
        }
    }
    return false;
}

bool Poldata::ptypeToBosque(const std::string &ptype,
                            std::string       &bosque) const
{
    for (const auto &i : _ptype)
    {
        if (ptype.compare(i.getType()) == 0)
        {
            bosque = i.getBosque();
            return true;
        }
    }
    return false;
}

bool Poldata::atypeToPtype(const std::string &atype,
                           std::string       &ptype) const
{
    if (atype.size() == 0)
    {
        return false;
    }
    auto ai = std::find_if(_alexandria.begin(), _alexandria.end(),
                           [atype](Ffatype const &fa)
                           { return fa.getType().compare(atype) == 0; });
    if (ai != _alexandria.end() && ai->getPtype().size() > 0)
    {
        ptype = ai->getPtype();
        return true;
    }
    return false;
}

bool Poldata::atypeToBtype(const std::string &atype,
                           std::string       &btype) const
{
    if (atype.size() == 0)
    {
        return false;
    }
    auto ai = std::find_if(_alexandria.begin(), _alexandria.end(),
                           [atype](Ffatype const &fa)
                           { return fa.getType().compare(atype) == 0; });
    if (ai != _alexandria.end() && ai->getBtype().size() > 0)
    {
        btype = ai->getBtype();
        return true;
    }
    return false;
}

int loGtbComp(GtBond *ba, GtBond *bb)
{
    std::string a1, a2, b1, b2;
    int         i;

    if (ba->getAtom1().compare(ba->getAtom2()) <= 0)
    {
        a1 = ba->getAtom1();
        a2 = ba->getAtom2();
    }
    else
    {
        a2 = ba->getAtom1();
        a1 = ba->getAtom2();
    }
    if (bb->getAtom1().compare(bb->getAtom2()) <= 0)
    {
        b1 = bb->getAtom1();
        b2 = bb->getAtom2();
    }
    else
    {
        b2 = bb->getAtom1();
        b1 = bb->getAtom2();
    }
    i = a1.compare(b1);
    if (0 == i)
    {
        i = a2.compare(b2);
    }

    return i;
}


int Poldata::gtbComp(const void *a, const void *b)
{
    GtBond    *ba = (GtBond *)a;
    GtBond    *bb = (GtBond *)b;
    int        i;

    i = loGtbComp(ba, bb);
    if ((0 == i) && ((0 != ba->getBondorder()) && (0 != bb->getBondorder())))
    {
        if (ba->getBondorder() < bb->getBondorder())
        {
            i = -1;
        }
        else if (ba->getBondorder() > bb->getBondorder())
        {
            i = 1;
        }
    }
    return i;
}

GtBondIterator Poldata::findBond(const std::string &atom1,
                                 const std::string &atom2,
                                 double bondorder)
{
    return std::find_if(_gtBond.begin(), _gtBond.end(),
                        [atom1,atom2,bondorder](const GtBond gtb)
                        { return ((((atom1.compare(gtb.getAtom1()) == 0) &&
                                    (atom2.compare(gtb.getAtom2()) == 0)) ||
                                   ((atom1.compare(gtb.getAtom2()) == 0) &&
                                    (atom2.compare(gtb.getAtom1()) == 0))) &&
                                  (bondorder == gtb.getBondorder()));
                        });
}

GtBondConstIterator Poldata::findBond(const std::string &atom1,
                                      const std::string &atom2,
                                      double bondorder) const
{
    GtBondConstIterator bb = _gtBond.begin(), be = _gtBond.end();
    return std::find_if(bb, be,
                        [atom1,atom2,bondorder](const GtBond gtb)
                        { return ((((atom1.compare(gtb.getAtom1()) == 0) &&
                                    (atom2.compare(gtb.getAtom2()) == 0)) ||
                                   ((atom1.compare(gtb.getAtom2()) == 0) &&
                                    (atom2.compare(gtb.getAtom1()) == 0))) &&
                                  (bondorder == gtb.getBondorder()));
                        });
}

void Poldata::addMiller(const std::string &miller,
                        int                atomnumber,
                        double             tauAhc,
                        double             alphaAhp,
                        const std::string &alexandria_equiv)
{
    Miller mil(miller, atomnumber, tauAhc, alphaAhp, alexandria_equiv);

    _miller.push_back(mil);
}

bool Poldata::getMillerPol(const std::string &miller,
                           int               *atomnumber,
                           double            *tauAhc,
                           double            *alphaAhp,
                           std::string       &alexandria_equiv) const
{
    MillerConstIterator mb  = _miller.begin(), me = _miller.end();
    MillerConstIterator mil = std::find_if(mb, me, [miller](Miller const &m)
                                           { return (miller.compare(m.getMiller()) == 0); });
    if (_miller.end() != mil)
    {
        *atomnumber = mil->getAtomnumber();
        *tauAhc     = mil->getTauAhc();
        *alphaAhp   = mil->getAlphaAhp();
        alexandria_equiv = mil->getAlexandriaEquiv();
        
        return true;
    }
    return false;
}

bool Poldata::getBosquePol(const std::string &bosque,
                           double            *polarizability) const
{
    BosqueConstIterator bb  = _bosque.begin(), be = _bosque.end();
    BosqueConstIterator bos = std::find_if(bb, be, [bosque](Bosque const &b)
                                           { return (bosque.compare(b.getBosque()) == 0); });
    if (_bosque.end() != bos)
    {
        *polarizability = bos->getPolarizability();
        
        return true;
    }
    return false;
}

/*
 * _gtBond stuff
 */
bool Poldata::setBondParams(const std::string &atom1, 
                            const std::string &atom2,
                            double length, 
                            double sigma, 
                            int ntrain,
                            double bondorder, 
                            const std::string &params)
{
    GtBondIterator gtB = findBond(atom1, atom2, bondorder);
    
    if (getBondEnd() != gtB)
    {
        gtB->setLength(length);
        gtB->setSigma(sigma);
        gtB->setNtrain(ntrain);
        gtB->setParams(params);
        
        return true;
    }
    return false;
}

bool Poldata::addBond(const std::string &atom1,
                      const std::string &atom2,
                      double length, 
                      double sigma,
                      int ntrain,
                      double bondorder, 
                      const std::string &params)
{
    if (setBondParams(atom1, atom2, length, sigma, ntrain, bondorder, params))
    {
        return true;
    }
    
    FfatypeIterator a1, a2;
    
    if (((a1 = btypeToAtype(atom1)) == _alexandria.end()) ||
        ((a2 = btypeToAtype(atom2)) == _alexandria.end()))
    {
        return false;
    }
    GtBond bond(atom1, atom2, params,
                a1->getElem(), a2->getElem(),
                length, sigma, bondorder, ntrain);
    
    _gtBond.push_back(bond);
    
    return true;
}

bool Poldata::searchBond(const std::string &atom1,
                         const std::string &atom2,
                         double *length,
                         double *sigma, 
                         int *ntrain,
                         double *bondorder, 
                         std::string &params) const
{
    int imin = 1, imax = 6;
    if (*bondorder > 0)
    {
        // only look for this particular bond order
        imin = imax = *bondorder;
    }
    for(int i = imin; (i <= imax); i++)
    {
        GtBondConstIterator gtB = findBond(atom1, atom2, i);
        if (_gtBond.end() != gtB)
        {
            *length = gtB->getLength();
            *sigma = gtB->getSigma();
            *ntrain = gtB->getNtrain();
            *bondorder = gtB->getBondorder();
            params.assign(gtB->getParams());

            return true;
        }
    }
    return false;
}

/*
 * gt_angle stuff
 */
bool Poldata::setAngleParams(const std::string &atom1,
                             const std::string &atom2,
                             const std::string &atom3, 
                             double angle, 
                             double sigma, 
                             int ntrain,
                             const std::string &params)
{
    GtAngleIterator gta = findAngle(atom1, atom2, atom3);
    
    if (_gtAngle.end() != gta)
    {
        gta->setAngle(angle);
        gta->setSigma(sigma);
        gta->setNtrain(ntrain);
        gta->setParams(params);
        
        return true;
    }
    return false;
}

bool Poldata::addAngle(const std::string &atom1, 
                       const std::string &atom2,
                       const std::string &atom3,
                       double angle,
                       double sigma,
                       int ntrain, 
                       const std::string &params)
{
    if (setAngleParams(atom1, atom2, atom3,
                       angle, sigma, ntrain, params))
    {
        return true;
    }
    
    FfatypeIterator a1, a2, a3;
    
    if (((a1 = btypeToAtype(atom1)) == _alexandria.end()) ||
        ((a2 = btypeToAtype(atom2)) == _alexandria.end()) ||
        ((a3 = btypeToAtype(atom3)) == _alexandria.end()))
    {
        return false;
    }

    GtAngle ang(atom1, atom2, atom3,
                params, angle, sigma, ntrain);

    _gtAngle.push_back(ang);
    
    return false;
}

GtAngleIterator Poldata::findAngle(const std::string &atom1,
                                   const std::string &atom2,
                                   const std::string &atom3)
{
    return std::find_if(_gtAngle.begin(), _gtAngle.end(),
                        [atom1,atom2,atom3](const GtAngle gta)
                        { return ((((atom1.compare(gta.getAtom1()) == 0) &&
                                    (atom3.compare(gta.getAtom3()) == 0)) ||
                                   ((atom1.compare(gta.getAtom3()) == 0) &&
                                    (atom3.compare(gta.getAtom1()) == 0))) &&
                                  (atom2.compare(gta.getAtom2()) == 0));
                        });
}

GtAngleConstIterator Poldata::findAngle(const std::string &atom1,
                                        const std::string &atom2,
                                        const std::string &atom3) const
{
    GtAngleConstIterator ba = _gtAngle.begin(), ea = _gtAngle.end();
    return std::find_if(ba, ea,
                        [atom1,atom2,atom3](const GtAngle gta)
                        { return ((((atom1.compare(gta.getAtom1()) == 0) &&
                                    (atom3.compare(gta.getAtom3()) == 0)) ||
                                   ((atom1.compare(gta.getAtom3()) == 0) &&
                                    (atom3.compare(gta.getAtom1()) == 0))) &&
                                  (atom2.compare(gta.getAtom2()) == 0));
                        });
}

bool Poldata::searchAngle(const std::string &atom1,
                          const std::string &atom2,
                          const std::string &atom3, 
                          double *angle, 
                          double *sigma,
                          int *ntrain, 
                          std::string &params) const
{
    GtAngleConstIterator gtA = findAngle(atom1, atom2, atom3);
    
    if (_gtAngle.end() != gtA)
    {
        *angle  = gtA->getAngle();
        *sigma  = gtA->getSigma();
        *ntrain = gtA->getNtrain();
        params.assign(gtA->getParams());
        
        return true;
    }
    return false;
}

/*
 * gt_dihedral stuff
 */
bool GtDihedral::compare(const GtDihedral &gtB) const
{
    return ((getAtom1().compare(gtB.getAtom1()) &&
             getAtom2().compare(gtB.getAtom2()) &&
             getAtom3().compare(gtB.getAtom3()) &&
             getAtom4().compare(gtB.getAtom4())) ||
            (getAtom1().compare(gtB.getAtom4()) &&
             getAtom2().compare(gtB.getAtom3()) &&
             getAtom3().compare(gtB.getAtom2()) &&
             getAtom4().compare(gtB.getAtom1())));
}

DihedralIterator Poldata::findDihedral(int egd,
                                       const std::string &atom1, 
                                       const std::string &atom2,
                                       const std::string &atom3, 
                                       const std::string &atom4)
{
    FfatypeIterator a1 = findAtype(atom1), a2 = findAtype(atom2), 
        a3 = findAtype(atom3), a4 = findAtype(atom4);
    if (a1 != _alexandria.end() && a2 != _alexandria.end() &&
        a3 != _alexandria.end() && a4 != _alexandria.end())
    {
        GtDihedral gtA(a1->getBtype(), a2->getBtype(),
                       a3->getBtype(), a4->getBtype(),
                       "", 0, 0, 0);
        return std::find_if(_gtDihedral[egd].begin(),
                            _gtDihedral[egd].end(),
                            [gtA](const GtDihedral &gtB)
                            {
                                return gtA.compare(gtB);
                            });
    }
    return _gtDihedral[egd].end();
}

DihedralConstIterator Poldata::findDihedral(int egd,
                                            const std::string &atom1, 
                                            const std::string &atom2,
                                            const std::string &atom3, 
                                            const std::string &atom4) const
{
    FfatypeConstIterator a1 = findAtype(atom1), a2 = findAtype(atom2), 
        a3 = findAtype(atom3), a4 = findAtype(atom4);
    if (a1 != _alexandria.end() && a2 != _alexandria.end() &&
        a3 != _alexandria.end() && a4 != _alexandria.end())
    {
        GtDihedral gtA(a1->getBtype(), a2->getBtype(),
                       a3->getBtype(), a4->getBtype(),
                       "", 0, 0, 0);
        DihedralConstIterator bd = _gtDihedral[egd].begin(), be = _gtDihedral[egd].end();
        return std::find_if(bd, be, [gtA](const GtDihedral &gtB)
                            {
                                return gtA.compare(gtB);
                            });
    }
    return _gtDihedral[egd].end();
}

bool Poldata::setDihedralParams(int egd,
                                const std::string &atom1, 
                                const std::string &atom2,
                                const std::string &atom3, 
                                const std::string &atom4,
                                double dihedral, 
                                double sigma,
                                int ntrain,
                                const std::string &params)
{
    DihedralIterator gtB = findDihedral(egd, atom1, atom2, atom3, atom4);
    if (_gtDihedral[egd].end() != gtB)
    {
        gtB->setDihedral(dihedral);
        gtB->setSigma(sigma);
        gtB->setNtrain(ntrain);
        gtB->setParams(params);
        
        return true;
    }
    return false;
}

bool Poldata::addDihedral(int egd,
                          const std::string &atom1, 
                          const std::string &atom2,
                          const std::string &atom3, 
                          const std::string &atom4, 
                          double dihedral,
                          double sigma, 
                          int ntrain, 
                          const std::string &params)
{
    if (setDihedralParams(egd, atom1, atom2,
                          atom3, atom4, dihedral,
                          sigma, ntrain, params))
    {
        return true;
    }
    FfatypeConstIterator a1 = findAtype(atom1), a2 = findAtype(atom2), 
        a3 = findAtype(atom3), a4 = findAtype(atom4);
    if (a1 != _alexandria.end() && a2 != _alexandria.end() &&
        a3 != _alexandria.end() && a4 != _alexandria.end())
    {
        GtDihedral dihed(atom1, atom2, atom3, atom4,
                         params, dihedral, sigma, ntrain);
    
        _gtDihedral[egd].push_back(dihed);
        return true;
    }
    return false;
}

bool Poldata::searchDihedral(int egd,
                             const std::string &atom1,
                             const std::string &atom2,
                             const std::string &atom3,
                             const std::string &atom4,
                             double *dihedral, 
                             double *sigma,
                             int *ntrain, 
                             std::string &params) const
{
    DihedralConstIterator gtRes = findDihedral(egd, atom1, atom2, atom3, atom4);
    if (_gtDihedral[egd].end() != gtRes)
    {
        *dihedral = gtRes->getDihedral();
        *sigma    = gtRes->getSigma();
        *ntrain   = gtRes->getNtrain();
        params    = gtRes->getParams();

        return true;
    }
    return false;
}

void Poldata::addSymcharges(const std::string &central,
                            const std::string &attached,
                            int numattach)
{
    unsigned int           i;
    Symcharges           * sc;
    for (i = 0; (i < _symcharges.size()); i++)
    {
        sc = &(_symcharges[i]);
        if ((strcasecmp(sc->getCentral().c_str(), central.c_str()) == 0) &&
            (strcasecmp(sc->getAttached().c_str(), attached.c_str()) == 0) &&
            (sc->getNumattach() == numattach))
        {
            break;
        }
    }
    if (i == _symcharges.size())
    {
        Symcharges symcharges(central, attached, numattach);

        _symcharges.push_back(symcharges);
    }
}

int Poldata::getNumprops(ChargeDistributionModel eqdModel) const
{
    unsigned int i, n = 0;

    for (i = 0; (i < _eep.size()); i++)
    {
        if (_eep[i].getEqdModel() == eqdModel)
        {
            n++;
        }
    }

    return n;
}

int Poldata::havePolSupport(const std::string &atype) const
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (atype.compare(_alexandria[i].getType()) == 0)
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
                        int zz) const
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
                     int zz) const
{
    EempropsConstIterator eer;

    if ((eer = findEem(eqdModel, name)) != EndEemprops())
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getQ(zz);
    }
    return -1;
}

double Poldata::getChi0(ChargeDistributionModel  eqdModel,
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
    for (auto &i : _epr)
    {
        if (i.getEqdModel() == eqdModel)
        {
            i.setEpref(epref);
            return;
        }
    }
    Epref epr(eqdModel, epref);
    _epr.push_back(epr);
}

const char *Poldata::getEpref(ChargeDistributionModel eqdModel) const
{
    for(auto &i : _epr)
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
            gmx_send_int(cr, i, _eep.size());
            gmx_send(cr, i, _eep.data(), _eep.size()*sizeof(_eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != _eep.size())
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        _eep.clear();
        gmx_recv(cr, 0, ep.data(), _eep.size()*sizeof(ep[0]));
        for (i = 0; (i < nep); i++)
        {
            _eep.push_back(ep[i]);
        }
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    getEemtypeName(_eep[i].getEqdModel()),
                    _eep[i].getName(), _eep[i].getChi0(),
                    _eep[i].getJ0());
            for (j = 0; ((int)j < _eep[i].getNzeta()); j++)
            {
                fprintf(debug, " %8.3f", _eep[i].getZeta(j));
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
            return _eep.end();
        }
    }
    return std::find_if(_eep.begin(), _eep.end(),
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
            return _eep.end();
        }
    }
    return std::find_if(_eep.begin(), _eep.end(),
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
    
    char buf[256];
    snprintf(buf, sizeof(buf), "More zeta values than q/row values for %s n = %d\n",
             getName(), static_cast<int>(nn));
    GMX_RELEASE_ASSERT(sz.size() <= nn, buf);
    snprintf(buf, sizeof(buf), "More q values than zeta/row values for %s n = %d\n",
             getName(), static_cast<int>(nn));
    GMX_RELEASE_ASSERT(sq.size() <= nn, buf);
    snprintf(buf, sizeof(buf), "More row values than q/zeta values for %s n = %d\n",
             getName(), static_cast<int>(nn));
    GMX_RELEASE_ASSERT(sr.size() <= nn, buf);

    for(size_t n = 0; n < nn; n++)
    {
        RowZetaQ rzq(atoi(sr[n].c_str()), atof(sz[n].c_str()), atof(sq[n].c_str()));
        rzq_.push_back(rzq);
    }
}
                   
Eemprops::Eemprops(ChargeDistributionModel eqdModel,
                   const std::string      &name,
                   const std::string      &rowstr,
                   const std::string      &zetastr,
                   const std::string      &qstr,
                   double                  J0,
                   double                  chi0) :
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

} // namespace alexandria
