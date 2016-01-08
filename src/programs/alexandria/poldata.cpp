/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "poldata.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "gromacs/gmxlib/ifunc.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "gmx_simple_comm.h"
#include "stringutil.h"

#define NOTSET -666

#define assignStr(dst, src) if (dst != NULL) {*dst = src; }

#define assignScal(dst, src) if (dst) *dst = src

#define vectorToArray(vec) (&(vec[0]))

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
    _nexcl          = NOTSET;
    _fudgeQQ        = NOTSET;
    _fudgeLJ        = NOTSET;
    _gtCombRule     = NOTSET;
    _gtVdwFtype     = NOTSET;
    _gtBondFtype    = NOTSET;
    _gtAngleFtype   = NOTSET;
    for (int i = 0; (i < egdNR); i++)
    {
        _gtDihedralFtype[i] = NOTSET;
    }

}

void Poldata::setFilename( std::string fn2)
{
    if ("" == fn2)
    {
        fprintf(stderr, "Trying to set Poldata filename to NULL\n");
        return;
    }
    if (0 != _filename.size())
    {
        fprintf(stderr, "Overwriting Poldata _filename from %s to %s\n",
                _filename.c_str(), fn2.c_str());

    }
    _filename = fn2;
}

void Poldata::setVdwFunction( const std::string func)
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


int Poldata::getVdwFtype()
{
    if (NOTSET == _gtVdwFtype)
    {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  _filename.c_str());
    }
    return _gtVdwFtype;
}

int Poldata::getNexcl()
{
    if (NOTSET == _nexcl)
    {
        gmx_fatal(FARGS, "Nexclusions not set in %s", _filename.c_str());
    }
    return _nexcl;
}


double Poldata::getFudgeQQ()
{
    if (NOTSET == _fudgeQQ)
    {
        gmx_fatal(FARGS, "_fudgeQQ not set in %s", _filename.c_str());
    }
    return _fudgeQQ;
}


double Poldata::getFudgeLJ()
{
    if (NOTSET == _fudgeLJ)
    {
        gmx_fatal(FARGS, "_fudgeLJ not set in %s", _filename.c_str());
    }
    return _fudgeLJ;
}

void Poldata::setCombinationRule( std::string func)
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

std::string Poldata::getCombinationRule()
{
    if (NOTSET == _gtVdwFtype)
    {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  _filename.c_str());
    }
    return _gtCombinationRule;
}

int Poldata::getCombRule()
{
    if (NOTSET == _gtCombRule)
    {
        gmx_fatal(FARGS, "No valid combination rule in %s", _filename.c_str());
    }
    return _gtCombRule;
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
                       double            refEnthalpy)
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

void Poldata::addBondingRule(std::string gtBrule, std::string atype,
                             std::string geometry, int numbonds,
                             double valence, int iAromatic,
                             std::string neighbors)
{
    unsigned int      i, j;

    for (j = 0; (j < _alexandria.size()); j++)
    {
        if (_alexandria[j].getType().compare(atype) == 0)
        {
            break;
        }
    }
    if (j < _alexandria.size())
    {
        for (i = 0; (i < _brule.size()); i++)
        {
            if (_brule[i].getRule().compare( gtBrule) == 0)
            {
                break;
            }
        }
        if (i == _brule.size())
        {

            Brule brule(_alexandria[j].getElem(), gtBrule, atype, neighbors, geometry,
                        numbonds, iAromatic, valence, split(neighbors, ' '));

            _brule.push_back(brule);
        }
        else
        {
            fprintf(stderr, "Bonding rule %s was already added to Poldata record\n", gtBrule.c_str());
        }
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Ignoring bonding rule involving unknown atom type %s\n",
                atype.c_str());
    }
}

void Poldata::setBondFunction( std::string fn)
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


void Poldata::setAngleFunction( std::string fn)
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


void Poldata::setDihedralFunction( int egd, std::string func)
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


void Poldata::setPtypePolarizability( const std::string ptype,
                                      double polarizability, double sigPol)
{
    Ptype            *sp;
    unsigned int      i;

    for (i = 0; (i < _ptype.size()); i++)
    {
        if (ptype.compare(_ptype[i].getType()) == 0)
        {
            sp                 = &(_ptype[i]);
            sp->setPolarizability(polarizability);
            sp->setSigPol(sigPol);
            break;
        }
    }
    if (i == _ptype.size())
    {
        fprintf(stderr, "No such ptype %s when trying to set the polarizability.\n", ptype.c_str());
    }
}


std::string Poldata::getLengthUnit()
{
    if (0 == _gtLengthUnit.size())
    {
        gmx_fatal(FARGS, "No length unit in %s",
                  (0 != _filename.size()) ? _filename.c_str() : "unknown");
    }

    return _gtLengthUnit;
}



std::string Poldata::getGeometry( std::string gtBrule)
{
    unsigned int i;

    if (gtBrule.size() != 0)
    {
        for (i = 0; (i < _brule.size()); i++)
        {
            if (_brule[i].getRule().compare(gtBrule) == 0)
            {
                return _brule[i].getGeometry();
            }
        }
    }

    return "";
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

int Poldata::countNeighbors(Brule *brule, int nbond, std::string nbhybrid[], int *score)
{
    int              j, ni = 0, IFound;
    std::vector<int> jj;

    *score = 0;
    jj.resize(nbond+1);
    for (unsigned int i = 0; (i < brule->getNb().size()); i++)
    {
        IFound = 0;
        for (j = 0; (j < nbond); j++)
        {
            if (
                (0 != nbhybrid[j].size()) &&
                (jj[j] == 0) &&
                (IFound == 0) &&
                strcasestrStart(brule->getNb()[i], nbhybrid[j])
                )
            {
                IFound  = 1;
                jj[j]   = j+1;
                ni++;
                (*score) += 1;
                if (nbhybrid[j].size() > 0)
                {
                    (*score) += 1;
                }
            }
        }
    }

    return ni;
}

int Poldata::bondingRuleValence( std::string gtBrule, double *valence)
{
    unsigned int i;

    for (i = 0; (i < _brule.size()); i++)
    {
        if (strcasecmp(gtBrule.c_str(), _brule[i].getRule().c_str()) == 0)
        {
            *valence = _brule[i].getValence();
            return 1;
        }
    }
    return 0;
}

int Poldata::getPtypePol( const std::string ptype,
                          double *polar, double *sigPol)
{
    unsigned int j;

    for (j = 0; (j < _ptype.size()); j++)
    {
        if (ptype.compare(_ptype[j].getType()) == 0)
        {
            if (NULL != polar)
            {
                *polar   = _ptype[j].getPolarizability();
            }

            if (NULL != sigPol)
            {
                *sigPol = _ptype[j].getSigPol();
            }
            return 1;
        }
    }
    return 0;
}

int Poldata::getAtypePol( const std::string atype,
                          double *polar, double *sigPol)
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (atype.compare(_alexandria[i].getType()) == 0)
        {
            return getPtypePol( _alexandria[i].getPtype(), polar, sigPol);
        }
    }
    return 0;
}


int Poldata::getAtypeRefEnthalpy( const std::string atype,
                                  double           *Href)
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (atype.compare(_alexandria[i].getType()) == 0)
        {
            *Href = _alexandria[i].getRefEnthalpy();
            return 1;
        }
    }
    return 0;
}

std::string Poldata::ptypeToMiller( const std::string ptype)
{
    unsigned int i;

    for (i = 0; (i < _ptype.size()); i++)
    {
        if (ptype.compare(_ptype[i].getType()) == 0)
        {
            return _ptype[i].getMiller();
        }
    }
    return "";
}

std::string Poldata::ptypeToBosque( const std::string ptype)
{
    unsigned int i;

    for (i = 0; (i < _ptype.size()); i++)
    {
        if (ptype.compare(_ptype[i].getType()) == 0)
        {
            return _ptype[i].getBosque();
        }
    }
    return "";
}

std::string Poldata::atypeToPtype( const std::string atype)
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (_alexandria[i].getType().compare(atype) == 0)
        {
            return _alexandria[i].getPtype();
        }
    }
    return "";
}

std::string Poldata::atypeToBtype( const std::string atype)
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (_alexandria[i].getType().compare(atype) == 0)
        {
            return _alexandria[i].getBtype();
        }
    }
    return "";
}

int Poldata::searchAtype(std::string         key,
                         Ffatype           * atype)
{
    unsigned int        i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (key.compare(_alexandria[i].getType()) == 0)
        {
            break;
        }
    }

    if (i < _alexandria.size())
    {
        *atype  =  _alexandria[i];
        return 1;
    }
    else
    {
        return 0;
    }
}

double Poldata::elemGetMaxValence( std::string elem)
{
    double          mv = 0;
    unsigned int    i;

    for (i = 0; (i < _brule.size()); i++)
    {
        if ((0 == gmx_strcasecmp(_brule[i].getElem().c_str(), elem.c_str())) &&
            (mv < _brule[i].getValence()))
        {
            mv = _brule[i].getValence();
        }
    }
    return mv;
}

double *Poldata::elemGetBondorders( std::string elem1, std::string elem2,
                                    double distance, double toler)
{
    double           dev, *bo = NULL;
    std::string      ba1, ba2;
    unsigned int     i, j, k, nbo;

    if ((0 == elem1.size()) || (0 == elem2.size()))
    {
        return 0;
    }
    nbo = 0;
    for (i = 0; (i < _gtBond.size()); i++)
    {
        if (0 == _gtBond[i].getElem1().size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].getType().compare(_gtBond[i].getAtom2()) == 0)
                {
                    _gtBond[i].setElem1(_alexandria[j].getElem());
                }
            }
        }
        if (0 == _gtBond[i].getElem2().size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].getType().compare(_gtBond[i].getAtom2()) == 0)
                {
                    _gtBond[i].setElem2(_alexandria[j].getElem());
                }
            }
        }
        ba1 = _gtBond[i].getElem1();
        ba2 = _gtBond[i].getElem2();
        if (((ba1.compare(elem1) == 0) && (ba2.compare(elem2) == 0)) ||
            ((ba1.compare(elem2) == 0) && (ba2.compare(elem1) == 0)))
        {
            dev = fabs((_gtBond[i].getLength() - distance)/_gtBond[i].getLength());
            if (dev < toler)
            {
                for (k = 0; (k < nbo); k++)
                {
                    if (_gtBond[i].getBondorder() == bo[k])
                    {
                        break;
                    }
                }
                if (k == nbo)
                {
                    srenew(bo, nbo+2);
                    bo[nbo]   = _gtBond[i].getBondorder();
                    bo[nbo+1] = 0;
                    nbo++;
                }
            }
        }
    }
    return bo;
}

int Poldata::elemIsBond( std::string elem1, std::string elem2,
                         double distance, double toler)
{
    std::string     ba1, ba2;
    double          dev, devBest = 100000;
    unsigned int    j, i;

    if ((0 == elem1.size()) || (0 == elem2.size()))
    {
        return 0;
    }
    for (i = 0; (i < _gtBond.size()); i++)
    {
        if (0 == _gtBond[i].getElem1().size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].getType().compare(_gtBond[i].getAtom2()) == 0)
                {
                    _gtBond[i].setElem1(_alexandria[j].getElem());
                }
            }
        }
        if (0 == _gtBond[i].getElem2().size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].getType().compare(_gtBond[i].getAtom2()) == 0)
                {
                    _gtBond[i].setElem2(_alexandria[j].getElem());
                }
            }
        }
        ba1 = _gtBond[i].getElem1();
        ba2 = _gtBond[i].getElem2();
        if (((ba1.compare(elem1) == 0) && (ba2.compare(elem2)) == 0) ||
            ((ba1.compare(elem2) == 0) && (ba2.compare(elem1) == 0)))
        {
            dev = fabs((_gtBond[i].getLength() - distance)/_gtBond[i].getLength());
            if (dev < devBest)
            {
                devBest = dev;
            }
        }
    }
    return (devBest < toler);
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

GtBond *Poldata::searchBond( std::string atom1, std::string atom2,
                             double bondorder)
{
    GtBond    key, *gtB;
    int       i;

    key.setAtom1(atom1);
    key.setAtom2(atom2);
    key.setBondorder(bondorder);

    gtB = (GtBond *) bsearch(&key, vectorToArray(_gtBond), _gtBond.size(), sizeof(key), gtbComp);
    if (NULL != gtB)
    {
        i = indexOfPointInVector(gtB, _gtBond);
        while ((i > 0) && (loGtbComp(&(_gtBond[i-1]), &(_gtBond[i])) == 0))
        {
            i--;
        }
        gtB = &(_gtBond[i]);
    }
    return gtB;
}


double Poldata::atypeBondorder( std::string atype1, std::string atype2,
                                double distance, double toler)
{
    double              dev, devBest = 100000;
    unsigned int        i;
    int                 iBest = -1;
    GtBond             *gtB;

    if (0 == atype1.size() || 0 == atype2.size())
    {
        return 0.0;
    }
    gtB = searchBond( atype1, atype2, 0);
    if (NULL != gtB)
    {
        i = indexOfPointInVector(gtB, _gtBond);
        do
        {
            dev = fabs(_gtBond[i].getLength() - distance);
            if (dev < devBest)
            {
                devBest = dev;
                iBest   = i;
            }
            i++;
        }
        while ((i < _gtBond.size()) &&
               (0 == loGtbComp(&(_gtBond[i]), &(_gtBond[i-1]))));
    }
    if (devBest < toler)
    {
        return _gtBond[iBest].getBondorder();
    }

    return 0.0;
}

void Poldata::addMiller(std::string   miller,
                        int           atomnumber,
                        double        tauAhc,
                        double        alphaAhp)
{
    Miller mil(miller, atomnumber, tauAhc, alphaAhp);

    _miller.push_back(mil);
}

void Poldata::setMillerUnits( std::string tauUnit, std::string ahpUnit)
{
    _millerTauUnit = tauUnit;
    _millerAhpUnit = ahpUnit;
}

void Poldata::getMillerUnits( std::string *tauUnit,
                              std::string *ahpUnit)
{
    assignStr(tauUnit, _millerTauUnit);
    assignStr(ahpUnit, _millerAhpUnit);
}

int Poldata::getMillerPol(
        std::string   miller,
        int          *atomnumber,
        double       *tauAhc,
        double       *alphaAhp)
{
    Miller            *mil;
    unsigned int       i;

    for (i = 0; (i < _miller.size()); i++)
    {
        if (miller.compare(_miller[i].getMiller()) == 0)
        {
            mil = &(_miller[i]);
            assignScal(atomnumber, mil->getAtomnumber());
            assignScal(tauAhc, mil->getTauAhc());
            assignScal(alphaAhp, mil->getAlphaAhp());

            return 1;
        }
    }

    return 0;
}

void Poldata::addBosque(
        std::string   bosque,
        double        polarizability)
{
    Bosque bos(bosque, polarizability);
    _bosque.push_back(bos);
}

int Poldata::getBosquePol(
        std::string   bosque,
        double       *polarizability)
{
    unsigned int i;

    for (i = 0; (i < _bosque.size()); i++)
    {
        if (strcasecmp(bosque.c_str(), _bosque[i].getBosque().c_str()) == 0)
        {
            *polarizability = _bosque[i].getPolarizability();
            return 1;
        }
    }
    return 0;
}


int Poldata::searchBondtype( std::string atom)
{
    for (size_t j = 0; (j < _btype.size()); j++)
    {
        if (_btype[j].compare(atom) == 0)
        {
            return j;
        }
    }
    return -1;
}

/*
 * _gtBond stuff
 */
int Poldata::setBondParams( std::string atom1, std::string atom2,
                            double length, double sigma, int ntrain,
                            double bondorder, std::string params)
{
    GtBond *gtB;
    size_t  i;
    for (i = 0; (i < _gtBond.size()); i++)
    {
        gtB = &(_gtBond[i]);
        if ((((gtB->getAtom1().compare(atom1) == 0) &&
              (gtB->getAtom2().compare(atom2) == 0)) ||
             ((gtB->getAtom1().compare(atom2) == 0) &&
              (gtB->getAtom2().compare(atom1) == 0))) &&
            ((bondorder == 0) || (gtB->getBondorder() == bondorder)))
        {
            break;
        }
    }
    if (i < _gtBond.size())
    {
        if (length > 0)
        {
            gtB->setLength(length);
        }
        if (sigma > 0)
        {
            gtB->setSigma(sigma);
        }
        if (ntrain > 0)
        {
            gtB->setNtrain(ntrain);
        }
        gtB->setParams(params);
        return 1;
    }
    return 0;
}

FfatypeIterator Poldata::searchType(const std::string &type)
{
    return std::find_if(_alexandria.begin(), _alexandria.end(), 
                        [type](Ffatype const &f) { return f.getType().compare(type); });
}

FfatypeIterator Poldata::searchBtype(const std::string &btype)
{
    return std::find_if(_alexandria.begin(), _alexandria.end(), 
                        [btype](Ffatype const &f) { return f.getBtype().compare(btype); });
}

FfatypeIterator Poldata::searchPtype(const std::string &ptype)
{
    return std::find_if(_alexandria.begin(), _alexandria.end(), 
                        [ptype](Ffatype const &f) { return f.getPtype().compare(ptype); });
}

int Poldata::addBond( std::string atom1, std::string atom2,
                      double length, double sigma, int ntrain,
                      double bondorder, std::string params)
{
    FfatypeIterator a1, a2;

    if (((a1 = searchBtype(atom1)) == _alexandria.end()) ||
        ((a2 = searchBtype(atom2)) == _alexandria.end()))
    {
        return 0;
    }
    if (setBondParams( atom1, atom2, length, sigma, ntrain,
                       bondorder, params) == 0)
    {
        GtBond bond(atom1, atom2, params,
                    a1->getElem(), a2->getElem(),
                    length, sigma, bondorder, ntrain);

        _gtBond.push_back(bond);
        qsort(vectorToArray(_gtBond), _gtBond.size(), sizeof(_gtBond[0]), gtbComp);
    }
    return 1;
}


int Poldata::searchBond( std::string atom1, std::string atom2,
                         double *length, double *sigma, int *ntrain,
                         double *bondorder, std::string *params)
{
    GtBond *gtB;

    if ((0 == atom1.size()) || (0 == atom2.size()))
    {
        return 0;
    }
    gtB = searchBond( atom1, atom2, 0);
    if (NULL != gtB)
    {
        if (((gtB->getAtom1().compare(atom1) == 0) &&
             (gtB->getAtom2().compare(atom2) == 0)) ||
            ((gtB->getAtom1().compare(atom2) == 0) &&
             (gtB->getAtom2().compare(atom1) == 0)))
        {
            assignScal(length, gtB->getLength());
            assignScal(sigma, gtB->getSigma());
            assignScal(ntrain, gtB->getNtrain());
            assignScal(bondorder, gtB->getBondorder());
            assignStr(params, gtB->getParams());

            return 1+indexOfPointInVector(gtB, _gtBond);
        }
    }
    return 0;
}

/*
 * gt_angle stuff
 */
int Poldata::setAngleParams(std::string atom1,
                            std::string atom2,
                            std::string atom3, 
                            double angle, 
                            double sigma, 
                            int ntrain,
                            std::string params)
{
    FfatypeIterator a1, a2, a3;

    if (((a1 = searchBtype(atom1)) == _alexandria.end()) ||
        ((a2 = searchBtype(atom2)) == _alexandria.end()) ||
        ((a3 = searchBtype(atom3)) == _alexandria.end()))
    {
        return 0;
    }

    size_t i = 0;
    
    for (auto gtB : _gtAngle)
    {
        if ((gtB.getAtom2().compare(a2->getBtype()) == 0) &&
            (((gtB.getAtom1().compare(a1->getBtype()) == 0) &&
              (gtB.getAtom3().compare(a3->getBtype()) == 0)) ||
             ((gtB.getAtom1().compare(a3->getBtype()) == 0) &&
              (gtB.getAtom3().compare(a1->getBtype()) == 0))))
        {
            break;
        }
        i++;
    }
    if (i < _gtAngle.size())
    {
        if (angle > 0)
        {
            _gtAngle[i].setAngle(angle);
        }
        if (sigma > 0)
        {
            _gtAngle[i].setSigma(sigma);
        }
        if (ntrain > 0)
        {
            _gtAngle[i].setNtrain(ntrain);
        }
        _gtAngle[i].setParams(params);
        
        return 1;
    }
    return 0;
}

int Poldata::addAngle(std::string atom1, std::string atom2,
                      std::string atom3, double angle, double sigma,
                      int ntrain, std::string params)
{
    if (0 == setAngleParams(atom1, atom2, atom3,
                            angle, sigma, ntrain, params))
    {
        GtAngle ang(atom1, atom2, atom3,
                    params, angle, sigma, ntrain);

        _gtAngle.push_back(ang);
    }
    return 1;
}

int Poldata::searchAngle( std::string atom1, std::string atom2,
                          std::string atom3, double *angle, double *sigma,
                          int *ntrain, std::string *params)
{
    GtAngle             *gtB;
    unsigned int         i;

    if ((0 == atom1.size()) || (0 == atom2.size()) || (0 == atom3.size()))
    {
        return 0;
    }
    for (i = 0; (i < _gtAngle.size()); i++)
    {
        gtB = &(_gtAngle[i]);
        if ((gtB->getAtom2().compare(atom2) == 0) &&
            (((gtB->getAtom1().compare(atom1) == 0) &&
              (gtB->getAtom3().compare(atom3) == 0)) ||
             ((gtB->getAtom1().compare(atom3) == 0) &&
              (gtB->getAtom3().compare(atom1) == 0))))
        {
            assignScal(angle, gtB->getAngle());
            assignScal(sigma, gtB->getSigma());
            assignScal(ntrain, gtB->getNtrain());
            assignStr(params, gtB->getParams());

            return i+1;
        }
    }
    return 0;
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

DihedralIterator Poldata::searchDihedral(int egd,
                                         const std::string &atom1, 
                                         const std::string &atom2,
                                         const std::string &atom3, 
                                         const std::string &atom4)
{
    FfatypeIterator a1, a2, a3, a4;

    if (((a1 = searchBtype(atom1)) == _alexandria.end()) ||
        ((a2 = searchBtype(atom2)) == _alexandria.end()) ||
        ((a3 = searchBtype(atom3)) == _alexandria.end()) ||
        ((a4 = searchBtype(atom4)) == _alexandria.end()))
    {
        return _gtDihedral[egd].end();
    }

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

int Poldata::setDihedralParams(int egd,
                               std::string atom1, std::string atom2,
                               std::string atom3, std::string atom4,
                               double dihedral, double sigma, int ntrain,
                               std::string params)
{
    DihedralIterator gtB = searchDihedral(egd, atom1, atom2, atom3, atom4);
    if (_gtDihedral[egd].end() != gtB)
    {
        gtB->setDihedral(dihedral);
        if (sigma > 0)
        {
            gtB->setSigma(sigma);
        }
        if (ntrain > 0)
        {
            gtB->setNtrain(ntrain);
        }
        gtB->setParams(params);
        return 1;
    }
    return 0;
}

int Poldata::addDihedral( int egd,
                          std::string atom1, std::string atom2,
                          std::string atom3, std::string atom4, double dihedral,
                          double sigma, int ntrain, std::string params)
{
    if (0 == Poldata::setDihedralParams( egd, atom1, atom2,
                                         atom3, atom4, dihedral,
                                         sigma, ntrain, params))
    {
        GtDihedral dihed(atom1, atom2, atom3, atom4,
                         params, dihedral, sigma, ntrain);

        _gtDihedral[egd].push_back(dihed);
        std::sort(_gtDihedral[egd].begin(),
                  _gtDihedral[egd].end(), 
                  [](const GtDihedral &a, const GtDihedral &b)
                  { return a.compare(b); });
    }
    return 1;
}

int Poldata::searchDihedral( int egd,
                             std::string atom1, std::string atom2,
                             std::string atom3, std::string atom4,
                             double *dihedral, double *sigma,
                             int *ntrain, std::string *params)
{
    DihedralIterator gtRes = searchDihedral(egd, atom1, atom2, atom3, atom4);
    if (_gtDihedral[egd].end() != gtRes)
    {
        assignScal(dihedral, gtRes->getDihedral());
        assignScal(sigma, gtRes->getSigma());
        assignScal(ntrain, gtRes->getNtrain());
        assignStr(params, gtRes->getParams());

        return 1 + (gtRes - _gtDihedral[egd].begin());
    }
    return 0;
}

void Poldata::addSymcharges( std::string central,
                             std::string attached, int numattach)
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

int Poldata::searchSymcharges( std::string central,
                               std::string attached, int numattach)
{
    Symcharges            *sc;
    unsigned int           i;

    for (i = 0; (i < _symcharges.size()); i++)
    {
        sc = &(_symcharges[i]);
        if ((strcasecmp(sc->getCentral().c_str(), central.c_str()) == 0) &&
            (strcasecmp(sc->getAttached().c_str(), attached.c_str()) == 0) &&
            (sc->getNumattach() == numattach))
        {
            return 1;
        }
    }

    return 0;
}

/* Electrostatics properties */
Eemprops *Poldata::getEep(ChargeDistributionModel eqdModel,
                          const std::string       name)
{
    unsigned int i;

    for (i = 0; (i < _eep.size()); i++)
    {
        if ((strcasecmp(_eep[i].getName().c_str(), name.c_str()) == 0) &&
            (_eep[i].getEqdModel() == eqdModel))
        {
            return &(_eep[i]);
        }
    }
    return NULL;
}

void Poldata::setEemprops(ChargeDistributionModel eqdModel, const  std::string name,
                          double J0, double chi0, const std::string zeta, const std::string q, const std::string row)
{

    Eemprops                *eep;
    std::vector<std::string> sz, sq, sr;

    eep = getEep(eqdModel, name);
    if (NULL == eep)
    {
        _eep.resize(_eep.size()+1 );
        eep = &(_eep[_eep.size()-1]);
    }
    eep->setEqdModel(eqdModel);
    eep->setName(name);
    eep->setJ0(J0);
    sz                      = split(zeta, ' ');
    sq                      = split(q, ' ');
    sr                      = split(row, ' ');
    eep->setZetastr(zeta);
    eep->setQstr(q);
    eep->setRowstr(row);
    unsigned int nn = std::min(sz.size(), std::min(sq.size(), sr.size()));
    unsigned int n;
    for (n = 0; (n < nn); n++)
    {
        if (n < MAXZETA)
        {
            eep->setZeta(n, atof(sz[n].c_str()));
            eep->setQ(n, atof(sq[n].c_str()));
            eep->setRow(n, atoi(sr[n].c_str()));
        }
    }
    if (sz.size() > nn)
    {
        fprintf(stderr, "Warning: more zeta values than q/row values for %s n = %d\n",
                name.c_str(), nn);
    }
    if (sq.size() > nn)
    {
        fprintf(stderr, "Warning: more q values than zeta/row values for %s n = %d\n",
                name.c_str(), nn);
    }
    if (sr.size() > nn)
    {
        fprintf(stderr, "Warning: more row values than q/zeta values for %s n = %d\n",
                name.c_str(), nn);
    }
    eep->setNzeta(nn);
    if (nn >= MAXZETA)
    {
        fprintf(stderr, "More than %d zeta and/or q values for %s\n", MAXZETA, eep->getName().c_str());
        eep->setNzeta(MAXZETA);
    }
    for (; (n < MAXZETA); n++)
    {
        eep->setZeta(n, 0);
        eep->setQ(n, 0);
        eep->setRow(n, 0);
    }
    eep->setChi0(chi0);
}

int Poldata::getNumprops( ChargeDistributionModel eqdModel)
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

int Poldata::havePolSupport( const std::string atype)
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

int Poldata::haveEemSupport( ChargeDistributionModel eqdModel,
                             const std::string       name,
                             gmx_bool                bAllowZeroParameters)
{

    Eemprops  *eep  = getEep(eqdModel, name);

    return (eep && (bAllowZeroParameters || ((eep->getJ0() > 0) && (eep->getChi0() > 0))));
}

double Poldata::getJ00( ChargeDistributionModel eqdModel, const std::string name)
{
    Eemprops  *eer;

    if ((eer = getEep(eqdModel, name)) != NULL)
    {
        return eer->getJ0();
    }
    else
    {
        gmx_fatal(FARGS, "No J0 data for eqdModel %d and name %s",
                  eqdModel, name.c_str());
    }
    return -1;
}

std::string Poldata::getQstr( ChargeDistributionModel eqdModel, std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->getQstr();
    }
    return "";
}

std::string Poldata::getRowstr( ChargeDistributionModel eqdModel, std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->getRowstr();
    }
    return "";
}

int Poldata::getRow( ChargeDistributionModel eqdModel, const std::string name, int zz)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getRow(zz);
    }
    return -1;
}

double Poldata::getZeta( ChargeDistributionModel eqdModel, const std::string name, int zz)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
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

int Poldata::getNzeta( ChargeDistributionModel eqdModel, const std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->getNzeta();
    }
    return 0;
}

double Poldata::getQ( ChargeDistributionModel eqdModel, const std::string name, int zz)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        range_check(zz, 0, eer->getNzeta());
        return eer->getQ(zz);
    }
    return -1;
}

double Poldata::getChi0( ChargeDistributionModel eqdModel, const  std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->getChi0();
    }
    else
    {
        gmx_fatal(FARGS, "No chi0 data for eqdModel %d and name %s", eqdModel, name.c_str());
    }
    return -1;
}

void Poldata::setEpref( ChargeDistributionModel eqdModel, std::string epref)
{
    unsigned int i;

    for (i = 0; (i < _epr.size()); i++)
    {
        if (_epr[i].getEqdModel() == eqdModel)
        {
            _epr[i].setEpref(epref);
            break;
        }
    }
    if (i == _epr.size())
    {
        Epref epr(eqdModel, epref);
        _epr.push_back(epr);
    }
}

std::string Poldata::getEpref( ChargeDistributionModel eqdModel)
{
    unsigned int i;

    for (i = 0; (i < _epr.size()); i++)
    {
        if (_epr[i].getEqdModel() == eqdModel)
        {
            return _epr[i].getEpref();
        }
    }
    return "";
}

void Poldata::commEemprops( t_commrec *cr)
{
    unsigned int          i, j, nep;
    std::vector<Eemprops> ep;

    if (NULL != debug)
    {
        fprintf(debug, "Going to update eemprops on node %d\n", cr->nodeid);
    }
    if (MASTER(cr))
    {
        for (i = 1; ((int)i < cr->nnodes); i++)
        {
            gmx_send_int(cr, i, _eep.size());
            gmx_send(cr, i, vectorToArray(_eep), _eep.size()*sizeof(_eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (nep != _eep.size())
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        ep.resize(_eep.size());
        gmx_recv(cr, 0, vectorToArray(ep), _eep.size()*sizeof(ep[0]));
        for (i = 0; (i < _eep.size()); i++)
        {
            _eep[i] = ep[i];
        }
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < nep); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    getEemtypeName(_eep[i].getEqdModel()).c_str(),
                    _eep[i].getName().c_str(), _eep[i].getChi0(),
                    _eep[i].getJ0());
            for (j = 0; ((int)j < _eep[i].getNzeta()); j++)
            {
                fprintf(debug, " %8.3f", _eep[i].getZeta(j));
            }
            fprintf(debug, "\n");
        }
    }
}

void Poldata::commForceParameters(t_commrec *cr)
{
    unsigned int          i, j, nep;
    std::vector<Eemprops> ep;

    if (NULL != debug)
    {
        fprintf(debug, "Going to update force parameters on node %d\n", cr->nodeid);
    }
    if (MASTER(cr))
    {
        for (i = 1; ((int)i < cr->nnodes); i++)
        {
            gmx_send_int(cr, i, _eep.size());
            gmx_send(cr, i, vectorToArray(_eep), _eep.size()*sizeof(_eep[0]));
        }
    }
    else
    {
        nep = gmx_recv_int(cr, 0);
        if (_eep.size() != nep)
        {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
        }
        ep.resize(nep);
        gmx_recv(cr, 0, vectorToArray(ep), _eep.size()*sizeof(ep[0]));
        for (i = 0; (i < _eep.size()); i++)
        {
            _eep[i] = ep[i];
        }
    }
    if (NULL != debug)
    {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < _eep.size()); i++)
        {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    getEemtypeName(_eep[i].getEqdModel()).c_str(),
                    _eep[i].getName().c_str(), _eep[i].getChi0(),
                    _eep[i].getJ0());
            for (j = 0; ((int)j < _eep[i].getNzeta()); j++)
            {
                fprintf(debug, " %8.3f", _eep[i].getZeta(j));
            }
            fprintf(debug, "\n");
        }
    }
}

typedef struct {
    ChargeDistributionModel eqd;
    const char             *name, *ref;
    gmx_bool                bWeight;
} t_eemtype_props;

t_eemtype_props eemtype_props[eqdNR] = {
    { eqdAXp,      "AXp",      "Maaren2014a",   FALSE },
    { eqdAXg,      "AXg",      "Maaren2014a",   TRUE },
    { eqdAXs,      "AXs",      "Maaren2014a",   TRUE },
    { eqdYang,     "Yang",     "Yang2006b",     TRUE },
    { eqdBultinck, "Bultinck", "Bultinck2002a", FALSE },
    { eqdRappe,    "Rappe",    "Rappe1991a",    TRUE }
};

ChargeDistributionModel Poldata::name2eemtype(const std::string name)
{
    unsigned int i;

    for (i = 0; (i < eqdNR); i++)
    {
        if (strcasecmp(name.c_str(), eemtype_props[i].name) == 0)
        {
            return eemtype_props[i].eqd;
        }
    }
    return eqdNR;
}

std::string Poldata::getEemtypeName(ChargeDistributionModel eem)
{
    unsigned int i;

    for (i = 0; (i < eqdNR); i++)
    {
        if (eem == eemtype_props[i].eqd)
        {
            return eemtype_props[i].name;
        }
    }

    return "";
}

}
