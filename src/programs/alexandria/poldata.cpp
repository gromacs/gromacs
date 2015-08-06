/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITN
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/cstringutil.h"
#include "poldata.h"
#include "gmx_simple_comm.h"
#include "stringutil.h"

#define assignStr(dst, src) if (dst != NULL) {*dst = src; }

#define assignScal(dst, src) if (dst) *dst = src

#define vectorToArray(vec) (&(vec[0]))

namespace alexandria
{

Poldata::Poldata()
    :
      _gtDihedralFunction(egdNR),
      _ngtDihedralC(egdNR),
      _gtDihedralFtype(egdNR),
      _gtDihedral(egdNR)
{


    _nptypeC      = 0;
    _nalexandriaC = 0;
    _nbruleC      = 0;
    _nexcl        = 0;
    _gtBondFtype  = 0;
    _gtAngleFtype = 0;
    _nmillerC     = 0;
    _nbosqueC     = 0;
    _nsymchargesC = 0;
    _nepC         = 0;
    _nerC         = 0;

    for (int x = 0; x < egdNR; x++)
    {
        _ngtDihedralC[x] = 0;
        // _gtDihedral[x]=NULL;
    }

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

void Poldata::addPtype(
        const std::string ptype,
        const std::string miller,
        const std::string bosque,
        double            polarizability,
        double            sigPol)
{

    unsigned int      i;

    for (i = 0; (i < _ptype.size()); i++)
    {
        if (_ptype[i].type.compare(ptype) == 0)
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

void Poldata::addBtype(
        const std::string btype)
{
    unsigned int i;

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

void Poldata::addAtype(
        const std::string elem,
        const std::string desc,
        const std::string atype,
        const std::string ptype,
        const std::string btype,
        const std::string vdwparams,
        double            refEnthalpy)
{

    unsigned int        i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (_alexandria[i].type.compare(atype) == 0)
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
        if (_alexandria[j].type.compare(atype) == 0)
        {
            break;
        }
    }
    if (j < _alexandria.size())
    {
        for (i = 0; (i < _brule.size()); i++)
        {
            if (_brule[i].rule.compare( gtBrule) == 0)
            {
                break;
            }
        }
        if (i == _brule.size())
        {

            Brule brule(_alexandria[j].elem, gtBrule, atype, neighbors, geometry,
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

int Poldata::getBondingRule(
        std::string *gt_brule, std::string *atype,
        std::string *geometry, int *numbonds,
        double *valence, int *iAromatic,
        std::string *neighbors)
{
    if (_nbruleC < _brule.size())
    {
        assignStr(gt_brule, _brule[_nbruleC].rule);
        assignStr(atype, _brule[_nbruleC].type);
        assignStr(geometry, _brule[_nbruleC].geometry);
        assignScal(numbonds, _brule[_nbruleC].numbonds);
        assignScal(valence, _brule[_nbruleC].valence);
        assignScal(iAromatic, _brule[_nbruleC].iAromatic);
        assignStr(neighbors, _brule[_nbruleC].neighbors);

        _nbruleC++;

        return 1;
    }
    return 0;
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
        if (ptype.compare(_ptype[i].type) == 0)
        {
            sp                 = &(_ptype[i]);
            sp->polarizability = polarizability;
            sp->sigPol         = sigPol;
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
            if (_brule[i].rule.compare(gtBrule) == 0)
            {
                return _brule[i].geometry;
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
            if (_alexandria[i].type.compare(atype) == 0)
            {
                return _alexandria[i].desc;
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
    for (unsigned int i = 0; (i < brule->nb.size()); i++)
    {
        IFound = 0;
        for (j = 0; (j < nbond); j++)
        {
            if (
                (0 != nbhybrid[j].size()) &&
                (jj[j] == 0) &&
                (IFound == 0) &&
                strcasestrStart(brule->nb[i], nbhybrid[j])
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

std::string *Poldata::getBondingRules( std::string elem,
                                       int nbond, std::string neighbors[],
                                       const std::string geometry,
                                       int iAromatic)
{
    unsigned int             nnb;
    unsigned int             i;
    int nptr                     = 0, best = -1, score;
    std::string             *ptr = NULL;

    for (i = 0; (i < _brule.size()); i++)
    {
        nnb = countNeighbors(&(_brule[i]), nbond, neighbors, &score);
        if ((_brule[i].elem.compare(elem) == 0) &&
            (strcasecmp(_brule[i].geometry.c_str(), geometry.c_str()) == 0) &&
            (nbond == _brule[i].numbonds) &&
            ((iAromatic >= 0 && (iAromatic == _brule[i].iAromatic)) ||
             (iAromatic < 0)) &&
            (nnb == _brule[i].nb.size()))
        {
            if (score > best)
            {
                if (NULL != ptr)
                {
                    sfree(ptr);
                    nptr = 0;
                    ptr  = NULL;
                }
            }
            if (score >= best)
            {
                srenew(ptr, ++nptr);
                ptr[nptr-1] = _brule[i].rule;
                best        = score;
            }
        }
    }
    srenew(ptr, ++nptr);
    ptr[nptr-1] = "";

    return ptr;
}

int Poldata::bondingRuleValence( std::string gtBrule, double *valence)
{
    unsigned int i;

    for (i = 0; (i < _brule.size()); i++)
    {
        if (strcasecmp(gtBrule.c_str(), _brule[i].rule.c_str()) == 0)
        {
            *valence = _brule[i].valence;
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
        if (ptype.compare(_ptype[j].type) == 0)
        {
            if (NULL != polar)
            {
                *polar   = _ptype[j].polarizability;
            }

            if (NULL != sigPol)
            {
                *sigPol = _ptype[j].sigPol;
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
        if (atype.compare(_alexandria[i].type) == 0)
        {
            return getPtypePol( _alexandria[i].ptype, polar, sigPol);
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
        if (atype.compare(_alexandria[i].type) == 0)
        {
            *Href = _alexandria[i].refEnthalpy;
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
        if (ptype.compare(_ptype[i].type) == 0)
        {
            return _ptype[i].miller;
        }
    }
    return "";
}

std::string Poldata::ptypeToBosque( const std::string ptype)
{
    unsigned int i;

    for (i = 0; (i < _ptype.size()); i++)
    {
        if (ptype.compare(_ptype[i].type) == 0)
        {
            return _ptype[i].bosque;
        }
    }
    return "";
}

int Poldata::getPtype(
        std::string        *ptype,
        std::string        *miller,
        std::string        *bosque,
        double             *polarizability,
        double             *sigPol)
{
    Ptype *sp;

    if (_nptypeC < _ptype.size())
    {
        sp = &(_ptype[_nptypeC]);
        assignScal(polarizability, sp->polarizability);
        assignScal(sigPol, sp->sigPol);
        assignStr(ptype, sp->type);
        assignStr(miller, sp->miller);
        assignStr(bosque, sp->bosque);
        _nptypeC++;
        return 1;
    }
    else
    {
        _nptypeC = 0;
    }

    return 0;
}

int Poldata::getAtype(
        std::string        *elem,
        std::string        *desc,
        std::string        *atype,
        std::string        *ptype,
        std::string        *btype,
        std::string        *vdwparams,
        double             *refEnthalpy)
{
    Ffatype *sp;

    if (_nalexandriaC < _alexandria.size())
    {
        sp = &(_alexandria[_nalexandriaC]);
        assignStr(elem, sp->elem);
        assignStr(desc, sp->desc);
        assignStr(atype, sp->type);
        assignStr(ptype, sp->ptype);
        assignStr(btype, sp->btype);
        assignStr(vdwparams, sp->vdwparams);
        *refEnthalpy = sp->refEnthalpy;
        _nalexandriaC++;
        return 1;
    }
    else
    {
        _nalexandriaC = 0;
    }

    return 0;
}

std::string Poldata::atypeToPtype( const std::string atype)
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (_alexandria[i].type.compare(atype) == 0)
        {
            return _alexandria[i].ptype;
        }
    }
    return "";
}

std::string Poldata::atypeToBtype( const std::string atype)
{
    unsigned int i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (_alexandria[i].type.compare(atype) == 0)
        {
            return _alexandria[i].btype;
        }
    }
    return "";
}

int Poldata::searchAtype(
        std::string         key,
        std::string        *elem,
        std::string        *desc,
        std::string        *atype,
        std::string        *ptype,
        std::string        *btype,
        std::string        *vdwparams)
{
    Ffatype            *sp;
    unsigned int        i;

    for (i = 0; (i < _alexandria.size()); i++)
    {
        if (key.compare(_alexandria[i].type) == 0)
        {
            break;
        }
    }

    if (i < _alexandria.size())
    {
        sp = &(_alexandria[i]);
        assignStr(elem, sp->elem);
        assignStr(desc, sp->desc);
        assignStr(atype, sp->type);
        assignStr(ptype, sp->ptype);
        assignStr(btype, sp->btype);
        assignStr(vdwparams, sp->vdwparams);

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
        if ((0 == gmx_strcasecmp(_brule[i].elem.c_str(), elem.c_str())) &&
            (mv < _brule[i].valence))
        {
            mv = _brule[i].valence;
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
        if (0 == _gtBond[i].elem1.size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].type.compare(_gtBond[i].atom2) == 0)
                {
                    _gtBond[i].elem1 = _alexandria[j].elem;
                }
            }
        }
        if (0 == _gtBond[i].elem2.size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].type.compare(_gtBond[i].atom2) == 0)
                {
                    _gtBond[i].elem2 = _alexandria[j].elem;
                }
            }
        }
        ba1 = _gtBond[i].elem1;
        ba2 = _gtBond[i].elem2;
        if (((ba1.compare(elem1) == 0) && (ba2.compare(elem2) == 0)) ||
            ((ba1.compare(elem2) == 0) && (ba2.compare(elem1) == 0)))
        {
            dev = fabs((_gtBond[i].length - distance)/_gtBond[i].length);
            if (dev < toler)
            {
                for (k = 0; (k < nbo); k++)
                {
                    if (_gtBond[i].bondorder == bo[k])
                    {
                        break;
                    }
                }
                if (k == nbo)
                {
                    srenew(bo, nbo+2);
                    bo[nbo]   = _gtBond[i].bondorder;
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
        if (0 == _gtBond[i].elem1.size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].type.compare(_gtBond[i].atom2) == 0)
                {
                    _gtBond[i].elem1 = _alexandria[j].elem;
                }
            }
        }
        if (0 == _gtBond[i].elem2.size())
        {
            for (j = 0; (j < _alexandria.size()); j++)
            {
                if (_alexandria[j].type.compare(_gtBond[i].atom2) == 0)
                {
                    _gtBond[i].elem2 =  _alexandria[j].elem;
                }
            }
        }
        ba1 = _gtBond[i].elem1;
        ba2 = _gtBond[i].elem2;
        if (((ba1.compare(elem1) == 0) && (ba2.compare(elem2)) == 0) ||
            ((ba1.compare(elem2) == 0) && (ba2.compare(elem1) == 0)))
        {
            dev = fabs((_gtBond[i].length - distance)/_gtBond[i].length);
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

    if (ba->atom1.compare(ba->atom2) <= 0)
    {
        a1 = ba->atom1;
        a2 = ba->atom2;
    }
    else
    {
        a2 = ba->atom1;
        a1 = ba->atom2;
    }
    if (bb->atom1.compare(bb->atom2) <= 0)
    {
        b1 = bb->atom1;
        b2 = bb->atom2;
    }
    else
    {
        b2 = bb->atom1;
        b1 = bb->atom2;
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
    if ((0 == i) && ((0 != ba->bondorder) && (0 != bb->bondorder)))
    {
        if (ba->bondorder < bb->bondorder)
        {
            i = -1;
        }
        else if (ba->bondorder > bb->bondorder)
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

    key.atom1     = atom1;
    key.atom2     = atom2;
    key.bondorder = bondorder;

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
            dev = fabs(_gtBond[i].length - distance);
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
        return _gtBond[iBest].bondorder;
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

int Poldata::getMiller(
        std::string  *miller,
        int          *atomnumber,
        double       *tauAhc,
        double       *alphaAhp)
{
    Miller            *mil;
    unsigned int       i;

    i = _nmillerC;

    if (i < _miller.size())
    {
        mil = &(_miller[i]);
        assignStr(miller, mil->miller);
        assignScal(atomnumber, mil->atomnumber);
        assignScal(tauAhc, mil->tauAhc);
        assignScal(alphaAhp, mil->alphaAhp);
        _nmillerC++;

        return 1;
    }
    else
    {
        _nmillerC = 0;
    }

    return 0;
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
        if (miller.compare(_miller[i].miller) == 0)
        {
            mil = &(_miller[i]);
            assignScal(atomnumber, mil->atomnumber);
            assignScal(tauAhc, mil->tauAhc);
            assignScal(alphaAhp, mil->alphaAhp);

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

int Poldata::getBosque(
        std::string  *bosque,
        double       *polarizability)
{
    if (_nbosqueC < _bosque.size())
    {
        assignStr(bosque, _bosque[_nbosqueC].bosque);
        assignScal(polarizability, _bosque[_nbosqueC].polarizability);
        _nbosqueC++;

        return 1;
    }
    else
    {
        _nbosqueC = 0;
    }

    return 0;
}

int Poldata::getBosquePol(
        std::string   bosque,
        double       *polarizability)
{
    unsigned int i;

    for (i = 0; (i < _bosque.size()); i++)
    {
        if (strcasecmp(bosque.c_str(), _bosque[i].bosque.c_str()) == 0)
        {
            *polarizability = _bosque[i].polarizability;
            return 1;
        }
    }
    return 0;
}


int Poldata::searchBondtype( std::string atom)
{
    unsigned int j;

    for (j = 0; (j < _btype.size()); j++)
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
    GtBond             *gtB;
    unsigned int        i;

    for (i = 0; (i < _gtBond.size()); i++)
    {
        gtB = &(_gtBond[i]);
        if ((((gtB->atom1.compare(atom1) == 0) &&
              (gtB->atom2.compare(atom2) == 0)) ||
             ((gtB->atom1.compare(atom2) == 0) &&
              (gtB->atom2.compare(atom1) == 0))) &&
            ((bondorder == 0) || (gtB->bondorder == bondorder)))
        {
            break;
        }
    }
    if (i < _gtBond.size())
    {
        if (length > 0)
        {
            gtB->length = length;
        }
        if (sigma > 0)
        {
            gtB->sigma = sigma;
        }
        if (ntrain > 0)
        {
            gtB->ntrain = ntrain;
        }
        gtB->params = params;
        return 1;
    }
    return 0;
}

int Poldata::addBond( std::string atom1, std::string atom2,
                      double length, double sigma, int ntrain,
                      double bondorder, std::string params)
{
    int        a1, a2;

    if (-1 == (a1 = searchBondtype( atom1)))
    {
        return 0;
    }
    if (-1 == (a2 = searchBondtype( atom2)))
    {
        return 0;
    }
    if (setBondParams( atom1, atom2, length, sigma, ntrain,
                       bondorder, params) == 0)
    {
        GtBond bond(atom1, atom2, params,
                    _alexandria[a1].elem, _alexandria[a2].elem,
                    length, sigma, bondorder, ntrain);

        _gtBond.push_back(bond);
        qsort(vectorToArray(_gtBond), _gtBond.size(), sizeof(_gtBond[0]), gtbComp);
    }
    return 1;
}

int Poldata::getBond( std::string *atom1, std::string *atom2,
                      double *length, double *sigma, int *ntrain,
                      double *bondorder, std::string *params)
{
    GtBond *gtB;

    if (_ngtBondC < _gtBond.size())
    {
        gtB = &(_gtBond[_ngtBondC]);
        assignStr(atom1, gtB->atom1);
        assignStr(atom2, gtB->atom2);
        assignScal(length, gtB->length);
        assignScal(sigma, gtB->sigma);
        assignScal(ntrain, gtB->ntrain);
        assignScal(bondorder, gtB->bondorder);
        assignStr(params, gtB->params);
        _ngtBondC++;

        return _ngtBondC;
    }
    _ngtBondC = 0;

    return 0;
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
        if (((gtB->atom1.compare(atom1) == 0) &&
             (gtB->atom2.compare(atom2) == 0)) ||
            ((gtB->atom1.compare(atom2) == 0) &&
             (gtB->atom2.compare(atom1) == 0)))
        {
            assignScal(length, gtB->length);
            assignScal(sigma, gtB->sigma);
            assignScal(ntrain, gtB->ntrain);
            assignScal(bondorder, gtB->bondorder);
            assignStr(params, gtB->params);

            return 1+indexOfPointInVector(gtB, _gtBond);
        }
    }
    return 0;
}

/*
 * gt_angle stuff
 */
int Poldata::setAngleParams( std::string atom1, std::string atom2,
                             std::string atom3, double angle, double sigma, int ntrain,
                             std::string params)
{
    GtAngle             *gtB;
    unsigned int         i;

    for (i = 0; (i < _gtAngle.size()); i++)
    {
        gtB = &(_gtAngle[i]);
        if ((gtB->atom2.compare(atom2) == 0) &&
            (((gtB->atom1.compare(atom1) == 0) &&
              (gtB->atom3.compare(atom3) == 0)) ||
             ((gtB->atom1.compare(atom3) == 0) &&
              (gtB->atom3.compare(atom1) == 0))))
        {
            break;
        }
    }
    if (i < _gtAngle.size())
    {
        if (angle > 0)
        {
            gtB->angle = angle;
        }
        if (sigma > 0)
        {
            gtB->sigma = sigma;
        }
        if (ntrain > 0)
        {
            gtB->ntrain = ntrain;
        }
        gtB->params = params;
        return 1;
    }
    return 0;
}

int Poldata::addAngle(std::string atom1, std::string atom2,
                      std::string atom3, double angle, double sigma,
                      int ntrain, std::string params)
{

    if ((-1 == searchBondtype( atom1)) ||
        (-1 == searchBondtype( atom2)) ||
        (-1 == searchBondtype( atom3)))
    {
        return 0;
    }

    if (0 == setAngleParams( atom1, atom2, atom3, angle, sigma, ntrain, params))
    {
        GtAngle ang(atom1, atom2, atom3,
                    params, angle, sigma, ntrain);

        _gtAngle.push_back(ang);

    }
    return 1;
}

int Poldata::getAngle( std::string *atom1, std::string *atom2,
                       std::string *atom3, double *angle, double *sigma,
                       int *ntrain, std::string *params)
{
    GtAngle *gtB;

    if (_ngtAngleC < _gtAngle.size())
    {
        gtB = &(_gtAngle[_ngtAngleC]);
        assignStr(atom1, gtB->atom1);
        assignStr(atom2, gtB->atom2);
        assignStr(atom3, gtB->atom3);
        assignScal(angle, gtB->angle);
        assignScal(sigma, gtB->sigma);
        assignScal(ntrain, gtB->ntrain);
        assignStr(params, gtB->params);
        _ngtAngleC++;

        return _ngtAngleC;
    }
    _ngtAngleC = 0;

    return 0;
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
    for (i = 0; (i < _ngtAngleC); i++)
    {
        gtB = &(_gtAngle[i]);
        if ((gtB->atom2.compare(atom2) == 0) &&
            (((gtB->atom1.compare(atom1) == 0) &&
              (gtB->atom3.compare(atom3) == 0)) ||
             ((gtB->atom1.compare(atom3) == 0) &&
              (gtB->atom3.compare(atom1) == 0))))
        {
            assignScal(angle, gtB->angle);
            assignScal(sigma, gtB->sigma);
            assignScal(ntrain, gtB->ntrain);
            assignStr(params, gtB->params);

            return i+1;
        }
    }
    return 0;
}


/*
 * gt_dihedral stuff
 */
int Poldata::gtdComp(const void *a, const void *b)
{
    GtDihedral    *gtA = (GtDihedral *)a;
    GtDihedral    *gtB = (GtDihedral *)b;
    int            n;

    if (0 == (n = gtA->atom1.compare(gtB->atom1)))
    {
        if (0 == (n = gtA->atom2.compare(gtB->atom2)))
        {
            if (0 == (n = gtA->atom3.compare(gtB->atom3)))
            {
                n = gtA->atom4.compare(gtB->atom4);
            }
        }
    }

    return n;
}

GtDihedral *Poldata::searchDihedral( int egd,
                                     std::string atom1, std::string atom2,
                                     std::string atom3, std::string atom4)
{
    GtDihedral    gtA, *gtRes, *gtDptr;
    int           nd;

    if ((0 == atom1.size()) || (0 == atom2.size()) || (0 == atom3.size()) || (0 == atom4.size()))
    {
        return NULL;
    }
    gtDptr     = vectorToArray(_gtDihedral[egd]);
    nd         = _gtDihedral[egd].size();
    gtA.atom1  = atom1;
    gtA.atom2  = atom2;
    gtA.atom3  = atom3;
    gtA.atom4  = atom4;
    gtRes      = (GtDihedral *) bsearch(&gtA, gtDptr, nd, sizeof(gtA), &gtdComp);
    if (NULL == gtRes)
    {
        gtA.atom1 = atom4;
        gtA.atom2 = atom3;
        gtA.atom3 = atom2;
        gtA.atom4 = atom1;
        gtRes     = (GtDihedral *) bsearch(&gtA, gtDptr, nd, sizeof(gtA), gtdComp);
    }
    return gtRes;
}

int Poldata::setDihedralParams( int egd,
                                std::string atom1, std::string atom2,
                                std::string atom3, std::string atom4,
                                double dihedral, double sigma, int ntrain,
                                std::string params)
{
    GtDihedral *gtB;


    gtB = searchDihedral( egd, atom1, atom2, atom3, atom4);
    if (NULL != gtB)
    {
        if (dihedral > 0)
        {
            gtB->dihedral = dihedral;
        }
        if (sigma > 0)
        {
            gtB->sigma = sigma;
        }
        if (ntrain > 0)
        {
            gtB->ntrain = ntrain;
        }
        gtB->params = params;
        return 1;
    }
    return 0;
}

int Poldata::addDihedral( int egd,
                          std::string atom1, std::string atom2,
                          std::string atom3, std::string atom4, double dihedral,
                          double sigma, int ntrain, std::string params)
{

    if ((-1 == searchBondtype( atom1)) ||
        (-1 == searchBondtype( atom2)) ||
        (-1 == searchBondtype( atom3)) ||
        (-1 == searchBondtype( atom4)))
    {
        return 0;
    }

    if (0 == Poldata::setDihedralParams( egd, atom1, atom2,
                                         atom3, atom4, dihedral,
                                         sigma, ntrain, params))
    {

        GtDihedral dihed(atom1, atom2, atom3, atom4,
                         params, dihedral, sigma, ntrain);

        _gtDihedral[egd].push_back(dihed);
        qsort(vectorToArray(_gtDihedral[egd]), _gtDihedral[egd].size(), sizeof(_gtDihedral[egd][0]),
              gtdComp);
    }
    return 1;
}

int Poldata::getDihedral( int egd,
                          std::string *atom1, std::string *atom2,
                          std::string *atom3, std::string *atom4, double *dihedral,
                          double *sigma, int *ntrain, std::string *params)
{
    GtDihedral *gtB;

    if (_ngtDihedralC[egd] < _gtDihedral[egd].size())
    {
        gtB = &(_gtDihedral[egd][_ngtDihedralC[egd]]);
        assignStr(atom1, gtB->atom1);
        assignStr(atom2, gtB->atom2);
        assignStr(atom3, gtB->atom3);
        assignStr(atom4, gtB->atom4);
        assignScal(dihedral, gtB->dihedral);
        assignScal(sigma, gtB->sigma);
        assignScal(ntrain, gtB->ntrain);
        assignStr(params, gtB->params);
        _ngtDihedralC[egd]++;

        return _ngtDihedralC[egd];
    }
    _ngtDihedralC[egd] = 0;

    return 0;
}

int Poldata::searchDihedral( int egd,
                             std::string atom1, std::string atom2,
                             std::string atom3, std::string atom4,
                             double *dihedral, double *sigma,
                             int *ntrain, std::string *params)
{
    GtDihedral *gtRes;

    gtRes = searchDihedral( egd, atom1, atom2, atom3, atom4);
    if (NULL != gtRes)
    {
        assignScal(dihedral, gtRes->dihedral);
        assignScal(sigma, gtRes->sigma);
        assignScal(ntrain, gtRes->ntrain);
        assignStr(params, gtRes->params);

        return 1 + (indexOfPointInVector(gtRes, _gtDihedral[egd]));
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
        if ((strcasecmp(sc->central.c_str(), central.c_str()) == 0) &&
            (strcasecmp(sc->attached.c_str(), attached.c_str()) == 0) &&
            (sc->numattach == numattach))
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

int Poldata::getSymcharges( std::string *central,
                            std::string *attached, int *numattach)
{
    Symcharges *sc;

    if (_nsymchargesC < _symcharges.size())
    {
        sc = &(_symcharges[_nsymchargesC]);
        assignStr(central, sc->central);
        assignStr(attached, sc->attached);
        assignScal(numattach, sc->numattach);
        _nsymchargesC++;

        return 1;
    }
    _nsymchargesC = 0;

    return 0;
}

int Poldata::searchSymcharges( std::string central,
                               std::string attached, int numattach)
{
    Symcharges            *sc;
    unsigned int           i;

    for (i = 0; (i < _symcharges.size()); i++)
    {
        sc = &(_symcharges[i]);
        if ((strcasecmp(sc->central.c_str(), central.c_str()) == 0) &&
            (strcasecmp(sc->attached.c_str(), attached.c_str()) == 0) &&
            (sc->numattach == numattach))
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
        if ((strcasecmp(_eep[i].name.c_str(), name.c_str()) == 0) &&
            (_eep[i].eqdModel == eqdModel))
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
    eep->eqdModel           = eqdModel;
    eep->name               = name;
    eep->J0                 = J0;
    sz                      = split(zeta, ' ');
    sq                      = split(q, ' ');
    sr                      = split(row, ' ');
    eep->zetastr            = zeta;
    eep->qstr               = q;
    eep->rowstr             = row;
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
    eep->nzeta = nn;
    if (nn >= MAXZETA)
    {
        fprintf(stderr, "More than %d zeta and/or q values for %s\n", MAXZETA, eep->name.c_str());
        eep->nzeta = MAXZETA;
    }
    for (; (n < MAXZETA); n++)
    {
        eep->setZeta(n, 0);
        eep->setQ(n, 0);
        eep->setRow(n, 0);
    }
    eep->chi0  = chi0;
}

int Poldata::getEemprops(
        ChargeDistributionModel *eqdModel, std::string *name,
        double *J0, double *chi0, std::string *zeta, std::string *q, std::string *row)
{
    if (_nepC < _eep.size())
    {
        assignScal(eqdModel, _eep[_nepC].eqdModel);
        assignStr(name, _eep[_nepC].name);
        assignScal(J0, _eep[_nepC].J0);
        assignStr(zeta, _eep[_nepC].zetastr);
        assignStr(q, _eep[_nepC].qstr);
        assignStr(row, _eep[_nepC].rowstr);
        assignScal(chi0, _eep[_nepC].chi0);
        _nepC++;
        return 1;
    }
    else
    {
        _nepC = 0;
        return 0;
    }
}

int Poldata::getNumprops( ChargeDistributionModel eqdModel)
{
    unsigned int i, n = 0;

    for (i = 0; (i < _eep.size()); i++)
    {
        if (_eep[i].eqdModel == eqdModel)
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
        if (atype.compare(_alexandria[i].type) == 0)
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

    return (eep && (bAllowZeroParameters || ((eep->J0 > 0) && (eep->chi0 > 0))));
}

double Poldata::getJ00( ChargeDistributionModel eqdModel, const std::string name)
{
    Eemprops  *eer;

    if ((eer = getEep(eqdModel, name)) != NULL)
    {
        return eer->J0;
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
        return eer->qstr;
    }
    return "";
}

std::string Poldata::getRowstr( ChargeDistributionModel eqdModel, std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->rowstr;
    }
    return "";
}

int Poldata::getRow( ChargeDistributionModel eqdModel, const std::string name, int zz)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        range_check(zz, 0, eer->nzeta);
        return eer->getRow(zz);
    }
    return -1;
}

double Poldata::getZeta( ChargeDistributionModel eqdModel, const std::string name, int zz)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        if ((zz < 0) || (zz >= eer->nzeta))
        {
            printf("Bleh\n");
        }
        range_check(zz, 0, eer->nzeta);
        return eer->getZeta(zz);
    }
    return -1;
}

int Poldata::getNzeta( ChargeDistributionModel eqdModel, const std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->nzeta;
    }
    return 0;
}

double Poldata::getQ( ChargeDistributionModel eqdModel, const std::string name, int zz)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        range_check(zz, 0, eer->nzeta);
        return eer->getQ(zz);
    }
    return -1;
}

double Poldata::getChi0( ChargeDistributionModel eqdModel, const  std::string name)
{
    Eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
    {
        return eer->chi0;
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
        if (_epr[i].eqdModel == eqdModel)
        {
            _epr[i].epref = epref;
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
        if (_epr[i].eqdModel == eqdModel)
        {
            return _epr[i].epref;
        }
    }
    return "";
}

int Poldata::listEpref( ChargeDistributionModel *eqdModel, std::string *epref)
{
    if (_nerC < _epr.size())
    {
        assignScal(eqdModel, _epr[_nerC].eqdModel);
        assignStr(epref, _epr[_nerC].epref);
        _nerC++;
        return 1;
    }
    _nerC = 0;

    return 0;
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
                    getEemtypeName(_eep[i].eqdModel).c_str(),
                    _eep[i].name.c_str(), _eep[i].chi0,
                    _eep[i].J0);
            for (j = 0; ((int)j < _eep[i].nzeta); j++)
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
                    getEemtypeName(_eep[i].eqdModel).c_str(),
                    _eep[i].name.c_str(), _eep[i].chi0,
                    _eep[i].J0);
            for (j = 0; ((int)j < _eep[i].nzeta); j++)
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
