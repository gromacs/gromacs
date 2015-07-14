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


//int gtb_comp(const void *a, const void *b);
//int gtd_comp(const void *a, const void *b);

#define assignStr(dst, src)  if (dst) { if (src) {*dst = strdup(src); }else{*dst = NULL; }}
#define assignScal(dst, src) if (dst) *dst = src

#define vectorToArray(vec) (&vec[0])

namespace alexandria
{

  Poldata::Poldata()
    :
    _gtDihedralFunction(egdNR),
   _ngtDihedral(egdNR),
    _ngtDihedralC(egdNR),
_gtDihedralFtype(egdNR),
 

    _gtDihedral(egdNR)
  {
 
    _nptype = 0;
 _nptypeC = 0;
    _nalexandria = 0;
 _nalexandriaC = 0;
    _nbtype = 0;
    _nbrule = 0;
 _nbruleC = 0;
    _nexcl = 0;
    _ngtBond = 0; 
_ngtBondC = 0;
 _gtBondFtype = 0;
    _ngtAngle = 0;
 _ngtAngleC = 0;
 _gtAngleFtype = 0;
      _nmiller = 0;
 _nmillerC = 0;  
    _nbosque = 0;
 _nbosqueC = 0;
    _nsymcharges = 0;
 _nsymchargesC = 0;
    _nep = 0; _nepC = 0;
    _ner = 0; _nerC = 0; 
    _filename = NULL;
    _alexandriaPolarUnit = NULL;
    _alexandriaPolarRef = NULL;
    _alexandriaForcefield = NULL;
    _gtVdwFunction  = NULL;
    _gtCombinationRule = NULL;
    _gtBondFunction = NULL;
    _gtLengthUnit = NULL;
    _gtAngleFunction = NULL;
    _gtAngleUnit = NULL;
    _millerAhpUnit = NULL;
    _bosquePolarUnit = NULL;
    


    for (int x =0; x < egdNR; x++){
      _gtDihedralFunction[x]=NULL;
      _ngtDihedral[x]=0;
      _ngtDihedralC[x]=0;
      // _gtDihedral[x]=NULL;
    }

    /* Initiate some crucial variables */
    _nexcl          = NOTSET;
    _fudgeQQ        = NOTSET;
    _fudgeLJ        = NOTSET;
    _gtCombRule   = NOTSET;
    _gtVdwFtype   = NOTSET;
    _gtBondFtype  = NOTSET;
    _gtAngleFtype = NOTSET;
    for (int i = 0; (i < egdNR); i++)
      {
        _gtDihedralFtype[i] = NOTSET;
      }

  }

  void Poldata::setFilename( char *fn2)
  {
    if (NULL == fn2)
      {
        fprintf(stderr, "Trying to set Poldata filename to NULL\n");
        return;
      }
    if (NULL != _filename)
      {
        fprintf(stderr, "Overwriting Poldata _filename from %s to %s\n",
                _filename, fn2);
        sfree(_filename);
      }
    _filename = strdup(fn2);
  }

  void Poldata::setVdwFunction( const char *func)
  {
    int i;

    if (NULL != _gtVdwFunction)
      {
        sfree(_gtVdwFunction);
      }
    for (i = 0; (i < F_NRE); i++)
      {
        if (strcasecmp(interaction_function[i].name, func) == 0)
	  {
            break;
	  }
      }
    if (i == F_NRE)
      {
        gmx_fatal(FARGS, "Van der Waals function '%s' does not exist in gromacs", func);
      }
    _gtVdwFtype    = i;
    _gtVdwFunction = strdup(func);
  }

  char *Poldata::getVdwFunction()
  {
    return _gtVdwFunction;
  }

  int Poldata::getVdwFtype()
  {
    if (NOTSET == _gtVdwFtype)
      {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  _filename);
      }
    return _gtVdwFtype;
  }

  void Poldata::setNexcl( int nexcl)
  {
    _nexcl = nexcl;
  }

  int Poldata::getNexcl()
  {
    if (NOTSET == _nexcl)
      {
        gmx_fatal(FARGS, "Nexclusions not set in %s", _filename);
      }
    return _nexcl;
  }

  void Poldata::setFudgeQQ( double fudgeQQ)
  {
    _fudgeQQ = fudgeQQ;
  }

  double Poldata::getFudgeQQ()
  {
    if (NOTSET == _fudgeQQ)
      {
        gmx_fatal(FARGS, "_fudgeQQ not set in %s", _filename);
      }
    return _fudgeQQ;
  }

  void Poldata::setFudgeLJ( double fudgeLJ)
  {
    _fudgeLJ = fudgeLJ;
  }

  double Poldata::getFudgeLJ()
  {
    if (NOTSET == _fudgeLJ)
      {
        gmx_fatal(FARGS, "_fudgeLJ not set in %s", _filename);
      }
    return _fudgeLJ;
  }

  void Poldata::setCombinationRule( char *func)
  {
    int i;

    if (NULL != _gtCombinationRule)
      {
        sfree(_gtCombinationRule);
      }
    for (i = 0; (i < eCOMB_NR); i++)
      {
        if (strcasecmp(ecomb_names[i], func) == 0)
	  {
            break;
	  }
      }
    if (i == eCOMB_NR)
      {
        gmx_fatal(FARGS, "Combination rule '%s' does not exist in gromacs", func);
      }
    _gtCombRule        = i;
    _gtCombinationRule = strdup(func);
  }

  char *Poldata::getCombinationRule()
  {
    if (NOTSET == _gtVdwFtype)
      {
        gmx_fatal(FARGS, "No valid function type for Van der Waals function in %s",
                  _filename);
      }
    return _gtCombinationRule;
  }

  int Poldata::getCombRule()
  {
    if (NOTSET == _gtCombRule)
      {
        gmx_fatal(FARGS, "No valid combination rule in %s", _filename);
      }
    return _gtCombRule;
  }

  void Poldata::addPtype(
			 const char   *ptype,
			 const char   *miller,
			 const char   *bosque,
			 double        polarizability,
			 double        sig_pol)
  {
    t_ptype *sp;
    int      i;

    for (i = 0; (i < _nptype); i++)
      {
        if (strcmp(_ptype[i].type, ptype) == 0)
	  {
            break;
	  }
      }
    if (i == _nptype)
      {
        _nptype++;
        _ptype.resize(_nptype);

        sp                 = &(_ptype[i]);
        sp->type           = strdup(ptype);
        sp->bosque         = strdup(bosque);
        sp->miller         = strdup(miller);
        sp->polarizability = polarizability;
        sp->sig_pol        = sig_pol;
      }
    else
      {
        fprintf(stderr, "Polarizability type %s was already added to Poldata record\n", ptype);
      }
  }

  void Poldata::addBtype(
			 const char   *btype)
  {
    int i;

    for (i = 0; (i < _nbtype); i++)
      {
        if (strcmp(btype, _btype[i]) == 0)
	  {
            break;
	  }
      }
    if (i == _nbtype)
      {
	_btype.resize(++_nbtype);
        _btype[i] = strdup(btype);
      }
  }

  void Poldata::addAtype(
			 const char   *elem,
			 const char   *desc,
			 const char   *atype,
			 const char   *ptype,
			 const char   *btype,
			 const char   *vdwparams,
			 double        refEnthalpy)
  {
    t_ffatype *sp;
    int        i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(_alexandria[i].type, atype) == 0)
	  {
            break;
	  }
      }
    if (i == _nalexandria)
      {
        _nalexandria++;
        _alexandria.resize(_nalexandria);

        sp                 = &(_alexandria[i]);
        sp->elem           = strdup(elem);
        sp->desc           = strdup(desc);
        sp->type           = strdup(atype);
        sp->ptype          = strdup(ptype);
        sp->btype          = strdup(btype);
        sp->vdwparams      = strdup(vdwparams);
        sp->ref_enthalpy   = refEnthalpy;
        addBtype(btype);
      }
    else
      {
        fprintf(stderr, "Atom type %s was already added to Poldata record\n", atype);
      }
  }

  void Poldata::addBondingRule(
			       char *gtBrule, char *atype,
			       char *geometry, int numbonds,
			       double valence, int iAromatic,
			       char *neighbors)
  {
    t_brule *sp;
    int      i, j;

    for (j = 0; (j < _nalexandria); j++)
      {
        if (strcmp(_alexandria[j].type, atype) == 0)
	  {
            break;
	  }
      }
    if (j < _nalexandria)
      {
        for (i = 0; (i < _nbrule); i++)
	  {
            if (strcmp(_brule[i].rule, gtBrule) == 0)
	      {
                break;
	      }
	  }
        if (i == _nbrule)
	  {
            _nbrule++;
            _brule.resize(_nbrule);

            sp                 = &(_brule[i]);
            sp->elem           = strdup(_alexandria[j].elem);
            sp->rule           = strdup(gtBrule);
            sp->type           = strdup(atype);
            sp->neighbors      = strdup(neighbors);
            sp->valence        = valence;
            sp->iAromatic      = iAromatic;
            sp->nb             = split(neighbors, ' ');
            sp->geometry       = strdup(geometry);
            sp->numbonds       = numbonds;
	  }
        else
	  {
            fprintf(stderr, "Bonding rule %s was already added to Poldata record\n", gtBrule);
	  }
      }
    else if (NULL != debug)
      {
        fprintf(debug, "Ignoring bonding rule involving unknown atom type %s\n",
                atype);
      }
  }

  int Poldata::getBondingRule(
			      char **gt_brule, char **atype,
			      char **geometry, int *numbonds,
			      double *valence, int *iAromatic,
			      char **neighbors)
  {
    if (_nbruleC < _nbrule)
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

  int Poldata::getNatypes()
  {
    return _nalexandria;
  }

  int Poldata::getNptypes()
  {
    return _nptype;
  }

  int Poldata::getNgtBond()
  {
    return _ngtBond;
  }

  int Poldata::getNgtAngle()
  {
    return _ngtAngle;
  }

  int Poldata::getNgtDihedral( int egd)
  {
    return _ngtDihedral[egd];
  }

  void Poldata::setBondFunction( char *fn)
  {
    int i;

    if (NULL != _gtBondFunction)
      {
        sfree(_gtBondFunction);
      }
    for (i = 0; (i < F_NRE); i++)
      {
        if (strcasecmp(interaction_function[i].name, fn) == 0)
	  {
            break;
	  }
      }
    if (i == F_NRE)
      {
        gmx_fatal(FARGS, "Bond function '%s' does not exist in gromacs", fn);
      }
    _gtBondFtype    = i;
    _gtBondFunction = strdup(fn);
  }

  char *Poldata::getBondFunction()
  {
    return _gtBondFunction;
  }

  void Poldata::setAngleFunction( char *fn)
  {
    int i;

    if (NULL != _gtAngleFunction)
      {
        sfree(_gtAngleFunction);
      }
    for (i = 0; (i < F_NRE); i++)
      {
        if (strcasecmp(interaction_function[i].name, fn) == 0)
	  {
            break;
	  }
      }
    if (i == F_NRE)
      {
        gmx_fatal(FARGS, "Angle function '%s' does not exist in gromacs", fn);
      }
    _gtAngleFtype    = i;
    _gtAngleFunction = strdup(fn);
  }

  char *Poldata::getAngleFunction()
  {
    return _gtAngleFunction;
  }

  int Poldata::getBondFtype()
  {
    return _gtBondFtype;
  }

  int Poldata::getAngleFtype()
  {
    return _gtAngleFtype;
  }

  int Poldata::getDihedralFtype( int egd)
  {
    return _gtDihedralFtype[egd];
  }

  void Poldata::setDihedralFunction( int egd, char *func)
  {
    int i;

    if (NULL != _gtDihedralFunction[egd])
      {
        sfree(_gtDihedralFunction[egd]);
      }
    for (i = 0; (i < F_NRE); i++)
      {
        if (strcasecmp(interaction_function[i].name, func) == 0)
	  {
            break;
	  }
      }
    if (i == F_NRE)
      {
        gmx_fatal(FARGS, "Dihedral function '%s' does not exist in gromacs", func);
      }
    _gtDihedralFtype[egd]    = i;
    _gtDihedralFunction[egd] = strdup(func);
  }

  char *Poldata::getDihedralFunction( int egd)
  {
    return _gtDihedralFunction[egd];
  }

  void Poldata::setPtypePolarizability( const char *ptype,
					double polarizability, double sig_pol)
  {
    t_ptype *sp;
    int      i;

    for (i = 0; (i < _nptype); i++)
      {
        if (strcmp(ptype, _ptype[i].type) == 0)
	  {
            sp                 = &(_ptype[i]);
            sp->polarizability = polarizability;
            sp->sig_pol        = sig_pol;
            break;
	  }
      }
    if (i == _nptype)
      {
        fprintf(stderr, "No such ptype %s when trying to set the polarizability.\n", ptype);
      }
  }

  char *Poldata::getPolarUnit()
  {
    return _alexandriaPolarUnit;
  }

  char *Poldata::getPolarRef()
  {
    return _alexandriaPolarRef;
  }

  char *Poldata::getLengthUnit()
  {
    if (NULL == _gtLengthUnit)
      {
        gmx_fatal(FARGS, "No length unit in %s",
                  (NULL != _filename) ? _filename : "unknown");
      }

    return _gtLengthUnit;
  }

  void Poldata::setPolarUnit( const char *polarUnit)
  {
    _alexandriaPolarUnit   = strdup(polarUnit);
  }

  void Poldata::setPolarRef( const char *polarRef)
  {
    _alexandriaPolarRef   = strdup(polarRef);
  }

  char *Poldata::getForceField()
  {
    return _alexandriaForcefield;
  }

  void Poldata::setForceField( const char *forcefield)
  {
    _alexandriaForcefield   = strdup(forcefield);
  }

  void Poldata::setLengthUnit( const char *length_unit)
  {
    _gtLengthUnit   = strdup(length_unit);
  }

  char *Poldata::getGeometry( char *gtBrule)
  {
    int i;

    if (gtBrule)
      {
        for (i = 0; (i < _nbrule); i++)
	  {
            if (strcmp(_brule[i].rule, gtBrule) == 0)
	      {
                return _brule[i].geometry;
	      }
	  }
      }

    return NULL;
  }

  char *Poldata::getDesc( char *atype)
  {
    int i;

    if (atype)
      {
        for (i = 0; (i < _nalexandria); i++)
	  {
            if (strcmp(_alexandria[i].type, atype) == 0)
	      {
                return _alexandria[i].desc;
	      }
	  }
      }

    return NULL;
  }

  gmx_bool Poldata::strcasestrStart(char *needle, char *haystack)
  {
    char *ptr;

    ptr = strcasestr(haystack, needle);

    return (ptr == haystack);
  }

  int Poldata::countNeighbors(t_brule *brule, int nbond, char *nbhybrid[], int *score)
  {
    int j, ni = 0, *jj, IFound;

    *score = 0;
    snew(jj, nbond+1);
    for (unsigned int i = 0; (i < brule->nb.size()); i++)
      {
        IFound = 0;
        for (j = 0; (j < nbond); j++)
	  {
            if (
                (NULL != nbhybrid[j]) &&
                (jj[j] == 0) &&
                (IFound == 0) &&
                strcasestrStart((char *)brule->nb[i].c_str(), nbhybrid[j])
                )
	      {
                IFound = 1;
                jj[j]   = j+1;
                ni++;
                (*score) += 1;
                if (strlen(nbhybrid[j]) > 0)
		  {
                    (*score) += 1;
		  }
	      }
	  }
      }
    sfree(jj);

    return ni;
  }

  char **Poldata::getBondingRules( char *elem,
				   int nbond, char *neighbors[],
				   const char *geometry,
				   int iAromatic)
  {
    unsigned int    nnb;
    int             i, nptr = 0, best = -1, score;
    char          **ptr = NULL;

    for (i = 0; (i < _nbrule); i++)
      {
        nnb = countNeighbors(&(_brule[i]), nbond, neighbors, &score);
        if ((strcmp(_brule[i].elem, elem) == 0) &&
            (strcasecmp(_brule[i].geometry, geometry) == 0) &&
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
    ptr[nptr-1] = NULL;

    return ptr;
  }

  int Poldata::bondingRuleValence( char *gtBrule, double *valence)
  {
    int i;

    for (i = 0; (i < _nbrule); i++)
      {
        if (strcasecmp(gtBrule, _brule[i].rule) == 0)
	  {
            *valence = _brule[i].valence;
            return 1;
	  }
      }
    return 0;
  }

  int Poldata::getPtypePol( const char *ptype,
			    double *polar, double *sig_pol)
  {
    int j;

    for (j = 0; (j < _nptype); j++)
      {
        if (strcmp(ptype, _ptype[j].type) == 0)
	  {
            if (NULL != polar)
	      {
                *polar   = _ptype[j].polarizability;
	      }

            if (NULL != sig_pol)
	      {
                *sig_pol = _ptype[j].sig_pol;
	      }
            return 1;
	  }
      }
    return 0;
  }

  int Poldata::getAtypePol( const char *atype,
			    double *polar, double *sigPol)
  {
    int i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(atype, _alexandria[i].type) == 0)
	  {
            return getPtypePol( _alexandria[i].ptype, polar, sigPol);
	  }
      }
    return 0;
  }


  int Poldata::getAtypeRefEnthalpy( const char *atype,
				    double *Href)
  {
    int i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(atype, _alexandria[i].type) == 0)
	  {
            *Href = _alexandria[i].ref_enthalpy;
            return 1;
	  }
      }
    return 0;
  }

  char *Poldata::ptypeToMiller( const char *ptype)
  {
    int i;

    for (i = 0; (i < _nptype); i++)
      {
        if (strcmp(ptype, _ptype[i].type) == 0)
	  {
            return _ptype[i].miller;
	  }
      }
    return NULL;
  }

  char *Poldata::ptypeToBosque( const char *ptype)
  {
    int i;

    for (i = 0; (i < _nptype); i++)
      {
        if (strcmp(ptype, _ptype[i].type) == 0)
	  {
            return _ptype[i].bosque;
	  }
      }
    return NULL;
  }

  int Poldata::getPtype(
			char        **ptype,
			char        **miller,
			char        **bosque,
			double       *polarizability,
			double       *sigPol)
  {
    t_ptype *sp;

    if (_nptypeC < _nptype)
      {
        sp = &(_ptype[_nptypeC]);
        assignScal(polarizability, sp->polarizability);
        assignScal(sigPol, sp->sig_pol);
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
			char        **elem,
			char        **desc,
			char        **atype,
			char        **ptype,
			char        **btype,
			char        **vdwparams,
			double       *refEnthalpy)
  {
    t_ffatype *sp;

    if (_nalexandriaC < _nalexandria)
      {
        sp = &(_alexandria[_nalexandriaC]);
        assignStr(elem, sp->elem);
        assignStr(desc, sp->desc);
        assignStr(atype, sp->type);
        assignStr(ptype, sp->ptype);
        assignStr(btype, sp->btype);
        assignStr(vdwparams, sp->vdwparams);
        *refEnthalpy = sp->ref_enthalpy;
        _nalexandriaC++;
        return 1;
      }
    else
      {
        _nalexandriaC = 0;
      }

    return 0;
  }

  const char *Poldata::atypeToPtype( const char *atype)
  {
    int i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(_alexandria[i].type, atype) == 0)
	  {
            return _alexandria[i].ptype;
	  }
      }
    return NULL;
  }

  const char *Poldata::atypeToBtype( const char *atype)
  {
    int i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(_alexandria[i].type, atype) == 0)
	  {
            return _alexandria[i].btype;
	  }
      }
    return NULL;
  }

  int Poldata::searchAtype(
			   char         *key,
			   char        **elem,
			   char        **desc,
			   char        **atype,
			   char        **ptype,
			   char        **btype,
			   char        **vdwparams)
  {
    t_ffatype *sp;
    int        i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(key, _alexandria[i].type) == 0)
	  {
            break;
	  }
      }

    if (i < _nalexandria)
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

  double Poldata::elemGetMaxValence( char *elem)
  {
    double mv = 0;
    int    i;

    for (i = 0; (i < _nbrule); i++)
      {
        if ((0 == gmx_strcasecmp(_brule[i].elem, elem)) &&
            (mv < _brule[i].valence))
	  {
            mv = _brule[i].valence;
	  }
      }
    return mv;
  }

  double *Poldata::elemGetBondorders( char *elem1, char *elem2,
				      double distance, double toler)
  {
    double  dev, *bo = NULL;
    char   *ba1, *ba2;
    int     i, j, k, nbo;

    if ((NULL == elem1) || (NULL == elem2))
      {
        return 0;
      }
    nbo = 0;
    for (i = 0; (i < _ngtBond); i++)
      {
        if (0 == strlen(_gtBond[i].elem1))
	  {
            for (j = 0; (j < _nalexandria); j++)
	      {
                if (strcmp(_alexandria[j].type, _gtBond[i].atom2) == 0)
		  {
                    strcpy(_gtBond[i].elem1, _alexandria[j].elem);
		  }
	      }
	  }
        if (0 == strlen(_gtBond[i].elem2))
	  {
            for (j = 0; (j < _nalexandria); j++)
	      {
                if (strcmp(_alexandria[j].type, _gtBond[i].atom2) == 0)
		  {
                    strcpy(_gtBond[i].elem2, _alexandria[j].elem);
		  }
	      }
	  }
        ba1 = _gtBond[i].elem1;
        ba2 = _gtBond[i].elem2;
        if (((strcmp(ba1, elem1) == 0) && (strcmp(ba2, elem2) == 0)) ||
            ((strcmp(ba1, elem2) == 0) && (strcmp(ba2, elem1) == 0)))
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

  int Poldata::elemIsBond( char *elem1, char *elem2,
			   double distance, double toler)
  {
    char  *ba1, *ba2;
    double dev, devBest = 100000;
    int    j, i;

    if ((NULL == elem1) || (NULL == elem2))
      {
        return 0;
      }
    for (i = 0; (i < _ngtBond); i++)
      {
        if (0 == strlen(_gtBond[i].elem1))
	  {
            for (j = 0; (j < _nalexandria); j++)
	      {
                if (strcmp(_alexandria[j].type, _gtBond[i].atom2) == 0)
		  {
                    strcpy(_gtBond[i].elem1, _alexandria[j].elem);
		  }
	      }
	  }
        if (0 == strlen(_gtBond[i].elem2))
	  {
            for (j = 0; (j < _nalexandria); j++)
	      {
                if (strcmp(_alexandria[j].type, _gtBond[i].atom2) == 0)
		  {
                    strcpy(_gtBond[i].elem2, _alexandria[j].elem);
		  }
	      }
	  }
        ba1 = _gtBond[i].elem1;
        ba2 = _gtBond[i].elem2;
        if (((strcmp(ba1, elem1) == 0) && (strcmp(ba2, elem2) == 0)) ||
            ((strcmp(ba1, elem2) == 0) && (strcmp(ba2, elem1) == 0)))
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

  int loGtbComp(t_gt_bond *ba, t_gt_bond *bb)
  {
    char *a1, *a2, *b1, *b2;
    int   i;

    if (strcmp(ba->atom1, ba->atom2) <= 0)
      {
        a1 = ba->atom1;
        a2 = ba->atom2;
      }
    else
      {
        a2 = ba->atom1;
        a1 = ba->atom2;
      }
    if (strcmp(bb->atom1, bb->atom2) <= 0)
      {
        b1 = bb->atom1;
        b2 = bb->atom2;
      }
    else
      {
        b2 = bb->atom1;
        b1 = bb->atom2;
      }
    i = strcmp(a1, b1);
    if (0 == i)
      {
        i = strcmp(a2, b2);
      }

    return i;
  }


  int Poldata::gtbComp(const void *a, const void *b)
  {
    t_gt_bond *ba = (t_gt_bond *)a;
    t_gt_bond *bb = (t_gt_bond *)b;
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

  t_gt_bond *Poldata::searchBond( char *atom1, char *atom2,
				  double bondorder)
  {
    t_gt_bond key, *gtB;
    int       i;

    key.atom1     = atom1;
    key.atom2     = atom2;
    key.bondorder = bondorder;

    gtB = (t_gt_bond *) bsearch(&key, vectorToArray(_gtBond), _ngtBond, sizeof(key), gtbComp);
    if (NULL != gtB)
      {
        i = indexOfPointInVector(gtB,_gtBond);
        while ((i > 0) && (loGtbComp(&(_gtBond[i-1]), &(_gtBond[i])) == 0))
	  {
            i--;
	  }
        gtB = &(_gtBond[i]);
      }
    return gtB;
  }

  double Poldata::atypeBondorder( char *atype1, char *atype2,
				  double distance, double toler)
  {
    double     dev, devBest = 100000;
    int        i, iBest = -1;
    t_gt_bond *gtB;

    if ((NULL == atype1) || (NULL == atype2))
      {
        return 0.0;
      }
    gtB = searchBond( atype1, atype2, 0);
    if (NULL != gtB)
      {
        i = indexOfPointInVector(gtB,_gtBond);
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
        while ((i < _ngtBond) &&
               (0 == loGtbComp(&(_gtBond[i]), &(_gtBond[i-1]))));
      }
    if (devBest < toler)
      {
        return _gtBond[iBest].bondorder;
      }

    return 0.0;
  }

  void Poldata::addMiller(
			  char         *miller,
			  int           atomnumber,
			  double        tauAhc,
			  double        alphaAhp)
  {
    t_miller *mil;

    _nmiller++;
    _miller.resize(_nmiller);
    mil             = &(_miller[_nmiller-1]);
    mil->miller     = strdup(miller);
    mil->atomnumber = atomnumber;
    mil->tau_ahc    = tauAhc;
    mil->alpha_ahp  = alphaAhp;
  }

  void Poldata::setMillerUnits( char *tauUnit, char *ahpUnit)
  {
    _millerTauUnit = strdup(tauUnit);
    _millerAhpUnit = strdup(ahpUnit);
  }

  void Poldata::getMillerUnits( char **tauUnit,
				char **ahpUnit)
  {
    assignStr(tauUnit, _millerTauUnit);
    assignStr(ahpUnit, _millerAhpUnit);
  }

  int Poldata::getMiller(
			 char        **miller,
			 int          *atomnumber,
			 double       *tauAhc,
			 double       *alphaAhp)
  {
    t_miller *mil;
    int       i;

    i = _nmillerC;

    if (i < _nmiller)
      {
        mil = &(_miller[i]);
        assignStr(miller, mil->miller);
        assignScal(atomnumber, mil->atomnumber);
        assignScal(tauAhc, mil->tau_ahc);
        assignScal(alphaAhp, mil->alpha_ahp);
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
			    char         *miller,
			    int          *atomnumber,
			    double       *tauAhc,
			    double       *alphaAhp)
  {
    t_miller *mil;
    int       i;

    for (i = 0; (i < _nmiller); i++)
      {
        if (strcmp(miller, _miller[i].miller) == 0)
	  {
            mil = &(_miller[i]);
            assignScal(atomnumber, mil->atomnumber);
            assignScal(tauAhc, mil->tau_ahc);
            assignScal(alphaAhp, mil->alpha_ahp);

            return 1;
	  }
      }

    return 0;
  }

  void Poldata::addBosque(
			  char         *bosque,
			  double        polarizability)
  {
    t_bosque *bs;

    _nbosque++;
    _bosque.resize(_nbosque);
    bs                 = &(_bosque[_nbosque-1]);
    bs->bosque         = strdup(bosque);
    bs->polarizability = polarizability;
  }

  int Poldata::getBosque(
			 char        **bosque,
			 double       *polarizability)
  {
    if (_nbosqueC < _nbosque)
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
			    char         *bosque,
			    double       *polarizability)
  {
    int i;

    for (i = 0; (i < _nbosque); i++)
      {
        if (strcasecmp(bosque, _bosque[i].bosque) == 0)
	  {
            *polarizability = _bosque[i].polarizability;
            return 1;
	  }
      }
    return 0;
  }

  void Poldata::setBosqueUnit( char *polarUnit)
  {
    _bosquePolarUnit   = strdup(polarUnit);
  }

  char *Poldata::getBosqueUnit()
  {
    return _bosquePolarUnit;
  }


  int Poldata::searchBondtype( char *atom)
  {
    int j;

    for (j = 0; (j < _nbtype); j++)
      {
        if (strcmp(_btype[j], atom) == 0)
	  {
            return j;
	  }
      }
    return -1;
  }

  /*
   * _gtBond stuff
   */
  int Poldata::setBondParams( char *atom1, char *atom2,
			      double length, double sigma, int ntrain,
			      double bondorder, char *params)
  {
    t_gt_bond *gtB;
    int        i;

    for (i = 0; (i < _ngtBond); i++)
      {
        gtB = &(_gtBond[i]);
        if (((((strcmp(gtB->atom1, atom1) == 0) &&
               (strcmp(gtB->atom2, atom2) == 0)) ||
              ((strcmp(gtB->atom1, atom2) == 0) &&
               (strcmp(gtB->atom2, atom1) == 0)))) &&
            ((bondorder == 0) || (gtB->bondorder == bondorder)))
	  {
            break;
	  }
      }
    if (i < _ngtBond)
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
        if (NULL != gtB->params)
	  {
            sfree(gtB->params);
	  }
        gtB->params = strdup(params);
        return 1;
      }
    return 0;
  }

  int Poldata::addBond( char *atom1, char *atom2,
			double length, double sigma, int ntrain,
			double bondorder, char *params)
  {
    t_gt_bond *gtB;
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
        _ngtBond++;
        _gtBond.resize(_ngtBond);
        gtB            = &(_gtBond[_ngtBond-1]);
        gtB->atom1     = strdup(atom1);
        strncpy(gtB->elem1, _alexandria[a1].elem, sizeof(gtB->elem1)-1);
        gtB->atom2     = strdup(atom2);
        strncpy(gtB->elem2, _alexandria[a2].elem, sizeof(gtB->elem2)-1);
        gtB->length    = length;
        gtB->sigma     = sigma;
        gtB->ntrain    = ntrain;
        gtB->bondorder = bondorder;
        gtB->params    = strdup(params);
        qsort(vectorToArray(_gtBond), _ngtBond, sizeof(_gtBond[0]), gtbComp);
      }
    return 1;
  }

  int Poldata::getBond( char **atom1, char **atom2,
			double *length, double *sigma, int *ntrain,
			double *bondorder, char **params)
  {
    t_gt_bond *gtB;

    if (_ngtBondC < _ngtBond)
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

  int Poldata::searchBond( char *atom1, char *atom2,
			   double *length, double *sigma, int *ntrain,
			   double *bondorder, char **params)
  {
    t_gt_bond *gtB;

    if ((NULL == atom1) || (NULL == atom2))
      {
        return 0;
      }
    gtB = searchBond( atom1, atom2, 0);
    if (NULL != gtB)
      {
        if (((strcmp(gtB->atom1, atom1) == 0) &&
             (strcmp(gtB->atom2, atom2) == 0)) ||
            ((strcmp(gtB->atom1, atom2) == 0) &&
             (strcmp(gtB->atom2, atom1) == 0)))
	  {
            assignScal(length, gtB->length);
            assignScal(sigma, gtB->sigma);
            assignScal(ntrain, gtB->ntrain);
            assignScal(bondorder, gtB->bondorder);
            assignStr(params, gtB->params);

            return 1+indexOfPointInVector(gtB,_gtBond);
	  }
      }
    return 0;
  }

  /*
   * gt_angle stuff
   */
  int Poldata::setAngleParams( char *atom1, char *atom2,
			       char *atom3, double angle, double sigma, int ntrain,
			       char *params)
  {
    t_gt_angle *gtB;
    int         i;

    for (i = 0; (i < _ngtAngle); i++)
      {
        gtB = &(_gtAngle[i]);
        if ((strcmp(gtB->atom2, atom2) == 0) &&
            (((strcmp(gtB->atom1, atom1) == 0) &&
              (strcmp(gtB->atom3, atom3) == 0)) ||
             ((strcmp(gtB->atom1, atom3) == 0) &&
              (strcmp(gtB->atom3, atom1) == 0))))
	  {
            break;
	  }
      }
    if (i < _ngtAngle)
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
        if (NULL != gtB->params)
	  {
            sfree(gtB->params);
	  }
        gtB->params = strdup(params);
        return 1;
      }
    return 0;
  }

  int Poldata::addAngle(
			char *atom1, char *atom2,
			char *atom3, double angle, double sigma,
			int ntrain, char *params)
  {
    t_gt_angle *gtB;

    if ((-1 == searchBondtype( atom1)) ||
        (-1 == searchBondtype( atom2)) ||
        (-1 == searchBondtype( atom3)))
      {
        return 0;
      }

    if (0 == setAngleParams( atom1, atom2, atom3, angle, sigma, ntrain, params))
      {
        _ngtAngle++;
        _gtAngle.resize(_ngtAngle);
        gtB          = &(_gtAngle[_ngtAngle-1]);
        gtB->atom1   = strdup(atom1);
        gtB->atom2   = strdup(atom2);
        gtB->atom3   = strdup(atom3);
        gtB->angle   = angle;
        gtB->sigma   = sigma;
        gtB->ntrain  = ntrain;
        gtB->params  = strdup(params);
      }
    return 1;
  }

  int Poldata::getAngle( char **atom1, char **atom2,
			 char **atom3, double *angle, double *sigma,
			 int *ntrain, char **params)
  {
    t_gt_angle *gtB;

    if (_ngtAngleC < _ngtAngleC)
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

  int Poldata::searchAngle( char *atom1, char *atom2,
			    char *atom3, double *angle, double *sigma,
			    int *ntrain, char **params)
  {
    t_gt_angle *gtB;
    int         i;

    if ((NULL == atom1) || (NULL == atom2) || (NULL == atom3))
      {
        return 0;
      }
    for (i = 0; (i < _ngtAngleC); i++)
      {
        gtB = &(_gtAngle[i]);
        if ((strcmp(gtB->atom2, atom2) == 0) &&
            (((strcmp(gtB->atom1, atom1) == 0) &&
              (strcmp(gtB->atom3, atom3) == 0)) ||
             ((strcmp(gtB->atom1, atom3) == 0) &&
              (strcmp(gtB->atom3, atom1) == 0))))
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

  void Poldata::setAngleUnit( char *angleUnit)
  {
    _gtAngleUnit   = strdup(angleUnit);
  }

  char *Poldata::getAngleUnit()
  {
    return _gtAngleUnit;
  }

  /*
   * gt_dihedral stuff
   */
  int Poldata::gtdComp(const void *a, const void *b)
  {
    t_gt_dihedral *gtA = (t_gt_dihedral *)a;
    t_gt_dihedral *gtB = (t_gt_dihedral *)b;
    int            n;

    if (0 == (n = strcmp(gtA->atom1, gtB->atom1)))
      {
        if (0 == (n = strcmp(gtA->atom2, gtB->atom2)))
	  {
            if (0 == (n = strcmp(gtA->atom3, gtB->atom3)))
	      {
                n = strcmp(gtA->atom4, gtB->atom4);
	      }
	  }
      }

    return n;
  }

  t_gt_dihedral *Poldata::searchDihedral( int egd,
					  char *atom1, char *atom2,
					  char *atom3, char *atom4)
  {
    t_gt_dihedral gtA, *gtRes, *gtDptr;
    int           nd;

    if ((NULL == atom1) || (NULL == atom2) || (NULL == atom3) || (NULL == atom4))
      {
        return NULL;
      }
    gtDptr    = vectorToArray(_gtDihedral[egd]);
    nd         = _ngtDihedral[egd];
    gtA.atom1 = atom1;
    gtA.atom2 = atom2;
    gtA.atom3 = atom3;
    gtA.atom4 = atom4;
    gtRes     = (t_gt_dihedral *) bsearch(&gtA, gtDptr, nd, sizeof(gtA), &gtdComp);
    if (NULL == gtRes)
      {
        gtA.atom1 = atom4;
        gtA.atom2 = atom3;
        gtA.atom3 = atom2;
        gtA.atom4 = atom1;
        gtRes     = (t_gt_dihedral *) bsearch(&gtA, gtDptr, nd, sizeof(gtA), gtdComp);
      }
    return gtRes;
  }

  int Poldata::setDihedralParams( int egd,
				  char *atom1, char *atom2,
				  char *atom3, char *atom4,
				  double dihedral, double sigma, int ntrain,
				  char *params)
  {
    t_gt_dihedral *gtB;


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
        if (NULL != gtB->params)
	  {
            sfree(gtB->params);
	  }
        gtB->params = strdup(params);
        return 1;
      }
    return 0;
  }

  int Poldata::addDihedral( int egd,
			    char *atom1, char *atom2,
			    char *atom3, char *atom4, double dihedral,
			    double sigma, int ntrain, char *params)
  {
    t_gt_dihedral *gtB;

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
        _ngtDihedral[egd]++;
        _gtDihedral[egd].resize(_ngtDihedral[egd]);
        gtB           = &(_gtDihedral[egd][_ngtDihedral[egd]-1]);
        gtB->atom1    = strdup(atom1);
        gtB->atom2    = strdup(atom2);
        gtB->atom3    = strdup(atom3);
        gtB->atom4    = strdup(atom4);
        gtB->dihedral = dihedral;
        gtB->sigma    = sigma;
        gtB->ntrain   = ntrain;
        gtB->params   = strdup(params);
        qsort(vectorToArray(_gtDihedral[egd]), _ngtDihedral[egd], sizeof(_gtDihedral[egd][0]),
              gtdComp);
      }
    return 1;
  }

  int Poldata::getDihedral( int egd,
			    char **atom1, char **atom2,
			    char **atom3, char **atom4, double *dihedral,
			    double *sigma, int *ntrain, char **params)
  {
    t_gt_dihedral *gtB;

    if (_ngtDihedralC[egd] < _ngtDihedral[egd])
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
			       char *atom1, char *atom2,
			       char *atom3, char *atom4,
			       double *dihedral, double *sigma,
			       int *ntrain, char **params)
  {
    t_gt_dihedral *gtRes;

    gtRes = searchDihedral( egd, atom1, atom2, atom3, atom4);
    if (NULL != gtRes)
      {
        assignScal(dihedral, gtRes->dihedral);
        assignScal(sigma, gtRes->sigma);
        assignScal(ntrain, gtRes->ntrain);
        assignStr(params, gtRes->params);

	return 1 + (indexOfPointInVector(gtRes,_gtDihedral[egd]));
      }
    return 0;
  }

  void Poldata::addSymcharges( char *central,
			       char *attached, int numattach)
  {
    t_symcharges *sc;
    int           i;

    for (i = 0; (i < _nsymcharges); i++)
      {
        sc = &(_symcharges[i]);
        if ((strcasecmp(sc->central, central) == 0) &&
            (strcasecmp(sc->attached, attached) == 0) &&
            (sc->numattach == numattach))
	  {
            break;
	  }
      }
    if (i == _nsymcharges)
      {
        _nsymcharges++;
        _symcharges.resize(_nsymcharges);
        sc              = &(_symcharges[i]);
        sc->central     = strdup(central);
        sc->attached    = strdup(attached);
        sc->numattach   = numattach;
      }
  }

  int Poldata::getSymcharges( char **central,
			      char **attached, int *numattach)
  {
    t_symcharges *sc;

    if (_nsymchargesC < _nsymcharges)
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

  int Poldata::searchSymcharges( char *central,
				 char *attached, int numattach)
  {
    t_symcharges *sc;
    int           i;

    for (i = 0; (i < _nsymcharges); i++)
      {
        sc = &(_symcharges[i]);
        if ((strcasecmp(sc->central, central) == 0) &&
            (strcasecmp(sc->attached, attached) == 0) &&
            (sc->numattach == numattach))
	  {
            return 1;
	  }
      }

    return 0;
  }

  /* Electrostatics properties */
  t_eemprops *Poldata::getEep(ChargeDistributionModel eqdModel,
			      const char *name)
  {
    int i;

    for (i = 0; (i < _nep); i++)
      {
        if ((strcasecmp(_eep[i].name, name) == 0) &&
            (_eep[i].eqd_model == eqdModel))
	  {
            return &(_eep[i]);
	  }
      }
    return NULL;
  }

  void Poldata::setEemprops(ChargeDistributionModel eqdModel, char *name,
			    double J0, double chi0, char *zeta, char *q, char *row)
  {
   
    t_eemprops              *eep;
    std::vector<std::string> sz, sq, sr;

    eep = getEep(eqdModel, name);
    if (NULL == eep)
      {
	_nep++;
	_eep.resize(_nep );
        eep = &(_eep[_nep-1]);
      }
    eep->eqd_model = eqdModel;
    strncpy(eep->name, name, EEMBUFSIZE-1);
    eep->name[EEMBUFSIZE-1] = '\0';
    eep->J0                 = J0;
    sz = split(zeta, ' ');
    sq = split(q, ' ');
    sr = split(row, ' ');
    strncpy(eep->zetastr, zeta, EEMBUFSIZE-1);
    strncpy(eep->qstr, q, EEMBUFSIZE-1);
    strncpy(eep->rowstr, row, EEMBUFSIZE-1);
    unsigned int nn = std::min(sz.size(), std::min(sq.size(), sr.size()));
    unsigned int n;
    for (n = 0; (n < nn); n++)
      {
        if (n < MAXZETA)
	  {
            eep->zeta[n] = atof(sz[n].c_str());
            eep->q[n]    = atof(sq[n].c_str());
            eep->row[n]  = atoi(sr[n].c_str());
	  }
      }
    if (sz.size() > nn)
      {
        fprintf(stderr, "Warning: more zeta values than q/row values for %s n = %d\n",
                name, nn);
      }
    if (sq.size() > nn)
      {
        fprintf(stderr, "Warning: more q values than zeta/row values for %s n = %d\n",
                name, nn);
      }
    if (sr.size() > nn)
      {
        fprintf(stderr, "Warning: more row values than q/zeta values for %s n = %d\n",
                name, nn);
      }
    eep->nzeta = nn;
    if (nn >= MAXZETA)
      {
        fprintf(stderr, "More than %d zeta and/or q values for %s\n", MAXZETA, eep->name);
        eep->nzeta = MAXZETA;
      }
    for (; (n < MAXZETA); n++)
      {
        eep->zeta[n] = 0;
        eep->q[n]    = 0;
        eep->row[n]  = 0;
      }
    eep->chi0  = chi0;
  }

  int Poldata::getEemprops(
			   ChargeDistributionModel *eqdModel, char **name,
			   double *J0, double *chi0, char **zeta, char **q, char **row)
  {
    if (_nepC < _nep)
      {
        assignScal(eqdModel, _eep[_nepC].eqd_model);
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
    int i, n = 0;

    for (i = 0; (i < _nep); i++)
      {
        if (_eep[i].eqd_model == eqdModel)
	  {
            n++;
	  }
      }

    return n;
  }

  int Poldata::havePolSupport( const char *atype)
  {
    int i;

    for (i = 0; (i < _nalexandria); i++)
      {
        if (strcmp(atype, _alexandria[i].type) == 0)
	  {
            return 1;
	  }
      }
    return 0;
  }

  int Poldata::haveEemSupport( ChargeDistributionModel eqdModel,
			       const char *name,
			       gmx_bool bAllowZeroParameters)
  {
   
    t_eemprops  *eep  = getEep(eqdModel, name);

    return (eep && (bAllowZeroParameters || ((eep->J0 > 0) && (eep->chi0 > 0))));
  }

  double Poldata::getJ00( ChargeDistributionModel eqdModel, char *name)
  {
    t_eemprops  *eer;

    if ((eer = getEep(eqdModel, name)) != NULL)
      {
        return eer->J0;
      }
    else
      {
        gmx_fatal(FARGS, "No J0 data for eqdModel %d and name %s",
                  eqdModel, name);
      }
    return -1;
  }

  char *Poldata::getQstr( ChargeDistributionModel eqdModel, char *name)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        return eer->qstr;
      }
    return NULL;
  }

  char *Poldata::getRowstr( ChargeDistributionModel eqdModel, char *name)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        return eer->rowstr;
      }
    return NULL;
  }

  int Poldata::getRow( ChargeDistributionModel eqdModel, char *name, int zz)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        range_check(zz, 0, eer->nzeta);
        return eer->row[zz];
      }
    return -1;
  }

  double Poldata::getZeta( ChargeDistributionModel eqdModel, char *name, int zz)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        if ((zz < 0) || (zz >= eer->nzeta))
	  {
            printf("Bleh\n");
	  }
        range_check(zz, 0, eer->nzeta);
        return eer->zeta[zz];
      }
    return -1;
  }

  int Poldata::getNzeta( ChargeDistributionModel eqdModel, char *name)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        return eer->nzeta;
      }
    return 0;
  }

  double Poldata::getQ( ChargeDistributionModel eqdModel, char *name, int zz)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        range_check(zz, 0, eer->nzeta);
        return eer->q[zz];
      }
    return -1;
  }

  double Poldata::getChi0( ChargeDistributionModel eqdModel, char *name)
  {
    t_eemprops *eer;

    if ((eer = getEep( eqdModel, name)) != NULL)
      {
        return eer->chi0;
      }
    else
      {
        gmx_fatal(FARGS, "No chi0 data for eqdModel %d and name %s", eqdModel, name);
      }
    return -1;
  }

  void Poldata::setEpref( ChargeDistributionModel eqdModel, char *epref)
  {
    int i;

    for (i = 0; (i < _ner); i++)
      {
        if (_epr[i].eqd_model == eqdModel)
	  {
            if (_epr[i].epref)
	      {
                sfree(_epr[i].epref);
	      }
            _epr[i].epref = strdup(epref);
            break;
	  }
      }
    if (i == _ner)
      {
	_epr.resize(++_ner);
        _epr[i].eqd_model = eqdModel;
        _epr[i].epref     = strdup(epref);
      }
  }

  char *Poldata::getEpref( ChargeDistributionModel eqdModel)
  {
    int i;

    for (i = 0; (i < _ner); i++)
      {
        if (_epr[i].eqd_model == eqdModel)
	  {
            return _epr[i].epref;
	  }
      }
    return NULL;
  }

  int Poldata::listEpref( ChargeDistributionModel *eqdModel, char **epref)
  {
    if (_nerC < _ner)
      {
	assignScal(eqdModel, _epr[_nerC].eqd_model);
        assignStr(epref, _epr[_nerC].epref);
        _nerC++;
        return 1;
      }
    _nerC = 0;

    return 0;
  }

  void Poldata::commEemprops( t_commrec *cr)
  {
    int         i, j, nep;
    t_eemprops *ep;

    if (NULL != debug)
      {
        fprintf(debug, "Going to update eemprops on node %d\n", cr->nodeid);
      }
    if (MASTER(cr))
      {
        for (i = 1; (i < cr->nnodes); i++)
	  {
            gmx_send_int(cr, i, _nep);
            gmx_send(cr, i, vectorToArray(_eep), _nep*sizeof(_eep[0]));
	  }
      }
    else
      {
        nep = gmx_recv_int(cr, 0);
        if (nep != _nep)
	  {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
	  }
        snew(ep, _nep);
        gmx_recv(cr, 0, ep, _nep*sizeof(ep[0]));
        for (i = 0; (i < _nep); i++)
	  {
            _eep[i] = ep[i];
	  }
        sfree(ep);
      }
    if (NULL != debug)
      {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < nep); i++)
	  {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    getEemtypeName(_eep[i].eqd_model),
                    _eep[i].name, _eep[i].chi0, _eep[i].J0);
            for (j = 0; (j < _eep[i].nzeta); j++)
	      {
                fprintf(debug, " %8.3f", _eep[i].zeta[j]);
	      }
            fprintf(debug, "\n");
	  }
      }
  }

  void Poldata::commForceParameters(t_commrec *cr)
  {
    int         i, j, nep;
    t_eemprops *ep;

    if (NULL != debug)
      {
        fprintf(debug, "Going to update force parameters on node %d\n", cr->nodeid);
      }
    if (MASTER(cr))
      {
        for (i = 1; (i < cr->nnodes); i++)
	  {
            gmx_send_int(cr, i, _nep);
            gmx_send(cr, i, vectorToArray(_eep), _nep*sizeof(_eep[0]));
	  }
      }
    else
      {
        nep = gmx_recv_int(cr, 0);
        if (nep != nep)
	  {
            gmx_fatal(FARGS, "Inconsistency in number of EEM parameters");
	  }
        snew(ep, nep);
        gmx_recv(cr, 0, ep, _nep*sizeof(ep[0]));
        for (i = 0; (i < _nep); i++)
	  {
            _eep[i] = ep[i];
	  }
        sfree(ep);
      }
    if (NULL != debug)
      {
        fprintf(debug, "  EEP  Atom      Chi      J00     Zeta\n");
        for (i = 0; (i < _nep); i++)
	  {
            fprintf(debug, "%5s %5s %8.3f %8.3f",
                    getEemtypeName(_eep[i].eqd_model),
                    _eep[i].name, _eep[i].chi0, _eep[i].J0);
            for (j = 0; (j < _eep[i].nzeta); j++)
	      {
                fprintf(debug, " %8.3f", _eep[i].zeta[j]);
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

  ChargeDistributionModel Poldata::name2eemtype(const char *name)
  {
    int i;

    for (i = 0; (i < eqdNR); i++)
      {
        if (strcasecmp(name, eemtype_props[i].name) == 0)
	  {
            return eemtype_props[i].eqd;
	  }
      }
    return eqdNR;
  }

  const char *Poldata::getEemtypeName(ChargeDistributionModel eem)
  {
    int i;

    for (i = 0; (i < eqdNR); i++)
      {
        if (eem == eemtype_props[i].eqd)
	  {
            return eemtype_props[i].name;
	  }
      }

    return NULL;
  }

}
