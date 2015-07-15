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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "poldata.h"
#include "gmx_resp.h"
#include "gentop_qgen.h"
#include "stringutil.h"
#include "gmx_ra.h"




namespace alexandria
{

  Ra::Ra(int atomnumber, int atype,
                 const char *atomtype, Poldata * pd,
	 ChargeDistributionModel iDistributionModel, std::vector<std::string> dzatoms)
  {
    int  k, zz;
    bool bRestr;

    bRestr = false;
    if (!dzatoms.empty())
      {
        k = 0;
        while (("" != dzatoms[k]) && !bRestr)
	  {
            bRestr = (strcasecmp(atomtype, dzatoms[k].c_str()) == 0);
            k++;
	  }
      }
    _atomnumber  = atomnumber;
    _atype       = atype;
    _atomtype    = (atomtype);
    _nZeta       = pd->getNzeta(iDistributionModel, _atomtype.c_str());
    if (_nZeta <= 0)
      {
        _bSetUpcorrectly = false;
	return;
      }

    _bRestrained = bRestr;

    _zeta.resize(_nZeta);
    _zetaRef.resize(_nZeta);
    _q.resize(_nZeta);
    _iz.resize(_nZeta);
    _iq.resize(_nZeta);
    _row.resize(_nZeta);

    for (zz = 0; (zz < _nZeta); zz++)
      {
        _iq[zz]       = -1;
        _q[zz]        = pd->getQ( iDistributionModel, _atomtype.c_str(), zz);
        _iz[zz]       = -1;
        _zetaRef[zz] =
	  _zeta[zz] = pd->getZeta( iDistributionModel, _atomtype.c_str(), zz);
        _row[zz]      = pd->getRow( iDistributionModel, _atomtype.c_str(), zz);
      }
    _bSetUpcorrectly = true;
  }


  bool Ra::setUpcorrectly(){
    return _bSetUpcorrectly;
  }

  Ra::~Ra()
  {
    /* sfree(_iz);
    sfree(_iq);
    sfree(_zeta);
    sfree(_q);
    sfree(_row);*/
    // sfree(_atomtype);
  }

  real Ra::getQ()
  {
    int  i;
    real q = 0;

    for (i = 0; (i < _nZeta); i++)
      {
        q += _q[i];
      }

    return q;
  }
}
