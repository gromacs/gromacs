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
                 ChargeDistributionModel iDistributionModel, char **dzatoms)
  {
    int  k, zz;
    bool bRestr;

    bRestr = false;
    if (NULL != dzatoms)
      {
        k = 0;
        while ((NULL != dzatoms[k]) && !bRestr)
	  {
            bRestr = (strcasecmp(atomtype, dzatoms[k]) == 0);
            k++;
	  }
      }
    this->atomnumber  = atomnumber;
    this->atype       = atype;
    this->atomtype    = strdup(atomtype);
    this->nZeta       = pd->getNzeta(iDistributionModel, this->atomtype);
    if (this->nZeta <= 0)
      {
        bSetUpcorrectly = false;
	return;
      }

    this->bRestrained = bRestr;

    snew(this->zeta, this->nZeta);
    snew(this->zeta_ref, this->nZeta);
    snew(this->q, this->nZeta);
    snew(this->iz, this->nZeta);
    snew(this->iq, this->nZeta);
    snew(this->row, this->nZeta);

    for (zz = 0; (zz < this->nZeta); zz++)
      {
        this->iq[zz]       = -1;
        this->q[zz]        = pd->getQ( iDistributionModel, this->atomtype, zz);
        this->iz[zz]       = -1;
        this->zeta_ref[zz] =
	  this->zeta[zz] = pd->getZeta( iDistributionModel, this->atomtype, zz);
        this->row[zz]      = pd->getRow( iDistributionModel, this->atomtype, zz);
      }
    bSetUpcorrectly = true;
  }


  bool Ra::setUpcorrectly(){
    return bSetUpcorrectly;
  }

  Ra::~Ra()
  {
    sfree(this->iz);
    sfree(this->iq);
    sfree(this->zeta);
    sfree(this->q);
    sfree(this->row);
    sfree(this->atomtype);
  }

  real Ra::get_q()
  {
    int  i;
    real q = 0;

    for (i = 0; (i < this->nZeta); i++)
      {
        q += this->q[i];
      }

    return q;
  }

}
