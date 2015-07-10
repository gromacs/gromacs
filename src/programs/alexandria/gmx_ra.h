/*
 * This source file is part of the Alexandria project.
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
#ifndef GMX_RA_H
#define GMX_RA_H

#include "poldata.h"
#include "gmx_resp.h"


namespace alexandria
{
class Ra
{
public:
Ra(int atomnumber, int atype,
                 const char *atomtype, Poldata * pd,
     ChargeDistributionModel iDistributionModel, char **dzatoms);
~Ra();

real get_q();

bool setUpcorrectly();


//private:

int  *row, atomnumber, atype;
bool  bRestrained;
char *atomtype;
int   nZeta;
real *q, *zeta, *zeta_ref;
int  *iq, *iz;
bool bSetUpcorrectly;


};
}
#endif
