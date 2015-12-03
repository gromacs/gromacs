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

#include "gmx_resp.h"
#include "poldata.h"


namespace alexandria
{
  class Ra
  {
  public:
    Ra(int atomnumber, int atype,
       const char *atomtype, Poldata * pd,
       ChargeDistributionModel iDistributionModel, std::vector<std::string>  dzatoms);
    ~Ra();

    real getQ();

    bool setUpcorrectly();

    real getZeta(int i){
      return _zeta[i];
    }

    void setZeta(int i, real value){
      _zeta[i] = value;
    }

    int getNZeta(){
      return _nZeta;
    }

    int getIq(int i){
      return _iq[i];
    }

    void setIq(int i,int value){
      _iq[i] = value;
    }

    int getIz(int i){
      return _iz[i];
    }

    void setIz(int i, int value){
      _iz[i] = value;
    }


    real getQ(int i){
      return _q[i];
    }

    void setQ(int i, real value){
      _q[i] = value;
    }

    int getRow(int i){
      return _row[i];
    }


    int getAtomnumber(){
      return _atomnumber;
    }

    void setAtomnumber(int value){
      _atomnumber = value;
    }

    real getZetaRef(int i){
      return _zetaRef[i];
    }

    std::string getAtomtype(){
      return _atomtype;
    }

    bool getBRestrained(){
      return _bRestrained;
    }

    int getAtype(){
      return _atype;
    }

    void setAtype(int i ){
      _atype = i;
    }

  private:

    std::vector<int>  _row; 
    int _atomnumber, _atype;
    bool  _bRestrained;
    std::string _atomtype;
    int   _nZeta;
    std::vector<real> _q, _zeta, _zetaRef;
    std::vector<int>  _iq, _iz;
    bool _bSetUpcorrectly;


  };
}
#endif
