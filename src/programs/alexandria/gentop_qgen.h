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
#ifndef GENTOP_QGEN_H
#define GENTOP_QGEN_H

#include <stdio.h>
#include "poldata.h"
#include "gmx_resp.h"

enum {
  eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_NOSUPPORT, eQGEN_ERROR, eQGEN_NR
};

enum ChargeGenerationAlgorithm {
  eqgNONE, eqgEEM, eqgESP, eqgRESP, eqgNR
};

namespace alexandria
{

  class GentopQgen
  {
  
  public:
  
    ~GentopQgen();


    GentopQgen(Poldata * pd, t_atoms *atoms,
	       gmx_atomprop_t aps,
	       rvec *x,
	       ChargeDistributionModel   iChargeDistributionModel,
	       ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
	       real hfac, int qtotal,
	       real epsr);
  
    // void done();
  
    int generateChargesSm(FILE *fp,
			    Poldata * pd, t_atoms *atoms,
			    real tol, int maxiter, gmx_atomprop_t aps,
			    real *chieq);

    int generateCharges(FILE *fp,
			 Resp * gr, const char *molname,
			 Poldata * pd,
			 t_atoms *atoms,
			 real tol, int maxiter, int maxcycle,
			 gmx_atomprop_t aps);
  
    void message(int len, char buf[], Resp * gr);
    gmx_bool SplitQ(ChargeDistributionModel iDistributionModel);
  
    /* The routines below return NOTSET if something is out of the ordinary */
    int getNzeta(int atom);
  
    int getRow(int atom, int z);
  
    double getQ(int atom, int z);

    void checkSupport(Poldata * pd, gmx_atomprop_t aps);
  
    void saveParams( Resp * gr);
  
    void getParams( Resp *  gr);
    double getZeta(int atom, int z);
  

    void print(FILE *fp, t_atoms *atoms);
  
    void debugFun(FILE *fp);

  private:
    gmx_bool                  _bWarned;
    ChargeDistributionModel   _iChargeDistributionModel;
    ChargeGenerationAlgorithm _iChargeGenerationAlgorithm;
    int                       _natom, _eQGEN;
    real                      _qtotal, _chieq, _hfac, _epsr;
    /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
    std::vector<std::string>                    _elem;
    std::vector<int>                      _atomnr;
    std::vector<real>                     _chi0, _rhs, _j00;
    rvec *                     _x;
    /* Jab is a matrix over atom pairs */
    std::vector<std::vector<real> >                    _Jab;
    /* For each atom i there are nZeta[i] row, q and zeta entries */
    std::vector<int>                      _nZeta;
    std::vector<std::vector<int> >                     _row;
    gmx_bool                  _bAllocSave;
    std::vector<std::vector<real> >                   _q, _zeta, _qsave, _zetasave;
  
  
    real calcJab(ChargeDistributionModel iChargeDistributionModel,
		  rvec xi, rvec xj,
		  int nZi, int nZj,
		  std::vector<real> zeta_i,
		  std::vector<real> zeta_j,
		  std::vector<int> rowi, 
		  std::vector<int> rowj);

    void calcJab();
  
    void solveQEem(FILE *fp,  real hardsness_factor);
  
    void updateJ00();
  
    real calcSij(int i, int j);

    void updatePd(t_atoms *atoms, Poldata * pd);
  
    int generateChargesBultinck(FILE *fp,
				  Poldata * pd, t_atoms *atoms,
				  gmx_atomprop_t aps);

    void calcRhs();  
  };
}
#endif
