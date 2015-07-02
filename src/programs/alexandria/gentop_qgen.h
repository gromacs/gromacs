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
  
 private:
  gmx_bool                  bWarned;
  ChargeDistributionModel   iChargeDistributionModel;
  ChargeGenerationAlgorithm iChargeGenerationAlgorithm;
  int                       natom, eQGEN;
  real                      qtotal, chieq, hfac, epsr;
  /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
  char                    **elem;
  int                      *atomnr;
  real                     *chi0, *rhs, *j00;
  rvec                     *x;
  /* Jab is a matrix over atom pairs */
  real                    **Jab;
  /* For each atom i there are nZeta[i] row, q and zeta entries */
  int                      *nZeta;
  int                     **row;
  gmx_bool                  bAllocSave;
  real                    **q, **zeta, **qsave, **zetasave;
  
  
  real calc_jab(ChargeDistributionModel iChargeDistributionModel,
		rvec xi, rvec xj,
		int nZi, int nZj,
		real *zeta_i, real *zeta_j,
		int *rowi, int *rowj);
  void calc_Jab();
  
  void solve_q_eem(FILE *fp,  real hardsness_factor);
  
  void update_J00();
  
  real calc_Sij(int i, int j);

  void update_pd(t_atoms *atoms, Poldata * pd);
  
  
  int generate_charges_bultinck(FILE *fp,
				Poldata * pd, t_atoms *atoms,
				gmx_atomprop_t aps);



  void calc_rhs();


 
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
  
  int generate_charges_sm(FILE *fp,
			Poldata * pd, t_atoms *atoms,
			real tol, int maxiter, gmx_atomprop_t aps,
                    real *chieq);

  int generate_charges(FILE *fp,
		       gmx_resp_t gr, const char *molname,
		       Poldata * pd,
		       t_atoms *atoms,
		       real tol, int maxiter, int maxcycle,
		       gmx_atomprop_t aps);
  
  void message(int len, char buf[], gmx_resp_t gr);
  gmx_bool SplitQ(ChargeDistributionModel iDistributionModel);
  
  /* The routines below return NOTSET if something is out of the ordinary */
  int get_nzeta(int atom);
  
  int get_row(int atom, int z);
  
  double get_q(int atom, int z);

  void check_support(Poldata * pd, gmx_atomprop_t aps);
  
  void save_params( gmx_resp_t gr);
  
  void get_params( gmx_resp_t gr);
  double get_zeta(int atom, int z);
  

void print(FILE *fp, t_atoms *atoms);
  
  void debugFun(FILE *fp);
  
};
}
#endif
