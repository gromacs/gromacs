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
#ifndef GMX_RESP_H
#define GMX_RESP_H

#include <stdio.h>
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "poldata.h"
#include "molprop.h"
#include "gmx_ra.h"





enum eParm {
  eparmQ, eparmZ, eparmNR
};


namespace alexandria
{
  class Resp
  {
  public:
    Resp(){}

    Resp(ChargeDistributionModel iDistributionModel,
	 bool bAXpRESP, real qfac, real b_hyper, real qtot,
	 real zmin, real zmax, real delta_z, bool bZatyp,
	 real watoms, real rDecrZeta,
	 bool bRandZeta, bool bRandQ, real penalty_fac, bool bFitZeta,
	 bool bEntropy, const char *dzatoms,
	 unsigned int seed);


    ~Resp();


    int atomicnumber2row(int elem);



    void statistics( int len, char buf[]);

    void summary(FILE *gp, 
		 std::vector<int> &symmetric_atoms);

    void updateAtomtypes( t_atoms *atoms);

    void fillZeta( alexandria::Poldata * pd);

    void fillQ( t_atoms *atoms);

    void addAtomCoords( rvec *x);

    bool addAtomInfo( t_atoms *atoms,
			alexandria::Poldata * pd);

    void getAtomInfo( t_atoms *atoms,
			t_symtab *symtab, rvec **x);

    const char *getStoichiometry();

    void addAtomSymmetry(std::vector<int> &symmetric_atoms);

    void addPoint( double x, double y,
		    double z, double V);

    void makeGrid( real spacing, matrix box, rvec x[]);

    void copyGrid(Resp * src);

    Resp * copy();

    void calcRms();

    double getRms( real *wtot);

    void potLsq( gmx_stats_t lsq);

    void calcRho();

    void calcPot();

    void readCube( const char *fn, bool bESPonly);

    void writeCube( const char *fn, char *title);

    void writeRho( const char *fn, char *title);

    void writeDiffCube(Resp * src,
			 const char *cube_fn, const char *hist_fn,
			 char *title, output_env_t oenv, int rho);

    void writeHisto(const char *fn,
		     char *title, output_env_t oenv);

    int  optimizeCharges(FILE *fp,  int maxiter,
			  real toler, real *rms);

    void potcomp(const char *potcomp,
		 const char *pdbdiff, output_env_t oenv);

    double getQtot( int atom);

    double getQ( int atom, int zz);

    double getZeta( int atom, int zz);

    void setQ( int atom, int zz, double q);

    void setZeta( int atom, int zz, double zeta);

  private:
    ChargeDistributionModel _iDistributionModel;
    int                     _nesp, _natom, _natype;
    double                  _qtot, _qsum, _watoms;
    double                  _rms, _rrms, _penalty, _pfac, _entropy, _wtot;
    dvec                    _origin, _space;
    bool                    _bZatype, _bFitZeta, _bEntropy;
    bool                    _bRandZeta, _bRandQ;
    bool                    _bAXpRESP;
    ivec                    _nxyz;
    real                    _qfac, _bHyper, _zmin, _zmax, _deltaZ, _qmin, _qmax, _rDecrZeta;
    unsigned int            _seed;
    int                     _nparam; /* Total number of parameters */
    std::vector<Ra *>                 _ra;
    std::vector<std::string>                  _dzatoms;
    const std::string             _stoichiometry;
    std::vector<double>                 _pot, _potCalc, _rho;
    rvec                   *_x, *_esp;

    
    void warning(const char *fn, int line);

    void getSetVector(bool         bSet,
			bool         bRandQ,
			bool         bRandZeta,
			unsigned int seed,
			double      *nmx);

    real myWeight(int iatom);

    void calcPenalty();

    void addParam( int atom, eParm eparm, int zz);

    static double chargeFunction(void * gr,double v[]);

};
}
#endif
