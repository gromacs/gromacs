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

    void update_atomtypes( t_atoms *atoms);

    void fill_zeta( alexandria::Poldata * pd);

    void fill_q( t_atoms *atoms);

    void add_atom_coords( rvec *x);

    bool add_atom_info( t_atoms *atoms,
			alexandria::Poldata * pd);

    void get_atom_info( t_atoms *atoms,
			t_symtab *symtab, rvec **x);

    const char *get_stoichiometry();

    void add_atom_symmetry(std::vector<int> &symmetric_atoms);

    void add_point( double x, double y,
		    double z, double V);

    void make_grid( real spacing, matrix box, rvec x[]);

    void copy_grid(Resp * src);

    Resp * copy();

    void calc_rms();

    double get_rms( real *wtot);

    void pot_lsq( gmx_stats_t lsq);

    void calc_rho();

    void calc_pot();

    void read_cube( const char *fn, bool bESPonly);

    void write_cube( const char *fn, char *title);

    void write_rho( const char *fn, char *title);

    void write_diff_cube(Resp * src,
			 const char *cube_fn, const char *hist_fn,
			 char *title, output_env_t oenv, int rho);

    void write_histo(const char *fn,
		     char *title, output_env_t oenv);

    int  optimize_charges(FILE *fp,  int maxiter,
			  real toler, real *rms);

    void potcomp(const char *potcomp,
		 const char *pdbdiff, output_env_t oenv);

    double get_qtot( int atom);

    double get_q( int atom, int zz);

    double get_zeta( int atom, int zz);

    void set_q( int atom, int zz, double q);

    void set_zeta( int atom, int zz, double zeta);

  private:
    ChargeDistributionModel iDistributionModel;
    int                     nesp, nrho, natom, natype, ngridp;
    double                  qtot, qsum, watoms;
    double                  rms, rrms, penalty, pfac, entropy, wtot;
    dvec                    origin, space;
    bool                    bZatype, bFitZeta, bEntropy;
    bool                    bRandZeta, bRandQ;
    bool                    bAXpRESP;
    ivec                    nxyz;
    real                    qfac, b_hyper, zmin, zmax, delta_z, qmin, qmax, rDecrZeta;
    unsigned int            seed;
    int                     nparam; /* Total number of parameters */
    Ra                 **ra;
    char                  **dzatoms;
    const char             *stoichiometry;
    double                 *pot, *pot_calc, *rho;
    rvec                   *x, *esp;

    
    void warning(const char *fn, int line);

    void get_set_vector(bool         bSet,
			bool         bRandQ,
			bool         bRandZeta,
			unsigned int seed,
			double      *nmx);

    real my_weight(int iatom);

    void calc_penalty();

    void add_param( int atom, eParm eparm, int zz);

    static double charge_function(void * gr,double v[]);

};
}
#endif
