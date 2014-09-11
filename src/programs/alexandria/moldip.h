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
#ifndef MOLDIP_H
#define MOLDIP_H

#include "gromacs/utility/real.h"
#include "mymol.h"

typedef struct {
    int       n, nopt, nconst, nopt_c;
    char    **name;
    int      *tot_count, *count;
    gmx_bool *bConst;
} t_index_count;

extern char *opt_index_count(t_index_count *ic);

enum {
    ermsBOUNDS, ermsMU, ermsQUAD, ermsCHARGE, ermsESP,
    ermsEPOT, ermsForce2, ermsTOT, ermsNR
};

namespace alexandria
{

class MolDip
{
    private:
    public:
        bool                           _bDone, _bFinal, _bGaussianBug;
        bool                           _bFitZeta;
        std::vector<alexandria::MyMol> _mymol;
        int                            _nmol_support;
        ChargeGenerationModel          _iModel;
        t_index_count                 *_ic;
        real                           _J0_0, _Chi0_0, _w_0, _J0_1, _Chi0_1, _w_1;
        real                           _hfac, _hfac0, _decrzeta, _epsr;
        real                           _ener[ermsNR], _fc[ermsNR];
        gmx_bool                       _bOptHfac, _bPol, _bQM;
        char                          *_fixchi;
        gmx_poldata_t                  _pd;
        gmx_atomprop_t                 _atomprop;
        t_commrec                     *_cr;

        //! Constructor
        MolDip();

        //! Destructor
        ~MolDip() {};

        void Init(t_commrec *cr, gmx_bool bQM, gmx_bool bGaussianBug,
                  ChargeGenerationModel iModel, real rDecrZeta, real epsr,
                  real J0_0, real Chi0_0, real w_0,
                  real J0_1, real Chi0_1, real w_1,
                  real fc_bound, real fc_mu, real fc_quad, real fc_charge,
                  real fc_esp, real fc_epot, real fc_force, char *fixchi,
                  gmx_bool bOptHfac, real hfac,
                  gmx_bool bPol, gmx_bool bFitZeta);
        void Read(FILE *fp, const char *fn, const char *pd_fn,
                  int minimum_data,
                  gmx_bool bZero,
                  char *opt_elem, char *const_elem,
                  char *lot,
                  output_env_t oenv, gmx_molselect_t gms,
                  real watoms, gmx_bool bCheckSupport, unsigned int seed);

        void CalcDeviation();
};

}

#endif
