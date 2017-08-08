/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author  Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author  David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef MYMOL_LOW_H
#define MYMOL_LOW_H

#include <assert.h>
#include <cstdio>
#include <cstring>

#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdtypes/forcerec.h"

#include "gentop_core.h"
#include "gentop_vsite.h"
#include "poldata.h"
#include "qgen_eem.h"
#include "qgen_resp.h"

namespace alexandria
{

class MyForceProvider : public IForceProvider
{
    private:
        std::vector<double> efield_;
        
    public:
    
        MyForceProvider();
        
        void calculateForces(const t_commrec  *cr,
                             const t_mdatoms  *mdatoms,
                             PaddedRVecVector *force,
                             double            t);
                         
        void setField(std::vector<double> field) { efield_ = field; }
};

bool is_planar(rvec xi, rvec xj, 
               rvec xk,  rvec xl, 
               t_pbc *pbc, real phi_toler);

bool is_linear(rvec xi, rvec xj, 
               rvec xk, t_pbc *pbc,
               real th_toler);

void add_excl(t_excls *excls, int e);

void add_excl_pair(t_excls excls[], int e1, int e2);

void remove_excl(t_excls *excls, int remove);

void prune_excl(t_excls excls[], t_atoms *atoms, gpp_atomtype_t atype);

void copy_atoms(t_atoms *src, t_atoms *dest);

void cp_plist(t_params                  *plist,
              int                        ftype,
              InteractionType            itype,
              std::vector<PlistWrapper> &plist_);
              
real calc_r13(const Poldata     &pd,
              const std::string  aai,
              const std::string  aaj,
              const std::string  aak,
              const real         angle);
              
real calc_relposition(const Poldata     &pd,
                      const std::string  aai,
                      const std::string  aaj,
                      const std::string  aak);  
                      
void updatePlist(const Poldata             &pd,
                 std::vector<PlistWrapper> &plist,
                 t_topology                *top);
                 
std::vector<double> getDoubles(const std::string &s);

void getLjParams(const Poldata     &pd,
                 const std::string &ai,
                 const std::string &aj,
                 double            *c6,
                 double            *cn);
                 
void getBhamParams(const Poldata     &pd,
                   const std::string &ai,
                   const std::string &aj,
                   double            *a,
                   double            *b,
                   double            *c);
                   
void plist_to_mtop(const Poldata             &pd,
                   std::vector<PlistWrapper>  plist,
                   gmx_mtop_t                *mtop_);
                   
void do_init_mtop(const Poldata            &pd,
                  gmx_mtop_t               *mtop,
                  char                    **molname,
                  t_atoms                  *atoms,
                  std::vector<PlistWrapper> plist,
                  t_inputrec               *ir,
                  t_symtab                 *symtab,
                  const char               *tabfn);
                  
void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka);

void mtop_update_cgs(gmx_mtop_t *mtop);

void put_in_box(int natom, matrix box, rvec x[], real dbox);

void write_zeta_q(FILE *fp, QgenEem * qgen,
                  t_atoms *atoms, ChargeDistributionModel iChargeDistributionModel);

void write_zeta_q2(QgenEem * qgen, gpp_atomtype_t atype,
                   t_atoms *atoms, ChargeDistributionModel iChargeDistributionModel);
                   
int get_subtype(directive d, int ftype);

void print_bondeds2(FILE                     *out,
                    directive                 d,
                    int                       plist_ftype,
                    int                       print_ftype,
                    std::vector<PlistWrapper> plist);

void write_top2(FILE *out, char *molname,
                t_atoms *at, gmx_bool bRTPresname,
                std::vector<PlistWrapper> plist_,
                t_excls excls[],
                gpp_atomtype_t atype, int *cgnr, int nrexcl,
                const Poldata &pd);
                
void print_top_header2(FILE *fp, const Poldata &pd,
                       gmx_atomprop_t aps, bool bPol,
                       std::vector<std::string> commercials,
                       bool bItp);

void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix);
                
}// namespace alexandria
#endif
