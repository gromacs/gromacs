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

#include <assert.h>
#include <cstdio>
#include <cstring>

#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

#include "mymol_low.h"
#include "poldata.h"

namespace alexandria
{

MyForceProvider::MyForceProvider() {};

void MyForceProvider::calculateForces(const t_commrec  *cr,
                                      const t_mdatoms  *mdatoms,
                                      PaddedRVecVector *force,
                                      double            t)
{
    rvec *f = as_rvec_array(force->data());
    for(int dim = 0; dim < DIM; dim++)
    {       
        double efield = FIELDFAC*efield_[dim]; 
        if (efield != 0)
        {      
            for(int i = 0; i < mdatoms->nr; ++i)
            {
                f[i][dim] += mdatoms->chargeA[i]*efield;
            }
        }
    }
    if (MASTER(cr) && nullptr != debug)
    {
        fprintf(debug, "Electric Field. t: %4g  Ex: %4g  Ey: %4g  Ez: %4g\n", 
                t, efield_[XX], efield_[YY], efield_[ZZ]);
    }
}

bool is_planar(rvec xi, rvec xj, rvec xk, 
               rvec xl, t_pbc *pbc,
               real phi_toler)
{
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real sign, phi;

    phi = RAD2DEG*dih_angle(xi, xj, xk, xl, pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);

    return (fabs(phi) < phi_toler);
}

bool is_linear(rvec xi, rvec xj, 
               rvec xk, t_pbc *pbc,
               real th_toler)
{
    int  t1, t2;
    rvec r_ij, r_kj;
    real costh, th;

    th = fabs(RAD2DEG*bond_angle(xi, xj, xk, pbc, r_ij, r_kj, &costh, &t1, &t2));
    if ((th > th_toler) || (th < 180-th_toler))
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Angle is %g, th_toler is %g\n", th, th_toler);
        }
        return true;
    }
    return false;
}

void add_excl(t_excls *excls, int e)
{
    int i;

    for (i = 0; (i < excls->nr); i++)
    {
        if (excls->e[i] == e)
        {
            return;
        }
    }
    srenew(excls->e, excls->nr+1);
    excls->e[excls->nr++] = e;
}

void add_excl_pair(t_excls excls[], int e1, int e2)
{
    if (e1 != e2)
    {
        add_excl(&excls[e1], e2);
        add_excl(&excls[e2], e1);
    }
}

void remove_excl(t_excls *excls, int remove)
{
    int i;

    for (i = remove+1; i < excls->nr; i++)
    {
        excls->e[i-1] = excls->e[i];
    }

    excls->nr--;
}

void let_shells_see_shells(t_excls excls[], t_atoms *atoms, gpp_atomtype_t atype)
{
    int i, k, ak;
    for (i = 0; i < atoms->nr; i++)
    {
        if (get_atomtype_ptype(atoms->atom[i].type, atype) == eptShell)
        {
            for (k = 0; k < excls[i].nr; )
            {
                ak = excls[i].e[k];
                if (get_atomtype_ptype(atoms->atom[ak].type, atype) == eptShell)
                {
                    remove_excl(&(excls[i]), k);
                }
                else
                {
                    k++;
                }
            }
        }
    }
}

void copy_atoms(t_atoms *src, t_atoms *dest)
{
    int i;

    if (dest->nr < src->nr)
    {
        srenew(dest->atom, src->nr);
        srenew(dest->atomname, src->nr);
        if (nullptr != src->atomtype)
        {
            srenew(dest->atomtype, src->nr);
        }
        else if (nullptr != dest->atomtype)
        {
            sfree(dest->atomtype);
            dest->atomtype = nullptr;
        }
        if (nullptr != src->atomtypeB)
        {
            srenew(dest->atomtypeB, src->nr);
        }
        else if (nullptr != dest->atomtypeB)
        {
            sfree(dest->atomtypeB);
            dest->atomtypeB = nullptr;
        }
    }
    dest->nr = src->nr;
    for (i = 0; (i < src->nr); i++)
    {
        dest->atom[i]      = src->atom[i];
        dest->atomname[i]  = src->atomname[i];
        if (nullptr != src->atomtype)
        {
            dest->atomtype[i]  = src->atomtype[i];
        }
        if (nullptr != src->atomtypeB)
        {
            dest->atomtypeB[i] = src->atomtypeB[i];
        }
    }
    if (dest->nres < src->nres)
    {
        srenew(dest->resinfo, src->nres);
    }

    if (nullptr != src->pdbinfo)
    {
        srenew(dest->pdbinfo, src->nres);
    }
    else if (nullptr != dest->pdbinfo)
    {
        sfree(dest->pdbinfo);
        dest->pdbinfo = nullptr;
    }
    dest->nres = src->nres;
    for (i = 0; (i < src->nres); i++)
    {
        dest->resinfo[i] = src->resinfo[i];
        if (nullptr != src->pdbinfo)
        {
            dest->pdbinfo[i] = src->pdbinfo[i];
        }
    }
}

void cp_plist(t_params                  *plist,
              int                        ftype,
              InteractionType            itype,
              std::vector<PlistWrapper> &plist_)
{
    if (plist->nr > 0)
    {
        PlistWrapper pw(itype, ftype);
        for (int i = 0; (i < plist->nr); i++)
        {
            for (int j = interaction_function[ftype].nratoms; j < MAXATOMLIST; j++)
            {
                plist->param[i].a[j] = 0;
            }
            for (int j = interaction_function[ftype].nrfpA; j < MAXFORCEPARAM; j++)
            {
                plist->param[i].c[j] = 0;
            }
            pw.addParam(plist->param[i]);
        }
        plist_.push_back(pw);
    }
}

real calc_r13(const Poldata     &pd,
              const std::string  aai,
              const std::string  aaj,
              const std::string  aak,
              const real         angle)
{
    std::string              params;
    size_t                   ntrain;
    double                   sigma;
    double                   rij = 0, rjk = 0;
    double                   r12 = 0, r23 = 0;
    double                   r13 = 0;

    std::vector<std::string> aij = {aai, aaj};
    std::vector<std::string> ajk = {aaj, aak};

    auto                     fs = pd.findForces(eitBONDS);
    auto                     lu = string2unit(fs->unit().c_str());

    pd.searchForce(aij, params, &rij, &sigma, &ntrain);
    pd.searchForce(ajk, params, &rjk, &sigma, &ntrain);

    r12 = convert2gmx(rij, lu);
    r23 = convert2gmx(rjk, lu);

    r13 = std::sqrt((r12*r12) + (r23*r23) - (2*r12*r23*std::cos(DEG2RAD*angle)));

    return r13;
}

real calc_relposition(const Poldata     &pd,
                      const std::string  aai,
                      const std::string  aaj,
                      const std::string  aak)
{
    std::string              params;
    size_t                   ntrain;
    double                   sigma;
    double                   rij               = 0, rjk = 0;
    double                   b0                = 0, b1 = 0;
    double                   relative_position = 0;

    std::vector<std::string> aij = {aai, aaj};
    std::vector<std::string> ajk = {aaj, aak};

    auto                     fs = pd.findForces(eitBONDS);
    auto                     lu = string2unit(fs->unit().c_str());

    pd.searchForce(aij, params, &rij, &sigma, &ntrain);
    pd.searchForce(ajk, params, &rjk, &sigma, &ntrain);

    b0 = convert2gmx(rij, lu);
    b1 = convert2gmx(rjk, lu);

    relative_position = (b1/(b0+b1));

    return relative_position;
}

immStatus updatePlist(const Poldata             &pd,
                      std::vector<PlistWrapper> &plist,
                      t_topology                *top,
                      bool                       bBASTAT)
{
    std::string              aai, aaj, aak, aal, params;
    std::vector<std::string> atoms, ptr;
    int                      lu, n;
    size_t                   ntrain;
    double                   value, sigma, r13 = 0;

    for (auto &pw : plist)
    {
        auto iType = pw.getItype();
        auto fs    = pd.findForces(iType);
        
        if (fs != pd.forcesEnd())
        {
            pw.setFtype(fs->fType());
            
            if (eitBONDS == iType)
            {
                lu = string2unit(fs->unit().c_str());
                for (auto b = pw.beginParam(); b < pw.endParam(); ++b)
                {
                    if (pd.atypeToBtype(*top->atoms.atomtype[b->a[0]], aai) &&
                        pd.atypeToBtype(*top->atoms.atomtype[b->a[1]], aaj))
                    {
                        atoms = {aai, aaj};
                        n     = 0;
                        if ((fs->searchForce(atoms, params, &value, &sigma, &ntrain)) != 0)
                        {
                            b->c[n++] = convert2gmx(value, lu);
                            ptr       = gmx::splitString(params);
                            for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                            {
                                b->c[n++] = gmx::doubleFromString(pi->c_str());
                            }
                        }
                        else if (!bBASTAT)
                        {
			    if (debug)
			    {
				fprintf(debug, "Could not find bond information for %s - %s\n",
					aai.c_str(), aaj.c_str());
			    }
                            return immNotSupportedBond;
                        }
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Unsupported atom types: %d, %d!\n",
                                  b->a[0], b->a[1]);
                    }
                }
            }
            else if (eitANGLES == iType ||
                     eitLINEAR_ANGLES == iType)
            {
                for (auto b = pw.beginParam(); b < pw.endParam(); ++b)
                {
                    if (pd.atypeToBtype(*top->atoms.atomtype[b->a[0]], aai) &&
                        pd.atypeToBtype(*top->atoms.atomtype[b->a[1]], aaj) &&
                        pd.atypeToBtype(*top->atoms.atomtype[b->a[2]], aak))
                    {
                        atoms = {aai, aaj, aak};
                        n     = 0;
                        if ((fs->searchForce(atoms, params, &value, &sigma, &ntrain)) != 0)
                        {
                            r13 = calc_r13(pd, aai, aaj, aak, value);
                            
                            b->c[n++] = value;
                            ptr       = gmx::splitString(params);
                            for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                            {
                                b->c[n++] = gmx::doubleFromString(pi->c_str());
                                if (n == 2)
                                {
                                    b->c[n++] = r13;
                                }
                            }
                        }
                        else if (!bBASTAT)
                        {
                            return immNotSupportedAngle;
                        }
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Unsuppotred atom types: %d, %d, %d!\n",
                                  b->a[0], b->a[1], b->a[2]);
                    }
                }
            }
            else if (eitPROPER_DIHEDRALS == iType ||
                     eitIMPROPER_DIHEDRALS == iType)
            {
                for (auto b = pw.beginParam(); b < pw.endParam(); ++b)
                {
                    if (pd.atypeToBtype(*top->atoms.atomtype[b->a[0]], aai) &&
                        pd.atypeToBtype(*top->atoms.atomtype[b->a[1]], aaj) &&
                        pd.atypeToBtype(*top->atoms.atomtype[b->a[2]], aak) &&
                        pd.atypeToBtype(*top->atoms.atomtype[b->a[3]], aal))
                    {
                        atoms = {aai, aaj, aak, aal};
                        n     = 0;
                        if ((fs->searchForce(atoms, params, &value, &sigma, &ntrain)) != 0)
                        {
                            b->c[n++] = value;
                            ptr       = gmx::splitString(params);
                            int n = 0;
                            for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                            {
                                if (n == 0)
                                {
                                    b->c[n++] = gmx::doubleFromString(pi->c_str());
                                }
                                else
                                {
                                    /*Multiplicity for Proper Dihedral must be integer
                                      This assumes that the second paramter is Multiplicity*/
                                    b->c[n++] = atoi(pi->c_str());
                                }
                            }
                        }
                        else if (!bBASTAT)
                        {
                            return immNotSupportedDihedral;
                        }
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Unsuppotred atom types: %d, %d, %d, %d!\n",
                                  b->a[0], b->a[1], b->a[2], b->a[3]);
                    }
                }
            }
        }
    }
    return immOK;
}

std::vector<double> getDoubles(const std::string &s)
{
    std::vector<double> d;

    for (auto &ss : gmx::splitString(s))
    {
        d.push_back(gmx::doubleFromString(ss.c_str()));
    }
    return d;
}

void getLjParams(const Poldata     &pd,
                 const std::string &ai,
                 const std::string &aj,
                 double            *c6,
                 double            *cn)
{
    std::vector<double> vdwi, vdwj;

    auto                fai = pd.findAtype(ai);
    if (fai != pd.getAtypeEnd())
    {
        vdwi  = getDoubles(fai->getVdwparams());
    }
    else
    {
        vdwi.resize(2, 0.0);
    }
    auto faj = pd.findAtype(aj);
    if (faj != pd.getAtypeEnd())
    {
        vdwj  = getDoubles(faj->getVdwparams());
    }
    else
    {
        vdwj.resize(2, 0.0);
    }

    switch (pd.getCombRule())
    {
        case eCOMB_GEOMETRIC:
            *c6 = std::sqrt((vdwi[0]) * (vdwj[0]));
            *cn = std::sqrt((vdwi[1]) * (vdwj[1]));
            break;
        case eCOMB_ARITHMETIC:
        {
            double sig  = 0.5 * ((vdwi[0]) + (vdwj[0]));
            double eps  = std::sqrt((vdwi[1]) + (vdwj[1]));
            double sig6 = std::pow(sig, 6.0);
            *c6 = 4*eps*sig6;
            *cn = *c6 * sig6;
        }
        break;
        case eCOMB_GEOM_SIG_EPS:
        {
            double sig  = std::sqrt((vdwi[0]) * (vdwj[0]));
            double eps  = std::sqrt((vdwi[1]) * (vdwj[1]));
            double sig6 = std::pow(sig, 6.0);
            *c6 = 4*eps*sig6;
            *cn = *c6 * sig6;
        }
        break;
        case eCOMB_NONE:
        case eCOMB_NR:
            gmx_fatal(FARGS, "Unsupported combination rule for Lennard Jones");
    }
}

void getBhamParams(const Poldata     &pd,
                   const std::string &ai,
                   const std::string &aj,
                   double            *a,
                   double            *b,
                   double            *c)
{
    std::vector<double> vdwi, vdwj;

    auto                fai = pd.findAtype(ai);
    if (fai != pd.getAtypeEnd())
    {
        vdwi  = getDoubles(fai->getVdwparams());
    }
    else
    {
        vdwi.resize(3, 0.0);
    }
    auto faj = pd.findAtype(aj);
    if (faj != pd.getAtypeEnd())
    {
        vdwj  = getDoubles(faj->getVdwparams());
    }
    else
    {
        vdwj.resize(3, 0.0);
    }

    switch (pd.getCombRule())
    {
        case eCOMB_GEOMETRIC:
            *a = std::sqrt((vdwi[0]) * (vdwj[0]));
            *b = std::sqrt((vdwi[1]) * (vdwj[1]));
            *c = std::sqrt((vdwi[2]) * (vdwj[2]));
            break;
        case eCOMB_ARITHMETIC:
            *a = 0.5 * ((vdwi[0]) + (vdwj[0]));
            *b = std::sqrt((vdwi[1]) * (vdwj[1]));
            *c = 0.5 * ((vdwi[2]) + (vdwj[2]));
            break;
        case eCOMB_KONG_MASON:
            *a = std::sqrt((vdwi[0]) * (vdwj[0]));
            *b = 2.0*(vdwi[1] * vdwj[1])/(vdwi[1] + vdwj[1]);
            *c = 0.25*((vdwi[2]/vdwi[0])+(vdwj[2]/vdwj[0]))*((vdwi[0]) + (vdwj[0]));
            break;
        case eCOMB_GEOM_SIG_EPS:
        case eCOMB_NONE:
        case eCOMB_NR:
            gmx_fatal(FARGS, "Unsupported combination rule for Buckingham");
    }
}

void plist_to_mtop(const Poldata             &pd,
                   std::vector<PlistWrapper>  plist,
                   gmx_mtop_t                *mtop_)
{
    double fudgeLJ;
    double reppow = 12.0;
    int    n      = 0;

    /* Generate pairs */
    fudgeLJ = pd.getFudgeLJ();

    int nfptot = mtop_->ffparams.ntypes;
    for (auto &pw : plist)
    {
        nfptot += pw.nParam()*NRFPA(pw.getFtype());
    }
    srenew(mtop_->ffparams.functype, nfptot);
    srenew(mtop_->ffparams.iparams, nfptot);
    for (int i = mtop_->ffparams.ntypes; i < nfptot; i++)
    {
        mtop_->ffparams.functype[i] = 0;
        memset(&mtop_->ffparams.iparams[i], 0, sizeof(mtop_->ffparams.iparams[i]));
    }

    for (auto &pw : plist)
    {
        int ftype  = pw.getFtype();
        int nra    = NRAL(ftype);
        int nrfp   = NRFPA(ftype);
        int nratot = pw.nParam()*(1+nra);
        if (nratot > 0 && debug)
        {
            fprintf(debug, "There are %d interactions of type %s\n", nratot/(nra+1),
                    interaction_function[ftype].name);
        }
        snew(mtop_->moltype[0].ilist[ftype].iatoms, nratot);
        int k = 0;
        for (auto j = pw.beginParam(); (j < pw.endParam()); ++j)
        {
            std::vector<real> c;
            c.resize(MAXFORCEPARAM, 0);
            int               l = 0;
            if (ftype == F_LJ14)
            {
                int ati = mtop_->moltype[0].atoms.atom[j->a[0]].type;
                int atj = mtop_->moltype[0].atoms.atom[j->a[1]].type;
                int tp  = ati*mtop_->ffparams.atnr+atj;
                c[l++] = mtop_->ffparams.iparams[tp].lj.c6*fudgeLJ;
                c[l++] = mtop_->ffparams.iparams[tp].lj.c12*fudgeLJ;
            }
            else
            {
                for (; (l < nrfp); l++)
                {
                    c[l] = j->c[l];
                    if (NOTSET == c[l])
                    {
                        c[l] = 0;
                    }
                }
            }
            for (; (l < MAXFORCEPARAM); l++)
            {
                c[l] = 0;
            }
            n = enter_params(&mtop_->ffparams, ftype, c.data(), 0, reppow, n, TRUE);
            mtop_->moltype[0].ilist[ftype].iatoms[k++] = n;
            for (l = 0; (l < nra); l++)
            {
                mtop_->moltype[0].ilist[ftype].iatoms[k++] = j->a[l];
            }
        }
        mtop_->moltype[0].ilist[ftype].nr = k;
    }
}

void do_init_mtop(const Poldata            &pd,
                  gmx_mtop_t               *mtop,
                  char                    **molname,
                  t_atoms                  *atoms,
                  std::vector<PlistWrapper> plist,
                  t_inputrec               *ir,
                  t_symtab                 *symtab,
                  const char               *tabfn)
{

    init_mtop(mtop);
    mtop->name     = molname;
    mtop->nmoltype = 1;
    snew(mtop->moltype, mtop->nmoltype);
    mtop->moltype[0].name = molname;
    mtop->nmolblock       = 1;
    snew(mtop->molblock, mtop->nmolblock);
    mtop->molblock[0].nmol        = 1;
    mtop->molblock[0].type        = 0;
    mtop->groups.grps[egcENER].nr = 1;
    mtop->molblock[0].natoms_mol  = atoms->nr;
    mtop->natoms                  = atoms->nr;
    init_t_atoms(&(mtop->moltype[0].atoms), atoms->nr, false);
    

    /*Count the number of atom types in the molecule*/
    int ntype      = 0;
    for (int i = 0; (i < atoms->nr); i++)
    {
        bool found = false;
        int  itp   = atoms->atom[i].type;
        mtop->moltype[0].atoms.atom[i] = atoms->atom[i];
        for (int j = 0; !found && (j < i); j++)
        {
            found = (itp == atoms->atom[j].type);
        }
        if (!found)
        {
            ntype++;
        }
    }

    snew(mtop->groups.grpname, ntype);
    snew(mtop->groups.grps[egcENER].nm_ind, ntype);

    int  ind = 0;
    for (int i = 0; i < atoms->nr; i++)
    {
        char  *atp   = *atoms->atomtype[i];
        bool   found = false;
        for (int j = 0; !found && (j < i); j++)
        {
            found = (strcmp(atp, *atoms->atomtype[j]) == 0);
        }
        if (!found)
        {
            mtop->groups.grpname[ind]              = put_symtab(symtab, atp);
            mtop->groups.grps[egcENER].nm_ind[ind] = ind;
            ind++;
        }
    }

    mtop->ffparams.atnr             = ntype;
    mtop->ffparams.ntypes           = ntype*ntype;   
    mtop->ffparams.reppow           = 12;  
    
    if (nullptr != tabfn)
    {
        mtop->groups.grps[egcENER].nr   = ntype;
        ir->opts.ngener                 = ntype;
        srenew(ir->opts.egp_flags, ntype*ntype);
        for (int k = 0; k < ntype; k++)
        {
            for (int m = k; m < ntype; m++)
            {
                ir->opts.egp_flags[k*ntype + m] |= EGP_TABLE;
            }
        }
    }
    
    int vdw_type = pd.getVdwFtype();
    snew(mtop->ffparams.functype, mtop->ffparams.ntypes);
    snew(mtop->ffparams.iparams, mtop->ffparams.ntypes);   
    for (int i = 0; (i < ntype); i++)
    {
        for (int j = 0; (j < ntype); j++)
        {
            int idx = ntype*i+j;
            mtop->ffparams.functype[idx] = vdw_type;
            switch (vdw_type)
            {
                case F_LJ:
                {
                    double c6  = 0;
                    double c12 = 0;
                    if (atoms->atom[i].ptype != eptShell &&
                        atoms->atom[i].ptype != eptVSite &&
                        atoms->atom[j].ptype != eptShell &&
                        atoms->atom[j].ptype != eptVSite)
                    {
                        getLjParams(pd,
                                    *(atoms->atomtype[i]),
                                    *(atoms->atomtype[j]),
                                    &c6, &c12);
                    }
                    mtop->ffparams.iparams[idx].lj.c6  = c6;
                    mtop->ffparams.iparams[idx].lj.c12 = c12;
                }
                break;
                case F_BHAM:
                {
                    double a = 0;
                    double b = 0;
                    double c = 0;
                    if (atoms->atom[i].ptype != eptShell &&
                        atoms->atom[i].ptype != eptVSite &&
                        atoms->atom[j].ptype != eptShell &&
                        atoms->atom[j].ptype != eptVSite)
                    {
                        getBhamParams(pd,
                                      *(atoms->atomtype[i]),
                                      *(atoms->atomtype[j]),
                                      &a, &b, &c);
                    }
                    mtop->ffparams.iparams[idx].bham.a = a;
                    mtop->ffparams.iparams[idx].bham.b = b;
                    mtop->ffparams.iparams[idx].bham.c = c;
                }
                break;
                default:
                    fprintf(stderr, "Invalid van der waals type %s\n",
                            pd.getVdwFunction().c_str());
            }
        }
    }   
    gmx_mtop_finalize(mtop);
    /* Create a charge group block */
    stupid_fill_block(&(mtop->moltype[0].cgs), atoms->nr, false);
    plist_to_mtop(pd, plist, mtop);
}

void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka)
{
    int i, j, k, nra;

    if (blocka->nr < natom)
    {
        srenew(blocka->index, natom+1);
        for (int i = blocka->nr; i < natom+1; i++)
        {
            blocka->index[i] = 0;
        }
    }
    nra = 0;
    for (i = 0; (i < natom); i++)
    {
        nra += excls_[i].nr;
    }
    snew(blocka->a, nra+1);
    nra = 0;
    for (i = j = 0; (i < natom); i++)
    {
        blocka->index[i] = nra;
        for (k = 0; (k < excls_[i].nr); k++)
        {
            blocka->a[j++] = excls_[i].e[k];
        }
        nra += excls_[i].nr;
    }
    blocka->index[natom] = nra;
    blocka->nr           = natom;
    blocka->nra          = nra;
}

void mtop_update_cgs(gmx_mtop_t *mtop)
{
    int i, j;

    for (i = 0; (i < mtop->nmoltype); i++)
    {
        if (mtop->moltype[i].atoms.nr > mtop->moltype[i].cgs.nr)
        {
            mtop->moltype[i].cgs.nr           = mtop->moltype[i].atoms.nr;
            mtop->moltype[i].cgs.nalloc_index = mtop->moltype[i].atoms.nr+1;
            srenew(mtop->moltype[i].cgs.index, mtop->moltype[i].cgs.nr+1);
            for (j = 0; (j <= mtop->moltype[i].cgs.nr); j++)
            {
                mtop->moltype[i].cgs.index[j] = j;
            }
        }
    }
}

void put_in_box(int natom, matrix box, rvec x[], real dbox)
{
    int  i, m;
    rvec xmin, xmax, xcom;

    clear_rvec(xcom);
    copy_rvec(x[0], xmin);
    copy_rvec(x[0], xmax);
    for (i = 0; (i < natom); i++)
    {
        rvec_inc(xcom, x[i]);
        for (m = 0; (m < DIM); m++)
        {
            if (xmin[m] > x[i][m])
            {
                xmin[m] = x[i][m];
            }
            else if (xmax[m] < x[i][m])
            {
                xmax[m] = x[i][m];
            }
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        xcom[m]  /= natom;
        box[m][m] = (dbox+xmax[m]-xmin[m]);
    }
}

void write_zeta_q(FILE                   *fp, 
                  QgenEem                *qgen,
                  t_atoms                *atoms, 
                  ChargeDistributionModel iChargeDistributionModel)
{
    int    i, ii, j, k, nz, row;
    double zeta, q;
    bool   bAtom, bTypeSet;

    if (nullptr == qgen)
    {
        return;
    }

    fprintf(fp, "[ charge_spreading ]\n");
    fprintf(fp, "; This section describes additional atom type properties.\n");
    fprintf(fp, "; Spreading type (stype) can be either Gaussian (AXg) or Slater (AXs).\n");
    fprintf(fp, "; The zeta are the same for atoms of the same type, and all but the last\n");
    fprintf(fp, "; charge as well. The final charge is different between atoms however,\n");
    fprintf(fp, "; and it is listed below in the [ atoms ] section.\n");
    fprintf(fp, "; atype stype  nq%s      zeta          q  ...\n",
            (iChargeDistributionModel == eqdAXs) ? "  row" : "");

    k = -1;
    for (i = 0; (i < atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
        {
            k++;
        }
        if (k == -1)
        {
            gmx_fatal(FARGS, "The first atom must be a real atom, not a shell");
        }
        nz = qgen->getNzeta( k);
        if (nz != NOTSET)
        {
            bTypeSet = false;
            for (ii = 0; !bTypeSet && (ii < i); ii++)
            {
                bTypeSet = (atoms->atom[ii].type == atoms->atom[i].type);
            }
            if (!bTypeSet)
            {
                fprintf(fp, "%5s %6s %3d",
                        *atoms->atomtype[i],
                        getEemtypeName(iChargeDistributionModel), (bAtom) ? nz : 1);
            }
            for (j = (bAtom ? 0 : nz); (j < (bAtom ? nz : nz)); j++)
            {
                row   = qgen->getRow( k, j);
                q     = qgen->getQ( k, j);
                zeta  = qgen->getZeta( k, j);
                if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET))
                {
                    if (j == nz-1)
                    {
                        atoms->atom[i].q      =
                            atoms->atom[i].qB = q;
                    }
                    if (!bTypeSet)
                    {
                        if (iChargeDistributionModel == eqdAXs)
                        {
                            fprintf(fp, "  %4d", row);
                        }
                        fprintf(fp, " %10f", zeta);
                        if (j < nz-1)
                        {
                            fprintf(fp, " %10f", q);
                        }
                    }
                }
            }
            if (!bTypeSet)
            {
                fprintf(fp, "\n");
            }
        }
    }
    fprintf(fp, "\n");
}

void write_zeta_q2(QgenEem                *qgen, 
                   gpp_atomtype_t          atype,
                   t_atoms                *atoms, 
                   ChargeDistributionModel iChargeDistributionModel)
{
    FILE      *fp;
    int        i, j, k, nz, row;
    double     zeta, q, qtot;
    gmx_bool   bAtom;

    if (nullptr == qgen)
    {
        return;
    }

    fp = fopen("zeta_q.txt", "w");
    k  = -1;
    for (i = 0; (i < atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
        {
            k++;
        }
        if (k == -1)
        {
            gmx_fatal(FARGS, "The first atom must be a real atom, not a shell");
        }
        nz = qgen->getNzeta( k);
        if (nz != NOTSET)
        {
            fprintf(fp, "%6s  %5s  %5d", getEemtypeName(iChargeDistributionModel),
                    get_atomtype_name(atoms->atom[i].type, atype),
                    (bAtom) ? nz-1 : 1);
            qtot = 0;
            for (j = (bAtom ? 0 : nz-1); (j < (bAtom ? nz-1 : nz)); j++)
            {
                row   = qgen->getRow( k, j);
                q     = qgen->getQ( k, j);
                zeta  = qgen->getZeta( k, j);
                if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET))
                {
                    qtot += q;
                    fprintf(fp, "%5d %10g %10g", row, zeta, q);
                }
            }
            atoms->atom[i].q = qtot;
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "\n");
    fclose(fp);
}

int get_subtype(directive d, int ftype)
{
    int i;
    for (i = 1; i < 20; i++)
    {
        if (ifunc_index(d, i) == ftype)
        {
            return i;
        }
    }
    return 1;
}

void print_bondeds(FILE                     *out,
                   directive                 d,
                   int                       plist_ftype,
                   int                       print_ftype,
                   std::vector<PlistWrapper> plist)
{
    auto p = SearchPlist(plist, plist_ftype);

    if (plist.end() == p || p->nParam() == 0)
    {
        return;
    }
    fprintf(out, "[ %s ]\n", dir2str(d));
    fprintf(out, ";atom i");
    for (int j = 1; (j < NRAL(print_ftype)); j++)
    {
        fprintf(out, "  %5c", j+'i');
    }
    fprintf(out, "   type  parameters\n");
    int subtype = get_subtype(d, print_ftype);
    for (auto i = p->beginParam(); (i < p->endParam()); ++i)
    {
        for (int j = 0; (j < NRAL(print_ftype)); j++)
        {
            fprintf(out, "  %5d", 1+i->a[j]);
        }
        fprintf(out, "  %5d", subtype);
        for (int j = 0; (j < NRFPA(print_ftype)); j++)
        {
            fprintf(out, "  %10g", i->c[j]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

void write_top(FILE                     *out, 
               char                     *molname,
               t_atoms                  *at, 
               gmx_bool                  bRTPresname,
               std::vector<PlistWrapper> plist_,
               t_excls                   excls[],
               gpp_atomtype_t            atype, 
               int                      *cgnr, 
               int                       nrexcl,
               const Poldata            &pd)
{
    if (at && atype && cgnr)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nrexcl);
        print_atoms(out, atype, at, cgnr, bRTPresname);
        for (auto fs = pd.forcesBegin(); fs != pd.forcesEnd(); fs++)
        {
            if (eitBONDS == fs->iType())
            {
                print_bondeds(out, d_bonds, fs->fType(), fs->fType(), plist_);
            }
            else if (eitANGLES == fs->iType() || eitLINEAR_ANGLES == fs->iType())
            {
                print_bondeds(out, d_angles, fs->fType(), fs->fType(), plist_);
            }
            else if (eitPROPER_DIHEDRALS == fs->iType() || eitIMPROPER_DIHEDRALS == fs->iType())
            {
                print_bondeds(out, d_dihedrals, fs->fType(), fs->fType(), plist_);
            }
        }
        print_bondeds(out, d_constraints, F_CONSTR, F_CONSTR, plist_);
        print_bondeds(out, d_constraints, F_CONSTRNC, F_CONSTRNC, plist_);
        print_bondeds(out, d_pairs, F_LJ14, F_LJ14, plist_);
        print_excl(out, at->nr, excls);
        print_bondeds(out, d_cmap, F_CMAP, F_CMAP, plist_);
        print_bondeds(out, d_polarization, F_POLARIZATION, F_POLARIZATION, plist_);
        print_bondeds(out, d_thole_polarization, F_THOLE_POL, F_THOLE_POL, plist_);
        print_bondeds(out, d_vsites2, F_VSITE2, F_VSITE2, plist_);
        print_bondeds(out, d_vsites3, F_VSITE3, F_VSITE3, plist_);
        print_bondeds(out, d_vsites3, F_VSITE3FD, F_VSITE3FD, plist_);
        print_bondeds(out, d_vsites3, F_VSITE3FAD, F_VSITE3FAD, plist_);
        print_bondeds(out, d_vsites3, F_VSITE3OUT, F_VSITE3OUT, plist_);
        print_bondeds(out, d_vsites4, F_VSITE4FD, F_VSITE4FD, plist_);
        print_bondeds(out, d_vsites4, F_VSITE4FDN, F_VSITE4FDN, plist_);
    }
}

void print_top_header(FILE                    *fp, 
                      const Poldata           &pd,
                      gmx_atomprop_t           aps, 
                      bool                     bPol,
                      ChargeDistributionModel  iChargeDistributionModel,
                      std::vector<std::string> commercials,
                      bool                     bItp)
{
    std::string   gt_old, gt_type;
    std::string   btype;
    int           atomnumber;
    real          mass;

    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
    fprintf(fp, "; Watch this space for information & commercials.\n");
    for (auto i = commercials.begin(); (i < commercials.end()); ++i)
    {
        fprintf(fp, "; %s\n", i->c_str());
    }
    fprintf(fp, ";\n");
    if (!bItp)
    {
        fprintf(fp, "[ defaults ]\n");
        fprintf(fp, "; nbfunc         comb-rule       gen-pairs       fudgeLJ     fudgeQQ\n");
        std::string ff = pd.getVdwFunction();
        if (strcasecmp(ff.c_str(), "LJ_SR") == 0)
        {
            ff = "LJ";
        }
        fprintf(fp, "%-15s  %-15s no           %10g  %10g\n\n",
                ff.c_str(),
                pd.getCombinationRule().c_str(),
                pd.getFudgeLJ(),
                pd.getFudgeQQ());

        fprintf(fp, "[ atomtypes ]\n");
        fprintf(fp, "%-7s%-6s  %6s  %11s  %10s  %5s %-s  %s\n",
                ";atype ", "btype", "at.num", "mass", "charge", "ptype",
                "Van_der_Waals", "Ref_Enthalpy");

        gt_old = "";

        for (auto aType = pd.getAtypeBegin(); aType != pd.getAtypeEnd(); aType++)
        {
            gt_type = aType->getType();
            btype   = aType->getBtype();
            if (gmx_atomprop_query(aps, epropMass, "", aType->getElem().c_str(), &mass))
            {
                atomnumber = gmx_atomprop_atomnumber(aps, aType->getElem().c_str());
                if ((0 ==  gt_old.size()) || (gt_old.compare(gt_type) != 0))
                {
                    char sgt_type[32];
                    snprintf(sgt_type, 32, "%s_s", gt_type.c_str());
                    if (0 == btype.size())
                    {
                        btype = gt_type;
                    }
                    fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  A     %-s  %s\n",
                            gt_type.c_str(), aType->getBtype().c_str(), atomnumber, mass, 0.0, aType->getVdwparams().c_str(),
                            aType->getRefEnthalpy().c_str());
                    if (bPol)
                    {
                        if (strcasecmp(ff.c_str(), "LJ") == 0)
                        {
                            fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0\n",
                                    sgt_type, sgt_type, 0, 0.0, 0.0);
                        }
                        else
                        {
                            fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0  0\n",
                                    sgt_type, sgt_type, 0, 0.0, 0.0);
                        }
                    }
                }
            }
            gt_old = gt_type;
        }
        fprintf(fp, "\n");
        if (iChargeDistributionModel == eqdAXpg ||
            iChargeDistributionModel == eqdAXps)
        {
            fprintf(fp, "[ distributed_charges ]\n");
            for (auto atype = pd.getAtypeBegin(); atype != pd.getAtypeEnd(); atype++)
            {
                auto eem = pd.findEem(iChargeDistributionModel, atype->getType());
                switch(iChargeDistributionModel)
                {
                case eqdAXpg:
                    {
                        fprintf(fp, "%-5s  1  %g\n",  atype->getType().c_str(),
                                eem->getZeta(0));
                        break;
                    }
                case eqdAXps:
                    {
                        fprintf(fp, "%-5s  2  %d  %g\n",  atype->getType().c_str(), 
                                eem->getRow(0), eem->getZeta(0));
                        break;
                    }
                default:
                    GMX_RELEASE_ASSERT(false, "Death horror");
                }
            }
            fprintf(fp, "\n");
        }
    }
}

void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix)
{
    rvec au = {0, 0 ,0};
    rvec bu = {0, 0, 0};

    svmul((1.0/norm(target_vec)), target_vec, au);
    svmul((1.0/norm(ref_vec)), ref_vec, bu);

    rotmatrix[0][0] = bu[0]*au[0];
    rotmatrix[0][1] = bu[0]*au[1];
    rotmatrix[0][2] = bu[0]*au[2];
    rotmatrix[1][0] = bu[1]*au[0];
    rotmatrix[1][1] = bu[1]*au[1];
    rotmatrix[1][2] = bu[1]*au[2];
    rotmatrix[2][0] = bu[2]*au[0];
    rotmatrix[2][1] = bu[2]*au[1];
    rotmatrix[2][2] = bu[2]*au[2];     
}
    

}// namespace alexandria
