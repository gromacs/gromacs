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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "mymol.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "gauss_io.h"
#include "gentop_core.h"
#include "gentop_vsite.h"
#include "gmx_simple_comm.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "stringutil.h"

#define STRLEN 256

namespace alexandria
{

const char *immsg(immStatus imm)
{
    static const char *msg[immNR] = {
        "Unknown status",
        "OK", "Zero Dipole", "No Quadrupole", "Charged",
        "Atom type problem", "Atom number problem", "Converting from molprop",
        "Determining bond order", "RESP Initialization",
        "Charge generation", "Requested level of theory missing",
        "QM Inconsistency (ESP dipole does not match Elec)",
        "Not in training set", "No experimental data",
        "Generating shells", "Generating bonds", "Communicating MolProp"
    };

    return msg[imm];
}

static bool is_planar(rvec xi, rvec xj, rvec xk, rvec xl, t_pbc *pbc,
                      real phi_toler)
{
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real sign, phi;

    phi = RAD2DEG*dih_angle(xi, xj, xk, xl, pbc, r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);

    return (fabs(phi) < phi_toler);
}

static bool is_linear(rvec xi, rvec xj, rvec xk, t_pbc *pbc,
                      real th_toler)
{
    int  t1, t2;
    rvec r_ij, r_kj;
    real costh, th;

    th = fabs(RAD2DEG*bond_angle(xi, xj, xk, pbc, r_ij, r_kj, &costh, &t1, &t2));
    if ((th > th_toler) || (th < 180-th_toler))
    {
        if (NULL != debug)
        {
            fprintf(debug, "Angle is %g, th_toler is %g\n", th, th_toler);
        }
        return true;
    }
    return false;
}

void MyMol::getForceConstants(const Poldata &pd)
{
    int         n;
    double      xx, sx, bo;
    std::string params;

    for (std::vector<PlistWrapper>::iterator pw = plist_.begin();
         (pw < plist_.end()); ++pw)
    {
        switch (pw->getFtype())
        {
            case F_BONDS:
                for (ParamIterator j = pw->beginParam(); (j < pw->endParam()); ++j)
                {
                    std::string cai, caj;
                    if (pd.atypeToBtype( *topology_->atoms.atomtype[j->a[0]], cai) &&
                        pd.atypeToBtype( *topology_->atoms.atomtype[j->a[1]], caj))
                    {
                        int ntrain;
                        if (!pd.searchBond(cai, caj,
                                           &xx, &sx, &ntrain, &bo, params))
                        {
                            j->c[0] = convert2gmx(xx, eg2cPm);
                            std::vector<std::string> ptr = gmx::splitString(params);
                            n = 0;
                            for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                            {
                                if ((pi->length() > 0) && (n < MAXFORCEPARAM-1))
                                {
                                    j->c[1+n] = atof(pi->c_str());
                                    n++;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Default bond parameters
                        j->c[0] = 0.15;
                        j->c[1] = 2e5;
                    }
                }
                break;
            case F_ANGLES:
                for (ParamIterator j = pw->beginParam(); (j < pw->endParam()); ++j)
                {
                    std::string cai, caj, cak;
                    if (pd.atypeToBtype( *topology_->atoms.atomtype[j->a[0]], cai) &&
                        pd.atypeToBtype( *topology_->atoms.atomtype[j->a[1]], caj) &&
                        pd.atypeToBtype( *topology_->atoms.atomtype[j->a[2]], cak))
                    {
                        if (pd.searchAngle(cai, caj, cak,
                                           &xx, &sx, NULL, params))
                        {
                            j->c[0] = xx;
                            std::vector<std::string> ptr = gmx::splitString(params);
                            n = 0;
                            for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                            {
                                if ((pi->length() > 0) && (n < MAXFORCEPARAM-1))
                                {
                                    j->c[1+n] = atof(pi->c_str());
                                    n++;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Default angle parameters
                        j->c[0] = 109;
                        j->c[1] = 400;
                    }
                }
                break;
            case F_PDIHS:
            case F_IDIHS:
                int egd = egdPDIHS;
                if (pw->getFtype() == F_IDIHS)
                {
                    egd = egdIDIHS;
                }
                for (ParamIterator j = pw->beginParam(); (j < pw->endParam()); ++j)
                {
                    std::string cai, caj, cak, cal;
                    if (pd.atypeToBtype( *topology_->atoms.atomtype[j->a[0]], cai) &&
                        pd.atypeToBtype( *topology_->atoms.atomtype[j->a[1]], caj) &&
                        pd.atypeToBtype( *topology_->atoms.atomtype[j->a[2]], cak) &&
                        pd.atypeToBtype( *topology_->atoms.atomtype[j->a[3]], cal))
                    {
                        if (pd.searchDihedral(egd, cai, caj, cak, cal,
                                              &xx, &sx, NULL, params))
                        {
                            j->c[0] = xx;
                            std::vector<std::string> ptr = gmx::splitString(params);
                            n = 0;
                            for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                            {
                                if ((pi->length() > 0) && (n < MAXFORCEPARAM-1))
                                {
                                    j->c[1+n] = atof(pi->c_str());
                                    n++;
                                }
                            }
                        }
                    }
                    else if (pw->getFtype() == F_PDIHS)
                    {
                        // Default dihedral parameters
                        j->c[0] = 0;
                        j->c[1] = 5;
                        j->c[2] = 0;
                    }
                    else
                    {
                        // Default improper dihedral parameters
                        j->c[0] = 0;
                        j->c[1] = 5;
                    }
                }
                break;
                //default:
                //break;
        }
    }
}

void MyMol::MakeSpecialInteractions(bool bUseVsites)
{
    std::vector < std::vector < unsigned int> > bonds;
    std::vector<int> nbonds;
    t_pbc            pbc;
    matrix           box;
    real             th_toler = 175;
    real             ph_toler = 5;

    clear_mat(box);
    set_pbc(&pbc, epbcNONE, box);

    bonds.resize(topology_->atoms.nr);
    for (alexandria::BondIterator bi = molProp()->BeginBond(); (bi < molProp()->EndBond()); bi++)
    {
        // Store bonds bidirectionally to get the number correct
        bonds[bi->getAi() - 1].push_back(bi->getAj() - 1);
        bonds[bi->getAj() - 1].push_back(bi->getAi() - 1);
    }
    nbonds.resize(topology_->atoms.nr);
    for (int i = 0; (i < topology_->atoms.nr); i++)
    {
        nbonds[i] = bonds[i].size();
    }
    for (int i = 0; (i < topology_->atoms.nr); i++)
    {
        /* Now test initial geometry */
        if ((bonds[i].size() == 2) &&
            is_linear(x_[i], x_[bonds[i][0]], x_[bonds[i][1]],
                      &pbc, th_toler))
        {
            if (NULL != debug)
            {
                fprintf(debug, "found linear angle %s-%s-%s in %s\n",
                        *topology_->atoms.atomtype[bonds[i][0]],
                        *topology_->atoms.atomtype[i],
                        *topology_->atoms.atomtype[bonds[i][1]],
                        molProp()->getMolname().c_str());
            }
            gvt_.addLinear(bonds[i][0], i, bonds[i][1]);
        }
        else if ((bonds[i].size() == 3) &&
                 is_planar(x_[i], x_[bonds[i][0]],
                           x_[bonds[i][1]], x_[bonds[i][2]],
                           &pbc, ph_toler))
        {
            if (NULL != debug)
            {
                fprintf(debug, "found planar group %s-%s-%s-%s in %s\n",
                        *topology_->atoms.atomtype[i],
                        *topology_->atoms.atomtype[bonds[i][0]],
                        *topology_->atoms.atomtype[bonds[i][1]],
                        *topology_->atoms.atomtype[bonds[i][2]],
                        molProp()->getMolname().c_str());
            }
            gvt_.addPlanar(i, bonds[i][0], bonds[i][1], bonds[i][2],
                           &nbonds[0]);
        }
    }
    int anr = topology_->atoms.nr;

    gvt_.generateSpecial(bUseVsites, &topology_->atoms, &x_,
                         plist_, symtab_, atype_, &excls_);
    bHaveVSites_ = (topology_->atoms.nr > anr);
}

static void cp_plist(t_params plist[], int ftype,
                     std::vector<PlistWrapper> &plist_)
{
    if (plist[ftype].nr > 0)
    {
        PlistWrapper pw(ftype);
        for (int i = 0; (i < plist[ftype].nr); i++)
        {
            pw.addParam(plist[ftype].param[i]);
        }
        plist_.push_back(pw);
    }
}

void MyMol::MakeAngles(bool bPairs, bool bDihs)
{
    t_nextnb nnb;
    t_restp  rtp;
    t_params plist[F_NRE];
    std::vector<PlistWrapper>::iterator pw;

    init_plist(plist);
    for (pw = plist_.begin(); (pw < plist_.end()); ++pw)
    {
        if (F_BONDS == pw->getFtype())
        {
            pr_alloc(pw->nParam(), &plist[F_BONDS]);
            int i = 0;
            for (ParamIterator pi = pw->beginParam();
                 (pi < pw->endParam()); ++pi)
            {
                t_param *src = &(*pi);
                cp_param(&(plist[F_BONDS].param[i]), src);
                i++;
            }
            plist[F_BONDS].nr = i;
            break;
        }
    }
    /* Make Angles and Dihedrals */
    snew(excls_, topology_->atoms.nr);
    init_nnb(&nnb, topology_->atoms.nr, nexcl_+2);
    gen_nnb(&nnb, plist);

    print_nnb(&nnb, "NNB");
    rtp.bKeepAllGeneratedDihedrals    = TRUE;
    rtp.bRemoveDihedralIfWithImproper = TRUE;
    rtp.bGenerateHH14Interactions     = TRUE;
    rtp.nrexcl = nexcl_;
    gen_pad(&nnb, &(topology_->atoms), &rtp, plist, excls_, NULL, FALSE);
    {
        t_blocka *EXCL;
        snew(EXCL, 1);
        generate_excl(nexcl_, topology_->atoms.nr, plist, &nnb, EXCL);
        for (int i = 0; (i < EXCL->nr); i++)
        {
            int ne = EXCL->index[i+1]-EXCL->index[i];
            srenew(excls_[i].e, ne);
            excls_[i].nr = 0;
            for (int j = EXCL->index[i]; (j < EXCL->index[i+1]); j++)
            {
                if (EXCL->a[j] != i)
                {
                    excls_[i].e[excls_[i].nr++] = EXCL->a[j];
                }
            }
        }
        done_blocka(EXCL);
        sfree(EXCL);
        if (NULL != debug)
        {
            for (int i = 0; (i < topology_->atoms.nr); i++)
            {
                fprintf(debug, "excl %d", i);
                for (int j = 0; (j < excls_[i].nr); j++)
                {
                    fprintf(debug, "  %2d", excls_[i].e[j]);
                }
                fprintf(debug, "\n");
            }
        }
    }
    done_nnb(&nnb);

    cp_plist(plist, F_ANGLES, plist_);
    if (bDihs)
    {
        cp_plist(plist, F_PDIHS, plist_);
    }
    cp_plist(plist, F_IDIHS, plist_);
    if (bPairs)
    {
        cp_plist(plist, F_LJ14, plist_);
    }
    for (int i = 0; (i < F_NRE); i++)
    {
        if (plist[i].nr > 0)
        {
            sfree(plist[i].param);
        }
    }
}

static void generate_nbparam(int ftype, int comb, double ci[], double cj[],
                             t_iparams *ip)
{
    double sig, eps;

    switch (ftype)
    {
        case F_LJ:
            switch (comb)
            {
                case eCOMB_GEOMETRIC:
                    /* Gromos rules */
                    ip->lj.c6  = sqrt(ci[0] * cj[0]);
                    ip->lj.c12 = sqrt(ci[1] * cj[1]);
                    break;

                case eCOMB_ARITHMETIC:
                    /* c0 and c1 are epsilon and sigma */
                    sig        = (ci[0]+cj[0])*0.5;
                    eps        = sqrt(ci[1]*cj[1]);
                    ip->lj.c6  = 4*eps*pow(sig, 6);
                    ip->lj.c12 = 4*eps*pow(sig, 12);

                    break;
                case eCOMB_GEOM_SIG_EPS:
                    /* c0 and c1 are epsilon and sigma */
                    sig        = sqrt(ci[0]*cj[0]);
                    eps        = sqrt(ci[1]*cj[1]);
                    ip->lj.c6  = 4*eps*pow(sig, 6);
                    ip->lj.c12 = 4*eps*pow(sig, 12);

                    break;
                default:
                    gmx_fatal(FARGS, "No such combination rule %d", comb);
            }
            break;
        default:
            gmx_fatal(FARGS, "No such function type supported %s",
                      interaction_function[ftype].name);
    }
}

static void do_init_mtop(const Poldata &pd,
                         gmx_mtop_t    *mtop_,
                         char         **molname,
                         t_atoms       *atoms)
{
    init_mtop(mtop_);
    mtop_->name     = molname;
    mtop_->nmoltype = 1;
    snew(mtop_->moltype, mtop_->nmoltype);
    mtop_->moltype[0].name = molname;
    mtop_->nmolblock       = 1;
    snew(mtop_->molblock, mtop_->nmolblock);
    mtop_->molblock[0].nmol        = 1;
    mtop_->molblock[0].type        = 0;
    mtop_->molblock[0].natoms_mol  = atoms->nr;
    mtop_->groups.grps[egcENER].nr = 1;

    //! Count the number of types in this molecule, at least 1 assuming there is one atom
    int ntype = 1;
    for (int i = 1; (i < atoms->nr); i++)
    {
        int  itp   = atoms->atom[i].type;
        bool found = false;
        for (int j = 0; !found && (j < i); j++)
        {
            found = (itp == atoms->atom[j].type);
        }
        if (!found)
        {
            ntype++;
        }
    }

    mtop_->ffparams.atnr   = ntype;
    mtop_->ffparams.ntypes = ntype*ntype;
    mtop_->ffparams.reppow = 12;

    int vdw_type = pd.getVdwFtype();

    snew(mtop_->ffparams.functype, mtop_->ffparams.ntypes);
    snew(mtop_->ffparams.iparams, mtop_->ffparams.ntypes);
    for (int i = 0; (i < ntype); i++)
    {
        for (int j = 0; (j < ntype); j++)
        {
            int idx = ntype*i+j;
            mtop_->ffparams.functype[idx] = vdw_type;
            switch (vdw_type)
            {
                case F_LJ:
                    //! NOTE  get the real parameters from the pd here
                    //! May need to set the atomtypes properly too.
                    mtop_->ffparams.iparams[idx].lj.c6  = 0;
                    mtop_->ffparams.iparams[idx].lj.c12 = 0;
                    break;
                case F_BHAM:
                    mtop_->ffparams.iparams[idx].bham.a = 0;
                    mtop_->ffparams.iparams[idx].bham.b = 0;
                    mtop_->ffparams.iparams[idx].bham.c = 0;
                    break;
                default:
                    fprintf(stderr, "Invalid van der waals type %s\n",
                            pd.getVdwFunction().c_str());
            }
        }
    }

    /* Create a charge group block */
    stupid_fill_block(&(mtop_->moltype[0].cgs), atoms->nr, FALSE);

    mtop_->natoms = atoms->nr;
    init_t_atoms(&(mtop_->moltype[0].atoms), atoms->nr, FALSE);
}

static void excls_to_blocka(int natom, t_excls excls_[], t_blocka *blocka)
{
    int i, j, k, nra;

    if (blocka->nr < natom)
    {
        srenew(blocka->index, natom+1);
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

static void plist_to_mtop(const Poldata             &pd,
                          std::vector<PlistWrapper>  plist,
                          gmx_mtop_t                *mtop_)
{
    double fudgeLJ;
    double reppow = 12.0;
    int    n      = 0;

    /* Generate pairs */
    fudgeLJ = pd.getFudgeLJ();

    int nfptot = mtop_->ffparams.ntypes;
    for (std::vector<PlistWrapper>::iterator pw = plist.begin();
         (pw < plist.end()); ++pw)
    {
        nfptot += pw->nParam()*NRFPA(pw->getFtype());
    }
    srenew(mtop_->ffparams.functype, nfptot);
    srenew(mtop_->ffparams.iparams, nfptot);

    for (std::vector<PlistWrapper>::iterator pw = plist.begin();
         (pw < plist.end()); ++pw)
    {
        int nra    = NRAL(pw->getFtype());
        int nrfp   = NRFPA(pw->getFtype());
        int nratot = pw->nParam()*(1+nra);
        snew(mtop_->moltype[0].ilist[pw->getFtype()].iatoms, nratot);
        int k = 0;
        for (ParamIterator j = pw->beginParam();
             (j < pw->endParam()); ++j)
        {
            real c[MAXFORCEPARAM];
            int  l = 0;
            if (pw->getFtype() == F_LJ14)
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
            n = enter_params(&mtop_->ffparams, pw->getFtype(), c, 0, reppow, n, TRUE);
            mtop_->moltype[0].ilist[pw->getFtype()].iatoms[k++] = n;
            for (l = 0; (l < nra); l++)
            {
                mtop_->moltype[0].ilist[pw->getFtype()].iatoms[k++] = j->a[l];
            }
        }
        mtop_->moltype[0].ilist[pw->getFtype()].nr = k;
    }
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

bool MyMol::IsSymmetric(real toler)
{
    int       i, j, m;
    real      mm, tm;
    rvec      com, test;
    gmx_bool *bSymm, bSymmAll;

    clear_rvec(com);
    tm = 0;
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        mm  = topology_->atoms.atom[i].m;
        tm += mm;
        for (m = 0; (m < DIM); m++)
        {
            com[m] += mm*x_[i][m];
        }
    }
    if (tm > 0)
    {
        for (m = 0; (m < DIM); m++)
        {
            com[m] /= tm;
        }
    }
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        rvec_dec(x_[i], com);
    }

    snew(bSymm, topology_->atoms.nr);
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        bSymm[i] = (norm(x_[i]) < toler);
        for (j = i+1; (j < topology_->atoms.nr) && !bSymm[i]; j++)
        {
            rvec_add(x_[i], x_[j], test);
            if (norm(test) < toler)
            {
                bSymm[i] = TRUE;
                bSymm[j] = TRUE;
            }
        }
    }
    bSymmAll = TRUE;
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        rvec_inc(x_[i], com);
    }

    return bSymmAll;
}

static void fill_inputrec(t_inputrec *ir)
{
    ir->cutoff_scheme = ecutsGROUP;
    ir->tabext        = 2; /* nm */
    ir->ePBC          = epbcNONE;
    ir->epsilon_r     = 1;
    ir->vdwtype       = evdwCUT;
    ir->coulombtype   = eelCUT;
    ir->eDispCorr     = edispcNO;
    snew(ir->opts.egp_flags, 1);
    snew(ir->fepvals, 1);
}

MyMol::MyMol() : gvt_(egvtALL)
{
    bHaveShells_       = false;
    bHaveVSites_       = false;
    cgnr_              = NULL;
    gr_                = NULL;
    immAtoms_          = immOK;
    immTopology_       = immOK;
    immCharges_        = immOK;
    shellfc_           = nullptr;
    snew(symtab_, 1);
    open_symtab(symtab_);
    atype_ = init_atomtype();
    clear_mat(box);
    mtop_  = NULL;
    ltop_  = NULL;
    md_    = NULL;
    mp_    = new MolProp;
    snew(state_, 1);
    snew(enerd_, 1);
    init_enerdata(1, 0, enerd_);

    /* Inputrec parameters */
    snew(inputrec_, 1);
    fill_inputrec(inputrec_);
}

MyMol::~MyMol()
{
    return;
    if (NULL != cgnr_)
    {
        sfree(cgnr_);
        cgnr_ = NULL;
    }
    if (NULL != atype_)
    {
        done_atomtype(atype_);
        atype_ = NULL;
    }
    if (NULL != inputrec_)
    {
        sfree(inputrec_);
        inputrec_ = NULL;
    }
    for (std::vector<PlistWrapper>::iterator pw = plist_.begin();
         (pw < plist_.end()); ++pw)
    {
        pw->eraseParams();
    }
    if (NULL != symtab_)
    {
        done_symtab(symtab_);
        sfree(symtab_);
        symtab_ = NULL;
    }
}

immStatus MyMol::GenerateAtoms(gmx_atomprop_t            ap,
                               const char               *lot,
                               ChargeDistributionModel   iChargeDistributionModel)
{
    int                 myunit;
    double              xx, yy, zz;
    int                 natom;
    immStatus           imm   = immOK;

    ExperimentIterator  ci = molProp()->getLot(lot);
    if (ci < molProp()->EndExperiment())
    {
        t_param nb;

        memset(&nb, 0, sizeof(nb));
        natom = 0;
        init_t_atoms(&(topology_->atoms), ci->NAtom(), FALSE);
        snew(x_, ci->NAtom());
        snew(topology_->atoms.atomtype, ci->NAtom());
        snew(topology_->atoms.atomtypeB, ci->NAtom());

        for (CalcAtomIterator cai = ci->BeginAtom(); (cai < ci->EndAtom()); cai++)
        {
            myunit = string2unit((char *)cai->getUnit().c_str());
            if (myunit == -1)
            {
                gmx_fatal(FARGS, "Unknown unit '%s' for atom coords",
                          cai->getUnit().c_str());
            }
            cai->getCoords(&xx, &yy, &zz);
            x_[natom][XX] = convert2gmx(xx, myunit);
            x_[natom][YY] = convert2gmx(yy, myunit);
            x_[natom][ZZ] = convert2gmx(zz, myunit);

            double q = 0;
            for (AtomicChargeIterator qi = cai->BeginQ(); (qi < cai->EndQ()); qi++)
            {
                // TODO Clean up this mess.
                if ((qi->getType().compare("ESP") == 0) ||
                    (name2eemtype(qi->getType()) == iChargeDistributionModel))
                {
                    myunit = string2unit((char *)qi->getUnit().c_str());
                    q      = convert2gmx(qi->getQ(), myunit);
                    break;
                }
            }
            topology_->atoms.atom[natom].q      =
                topology_->atoms.atom[natom].qB = q;

            t_atoms_set_resinfo(&(topology_->atoms), natom, symtab_, molProp()->getMolname().c_str(), 1, ' ', 1, ' ');
            topology_->atoms.atomname[natom]        = put_symtab(symtab_, cai->getName().c_str());
            topology_->atoms.atom[natom].atomnumber = gmx_atomprop_atomnumber(ap, cai->getName().c_str());

            real mass = 0;
            if (!gmx_atomprop_query(ap, epropMass, "???", cai->getName().c_str(), &mass))
            {
                fprintf(stderr, "Could not find mass for %s\n", cai->getName().c_str());
            }
            topology_->atoms.atom[natom].m      =
                topology_->atoms.atom[natom].mB = mass;

            strcpy(topology_->atoms.atom[natom].elem, gmx_atomprop_element(ap, topology_->atoms.atom[natom].atomnumber));

            topology_->atoms.atom[natom].resind = 0;
            // First set the atomtype
            topology_->atoms.atomtype[natom]      =
                topology_->atoms.atomtypeB[natom] = put_symtab(symtab_, cai->getObtype().c_str());

            natom++;
        }
        for (int i = 0; (i < natom); i++)
        {
            topology_->atoms.atom[i].type      =
                topology_->atoms.atom[i].typeB = add_atomtype(atype_, symtab_,
                                                              &(topology_->atoms.atom[i]),
                                                              *topology_->atoms.atomtype[i],
                                                              &nb,
                                                              0, 0.0, 0.0, 0.0,
                                                              topology_->atoms.atom[i].atomnumber,
                                                              0.0, 0.0);
        }
        topology_->atoms.nr   = natom;
        topology_->atoms.nres = 1;
    }
    else
    {
        imm = immLOT;
    }
    if (NULL != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molProp()->getMolname().c_str(), lot, natom);
    }

    return imm;
}

immStatus MyMol::GenerateTopology(gmx_atomprop_t          ap,
                                  const Poldata          &pd,
                                  const char             *lot,
                                  ChargeDistributionModel iChargeDistributionModel,
                                  int                     nexcl,
                                  bool                    bUseVsites,
                                  bool                    bPairs,
                                  bool                    bDih)
{
    immStatus                imm = immOK;
    int                      ftb;
    t_param                  b;

    if (NULL != debug)
    {
        fprintf(debug, "Generating topology_ for %s\n", molProp()->getMolname().c_str());
    }

    nexcl_ = nexcl;
    molProp()->GenerateComposition(pd);
    if (molProp()->NAtom() <= 0)
    {
        imm = immAtomTypes;
    }
    if (immOK == imm)
    {
        snew(topology_, 1);
        init_top(topology_);
        /* get atoms */
        imm = GenerateAtoms(ap, lot, iChargeDistributionModel);
    }
    /* Store bonds in harmonic potential list first, update type later */
    ftb = F_BONDS;
    if (immOK == imm)
    {
        memset(&b, 0, sizeof(b));
        for (alexandria::BondIterator bi = molProp()->BeginBond(); (bi < molProp()->EndBond()); bi++)
        {
            b.a[0] = bi->getAi() - 1;
            b.a[1] = bi->getAj() - 1;
            add_param_to_plist(plist_, ftb, b);
        }
        if (molProp()->NBond() == 0)
        {
            imm = immGenBonds;
        }
    }
    if (immOK == imm)
    {
        /* Make Angles and Dihedrals. This needs the bonds to be F_BONDS. */
        MakeAngles(bPairs, bDih);

        /* Linear angles and or vsites etc. */
        MakeSpecialInteractions(bUseVsites);

        /* Move the plist_ to the correct function */
        //mv_plists(pd, plist_, true);

        getForceConstants(pd);

        char **molnameptr = put_symtab(symtab_, molProp()->getMolname().c_str());
        snew(mtop_, 1);
        do_init_mtop(pd, mtop_, molnameptr, &topology_->atoms);

        plist_to_mtop(pd, plist_, mtop_);
        excls_to_blocka(topology_->atoms.nr, excls_,
                        &(mtop_->moltype[0].excls));

        ltop_ = gmx_mtop_generate_local_top(mtop_, inputrec_);
    }

    return imm;
}

void MyMol::CalcMultipoles()
{
    int                     i, m;
    rvec                    mu, mm;
    real                    r2, dfac, q;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    int                     at_global;

    clear_rvec(mu);
    aloop = gmx_mtop_atomloop_all_init(mtop_);
    i     = 0;
    clear_mat(Q_calc);
    clear_rvec(coq);
    while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
    {
        q = atom->q;
        svmul(ENM2DEBYE*q, x_[i], mm);
        rvec_inc(mu, mm);

        dfac = q*0.5*10*ENM2DEBYE;
        r2   = iprod(x_[i], x_[i]);
        for (m = 0; (m < DIM); m++)
        {
            Q_calc[m][m] += dfac*(3*gmx::square(x_[i][m]) - r2);
        }
        Q_calc[XX][YY] += dfac*3*(x_[i][XX]+coq[XX])*(x_[i][YY]+coq[YY]);
        Q_calc[XX][ZZ] += dfac*3*(x_[i][XX]+coq[XX])*(x_[i][ZZ]+coq[ZZ]);
        Q_calc[YY][ZZ] += dfac*3*(x_[i][YY]+coq[YY])*(x_[i][ZZ]+coq[ZZ]);

        i++;
    }
    GMX_RELEASE_ASSERT(i == topology_->atoms.nr, "Inconsistency 1 in mymol.cpp");
    copy_rvec(mu, mu_calc);
    dip_calc = norm(mu);
}

immStatus MyMol::GenerateCharges(const Poldata &pd,
                                 gmx_atomprop_t ap,
                                 ChargeDistributionModel iChargeDistributionModel,
                                 ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                                 real hfac, real epsr,
                                 const char *lot,
                                 bool bSymmetricCharges,
                                 const char *symm_string)
{
    int       i, eQGEN;
    char      qgen_msg[STRLEN];
    immStatus imm = immOK;

    qgen_ = new GentopQgen(pd, &topology_->atoms, ap, x_,
                           iChargeDistributionModel,
                           iChargeGenerationAlgorithm,
                           hfac, molProp()->getCharge(), epsr);
    if (NULL == qgen_)
    {
        return immChargeGeneration;
    }
    if (iChargeGenerationAlgorithm == eqgNONE)
    {
        return imm;
    }
    if (immOK == imm)
    {
        if (bSymmetricCharges)
        {
            std::vector<PlistWrapper>::iterator pw = SearchPlist(plist_, F_BONDS);
            if (plist_.end() != pw)
            {
                symmetrize_charges(bSymmetricCharges,
                                   &topology_->atoms,
                                   pw,
                                   pd, ap, symm_string, symmetric_charges_);
            }
        }
    }

    if (immOK == imm)
    {
        eQGEN = eQGEN_OK;
        switch (iChargeGenerationAlgorithm)
        {
            case eqgRESP:
                if (gr_->addAtomInfo(&topology_->atoms, pd))
                {
                    gr_->addAtomSymmetry(symmetric_charges_);
                    gr_->updateAtomtypes(&(topology_->atoms));
                    gr_->summary(debug, symmetric_charges_);
                    gr_->addAtomCoords(x_);
                    /* Even if we get the right LoT it may still not have
                     * the ESP
                     */
                    ExperimentIterator ci = molProp()->getLotPropType(lot,
                                                                      MPO_POTENTIAL,
                                                                      NULL);
                    if (ci != molProp()->EndExperiment())
                    {
                        //printf("There are %d potential points\n",ci->NPotential());
                        for (ElectrostaticPotentialIterator epi = ci->BeginPotential(); (epi < ci->EndPotential()); ++epi)
                        {
                            /* Maybe not convert to gmx ? */
                            int xu = string2unit(epi->getXYZunit().c_str());
                            int vu = string2unit(epi->getVunit().c_str());
                            if (-1 == xu)
                            {
                                xu = eg2cAngstrom;
                            }
                            if (-1 == vu)
                            {
                                vu = eg2cHartree_e;
                            }
                            gr_->addPoint(convert2gmx(epi->getX(), xu),
                                          convert2gmx(epi->getY(), xu),
                                          convert2gmx(epi->getZ(), xu),
                                          convert2gmx(epi->getV(), vu));
                        }
                    }
                    eQGEN = qgen_->generateCharges(NULL,
                                                   gr_, molProp()->getMolname().c_str(),
                                                   pd, &topology_->atoms, 0.0001,
                                                   10000, 1, ap);
                }
                break;
            case eqgESP:
                break;
            case eqgNONE:
                /* Check which algorithm to use for charge generation */
                strcpy(qgen_msg, "");
                printf("Using zero charges!\n");
                for (i = 0; (i < topology_->atoms.nr); i++)
                {
                    topology_->atoms.atom[i].q  = topology_->atoms.atom[i].qB = 0;
                }
                eQGEN = eQGEN_OK;
                break;
            default:
                if (NULL == qgen_)
                {
                    gmx_fatal(FARGS, "Can not generate charges for %s. Probably due to issues with atomtype detection or support.\n", molProp()->getMolname().c_str());
                }
                eQGEN = qgen_->generateCharges(NULL,
                                               NULL, molProp()->getMolname().c_str(),
                                               pd, &topology_->atoms, 0.0001,
                                               10000, 1, ap);
                break;
        }
        qgen_->message(sizeof(qgen_msg), qgen_msg, gr_);
        if (eQGEN_OK != eQGEN)
        {
            imm = immChargeGeneration;
        }
    }
    return imm;
}

immStatus MyMol::GenerateGromacs(t_commrec *cr)
{
    int nalloc = 2 * topology_->atoms.nr;

    snew(f_, nalloc);
    fr_ = mk_forcerec();
    init_forcerec(NULL, fr_, NULL, inputrec_, mtop_, cr,
                  box, NULL, NULL, NULL, NULL, TRUE, -1);

    init_state(state_, topology_->atoms.nr, 1, 1, 1, 0);
    ltop_ = gmx_mtop_generate_local_top(mtop_, inputrec_);
    md_   = init_mdatoms(NULL, mtop_, FALSE);
    for (int i = 0; (i < topology_->atoms.nr); i++)
    {
        copy_rvec(x_[i], state_->x[i]);
    }
    return immOK;
}

static void put_in_box(int natom, matrix box, rvec x[], real dbox)
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

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];

    put_in_box(topology_->atoms.nr, box, x_, 0.3);
    sprintf(title, "%s processed by alexandria", molProp()->getMolname().c_str());
    write_sto_conf(fn, title, &topology_->atoms, x_, NULL, epbcNONE, box);
}

static void write_zeta_q(FILE *fp, GentopQgen * qgen,
                         t_atoms *atoms, ChargeDistributionModel iChargeDistributionModel)
{
    int    i, ii, j, k, nz, row;
    double zeta, q;
    bool   bAtom, bTypeSet;

    if (NULL == qgen)
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

static void write_zeta_q2(GentopQgen * qgen, gpp_atomtype_t atype,
                          t_atoms *atoms, ChargeDistributionModel iChargeDistributionModel)
{
    FILE      *fp;
    int        i, j, k, nz, row;
    double     zeta, q, qtot;
    gmx_bool   bAtom;

    if (NULL == qgen)
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

static int get_subtype(directive d, int ftype)
{
    int i;
    for (i = 1; (i < 20); i++)
    {
        if (ifunc_index(d, i) == ftype)
        {
            return i;
        }
    }
    return 1;
}

static void print_bondeds2(FILE                     *out,
                           directive                 d,
                           int                       plist_ftype,
                           int                       print_ftype,
                           std::vector<PlistWrapper> plist)
{
    std::vector<PlistWrapper>::iterator p = SearchPlist(plist, plist_ftype);

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
    for (ParamIterator i = p->beginParam(); (i < p->endParam()); ++i)
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

static void write_top2(FILE *out, char *molname,
                       t_atoms *at, gmx_bool bRTPresname,
                       std::vector<PlistWrapper> plist_,
                       t_excls excls[],
                       gpp_atomtype_t atype, int *cgnr, int nrexcl,
                       const Poldata &pd)
/* NOTE: nrexcl is not the size of *excl! */
{
    if (at && atype && cgnr)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nrexcl);

        print_atoms(out, atype, at, cgnr, bRTPresname);
        print_bondeds2(out, d_bonds, F_BONDS,
                       pd.getBondFtype(),
                       plist_);
        print_bondeds2(out, d_constraints, F_CONSTR, F_CONSTR, plist_);
        print_bondeds2(out, d_constraints, F_CONSTRNC, F_CONSTRNC, plist_);
        print_bondeds2(out, d_pairs, F_LJ14, F_LJ14, plist_);
        print_excl(out, at->nr, excls);
        print_bondeds2(out, d_angles, F_ANGLES,
                       pd.getAngleFtype(), plist_);
        print_bondeds2(out, d_angles, F_LINEAR_ANGLES, F_LINEAR_ANGLES,
                       plist_);
        print_bondeds2(out, d_dihedrals, F_PDIHS,
                       pd.getDihedralFtype( egdPDIHS), plist_);
        /* Check whether the dihedrals use the same function */
        print_bondeds2(out, d_dihedrals, F_IDIHS,
                       pd.getDihedralFtype( egdIDIHS), plist_);
        print_bondeds2(out, d_cmap, F_CMAP, F_CMAP, plist_);
        print_bondeds2(out, d_polarization, F_POLARIZATION, F_POLARIZATION,
                       plist_);
        print_bondeds2(out, d_thole_polarization, F_THOLE_POL, F_THOLE_POL, plist_);
        print_bondeds2(out, d_vsites2, F_VSITE2, F_VSITE2, plist_);
        print_bondeds2(out, d_vsites3, F_VSITE3, F_VSITE3, plist_);
        print_bondeds2(out, d_vsites3, F_VSITE3FD, F_VSITE3FD, plist_);
        print_bondeds2(out, d_vsites3, F_VSITE3FAD, F_VSITE3FAD, plist_);
        print_bondeds2(out, d_vsites3, F_VSITE3OUT, F_VSITE3OUT, plist_);
        print_bondeds2(out, d_vsites4, F_VSITE4FD, F_VSITE4FD, plist_);
        print_bondeds2(out, d_vsites4, F_VSITE4FDN, F_VSITE4FDN, plist_);
    }
}


static void print_top_header2(FILE *fp, const Poldata &pd,
                              gmx_atomprop_t aps, bool bPol,
                              std::vector<std::string> commercials,
                              bool bItp)
{
    std::string   gt_old, gt_type;
    std::string   btype;
    int           atomnumber;
    real          mass;

    fprintf(fp, ";\n");
    fprintf(fp, "; Topology generated by alexandria gentop.\n");
    fprintf(fp, "; Watch this space for information & commercials.\n");
    for (std::vector<std::string>::iterator i = commercials.begin(); (i < commercials.end()); ++i)
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

        for (FfatypeConstIterator aType = pd.getAtypeBegin();
             aType != pd.getAtypeEnd(); aType++)
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
                    fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  A     %-s  %f\n",
                            gt_type.c_str(), aType->getBtype().c_str(), atomnumber, mass, 0.0, aType->getVdwparams().c_str(),
                            aType->getRefEnthalpy());
                    if (bPol)
                    {
                        fprintf(fp, "%-6s %-6s %6d  %12.6f  %10.4f  S     0  0\n",
                                sgt_type, sgt_type, 0, 0.0, 0.0);
                    }
                }
            }
            gt_old = gt_type;
        }
        fprintf(fp, "\n");
    }
}

void MyMol::PrintTopology(const char             *fn,
                          ChargeDistributionModel iChargeDistributionModel,
                          bool                    bVerbose,
                          const Poldata          &pd,
                          gmx_atomprop_t          aps)
{
    FILE                    *fp;
    t_mols                   printmol;
    bool                     bITP;
    std::vector<std::string> commercials;
    char                     buf[256];

    CalcQPol(pd);

    if (molProp()->getMolname().size() > 0)
    {
        printmol.name = strdup(molProp()->getMolname().c_str());
    }
    else if (molProp()->formula().size() > 0)
    {
        printmol.name = strdup(molProp()->formula().c_str());
    }
    else
    {
        printmol.name = strdup("Onbekend");
    }
    printmol.nr   = 1;

    /* Write topology_ file */
    bITP = (fn2ftp(fn) == efITP);
    fp   = gmx_ffopen(fn, "w");

    snprintf(buf, sizeof(buf), "ref_enthalpy   = %.3f kJ/mol", ref_enthalpy_);
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "polarizability = %.3f +/- %.3f A^3",
             polarizability_, sig_pol_);
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "total charge   = %d e", molProp()->getCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "total mass     = %.3f Da", molProp()->getMass());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "total dipole   = %.3f D", mutot_);
    commercials.push_back(buf);
    print_top_header2(fp, pd, aps, bHaveShells_, commercials, bITP);

    if (bHaveShells_ || (iChargeDistributionModel == eqdAXg) || (iChargeDistributionModel == eqdAXs))
    {
        //write_zeta_q(fp, qgen_, &topology_->atoms, iChargeDistributionModel);
        //write_zeta_q2(qgen,atype,&topology_->atoms,pd,iChargeDistributionModel);
    }
    write_top2(fp, printmol.name, &topology_->atoms, FALSE,
               plist_, excls_, atype_, cgnr_, nexcl_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, getForceField().c_str(), NULL, 0, NULL, 1, &printmol);
    }

    if (bVerbose)
    {
        printf("There are %4d proper dihedrals, %4d impropers\n"
               "          %4d angles, %4d linear angles\n"
               "          %4d pairs, %4d bonds, %4d atoms\n"
               "          %4d polarizations\n",
               CountPlist(plist_, F_PDIHS),
               CountPlist(plist_, F_IDIHS),
               CountPlist(plist_, F_ANGLES),
               CountPlist(plist_, F_LINEAR_ANGLES),
               CountPlist(plist_, F_LJ14),
               CountPlist(plist_, F_BONDS),
               topology_->atoms.nr,
               CountPlist(plist_, F_POLARIZATION));
        for (std::vector<std::string>::iterator i = commercials.begin(); (i < commercials.end()); ++i)
        {
            printf("%s\n", i->c_str());
        }
    }

    fclose(fp);

}

static void add_excl(t_excls *excls, int e)
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

static void add_excl_pair(t_excls excls[], int e1, int e2)
{
    if (e1 != e2)
    {
        add_excl(&excls[e1], e2);
        add_excl(&excls[e2], e1);
    }
}

static void remove_excl(t_excls *excls, int remove)
{
    int i;

    for (i = remove+1; i < excls->nr; i++)
    {
        excls->e[i-1] = excls->e[i];
    }

    excls->nr--;
}

static void prune_excl(t_excls excls[], t_atoms *atoms, gpp_atomtype_t atype)
{
    int i, k, ak;

    for (i = 0; (i < atoms->nr); i++)
    {
        if (get_atomtype_ptype(atoms->atom[i].type, atype) != eptShell)
        {
            for (k = 0; (k < excls[i].nr); )
            {
                ak = excls[i].e[k];
                if (get_atomtype_ptype(atoms->atom[ak].type, atype) != eptShell)
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

static void copy_atoms(t_atoms *src, t_atoms *dest)
{
    int i;

    if (dest->nr < src->nr)
    {
        srenew(dest->atom, src->nr);
        srenew(dest->atomname, src->nr);
        if (NULL != src->atomtype)
        {
            srenew(dest->atomtype, src->nr);
        }
        else if (NULL != dest->atomtype)
        {
            sfree(dest->atomtype);
            dest->atomtype = NULL;
        }
        if (NULL != src->atomtypeB)
        {
            srenew(dest->atomtypeB, src->nr);
        }
        else if (NULL != dest->atomtypeB)
        {
            sfree(dest->atomtypeB);
            dest->atomtypeB = NULL;
        }
    }
    dest->nr = src->nr;
    for (i = 0; (i < src->nr); i++)
    {
        dest->atom[i]      = src->atom[i];
        dest->atomname[i]  = src->atomname[i];
        if (NULL != src->atomtype)
        {
            dest->atomtype[i]  = src->atomtype[i];
        }
        if (NULL != src->atomtypeB)
        {
            dest->atomtypeB[i] = src->atomtypeB[i];
        }
    }
    if (dest->nres < src->nres)
    {
        srenew(dest->resinfo, src->nres);
    }

    if (NULL != src->pdbinfo)
    {
        srenew(dest->pdbinfo, src->nres);
    }
    else if (NULL != dest->pdbinfo)
    {
        sfree(dest->pdbinfo);
        dest->pdbinfo = NULL;
    }
    dest->nres = src->nres;
    for (i = 0; (i < src->nres); i++)
    {
        dest->resinfo[i] = src->resinfo[i];
        if (NULL != src->pdbinfo)
        {
            dest->pdbinfo[i] = src->pdbinfo[i];
        }
    }
}

void MyMol::AddShells(const Poldata &pd, ePolar epol, ChargeDistributionModel iModel)
{
    int              i, j, k, iat, shell, ns = 0;
    std::vector<int> renum, inv_renum;
    char             buf[32], **newname;
    t_param          p;
    t_atom          *shell_atom;
    t_atoms         *newa;
    t_excls         *newexcls;
    rvec            *newx;
    double           pol, sigpol;

    if (epol == epolNo)
    {
        return;
    }
    int      maxatom = topology_->atoms.nr*2+2;
    srenew(x_, maxatom);
    snew(shell_atom, 1);
    shell_atom->ptype = eptShell;
    memset(&p, 0, sizeof(p));
    inv_renum.reserve(topology_->atoms.nr*2);
    for (i = 0; (i < topology_->atoms.nr*2); i++)
    {
        inv_renum[i] = -1;
    }
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        if (pd.getAtypePol(*topology_->atoms.atomtype[i],
                           &pol, &sigpol))
        {
            renum.push_back(i+ns);
            inv_renum[i+ns] = i;
            ns++;
            p.a[0] = renum[i];
            p.a[1] = renum[i]+1;
            p.c[0] = 0.001*pol;
            add_param_to_plist(plist_, F_POLARIZATION, p);
        }
    }
    renum[topology_->atoms.nr] = topology_->atoms.nr + ns;
    printf("added %d shells\n", ns);
    if (ns > 0)
    {
        /* Make new atoms and x arrays */
        snew(newa, 1);
        init_t_atoms(newa, topology_->atoms.nr+ns, TRUE);
        snew(newa->atomtype, topology_->atoms.nr+ns);
        snew(newa->atomtypeB, topology_->atoms.nr+ns);
        newa->nres = topology_->atoms.nres;
        snew(newx, newa->nr);
        snew(newname, newa->nr);

        /* Make new exclusion array, and put the shells in it */
        snew(newexcls, newa->nr);
        /* TODO: other polarization types */
        std::vector<PlistWrapper>::iterator pw = SearchPlist(plist_, F_POLARIZATION);
        if (plist_.end() != pw)
        {
            for (ParamIterator j = pw->beginParam();
                 (j < pw->endParam()); ++j)
            {
                // Exclude nucleus and shell from each other
                add_excl_pair(newexcls, j->a[0], j->a[1]);
            }
            for (ParamIterator j = pw->beginParam();
                 (j < pw->endParam()); ++j)
            {
                // Now add the exclusions from the nucleus to the shell.
                // We know that the nuclues is 0 since we just made the list
                int i0 = inv_renum[j->a[0]];
                for (int j0 = 0; (j0 < excls_[i0].nr); j0++)
                {
                    add_excl_pair(newexcls, j->a[0], renum[excls_[i0].e[j0]]);
                    add_excl_pair(newexcls, j->a[1], renum[excls_[i0].e[j0]]);
                }
            }
            for (ParamIterator j = pw->beginParam();
                 (j < pw->endParam()); ++j)
            {
                for (int j0 = 0; (j0 < newexcls[j->a[0]].nr); j0++)
                {
                    add_excl_pair(newexcls, j->a[1], newexcls[j->a[0]].e[j0]);
                }
            }
        }
        // Now copy the old atoms to the new structures
        for (i = 0; (i < topology_->atoms.nr); i++)
        {
            newa->atom[renum[i]]      = topology_->atoms.atom[i];
            newa->atomname[renum[i]]  = put_symtab(symtab_, *topology_->atoms.atomname[i]);
            newa->atomtype[renum[i]]  = put_symtab(symtab_, *topology_->atoms.atomtype[i]);
            newa->atomtypeB[renum[i]] = put_symtab(symtab_, *topology_->atoms.atomtypeB[i]);
            copy_rvec(x_[i], newx[renum[i]]);
            newname[renum[i]] = *topology_->atoms.atomtype[i];
            t_atoms_set_resinfo(newa, renum[i], symtab_,
                                *topology_->atoms.resinfo[topology_->atoms.atom[i].resind].name,
                                topology_->atoms.atom[i].resind, ' ', 1, ' ');
        }
        // Now insert the shell particles
        for (i = 0; (i < topology_->atoms.nr); i++)
        {
            iat = renum[i];
            for (j = iat+1; (j < renum[i+1]); j++)
            {
                newa->atom[j]            = topology_->atoms.atom[i];
                newa->atom[iat].q        = pd.getQ(iModel, *topology_->atoms.atomtype[i], 0);
                newa->atom[iat].qB       = pd.getQ(iModel, *topology_->atoms.atomtype[i], 0);
                newa->atom[j].m          = 0;
                newa->atom[j].mB         = 0;
                newa->atom[j].atomnumber = 0;
                sprintf(buf, "%s_s", get_atomtype_name(topology_->atoms.atom[i].type,
                                                       atype_));
                newname[j] = strdup(buf);
                shell      = add_atomtype(atype_, symtab_, shell_atom, buf, &p,
                                          0, 0, 0, 0, 0, 0, 0);
                newa->atom[j].type          = shell;
                newa->atom[j].typeB         = shell;
                newa->atomtype[j]           =
                    newa->atomtypeB[j]      = put_symtab(symtab_, buf);
                newa->atom[j].ptype         = eptShell;
                newa->atom[j].q             = newa->atom[j].q - newa->atom[iat].q;
                newa->atom[j].qB            = 0;
                newa->atom[j].resind        = topology_->atoms.atom[i].resind;
                sprintf(buf, "%ss", *(topology_->atoms.atomname[i]));
                newa->atomname[j] = put_symtab(symtab_, buf);
                copy_rvec(x_[i], newx[j]);
            }
        }
        /* Copy newa to atoms */
        copy_atoms(newa, &topology_->atoms);
        /* Copy coordinates and smnames */
        for (i = 0; (i < newa->nr); i++)
        {
            copy_rvec(newx[i], x_[i]);
            topology_->atoms.atomtype[i] = put_symtab(symtab_, newname[i]);
        }
        sfree(newx);
        sfree(newname);
        /* Copy exclusions, may need to empty the original first */
        sfree(excls_);
        excls_ = newexcls;

        for (PlistWrapperIterator i = plist_.begin();
             (i < plist_.end()); ++i)
        {
            if (i->getFtype() != F_POLARIZATION)
            {
                for (ParamIterator j = i->beginParam();
                     (j < i->endParam()); ++j)
                {
                    for (k = 0; (k < NRAL(i->getFtype())); k++)
                    {
                        j->a[k] = renum[j->a[k]];
                    }
                }
            }
        }
        bHaveShells_ = true;
    }
    sfree(shell_atom);
}

immStatus MyMol::GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge)
{
    real qtot, mtot;

    if ((cgnr_ = generate_charge_groups(ecg, &topology_->atoms,
                                        plist_,
                                        bUsePDBcharge,
                                        &qtot, &mtot)) == NULL)
    {
        return immChargeGeneration;
    }

    if (ecg != ecgAtom)
    {
        //sort_on_charge_groups(cgnr_, &topology_->atoms,
        //                    plist_, x_, excls_, ndxfn, nmol);
    }
    return immOK;
}

void MyMol::GenerateCube(ChargeDistributionModel iChargeDistributionModel,
                         const Poldata          &pd,
                         real                    spacing,
                         const char             *reffn,
                         const char             *pcfn,
                         const char             *pdbdifffn,
                         const char             *potfn,
                         const char             *rhofn,
                         const char             *hisfn,
                         const char             *difffn,
                         const char             *diffhistfn,
                         const gmx_output_env_t *oenv)
{
    char       *gentop_version = (char *)"v0.99b";
    Resp     *  grref;

    if (NULL != gr_)
    {
        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        grref = gr_->copy();
        gr_->potcomp(pcfn, pdbdifffn, oenv);
        if ((NULL != potfn) || (NULL != hisfn) || (NULL != rhofn) ||
            ((NULL != difffn) && (NULL != reffn)))
        {
            char buf[256];

            sprintf(buf, "Potential generated by %s based on %s charges",
                    gentop_version,
                    getEemtypeName(iChargeDistributionModel));

            if (NULL != difffn)
            {
                grref->addAtomInfo(&topology_->atoms, pd);
                grref->addAtomSymmetry(symmetric_charges_);
                grref->readCube(reffn, FALSE);
                gr_->copyGrid(grref);
            }
            else
            {
                gr_->makeGrid(spacing, box, x_);
            }
            if (NULL != rhofn)
            {
                sprintf(buf, "Electron density generated by %s based on %s charges",
                        gentop_version, getEemtypeName(iChargeDistributionModel));
                gr_->calcRho();
                gr_->writeRho(rhofn, buf);
            }
            sprintf(buf, "Potential generated by %s based on %s charges",
                    gentop_version, getEemtypeName(iChargeDistributionModel));
            if (NULL != potfn)
            {
                gr_->calcPot();
                gr_->writeCube(potfn, buf);
            }
            if (NULL != hisfn)
            {
                gr_->writeHisto(hisfn, buf, oenv);
            }
            if ((NULL != difffn) || (NULL != diffhistfn))
            {
                sprintf(buf, "Potential difference generated by %s based on %s charges",
                        gentop_version,
                        getEemtypeName(iChargeDistributionModel));

                gr_->writeDiffCube(grref, difffn, diffhistfn, buf, oenv, 0);
                delete grref;
            }
        }
        delete grref;
    }
}

immStatus MyMol::getExpProps(gmx_bool bQM, gmx_bool bZero, char *lot,
                             alexandria::GaussAtomProp &gap)
{
    immStatus    imm = immOK;
    unsigned int m, nwarn = 0;
    double       value, dv0, dv298, T = -1, error, vec[3];
    tensor       quadrupole;
    std::string  myref, mylot;
    int          ia;

    if (molProp()->getPropRef(MPO_DIPOLE, (bQM ? iqmQM : iqmBoth),
                              lot, NULL, (char *)"elec",
                              &value, &error, &T, myref, mylot,
                              vec, quadrupole))
    {
        if (!bZero)
        {
            imm = immZeroDip;
        }
    }
    else
    {
        dip_exp  = value;
        dip_err  = error;
        for (m = 0; (m < DIM); m++)
        {
            mu_exp[m] = vec[m];
        }
        mu_exp2 = gmx::square(value);
        if (error <= 0)
        {
            if (debug)
            {
                fprintf(debug, "WARNING: Error for %s is %g, assuming it is 10%%.\n",
                        molProp()->getMolname().c_str(), error);
            }
            nwarn++;
            error = 0.1*value;
        }
        dip_weight = gmx::square(1.0/error);
    }
    /* Check handling of LOT */
    if (molProp()->getPropRef(MPO_DIPOLE, iqmQM,
                              (char *)mylot.c_str(), NULL, (char *)"ESP", &value, &error, &T,
                              myref, mylot, vec, quadrupole))
    {
        for (m = 0; (m < DIM); m++)
        {
            mu_esp[m] = vec[m];
        }
    }
    if (molProp()->getProp(MPO_ENERGY, (bQM ? iqmQM : iqmBoth),
                           lot, NULL, (char *)"DeltaHform", &T, &value, NULL))
    {
        gmx_fatal(FARGS, "Please rewrite code");
        Hform = value;
        Emol  = value;
        for (ia = 0; (ia < topology_->atoms.nr); ia++)
        {
            if (gap.getValue(*topology_->atoms.atomname[ia],
                             (char *)"exp", (char *)"DHf(0K)", 0, &dv0) &&
                gap.getValue(*topology_->atoms.atomname[ia],
                             (char *)"exp", (char *)"H(0K)-H(298.15K)",
                             298.15, &dv298))
            {
                Emol -= convert2gmx(dv0+dv298, eg2cHartree);
            }
            else
            {
                Emol = 0;
                break;
            }
        }
        if (ia < topology_->atoms.nr)
        {
            imm = immNoData;
        }
    }
    else
    {
        imm = immNoData;
    }
    return imm;
}

void MyMol::CalcQPol(const Poldata &pd)

{
    int     i, m, np;
    double  poltot, pol, sigpol, sptot, ereftot, eref;
    rvec    mu;

    poltot  = 0;
    sptot   = 0;
    ereftot = 0;
    np      = 0;
    clear_rvec(mu);
    for (i = 0; (i < topology_->atoms.nr); i++)
    {
        if (pd.getAtypePol(*topology_->atoms.atomtype[i], &pol, &sigpol))
        {
            np++;
            poltot += pol;
            sptot  += gmx::square(sigpol);
        }
        if (1 ==
            pd.getAtypeRefEnthalpy(*topology_->atoms.atomtype[i], &eref))
        {
            ereftot += eref;
        }
        for (m = 0; (m < DIM); m++)
        {
            mu[m] += x_[i][m]*topology_->atoms.atom[i].q;
        }
    }
    mutot_          = ENM2DEBYE*norm(mu);
    ref_enthalpy_   = ereftot;
    polarizability_ = poltot;
    sig_pol_        = sqrt(sptot/topology_->atoms.nr);
}

void MyMol::UpdateIdef(const Poldata &pd,
                       int            bondtype)
{
    int          gt;
    std::string  aai, aaj, aak, aal, params;
    int          lu;
    double       value;

    lu = string2unit(pd.getLengthUnit().c_str());
    switch (bondtype)
    {
        case ebtsBONDS:
        {
            int ftb = pd.getBondFtype();
            for (int i = 0; (i < ltop_->idef.il[ftb].nr); i += interaction_function[ftb].nratoms+1)
            {
                int         tp  = ltop_->idef.il[ftb].iatoms[i];
                int         ai  = ltop_->idef.il[ftb].iatoms[i+1];
                int         aj  = ltop_->idef.il[ftb].iatoms[i+2];
                std::string aai, aaj;
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj))
                {
                    /* Here unfortunately we need a case statement for the types */
                    double sigma, bondorder;
                    int    ntrain;
                    if (pd.searchBond(aai, aaj, &value, &sigma, &ntrain, &bondorder, params))
                    {
                        mtop_->ffparams.iparams[tp].morse.b0A = convert2gmx(value, lu);

                        std::vector<std::string> ptr = splitString(params);
                        int n = 0;
                        for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                if (n == 0)
                                {
                                    mtop_->ffparams.iparams[tp].morse.cbA = atof(pi->c_str());
                                }
                                else
                                {
                                    mtop_->ffparams.iparams[tp].morse.betaA = atof(pi->c_str());
                                }
                                n++;
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "There are no parameters for bond %s-%s in the force field", aai.c_str(), aaj.c_str());
                }
            }
        }
        case ebtsANGLES:
        {
            int fta = pd.getAngleFtype();
            for (int i = 0; (i < ltop_->idef.il[fta].nr); i += interaction_function[fta].nratoms+1)
            {
                int         tp  = ltop_->idef.il[fta].iatoms[i];
                int         ai  = ltop_->idef.il[fta].iatoms[i+1];
                int         aj  = ltop_->idef.il[fta].iatoms[i+2];
                int         ak  = ltop_->idef.il[fta].iatoms[i+3];
                std::string aai, aaj, aak;
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[ak], aak))
                {
                    double sigma;
                    int    ntrain;
                    if (pd.searchAngle( aai, aaj, aak, &value,
                                        &sigma, &ntrain, params))
                    {
                        mtop_->ffparams.iparams[tp].harmonic.rA     =
                            mtop_->ffparams.iparams[tp].harmonic.rB = value;
                        std::vector<std::string> ptr = splitString(params);
                        for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                mtop_->ffparams.iparams[tp].harmonic.krA     =
                                    mtop_->ffparams.iparams[tp].harmonic.krB = atof(pi->c_str());
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "here are no parameters for angle %s-%s-%s in the force field", aai.c_str(), aaj.c_str(), aak.c_str());
                }
            }
        }
        case ebtsPDIHS:
        case ebtsIDIHS:
        {
            for (int dd = 0; (dd < egdNR); dd++)
            {
                int ftd  = pd.getDihedralFtype(dd);
                for (int i = 0; (i < ltop_->idef.il[ftd].nr); i += interaction_function[ftd].nratoms+1)
                {
                    int         tp  = ltop_->idef.il[ftd].iatoms[i];
                    int         ai  = ltop_->idef.il[ftd].iatoms[i+1];
                    int         aj  = ltop_->idef.il[ftd].iatoms[i+2];
                    int         ak  = ltop_->idef.il[ftd].iatoms[i+3];
                    int         al  = ltop_->idef.il[ftd].iatoms[i+4];
                    std::string aai, aaj, aak;
                    if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                        pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj) &&
                        pd.atypeToBtype(*topology_->atoms.atomtype[ak], aak) &&
                        pd.atypeToBtype(*topology_->atoms.atomtype[al], aal))
                    {
                        double sigma;
                        int    ntrain;
                        if ((gt = pd.searchDihedral(dd, aai, aaj, aak, aal,
                                                    &value, &sigma, &ntrain, params)) != 0)
                        {
                            std::vector<std::string> ptr = splitString(params);
                            if (dd == egdPDIHS)
                            {
                                mtop_->ffparams.iparams[tp].pdihs.phiA = value;

                                if (ptr[0].length() > 0)
                                {
                                    mtop_->ffparams.iparams[tp].pdihs.cpA     =
                                        mtop_->ffparams.iparams[tp].pdihs.cpB =
                                            atof(ptr[0].c_str());
                                }
                                if (ptr[1].length() > 0)
                                {
                                    mtop_->ffparams.iparams[tp].pdihs.mult = atof(ptr[1].c_str());
                                }
                            }
                            else
                            {
                                mtop_->ffparams.iparams[tp].harmonic.rA     =
                                    mtop_->ffparams.iparams[tp].harmonic.rB = value;
                                if (ptr[0].length() > 0)
                                {
                                    mtop_->ffparams.iparams[tp].harmonic.krA     =
                                        mtop_->ffparams.iparams[tp].harmonic.krB = atof(ptr[0].c_str());
                                }
                            }
                        }
                    }
                    else
                    {
                        gmx_fatal(FARGS, "There are no parameters for angle %s-%s-%s in the force field", aai.c_str(), aaj.c_str(), aak.c_str());
                    }
                }
            }
        }
    }
}

}
