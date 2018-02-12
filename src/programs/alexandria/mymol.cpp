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

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"

#include "mymol.h"
#include "mymol_low.h"
      
namespace alexandria
{

static const double A2CM = E_CHARGE*1.0e-10;        /* e Angstrom to Coulomb meter */

static const double CM2D = SPEED_OF_LIGHT*1.0e+24;  /* Coulomb meter to Debye */

static inline int delta(int a, int b) { return ( a == b ) ? 1 : 0;}

static inline double e2d(double a) {return a*ENM2DEBYE;}

static void vsiteType_to_atomType(const std::string &vsiteType, std::string *atomType)  
{
    std::size_t pos = vsiteType.find("L");
   *atomType        = vsiteType.substr (0, pos);
}

const char *qTypeName(qType qt)
{
    switch(qt)
    {
    case qtESP:       return "ESP";
    case qtMulliken:  return "Mulliken";
    case qtHirshfeld: return "Hirshfeld";
    case qtCM5:       return "CM5";
    case qtCalc:      return "Calculated";
    case qtElec:      return "Electronic";
    default:
        return "Unknown charge type";
    }
}

const char *immsg(immStatus imm)
{
    static const char *msg[immNR] = {
        "Unknown status",
        "OK", 
        "Zero Dipole", 
        "No Quadrupole", 
        "Charged",
        "Atom type problem", 
        "Atom number problem", 
        "Converting from molprop",
        "Determining bond order", 
        "RESP Initialization",
        "Charge generation", 
        "Requested level of theory missing",
        "QM Inconsistency (ESP dipole does not match Elec)",
        "Not in training set", 
        "No experimental data",
        "Generating shells", 
        "Generating bonds", 
        "Communicating MolProp",
        "Zeta is zero", 
        "The number of data is lower than mindata", 
        "No Dipole moment",
        "NotSupportedBond", 
        "NotSupportedAngle", 
        "NotSupportedDihedral"
    };
    return msg[imm];
} 

MyMol::MyMol() : gvt_(evtALL)
{
    bHaveShells_       = false;
    bHaveVSites_       = false;
    bNeedVsites_       = false;
    cgnr_              = nullptr;
    immAtoms_          = immOK;
    immTopology_       = immOK;
    immCharges_        = immOK;
    shellfc_           = nullptr;
    vsite_             = nullptr;
    snew(symtab_, 1);
    open_symtab(symtab_);
    atype_ = init_atomtype();
    clear_mat(box_);
    mtop_          = nullptr;
    fr_            = nullptr;
    ltop_          = nullptr;
    mdatoms_       = nullptr;
    mp_            = new MolProp;
    state_         = new t_state;
    state_->flags |= (1<<estX);
    state_->flags |= (1<<estV);
    state_->flags |= (1<<estCGP);
    snew(enerd_, 1);
    snew(fcd_, 1);
    init_enerdata(1, 0, enerd_);
    init_nrnb(&nrnb_);
    for(int j=0; j < qtNR; j++)
    {
        clear_rvec(mu_qm_[j]);
        clear_mat(Q_qm_[j]);
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
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        mm  = topology_->atoms.atom[i].m;
        tm += mm;
        for (m = 0; (m < DIM); m++)
        {
            com[m] += mm*state_->x[i][m];
        }
    }
    if (tm > 0)
    {
        for (m = 0; m < DIM; m++)
        {
            com[m] /= tm;
        }
    }
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_dec(state_->x[i], com);
    }

    snew(bSymm, topology_->atoms.nr);
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        bSymm[i] = (norm(state_->x[i]) < toler);
        for (j = i+1; (j < topology_->atoms.nr) && !bSymm[i]; j++)
        {
            rvec_add(state_->x[i], state_->x[j], test);
            if (norm(test) < toler)
            {
                bSymm[i] = true;
                bSymm[j] = true;
            }
        }
    }
    bSymmAll = true;
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_inc(state_->x[i], com);
    }

    return bSymmAll;
}

void MyMol::findInPlaneAtoms(int ca, std::vector<int> &atoms)
{
    int bca = 0;
    /*First try to find the atom bound to ca.*/
    for (auto bi = molProp()->BeginBond(); 
         bi < molProp()->EndBond(); bi++)
    {
        if ((ca == (bi->getAj() - 1) || 
             ca == (bi->getAi() - 1)))
        {
            if (ca == (bi->getAi() - 1))
            {
                bca = (bi->getAj() - 1);
                atoms.push_back(bca);
            }
            else
            {
                bca = (bi->getAi() - 1);
                atoms.push_back(bca);
            }
        }
    }   
    /*Now try to find atoms bound to bca, except ca.*/
    for (auto bi = molProp()->BeginBond(); 
         bi < molProp()->EndBond(); bi++)
    {
        if ((ca != (bi->getAj() - 1)   && 
             ca != (bi->getAi() - 1))  &&
            (bca == (bi->getAj() - 1)  || 
             bca == (bi->getAi() - 1)))
        {
            if (bca == (bi->getAi() - 1))
            {
                atoms.push_back(bi->getAj() - 1);
            }
            else
            {
                atoms.push_back(bi->getAi() - 1);
            }
        }
    }   
}

void MyMol::findOutPlaneAtoms(int ca, std::vector<int> &atoms)
{
    for (auto bi = molProp()->BeginBond(); 
         bi < molProp()->EndBond(); bi++)
    {
        if (bi->getBondOrder() == 1  && 
            (ca == (bi->getAj() - 1) || 
             ca == (bi->getAi() - 1)))
        {
            if (ca == (bi->getAi() - 1))
            {
                atoms.push_back(bi->getAj() - 1);
            }
            else
            {
                atoms.push_back(bi->getAi() - 1);
            }
        }
    }
}

bool MyMol::IsVsiteNeeded(std::string        atype,
                          const Poldata      &pd)
{    
    auto vsite = pd.findVsite(atype);
    if (vsite != pd.getVsiteEnd())
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*
 * Make Linear Angles, Improper Dihedrals, and Virtual Sites
 */
void MyMol::MakeSpecialInteractions(const Poldata &pd,
                                    bool           bUseVsites)
{
    std::vector<std::vector<unsigned int>> bonds;
    std::vector<int> nbonds;
    t_pbc            pbc;
    real             th_toler = 175;
    real             ph_toler = 5;

    rvec *x = as_rvec_array(state_->x.data());
    
    clear_mat(box_);
    set_pbc(&pbc, epbcNONE, box_);

    bonds.resize(topology_->atoms.nr);
    for (auto bi = molProp()->BeginBond(); (bi < molProp()->EndBond()); bi++)
    {
        // Store bonds bidirectionally to get the number correct
        bonds[bi->getAi() - 1].push_back(bi->getAj() - 1);
        bonds[bi->getAj() - 1].push_back(bi->getAi() - 1);
    }
    nbonds.resize(topology_->atoms.nr);
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        nbonds[i] = bonds[i].size();
    }
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        /* Now test initial geometry */
        if ((bonds[i].size() == 2) &&
            is_linear(x[i], x[bonds[i][0]], x[bonds[i][1]],
                      &pbc, th_toler))
        {
            if (nullptr != debug)
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
                 is_planar(x[i], x[bonds[i][0]],
                           x[bonds[i][1]], x[bonds[i][2]],
                           &pbc, ph_toler))
        {
            if (nullptr != debug)
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
        if (bUseVsites)
        {
            const auto atype(*topology_->atoms.atomtype[i]);
            if (IsVsiteNeeded(atype, pd))
            {
                std::vector<int> atoms;
                auto vsite = pd.findVsite(atype);
                if(vsite->type() == evtIN_PLANE)
                {
                    atoms.push_back(i);
                    findInPlaneAtoms(i, atoms);
                    if (vsite->ncontrolatoms() == static_cast<int>(atoms.size()))
                    {
                        gvt_.addInPlane(vsite->ncontrolatoms(), 
                                        vsite->nvsite(),
                                        atoms[0], atoms[1], 
                                        atoms[2], atoms[3]);
                    }
                }
                else if (vsite->type() == evtOUT_OF_PLANE)
                { 
                    atoms.push_back(i);
                    findOutPlaneAtoms(i, atoms);
                    if (vsite->ncontrolatoms() == static_cast<int>(atoms.size()))
                    {
                        gvt_.addOutPlane(vsite->ncontrolatoms(),
                                         vsite->nvsite(),
                                         atoms[0], atoms[1], 
                                         atoms[2]);
                    }
                }
                if(!bNeedVsites_)
                {
                    bNeedVsites_ = true;
                }
            }
        }
    }
    auto anr = topology_->atoms.nr;
    gvt_.generateSpecial(pd, bUseVsites, &topology_->atoms, 
                         &x, plist_, symtab_, atype_, &excls_, state_);
    bHaveVSites_ = (topology_->atoms.nr > anr);
}

/*
 * Make Harmonic Angles, Proper Dihedrals, and 14 Pairs.
 * This needs the bonds to be F_BONDS.
 */
void MyMol::MakeAngles(bool bPairs,
                       bool bDihs)
{
    t_nextnb                            nnb;
    t_restp                             rtp;
    t_params                            plist[F_NRE];

    init_plist(plist);
    for (auto &pw : plist_)
    {
        if (F_BONDS == pw.getFtype())
        {
            pr_alloc(pw.nParam(), &(plist[F_BONDS]));
            auto i = 0;
            for (auto pi = pw.beginParam(); (pi < pw.endParam()); ++pi)
            {
                for (auto j = 0; j < MAXATOMLIST; j++)
                {
                    plist[F_BONDS].param[i].a[j] = pi->a[j];
                }
                for (auto j = 0; j < MAXFORCEPARAM; j++)
                {
                    plist[F_BONDS].param[i].c[j] = pi->c[j];
                }
                i++;
            }
            plist[F_BONDS].nr = i;
            break;
        }
    }
    /* Make Harmonic Angles and Proper Dihedrals */
    snew(excls_, topology_->atoms.nr);
    init_nnb(&nnb, topology_->atoms.nr, nexcl_); //nexcl_+2
    gen_nnb(&nnb, plist);

    print_nnb(&nnb, "NNB");
    rtp.bKeepAllGeneratedDihedrals    = bDihs;
    rtp.bRemoveDihedralIfWithImproper = bDihs;
    rtp.bGenerateHH14Interactions     = bPairs;
    rtp.nrexcl                        = nexcl_;

    gen_pad(&nnb, &(topology_->atoms), &rtp, plist, excls_, nullptr, false);

    t_blocka *EXCL;
    snew(EXCL, 1);
    generate_excl(nexcl_, topology_->atoms.nr, plist, &nnb, EXCL);
    for (int i = 0; i < EXCL->nr; i++)
    {
        int ne = EXCL->index[i+1]-EXCL->index[i];
        srenew(excls_[i].e, ne);
        excls_[i].nr = 0;
        for (auto j = EXCL->index[i]; j < EXCL->index[i+1]; j++)
        {
            if (EXCL->a[j] != i)
            {
                excls_[i].e[excls_[i].nr++] = EXCL->a[j];
            }
        }
        // Set the rest of the memory to zero
        for (auto j = excls_[i].nr; j < ne; j++)
        {
            excls_[i].e[j] = 0;
        }
    }
    done_blocka(EXCL);
    sfree(EXCL);
    if (nullptr != debug)
    {
        for (auto i = 0; i < topology_->atoms.nr; i++)
        {
            fprintf(debug, "excl %d", i);
            for (auto j = 0; j < excls_[i].nr; j++)
            {
                fprintf(debug, "  %2d", excls_[i].e[j]);
            }
            fprintf(debug, "\n");
        }
    }
    done_nnb(&nnb);

    cp_plist(&plist[F_ANGLES], F_ANGLES, eitANGLES, plist_);

    if (bDihs)
    {
        cp_plist(&plist[F_PDIHS], F_PDIHS, eitPROPER_DIHEDRALS, plist_);
    }
    if (bPairs)
    {
        cp_plist(&plist[F_LJ14], F_LJ14, eitLJ14, plist_);
    }
    for (auto i = 0; i < F_NRE; i++)
    {
        if (plist[i].nr > 0)
        {
            sfree(plist[i].param);
        }
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
        init_t_atoms(&(topology_->atoms), ci->NAtom(), false);
        snew(topology_->atoms.atomtype, ci->NAtom());
        snew(topology_->atoms.atomtypeB, ci->NAtom());

        for (auto cai = ci->BeginAtom(); (cai < ci->EndAtom()); cai++)
        {
            myunit = string2unit((char *)cai->getUnit().c_str());
            if (myunit == -1)
            {
                gmx_fatal(FARGS, "Unknown length unit '%s' for atom coordinates",
                          cai->getUnit().c_str());
            }
            cai->getCoords(&xx, &yy, &zz); 
                       
            state_->x[natom][XX] = convert2gmx(xx, myunit);
            state_->x[natom][YY] = convert2gmx(yy, myunit);
            state_->x[natom][ZZ] = convert2gmx(zz, myunit);

            double q = 0;
            for (auto qi = cai->BeginQ(); (qi < cai->EndQ()); qi++)
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
        for (auto i = 0; i < natom; i++)
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
    if (nullptr != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molProp()->getMolname().c_str(), lot, natom);
    }

    return imm;
}

immStatus MyMol::checkAtoms(const Poldata &pd)
{
    auto nmissing = 0;
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        const auto atype(*topology_->atoms.atomtype[i]);
        auto       fa = pd.findAtype(atype);
        if (fa == pd.getAtypeEnd())
        {
            printf("Could not find a force field entry for atomtype %s atom %d\n",
                   *topology_->atoms.atomtype[i], i+1);
            nmissing++;
        }
    }
    if (nmissing > 0)
    {
        return immAtomTypes;
    }
    return immOK;
}

immStatus MyMol::zeta2atoms(ChargeDistributionModel eqdModel,
                            const Poldata          &pd)
{
    /*Here, we add zeta for the core. addShell will 
      take care of the zeta for the shells later*/
    auto zeta = 0.0;
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        zeta = pd.getZeta(eqdModel, *topology_->atoms.atomtype[i], 0);       
        if (zeta == 0 && (eqdModel != eqdAXp && 
                          eqdModel != eqdAXpp && 
                          eqdModel != eqdBultinck))
        {
            printf("Zeta is zero for %s atom\n", *topology_->atoms.atomtype[i]);
            return immZeroZeta;
        }
        topology_->atoms.atom[i].zetaA = 
            topology_->atoms.atom[i].zetaB = zeta;
    }    
    return immOK;
}

immStatus MyMol::GenerateTopology(gmx_atomprop_t          ap,
                                  const Poldata          &pd,
                                  const char             *lot,
                                  ChargeDistributionModel iChargeDistributionModel,
                                  bool                    bUseVsites,
                                  bool                    bPairs,
                                  bool                    bDih,
                                  bool                    bAddShells,
                                  bool                    bBASTAT,
                                  const char             *tabfn)
{
    immStatus   imm = immOK;
    int         ftb;
    t_param     b;
    std::string btype1, btype2;
    
    if (nullptr != debug)
    {
        fprintf(debug, "Generating topology_ for %s\n", molProp()->getMolname().c_str());
    }
    nexcl_ = pd.getNexcl();
    molProp()->GenerateComposition(pd);
    if (molProp()->NAtom() <= 0)
    {
        imm = immAtomTypes;
    }
    if (immOK == imm)
    {
        snew(topology_, 1);
        init_top(topology_);
        state_change_natoms(state_, molProp()->NAtom());
        imm = GenerateAtoms(ap, lot, iChargeDistributionModel);    
    }
    if (immOK == imm)
    {
        imm = checkAtoms(pd);
    }
    if (immOK == imm)
    {
        imm = zeta2atoms(iChargeDistributionModel, pd);
    }   
    /* Store bonds in harmonic potential list first, update type later */
    ftb = F_BONDS;
    if (immOK == imm)
    {
        for (auto fs = pd.forcesBegin(); fs != pd.forcesEnd(); fs++)
        {
            if (eitBONDS == fs->iType())
            {
                ListedForceConstIterator force;
                auto lengthUnit = string2unit(fs->unit().c_str());
                if (-1 == lengthUnit)
                {
                    gmx_fatal(FARGS, "No such length unit '%s' for bonds", fs->unit().c_str());
                }
                memset(&b, 0, sizeof(b));
                for (auto bi = molProp()->BeginBond(); bi < molProp()->EndBond(); bi++)
                {
                    b.a[0] = bi->getAi() - 1;
                    b.a[1] = bi->getAj() - 1;
                    pd.atypeToBtype(*topology_->atoms.atomtype[b.a[0]], btype1);
                    pd.atypeToBtype(*topology_->atoms.atomtype[b.a[1]], btype2);
                    std::vector<std::string> atoms = {btype1, btype2};
                    if (pd.findForce(atoms, &force))
                    {
                        std::string         pp = force->params();
                        std::vector<double> dd = getDoubles(pp);
                        int                 ii = 0;
                        b.c[ii++] = convert2gmx(force->refValue(), lengthUnit);
                        for (auto &d : dd)
                        {
                            b.c[ii++] = d;
                        }
                        add_param_to_plist(plist_, ftb, eitBONDS, b);
                    }
                    else
                    {
                        // Insert dummy bond to be replaced later
                        for (auto i = 0; i < MAXFORCEPARAM; i++)
                        {
                            b.c[i] = 0;
                        }
                        add_param_to_plist(plist_, ftb, eitBONDS, b);
                    }
                }
            }
        }
        auto pw = SearchPlist(plist_, ftb);
        if (plist_.end() == pw || pw->nParam() == 0)
        {
            imm = immGenBonds;
        }
    }
    if (immOK == imm)
    {
        MakeAngles(bPairs, bDih);

        MakeSpecialInteractions(pd, bUseVsites);

        imm = updatePlist(pd, plist_, topology_, bBASTAT);
    }
    if (immOK == imm)
    {
        auto atntot = 0;
        for (auto i = 0; i < topology_->atoms.nr; i++)
        {
            auto atn = topology_->atoms.atom[i].atomnumber;
            atntot  += atn;
            for (auto m = 0; m < DIM; m++)
            {
                coc_[m] += state_->x[i][m]*atn;
            }
        }
        svmul((1.0/atntot), coc_, coc_);
    }   
    if (bAddShells && imm == immOK)
    {
        addShells(pd, iChargeDistributionModel);
    }
    if (imm == immOK)
    {
        snew(mtop_, 1);
        
        char **molnameptr = put_symtab(symtab_, molProp()->getMolname().c_str());

        do_init_mtop(pd, mtop_, molnameptr, &topology_->atoms, plist_, inputrec_, symtab_, tabfn);

        excls_to_blocka(topology_->atoms.nr, excls_, &(mtop_->moltype[0].excls));
    }
    if (bAddShells && imm == immOK)
    {
        shellfc_ = init_shell_flexcon(debug, mtop_, 0, 1, false);
    }
    if (nullptr == ltop_ && imm == immOK)
    {
        ltop_ = gmx_mtop_generate_local_top(mtop_, false);
        ltop_->idef.nthreads = 1;
    }
    return imm;
}

void MyMol::addShells(const Poldata          &pd,
                      ChargeDistributionModel iModel)
{

    int                 shell  = 0;
    int                 nshell = 0;
    double              pol    = 0; 
    double              sigpol = 0; 
      
    std::vector<int>    renum;
    std::vector<int>    inv_renum;

    char                buf[32];
    char              **newname;
    t_atoms            *newatoms;
    t_excls            *newexcls;
    PaddedRVecVector    newx;    
    t_param             p;
    
    auto polarUnit = string2unit(pd.getPolarUnit().c_str());
    if (-1 == polarUnit)
    {
        gmx_fatal(FARGS, "No such polarizability unit '%s'", pd.getPolarUnit().c_str());
    }
    
    /*Calculate the total number of particles.*/
    auto nParticles = topology_->atoms.nr;    
    for(int i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            nParticles++; // We add 1 particle as shell per Atom and Vsite
        }
    }               
    state_change_natoms(state_, nParticles);
    inv_renum.resize(nParticles, -1);
    renum.resize(topology_->atoms.nr, 0);
        
    /*Renumber the atoms.*/
    for(int i = 0; i < topology_->atoms.nr; i++)
    {
        renum[i] = i + nshell;
        inv_renum[i + nshell] = i;
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            nshell++;
        }
    }
    
    /*Add Polarization to the plist.*/
    memset(&p, 0, sizeof(p));
    for(int i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            std::string atomtype;
            vsiteType_to_atomType(*topology_->atoms.atomtype[i], &atomtype);
            auto fa    = pd.findAtype(atomtype);
            if (pd.getAtypeEnd() != fa)
            {
                auto ztype = fa->getZtype();
                if (pd.getAtypePol(atomtype, &pol, &sigpol) && (pol > 0) &&
                    (pd.getNzeta(iModel, ztype) == 2))
                {
                    p.a[0] = renum[i];
                    p.a[1] = renum[i]+1;
                    if(bHaveVSites_)
                    {
                        auto vsite = pd.findVsite(atomtype);
                        if (vsite != pd.getVsiteEnd())
                        {
                            pol /= vsite->nvsite();
                        }
                    }
                    p.c[0] = convert2gmx(pol, polarUnit);
                    add_param_to_plist(plist_, F_POLARIZATION, eitPOLARIZATION, p);
                }
                else
                {
                    gmx_fatal(FARGS, "Polarizability is %f for %s ztype %s btype %s ptype %s.\n", 
                              pol, *topology_->atoms.atomtype[i], ztype.c_str(),
                              fa->getBtype().c_str(), fa->getPtype().c_str());
                }
            }
            else
            {
                printf("Can not find atomtype %s in poldata\n", atomtype.c_str());
            }
        }
    }
    
    t_atom *shell_atom;
    snew(shell_atom, 1);
    shell_atom->ptype = eptShell;

    /* Make new atoms and x arrays. */
    snew(newatoms, 1);
    init_t_atoms(newatoms, nParticles, true);
    snew(newatoms->atomtype, nParticles);
    snew(newatoms->atomtypeB, nParticles);
    newatoms->nres = topology_->atoms.nres;
    newx.resize(newatoms->nr);
    snew(newname, newatoms->nr);
    
    /* Make a new exclusion array and put the shells in it. */
    snew(newexcls, newatoms->nr);
    
    /* Add exclusion for F_POLARIZATION. */
    auto pw = SearchPlist(plist_, F_POLARIZATION);
    if (plist_.end() != pw)
    {
        // Exclude the vsites and the atoms from their own shell.
        for (auto j = pw->beginParam(); (j < pw->endParam()); ++j)
        {
            add_excl_pair(newexcls, j->a[0], j->a[1]);
        }
        // Make a copy of the exclusions of the Atom or Vsite for the shell.
        for (auto j = pw->beginParam(); (j < pw->endParam()); ++j)
        {
            // We know that the Atom or Vsite is 0 as we added it to plist as such.
            int  i0 = inv_renum[j->a[0]];
            char buf[256];
            snprintf(buf, sizeof(buf), "Uninitialized inv_renum entry for atom %d (%d) shell %d (%d)",
                     j->a[0], inv_renum[j->a[0]],
                     j->a[1], inv_renum[j->a[1]]);
            GMX_RELEASE_ASSERT(i0 >= 0, buf);
            for (auto j0 = 0; j0 < excls_[i0].nr; j0++)
            {
                add_excl_pair(newexcls, j->a[0], renum[excls_[i0].e[j0]]);
                add_excl_pair(newexcls, j->a[1], renum[excls_[i0].e[j0]]);
            }
        }
        for (auto j = pw->beginParam(); j < pw->endParam(); ++j)
        {
            for (auto j0 = 0; j0 < newexcls[j->a[0]].nr; j0++)
            {
                add_excl_pair(newexcls, j->a[1], newexcls[j->a[0]].e[j0]);
            }
        }
    }
    
     /* Copy the old atoms to the new structures. */
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        newatoms->atom[renum[i]]      = topology_->atoms.atom[i];
        newatoms->atomname[renum[i]]  = put_symtab(symtab_, *topology_->atoms.atomname[i]);
        newatoms->atomtype[renum[i]]  = put_symtab(symtab_, *topology_->atoms.atomtype[i]);
        newatoms->atomtypeB[renum[i]] = put_symtab(symtab_, *topology_->atoms.atomtypeB[i]);
        copy_rvec(state_->x[i], newx[renum[i]]);
        newname[renum[i]] = *topology_->atoms.atomtype[i];
        t_atoms_set_resinfo(newatoms, renum[i], symtab_,
                            *topology_->atoms.resinfo[topology_->atoms.atom[i].resind].name,
                            topology_->atoms.atom[i].resind, ' ', 1, ' ');
    }
    
    for (int i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom ||
            topology_->atoms.atom[i].ptype == eptVSite)
        {
            std::string atomtype;
            auto iat = renum[i];  // Atom or Vsite
            auto j   = iat + 1;   // Shell sits next to the Atom or Vsite
            
            auto atomtypeName = get_atomtype_name(topology_->atoms.atom[i].type, atype_);
            vsiteType_to_atomType(atomtypeName, &atomtype);
            
            newatoms->atom[j]               = topology_->atoms.atom[i];
            newatoms->atom[j].m             = 0;
            newatoms->atom[j].mB            = 0;
            newatoms->atom[j].atomnumber    = topology_->atoms.atom[i].atomnumber;
            sprintf(buf, "%s_s", topology_->atoms.atom[i].elem);
            newatoms->atomname[j]           = put_symtab(symtab_, buf);
            sprintf(buf, "%s_s", atomtype.c_str());
            newname[j]                      = strdup(buf);           
            shell                           = add_atomtype(atype_, symtab_, shell_atom, buf, &p, 0, 0, 0, 0, 0, 0, 0);
            newatoms->atom[j].type          = shell;
            newatoms->atom[j].typeB         = shell;
            newatoms->atomtype[j]           = put_symtab(symtab_, buf);
            newatoms->atomtypeB[j]          = put_symtab(symtab_, buf);
            newatoms->atom[j].ptype         = eptShell;
            newatoms->atom[j].zetaA         = pd.getZeta(iModel, atomtype, 1);
            newatoms->atom[j].zetaB         = newatoms->atom[j].zetaA;
            newatoms->atom[j].resind        = topology_->atoms.atom[i].resind;            
            copy_rvec(state_->x[i], newx[j]);   
                     
            newatoms->atom[j].q = 
                newatoms->atom[j].qB = pd.getQ(iModel, atomtype, 1);
            if(bHaveVSites_)
            {
                auto vsite = pd.findVsite(atomtype);
                if (vsite != pd.getVsiteEnd())
                {
                    newatoms->atom[j].q /= vsite->nvsite();
                    newatoms->atom[j].qB = newatoms->atom[j].q;
                }
            }
        }   
    }
    
    /* Copy newatoms to atoms */
    copy_atoms(newatoms, &topology_->atoms);

    for (int i = 0; i < newatoms->nr; i++)
    {
        copy_rvec(newx[i], state_->x[i]);
        topology_->atoms.atomtype[i] = put_symtab(symtab_, newname[i]);
    }   
    sfree(newname);
    
    /* Copy exclusions, empty the original first */
    sfree(excls_);
    excls_ = newexcls;
    
    //let_shells_see_shells(excls_, &topology_->atoms, atype_);
    
    /*Now renumber atoms in all other plist interaction types */
    for (auto i = plist_.begin(); i < plist_.end(); ++i)
    {
        if (i->getFtype() != F_POLARIZATION)
        {
            for (auto j = i->beginParam(); j < i->endParam(); ++j)
            {
                for (int k = 0; k < NRAL(i->getFtype()); k++)
                {
                    j->a[k] = renum[j->a[k]];
                }
            }
        }
    }
    bHaveShells_ = true;
    sfree(shell_atom);
}

immStatus MyMol::GenerateCharges(const Poldata             &pd,
                                 const gmx::MDLogger       &mdlog,
                                 gmx_atomprop_t             ap,
                                 ChargeDistributionModel    iChargeDistributionModel,
                                 ChargeGenerationAlgorithm  iChargeGenerationAlgorithm,
                                 real                       watoms,
                                 real                       hfac,
                                 const char                *lot,
                                 bool                       bSymmetricCharges,
                                 const char                *symm_string,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 int                        maxiter,
                                 int                        maxESP,
                                 real                       tolerance,
                                 const gmx_output_env_t    *oenv)
{
    std::vector<double> qq;
    immStatus imm         = immOK;
    bool      converged   = false;
    int       iter        = 0;
    
    gmx_omp_nthreads_init(mdlog, cr, 1, 1, 0, false, false);
    GenerateGromacs(mdlog, cr, tabfn, hwinfo, iChargeDistributionModel);
    if (bSymmetricCharges)
    {
        auto bonds = SearchPlist(plist_, eitBONDS);
        if (plist_.end() != bonds)
        {
            symmetrize_charges(bSymmetricCharges, &topology_->atoms, bonds,
                               pd, ap, symm_string, symmetric_charges_);
        }
    }
    else
    {
        for (auto i = 0; i < topology_->atoms.nr; i++)
        {
            symmetric_charges_.push_back(i);
        }
    }
    switch (iChargeGenerationAlgorithm)
    {
        case eqgNONE:
            if (debug)
            {
                fprintf(debug, "WARNING! Using zero charges for %s!\n",
                        molProp()->getMolname().c_str());
            }
            for (auto i = 0; i < topology_->atoms.nr; i++)
            {
                topology_->atoms.atom[i].q  = topology_->atoms.atom[i].qB = 0;
            }
            return immOK;
        case eqgESP:
        {          
            Qgresp_.setChargeDistributionModel(iChargeDistributionModel);
            Qgresp_.setAtomWeight(watoms);
            Qgresp_.setAtomInfo(&topology_->atoms, pd, state_->x, molProp()->getCharge());
            Qgresp_.setAtomSymmetry(symmetric_charges_);
            Qgresp_.setMolecularCharge(molProp()->getCharge());
            Qgresp_.summary(debug);
            
            auto ci = molProp()->getLotPropType(lot, MPO_POTENTIAL, nullptr);
            if (ci != molProp()->EndExperiment())
            {
                int mod  = 100/maxESP;                
                int iesp = 0;
                for (auto epi = ci->BeginPotential(); epi < ci->EndPotential(); ++epi, ++iesp)
                {
                    if (Qgresp_.myWeight(iesp) == 0 || ((iesp-ci->NAtom()) % mod) != 0)
                    {
                        continue;
                    }
                    auto xu = string2unit(epi->getXYZunit().c_str());
                    auto vu = string2unit(epi->getVunit().c_str());
                    if (-1 == xu)
                    {
                        gmx_fatal(FARGS, "No such length unit '%s' for potential",
                                  epi->getXYZunit().c_str());
                    }
                    if (-1 == vu)
                    {
                        gmx_fatal(FARGS, "No such potential unit '%s' for potential",
                                  epi->getVunit().c_str());
                    }
                    Qgresp_.addEspPoint(convert2gmx(epi->getX(), xu),
                                        convert2gmx(epi->getY(), xu),
                                        convert2gmx(epi->getZ(), xu),
                                        convert2gmx(epi->getV(), vu));
                }
                if (debug)
                {
                    fprintf(debug, "Added %zu ESP points to the RESP structure.\n", Qgresp_.nEsp());
                }
            }
            double chi2[2]   = {1e8, 1e8};
            real   rrms      = 0;
            real   wtot      = 0;
            int    cur       = 0; 
            EspRms_          = 0;
            iter             = 0;

            Qgresp_.optimizeCharges();
            Qgresp_.calcPot();
            EspRms_ = chi2[cur] = Qgresp_.getRms(&wtot, &rrms);
            if (debug)
            {
                fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
            }
            do
            {
                Qgresp_.optimizeCharges();
                for (auto i = 0; i < topology_->atoms.nr; i++)
                {
                    mtop_->moltype[0].atoms.atom[i].q      =
                        mtop_->moltype[0].atoms.atom[i].qB = Qgresp_.getAtomCharge(i);
                }
                if (nullptr != shellfc_)
                {
                    computeForces(nullptr, cr);
                }
                Qgresp_.updateAtomCoords(state_->x);
                Qgresp_.calcPot();
                EspRms_ = chi2[cur] = Qgresp_.getRms(&wtot, &rrms);
                if (debug)
                {
                    fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
                }
                converged = (fabs(chi2[cur] - chi2[1-cur]) < tolerance) || (nullptr == shellfc_);
                cur       = 1-cur;
                iter++;
            }
            while ((!converged) && (iter < maxiter));
            for (auto i = 0; i < topology_->atoms.nr; i++)
            {
                topology_->atoms.atom[i].q      =
                    topology_->atoms.atom[i].qB = Qgresp_.getAtomCharge(i);
            }
            if(nullptr != oenv)
            {
                Qgresp_.plotLsq(oenv);
            }
        }
        break;
        case eqgEEM:
        {
            Qgeem_.setInfo(pd, &topology_->atoms,
                           iChargeDistributionModel,
                           hfac, molProp()->getCharge(),
                           bHaveShells_);
                                       
            auto q     = Qgeem_.q();
            auto natom = Qgeem_.natom();
            
            qq.resize(natom + 1);
            for (auto i = 0; i < natom + 1; i++)
            {           
                qq[i] = q[i][0];
            }               
            iter = 0;
            do
            {                
                if (eQGEN_OK == Qgeem_.generateCharges(nullptr,
                                                       molProp()->getMolname().c_str(),
                                                       pd, &topology_->atoms,
                                                       state_->x))
                {   
                    for (auto i = 0; i < mtop_->natoms; i++)
                    {
                        mtop_->moltype[0].atoms.atom[i].q = 
                            mtop_->moltype[0].atoms.atom[i].qB = topology_->atoms.atom[i].q;     
                        
                    }                                                        
                    if (nullptr != shellfc_)
                    {
                        computeForces(nullptr, cr);
                    }                    
                    q       = Qgeem_.q(); 
                    EemRms_ = 0;                  
                    for (auto i = 0; i < natom + 1; i++)
                    {
                        EemRms_  += gmx::square(qq[i] - q[i][0]);
                        qq[i]     = q[i][0];
                    }                    
                    EemRms_  /= natom;
                    converged = (EemRms_ < tolerance) || (nullptr == shellfc_);
                    iter++;
                }
                else
                {
                    imm = immChargeGeneration;
                }           
            }
            while(imm == immOK && (!converged) && (iter < maxiter));
            for (auto i = 0; i < mtop_->natoms; i++)
            {
                mtop_->moltype[0].atoms.atom[i].q = 
                    mtop_->moltype[0].atoms.atom[i].qB = topology_->atoms.atom[i].q;     
                
            }             
            if (!converged)
            {
                printf("EEM did not converge to %g. rms: %g\n", tolerance, EemRms_);
            }
        }
        break;
        default:
            gmx_fatal(FARGS, "Not implemented");
            break;
    }
    return imm;
}

immStatus MyMol::GenerateGromacs(const gmx::MDLogger       &mdlog,
                                 t_commrec                 *cr,
                                 const char                *tabfn,
                                 gmx_hw_info_t             *hwinfo,
                                 ChargeDistributionModel    ieqd)
{
    GMX_RELEASE_ASSERT(nullptr != mtop_, "mtop_ == nullptr. You forgot to call GenerateTopology");
    auto nalloc = 2*topology_->atoms.nr + 1;

    if (nullptr == fr_)
    {
        fr_ = mk_forcerec();
    }
    if (nullptr != tabfn)
    {
        inputrec_->coulombtype = eelUSER;
    }

    f_.resize(nalloc);
    optf_.resize(nalloc);
   
    fr_->hwinfo = hwinfo;
    init_forcerec(nullptr, mdlog, fr_, nullptr, inputrec_, mtop_, cr,
                  box_, tabfn, tabfn, nullptr, nullptr, true, -1);
    setup_bonded_threading(fr_, &ltop_->idef);
    wcycle_    = wallcycle_init(debug, 0, cr);
    mdatoms_   = init_mdatoms(nullptr, mtop_, false);
    atoms2md(mtop_, inputrec_, -1, nullptr, topology_->atoms.nr, mdatoms_);   
    
    if (nullptr != shellfc_)
    {
        make_local_shells(cr, mdatoms_, shellfc_);
    }
    if (ieqd != eqdAXs && ieqd != eqdAXps)
    {
        for (auto i = 0; i < mtop_->natoms; i++)
        {
            mdatoms_->row[i] = 0;
        }
    }
    return immOK;
}

void MyMol::computeForces(FILE *fplog, t_commrec *cr)
{
    tensor          force_vir;
    unsigned long   force_flags = ~0;
    double          t           = 0;
       
    clear_mat (force_vir);
    for (auto i = 0; i < mtop_->natoms; i++)
    {
        mdatoms_->chargeA[i] = mtop_->moltype[0].atoms.atom[i].q;  
        if (nullptr != debug)
        {
            fprintf(debug, "QQQ Setting q[%d] to %g\n", i, mdatoms_->chargeA[i]);
        }
    }
    if (nullptr != shellfc_)
    {
        auto nnodes = cr->nnodes;
        cr->nnodes  = 1;
        if (bHaveVSites_)
        {
            vsite_ = init_vsite(mtop_, cr, false);
        }
        relax_shell_flexcon(fplog, cr, false, 0, inputrec_, 
                            true, force_flags, ltop_, nullptr, 
                            enerd_, fcd_, state_,
                            &f_, force_vir, mdatoms_,
                            &nrnb_, wcycle_, nullptr,
                            &(mtop_->groups), shellfc_, 
                            fr_, false, t, nullptr,
                            vsite_);
        cr->nnodes = nnodes;     
    }
    else
    {   
        if (bHaveVSites_)
        {
            vsite_ = init_vsite(mtop_, cr, false);
        }
        do_force(fplog, cr, inputrec_, 0,
                 &nrnb_, wcycle_, ltop_,
                 &(mtop_->groups),
                 box_, &state_->x, nullptr,
                 &f_, force_vir, mdatoms_,
                 enerd_, fcd_,
                 state_->lambda, nullptr,
                 fr_, vsite_, nullptr, t,
                 nullptr, false,
                 force_flags);
    }
}

void MyMol::changeCoordinate(ExperimentIterator ei)
{
    double  xx, yy, zz;
    int     unit, natom = 0;

    for (auto eia = ei->BeginAtom(); eia < ei->EndAtom(); eia++)
    {
        unit = string2unit((char *)eia->getUnit().c_str());
        eia->getCoords(&xx, &yy, &zz);         
        state_->x[natom][XX] = convert2gmx(xx, unit);
        state_->x[natom][YY] = convert2gmx(yy, unit);
        state_->x[natom][ZZ] = convert2gmx(zz, unit);

        natom++;
    }
}

bool MyMol::getOptimizedGeometry(rvec *x)
{
    double  xx, yy, zz;
    int     unit, natom = 0;
    bool    opt = false;
    
    for (auto ei = molProp()->BeginExperiment(); 
         (!opt) && (ei < molProp()->EndExperiment()); ++ei)
    {
        if (JOB_OPT == ei->getJobtype())
        {
            for (auto eia = ei->BeginAtom(); eia < ei->EndAtom(); eia++)
            {
                unit = string2unit((char *)eia->getUnit().c_str());
                eia->getCoords(&xx, &yy, &zz);           
                x[natom][XX] = convert2gmx(xx, unit);
                x[natom][YY] = convert2gmx(yy, unit);
                x[natom][ZZ] = convert2gmx(zz, unit);           
                natom++;
            }
            opt = true;
        }
    }
    return opt;
}

void MyMol::CalcDipole()
{
    rvec mu;
    CalcDipole(mu);
    set_muQM(qtCalc, mu);
}

void MyMol::CalcDipole(rvec mu)
{
    clear_rvec(mu);
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        auto q = e2d(topology_->atoms.atom[i].q);
        for (auto m = 0; m < DIM; m++)
        {
            mu[m] += state_->x[i][m]*q;
        }
    }
}

void MyMol::CalcQuadrupole()
{
    real   r2, q;   
    rvec   r; /* distance of atoms to center of charge */
    tensor Q;
    clear_mat(Q);   
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        rvec_sub(state_->x[i], coc_, r);
        r2   = iprod(r, r);
        q    = topology_->atoms.atom[i].q;
        for (auto m = 0; m < DIM; m++)
        {
            for (auto n = 0; n < DIM; n++) 
            {
                Q[m][n] += 0.5*q*(3.0*r[m]*r[n] - r2*delta(m, n))*NM2A*A2CM*CM2D*10;
            }
        }
    }
    set_QQM(qtCalc, Q);
}

void MyMol::CalcQMbasedMoments(double *q, rvec mu, tensor Q)
{
    int   i, j;
    real  r2;
    rvec  r;  /* distance of atoms to center of charge */
    
    clear_rvec(mu);
    clear_mat(Q); 
    for (i = j = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom || 
            topology_->atoms.atom[i].ptype == eptNucleus)
        {
            rvec_sub(state_->x[i], coc_, r);
            r2       = iprod(r, r);
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += (state_->x[i][m]*e2d(q[j]));
                for (auto n = 0; n < DIM; n++)
                {
                    Q[m][n] += 0.5*q[j]*(3.0*r[m]*r[n] - r2*delta(m, n))*NM2A*A2CM*CM2D*10;
                }
            }
            j++;
        }
    }       
}

void MyMol::CalcQPol(const Poldata &pd, rvec mu)

{
    int     i, np;
    double  poltot, pol, sigpol, sptot, ereftot, eref;

    poltot  = 0;
    sptot   = 0;
    ereftot = 0;
    np      = 0;
    for (i = 0; i < topology_->atoms.nr; i++)
    {
        if (pd.getAtypePol(*topology_->atoms.atomtype[i], &pol, &sigpol))
        {
            np++;
            poltot += pol;
            sptot  += gmx::square(sigpol);
        }
        if (pd.getAtypeRefEnthalpy(*topology_->atoms.atomtype[i], &eref))
        {
            ereftot += eref;
        }
    }
    CalcDipole(mu);
    CalcQuadrupole();
    ref_enthalpy_   = ereftot;
    polarizability_ = poltot;
    sig_pol_        = sqrt(sptot/topology_->atoms.nr);
}

void MyMol::CalcAnisoPolarizability(tensor polar, double *anisoPol)
{
    auto a = gmx::square(polar[XX][XX] - polar[YY][YY]);
    auto b = gmx::square(polar[XX][XX] - polar[ZZ][ZZ]);
    auto c = gmx::square(polar[ZZ][ZZ] - polar[YY][YY]);
    auto d = 6 * (gmx::square(polar[XX][YY]) + gmx::square(polar[XX][ZZ]) + gmx::square(polar[ZZ][YY]));
    
    *anisoPol = sqrt(1/2.0) * sqrt(a + b + c + d);
}

void MyMol::CalcPolarizability(double     efield,
                               t_commrec *cr,
                               FILE      *fplog)
{
    const double        POLFAC = 29.957004; /* C.m**2.V*-1 to Å**3 */
    std::vector<double> field(DIM, 0);
    MyForceProvider    *myforce;
    rvec                mu_ref, mu_tot;
    
    myforce                            = new MyForceProvider;
    fr_->forceBufferNoVirialSummation  = new PaddedRVecVector; 
    fr_->forceBufferNoVirialSummation->resize(f_.size());     
    fr_->efield      = myforce; 
    fr_->bF_NoVirSum = true;  
    myforce->setField(field); 
    CalcDipole(mu_ref);       
    computeForces(fplog, cr);
    for (auto m = 0; m < DIM; m++)
    {
        field[m] = efield;
        myforce->setField(field);
        computeForces(fplog, cr);
        CalcDipole(mu_tot);
        for (auto n = 0; n < DIM; n++)
        {
            alpha_calc_[n][m] = ((mu_tot[n]-mu_ref[n])/efield)*(POLFAC);
        }
        isoPol_calc_     += alpha_calc_[m][m];
        field[m] = 0.0;
    }
    isoPol_calc_ /= DIM;
    CalcAnisoPolarizability(alpha_calc_, &anisoPol_calc_);
}

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];

    put_in_box(topology_->atoms.nr, box_, as_rvec_array(state_->x.data()), 0.3);
    sprintf(title, "%s processed by alexandria", molProp()->getMolname().c_str());
    write_sto_conf(fn, title, &topology_->atoms, as_rvec_array(state_->x.data()), nullptr, epbcNONE, box_);
}

void MyMol::PrintTopology(const char             *fn,
                          ChargeDistributionModel iChargeDistributionModel,
                          bool                    bVerbose,
                          const Poldata          &pd,
                          gmx_atomprop_t          aps,
                          t_commrec              *cr,
                          double                  efield,
                          const char             *lot)
{
    FILE                    *fp;
    bool                     bITP;

    bITP = (fn2ftp(fn) == efITP);
    fp   = gmx_ffopen(fn, "w");

    PrintTopology(fp, iChargeDistributionModel,
                  bVerbose, pd, aps, bITP, cr, efield, lot);

    fclose(fp);
}

void add_tensor(std::vector<std::string> *commercials,
                const char *title, const tensor &Q)
{
    char buf[256];
    snprintf(buf, sizeof(buf), "%s:\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; ( %6.2f %6.2f %6.2f )\n", 
             title,
             Q[XX][XX], Q[XX][YY], Q[XX][ZZ], 
             Q[YY][XX], Q[YY][YY], Q[YY][ZZ],
             Q[ZZ][XX], Q[ZZ][YY], Q[ZZ][ZZ]);
    commercials->push_back(buf);
}

void rotate_tensor(tensor Q, tensor Qreference)
{
    matrix rotmatrix;
    rvec   tmpvec;
    // TODO: this code is not correct!
    // the whole tensor should be taken into account, not
    // just the components. All vectors should be transformed
    // by the same matrix.
    for(int m = 0; m < DIM; m++)
    {
        if (norm(Q[m]) > 0 && norm(Qreference[m]) > 0)
        {
            calc_rotmatrix(Q[m], Qreference[m], rotmatrix);
            mvmul(rotmatrix, Q[m], tmpvec);
            copy_rvec(tmpvec, Q[m]);
        }
    }
}

void MyMol::PrintTopology(FILE                   *fp,
                          ChargeDistributionModel iChargeDistributionModel,
                          bool                    bVerbose,
                          const Poldata          &pd,
                          gmx_atomprop_t          aps,
                          bool                    bITP,
                          t_commrec              *cr,
                          double                  efield,
                          const char             *lot)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    double                   vec[DIM];
    double                   value, error, T;    
    std::string              myref, mylot;
    rvec                     mu;

    if (fp == nullptr)
    {
        return;
    }

    CalcQPol(pd, mu);
    
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
        printmol.name = strdup("Unknown");
    }
    
    printmol.nr = 1;
    
    snprintf(buf, sizeof(buf), "Total Mass = %.3f (Da)", molProp()->getMass());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Reference_Enthalpy = %.3f (kJ/mol)", ref_enthalpy_);
    commercials.push_back(buf);      
    snprintf(buf, sizeof(buf), "Total Charge = %d (e)", molProp()->getCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Charge Type  = %s\n", getEemtypeName(iChargeDistributionModel));
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Alexandria Dipole Moment (Debye):\n"
             "; ( %.2f %6.2f %6.2f ) Total= %.2f\n", 
             mu[XX], mu[YY], mu[ZZ], 
             norm(mu));
    commercials.push_back(buf);
    
    tensor myQ;
    if (molProp()->getPropRef(MPO_DIPOLE, iqmBoth, lot, "",
                              (char *)"electronic", &value, &error,
                              &T, myref, mylot, vec, myQ))
    {
        set_muQM(qtElec, vec);
        if (value > 0)
        {
            rotateDipole(mu_qm_[qtElec], mu);
        }
        snprintf(buf, sizeof(buf), "%s Dipole Moment (Debye):\n"
                 "; ( %.2f %6.2f %6.2f ) Total= %.2f\n", 
                 lot, 
                 mu_qm_[qtElec][XX], mu_qm_[qtElec][YY], mu_qm_[qtElec][ZZ], 
                 norm(mu_qm_[qtElec]));
        commercials.push_back(buf);
    }
    
    add_tensor(&commercials, 
               "Alexandria Traceless Quadrupole Moments (Buckingham)",
               QQM(qtCalc));
    
    if (molProp()->getPropRef(MPO_QUADRUPOLE, iqmBoth, lot, "",
                              (char *)"electronic", &value, &error,
                              &T, myref, mylot, vec, myQ))
    {
        set_QQM(qtElec, myQ);
        rotate_tensor(Q_qm_[qtElec], Q_qm_[qtCalc]);
        snprintf(buf, sizeof(buf), "%s Traceless Quadrupole Moments (Buckingham)", lot);
        add_tensor(&commercials, buf , Q_qm_[qtElec]);
    }
    
    snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Additive Law): %.2f +/- %.2f (A^3)\n", polarizability_, sig_pol_);
    commercials.push_back(buf);
       
    if (efield > 0 && nullptr != cr)
    {    
        CalcPolarizability(efield, cr, debug);
        add_tensor(&commercials, "Alexandria Polarizability components (A^3)", alpha_calc_);
        
        snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability (Interactive): %.2f (A^3)\n", isoPol_calc_);
        commercials.push_back(buf);
        
        snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n", anisoPol_calc_);
        commercials.push_back(buf);
    
        if (molProp()->getPropRef(MPO_POLARIZABILITY, iqmBoth, lot, "",
                                  (char *)"electronic", &isoPol_elec_, &error,
                                  &T, myref, mylot, vec, alpha_elec_))
        {
            CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
            snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", lot);
            add_tensor(&commercials, buf, alpha_elec_);
       
            snprintf(buf, sizeof(buf), "%s Isotropic Polarizability: %.2f (A^3)\n", lot, isoPol_elec_);
            commercials.push_back(buf);    
            snprintf(buf, sizeof(buf), "%s Anisotropic Polarizability: %.2f (A^3)\n", lot, anisoPol_elec_);
            commercials.push_back(buf);
        }        
    }
          
    print_top_header(fp, pd, aps, bHaveShells_, iChargeDistributionModel, commercials, bITP);
    write_top(fp, printmol.name, &topology_->atoms, false,
              plist_, excls_, atype_, cgnr_, nexcl_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, getForceField().c_str(), nullptr, 0, nullptr, 1, &printmol);
    }
    if (bVerbose)
    {
        for (auto &p : plist_)
        {
            if (p.nParam() > 0)
            {
                printf("There are %4d %s interactions\n", p.nParam(),
                       interaction_function[p.getFtype()].name);
            }
        }
        for (auto i = commercials.begin(); (i < commercials.end()); ++i)
        {
            printf("%s\n", i->c_str());
        }
    }

    sfree(printmol.name);
}

immStatus MyMol::GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge)
{
    real qtot, mtot;
    if ((cgnr_ = generate_charge_groups(ecg, &topology_->atoms,
                                        plist_,
                                        bUsePDBcharge,
                                        &qtot, &mtot)) == nullptr)
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
    if ((nullptr  != potfn) || (nullptr != hisfn) || (nullptr != rhofn) ||
        ((nullptr != difffn) && (nullptr != reffn)))
    {
        char      buf[256];
        char     *gentop_version = (char *)"v0.99b";
        QgenResp  grref;

        Qgresp_.potcomp(pcfn, pdbdifffn, oenv);

        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        grref = Qgresp_;

        sprintf(buf, "Potential generated by %s based on %s charges",
                gentop_version,
                getEemtypeName(iChargeDistributionModel));

        if (nullptr != difffn)
        {
            grref.setAtomInfo(&topology_->atoms, pd, state_->x, molProp()->getCharge());
            grref.setAtomSymmetry(symmetric_charges_);
            grref.readCube(reffn, FALSE);
            Qgresp_ = grref;
        }
        else
        {
            Qgresp_.makeGrid(spacing, box_, as_rvec_array(state_->x.data()));
        }
        if (nullptr != rhofn)
        {
            sprintf(buf, "Electron density generated by %s based on %s charges",
                    gentop_version, getEemtypeName(iChargeDistributionModel));
            Qgresp_.calcRho();
            Qgresp_.writeRho(rhofn, buf);
        }
        sprintf(buf, "Potential generated by %s based on %s charges",
                gentop_version, getEemtypeName(iChargeDistributionModel));
        if (nullptr != potfn)
        {
            Qgresp_.calcPot();
            Qgresp_.writeCube(potfn, buf);
        }
        if (nullptr != hisfn)
        {
            Qgresp_.writeHisto(hisfn, buf, oenv);
        }
        if ((nullptr != difffn) || (nullptr != diffhistfn))
        {
            sprintf(buf, "Potential difference generated by %s based on %s charges",
                    gentop_version,
                    getEemtypeName(iChargeDistributionModel));

            Qgresp_.writeDiffCube(grref, difffn, diffhistfn, buf, oenv, 0);
        }
    }
}

void MyMol::rotateDipole(rvec mu, rvec muReference)
{
    matrix rotmatrix;
    rvec   tmpvec;
    calc_rotmatrix(mu, muReference, rotmatrix);
    mvmul(rotmatrix, mu, tmpvec);
    copy_rvec(tmpvec, mu);
}

void MyMol::setQandMoments(qType qt, int natom, double q[])
{
    int i, j;
    
    if (natom > 0)
    {
        charge_QM_[qt].resize(natom);
        for (i = j = 0; i < topology_->atoms.nr; i++)
        {
            if (topology_->atoms.atom[i].ptype == eptAtom || 
                topology_->atoms.atom[i].ptype == eptNucleus)
            {
                charge_QM_[qt][j] = q[j];
                j++;
            }
        }
        CalcQMbasedMoments(q, mu_qm_[qt], Q_qm_[qt]);
    }
}

immStatus MyMol::getExpProps(gmx_bool bQM, gmx_bool bZero,
                             gmx_bool bZPE, const char *lot,
                             const Poldata &pd)
{
    int          ia    = 0;
    int          natom = 0;
    unsigned int nwarn = 0;
    double       value = 0;
    double       Hatom = 0;
    double       ZPE   = 0;
    double       error = 0; 
    double       T     = -1; 
    double       vec[DIM];
    tensor       quadrupole = {{0,0,0},{0,0,0},{0,0,0}};
    tensor       polar      = {{0,0,0},{0,0,0},{0,0,0}};
    std::string  myref;
    std::string  mylot;
    bool         esp_dipole_found  = false;
    
    
    for (auto i = 0; i < topology_->atoms.nr; i++)
    {
        if (topology_->atoms.atom[i].ptype == eptAtom || 
            topology_->atoms.atom[i].ptype == eptNucleus)
        {
            natom++;
        }
    }   
    double q[natom];
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "", 
                              (char *)"ESP charges",
                              &value, &error, &T,
                              myref, mylot, q, quadrupole))
    {
        setQandMoments(qtESP, natom, q);
        esp_dipole_found = true;
    }
    T = -1;
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "", 
                              (char *)"Mulliken charges",
                              &value, &error, &T,
                              myref, mylot, q, quadrupole))
    {
        setQandMoments(qtMulliken, natom, q);
        if (esp_dipole_found && dipQM(qtMulliken) > 0)
        {
            rotate_tensor(Q_qm_[qtMulliken], Q_qm_[qtESP]);
        }
    }
    T = -1;
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "", 
                              (char *)"Hirshfeld charges",
                              &value, &error, &T,
                              myref, mylot, q, quadrupole))
    {
        setQandMoments(qtHirshfeld, natom, q);
        if (esp_dipole_found && dipQM(qtHirshfeld) > 0)
        {
            rotate_tensor(Q_qm_[qtHirshfeld], Q_qm_[qtESP]);
        }

    }
    T = -1;
    if (molProp()->getPropRef(MPO_CHARGE, iqmQM,
                              (char *)mylot.c_str(), "", 
                              (char *)"CM5 charges",
                              &value, &error, &T,
                              myref, mylot, q, quadrupole))
    {
        setQandMoments(qtCM5, natom, q);
        if (esp_dipole_found && dipQM(qtCM5) > 0)
        {
            rotate_tensor(Q_qm_[qtCM5], Q_qm_[qtESP]);
        }

    } 
    T = -1;
    immStatus imm = immOK;
    if (molProp()->getProp(MPO_ENERGY, (bQM ? iqmQM : iqmBoth),
                           lot, "", (char *)"DeltaHform", &value, &error, &T))
    {
        Hform_ = value;
        Emol_  = value;
        for (ia = 0; ia < topology_->atoms.nr; ia++)
        {
            if (topology_->atoms.atom[ia].ptype == eptAtom ||
                topology_->atoms.atom[ia].ptype == eptNucleus)
            {
                if (pd.getAtypeRefEnthalpy(*topology_->atoms.atomtype[ia], &Hatom))
                {
                    Emol_ -= Hatom;
                }
                else
                {
                    if (debug)
                    {
                        fprintf(debug, "WARNING: NO reference enthalpy for molecule %s.\n",
                                molProp()->getMolname().c_str());
                    }
                    Emol_ = 0;
                    break;
                }
            }
        }
        if (bZPE)
        {

            if (molProp()->getProp(MPO_ENERGY, iqmBoth, lot, "",
                                   (char *)"ZPE", &ZPE, &error, &T))
            {
                Emol_ -= ZPE;
            }
            else
            {
                if (debug)
                {
                    fprintf(debug, "WARNING: NO Zero-point energy for molecule %s.\n",
                            molProp()->getMolname().c_str());
                }
            }
        }
        if (ia < topology_->atoms.nr)
        {
            imm = immNoData;
        }
    }
    T = -1;
    if (molProp()->getPropRef(MPO_DIPOLE, iqmQM , lot, "", 
                              (char *)"electronic",
                              &value, &error, &T, myref, mylot,
                              vec, quadrupole))
    {       
        dip_exp_  = value;
        dip_err_  = error;
        set_muQM(qtElec, vec);
        
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
        dip_weight_ = gmx::square(1.0/error);
        
        if (!bZero && (dipQM(qtElec) - 0.0) < 1e-2)
        {
          imm = immZeroDip;
        }
        if (immOK == imm && esp_dipole_found)
        {
            rotateDipole(mu_qm_[qtElec], mu_qm_[qtESP]);
        }
    }
    else
    {
      imm = immNoDipole;
    }
    T = -1;
    if (molProp()->getPropRef(MPO_QUADRUPOLE, iqmQM,
                              lot, "", (char *)"electronic",
                              &value, &error, &T, myref, mylot,
                              vec, quadrupole))
    {
        set_QQM(qtElec, quadrupole);
        if (immOK == imm && esp_dipole_found && norm(mu_qm_[qtElec]) > 0)
        {
            rotate_tensor(Q_qm_[qtElec], Q_qm_[qtESP]);
        }
    }
    T = -1;  
    if (molProp()->getPropRef(MPO_POLARIZABILITY, iqmQM, 
                              lot, "", (char *)"electronic", 
                              &isoPol_elec_, &error, &T, 
                              myref, mylot, vec, polar))
    {        
        copy_mat(polar, alpha_elec_);
        CalcAnisoPolarizability(alpha_elec_, &anisoPol_elec_);
    } 
    return imm;
}

void MyMol::UpdateIdef(const Poldata   &pd,
                       InteractionType  iType)
{
    std::string              aai, aaj, aak, aal, params;
    std::vector<std::string> atoms, ptr;
    int                      lu, n;
    size_t                   ntrain            = 0;
    double                   value             = 0;
    double                   sigma             = 0;
    double                   r13               = 0;
    double                   relative_position = 0;

    switch (iType)
    {
        case eitBONDS:
        {
            auto fs = pd.findForces(iType);
            lu = string2unit(fs->unit().c_str());
            if (-1 == lu)
            {
                gmx_fatal(FARGS, "Unknown length unit '%s' for bonds",
                          fs->unit().c_str());
            }
            auto ftb = fs->fType();
            for (auto i = 0; i < ltop_->idef.il[ftb].nr; i += interaction_function[ftb].nratoms+1)
            {
                auto  tp  = ltop_->idef.il[ftb].iatoms[i];
                auto  ai  = ltop_->idef.il[ftb].iatoms[i+1];
                auto  aj  = ltop_->idef.il[ftb].iatoms[i+2];
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj))
                {
                    atoms = {aai, aaj};
                    if (pd.searchForce(atoms, params, &value, &sigma, &ntrain, iType))
                    {
                        mtop_->ffparams.iparams[tp].morse.b0A     =
                            mtop_->ffparams.iparams[tp].morse.b0B = convert2gmx(value, lu);

                        ltop_->idef.iparams[tp].morse.b0A     =
                            ltop_->idef.iparams[tp].morse.b0B = convert2gmx(value, lu);

                        ptr = gmx::splitString(params);
                        n   = 0;
                        for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                if (n == 0)
                                {
                                    mtop_->ffparams.iparams[tp].morse.cbA     =
                                        mtop_->ffparams.iparams[tp].morse.cbB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].morse.cbA     =
                                        ltop_->idef.iparams[tp].morse.cbB = gmx::doubleFromString(pi->c_str());
                                }
                                else
                                {
                                    mtop_->ffparams.iparams[tp].morse.betaA     =
                                        mtop_->ffparams.iparams[tp].morse.betaB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].morse.betaA     =
                                        ltop_->idef.iparams[tp].morse.betaB = gmx::doubleFromString(pi->c_str());
                                }
                                n++;
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "There are no parameters for bond %s-%s in the force field",
                              aai.c_str(), aaj.c_str());
                }
            }
        }
        break;
        case eitANGLES:
        {
            auto fs  = pd.findForces(iType);
            auto fta = fs->fType();
            for (auto i = 0; i < ltop_->idef.il[fta].nr; i += interaction_function[fta].nratoms+1)
            {
                auto  tp  = ltop_->idef.il[fta].iatoms[i];
                auto  ai  = ltop_->idef.il[fta].iatoms[i+1];
                auto  aj  = ltop_->idef.il[fta].iatoms[i+2];
                auto  ak  = ltop_->idef.il[fta].iatoms[i+3];
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[ak], aak))
                {
                    atoms = {aai, aaj, aak};
                    if (pd.searchForce(atoms, params, &value, &sigma, &ntrain, iType))
                    {

                        r13 = calc_r13(pd, aai, aaj, aak, value);

                        mtop_->ffparams.iparams[tp].u_b.thetaA     =
                            mtop_->ffparams.iparams[tp].u_b.thetaB = value;

                        ltop_->idef.iparams[tp].u_b.thetaA     =
                            ltop_->idef.iparams[tp].u_b.thetaB = value;

                        mtop_->ffparams.iparams[tp].u_b.r13A     =
                            mtop_->ffparams.iparams[tp].u_b.r13B = r13;

                        ltop_->idef.iparams[tp].u_b.r13A     =
                            ltop_->idef.iparams[tp].u_b.r13B = r13;

                        ptr = gmx::splitString(params);
                        n   = 0;
                        for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                if (n == 0)
                                {
                                    mtop_->ffparams.iparams[tp].u_b.kthetaA     =
                                        mtop_->ffparams.iparams[tp].u_b.kthetaB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].u_b.kthetaA     =
                                        ltop_->idef.iparams[tp].u_b.kthetaB = gmx::doubleFromString(pi->c_str());
                                }
                                else
                                {
                                    mtop_->ffparams.iparams[tp].u_b.kUBA     =
                                        mtop_->ffparams.iparams[tp].u_b.kUBB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].u_b.kUBA     =
                                        ltop_->idef.iparams[tp].u_b.kUBB = gmx::doubleFromString(pi->c_str());
                                }
                                n++;
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "There are no parameters for angle %s-%s-%s in the force field",
                              aai.c_str(), aaj.c_str(), aak.c_str());
                }
            }
        }
        break;
        case eitLINEAR_ANGLES:
        {
            auto fs  = pd.findForces(iType);
            auto fta = fs->fType();
            for (auto i = 0; i < ltop_->idef.il[fta].nr; i += interaction_function[fta].nratoms+1)
            {
                auto  tp  = ltop_->idef.il[fta].iatoms[i];
                auto  ai  = ltop_->idef.il[fta].iatoms[i+1];
                auto  aj  = ltop_->idef.il[fta].iatoms[i+2];
                auto  ak  = ltop_->idef.il[fta].iatoms[i+3];
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[ak], aak))
                {
                    atoms = {aai, aaj, aak};
                    if (pd.searchForce(atoms, params, &value, &sigma, &ntrain, iType))
                    {
                        r13 = calc_r13(pd, aai, aaj, aak, value);

                        relative_position = calc_relposition(pd, aai, aaj, aak);

                        mtop_->ffparams.iparams[tp].linangle.aA     =
                            mtop_->ffparams.iparams[tp].linangle.aB = relative_position;

                        ltop_->idef.iparams[tp].linangle.aA     =
                            ltop_->idef.iparams[tp].linangle.aB = relative_position;

                        mtop_->ffparams.iparams[tp].linangle.r13A     =
                            mtop_->ffparams.iparams[tp].linangle.r13B = r13;

                        ltop_->idef.iparams[tp].linangle.r13A     =
                            ltop_->idef.iparams[tp].linangle.r13B = r13;

                        ptr = gmx::splitString(params);
                        n   = 0;
                        for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                if (n == 0)
                                {
                                    mtop_->ffparams.iparams[tp].linangle.klinA     =
                                        mtop_->ffparams.iparams[tp].linangle.klinB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].linangle.klinA     =
                                        ltop_->idef.iparams[tp].linangle.klinB = gmx::doubleFromString(pi->c_str());
                                }
                                else
                                {
                                    mtop_->ffparams.iparams[tp].linangle.kUBA     =
                                        mtop_->ffparams.iparams[tp].linangle.kUBB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].linangle.kUBA     =
                                        ltop_->idef.iparams[tp].linangle.kUBB = gmx::doubleFromString(pi->c_str());
                                }
                                n++;
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "There are no parameters for linear angle %s-%s-%s in the force field",
                              aai.c_str(), aaj.c_str(), aak.c_str());
                }
            }
        }
        break;
        case eitPROPER_DIHEDRALS:
        {
            auto fs  = pd.findForces(iType);
            auto ftd = fs->fType();
            for (auto i = 0; i < ltop_->idef.il[ftd].nr; i += interaction_function[ftd].nratoms+1)
            {
                auto  tp  = ltop_->idef.il[ftd].iatoms[i];
                auto  ai  = ltop_->idef.il[ftd].iatoms[i+1];
                auto  aj  = ltop_->idef.il[ftd].iatoms[i+2];
                auto  ak  = ltop_->idef.il[ftd].iatoms[i+3];
                auto  al  = ltop_->idef.il[ftd].iatoms[i+4];
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[ak], aak) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[al], aal))
                {
                    atoms = {aai, aaj, aak, aal};
                    if (pd.searchForce(atoms, params, &value, &sigma, &ntrain, iType))
                    {
                        mtop_->ffparams.iparams[tp].pdihs.phiA     =
                            mtop_->ffparams.iparams[tp].pdihs.phiB = value;

                        ltop_->idef.iparams[tp].pdihs.phiA     =
                            ltop_->idef.iparams[tp].pdihs.phiB = value;

                        ptr = gmx::splitString(params);
                        n   = 0;
                        for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                if (n == 0)
                                {
                                    mtop_->ffparams.iparams[tp].pdihs.cpA     =
                                        mtop_->ffparams.iparams[tp].pdihs.cpB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].pdihs.cpA     =
                                        ltop_->idef.iparams[tp].pdihs.cpB = gmx::doubleFromString(pi->c_str());
                                }
                                else
                                {
                                    /*Multiplicity for Proper Dihedral must be integer
                                       This assumes that the second paramter is Multiplicity*/
                                    mtop_->ffparams.iparams[tp].pdihs.mult = atoi(pi->c_str());

                                    ltop_->idef.iparams[tp].pdihs.mult = atoi(pi->c_str());
                                }
                                n++;
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "There are no parameters for proper dihedral %s-%s-%s-%s in the force field",
                              aai.c_str(), aaj.c_str(), aak.c_str(), aal.c_str());
                }
            }
        }
        break;
        case eitIMPROPER_DIHEDRALS:
        {
            auto fs  = pd.findForces(iType);
            auto ftd = fs->fType();
            for (auto i = 0; i < ltop_->idef.il[ftd].nr; i += interaction_function[ftd].nratoms+1)
            {
                auto  tp  = ltop_->idef.il[ftd].iatoms[i];
                auto  ai  = ltop_->idef.il[ftd].iatoms[i+1];
                auto  aj  = ltop_->idef.il[ftd].iatoms[i+2];
                auto  ak  = ltop_->idef.il[ftd].iatoms[i+3];
                auto  al  = ltop_->idef.il[ftd].iatoms[i+4];
                if (pd.atypeToBtype(*topology_->atoms.atomtype[ai], aai) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[aj], aaj) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[ak], aak) &&
                    pd.atypeToBtype(*topology_->atoms.atomtype[al], aal))
                {
                    atoms = {aai, aaj, aak, aal};
                    if (pd.searchForce(atoms, params, &value, &sigma, &ntrain, iType))
                    {
                        mtop_->ffparams.iparams[tp].harmonic.rA     =
                            mtop_->ffparams.iparams[tp].harmonic.rB = value;

                        ltop_->idef.iparams[tp].harmonic.rA     =
                            ltop_->idef.iparams[tp].harmonic.rB = value;

                        ptr = gmx::splitString(params);
                        n   = 0;
                        for (auto pi = ptr.begin(); pi < ptr.end(); ++pi)
                        {
                            if (pi->length() > 0)
                            {
                                if (n == 0)
                                {
                                    mtop_->ffparams.iparams[tp].harmonic.krA     =
                                        mtop_->ffparams.iparams[tp].harmonic.krB = gmx::doubleFromString(pi->c_str());

                                    ltop_->idef.iparams[tp].harmonic.krA     =
                                        ltop_->idef.iparams[tp].harmonic.krB = gmx::doubleFromString(pi->c_str());
                                }
                                n++;
                            }
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "There are no parameters for imporper proper dihedral %s-%s-%s-%s in the force field",
                              aai.c_str(), aaj.c_str(), aak.c_str(), aal.c_str());
                }
            }
        }
        break;
        case eitVDW:
        {
            auto ftv   = pd.getVdwFtype();
            auto ntype = mtop_->ffparams.atnr;
            for (auto i = 0; (i < ntype); i++)
            {
                for (auto j = 0; (j < ntype); j++)
                {
                    auto idx = ntype*i+j;
                    switch (ftv)
                    {
                        case F_LJ:
                        {
                            double c6, c12;
                            getLjParams(pd,
                                        *(topology_->atoms.atomtype[i]),
                                        *(topology_->atoms.atomtype[j]),
                                        &c6, &c12);
                            mtop_->ffparams.iparams[idx].lj.c6  = c6;
                            mtop_->ffparams.iparams[idx].lj.c12 = c12;
                        }
                        break;
                        case F_BHAM:
                        {
                            double a, b, c;
                            getBhamParams(pd,
                                          *(topology_->atoms.atomtype[i]),
                                          *(topology_->atoms.atomtype[j]),
                                          &a, &b, &c);
                            mtop_->ffparams.iparams[idx].bham.a = a;
                            mtop_->ffparams.iparams[idx].bham.b = b;
                            mtop_->ffparams.iparams[idx].bham.c = c;
                        }
                        break;
                        default:
                            fprintf(stderr, "Invalid van der waals type %s\n",
                                    pd.getVdwFunction().c_str());
                    }
                }
            }
        }
        break;
        case eitLJ14:
        case eitPOLARIZATION:
        {
            auto pw = SearchPlist(plist_, iType);
            if (plist_.end() != pw)
            {
                auto ft = pw->getFtype();
                for (auto i = 0; i < ltop_->idef.il[ft].nr; i += interaction_function[ft].nratoms+1)
                {
                    auto tp  = ltop_->idef.il[ft].iatoms[i];
                    auto ai  = ltop_->idef.il[ft].iatoms[i+1];
                    if (pd.getAtypePol(*topology_->atoms.atomtype[ai], &value, &sigma))
                    {
                        mtop_->ffparams.iparams[tp].polarize.alpha = value;
                        ltop_->idef.iparams[tp].polarize.alpha     = value;
                    }
                }
            } 
        }
        break;
        case eitVSITE2:
        case eitVSITE3FAD:
        case eitVSITE3OUT:
        case eitCONSTR:
        case eitNR:
            break;
    }
}
}// namespace alexandria
