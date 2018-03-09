/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "gentop_vsite.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <map>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/smalloc.h"

#include "mymol_low.h"
#include "plistwrapper.h"
#include "poldata.h"

namespace alexandria
{

gv_inplane::gv_inplane(int natom, int nvsite, int ca, 
                       int bca,   int bbca1,  int bbca2)
    :
      natom_(natom),
      nvsite_(nvsite),
      ca_(ca),
      bca_(bca),
      bbca1_(bbca1),
      bbca2_(bbca2)
{}

gv_inplaneIterator GentopVsites::findInPlane(int natom, int nvsite, int ca, 
                                             int bca,   int bbca1,  int bbca2)
{
    gv_inplaneIterator b = inplaneBegin(), e = inplaneEnd();
    return std::find_if(b, e, [natom, nvsite, ca, bca, bbca1, bbca2](const gv_inplane &inplane)
                        {
                            return (natom   == inplane.natom()  && 
                                    nvsite  == inplane.nvsite() &&
                                    ca      == inplane.ca()     &&
                                    bca     == inplane.bca()    &&
                                    ((bbca1 == inplane.bbca1()  && 
                                      bbca2 == inplane.bbca2()) ||
                                     (bbca1 == inplane.bbca2()  && 
                                      bbca2 == inplane.bbca1())));
                        });
}
        
gv_inplaneConstIterator GentopVsites::findInPlane(int natom, int nvsite, int ca, 
                                                  int bca,   int bbca1,  int bbca2) const
{
    gv_inplaneConstIterator b = inplaneBegin(), e = inplaneEnd();
    return std::find_if(b, e, [natom, nvsite, ca, bca, bbca1, bbca2](const gv_inplane &inplane)
                        {
                            return (natom   == inplane.natom()  && 
                                    nvsite  == inplane.nvsite() &&
                                    ca      == inplane.ca()     &&
                                    bca     == inplane.bca()    &&
                                    ((bbca1 == inplane.bbca1()  && 
                                      bbca2 == inplane.bbca2()) ||
                                     (bbca1 == inplane.bbca2()  && 
                                      bbca2 == inplane.bbca1())));
                        });
}

gv_inplaneIterator GentopVsites::findInPlane(int nvsite, int ca)
{
    gv_inplaneIterator b = inplaneBegin(), e = inplaneEnd();
    return std::find_if(b, e, [nvsite, ca](const gv_inplane &inplane)
                        {
                            return (nvsite  == inplane.nvsite() && ca == inplane.ca());
                        });
}

gv_inplaneConstIterator GentopVsites::findInPlane(int nvsite, int ca) const
{
    gv_inplaneConstIterator b = inplaneBegin(), e = inplaneEnd();
    return std::find_if(b, e, [nvsite, ca](const gv_inplane &inplane)
                        {
                            return (nvsite  == inplane.nvsite() && ca == inplane.ca());
                        });
}

void GentopVsites::addInPlane(int natom, int nvsite, int ca, 
                              int bca,   int bbca1,  int bbca2)
{
    auto inplane = findInPlane(natom, nvsite, ca, bca, bbca1, bbca2);
    if(inplane != inplaneEnd())
    {
        return;
    }
    gv_inplane p(natom, nvsite, ca, bca, bbca1, bbca2);
    inplane_.push_back(std::move(p));
}

gv_outplane::gv_outplane(int natom, int nvsite,int ca, int bca1, int bca2)
    :
      natom_(natom),
      nvsite_(nvsite),
      ca_(ca),
      bca1_(bca1),
      bca2_(bca2)
{}

gv_outplaneIterator GentopVsites::findOutPlane(int natom, int nvsite, int ca, int bca1, int bca2)
{
    gv_outplaneIterator b = outplaneBegin(), e = outplaneEnd();
    return std::find_if(b, e, [natom, nvsite, ca, bca1, bca2](const gv_outplane &outplane)
                        {
                            return (natom  == outplane.natom()  && 
                                    ca     == outplane.ca()     &&
                                    nvsite == outplane.nvsite() &&
                                    ((bca1 == outplane.bca1() && bca2 == outplane.bca2()) ||
                                     (bca1 == outplane.bca2() && bca2 == outplane.bca1())));
                        });
}
        
gv_outplaneConstIterator GentopVsites::findOutPlane(int natom, int nvsite, int ca, int bca1, int bca2) const
{
    gv_outplaneConstIterator b = outplaneBegin(), e = outplaneEnd();
    return std::find_if(b, e, [natom, nvsite, ca, bca1, bca2](const gv_outplane &outplane)
                        {
                            return (natom  == outplane.natom()  && 
                                    ca     == outplane.ca()     &&
                                    nvsite == outplane.nvsite() &&
                                    ((bca1 == outplane.bca1() && bca2 == outplane.bca2()) ||
                                     (bca1 == outplane.bca2() && bca2 == outplane.bca1())));
                        });

}

gv_outplaneIterator GentopVsites::findOutPlane(int nvsite, int ca)
{
    gv_outplaneIterator b = outplaneBegin(), e = outplaneEnd();
    return std::find_if(b, e, [nvsite, ca](const gv_outplane &outplane)
                        {
                            return (ca == outplane.ca() && nvsite == outplane.nvsite());
                        });
}
        
gv_outplaneConstIterator GentopVsites::findOutPlane(int nvsite, int ca) const
{
    gv_outplaneConstIterator b = outplaneBegin(), e = outplaneEnd();
    return std::find_if(b, e, [nvsite, ca](const gv_outplane &outplane)
                        {
                            return (ca == outplane.ca() && nvsite == outplane.nvsite());
                        });

}

void GentopVsites::addOutPlane(int natom, int nvsite, int ca, int bca1, int bca2)
{
    auto outplane = findOutPlane(natom, nvsite, ca, bca1, bca2);    
    if (outplane != outplaneEnd())
    {
        return;
    }
    gv_outplane p(natom, nvsite, ca, bca1, bca2);
    outplane_.push_back(std::move(p));
}


int GentopVsites::nVsites()
{
    int n = 0;
    for (const auto &oplane : outplane_)
    {
        n += oplane.nvsite();
    }
    for (const auto &iplane : inplane_)
    {
        n += iplane.nvsite();
    }
    return n;
}

void GentopVsites::addLinear(int ai, int aj, int ak)
{
    unsigned int i;

    for (i = 0; (i < linear_.size()); i++)
    {
        if (linear_[i].a[1] == aj)
        {
            if (((linear_[i].a[0] == ai)  && (linear_[i].a[2] == ak)) ||
                ((linear_[i].a[0] == aj)  && (linear_[i].a[2] == ai)))
            {
                break;
            }
        }
    }
    if (i == linear_.size())
    {
        gv_linear l;
        l.nline = -1;
        l.a[0]  = ai;
        l.a[1]  = aj;
        l.a[2]  = ak;
        linear_.push_back(l);
    }
}

void GentopVsites::addPlanar(int ai, int aj, int ak, int al, int nbonds[])
{
    unsigned int i;

    for (i = 0; (i < planar_.size()); i++)
    {
        if (((planar_[i].a[0] == ai) && (planar_[i].a[1] == aj)  &&
             (planar_[i].a[2] == ak) && (planar_[i].a[3] == al)) ||
            ((planar_[i].a[0] == al) && (planar_[i].a[1] == ak)  &&
             (planar_[i].a[2] == aj) && (planar_[i].a[3] == ai)))
        {
            break;
        }
    }
    if (i == planar_.size())
    {
        gv_planar p;

        p.a[0] = ai;
        p.a[1] = aj;
        p.a[2] = ak;
        p.a[3] = al;
        for (unsigned int j = 0; (j < 4); j++)
        {
            p.nb[j] = nbonds[p.a[j]];
        }
        planar_.push_back(p);
    }
}

void GentopVsites::addRingPlanar(int natom, int aa[], int nbonds[])
{
    unsigned int i;

    for (i = 0; (i < ringplanar_.size()); i++)
    {
        bool bSame = true;
        for (int j = 0; (j < natom); j++)
        {
            bSame = bSame && (ringplanar_[i].a[j] == aa[j]);
        }
        if (bSame)
        {
            break;
        }
    }
    if (i == ringplanar_.size())
    {
        gv_ringplanar rp;
        for (int j = 0; (j < natom); j++)
        {
            rp.a[j]  = aa[j];
            rp.nb[j] = nbonds[rp.a[j]];
        }
        ringplanar_.push_back(rp);
    }
}

static void calc_vsite2parm(t_atoms                   *atoms,
                            std::vector<PlistWrapper> &plist,
                            rvec                      **x,
                            gv_linear                 *gvl,
                            t_symtab                  *symtab,
                            gpp_atomtype_t             atype)
{
    int              i, j, natoms, mt;
    const   char    *ml = "ML";
    double           mI, mJ, mK, mL, mT, com, I;
    double           rB, rC, rD, rVV, mV, ac[4];
    rvec             dx;
    t_param          pp, nbml;
    t_atom           aml;

    if (gvl->nline <= 0)
    {
        return;
    }
    memset(&nbml, 0, sizeof(nbml));
    rvec_sub((*x)[gvl->a[0]], (*x)[gvl->a[1]], dx);
    rB    = norm(dx);
    rvec_sub((*x)[gvl->a[1]], (*x)[gvl->a[2]], dx);
    rC    = rB+norm(dx);
    mI    = atoms->atom[gvl->a[0]].m;
    mJ    = atoms->atom[gvl->a[1]].m;
    mK    = atoms->atom[gvl->a[2]].m;
    if (gvl->nline == 4)
    {
        rvec_sub((*x)[gvl->a[2]], (*x)[gvl->a[3]], dx);
        rD = rC+norm(dx);
        mL = atoms->atom[gvl->a[3]].m;
    }
    else
    {
        mL = 0;
        rD = 0;
    }
    mT    = mI+mJ+mK+mL;
    /* We need to keep the COM at the same position and the moment of inertia.
     * In order to do this we have two variables, the position of the dummy
     * and the relative masses (we also need to keep the total mass constant).
     * The atom I should be the one connecting to the remainder of the molecule.
     * We put the first atom I at coordinate 0.
     */
    com   = (mJ*rB+mK*rC+mL*rD)/(mT);
    I     = (mI*gmx::square(com) + mJ*gmx::square(rB-com) +
             mK*gmx::square(rC-com) + mL*gmx::square(rD-com));
    rVV   = com+I/(com*mT);
    mV    = com*mT/rVV;
    if (nullptr != debug)
    {
        fprintf(debug, "com = %g, I = %g, rVV = %g mV = %g rB = %g rC = %g rD = %g\n",
                com, I, rVV, mV, rB, rC, rD);
    }
    mI    = (mJ+mK+mL-mV);
    if (mI <= 0)
    {
        gmx_fatal(FARGS, "Zero or negative mass %f in virtual site construction", mI);
    }
    ac[0] = 0;
    ac[1] = (rB/rVV);
    ac[2] = (rC/rVV);
    ac[3] = (rD/rVV);

    natoms = atoms->nr;
    add_t_atoms(atoms, 1, 0);
    srenew(*x, natoms+1);

    /* Add coordinates for mass-particles */
    clear_rvec((*x)[natoms]);
    for (i = 1; (i < 3); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            (*x)[natoms][j] += (atoms->atom[gvl->a[i]].m/mV)*(*x)[gvl->a[i]][j];
        }
    }
    /* Update t_atoms for atoms that change/lose their masses */
    atoms->atom[gvl->a[0]].m  = mI;
    atoms->atom[gvl->a[0]].mB = mI;
    for (i = 1; (i < 3); i++)
    {
        atoms->atom[gvl->a[i]].m  = 0;
        atoms->atom[gvl->a[i]].mB = 0;
    }
    /* Set information in atoms structure for the new mass-particles */
    memset(&aml, 0, sizeof(aml));
    mt = add_atomtype(atype, symtab, &aml, ml, &nbml, -1, 0, 0, 0, 0, 0, 0);
    for (i = natoms; (i <= natoms); i++)
    {
        atoms->atom[i].m          = mV;
        atoms->atom[i].mB         = mV;
        atoms->atom[i].atomnumber = 0;
        if (nullptr != atoms->atomname)
        {
            atoms->atomname[i]        = put_symtab(symtab, ml);
        }
        if (nullptr != atoms->atomtype)
        {
            atoms->atomtype[i]        = put_symtab(symtab, ml);
        }
        if (nullptr != atoms->atomtypeB)
        {
            atoms->atomtypeB[i]       = put_symtab(symtab, ml);
        }
        atoms->atom[i].type       = mt;
        atoms->atom[i].typeB      = mt;
    }

    /* Add constraint between masses */
    memset(&pp, 0, sizeof(pp));
    pp.a[0] = gvl->a[0];
    pp.a[1] = natoms;
    pp.c[0] = rVV;
    add_param_to_plist(plist, F_CONSTR, eitCONSTR, pp);

    /* Add vsites */
    for (i = 1; (i < gvl->nline); i++)
    {
        memset(&pp, 0, sizeof(pp));
        pp.a[0] = gvl->a[i];
        pp.a[1] = gvl->a[0];
        pp.a[2] = natoms;
        pp.c[0] = ac[i];
        add_param_to_plist(plist, F_VSITE2, eitVSITE2, pp);
    }
}

void GentopVsites::mergeLinear(bool bGenVsites)
{
    //int k, l, ai, aj, ndbl;

    for (size_t i = 0; i < linear_.size(); i++)
    {
        linear_[i].nline = 3;
    }


    if (!bGenVsites)
    {
        return;
    }

    /*
    for (unsigned int i = 0; (i < linear_.size()); i++)
    {
        for (unsigned int j = i+1; (j < linear_.size()); j++)
        {
            ndbl = 0;
            for (k = 0; (k < linear_[i].nline); k++)
            {
                ai = linear_[i].a[k];
                for (l = 0; (l < linear_[j].nline); l++)
                {
                    aj = linear_[j].a[l];
                    if (ai == aj)
                    {
                        ndbl++;
                    }
                }
            }

            ndbl /= 2;
            if (ndbl > 0)
            {
                fprintf(stderr, "WARNING: merging two linear vsites into one. Please check result.\n");
                if (nullptr != debug)
                {
                    fprintf(debug, "Linear group j");
                    for (l = 0; (l < linear_[j].nline); l++)
                    {
                        fprintf(debug, " %d", linear_[j].a[l]);
                    }
                    fprintf(debug, "\n");
                    fprintf(debug, "Linear group i");
                    for (k = 0; (k < linear_[i].nline); k++)
                    {
                        fprintf(debug, " %d", linear_[i].a[k]);
                    }
                    fprintf(debug, "\n");
                }
                if ((linear_[j].a[0] == linear_[i].a[1]) &&
                    (linear_[j].a[1] == linear_[i].a[2]))
                {
                    linear_[i].a[3]  = linear_[j].a[2];
                    linear_[i].nline = 4;
                    linear_[j].nline = 0;
                }
                else if ((linear_[i].a[0] == linear_[j].a[1]) &&
                         (linear_[i].a[1] == linear_[j].a[2]))
                {
                    linear_[j].a[3]  = linear_[i].a[2];
                    linear_[j].nline = 4;
                    linear_[i].nline = 0;
                }
                else
                {
                    gmx_fatal(FARGS, "Atoms in strange order in linear vsites. Check debug file.");
                }
                if (nullptr != debug)
                {
                    fprintf(debug, "Linear group j");
                    for (l = 0; (l < linear_[j].nline); l++)
                    {
                        fprintf(debug, " %d", linear_[j].a[l]);
                    }
                    fprintf(debug, "\n");
                    fprintf(debug, "Linear group i");
                    for (k = 0; (k < linear_[i].nline); k++)
                    {
                        fprintf(debug, " %d", linear_[i].a[k]);
                    }
                    fprintf(debug, "\n");
                }
            }
        }
    }
    */
}

static void set_linear_angle_params(const int                  atoms[],
                                    std::vector<PlistWrapper> &plist)
{
    t_param pp;
    bool    found = false;

    auto    pangle = SearchPlist(plist, eitANGLES);

    if (plist.end() == pangle)
    {
        fprintf(stderr, "Cannot find either the angles in the plist to set the linear angle params.");
        return;
    }

    for (auto ang = pangle->beginParam(); ang < pangle->endParam(); ++ang)
    {
        if (((ang->a[0] == atoms[0]) && (ang->a[2] == atoms[2])) ||
            ((ang->a[2] == atoms[0]) && (ang->a[0] == atoms[2])))
        {
            delete_params(plist, pangle->getFtype(), atoms);
            memset(&pp, 0, sizeof(pp));
            for (int i = 0; (i < 3); i++)
            {
                pp.a[i] = atoms[i];
            }
            add_param_to_plist(plist, F_LINEAR_ANGLES, eitLINEAR_ANGLES, pp);
            found = true;
            break;
        }
    }

    if (!found && nullptr != debug)
    {
        fprintf(debug, "Cannot find bond types in linear angle %d-%d-%d\n",
                atoms[0], atoms[1], atoms[2]);
    }
}

void GentopVsites::gen_Vsites(const Poldata             &pd,
                              t_atoms                   *atoms,
                              std::vector<PlistWrapper> &plist,
                              gpp_atomtype              *atype,
                              t_symtab                  *symtab,
                              t_excls                   **excls,
                              t_state                    *state)
{
    int                      nvsite = 0;
    size_t                   ntrain = 0;
    double                   sigma  = 0;       
    double                   akjl   = 0;
    double                   bij    = 0;
    double                   bjk    = 0;
    double                   bjl    = 0;
    double                   aijk   = 0;
    double                   aijl   = 0;
    
    std::string              j, k, l, m;
    std::string              params;    
    std::vector<int>         renum;
    std::vector<int>         inv_renum;
    std::map<int , int>      nuclei;    // std::map<original atom index, number of vsites>
    std::vector<std::string> Akjl;
    std::vector<std::string> Bjk, Bjl;
        
    char                     buf[32];
    char                   **newname;
    t_atoms                 *newatoms;
    t_excls                 *newexcls;
    PaddedRVecVector         newx;    
    t_param                  vs;
    
    memset(&vs, 0, sizeof(vs));    
    int nParticles = atoms->nr + nVsites();
    state_change_natoms(state, nParticles);
    renum.resize(atoms->nr + 1, 0);
    inv_renum.resize(nParticles, -1);
    
    /*Renumber the atoms*/
    for (int i = 0; i < atoms->nr; i++)
    {
        renum[i] = i+nvsite;
        inv_renum[i+nvsite] = i;
        const auto atype(*atoms->atomtype[i]);
        auto vsite = pd.findVsite(atype);
        if (vsite != pd.getVsiteEnd())
        {
            nvsite += vsite->nvsite();
            nuclei.insert(std::pair<int,int>(i, vsite->nvsite()));            
        }
    }
    renum[atoms->nr] = nParticles;
    
    /*Add the virtual sites to the plist*/
    for (int i = 0; i < atoms->nr; i++)
    {
        const auto atype(*atoms->atomtype[i]);
        auto vsite = pd.findVsite(atype);
        if (vsite != pd.getVsiteEnd())
        {
            auto lengthUnit = string2unit(pd.getVsite_length_unit().c_str());
            if (-1 == lengthUnit)
            {
                gmx_fatal(FARGS, "No such length unit '%s'", pd.getVsite_length_unit().c_str());
            }
            if(vsite->type() == evtIN_PLANE)
            {
                auto inplane = findInPlane(vsite->nvsite(), i);
                if (inplane != inplaneEnd())
                {
                    bij        = convert2gmx(vsite->distance(), lengthUnit);
                    aijk       = vsite->angle();
                    aijl       = aijk; // Must be in Degree here, will be converetd to RAD inside GROMACS.                 
                    for (int n = 1; n <= vsite->nvsite(); n++)
                    {
                        vs.a[0] = renum[inplane->ca()] + n; /*vsite    i   */
                        vs.a[1] = renum[inplane->ca()];     /*nucleus  j   */
                        vs.a[2] = renum[inplane->bca()];    /*         m   */
                        vs.a[3] = renum[inplane->bbca1()];  /*        / \  */                         
                        vs.c[0] = aijk;                     /*       k   l */
                        vs.c[1] = bij;  
                        if (n == vsite->nvsite())
                        {
                            vs.a[3] = renum[inplane->bbca2()];                            
                            vs.c[0] = aijl;
                        }
                        add_param_to_plist(plist, F_VSITE3FAD, eitVSITE3FAD, vs);
                    }                                                
                }
            }
            else if (vsite->type() == evtOUT_OF_PLANE)
            {
                auto outplane = findOutPlane(vsite->nvsite(), i);
                if (outplane != outplaneEnd())
                {
                    if (pd.atypeToBtype(*atoms->atomtype[outplane->bca1()], k) &&
                        pd.atypeToBtype(*atoms->atomtype[outplane->ca()],   j) &&
                        pd.atypeToBtype(*atoms->atomtype[outplane->bca2()], l))
                    {
                        bij      = convert2gmx(vsite->distance(), lengthUnit);
                        aijk     = 
                            aijl = DEG2RAD*vsite->angle(); // Must be in RAD here.
                        Akjl     = {k, j, l};
                        Bjk      = {j, k};
                        Bjl      = {j, l};
                        if (pd.searchForce(Akjl, params, &akjl, &sigma, &ntrain, eitANGLES) &&
                            pd.searchForce(Bjk,  params, &bjk,  &sigma, &ntrain, eitBONDS)  &&
                            pd.searchForce(Bjl,  params, &bjl,  &sigma, &ntrain, eitBONDS))
                        {   
                            bjk  = convert2gmx(bjk, lengthUnit);
                            bjl  = convert2gmx(bjl, lengthUnit);
                            akjl = DEG2RAD*akjl;
                            
                            auto pijk = std::cos(aijk)*bij;
                            auto pijl = std::cos(aijl)*bij;
                            auto a    = ( pijk + (pijk*std::cos(akjl)-pijl) * std::cos(akjl) / gmx::square(std::sin(akjl)) ) / bjk;
                            auto b    = ( pijl + (pijl*std::cos(akjl)-pijk) * std::cos(akjl) / gmx::square(std::sin(akjl)) ) / bjl;
                            auto c    = std::sqrt( gmx::square(bij) -
                                                   ( gmx::square(pijk) - 2*pijk*pijl*std::cos(akjl) + gmx::square(pijl) )
                                                   / gmx::square(std::sin(akjl)) ) / ( bjk*bjl*std::sin(akjl) );
                            
                            for (int n = 1; n <= vsite->nvsite(); n++)
                            {
                                vs.a[0] = renum[outplane->ca()] + n; /* vsite     i   */
                                vs.a[1] = renum[outplane->ca()];     /* nucleus   j   */
                                vs.a[2] = renum[outplane->bca1()];   /*          / \  */
                                vs.a[3] = renum[outplane->bca2()];   /*         k   l */            
                                vs.c[0] = a;
                                vs.c[1] = b;  
                                vs.c[2] = c;
                                if (n == vsite->nvsite())
                                {
                                    vs.c[2] = (-1)*c;
                                }
                                add_param_to_plist(plist, F_VSITE3OUT, eitVSITE3OUT, vs);
                            }                            
                        }
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "Undefined virtual site type!\n");
            }
        }        
    }
        
    if (nvsite == nVsites())
    {
        t_atom          *vsite_atom;
        snew(vsite_atom, 1);
        vsite_atom->ptype = eptVSite;

        /* Make new atoms and x arrays */
        snew(newatoms, 1);
        init_t_atoms(newatoms, nParticles, true);
        snew(newatoms->atomtype, nParticles);
        snew(newatoms->atomtypeB, nParticles);
        newatoms->nres = atoms->nres;
        newx.resize(newatoms->nr);
        snew(newname, newatoms->nr);
        
        /* Make a new exclusion array and put the vsites in it */
        snew(newexcls, newatoms->nr);
        
        /* Add exclusion for F_VSITE3OUT virtual type */
        auto pl1 = SearchPlist(plist, F_VSITE3OUT);
        if (plist.end() != pl1)
        {
            // Exclude vsite and nucleus from each other.
            for (auto j = pl1->beginParam(); j < pl1->endParam(); ++j)
            {
                add_excl_pair(newexcls, j->a[0], j->a[1]);
            }
            // Make a copy of the exclusions of the nucleus for the vsite.
            for (auto j = pl1->beginParam(); j < pl1->endParam(); ++j)
            {                
                // We know that the nuclues is 1 as we added it to plist as such.
                int  i0 = inv_renum[j->a[1]];
                for (auto j0 = 0; j0 < (*excls)[i0].nr; j0++)
                {
                    add_excl_pair(newexcls, j->a[1], renum[(*excls)[i0].e[j0]]);
                    add_excl_pair(newexcls, j->a[0], renum[(*excls)[i0].e[j0]]);
                }
            }
            for (auto j = pl1->beginParam(); j < pl1->endParam(); ++j)
            {
                for (auto j0 = 0; j0 < newexcls[j->a[1]].nr; j0++)
                {
                    add_excl_pair(newexcls, j->a[0], newexcls[j->a[1]].e[j0]);
                }
            }
        }
        /* Add exclusion for F_VSITE3FAD virtual type */
        auto pl2 = SearchPlist(plist, F_VSITE3FAD);
        if (plist.end() != pl2)
        {
            for (auto j = pl2->beginParam(); j < pl2->endParam(); ++j)
            {
                add_excl_pair(newexcls, j->a[0], j->a[1]);
            }
            for (auto j = pl2->beginParam(); j < pl2->endParam(); ++j)
            {
                int  i0 = inv_renum[j->a[1]];
                char buf[256];
                snprintf(buf, sizeof(buf), "Uninitialized inv_renum entry for atom %d (%d) vsite %d (%d)",
                         j->a[1], inv_renum[j->a[1]],
                         j->a[0], inv_renum[j->a[0]]);
                GMX_RELEASE_ASSERT(i0 >= 0, buf);
                for (auto j0 = 0; j0 < (*excls)[i0].nr; j0++)
                {
                    add_excl_pair(newexcls, j->a[1], renum[(*excls)[i0].e[j0]]);
                    add_excl_pair(newexcls, j->a[0], renum[(*excls)[i0].e[j0]]);
                }
            }
            for (auto j = pl2->beginParam(); j < pl2->endParam(); ++j)
            {
                for (auto j0 = 0; j0 < newexcls[j->a[1]].nr; j0++)
                {
                    add_excl_pair(newexcls, j->a[0], newexcls[j->a[1]].e[j0]);
                }
            }
        }
        
        /* Set the particle type for atoms having vsites to eptNucleus. */
        for (const auto& nucleus : nuclei)
        {
            if (nucleus.second > 0)
            {
                atoms->atom[nucleus.first].ptype = eptNucleus;
            }
            else
            {
                gmx_fatal(FARGS, "Nucleus %d is supposed to have virtual site!\n", nucleus.first);
            }
        }
        
        /* Copy the old atoms to the new structures. */
        for (int i = 0; i < atoms->nr; i++)
        {
            newatoms->atom[renum[i]]      = atoms->atom[i];
            newatoms->atomname[renum[i]]  = put_symtab(symtab, *atoms->atomname[i]);
            newatoms->atomtype[renum[i]]  = put_symtab(symtab, *atoms->atomtype[i]);
            newatoms->atomtypeB[renum[i]] = put_symtab(symtab, *atoms->atomtypeB[i]);
            newname[renum[i]]             = *atoms->atomtype[i];            
            copy_rvec(state->x[i], newx[renum[i]]);
            t_atoms_set_resinfo(newatoms, renum[i], symtab, 
                                *atoms->resinfo[atoms->atom[i].resind].name,
                                atoms->atom[i].resind, ' ', 1, ' ');
        }
        
        /* Now insert the virtual particles. */
        for (int i = 0; i < atoms->nr; i++)
        {
            auto nucleus = nuclei.find(i);
            if (nucleus != nuclei.end())
            {
                auto iat    = renum[nucleus->first];
                auto nvsite = nucleus->second;
                for(int j = 1; j <= nvsite; j++)
                {
                    newatoms->atom[iat + j]               = atoms->atom[i];
                    newatoms->atom[iat + j].m             = 0;
                    newatoms->atom[iat + j].mB            = 0;
                    newatoms->atom[iat + j].q             = 0;
                    newatoms->atom[iat + j].qB            = 0;
                    newatoms->atom[iat + j].zetaA         = 0;
                    newatoms->atom[iat + j].zetaB         = 0;
                    newatoms->atom[iat + j].atomnumber    = atoms->atom[i].atomnumber;
                    sprintf(buf, "%sL%d", get_atomtype_name(atoms->atom[i].type, atype), j);
                    newname[iat + j] = strdup(buf);
                    auto vsite                            = add_atomtype(atype, symtab, vsite_atom, buf, &vs, 0, 0, 0, 0, 0, 0, 0);
                    newatoms->atom[iat + j].type          = vsite;
                    newatoms->atom[iat + j].typeB         = vsite;
                    newatoms->atomtype[iat + j]           = put_symtab(symtab, buf);
                    newatoms->atomtypeB[iat + j]          = put_symtab(symtab, buf);
                    newatoms->atom[iat + j].ptype         = eptVSite;
                    newatoms->atom[iat + j].resind        = atoms->atom[i].resind;
                    sprintf(buf, "%sL%d", *(atoms->atomname[i]), j);
                    newatoms->atomname[iat + j] = put_symtab(symtab, buf);
                    copy_rvec(state->x[i], newx[iat + j]);
                }
            }
        }
        
        /* Copy newatoms to atoms */
        copy_atoms(newatoms, atoms);
        
        /* Copy coordinates and names */
        for (int i = 0; i < newatoms->nr; i++)
        {
            copy_rvec(newx[i], state->x[i]);
            atoms->atomtype[i] = put_symtab(symtab, newname[i]);
        }
        sfree(newname);
        
        /* Copy exclusions, empty the original first */
        sfree((*excls));
        (*excls) = newexcls;
                
        /*Now renumber atoms in all other plist interaction types */
        for (auto pw = plist.begin(); pw < plist.end(); ++pw)
        {
            if (pw->getFtype() != F_VSITE3FAD && pw->getFtype() != F_VSITE3OUT)
            {
                for (auto j = pw->beginParam(); j < pw->endParam(); ++j)
                {
                    for (int k = 0; k < NRAL(pw->getFtype()); k++)
                    {
                        j->a[k] = renum[j->a[k]];
                    }
                }
            }
        }
        sfree(vsite_atom);
    }
    else
    {
        gmx_fatal(FARGS, 
                  "Number of vistes counted here does not match the number of vsites needed: %d vs %d",
                  nvsite, nVsites());
    }
}

void GentopVsites::generateSpecial(const Poldata              &pd,
                                   bool                        bUseVsites,
                                   t_atoms                    *atoms,
                                   rvec                      **x,
                                   std::vector<PlistWrapper>  &plist,
                                   t_symtab                   *symtab,
                                   gpp_atomtype_t              atype,
                                   t_excls                   **excls,
                                   t_state                    *state)
{
    int          j, nlin_at;
    int          a[MAXATOMLIST], aa[2];
    t_param      pp;
    unsigned int ftb, fta, ftl, ftp, fti;

    mergeLinear(bUseVsites);

    ftb     = pd.findForces(eitBONDS)->fType();
    fta     = pd.findForces(eitANGLES)->fType();
    ftl     = pd.findForces(eitLINEAR_ANGLES)->fType();
    ftp     = pd.findForces(eitPROPER_DIHEDRALS)->fType();
    fti     = pd.findForces(eitIMPROPER_DIHEDRALS)->fType();

    nlin_at = 0;
    for (unsigned int i = 0; (i < linear_.size()); i++)
    {
        nlin_at += linear_[i].nline;
    }

    if (nullptr != debug)
    {
        fprintf(debug, "Generating %d linear ", nlin_at);
        if (bUseVsites)
        {
            fprintf(debug, "vsites");
        }
        else
        {
            fprintf(debug, "angles");
        }
        fprintf(debug, " and %d impropers\n", (int)planar_.size());
    }
    if ((gvt_ == evtLINEAR) || (gvt_ == evtALL))
    {
        /* If we use vsites (discouraged) each triplet of atoms in a linear arrangement
         * is described by
         * two dummy masses connected by a constraint, that maintain the
         * moment of inertia. The three atoms are then generated from these two
         * positions using a virtual_site2 construction. We need to add 2 extra
         * particles for each linear group.
         * In case we use the special linear angle terms, life gets a lot easier!
         */
        
        for (size_t i = 0; i < linear_.size(); i++)
        {
            for (j = 0; (j < linear_[i].nline); j++)
            {
                a[j] = linear_[i].a[j];
            }
            for (; (j < MAXATOMLIST); j++)
            {
                a[j] = -1;
            }
            delete_params(plist, fta, a);
            delete_params(plist, ftl, a);
            delete_params(plist, ftp, a);
            delete_params(plist, fti, a);
            
            //if (bUseVsites) adding vsite for linears need to be fixed
            if (0)
            {                
                /* Complicated algorithm, watch out */
                for (j = 0; (j < linear_[i].nline-1); j++)
                {
                    aa[0] = a[j+0];
                    aa[1] = a[j+1];
                    delete_params(plist, ftb, aa);
                }

                /* Compute details for the new masses and vsites,
                 * and update everything
                 */
                calc_vsite2parm(atoms, plist, x, &linear_[i], symtab, atype);
                srenew((*excls), atoms->nr);
                for (j = atoms->nr-2; (j <= atoms->nr-1); j++)
                {
                    (*excls)[j].nr = 1;
                    snew((*excls)[j].e, 1);
                    (*excls)[j].e[0] = j;
                }                
            }
            else
            {
                set_linear_angle_params(a, plist);
            }
        }
    }
    if ((gvt_ == evtPLANAR) || (gvt_ == evtALL))
    {
        for (unsigned int i = 0; (i < planar_.size()); i++)
        {
            /* First delete superfluous dihedrals */
            a[1] = planar_[i].a[0];
            a[3] = -1;
            for (j = 1; (j < 4); j++)
            {
                if (planar_[i].nb[j] == 1)
                {
                    if (j == 1)
                    {
                        a[0] = planar_[i].a[1];
                        a[2] = planar_[i].a[2];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                        a[2] = planar_[i].a[3];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                    }
                    else if (j == 2)
                    {
                        a[0] = planar_[i].a[2];
                        a[2] = planar_[i].a[1];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                        a[2] = planar_[i].a[3];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                    }
                    else if (j == 3)
                    {
                        a[0] = planar_[i].a[3];
                        a[2] = planar_[i].a[1];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                        a[2] = planar_[i].a[2];
                        delete_params(plist, ftp, a);
                        delete_params(plist, fti, a);
                    }
                }
            }
            /* Now add impropers! */
            memset(&pp, 0, sizeof(pp));
            for (j = 0; (j < 4); j++)
            {
                pp.a[j] = planar_[i].a[j];
            }
            add_param_to_plist(plist, F_IDIHS, eitIMPROPER_DIHEDRALS, pp);
        }
    }
    if (bUseVsites && (inplane_.size() > 0 || outplane_.size() > 0))
    {
        gen_Vsites(pd, atoms, plist, atype, symtab, excls, state);
    }
}

}
