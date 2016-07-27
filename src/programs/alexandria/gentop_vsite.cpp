/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#include "gentop_vsite.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "plistwrapper.h"
#include "poldata.h"

namespace alexandria
{

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

void GentopVsites::addOutOfPlane(int ai, int aj, int ak, int al,
                                 int nbonds[])
{
    unsigned int i;

    for (i = 0; (i < outofplane_.size()); i++)
    {
        if (((outofplane_[i].a[0] == ai) && (outofplane_[i].a[1] == aj)  &&
             (outofplane_[i].a[2] == ak) && (outofplane_[i].a[3] == al)) ||
            ((outofplane_[i].a[0] == al) && (outofplane_[i].a[1] == ak)  &&
             (outofplane_[i].a[2] == aj) && (outofplane_[i].a[3] == ai)))
        {
            break;
        }
    }
    if (i == outofplane_.size())
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
        outofplane_.push_back(p);
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

static void calc_vsite2parm(t_atoms *atoms,
                            std::vector<PlistWrapper> &plist,
                            rvec **x,
                            gv_linear *gvl, t_symtab *symtab,
                            gpp_atomtype_t atype)
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
    if (NULL != debug)
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
        if (NULL != atoms->atomname)
        {
            atoms->atomname[i]        = put_symtab(symtab, ml);
        }
        if (NULL != atoms->atomtype)
        {
            atoms->atomtype[i]        = put_symtab(symtab, ml);
        }
        if (NULL != atoms->atomtypeB)
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
    int k, l, ai, aj, ndbl;

    for (unsigned int i = 0; (i < linear_.size()); i++)
    {
        linear_[i].nline = 3;
    }

    if (!bGenVsites)
    {
        return;
    }

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
                if (NULL != debug)
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
                if (NULL != debug)
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
}

static void set_linear_angle_params(const int                  atoms[],
                                    std::vector<PlistWrapper> &plist,
                                    real                       klin)
{
    t_param pp;
    real    b0 = 0;
    real    b1 = 0;

    auto    pangle = SearchPlist(plist, eitANGLES);
    auto    pbond  = SearchPlist(plist, eitBONDS);

    if (plist.end() == pangle || plist.end() == pbond)
    {
        fprintf(stderr, "Cannot find either the angles or the bonds in the plist to set the linear angle params.");
        return;
    }

    for (auto ang = pangle->beginParam(); ang < pangle->endParam(); ++ang)
    {
        if (((ang->a[0] == atoms[0]) && (ang->a[2] == atoms[2])) ||
            ((ang->a[2] == atoms[0]) && (ang->a[0] == atoms[2])))
        {
	    for (auto b = pbond->beginParam(); b < pbond->endParam(); ++b)
	    {
	        if (((b->a[0] == atoms[0]) && (b->a[1] == atoms[1])) ||
		    ((b->a[0] == atoms[1]) && (b->a[1] == atoms[0])))
		{
		    b0 = b->c[0];
		}
		else if (((b->a[0] == atoms[2]) && (b->a[1] == atoms[1])) ||
			 ((b->a[0] == atoms[1]) && (b->a[1] == atoms[2])))
		{
		    b1 = b->c[0];
		}
	    }
        }
    }

    if (b0 > 0 && b1 > 0)
    {
        delete_params(plist, pangle->getFtype(), atoms);
        memset(&pp, 0, sizeof(pp));
        for (int i = 0; (i < 3); i++)
        {
            pp.a[i] = atoms[i];
        }
        pp.c[0] = klin;
	pp.c[1] = (b1/(b1+b0));
        add_param_to_plist(plist, F_LINEAR_ANGLES, eitLINEAR_ANGLES, pp);
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Cannot find bonds in linear angle %d-%d-%d b0: %g and b1: %g\n",
                atoms[0], atoms[1], atoms[2], b0, b1);
    }
}

void GentopVsites::generateSpecial(const Poldata              &pd,
                                   bool                        bUseVsites,
                                   t_atoms                    *atoms,
                                   rvec                      **x,
                                   std::vector<PlistWrapper>  &plist,
                                   t_symtab                   *symtab,
                                   gpp_atomtype_t              atype,
                                   t_excls                   **excls)
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

    if (NULL != debug)
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
    if ((gvt_ == egvtLINEAR) || (gvt_ == egvtALL))
    {
        /* If we use vsites (discouraged) each triplet of atoms in a linear arrangement
         * is described by
         * two dummy masses connected by a constraint, that maintain the
         * moment of inertia. The three atoms are then generated from these two
         * positions using a virtual_site2 construction. We need to add 2 extra
         * particles for each linear group.
         * In case we use the special linear angle terms, life gets a lot easier!
         */
        for (unsigned int i = 0; (i < linear_.size()); i++)
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

            if (bUseVsites)
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
                calc_vsite2parm(atoms, plist,
                                x, &linear_[i], symtab, atype);
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
                /* Need to add parameters here! */
                real kth = 400;
                set_linear_angle_params(a, plist, kth);
            }
        }
    }
    if ((gvt_ == egvtPLANAR) || (gvt_ == egvtALL))
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
}

}
