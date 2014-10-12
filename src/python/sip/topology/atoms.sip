/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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

struct t_atom {

%TypeHeaderCode
#include <gromacs/topology/atoms.h>
%End

    double         m; double q;        /* Mass and charge                      */
    double         mB; double qB;      /* Mass and charge for Free Energy calc */
    unsigned short type;        /* Atom type                            */
    unsigned short typeB;       /* Atom type for Free Energy calc       */
    int            ptype;       /* Particle type                        */
    int            resind;      /* Index into resinfo (in t_atoms)      */
    int            atomnumber;  /* Atomic Number or NOTSET              */
    char           *elem /NoSetter/;       /* Element name                         */
};

struct t_resinfo {

%TypeHeaderCode
#include <gromacs/topology/atoms.h>
%End

    char          **name {      /* Pointer to the residue name          */
    %GetCode
        sipPy = sipBuildResult(NULL, "s", *(sipCpp->name));
    %End
    %SetCode

    %End
    };
    int             nr;         /* Residue number                       */
    unsigned char   ic;         /* Code for insertion of residues       */
    int             chainnum;   /* Iincremented at TER or new chain id  */
    char            chainid;    /* Chain identifier written/read to pdb */
    char          **rtp {       /* rtp building block name (optional)   */
    %GetCode
        sipPy = sipBuildResult(NULL, "s", *(sipCpp->rtp));
    %End
    %SetCode

    %End
    };
};

struct t_pdbinfo {

%TypeHeaderCode
#include <gromacs/topology/atoms.h>
%End

    int      type;              /* PDB record name                      */
    int      atomnr;            /* PDB atom number                      */
    char     altloc;            /* Alternate location indicator         */
    char    *atomnm /NoSetter/; /* True atom name including leading spaces */
    double   occup;             /* Occupancy                            */
    double   bfac;              /* B-factor                             */
    bool bAnisotropic;          /* (an)isotropic switch                 */
    int     *uij {
    %GetCode
        sipPy = sipConvertToArray(sipCpp->uij, "i", 6, SIP_READ_ONLY);
    %End
    %SetCode

    %End
    };
};

struct t_grps {

%TypeHeaderCode
#include <gromacs/topology/atoms.h>
%End

    int   nr;                   /* Number of different groups           */
    int  *nm_ind {
    %GetCode
        sipPy = sipConvertToArray(sipCpp->nm_ind, "i", sipCpp->nr, SIP_READ_ONLY);
    %End
    %SetCode

    %End
    };
};

struct t_atoms {

%TypeHeaderCode
#include <gromacs/topology/atoms.h>
// Sip creates 'copy_[typename]' functions. In this case it coincides with function from atoms.h, so we need to rename it.
#define copy_t_atoms copy_t_atoms_sip
%End

    int            nr;          /* Nr of atoms                          */
    t_atom        *atom {
    %GetCode
        sipPy = sipConvertToTypedArray(sipCpp->atom, sipType_t_atom, NULL, sizeof(t_atom), sipCpp->nr, SIP_READ_ONLY);
    %End
    %SetCode

    %End
    };
                                /*[> The following entries will not       <]
                                [> always be used (nres==0)             <]
    char          ***atomname;  [> Array of pointers to atom name       <]
                                [> use: (*(atomname[i]))                <]
    char          ***atomtype;  [> Array of pointers to atom types      <]
                                [> use: (*(atomtype[i]))                <]
    char          ***atomtypeB; [> Array of pointers to B atom types    <]
                                [> use: (*(atomtypeB[i]))               <]*/
    int              nres;
    t_resinfo       *resinfo {
    %GetCode
        sipPy = sipConvertToTypedArray(sipCpp->resinfo, sipType_t_resinfo, NULL, sizeof(t_resinfo), sipCpp->nres, SIP_READ_ONLY);
    %End
    %SetCode

    %End
    };
    // FIXME: how many pdbinfo's are there?
    t_pdbinfo       *pdbinfo {
    %GetCode
        sipPy = sipConvertToTypedArray(sipCpp->pdbinfo, sipType_t_pdbinfo, NULL, sizeof(t_pdbinfo), sipCpp->nres, SIP_READ_ONLY);
    %End
    %SetCode

    %End
    };
};
