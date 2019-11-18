/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "ifunc.h"

#define def_bonded(str, lstr, nra, nrpa, nrpb)    \
    {                                             \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND \
    }

#define def_bondedz(str, lstr, nra, nrpa, nrpb)                \
    {                                                          \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_LIMZERO \
    }

#define def_bondedt(str, lstr, nra, nrpa, nrpb)                  \
    {                                                            \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_TABULATED \
    }

#define def_bondedtz(str, lstr, nra, nrpa, nrpb)                              \
    {                                                                         \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_TABULATED | IF_LIMZERO \
    }

#define def_angle(str, lstr, nra, nrpa, nrpb)                \
    {                                                        \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_ATYPE \
    }

#define def_bond(str, lstr, nra, nrpa, nrpb)                               \
    {                                                                      \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_CHEMBOND | IF_BTYPE \
    }

#define def_bondt(str, lstr, nra, nrpa, nrpb)                                  \
    {                                                                          \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_CHEMBOND | IF_TABULATED \
    }

#define def_bondnb(str, lstr, nra, nrpa, nrpb)                  \
    {                                                           \
        str, lstr, (nra), (nrpa), (nrpb), IF_BOND | IF_CHEMBOND \
    }

#define def_vsite(str, lstr, nra, nrpa)       \
    {                                         \
        str, lstr, (nra), (nrpa), 0, IF_VSITE \
    }

#define def_shk(str, lstr, nra, nrpa, nrpb)             \
    {                                                   \
        str, lstr, (nra), (nrpa), (nrpb), IF_CONSTRAINT \
    }

#define def_shkcb(str, lstr, nra, nrpa, nrpb)                         \
    {                                                                 \
        str, lstr, (nra), (nrpa), (nrpb), IF_CONSTRAINT | IF_CHEMBOND \
    }

#define def_nb(str, lstr, nra, nrp)         \
    {                                       \
        str, lstr, (nra), (nrp), 0, IF_NULL \
    }

#define def_nofc(str, lstr)         \
    {                               \
        str, lstr, 0, 0, 0, IF_NULL \
    }

/* this MUST correspond to the enum in src/gromacs/topology/ifunc.h */
const t_interaction_function interaction_function[F_NRE] = {
    def_bond("BONDS", "Bond", 2, 2, 2), def_bond("G96BONDS", "G96Bond", 2, 2, 2),
    def_bond("MORSE", "Morse", 2, 3, 3), def_bond("CUBICBONDS", "Cubic Bonds", 2, 3, 0),
    def_bondnb("CONNBONDS", "Connect Bonds", 2, 0, 0),
    def_bonded("HARMONIC", "Harmonic Pot.", 2, 2, 2), def_bondnb("FENEBONDS", "FENE Bonds", 2, 2, 0),
    def_bondt("TABBONDS", "Tab. Bonds", 2, 2, 2), def_bondedtz("TABBONDSNC", "Tab. Bonds NC", 2, 2, 2),
    def_bonded("RESTRAINTPOT", "Restraint Pot.", 2, 4, 4), def_angle("ANGLES", "Angle", 3, 2, 2),
    def_angle("G96ANGLES", "G96Angle", 3, 2, 2), def_angle("RESTRANGLES", "Restricted Angles", 3, 2, 2),
    def_angle("LINEAR_ANGLES", "Lin. Angle", 3, 2, 2),
    def_bonded("CROSS_BOND_BOND", "Bond-Cross", 3, 3, 0),
    def_bonded("CROSS_BOND_ANGLE", "BA-Cross", 3, 4, 0), def_angle("UREY_BRADLEY", "U-B", 3, 4, 4),
    def_angle("QANGLES", "Quartic Angles", 3, 6, 0), def_bondedt("TABANGLES", "Tab. Angles", 3, 2, 2),
    def_bonded("PDIHS", "Proper Dih.", 4, 3, 3), def_bonded("RBDIHS", "Ryckaert-Bell.", 4, 6, 6),
    def_bonded("RESTRDIHS", "Restricted Dih.", 4, 2, 2), def_bonded("CBTDIHS", "CBT Dih.", 4, 6, 6),
    def_bonded("FOURDIHS", "Fourier Dih.", 4, 4, 4), def_bonded("IDIHS", "Improper Dih.", 4, 2, 2),
    def_bonded("PIDIHS", "Improper Dih.", 4, 3, 3), def_bondedt("TABDIHS", "Tab. Dih.", 4, 2, 2),
    def_bonded("CMAP", "CMAP Dih.", 5, -1, -1), def_nofc("GB12", "GB 1-2 Pol. (unused)"),
    def_nofc("GB13", "GB 1-3 Pol. (unused)"), def_nofc("GB14", "GB 1-4 Pol. (unused)"),
    def_nofc("GBPOL", "GB Polarization (unused)"), def_nofc("NPSOLVATION", "Nonpolar Sol. (unused)"),
    def_bondedz("LJ14", "LJ-14", 2, 2, 2), def_nofc("COUL14", "Coulomb-14"),
    def_bondedz("LJC14_Q", "LJC-14 q", 2, 5, 0), def_bondedz("LJC_NB", "LJC Pairs NB", 2, 4, 0),
    def_nb("LJ_SR", "LJ (SR)", 2, 2), def_nb("BHAM", "Buck.ham (SR)", 2, 3),
    def_nofc("LJ_LR", "LJ (unused)"), def_nofc("BHAM_LR", "B.ham (unused)"),
    def_nofc("DISPCORR", "Disper. corr."), def_nofc("COUL_SR", "Coulomb (SR)"),
    def_nofc("COUL_LR", "Coul (unused)"), def_nofc("RF_EXCL", "RF excl."),
    def_nofc("COUL_RECIP", "Coul. recip."), def_nofc("LJ_RECIP", "LJ recip."),
    def_nofc("DPD", "DPD"), def_bondnb("POLARIZATION", "Polarization", 2, 1, 0),
    def_bonded("WATERPOL", "Water Pol.", 5, 6, 0), def_bonded("THOLE", "Thole Pol.", 4, 3, 0),
    def_bondnb("ANHARM_POL", "Anharm. Pol.", 2, 3, 0), def_bonded("POSRES", "Position Rest.", 1, 3, 3),
    def_bonded("FBPOSRES", "Flat-bottom posres", 1, 3, 0), def_bonded("DISRES", "Dis. Rest.", 2, 6, 0),
    def_nofc("DISRESVIOL", "D.R.Viol. (nm)"), def_bonded("ORIRES", "Orient. Rest.", 2, 6, 0),
    def_nofc("ORDEV", "Ori. R. RMSD"), def_bonded("ANGRES", "Angle Rest.", 4, 3, 3),
    def_bonded("ANGRESZ", "Angle Rest. Z", 2, 3, 3), def_bonded("DIHRES", "Dih. Rest.", 4, 3, 3),
    def_nofc("DIHRESVIOL", "Dih. Rest. Viol."), /* obsolete */
    def_shkcb("CONSTR", "Constraint", 2, 1, 1), def_shk("CONSTRNC", "Constr. No Conn.", 2, 1, 1),
    def_shkcb("SETTLE", "Settle", 3, 2, 0), def_vsite("VSITE2", "Virtual site 2", 3, 1),
    def_vsite("VSITE2FD", "Virtual site 2fd", 3, 1), def_vsite("VSITE3", "Virtual site 3", 4, 2),
    def_vsite("VSITE3FD", "Virtual site 3fd", 4, 2), def_vsite("VSITE3FAD", "Virtual site 3fad", 4, 2),
    def_vsite("VSITE3OUT", "Virtual site 3out", 4, 3), def_vsite("VSITE4FD", "Virtual site 4fd", 5, 3),
    def_vsite("VSITE4FDN", "Virtual site 4fdn", 5, 3), def_vsite("VSITEN", "Virtual site N", 2, 2),
    def_nofc("COM_PULL", "COM Pull En."), def_nofc("DENSITYFIT", "Density fitting"),
    def_nofc("EQM", "Quantum En."), def_nofc("EPOT", "Potential"), def_nofc("EKIN", "Kinetic En."),
    def_nofc("ETOT", "Total Energy"), def_nofc("ECONS", "Conserved En."),
    def_nofc("TEMP", "Temperature"), def_nofc("VTEMP", "Vir. Temp. (not used)"),
    /* Note that pressure names can not be more than 8 char's,
     * because " (bar)" is appended to them.
     */
    def_nofc("PDISPCORR", "Pres. DC"), def_nofc("PRES", "Pressure"),
    def_nofc("DH/DL_CON", "dH/dl constr."), /* obsolete */
    def_nofc("DV/DL", "dVremain/dl"), def_nofc("DK/DL", "dEkin/dl"), def_nofc("DVC/DL", "dVcoul/dl"),
    def_nofc("DVV/DL", "dVvdw/dl"), def_nofc("DVB/DL", "dVbonded/dl"),
    def_nofc("DVR/DL", "dVrestraint/dl"), def_nofc("DVT/DL", "dVtemperature/dl")
};
