/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \file
 * \brief
 * Defines interaction functions.
 *
 * \ingroup module_topology
 * \inlibraryapi
 */
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "gromacs/topology/ifunc.h"

/*! \brief Interaction function definitions.
 *
 * This MUST correspond to the enum in api/legacy/include/gromacs/topology/ifunc.h.
 *
 * Note also that the longname field of the interaction is used for
 * printing e.g. the mdrun log file in a columnar style, and pr_ebin
 * makes available only 14 printing characters (ie not including the
 * terminating '\0'). So please abbreviate accordingly,
 * e.g. "Conserved En."
 */
// clang-format off
constexpr gmx::EnumerationArray<InteractionFunction, t_interaction_function> interaction_function = {
//    Name                Long name         Number of atoms
//                                             Number of A-state parameters
//                                                Number of B-state parameters
//                                                   Interaction function categories
    {{"BONDS",            "Bond",           2, 2, 2, IF_BOND | IF_CHEMBOND | IF_BTYPE },
     {"G96BONDS",         "G96Bond",        2, 2, 2, IF_BOND | IF_CHEMBOND | IF_BTYPE },
     {"MORSE",            "Morse",          2, 3, 3, IF_BOND | IF_CHEMBOND | IF_BTYPE },
     {"CUBICBONDS",       "Cubic Bonds",    2, 3, 0, IF_BOND | IF_CHEMBOND | IF_BTYPE },
     {"CONNBONDS",        "Connect Bonds",  2, 0, 0, IF_BOND | IF_CHEMBOND },
     {"HARMONIC",         "Harmonic Pot.",  2, 2, 2, IF_BOND },
     {"FENEBONDS",        "FENE Bonds",     2, 2, 0, IF_BOND | IF_CHEMBOND },
     {"TABBONDS",         "Tab. Bonds",     2, 2, 2, IF_BOND | IF_CHEMBOND | IF_BTYPE | IF_TABULATED },
     {"TABBONDSNC",       "Tab. Bonds NC",  2, 2, 2, IF_BOND | IF_TABULATED | IF_LIMZERO },
     {"RESTRAINTPOT",     "Restraint Pot.", 2, 4, 4, IF_BOND },
     {"ANGLES",           "Angle",          3, 2, 2, IF_BOND | IF_ATYPE },
     {"G96ANGLES",        "G96Angle",       3, 2, 2, IF_BOND | IF_ATYPE },
     {"RESTRANGLES",      "Restr. Angles",  3, 2, 2, IF_BOND | IF_ATYPE },
     {"LINEAR_ANGLES",    "Lin. Angle",     3, 2, 2, IF_BOND | IF_ATYPE },
     {"CROSS_BOND_BOND",  "Bond-Cross",     3, 3, 0, IF_BOND },
     {"CROSS_BOND_ANGLE", "BA-Cross",       3, 4, 0, IF_BOND },
     {"UREY_BRADLEY",     "U-B",            3, 4, 4, IF_BOND | IF_ATYPE },
     {"QANGLES",          "Quartic Angles", 3, 6, 0, IF_BOND | IF_ATYPE },
     // TODO should TABANGLES include IF_ATYPE?
     {"TABANGLES",        "Tab. Angles",    3, 2, 2, IF_BOND | IF_TABULATED },
     {"PDIHS",            "Proper Dih.",    4, 3, 3, IF_BOND | IF_DIHEDRAL },
     {"RBDIHS",           "Ryckaert-Bell.", 4, 6, 6, IF_BOND | IF_DIHEDRAL },
     {"RESTRDIHS",        "Restr. Dih.",    4, 2, 2, IF_BOND | IF_DIHEDRAL },
     {"CBTDIHS",          "CBT Dih.",       4, 6, 6, IF_BOND | IF_DIHEDRAL },
     {"FOURDIHS",         "Fourier Dih.",   4, 4, 4, IF_BOND | IF_DIHEDRAL },
     {"IDIHS",            "Improper Dih.",  4, 2, 2, IF_BOND | IF_DIHEDRAL },
     {"PIDIHS",           "Per. Imp. Dih.", 4, 3, 3, IF_BOND | IF_DIHEDRAL },
     {"TABDIHS",          "Tab. Dih.",      4, 2, 2, IF_BOND | IF_DIHEDRAL | IF_TABULATED },
     {"CMAP",             "CMAP Dih.",      5, 1, 1, IF_BOND | IF_DIHEDRAL },
     {"GB12",             "GB 1-2 Pol.",    0, 0, 0, IF_NULL }, /* unused */
     {"GB13",             "GB 1-3 Pol.",    0, 0, 0, IF_NULL }, /* unused */
     {"GB14",             "GB 1-4 Pol.",    0, 0, 0, IF_NULL }, /* unused */
     {"GBPOL",            "GB Polariz.",    0, 0, 0, IF_NULL }, /* unused */
     {"NPSOLVATION",      "Nonpolar Sol.",  0, 0, 0, IF_NULL }, /* unused */
     {"LJ14",             "LJ-14",          2, 2, 2, IF_BOND | IF_PAIR | IF_LIMZERO },
     {"COUL14",           "Coulomb-14",     0, 0, 0, IF_NULL },
     {"LJC14_Q",          "LJC-14 q",       2, 5, 0, IF_BOND | IF_PAIR | IF_LIMZERO },
     {"LJC_NB",           "LJC Pairs NB",   2, 4, 0, IF_BOND | IF_PAIR | IF_LIMZERO },
     {"LJ_SR",            "LJ (SR)",        2, 2, 0, IF_NULL },
     {"BHAM",             "Buck.ham (SR)",  2, 3, 0, IF_NULL },
     {"LJ_LR",            "LJ",             0, 0, 0, IF_NULL }, /* unused */
     {"BHAM_LR",          "B.ham",          0, 0, 0, IF_NULL }, /* unused */
     {"DISPCORR",         "Disper. corr.",  0, 0, 0, IF_NULL },
     {"COUL_SR",          "Coulomb (SR)",   0, 0, 0, IF_NULL },
     {"COUL_LR",          "Coul",           0, 0, 0, IF_NULL }, /* unused */
     {"RF_EXCL",          "RF excl.",       0, 0, 0, IF_NULL },
     {"COUL_RECIP",       "Coul. recip.",   0, 0, 0, IF_NULL },
     {"LJ_RECIP",         "LJ recip.",      0, 0, 0, IF_NULL },
     {"DPD",              "DPD",            0, 0, 0, IF_NULL },
     {"POLARIZATION",     "Polarization",   2, 1, 0, IF_BOND | IF_CHEMBOND },
     {"WATERPOL",         "Water Pol.",     5, 6, 0, IF_BOND },
     {"THOLE",            "Thole Pol.",     4, 3, 0, IF_BOND },
     {"ANHARM_POL",       "Anharm. Pol.",   2, 3, 0, IF_BOND | IF_CHEMBOND },
     {"POSRES",           "Position Rest.", 1, 3, 3, IF_BOND },
     {"FBPOSRES",         "Flat-b. P-R.",   1, 3, 0, IF_BOND },
     {"DISRES",           "Dis. Rest.",     2, 6, 0, IF_BOND },
     {"DISRESVIOL",       "D.R.Viol. (nm)", 0, 0, 0, IF_NULL },
     {"ORIRES",           "Orient. Rest.",  2, 6, 0, IF_BOND },
     {"ORDEV",            "Ori. R. RMSD",   0, 0, 0, IF_NULL },
     {"ANGRES",           "Angle Rest.",    4, 3, 3, IF_BOND },
     {"ANGRESZ",          "Angle Rest. Z",  2, 3, 3, IF_BOND },
     {"DIHRES",           "Dih. Rest.",     4, 3, 3, IF_BOND },
     {"DIHRESVIOL",       "Dih. Rest. Vi.", 0, 0, 0, IF_NULL }, /* obsolete */
     {"CONSTR",           "Constraint",     2, 1, 1, IF_CONSTRAINT | IF_CHEMBOND },
     {"CONSTRNC",         "Constr. No Co.", 2, 1, 1, IF_CONSTRAINT },
     {"SETTLE",           "Settle",         3, 2, 0, IF_CONSTRAINT | IF_CHEMBOND },
     {"VSITE1",           "Virtual site 1", 2, 0, 0, IF_VSITE },
     {"VSITE2",           "Virtual site 2", 3, 1, 0, IF_VSITE },
     {"VSITE2FD",         "Virt. site 2fd", 3, 1, 0, IF_VSITE },
     {"VSITE3",           "Virtual site 3", 4, 2, 0, IF_VSITE },
     {"VSITE3FD",         "Virt. site 3fd", 4, 2, 0, IF_VSITE },
     {"VSITE3FAD",        "Vir. site 3fad", 4, 2, 0, IF_VSITE },
     {"VSITE3OUT",        "Vir. site 3out", 4, 3, 0, IF_VSITE },
     {"VSITE4FD",         "Virt. site 4fd", 5, 3, 0, IF_VSITE },
     {"VSITE4FDN",        "Vir. site 4fdn", 5, 3, 0, IF_VSITE },
     {"VSITEN",           "Virtual site N", 2, 2, 0, IF_VSITE },
     {"COM_PULL",         "COM Pull En.",   0, 0, 0, IF_NULL },
     {"DENSITYFIT",       "Dens. fitting",  0, 0, 0, IF_NULL },
     {"EQM",              "Quantum En.",    0, 0, 0, IF_NULL },
     {"ENNPOT",           "NN Potential",   0, 0, 0, IF_NULL },
     {"EPOT",             "Potential",      0, 0, 0, IF_NULL },
     {"EKIN",             "Kinetic En.",    0, 0, 0, IF_NULL },
     {"ETOT",             "Total Energy",   0, 0, 0, IF_NULL },
     {"ECONS",            "Conserved En.",  0, 0, 0, IF_NULL },
     {"TEMP",             "Temperature",    0, 0, 0, IF_NULL },
     {"VTEMP",            "Vir. Temp.",     0, 0, 0, IF_NULL }, /* unused */
    /* Note that pressure names can not be more than 8 characters,
     * because " (bar)" is appended to them. */
     {"PDISPCORR",        "Pres. DC",       0, 0, 0, IF_NULL },
     {"PRES",             "Pressure",       0, 0, 0, IF_NULL },
     {"DH/DL_CON",        "dH/dl constr.",  0, 0, 0, IF_NULL }, /* obsolete */
     {"DV/DL",            "dVremain/dl",    0, 0, 0, IF_NULL },
     {"DK/DL",            "dEkin/dl",       0, 0, 0, IF_NULL },
     {"DVC/DL",           "dVcoul/dl",      0, 0, 0, IF_NULL },
     {"DVV/DL",           "dVvdw/dl",       0, 0, 0, IF_NULL },
     {"DVB/DL",           "dVbonded/dl",    0, 0, 0, IF_NULL },
     {"DVR/DL",           "dVrestraint/dl", 0, 0, 0, IF_NULL },
     {"DVT/DL",           "dVtemp/dl",      0, 0, 0, IF_NULL }}
};
// clang-format on
