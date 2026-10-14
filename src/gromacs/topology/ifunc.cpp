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
    {{"BONDS",            "Bond",           2, 2, 2, sc_ifcBondedForce | sc_ifcChemicalBond | sc_ifcBType },
     {"G96BONDS",         "G96Bond",        2, 2, 2, sc_ifcBondedForce | sc_ifcChemicalBond | sc_ifcBType },
     {"MORSE",            "Morse",          2, 3, 3, sc_ifcBondedForce | sc_ifcChemicalBond | sc_ifcBType },
     {"CUBICBONDS",       "Cubic Bonds",    2, 3, 0, sc_ifcBondedForce | sc_ifcChemicalBond | sc_ifcBType },
     {"CONNBONDS",        "Connect Bonds",  2, 0, 0, sc_ifcBondedForce | sc_ifcChemicalBond },
     {"HARMONIC",         "Harmonic Pot.",  2, 2, 2, sc_ifcBondedForce },
     {"FENEBONDS",        "FENE Bonds",     2, 2, 0, sc_ifcBondedForce | sc_ifcChemicalBond },
     {"TABBONDS",         "Tab. Bonds",     2, 2, 2, sc_ifcBondedForce | sc_ifcChemicalBond | sc_ifcBType | sc_ifcTabulated },
     {"TABBONDSNC",       "Tab. Bonds NC",  2, 2, 2, sc_ifcBondedForce | sc_ifcTabulated | sc_ifcZeroAtLargeDistances },
     {"RESTRAINTPOT",     "Restraint Pot.", 2, 4, 4, sc_ifcBondedForce },
     {"ANGLES",           "Angle",          3, 2, 2, sc_ifcBondedForce | sc_ifcAngle },
     {"G96ANGLES",        "G96Angle",       3, 2, 2, sc_ifcBondedForce | sc_ifcAngle },
     {"RESTRANGLES",      "Restr. Angles",  3, 2, 2, sc_ifcBondedForce | sc_ifcAngle },
     {"LINEAR_ANGLES",    "Lin. Angle",     3, 2, 2, sc_ifcBondedForce | sc_ifcAngle },
     {"CROSS_BOND_BOND",  "Bond-Cross",     3, 3, 0, sc_ifcBondedForce },
     {"CROSS_BOND_ANGLE", "BA-Cross",       3, 4, 0, sc_ifcBondedForce },
     {"UREY_BRADLEY",     "U-B",            3, 4, 4, sc_ifcBondedForce | sc_ifcAngle },
     {"QANGLES",          "Quartic Angles", 3, 6, 0, sc_ifcBondedForce | sc_ifcAngle },
     {"TABANGLES",        "Tab. Angles",    3, 2, 2, sc_ifcBondedForce | sc_ifcTabulated | sc_ifcAngle },
     {"PDIHS",            "Proper Dih.",    4, 3, 3, sc_ifcBondedForce | sc_ifcDihedral },
     {"RBDIHS",           "Ryckaert-Bell.", 4, 6, 6, sc_ifcBondedForce | sc_ifcDihedral },
     {"RESTRDIHS",        "Restr. Dih.",    4, 2, 2, sc_ifcBondedForce | sc_ifcDihedral },
     {"CBTDIHS",          "CBT Dih.",       4, 6, 6, sc_ifcBondedForce | sc_ifcDihedral },
     {"FOURDIHS",         "Fourier Dih.",   4, 4, 4, sc_ifcBondedForce | sc_ifcDihedral },
     {"IDIHS",            "Improper Dih.",  4, 2, 2, sc_ifcBondedForce | sc_ifcDihedral },
     {"PIDIHS",           "Per. Imp. Dih.", 4, 3, 3, sc_ifcBondedForce | sc_ifcDihedral },
     {"TABDIHS",          "Tab. Dih.",      4, 2, 2, sc_ifcBondedForce | sc_ifcDihedral | sc_ifcTabulated },
     {"CMAP",             "CMAP Dih.",      5, 1, 1, sc_ifcBondedForce | sc_ifcDihedral },
     {"GB12",             "GB 1-2 Pol.",    0, 0, 0, sc_ifcNull }, /* unused */
     {"GB13",             "GB 1-3 Pol.",    0, 0, 0, sc_ifcNull }, /* unused */
     {"GB14",             "GB 1-4 Pol.",    0, 0, 0, sc_ifcNull }, /* unused */
     {"GBPOL",            "GB Polariz.",    0, 0, 0, sc_ifcNull }, /* unused */
     {"NPSOLVATION",      "Nonpolar Sol.",  0, 0, 0, sc_ifcNull }, /* unused */
     {"LJ14",             "LJ-14",          2, 2, 2, sc_ifcBondedForce | sc_ifcPair | sc_ifcZeroAtLargeDistances },
     {"COUL14",           "Coulomb-14",     0, 0, 0, sc_ifcNull },
     {"LJC14_Q",          "LJC-14 q",       2, 5, 0, sc_ifcBondedForce | sc_ifcPair | sc_ifcZeroAtLargeDistances },
     {"LJC_NB",           "LJC Pairs NB",   2, 4, 0, sc_ifcBondedForce | sc_ifcPair | sc_ifcZeroAtLargeDistances },
     {"LJ_SR",            "LJ (SR)",        2, 2, 0, sc_ifcNull },
     {"BHAM",             "Buck.ham (SR)",  2, 3, 0, sc_ifcNull },
     {"LJ_LR",            "LJ",             0, 0, 0, sc_ifcNull }, /* unused */
     {"BHAM_LR",          "B.ham",          0, 0, 0, sc_ifcNull }, /* unused */
     {"DISPCORR",         "Disper. corr.",  0, 0, 0, sc_ifcNull },
     {"COUL_SR",          "Coulomb (SR)",   0, 0, 0, sc_ifcNull },
     {"COUL_LR",          "Coul",           0, 0, 0, sc_ifcNull }, /* unused */
     {"RF_EXCL",          "RF excl.",       0, 0, 0, sc_ifcNull },
     {"COUL_RECIP",       "Coul. recip.",   0, 0, 0, sc_ifcNull },
     {"LJ_RECIP",         "LJ recip.",      0, 0, 0, sc_ifcNull },
     {"DPD",              "DPD",            0, 0, 0, sc_ifcNull },
     {"POLARIZATION",     "Polarization",   2, 1, 0, sc_ifcBondedForce | sc_ifcChemicalBond },
     {"WATERPOL",         "Water Pol.",     5, 6, 0, sc_ifcBondedForce },
     {"THOLE",            "Thole Pol.",     4, 3, 0, sc_ifcBondedForce },
     {"ANHARM_POL",       "Anharm. Pol.",   2, 3, 0, sc_ifcBondedForce | sc_ifcChemicalBond },
     {"POSRES",           "Position Rest.", 1, 3, 3, sc_ifcBondedForce },
     {"FBPOSRES",         "Flat-b. P-R.",   1, 3, 0, sc_ifcBondedForce },
     {"DISRES",           "Dis. Rest.",     2, 6, 0, sc_ifcBondedForce },
     {"DISRESVIOL",       "D.R.Viol. (nm)", 0, 0, 0, sc_ifcNull },
     {"ORIRES",           "Orient. Rest.",  2, 6, 0, sc_ifcBondedForce },
     {"ORDEV",            "Ori. R. RMSD",   0, 0, 0, sc_ifcNull },
     {"ANGRES",           "Angle Rest.",    4, 3, 3, sc_ifcBondedForce },
     {"ANGRESZ",          "Angle Rest. Z",  2, 3, 3, sc_ifcBondedForce },
     {"DIHRES",           "Dih. Rest.",     4, 3, 3, sc_ifcBondedForce },
     {"DIHRESVIOL",       "Dih. Rest. Vi.", 0, 0, 0, sc_ifcNull }, /* obsolete */
     {"CONSTR",           "Constraint",     2, 1, 1, sc_ifcConstraint | sc_ifcChemicalBond },
     {"CONSTRNC",         "Constr. No Co.", 2, 1, 1, sc_ifcConstraint },
     {"SETTLE",           "Settle",         3, 2, 0, sc_ifcConstraint | sc_ifcChemicalBond },
     {"VSITE1",           "Virtual site 1", 2, 0, 0, sc_ifcVirtualSite },
     {"VSITE2",           "Virtual site 2", 3, 1, 0, sc_ifcVirtualSite },
     {"VSITE2FD",         "Virt. site 2fd", 3, 1, 0, sc_ifcVirtualSite },
     {"VSITE3",           "Virtual site 3", 4, 2, 0, sc_ifcVirtualSite },
     {"VSITE3FD",         "Virt. site 3fd", 4, 2, 0, sc_ifcVirtualSite },
     {"VSITE3FAD",        "Vir. site 3fad", 4, 2, 0, sc_ifcVirtualSite },
     {"VSITE3OUT",        "Vir. site 3out", 4, 3, 0, sc_ifcVirtualSite },
     {"VSITE4FD",         "Virt. site 4fd", 5, 3, 0, sc_ifcVirtualSite },
     {"VSITE4FDN",        "Vir. site 4fdn", 5, 3, 0, sc_ifcVirtualSite },
     {"VSITEN",           "Virtual site N", 2, 2, 0, sc_ifcVirtualSite },
     {"COM_PULL",         "COM Pull En.",   0, 0, 0, sc_ifcNull },
     {"DENSITYFIT",       "Dens. fitting",  0, 0, 0, sc_ifcNull },
     {"EQM",              "Quantum En.",    0, 0, 0, sc_ifcNull },
     {"ENNPOT",           "NN Potential",   0, 0, 0, sc_ifcNull },
     {"EPOT",             "Potential",      0, 0, 0, sc_ifcNull },
     {"EKIN",             "Kinetic En.",    0, 0, 0, sc_ifcNull },
     {"ETOT",             "Total Energy",   0, 0, 0, sc_ifcNull },
     {"ECONS",            "Conserved En.",  0, 0, 0, sc_ifcNull },
     {"TEMP",             "Temperature",    0, 0, 0, sc_ifcNull },
     {"VTEMP",            "Vir. Temp.",     0, 0, 0, sc_ifcNull }, /* unused */
    /* Note that pressure names can not be more than 8 characters,
     * because " (bar)" is appended to them. */
     {"PDISPCORR",        "Pres. DC",       0, 0, 0, sc_ifcNull },
     {"PRES",             "Pressure",       0, 0, 0, sc_ifcNull },
     {"DH/DL_CON",        "dH/dl constr.",  0, 0, 0, sc_ifcNull }, /* obsolete */
     {"DV/DL",            "dVremain/dl",    0, 0, 0, sc_ifcNull },
     {"DK/DL",            "dEkin/dl",       0, 0, 0, sc_ifcNull },
     {"DVC/DL",           "dVcoul/dl",      0, 0, 0, sc_ifcNull },
     {"DVV/DL",           "dVvdw/dl",       0, 0, 0, sc_ifcNull },
     {"DVB/DL",           "dVbonded/dl",    0, 0, 0, sc_ifcNull },
     {"DVR/DL",           "dVrestraint/dl", 0, 0, 0, sc_ifcNull },
     {"DVT/DL",           "dVtemp/dl",      0, 0, 0, sc_ifcNull }}
};
// clang-format on
