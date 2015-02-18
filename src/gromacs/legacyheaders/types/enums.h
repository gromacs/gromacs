/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#ifndef ENUMS_H_
#define ENUMS_H_

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/* note: these enums should correspond to the names in gmxlib/names.c */

enum {
    epbcXYZ, epbcNONE, epbcXY, epbcSCREW, epbcNR
};

enum {
    etcNO, etcBERENDSEN, etcNOSEHOOVER, etcYES, etcANDERSEN, etcANDERSENMASSIVE, etcVRESCALE, etcNR
}; /* yes is an alias for berendsen */

#define ETC_ANDERSEN(e) (((e) == etcANDERSENMASSIVE) || ((e) == etcANDERSEN))

enum {
    epcNO, epcBERENDSEN, epcPARRINELLORAHMAN, epcISOTROPIC, epcMTTK, epcNR
}; /* isotropic is an alias for berendsen */

/* trotter decomposition extended variable parts */
enum {
    etrtNONE, etrtNHC, etrtBAROV, etrtBARONHC, etrtNHC2, etrtBAROV2, etrtBARONHC2,
    etrtVELOCITY1, etrtVELOCITY2, etrtPOSITION, etrtSKIPALL, etrtNR
};

/* sequenced parts of the trotter decomposition */
enum {
    ettTSEQ0,  ettTSEQ1,  ettTSEQ2,  ettTSEQ3,  ettTSEQ4, ettTSEQMAX
};

enum {
    epctISOTROPIC, epctSEMIISOTROPIC, epctANISOTROPIC,
    epctSURFACETENSION, epctNR
};

enum {
    erscNO, erscALL, erscCOM, erscNR
};

enum {
    ecutsVERLET, ecutsGROUP, ecutsNR
};

/* Coulomb / VdW interaction modifiers.
 * grompp replaces eintmodPOTSHIFT_VERLET by eintmodPOTSHIFT or eintmodNONE.
 * Exactcutoff is only used by Reaction-field-zero, and is not user-selectable.
 */
enum eintmod {
    eintmodPOTSHIFT_VERLET, eintmodPOTSHIFT, eintmodNONE, eintmodPOTSWITCH, eintmodEXACTCUTOFF, eintmodFORCESWITCH, eintmodNR
};

/*
 * eelNOTUSED1 used to be GB, but to enable generalized born with different
 * forms of electrostatics (RF, switch, etc.) in the future it is now selected
 * separately (through the implicit_solvent option).
 */
enum {
    eelCUT,     eelRF,     eelGRF,   eelPME,  eelEWALD,  eelP3M_AD,
    eelPOISSON, eelSWITCH, eelSHIFT, eelUSER, eelGB_NOTUSED, eelRF_NEC, eelENCADSHIFT,
    eelPMEUSER, eelPMESWITCH, eelPMEUSERSWITCH, eelRF_ZERO, eelNR
};

/* Ewald geometry */
enum {
    eewg3D, eewg3DC, eewgNR
};

#define EEL_RF(e) ((e) == eelRF || (e) == eelGRF || (e) == eelRF_NEC || (e) == eelRF_ZERO )

#define EEL_PME(e)  ((e) == eelPME || (e) == eelPMESWITCH || (e) == eelPMEUSER || (e) == eelPMEUSERSWITCH || (e) == eelP3M_AD)
#define EEL_PME_EWALD(e) (EEL_PME(e) || (e) == eelEWALD)
#define EEL_FULL(e) (EEL_PME_EWALD(e) || (e) == eelPOISSON)

#define EEL_USER(e) ((e) == eelUSER || (e) == eelPMEUSER || (e) == (eelPMEUSERSWITCH))

enum {
    evdwCUT, evdwSWITCH, evdwSHIFT, evdwUSER, evdwENCADSHIFT,
    evdwPME, evdwNR
};

enum {
    eljpmeGEOM, eljpmeLB, eljpmeNR
};

#define EVDW_PME(e) ((e) == evdwPME)

enum {
    ensGRID, ensSIMPLE, ensNR
};

/* eiVV is normal velocity verlet -- eiVVAK uses 1/2*(KE(t-dt/2)+KE(t+dt/2)) as the kinetic energy, and the half step kinetic
   energy for temperature control */

enum {
    eiMD, eiSteep, eiCG, eiBD, eiSD2, eiNM, eiLBFGS, eiTPI, eiTPIC, eiSD1, eiVV, eiVVAK, eiNR
};
#define EI_VV(e) ((e) == eiVV || (e) == eiVVAK)
#define EI_MD(e) ((e) == eiMD || EI_VV(e))
#define EI_SD(e) ((e) == eiSD1 || (e) == eiSD2)
#define EI_RANDOM(e) (EI_SD(e) || (e) == eiBD)
/*above integrators may not conserve momenta*/
#define EI_DYNAMICS(e) (EI_MD(e) || EI_SD(e) || (e) == eiBD)
#define EI_ENERGY_MINIMIZATION(e) ((e) == eiSteep || (e) == eiCG || (e) == eiLBFGS)
#define EI_TPI(e) ((e) == eiTPI || (e) == eiTPIC)

#define EI_STATE_VELOCITY(e) (EI_MD(e) || EI_SD(e))

enum {
    econtLINCS, econtSHAKE, econtNR
};

enum {
    edrNone, edrSimple, edrEnsemble, edrNR
};

enum {
    edrwConservative, edrwEqual, edrwNR
};

/* Combination rule things */
enum {
    eCOMB_NONE, eCOMB_GEOMETRIC, eCOMB_ARITHMETIC, eCOMB_GEOM_SIG_EPS, eCOMB_NR
};

/* NBF selection */
enum {
    eNBF_NONE, eNBF_LJ, eNBF_BHAM, eNBF_NR
};

/* simulated tempering methods */
enum {
    esimtempGEOMETRIC, esimtempEXPONENTIAL, esimtempLINEAR, esimtempNR
};
/* FEP selection */
enum {
    efepNO, efepYES, efepSTATIC, efepSLOWGROWTH, efepEXPANDED, efepNR
};
/* if efepNO, there are no evaluations at other states.
   if efepYES, treated equivalently to efepSTATIC.
   if efepSTATIC, then lambdas do not change during the simulation.
   if efepSLOWGROWTH, then the states change monotonically throughout the simulation.
   if efepEXPANDED, then expanded ensemble simulations are occuring.
 */

/* FEP coupling types */
enum {
    efptFEP, efptMASS, efptCOUL, efptVDW, efptBONDED, efptRESTRAINT, efptTEMPERATURE, efptNR
};

/* Printing the energy to the free energy dhdl file. YES is an alias to TOTAL, and
 * will be converted in readir, so we never have to account for it in code.
 */
enum {
    edHdLPrintEnergyNO, edHdLPrintEnergyTOTAL, edHdLPrintEnergyPOTENTIAL, edHdLPrintEnergyYES, edHdLPrintEnergyNR
};

/* How the lambda weights are calculated:
   elamstatsMETROPOLIS = using the metropolis criteria
   elamstatsBARKER = using the Barker critera for transition weights - also called unoptimized Bennett
   elamstatsMINVAR = using Barker + minimum variance for weights
   elamstatsWL = Wang-Landu (using visitation counts)
   elamstatsWWL = Weighted Wang-Landau (using optimized gibbs weighted visitation counts)
 */
enum {
    elamstatsNO, elamstatsMETROPOLIS, elamstatsBARKER, elamstatsMINVAR, elamstatsWL, elamstatsWWL, elamstatsNR
};

#define ELAMSTATS_EXPANDED(e) ((e) > elamstatsNO)

#define EWL(e) ((e) == elamstatsWL || (e) == elamstatsWWL)

/* How moves in lambda are calculated:
   elmovemcMETROPOLIS - using the Metropolis criteria, and 50% up and down
   elmovemcBARKER - using the Barker criteria, and 50% up and down
   elmovemcGIBBS - computing the transition using the marginalized probabilities of the lambdas
   elmovemcMETGIBBS - computing the transition using the metropolized version of Gibbs (Monte Carlo Strategies in Scientific computing, Liu, p. 134)
 */
enum {
    elmcmoveNO, elmcmoveMETROPOLIS, elmcmoveBARKER, elmcmoveGIBBS, elmcmoveMETGIBBS, elmcmoveNR
};

/* how we decide whether weights have reached equilibrium
   elmceqNO - never stop, weights keep going
   elmceqYES - fix the weights from the beginning; no movement
   elmceqWLDELTA - stop when the WL-delta falls below a certain level
   elmceqNUMATLAM - stop when we have a certain number of samples at every step
   elmceqSTEPS - stop when we've run a certain total number of steps
   elmceqSAMPLES - stop when we've run a certain total number of samples
   elmceqRATIO - stop when the ratio of samples (lowest to highest) is sufficiently large
 */
enum {
    elmceqNO, elmceqYES, elmceqWLDELTA, elmceqNUMATLAM, elmceqSTEPS, elmceqSAMPLES, elmceqRATIO, elmceqNR
};

/* separate_dhdl_file selection */
enum
{
    /* NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool */
    esepdhdlfileYES, esepdhdlfileNO, esepdhdlfileNR
};

/* dhdl_derivatives selection */
enum
{
    /* NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool */
    edhdlderivativesYES, edhdlderivativesNO, edhdlderivativesNR
};

/* Solvent model */
enum {
    esolNO, esolSPC, esolTIP4P, esolNR
};

/* Dispersion correction */
enum {
    edispcNO, edispcEnerPres, edispcEner, edispcAllEnerPres, edispcAllEner, edispcNR
};

/* Center of mass motion selection */
enum {
    ecmLINEAR, ecmANGULAR, ecmNO, ecmNR
};

/* New version of simulated annealing */
enum {
    eannNO, eannSINGLE, eannPERIODIC, eannNR
};

/* Implicit solvent algorithms */
enum {
    eisNO, eisGBSA, eisNR
};

/* Algorithms for calculating GB radii */
enum {
    egbSTILL, egbHCT, egbOBC, egbNR
};

enum {
    esaAPPROX, esaNO, esaSTILL, esaNR
};

/* Wall types */
enum {
    ewt93, ewt104, ewtTABLE, ewt126, ewtNR
};

/* Pull stuff */
enum {
    epullUMBRELLA, epullCONSTRAINT, epullCONST_F, epullFLATBOTTOM, epullNR
};

enum {
    epullgDIST, epullgDIR, epullgCYL, epullgDIRPBC, epullgDIRRELATIVE, epullgNR
};

/* Enforced rotation groups */
enum {
    erotgISO, erotgISOPF,
    erotgPM, erotgPMPF,
    erotgRM, erotgRMPF,
    erotgRM2, erotgRM2PF,
    erotgFLEX, erotgFLEXT,
    erotgFLEX2, erotgFLEX2T,
    erotgNR
};

enum {
    erotgFitRMSD, erotgFitNORM, erotgFitPOT, erotgFitNR
};

/* Direction along which ion/water swaps happen in "Computational
 * Electrophysiology" (CompEL) setups */
enum eSwaptype {
    eswapNO, eswapX, eswapY, eswapZ, eSwapTypesNR
};

/* QMMM */
enum {
    eQMmethodAM1, eQMmethodPM3, eQMmethodRHF,
    eQMmethodUHF, eQMmethodDFT, eQMmethodB3LYP, eQMmethodMP2, eQMmethodCASSCF, eQMmethodB3LYPLAN,
    eQMmethodDIRECT, eQMmethodNR
};

enum {
    eQMbasisSTO3G, eQMbasisSTO3G2, eQMbasis321G,
    eQMbasis321Gp, eQMbasis321dGp, eQMbasis621G,
    eQMbasis631G, eQMbasis631Gp, eQMbasis631dGp,
    eQMbasis6311G, eQMbasisNR
};

enum {
    eQMMMschemenormal, eQMMMschemeoniom, eQMMMschemeNR
};

enum {
    eMultentOptName, eMultentOptNo, eMultentOptLast, eMultentOptNR
};

/* flat-bottom posres geometries */
enum {
    efbposresZERO, efbposresSPHERE, efbposresCYLINDER, efbposresX, efbposresY, efbposresZ,
    efbposresCYLINDERX, efbposresCYLINDERY, efbposresCYLINDERZ, efbposresNR
};

enum {
    eAdressOff, eAdressConst, eAdressXSplit, eAdressSphere, eAdressNR
};

enum {
    eAdressICOff, eAdressICThermoForce, eAdressICNR
};

enum {
    eAdressSITEcom, eAdressSITEcog, eAdressSITEatom, eAdressSITEatomatom, eAdressSITENR
};


/* The interactions contained in a (possibly merged) table
 * for computing electrostatic, VDW repulsion and/or VDW dispersion
 * contributions.
 */
enum gmx_table_interaction
{
    GMX_TABLE_INTERACTION_ELEC,
    GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWDISP,
    GMX_TABLE_INTERACTION_NR
};

/* Different formats for table data. Cubic spline tables are typically stored
 * with the four Y,F,G,H intermediate values (check tables.c for format), which
 * makes it easy to load with a single 4-way SIMD instruction too.
 * Linear tables only need one value per table point, or two if both V and F
 * are calculated. However, with SIMD instructions this makes the loads unaligned,
 * and in that case we store the data as F, D=F(i+1)-F(i), V, and then a blank value,
 * which again makes it possible to load as a single instruction.
 */
enum gmx_table_format
{
    GMX_TABLE_FORMAT_CUBICSPLINE_YFGH,
    GMX_TABLE_FORMAT_LINEAR_VF,
    GMX_TABLE_FORMAT_LINEAR_V,
    GMX_TABLE_FORMAT_LINEAR_F,
    GMX_TABLE_FORMAT_LINEAR_FDV0,
    GMX_TABLE_FORMAT_NR
};

/* Neighborlist geometry type.
 * Kernels will compute interactions between two particles,
 * 3-center water, 4-center water or coarse-grained beads.
 */
enum gmx_nblist_kernel_geometry
{
    GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE,
    GMX_NBLIST_GEOMETRY_WATER3_PARTICLE,
    GMX_NBLIST_GEOMETRY_WATER3_WATER3,
    GMX_NBLIST_GEOMETRY_WATER4_PARTICLE,
    GMX_NBLIST_GEOMETRY_WATER4_WATER4,
    GMX_NBLIST_GEOMETRY_CG_CG,
    GMX_NBLIST_GEOMETRY_NR
};

/* Types of electrostatics calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections in the MDP file;
 * many interactions for instance map to tabulated kernels.
 */
enum gmx_nbkernel_elec
{
    GMX_NBKERNEL_ELEC_NONE,
    GMX_NBKERNEL_ELEC_COULOMB,
    GMX_NBKERNEL_ELEC_REACTIONFIELD,
    GMX_NBKERNEL_ELEC_CUBICSPLINETABLE,
    GMX_NBKERNEL_ELEC_GENERALIZEDBORN,
    GMX_NBKERNEL_ELEC_EWALD,
    GMX_NBKERNEL_ELEC_NR
};

/* Types of vdw calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections in the MDP file;
 * many interactions for instance map to tabulated kernels.
 */
enum gmx_nbkernel_vdw
{
    GMX_NBKERNEL_VDW_NONE,
    GMX_NBKERNEL_VDW_LENNARDJONES,
    GMX_NBKERNEL_VDW_BUCKINGHAM,
    GMX_NBKERNEL_VDW_CUBICSPLINETABLE,
    GMX_NBKERNEL_VDW_LJEWALD,
    GMX_NBKERNEL_VDW_NR
};
/* Types of interactions inside the neighborlist
 */
enum gmx_nblist_interaction_type
{
    GMX_NBLIST_INTERACTION_STANDARD,
    GMX_NBLIST_INTERACTION_FREE_ENERGY,
    GMX_NBLIST_INTERACTION_ADRESS,
    GMX_NBLIST_INTERACTION_NR
};

#ifdef __cplusplus
}
#endif

#endif /* ENUMS_H_ */
