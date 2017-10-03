/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares enumerated types used throughout the code.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_MD_ENUMS_H
#define GMX_MDTYPES_MD_ENUMS_H

#include "gromacs/utility/basedefinitions.h"

/*! \brief Return a string from a list of strings
 *
 * If index if within 0 .. max_index-1 returns the corresponding string
 * or "no name defined" otherwise, in other words this is a range-check that does
 * not crash.
 * \param[in] index     The index in the array
 * \param[in] max_index The length of the array
 * \param[in] names     The array
 * \return the correct string or "no name defined"
 */
const char *enum_name(int index, int max_index, const char *names[]);

//! Boolean strings no or yes
extern const char *yesno_names[BOOL_NR+1];

//! \brief The two compartments for CompEL setups.
enum eCompartment {
    eCompA, eCompB, eCompNR
};

/*! \brief The channels that define with their COM the compartment boundaries in CompEL setups.
 *
 * In principle one could also use modified setups with more than two channels.
 */
enum eChannel {
    eChan0, eChan1, eChanNR
};

/*! \brief Temperature coupling type
 *
 * yes is an alias for berendsen
 */
enum {
    etcNO, etcBERENDSEN, etcNOSEHOOVER, etcYES, etcANDERSEN, etcANDERSENMASSIVE, etcVRESCALE, etcNR
};
//! Strings corresponding to temperatyre coupling types
extern const char *etcoupl_names[etcNR+1];
//! Macro for selecting t coupling string
#define ETCOUPLTYPE(e) enum_name(e, etcNR, etcoupl_names)
//! Return whether this is andersen coupling
#define ETC_ANDERSEN(e) (((e) == etcANDERSENMASSIVE) || ((e) == etcANDERSEN))

/*! \brief Pressure coupling types
 *
 * isotropic is an alias for berendsen
 */
enum {
    epcNO, epcBERENDSEN, epcPARRINELLORAHMAN, epcISOTROPIC, epcMTTK, epcNR
};
//! String corresponding to pressure coupling algorithm
extern const char *epcoupl_names[epcNR+1];
//! Macro to return the correct pcoupling string
#define EPCOUPLTYPE(e) enum_name(e, epcNR, epcoupl_names)

//! Flat-bottom posres geometries
enum {
    efbposresZERO, efbposresSPHERE, efbposresCYLINDER, efbposresX, efbposresY, efbposresZ,
    efbposresCYLINDERX, efbposresCYLINDERY, efbposresCYLINDERZ, efbposresNR
};

//! Relative coordinate scaling type for position restraints.
enum {
    erscNO, erscALL, erscCOM, erscNR
};
//! String corresponding to relativ coordinate scaling.
extern const char *erefscaling_names[erscNR+1];
//! Macro to select correct coordinate scaling string.
#define EREFSCALINGTYPE(e) enum_name(e, erscNR, erefscaling_names)

//! Trotter decomposition extended variable parts.
enum {
    etrtNONE, etrtNHC, etrtBAROV, etrtBARONHC, etrtNHC2, etrtBAROV2, etrtBARONHC2,
    etrtVELOCITY1, etrtVELOCITY2, etrtPOSITION, etrtSKIPALL, etrtNR
};

//! Sequenced parts of the trotter decomposition.
enum {
    ettTSEQ0,  ettTSEQ1,  ettTSEQ2,  ettTSEQ3,  ettTSEQ4, ettTSEQMAX
};

//! Pressure coupling type
enum {
    epctISOTROPIC, epctSEMIISOTROPIC, epctANISOTROPIC,
    epctSURFACETENSION, epctNR
};
//! String corresponding to pressure coupling type
extern const char *epcoupltype_names[epctNR+1];
//! Macro to select the right string for pcoupl type
#define EPCOUPLTYPETYPE(e) enum_name(e, epctNR, epcoupltype_names)

//! \\brief Cutoff scheme
enum {
    ecutsVERLET, ecutsGROUP, ecutsNR
};
//! String corresponding to cutoff scheme
extern const char *ecutscheme_names[ecutsNR+1];
//! Macro to select the right string for cutoff scheme
#define ECUTSCHEME(e)  enum_name(e, ecutsNR, ecutscheme_names)

/*! \brief Coulomb / VdW interaction modifiers.
 *
 * grompp replaces eintmodPOTSHIFT_VERLET by eintmodPOTSHIFT or eintmodNONE.
 * Exactcutoff is only used by Reaction-field-zero, and is not user-selectable.
 */
enum eintmod {
    eintmodPOTSHIFT_VERLET, eintmodPOTSHIFT, eintmodNONE, eintmodPOTSWITCH, eintmodEXACTCUTOFF, eintmodFORCESWITCH, eintmodNR
};
//! String corresponding to interaction modifiers
extern const char *eintmod_names[eintmodNR+1];
//! Macro to select the correct string for modifiers
#define INTMODIFIER(e) enum_name(e, eintmodNR, eintmod_names)

/*! \brief Cut-off treatment for Coulomb
 *
 * eelNOTUSED1 used to be GB, but to enable generalized born with different
 * forms of electrostatics (RF, switch, etc.) in the future it is now selected
 * separately (through the implicit_solvent option).
 */
enum {
    eelCUT,     eelRF,     eelGRF,   eelPME,  eelEWALD,  eelP3M_AD,
    eelPOISSON, eelSWITCH, eelSHIFT, eelUSER, eelGB_NOTUSED, eelRF_NEC_UNSUPPORTED, eelENCADSHIFT,
    eelPMEUSER, eelPMESWITCH, eelPMEUSERSWITCH, eelRF_ZERO, eelNR
};
//! String corresponding to Coulomb treatment
extern const char *eel_names[eelNR+1];
//! Macro for correct string for Coulomb treatment
#define EELTYPE(e)     enum_name(e, eelNR, eel_names)

//! Ewald geometry.
enum {
    eewg3D, eewg3DC, eewgNR
};
//! String corresponding to Ewald geometry
extern const char *eewg_names[eewgNR+1];

//! Macro telling us whether we use reaction field
#define EEL_RF(e) ((e) == eelRF || (e) == eelGRF || (e) == eelRF_NEC_UNSUPPORTED || (e) == eelRF_ZERO )

//! Macro telling us whether we use PME
#define EEL_PME(e)  ((e) == eelPME || (e) == eelPMESWITCH || (e) == eelPMEUSER || (e) == eelPMEUSERSWITCH || (e) == eelP3M_AD)
//! Macro telling us whether we use PME or full Ewald
#define EEL_PME_EWALD(e) (EEL_PME(e) || (e) == eelEWALD)
//! Macro telling us whether we use full electrostatics of any sort
#define EEL_FULL(e) (EEL_PME_EWALD(e) || (e) == eelPOISSON)
//! Macro telling us whether we use user defined electrostatics
#define EEL_USER(e) ((e) == eelUSER || (e) == eelPMEUSER || (e) == (eelPMEUSERSWITCH))

//! Van der Waals interaction treatment
enum {
    evdwCUT, evdwSWITCH, evdwSHIFT, evdwUSER, evdwENCADSHIFT,
    evdwPME, evdwNR
};
//! String corresponding to Van der Waals treatment
extern const char *evdw_names[evdwNR+1];
//! Macro for selecting correct string for VdW treatment
#define EVDWTYPE(e)    enum_name(e, evdwNR, evdw_names)

//! Type of long-range VdW treatment of combination rules
enum {
    eljpmeGEOM, eljpmeLB, eljpmeNR
};
//! String for LJPME combination rule treatment
extern const char *eljpme_names[eljpmeNR+1];
//! Macro for correct LJPME comb rule name
#define ELJPMECOMBNAMES(e) enum_name(e, eljpmeNR, eljpme_names)

//! Macro to tell us whether we use LJPME
#define EVDW_PME(e) ((e) == evdwPME)

//! Neighborsearching algorithm
enum {
    ensGRID, ensSIMPLE, ensNR
};
//! String corresponding to neighborsearching
extern const char *ens_names[ensNR+1];
//! Macro for correct NS algorithm
#define ENS(e)         enum_name(e, ensNR, ens_names)

/*! \brief Integrator algorithm
 *
 * eiSD2 has been removed, but we keep a renamed enum entry,
 * so we can refuse to do MD with such .tpr files.
 * eiVV is normal velocity verlet
 * eiVVAK uses 1/2*(KE(t-dt/2)+KE(t+dt/2)) as the kinetic energy,
 * and the half step kinetic energy for temperature control
 */
enum {
    eiMD, eiSteep, eiCG, eiBD, eiSD2_REMOVED, eiNM, eiLBFGS, eiTPI, eiTPIC, eiSD1, eiVV, eiVVAK, eiNR
};
//! Name of the integrator algorithm
extern const char *ei_names[eiNR+1];
//! Macro returning integrator string
#define EI(e)          enum_name(e, eiNR, ei_names)
//! Do we use velocity Verlet
#define EI_VV(e) ((e) == eiVV || (e) == eiVVAK)
//! Do we use molecular dynamics
#define EI_MD(e) ((e) == eiMD || EI_VV(e))
//! Do we use stochastic dynamics
#define EI_SD(e) ((e) == eiSD1)
//! Do we use any stochastic integrator
#define EI_RANDOM(e) (EI_SD(e) || (e) == eiBD)
/*above integrators may not conserve momenta*/
//! Do we use any type of dynamics
#define EI_DYNAMICS(e) (EI_MD(e) || EI_RANDOM(e))
//! Or do we use minimization
#define EI_ENERGY_MINIMIZATION(e) ((e) == eiSteep || (e) == eiCG || (e) == eiLBFGS)
//! Do we apply test particle insertion
#define EI_TPI(e) ((e) == eiTPI || (e) == eiTPIC)
//! Do we deal with particle velocities
#define EI_STATE_VELOCITY(e) (EI_MD(e) || EI_SD(e))

//! Constraint algorithm
enum {
    econtLINCS, econtSHAKE, econtNR
};
//! String corresponding to constraint algorithm
extern const char *econstr_names[econtNR+1];
//! Macro to select the correct string
#define ECONSTRTYPE(e) enum_name(e, econtNR, econstr_names)

//! Distance restraint refinement algorithm
enum {
    edrNone, edrSimple, edrEnsemble, edrNR
};
//! String corresponding to distance restraint algorithm
extern const char *edisre_names[edrNR+1];
//! Macro to select the right disre algorithm string
#define EDISRETYPE(e)  enum_name(e, edrNR, edisre_names)

//! Distance restraints weighting type
enum {
    edrwConservative, edrwEqual, edrwNR
};
//! String corresponding to distance restraint weighting
extern const char *edisreweighting_names[edrwNR+1];
//! Macro corresponding to dr weighting
#define EDISREWEIGHTING(e)  enum_name(e, edrwNR, edisreweighting_names)

//! Combination rule algorithm.
enum {
    eCOMB_NONE, eCOMB_GEOMETRIC, eCOMB_ARITHMETIC, eCOMB_GEOM_SIG_EPS, eCOMB_NR
};
//! String for combination rule algorithm
extern const char *ecomb_names[eCOMB_NR+1];
//! Macro to select the comb rule string
#define ECOMBNAME(e)   enum_name(e, eCOMB_NR, ecomb_names)

//! Van der Waals potential.
enum {
    eNBF_NONE, eNBF_LJ, eNBF_BHAM, eNBF_NR
};
//! String corresponding to Van der Waals potential
extern const char *enbf_names[eNBF_NR+1];
//! Macro for correct VdW potential string
#define ENBFNAME(e)    enum_name(e, eNBF_NR, enbf_names)

//! Simulated tempering methods.
enum {
    esimtempGEOMETRIC, esimtempEXPONENTIAL, esimtempLINEAR, esimtempNR
};
//! String corresponding to simulated tempering
extern const char *esimtemp_names[esimtempNR+1];
//! Macro for correct tempering string
#define ESIMTEMP(e)    enum_name(e, esimtempNR, esimtemp_names)

/*! \brief Free energy perturbation type
 *
 * efepNO, there are no evaluations at other states.
 * efepYES, treated equivalently to efepSTATIC.
 * efepSTATIC, then lambdas do not change during the simulation.
 * efepSLOWGROWTH, then the states change monotonically
 * throughout the simulation.
 * efepEXPANDED, then expanded ensemble simulations are occuring.
 */
enum {
    efepNO, efepYES, efepSTATIC, efepSLOWGROWTH, efepEXPANDED, efepNR
};
//! String corresponding to FEP type.
extern const char *efep_names[efepNR+1];
//! Macro corresponding to FEP string.
#define EFEPTYPE(e)    enum_name(e, efepNR, efep_names)

//! Free energy pertubation coupling types.
enum {
    efptFEP, efptMASS, efptCOUL, efptVDW, efptBONDED, efptRESTRAINT, efptTEMPERATURE, efptNR
};
//! String for FEP coupling type
extern const char *efpt_names[efptNR+1];
//! Long names for FEP coupling type
extern const char *efpt_singular_names[efptNR+1];

/*! \brief What to print for free energy calculations
 *
 * Printing the energy to the free energy dhdl file.
 * YES is an alias to TOTAL, and
 * will be converted in readir, so we never have to account for it in code.
 */
enum {
    edHdLPrintEnergyNO, edHdLPrintEnergyTOTAL, edHdLPrintEnergyPOTENTIAL, edHdLPrintEnergyYES, edHdLPrintEnergyNR
};
//! String corresponding to printing of free energy
extern const char *edHdLPrintEnergy_names[edHdLPrintEnergyNR+1];

/*! \brief How the lambda weights are calculated
 *
 * elamstatsMETROPOLIS - using the metropolis criteria
 * elamstatsBARKER     - using the Barker critera for transition weights,
 *                       also called unoptimized Bennett
 * elamstatsMINVAR     - using Barker + minimum variance for weights
 * elamstatsWL         - Wang-Landu (using visitation counts)
 * elamstatsWWL        - Weighted Wang-Landau (using optimized Gibbs
 *                       weighted visitation counts)
 */
enum {
    elamstatsNO, elamstatsMETROPOLIS, elamstatsBARKER, elamstatsMINVAR, elamstatsWL, elamstatsWWL, elamstatsNR
};
//! String corresponding to lambda weights
extern const char *elamstats_names[elamstatsNR+1];
//! Macro telling us whether we use expanded ensemble
#define ELAMSTATS_EXPANDED(e) ((e) > elamstatsNO)
//! Macro telling us whether we use some kind of Wang-Landau
#define EWL(e) ((e) == elamstatsWL || (e) == elamstatsWWL)

/*! \brief How moves in lambda are calculated
 *
 * elmovemcMETROPOLIS - using the Metropolis criteria, and 50% up and down
 * elmovemcBARKER     - using the Barker criteria, and 50% up and down
 * elmovemcGIBBS      - computing the transition using the marginalized
 *                      probabilities of the lambdas
 * elmovemcMETGIBBS   - computing the transition using the metropolized
 *                      version of Gibbs (Monte Carlo Strategies in
 *                      Scientific computing, Liu, p. 134)
 */
enum {
    elmcmoveNO, elmcmoveMETROPOLIS, elmcmoveBARKER, elmcmoveGIBBS, elmcmoveMETGIBBS, elmcmoveNR
};
//! String corresponding to lambda moves
extern const char *elmcmove_names[elmcmoveNR+1];

/*! \brief How we decide whether weights have reached equilibrium
 *
 * elmceqNO       - never stop, weights keep going
 * elmceqYES      - fix the weights from the beginning; no movement
 * elmceqWLDELTA  - stop when the WL-delta falls below a certain level
 * elmceqNUMATLAM - stop when we have a certain number of samples at
 *                  every step
 * elmceqSTEPS    - stop when we've run a certain total number of steps
 * elmceqSAMPLES  - stop when we've run a certain total number of samples
 * elmceqRATIO    - stop when the ratio of samples (lowest to highest)
 *                  is sufficiently large
 */
enum {
    elmceqNO, elmceqYES, elmceqWLDELTA, elmceqNUMATLAM, elmceqSTEPS, elmceqSAMPLES, elmceqRATIO, elmceqNR
};
//! String corresponding to equilibrium algorithm
extern const char *elmceq_names[elmceqNR+1];

/*! \brief separate_dhdl_file selection
 *
 * NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool
 */
enum
{
    esepdhdlfileYES, esepdhdlfileNO, esepdhdlfileNR
};
//! String corresponding to separate DHDL file selection
extern const char *separate_dhdl_file_names[esepdhdlfileNR+1];
//! Monster macro for DHDL file selection
#define SEPDHDLFILETYPE(e) enum_name(e, esepdhdlfileNR, separate_dhdl_file_names)

/*! \brief dhdl_derivatives selection \
 *
 * NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool
 */
enum
{
    edhdlderivativesYES, edhdlderivativesNO, edhdlderivativesNR
};
//! String for DHDL derivatives
extern const char *dhdl_derivatives_names[edhdlderivativesNR+1];
//! YAMM (Yet another monster macro)
#define DHDLDERIVATIVESTYPE(e) enum_name(e, edhdlderivativesNR, dhdl_derivatives_names)

/*! \brief Solvent model
 *
 * Distinguishes classical water types with 3 or 4 particles
 */
enum {
    esolNO, esolSPC, esolTIP4P, esolNR
};
//! String corresponding to solvent type
extern const char *esol_names[esolNR+1];
//! Macro lest we print the wrong solvent model string
#define ESOLTYPE(e)    enum_name(e, esolNR, esol_names)

//! Dispersion correction.
enum {
    edispcNO, edispcEnerPres, edispcEner, edispcAllEnerPres, edispcAllEner, edispcNR
};
//! String corresponding to dispersion corrections
extern const char *edispc_names[edispcNR+1];
//! Macro for dispcorr string
#define EDISPCORR(e)   enum_name(e, edispcNR, edispc_names)

//! Center of mass motion removal algorithm.
enum {
    ecmLINEAR, ecmANGULAR, ecmNO, ecmLINEAR_ACCELERATION_CORRECTION, ecmNR
};
//! String corresponding to COM removal
extern const char *ecm_names[ecmNR+1];
//! Macro for COM removal string
#define ECOM(e)        enum_name(e, ecmNR, ecm_names)

//! Algorithm for simulated annealing.
enum {
    eannNO, eannSINGLE, eannPERIODIC, eannNR
};
//! String for simulated annealing
extern const char *eann_names[eannNR+1];
//! And macro for simulated annealing string
#define EANNEAL(e)      enum_name(e, eannNR, eann_names)

//! Implicit solvent algorithms.
enum {
    eisNO, eisGBSA, eisNR
};
//! String corresponding to implicit solvent.
extern const char *eis_names[eisNR+1];
//! Macro for implicit solvent string.
#define EIMPLICITSOL(e) enum_name(e, eisNR, eis_names)

//! Algorithms for calculating GB radii.
enum {
    egbSTILL, egbHCT, egbOBC, egbNR
};
//! String for GB algorithm name.
extern const char *egb_names[egbNR+1];
//! Macro for GB string.
#define EGBALGORITHM(e) enum_name(e, egbNR, egb_names)

//! Surface area algorithm for implicit solvent.
enum {
    esaAPPROX, esaNO, esaSTILL, esaNR
};
//! String corresponding to surface area algorithm.
extern const char *esa_names[esaNR+1];
//! brief Macro for SA algorithm string.
#define ESAALGORITHM(e) enum_name(e, esaNR, esa_names)

//! Wall types.
enum {
    ewt93, ewt104, ewtTABLE, ewt126, ewtNR
};
//! String corresponding to wall type
extern const char *ewt_names[ewtNR+1];
//! Macro for wall type string
#define EWALLTYPE(e)   enum_name(e, ewtNR, ewt_names)

//! Pulling algorithm.
enum {
    epullUMBRELLA, epullCONSTRAINT, epullCONST_F, epullFLATBOTTOM, epullFLATBOTTOMHIGH, epullEXTERNAL, epullNR
};
//! String for pulling algorithm
extern const char *epull_names[epullNR+1];
//! Macro for pulling string
#define EPULLTYPE(e)   enum_name(e, epullNR, epull_names)

//! Control of pull groups
enum {
    epullgDIST, epullgDIR, epullgCYL, epullgDIRPBC, epullgDIRRELATIVE, epullgANGLE, epullgDIHEDRAL, epullgANGLEAXIS, epullgNR
};
//! String for pull groups
extern const char *epullg_names[epullgNR+1];
//! Macro for pull group string
#define EPULLGEOM(e)   enum_name(e, epullgNR, epullg_names)

//! Enforced rotation groups.
enum {
    erotgISO, erotgISOPF,
    erotgPM, erotgPMPF,
    erotgRM, erotgRMPF,
    erotgRM2, erotgRM2PF,
    erotgFLEX, erotgFLEXT,
    erotgFLEX2, erotgFLEX2T,
    erotgNR
};
//! Rotation group names
extern const char *erotg_names[erotgNR+1];
//! Macro for rot group names
#define EROTGEOM(e)    enum_name(e, erotgNR, erotg_names)
//! String for rotation group origin names
extern const char *erotg_originnames[erotgNR+1];
//! Macro for rot group origin names
#define EROTORIGIN(e)  enum_name(e, erotgOriginNR, erotg_originnames)

//! Rotation group fitting type
enum {
    erotgFitRMSD, erotgFitNORM, erotgFitPOT, erotgFitNR
};
//! String corresponding to rotation group fitting
extern const char *erotg_fitnames[erotgFitNR+1];
//! Macro for rot group fit names
#define EROTFIT(e)     enum_name(e, erotgFitNR, erotg_fitnames)

/*! \brief Direction along which ion/water swaps happen
 *
 * Part of "Computational Electrophysiology" (CompEL) setups
 */
enum eSwaptype {
    eswapNO, eswapX, eswapY, eswapZ, eSwapTypesNR
};
//! Names for swapping
extern const char *eSwapTypes_names[eSwapTypesNR+1];
//! Macro for swapping string
#define ESWAPTYPE(e)   enum_name(e, eSwapTypesNR, eSwapTypes_names)

/*! \brief Swap group splitting type
 *
 * These are just the fixed groups we need for any setup. In t_swap's grp
 * entry after that follows the variable number of swap groups.
 */
enum {
    eGrpSplit0, eGrpSplit1, eGrpSolvent, eSwapFixedGrpNR
};
//! String for swap group splitting
extern const char *eSwapFixedGrp_names[eSwapFixedGrpNR+1];

//! QMMM methods.
enum {
    eQMmethodAM1, eQMmethodPM3, eQMmethodRHF,
    eQMmethodUHF, eQMmethodDFT, eQMmethodB3LYP, eQMmethodMP2, eQMmethodCASSCF, eQMmethodB3LYPLAN,
    eQMmethodDIRECT, eQMmethodNR
};
//! String corresponding to QMMM methods
extern const char *eQMmethod_names[eQMmethodNR+1];
//! Macro to pick QMMM method name
#define EQMMETHOD(e)   enum_name(e, eQMmethodNR, eQMmethod_names)

//! QMMM basis function for QM part
enum {
    eQMbasisSTO3G, eQMbasisSTO3G2, eQMbasis321G,
    eQMbasis321Gp, eQMbasis321dGp, eQMbasis621G,
    eQMbasis631G, eQMbasis631Gp, eQMbasis631dGp,
    eQMbasis6311G, eQMbasisNR
};
//! Name for QMMM basis function
extern const char *eQMbasis_names[eQMbasisNR+1];
//! Macro to pick right basis function string
#define EQMBASIS(e)    enum_name(e, eQMbasisNR, eQMbasis_names)

//! QMMM scheme
enum {
    eQMMMschemenormal, eQMMMschemeoniom, eQMMMschemeNR
};
//! QMMMM scheme names
extern const char *eQMMMscheme_names[eQMMMschemeNR+1];
//! Macro to pick QMMMM scheme name
#define EQMMMSCHEME(e) enum_name(e, eQMMMschemeNR, eQMMMscheme_names)

/*! \brief Neighborlist geometry type.
 *
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
//! String corresponding to nblist geometry names
extern const char *gmx_nblist_geometry_names[GMX_NBLIST_GEOMETRY_NR+1];

/*! \brief Types of electrostatics calculations
 *
 * Types of electrostatics calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections
 * in the MDP file; many interactions for instance map to tabulated kernels.
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
//! String corresponding to electrostatics kernels
extern const char *gmx_nbkernel_elec_names[GMX_NBKERNEL_ELEC_NR+1];

/*! \brief Types of vdw calculations available
 *
 * Types of vdw calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections
 * in the MDP file; many interactions for instance map to tabulated kernels.
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
//! String corresponding to VdW kernels
extern const char *gmx_nbkernel_vdw_names[GMX_NBKERNEL_VDW_NR+1];

//! \brief Types of interactions inside the neighborlist
enum gmx_nblist_interaction_type
{
    GMX_NBLIST_INTERACTION_STANDARD,
    GMX_NBLIST_INTERACTION_FREE_ENERGY,
    GMX_NBLIST_INTERACTION_NR
};
//! String corresponding to interactions in neighborlist code
extern const char *gmx_nblist_interaction_names[GMX_NBLIST_INTERACTION_NR+1];

#endif /* GMX_MDTYPES_MD_ENUMS_H */
