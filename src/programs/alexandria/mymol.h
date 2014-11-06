/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MYMOL_H
#define MYMOL_H

#include "gromacs/utility/real.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/topology/atomprop.h"
#include "gmx_resp.h"
#include "gentop_qgen.h"
#include "gentop_vsite.h"
#include "gentop_core.h"
#include "molprop.h"
#include "molselect.h"
#include "poldata.h"
#include "gauss_io.h"

enum immStatus {
    immUnknown,
    immOK, immZeroDip, immNoQuad, immCharged,
    immAtomTypes, immAtomNumber, immMolpropConv, immBondOrder, immRespInit,
    immChargeGeneration, immLOT,
    immQMInconsistency, immTest, immNoData,
    immGenShells, immGenBonds, immCommProblem, immNR
};

enum eDih {
    edihNo, edihOne, edihAll, edihNR
};

enum eSupport {
    eSupportNo, eSupportLocal, eSupportRemote, eSupportNR
};

namespace alexandria
{

/*! \brief
 * Contains molecular properties from a range of sources.
 * Overloads the regular molprop and adds a lot of functionality.
 * For one thing, it can generate molprop contents from a coordinate
 * file if needed.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MyMol : public MolProp
{
    private:
        //! Gromacs structures
        int             nexcl_;
        int            *symmetric_charges_;
        int            *cgnr_;
        t_excls        *excls;
        GentopVsites    gvt;
        immStatus       immAtoms_, immCharges_, immTopology_;
        std::string     forcefield_;
        bool            bHaveShells_, bHaveVSites_;

        //! Determine whether a molecule has symmetry (within a certain tolerance)
        bool IsSymmetric(real toler);

        //! Generate Atoms based on quantum calculation with specified level of theory
        immStatus GenerateAtoms(gmx_atomprop_t        ap,
                                const char           *lot,
                                ChargeGenerationModel iModel);

        //! Read atoms?

        //! Generate bonds between atoms
        /*int MakeBonds(gmx_poldata_t pd,
                      gmx_conect gc,
                      gmx_bool bPBC,matrix box,gmx_atomprop_t aps,real tol,
                      gmx_bool bMovePlists);
         */
        //! Generate angles, dihedrals, exclusions etc.
        void MakeAngles();

        //! Generate virtual sites
        void MakeVsites(unsigned int flags);
        
    public:
        rvec           *x_, *f_, *buf, mu_exp, mu_calc, mu_esp, coq;
        matrix          box;
        real            dip_exp, mu_exp2, dip_err, dip_weight, dip_calc, chieq, Hform, Emol, Ecalc, Force2;
        real           *qESP;
        tensor          Q_exp, Q_calc, Q_esp;
        int             bts[ebtsNR];
        eSupport        eSupp;
        t_state         state_;
        t_forcerec     *fr_;
        t_params       *plist_;
        gmx_mtop_t     *mtop_;
        gmx_localtop_t *ltop_;
        gpp_atomtype_t  atype_;
        gentop_qgen_t   qgen_;
        t_symtab       *symtab_;
        t_inputrec     *inputrec_;
        gmx_shellfc_t   shell_;
        gmx_enerdata_t  enerd_;
        gmx_resp_t      gr_;
        t_mdatoms      *md_;
        t_topology     *topology_;

        //! Constructor
        MyMol();

        //! Destructor
        ~MyMol();

        //! Generate the topology structure
        immStatus GenerateTopology(gmx_atomprop_t        ap,
                                   gmx_poldata_t         pd,
                                   const char           *lot,
                                   ChargeGenerationModel iModel,
                                   bool                  bPol,
                                   int                   nexcl);
        //! Generate Charges
        immStatus GenerateCharges(gmx_poldata_t pd, gmx_atomprop_t ap,
                                  ChargeGenerationModel iModel, real hfac, real epsr,
                                  const char *lot,
                                  bool bSymmetricCharges,
                                  const char *symm_string);

        // Collect the experimental properties
        immStatus GetExpProps(gmx_bool bQM, gmx_bool bZero, char *lot,
                              alexandria::GaussAtomProp &gap);

        //! Print the topology that was generated previously in GROMACS format.
        //! fp is a File pointer opened previously.
        void PrintTopology(const char           *fn,
                           ChargeGenerationModel iModel,
                           bool                  bVerbose);

        //! Print some info about the molecule to a file
        void PrintQPol(FILE *fp, gmx_poldata_t pd);

        //! Print a rtp entry
        void PrintRTPEntry(const char *fn);

        //! Set the force field
        void SetForceField(const char *ff) { forcefield_.assign(ff); }

        //! Updated internal structures due to changes in pd
        void UpdateIdef(gmx_poldata_t pd, bool bOpt[]);

        //! Get the force field
        std::string GetForceField() { return forcefield_; }

        void CalcMultipoles();

        void GenerateVsitesShells(gmx_poldata_t pd, bool bGenVsites, bool bAddShells,
                                  bool bPairs, eDih edih);

        immStatus GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge,
                                       const char *ndxfn, int nmol);

        immStatus GenerateGromacs(output_env_t oenv, t_commrec *cr);

        void GenerateCube(ChargeGenerationModel iModel,
                          gmx_poldata_t         pd,
                          real                  spacing,
                          const char           *reffn,
                          const char           *pcfn,
                          const char           *pdbdifffn,
                          const char           *potfn,
                          const char           *rhofn,
                          const char           *hisfn,
                          const char           *difffn,
                          const char           *diffhistfn,
                          output_env_t          oenv);

        //! Print the coordinates corresponding to topology after adding shell particles
        //! and/or vsites. fp is a File pointer opened previously.
        void PrintConformation(const char *fn);

        //! Routine initiating the internal GROMACS structures
        /*immStatus Initxx(FILE *fp,GaussAtomProp &gap,
                       gmx_bool bQM,char *lot,gmx_bool bZero,
                       gmx_poldata_t pd,gmx_atomprop_t aps,
                       ChargeGenerationModel iModel,t_commrec *cr,int *nwarn,
                       gmx_bool bCharged,const output_env_t oenv,
                       real th_toler,real ph_toler,
                       real dip_toler,real hfac,gmx_bool bH14,
                       gmx_bool bAllDihedrals,gmx_bool bRemoveDoubleDihedrals,
                       int nexcl,gmx_bool bESP,
                       real watoms,real rDecrZeta,gmx_bool bPol,gmx_bool bFitZeta);*/

};

const char *immsg(immStatus imm);

}

void mv_plists(gmx_poldata_t pd, t_params plist[], gmx_bool bForward);

#define gmx_assert(n, m) if (n != m) { gmx_fatal(FARGS, "Variable %s = %d, should have been %d",#n, n, m); }

#endif
