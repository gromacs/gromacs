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

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/real.h"

#include "gauss_io.h"
#include "gentop_core.h"
#include "gentop_qgen.h"
#include "gentop_vsite.h"
#include "gmx_resp.h"
#include "molprop.h"
#include "molselect.h"
#include "poldata.h"

struct gmx_shellfc_t;
struct t_state;
struct t_topology;

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

enum ePolar {
    epolNo, epolAllAtom, epolUnited, epolNR
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
class MyMol
{
    private:
        //! The molprop
        MolProp         *mp_;
        //! Gromacs structures
        int              nexcl_;
        //! List of symmetric charges
        std::vector<int> symmetric_charges_;
        int             *cgnr_;
        t_excls         *excls_;
        GentopVsites     gvt_;
        immStatus        immAtoms_, immCharges_, immTopology_;
        std::string      forcefield_;
        bool             bHaveShells_, bHaveVSites_;
        double           ref_enthalpy_, mutot_;
        double           polarizability_, sig_pol_;
        //! Determine whether a molecule has symmetry (within a certain tolerance)
        bool IsSymmetric(real toler);

        //! Generate Atoms based on quantum calculation with specified level of theory
        immStatus GenerateAtoms(gmx_atomprop_t            ap,
                                const char               *lot,
                                ChargeDistributionModel   iModel);

        //! Generate angles, dihedrals, exclusions etc.
        void MakeAngles(bool bPairs, bool bDihs);

        //! Generate virtual sites or linear angles
        void MakeSpecialInteractions(bool bUseVsites);

        //! Fetch the force constants
        void getForceConstants(Poldata * pd);
    public:
        rvec                     *x_, *f_, *buf, mu_exp, mu_calc, mu_esp, coq;
        matrix                    box;
        real                      dip_exp, mu_exp2, dip_err, dip_weight, dip_calc, chieq, Hform, Emol, Ecalc, Force2;
        real                     *qESP;
        tensor                    Q_exp, Q_calc, Q_esp;
        eSupport                  eSupp;
        t_state                   state_;
        t_forcerec               *fr_;

        std::vector<PlistWrapper> plist_;

        gmx_mtop_t               *mtop_;
        gmx_localtop_t           *ltop_;
        gpp_atomtype_t            atype_;
        gmx_shellfc_t            *shellfc_;
        GentopQgen               *qgen_;
        t_symtab                 *symtab_;
        t_inputrec               *inputrec_;
        gmx_enerdata_t            enerd_;
        Resp                     *gr_;
        t_mdatoms                *md_;
        t_topology               *topology_;

        //! Constructor
        MyMol();

        //! Destructor
        ~MyMol();

        //! Return my inner molprop
        MolProp *molProp() { return mp_; }
        //! Generate the topology structure
        immStatus GenerateTopology(gmx_atomprop_t          ap,
                                   Poldata                *pd,
                                   const char             *lot,
                                   ChargeDistributionModel iModel,
                                   int                     nexcl,
                                   bool                    bUseVsites,
                                   bool                    bPairs,
                                   bool                    bDih);
        //! Generate Charges
        immStatus GenerateCharges(Poldata *pd, gmx_atomprop_t ap,
                                  ChargeDistributionModel iModel,
                                  ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                                  real hfac, real epsr,
                                  const char *lot,
                                  bool bSymmetricCharges,
                                  const char *symm_string);

        // Collect the experimental properties
        immStatus getExpProps(gmx_bool bQM, gmx_bool bZero, char *lot,
                              alexandria::GaussAtomProp &gap);

        //! Print the topology that was generated previously in GROMACS format.
        //! fp is a File pointer opened previously.
        void PrintTopology(const char             *fn,
                           ChargeDistributionModel iModel,
                           bool                    bVerbose,
                           Poldata           *pd,
                           gmx_atomprop_t          aps);



        //! Compute/derive global info about the molecule
        void CalcQPol(Poldata * pd);


        //! Set the force field
        void SetForceField(const char *ff) { forcefield_.assign(ff); }

        //! Updated internal structures due to changes in pd
        void UpdateIdef(Poldata *pd, bool bOpt[]);

        //! Get the force field
        std::string getForceField() { return forcefield_; }

        void CalcMultipoles();

        void AddShells(Poldata *pd, ePolar epol, ChargeDistributionModel iModel);

        immStatus GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge);

        immStatus GenerateGromacs(t_commrec *cr);

        void GenerateCube(ChargeDistributionModel iModel,
                          Poldata                *pd,
                          real                    spacing,
                          const char             *reffn,
                          const char             *pcfn,
                          const char             *pdbdifffn,
                          const char             *potfn,
                          const char             *rhofn,
                          const char             *hisfn,
                          const char             *difffn,
                          const char             *diffhistfn,
                          const gmx_output_env_t *oenv);

        //! Print the coordinates corresponding to topology after adding shell particles
        //! and/or vsites. fp is a File pointer opened previously.
        void PrintConformation(const char *fn);
};

const char *immsg(immStatus imm);

}

#define gmx_assert(n, m) if (n != m) { gmx_fatal(FARGS, "Variable %s = %d, should have been %d",#n, n, m); }

#endif
