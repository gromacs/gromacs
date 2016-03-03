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
#include "gentop_vsite.h"
#include "molprop.h"
#include "molselect.h"
#include "poldata.h"
#include "qgen_eem.h"
#include "qgen_resp.h"

struct gmx_enerdata_t;
struct gmx_shellfc_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
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
        QgenResp         gr_;
        immStatus        immAtoms_, immCharges_, immTopology_;
        std::string      forcefield_;
        bool             bHaveShells_, bHaveVSites_;
        double           ref_enthalpy_, mutot_;
        double           polarizability_, sig_pol_;
        double           EspRms_;
        //! Determine whether a molecule has symmetry (within a certain tolerance)
        bool IsSymmetric(real toler);

        //! Generate Atoms based on quantum calculation with specified level of theory
        immStatus GenerateAtoms(gmx_atomprop_t            ap,
                                const char               *lot,
                                ChargeDistributionModel   iModel);

        //! Generate angles, dihedrals, exclusions etc.
        void MakeAngles(bool bPairs, 
                        bool bDihs);

        //! Generate virtual sites or linear angles
        void MakeSpecialInteractions(bool bUseVsites);

        //! Add shell particles
        void addShells(const Poldata &pd, ChargeDistributionModel iModel);

        //! Check whether atom types exist in the force field
        immStatus checkAtoms(const Poldata &pd);

        //! Fetch the force constants
        void getForceConstants(const Poldata &pd);
    public:
        rvec                     *x_, *f_, *buf, mu_exp, mu_calc, mu_esp, coq;
        matrix                    box_;
        real                      dip_exp, mu_exp2, dip_err, dip_weight, dip_calc, chieq, Hform, Emol, Ecalc, Force2;
        real                     *qESP;
        tensor                    Q_exp, Q_calc, Q_esp;
        eSupport                  eSupp;
        t_state                  *state_;
        t_forcerec               *fr_;

        std::vector<PlistWrapper> plist_;

        gmx_mtop_t               *mtop_;
        gmx_localtop_t           *ltop_;
        gpp_atomtype_t            atype_;
        gmx_shellfc_t            *shellfc_;
        t_symtab                 *symtab_;
        t_inputrec               *inputrec_;
        gmx_enerdata_t           *enerd_;
        t_mdatoms                *mdatoms_;
        t_topology               *topology_;

        //! Constructor
        MyMol();

        //! Destructor
        ~MyMol();

        //! Return my inner molprop
        MolProp *molProp() const { return mp_; }

        //! Generate the topology structure
        immStatus GenerateTopology(gmx_atomprop_t            ap,
                                   const Poldata            &pd,
                                   const char               *lot,
                                   ChargeDistributionModel   iModel,
                                   int                       nexcl,
                                   bool                      bUseVsites,
                                   bool                      bPairs,
                                   bool                      bDih,
                                   bool                      bAddShells);
        //! Generate Charges
        immStatus GenerateCharges(const Poldata             &pd, 
                                  gmx_atomprop_t             ap,
                                  ChargeDistributionModel    iModel,
                                  ChargeGenerationAlgorithm  iChargeGenerationAlgorithm,
                                  real                       watoms,
                                  real                       hfac,
                                  const char                *lot,
                                  bool                       bSymmetricCharges,
                                  const char                *symm_string,
                                  t_commrec                 *cr,
                                  const char                *tabfn);

        double espRms() const { return EspRms_; }
        
        // Collect the experimental properties
        immStatus getExpProps(gmx_bool bQM, gmx_bool bZero, char *lot,
                              alexandria::GaussAtomProp &gap);

        //! Print the topology that was generated previously in GROMACS format.
        //! fp is a File pointer opened previously.
        void PrintTopology(const char             *fn,
                           ChargeDistributionModel iModel,
                           bool                    bVerbose,
                           const Poldata          &pd,
                           gmx_atomprop_t          aps);

        //! Compute/derive global info about the molecule
        void CalcQPol(const Poldata &pd);

        //! Relax the shells (if any) in the molecule
        void relaxShells(t_commrec *cr);

        //! Set the force field
        void SetForceField(const char *ff) { forcefield_.assign(ff); }

        //! Update internal structures for bondtype due to changes in pd
        void UpdateIdef(const Poldata &pd,
                        int            bondtype);

        //! Get the force field
        std::string getForceField() { return forcefield_; }

        void CalcMultipoles();

        immStatus GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge);

        immStatus GenerateGromacs(t_commrec *cr,
                                  const char *tabfn);

        void GenerateCube(ChargeDistributionModel iModel,
                          const Poldata          &pd,
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

#endif
