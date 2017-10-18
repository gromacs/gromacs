/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012-2016, by the GROMACS development team, led by
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
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "mymol_low.h"
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
struct t_topology;

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
        /*! \brief
         * The molprop
         */
        MolProp         *mp_;
        /*! \brief
         * Gromacs structures
         */
        int              nexcl_;
        int             *cgnr_;
        bool             bHaveShells_;
        bool             bHaveVSites_;
        bool             bNeedVsites_;
        double           ref_enthalpy_;
        double           polarizability_;
        double           sig_pol_;
        double           EspRms_;
        double           EemRms_;
        t_excls         *excls_;
        immStatus        immAtoms_;
        immStatus        immCharges_;
        immStatus        immTopology_;
        gmx_vsite_t     *vsite_;
        GentopVsites     gvt_;
        std::string      forcefield_;


        bool             IsSymmetric(real toler);

        /*! \brief
         * Generate Atoms based on quantum calculation with specified level of theory
         *
         * \param[in] ap      Gromacs atom properties
         * \param[in] lot     Level of theory used for QM calculations
         * \param[in] iModel  The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         */
        immStatus GenerateAtoms(gmx_atomprop_t            ap,
                                const char               *lot,
                                ChargeDistributionModel   iModel);

        /*! \brief
         * Generate angles, dihedrals, exclusions etc.
         *
         * \param[in] bPairs
         * \param[in] bDihs
         */
        void MakeAngles(bool bPairs,
                        bool bDihs);

        /*! \brief
         * Generate virtual sites or linear angles
         *
         * \param[in] bUseVsites
         */
        void MakeSpecialInteractions(const Poldata &pd,
                                     bool           bUseVsites);

        /*! \brief
         * Add shell particles
         *
         * \param[in] pd       Data structure containing atomic properties
         * \param[in] iModel   The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         */
        void addShells(const Poldata &pd, ChargeDistributionModel iModel);

        /*! \brief
         * Check whether atom types exist in the force field
         *
         * \param[in] pd
         */
        immStatus checkAtoms(const Poldata &pd);
        
        
        immStatus zeta2atoms(ChargeDistributionModel iChargeDistributionModel,
                             const Poldata &pd);

        bool IsVsiteNeeded(std::string       atype,
                           const Poldata     &pd);
                           
        void findInPlaneAtoms(int ca, std::vector<int> &atoms);
        
        void findOutPlaneAtoms(int ca, std::vector<int> &atoms);
       
    public:
        rvec                     *buf_;
        rvec                      coc_           = {0,0,0};
        rvec                      mu_elec_       = {0,0,0};
        rvec                      mu_calc_       = {0,0,0};
        rvec                      mu_esp_        = {0,0,0};
        rvec                      mu_mulliken_   = {0,0,0};
        rvec                      mu_hirshfeld_  = {0,0,0};
        rvec                      mu_cm5_        = {0,0,0};
        double                    dip_exp_       = 0;
        double                    dip_elec_      = 0;
        double                    dip_calc_      = 0;
        double                    dip_esp_       = 0;
        double                    dip_mulliken_  = 0;
        double                    dip_hirshfeld_ = 0;
        double                    dip_cm5_       = 0;
        double                   *qESP_;
        double                   *qMulliken_;
        double                   *qCM5_;
        double                   *qHirshfeld_;       
        double                    mu_elec2_      = 0;
        double                    dip_err_       = 0;
        double                    dip_weight_    = 0;
        double                    chieq_         = 0;
        double                    Hform_         = 0;
        double                    Emol_          = 0;
        double                    Ecalc_         = 0;
        double                    OptEcalc_      = 0;
        double                    Force2_        = 0;
        double                    OptForce2_     = 0;
        double                    isoPol_elec_   = 0;
        double                    isoPol_calc_   = 0;
        double                    anisoPol_elec_ = 0;
        double                    anisoPol_calc_ = 0;
        matrix                    box_;
        tensor                    Q_elec_        = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    Q_calc_        = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    Q_esp_         = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    Q_mulliken_    = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    Q_hirshfeld_   = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    Q_cm5_         = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    alpha_elec_    = {{0,0,0},{0,0,0},{0,0,0}};
        tensor                    alpha_calc_    = {{0,0,0},{0,0,0},{0,0,0}};
        eSupport                  eSupp_;
        t_state                  *state_;
        t_forcerec               *fr_;
        PaddedRVecVector          f_;
        PaddedRVecVector          optf_;
        std::vector<int>          symmetric_charges_;
        QgenResp                  Qgresp_;
        QgenEem                   Qgeem_;
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
        t_fcdata                 *fcd_;
        t_nrnb                   nrnb_;
        gmx_wallcycle_t          wcycle_;

        /*! \brief
         * Constructor
         */
        MyMol();

        /*! \brief
         * Return my inner molprop
         */
        MolProp *molProp() const { return mp_; }

        /*! \brief
         * It generates the topology structure which will be used to print the topology file.
         *
         * \param[in] ap          Gromacs atom properties
         * \param[in] pd          Data structure containing atomic properties
         * \param[in] lot         The level of theory used for QM calculation
         * \param[in] iModel      The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] nexcl       Number of Exclusions
         * \param[in] bUseVsites  Add virtual sites to the topology structure
         * \param[in] bPairs      Add pairs to the topology structure
         * \param[in] bDih        Add dihedrals to the topology structure
         * \param[in] bAddShells  Add shells to the topology structure
         */
        immStatus GenerateTopology(gmx_atomprop_t            ap,
                                   const Poldata            &pd,
                                   const char               *lot,
                                   ChargeDistributionModel   iModel,
                                   bool                      bUseVsites,
                                   bool                      bPairs,
                                   bool                      bDih,
                                   bool                      bAddShells,
                                   const char               *tabfn);
        /*! \brief
         *  Computes isotropic polarizability at the presence of external
         *  electric field (under construction!!!)
         *
         * \param[in]  efield   Strenght of the external electric field
         * \param[in]  fplog
         * \param[in]  cr
         * \param[out] isoPol   Isotropic polarizability
         */
        void CalcPolarizability(double efield, t_commrec *cr, FILE *fplog);

        /*! \brief
         * Generate atomic partial charges
         *
         * \param[in] pd                             Data structure containing atomic properties
         * \param[in] ap                             Gromacs atom properties
         * \param[in] iModel                         The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] iChargeGenerationAlgorithm     The algorithm calculating the partial charge (e.x. ESP, RESP)
         * \param[in] watoms
         * \param[in] hfac
         * \param[in] lot                            The level of theory used for QM calculation
         * \param[in] bSymmetricCharges              Consider molecular symmetry to calculate partial charge
         * \param[in] symm_string                    The type of molecular symmetry
         * \param[in] cr
         * \param[in] tabfn
         */
        immStatus GenerateCharges(const Poldata             &pd,
                                  const gmx::MDLogger       &fplog,
                                  gmx_atomprop_t             ap,
                                  ChargeDistributionModel    iModel,
                                  ChargeGenerationAlgorithm  iChargeGenerationAlgorithm,
                                  real                       watoms,
                                  real                       hfac,
                                  const char                *lot,
                                  bool                       bSymmetricCharges,
                                  const char                *symm_string,
                                  t_commrec                 *cr,
                                  const char                *tabfn,
                                  gmx_hw_info_t             *hwinfo,
                                  int                        qcycle,
                                  real                       qtol,
                                  const gmx_output_env_t    *oenv);

        /*! \brief
         * Return the root-mean square deviation of
         * the generated ESP from the QM ESP.
         *
         */
        double espRms() const { return EspRms_; }

        /*! \brief
         * Collect the experimental properties
         *
         * \param[in] bQM
         * \param[in] bZero
         * \param[in] lot      The level of theory used for QM calculation
         * \param[in] gap      Gaussian atom property
         */
        immStatus getExpProps(gmx_bool bQM, gmx_bool bZero,
                              gmx_bool bZPE, const char *lot,
                              const Poldata &pd);
                              
        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fn        A File pointer opened previously.
         * \param[in] iModel    The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] bVerbose  Verobse
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] aps       Gromacs atom properties
         */
        void PrintTopology(const char             *fn,
                           ChargeDistributionModel iModel,
                           bool                    bVerbose,
                           const Poldata          &pd,
                           gmx_atomprop_t          aps,
                           t_commrec              *cr,
                           double                  efield,
                           const char             *lot);

        /*! \brief
         * Print the topology that was generated previously in GROMACS format.
         *
         * \param[in] fn        A File pointer opened previously.
         * \param[in] iModel    The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] bVerbose  Verobse
         * \param[in] pd        Data structure containing atomic properties
         * \param[in] aps       Gromacs atom properties
         */
        void PrintTopology(FILE                    *fp,
                           ChargeDistributionModel  iModel,
                           bool                     bVerbose,
                           const Poldata           &pd,
                           gmx_atomprop_t           aps,
                           bool                     bITP,
                           t_commrec               *cr,
                           double                   efield,
                           const char              *lot);

        /*! \brief
         *  Compute or derive global info about the molecule
         *
         * \param[in] pd   Data structure containing atomic properties
         */
        void CalcQPol(const Poldata &pd, rvec mu);

        void CalcDipole();
        
        void CalcDipole(rvec mu);
        
        void CalcAnisoPolarizability(tensor polar, double *anisoPol);

        /*! \brief
         * Relax the shells (if any) or compute the forces in the molecule
         *
         * \param[in] fplog
         * \param[in] cr
         * \param[in] mu_tot
         */
        void computeForces(FILE *fplog, t_commrec *cr);

        /*! \brief
         * Set the force field
         *
         * \param[in] ff   Force field
         */

        /*! \brief
         * Change the coordinate of the molecule based
         * on the coordinate of the conformation stored
         * in molprop experiment class.
         *
         * \param[in] ei   ExperimentIterator
         */
        void changeCoordinate(ExperimentIterator ei);
        
        bool getOptimizedGeometry(rvec *x);

        void SetForceField(const char *ff) { forcefield_.assign(ff); }

        /*! \brief
         * Update internal structures for bondtype due to changes in pd
         *
         * \param[in] pd      Data structure containing atomic properties
         * \param[in] iType   Interaction type
         */
        void UpdateIdef(const Poldata   &pd,
                        InteractionType  iType);

        /*! \brief
         * Get the force field
         *
         * \param[out] forcefield_     Force field
         */
        std::string getForceField() { return forcefield_; }

        /*! \brief
         * Calculate quadrupole tensor
         */
        void CalcQuadrupole();

        /*! \brief
          CalcQMbasedMoments calculates total dipole moment,
          dipole components, and quadrupoles using QM-based charges like
          Mulliken, Hirshfeld, CM5, etc. Since there is no Shell particle in 
          QM calculations, it loops over eptAtoms, only. 
        */
        void CalcQMbasedMoments(double *q, double *dip, rvec mu, tensor Q);

        /*! \brief
         * Generate Charge Groups
         *
         * \param[in] ecg
         * \param[in] bUsePDBcharge
         */
        immStatus GenerateChargeGroups(eChargeGroup ecg, bool bUsePDBcharge);

        immStatus GenerateGromacs(const gmx::MDLogger &mdlog,
                                  t_commrec           *cr,
                                  const char          *tabfn,
                                  gmx_hw_info_t       *hwinfo);

        /*! \brief
         * Generate cube
         *
         * \param[in] iModel      The distrbution model of charge (e.x. point charge, gaussian, and slater models)
         * \param[in] pd          Data structure containing atomic properties
         * \param[in] spacing     The grid space
         * \param[in] reffn
         * \param[in] pcfn
         * \param[in] pdbdifffn
         * \param[in] potfn
         * \param[in] rhofn
         * \param[in] hisfn
         * \param[in] difffn
         * \param[in] diffhistfn
         * \param[in] oenv
         */
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

        /*! \brief
         * Print the coordinates corresponding to topology after adding shell particles and/or vsites.
         *
         * \param[in] fn A File pointer opened previously.
         */
        void PrintConformation(const char *fn);

        void setInputrec(t_inputrec  *ir)
        {
            inputrec_ = ir;
        }
};

const char *immsg(immStatus imm);

}

#endif
