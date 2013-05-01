#ifndef _MYMOL2_HPP
#define _MYMOL2_HPP

#include "typedefs.h"
#include "vsite.h"
#include "gpp_atomtype.h"
#include "pdb2top.h"
#include "atomprop.h"
#include "gmx_resp.hpp"
#include "gentop_qgen.hpp"
#include "gentop_vsite.hpp"
#include "gentop_core.hpp"
#include "molprop.hpp"
#include "molselect.h"
#include "poldata.h"
#include "gauss_io.hpp"

enum immStatus { 
    immUnknown,
    immOK, immZeroDip, immNoQuad, immCharged, 
    immAtomTypes, immAtomNumber, immMolpropConv, immBondOrder, immRespInit,
    immChargeGeneration,
    immQMInconsistency, immTest, immNoData, immNR 
};

enum eDih { edihNo, edihOne, edihAll, edihNR };

enum eSupport { eSupportNo, eSupportLocal, eSupportRemote, eSupportNR };

namespace alexandria {

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
    int            nshell;
    gmx_bool       *bRing;
    char           **smnames;
    //! Gromacs structures
    int            nexcl_;
    t_excls        *excls;
    int            *symmetric_charges_;
    int            *cgnr_;
    gentop_vsite_t gvt;
    immStatus      immAtoms,immCharges,immTopology;
    std::string    forcefield_;
    bool           bHaveShells_,bHaveVSites_;
    
    //! Determine whether a molecule has symmetry (within a certain tolerance)
    bool IsSymmetric(real toler);
    
    //! Generate Atoms based on quantum calculation with specified level of theory
    immStatus GenerateAtoms(gmx_atomprop_t ap,
                            gmx_poldata_t pd,
                            const char *lot,
                            const char *q_algorithm);

    //! Read atoms?
                            
    //! Generate bonds between atoms
    int MakeBonds(gmx_poldata_t pd,
                  gmx_conect gc,t_params plist[],int nbond[],
                  gmx_bool bH14,gmx_bool bAllDihedrals,
                  gmx_bool bRemoveDoubleDihedrals,
                  gmx_bool bPBC,matrix box,gmx_atomprop_t aps,real tol,
                  gmx_bool bMovePlists);
                            
public:
    rvec           *x,*f,*buf,mu_exp,mu_calc,mu_esp,coq;
    matrix         box;
    real           dip_exp,mu_exp2,dip_err,dip_weight,dip_calc,chieq,Hform,Emol,Ecalc,Force2;
    real           *qESP;
    tensor         Q_exp,Q_calc,Q_esp;
    int            bts[ebtsNR];
    eSupport       eSupp;
    t_state        state;
    t_forcerec     *fr;
    gmx_mtop_t     mtop;
    gmx_localtop_t *ltop;
    gpp_atomtype_t atype;
    gentop_qgen_t  qgen;
    t_symtab       symtab;
    t_inputrec     ir;
    gmx_shellfc_t  shell;
    gmx_enerdata_t enerd;
    gmx_resp_t     gr;
    t_mdatoms      *md;
    t_topology     *topology;
    t_params       plist[F_NRE];
                               
    //! Constructor
    MyMol();
    
    //! Destructor
    ~MyMol();
    
    //! Generate the topology structure
    immStatus GenerateTopology(gmx_atomprop_t ap,
                               gmx_poldata_t pd,
                               const char *lot,
                               const char *q_algorithm,
                               bool bPol,
                               int nexcl);
    //! Generate Charges
    immStatus GenerateCharges(gmx_poldata_t pd,gmx_atomprop_t ap,
                              int iModel,real hfac,real epsr,
                              real qtol,int maxiter,int maxcycle,
                              const char *lot,
                              bool bSymmetricCharges,
                              const char *symm_string);
    
    //! Print the topology that was generated previously in GROMACS format.
    //! fp is a File pointer opened previously.
    void PrintTopology(const char *fn,gmx_poldata_t pd,int iModel,
                       const char *forcefield,bool bVerbose);
    
    //! Print a rtp entry
    void PrintRTPEntry(const char *fn);

    //! Set the force field
    void SetForceField(const char *ff) { forcefield_.assign(ff); }
    
    //! Get the force field
    std::string GetForceField() { return forcefield_; }
    
    void CalcMultipoles();

    void GenerateVsitesShells(gmx_poldata_t pd,bool bGenVsites,bool bAddShells,
                              bool bPairs,eDih edih);
    
    void GenerateChargeGroups(eChargeGroup ecg,bool bUsePDBcharge,
                              const char *ndxfn,int nmol);
    
    void GenerateCube(int iModel,
                      gmx_poldata_t pd,
                      real spacing,
                      const char *reffn,
                      const char *pcfn,
                      const char *pdbdifffn,
                      const char *potfn,
                      const char *rhofn,
                      const char *hisfn,
                      const char *difffn,
                      const char *diffhistfn,
                      output_env_t oenv);
    
    //! Print the coordinates corresponding to topology after adding shell particles
    //! and/or vsites. fp is a File pointer opened previously.
    void PrintConformation(const char *fn);
    
    //! Routine initiating the internal GROMACS structures
    immStatus Initxx(FILE *fp,GaussAtomProp &gap,
                   gmx_bool bQM,char *lot,gmx_bool bZero,
                   gmx_poldata_t pd,gmx_atomprop_t aps,
                   int  iModel,t_commrec *cr,int *nwarn,
                   gmx_bool bCharged,const output_env_t oenv,
                   real th_toler,real ph_toler,
                   real dip_toler,real hfac,gmx_bool bH14,
                   gmx_bool bAllDihedrals,gmx_bool bRemoveDoubleDihedrals,
                   int nexcl,gmx_bool bESP,
                   real watoms,real rDecrZeta,gmx_bool bPol,gmx_bool bFitZeta);
             
};

const char *immsg(immStatus imm);

}

void mv_plists(gmx_poldata_t pd,t_params plist[],gmx_bool bForward);

#define gmx_assert(n,m) if (n != m) { gmx_fatal(FARGS,"Variable %s = %d, should have been %d",#n,n,m); }

#endif
