/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_RESP_ATOM_H
#define GMX_RESP_ATOM_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

#include "poldata.h"

namespace alexandria
{
class RespZeta
{
 public:
    RespZeta(int row, real q, real zeta) : row_(row), q_(q), zeta_(zeta), 
        zetaRef_(zeta) { zindex_ = -1; }
    
    int row() const { return row_; };
    
    real q() const { return q_; }
    
    void setQ(real q) { q_ = q; }
    
    real zeta() const { return zeta_; } 
    
    void setZeta(double z) { zeta_ = z; }
    
    real zetaRef() const { return zetaRef_; } 
    
    void setZetaRef(double z) { zetaRef_ = z; }
    
    int zIndex() const { return zindex_; }
    
    void setZindex(int zi) { zindex_ = zi; }
    
 private:
    //! The row in the periodic table for each of the charge components
    int  row_;
    //! Charge of each of the components
    real q_;
    //! Inverse screening length of each of the components
    real zeta_;
    //! Reference (starting) value for zeta
    real zetaRef_;
    //! Parameter optimization index
    int  zindex_;
};
//! Loop over RespZeta
typedef std::vector<RespZeta>::iterator RespZetaIterator;

//! Loop over RespZeta
typedef std::vector<RespZeta>::const_iterator RespZetaConstIterator;

class RespAtomType
{
    public:
        RespAtomType(int atype,
                     const char *atomtype, 
                     const Poldata &pd,
                     ChargeDistributionModel iDistributionModel,
                     const std::vector<std::string> &dzatoms);

        size_t getNZeta() const { return rz_.size(); }

        const std::string &getAtomtype() const { return atomtype_; }

        void setAtomtype(const std::string &type) { atomtype_ = type; }
        
        bool getBRestrained() const { return bRestrained_; }

        int getAtype() const { return atype_; }

        void setAtype(int i) { atype_ = i; }
        
        RespZetaIterator beginRZ() { return rz_.begin(); }

        RespZetaIterator endRZ() { return rz_.end(); }

        RespZetaConstIterator beginRZ() const { return rz_.begin(); }

        RespZetaConstIterator endRZ() const { return rz_.end(); }

 private:
        //! Atom type index
        int                   atype_;
        //! Signify whether this charge should be restrained during fitting
        bool                  bRestrained_;
        //! String corresponding to atom type
        std::string           atomtype_;
        //! Arrays of charge components
        std::vector<RespZeta> rz_;
};

//! Iterator over RespAtomType vector
typedef std::vector<RespAtomType>::iterator RespAtomTypeIterator;

//! Const Iterator over RespAtomType vector
typedef std::vector<RespAtomType>::const_iterator RespAtomTypeConstIterator;

class RespAtom
{
    public:
    RespAtom(int atomnumber, int atype, real q, gmx::RVec x) 
     : atomnumber_(atomnumber), atype_(atype), q_(q), x_(x) 
    { qindex_ = -1; }

        int atype() const { return atype_; }

        void setAtype(int i) { atype_ = i; }

        int atomnumber() const { return atomnumber_; }

        void setAtomnumber(int i) { atomnumber_ = i; }

        real charge() const { return q_; }

        void setCharge(real value) { q_ = value; }
        
        const gmx::RVec &x() const { return x_; }
    
        void setX(const gmx::RVec &x) { x_ = x; }
        
        int qIndex() const { return qindex_; }
        
        void setQindex(int qi) { qindex_ = qi; }
    private:
        //! Atomic number
        int       atomnumber_;
        //! Atom type
        int       atype_;
        //! Total charge of the atom (which is optimized in resp)
        real      q_;
        //! Coordinates for this atom
        gmx::RVec x_;
        //! Index in parameterization array
        int       qindex_;
};
//! Let's loop over RespAtoms, shall we?
typedef std::vector<RespAtom>::iterator RespAtomIterator;

//! Let's loop over RespAtoms, shall we?
typedef std::vector<RespAtom>::const_iterator RespAtomConstIterator;

//! Optimizable entities
enum eParm {
    //! Optimize charge
    eparmQ, 
    //! Optimize zeta
    eparmZ 
};

class RespParam
{
 public:
    RespParam(eParm eparm, int aindex, int zindex) : 
    eparm_(eparm), aindex_(aindex), zindex_(zindex) {}
 
    eParm eParam() const { return eparm_; }
    
    int aIndex() const { return aindex_; }
    
    int zIndex() const { return zindex_; }
    
 private:
    //! Type of parameter
    eParm eparm_;
    //! Atom index (for charges) or atype index (for zeta)
    int   aindex_;
    //! Zeta index (in the RespAtomType) in case we're optimizing zeta
    int   zindex_;
};

} // namespace

#endif
