/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_RESP_ATOM_H
#define GMX_RESP_ATOM_H

#include <vector>

#include "gromacs/math/vectypes.h"

#include "poldata.h"

namespace alexandria
{
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
        
        RowZetaQIterator beginRZ() { return rz_.begin(); }

        RowZetaQIterator endRZ() { return rz_.end(); }

        RowZetaQConstIterator beginRZ() const { return rz_.begin(); }

        RowZetaQConstIterator endRZ() const { return rz_.end(); }

 private:
        //! Atom type index
        int                   atype_;
        //! Signify whether this charge should be restrained during fitting
        bool                  bRestrained_;
        //! String corresponding to atom type
        std::string           atomtype_;
        //! Arrays of charge components
        std::vector<RowZetaQ> rz_;
};

//! Iterator over RespAtomType vector
typedef std::vector<RespAtomType>::iterator RespAtomTypeIterator;

//! Const Iterator over RespAtomType vector
typedef std::vector<RespAtomType>::const_iterator RespAtomTypeConstIterator;

class RespAtom
{
    public:
    RespAtom(int atomnumber, int atype, double q, gmx::RVec x) 
     : atomnumber_(atomnumber), atype_(atype), q_(q), x_(x) 
    { qindex_ = -1; }

        int atype() const { return atype_; }

        void setAtype(int i) { atype_ = i; }

        int atomnumber() const { return atomnumber_; }

        void setAtomnumber(int i) { atomnumber_ = i; }

        double charge() const { return q_; }

        void setCharge(double value) { q_ = value; }
        
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
        double    q_;
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
    RespParam(eParm eparm, size_t aindex, size_t zindex) : 
    eparm_(eparm), aindex_(aindex), zindex_(zindex) {}
 
    eParm eParam() const { return eparm_; }
    
    size_t aIndex() const { return aindex_; }
    
    size_t zIndex() const { return zindex_; }
    
 private:
    //! Type of parameter
    eParm  eparm_;
    //! Atom index (for charges) or atype index (for zeta)
    size_t aindex_;
    //! Zeta index (in the RespAtomType) in case we're optimizing zeta
    size_t zindex_;
};

//! Looking for the right Resp parameters? Here to help.
typedef std::vector<RespParam>::iterator RespParamIterator;

} // namespace

#endif
