/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef PLISTWRAPPER_H
#define PLISTWRAPPER_H

//#include <algorithm>
#include <vector>

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"

namespace alexandria
{
//! Interaction type
enum InteractionType 
{
    InteractionType_BONDS  = ebtsBONDS,
    InteractionType_ANGLES = ebtsANGLES,
    InteractionType_PDIHS  = ebtsPDIHS,
    InteractionType_IDIHS  = ebtsIDIHS,
    InteractionType_LINEAR_ANGLES,
    InteractionType_LJ14,
    InteractionType_Polarization,
    InteractionType_CONSTR,
    InteractionType_VSITE2,
};

//! Utility typedef
using ParamIterator = typename std::vector<t_param>::iterator;

//! Cleaner version of plist array
class PlistWrapper
{
    private:
        //! Function type
        int                  ftype_;
        //! Interaction type
        InteractionType      itype_;
        //! Array of parameters
        std::vector<t_param> p_;
    public:
        //! Constructor
        PlistWrapper(InteractionType itype,
                     int             ftype) :ftype_(ftype), itype_(itype) {}
       
        //! Add one parameter
        void addParam(t_param p) { p_.push_back(p); }
        
        //! Return the function type
        int getFtype() const { return ftype_; }

        //! Update the function type
        void setFtype(int ftype) { ftype_ = ftype; }
        
        //! Return the interaction type
        InteractionType interactionType() const { return itype_; }
        
        //! Loop over parameters        
        ParamIterator beginParam() { return p_.begin(); }
        
        //! Loop over parameters        
        ParamIterator endParam() { return p_.end(); }
        
        //! Remove one parameter from the array and return array for next
        ParamIterator eraseParam(ParamIterator p) { return p_.erase(p); }

        //! Remove all parameters
        void eraseParams() { p_.clear(); }
        
        //! Return number of parameters
        unsigned int nParam() { return p_.size(); }
};

//! Another utility typedef for a looper
using  PlistWrapperIterator = typename std::vector<PlistWrapper>::iterator;

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype);

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, InteractionType itype);

unsigned int CountPlist(std::vector<PlistWrapper> &plist, int ftype);

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int                  ftype,
                        const t_param             &p);

void delete_params(std::vector<PlistWrapper> &plist_,
                   const int                  ftype,
                   const int                  alist[]);

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        int                        ftype,
                        InteractionType            itype,
                        const t_param             &p);
}

#endif
