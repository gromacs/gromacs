/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef PLISTWRAPPER_H
#define PLISTWRAPPER_H

//#include <algorithm>
#include <vector>

#include "gromacs/gmxpreprocess/grompp-impl.h"

namespace alexandria
{
//! Utility typedef
typedef std::vector<t_param>::iterator ParamIterator;
//! Cleaner version of plist array
class PlistWrapper
{
    private:
    //! Function type
        int                  ftype_;
        //! Array of parameters
        std::vector<t_param> p_;
    public:
        //! Constructor
        PlistWrapper(int ftype) :ftype_(ftype) {}
       
        //! Add one parameter
        void addParam(t_param p) { p_.push_back(p); }
        
        //! Return the function type
        int getFtype() { return ftype_; }

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
typedef std::vector<PlistWrapper>::iterator PlistWrapperIterator;

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype);

unsigned int CountPlist(std::vector<PlistWrapper> &plist, int ftype);

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int                  ftype,
                        const t_param             &p);

void delete_params(std::vector<PlistWrapper> &plist_,
                   const int                  ftype,
                   const int                  alist[]);

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int                  ftype,
                        const t_param             &p);
}

#endif
