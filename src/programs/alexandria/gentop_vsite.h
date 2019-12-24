/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#ifndef GENTOP_VSITE_H
#define GENTOP_VSITE_H

#include <vector>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/state.h"

#include "plistwrapper.h"
#include "poldata.h"

struct gpp_atomtype;

namespace alexandria
{
typedef struct {
    int nline; /* Must be 3 or 4 */
    int a[4];
} gv_linear;

typedef struct {
    int a[4];
    int nb[4];
} gv_planar;

typedef struct {
    int natom;
    int a[6];
    int nb[6];
} gv_ringplanar;


/*
             ca                      o
             |                       ||
            bca           e.x.       c
          /     \                  /   \
      bbca1     bbca2             ch3   ch3
 */
 
class gv_inplane
{
    public:
      gv_inplane () {};  
      
      gv_inplane (int natom, int nvsite, int ca, 
                  int bca,   int bbca1,  int bbca2);
      
      int natom() const {return natom_;}
      
      int ca() const {return ca_;}
      
      int bca() const {return bca_;}
      
      int bbca1()  const {return bbca1_;}
      
      int bbca2() const {return bbca2_;}
      
      int nvsite() const {return nvsite_;}
      
    private:    
      //! Number of atoms needed to place virtual sites. (must be 4)
      int natom_;
      //! Number of virtual sites.
      int nvsite_;      
      //! Central atom on which virtual sites are placed.
      int ca_;
      //! Atom bound to the central atom.
      int bca_;
      //! Atom 1 bound to bca_.
      int bbca1_;
      //! Atom 2 bound to bca_.
      int bbca2_;
};

using gv_inplaneIterator      = typename std::vector<gv_inplane>::iterator;
using gv_inplaneConstIterator = typename std::vector<gv_inplane>::const_iterator;
 

/*
           ca                    o
         /    \                /    \
       ca1    ca2     e.x.    c      c
         \ __ /               \\    //
                               c __ c
 */

class gv_outplane
{
    public:
      gv_outplane () {};  
      
      gv_outplane (int natom, int nvsite, int ca, int bca1, int bca2);
      
      int natom() const {return natom_;}
      
      int ca() const {return ca_;}
      
      int bca1() const {return bca1_;}
      
      int bca2() const {return bca2_;}
      
      int nvsite() const {return nvsite_;}
      
    private:    
      //! Number of atoms needed to place virtual sites. (must be 4)
      int natom_;     
      //! Number of virtual sites.
      int nvsite_;      
      //! Central atom on which virtual sites are placed.
      int ca_;
      //! Atom 1 bound to ca_.
      int bca1_;
      //! Atom 2 bound to ca_.
      int bca2_;
};

using gv_outplaneIterator      = typename std::vector<gv_outplane>::iterator;
using gv_outplaneConstIterator = typename std::vector<gv_outplane>::const_iterator;

class GentopVsites
{
    private:
        //! The type of vsites we are storing
        unsigned int               gvt_;
        //! The linear vsites
        std::vector<gv_linear>     linear_;
        //! The planar vsites
        std::vector<gv_planar>     planar_;
        //! The out-of-plane vsites
        std::vector<gv_outplane>   outplane_;
        //! The out-of-plane vsites
        std::vector<gv_inplane>    inplane_;
        //! The ring-planar vsites
        std::vector<gv_ringplanar> ringplanar_;

    public:
        //! Default constructor
        GentopVsites(unsigned int gvt) { gvt_ = gvt; }

        unsigned int getVsitesType() { return gvt_; }

        bool bHaveVsites()
        {
            return ((inplane_.size()  > 0) ||
                    (outplane_.size() > 0));
        }

        /*! \brief Add a linear vsite
         *
         * param[in] a1 first atom
         * param[in] a2 second atom
         * param[in] a3 atom number of the vsite
         */
        void addLinear(int a1, int a2, int a3);

        /*! \brief Merge linear sites
         *
         * \param[in] bGenVsites
         */
        void mergeLinear(bool bGenVsites);

        /*! \brief Add a planar vsite
         *
         * param[in] a1 first atom
         * param[in] a2 second atom
         * param[in] a3 third atom
         * param[in] a4 atom number of the vsite
         * param[in] nbonds number of bonds for each of the atoms
         */
        void addPlanar(int a1, int a2, int a3, int a4, int nbonds[]);
        
        /*! \brief Add an out of plane vsite
         *
         * param[in] a1 first atom
         * param[in] a2 second atom
         * param[in] a3 third atom
         * param[in] a4 atom number of the vsite
         * param[in] nbonds number of bonds for each of the atoms
         */
        void addOutPlane(int natom, int nvsite, int ca, int bca1, int bca2);
        
        gv_outplaneIterator outplaneBegin() {return outplane_.begin(); }

        gv_outplaneConstIterator outplaneBegin() const {return outplane_.begin(); }

        gv_outplaneIterator outplaneEnd() {return outplane_.end(); }

        gv_outplaneConstIterator outplaneEnd() const {return outplane_.end(); }
        
        gv_outplaneIterator findOutPlane(int natom, int nvsite, int ca, int bca1, int bca2);
        
        gv_outplaneConstIterator findOutPlane(int natom, int nvsite, 
                                              int ca,    int bca1, int bca2) const;
                                              
        gv_outplaneIterator findOutPlane(int nvsite, int ca);
        
        gv_outplaneConstIterator findOutPlane(int nvsite, int ca) const;
        
        void addInPlane(int natom, int nvsite, int ca, 
                        int bca,   int bbca1, int bbca2);
        
        gv_inplaneIterator inplaneBegin() {return inplane_.begin(); }

        gv_inplaneConstIterator inplaneBegin() const {return inplane_.begin(); }

        gv_inplaneIterator inplaneEnd() {return inplane_.end(); }

        gv_inplaneConstIterator inplaneEnd() const {return inplane_.end(); }
        
        gv_inplaneIterator findInPlane(int natom, int nvsite, int ca, 
                                       int bca,   int bbca1,  int bbca2);
        
        gv_inplaneConstIterator findInPlane(int natom, int nvsite, int ca, 
                                            int bca,   int bbca1,  int bbca2) const;
                                            
        gv_inplaneIterator findInPlane(int nvsite, int ca);
        
        gv_inplaneConstIterator findInPlane(int nvsite, int ca) const;
        
        int nVsites();
        
        void gen_Vsites(const Poldata             *pd,
                        t_atoms                   *atoms,
                        std::vector<PlistWrapper> &plist,
                        gpp_atomtype              *atype,
                        t_symtab                  *symtab,
                        t_excls                   **excls,
                        t_state                    *state);

        /*! \brief Add a ring-planar vsite
         *
         * param[in] natom number of atoms in the ring
         * param[in] a the ring atoms
         * param[in] nbonds number of bonds for each of the atoms
         */
        void addRingPlanar(int natom, int a[], int nbonds[]);

        /*! \brief Does checks on vsites and more
         *
         * Generate linear angles and merges linear vsites in case
         * there are more than 1 in a row.
         */
        void generateSpecial(const Poldata             *pd,
                             bool                       bUseVsites,
                             t_atoms                   *atoms,
                             rvec                     **x,
                             std::vector<PlistWrapper> &plist,
                             t_symtab                  *symtab,
                             gpp_atomtype              *atype,
                             t_excls                  **excls,
                             t_state                   *state);
};

}

#endif
