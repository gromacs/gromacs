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
#ifndef GENTOP_VSITE_H
#define GENTOP_VSITE_H

#include <algorithm>

#include "gromacs/gmxpreprocess/grompp.h"
#include "poldata.h"

#define egvtNO           0
#define egvtLINEAR       1
#define egvtPLANAR       2
#define egvtRING_PLANAR  4
#define egvtOUT_OF_PLANE 8
#define egvtALL (egvtLINEAR | egvtPLANAR | egvtRING_PLANAR | egvtOUT_OF_PLANE)

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
        std::vector<gv_planar>     outofplane_;
        //! The ring-planar vsites
        std::vector<gv_ringplanar> ringplanar_;

    public:
        //! Default constructor
        GentopVsites() { gvt_ = egvtNO; };

        //! Default destructor
        ~GentopVsites() {};

        void setVsitesType(unsigned int gvt) { gvt_ = gvt; }

        unsigned int getVsitesType() { return gvt_; }

        bool bHaveVsites()
        {
            return ((linear_.size() > 0) ||
                    (planar_.size() > 0) ||
                    (outofplane_.size() > 0) ||
                    (ringplanar_.size() > 0));
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
        void addOutOfPlane(int a1, int a2, int a3, int a4, int nbonds[]);

        /*! \brief Add a ring-planar vsite
         *
         * param[in] natom number of atoms in the ring
         * param[in] a the ring atoms
         * param[in] nbonds number of bonds for each of the atoms
         */
        void addRingPlanar(int natom, int a[], int nbonds[]);

        //! The fun stuff
        void generateSpecial(bool bGenVsites,
                             t_atoms *atoms, rvec **x,
                             t_params plist[],
                             t_symtab *symtab, gpp_atomtype_t atype,
                             t_excls **excls, gmx_poldata_t pd);
};

}

#endif
