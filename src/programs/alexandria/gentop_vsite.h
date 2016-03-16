/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GENTOP_VSITE_H
#define GENTOP_VSITE_H

#include <vector>

#include "gromacs/gmxpreprocess/grompp.h"

#include "plistwrapper.h"
#include "poldata.h"

struct gpp_atomtype;

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
        GentopVsites(unsigned int gvt) { gvt_ = gvt; }

        //! Default destructor
        ~GentopVsites() {};

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

        /*! \brief Does checks on vsites and more
         *
         * Generate linear angles and merges linear vsites in case
         * there are more than 1 in a row.
         */
        void generateSpecial(const Poldata             &pd,
                             bool                       bUseVsites,
                             t_atoms                   *atoms,
                             rvec                     **x,
                             std::vector<PlistWrapper> &plist,
                             t_symtab                  *symtab,
                             gpp_atomtype              *atype,
                             t_excls                  **excls);
};

}

#endif
