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
/*! Processing of input files and force fields
 * \ingroup group_preprocessing
 * \brief
 * Provides tools for input processing based on the alexandria force field
 *
 * The tools in the alexandria directory are under heavy development still.
 * Once finished they will provide processing of small molecules files based
 * on quantum chemistry or other input source through the OpenBabel library.
 * Assigning of atom types, derivation of charges, and generation of
 * molecular topology files (or objects) is implemented.
 *
 * \author David van der Spoel <david.vanderspoel@gmail.com>
 * \ingroup module_alexandria
 */
#ifndef MOLPROP_H
#define MOLPROP_H

#include <string.h>

#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/real.h"

#include "phase.h"
#include "poldata.h"

/*! \brief
 * Enumerated type holding the types of observables stored in MolProp
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum MolPropObservable {
    MPO_POTENTIAL,
    MPO_DIPOLE,
    MPO_QUADRUPOLE,
    MPO_POLARIZABILITY,
    MPO_ENERGY,
    MPO_ENTROPY,
    MPO_NR
};

//! Enum to select either QM or Experimental data or either
enum iqmType {
    iqmExp, iqmBoth, iqmQM, iqmNR
};

//! Strings describing the MolPropObservable enum elements
extern const char *mpo_name[MPO_NR];

//! Strings describing the MolPropObservable enum units
extern const char *mpo_unit[MPO_NR];

/*! \brief
 * Enumerated type holding the result status of communication operations
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum CommunicationStatus {
    CS_OK        = 6666,
    CS_ERROR     = 7777,
    CS_SEND_DATA = 8888,
    CS_RECV_DATA = 9999
};

//! String describing the CommunicationStatus enum elements
extern const char *cs_name(CommunicationStatus cs);

/*! \brief
 * Contains all classes related to alexandria force field tools
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
namespace alexandria
{
/*! \brief
 * Specifies the name of an atom type and the number in a molecular composition
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class AtomNum
{
    private:
        std::string _catom;
        int         _cnumber;
    public:
        //! Empty constructor
        AtomNum() {};

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(const char *catom, int cnumber) { SetAtom(catom); SetNumber(cnumber); }

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(std::string catom, int cnumber) { SetAtom(catom); SetNumber(cnumber); }

        //! Return the name of the atom for this AtomNum
        const std::string &getAtom() const { return _catom; }

        //! Set the name of the atom for this AtomNum
        void SetAtom(std::string catom) { _catom = catom; }

        //! Set the name of the atom for this AtomNum
        void SetAtom(const char *catom) { _catom.assign(catom); }

        //! Return the number of atoms for this AtomNum
        int getNumber() const { return _cnumber; }

        //! Set the number of atoms for this AtomNum
        void SetNumber(int cnumber) { _cnumber = cnumber; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over a vector of AtomNum
typedef std::vector<AtomNum>::iterator AtomNumIterator;

/*! \brief
 * Contains the molecular composition in terms of atoms
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularComposition
{
    private:
        std::string          _compname;
        std::vector<AtomNum> _atomnum;
    public:
        //! Empty constructor
        MolecularComposition() {}

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(const char *compname) { _compname.assign(compname); }

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(std::string compname) { _compname = compname; }

        //! Return the composition name
        const std::string &getCompName() const { return _compname; }

        //! Set the composition name
        void SetCompName(std::string compname) { _compname = compname; }

        //! Set the composition name
        void SetCompName(char *compname) { _compname.assign(compname); }

        //! Add an AtomNum struct to the composition
        void AddAtom(AtomNum an);

        //! Remove the atom with name catom from the composition
        void DeleteAtom(const char *catom)
        {
            std::string _str(catom); DeleteAtom(_str);
        }

        //! Remove the atom with name catom from the composition
        void DeleteAtom(std::string catom);

        //! Replace the oldatom by newatom
        void ReplaceAtom(const char *oldatom, const char *newatom)
        {
            std::string so(oldatom), sn(newatom); ReplaceAtom(so, sn);
        }
        //! Replace the oldatom by newatom
        void ReplaceAtom(std::string oldatom, std::string newatom);

        //! Return iterator to begin looping over AtomNum
        AtomNumIterator BeginAtomNum() { return _atomnum.begin(); }

        //! Return iterator to end looping over AtomNum
        AtomNumIterator EndAtomNum() { return _atomnum.end(); }

        //! Return iterator pointing to a specific atom or EndAtomNum if not found
        AtomNumIterator SearchAtom(std::string an);

        //! Return the number of atoms of a certain type
        int CountAtoms(const char *atom);

        //! Return the number of atoms of a certain type
        int CountAtoms(std::string atom);

        //! Return the total number of atoms
        int CountAtoms();

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularComposition items
typedef std::vector<MolecularComposition>::iterator MolecularCompositionIterator;

/*! \brief
 * Generic molecular property base clase
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class GenericProperty
{
    private:
        std::string type_, unit_;
        ePhase      eP_;
        double      T_;
    public:
        //! Empty constructor
        GenericProperty() { T_ = 0; eP_ = epNR; };

        /*! \brief
         * Creates a new GenericProperty object.
         *
         * \param[in] type  Type of the property
         * \param[in] unit  Unit of the property
         * \param[in] T     Temperature
         */
        GenericProperty(std::string type, std::string unit, double T, ePhase ep)
        { SetType(type); SetUnit(unit); setTemperature(T); setPhase(ep); }

        //! Return the property type
        const std::string &getType() const { return type_; }

        //! Return the unit of the property
        const std::string &getUnit() const { return unit_; }

        //! Return the temperature
        double getTemperature() const { return T_; }

        //! Return the phase
        ePhase getPhase() const { return eP_; }

        //! Set the type of the property
        void SetType(std::string type);

        //! Set the unit of the property
        void SetUnit(std::string unit);

        //! Set the temperature of the property
        void setTemperature(double T) { T_ = T; }

        void setPhase(ePhase ep) { eP_ = ep; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};

/*! \brief
 * Contains the elements of the molecular quadrupole
 *
 * The six elements of the upper diagonal of a quadrupole tensor are stored.
 * The values are dependent on the orientation of the molecule and are relative
 * to the center of charge in case that the total charge or total dipole of
 * the molecules are non-zero.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularQuadrupole : public GenericProperty
{
    private:
        double xx_, yy_, zz_, xy_, xz_, yz_;
    public:
        //! Empty constructor
        MolecularQuadrupole() {}

        //! Constructor initiating all elements of the quadrupole tensor
        MolecularQuadrupole(std::string type, std::string unit, double T,
                            double xx, double yy, double zz,
                            double xy, double xz, double yz) : GenericProperty(type, unit, T, epGAS) { Set(xx, yy, zz, xy, xz, yz); };

        //! Set all the elements of the qudrupole tensor
        void Set(double xx, double yy, double zz, double xy, double xz, double yz) { xx_ = xx; yy_ = yy; zz_ = zz; xy_ = xy; xz_ = xz; yz_ = yz; };

        //! get all the elements of the qudrupole tensor
        void get(double *xx, double *yy, double *zz, double *xy, double *xz, double *yz) { *xx = xx_; *yy = yy_; *zz = zz_; *xy = xy_; *xz = xz_; *yz = yz_; };

        //! Return the XX component of the quadrupole tensor
        double getXX() const { return xx_; }

        //! Return the YY component of the quadrupole tensor
        double getYY() const { return yy_; }

        //! Return the ZZ component of the quadrupole tensor
        double getZZ() const { return zz_; }

        //! Return the XY component of the quadrupole tensor
        double getXY() const { return xy_; }

        //! Return the XZ component of the quadrupole tensor
        double getXZ() const { return xz_; }

        //! Return the YZ component of the quadrupole tensor
        double getYZ() const { return yz_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularQuadrupole items
typedef std::vector<MolecularQuadrupole>::iterator MolecularQuadrupoleIterator;

/*! \brief
 * Contains the elements of the molecular polarizability tensor
 *
 * The six elements of the upper diagonal of a polarizability tensor are stored
 * along with the average molecular polarizability and the error if known.
 * The values are dependent on the orientation of the molecule.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularPolarizability : public GenericProperty
{
    private:
        double xx_, yy_, zz_, xy_, xz_, yz_, _average, _error;
    public:
        //! Empty constructor
        MolecularPolarizability() {}

        //! Constructor initiating all elements of the quadrupole tensor
        MolecularPolarizability(std::string type, std::string unit, double T,
                                double xx, double yy, double zz,
                                double xy, double xz, double yz,
                                double average, double error) : GenericProperty(type, unit, T, epGAS) { Set(xx, yy, zz, xy, xz, yz, average, error); };

        //! Set all the elements of the polarizability tensor
        void Set(double xx, double yy, double zz,
                 double xy, double xz, double yz,
                 double average, double error);

        //! get all the elements of the polarizability tensor
        void get(double *xx, double *yy, double *zz,
                 double *xy, double *xz, double *yz,
                 double *average, double *error) const
        {
            *xx      = xx_; *yy = yy_; *zz = zz_;
            *xy      = xy_; *xz = xz_; *yz = yz_;
            *average = _average; *error = _error;
        }

        //! Return the XX component of the polarizability tensor
        double getXX() const { return xx_; }

        //! Return the YY component of the polarizability tensor
        double getYY() const { return yy_; }

        //! Return the ZZ component of the polarizability tensor
        double getZZ() const { return zz_; }

        //! Return the XY component of the polarizability tensor
        double getXY() const { return xy_; }

        //! Return the XZ component of the polarizability tensor
        double getXZ() const { return xz_; }

        //! Return the YZ component of the polarizability tensor
        double getYZ() const { return yz_; }

        //! Return the average of the polarizability tensor
        double getAverage() const { return _average; }

        //! Return the error in the polarizability tensor
        double getError() const { return _error; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularPolarizability items
typedef std::vector<MolecularPolarizability>::iterator MolecularPolarizabilityIterator;

/*! \brief
 * Contains a molecular energy
 *
 * Different energy terms associated with a molecule can be stored based
 * on quantum chemistry calculations or experimental data.
 * For example the Enthalpy of formation at different temperatures.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularEnergy : public GenericProperty
{
    private:
        double _value, _error;
    public:
        //! Empty constructor
        MolecularEnergy() {};

        //! Constructor storing all properties related to this energy term
        MolecularEnergy(std::string type, std::string unit, double T, ePhase ep, double value, double error) : GenericProperty(type, unit, T, ep) { Set(value, error); };

        //! Set the value and error for the energy
        void Set(double value, double error) { _value = value; _error = error; };

        //! get the value and error for this energy
        void get(double *value, double *error) const { *value = _value; *error = _error; };

        //! Set the value for the energy
        void SetValue(double value) { _value = value; };

        //! Return the energy value
        double getValue() const { return _value; };

        //! Set the error in the energy value
        void SetError(double error) { _value = error; };

        //! Return the error in the energy
        double getError() const { return _error; };

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularEnergy items
typedef std::vector<MolecularEnergy>::iterator MolecularEnergyIterator;

/*! \brief
 * Contains the dipole vector
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularDipole : public GenericProperty
{
    private:
        double _x, _y, _z;
        double _aver, _error;
    public:
        //! Empty constructor
        MolecularDipole() {}

        //! Constructor storing all properties related to this dipole
        MolecularDipole(std::string type, std::string unit, double T,
                        double x, double y, double z, double aver, double error) : GenericProperty(type, unit, T, epGAS) { Set(x, y, z, aver, error); }

        //! Set all properties related to this dipole
        void Set(double x, double y, double z, double aver, double error) { _x = x; _y = y; _z = z; _aver = aver; _error = error; };

        //! Return all properties of this dipole
        void get(double *x, double *y, double *z, double *aver, double *error) const
        { *x = _x; *y = _y; *z = _z; *aver = _aver; *error = _error; };

        //! Return the average dipole value
        double getAver() const { return _aver; }

        //! Return the error in the average dipole
        double getError() const { return _error; }

        //! Return the X component of the dipole
        double getX() const { return _x; }

        //! Return the Y component of the dipole
        double getY() const { return _y; }

        //! Return the Z component of the dipole
        double getZ() const { return _z; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over a vector of MolecularDipole
typedef std::vector<MolecularDipole>::iterator MolecularDipoleIterator;

/*! \brief
 * Contains the electrostatic potential in a coordinate close to a molecule.
 *
 * The electrostatic potential (ESP) can be computed using quantum chemistry and
 * stored in this class.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ElectrostaticPotential
{
    private:
        std::string xyzUnit_, vUnit_;
        int         espID_;
        double      x_, y_, z_, V_;
    public:
        //! Empty constructor
        ElectrostaticPotential() {}

        //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        ElectrostaticPotential(std::string xyz_unit, std::string V_unit, int espid, double x, double y, double z, double V) { Set(xyz_unit, V_unit, espid, x, y, z, V); };

        //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        ElectrostaticPotential(const char *xyz_unit, const char *V_unit, int espid, double x, double y, double z, double V) { Set(xyz_unit, V_unit, espid, x, y, z, V); };

        //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Set(std::string xyz_unit, std::string V_unit, int espid, double x, double y, double z, double V) { xyzUnit_ = xyz_unit; vUnit_ = V_unit; espID_ = espid; x_ = x; y_ = y; z_ = z; V_ = V; };

        //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Set(const char *xyz_unit, const char *V_unit, int espid, double x, double y, double z, double V)
        {
            std::string x_(xyz_unit), V_(V_unit); Set(x_, V_, espid, x, y, z, V);
        }

        //! Return the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void get(char **xyz_unit, char **V_unit, int *espid, double *x, double *y, double *z, double *V) const
        {
            *xyz_unit = strdup(xyzUnit_.c_str()); *V_unit = strdup(vUnit_.c_str()); *espid = espID_; *x = x_; *y = y_; *z = z_; *V = V_;
        };

        //! Return the unit of the coordinates
        const std::string &getXYZunit() const { return xyzUnit_; }

        //! Return the unit of the potential
        const std::string &getVunit() const { return vUnit_; }

        //! Return the ESP id from the original calculation
        int getEspid() const { return espID_; }

        //! Return the X coordinate of the ESP point
        double getX() const { return x_; }

        //! Return the Y coordinate of the ESP point
        double getY() const { return y_; }

        //! Return the Z coordinate of the ESP point
        double getZ() const { return z_; }

        //! Return the electrostatic potential at this point in space
        double getV() const { return V_; }

        //! Set the potential for this instance of the class
        void SetV(double V) { V_ = V; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over ElectrostaticPotential items
typedef std::vector<ElectrostaticPotential>::iterator ElectrostaticPotentialIterator;

/*! \brief
 * Chemical bond in a molecule with associated bond order.
 *
 * The chemical bonds in a molecule are stored here along with the bond order
 * which together can be used for determining the atom types in a force field
 * calculation. In the present implementation this data is generated by OpenBabel.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Bond
{
    private:
        int _ai, _aj, _bondorder;
    public:
        //! Empty constructor
        Bond() {}

        //! Constructor setting the ids of the atoms and the bondorder
        Bond(int ai, int aj, int bondorder) { Set(ai, aj, bondorder); }

        //! Sets the ids of the atoms and the bondorder
        void Set(int ai, int aj, int bondorder) {_ai = ai; _aj = aj; _bondorder = bondorder; };

        //! Returns the ids of the atoms and the bondorder
        void get(int *ai, int *aj, int *bondorder) const
        { *ai = _ai; *aj = _aj; *bondorder = _bondorder; };

        //! Returns the first atom id
        int getAi() const { return _ai; }

        //! Returns the second atom id
        int getAj() const { return _aj; }

        //! Returns the bondorder
        int getBondOrder() const { return _bondorder; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over Bond items
typedef std::vector<Bond>::iterator BondIterator;

/*! \brief
 * Contains the charge of an atom
 *
 * The charge of an atom can be derived from quantum chemistry codes in different
 * ways, despite it not being a physical observable.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class AtomicCharge : public GenericProperty
{
    private:
        double _q;
    public:
        //! Empty constructor
        AtomicCharge() {}

        //! Constructor setting type, unit and charge itself
        AtomicCharge(std::string type, std::string unit, double T, double q) : GenericProperty(type, unit, T, epGAS) { SetQ(q); };

        //! Set the charge to q
        void SetQ(double q) { _q = q; };

        //! Return the charge
        double getQ() const { return _q; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over AtomicCharge items
typedef std::vector<AtomicCharge>::iterator AtomicChargeIterator;

/*! \brief
 * Contains data on an atom based on a calculation.
 *
 * This class coordinates, name, atom type and an array of
 * AtomicCharge values based on different methods.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class CalcAtom
{
    private:
        std::string               name_, obType_, unit_;
        double                    x_, y_, z_;
        int                       atomID_;
        std::vector<AtomicCharge> q_;
    public:
        //! Empty constructor
        CalcAtom() {}

        //! Constructor initiating the name, type and atomid
        CalcAtom(const char *name, const char *obtype, int atomid)
        {
            name_.assign(name); obType_.assign(obtype); atomID_ = atomid;
        };

        //! Constructor initiating the name, type and atomid
        CalcAtom(std::string name, std::string obtype, int atomid)
        {
            name_ = name; obType_ = obtype; atomID_ = atomid;
        };


        //! Function returning true if the two atoms are equal
        bool Equal(CalcAtom ca);

        //! Add an AtomicCharge element to the atom
        void AddCharge(AtomicCharge aq);

        //! Begin Iterator over AtomicCharge items
        AtomicChargeIterator BeginQ() { return q_.begin(); }

        //! End Iterator over AtomicCharge items
        AtomicChargeIterator EndQ() { return q_.end(); }

        //! Return the atom id of the atom
        int getAtomid() const { return atomID_; }

        //! Return the name of the atom
        const std::string &getName() const { return name_; }

        //! Return the OpenBabel type of the atom
        const std::string &getObtype() const { return obType_; }

        //! Return the unit of the coordinates of the atom
        const std::string &getUnit() const { return unit_; }

        //! Set the unit of the coordinates of the atom
        void SetUnit(std::string unit);

        //! Set the unit of the coordinates of the atom
        void SetUnit(const char *unit) { std::string s(unit); SetUnit(s); }

        //! Set the coordinates of the atom
        void SetCoords(double x, double y, double z) { x_ = x; y_ = y; z_ = z; }

        //! Return all the coordinates of the atom
        void getCoords(double *x, double *y, double *z) const
        { *x = x_; *y = y_; *z = z_; }

        //! Return the X coordinate of the atom
        double getX() const { return x_; }

        //! Return the Y coordinate of the atom
        double getY() const { return y_; }

        //! Return the Z coordinate of the atom
        double getZ() const { return z_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over CalcAtom items
typedef std::vector<CalcAtom>::iterator CalcAtomIterator;

enum DataSource {
    dsExperiment, dsTheory
};

const char *dataSourceName(DataSource ds);

DataSource dataSourceFromName(const std::string &name);

/*! \brief
 * Contains molecular data based on experiments
 *
 * This is a composite class holding the results from experiment, either
 * the dipole, polarizability, energy or the quadrupole. A reference to the
 * publication (or handbook) containing the data is stored, as well as the
 * conformation of the molecule (if known).
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Experiment
{
    private:
        DataSource                           dataSource_;
        std::string                          reference_, conformation_, jobtype_;
        std::string                          _program, _method, _basisset, _datafile;
        std::vector<CalcAtom>                _catom;
        std::vector<ElectrostaticPotential>  _potential;
        std::vector<MolecularDipole>         dipole_;
        std::vector<MolecularEnergy>         energy_;
        std::vector<MolecularQuadrupole>     quadrupole_;
        std::vector<MolecularPolarizability> polar_;

    public:
        //! Empty constructor
        Experiment() { }

        //! Constructor initiating an Experiment with reference and conformation
        Experiment(std::string reference, std::string conformation) :
            dataSource_(dsExperiment), reference_(reference),
            conformation_(conformation)
        {}

        //! Constructor initiating a Calculation
        Experiment(std::string program, std::string method,
                   std::string basisset, std::string reference,
                   std::string conformation, std::string datafile,
                   std::string jobtype) :
            dataSource_(dsTheory),
            reference_(reference), conformation_(conformation), jobtype_(jobtype),
            _program(program), _method(method), _basisset(basisset),
            _datafile(datafile)
        {}

        //! Return the type of data
        DataSource dataSource() const { return dataSource_; }

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

        //! Add a MolecularDipole element containing polarizability data
        void AddPolar(MolecularPolarizability mdp) { polar_.push_back(mdp); }

        //! Return Begin Iterator over polarizability elements
        MolecularPolarizabilityIterator BeginPolar() { return polar_.begin(); }

        //! Return End Iterator over polarizability elements
        MolecularPolarizabilityIterator EndPolar()   { return polar_.end(); }

        //! Number of polarizability values
        int NPolar() { return polar_.size(); }

        //! Number of dipole values
        int NDipole() { return dipole_.size(); }

        //! Add a MolecularDipole element containing dipole data
        void AddDipole(MolecularDipole mdp) { dipole_.push_back(mdp); }

        //! Return Begin Iterator over dipole elements
        MolecularDipoleIterator BeginDipole() { return dipole_.begin(); }

        //! Return End Iterator over dipole elements
        MolecularDipoleIterator EndDipole()   { return dipole_.end(); }

        //! Number of quadrupole values
        int NQuadrupole() { return quadrupole_.size(); }

        //! Add a MolecularQuadrupole element
        void AddQuadrupole(MolecularQuadrupole mq) { quadrupole_.push_back(mq); }

        //! Return Begin Iterator over quadrupole elements
        MolecularQuadrupoleIterator BeginQuadrupole() { return quadrupole_.begin(); }

        //! Return End Iterator over quadrupole elements
        MolecularQuadrupoleIterator EndQuadrupole() { return quadrupole_.end(); }

        //! Add a MolecularEnergy element
        void AddEnergy(MolecularEnergy me) { energy_.push_back(me); }

        //! Return Begin Iterator over energy elements
        MolecularEnergyIterator BeginEnergy() { return energy_.begin(); }

        //! Return End Iterator over energy elements
        MolecularEnergyIterator EndEnergy()   { return energy_.end(); }

        //! Return the molecular conformation
        const std::string &getConformation() const { return conformation_; }

        //! Return the literature reference
        const std::string &getReference() const { return reference_; }

        //! Return the type of calculation
        const std::string &getJobtype() const { return jobtype_; }

        //! Add a CalcAtom object to the list of atoms
        void AddAtom(CalcAtom ca);

        //! Return the number of atoms
        int NAtom() { return _catom.size(); }

        //! Iterator Begin over CalcAtom objects
        CalcAtomIterator BeginAtom() { return _catom.begin(); }

        //! Iterator End over CalcAtom objects
        CalcAtomIterator EndAtom() { return _catom.end(); }

        //! Return iterator pointint to this particular atom or EndAtom() if not found
        CalcAtomIterator SearchAtom(CalcAtom ca);

        //! Add ElectrostaticPotential element to the array
        void AddPotential(ElectrostaticPotential ep) { _potential.push_back(ep); }

        //! Return the number of potential points
        int NPotential() { return _potential.size(); };

        //! Iterator Begin over ElectrostaticPotential objects
        ElectrostaticPotentialIterator BeginPotential() { return _potential.begin(); }

        //! Iterator End over ElectrostaticPotential objects
        ElectrostaticPotentialIterator EndPotential() { return _potential.end(); }

        //! Return the program used to perform the calculation
        const std::string &getProgram() const { return _program; }

        //! Return the basis set used to perform the calculation
        const std::string &getBasisset() const { return _basisset; }

        //! Return the method used to perform the calculation
        const std::string &getMethod() const { return _method; }

        //! Return the datafile from which the calculation output was extracted
        const std::string &getDatafile() const { return _datafile; }

        /*! \brief
         * Convenience function that fetches a value from this experiment
         *
         * \param[in]  type   The type of the data (dependent on whether it is dipole, energy etc.)
         * \param[in]  mpo    Enum selecting the type of data to fetch
         * \param[out] value  The value of e.g. the energy
         * \param[out] error  The error in the value of e.g. the energy
         * \param[out] T      Temperature
         * \param[out] vec    Vector data to be output, dipole or diagonal element of polarizability
         * \param[out] quadrupole The quadrupole tensor
         * \return true on success
         */
        bool getVal(const char *type, MolPropObservable mpo,
                    double *value, double *error, double *T,
                    double vec[3], tensor quadrupole);

        bool getHF(double *value);

        //! Merge in another object. Return number of warnings.
        int Merge(std::vector<Experiment>::iterator src);

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over Experiment items
typedef std::vector<Experiment>::iterator ExperimentIterator;

/*! \brief
 * Contains molecular properties from a range of sources.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolProp
{
    private:
        int                               _index;
        double                            _mass;
        int                               _charge, _multiplicity;
        std::string                       _formula, _texform, _molname, _iupac, _cas, _cid, _inchi;
        std::vector<std::string>          category_;
        std::vector<MolecularComposition> _mol_comp;
        std::vector<Experiment>           _exper;
        std::vector<Bond>                 _bond;
    public:
        //! Construct a number MolProp object
        MolProp() { _index = -1; _mass = 0; _charge = 0; _multiplicity = 0; }

        /*! \brief
         * Check the internal consistency of this object
         *
         * \todo Implement this
         */
        void CheckConsistency();

        //! Set the index number for sorting
        void SetIndex(int index) { _index = index; }

        //! Return the index number for sorting
        int getIndex() const { return _index; }

        //! Set the molecular mass
        void SetMass(double mass) { _mass = mass; }

        //! Return the molecular mass
        double getMass() const { return _mass; }

        //! Set the total charge of the molecule
        void SetCharge(double charge) { _charge = charge; }

        //! Return the total charge of the molecule
        int getCharge() const { return _charge; }

        //! Set the multiplicity of the molecule
        void SetMultiplicity(int multiplicity) { _multiplicity = multiplicity; }

        //! Return the multiplicity  of the molecule
        int getMultiplicity() const { return _multiplicity; }

        /*! \brief
         * Merge the content of another MolProp into this one
         *
         * \param[in] mpi The object to be merged into the present one
         * \return Number of warnings
         * \todo Check and double check
         */
        int Merge(std::vector<MolProp>::iterator mpi);

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

        //! Set the LaTeX formula
        void SetTexFormula(const std::string &formula) { _texform.assign(formula); }

        //! Set the formula
        void SetFormula(const std::string &formula) { _formula.assign(formula); }

        //! Return the formula
        const std::string &formula() const { return _formula; }

        //! Return the LaTeX formula
        const std::string &getTexFormula() const;

        /*! \brief
         * Generate the chemical formula for this molecule based on atoms
         * present in a calculation
         *
         * \param[in] ap Data structure containing information about atoms
         * \todo Check and double check. If there is no calculation data no
         * formula can be generated
         */
        bool GenerateFormula(gmx_atomprop_t ap);

        //! Set the molname
        void SetMolname(const char *molname) { _molname.assign(molname); }

        //! Set the molname
        void SetMolname(std::string molname) { _molname = molname; }

        //! Return the molname
        const std::string &getMolname() const { return _molname; }

        //! Set the IUPAC name
        void SetIupac(const char *iupac) { _iupac.assign(iupac); }

        //! Set the IUPAC name
        void SetIupac(std::string iupac) { _iupac = iupac; }

        //! Return IUPAC name or, if not found, the molname
        const std::string &getIupac() const
        {
            if (_iupac.size() > 0)
            {
                return _iupac;
            }
            else
            {
                return _molname;
            }
        }

        //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
        void SetCas(const char *cas) { _cas.assign(cas); }

        //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
        void SetCas(std::string cas) { _cas.assign(cas); }

        //! Return the CAS (Chemical Abstract Service) identifier, see http:://www.cas.org
        const std::string &getCas() const { return _cas; }

        //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
        void SetCid(const char *cid) { _cid.assign(cid); }

        //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
        void SetCid(std::string cid) { _cid = cid; }

        //! Return the CID (Chemspider identifier) see http:://www.chemspider.com
        const std::string &getCid() const { return _cid; }

        //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        void SetInchi(const char *inchi) { _inchi.assign(inchi); }

        //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        void SetInchi(std::string inchi) { _inchi = inchi; }

        //! Return the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        const std::string &getInchi() const { return _inchi; }

        //! Convenience function
        bool getPropRef(MolPropObservable mpo, iqmType iQM, char *lot,
                        const char *conf, const char *type,
                        double *value, double *error, double *T,
                        std::string &ref, std::string &mylot,
                        double vec[3], tensor quadrupole);

        //! And another one
        bool getProp(MolPropObservable mpo, iqmType iQM, char *lot,
                     char *conf, char *type,
                     double *value, double *error, double *T);

        //! Returns true if the HF energy of the optimized geometry exists and returns the HF
        bool getOptHF(double *value);

        //! Add a classification category for this molecule
        void AddCategory(const std::string &category)
        {
            if (!SearchCategory(category)) { category_.push_back(category); }
        }

        //! Return the number of categories
        int NCategory() { return category_.size(); }

        //! Return iterator Begin over categories
        std::vector<std::string>::iterator BeginCategory() { return category_.begin(); }

        //! Return iterator End over categories
        std::vector<std::string>::iterator EndCategory() { return category_.end(); }

        //! Return true if catname is an existing category
        bool SearchCategory(const std::string &catname) const;

        //! Delete a composition type if present
        void DeleteComposition(const std::string &compname);

        //! Add a composition entry
        void AddComposition(MolecularComposition mc);

        std::vector<MolecularComposition> &MolComp()
        { return _mol_comp; }

        const std::vector<MolecularComposition> &MolComp() const
        { return _mol_comp; }

        //! Begin Iterator over MolecularCompostion items
        MolecularCompositionIterator BeginMolecularComposition() { return _mol_comp.begin(); }

        //! End Iterator over MolecularCompostion items
        MolecularCompositionIterator EndMolecularComposition()   { return _mol_comp.end(); }

        //! Last Iterator over MolecularCompostion items
        MolecularComposition *LastMolecularComposition()   { return &(_mol_comp.back()); }

        //! Search for particular MolecularCompostion item or return EndMolecularComposition if not found
        MolecularCompositionIterator SearchMolecularComposition(std::string str);

        //! Return number of atoms in the first composition if present, or 0 otherwise
        int NAtom();

        //! Routine to generate compositions based on calculation data
        bool GenerateComposition(const Poldata &pd);

        //! Returns boolean stating whether a particular composition is present
        bool HasComposition(char *composition) { std::string _str(composition); return HasComposition(_str); }

        //! Returns boolean stating whether a particular composition is present
        bool HasComposition(std::string composition);

        //! Add a Bond element
        void AddBond(Bond b);

        //! Check whether a Bond element is present already
        bool BondExists(Bond b);

        //! Return the number of Bond elements
        int NBond() const { return _bond.size(); }

        //! Begin Iterator over Bond elements
        BondIterator BeginBond() { return _bond.begin(); }

        //! End Iterator over Bond elements
        BondIterator EndBond() { return _bond.end(); }

        //! Add an experiment
        void AddExperiment(Experiment myexp) { _exper.push_back(myexp); }

        void Stats()
        {
            printf("%s - %s - %d experiments\n",
                   _molname.c_str(), _formula.c_str(),
                   (int)_exper.size());
        }
        //! Return the number of experiments
        int NExperiment() const { return _exper.size(); }

        //! Iterator Begin over experiments
        ExperimentIterator BeginExperiment() { return _exper.begin(); }

        //! Iterator End over experiments
        ExperimentIterator EndExperiment() { return _exper.end(); }

        //! Return pointer to the last inserted experiment or NULL if the number of experiments is zero
        Experiment *LastExperiment()
        {
            if (NExperiment() > 0)
            {
                return &(_exper.back());
            }
            else
            {
                return nullptr;
            }
        }

        //! Return a calculation iterator corresponding to the level of theory (lot) parameter, or EndExperiment in case it is not found
        ExperimentIterator getLot(const char *lot);

        //! Return a calculation iterator corresponding to the level of theory (lot) parameter, or EndCalculation in case it is not found
        //! The operator should hold the requested observable of the type (can be NULL)
        ExperimentIterator getLotPropType(const char       *lot,
                                          MolPropObservable mpo,
                                          const char       *type);
/*! \brief
 * Sends this object over an MPI connection
 *
 * \param[in] commrec   GROMACS data structure for MPI communication
 * \param[in] dest      Destination processor
 * \return the CommunicationStatus of the operation
 */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolProp items
typedef std::vector<MolProp>::iterator MolPropIterator;
typedef std::vector<MolProp>::const_iterator MolPropConstIterator;

/*! \brief Utility to compare temperatures
 *
 * Compares two temperatures
 * \param[in] Tref The reference temperature
 * \param[in] T    The temperature to be tested
 * \return true if the reference T < 0, or the difference between the two
 *              T is negligable.
 */
bool bCheckTemperature(double Tref, double T);

}

#endif
