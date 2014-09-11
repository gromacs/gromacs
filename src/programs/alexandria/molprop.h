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
/*! \defgroup module_alexandria Processing of input files and force fields
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
 * \inpublicapi
 * \ingroup module_alexandria
 */
#ifndef MOLPROP_H
#define MOLPROP_H

#include <string>
#include <vector>
#include <string.h>
#include "gromacs/utility/real.h"
#include "gromacs/topology/atomprop.h"
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

        //! Destructor
        ~AtomNum() {};

        //! Return the name of the atom for this AtomNum
        std::string GetAtom() { return _catom; }

        //! Set the name of the atom for this AtomNum
        void SetAtom(std::string catom) { _catom = catom; }

        //! Set the name of the atom for this AtomNum
        void SetAtom(const char *catom) { _catom.assign(catom); }

        //! Return the number of atoms for this AtomNum
        int GetNumber() { return _cnumber; }

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

        //! Destructor
        ~MolecularComposition() {}

        //! Return the composition name
        std::string GetCompName() { return _compname; }

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
        std::string _type, _unit;
    public:
        //! Empty constructor
        GenericProperty() {}

        /*! \brief
         * Creates a new GenericProperty object.
         *
         * \param[in] type  Type of the property
         * \param[in] unit  Unit of the property
         */
        GenericProperty(char *type, char *unit) { SetType(type); SetUnit(unit); };

        /*! \brief
         * Creates a new GenericProperty object.
         *
         * \param[in] type  Type of the property
         * \param[in] unit  Unit of the property
         */
        GenericProperty(std::string type, std::string unit) { SetType(type); SetUnit(unit); };

        //! Destructor
        ~GenericProperty() {};

        //! Return the property type
        std::string GetType() { return _type; }

        //! Return the unit of the property
        std::string GetUnit() { return _unit; }

        //! Set the type of the property
        void SetType(std::string type);

        //! Set the unit of the property
        void SetUnit(std::string unit);

        //! Set the type of the property
        void SetType(const char *type) { std::string s(type); SetType(s); }

        //! Set the unit of the property
        void SetUnit(const char *unit) { std::string s(unit); SetUnit(s); }

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
        double _xx, _yy, _zz, _xy, _xz, _yz;
    public:
        //! Empty constructor
        MolecularQuadrupole() {}

        //! Constructor initiating all elements of the quadrupole tensor
        MolecularQuadrupole(char *type, char *unit, double xx, double yy, double zz, double xy, double xz, double yz) : GenericProperty(type, unit) { Set(xx, yy, zz, xy, xz, yz); };

        //! Constructor initiating all elements of the quadrupole tensor
        MolecularQuadrupole(std::string type, std::string unit, double xx, double yy, double zz, double xy, double xz, double yz) : GenericProperty(type, unit) { Set(xx, yy, zz, xy, xz, yz); };

        //! Destructor
        ~MolecularQuadrupole() {};

        //! Set all the elements of the qudrupole tensor
        void Set(double xx, double yy, double zz, double xy, double xz, double yz) { _xx = xx; _yy = yy; _zz = zz; _xy = xy; _xz = xz; _yz = yz; };

        //! Get all the elements of the qudrupole tensor
        void Get(double *xx, double *yy, double *zz, double *xy, double *xz, double *yz) { *xx = _xx; *yy = _yy; *zz = _zz; *xy = _xy; *xz = _xz; *yz = _yz; };

        //! Return the XX component of the quadrupole tensor
        double GetXX() { return _xx; }

        //! Return the YY component of the quadrupole tensor
        double GetYY() { return _yy; }

        //! Return the ZZ component of the quadrupole tensor
        double GetZZ() { return _zz; }

        //! Return the XY component of the quadrupole tensor
        double GetXY() { return _xy; }

        //! Return the XZ component of the quadrupole tensor
        double GetXZ() { return _xz; }

        //! Return the YZ component of the quadrupole tensor
        double GetYZ() { return _yz; }

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
        MolecularEnergy(char *type, char *unit, double value, double error) : GenericProperty(type, unit) { Set(value, error); };

        //! Constructor storing all properties related to this energy term
        MolecularEnergy(std::string type, std::string unit, double value, double error) : GenericProperty(type, unit) { Set(value, error); };

        //! Destructor
        ~MolecularEnergy() {};

        //! Set the value and error for the energy
        void Set(double value, double error) { _value = value; _error = error; };

        //! Get the value and error for this energy
        void Get(double *value, double *error) { *value = _value; *error = _error; };

        //! Set the value for the energy
        void SetValue(double value) { _value = value; };

        //! Return the energy value
        double GetValue() { return _value; };

        //! Set the error in the energy value
        void SetError(double error) { _value = error; };

        //! Return the error in the energy
        double GetError() { return _error; };

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
 * Contains either a dipole vector, or the diagonal of the polarizability tensor.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularDipPolar : public GenericProperty
{
    private:
        double _x, _y, _z;
        double _aver, _error;
    public:
        //! Empty constructor
        MolecularDipPolar() {}

        //! Constructor storing all properties related to this dipole/polarizability
        MolecularDipPolar(char *type, char *unit, double x, double y, double z, double aver, double error) : GenericProperty(type, unit) { Set(x, y, z, aver, error); }

        //! Constructor storing all properties related to this dipole/polarizability
        MolecularDipPolar(std::string type, std::string unit, double x, double y, double z, double aver, double error) : GenericProperty(type, unit) { Set(x, y, z, aver, error); }

        //! Destructor
        ~MolecularDipPolar() {};

        //! Set all properties related to this dipole/polarizability
        void Set(double x, double y, double z, double aver, double error) { _x = x; _y = y; _z = z; _aver = aver; _error = error; };

        //! Return all properties of this dipole/polarizability
        void Get(double *x, double *y, double *z, double *aver, double *error) { *x = _x; *y = _y; *z = _z; *aver = _aver; *error = _error; };

        //! Return the average dipole/polarizability value
        double GetAver() { return _aver; }

        //! Return the error in the average dipole/polarizability
        double GetError() { return _error; }

        //! Return the X component of the dipole/polarizability
        double GetX() { return _x; }

        //! Return the Y component of the dipole/polarizability
        double GetY() { return _y; }

        //! Return the Z component of the dipole/polarizability
        double GetZ() { return _z; }

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
//! Iterates over a vector of MolecularDipPolar
typedef std::vector<MolecularDipPolar>::iterator MolecularDipPolarIterator;

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
        std::string _xyz_unit, _V_unit;
        int         _espid;
        double      _x, _y, _z, _V;
    public:
        //! Empty constructor
        ElectrostaticPotential() {}

        //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        ElectrostaticPotential(std::string xyz_unit, std::string V_unit, int espid, double x, double y, double z, double V) { Set(xyz_unit, V_unit, espid, x, y, z, V); };

        //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        ElectrostaticPotential(const char *xyz_unit, const char *V_unit, int espid, double x, double y, double z, double V) { Set(xyz_unit, V_unit, espid, x, y, z, V); };

        //! Destructor
        ~ElectrostaticPotential() {};

        //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Set(std::string xyz_unit, std::string V_unit, int espid, double x, double y, double z, double V) { _xyz_unit = xyz_unit; _V_unit = V_unit; _espid = espid; _x = x; _y = y; _z = z; _V = V; };

        //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Set(const char *xyz_unit, const char *V_unit, int espid, double x, double y, double z, double V)
        {
            std::string _x(xyz_unit), _V(V_unit); Set(_x, _V, espid, x, y, z, V);
        }

        //! Return the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Get(char **xyz_unit, char **V_unit, int *espid, double *x, double *y, double *z, double *V)
        {
            *xyz_unit = strdup(_xyz_unit.c_str()); *V_unit = strdup(_V_unit.c_str()); *espid = _espid; *x = _x; *y = _y; *z = _z; *V = _V;
        };

        //! Return the unit of the coordinates
        std::string GetXYZunit() { return _xyz_unit; }

        //! Return the unit of the potential
        std::string GetVunit() { return _V_unit; }

        //! Return the ESP id from the original calculation
        int GetEspid() { return _espid; }

        //! Return the X coordinate of the ESP point
        double GetX() { return _x; }

        //! Return the Y coordinate of the ESP point
        double GetY() { return _y; }

        //! Return the Z coordinate of the ESP point
        double GetZ() { return _z; }

        //! Return the electrostatic potential at this point in space
        double GetV() { return _V; }

        //! Set the potential for this instance of the class
        void SetV(double V) { _V = V; }

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

        //! Destructor
        ~Bond() {};

        //! Sets the ids of the atoms and the bondorder
        void Set(int ai, int aj, int bondorder) {_ai = ai; _aj = aj; _bondorder = bondorder; };

        //! Returns the ids of the atoms and the bondorder
        void Get(int *ai, int *aj, int *bondorder) { *ai = _ai; *aj = _aj; *bondorder = _bondorder; };

        //! Returns the first atom id
        int GetAi() { return _ai; }

        //! Returns the second atom id
        int GetAj() { return _aj; }

        //! Returns the bondorder
        int GetBondOrder() { return _bondorder; }

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
        AtomicCharge(const char *type, const char *unit, double q) : GenericProperty(type, unit) { SetQ(q); };

        //! Constructor setting type, unit and charge itself
        AtomicCharge(std::string type, std::string unit, double q) : GenericProperty(type, unit) { SetQ(q); };

        //! Destructor
        ~AtomicCharge() {};

        //! Set the charge to q
        void SetQ(double q) { _q = q; };

        //! Return the charge
        double GetQ() { return _q; }

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
        std::string               _name, _obtype, _unit;
        double                    _x, _y, _z;
        int                       _atomid;
        std::vector<AtomicCharge> _q;
    public:
        //! Empty constructor
        CalcAtom() {}

        //! Constructor initiating the name, type and atomid
        CalcAtom(const char *name, const char *obtype, int atomid)
        {
            _name.assign(name); _obtype.assign(obtype); _atomid = atomid;
        };

        //! Constructor initiating the name, type and atomid
        CalcAtom(std::string name, std::string obtype, int atomid)
        {
            _name = name; _obtype = obtype; _atomid = atomid;
        };

        //! Destructur
        ~CalcAtom() {};


        //! Function returning true if the two atoms are equal
        bool Equal(CalcAtom ca);

        //! Add an AtomicCharge element to the atom
        void AddCharge(AtomicCharge aq);

        //! Begin Iterator over AtomicCharge items
        AtomicChargeIterator BeginQ() { return _q.begin(); }

        //! End Iterator over AtomicCharge items
        AtomicChargeIterator EndQ() { return _q.end(); }

        //! Return the atom id of the atom
        int GetAtomid() { return _atomid; }

        //! Return the name of the atom
        std::string GetName() { return _name; }

        //! Return the OpenBabel type of the atom
        std::string GetObtype() { return _obtype; }

        //! Return the unit of the coordinates of the atom
        std::string GetUnit() { return _unit; }

        //! Set the unit of the coordinates of the atom
        void SetUnit(std::string unit);

        //! Set the unit of the coordinates of the atom
        void SetUnit(const char *unit) { std::string s(unit); SetUnit(s); }

        //! Set the coordinates of the atom
        void SetCoords(double x, double y, double z) { _x = x; _y = y; _z = z; }

        //! Return all the coordinates of the atom
        void GetCoords(double *x, double *y, double *z) { *x = _x; *y = _y; *z = _z; }

        //! Return the X coordinate of the atom
        double GetX() { return _x; }

        //! Return the Y coordinate of the atom
        double GetY() { return _y; }

        //! Return the Z coordinate of the atom
        double GetZ() { return _z; }

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
        std::string                      _reference, _conformation;
        std::vector<MolecularDipPolar>   _polar, _dipole;
        std::vector<MolecularEnergy>     _energy;
        std::vector<MolecularQuadrupole> _quadrupole;

    public:
        //! Empty constructor
        Experiment() { }

        //! Constructor initiating an Experiment with reference and conformation
        Experiment(std::string reference, std::string conformation)
        {
            _reference = reference; _conformation = conformation;
        };

        //! Constructor initiating an Experiment with reference and conformation
        Experiment(const char *reference, const char *conformation)
        {
            _reference.assign(reference); _conformation.assign(conformation);
        };

        //! Default destructor
        ~Experiment() {};

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

        //! Add a MolecularDipPolar element containing polarizability data
        void AddPolar(MolecularDipPolar mdp) { _polar.push_back(mdp); }

        //! Return Begin Iterator over polarizability elements
        MolecularDipPolarIterator BeginPolar() { return _polar.begin(); }

        //! Return End Iterator over polarizability elements
        MolecularDipPolarIterator EndPolar()   { return _polar.end(); }

        //! Number of polarizability values
        int NPolar() { return _polar.size(); }

        //! Number of dipole values
        int NDipole() { return _dipole.size(); }

        //! Add a MolecularDipPolar element containing dipole data
        void AddDipole(MolecularDipPolar mdp) { _dipole.push_back(mdp); }

        //! Return Begin Iterator over dipole elements
        MolecularDipPolarIterator BeginDipole() { return _dipole.begin(); }

        //! Return End Iterator over dipole elements
        MolecularDipPolarIterator EndDipole()   { return _dipole.end(); }

        //! Number of quadrupole values
        int NQuadrupole() { return _quadrupole.size(); }

        //! Add a MolecularQuadrupole element
        void AddQuadrupole(MolecularQuadrupole mq) { _quadrupole.push_back(mq); }

        //! Return Begin Iterator over quadrupole elements
        MolecularQuadrupoleIterator BeginQuadrupole() { return _quadrupole.begin(); }

        //! Return End Iterator over quadrupole elements
        MolecularQuadrupoleIterator EndQuadrupole() { return _quadrupole.end(); }

        //! Add a MolecularEnergy element
        void AddEnergy(MolecularEnergy me) { _energy.push_back(me); }

        //! Return Begin Iterator over energy elements
        MolecularEnergyIterator BeginEnergy() { return _energy.begin(); }

        //! Return End Iterator over energy elements
        MolecularEnergyIterator EndEnergy()   { return _energy.end(); }

        //! Return the molecular conformation
        std::string GetConformation() { return _conformation; }

        //! Return the literature reference
        std::string GetReference() { return _reference; }

        /*! \brief
         * Convenience function that fetches a value from this experiment
         *
         * \param[in]  type   The type of the data (dependent on whether it is dipole, energy etc.)
         * \param[in]  mpo    Enum selecting the type of data to fetch
         * \param[out] value  The value of e.g. the energy
         * \param[out] error  The error in the value of e.g. the energy
         * \param[out] vec    Vector data to be output, dipole or diagonal element of polarizability
         * \param[out] quadrupole The quadrupole tensor
         * \return true on success
         */
        bool GetVal(const char *type, MolPropObservable mpo,
                    double *value, double *error, double vec[3],
                    tensor quadrupole);

        //! Merge in another object - Low level function
        void MergeLow(Experiment *src);

        //! Merge in another object
        void Merge(Experiment &src);

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
 * Contains data on a molecule based on a calculation.
 *
 * This is a composite class holding the results from a calculation, either
 * the dipole, polarizability, energy, the quadrupole, or the electrostatic
 * potential. The type of method used and basisset where appropriate are stored,
 * along with the program used and the datafile where the data came from.
 * The conformation of the molecule (if known) is stored.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Calculation : public Experiment
{
    private:
        std::string                         _program, _method, _basisset, _datafile;
        std::vector<CalcAtom>               _catom;
        std::vector<ElectrostaticPotential> _potential;
    public:
        //! Empty constructor
        Calculation() { }

        //! Constructor for calculations with program, method, basisset, reference, conformation and datafile
        Calculation(std::string program, std::string method,
                    std::string basisset, std::string reference,
                    std::string conformation, std::string datafile) : Experiment(reference, conformation)
        {
            _program  = program; _method = method;
            _basisset = basisset; _datafile = datafile;
        };

        //! Constructor for calculations with program, method, basisset, reference, conformation and datafile
        Calculation(const char *program, const char *method,
                    const char *basisset, const char *reference,
                    const char *conformation, const char *datafile) : Experiment(reference, conformation)
        {
            _program.assign(program); _method.assign(method);
            _basisset.assign(basisset); _datafile.assign(datafile);
        };

        //! Destructor
        ~Calculation() {};

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

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
        std::string GetProgram() { return _program; }

        //! Return the basis set used to perform the calculation
        std::string GetBasisset() { return _basisset; }

        //! Return the method used to perform the calculation
        std::string GetMethod() { return _method; }

        //! Return the datafile from which the calculation output was extracted
        std::string GetDatafile() { return _datafile; }

        //! Merge in another object
        void Merge(Calculation &src);

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
//! Iterates over Calculation items
typedef std::vector<Calculation>::iterator CalculationIterator;
//typedef ExperimentIterator CalculationIterator;

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
        std::vector<std::string>          _category;
        std::vector<MolecularComposition> _mol_comp;
        std::vector<Calculation>          _calc;
        std::vector<Experiment>           _exper;
        std::vector<Bond>                 _bond;
    public:
        //! Construct a number MolProp object
        MolProp() { _index = NOTSET; _mass = 0; _charge = 0; _multiplicity = 0; }

        //! Destructor
        ~MolProp() {}

        /*! \brief
         * Check the internal consistency of this object
         *
         * \todo Implement this
         */
        void CheckConsistency();

        //! Set the index number for sorting
        void SetIndex(int index) { _index = index; }

        //! Return the index number for sorting
        int GetIndex() { return _index; }

        //! Set the molecular mass
        void SetMass(double mass) { _mass = mass; }

        //! Return the molecular mass
        double GetMass() { return _mass; }

        //! Set the total charge of the molecule
        void SetCharge(double charge) { _charge = charge; }

        //! Return the total charge of the molecule
        int GetCharge() { return _charge; }

        //! Set the multiplicity of the molecule
        void SetMultiplicity(int multiplicity) { _multiplicity = multiplicity; }

        //! Return the multiplicity  of the molecule
        int GetMultiplicity() { return _multiplicity; }

        /*! \brief
         * Merge the content of another MolProp into this one
         *
         * \param[in] mpi The object to be merged into the present one
         * \todo Check and double check
         */
        void Merge(MolProp &mpi);

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

        //! Set the LaTeX formula
        void SetTexFormula(const char *formula) { _texform.assign(formula); }

        //! Set the formula
        void SetFormula(const char *formula) { _formula.assign(formula); }

        //! Set the formula
        void SetFormula(std::string formula) { _formula = formula; }

        //! Return the formula
        std::string GetFormula() { return _formula; }

        //! Return the LaTeX formula
        std::string GetTexFormula();

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
        std::string GetMolname() { return _molname; }

        //! Set the IUPAC name
        void SetIupac(const char *iupac) { _iupac.assign(iupac); }

        //! Set the IUPAC name
        void SetIupac(std::string iupac) { _iupac = iupac; }

        //! Return IUPAC name or, if not found, the molname
        std::string GetIupac() { if (_iupac.size() > 0) { return _iupac; } else{ return _molname; } }

        //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
        void SetCas(const char *cas) { _cas.assign(cas); }

        //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
        void SetCas(std::string cas) { _cas.assign(cas); }

        //! Return the CAS (Chemical Abstract Service) identifier, see http:://www.cas.org
        std::string GetCas() { return _cas; }

        //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
        void SetCid(const char *cid) { _cid.assign(cid); }

        //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
        void SetCid(std::string cid) { _cid = cid; }

        //! Return the CID (Chemspider identifier) see http:://www.chemspider.com
        std::string GetCid() { return _cid; }

        //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        void SetInchi(const char *inchi) { _inchi.assign(inchi); }

        //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        void SetInchi(std::string inchi) { _inchi = inchi; }

        //! Return the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        std::string GetInchi() { return _inchi; }

        //! Convenience function
        bool GetPropRef(MolPropObservable mpo, iqmType iQM, char *lot,
                        const char *conf, const char *type, double *value, double *error,
                        char **ref, char **mylot,
                        double vec[3], tensor quadrupole);

        //! And another one
        bool GetProp(MolPropObservable mpo, iqmType iQM, char *lot,
                     char *conf, char *type, double *value);

        //! Add a classification category for this molecule
        void AddCategory(const char *category)
        {
            std::string _str(category); AddCategory(_str);
        }

        //! Add a classification category for this molecule
        void AddCategory(std::string category)
        {
            if (SearchCategory(category) == 0) { _category.push_back(category); }
        }

        //! Return the number of categories
        int NCategory() { return _category.size(); }

        //! Return iterator Begin over categories
        std::vector<std::string>::iterator BeginCategory() { return _category.begin(); }

        //! Return iterator End over categories
        std::vector<std::string>::iterator EndCategory() { return _category.end(); }

        //! Return true if catname is an existing category
        bool SearchCategory(const char *catname)
        {
            std::string _str(catname); return SearchCategory(_str);
        }

        //! Return true if catname is an existing category
        bool SearchCategory(std::string catname);

        //! Delete a composition type if present
        void DeleteComposition(std::string compname);

        //! Delete a composition type if present
        void DeleteComposition(const char *compname)
        {
            std::string _str(compname); DeleteComposition(_str);
        }

        //! Add a composition entry
        void AddComposition(MolecularComposition mc);

        //! Begin Iterator over MolecularCompostion items
        MolecularCompositionIterator BeginMolecularComposition() { return _mol_comp.begin(); }

        //! End Iterator over MolecularCompostion items
        MolecularCompositionIterator EndMolecularComposition()   { return _mol_comp.end(); }

        //! Last Iterator over MolecularCompostion items
        MolecularComposition *LastMolecularComposition()   { return &(_mol_comp.back()); }

        //! Search for particular MolecularCompostion item or return EndMolecularComposition if not found
        MolecularCompositionIterator SearchMolecularComposition(const char *str)
        {
            std::string _str(str); return SearchMolecularComposition(_str);
        }

        //! Search for particular MolecularCompostion item or return EndMolecularComposition if not found
        MolecularCompositionIterator SearchMolecularComposition(std::string str);

        //! Return number of atoms in the first composition if present, or 0 otherwise
        int NAtom();

        //! Routine to generate compositions based on calculation data
        bool GenerateComposition(gmx_poldata_t pd);

        //! Returns boolean stating whether a particular composition is present
        bool HasComposition(char *composition) { std::string _str(composition); return HasComposition(_str); }

        //! Returns boolean stating whether a particular composition is present
        bool HasComposition(std::string composition);

        //! Add a Bond element
        void AddBond(Bond b);

        //! Check whether a Bond element is present already
        bool BondExists(Bond b);

        //! Return the number of Bond elements
        int NBond() { return _bond.size(); }

        //! Begin Iterator over Bond elements
        BondIterator BeginBond() { return _bond.begin(); }

        //! End Iterator over Bond elements
        BondIterator EndBond() { return _bond.end(); }

        //! Add an experiment
        void AddExperiment(Experiment myexp) { _exper.push_back(myexp); }

        void Stats()
        {
            printf("%s - %s - %d experiments - %d calculations\n",
                   _molname.c_str(), _formula.c_str(),
                   (int)_exper.size(), (int)_calc.size());
        }
        //! Return the number of experiments
        int NExperiment() { return _exper.size(); }

        //! Iterator Begin over experiments
        ExperimentIterator BeginExperiment() { return _exper.begin(); }

        //! Iterator End over experiments
        ExperimentIterator EndExperiment() { return _exper.end(); }

        //! Return pointer to the last inserted experiment or NULL if the number of experiments is zero
        Experiment *LastExperiment()
        {
            if (NExperiment() > 0) { return &(_exper.back()); } else{ return NULL; }
        }

        //! Add a calculation
        void AddCalculation(Calculation calc) { _calc.push_back(calc); }

        //! Return the number of calculations
        int NCalculation() { return _calc.size(); }

        //! Iterator Begin over calculations
        CalculationIterator BeginCalculation() { return _calc.begin(); }

        //! Iterator End over calculations
        CalculationIterator EndCalculation() { return _calc.end(); }

        //! Return pointer to the last inserted calculation or NULL if the number of calculations is zero
        Calculation *LastCalculation()
        {
            if (NCalculation() > 0) { return &(_calc.back()); } else{ return NULL; }
        }

        //! Return a calculation iterator corresponding to the level of theory (lot) parameter, or EndCalculation in case it is not found
        CalculationIterator GetLot(const char *lot);

        //! Return a calculation iterator corresponding to the level of theory (lot) parameter, or EndCalculation in case it is not found
        //! The operator should hold the requested observable of the type (can be NULL)
        CalculationIterator GetLotPropType(const char       *lot,
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
}

#endif
