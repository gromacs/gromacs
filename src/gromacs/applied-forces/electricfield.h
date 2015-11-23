/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_APPLIED_FORCES_ELECTRICFIELD_H
#define GMX_APPLIED_FORCES_ELECTRICFIELD_H

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/real.h"

struct t_commrec;
struct t_inpfile;
struct warninp;

/*! \libinternal
 * \brief Declaration of storage unit for fields
 *
 */
class ElectricFieldData
{
    public:
        ElectricFieldData() { a_ = 0; omega_ = 0; t0_ = 0; sigma_ = 0; };

        /*! \brief Initiate the field values
         *
         * \param[in] a     Amplitude
         * \param[in] omega Frequency
         * \param[in] t0    Peak of the pulse
         * \param[in] sigma Width of the pulse
         */
        void setField(real a, real omega, real t0, real sigma)
        {
            a_     = a;
            omega_ = omega;
            t0_    = t0;
            sigma_ = sigma;
        }

        //! Return the amplitude
        real a()     const { return a_; }
        //! Return the frequency
        real omega() const { return omega_; }
        //! Return the time for the peak of the pulse
        real t0()    const { return t0_; }
        //! Return the width of the pulse (0 means inifinite)
        real sigma() const { return sigma_; }
    private:
        //! Coeffient (V / nm)
        real a_;
        //! Frequency (1/ps)
        real omega_;
        //! Central time point (ps) for pulse
        real t0_;
        //! Width of pulse (ps, if zero there is no pulse)
        real sigma_;
};

/*! \libinternal
 * \brief Describe time dependent electric field
 *
 * This class describes the time dependent electric field that can
 * be applied to all charges in a simulation. The field is described
 * by the following:
 *     E(t) = A cos(omega*(t-t0))*exp(-sqr(t-t0)/(2.0*sqr(sigma)));
 * If sigma = 0 there is no pulse and we have instead
 *     E(t) = A cos(omega*t)
 *
 * force is kJ mol^-1 nm^-1 = e * kJ mol^-1 nm^-1 / e
 *
 * WARNING:
 * There can be problems with the virial.
 * Since the field is not self-consistent this is unavoidable.
 * For neutral molecules the virial is correct within this approximation.
 * For neutral systems with many charged molecules the error is small.
 * But for systems with a net charge or a few charged molecules
 * the error can be significant when the field is high.
 * Solution: implement a self-consistent electric field into PME.
 */
class ElectricField : public gmx::IInputRecExtension, public IForceProvider
{
    public:
        ElectricField() : isSet_(false), fp_(nullptr) {}

        /*! \brief Read or write the data from a tpx file
         *
         * \param[in] fio   The file input/output handle
         * \param[in] bRead true if we should read, false otherwise
         */
        virtual void doTpxIO(t_fileio *fio, bool bRead);
        virtual void readMdp(int *ninp_p, t_inpfile **inp_p, warninp *wi);

        /*! \brief Extract relevant fields from an mdp file
         *
         * \param[in] dim          The direction XX, YY, or ZZ
         * \param[in] staticField  Static components of the field
         * \param[in] dynamicField Dynamic components of the field
         * \param[in] wi           Warning control structure
         */
        void decodeMdp(int         dim,
                       const char *staticField,
                       const char *dynamicField,
                       warninp    *wi);

        /*! \brief Broadcast the contents over processors
         *
         * \param[in] cr The communication record
         */
        virtual void broadCast(const t_commrec *cr);

        /*! \brief Initiate output parameters
         *
         * \param[in] fplog File pointer for log messages
         * \param[in] nfile Number of files
         * \param[in] fnm   Array of filenames and properties
         * \param[in] bAppendFiles Whether or not we should append to files
         * \param[in] oenv  The output environment for xvg files
         */
        virtual void initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                                bool bAppendFiles, const gmx_output_env_t *oenv);
        //! Finalize output
        virtual void finishOutput();

        /*! \brief Set/initiate relevant options in the forcerec structure
         *
         * \param[inout] fr The forcerec structure
         */
        virtual void initForcerec(t_forcerec *fr);

        /*! \brief Add a component to the electric field
         *
         * The electric field has three spatial dimensions that are
         * added to the data structure one at a time.
         * \param[in] dim   Dimension, XX, YY, ZZ (0, 1, 2)
         * \param[in] a     Amplitude of the field in V/nm
         * \param[in] omega Frequency (1/ps)
         * \param[in] t0    Time of pulse peak (ps)
         * \param[in] sigma Width of peak (ps)
         */
        void setFieldTerm(int dim, real a, real omega, real t0, real sigma);

        /*! \brief Return the field strength
         *
         * \param[in] dim The spatial direction
         * \param[in] t   The time (ps)
         * \return The field strength in V/nm units
         */
        real field(int dim, real t) const;

        /*! \brief Return amplitude of field
         *
         * \param[in] dim Direction of the field (XX, YY, ZZ)
         * \return Amplitude of the field
         */
        real a(int dim)     const { return efield_[dim].a(); }
        /*! \brief Return frequency of field (1/ps)
         *
         * \param[in] dim Direction of the field (XX, YY, ZZ)
         * \return Frequency of the field
         */
        real omega(int dim) const { return efield_[dim].omega(); }
        /*! \brief Return time of pulse peak
         *
         * \param[in] dim Direction of the field (XX, YY, ZZ)
         * \return Time of pulse peak
         */
        real t0(int dim) const { return efield_[dim].t0(); }
        /*! \brief Return width of the pulse
         *
         * \param[in] dim Direction of the field (XX, YY, ZZ)
         * \return Width of the pulse
         */
        real sigma(int dim) const { return efield_[dim].sigma(); }

        //! Return whether or not to apply a field
        bool applyField() const;

        /*! \brief Dump the parameters to a text file
         *
         * \param[in] fp     File pointer to write to
         * \param[in] indent Level of indentation to start printing
         */
        virtual void printParameters(FILE *fp, int indent);

        /*! \brief Print the field components to a file
         *
         * \param[in] fp  File pointer, must be open
         * \param[in] t   The time
         * Will throw and exit with fatal error if file is not open.
         */
        void printComponents(FILE *fp, double t) const;

        /*! \brief Compare the parameters of two fields
         *
         * \param[in] fp     File pointer to write the result of the comparison
         * \param[in] field2 The other electric field structure
         * \param[in] reltol The relative tolerance when comparing reals
         * \param[in] abstol The absolute tolerance when comparing reals
         */
        virtual void compare(FILE                          *fp,
                             const gmx::IInputRecExtension *field2,
                             real                           reltol,
                             real                           abstol);

        // From IForceProvider
        virtual void calculateForces(const t_commrec *cr,
                                     const t_mdatoms *atoms,
                                     rvec f[], double t);
    private:
        //! Has the field been initiated
        bool              isSet_;
        //! The field strength in each dimension
        ElectricFieldData efield_[DIM];
        //! Log file pointer
        FILE             *fp_;
};
#endif
