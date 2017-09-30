/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares data structure and utilities for electric fields
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "electricfield.h"

#include <cmath>

#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Declaration of storage unit for fields
 */
class ElectricFieldData
{
    public:
        ElectricFieldData() : a_(0), omega_(0), t0_(0), sigma_(0)
        {
        }

        /*! \brief
         * Adds an option section to specify parameters for this field component.
         */
        void initMdpOptions(IOptionsContainerWithSections *options, const char *sectionName)
        {
            auto section = options->addSection(OptionSection(sectionName));
            section.addOption(RealOption("E0").store(&a_));
            section.addOption(RealOption("omega").store(&omega_));
            section.addOption(RealOption("t0").store(&t0_));
            section.addOption(RealOption("sigma").store(&sigma_));
        }
        /*! \brief
         * Creates mdp parameters for this field component.
         */
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder, const std::string &name) const
        {
            builder->addUniformArray<real>("E-" + name, {1, a_, -1});
            builder->addUniformArray<real>("E-" + name + "t", {omega_, t0_, sigma_});
        }

        /*! \brief Evaluates this field component at given time.
         *
         * \param[in] t The time to evualate at
         * \return The electric field
         */
        real evaluate(real t) const
        {
            if (sigma_ > 0)
            {
                return a_ * (std::cos(omega_*(t-t0_))
                             * std::exp(-square(t-t0_)/(2.0*square(sigma_))));
            }
            else
            {
                return a_ * std::cos(omega_*t);
            }
        }

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

/*! \internal
 * \brief Describe time dependent electric field
 *
 * Class that implements a force to be evaluated in mdrun.
 * The electric field can be pulsed and oscillating, simply
 * oscillating, or static, in each of X,Y,Z directions.
 */
class ElectricField final : public IMDModule,
                            public IMdpOptionProvider, public IMDOutputProvider,
                            public IForceProvider
{
    public:
        ElectricField() : fpField_(nullptr) {}

        // From IMDModule
        IMdpOptionProvider *mdpOptionProvider() override { return this; }
        IMDOutputProvider *outputProvider() override { return this; }
        void initForceProviders(ForceProviders *forceProviders) override
        {
            if (isActive())
            {
                forceProviders->addForceProvider(this);
            }
        }

        // From IMdpOptionProvider
        void initMdpTransform(IKeyValueTreeTransformRules *transform) override;
        void initMdpOptions(IOptionsContainerWithSections *options) override;
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const override;

        // From IMDOutputProvider
        void initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                        bool bAppendFiles, const gmx_output_env_t *oenv) override;
        void finishOutput() override;

        // From IForceProvider
        //! \copydoc IForceProvider::calculateForces()
        void calculateForces(const t_commrec       *cr,
                             const t_mdatoms       *mdatoms,
                             const matrix           box,
                             double                 t,
                             const rvec            *x,
                             gmx::ForceWithVirial  *forceWithVirial) override;

    private:
        //! Return whether or not to apply a field
        bool isActive() const;

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

        /*! \brief Print the field components to a file
         *
         * \param[in] t   The time
         * Will throw and exit with fatal error if file is not open.
         */
        void printComponents(double t) const;

        //! The field strength in each dimension
        ElectricFieldData efield_[DIM];
        //! File pointer for electric field
        FILE             *fpField_;
};

//! Converts static parameters from mdp format to E0.
real convertStaticParameters(const std::string &value)
{
    // TODO: Better context for the exceptions here (possibly
    // also convert them to warning_errors or such).
    const std::vector<std::string> sx = splitString(value);
    if (sx.empty())
    {
        return 0.0;
    }
    const int n = fromString<int>(sx[0]);
    if (n <= 0)
    {
        return 0.0;
    }
    if (n != 1)
    {
        GMX_THROW(InvalidInputError("Only one electric field term supported for each dimension"));
    }
    if (sx.size() != 3)
    {
        GMX_THROW(InvalidInputError("Expected exactly one electric field amplitude value"));
    }
    return fromString<real>(sx[1]);
}

//! Converts dynamic parameters from mdp format to (omega, t0, sigma).
void convertDynamicParameters(gmx::KeyValueTreeObjectBuilder *builder,
                              const std::string              &value)
{
    const std::vector<std::string> sxt = splitString(value);
    if (sxt.empty())
    {
        return;
    }
    const int n = fromString<int>(sxt[0]);
    switch (n)
    {
        case 1:
            if (sxt.size() != 3)
            {
                GMX_THROW(InvalidInputError("Please specify 1 omega 0 for non-pulsed fields"));
            }
            builder->addValue<real>("omega", fromString<real>(sxt[1]));
            break;
        case 3:
            if (sxt.size() != 7)
            {
                GMX_THROW(InvalidInputError("Please specify 1 omega 0 t0 0 sigma 0 for pulsed fields"));
            }
            builder->addValue<real>("omega", fromString<real>(sxt[1]));
            builder->addValue<real>("t0", fromString<real>(sxt[3]));
            builder->addValue<real>("sigma", fromString<real>(sxt[5]));
            break;
        default:
            GMX_THROW(InvalidInputError("Incomprehensible input for electric field"));
    }
}

void ElectricField::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    rules->addRule().from<std::string>("/E-x").to<real>("/electric-field/x/E0")
        .transformWith(&convertStaticParameters);
    rules->addRule().from<std::string>("/E-xt").toObject("/electric-field/x")
        .transformWith(&convertDynamicParameters);
    rules->addRule().from<std::string>("/E-y").to<real>("/electric-field/y/E0")
        .transformWith(&convertStaticParameters);
    rules->addRule().from<std::string>("/E-yt").toObject("/electric-field/y")
        .transformWith(&convertDynamicParameters);
    rules->addRule().from<std::string>("/E-z").to<real>("/electric-field/z/E0")
        .transformWith(&convertStaticParameters);
    rules->addRule().from<std::string>("/E-zt").toObject("/electric-field/z")
        .transformWith(&convertDynamicParameters);
}

void ElectricField::initMdpOptions(IOptionsContainerWithSections *options)
{
    auto section = options->addSection(OptionSection("electric-field"));
    efield_[XX].initMdpOptions(&section, "x");
    efield_[YY].initMdpOptions(&section, "y");
    efield_[ZZ].initMdpOptions(&section, "z");
}

void ElectricField::buildMdpOutput(KeyValueTreeObjectBuilder *builder) const
{
    const char *const comment[] = {
        "; Electric fields",
        "; Format for E-x, etc. is: number of cosines (int; only 1 is supported),",
        "; amplitude (real; V/nm), and phase (real; value is meaningless",
        "; for a cosine of frequency 0.",
        "; Format for E-xt, etc. is: omega (1/ps), time for the pulse peak (ps),",
        "; and sigma (ps) width of the pulse. Sigma = 0 removes the pulse,",
        "; leaving the field to be a cosine function."
    };
    builder->addValue<std::string>("comment-electric-field", joinStrings(comment, "\n"));
    efield_[XX].buildMdpOutput(builder, "x");
    efield_[YY].buildMdpOutput(builder, "y");
    efield_[ZZ].buildMdpOutput(builder, "z");
}

void ElectricField::initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                               bool bAppendFiles, const gmx_output_env_t *oenv)
{
    if (isActive())
    {
        please_cite(fplog, "Caleman2008a");

        // Optional outpuf file showing the field, see manual.
        if (opt2bSet("-field", nfile, fnm))
        {
            if (bAppendFiles)
            {
                fpField_ = gmx_fio_fopen(opt2fn("-field", nfile, fnm), "a+");
            }
            else
            {
                fpField_ = xvgropen(opt2fn("-field", nfile, fnm),
                                    "Applied electric field", "Time (ps)",
                                    "E (V/nm)", oenv);
            }
        }
    }
}

void ElectricField::finishOutput()
{
    if (fpField_ != nullptr)
    {
        // This is opened sometimes with xvgropen, sometimes with
        // gmx_fio_fopen, so we use the least common denominator for closing.
        gmx_fio_fclose(fpField_);
        fpField_ = nullptr;
    }
}

real ElectricField::field(int dim, real t) const
{
    return efield_[dim].evaluate(t);
}

bool ElectricField::isActive() const
{
    return (efield_[XX].a() != 0 ||
            efield_[YY].a() != 0 ||
            efield_[ZZ].a() != 0);
}

void ElectricField::printComponents(double t) const
{
    fprintf(fpField_, "%10g  %10g  %10g  %10g\n", t,
            field(XX, t), field(YY, t), field(ZZ, t));
}

void ElectricField::calculateForces(const t_commrec       *cr,
                                    const t_mdatoms       *mdatoms,
                                    const matrix           /* box */,
                                    double                 t,
                                    const rvec             * /* x */,
                                    gmx::ForceWithVirial  *forceWithVirial)
{
    if (isActive())
    {
        // NOTE: The non-conservative electric field does not have a virial
        rvec *f = as_rvec_array(forceWithVirial->force_.data());

        for (int m = 0; (m < DIM); m++)
        {
            real Ext = FIELDFAC*field(m, t);

            if (Ext != 0)
            {
                // TODO: Check parallellism
                for (int i = 0; i < mdatoms->homenr; ++i)
                {
                    // NOTE: Not correct with perturbed charges
                    f[i][m] += mdatoms->chargeA[i]*Ext;
                }
            }
        }
        if (MASTER(cr) && fpField_ != nullptr)
        {
            printComponents(t);
        }
    }
}

}   // namespace

std::unique_ptr<IMDModule> createElectricFieldModule()
{
    return std::unique_ptr<IMDModule>(new ElectricField());
}

} // namespace gmx
