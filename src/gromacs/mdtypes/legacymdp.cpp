/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Defines data structure for inputrec fields that currently
 * simply map to legacy .mdp fields.
 *
 * \todo Eliminate the need for this handling and data structure.
 *
 * \ingroup module_mdtypes
 */

#include "legacymdp.h"

#include <functional>
#include <memory>
#include <string>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

class LegacyMdp::Impl
{
    public:
        Impl(t_inputrec *ir) : ir_(ir), is_(nullptr), inputrecStrings_(new InputrecStrings),
                               opts_(nullptr), gromppOptions_(new GromppOptions)
        {
            t_inputrec_strings *is;
            snew(is, 1);
            is_.reset(is);
            t_gromppopts *opts;
            snew(opts, 1);
            opts_.reset(opts);
        }
        //! Handle to associated inputrec.
        t_inputrec                      *ir_;
        //! Struct for inputrec strings.
        unique_cptr<t_inputrec_strings>  is_;
        //! Struct for inputrec strings.
        std::unique_ptr<InputrecStrings> inputrecStrings_;
        //! Struct for grompp options.
        unique_cptr<t_gromppopts>        opts_;
        //! Struct for grompp options.
        std::unique_ptr<GromppOptions>   gromppOptions_;
};

LegacyMdp::LegacyMdp(t_inputrec *ir)
    : impl_(new Impl(ir)), is_(impl_->is_.get()), inputrecStrings_(impl_->inputrecStrings_.get()),
      opts_(impl_->opts_.get()), gromppOptions_(impl_->gromppOptions_.get())
{
}

LegacyMdp::~LegacyMdp()
{
}

namespace
{

std::string stringToString(const std::string &value)
{
    return value;
}

using GetEnumValueFromString = std::function<int(const std::string &)>;

// must be null terminated
GetEnumValueFromString stringToEnum(const char **possibleValues)
{
    return [possibleValues](const std::string &s) {
               for (int i = 0; possibleValues[i] != nullptr; ++i)
               {
                   if (gmx_strcasecmp_min(possibleValues[i], s.c_str()) == 0)
                   {
                       return i;
                   }
               }
               GMX_THROW(InvalidInputError(formatString("Invalid value '%s'", s.c_str())));
    };
}

using GetStringFromEnumValue = std::function<std::string(const int &)>;
GetStringFromEnumValue enumToString(const char **possibleValues)
{
    return [possibleValues](const int &value) {
               return possibleValues[value];
    };
}

}   // namespace

void LegacyMdp::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    rules->addRule().from<std::string>("/include").to<std::string>("/legacy/include").transformWith(fromStdString<std::string>);
    rules->addRule().from<std::string>("/define").to<std::string>("/legacy/define").transformWith(fromStdString<std::string>);
    rules->addRule().from<std::string>("/integrator").to<int>("/legacy/integrator").transformWith(stringToEnum(ei_names));
    rules->addRule().from<std::string>("/tinit").to<double>("/legacy/tinit").transformWith(fromStdString<double>);
    rules->addRule().from<std::string>("/init-step").to<gmx_int64_t>("/legacy/init-step").transformWith(fromStdString<gmx_int64_t>);
    rules->addRule().from<std::string>("/nstcomm").to<int>("/legacy/nstcomm").transformWith(fromStdString<int>);
    rules->addRule().from<std::string>("/compressed-x-grps").to<std::string>("/legacy/compressed-x-grps").transformWith(fromStdString<std::string>);
}

void LegacyMdp::initMdpOptions(IOptionsContainerWithSections *options)
{
    options->addOption(StringOption("include").defaultValue("").store(&impl_->gromppOptions_->include));
    options->addOption(StringOption("define").defaultValue("").store(&impl_->gromppOptions_->define));
    options->addOption(EnumOption<int>("integrator").enumValueFromNullTerminatedArray(ei_names).defaultValue(0).store(&impl_->ir_->eI));
    options->addOption(DoubleOption("tinit").defaultValue(0).store(&impl_->ir_->init_t));
    options->addOption(Int64Option("init-step").defaultValue(0).store(&impl_->ir_->init_step));
    options->addOption(IntegerOption("nstcomm").defaultValue(100).store(&impl_->ir_->nstcomm));
    options->addOption(StringOption("compressed-x-grps").defaultValue("").store(&impl_->inputrecStrings_->x_compressed_groups));
}

namespace
{

/*! \brief Make a rule transforming a key to an object that contains a
    comment entry, and the real key-value entry. */
template <typename FromType>
void makeEntry(KeyValueTreeObjectBuilder *builder,
               const KeyValueTreeObject  &object,
               const std::string         &key,
               const std::string         &comment)
{
    if (!comment.empty())
    {
        builder->addValue<std::string>("comment-" + key, comment);
    }
    auto value = object[key].cast<FromType>();
    builder->addValue<std::string>(key, toString<FromType>(value));
}

/*! \brief \copydoc makeEntry()
 *
 * Enums (whether naked ints or enum classes) will continue to need
 * special handling when mapping back to mdp strings. */
void makeEntryForEnum(KeyValueTreeObjectBuilder *builder,
                      const KeyValueTreeObject  &object,
                      const std::string         &key,
                      const char               **possibleValues,
                      const std::string         &comment)
{
    if (!comment.empty())
    {
        builder->addValue<std::string>("comment-" + key, comment);
    }
    auto value = object[key].cast<int>();
    builder->addValue<std::string>(key, possibleValues[value]);
}

void makeMdpEntries(KeyValueTreeObjectBuilder *builder,
                    const KeyValueTreeObject  &legacyObject)
{
    // Unpack the transformed object
    auto object = legacyObject;

    /* TODO Note that until "/legacy" gets some more structure, there
     * is no good way to separate comments that are section headers
     * from the comments for the first key-value pair in that
     * section. But perhaps we'll evolve modules before structuring
     * /legacy, and now the problem might resolve naturally. */
    makeEntry<std::string>(builder, object, "include",
                           "; VARIOUS PREPROCESSING OPTIONS\n"
                           "; Preprocessor information: use cpp syntax.\n"
                           "; e.g.: -I/home/joe/doe -I/home/mary/roe");
    makeEntry<std::string>(builder, object, "define",
                           "e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)");
    makeEntryForEnum      (builder, object, "integrator", ei_names,
                           "; RUN CONTROL PARAMETERS");
    makeEntry<double>     (builder, object, "tinit",
                           "; Start time and timestep in ps");
    makeEntry<gmx_int64_t>(builder, object, "init-step",
                           "; For exact run continuation or redoing part of a run");
    makeEntry<int>        (builder, object, "nstcomm",
                           "; number of steps for center of mass motion removal");
    makeEntry<std::string>(builder, object, "compressed-x-grps",
                           "; This selects the subset of atoms for the compressed\n"
                           "; trajectory file. You can select multiple groups. By\n"
                           "; default, all atoms will be written.");
}

}   // namespace

void LegacyMdp::initMdpBackTransform(IKeyValueTreeTransformRules *rules) const
{
    rules->addRule()
        .from<KeyValueTreeObject>("/legacy")
        .toObject("/legacy")
        .transformWith(&makeMdpEntries);
}

void LegacyMdp::initOutput(FILE */*fplog*/, int /*nfile*/, const t_filenm /*fnm*/[],
                           bool /*bAppendFiles*/, const gmx_output_env_t */*oenv*/)
{
}

void LegacyMdp::finishOutput()
{
}

void LegacyMdp::initForcerec(t_forcerec */*fr*/)
{
}

std::unique_ptr<LegacyMdp> createLegacyMdpModule(t_inputrec *ir)
{
    return std::unique_ptr<LegacyMdp>(new LegacyMdp(ir));
}

} // namespace
