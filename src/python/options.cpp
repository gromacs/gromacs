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

#include "gmxpre.h"

#include "options.h"

#include <map>
#include <string>

#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace pyapi
{

PyOptions::PyOptions() :
    options_ {}
{
    //print_options(*this);
}

PyOptions::PyOptions(std::string filename) :
    options_ {},
filename_ {
    filename
}
{}

PyOptions::~PyOptions()
{}

gmx::Options* PyOptions::data()
{
    // Return a pointer to the Options object
    return &options_;
}


bool PyOptions::parse()
{
    /// Helper class for RAII wrapping of gmx::OptionsAssigner
    class Assigner
    {
        public:
            Assigner(gmx::Options &options) : assigner_ {&options}
            {
                assigner_.start();
            };
            ~Assigner()
            {
                assigner_.finish();
            };

            bool startOption(const std::string &name)
            {
                assigner_.startOption(name.c_str());
                return true;
            };

            bool addSingleValue(const std::string &value)
            {
                try
                {
                    assigner_.appendValue(value.c_str());
                }
                catch (GromacsException &ex)
                {
                    assigner_.finishOption();
                    return false;
                }
                assigner_.finishOption();
                return true;
            };

        private:
            gmx::OptionsAssigner assigner_;
    };
    // use gmx::OptionsAssigner to set parsed options, but refer to
    // CommandLineParser().impl_->assigner for example helper functions.
    // Note that toOptionName just identifies command-line flags and strips
    // them into the form of an option name. We can do this with py dict...
    // E.g.
    //    assigner.start();
    //    // loop over args
    //    const char* const optionName = toOptionName(arg);
    //    // note, arg sets optionName = nullptr implies append to previous option
    //    try
    //        is_valid_name = assigner.tryStartOption(optionName);
    //        // handle invalid options...
    //    catch (UserInputError &ex)
    //        // badly formed option argument
    //    // If instead we are appending:
    //    assigner.appendValue(arg) // throws GromacsException
    //    // conclude single or multivalue options. Multivalue options receive
    //    // additional checking and get another chance to throw.
    //    assigner.finishOption() // throws UserInputError
    //    // end loop
    //    assigner.finish();
    //
    // In the longer term, the options object has a parser available after
    // interacting with the module runner, and can then
    // provide useful user help. It probably doesn't make sense to allow kwargs
    // in a generic and configurable Python class, but maybe the various
    // components can be independently configured or accessed statically to
    // build an argparse-like object. I suppose that's what the help generator
    // does...

    // Scope the lifetime of assigner
    {
        Assigner assigner {
            options_
        };
        // TrajectoryRunnerCommon names the filename option "f"

        const std::string name {
            "f"
        };

        try     // assign a new option
        {
            if (assigner.startOption(name.c_str()))
            {
                // assign (all of) the value(s) for name
                if (!assigner.addSingleValue(filename_.c_str()))
                {
                    // TODO: replace with namespace-scoped API error
                    throw(InvalidInputError("bad option value"));
                }
                // else value successfully appended
            }
            else
            {
                // TODO: replace with namespace-scoped API error
                throw(InvalidInputError("bad option name"));
            };
        }
        catch (InvalidInputError &ex)
        {
            // InvalidInputError is thrown if option is not recognized or inappropriate (e.g. specified more than once)
            throw(ex);
        }
    }
    options_.finish();
    //print_options(*this);
    return true;
}

void PyOptions::view_traverse(gmx::OptionsVisitor &&visitor) const
{
    visitor.visitSection(options_.rootSection());
}

void print_options(const PyOptions &pyoptions)
{
    /// Provide a way to explore the PyOptions object
    /*! gmx::OptionsVisitor defines an interface for visitor classes.
     *  gmx::OptionsIterator provides a means to traverse an Options collection
     *  as a non-member from arbitrary calling code, rather than as a member function of the collection,
     *  which would be more like the Visitor pattern I'm used to.
     *  gmx::OptionsIterator decorates an Options object to provide
     *  acceptSection(OptionsVisitor*) and acceptOptions(OptionsVisitor*)
     *  so that a visitor object should have its visit methods
     *  called directly. E.g. Visitor().visitSection(options.rootSection()) calls
     *  OptionsIteratory(options.rootSection()).acceptSections(visitor) and again
     *  for acceptOptions(visitor). It is not documented whether OptionsIterator.acceptSections(visitor) is made recursive through the Visitor's implementation of visitSection.
     */
    // TODO: IOptionsContainer can only throw APIError and has no way to tell
    // the caller that an option is already defined. However, it is implemented
    // in gmx::Options, so the calling code can get there in a roundabout way.
    // TODO: There is no const version of gmx::OptionsVisitor
    class Visitor : public gmx::OptionsVisitor
    {
        public:
            Visitor() :                 gmx::OptionsVisitor(), observations_ {}
            {};
            std::map<std::string, bool> observations_;
            virtual void visitSection(const gmx::OptionSectionInfo &section)
            {
                // note hierarchy...
                // const std::string name = section.name()
                gmx::OptionsIterator iterator(section);
                iterator.acceptSections(this);
                iterator.acceptOptions(this);
            }
            virtual void visitOption(const OptionInfo &option)
            {
                // Do something...
//            const std::string name   = option.name();
//            const bool        is_set = option.isSet();
                observations_.emplace(option.name(), option.isSet());

                // Can't see values? OptionInfo::option() returns a AbstractOptionStorage& but is a protected function.
                // There does not appear to be a way to get the OptionType (template
                // parameter) settings object used in addOption. OptionType is
                // derived from AbstractOption, I think. Unless there is a default
                // value, there is no option value until the Assigner runs, which
                // operates on a full gmx::Options object. Then the value is owned
                // by the caller of IOptionsContainer.addOption()
                // Where does the reference in AbstractOption.store(T*) end up?
                // Options have a T* store_ that points to the storage defined
                // in the object passed to addOption().
                // OptionSectionImpl::Group::addOptionImpl(AbstractOption& settings) calls
                // settings.createStorage() to get a OptionSectionImpl::AbstractOptionStoragePointer object.
                // When the Assigner gets to appendValue, there is ultimately a
                // commitValues() template method that calls OptionStorageTemplate<T>::store_->append(value). Here, store_ was
                // initialized in the constructor and is a
                // std::unique_ptr<IOptionValueStore<T> >
                // IOptionValueStore<T> is a wrapper that is implemented by various
                // basic ValueStore types to
                // provide some methods and a wrapped ArrayRef<T>.
                // Assigner::Impl uses Options::impl_->rootSection_ to get a
                // internal::OptionSectionImpl which provides findOption() to get
                // an AbstractOptionStorage*. AbstractOptionsStorage->appendValue()
                // results in the values actually being set through virtual functions
                // in the inaccessible OptionStorage object.
                // There appears to be no way to get from either an Options or
                // OptionInfo object to the OptionStorageTemplate object that can
                // see the actual storage. I would need to implement an alternative
                // IOptionsContainer or cause the modules to provide some interface.
                // Also, the parsed option values appear temporarily in the machinery
                // but I think are gone after assignment completes.
                //
                // In short, if I want to see the values being passed now, I can
                // look at the raw memory for the storage destinations, the values
                // being handled by the Assigner, or add printf debugging into the
                // core options handling code...
            }
    };
    pyoptions.view_traverse(Visitor());
}


} // end namespace pyapi
} // end namespace gmx
