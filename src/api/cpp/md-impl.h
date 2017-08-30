#ifndef GMXAPI_MD_IMPL_H
#define GMXAPI_MD_IMPL_H
/*! \file
 * \brief Declarations for molecular dynamics API implementation details.
 *
 * \ingroup gmxapi
 */

#include "gmxapi/md.h"
#include <memory>
#include <string>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

// Container for MD engine input parameters.
class MDInput
{
    public:
        MDInput();

        /// Take over and wrap input data structures.
        MDInput(std::unique_ptr<t_inputrec> &&inputRecord,
                std::unique_ptr<t_state>    &&state,
                std::unique_ptr<gmx_mtop_t> &&topology);

        /// Get input record and state from TPR file.
        static std::unique_ptr<MDInput> from_tpr_file(std::string filename);

        int nAtoms() const;

        /// Return a copy of the KeyValueTreeObject
        gmx::KeyValueTreeObject params() const;

        const t_state* state() const;

        std::unique_ptr<t_inputrec> inputRecord_;
        std::unique_ptr<t_state>    state_;
        std::unique_ptr<gmx_mtop_t> topology_;
};

class EmptyMD : public MDEngine
{

};

/*! \brief Data based on C structs.
 *
 * Provide a ModuleMD implementation that can be created from a MDInput struct
 * of the sort produced by MDInput::from_tpr_file().
 */
class MDStateFromMDInput : public MDEngine
{
    private:
        std::unique_ptr<MDInput> input_;
        std::string              metadata_;

    public:
        MDStateFromMDInput();
        virtual ~MDStateFromMDInput() final;

        explicit MDStateFromMDInput(std::unique_ptr<MDInput> input);
        MDStateFromMDInput(std::unique_ptr<MDInput> input, std::string metadata);
        virtual std::unique_ptr<MDBuilder> builder() override;
        virtual const std::string info() const override;
};

/*!
 * \brief Thin implementation just for holding a tpr filename
 */
class MDStatePlaceholder : public MDEngine
{
    public:
        std::string filename_;
        explicit MDStatePlaceholder(const std::string &filename) : filename_ {filename}
        {};

        const std::string info() const override
        {
            std::string output("MDStatePlaceholder initialized");
            output.append(" with filename: \"");
            output.append(filename_);
            output.append("\"\n");
            return output;
        }

        std::unique_ptr<MDBuilder> builder() override
        {
            class NonBuilder : public MDBuilder
            {
                public:
                    std::string filename_;
                    explicit NonBuilder(const std::string &filename) : filename_ {filename}
                    {};
                    std::unique_ptr<MDEngine> build() override
                    {
                        return nullptr;
                    }

                    std::string inputAsTprFilename() override
                    {
                        return filename_;
                    }
            };
            std::unique_ptr<MDBuilder> builder = gmx::compat::make_unique<NonBuilder>(filename_);
            return builder;
        }
};

}      // namespace gmxapi

#endif // header guard
