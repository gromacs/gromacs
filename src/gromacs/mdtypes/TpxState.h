//
// Created by Eric Irrgang on 8/27/17.
//

#ifndef GMX_MDTYPES_TPXSTATE_H
#define GMX_MDTYPES_TPXSTATE_H

#include <atomic>
#include <memory>
#include <mutex>
#include <string>

struct t_inputrec;
class t_state;
struct gmx_mtop_t;

namespace gmx
{

class TpxState final
{
    private:
        // TpxState is currently always file-backed.
        std::string filename_;
        std::shared_ptr<t_inputrec>                     inputrecInstance_;
        std::shared_ptr<t_state>                        stateInstance_;
        gmx_mtop_t                                     *mtop_;
        std::atomic<bool>                               initialized_;
        std::atomic<bool>                               dirty_;
        mutable std::mutex                              exclusive_;
    public:
        TpxState();
        ~TpxState();

        // Copy semantics TBD. In addition to unclear copy semantics of members, probably need to use setters to allow
        // for notifications of data changs.
        TpxState(const TpxState&) = delete;
        TpxState               &operator=(const TpxState &) = delete;

        // Move should be okay
        TpxState(TpxState && old) noexcept;
        TpxState &operator=(TpxState &&) noexcept;

        static std::unique_ptr<TpxState> initializeFromFile(const char* filename);
        static std::unique_ptr<TpxState> initializeFromFile(const std::string &filename);

        // Takes ownership of arguments to be members of new object.
        static std::unique_ptr<TpxState>
        initializeFromWrappers(std::unique_ptr<t_inputrec> inputRecord,
                               std::unique_ptr<t_state>    state,
                               std::unique_ptr<gmx_mtop_t> mtop);

        t_inputrec* getRawInputrec();
        gmx_mtop_t* getRawMtop();
        t_state*    getRawState();

        /// \returns Whether data has been loaded into the object
        bool isInitialized() const;

        /// \returns true if we do not have a guarantee that the object is in a self-consistent state
        bool isDirty() const;

        /// \brief Allow caller to assert the validity of an instance.
        void markClean();

        const char* filename() const;
};


}      // end namespace gmx

#endif //GROMACS_TPXSTATE_H
