#ifndef FRAME_SET_HPP
#define FRAME_SET_HPP

#include <vector>
#include "utils.hpp"
#include "gromacs/fileio/trxio.h"

// TODO
struct gmx_output_env_t;

namespace FCA{


class FrameSet {
    std::vector< std::unique_ptr< rvec[] > > frames;
    std::vector<real> t_frames;
    matrix box;
    int nbAtoms;

public:
    explicit FrameSet(const char* trajfile,
                      const gmx_output_env_t* oenv,
                      const int load_step = 1)
            : t_frames(0), nbAtoms(0){
        if(trajfile == nullptr){
            return;
        }

        t_trxstatus* status;
        rvec* x = nullptr;
        real t_frame;
        nbAtoms = read_first_x(oenv, &status, trajfile, &t_frame, &x, box);

        int idxFrame = 0;
        do {
            if(idxFrame % load_step == 0) {
                std::unique_ptr< rvec[] > currentFrame(new rvec[nbAtoms]);
                copy_rvecn(x, currentFrame.get(), 0, nbAtoms);
                frames.emplace_back(std::move(currentFrame));
                t_frames.emplace_back(t_frame);
            }
            idxFrame += 1;
        } while(read_next_x(oenv, status, &t_frame, /*nbAtoms,*/ x, box));

        close_trx(status);
    }

    int getNbAtoms() const{
        return nbAtoms;
    }

    const matrix& getBox() const{
        return box;
    }

    real getTimeFrame() const{
        return t_frames.front();
    }

    const std::vector< std::unique_ptr< rvec[] > >& getFrames() const{
        return frames;
    }

    std::vector< std::unique_ptr< rvec[] > > releaseFrames() {
        return std::move(frames);
    }
};

}

#endif
