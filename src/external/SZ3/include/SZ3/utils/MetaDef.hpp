#ifndef _meta_def_hpp
#define _meta_def_hpp

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstddef>

namespace SZMETA {

//    using namespace std;

#define MAX(a, b) ((a<b)?(a):(b))
#define MIN(a, b) ((a<b)?(a):(b))

#define SELECTOR_RADIUS 8
#define SELECTOR_LORENZO 0
#define SELECTOR_REGRESSION 1
#define SELECTOR_LORENZO_2LAYER 2
#define RegCoeffNum3d 4
#define RegErrThreshold 0.1
#define RegCoeffRadius 32768
#define RegCoeffCapacity 65536
#define LorenzeNoise1d 0.5
#define LorenzeNoise2d 0.81
#define LorenzeNoise3d 1.22
#define Lorenze2LayerNoise1d 1.08
#define Lorenze2LayerNoise2d 2.76
#define Lorenze2LayerNoise3d 6.8

    struct meta_params {
        int block_size;
        int prediction_dim;
        bool use_lorenzo;
        bool use_lorenzo_2layer;
        bool use_regression_linear;
        int lorenzo_padding_layer;
        int capacity;
        float regression_param_eb_independent;
        float regression_param_eb_linear;
        float reg_eb_base, reg_eb_1;
        float sample_ratio;
        bool lossless;

        meta_params(bool bi = false, int bs = 6, int pd = 3, int iqi = 0, bool lo = true, bool lo2 = false,
                    bool url = true,
                    float precision = 0.0,
                    float _reg_eb_base = RegErrThreshold, float _reg_eb_1 = -1,
                    int cp = 0,
                    float sr = 1.0, bool ll = true) :
                block_size(bs), prediction_dim(pd),
                use_lorenzo(lo), use_lorenzo_2layer(lo2), use_regression_linear(url),
                capacity(cp), sample_ratio(sr), lossless(ll) {
            lorenzo_padding_layer = 2;
            reg_eb_base = _reg_eb_base;
            reg_eb_1 = _reg_eb_1;
            if (_reg_eb_1 == -1) {
                reg_eb_1 = block_size;
            }
            regression_param_eb_linear = reg_eb_base * precision / RegCoeffNum3d / (float) block_size;
            regression_param_eb_independent = reg_eb_1 * regression_param_eb_linear;

        }
    };


    struct DSize_3d {
        size_t d1;
        size_t d2;
        size_t d3;
        size_t num_elements;
        int block_size;
        int max_num_block_elements;
        size_t num_x;
        size_t num_y;
        size_t num_z;
        size_t num_blocks;
        size_t dim0_offset;
        size_t dim1_offset;
        int sample_nx;
        int sample_ny;
        int sample_nz;
        size_t sample_distance;

        DSize_3d() {}

        DSize_3d(size_t r1, size_t r2, size_t r3, int bs) {
            d1 = r1, d2 = r2, d3 = r3;
            num_elements = r1 * r2 * r3;
            block_size = bs;
            max_num_block_elements = bs * bs * bs;
            num_x = (r1 - 1) / block_size + 1;
            num_y = (r2 - 1) / block_size + 1;
            num_z = (r3 - 1) / block_size + 1;
            num_blocks = num_x * num_y * num_z;
            dim0_offset = r2 * r3;
            dim1_offset = r3;
        }
    };

    template<typename T>
    struct meanInfo {
        bool use_mean;
        T mean;

        meanInfo(bool use = false, T mean_ = 0) {
            use_mean = use;
            mean = mean_;
        }
    };

}
#endif

