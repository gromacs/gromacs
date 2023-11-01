//
// Created by Kai Zhao on 1/29/21.
//

#ifndef SZ3_PREFILTER_HPP
#define SZ3_PREFILTER_HPP

#include "SZ3/preprocessor/PreProcessor.hpp"

namespace SZ3 {
    template<class T, uint N>

    class PreFilter : public concepts::PreprocessorInterface<T, N> {

        void preprocess(T *data, std::array<size_t, N> dims, std::pair<T, T> range, T defaultValue) {
            for (T &d : data) {
                if (d > range.second || d < range.first) {
                    d = defaultValue;
                }
            }
        }
    };
}
#endif //SZ3_PRETRANSPOSE_H
