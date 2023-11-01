//
// Created by Kai Zhao on 1/29/21.
//

#ifndef SZ3_TRANSPOSE_H
#define SZ3_TRANSPOSE_H

#include "SZ3/preprocessor/PreProcessor.hpp"

namespace SZ3 {
    template<class T, uint N>

    class Transpose : public concepts::PreprocessorInterface<T, N> {
    public:
        void preprocess(T *data, std::array<size_t, N> dims, std::array<size_t, N> axes) {
            static_assert(N < 5, "Data in 5D and above is not supported yet.");
            if (N == 1) {
                return;
            }
            size_t num = dims[N - 1];
            std::array<size_t, N> stride;
            stride[N - 1] = 1;
            for (int i = N - 2; i >= 0; i--) {
                stride[i] = stride[i + 1] * dims[i + 1];
                num *= dims[i];
            }
            std::vector<T> ori(data, data + num);
            transpose(data, dims, axes, ori.data(), stride);
        }

    private:
        template<uint NN = N>
        inline typename std::enable_if<NN == 2, void>::type
        transpose(T *data, std::array<size_t, N> dims, std::array<size_t, N> axes, T *ori,
                  std::array<size_t, N> stride) {
            size_t idx = 0;
            for (size_t i = 0; i < dims[axes[0]]; i++) {
                for (size_t j = 0; j < dims[axes[1]]; j++) {
                    data[idx++] = ori[i * stride[axes[0]]
                                      + j * stride[axes[1]]];
                }
            }
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 3, void>::type
        transpose(T *data, std::array<size_t, N> dims, std::array<size_t, N> axes, T *ori,
                  std::array<size_t, N> stride) {
            size_t idx = 0;
            for (size_t i = 0; i < dims[axes[0]]; i++) {
                for (size_t j = 0; j < dims[axes[1]]; j++) {
                    for (size_t k = 0; k < dims[axes[2]]; k++) {
                        data[idx++] = ori[i * stride[axes[0]]
                                          + j * stride[axes[1]]
                                          + k * stride[axes[2]]];
                    }
                }
            }
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 4, void>::type
        transpose(T *data, std::array<size_t, N> dims, std::array<size_t, N> axes, T *ori,
                  std::array<size_t, N> stride) {
            size_t idx = 0;
            for (size_t i = 0; i < dims[axes[0]]; i++) {
                for (size_t j = 0; j < dims[axes[1]]; j++) {
                    for (size_t k = 0; k < dims[axes[2]]; k++) {
                        for (size_t t = 0; t < dims[axes[3]]; t++) {
                            data[idx++] = ori[i * stride[axes[0]]
                                              + j * stride[axes[1]]
                                              + k * stride[axes[2]]
                                              + t * stride[axes[3]]];
                        }
                    }
                }
            }
        }


    };
}
#endif //SZ3_TRANSPOSE_H
