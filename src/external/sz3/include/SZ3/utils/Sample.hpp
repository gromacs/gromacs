#ifndef SZ3_SAMPLE_HPP
#define SZ3_SAMPLE_HPP

#include <cassert>
#include <vector>

#include "SZ3/def.hpp"

namespace SZ3 {
template <class T, uint N>
inline void profiling_block(T* data, std::vector<size_t>& dims, std::vector<std::vector<size_t>>& starts,
                            size_t block_size, double abseb, size_t stride = 4) {
    assert(dims.size() == N);
    if (stride == 0) stride = block_size;
    if constexpr (N == 4) {
        size_t dimx = dims[0], dimy = dims[1], dimz = dims[2], dimw = dims[3], dimyzw = dimy * dimz * dimw,
               dimzw = dimz * dimw;
        if (dimx < block_size || dimy < block_size || dimz < block_size || dimw < block_size) {
            return;
        }
        for (size_t i = 0; i < dimx - block_size; i += block_size) {
            for (size_t j = 0; j < dimy - block_size; j += block_size) {
                for (size_t k = 0; k < dimz - block_size; k += block_size) {
                    for (size_t l = 0; l < dimw - block_size; l += block_size) {
                        size_t start_idx = i * dimyzw + j * dimzw + k * dimw + l;
                        T min = data[start_idx];
                        T max = data[start_idx];
                        for (size_t ii = 0; ii <= block_size; ii += stride) {
                            for (size_t jj = 0; jj <= block_size; jj += stride) {
                                for (size_t kk = 0; kk <= block_size; kk += stride) {
                                    for (size_t ll = 0; ll <= block_size; ll += stride) {
                                        size_t cur_idx = start_idx + ii * dimyzw + jj * dimzw + kk * dimw + ll;
                                        T cur_value = data[cur_idx];
                                        if (cur_value < min)
                                            min = cur_value;
                                        else if (cur_value > max)
                                            max = cur_value;
                                    }
                                }
                            }
                        }
                        if (max - min > abseb) {
                            size_t a[4] = {i, j, k, l};
                            starts.push_back(std::vector<size_t>(a, a + 4));
                        }
                    }
                }
            }
        }
    } else if constexpr (N == 3) {
        size_t dimx = dims[0], dimy = dims[1], dimz = dims[2], dimyz = dimy * dimz;
        if (dimx < block_size || dimy < block_size || dimz < block_size) {
            return;
        }
        for (size_t i = 0; i < dimx - block_size; i += block_size) {
            for (size_t j = 0; j < dimy - block_size; j += block_size) {
                for (size_t k = 0; k < dimz - block_size; k += block_size) {
                    size_t start_idx = i * dimyz + j * dimz + k;
                    T min = data[start_idx];
                    T max = data[start_idx];
                    for (size_t ii = 0; ii <= block_size; ii += stride) {
                        for (size_t jj = 0; jj <= block_size; jj += stride) {
                            for (size_t kk = 0; kk <= block_size; kk += stride) {
                                size_t cur_idx = start_idx + ii * dimyz + jj * dimz + kk;
                                T cur_value = data[cur_idx];
                                if (cur_value < min)
                                    min = cur_value;
                                else if (cur_value > max)
                                    max = cur_value;
                            }
                        }
                    }
                    if (max - min > abseb) {
                        size_t a[3] = {i, j, k};
                        starts.push_back(std::vector<size_t>(a, a + 3));
                    }
                }
            }
        }
    } else if constexpr (N == 2) {
        size_t dimx = dims[0], dimy = dims[1];
        if (dimx < block_size || dimy < block_size) {
            return;
        }
        for (size_t i = 0; i < dimx - block_size; i += block_size) {
            for (size_t j = 0; j < dimy - block_size; j += block_size) {
                size_t start_idx = i * dimy + j;
                T min = data[start_idx];
                T max = data[start_idx];
                for (size_t ii = 0; ii <= block_size; ii += stride) {
                    for (size_t jj = 0; jj <= block_size; jj += stride) {
                        size_t cur_idx = start_idx + ii * dimy + jj;
                        T cur_value = data[cur_idx];
                        if (cur_value < min)
                            min = cur_value;
                        else if (cur_value > max)
                            max = cur_value;
                    }
                }
                if (max - min > abseb) {
                    size_t a[2] = {i, j};
                    starts.push_back(std::vector<size_t>(a, a + 2));
                }
            }
        }
    } else {
        size_t dimx = dims[0];
        if (dimx < block_size) {
            return;
        }
        for (size_t i = 0; i < dimx - block_size; i += block_size) {
            size_t start_idx = i;
            T min = data[start_idx];
            T max = data[start_idx];
            for (size_t ii = 0; ii <= block_size; ii += stride) {
                size_t cur_idx = start_idx + ii;
                T cur_value = data[cur_idx];
                if (cur_value < min)
                    min = cur_value;
                else if (cur_value > max)
                    max = cur_value;
            }
            if (max - min > abseb) {
                size_t a[1] = {i};
                starts.push_back(std::vector<size_t>(a, a + 1));
            }
        }
    }
}

template <class T, uint N>
inline void sample_blocks(T* data, std::vector<T>& sampling_data, std::vector<size_t>& dims,
                          std::vector<size_t>& starts, size_t block_size) {
    assert(dims.size() == N);
    assert(starts.size() == N);
    if constexpr (N == 4) {
        if (dims[0] < block_size || dims[1] < block_size || dims[2] < block_size || dims[3] < block_size) {
            return;
        }
        size_t sample_num = block_size * block_size * block_size * block_size;
        sampling_data.resize(sample_num, 0);
        size_t startx = starts[0], starty = starts[1], startz = starts[2], startw = starts[3], dimy = dims[1],
               dimz = dims[2], dimw = dims[3];
        size_t cubic_block_size = block_size * block_size * block_size, square_block_size = block_size * block_size,
               dimyzw = dimy * dimz * dimw, dimzw = dimz * dimw;
        for (size_t i = 0; i < block_size; i++) {
            for (size_t j = 0; j < block_size; j++) {
                for (size_t k = 0; k < block_size; k++) {
                    for (size_t l = 0; l < block_size; l++) {
                        size_t sample_idx = i * cubic_block_size + j * square_block_size + k * block_size + l;
                        size_t idx = (i + startx) * dimyzw + (j + starty) * dimzw + (k + startz) * dimw + (l + startw);
                        sampling_data[sample_idx] = data[idx];
                    }
                }
            }
        }
    } else if constexpr (N == 3) {
        if (dims[0] < block_size || dims[1] < block_size || dims[2] < block_size) {
            return;
        }
        size_t sample_num = block_size * block_size * block_size;
        sampling_data.resize(sample_num, 0);
        size_t startx = starts[0], starty = starts[1], startz = starts[2], dimy = dims[1], dimz = dims[2];
        size_t square_block_size = block_size * block_size, dimyz = dimy * dimz;
        for (size_t i = 0; i < block_size; i++) {
            for (size_t j = 0; j < block_size; j++) {
                for (size_t k = 0; k < block_size; k++) {
                    size_t sample_idx = i * square_block_size + j * block_size + k;
                    size_t idx = (i + startx) * dimyz + (j + starty) * dimz + k + startz;
                    sampling_data[sample_idx] = data[idx];
                }
            }
        }
    } else if constexpr (N == 2) {
        if (dims[0] < block_size || dims[1] < block_size) {
            return;
        }
        size_t sample_num = block_size * block_size;
        sampling_data.resize(sample_num, 0);
        size_t startx = starts[0], starty = starts[1], dimy = dims[1];
        for (size_t i = 0; i < block_size; i++) {
            for (size_t j = 0; j < block_size; j++) {
                size_t sample_idx = i * block_size + j;
                size_t idx = (i + startx) * dimy + (j + starty);
                sampling_data[sample_idx] = data[idx];
            }
        }
    } else if constexpr (N == 1) {
        if (dims[0] < block_size) {
            return;
        }
        size_t sample_num = block_size;
        sampling_data.resize(sample_num, 0);
        size_t startx = starts[0];
        for (size_t i = 0; i < block_size; i++) {
            size_t sample_idx = i;
            size_t idx = (i + startx);
            sampling_data[sample_idx] = data[idx];
        }
    }
}

template <class T, uint N>
void sampleBlocks(T* data, std::vector<size_t>& dims, size_t sampleBlockSize,
                  std::vector<std::vector<T>>& sampled_blocks, double sample_rate, int profiling,
                  std::vector<std::vector<size_t>>& starts) {
    for (uint i = 0; i < N; i++) {
        if (dims[i] < sampleBlockSize) {
            return;
        }
    }
    for (uint i = 0; i < sampled_blocks.size(); i++) {
        std::vector<T>().swap(sampled_blocks[i]);
    }
    std::vector<std::vector<T>>().swap(sampled_blocks);
    for (uint i = 0; i < sampled_blocks.size(); i++) {
        std::vector<T>().swap(sampled_blocks[i]);
    }
    std::vector<std::vector<T>>().swap(sampled_blocks);
    size_t totalblock_num = 1;
    for (uint i = 0; i < N; i++) {
        totalblock_num *= static_cast<int>((dims[i] - 1) / sampleBlockSize);
    }
    size_t idx = 0;
    if (profiling) {
        size_t num_filtered_blocks = starts.size();
        size_t sample_stride = static_cast<size_t>(num_filtered_blocks / (totalblock_num * sample_rate));
        if (sample_stride <= 0) sample_stride = 1;
        for (size_t i = 0; i < num_filtered_blocks; i += sample_stride) {
            std::vector<T> s_block;
            sample_blocks<T, N>(data, s_block, dims, starts[i], sampleBlockSize + 1);
            sampled_blocks.push_back(s_block);
        }
    } else {
        size_t sample_stride = static_cast<size_t>(1.0 / sample_rate);
        if (sample_stride <= 0) sample_stride = 1;
        if constexpr (N == 1) {
            for (size_t x_start = 0; x_start < dims[0] - sampleBlockSize; x_start += sampleBlockSize) {
                if (idx % sample_stride == 0) {
                    std::vector<size_t> block_start{x_start};
                    std::vector<T> s_block;
                    sample_blocks<T, N>(data, s_block, dims, block_start, sampleBlockSize + 1);
                    sampled_blocks.push_back(s_block);
                }
                idx += 1;
            }
        } else if constexpr (N == 2) {
            for (size_t x_start = 0; x_start < dims[0] - sampleBlockSize; x_start += sampleBlockSize) {
                for (size_t y_start = 0; y_start < dims[1] - sampleBlockSize; y_start += sampleBlockSize) {
                    if (idx % sample_stride == 0) {
                        std::vector<size_t> block_start{x_start, y_start};
                        std::vector<T> s_block;
                        sample_blocks<T, N>(data, s_block, dims, block_start, sampleBlockSize + 1);
                        sampled_blocks.push_back(s_block);
                    }
                    idx += 1;
                }
            }
        } else if constexpr (N == 3) {
            for (size_t x_start = 0; x_start < dims[0] - sampleBlockSize; x_start += sampleBlockSize) {
                for (size_t y_start = 0; y_start < dims[1] - sampleBlockSize; y_start += sampleBlockSize) {
                    for (size_t z_start = 0; z_start < dims[2] - sampleBlockSize; z_start += sampleBlockSize) {
                        if (idx % sample_stride == 0) {
                            std::vector<size_t> block_start{x_start, y_start, z_start};
                            std::vector<T> s_block;
                            sample_blocks<T, N>(data, s_block, dims, block_start, sampleBlockSize + 1);
                            sampled_blocks.push_back(s_block);
                        }
                        idx += 1;
                    }
                }
            }
        } else if constexpr (N == 4) {
            for (size_t x_start = 0; x_start < dims[0] - sampleBlockSize; x_start += sampleBlockSize) {
                for (size_t y_start = 0; y_start < dims[1] - sampleBlockSize; y_start += sampleBlockSize) {
                    for (size_t z_start = 0; z_start < dims[2] - sampleBlockSize; z_start += sampleBlockSize) {
                        for (size_t w_start = 0; w_start < dims[3] - sampleBlockSize; w_start += sampleBlockSize) {
                            if (idx % sample_stride == 0) {
                                std::vector<size_t> block_start{x_start, y_start, z_start, w_start};
                                std::vector<T> s_block;
                                sample_blocks<T, N>(data, s_block, dims, block_start, sampleBlockSize + 1);
                                sampled_blocks.push_back(s_block);
                            }
                            idx += 1;
                        }
                    }
                }
            }
        }
    }
}
} // namespace SZ3

#endif
