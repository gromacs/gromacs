#ifndef _meta_lorenzo_hpp
#define _meta_lorenzo_hpp

#include "SZ3/utils/MetaDef.hpp"

namespace SZMETA {

    template<typename T>
    inline T
    lorenzo_predict_1d(const T *data_pos, size_t dim0_offset) {
        return data_pos[-1];
    }

    template<typename T>
    inline T
    lorenzo_predict_1d_2layer(const T *data_pos, size_t dim0_offset) {
        return 2 * data_pos[-1] - data_pos[-2];
    }


    template<typename T>
    inline T
    lorenzo_predict_2d(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return data_pos[-1] + data_pos[-dim0_offset] - data_pos[-1 - dim0_offset];
    }

    template<typename T>
    inline T
    lorenzo_predict_2d_2layer(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return 2 * data_pos[-dim0_offset]
               - data_pos[-2 * dim0_offset]
               + 2 * data_pos[-1]
               - 4 * data_pos[-1 - dim0_offset]
               + 2 * data_pos[-1 - 2 * dim0_offset]
               - data_pos[-2]
               + 2 * data_pos[-2 - dim0_offset]
               - data_pos[-2 - 2 * dim0_offset];
    }

    template<typename T>
    inline T
    lorenzo_predict_3d(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return data_pos[-1] + data_pos[-dim1_offset] + data_pos[-dim0_offset]
               - data_pos[-dim1_offset - 1] - data_pos[-dim0_offset - 1]
               - data_pos[-dim0_offset - dim1_offset] + data_pos[-dim0_offset - dim1_offset - 1];
    }

    template<typename T>
    inline T
    lorenzo_predict_3d_2layer(const T *data_pos, size_t dim0_offset, size_t dim1_offset) {
        return 2 * data_pos[-1]
               - data_pos[-2]
               + 2 * data_pos[-dim1_offset]
               - 4 * data_pos[-dim1_offset - 1]
               + 2 * data_pos[-dim1_offset - 2]
               - data_pos[-2 * dim1_offset]
               + 2 * data_pos[-2 * dim1_offset - 1]
               - data_pos[-2 * dim1_offset - 2]
               + 2 * data_pos[-dim0_offset]
               - 4 * data_pos[-dim0_offset - 1]
               + 2 * data_pos[-dim0_offset - 2]
               - 4 * data_pos[-dim0_offset - dim1_offset]
               + 8 * data_pos[-dim0_offset - dim1_offset - 1]
               - 4 * data_pos[-dim0_offset - dim1_offset - 2]
               + 2 * data_pos[-dim0_offset - 2 * dim1_offset]
               - 4 * data_pos[-dim0_offset - 2 * dim1_offset - 1]
               + 2 * data_pos[-dim0_offset - 2 * dim1_offset - 2]
               - data_pos[-2 * dim0_offset]
               + 2 * data_pos[-2 * dim0_offset - 1]
               - data_pos[-2 * dim0_offset - 2]
               + 2 * data_pos[-2 * dim0_offset - dim1_offset]
               - 4 * data_pos[-2 * dim0_offset - dim1_offset - 1]
               + 2 * data_pos[-2 * dim0_offset - dim1_offset - 2]
               - data_pos[-2 * dim0_offset - 2 * dim1_offset]
               + 2 * data_pos[-2 * dim0_offset - 2 * dim1_offset - 1]
               - data_pos[-2 * dim0_offset - 2 * dim1_offset - 2];
    }


    template<typename T, class Quantizer>
    inline void
    lorenzo_predict_quantize_3d(const meanInfo<T> &mean_info, const T *data_pos, T *buffer, T precision,
                                T recip_precision, int capacity, int intv_radius,
                                int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                size_t buffer_dim1_offset,
                                size_t dim0_offset, size_t dim1_offset, int *&type_pos, int *unpred_count_buffer,
                                T *unpred_buffer, size_t offset, int padding_layer,
                                bool use_2layer, Quantizer &quantizer, int pred_dim) {

        const T *cur_data_pos = data_pos;
        T *buffer_pos = buffer + padding_layer * (buffer_dim0_offset + buffer_dim1_offset + 1);
        int radius = quantizer.get_radius();
        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {

                    T *cur_buffer_pos = buffer_pos + k;
                    T cur_data = cur_data_pos[k];
                    T pred;
                    if (mean_info.use_mean && fabs(cur_data - mean_info.mean) <= precision) {
                        type_pos[k] = radius;
                        *cur_buffer_pos = mean_info.mean;
                    } else {
                        if (use_2layer) {
                            if (pred_dim == 3) {
                                pred = lorenzo_predict_3d_2layer(cur_buffer_pos, buffer_dim0_offset, buffer_dim1_offset);
                            } else if (pred_dim == 2) {
                                pred = lorenzo_predict_2d_2layer(cur_buffer_pos, buffer_dim0_offset, buffer_dim1_offset);
                            } else {
                                pred = lorenzo_predict_1d_2layer(cur_buffer_pos, buffer_dim0_offset);
                            }
                        } else {
                            if (pred_dim == 3) {
                                pred = lorenzo_predict_3d(cur_buffer_pos, buffer_dim0_offset, buffer_dim1_offset);
                            } else if (pred_dim == 2) {
                                pred = lorenzo_predict_2d(cur_buffer_pos, buffer_dim0_offset, buffer_dim1_offset);
                            } else {
                                pred = lorenzo_predict_1d(cur_buffer_pos, buffer_dim0_offset);
                            }
                        }
//                    *cur_buffer_pos = cur_data;
//                    type_pos[k] = quantizer.quantize_and_overwrite(*cur_buffer_pos, pred);
                        type_pos[k] = quantizer.quantize_and_overwrite(cur_data, pred, *cur_buffer_pos);
                        if (mean_info.use_mean && type_pos[k] >= radius) {
                            type_pos[k] += 1;
                        }
                    }
                }
                type_pos += size_z;
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }

    template<typename T, class Quantizer>
    inline void
    lorenzo_predict_recover_3d(const meanInfo<T> &mean_info, T *buffer, T precision, int intv_radius,
                               int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                               size_t buffer_dim1_offset,
                               size_t dim0_offset, size_t dim1_offset,
                               const int *&type_pos, int *unpred_count_buffer, const T *unpred_data_buffer,
                               const int offset,
                               T *dec_data_pos, const int layer,
                               bool use_2layer, Quantizer &quantizer, int pred_dim) {
        T *cur_data_pos = dec_data_pos;
        T *buffer_pos = buffer + layer * (buffer_dim0_offset + buffer_dim1_offset + 1);
        int radius = quantizer.get_radius();
        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    int index = j * size_z + k;
                    int type_val = type_pos[index];
                    T *cur_buffer_pos = buffer_pos + k;
                    if (type_val == 0) {
                        cur_data_pos[k] = *cur_buffer_pos = quantizer.recover_unpred();
                    } else if (mean_info.use_mean && type_val == radius) {
                        cur_data_pos[k] = *cur_buffer_pos = mean_info.mean;
                    } else {
                        T pred;
//                        pred = predict(cur_buffer_pos, buffer_dim1_offset, buffer_dim0_offset);
                        if (use_2layer) {
                            if (pred_dim == 3) {
                                pred = lorenzo_predict_3d_2layer(cur_buffer_pos, buffer_dim0_offset,
                                                                 buffer_dim1_offset);
                            } else if (pred_dim == 2) {
                                pred = lorenzo_predict_2d_2layer(cur_buffer_pos, buffer_dim0_offset,
                                                                 buffer_dim1_offset);
                            } else {
                                pred = lorenzo_predict_1d_2layer(cur_buffer_pos, buffer_dim0_offset);
                            }
                        } else {
                            if (pred_dim == 3) {
                                pred = lorenzo_predict_3d(cur_buffer_pos, buffer_dim0_offset, buffer_dim1_offset);
                            } else if (pred_dim == 2) {
                                pred = lorenzo_predict_2d(cur_buffer_pos, buffer_dim0_offset, buffer_dim1_offset);
                            } else {
                                pred = lorenzo_predict_1d(cur_buffer_pos, buffer_dim0_offset);
                            }
                        }
                        if (mean_info.use_mean && type_val > radius) {
                            type_val -= 1;
                        }
                        cur_data_pos[k] = *cur_buffer_pos = quantizer.recover_pred(pred, type_val);
                    }

                }
                buffer_pos += buffer_dim1_offset;
                cur_data_pos += dim1_offset;
            }
            type_pos += size_y * size_z;
            buffer_pos += buffer_dim0_offset - size_y * buffer_dim1_offset;
            cur_data_pos += dim0_offset - size_y * dim1_offset;
        }
    }

}
#endif
