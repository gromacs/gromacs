#ifndef _meta_regression_hpp
#define _meta_regression_hpp

#include "SZ3/utils/MetaDef.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZMETA {

    template<typename T>
    inline int
    quantize_reg_coeff(float pred, T cur_data, double precision, double recip_precision, int capacity, int intv_radius,
                       T *&unpredictable_data_pos, T *decompressed) {
        double diff = cur_data - pred;
        double quant_diff = fabs(diff) * recip_precision + 1;
        if (quant_diff < capacity) {
            quant_diff = (diff > 0) ? quant_diff : -quant_diff;
            int quant_index = (int) (quant_diff / 2) + intv_radius;
            T decompressed_data = pred + 2 * (quant_index - intv_radius) * precision;
            *decompressed = decompressed_data;
            if (fabs(decompressed_data - cur_data) <= precision) return quant_index;
        }
        *decompressed = cur_data;
        *(unpredictable_data_pos++) = cur_data;
        return 0;
    }

    template<typename T>
    inline T
    recover_reg_coeff(float pred, double precision, int type_val, int intv_radius, const T *&unpredictable_data_pos) {
        if (type_val == 0) {
            return *(unpredictable_data_pos++);
        } else {
            return pred + 2 * (type_val - intv_radius) * precision;
        }
    }

    template<typename T>
    inline void
    compress_regression_coefficient_3d(const int coeff_num, const T *precisions, const T *recip_precisions,
                                       float *reg_params_pos,
                                       int *reg_params_type_pos, float *&reg_unpredictable_data_pos) {
        float *prev_reg_params = reg_params_pos - coeff_num;
        for (int i = 0; i < coeff_num; i++) {
            *(reg_params_type_pos++) = quantize_reg_coeff(*prev_reg_params, *reg_params_pos, precisions[i],
                                                          recip_precisions[i],
                                                          RegCoeffCapacity, RegCoeffRadius, reg_unpredictable_data_pos,
                                                          reg_params_pos);
            prev_reg_params++, reg_params_pos++;
        }
    }

    template<typename T>
    float *
    decode_regression_coefficients(const unsigned char *&compressed_pos, size_t reg_count, int block_size, T precision,
                                   const meta_params &params) {
        size_t reg_unpredictable_count = 0;
        size_t remaining_length = RegCoeffNum3d * reg_count;//fake, has no meaning
        SZ3::read(reg_unpredictable_count, compressed_pos, remaining_length);
        const float *reg_unpredictable_data_pos = (const float *) compressed_pos;
        compressed_pos += reg_unpredictable_count * sizeof(float);

//        int *reg_type = Huffman_decode_tree_and_data(2 * RegCoeffCapacity, RegCoeffNum3d * reg_count, compressed_pos);
        SZ3::HuffmanEncoder<int> selector_encoder = SZ3::HuffmanEncoder<int>();
        selector_encoder.load(compressed_pos, remaining_length);
        auto reg_vector = selector_encoder.decode(compressed_pos, RegCoeffNum3d * reg_count);
        selector_encoder.postprocess_decode();
        int *reg_type = reg_vector.data();

        float *reg_params = (float *) malloc(RegCoeffNum3d * (reg_count + 1) * sizeof(float));
        for (int i = 0; i < RegCoeffNum3d; i++)
            reg_params[i] = 0;
        T reg_precisions[RegCoeffNum3d];
        for (int i = 0; i < RegCoeffNum3d - 1; i++) {
            reg_precisions[i] = params.regression_param_eb_linear;
        }
        reg_precisions[RegCoeffNum3d - 1] = params.regression_param_eb_independent;
        float *prev_reg_params = reg_params;
        float *reg_params_pos = reg_params + RegCoeffNum3d;
        const int *type_pos = (const int *) reg_type;
        for (int i = 0; i < reg_count; i++) {
            for (int j = 0; j < RegCoeffNum3d; j++) {
                *reg_params_pos = recover_reg_coeff(*prev_reg_params, reg_precisions[j], *(type_pos++), RegCoeffRadius,
                                                    reg_unpredictable_data_pos);
                prev_reg_params++, reg_params_pos++;
            }
        }
//        free(reg_type);
        return reg_params;
    }

    void
    encode_regression_coefficients(const int *reg_params_type, const float *reg_unpredictable_data, size_t reg_count,
                                   size_t reg_unpredictable_count, SZ3::HuffmanEncoder<int> &reg_huffman, unsigned char *&compressed_pos) {
        SZ3::write(reg_unpredictable_count, compressed_pos);
        SZ3::write(reg_unpredictable_data, reg_unpredictable_count, compressed_pos);

//        reg_huffman.preprocess_encode(reg_params_type, reg_count, 0);
        reg_huffman.save(compressed_pos);
        reg_huffman.encode(reg_params_type, reg_count, compressed_pos);
        reg_huffman.postprocess_encode();
//        Huffman_encode_tree_and_data(2 * RegCoeffCapacity, reg_params_type, reg_count, compressed_pos);
    }

    template<typename T>
    inline void
    compute_regression_coeffcients_3d(const T *data_pos, int size_x, int size_y, int size_z, size_t dim0_offset,
                                      size_t dim1_offset, float *reg_params_pos) {
        /*Calculate regression coefficients*/
        const T *cur_data_pos = data_pos;
        float fx = 0.0;
        float fy = 0.0;
        float fz = 0.0;
        float f = 0;
        float sum_x, sum_y;
        T curData;
        for (int i = 0; i < size_x; i++) {
            sum_x = 0;
            for (int j = 0; j < size_y; j++) {
                sum_y = 0;
                for (int k = 0; k < size_z; k++) {
                    curData = *cur_data_pos;
                    sum_y += curData;
                    fz += curData * k;
                    cur_data_pos++;
                }
                fy += sum_y * j;
                sum_x += sum_y;
                cur_data_pos += (dim1_offset - size_z);
            }
            fx += sum_x * i;
            f += sum_x;
            cur_data_pos += (dim0_offset - size_y * dim1_offset);
        }
        float coeff = 1.0 / (size_x * size_y * size_z);
        reg_params_pos[0] = (2 * fx / (size_x - 1) - f) * 6 * coeff / (size_x + 1);
        reg_params_pos[1] = (2 * fy / (size_y - 1) - f) * 6 * coeff / (size_y + 1);
        reg_params_pos[2] = (2 * fz / (size_z - 1) - f) * 6 * coeff / (size_z + 1);
        reg_params_pos[3] = f * coeff - ((size_x - 1) * reg_params_pos[0] / 2 + (size_y - 1) * reg_params_pos[1] / 2 +
                                         (size_z - 1) * reg_params_pos[2] / 2);
    }


    template<typename T>
    inline T
    regression_predict_3d(const float *reg_params_pos, int x, int y, int z) {
        return reg_params_pos[0] * x + reg_params_pos[1] * y + reg_params_pos[2] * z + reg_params_pos[3];
    }

    template<typename T, class Quantizer>
    inline void
    regression_predict_quantize_3d(const T *data_pos, const float *reg_params_pos, T *buffer, T precision,
                                   T recip_precision, int capacity,
                                   int intv_radius, int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                   size_t buffer_dim1_offset,
                                   size_t dim0_offset, size_t dim1_offset, int *&type_pos,
                                   int *unpred_count_buffer, T *unpred_buffer, size_t offset, int lorenzo_layer,
                                   Quantizer &quantizer) {
        for (int i = 0; i < size_x; i++) {
//        const T *cur_data_pos = data_pos + i * dim0_offset;
            // T * buffer_pos = buffer + (i+1)*buffer_dim0_offset + buffer_dim1_offset + 1;
            T *buffer_pos =
                    buffer + (i + lorenzo_layer) * buffer_dim0_offset + lorenzo_layer * buffer_dim1_offset +
                    lorenzo_layer;
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    T cur_data = data_pos[i * dim0_offset + j * dim1_offset + k];
                    T pred = (T) (reg_params_pos[0] * (float) i + reg_params_pos[1] * (float) j +
                                  reg_params_pos[2] * (float) k +
                                  reg_params_pos[3]);

//                    buffer_pos[j * buffer_dim1_offset + k] = cur_data;
//                    type_pos[j * size_z + k] = quantizer.quantize_and_overwrite(buffer_pos[j * buffer_dim1_offset + k], pred);

                    type_pos[j * size_z + k] = quantizer.quantize_and_overwrite(cur_data, pred,
                                                                                buffer_pos[j * buffer_dim1_offset +
                                                                                           k]);

//                    T diff = cur_data - pred;
//                    int quant_index = (int) (fabs(diff) * recip_precision) + 1;
//                    if (quant_index < capacity) {
//                        quant_index >>= 1;
//                        int half_index = quant_index;
//                        quant_index <<= 1;
//                        if (diff < 0) {
//                            quant_index = -quant_index;
//                            type_pos[j * size_z + k] = intv_radius - half_index;
//                        } else type_pos[j * size_z + k] = intv_radius + half_index;
//                        T decompressed_data = pred + (T) quant_index * precision;
//                        if (fabs(cur_data - decompressed_data) > precision) {
//                            int index = j * size_z + k;
//                            type_pos[index] = 0;
//                            unpred_buffer[index * offset + unpred_count_buffer[index]] = buffer_pos[
//                                    j * buffer_dim1_offset +
//                                    k] = cur_data;
//                            unpred_count_buffer[index]++;
//                        } else buffer_pos[j * buffer_dim1_offset + k] = decompressed_data;
//                    } else {
//                        int index = j * size_z + k;
//                        type_pos[index] = 0;
//                        unpred_buffer[index * offset + unpred_count_buffer[index]] = buffer_pos[j * buffer_dim1_offset +
//                                                                                                k] = cur_data;
//                        unpred_count_buffer[index]++;
//                    }
                }
            }
            type_pos += size_y * size_z;
        }
//    exit(0);

    }

    template<typename T, class Quantizer>
    inline void
    regression_predict_recover_3d(const float *reg_params_pos, T *buffer, T precision, int intv_radius,
                                  int size_x, int size_y, int size_z, size_t buffer_dim0_offset,
                                  size_t buffer_dim1_offset,
                                  size_t dim0_offset, size_t dim1_offset, const int *&type_pos,
                                  int *unpred_count_buffer,
                                  const T *unpred_data_buffer, const int offset, T *dec_data_pos,
                                  int lorenzo_layer, Quantizer &quantizer) {
        T *cur_data_pos = dec_data_pos;
        T *buffer_pos = buffer + lorenzo_layer * (buffer_dim0_offset + buffer_dim1_offset + 1);
        for (int i = 0; i < size_x; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_z; k++) {
                    int index = j * size_z + k;
                    int type_val = type_pos[index];
                    if (type_val == 0) {
                        cur_data_pos[j * dim1_offset + k] = buffer_pos[j * buffer_dim1_offset + k]
                                = quantizer.recover_unpred();
//                        cur_data_pos[j * dim1_offset + k] = buffer_pos[j * buffer_dim1_offset + k] = unpred_data_buffer[
//                                index * offset + unpred_count_buffer[index]];
//                        unpred_count_buffer[index]++;
                    } else {

                        T pred = (T) (reg_params_pos[0] * (float) i + reg_params_pos[1] * (float) j +
                                      reg_params_pos[2] * (float) k +
                                      reg_params_pos[3]);
                        cur_data_pos[j * dim1_offset + k] = buffer_pos[j * buffer_dim1_offset +
                                                                       k] = quantizer.recover_pred(pred, type_val);
//                        cur_data_pos[j * dim1_offset + k] = buffer_pos[j * buffer_dim1_offset + k] =
//                                pred + (T) (2 * (type_val - intv_radius)) * precision;
                    }
                }
            }
            type_pos += size_y * size_z;
            cur_data_pos += dim0_offset;
            buffer_pos += buffer_dim0_offset;
        }
    }
}
#endif
