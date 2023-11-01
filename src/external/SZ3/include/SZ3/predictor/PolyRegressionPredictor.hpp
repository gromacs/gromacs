#ifndef _SZ_POLY_REGRESSION_PREDICTOR_HPP
#define _SZ_POLY_REGRESSION_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "PolyRegressionCoeffAux.hpp"
#include <cstring>
#include <iostream>

namespace SZ3 {

    // N-d regression predictor
    template<class T, uint N, uint M = (N + 1) * (N + 2) / 2>
    class PolyRegressionPredictor : public concepts::PredictorInterface<T, N> {
    public:
        static const uint8_t predictor_id = 0b00000011;

//        PolyRegressionPredictor() : quantizer_independent(0), quantizer_liner(0), quantizer_poly(0), current_coeffs{0} {
//            init_poly();
//        }

        PolyRegressionPredictor(uint block_size, T eb) : quantizer_independent(eb / 5 / block_size),
                                                         quantizer_liner(eb / 20 / block_size),
                                                         quantizer_poly(eb / 100 / block_size),
                                                         prev_coeffs{0}, current_coeffs{0} {
            init_poly(block_size);
        }

//        PolyRegressionPredictor(T eb1, T eb2, T eb3) : quantizer_independent(eb1),
//                                                       quantizer_liner(eb2),
//                                                       quantizer_poly(eb3),
//                                                       prev_coeffs{0}, current_coeffs{0} {
//            init_poly();
//        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        void precompress_data(const iterator &) const noexcept {}

        void postcompress_data(const iterator &) const noexcept {}

        void predecompress_data(const iterator &) const noexcept {}

        void postdecompress_data(const iterator &) const noexcept {}

        inline T estimate_error(const iterator &iter) const noexcept {
            return fabs(*iter - predict(iter));
        }

        bool precompress_block(const std::shared_ptr<Range> &range) noexcept {
            auto dims = range->get_dimensions();
            for (const auto &dim: dims) {
                if (dim <= 2) {
                    return false;
                }
            }
            std::array<double, M> sum{0};
            {
                for (auto iter = range->begin(); iter != range->end(); ++iter) {
                    T data = *iter;
                    auto poly_index = get_poly_index<N>(iter);
                    for (int i = 0; i < M; i++) {
                        sum[i] += poly_index[i] * data;
                    }
                }
            }
//            std::array<double, M> current_coeffs{0};
            std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
            auto coef_aux = coef_aux_list[get_coef_aux_list_idx(dims)];

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    current_coeffs[i] += coef_aux[i * M + j] * sum[j];
                }
            }
            return true;
        }

        void precompress_block_commit() noexcept {
            pred_and_quantize_coefficients();
            std::copy(current_coeffs.begin(), current_coeffs.end(), prev_coeffs.begin());
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 1, std::array<double, M>>::type get_poly_index(const iterator &iter) const {
            double i = iter.get_local_index(0);

            return std::array<double, M>{1, i, i * i};
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 2, std::array<double, M>>::type get_poly_index(const iterator &iter) const {
            double i = iter.get_local_index(0);
            double j = iter.get_local_index(1);

            return std::array<double, M>{1, i, j, i * i, i * j, j * j};
        }

        template<uint NN = N>
        inline typename std::enable_if<NN != 1 && NN != 2, std::array<double, M>>::type
        get_poly_index(const iterator &iter) const {
            double i = iter.get_local_index(0);
            double j = iter.get_local_index(1);
            double k = iter.get_local_index(2);

            return std::array<double, M>{1, i, j, k, i * i, i * j, i * k, j * j, j * k, k * k};
        }

        inline T predict(const iterator &iter) const noexcept {
            T pred = 0;
            auto poly_index = get_poly_index<N>(iter);
            for (int i = 0; i < M; i++) {
                pred += poly_index[i] * current_coeffs[i];
            }
            return pred;
        }

        void save(uchar *&c) const {
            c[0] = predictor_id;
            c += 1;
            *reinterpret_cast<size_t *>(c) = regression_coeff_quant_inds.size();
            c += sizeof(size_t);
            if (!regression_coeff_quant_inds.empty()) {
                quantizer_independent.save(c);
                quantizer_liner.save(c);
                quantizer_poly.save(c);
                HuffmanEncoder<int> encoder = HuffmanEncoder<int>();

                encoder.preprocess_encode(regression_coeff_quant_inds,
                                          2 * std::max(std::max(quantizer_independent.get_radius(), quantizer_liner.get_radius()),
                                                       quantizer_poly.get_radius()));
                encoder.save(c);
                encoder.encode(regression_coeff_quant_inds, c);
                encoder.postprocess_encode();
            }
        }

        bool predecompress_block(const std::shared_ptr<Range> &range) noexcept {
            for (const auto &dim: range->get_dimensions()) {
                if (dim <= 2) {
                    return false;
                }
            }
            pred_and_recover_coefficients();
            return true;
        }

        void load(const uchar *&c, size_t &remaining_length) {
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            size_t coeff_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            remaining_length -= sizeof(size_t);
            if (coeff_size != 0) {
                quantizer_independent.load(c, remaining_length);
                quantizer_liner.load(c, remaining_length);
                quantizer_poly.load(c, remaining_length);
                HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
                encoder.load(c, remaining_length);
                regression_coeff_quant_inds = encoder.decode(c, coeff_size);
                encoder.postprocess_decode();
            }
            std::fill(current_coeffs.begin(), current_coeffs.end(), 0);
            regression_coeff_index = 0;
        }

        void print() const {
            std::cout << "2-Layer Regression predictor, indendent term eb = " << quantizer_independent.get_eb() << "\n";
            std::cout << "2-Layer Regression predictor, linear term eb = " << quantizer_liner.get_eb() << "\n";
            std::cout << "2-Layer Regression predictor, poly term eb = " << quantizer_poly.get_eb() << "\n";
        }

        void clear() {
            quantizer_independent.clear();
            quantizer_liner.clear();
            quantizer_poly.clear();
            regression_coeff_quant_inds.clear();
            regression_coeff_index = 0;
            current_coeffs = {0};
            prev_coeffs = {0};
        }

    private:
        LinearQuantizer<T> quantizer_independent, quantizer_liner, quantizer_poly;
        std::vector<int> regression_coeff_quant_inds;
        size_t regression_coeff_index = 0;
        std::array<T, M> current_coeffs;
        std::array<T, M> prev_coeffs;
        std::vector<std::array<T, M * M>> coef_aux_list;
        const std::vector<int> COEF_AUX_MAX_BLOCK = {5000, 4096, 64, 16};//The maximum block size supported by PolyReg.

        void init_poly(size_t block_size) {
            float *data;
            size_t num;
            if (N == 1) {
                data = COEFF_1D;
                num = sizeof(COEFF_1D) / sizeof(float);
            } else if (N == 2) {
                data = COEFF_2D;
                num = sizeof(COEFF_2D) / sizeof(float);
            } else if (N == 3) {
                data = COEFF_3D;
                num = sizeof(COEFF_3D) / sizeof(float);
            } else {
                printf("Poly regression only supports 1D, 2D, and 3D datasets.\n");
                exit(1);
            }
            if (block_size > COEF_AUX_MAX_BLOCK[N]) {
                printf("%dD Poly regression supports block size upto %d\n.", N, COEF_AUX_MAX_BLOCK[N]);
                exit(1);
            }

            auto coef_aux_p = &data[0];
            coef_aux_list = std::vector<std::array<T, M * M>>(COEF_AUX_MAX_BLOCK[0], {0});
            while (coef_aux_p < &data[0] + num) {
                std::array<size_t, N> dims;
                for (auto &idx: dims) {
                    idx = *coef_aux_p++;
                }
                std::copy_n(coef_aux_p, M * M, coef_aux_list[get_coef_aux_list_idx(dims)].begin());
                coef_aux_p += M * M;
            }
//            display_coef_aux(coef_aux_list[get_coef_aux_list_idx( std::array<size_t, N>{6, 6, 6})]);
        }

//        std::array<T, M>
//        compute_regression_coefficients(const std::shared_ptr<Range> &range, const std::array<size_t, N> &dims) const {
//            std::array<double, M> sum{0};
//            {
//                for (auto iter = range->begin(); iter != range->end(); ++iter) {
//                    T data = *iter;
//                    auto poly_index = get_poly_index<N>(iter);
//                    for (int i = 0; i < M; i++) {
//                        sum[i] += poly_index[i] * data;
//                    }
//                }
//            }
//            std::array<double, M> coeffs{0};
//            auto coef_aux_idx = coef_aux_list[get_coef_aux_list_idx(dims)];
//
//            for (int i = 0; i < M; i++) {
//                for (int j = 0; j < M; j++) {
//                    coeffs[i] += coef_aux_idx[i * M + j] * sum[j];
//                }
//            }
//
//            std::array<T, M> coeffsT;
//            for (int i = 0; i < M; i++) {
//                coeffsT[i] = coeffs[i];
//            }
//            return coeffsT;
//        }

        void pred_and_quantize_coefficients() {
            regression_coeff_quant_inds.push_back(
                    quantizer_independent.quantize_and_overwrite(current_coeffs[0], prev_coeffs[0]));
            for (int i = 1; i < N + 1; i++) {
                regression_coeff_quant_inds.push_back(
                        quantizer_liner.quantize_and_overwrite(current_coeffs[i], prev_coeffs[i]));
            }
            for (int i = N + 1; i < M; i++) {
                regression_coeff_quant_inds.push_back(
                        quantizer_poly.quantize_and_overwrite(current_coeffs[i], prev_coeffs[i]));
            }
        }

        void pred_and_recover_coefficients() {
            current_coeffs[0] = quantizer_independent.recover(current_coeffs[0],
                                                              regression_coeff_quant_inds[regression_coeff_index++]);
            for (int i = 1; i < N + 1; i++) {
                current_coeffs[i] = quantizer_liner.recover(current_coeffs[i],
                                                            regression_coeff_quant_inds[regression_coeff_index++]);
            }
            for (int i = N + 1; i < M; i++) {
                current_coeffs[i] = quantizer_poly.recover(current_coeffs[i],
                                                           regression_coeff_quant_inds[regression_coeff_index++]);
            }
        }

        void display_coef_aux(std::array<T, M * M> aux) {
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++)
                    std::cout << std::setw(10) << std::setprecision(6) << aux[i * M + j] << " ";
                std::cout << std::endl;
            }
        }

        inline size_t get_coef_aux_list_idx(const std::array<size_t, N> &dims) const {
            auto coef_aux_index = 0;
            for (auto &dim: dims) {
                coef_aux_index = coef_aux_index * COEF_AUX_MAX_BLOCK[N] + dim;
            }
            return coef_aux_index;
        }
    };

}

#endif
