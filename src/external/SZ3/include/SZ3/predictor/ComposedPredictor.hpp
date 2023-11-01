#ifndef _SZ_COMPOSED_PREDICTOR_HPP
#define _SZ_COMPOSED_PREDICTOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include <cassert>
#include <iostream>
#include <memory>

namespace SZ3 {

    template<class T, uint N>
    class ComposedPredictor : public concepts::PredictorInterface<T, N> {
    public:
        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        void precompress_data(const iterator &iter) const noexcept {
            for (const auto &p:predictors) {
                p->precompress_data(iter);
            }
        }

        void postcompress_data(const iterator &iter) const noexcept {
            for (const auto &p:predictors) {
                p->postcompress_data(iter);
            }
        }

        void predecompress_data(const iterator &iter) const noexcept {
            for (const auto &p:predictors) {
                p->predecompress_data(iter);
            }
        }

        void postdecompress_data(const iterator &iter) const noexcept {
            for (const auto &p:predictors) {
                p->postdecompress_data(iter);
            }
        }


        bool precompress_block(const std::shared_ptr<Range> &range) {
            std::vector<bool> precompress_block_result;
            for (const auto &p:predictors) {
                precompress_block_result.push_back(p->precompress_block(range));
            }
            const auto &dims = range->get_dimensions();
            int min_dimension = *std::min_element(dims.begin(), dims.end());

            do_estimate_error(range->begin(), min_dimension);

            sid = std::distance(predict_error.begin(), std::min_element(predict_error.begin(), predict_error.end()));
            return precompress_block_result[sid];
        }

        void precompress_block_commit() {
            selection.push_back(sid);
            predictors[sid]->precompress_block_commit();
        }

        bool predecompress_block(const std::shared_ptr<Range> &range) {
            sid = selection[current_index++];
            return predictors[sid]->predecompress_block(range);
        }

        void save(uchar *&c) const {
            auto tmp = c;
            for (const auto &p:predictors) {
                p->save(c);
            }
            // store selection

            *reinterpret_cast<size_t *>(c) = (size_t) selection.size();
            c += sizeof(size_t);
            if (selection.size()) {
                HuffmanEncoder<int> selection_encoder;
                selection_encoder.preprocess_encode(selection, predict_error.size());
                selection_encoder.save(c);
                selection_encoder.encode(selection, c);
                selection_encoder.postprocess_encode();
            }
//            *reinterpret_cast<size_t *>(c) = (size_t) selection.size();
//            c += sizeof(size_t);
//            memcpy(c, selection.data(), selection.size() * sizeof(int));
//            c += selection.size() * sizeof(int);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            auto tmp = c;
            for (const auto &p:predictors) {
                p->load(c, remaining_length);
            }

            // load selection
            // TODO: check correctness
            size_t selection_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            if (selection_size > 0) {
                remaining_length -= sizeof(size_t);
                HuffmanEncoder<int> selection_encoder;
                selection_encoder.load(c, remaining_length);
                this->selection = selection_encoder.decode(c, selection_size);
                selection_encoder.postprocess_decode();
            }
//            size_t selection_size = *reinterpret_cast<const size_t *>(c);
//            c += sizeof(size_t);
//            this->selection = std::vector<int>(reinterpret_cast<const int *>(c),
//                                               reinterpret_cast<const int *>(c) + selection_size);
//            c += selection_size * sizeof(int);
        }

        inline T predict(const iterator &iter) const noexcept {
            return predictors[sid]->predict(iter);
        }

        int get_sid() const { return sid; }

        void set_sid(int _sid) {
            sid = _sid;
        }

        T estimate_error(const iterator &iter) const noexcept {
            return predictors[sid]->estimate_error(iter);
        }

        void print() const {
            std::vector<size_t> cnt(predictors.size(), 0);
            size_t cnt_total = 0;
            for (auto &sel:selection) {
                cnt[sel]++;
                cnt_total++;
            }
            for (int i = 0; i < predictors.size(); i++) {
//                predictors[i]->print();
                printf("Blocks:%ld, Percentage:%.2f\n", cnt[i], 1.0 * cnt[i] / cnt_total);
            }
        }

//        template<typename P1>
//        void instantiate(P1 p1) {
//            predictors.push_back(std::move(p1));
//        }
//
//        template<typename P1>
//        void unpack(P1 p1) {
//            instantiate(p1);
//        }
//
//        template<typename P1, typename... Rest>
//        void unpack(P1 p1, Rest... Rs) {
//            instantiate<P1>(p1);
//            unpack(Rs...);
//        }
//
//        template<class... Predictors>
//        ComposedPredictor(Predictors &&... Ps) {
//            unpack(Ps...);
//        }

        ComposedPredictor(std::vector<std::shared_ptr<concepts::PredictorInterface < T, N>>

        > predictors) {
            this->predictors = predictors;
            predict_error.resize(predictors.size());
        }

        void clear() {
            for (auto &pred:predictors) {
                pred->clear();
            }
            selection.clear();
        }

    private:
        std::vector<std::shared_ptr<concepts::PredictorInterface < T, N>>>
        predictors;
        std::vector<int> selection;
        int sid = 0;                            // selected index
        size_t current_index = 0;            // for decompression only
        std::vector<double> predict_error;

        template<uint NN = N>
        inline typename std::enable_if<NN == 1, void>::type
        do_estimate_error(const iterator &iter, int min_dimension) {
            std::fill(predict_error.begin(), predict_error.end(), 0);
            auto iter1 = iter;
            iter1.move(min_dimension - 1);
            for (int p = 0; p < predictors.size(); p++) {
                predict_error[p] += predictors[p]->estimate_error(iter);
                predict_error[p] += predictors[p]->estimate_error(iter1);
            }
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 2, void>::type
        do_estimate_error(const iterator &iter, int min_dimension) {
            std::fill(predict_error.begin(), predict_error.end(), 0);
            auto iter1 = iter, iter2 = iter;
            iter2.move(0, min_dimension - 1);
            for (int i = 2; i < min_dimension; i++) {
                for (int p = 0; p < predictors.size(); p++) {
                    predict_error[p] += predictors[p]->estimate_error(iter1);
                    predict_error[p] += predictors[p]->estimate_error(iter2);
                }
                iter1.move(1, 1);
                iter2.move(1, -1);
            }
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 3, void>::type
        do_estimate_error(const iterator &iter, int min_dimension) {
            std::fill(predict_error.begin(), predict_error.end(), 0);
//            std::vector<double> err(predictors.size(), 0);
            auto iter1 = iter, iter2 = iter, iter3 = iter, iter4 = iter;
            iter2.move(0, 0, min_dimension - 1);
            iter3.move(0, min_dimension - 1, 0);
            iter4.move(0, min_dimension - 1, min_dimension - 1);
            for (int i = 2; i < min_dimension; i++) {
                for (int p = 0; p < predictors.size(); p++) {
                    predict_error[p] += predictors[p]->estimate_error(iter1);
                    predict_error[p] += predictors[p]->estimate_error(iter2);
                    predict_error[p] += predictors[p]->estimate_error(iter3);
                    predict_error[p] += predictors[p]->estimate_error(iter4);
                }
                iter1.move(1, 1, 1);
                iter2.move(1, 1, -1);
                iter3.move(1, -1, 1);
                iter4.move(1, -1, -1);
            }
        }

        template<uint NN = N>
        inline typename std::enable_if<NN >= 4, void>::type
        do_estimate_error(const iterator &iter, int min_dimension) {
            std::fill(predict_error.begin(), predict_error.end(), 0);
//            std::vector<double> err(predictors.size(), 0);
            auto iter1 = iter, iter2 = iter, iter3 = iter, iter4 = iter,
                    iter5 = iter, iter6 = iter, iter7 = iter, iter8 = iter;;
            iter2.move(0, 0, 0, min_dimension - 1);
            iter3.move(0, 0, min_dimension - 1, 0);
            iter4.move(0, 0, min_dimension - 1, min_dimension - 1);
            iter5.move(0, min_dimension - 1, 0, 0);
            iter6.move(0, min_dimension - 1, 0, min_dimension - 1);
            iter7.move(0, min_dimension - 1, min_dimension - 1, 0);
            iter8.move(0, min_dimension - 1, min_dimension - 1, min_dimension - 1);
            for (int i = 2; i < min_dimension; i++) {
                for (int p = 0; p < predictors.size(); p++) {
                    predict_error[p] += predictors[p]->estimate_error(iter1);
                    predict_error[p] += predictors[p]->estimate_error(iter2);
                    predict_error[p] += predictors[p]->estimate_error(iter3);
                    predict_error[p] += predictors[p]->estimate_error(iter4);
                    predict_error[p] += predictors[p]->estimate_error(iter5);
                    predict_error[p] += predictors[p]->estimate_error(iter6);
                    predict_error[p] += predictors[p]->estimate_error(iter7);
                    predict_error[p] += predictors[p]->estimate_error(iter8);
                }
                iter1.move(1, 1, 1, 1);
                iter2.move(1, 1, 1, -1);
                iter3.move(1, 1, -1, 1);
                iter4.move(1, 1, -1, -1);
                iter5.move(1, -1, 1, 1);
                iter6.move(1, -1, 1, -1);
                iter7.move(1, -1, -1, 1);
                iter8.move(1, -1, -1, -1);
            }
        }
    };

}


#endif
