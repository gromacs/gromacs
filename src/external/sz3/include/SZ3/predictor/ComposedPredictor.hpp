#ifndef SZ3_COMPOSED_PREDICTOR_HPP
#define SZ3_COMPOSED_PREDICTOR_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/predictor/Predictor.hpp"

namespace SZ3 {

template <class T, uint N, class Encoder = HuffmanEncoder<int>>
class ComposedPredictor : public concepts::PredictorInterface<T, N> {
   public:
    using block_iter = typename block_data<T, N>::block_iterator;

    ComposedPredictor(std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors_)
        : predictors(predictors_) {
        if (predictors.empty()) {
            throw std::invalid_argument("Empty predictor list for ComposedPredictor.");
        }
    }

    bool precompress(const block_iter &block) override {
        std::vector<double> predict_error(predictors.size(), 0);
        std::vector<bool> isvalid(predictors.size());
        for (size_t i = 0; i < predictors.size(); i++) {
            isvalid[i] = predictors[i]->precompress(block);
            if (isvalid[i]) {
                block_iter::foreach_sampling(block, [&](T *c, const std::array<size_t, N> &index) {
                    predict_error[i] += predictors[i]->estimate_error(block, c, index);
                });
            } else {
                predict_error[i] = std::numeric_limits<double>::max();
            }
        }
        sid = std::distance(predict_error.begin(), std::min_element(predict_error.begin(), predict_error.end()));
        return isvalid[sid];
    }

    void precompress_block_commit() override {
        selection.push_back(sid);
        predictors[sid]->precompress_block_commit();
    }

    bool predecompress(const block_iter &block) override {
        // selection and its entries are untrusted, and predict()/estimate_error() reuse the sid set here.
        if (current_index >= selection.size())
            throw std::out_of_range("SZ3: ran out of predictor selections while decompressing");
        sid = selection[current_index++];
        if (sid < 0 || static_cast<size_t>(sid) >= predictors.size())
            throw std::out_of_range("SZ3: predictor selection index is out of range");
        return predictors[sid]->predecompress(block);
    }

    size_t size_est() override {
        size_t total = sizeof(uint64_t);
        for (const auto &p : predictors) {
            total += p->size_est();
        }
        if (!selection.empty()) {
            total += Encoder::size_bound(selection.size(), predictors.size());
        }
        return total;
    }

    void save(uchar *&c) override {
        for (const auto &p : predictors) {
            p->save(c);
        }
        write<uint64_t>(selection.size(), c);
        if (selection.size() > 0) {
            Encoder selection_encoder;
            selection_encoder.preprocess_encode(selection, predictors.size());
            selection_encoder.save(c);
            selection_encoder.encode(selection, c);
            selection_encoder.postprocess_encode();
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        for (const auto &p : predictors) {
            p->load(c, remaining_length);
        }
        uint64_t selection_size = 0;
        read(selection_size, c, remaining_length);
        if (selection_size > 0) {
            Encoder selection_encoder;
            selection_encoder.load(c, remaining_length);
            this->selection = selection_encoder.decode(c, selection_size, remaining_length);
            selection_encoder.postprocess_decode();
        }
    }

    T predict(const block_iter &block, T *d, const std::array<size_t, N> &index) override {
        return predictors[sid]->predict(block, d, index);
    }

    T estimate_error(const block_iter &block, T *d, const std::array<size_t, N> &index) override {
        return predictors[sid]->estimate_error(block, d, index);
    }

    void print() const override {
        std::vector<size_t> cnt(predictors.size(), 0);
        size_t cnt_total = 0;
        for (auto &sel : selection) {
            cnt[sel]++;
            cnt_total++;
        }
        for (size_t i = 0; i < predictors.size(); i++) {
            //                predictors[i]->print();
            printf("Blocks:%zu, Percentage:%.2f\n", cnt[i], 1.0 * cnt[i] / cnt_total);
        }
    }

    size_t get_padding() override {
        size_t max_padding = 0;
        for (const auto &p : predictors) {
            max_padding = std::max(max_padding, p->get_padding());
        }
        return max_padding;
    }

   private:
    std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors;
    std::vector<int> selection;
    int sid = 0;               // selected index
    size_t current_index = 0;  // for decompression only
};

}  // namespace SZ3

#endif
