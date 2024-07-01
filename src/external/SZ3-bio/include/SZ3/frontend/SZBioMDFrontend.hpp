#ifndef SZ3_SZBIOMD_FRONTEND
#define SZ3_SZBIOMD_FRONTEND

/**
 */

#include "Frontend.hpp"
//#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include <list>

namespace SZ3 {
    
    template<class T, uint N, class Quantizer>
    class SZBioMDFrontend : public concepts::FrontendInterface<T, N> {
     public:
        SZBioMDFrontend(const Config &conf, Quantizer quantizer) :
            quantizer(quantizer),
            conf(conf) {
            if (N != 1 && N != 2 && N != 3) {
                throw std::invalid_argument("SZBioFront only support 1D, 2D or 3D data");
            }
        }
        
        ~SZBioMDFrontend() {
            clear();
        }
        
        void print() {};
        
        std::vector<int> compress(T *data) {
            if (N == 1) {
                return compress_1d(data);
            } else if (N == 2) {
                return compress_2d(data);
            } else {
                return compress_3d(data);
            }
        };
        
        T *decompress(std::vector<int> &quant_inds, T *dec_data) {
            if (N == 1) {
                return decompress_1d(quant_inds, dec_data);
            } else if (N == 2) {
                return decompress_2d(quant_inds, dec_data);
            } else {
                return decompress_3d(quant_inds, dec_data);
            }
        };
        
        void save(uchar *&c) {
            write(site, c);
            write(firstFillFrame_, c);
            write(fillValue_, c);
            quantizer.save(c);
        }
        
        void load(const uchar *&c, size_t &remaining_length) {
            clear();
            const uchar *c_pos = c;
            read(site, c, remaining_length);
            read(firstFillFrame_, c, remaining_length);
            read(fillValue_, c, remaining_length);
            quantizer.load(c, remaining_length);
            remaining_length -= c_pos - c;
        }
        
        void clear() {
            quantizer.clear();
        }
        
        size_t size_est() {
            return quantizer.size_est(); //unpred
        }
        
        int get_radius() const {
            return quantizer.get_radius();
        }
        
        size_t get_num_elements() const {
            return conf.num;
        };
     
     private:
        std::vector<int> compress_1d(T *data) {
            std::vector<int> quant_bins(conf.num);
            quant_bins[0] = quantizer.quantize_and_overwrite(data[0], 0);
            for (size_t i = 1; i < conf.num; i++) {
                quant_bins[i] = quantizer.quantize_and_overwrite(data[i], data[i - 1]);
            }
            return quant_bins;
        }
        
        T *decompress_1d(std::vector<int> &quant_inds, T *dec_data) {
            dec_data[0] = quantizer.recover(0, quant_inds[0]);
            for (size_t i = 1; i < conf.num; i++) {
                dec_data[i] = quantizer.recover(dec_data[i - 1], quant_inds[i]);
            }
            return dec_data;
        }
        
        int cal_site(T *data, std::vector<size_t> dims) {
            size_t numDims = dims.size();
            if (numDims < 2) {
                return 0;
            }
            std::vector<int> sites;
            for (int j = 0; j < std::min<size_t>(dims[numDims - 1], 5); j++) {
                size_t lprev = 0, lavg = 0, lcnt = 0;
                for (size_t i = 1; i < std::min<size_t>(dims[numDims - 2], 100); i++) {
                    auto c = data[i * dims[numDims - 1] + j], p = data[(i - 1) * dims[numDims - 1] + j];
                    if (fabs(c - p) / c > 0.5) {
                        sites.push_back(i - lprev);
//                        printf("%d %d\n", i, i - lprev);
                        lprev = i;
                    }
                }
            }
            ska::unordered_map<int, size_t> frequency;
            for (size_t i = 0; i < sites.size(); i++) {
                frequency[sites[i]]++;
            }
            int maxCount = 0, res = 0;
            for (const auto &kv : frequency) {
                auto k = kv.first;
                auto f = kv.second;
//                printf("k %d f %d\n", k ,f);
                if (maxCount < f) {
                    res = k;
                    maxCount = f;
                }
            }
            return (res <= 2 || res > 10) ? 0 : res;
        }
        
        /* Start from the last frame and look for frames filled with the same value.
         * Those frames do not need to be compressed - they can all just be filled. */
        std::tuple<size_t, T> findFillValueAndFirstFilledFrame(T *data, std::vector<size_t> dims) {
            size_t numDims = dims.size();
            if (numDims < 3) {
                return std::make_tuple(dims[0], 0);
            }
            size_t frameStride = dims[1] * dims[2];
            size_t firstFillFrame = dims[0];
            
            if (firstFillFrame == 0) {
                return std::make_tuple(firstFillFrame, 0);
            }
            /* Assume that the first value of the last frame is the potential fill value */
            T fillFrameValue = data[(dims[0] - 1) * frameStride];
            /* To simplify compression/decompression below, assume that the first frame needs
             * to be compressed. */
            for (size_t i = dims[0] - 1; i > 0; i--) {
                size_t idx = i * frameStride;
                bool allFrameValuesAreFillValue = true;
                
                for (size_t j = 0; j < dims[1] * dims[2]; j++) {
                    size_t idy = idx + j;
                    if (data[idy] != fillFrameValue) {
                        allFrameValuesAreFillValue = false;
                        break;
                    }
                }
                if (allFrameValuesAreFillValue) {
                    firstFillFrame = i;
                } else {
                    break;
                }
            }
            return std::make_tuple(firstFillFrame, fillFrameValue);
        }
        
        std::vector<int> compress_2d(T *data) {
            std::vector<int> quant_bins(conf.num);
            auto dims = conf.dims;
            site = cal_site(data, conf.dims);
            printf("# of site in the MD simulation guessed by SZ3 = %d\n", site);
            
            firstFillFrame_ = dims[0];
            fillValue_ = 0;
            
            //i==0 & j==0
            for (size_t k = 0; k < dims[1]; k++) { //xyz
                size_t idx = k;
                quant_bins[idx] = quantizer.quantize_and_overwrite(data[idx], 0);
            }
            
            //i==0
            for (size_t j = 1; j < dims[0]; j++) { //atoms
                /* Assume that we are compressing water coordinates. Water models are constructed the following way:
                 * 3-site: O-H-H
                 * 4-site: O-H-H-VS
                 * 5-site: O-H-H-VS-VS
                 * The oxygen atom is always most centrally located (on average closest to the others).
                 * Use the oxygen atom as reference.
                 */
                size_t site_reference_offset = site != 0 ? std::max<size_t>(1, j % site) : 1;
                for (size_t k = 0; k < dims[1]; k++) { //xyz
                    size_t idx = j * dims[1] + k;
                    size_t idx1 = (j - site_reference_offset) * dims[1] + k;
                    quant_bins[idx] = quantizer.quantize_and_overwrite(data[idx], data[idx1]);
                }
            }
            
            return quant_bins;
        }
        
        T *decompress_2d(std::vector<int> &quant_inds, T *dec_data) {
            
            auto dims = conf.dims;
            
            //i==0 & j==0
            for (size_t k = 0; k < dims[1]; k++) { //xyz
                size_t idx = k;
                dec_data[idx] = quantizer.recover(0, quant_inds[idx]);
            }
            
            //i==0
            for (size_t j = 1; j < dims[0]; j++) { //atoms
                /* Assume that we are decompressing water coordinates. Water models are constructed the following way:
                 * 3-site: O-H-H
                 * 4-site: O-H-H-VS
                 * 5-site: O-H-H-VS-VS
                 * The oxygen atom is always most centrally located (on average closest to the others).
                 * Use the oxygen atom as reference.
                 */
                size_t site_reference_offset = site != 0 ? std::max<size_t>(1, j % site) : 1;
                for (size_t k = 0; k < dims[1]; k++) { //xyz
                    size_t idx = j * dims[1] + k;
                    size_t idx1 = (j - site_reference_offset) * dims[1] + k;
                    dec_data[idx] = quantizer.recover(dec_data[idx1], quant_inds[idx]);
                }
            }
            
            return dec_data;
        }
        
        std::vector<int> compress_3d(T *data) {
            std::vector<int> quant_bins(conf.num);
            auto dims = conf.dims;
            std::vector<size_t> stride({dims[1] * dims[2], dims[2], 1});
            site = cal_site(data + stride[0], conf.dims);
            printf("# of site in the MD simulation guessed by SZ3 = %d\n", site);
            
            /* Find out if the last frames are all filled with the same value. */
            std::tuple<size_t, T> fillValueSettings = findFillValueAndFirstFilledFrame(data, dims);
            firstFillFrame_ = std::get<0>(fillValueSettings);
            fillValue_ = std::get<1>(fillValueSettings);
            size_t lastFrame = std::min(dims[0], firstFillFrame_);
            
            //TODO determine the # of system
            //i==0 & j==0
            for (size_t k = 0; k < dims[2]; k++) { //xyz
                size_t idx = k;
                quant_bins[idx] = quantizer.quantize_and_overwrite(data[idx], 0);
            }
            
            //i==0
            for (size_t j = 1; j < dims[1]; j++) { //atoms
                /* Assume that we are compressing water coordinates. Water models are constructed the following way:
                 * 3-site: O-H-H
                 * 4-site: O-H-H-VS
                 * 5-site: O-H-H-VS-VS
                 * The oxygen atom is always most centrally located (on average closest to the others).
                 * Use the oxygen atom as reference.
                 */
                size_t site_reference_offset = site != 0 ? std::max<size_t>(1, j % site) : 1;
                for (size_t k = 0; k < dims[2]; k++) { //xyz
                    size_t idx = j * stride[1] + k;
                    size_t idx1 = (j - site_reference_offset) * stride[1] + k;
                    quant_bins[idx] = quantizer.quantize_and_overwrite(data[idx], data[idx1]);
                }
            }
            
            for (size_t i = 1; i < lastFrame; i++) {//time
                for (size_t j = 0; j < dims[1]; j++) { //atoms
                    size_t site_reference_offset = site != 0 ? j % site : 1;
                    for (size_t k = 0; k < dims[2]; k++) { //xyz
                        size_t idx = i * stride[0] + j * stride[1] + k;
                        size_t idx1 = (i - 1) * stride[0] + j * stride[1] + k;
                        if (j == 0 || (site != 0 && j % site == 0)) {// time -1
                            quant_bins[idx] =
                                quantizer.quantize_and_overwrite(data[idx], data[idx1]);
                        } else { // time -1 & atom -1
                            size_t idx2 = i * stride[0] + (j - site_reference_offset) * stride[1] + k;
                            size_t idx3 = (i - 1) * stride[0] + (j - site_reference_offset) * stride[1] + k;
                            quant_bins[idx] =
                                quantizer.quantize_and_overwrite(data[idx], data[idx1] + data[idx2] - data[idx3]);
                        }
                    }
                }
            }
            
            return quant_bins;
        }
        
        T *decompress_3d(std::vector<int> &quant_inds, T *dec_data) {
            
            auto dims = conf.dims;
            std::vector<size_t> stride({dims[1] * dims[2], dims[2], 1});
//            quant_bins[0] = quantizer.quantize_and_overwrite(data[0], 0);
            
            size_t lastFrame = std::min(dims[0], firstFillFrame_);
            
            //i==0 & j==0
            for (size_t k = 0; k < dims[2]; k++) { //xyz
                size_t idx = k;
                dec_data[idx] = quantizer.recover(0, quant_inds[idx]);
            }
            
            //i==0
            for (size_t j = 1; j < dims[1]; j++) { //atoms
                /* Assume that we are decompressing water coordinates. Water models are constructed the following way:
                 * 3-site: O-H-H
                 * 4-site: O-H-H-VS
                 * 5-site: O-H-H-VS-VS
                 * The oxygen atom is always most centrally located (on average closest to the others).
                 * Use the oxygen atom as reference.
                 */
                size_t site_reference_offset = site != 0 ? std::max<size_t>(1, j % site) : 1;
                for (size_t k = 0; k < dims[2]; k++) { //xyz
                    size_t idx = j * stride[1] + k;
                    size_t idx1 = (j - site_reference_offset) * stride[1] + k;
                    dec_data[idx] = quantizer.recover(dec_data[idx1], quant_inds[idx]);
                }
            }
            
            for (size_t i = 1; i < lastFrame; i++) {//time
                for (size_t j = 0; j < dims[1]; j++) { //atoms
                    size_t site_reference_offset = site != 0 ? j % site : 1;
                    for (size_t k = 0; k < dims[2]; k++) { //xyz
                        size_t idx = i * stride[0] + j * stride[1] + k;
                        size_t idx1 = (i - 1) * stride[0] + j * stride[1] + k;
                        if (j == 0 || (site != 0 && j % site == 0)) {// time -1
                            dec_data[idx] = quantizer.recover(dec_data[idx1], quant_inds[idx]);
                        } else { // time -1 & atom -1
                            size_t idx2 = i * stride[0] + (j - site_reference_offset) * stride[1] + k;
                            size_t idx3 = (i - 1) * stride[0] + (j - site_reference_offset) * stride[1] + k;
                            dec_data[idx] = quantizer.recover(dec_data[idx1] + dec_data[idx2] - dec_data[idx3], quant_inds[idx]);
                        }
                    }
                }
            }
            
            /* Fill frames at the end with the fill value. */
            for (size_t i = firstFillFrame_; i < dims[0]; i++) {
                size_t idx = i * stride[0];
                for (size_t j = 0; j < dims[1] * dims[2]; j++) {
                    size_t idy = idx + j;
                    dec_data[idy] = fillValue_;
                }
            }
            return dec_data;
        }
        
        Quantizer quantizer;
        Config conf;
        int site = 0;
        size_t firstFillFrame_;
        T fillValue_;
        
    };
    
    template<class T, uint N, class Predictor>
    SZBioMDFrontend<T, N, Predictor>
    make_sz_bio_frontend(const Config &conf, Predictor predictor) {
        return SZBioMDFrontend<T, N, Predictor>(conf, predictor);
    }
}

#endif
