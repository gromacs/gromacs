#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZBioMDFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/PolyRegressionPredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/def.hpp"
#include <cmath>
#include <memory>

namespace SZ3 {


    template<class T, uint N>
    char *SZ_compress_bioMD(Config &conf, T *data, size_t &outSize) {

        assert(N == conf.N);
        assert(conf.cmprAlgo == ALGO_BIOMD);
        calAbsErrorBound(conf, data);

        char *cmpData;
        auto quantizer = LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
        auto sz = make_sz_general_compressor<T, N>(make_sz_bio_frontend<T, N>(conf, quantizer), HuffmanEncoder<int>(),
                                                   Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
        return cmpData;
    }


    template<class T, uint N>
    void SZ_decompress_bioMD(const Config &conf, char *cmpData, size_t cmpSize, T *decData) {
        assert(conf.cmprAlgo == ALGO_BIOMD);

        uchar const *cmpDataPos = (uchar *) cmpData;
        LinearQuantizer<T> quantizer;
        auto sz = make_sz_general_compressor<T, N>(make_sz_bio_frontend<T, N>(conf, quantizer),
                                                   HuffmanEncoder<int>(), Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);

    }
}
#endif