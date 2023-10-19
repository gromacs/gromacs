//
// Created by Kai Zhao on 7/1/21.
//

#ifndef SZ_MDZ_H
#define SZ_MDZ_H

#include <SZ3/compressor/SZGeneralCompressor.hpp>
#include <SZ3/predictor/Predictor.hpp>
#include <SZ3/predictor/LorenzoPredictor.hpp>
#include <SZ3/predictor/RegressionPredictor.hpp>
#include <SZ3/predictor/ComposedPredictor.hpp>
#include <SZ3/utils/FileUtil.hpp>
#include <SZ3/utils/Timer.hpp>
#include <SZ3/utils/Statistic.hpp>
#include <SZ3/def.hpp>
#include <SZ3/quantizer/IntegerQuantizer.hpp>
#include <SZ3/lossless/Lossless_zstd.hpp>
#include <SZ3/encoder/HuffmanEncoder.hpp>
#include <SZ3/frontend/SZGeneralFrontend.hpp>
#include <SZ3/utils/QuantOptimizatioin.hpp>
#include "KmeansUtil.hpp"
#include "TimeBasedFrontend.hpp"
#include "ExaaltCompressor.hpp"

double total_compress_time = 0;
double total_decompress_time = 0;
const char *compressor_names[] = {"VQ", "VQT", "MT", "LR", "TS"};


template<typename T, uint N, class Predictor>
SZ::concepts::CompressorInterface<T> *
make_sz_timebased2(const SZ::Config &conf, Predictor predictor, T *data_ts0) {
    return new SZ::SZGeneralCompressor<T, N, SZ::TimeBasedFrontend<T, N, Predictor, SZ::LinearQuantizer<T>>,
            SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::TimeBasedFrontend<T, N, Predictor, SZ::LinearQuantizer<T>>(conf, predictor,
                                                                           SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
                                                                           data_ts0),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
}

template<typename T, uint N>
SZ::concepts::CompressorInterface<T> *
make_sz_timebased(const SZ::Config &conf, T *data_ts0) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N - 1>>> predictors;

    int use_single_predictor =
            (conf.lorenzo + conf.lorenzo2 + conf.regression) == 1;
    if (conf.lorenzo) {
        if (use_single_predictor) {
            return make_sz_timebased2<T, N>(conf, SZ::LorenzoPredictor<T, N - 1, 1>(conf.absErrorBound), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 1>>(conf.absErrorBound));
        }
    }
    if (conf.lorenzo2) {
        if (use_single_predictor) {
            return make_sz_timebased2<T, N>(conf, SZ::LorenzoPredictor<T, N - 1, 2>(conf.absErrorBound), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N - 1, 2>>(conf.absErrorBound));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return make_sz_timebased2<T, N>(conf, SZ::RegressionPredictor<T, N - 1>(conf.blockSize, conf.absErrorBound), data_ts0);
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N - 1>>(conf.blockSize, conf.absErrorBound));
        }
    }

    return make_sz_timebased2<T, N>(conf, SZ::ComposedPredictor<T, N - 1>(predictors), data_ts0);
}

template<typename T, uint N, class Predictor>
SZ::concepts::CompressorInterface<T> *
make_sz2(const SZ::Config &conf, Predictor predictor) {

    return new SZ::SZGeneralCompressor<T, N, SZ::SZGeneralFrontend<T, N, Predictor, SZ::LinearQuantizer<T>>,
            SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
            SZ::SZGeneralFrontend<T, N, Predictor, SZ::LinearQuantizer<T>>(conf, predictor,
                                                                           SZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd());
}

template<typename T, uint N>
SZ::concepts::CompressorInterface<T> *
make_sz(const SZ::Config &conf) {
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;

    int use_single_predictor =
            (conf.lorenzo + conf.lorenzo2 + conf.regression) == 1;
    if (conf.lorenzo) {
        if (use_single_predictor) {
            return make_sz2<T, N>(conf, SZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        }
    }
    if (conf.lorenzo2) {
        if (use_single_predictor) {
            return make_sz2<T, N>(conf, SZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound));
        } else {
            predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return make_sz2<T, N>(conf, SZ::RegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound));
        } else {
            predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N >>(conf.blockSize, conf.absErrorBound));
        }
    }

    return make_sz2<T, N>(conf, SZ::ComposedPredictor<T, N>(predictors));
}


template<typename T, uint N>
float *
VQ(SZ::Config conf, size_t ts, T *data, size_t &compressed_size, bool decom,
   int method, float level_start, float level_offset, int level_num) {
    if (level_num == 0) {
        printf("VQ/VQT not availble on current dataset, please use ADP or MT\n");
        exit(0);
    }

    auto sz = SZ::SZ_Exaalt_Compressor<float, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>,
            SZ::Lossless_zstd>(conf, SZ::LinearQuantizer<float>(conf.absErrorBound, conf.quantbinCnt / 2),
                               SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), method);
    sz.set_level(level_start, level_offset, level_num);

    SZ::Timer timer(true);
    SZ::uchar *compressed;
    compressed = sz.compress(data, compressed_size);
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        delete[]compressed;
        return nullptr;
    }
    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    timer.start();
    auto ts_dec_data = sz.decompress(compressed, compressed_size);
    total_decompress_time += timer.stop("Decompression");

    delete[]compressed;
    return ts_dec_data;
}


template<typename T, uint N>
float *
MT(SZ::Config conf, size_t ts, T *data, size_t &compressed_size, bool decom, T *ts0) {
//    printf("eb=%.8f\n", conf.eb);
    auto sz = make_sz_timebased<T, N>(conf, ts0);

    SZ::uchar *compressed;
    SZ::Timer timer(true);
    compressed = sz->compress(conf, data, compressed_size);
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        delete[] compressed;
        return nullptr;
    }

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;


    timer.start();
    auto ts_dec_data = sz->decompress(compressed, compressed_size, conf.num);
    total_decompress_time += timer.stop("Decompression");

    delete sz;
    delete[] compressed;
    return ts_dec_data;
}

template<typename T, uint N>
float *
SZ2(SZ::Config conf, size_t ts, T *data, size_t &compressed_size, bool decom) {
    auto sz = make_sz<T, N>(conf);
    SZ::uchar *compressed;

    SZ::Timer timer(true);
    compressed = sz->compress(conf, data, compressed_size);
    total_compress_time += timer.stop("Compression");
    if (!decom) {
        delete[] compressed;
        return nullptr;
    }

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    timer.start();
    auto ts_dec_data = sz->decompress(compressed, compressed_size, conf.num);
    total_decompress_time += timer.stop("Decompression");

    delete sz;
    delete[] compressed;
    return ts_dec_data;
}


template<typename T, uint N>
void select(SZ::Config conf, int &method, size_t ts, T *data_all,
            float level_start, float level_offset, int level_num, T *data_ts0, size_t timestep_batch) {
//        && (ts_last_select == -1 || t - ts_last_select >= conf.timestep_batch * 10)) {
    std::cout << "****************** BEGIN Selection ****************" << std::endl;
//        ts_last_select = ts;
    std::vector<size_t> compressed_size(10, std::numeric_limits<size_t>::max());
    std::vector<T> data1;
    size_t t = ts;
    if (ts == 0) {
        t = conf.dims[0] / 2;
        conf.dims[0] /= 2;
    }
    if (timestep_batch > 10) {
        conf.dims[0] = 10;
    }
    conf.num = conf.dims[0] * conf.dims[1];

    std::cout << conf.dims[0] << " " << conf.dims[1] << " " << t << std::endl;

    if (level_num > 0) {
        data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
        VQ<T, N>(conf, t, data1.data(), compressed_size[0], false, 0, level_start, level_offset, level_num);

        data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
        VQ<T, N>(conf, t, data1.data(), compressed_size[1], false, 1, level_start, level_offset, level_num);
    } else {
        data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
        SZ2<T, N>(conf, t, data1.data(), compressed_size[3], false);
    }

    data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
    MT<T, N>(conf, t, data1.data(), compressed_size[2], false, data_ts0);

//    data1 = std::vector(&data_all[t * conf.dims[1]], &data_all[t * conf.dims[1]] + conf.num);
//    MT(conf, t, data1.data(), compressed_size[4], false, (T *) nullptr);

    method = std::distance(compressed_size.begin(),
                           std::min_element(compressed_size.begin(), compressed_size.end()));
    printf("Select %s as Compressor, timestep=%lu, method=%d\n",
           compressor_names[method],
           ts, method);
    std::cout << "****************** END Selection ****************" << std::endl;
}

template<typename Type>
std::unique_ptr<Type[]> readfile(const char *file, size_t start, size_t num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        return 0;
    }
    fin.seekg(0, std::ios::end);
    const size_t total_num_elements = fin.tellg() / sizeof(Type);
    assert(start + num <= total_num_elements);
    fin.seekg(start * sizeof(Type), std::ios::beg);
    auto data = std::make_unique<Type[]>(num);
    fin.read(reinterpret_cast<char *>(&data[0]), num * sizeof(Type));
    fin.close();
    return data;
}


template<typename T, uint N>
SZ::uchar *LAMMPS_compress(SZ::Config conf, T *data, int method, size_t &compressed_size,
                           float level_start, float level_offset, int level_num, T *ts0) {
    if ((method == 0 || method == 1) && level_num == 0) {
        printf("VQ/VQT not available on current dataset, please use ADP or MT\n");
        exit(0);
    }
    SZ::uchar *compressed_data;
    if (method == 0 || method == 1) {
        auto sz = SZ::SZ_Exaalt_Compressor<T, N, SZ::LinearQuantizer<float>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd>(
                conf, SZ::LinearQuantizer<float>(conf.absErrorBound, conf.quantbinCnt / 2),
                SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd(), method);
        sz.set_level(level_start, level_offset, level_num);
        compressed_data = sz.compress(data, compressed_size);
    } else if (method == 2 || method == 4) {
        auto sz = make_sz_timebased<T, N>(conf, ts0);
        compressed_data = sz->compress(conf, data, compressed_size);
    } else {
        auto sz = make_sz<T, N>(conf);
        compressed_data = sz->compress(conf, data, compressed_size);
    }
//    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
//    std::cout << "Compression Ratio = " << ratio << std::endl;
//    std::cout << "Compressed size = " << compressed_size << std::endl;
    return compressed_data;
}

template<typename T, uint N>
int LAMMPS_select_compressor(SZ::Config conf, T *data, bool firsttime,
                             float level_start, float level_offset, int level_num, T *data_ts0) {
    std::cout << "****************** BEGIN Selection ****************" << std::endl;

    std::vector<size_t> compressed_size(10, std::numeric_limits<size_t>::max());

    if (firsttime) {
        conf.dims[0] /= 2;
        conf.num = conf.dims[0] * conf.dims[1];
        data += conf.num;
    }
    if (conf.dims[0] > 10) {
        conf.dims[0] = 10;
        conf.num = conf.dims[0] * conf.dims[1];
    }

    std::vector<T> data1;
    SZ::uchar *cmpr;
    if (level_num > 0) {

        data1 = std::vector<T>(data, data + conf.num);
        cmpr = LAMMPS_compress<T, N>(conf, data1.data(), 0, compressed_size[0], level_start, level_offset, level_num, data_ts0);
        delete[] cmpr;

        data1 = std::vector<T>(data, data + conf.num);
        cmpr = LAMMPS_compress<T, N>(conf, data1.data(), 1, compressed_size[1], level_start, level_offset, level_num, data_ts0);
        delete[] cmpr;
    } else {
        data1 = std::vector<T>(data, data + conf.num);
        cmpr = LAMMPS_compress<T, N>(conf, data1.data(), 3, compressed_size[3], level_start, level_offset, level_num, data_ts0);
        delete[] cmpr;
    }

    data1 = std::vector<T>(data, data + conf.num);
    cmpr = LAMMPS_compress<T, N>(conf, data1.data(), 2, compressed_size[2], level_start, level_offset, level_num, data_ts0);
    delete[] cmpr;

//    data1 = std::vector(&data[t * conf.dims[1]], &data[t * conf.dims[1]] + conf.num);
//    MT(conf, t, data1.data(), compressed_size[4], false, (T *) nullptr);

    int method = std::distance(compressed_size.begin(),
                               std::min_element(compressed_size.begin(), compressed_size.end()));
    printf("Select %s as Compressor, method=%d\n", compressor_names[method], method);
    std::cout << "****************** END Selection ****************" << std::endl;
    return method;
}

template<typename T, uint N>
size_t MDZ_Compress(SZ::Config conf, T *input_data, T *dec_data, size_t batch_size, int method = -1) {
    if (N != 2) {
        throw std::invalid_argument("dimension should be 2");
    }
    if (batch_size == 0) {
        batch_size = conf.dims[0];
    }
    int method_batch = 0;
    if (method == -1) {
        method_batch = 50;
    }
    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.absErrorBound
              << ", method = " << method
              << ", method_update_batch = " << method_batch
              << ", batch_size = " << batch_size
              << ", quan_state_num = " << conf.quantbinCnt
              //              << ", encoder = " << conf.encoder_op
              //              << ", lossless = " << conf.lossless_op
              << std::endl;

//    auto data_all = readfile<T>(input_path.data(), 0, conf.num);
//    auto data_all = input_data;
    auto data_ts0 = std::vector<T>(input_data, input_data + conf.dims[1]);

    float level_start, level_offset;
    int level_num = 0;
    if (method != 2 && method != 3 && method != 4) {
        size_t sample_num = 0.1 * conf.dims[1];
        sample_num = std::min(sample_num, (size_t) 20000);
        sample_num = std::max(sample_num, std::min((size_t) 5000, conf.dims[1]));
        SZ::get_cluster(input_data, conf.dims[1], level_start, level_offset, level_num,
                        sample_num);
        if (level_num > conf.dims[1] * 0.25) {
            level_num = 0;
        }
        if (level_num != 0) {
//            printf("start = %.3f , level_offset = %.3f, nlevel=%d\n", level_start, level_offset, level_num);
        }
    }

    auto dims = conf.dims;
    auto total_num = conf.num;
//    std::vector<T> dec_data(total_num);
    size_t total_compressed_size = 0;
    double compressed_size_pre = total_num * sizeof(T);
    int current_method = method;

    for (size_t ts = 0; ts < dims[0]; ts += batch_size) {
        conf.dims[0] = (ts + batch_size > dims[0] ? dims[0] - ts : batch_size);
        conf.num = conf.dims[0] * conf.dims[1];

//        auto data = SZ::readfile<T>(conf.input_path.data(), ts * conf.dims[1], conf.num);
        T *data = &input_data[ts * conf.dims[1]];

        T max = *std::max_element(data, data + conf.num);
        T min = *std::min_element(data, data + conf.num);
        if (conf.errorBoundMode == SZ::EB_ABS) {
            conf.relErrorBound = conf.absErrorBound / (max - min);
        } else if (conf.errorBoundMode == SZ::EB_REL) {
            conf.absErrorBound = conf.relErrorBound * (max - min);
        }

        std::cout << "****************** Compression From " << ts << " to " << ts + conf.dims[0] - 1
                  << " ******************" << std::endl;
//        std::cout<<method_batch<<" "<<ts<<" "<<conf.batch_size<<" "<<method_batch<<std::endl;
        if (method_batch > 0 && ts / batch_size % method_batch == 0) {
            select<T, N>(conf, current_method, ts, input_data, level_start, level_offset, level_num, data_ts0.data(), batch_size);
        }
        printf("Compressor = %s\n", compressor_names[current_method]);

        T *ts_dec_data;
        size_t compressed_size;
        if (current_method == 0) {
            ts_dec_data = VQ<T, N>(conf, ts, data, compressed_size, true, current_method, level_start, level_offset,
                                   level_num);
        } else if (current_method == 1) {
            ts_dec_data = VQ<T, N>(conf, ts, data, compressed_size, true, current_method, level_start, level_offset,
                                   level_num);
        } else if (current_method == 2) {
            ts_dec_data = MT<T, N>(conf, ts, data, compressed_size, true, data_ts0.data());
        } else if (current_method == 4) {
            ts_dec_data = MT<T, N>(conf, ts, data, compressed_size, true, (T *) nullptr);
        } else {
            ts_dec_data = SZ2<T, N>(conf, ts, data, compressed_size, true);
        }
        total_compressed_size += compressed_size;
//        if (compressed_size > 4.0 * compressed_size_pre) {
//            select(conf, current_method, ts, data_all.get(), level_start, level_offset, level_num, data_ts0.get());
//        }
        compressed_size_pre = compressed_size;
        memcpy(&dec_data[ts * conf.dims[1]], ts_dec_data, conf.num * sizeof(T));
    }


    return total_compressed_size;
}

#endif //SZ3_MDZ_H
