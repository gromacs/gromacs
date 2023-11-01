#include "SZ3/compressor/SZTruncateCompressor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <random>
#include <array>

std::string src_file_name;
float relative_error_bound = 0;
float compression_time = 0;
int byteLen = 2;

template<typename T, uint N, class Lossless>
float SZ_compress(std::unique_ptr<T[]> const &data,
                  const SZ3::Config &conf, Lossless lossless) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.absErrorBound
              << ", byteLen = " << byteLen
              << ", lossless = " << conf.lossless
              << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    auto sz = SZ3::make_sz_truncate_compressor<T, N>(conf, lossless, byteLen);

    SZ3::Timer timer;
    timer.start();
    std::cout << "****************** Compression ******************" << std::endl;


    size_t compressed_size = 0;
    std::unique_ptr<SZ3::uchar[]> compressed;
    compressed.reset(sz.compress(conf,data.get(), compressed_size));

    compression_time = timer.stop("Compression");

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 10000);
    std::stringstream ss;
    ss << src_file_name.substr(src_file_name.rfind('/') + 1)
       << "." << relative_error_bound << "." << dis(gen) << ".sz3";
    auto compressed_file_name = ss.str();
    SZ3::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ3::readfile<SZ3::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    T *dec_data = sz.decompress(compressed.get(), compressed_size, conf.num);
    timer.stop("Decompression");

    SZ3::verify<T>(data_.data(), dec_data, conf.num);

    remove(compressed_file_name.c_str());
//    auto decompressed_file_name = compressed_file_name + ".out";
//    SZ3::writefile(decompressed_file_name.c_str(), dec_data, conf.num);
//    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    delete[] dec_data;
    return ratio;
}

template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float eb,
                             std::array<size_t, N> dims) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.absErrorBound = eb;
    if (argp < argc) {
        conf.lossless = atoi(argv[argp++]);
    }
    if (conf.lossless > 0) {
        return SZ_compress<T, N>(data, conf, SZ3::Lossless_zstd(conf.lossless));
    } else {
        return SZ_compress<T, N>(data, conf, SZ3::Lossless_bypass());
    }
}


int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "usage: " << argv[0] << " data_file -num_dim dim0 .. dimn byteLens " << std::endl;
        std::cout << "example: " << argv[0] << " qmcpack.dat -3 33120 69 69 2" << std::endl;
        return 0;
    }

    size_t num = 0;
    auto data = SZ3::readfile<float>(argv[1], num);
    src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    if (argp < argc) {
        byteLen = atoi(argv[argp++]);
    }

    double eb = 0;


    if (dim == 1) {
        SZ_compress_parse_args<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_compress_parse_args<float, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_compress_parse_args<float, 3>(argc, argv, argp, data, eb,
                                         std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_compress_parse_args<float, 4>(argc, argv, argp, data, eb,
                                         std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }


    return 0;
}
