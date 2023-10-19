//
// Created by Kai Zhao on 11/10/22.
//

#include <SZ3/api/sz.hpp>


int main(int argc, char **argv) {


    std::vector<size_t> dims({100, 200, 300});
    SZ::Config conf({dims[0], dims[1], dims[2]});
    conf.cmprAlgo = SZ::ALGO_INTERP_LORENZO;
    conf.errorBoundMode = SZ::EB_ABS; // refer to def.hpp for all supported error bound mode
    conf.absErrorBound = 1E-3; // absolute error bound 1e-3


    std::vector<float> input_data(conf.num);
    std::vector<float> dec_data(conf.num);

    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
            for (size_t k = 0; k < dims[2]; ++k) {
                double x = static_cast<double>(i) - static_cast<double>(dims[0]) / 2.0;
                double y = static_cast<double>(j) - static_cast<double>(dims[1]) / 2.0;
                double z = static_cast<double>(k) - static_cast<double>(dims[2]) / 2.0;
                input_data[i * dims[1] + j] = static_cast<float>(.0001 * y * sin(y) + .0005 * cos(pow(x, 2) + x) + z);
            }
        }
    }

    size_t cmpSize;
    char *cmpData = SZ_compress(conf, input_data.data(), cmpSize);
    auto dec_data_p = dec_data.data();
    SZ_decompress(conf, cmpData, cmpSize, dec_data_p);
    printf("%lu ", conf.num);
    return 0;


}