//
// Created by Kai Zhao on 11/10/22.
//

#include <mdz.hpp>


int main(int argc, char **argv) {


    std::vector<size_t> dims({100, 200});
    SZ3::Config conf({dims[0], dims[1]});

    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-6;
//    conf.blockSize = 128;
//    conf.stride = 128;
    conf.quantbinCnt = 1024;

    std::vector<float> input_data(conf.num);
    std::vector<float> dec_data(conf.num);

    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
            double x = static_cast<double>(i) - static_cast<double>(dims[0]) / 2.0;
            double y = static_cast<double>(j) - static_cast<double>(dims[1]) / 2.0;
            input_data[i * dims[1] + j] = static_cast<float>(.0001 * y * sin(y) + .0005 * cos(pow(x, 2) + x));
        }
    }

    size_t compressed_size = MDZ_Compress<float, 2>(conf, input_data.data(), dec_data.data(), 10);
    printf("compression ratio = %.5f", dims[0] * dims[1] * sizeof(float) * 1.0 / compressed_size);
    return 0;

}