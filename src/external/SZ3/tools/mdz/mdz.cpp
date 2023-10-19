#include <vector>

#include <mdz.hpp>


int main(int argc, char **argv) {


    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 2);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    SZ::Config conf({1, dims[0]});
    if (dim == 2) {
        conf = SZ::Config({dims[0], dims[1]});
    }
    std::string input_path = argv[1];

    char *eb_op = argv[argp++] + 1;
    if (*eb_op == 'a') {
        conf.errorBoundMode = SZ::EB_ABS;
        conf.absErrorBound = atof(argv[argp++]);
    } else {
        conf.errorBoundMode = SZ::EB_REL;
        conf.relErrorBound = atof(argv[argp++]);
    }
    size_t batch_size = 0;
    if (argp < argc) {
        batch_size = atoi(argv[argp++]);
    }
    int method = -1;
    if (argp < argc) {
        method = atoi(argv[argp++]);
    }

    conf.blockSize = 128;
    conf.stride = 128;
    conf.quantbinCnt = 1024;
//    conf.enable_regression = false;
//    conf.quant_state_num = 4096;
    if (argp < argc) {
        conf.quantbinCnt = atoi(argv[argp++]);
    }
    auto input_data = readfile<float>(input_path.data(), 0, conf.num);
    std::vector<float> input2(input_data.get(), input_data.get() + conf.num);
    std::vector<float> dec_data(conf.num);

    size_t compressed_size = MDZ_Compress<float, 2>(conf, input_data.get(), dec_data.data(), batch_size, method);
    float ratio = conf.num * 1.0 * sizeof(float) / compressed_size;

//    std::stringstream ss;
//    ss << input_path.substr(input_path.rfind('/') + 1)
//       << ".b" << batch_size
//       << "." << conf.relErrorBound << ".md-" << method << ".out";
//    std::cout << "Decompressed file = " << ss.str() << std::endl;
//    SZ::writefile(ss.str().data(), dec_data.data(), total_num);

    double max_diff, psnr, nrmse;
    SZ::verify<float>(input2.data(), dec_data.data(), conf.num, psnr, nrmse, max_diff);
    std::cout << "****************** Final ****************" << std::endl;
    printf("method=md, file=%s, block=%lu, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f, timestep_op=%d\n",
           input_path.data(), batch_size,
           ratio,
           conf.relErrorBound,
           max_diff, psnr, nrmse,
           total_compress_time, total_decompress_time,
           method);
}