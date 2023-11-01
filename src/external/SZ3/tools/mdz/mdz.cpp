#include <vector>

#include <mdz.hpp>


void usage() {
    printf("Usage: \n");
    printf("For 1D input:   mdz file_path -1 n_atoms                 -r reb\n");
    printf("For 2D input:   mdz file_path -2 n_frames n_atoms        -r reb\n");
    printf("For 3D input:   mdz file_path -3 n_frames n_atoms n_dims(x,y,z) -r reb\n");
    exit(0);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage();
    }

    int dim = atoi(argv[2] + 1);
    if (dim > 3) {
        usage();
    }
    assert(1 <= dim && dim <= 2);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    SZ3::Config conf({1, dims[0]});
    if (dim == 2) {
        conf = SZ3::Config({dims[0], dims[1]});
    } else if (dim == 3) {
        conf = SZ3::Config({dims[0], dims[1], dims[2]});
    }

    std::string input_path = argv[1];

    if (argp >= argc) {
        usage();
    }
    char *eb_op = argv[argp++] + 1;
    if (*eb_op == 'a') {
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = atof(argv[argp++]);
    } else {
        conf.errorBoundMode = SZ3::EB_REL;
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

    size_t compressed_size;
    if (dim == 2) {
        compressed_size = MDZ_Compress<float, 2>(conf, input_data.get(), dec_data.data(), batch_size, method);
    } else if (dim == 3) {
        compressed_size = MDZ_Compress<float, 3>(conf, input_data.get(), dec_data.data(), batch_size, method);
    }
    float ratio = conf.num * 1.0 * sizeof(float) / compressed_size;

    double max_diff, psnr, nrmse;
    printf("\nBatch=%lu\nCompression ratio=%.3f\nCompression time=%.3f\nDecompression time=%.3f\n",
           (batch_size == 0 ? dims[0] : batch_size), ratio,
           total_compress_time, total_decompress_time);

    SZ3::verify<float>(input2.data(), dec_data.data(), conf.num, psnr, nrmse, max_diff);
//    std::cout << "****************** Final ****************" << std::endl;
//    printf("method=md, file=%s, block=%lu, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f, timestep_op=%d\n",
//           input_path.data(), batch_size,
//           ratio,
//           conf.relErrorBound,
//           max_diff, psnr, nrmse,
//           total_compress_time, total_decompress_time,
//           method);
}