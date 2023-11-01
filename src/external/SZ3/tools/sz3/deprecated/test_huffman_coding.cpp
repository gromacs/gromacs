#include "SZ3/encoder/HuffmanEncoder.hpp"
#include <vector>
#include <iostream>

using namespace std;

int main() {
    const int N = 65536;
    const int capacity = 1024;
    std::vector<int> type(N);
    for (int i = 0; i < N; i++) {
        type[i] = rand() % capacity;
    }

    unsigned char *compressed = (unsigned char *) malloc(N * sizeof(int));
    {
        SZ3::HuffmanEncoder<int> encoder;
        encoder.preprocess_encode(type, capacity);
        unsigned char *compressed_pos = compressed;
        cout << "save encoder" << endl;
        encoder.save(compressed_pos);
        cout << "tree size = ";
        cout << compressed_pos - compressed << endl;
        auto size = encoder.encode(type, compressed_pos);
        // auto size = encoder.encode_overall(type, compressed_pos);
        cout << N * sizeof(int) << " " << size << endl;
        encoder.postprocess_encode();
        // const unsigned char * compressed_pos_2 = compressed;
        // auto dec_type = encoder.decode(compressed_pos_2, N);
        // for(int i=0; i<N; i++){
        // 	if(type[i] != dec_type[i]){
        // 		cout << "decompressed type is not correct\n";
        // 		exit(0);
        // 	}
        // }
    }
    {
        SZ3::HuffmanEncoder<int> encoder;
        const unsigned char *compressed_pos = compressed;
        size_t length = sizeof(int);
        cout << "load" << endl;
        encoder.load(compressed_pos, length);
        cout << compressed_pos - compressed << endl;
        auto dec_type = encoder.decode(compressed_pos, N);
        for (int i = 0; i < N; i++) {
            if (type[i] != dec_type[i]) {
                cout << "2 decompressed type is not correct\n";
                exit(0);
            }
        }
        encoder.postprocess_decode();
    }
    free(compressed);

}
