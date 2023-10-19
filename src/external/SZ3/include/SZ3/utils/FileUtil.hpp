//
// Created by Kai Zhao on 12/9/19.
//

#ifndef _SZ_FILE_UTIL
#define _SZ_FILE_UTIL

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <random>
#include <sstream>
#include <memory>

namespace SZ {

    template<typename Type>
    void readfile(const char *file, const size_t num, Type *data) {
        std::ifstream fin(file, std::ios::binary);
        if (!fin) {
            std::cout << " Error, Couldn't find the file: " << file << "\n";
            exit(0);
        }
        fin.seekg(0, std::ios::end);
        const size_t num_elements = fin.tellg() / sizeof(Type);
        assert(num_elements == num && "File size is not equals to the input setting");
        fin.seekg(0, std::ios::beg);
        fin.read(reinterpret_cast<char *>(data), num_elements * sizeof(Type));
        fin.close();
    }

    template<typename Type>
    std::unique_ptr<Type[]> readfile(const char *file, size_t &num) {
        std::ifstream fin(file, std::ios::binary);
        if (!fin) {
            std::cout << " Error, Couldn't find the file: " << file << std::endl;
            exit(0);
        }
        fin.seekg(0, std::ios::end);
        const size_t num_elements = fin.tellg() / sizeof(Type);
        fin.seekg(0, std::ios::beg);
//        auto data = SZ::compat::make_unique<Type[]>(num_elements);
        auto data = std::make_unique<Type[]>(num_elements);
        fin.read(reinterpret_cast<char *>(&data[0]), num_elements * sizeof(Type));
        fin.close();
        num = num_elements;
        return data;
    }

    template<typename Type>
    void writefile(const char *file, Type *data, size_t num_elements) {
        std::ofstream fout(file, std::ios::binary);
        fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
        fout.close();
    }

    template<typename Type>
    void writeTextFile(const char *file, Type *data, size_t num_elements) {
        std::ofstream fout(file);
        if (fout.is_open()) {
            for (size_t i = 0; i < num_elements; i++) {
                fout << data[i] << std::endl;
            }
            fout.close();
        } else {
            std::cout << "Error, unable to open file for output: " << file << std::endl;
            exit(0);
        }
    }

}

#endif
