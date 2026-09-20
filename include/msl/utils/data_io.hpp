/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\utils\data_io.hpp
** -----
** File Created: Thursday, 23rd October 2025 15:31:49
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 23rd October 2025 15:31:56
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

// Utility functions for data I/O

#ifndef MSL_UTILS_DATA_IO_HPP
#define MSL_UTILS_DATA_IO_HPP

#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

namespace msl::utils {
inline std::vector<double>
ReadData(const std::string &file_name, int col, int row) {
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: file open failed." << std::endl;
        return {};
    }

    // 读取矩阵数据(按列)
    std::vector<double> data(col * row);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            file >> data[i + j * row];
        }
    }

    return data;
}

inline void WriteData(const std::string &file_name,
                      const std::vector<double> &data,
                      int col,
                      int row) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: file open failed." << std::endl;
        return;
    }

    // 写入矩阵数据(按列)
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            file << data[i + j * row] << " ";
        }
        file << std::endl;
    }
}

inline void WriteComplexData(const std::string &file_name,
                             std::vector<std::complex<double>> &data,
                             int col,
                             int row) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: file open failed." << std::endl;
        return;
    }

    // 写入矩阵数据(按列)
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            file << data[i + j * row] << " ";
        }
        file << std::endl;
    }
}

inline std::string ReadJson(const std::string &file_name) {
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: file open failed." << std::endl;
        return "";
    }

    std::string json_str((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    return json_str;
}
} // namespace msl::utils

#endif // MSL_UTILS_DATA_IO_HPP