/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\utils\msl_assert.hpp
** -----
** File Created: Thursday, 9th October 2025 22:33:56
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 9th October 2025 22:34:02
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

// Description: Custom assert handler for MSL

#ifndef MSL_ASSERT_HPP
#define MSL_ASSERT_HPP

#include <fstream>
#include <iostream>

namespace msl {
inline void handle_assert(const char *expr,
                          const char *msg,
                          const char *file,
                          int line,
                          const char *func) {
    std::string text = std::string("[MSL Assert] ") + msg + "\n"
                       + "  Expr: " + expr + "\n" + "  Location: " + file + ":"
                       + std::to_string(line) + " (" + func + ")\n";

    // --- 打印到标准错误流 ---
    std::cerr << text;

    // --- 同时写入日志文件 ---
    // TODO: 可以改成更灵活的日志系统
    std::ofstream log("msl_error.log", std::ios::app);
    if (log.is_open()) {
        log << text;
    }

    std::abort();
}
} // namespace msl

#endif // MSL_ASSERT_HPP