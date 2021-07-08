//
// Created by L.Nagy on 24/05/2021.
//

#ifndef SD_COOLING_DEBUG_HPP
#define SD_COOLING_DEBUG_HPP

#ifndef NDEBUG
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#endif

#ifndef NDEBUG

#ifndef DEBUG_NO_COLOUR
#define BOLD_RED  "\033[1;31m"
#define BOLD_CYAN "\033[1;36m"
#define NORMAL    "\033[0m"
#define DEBUG_FORMAT BOLD_CYAN << "DEBUG MESSAGE: " << NORMAL << __FILE__ << ":" << __LINE__ << " " << BOLD_RED       \
                               << __func__ << "(): " << NORMAL
#else
#define DEBUG_FORMAT "DEBUG MESSAGE: " << __FILE__ << ":" << __LINE__ << " " << __func__ << "(): "
#endif

#define DEBUG_INDENT "    "
#define DEBUG_NUMBER_FORMAT_WIDTH 30
#define DEBUG_NUMBER_FORMAT_PRECISION 15
#define DEBUG_NUMBER_FORMAT std::right << std::scientific << std::setw(DEBUG_NUMBER_FORMAT_WIDTH)                     \
                                       << std::setprecision(DEBUG_NUMBER_FORMAT_PRECISION)

#define DEBUG_MSG(x) do {                                                                                             \
                       std::cout << DEBUG_FORMAT << x << std::endl;                                                   \
                     } while (0)

#define DEBUG_MSG_SIMPLE_VAR(x) do {                                                                                  \
                                  std::cout << DEBUG_FORMAT << #x << " = " << x << std::endl;                         \
                                } while (0)

#define DEBUG_MSG_STD_VECTOR(x) do {                                                                                  \
                                  std::cout << DEBUG_FORMAT << "contents of '" << #x << "'" << std::endl;             \
                                  for (size_t x03ddj2d = 0; x03ddj2d < x.size(); ++x03ddj2d) {                        \
                                    std::cout << DEBUG_FORMAT << DEBUG_INDENT << #x"[" << x03ddj2d << "] = "          \
                                              << x[x03ddj2d] << std::endl;                                            \
                                  }                                                                                   \
                                } while (0)

#define DEBUG_MSG_STD_VECTOR_OF_VECTORS(x) do {                                                                       \
                                             std::cout << DEBUG_FORMAT << "contents of '" << #x << "'" << std::endl;  \
                                             for (size_t x03ddj2d = 0; x03ddj2d < x.size(); ++x03ddj2d) {             \
                                               for (size_t zb8ff23d = 0; zb8ff23d < x[x03ddj2d].size(); ++zb8ff23d) { \
                                                 std::cout << DEBUG_FORMAT << DEBUG_INDENT << #x << "[" << x03ddj2d   \
                                                           << "]" << "[" << zb8ff23d << "] = " << DEBUG_NUMBER_FORMAT \
                                                           << x[x03ddj2d][zb8ff23d] << std::endl;                     \
                                               }                                                                      \
                                             }                                                                        \
                                           } while (0)

#define DEBUG_MSG_STD_VECTOR_OF_MAPS(x) do {                                                                          \
                                          std::cout << DEBUG_FORMAT << "contents of '" << #x << "'" << std::endl;     \
                                          for (size_t x03ddj2d = 0; x03ddj2d < x.size(); ++x03ddj2d) {                \
                                            for (const auto& key_val : x[x03ddj2d]) {                                 \
                                              std::cout << DEBUG_FORMAT << DEBUG_INDENT << "[" << x03ddj2d            \
                                                        << "] --> [" << key_val.first << "] = " << key_val.second     \
                                                        << std::endl;                                                 \
                                            }                                                                         \
                                          }                                                                           \
                                        } while (0)

#define DEBUG_MSG_EIGEN_MATRIX(x) do {                                                                                \
                                    std::cout << DEBUG_FORMAT << "contents of '" << #x << "'" << std::endl;           \
                                    for (size_t x03ddj2d = 0; x03ddj2d < x.rows(); ++x03ddj2d) {                      \
                                      std::cout << DEBUG_FORMAT << DEBUG_INDENT;                                      \
                                      for (size_t zb8ff23d = 0; zb8ff23d < x.cols(); ++zb8ff23d) {                    \
                                        std::cout << DEBUG_NUMBER_FORMAT << x(x03ddj2d, zb8ff23d);                    \
                                        if (zb8ff23d != x.cols() - 1) {                                               \
                                          std::cout << DEBUG_INDENT;                                                  \
                                        }                                                                             \
                                      }                                                                               \
                                      std::cout << std::endl;                                                         \
                                    }                                                                                 \
                                  } while (0)

#else

#define DEBUG_MSG(x) do {} while (0)
#define DEBUG_MSG_SIMPLE_VAR(x) do {} while (0)
#define DEBUG_MSG_STD_VECTOR(x) do {} while (0)
#define DEBUG_MSG_STD_VECTOR_OF_VECTORS(x) do {} while (0)
#define DEBUG_MSG_STD_VECTOR_OF_MAPS(x) do {} while (0)
#define DEBUG_MSG_EIGEN_MATRIX(x) do {} while (0)

#endif

#endif //SD_COOLING_DEBUG_HPP
