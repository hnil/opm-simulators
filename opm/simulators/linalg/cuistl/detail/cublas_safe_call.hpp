/*
  Copyright 2022-2023 SINTEF AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CUBLAS_SAFE_CALL_HPP
#define CUBLAS_SAFE_CALL_HPP
#include <cublas_v2.h>
#include <exception>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#define CHECK_CUBLAS_ERROR_TYPE(code, x)                                                                               \
    if (code == x) {                                                                                                   \
        return #x;                                                                                                     \
    }
namespace
{
/**
 * @brief getCublasErrorMessage Converts an error code returned from a cublas function a human readable string.
 * @param code an error code from a cublas routine
 * @return a human readable string.
 */
inline std::string
getCublasErrorMessage(int code)
{
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_SUCCESS);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_NOT_INITIALIZED);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_ALLOC_FAILED);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_INVALID_VALUE);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_ARCH_MISMATCH);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_MAPPING_ERROR);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_EXECUTION_FAILED);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_INTERNAL_ERROR);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_NOT_SUPPORTED);
    CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_LICENSE_ERROR);

    return fmt::format("UNKNOWN CUBLAS ERROR {}.", code);
}
} // namespace
#undef CHECK_CUSPRASE_ERROR_TYPE

/**
 * @brief OPM_CUBLAS_SAFE_CALL checks the return type of the cublas expression (function call) and throws an exception
 * if it does not equal CUBLAS_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <cublas_v2.h>
 * #include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
 *
 * void some_function() {
 *     cublasHandle_t cublasHandle;
 *     OPM_CUBLAS_SAFE_CALL(cublasCreate(&cublasHandle));
 * }
 * @endcode
 *
 * @note This should be used for any call to cuBlas unless you have a good reason not to.
 */
#define OPM_CUBLAS_SAFE_CALL(expression)                                                                               \
    do {                                                                                                               \
        cublasStatus_t error = expression;                                                                             \
        if (error != CUBLAS_STATUS_SUCCESS) {                                                                          \
            OPM_THROW(std::runtime_error,                                                                              \
                      fmt::format("cuBLAS expression did not execute correctly. Expression was: \n\n"                  \
                                  "    {}\n\n"                                                                         \
                                  "in function {}, in {}, at line {}.\n"                                               \
                                  "CuBLAS error code was: {}\n",                                                       \
                                  #expression,                                                                         \
                                  __func__,                                                                            \
                                  __FILE__,                                                                            \
                                  __LINE__,                                                                            \
                                  getCublasErrorMessage(error)));                                                      \
        }                                                                                                              \
    } while (false)
#endif // CUBLAS_SAFE_CALL_HPP
