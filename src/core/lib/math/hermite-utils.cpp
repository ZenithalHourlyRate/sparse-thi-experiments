//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2023, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  This code provides Chebyshev approximation utilities
 */

#include "math/hermite.h"
#include "utils/exception.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace lbcrypto {

namespace {
double gamma_weight(int n, int k, int p) {
    double kp = static_cast<double>(k) / p;
    return std::tgamma(n + 1 + kp) / (std::tgamma(n + 2) * std::tgamma(kp));
}
};  // namespace

double GetHermiteTrigAKPThreshold(std::function<double(int64_t)> func, uint32_t p, size_t order, double scale) {
    using namespace std::complex_literals;

    auto coeffs = GetHermiteTrigCoefficients(func, p, order, scale);
    auto omega  = std::exp(2i * M_PI / double(p));

    double max_value = -1.0;

    // factorial(order + 1)
    double factorial = std::tgamma(double(order) + 2.0);

    for (uint32_t x = 0; x < p; ++x) {
        std::complex<double> series(0.0, 0.0);

        for (uint32_t k = 1; k < p; ++k) {
            std::complex<double> derivative_factor = std::pow(2.0 * M_PI * 1i * double(k), double(order + 1));

            std::complex<double> omega_pow = std::pow(omega, double(k) * double(x));

            series += derivative_factor * coeffs[k] * omega_pow;
        }

        series /= factorial;

        double magnitude = std::abs(series.real());
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }

    // 1.0 for 2Re[]
    double t_bits = -(1.0 + std::log2(max_value)) / double(order);
    return t_bits;
}

double GetHermiteTrigSparseTHIThreshold(std::function<double(int64_t)> func, uint32_t p, size_t order, double scale) {
    using namespace std::complex_literals;
    auto omega = std::exp(2i * M_PI / double(p));

    // Compute IDFT
    std::vector<std::complex<double>> idft;
    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            // the extra /scale for postScaling/scaleTHI
            ret += (double(func(ell))) / double(scale) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        idft.push_back(ret);
    }
    // now the effective polynomial
    std::vector<std::complex<double>> e(p / 2 + 1);
    // extra /2 for computing Re
    e[0] = idft[0] / 2.0;  // actually not used
    for (size_t m = 1; m < p / 2; ++m) {
        e[m] = idft[m];
    }
    e[p / 2] = idft[p / 2] / 2.0;

    // Compute the gammas
    std::vector<double> weights(p / 2 + 1, 0.0);
    for (int k = 1; k <= p / 2; ++k) {
        weights[k] = gamma_weight(order, k, p);
    }

    double max_value = -1.0;

    for (int x = 0; x < p; ++x) {
        std::complex<double> series(0.0, 0.0);

        for (int k = 1; k <= p / 2; ++k) {
            // exp(2j * pi * k * x / p)
            double angle = 2.0 * M_PI * k * x / p;
            std::complex<double> omega_pow(std::cos(angle), std::sin(angle));

            series += e[k] * omega_pow * weights[k];
        }

        double magnitude = std::abs(series);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }

    // t_bits = -(1 + 1/n) * log2(2 * pi * p) - (log2(max_value) + 1) / n
    double log2_2pip = std::log2(2.0 * M_PI * p);
    double t_bits    = -(1.0 + 1.0 / order) * log2_2pip - (std::log2(max_value) + 1.0) / order;

    return t_bits;
}

double GetHermiteTrigFullTHIThreshold(std::function<std::complex<double>(int64_t)> func, uint32_t p, size_t order,
                                      double scale) {
    using namespace std::complex_literals;
    auto omega = std::exp(2i * M_PI / double(p));

    // Compute IDFT
    std::vector<std::complex<double>> idft;
    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            // the extra /scale for postScaling/scaleTHI
            ret += func(ell) / double(scale) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        idft.push_back(ret);
    }

    // Compute the gammas
    std::vector<double> weights(p, 0.0);
    for (int k = 1; k < p; ++k) {
        weights[k] = gamma_weight(order, k, p);
    }

    double max_value = -1.0;

    for (int x = 0; x < p; ++x) {
        std::complex<double> series(0.0, 0.0);

        for (int k = 1; k < p; ++k) {
            // exp(2j * pi * k * x / p)
            double angle = 2.0 * M_PI * k * x / p;
            std::complex<double> omega_pow(std::cos(angle), std::sin(angle));

            series += idft[k] * omega_pow * weights[k];
        }

        double magnitude = std::abs(series);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }

    // t_bits = -(1 + 1/n) * log2(2 * pi * p) - (log2(max_value)) / n
    double log2_2pip = std::log2(2.0 * M_PI * p);
    double t_bits    = -(1.0 + 1.0 / order) * log2_2pip - (std::log2(max_value)) / order;

    return t_bits;
}

double GetHermiteTrigSparseTHIThresholdForComplexLUT(std::function<std::complex<double>(int64_t)> func, uint32_t p,
                                                     size_t order, double scale) {
    auto real_func = [&func](int64_t x) {
        std::complex<double> val = func(x);
        return (val.real());
    };
    auto imag_func = [&func](int64_t x) {
        std::complex<double> val = func(x);
        return (val.imag());
    };
    auto t_bits_real = GetHermiteTrigSparseTHIThreshold(real_func, p, order, scale);
    auto t_bits_imag = GetHermiteTrigSparseTHIThreshold(imag_func, p, order, scale);
    return std::min(t_bits_real, t_bits_imag) - 0.5 / double(order);
}

double GetHermiteTrigAKPThresholdForComplexLUT(std::function<std::complex<double>(int64_t)> func, uint32_t p,
                                               size_t order, double scale) {
    auto real_func = [&func](int64_t x) {
        std::complex<double> val = func(x);
        return (val.real());
    };
    auto imag_func = [&func](int64_t x) {
        std::complex<double> val = func(x);
        return (val.imag());
    };
    auto t_bits_real = GetHermiteTrigAKPThreshold(real_func, p, order, scale);
    auto t_bits_imag = GetHermiteTrigAKPThreshold(imag_func, p, order, scale);
    return std::min(t_bits_real, t_bits_imag) - 0.5 / double(order);
}

double PredictHermiteTrigNoise(double t_bits, size_t order, double inputPrec) {
    return order * (inputPrec - t_bits) + inputPrec;
}

}  // namespace lbcrypto
