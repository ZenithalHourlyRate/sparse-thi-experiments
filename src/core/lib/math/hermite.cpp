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

static bool IsNotEqualZero(std::complex<double> v) {
    // TODO: tune this delta value during the fbt refactor
    constexpr double delta = 0x1p-32;  // 2**-32
    return (std::fabs(v.real()) >= delta) || (std::fabs(v.imag()) >= delta);
}
namespace lbcrypto {

std::ostream& operator<<(std::ostream& s, InterpolationMethod m) {
    switch (m) {
        case HERMITE_AKP25:
            s << "AKP25";
            break;
        case HERMITE_BKSS24:
            s << "BKSS24";
            break;
        case HERMITE_FULL_THI:
            s << "FULL THI";
            break;
        case HERMITE_SPARSE_THI:
            s << "Sparse THI";
            break;
        case HERMITE_BKSS24_NEW:
            s << "BKSS24_NEW";
            break;
        default:
            s << "UNKNOWN";
            break;
    }
    return s;
}

//=============================================================================
// AKP
//=============================================================================

std::vector<std::complex<double>> internalGetHermiteTrigCoefficients(std::function<double(int64_t)> func, uint32_t p,
                                                                     size_t order, double scale) {
    using namespace std::complex_literals;
    if (p == 0)
        OPENFHE_THROW("The degree of approximation can not be zero");

    switch (order) {
        case 1: {
            uint32_t degree = 0;
            std::vector<std::complex<double>> coeffs(p);

            for (uint32_t i = 0; i < p; ++i) {
                for (uint32_t j = 0; j < p; ++j)
                    coeffs[i] += static_cast<double>(func(j)) * std::exp((-2. * M_PI * i * j / p) * 1i);
                // No multiplication by 2 is to account for taking the real part
                coeffs[i] *= static_cast<double>(p - i) / static_cast<double>(p * p) / scale;
                if (IsNotEqualZero(coeffs[i]))
                    degree = i;
            }
            coeffs[0] /= 2.0;
            coeffs.resize(degree + 1);
            return coeffs;
        } break;
        case 2: {
            uint32_t pby2{p >> 1};
            uint32_t coeffTotal{p + pby2 + 1};
            std::vector<std::complex<double>> coeffs(coeffTotal);
            std::vector<std::complex<double>> alpha(p);
            std::vector<std::complex<double>> beta(pby2);
            std::vector<double> gamma(pby2);
            std::vector<std::complex<double>> delta(pby2);
            std::vector<std::complex<double>> omega(pby2);

            for (uint32_t i = 0; i < p; ++i) {
                for (uint32_t j = 0; j < p; ++j)
                    alpha[i] += static_cast<double>(func(j)) * std::exp((-2. * M_PI * i * j / p) * 1i);
                // The last /2 is to account for taking the real part
                alpha[i] *= 2. * static_cast<double>(p - i) / static_cast<double>(p * p) / 2. / scale;
            }
            alpha[0] /= 2.0;

            if ((p & 1) == 0)
                gamma.back() = 1.0;

            double factor = 1.0;
            for (uint32_t i = 1; i <= pby2; ++i) {
                for (uint32_t j = 0; j < p; ++j) {
                    auto y = static_cast<double>(func(j));
                    beta[i - 1] += y * std::exp((-2. * M_PI * i * j / p) * 1i);
                    delta[i - 1] += y * std::exp((-2. * M_PI * (p + i) * j / p) * 1i);
                    omega[i - 1] += y * std::exp((-2. * M_PI * (p - i) * j / p) * 1i);
                }
                // The last /2 is to account for taking the real part
                // factor = (2. - gamma[i - 1]) * i * static_cast<double>(p - i) / static_cast<double>(p * p * p) / 2. / scale;
                factor = (2. - gamma[i - 1]) * i * static_cast<double>(p - i) / static_cast<double>(p * p) /
                         static_cast<double>(p) / 2. /
                         scale;  // for large p, p*p*p overflows, so we separate the division
                beta[i - 1] *= factor;
                delta[i - 1] *= factor / 2.;
                omega[i - 1] *= factor / 2.;
            }

            uint32_t degree = 0;
            coeffs[0]       = alpha[0];
            for (uint32_t i = 1; i < coeffTotal; ++i) {
                if (i < p)
                    coeffs[i] = alpha[i];
                if (i <= pby2)
                    coeffs[i] += beta[i - 1];
                if (pby2 <= i && i < p)
                    coeffs[i] -= omega[p - i - 1];
                if (i > p)
                    coeffs[i] -= delta[i - p - 1];
                if (IsNotEqualZero(coeffs[i]))
                    degree = i;
            }
            coeffs.resize(degree + 1);
            return coeffs;
        } break;
        case 3: {
            uint32_t coeffTotal{p + p};
            std::vector<std::complex<double>> coeffs(coeffTotal);
            std::vector<std::complex<double>> alpha(p);
            std::vector<std::complex<double>> beta(p - 1);
            std::vector<std::complex<double>> delta(p - 1);
            std::vector<std::complex<double>> omega(p - 1);

            for (uint32_t i = 0; i < p; ++i) {
                for (uint32_t j = 0; j < p; ++j)
                    alpha[i] += static_cast<double>(func(j)) * std::exp((-2. * M_PI * i * j / p) * 1i);
                // The last /2 is to account for taking the real part
                alpha[i] *= 2. * static_cast<double>(p - i) / static_cast<double>(p * p) / 2. / scale;
            }
            alpha[0] /= 2.0;

            double factor = 1.0;
            for (uint32_t i = 1; i <= p - 1; ++i) {
                for (uint32_t j = 0; j < p; ++j) {
                    auto y = static_cast<double>(func(j));
                    beta[i - 1] += y * std::exp((-2. * M_PI * i * j / p) * 1i);
                    delta[i - 1] += y * std::exp((-2. * M_PI * (p + i) * j / p) * 1i);
                    omega[i - 1] += y * std::exp((-2. * M_PI * (p - i) * j / p) * 1i);
                }
                // The last /2 is to account for taking the real part
                factor = 2. * i * static_cast<double>(p - i) * static_cast<double>(2. * p - i) / 3. /
                         static_cast<double>(p * p) / static_cast<double>(p * p) / 2. /
                         scale;  // for large p, p*p*p*p overflows, so we separate the division
                beta[i - 1] *= factor;
                delta[i - 1] *= factor / 2.;
                omega[i - 1] *= factor / 2.;
            }

            uint32_t degree = 0;
            coeffs[0]       = alpha[0];
            for (uint32_t i = 1; i < coeffTotal; ++i) {
                if (i < p)
                    coeffs[i] = alpha[i];
                if (i <= p - 1)
                    coeffs[i] += beta[i - 1];
                if (1 <= i && i < p)
                    coeffs[i] -= omega[p - i - 1];
                if (i > p)
                    coeffs[i] -= delta[i - p - 1];
                if (IsNotEqualZero(coeffs[i]))
                    degree = i;
            }
            coeffs.resize(degree + 1);
            return coeffs;
        } break;
        default:
            OPENFHE_THROW("Order must be 1, 2, or 3");
    }
}

std::vector<std::complex<double>> GetHermiteTrigCoefficients(std::function<int64_t(int64_t)> func, uint32_t p,
                                                             size_t order, double scale) {
    auto funcDouble = [&](int64_t x) -> double {
        return static_cast<double>(func(x));
    };
    return internalGetHermiteTrigCoefficients(funcDouble, p, order, scale);
}

std::vector<std::complex<double>> GetHermiteTrigCoefficientsForComplexLUT(
    std::function<std::complex<double>(int64_t)> func, uint32_t p, size_t order, double scale) {
    auto funcReal = [&](int64_t x) -> double {
        return func(x).real();
    };
    auto funcImag = [&](int64_t x) -> double {
        return func(x).imag();
    };
    auto coeffsReal = internalGetHermiteTrigCoefficients(funcReal, p, order, scale);
    auto coeffsImag = internalGetHermiteTrigCoefficients(funcImag, p, order, scale);
    // append coeffsImag to coeffsReal
    coeffsReal.insert(coeffsReal.end(), coeffsImag.begin(), coeffsImag.end());
    return coeffsReal;
}

//=============================================================================
// FULL THI
//=============================================================================

std::vector<std::complex<double>> GetHermiteTrigCoefficientsFullTHIForComplexLUT(
    std::function<std::complex<double>(int64_t)> func, uint32_t p, size_t order, double scale) {
    using namespace std::complex_literals;
    auto omega = std::exp(2i * M_PI / double(p));

    // Compute IDFT
    std::vector<std::complex<double>> idft;
    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            ret += func(ell) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        idft.push_back(ret);
    }

    // Compute coeffs
    std::vector<std::complex<double>> coeffs;
    if (order == 1) {
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back((1.0 + double(m) / p) * idft[m]);
        }
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back((-double(m) / p) * idft[m]);
        }
    }
    if (order == 2) {
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back((1.0 + double(m) * (double(m) + 3 * p) / (2.0 * p * p)) * idft[m]);
        }
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back((-double(m) * (double(m) + 2 * p) / (double(p) * p)) * idft[m]);
        }
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back((double(m) * (double(m) + p) / (2.0 * double(p) * p)) * idft[m]);
        }
    }
    if (order == 3) {
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back(((double(m) + p) * (double(m) + 2 * p) * (double(m) + 3 * p) / (6.0 * p * p * p)) *
                             idft[m]);
        }
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back(-((double(m)) * (double(m) + 2 * p) * (double(m) + 3 * p) / (2.0 * p * p * p)) * idft[m]);
        }
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back(((double(m)) * (double(m) + p) * (double(m) + 3 * p) / (2.0 * p * p * p)) * idft[m]);
        }
        for (size_t m = 0; m != p; ++m) {
            coeffs.push_back(-(double(m) * (double(m) + p) * (double(m) + 2 * p) / (6.0 * p * p * p)) * idft[m]);
        }
    }

    // normalization
    for (size_t m = 0; m != coeffs.size(); ++m) {
        coeffs[m] /= scale;
    }

    return coeffs;
}

std::vector<std::complex<double>> GetHermiteTrigCoefficientsFullTHI(std::function<int64_t(int64_t)> func, uint32_t p,
                                                                    size_t order, double scale) {
    auto funcComplex = [&](int64_t x) -> std::complex<double> {
        return std::complex<double>(static_cast<double>(func(x)), 0.0);
    };
    return GetHermiteTrigCoefficientsFullTHIForComplexLUT(funcComplex, p, order, scale);
}

//=============================================================================
// SparseTHI
//=============================================================================

std::vector<std::complex<double>> internalGetHermiteTrigCoefficientsSparseTHI(std::function<double(int64_t)> func,
                                                                              uint32_t p, size_t order, double scale) {
    using namespace std::complex_literals;
    auto omega = std::exp(2i * M_PI / double(p));

    // Compute IDFT
    std::vector<std::complex<double>> idft;
    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            ret += func(ell) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        idft.push_back(ret);
    }

    // MixedConstraints / HYBRID for arbitrary order (Revisit.pdf Theorem 6)
    std::vector<std::complex<double>> coeffs(((2 * order + 1) * p) / 2 + 1, 0);
    coeffs[0] = idft[0];

    auto factorial = [&](size_t n) {
        double out = 1.0;
        for (size_t i = 2; i <= n; ++i)
            out *= static_cast<double>(i);
        return out;
    };
    auto binomial = [&](size_t n, size_t k) {
        if (k > n)
            return 0.0;
        if (k > n - k)
            k = n - k;
        double out = 1.0;
        for (size_t i = 1; i <= k; ++i) {
            out *= static_cast<double>(n - k + i);
            out /= static_cast<double>(i);
        }
        return out;
    };

    const double orderFactorial = factorial(order);
    const double pPower         = std::pow(static_cast<double>(p), static_cast<double>(order));

    for (size_t ell = 0; ell <= order; ++ell) {
        const double sign   = (ell % 2 == 0) ? 1.0 : -1.0;
        const double choose = binomial(order, ell);

        double nyquistProd = 1.0;
        for (size_t j = 0; j <= order; ++j) {
            if (j != ell)
                nyquistProd *= static_cast<double>(2 * j + 1);
        }
        coeffs[ell * p + p / 2] = sign * choose * nyquistProd * idft[p / 2].real() /
                                  (std::pow(2.0, static_cast<double>(order)) * orderFactorial);

        for (size_t m = 1; m < p / 2; ++m) {
            double prod = 1.0;
            for (size_t j = 0; j <= order; ++j) {
                if (j != ell)
                    prod *= static_cast<double>(j * p + m);
            }
            coeffs[ell * p + m] = sign * choose * 2.0 * prod * idft[m] / (orderFactorial * pPower);
        }
    }

    // normalization
    // The last /2 is to account for taking the real part
    for (size_t m = 0; m != coeffs.size(); ++m) {
        coeffs[m] /= scale * 2;
    }

    return coeffs;
}

std::vector<std::complex<double>> GetHermiteTrigCoefficientsSparseTHI(std::function<int64_t(int64_t)> func, uint32_t p,
                                                                      size_t order, double scale) {
    auto funcDouble = [&](int64_t x) -> double {
        return static_cast<double>(func(x));
    };
    return internalGetHermiteTrigCoefficientsSparseTHI(funcDouble, p, order, scale);
}

std::vector<std::complex<double>> GetHermiteTrigCoefficientsSparseTHIForComplexLUT(
    std::function<std::complex<double>(int64_t)> func, uint32_t p, size_t order, double scale) {
    auto funcReal = [&](int64_t x) -> double {
        return func(x).real();
    };
    auto funcImag = [&](int64_t x) -> double {
        return func(x).imag();
    };
    auto coeffsReal = internalGetHermiteTrigCoefficientsSparseTHI(funcReal, p, order, scale);
    auto coeffsImag = internalGetHermiteTrigCoefficientsSparseTHI(funcImag, p, order, scale);
    // append coeffsImag to coeffsReal
    coeffsReal.insert(coeffsReal.end(), coeffsImag.begin(), coeffsImag.end());
    return coeffsReal;
}

//=============================================================================
// BKSS24
//=============================================================================

std::vector<std::complex<double>> GetHermiteTrigCoefficientsBKSSForComplexLUT(
    std::function<std::complex<double>(int64_t)> func, uint32_t p, double scale) {
    using namespace std::complex_literals;
    auto omega = std::exp(2i * M_PI / double(p));

    // Compute IDFT
    std::vector<std::complex<double>> idft;
    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            ret += func(ell) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        idft.push_back(ret);
    }

    std::vector<std::complex<double>> f0 = {idft[0]};
    std::vector<std::complex<double>> Pf(p, 0);
    std::vector<std::complex<double>> Qf(p, 0);
    std::vector<std::complex<double>> Pfr(p, 0);
    std::vector<std::complex<double>> Qfr(p, 0);

    for (size_t k = 1; k <= p / 2; ++k) {
        Pf[k] = idft[k] * (double(p) - k) * (double(k) + 1) / double(p);
        Qf[k] = idft[k] * double(k) * (double(p) - k) / double(p);
    }
    for (size_t k = p / 2 + 1; k < p; ++k) {
        Pf[k] = idft[k] * (double(p) - k) / double(p);
    }
    // Their formula in paper for Pfr, Qfr is wrong...
    for (size_t k = 1; k <= p / 2 - 1; ++k) {
        Pfr[k] = idft[p - k] * (double(p) - k) * (double(k) + 1) / double(p);
        Qfr[k] = idft[p - k] * double(k) * (double(p) - k) / double(p);
    }
    for (size_t k = p / 2; k < p; ++k) {
        Pfr[k] = idft[p - k] * (double(p) - k) / double(p);
    }

    // normalization
    f0[0] /= scale;
    for (size_t m = 0; m != p; ++m) {
        Pf[m] /= scale;
        Qf[m] /= scale;
        Pfr[m] /= scale;
        Qfr[m] /= scale;
    }

    // Pack them for convenience of APIs...
    // {Pf, Pfr, Qf, Qfr, f0}
    std::vector<std::complex<double>> coeffs;
    coeffs.reserve(1 + 4 * p);
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(Pf[k]);
    }
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(Pfr[k]);
    }
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(Qf[k]);
    }
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(Qfr[k]);
    }
    coeffs.push_back(f0[0]);
    return coeffs;
}

std::vector<std::complex<double>> GetHermiteTrigCoefficientsBKSS(std::function<int64_t(int64_t)> func, uint32_t p,
                                                                 double scale) {
    auto funcComplex = [&](int64_t x) -> std::complex<double> {
        return std::complex<double>(static_cast<double>(func(x)), 0.0);
    };
    return GetHermiteTrigCoefficientsBKSSForComplexLUT(funcComplex, p, scale);
}

std::vector<std::complex<double>> GetHermiteTrigCoefficientsBKSSNew(std::function<int64_t(int64_t)> func, uint32_t p,
                                                                    double scale) {
    using namespace std::complex_literals;
    auto omega = std::exp(2i * M_PI / double(p));

    // Compute IDFT
    std::vector<std::complex<double>> idft;
    for (size_t m = 0; m != p; ++m) {
        std::complex<double> ret = 0;
        for (size_t ell = 0; ell != p; ++ell) {
            ret += static_cast<double>(func(ell)) * std::pow(omega, -double(ell) * m);
        }
        ret /= double(p);
        idft.push_back(ret);
    }

    std::vector<std::complex<double>> f0 = {idft[0]};
    std::vector<std::complex<double>> Pf(p, 0);
    std::vector<std::complex<double>> Qf(p, 0);
    std::complex<double> rem;

    for (size_t k = 1; k <= p / 2 - 1; ++k) {
        Pf[k] = idft[k] * (double(p) - k) * (double(k) + 1) / double(p);
        Qf[k] = idft[k] * double(k) * (double(p) - k) / double(p);
    }
    for (size_t k = p / 2 + 1; k < p; ++k) {
        Pf[k] = idft[k] * (double(p) - k) / double(p);
    }
    Pf[p / 2] = idft[p / 2] * (double(p) / 2.0) / double(p);
    rem       = idft[p / 2] * (double(p) / 2.0) * (double(p) / 2.0) / double(p);

    // 2Re[Pf(z)] + 2Re[Qf(z)] + L_{p/2} / p * (p/2) * (p/2) * (z^{p/2}) * (1-z*conj(z))

    // normalization
    f0[0] /= scale;
    rem /= scale;
    for (size_t m = 0; m != p; ++m) {
        Pf[m] /= scale;
        Qf[m] /= scale;
    }

    // Pack them for convenience of APIs...
    // {Pf, Pfr (0), Qf, Qfr (0), f0}
    std::vector<std::complex<double>> coeffs;
    coeffs.reserve(1 + 4 * p);
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(Pf[k]);
    }
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(0.0);
    }
    coeffs[p] = rem;
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(Qf[k]);
    }
    for (size_t k = 0; k < p; ++k) {
        coeffs.push_back(0.0);
    }
    coeffs.push_back(f0[0]);
    return coeffs;
}

}  // namespace lbcrypto
