
#include "math/hermite.h"
#include "openfhe.h"

#include <functional>

using namespace lbcrypto;
using CiphertextT        = ConstCiphertext<DCRTPoly>;
using MutableCiphertextT = Ciphertext<DCRTPoly>;
using CCParamsT          = CCParams<CryptoContextBFVRNS>;
using CryptoContextT     = CryptoContext<DCRTPoly>;
using EvalKeyT           = EvalKey<DCRTPoly>;
using PlaintextT         = Plaintext;
using PrivateKeyT        = PrivateKey<DCRTPoly>;
using PublicKeyT         = PublicKey<DCRTPoly>;

// DecryptCore not accessible from CryptoContext
// so copy from @openfhe//src/pke/lib/schemerns/rns-pke.cpp
DCRTPoly DecryptCore(const std::vector<DCRTPoly>& cv, const PrivateKey<DCRTPoly> privateKey) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = cv[0].GetParams()->GetParams().size();

    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);

    DCRTPoly b(cv[0]);
    b.SetFormat(Format::EVALUATION);

    DCRTPoly ci;
    for (size_t i = 1; i < cv.size(); i++) {
        ci = cv[i];
        ci.SetFormat(Format::EVALUATION);

        b += sPower * ci;
        sPower *= scopy;
    }
    return b;
}

void unused_func() {
    // auto overflow_and_middle = x >> (bigq.GetMSB() - PInput_global.GetMSB());
    // auto overflow            = x >> (bigq.GetMSB() - 1);
    // auto overflowBig         = overflow << (PInput_global.GetMSB() - 1);
    // auto middle              = overflow_and_middle - overflowBig;
    // auto noise               = x - (overflow_and_middle << (bigq.GetMSB() - PInput_global.GetMSB()));

    //if (msg == "PartialSum") {
    //    std::cout << msg << " coeff[" << i << "] = " << (neg ? "-" : "") << x
    //              << " Middle: " << (neg ? PInput_global - middle : middle) << " Overflow: " << overflow
    //              << " Value: " << (neg ? "-" : "") << std::setprecision(20)
    //              << (overflow_and_middle.ConvertToDouble() / PInput_global.ConvertToDouble())
    //              << " Raw: " << (x.ConvertToDouble() / bigq.ConvertToDouble()) << std::endl;
    //}

    // for (size_t i = 0; i < 16; i++) {
    //     std::cout << std::defaultfloat << msg << "  Slot FFT[" << i << "] = " << std::setprecision(20)
    //               << result[i].real() << std::defaultfloat << " diff " << raws[i] - (result[i].real())
    //               << " diff log2 " << std::log2(std::abs(raws[i] - (result[i].real()))) << std::endl;
    // }

    // bit reverse the result
    // auto BitReverse = [](std::vector<std::complex<double>>& vals) {
    //     uint32_t size = vals.size();
    //     for (size_t i = 1, j = 0; i < size; ++i) {
    //         size_t bit = size >> 1;
    //         for (; j >= bit; bit >>= 1) {
    //             j -= bit;
    //         }
    //         j += bit;
    //         if (i < j) {
    //             swap(vals[i], vals[j]);
    //         }
    //     }
    // };
    // BitReverse(result);

    //std::cout << "  After CoeffsToSlots BitReverse:" << std::endl;
    //for (size_t i = 0; i < 16; i++) {
    //    std::cout << "    Slot FFT[" << i << "] = " << std::setprecision(20) << result[i].real() * 25 << std::endl;
    //}
}

std::vector<std::complex<double>> multiplyByUComplex(std::vector<double> input) {
    // Original computes:
    // result[i] = sum_{j=0}^{2*slots-1} input[j] * exp(2*pi*i*j*rotGroup[i]/m)
    //
    // This is an unnormalized inverse DFT of length m = 4*slots,
    // where input is zero-padded from length 2*slots to length m.

    if (input.size() % 2 != 0) {
        throw std::invalid_argument("input.size() must be even");
    }

    const size_t slots = input.size() / 2;
    const size_t m     = 4 * slots;

    if (m == 0 || (m & (m - 1)) != 0) {
        throw std::invalid_argument("m = 4 * slots must be a power of two");
    }

    std::vector<std::complex<double>> a(m, std::complex<double>(0.0, 0.0));
    for (size_t i = 0; i < input.size(); ++i) {
        a[i] = std::complex<double>(input[i], 0.0);
    }

    // Bit-reversal permutation.
    for (size_t i = 1, j = 0; i < m; ++i) {
        size_t bit = m >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;

        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // Unnormalized inverse FFT:
    // a[k] = sum_j input[j] * exp(+2*pi*i*j*k/m)
    const double pi = std::acos(-1.0);

    for (size_t len = 2; len <= m; len <<= 1) {
        const double angle = 2.0 * pi / static_cast<double>(len);
        const std::complex<double> wLen(std::cos(angle), std::sin(angle));

        for (size_t start = 0; start < m; start += len) {
            std::complex<double> w(1.0, 0.0);

            for (size_t j = 0; j < len / 2; ++j) {
                const std::complex<double> u = a[start + j];
                const std::complex<double> v = a[start + j + len / 2] * w;

                a[start + j]           = u + v;
                a[start + j + len / 2] = u - v;

                w *= wLen;
            }
        }
    }

    // Gather primitive-root indices: 1, 5, 5^2, ... mod m.
    std::vector<std::complex<double>> result(slots);

    const size_t mmask = m - 1;
    size_t fivePows    = 1;

    for (size_t i = 0; i < slots; ++i) {
        result[i] = a[fivePows & mmask];
        fivePows  = (fivePows * 5) & mmask;
    }

    return result;
}

std::vector<double> multiplyByU(std::vector<double> input) {
    auto result = multiplyByUComplex(input);
    std::vector<double> realResult;
    for (size_t i = 0; i < result.size(); ++i) {
        realResult.push_back(result[i].real());
    }
    return realResult;
}

struct seriesPowersCustom {
    std::vector<std::complex<double>> powers;
    std::vector<std::complex<double>> powers2;
    std::complex<double> power2km1;
    int k;
    int m;
};

seriesPowersCustom evalPowersPS(std::complex<double> x, const std::vector<std::complex<double>>& coefficients) {
    seriesPowersCustom ret;

    auto n     = Degree(coefficients);
    auto degs  = ComputeDegreesPS(n);
    uint32_t k = degs[0];
    uint32_t m = degs[1];
    ret.k      = k;
    ret.m      = m;

    std::vector<std::complex<double>> powers;
    powers.reserve(k);
    powers.push_back(x);

    // computes all powers up to k for x
    uint32_t powerOf2 = 2;
    uint32_t rem      = 0;
    for (uint32_t i = 2; i <= k; i++) {
        if (rem == 0) {
            powers.push_back(powers[(powerOf2 >> 1) - 1] * powers[(powerOf2 >> 1) - 1]);
        }
        else {
            powers.push_back(powers[powerOf2 - 1] * powers[rem - 1]);
        }
        if (++rem == powerOf2) {
            powerOf2 <<= 1;
            rem = 0;
        }
    }

    // computes powers of form k*2^i for x and the product of the powers in power2, that yield x^{k(2*m - 1)}
    std::vector<std::complex<double>> powers2;
    powers2.reserve(m);
    powers2.push_back(powers.back());
    auto power2km1 = powers.back();

    for (uint32_t i = 1; i < m; i++) {
        powers2.push_back(powers2[i - 1] * powers2[i - 1]);
        power2km1 = power2km1 * powers2.back();
    }

    ret.powers    = powers;
    ret.powers2   = powers2;
    ret.power2km1 = power2km1;
    return ret;
}

static std::complex<double> InnerEvalPolyPSCustom(std::complex<double> x,
                                                  const std::vector<std::complex<double>>& coefficients, uint32_t k,
                                                  uint32_t m, std::vector<std::complex<double>>& powers,
                                                  std::vector<std::complex<double>>& powers2) {
    // Compute k*2^m because we use it often
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Divide coefficients by x^{k*2^{m-1}}
    std::vector<std::complex<double>> xkm(static_cast<int32_t>(k2m2k + k) + 1, 0.0);
    xkm.back() = 1;

    auto divqr = LongDivisionPoly(coefficients, xkm);

    // Subtract x^{k(2^{m-1} - 1)} from r
    auto r2 = divqr->r;
    if (static_cast<int32_t>(k2m2k - Degree(divqr->r)) <= 0) {
        r2[static_cast<int32_t>(k2m2k)] -= 1;
        r2.resize(Degree(r2) + 1);
    }
    else {
        r2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
        r2.back() = -1;
    }

    // Divide r2 by q
    auto divcs = LongDivisionPoly(r2, divqr->q);

    // Add x^{k(2^{m-1} - 1)} to s
    auto s2 = divcs->r;
    s2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
    s2.back() = 1;

    std::complex<double> cu;
    uint32_t dc = Degree(divcs->q);
    bool flag_c = false;

    if (dc >= 1) {
        if (dc == 1) {
            if (IsNotEqualOne(divcs->q[1])) {
                cu = powers.front() * divcs->q[1];
            }
            else {
                cu = powers.front();
            }
        }
        else {
            std::vector<std::complex<double>> ctxs(dc);
            std::vector<std::complex<double>> weights(dc);

            cu = 0;

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i]    = powers[i];
                weights[i] = divcs->q[i + 1];
                cu += ctxs[i] * weights[i];
            }
        }

        // adds the free term (at x^0)
        cu += divcs->q.front();
        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    std::complex<double> qu;

    if (Degree(divqr->q) > k) {
        qu = InnerEvalPolyPSCustom(x, divqr->q, k, m - 1, powers, powers2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto qcopy = divqr->q;
        qcopy.resize(k);
        if (Degree(qcopy) > 0) {
            std::vector<std::complex<double>> ctxs(Degree(qcopy));
            std::vector<std::complex<double>> weights(Degree(qcopy));

            qu = 0;

            for (uint32_t i = 0; i < Degree(qcopy); i++) {
                ctxs[i]    = powers[i];
                weights[i] = divqr->q[i + 1];
                qu += ctxs[i] * weights[i];
            }

            // the highest order term will always be 1 because q is monic
            qu += powers[k - 1];
        }
        else {
            qu = powers[k - 1];
        }
        // adds the free term (at x^0)
        qu += divqr->q.front();
    }

    uint32_t ds = Degree(s2);
    std::complex<double> su;

    if (std::equal(s2.begin(), s2.end(), divqr->q.begin())) {
        su = qu;
    }
    else {
        if (ds > k) {
            su = InnerEvalPolyPSCustom(x, s2, k, m - 1, powers, powers2);
        }
        else {
            // ds = k from construction
            // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
            auto scopy = s2;
            scopy.resize(k);
            if (Degree(scopy) > 0) {
                std::vector<std::complex<double>> ctxs(Degree(scopy));
                std::vector<std::complex<double>> weights(Degree(scopy));

                su = 0;

                for (uint32_t i = 0; i < Degree(scopy); ++i) {
                    ctxs[i]    = powers[i];
                    weights[i] = s2[i + 1];
                    su += ctxs[i] * weights[i];
                }

                // the highest order term will always be 1 because q is monic
                su += powers[k - 1];
            }
            else {
                su = powers[k - 1];
            }
            // adds the free term (at x^0)
            su += s2.front();
        }
    }

    std::complex<double> result;

    if (flag_c) {
        result = powers2[m - 1] + cu;
    }
    else {
        result = powers2[m - 1] + divcs->q.front();
    }

    result = result * qu;
    result += su;

    return result;
}

std::complex<double> evalPolyPSWithPrecompCustom(seriesPowersCustom ctxtPowers,
                                                 const std::vector<std::complex<double>>& coefficients) {
    auto f2 = coefficients;
    auto n  = Degree(f2);
    f2.resize(n + 1);

    auto powers    = ctxtPowers.powers;
    auto powers2   = ctxtPowers.powers2;
    auto power2km1 = ctxtPowers.power2km1;
    unsigned k     = ctxtPowers.k;
    auto m         = ctxtPowers.m;

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Add x^{k(2^m - 1)} to the polynomial that has to be evaluated
    // std::vector<double> f2 = coefficients;
    f2.resize(2 * k2m2k + k + 1, 0.0);
    f2.back() = 1;

    // Divide f2 by x^{k*2^{m-1}}
    std::vector<std::complex<double>> xkm(static_cast<int32_t>(k2m2k + k) + 1);
    xkm.back() = 1;
    auto divqr = LongDivisionPoly(f2, xkm);

    // Subtract x^{k(2^{m-1} - 1)} from r
    auto r2 = divqr->r;
    if (static_cast<int32_t>(k2m2k - Degree(divqr->r)) <= 0) {
        r2[static_cast<int32_t>(k2m2k)] -= 1;
        r2.resize(Degree(r2) + 1);
    }
    else {
        r2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
        r2.back() = -1;
    }

    // Divide r2 by q
    auto divcs = LongDivisionPoly(r2, divqr->q);

    // Add x^{k(2^{m-1} - 1)} to s
    auto s2 = divcs->r;
    s2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
    s2.back() = 1;

    uint32_t dc = Degree(divcs->q);
    bool flag_c = false;

    std::complex<double> cu;

    if (dc >= 1) {
        if (dc == 1) {
            if (IsNotEqualOne(divcs->q[1])) {
                cu = powers.front() * divcs->q[1];
            }
            else {
                cu = powers.front();
            }
        }
        else {
            std::vector<std::complex<double>> ctxs(dc);
            std::vector<std::complex<double>> weights(dc);

            cu = 0;

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i]    = powers[i];
                weights[i] = divcs->q[i + 1];
                cu += ctxs[i] * weights[i];
            }
        }

        // adds the free term (at x^0)
        cu += divcs->q.front();
        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    std::complex<double> qu;

    if (Degree(divqr->q) > k) {
        qu = InnerEvalPolyPSCustom(powers[0], divqr->q, k, m - 1, powers, powers2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto qcopy = divqr->q;
        qcopy.resize(k);
        if (Degree(qcopy) > 0) {
            std::vector<std::complex<double>> ctxs(Degree(qcopy));
            std::vector<std::complex<double>> weights(Degree(qcopy));

            qu = 0;

            for (uint32_t i = 0; i < Degree(qcopy); i++) {
                ctxs[i]    = powers[i];
                weights[i] = divqr->q[i + 1];
                qu += powers[i] * weights[i];
            }

            // the highest order term will always be 1 because q is monic
            qu += powers[k - 1];
        }
        else {
            qu = powers[k - 1];
        }
        // adds the free term (at x^0)
        qu += divqr->q.front();
    }

    uint32_t ds = Degree(s2);
    std::complex<double> su;

    if (std::equal(s2.begin(), s2.end(), divqr->q.begin())) {
        su = qu;
    }
    else {
        if (ds > k) {
            su = InnerEvalPolyPSCustom(powers[0], s2, k, m - 1, powers, powers2);
        }
        else {
            // ds = k from construction
            // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
            auto scopy = s2;
            scopy.resize(k);
            if (Degree(scopy) > 0) {
                std::vector<std::complex<double>> ctxs(Degree(scopy));
                std::vector<std::complex<double>> weights(Degree(scopy));

                su = 0;

                for (uint32_t i = 0; i < Degree(scopy); i++) {
                    ctxs[i]    = powers[i];
                    weights[i] = s2[i + 1];
                    su += ctxs[i] * weights[i];
                }

                // the highest order term will always be 1 because q is monic
                su += powers[k - 1];
            }
            else {
                su = powers[k - 1];
            }
            // adds the free term (at x^0)
            su += s2.front();
        }
    }

    std::complex<double> result;

    if (flag_c) {
        result = powers2[m - 1] + cu;
    }
    else {
        result = powers2[m - 1] + divcs->q.front();
    }

    result = result * qu;
    result += su;
    result -= power2km1;

    return result;
}

// Cheby PS...

seriesPowersCustom internalEvalChebyPolysPS(const std::complex<double>& x,
                                            const std::vector<std::complex<double>>& coefficients, double a, double b) {
    seriesPowersCustom ret;
    auto n     = Degree(coefficients);
    auto degs  = ComputeDegreesPS(n);
    uint32_t k = degs[0];
    uint32_t m = degs[1];

    // computes linear transformation y = -1 + 2 (x-a)/(b-a)
    // consumes one level when a <> -1 && b <> 1
    std::vector<std::complex<double>> T(k);
    if ((a - std::round(a) < 1e-10) && (b - std::round(b) < 1e-10) && (std::round(a) == -1.0) &&
        (std::round(b) == 1.0)) {
        // no linear transformation is needed if a = -1, b = 1
        // T_1(y) = y
        T[0] = x;
    }
    else {
        // linear transformation is needed
        double alpha = 2 / (b - a);
        double beta  = 2 * a / (b - a);

        T[0] = x * alpha;
        T[0] += -1.0 - beta;
    }

    std::complex<double> y = T[0];

    // Computes Chebyshev polynomials up to degree k
    // for y: T_1(y) = y, T_2(y), ... , T_k(y)
    // uses binary tree multiplication
    for (uint32_t i = 2; i <= k; ++i) {
        // if i is a power of two
        if (!(i & (i - 1))) {
            // compute T_{2i}(y) = 2*T_i(y)^2 - 1
            auto square = T[i / 2 - 1] * T[i / 2 - 1];
            T[i - 1]    = square + square;
            T[i - 1] += -1.0;
        }
        else {
            // non-power of 2
            if (i % 2 == 1) {
                // if i is odd
                // compute T_{2i+1}(y) = 2*T_i(y)*T_{i+1}(y) - y
                auto prod = T[i / 2 - 1] * T[i / 2];
                T[i - 1]  = prod + prod;

                T[i - 1] -= y;
            }
            else {
                // i is even but not power of 2
                // compute T_{2i}(y) = 2*T_i(y)^2 - 1
                auto square = T[i / 2 - 1] * T[i / 2 - 1];
                T[i - 1]    = square + square;
                T[i - 1] += -1.0;
            }
        }
    }

    std::vector<std::complex<double>> T2(m);
    // Compute the Chebyshev polynomials T_k(y), T_{2k}(y), T_{4k}(y), ... , T_{2^{m-1}k}(y)
    // T2[0] is used as a placeholder
    T2.front() = T.back();
    for (uint32_t i = 1; i < m; i++) {
        auto square = T2[i - 1] * T2[i - 1];
        T2[i]       = square + square;
        T2[i] += -1.0;
    }

    // computes T_{k(2*m - 1)}(y)
    auto T2km1 = T2.front();
    for (uint32_t i = 1; i < m; i++) {
        // compute T_{k(2*m - 1)} = 2*T_{k(2^{m-1}-1)}(y)*T_{k*2^{m-1}}(y) - T_k(y)
        auto prod = T2km1 * T2[i];
        T2km1     = prod + prod;
        T2km1 -= T2.front();
    }

    // We also need to reduce the number of levels of T[k-1] and of T2[0] by another level.
    //  cc->LevelReduceInPlace(T[k-1], nullptr);
    //  cc->LevelReduceInPlace(T2.front(), nullptr);

    ret.powers    = T;
    ret.powers2   = T2;
    ret.power2km1 = T2km1;
    ret.k         = k;
    ret.m         = m;

    return ret;
}

static std::complex<double> InnerEvalChebyshevPS(std::complex<double>& x,
                                                 const std::vector<std::complex<double>>& coefficients, uint32_t k,
                                                 uint32_t m, std::vector<std::complex<double>>& T,
                                                 std::vector<std::complex<double>>& T2) {
    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Divide coefficients by T^{k*2^{m-1}}
    std::vector<std::complex<double>> Tkm(static_cast<int32_t>(k2m2k + k) + 1);
    Tkm.back() = 1;
    auto divqr = LongDivisionChebyshev(coefficients, Tkm);

    // Subtract x^{k(2^{m-1} - 1)} from r
    auto r2 = divqr->r;
    if (static_cast<int32_t>(k2m2k - Degree(divqr->r)) <= 0) {
        r2[static_cast<int32_t>(k2m2k)] -= 1;
        r2.resize(Degree(r2) + 1);
    }
    else {
        r2.resize(static_cast<int32_t>(k2m2k + 1));
        r2.back() = -1;
    }

    // Divide r2 by q
    auto divcs = LongDivisionChebyshev(r2, divqr->q);

    // Add x^{k(2^{m-1} - 1)} to s
    auto s2 = divcs->r;
    s2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
    s2.back() = 1;

    // Evaluate c at u
    std::complex<double> cu;
    uint32_t dc = Degree(divcs->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            if (IsNotEqualOne(divcs->q[1])) {
                cu = T.front() * divcs->q[1];
            }
            else {
                cu = T.front();
            }
        }
        else {
            std::vector<std::complex<double>> ctxs(dc);
            std::vector<std::complex<double>> weights(dc);

            cu = 0;
            for (uint32_t i = 0; i < dc; ++i) {
                ctxs[i]    = T[i];
                weights[i] = divcs->q[i + 1];
                cu += ctxs[i] * weights[i];
            }
        }

        // adds the free term (at x^0)
        cu += divcs->q.front() / 2.0;
        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    std::complex<double> qu;

    if (Degree(divqr->q) > k) {
        qu = InnerEvalChebyshevPS(x, divqr->q, k, m - 1, T, T2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto qcopy = divqr->q;
        qcopy.resize(k);
        if (Degree(qcopy) > 0) {
            std::vector<std::complex<double>> ctxs(Degree(qcopy));
            std::vector<std::complex<double>> weights(Degree(qcopy));

            qu = 0;

            for (uint32_t i = 0; i < Degree(qcopy); i++) {
                ctxs[i]    = T[i];
                weights[i] = divqr->q[i + 1];
                qu += ctxs[i] * weights[i];
            }
            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            std::complex<double> sum = T[k - 1];
            uint32_t limit           = log2(ToReal(divqr->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = sum + sum;
            }
            qu += sum;
        }
        else {
            std::complex<double> sum = T[k - 1];
            uint32_t limit           = log2(ToReal(divqr->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = sum + sum;
            }
            qu = sum;
        }

        // adds the free term (at x^0)
        qu += divqr->q.front() / 2.0;
        // The number of levels of qu is the same as the number of levels of T[k-1] or T[k-1] + 1.
        // No need to reduce it to T2[m-1] because it only reaches here when m = 2.
    }

    std::complex<double> su;

    if (Degree(s2) > k) {
        su = InnerEvalChebyshevPS(x, s2, k, m - 1, T, T2);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto scopy = s2;
        scopy.resize(k);
        if (Degree(scopy) > 0) {
            std::vector<std::complex<double>> ctxs(Degree(scopy));
            std::vector<std::complex<double>> weights(Degree(scopy));

            su = 0;

            for (uint32_t i = 0; i < Degree(scopy); i++) {
                ctxs[i]    = T[i];
                weights[i] = s2[i + 1];
                su += ctxs[i] * weights[i];
            }

            // the highest order coefficient will always be 1 because s2 is monic.
            su += T[k - 1];
        }
        else {
            su = T[k - 1];
        }

        // adds the free term (at x^0)
        su += s2.front() / 2.0;
    }

    std::complex<double> result;

    if (flag_c) {
        result = T2[m - 1] + cu;
    }
    else {
        result = T2[m - 1] + divcs->q.front() / 2.0;
    }

    result = result * qu;
    result += su;

    return result;
}

std::complex<double> evalChebyshevSeriesPSWithPrecompCustom(seriesPowersCustom ctxtPolys,
                                                            const std::vector<std::complex<double>>& coefficients) {
    auto f2 = coefficients;
    auto n  = Degree(f2);
    f2.resize(n + 1);

    auto T     = ctxtPolys.powers;
    auto T2    = ctxtPolys.powers2;
    auto T2km1 = ctxtPolys.power2km1;
    unsigned k = ctxtPolys.k;
    unsigned m = ctxtPolys.m;

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Add T^{k(2^m - 1)}(y) to the polynomial that has to be evaluated
    f2.resize(2 * k2m2k + k + 1, 0.0);
    f2.back() = 1;

    // Divide f2 by T^{k*2^{m-1}}
    std::vector<std::complex<double>> Tkm(k2m2k + k + 1);
    Tkm.back() = 1;
    auto divqr = LongDivisionChebyshev(f2, Tkm);

    // Subtract x^{k(2^{m-1} - 1)} from r
    auto r2 = divqr->r;
    if (static_cast<int32_t>(k2m2k - Degree(r2)) <= 0) {
        r2[static_cast<int32_t>(k2m2k)] -= 1;
        r2.resize(Degree(r2) + 1);
    }
    else {
        r2.resize(static_cast<int32_t>(k2m2k + 1));
        r2.back() = -1;
    }

    // Divide r2 by q
    auto divcs = LongDivisionChebyshev(r2, divqr->q);

    // Add x^{k(2^{m-1} - 1)} to s
    auto s2 = divcs->r;
    s2.resize(k2m2k + 1);
    s2.back() = 1;

    // Evaluate c at u
    std::complex<double> cu;
    uint32_t dc = Degree(divcs->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            if (IsNotEqualOne(divcs->q[1])) {
                cu = T.front() * divcs->q[1];
            }
            else {
                cu = T.front();
            }
        }
        else {
            std::vector<std::complex<double>> ctxs(dc);
            std::vector<std::complex<double>> weights(dc);

            cu = 0;

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i]    = T[i];
                weights[i] = divcs->q[i + 1];
                cu += ctxs[i] * weights[i];
            }
        }

        // adds the free term (at x^0)
        cu += divcs->q.front() / 2.0;
        // TODO : Andrey why not T2[m-1]->GetLevel() instead?
        // Need to reduce levels to the level of T2[m-1].
        //    uint32_t levelDiff = y->GetLevel() - cu->GetLevel() + ceil(log2(k)) + m - 1;
        //    cc->LevelReduceInPlace(cu, nullptr, levelDiff);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    std::complex<double> qu;

    if (Degree(divqr->q) > k) {
        qu = InnerEvalChebyshevPS(T[0], divqr->q, k, m - 1, T, T2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto qcopy = divqr->q;
        qcopy.resize(k);
        if (Degree(qcopy) > 0) {
            std::vector<std::complex<double>> ctxs(Degree(qcopy));
            std::vector<std::complex<double>> weights(Degree(qcopy));

            qu = 0;

            for (uint32_t i = 0; i < Degree(qcopy); ++i) {
                ctxs[i]    = T[i];
                weights[i] = divqr->q[i + 1];
                qu += ctxs[i] * weights[i];
            }
            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            std::complex<double> sum = T[k - 1];
            uint32_t limit           = log2(ToReal(divqr->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = sum + sum;
            }
            qu += sum;
        }
        else {
            std::complex<double> sum = T[k - 1];
            uint32_t limit           = log2(ToReal(divqr->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = sum + sum;
            }
            qu = sum;
        }

        // adds the free term (at x^0)
        qu += divqr->q.front() / 2.0;
        // The number of levels of qu is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so the number of levels of qu and T2[m-1] will be the same.
    }

    std::complex<double> su;

    if (Degree(s2) > k) {
        su = InnerEvalChebyshevPS(T[0], s2, k, m - 1, T, T2);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto scopy = s2;
        scopy.resize(k);
        if (Degree(scopy) > 0) {
            std::vector<std::complex<double>> ctxs(Degree(scopy));
            std::vector<std::complex<double>> weights(Degree(scopy));

            su = 0;

            for (uint32_t i = 0; i < Degree(scopy); ++i) {
                ctxs[i]    = T[i];
                weights[i] = s2[i + 1];
                su += ctxs[i] * weights[i];
            }

            // the highest order coefficient will always be 1 because s2 is monic.
            su += T[k - 1];
        }
        else {
            su = T[k - 1];
        }

        // adds the free term (at x^0)
        su += s2.front() / 2.0;
        // The number of levels of su is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so need to reduce the number of levels by 1.
    }

    // TODO : Andrey : here is different from 895 line
    // Reduce number of levels of su to number of levels of T2km1.
    //  cc->LevelReduceInPlace(su, nullptr);

    std::complex<double> result;

    if (flag_c) {
        result = T2[m - 1] + cu;
    }
    else {
        result = T2[m - 1] + divcs->q.front() / 2.0;
    }

    result = result * qu;

    result += su;
    result -= T2km1;

    return result;
}

// Keep the example-side cleartext/debug path synchronized with the library tables.
static const inline std::vector<std::complex<double>>& coeff_exp_16_double_46 =
    lbcrypto::FHECKKSRNS::GetCoeffExp16Double46();
static const inline std::vector<std::complex<double>>& coeff_exp_16_double_58 =
    lbcrypto::FHECKKSRNS::GetCoeffExp16Double58();
static const inline std::vector<std::complex<double>>& coeff_exp_25_double_58 =
    lbcrypto::FHECKKSRNS::GetCoeffExp25Double58();
static const inline std::vector<std::complex<double>>& coeff_exp_25_double_66 =
    lbcrypto::FHECKKSRNS::GetCoeffExp25Double66();

// Chebyshev series coefficients for the SPARSE case (degree 44)
static const inline std::vector<std::complex<double>> g_coefficientsSparse = [] {
    const auto& coeffs = lbcrypto::FHECKKSRNS::GetGCoefficientsSparse();
    return std::vector<std::complex<double>>(coeffs.begin(), coeffs.end());
}();

// Chebyshev series coefficients for the SPARSE ENCAPSULATED case (degree 32)
static const inline std::vector<std::complex<double>> g_coefficientsSparseEncapsulated = [] {
    const auto& coeffs = lbcrypto::FHECKKSRNS::GetGCoefficientsSparseEncapsulated();
    return std::vector<std::complex<double>>(coeffs.begin(), coeffs.end());
}();

const uint32_t R_SPARSE = 3;

void ApplyDoubleAngleIterations(std::complex<double>& input, uint32_t numIter) {
    constexpr double twoPi = 2.0 * M_PI;

    for (int32_t i = 1 - numIter; i <= 0; ++i) {
        double scalar = -std::pow(twoPi, -std::pow(2.0, i));
        input *= input;
        input += input + scalar;
    }
}

// void disturb() {
//                 if (i == 0) {
//                     auto t            = c2s[0];
//                     auto tReal        = t.real();
//                     auto tRealRounded = std::complex<double>(std::round(tReal * p_in_global) / p_in_global, 0) +
//                                         std::complex<double>(std::exp2(-30), std::exp2(-30));
//                     std::cout << msg << std::setprecision(20) << " c2s[0]: " << t << " Rounded: " << tRealRounded
//                               << std::endl;
//                     auto expVal = std::exp(std::complex<double>(0, 2 * M_PI) * tRealRounded);
//                     std::cout << msg << std::setprecision(20) << " exp tRealRounded: " << expVal << std::endl;
//
//                     std::complex<double> sum = 0;
//                     auto power               = std::complex(1.0, 0.0);
//                     // do a poly eval using coeffcomp
//                     for (size_t j = 0; j != coeffcomp_global.size(); ++j) {
//                         sum += coeffcomp_global[j] * power;
//                         power *= expVal;
//                     }
//                     sum *= 4;  // related to scaleTHI
//                     std::cout << msg << std::setprecision(20) << " evalPolyNaive at exp tRealRounded: " << sum
//                               << std::endl;
//
//                     auto evalPoly =
//                         evalPolyPSWithPrecompCustom(evalPowersPS(expVal, coeffcomp_global), coeffcomp_global);
//                     evalPoly *= 4;
//                     std::cout << msg << std::setprecision(20) << " evalPoly at exp tRealRounded: " << evalPoly
//                               << std::endl;
//                     auto noise = extract_noise(evalPoly.real(), p_out_global);
//                     std::cout << msg << std::setprecision(20) << " noise at real: 2^" << noise << std::endl;
//                 }
// }

// void unused() {
//     // print the a
//     const DCRTPoly& s = keyPair.secretKey->GetPrivateElement();
//     size_t sizeQ      = s.GetParams()->GetParams().size();
//
//     size_t sizeQl = cv[0].GetParams()->GetParams().size();
//
//     size_t diffQl = sizeQ - sizeQl;
//
//     auto scopy(s);
//     scopy.DropLastElements(diffQl);
//
//     DCRTPoly sPower(scopy);
//
//     //DCRTPoly b(cv[0], 0);
//     //b.SetFormat(Format::EVALUATION);
//
//     const auto& q{cv[0].GetParams()->GetModulus()};
//     const auto& half{q >> 1};
//     const auto qDividedBy6 = q / 6;
//     const auto qDividedBy4 = q / 4;
//
//     std::cout << " q : " << q << std::endl;
//
//     auto res(cv[0]);
//     res.SetFormat(Format::COEFFICIENT);
//     res.SetValuesToZero();
//     DCRTPoly ci;
//
//     for (size_t i = 0; i < cv.size(); i++) {
//         ci = cv[i];
//         ci.SetFormat(Format::COEFFICIENT);
//         auto bigCi      = ci.CRTInterpolate();
//         auto& bigCoeffs = bigCi.GetValues();
//
//         auto ciCopy(cv[i]);
//         ciCopy.SetFormat(Format::COEFFICIENT);
//         ciCopy.SetValuesToZero();
//         for (size_t i = 0; i != ci.GetRingDimension(); ++i) {
//             auto x   = bigCoeffs[i];
//             bool neg = false;
//             if (x > half) {
//                 x   = q - x;
//                 neg = true;
//             }
//             BigInteger flag = 0;
//             if (x < qDividedBy4)
//                 if (neg)
//                     flag = q - BigInteger(1);
//                 else
//                     flag = 1;
//             else if (x > qDividedBy4 && !neg)
//                 flag = 3;
//             else
//                 flag = q - BigInteger(3);
//             ciCopy.GetAllElements()[0][i] = BigInteger(flag);
//         }
//
//         ciCopy.SetFormat(Format::EVALUATION);
//
//         if (i == 1)
//             res += sPower * ciCopy;
//         else
//             res = ciCopy;
//         //sPower *= scopy;
//     }
//     res.SetFormat(Format::COEFFICIENT);
//     for (size_t j = 0; j != 16; ++j) {
//         auto y   = res.GetAllElements()[0][j * (4096 / 16)];
//         bool neg = false;
//         if (y > half.ConvertToInt()) {
//             y   = q.ConvertToInt() - y;
//             neg = true;
//         }
//         std::cout << "res " << j << " " << (neg ? "-" : "") << y.ConvertToDouble() / 8 << " " << std::endl;
//     }
// }

// A cleartext version of CryptoContext<...>::EvalChebyshevFunction(...)
std::complex<double> customEvalChebyshevFunctionPtxt(std::vector<std::complex<double>> coeffs, std::complex<double> x) {
    // The standard practice is to halve the 1st coefficient.
    // See, for example, the Chebyshev Series section at
    // https://www.cfm.brown.edu/people/dobrush/am34/Mathematica/ch5/chebyshev.html
    // and derivation of Eq. (6) in https://arxiv.org/pdf/1810.04282.
    // The halving requirement follows from the discrete orthogonality relation for Chebyshev polynomials,
    // i.e., Eq. (4) in https://arxiv.org/pdf/1810.04282.
    coeffs[0] /= 2.0;

    auto x2 = 2.0 * x;

    std::complex<double> t_prev = 1.0;  // T0(x) = 1
    auto t_j                    = x;    // T1(x) = x
    auto y                      = coeffs[0] + coeffs[1] * x;
    // Use the recursive formula T_{i+1}(X) = 2x T_i(x) - T_{i-1}(x)
    for (size_t j = 2; j < coeffs.size(); j++) {
        // Compute T_j(x) and add it to the approximation
        auto t_next = x2 * t_j - t_prev;
        t_prev      = t_j;
        t_j         = t_next;
        y += coeffs[j] * t_next;
    }
    return y;
}

// Stable Chebyshev series evaluation using Clenshaw algorithm.
// Computes: f(x) = sum_{k=0}^{n-1} c[k] * T_k(x)
// where c[0] is the *unhalved* Chebyshev coefficient (i.e., standard convention where
// the series is (c0/2) + sum_{k>=1} c_k T_k(x)).
std::complex<double> evalChebyshevClenshaw(const std::vector<std::complex<double>>& c, const std::complex<double>& x) {
    const size_t n = c.size();
    if (n == 0)
        return {0.0, 0.0};
    if (n == 1)
        return c[0] / 2.0;

    std::complex<double> b_kp1 = 0.0;  // b_{k+1}
    std::complex<double> b_kp2 = 0.0;  // b_{k+2}

    // for k = n-1 down to 1:
    // b_k = 2x b_{k+1} - b_{k+2} + c_k
    const std::complex<double> two_x = 2.0 * x;
    for (size_t kk = n - 1; kk >= 1; --kk) {
        std::complex<double> b_k = two_x * b_kp1 - b_kp2 + c[kk];
        b_kp2                    = b_kp1;
        b_kp1                    = b_k;
        if (kk == 1)
            break;  // avoid size_t underflow
    }

    // result = x b_1 - b_2 + c0/2
    return x * b_kp1 - b_kp2 + c[0] / 2.0;
}
