#pragma once

#include "math/hermite.h"
#include "openfhe.h"
#include "schemelet/rlwe-mp.h"
#include "fbt-test-utils.h"
#include "fbt-benchmark-utils.h"

#include <algorithm>
#include <chrono>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <ratio>
#include <string>
#include <vector>

// #define TRACE
#ifndef TRACE
    #define BENCHMARK
#endif

using namespace lbcrypto;

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------
using time_point_t =
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<long long, std::ratio<1, 1000000000>>>;

// ---------------------------------------------------------------------------
// Structs
// ---------------------------------------------------------------------------

// Plaintext values captured at each FBT stage, used for precision analysis.
struct PtValues {
    std::vector<double> input;
    std::vector<double> modraise;
    std::vector<std::complex<double>> c2s;
    std::vector<std::complex<double>> evalExp;
    std::vector<std::complex<double>> lut;
    std::vector<double> output;
};

// Timing state — always active.
struct TimingState {
    std::map<std::string, time_point_t> times;
    uint32_t keySwitchCounterAtEvalExp = 0;
    uint32_t keySwitchCounterAtLUT     = 0;
};

// Debug context — only meaningful when BENCHMARK is not defined.
struct DebugContext {
    CryptoContextT cc;
    PrivateKeyT sk;
    BigInteger pInput;
    SecretKeyDist skDist;
    bool complexLUT = false;
    size_t order    = 1;
    size_t slots    = 0;
    size_t pIn      = 0;
    size_t pOut     = 0;
    std::vector<std::complex<double>> coeffComp;
    PtValues precise;
    PtValues encrypted;
};

struct ArbitraryLUTRunRequest {
    int runNumber;
    uint32_t addNoiseBase;
    bool warmup;
};

struct ArbitraryLUTRunResult {
    int runNumber;
    uint32_t addNoiseBase;
    bool warmup;
    uint32_t maxError;
    double overallMs;
    double lutMs;
    uint32_t evalExpKeySwitchCount;
    uint32_t lutKeySwitchCount;
    PrecisionState precision;
};

// ---------------------------------------------------------------------------
// Globals (defined in fbt-benchmark.cpp)
// ---------------------------------------------------------------------------
extern const BigInteger QBFVINIT;
extern const BigInteger QBFVINITLARGE;

// Referenced externally by ckksrns-fhe.cpp — do not rename.
extern lbcrypto::InterpolationMethod __interpolation_method_global;
extern uint32_t __eval_exp_degree_global;

extern TimingState g_timing;
extern PrecisionState g_precision;
extern RunMode g_runMode;
extern DebugContext g_debug;

inline double extract_noise(double val, size_t p) {
    return val - std::round(val * p) / p;
}

inline std::complex<double> extract_noise_complex(std::complex<double> val, size_t p) {
    return val - std::complex<double>(std::round(val.real() * p) / p, std::round(val.imag() * p) / p);
}

inline double extract_precise(double val, size_t p) {
    return std::round(val * p) / p;
}

inline std::complex<double> extract_precise_complex(std::complex<double> val, size_t p) {
    return {std::round(val.real() * p) / p, std::round(val.imag() * p) / p};
}

inline std::vector<double> get_precise_array(const std::vector<double>& vals, size_t p) {
    std::vector<double> out(vals.size());
    std::transform(vals.begin(), vals.end(), out.begin(), [p](double v) { return extract_precise(v, p); });
    return out;
}

inline std::vector<std::complex<double>> get_precise_array_complex(const std::vector<std::complex<double>>& vals,
                                                                   size_t p) {
    std::vector<std::complex<double>> out(vals.size());
    std::transform(vals.begin(), vals.end(), out.begin(), [p](auto v) { return extract_precise_complex(v, p); });
    return out;
}

inline double get_precision_array(const std::vector<double>& vals, size_t p) {
    double worst = -std::numeric_limits<double>::infinity();
    for (auto v : vals)
        worst = std::max(worst, std::log2(std::abs(extract_noise(v, p))));
    return worst;
}

inline double get_precision_array_complex(const std::vector<std::complex<double>>& vals, size_t p) {
    double worst = -std::numeric_limits<double>::infinity();
    for (auto v : vals)
        worst = std::max(worst, std::log2(std::abs(extract_noise_complex(v, p))));
    return worst;
}

inline std::vector<std::complex<double>> complex_array_diff(const std::vector<std::complex<double>>& a,
                                                            const std::vector<std::complex<double>>& b) {
    std::vector<std::complex<double>> out(a.size());
    std::transform(a.begin(), a.end(), b.begin(), out.begin(), std::minus<std::complex<double>>{});
    return out;
}

inline std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>> complex_exp_array_diff(
    const std::vector<std::complex<double>>& base, const std::vector<std::complex<double>>& off) {
    std::vector<std::complex<double>> inCircle, offCircle;
    for (size_t i = 0; i < base.size(); ++i) {
        auto diff       = off[i] - base[i];
        auto projLength = std::real(diff * std::conj(base[i])) / std::norm(base[i]);
        auto proj       = projLength * base[i];
        inCircle.push_back(diff - proj);
        offCircle.push_back(proj);
    }
    return {inCircle, offCircle};
}

const std::array<uint8_t, 256> aesSBox = {
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76, 0xCA, 0x82, 0xC9,
    0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0, 0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F,
    0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15, 0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07,
    0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75, 0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3,
    0x29, 0xE3, 0x2F, 0x84, 0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58,
    0xCF, 0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8, 0x51, 0xA3,
    0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2, 0xCD, 0x0C, 0x13, 0xEC, 0x5F,
    0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73, 0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88,
    0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB, 0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC,
    0x62, 0x91, 0x95, 0xE4, 0x79, 0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A,
    0xAE, 0x08, 0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A, 0x70,
    0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E, 0xE1, 0xF8, 0x98, 0x11,
    0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF, 0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42,
    0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16};

// ---------------------------------------------------------------------------
// Function declarations
// ---------------------------------------------------------------------------
void print_times();

double __debug(CiphertextT ct, std::string msg);

std::vector<ArbitraryLUTRunResult> ArbitraryLUT(
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput, BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, bool allInputs, NoiseInjectionStage noiseInjectionStage,
    lbcrypto::ScalingTechnique scalingTechnique, lbcrypto::SecretKeyDist secretKeyDist,
    const std::vector<uint32_t>& levelBudget, std::function<int64_t(int64_t)> func,
    const std::vector<ArbitraryLUTRunRequest>& runRequests, InterpolationMethod method = HERMITE_AKP25,
    std::function<std::complex<double>(int64_t)> complexFunc = nullptr, bool complexLUT = false);
