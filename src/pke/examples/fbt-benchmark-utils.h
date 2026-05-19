#pragma once

#include "math/hermite.h"
#include "openfhe.h"

#include <limits>
#include <cstdint>
#include <string>
#include <map>
#include <vector>

struct PrecisionState {
    double inputPrecision           = std::numeric_limits<double>::quiet_NaN();
    double lutPrecision             = std::numeric_limits<double>::quiet_NaN();
    double outputPrecision          = std::numeric_limits<double>::quiet_NaN();
    double predictedOutputPrecision = std::numeric_limits<double>::quiet_NaN();
    double predictionGapBits        = std::numeric_limits<double>::quiet_NaN();
    double lutMarginBits            = std::numeric_limits<double>::quiet_NaN();
    double outputMarginBits         = std::numeric_limits<double>::quiet_NaN();
};

enum class RunMode {
    Benchmark,
    Precision,
    Verify,
};

enum class NoiseInjectionStage {
    InputSlots,
    PreEvalExp,
};

struct BenchmarkOptions {
    RunMode runMode                             = RunMode::Benchmark;
    bool emitCsv                                = false;
    bool allInputs                              = false;
    bool complexLUT                             = false;
    bool unsafe                                 = false;
    uint32_t scalingFactor                      = 50;
    uint32_t ringDim                            = 1u << 8;
    uint32_t slots                              = 32;
    uint32_t evalExpDegree                      = 58;
    uint32_t order                              = 1;
    uint32_t pInput                             = 16;
    uint32_t pOutput                            = 16;
    NoiseInjectionStage noiseInjectionStage     = NoiseInjectionStage::PreEvalExp;
    lbcrypto::ScalingTechnique scalingTechnique = lbcrypto::FIXEDMANUAL;
    lbcrypto::SecretKeyDist secretKeyDist       = lbcrypto::SPARSE_TERNARY;
    std::vector<uint32_t> levelBudget           = {1, 1};
    std::vector<uint32_t> addNoiseBases;
    std::vector<lbcrypto::InterpolationMethod> methods;
};

struct MethodResult {
    std::string lutName;
    lbcrypto::InterpolationMethod method;
    double avgOverallMs;
    double maxOverallMs;
    double avgLutMs;
    double maxLutMs;
    uint32_t evalExpKeySwitchCount;
    uint32_t lutKeySwitchCount;
};

struct RunResult {
    std::string lutName;
    lbcrypto::InterpolationMethod method;
    int runNumber;
    uint32_t addNoiseBase;
    uint32_t maxError;
    double lutMs;
    uint32_t evalExpKeySwitchCount;
    uint32_t lutKeySwitchCount;
    double inputPrecision;
    double lutPrecision;
    double lutMarginBits;
    double outputPrecision;
    double predictedOutputPrecision;
    double predictionGapBits;
    double outputMarginBits;
};

BenchmarkOptions ParseBenchmarkOptions(int argc, char* argv[]);

const char* RunModeName(RunMode mode);
const char* NoiseInjectionStageName(NoiseInjectionStage stage);

const char* ScalingTechniqueName(lbcrypto::ScalingTechnique technique);

const char* SecretKeyDistName(lbcrypto::SecretKeyDist keyDist);

void PrintBenchmarkHeader(const BenchmarkOptions& options, lbcrypto::SecretKeyDist keyDist, uint32_t slots,
                          uint32_t mulDepth);

void PrintPrecisionCsv(const std::vector<RunResult>& runResults);

void PrintVerifyLine(const std::string& lutName, lbcrypto::InterpolationMethod method, uint32_t maxError);

void PrintVerifySummary(size_t passed, size_t total);

void PrintBenchmarkRunLine(const std::string& lutName, lbcrypto::InterpolationMethod method, int runNumber,
                           uint32_t addNoiseBase, uint32_t maxError, double overallMs, double lutMs,
                           uint32_t evalExpKeySwitchCount, uint32_t lutKeySwitchCount);

void PrintBenchmarkSummaryLine(const MethodResult& result);

void PrintPrecisionRunLine(const std::string& lutName, lbcrypto::InterpolationMethod method, uint32_t addNoiseBase,
                           uint32_t maxError, const PrecisionState& precision, uint32_t evalExpKeySwitchCount,
                           uint32_t lutKeySwitchCount);
