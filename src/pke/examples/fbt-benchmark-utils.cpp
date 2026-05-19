#include "fbt-benchmark-utils.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace {

const std::vector<lbcrypto::InterpolationMethod> kDefaultMethods = {
    lbcrypto::HERMITE_AKP25,
    lbcrypto::HERMITE_SPARSE_THI,
    lbcrypto::HERMITE_BKSS24,
};

lbcrypto::InterpolationMethod ParseMethodName(const std::string& value) {
    static const std::map<std::string, lbcrypto::InterpolationMethod> kMethodMap = {
        {"AKP25", lbcrypto::HERMITE_AKP25},           {"MIXEDCONSTRAINTS", lbcrypto::HERMITE_SPARSE_THI},
        {"HYBRID", lbcrypto::HERMITE_SPARSE_THI},     {"SPARSETHI", lbcrypto::HERMITE_SPARSE_THI},
        {"SPARSE-THI", lbcrypto::HERMITE_SPARSE_THI}, {"SPARSE_THI", lbcrypto::HERMITE_SPARSE_THI},
        {"BKSS24", lbcrypto::HERMITE_BKSS24},         {"FULLCOMPLEX", lbcrypto::HERMITE_FULL_THI},
        {"FULL_COMPLEX", lbcrypto::HERMITE_FULL_THI}, {"FULL_THI", lbcrypto::HERMITE_FULL_THI},
        {"FULLTHI", lbcrypto::HERMITE_FULL_THI},
    };
    auto it = kMethodMap.find(value);
    if (it == kMethodMap.end())
        throw std::invalid_argument("unsupported method: " + value);
    return it->second;
}

std::vector<uint32_t> ParseNoiseList(const std::string& value) {
    std::vector<uint32_t> out;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty())
            continue;
        out.push_back(static_cast<uint32_t>(std::stoul(item)));
    }
    if (out.empty())
        throw std::invalid_argument("empty noise-base list");
    return out;
}

std::vector<uint32_t> ParseLevelBudget(const std::string& value) {
    std::vector<uint32_t> out;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty())
            continue;
        out.push_back(static_cast<uint32_t>(std::stoul(item)));
    }
    if (out.size() != 2)
        throw std::invalid_argument("level-budget must have exactly two comma-separated integers");
    return out;
}

std::vector<lbcrypto::InterpolationMethod> ParseMethodList(const std::string& value) {
    std::vector<lbcrypto::InterpolationMethod> out;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty())
            continue;
        out.push_back(ParseMethodName(item));
    }
    if (out.empty())
        throw std::invalid_argument("empty methods list");
    return out;
}

NoiseInjectionStage ParseNoiseInjectionStageName(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (value == "INPUT" || value == "INPUT_SLOTS" || value == "INPUT-SLOTS")
        return NoiseInjectionStage::InputSlots;
    if (value == "PRE_EVALEXP" || value == "PRE-EVALEXP" || value == "PREEVALEXP")
        return NoiseInjectionStage::PreEvalExp;
    throw std::invalid_argument("unsupported noise injection stage: " + value);
}

lbcrypto::ScalingTechnique ParseScalingTechniqueName(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (value == "FIXEDMANUAL")
        return lbcrypto::FIXEDMANUAL;
    if (value == "FLEXIBLEMANUAL")
        return lbcrypto::FLEXIBLEMANUAL;
    throw std::invalid_argument("unsupported scaling technique: " + value);
}

lbcrypto::SecretKeyDist ParseSecretKeyDistName(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (value == "SPARSE_TERNARY" || value == "SPARSE-TERNARY")
        return lbcrypto::SPARSE_TERNARY;
    if (value == "SPARSE_ENCAPSULATED" || value == "SPARSE-ENCAPSULATED")
        return lbcrypto::SPARSE_ENCAPSULATED;
    throw std::invalid_argument("unsupported key distribution: " + value);
}

std::string NoiseBasesString(const std::vector<uint32_t>& values) {
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0)
            os << ",";
        os << values[i];
    }
    return os.str();
}

}  // namespace

BenchmarkOptions ParseBenchmarkOptions(int argc, char* argv[]) {
    BenchmarkOptions options;
    options.methods      = kDefaultMethods;
    bool noiseOverride   = false;
    bool methodsOverride = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mode=precision" || arg == "--precision")
            options.runMode = RunMode::Precision;
        else if (arg == "--mode=benchmark" || arg == "--benchmark")
            options.runMode = RunMode::Benchmark;
        else if (arg == "--mode=verify" || arg == "--verify")
            options.runMode = RunMode::Verify;
        else if (arg == "--csv")
            options.emitCsv = true;
        else if (arg == "--all-inputs")
            options.allInputs = true;
        else if (arg == "--complex-lut")
            options.complexLUT = true;
        else if (arg == "--unsafe")
            options.unsafe = true;
        else if (arg.rfind("--ring-dim=", 0) == 0)
            options.ringDim = static_cast<uint32_t>(std::stoul(arg.substr(11)));
        else if (arg.rfind("--sf=", 0) == 0)
            options.scalingFactor = static_cast<uint32_t>(std::stoul(arg.substr(5)));
        else if (arg.rfind("--slots=", 0) == 0)
            options.slots = static_cast<uint32_t>(std::stoul(arg.substr(8)));
        else if (arg.rfind("--eval-exp-degree=", 0) == 0)
            options.evalExpDegree = static_cast<uint32_t>(std::stoul(arg.substr(18)));
        else if (arg.rfind("--order=", 0) == 0)
            options.order = static_cast<uint32_t>(std::stoul(arg.substr(8)));
        else if (arg.rfind("--methods=", 0) == 0) {
            options.methods = ParseMethodList(arg.substr(10));
            methodsOverride = true;
        }
        else if (arg.rfind("--p-input=", 0) == 0)
            options.pInput = static_cast<uint32_t>(std::stoul(arg.substr(10)));
        else if (arg.rfind("--p-output=", 0) == 0)
            options.pOutput = static_cast<uint32_t>(std::stoul(arg.substr(11)));
        else if (arg.rfind("--noise-stage=", 0) == 0)
            options.noiseInjectionStage = ParseNoiseInjectionStageName(arg.substr(14));
        else if (arg.rfind("--scaling-technique=", 0) == 0)
            options.scalingTechnique = ParseScalingTechniqueName(arg.substr(20));
        else if (arg.rfind("--key-dist=", 0) == 0)
            options.secretKeyDist = ParseSecretKeyDistName(arg.substr(11));
        else if (arg.rfind("--sk-dist=", 0) == 0)
            options.secretKeyDist = ParseSecretKeyDistName(arg.substr(10));
        else if (arg.rfind("--secret-key-dist=", 0) == 0)
            options.secretKeyDist = ParseSecretKeyDistName(arg.substr(18));
        else if (arg.rfind("--level-budget=", 0) == 0)
            options.levelBudget = ParseLevelBudget(arg.substr(15));
        else if (arg.rfind("--noise-base=", 0) == 0) {
            options.addNoiseBases = ParseNoiseList(arg.substr(13));
            noiseOverride         = true;
        }
        else {
            throw std::invalid_argument("Usage: " + std::string(argv[0]) +
                                        " [--mode=benchmark|precision|verify] [--csv] "
                                        "[--all-inputs] [--complex-lut] "
                                        "[--ring-dim=N] [--sf=N] [--slots=N] [--eval-exp-degree=N] [--order=N] "
                                        "[--methods=AKP25,MIXEDCONSTRAINTS,BKSS24] "
                                        "[--p-input=N] [--p-output=N] "
                                        "[--noise-stage=INPUT_SLOTS|PRE_EVALEXP] "
                                        "[--scaling-technique=FIXEDMANUAL|FLEXIBLEMANUAL] "
                                        "[--key-dist=SPARSE_TERNARY|SPARSE_ENCAPSULATED] "
                                        "[--level-budget=A,B] "
                                        "[--noise-base=N[,M...]]");
        }
    }

    if (!noiseOverride) {
        options.addNoiseBases = (options.runMode == RunMode::Precision) ?
                                    std::vector<uint32_t>{10, 13, 16, 19, 22, 25, 28, 31, 34} :
                                    std::vector<uint32_t>{10};
    }
    if (options.complexLUT && !methodsOverride) {
        options.methods = {lbcrypto::HERMITE_AKP25, lbcrypto::HERMITE_SPARSE_THI};
    }
    return options;
}

const char* RunModeName(RunMode mode) {
    switch (mode) {
        case RunMode::Benchmark:
            return "benchmark";
        case RunMode::Precision:
            return "precision";
        case RunMode::Verify:
            return "verify";
    }
    return "unknown";
}

const char* NoiseInjectionStageName(NoiseInjectionStage stage) {
    switch (stage) {
        case NoiseInjectionStage::InputSlots:
            return "INPUT_SLOTS";
        case NoiseInjectionStage::PreEvalExp:
            return "PRE_EVALEXP";
    }
    return "UNKNOWN";
}

const char* ScalingTechniqueName(lbcrypto::ScalingTechnique technique) {
    switch (technique) {
        case lbcrypto::FIXEDMANUAL:
            return "FIXEDMANUAL";
        case lbcrypto::FLEXIBLEMANUAL:
            return "FLEXIBLEMANUAL";
        default:
            return "UNKNOWN";
    }
}

const char* SecretKeyDistName(lbcrypto::SecretKeyDist keyDist) {
    switch (keyDist) {
        case lbcrypto::SPARSE_TERNARY:
            return "SPARSE_TERNARY";
        case lbcrypto::SPARSE_ENCAPSULATED:
            return "SPARSE_ENCAPSULATED";
        default:
            return "UNKNOWN";
    }
}

void PrintBenchmarkHeader(const BenchmarkOptions& options, lbcrypto::SecretKeyDist keyDist, uint32_t slots,
                          uint32_t mulDepth) {
    std::cout << "Mode=" << RunModeName(options.runMode) << " ring_dim=" << options.ringDim
              << " key_dist=" << SecretKeyDistName(keyDist) << " slots=" << slots << " sf=" << options.scalingFactor
              << " eval_exp_degree=" << options.evalExpDegree << " order=" << options.order
              << " noise_stage=" << NoiseInjectionStageName(options.noiseInjectionStage)
              << " scaling=" << ScalingTechniqueName(options.scalingTechnique) << " mul_depth=" << mulDepth
              << " p_input=" << options.pInput << " p_output=" << options.pOutput
              << " all_inputs=" << (options.allInputs ? 1 : 0) << " complex_lut=" << (options.complexLUT ? 1 : 0)
              << " level_budget=" << options.levelBudget[0] << "," << options.levelBudget[1] << "\n";
    if (options.runMode == RunMode::Benchmark) {
        std::cout << "Measured runs=5 warmup_runs=1 noise_base=" << NoiseBasesString(options.addNoiseBases) << "\n";
    }
    else if (options.runMode == RunMode::Precision) {
        std::cout << "Noise bases=" << NoiseBasesString(options.addNoiseBases) << "\n";
    }
    else {
        std::cout << "Noise base=" << NoiseBasesString(options.addNoiseBases) << "\n";
    }
    std::cout << "\n";
}

void PrintPrecisionCsv(const std::vector<RunResult>& runResults) {
    std::cout << "CSV_BEGIN\n";
    std::cout
        << "lut,method,run,noise_base,max_error,input_precision,lut_precision,lut_margin,output_precision,predicted_output_precision,prediction_gap,output_margin,eval_exp_key_switch_count,key_switch_count\n";
    for (const auto& r : runResults) {
        std::cout << r.lutName << "," << r.method << "," << r.runNumber << "," << r.addNoiseBase << "," << r.maxError
                  << "," << r.inputPrecision << "," << r.lutPrecision << "," << r.lutMarginBits << ","
                  << r.outputPrecision << "," << r.predictedOutputPrecision << "," << r.predictionGapBits << ","
                  << r.outputMarginBits << "," << r.evalExpKeySwitchCount << "," << r.lutKeySwitchCount << "\n";
    }
    std::cout << "CSV_END\n";
}

void PrintVerifyLine(const std::string& lutName, lbcrypto::InterpolationMethod method, uint32_t maxError) {
    std::cout << std::left << std::setw(12) << lutName << std::left << std::setw(28) << method << " "
              << (maxError == 0 ? "PASS" : "FAIL") << " max_error=" << maxError << "\n";
}

void PrintVerifySummary(size_t passed, size_t total) {
    std::cout << "\nVerification: " << (passed == total ? "PASS" : "FAIL") << "  (" << passed << "/" << total
              << " correct)\n";
}

void PrintBenchmarkRunLine(const std::string& lutName, lbcrypto::InterpolationMethod method, int runNumber,
                           uint32_t addNoiseBase, uint32_t maxError, double overallMs, double lutMs,
                           uint32_t evalExpKeySwitchCount, uint32_t lutKeySwitchCount) {
    std::cout << std::left << std::setw(12) << lutName << std::left << std::setw(28) << method << " run=" << runNumber
              << " noise=" << addNoiseBase << " " << (maxError == 0 ? "PASS" : "FAIL") << " max_error=" << maxError
              << " overall_ms=" << std::fixed << std::setprecision(2) << overallMs << " lut_ms="
              << lutMs
              //<< " exp_ks=" << evalExpKeySwitchCount
              << " lut_ks=" << lutKeySwitchCount << "\n";
}

void PrintBenchmarkSummaryLine(const MethodResult& result) {
    std::cout << std::left << std::setw(12) << result.lutName << std::left << std::setw(28) << result.method
              << " summary avg_overall_ms=" << std::fixed << std::setprecision(2)
              << result.avgOverallMs
              //<< " max_overall_ms=" << result.maxOverallMs
              << " avg_lut_ms="
              << result.avgLutMs
              // << " max_lut_ms=" << result.maxLutMs
              //<< " exp_ks=" << result.evalExpKeySwitchCount
              << " lut_ks=" << result.lutKeySwitchCount << "\n";
}

void PrintPrecisionRunLine(const std::string& lutName, lbcrypto::InterpolationMethod method, uint32_t addNoiseBase,
                           uint32_t maxError, const PrecisionState& precision, uint32_t evalExpKeySwitchCount,
                           uint32_t lutKeySwitchCount) {
    std::cout << std::left << std::setw(12) << lutName << std::left << std::setw(28) << method
              << " noise=" << addNoiseBase << " " << (maxError == 0 ? "PASS" : "FAIL") << " max_error=" << maxError
              << " in=" << std::fixed << std::setprecision(3) << precision.inputPrecision << " lut="
              << precision.lutPrecision
              //<< " out=" << precision.outputPrecision
              //<< " pred=" << precision.predictedOutputPrecision << " pred_gap=" << precision.predictionGapBits
              //<< " exp_ks=" << evalExpKeySwitchCount
              << " lut_ks=" << lutKeySwitchCount;
    //if (precision.predictionGapBits > 0.1)
    //    std::cout << " pred_status=optimistic";
    std::cout << "\n";
}
