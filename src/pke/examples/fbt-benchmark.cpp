#include "fbt-benchmark.h"
#include "fbt-benchmark-utils.h"

#include <array>
#include <iterator>
#include <numeric>
#include <random>

// ---------------------------------------------------------------------------
// Global definitions
// ---------------------------------------------------------------------------
const BigInteger QBFVINIT(BigInteger(1) << 60);
const BigInteger QBFVINITLARGE(BigInteger(1) << 80);

// Referenced externally by ckksrns-fhe.cpp — do not rename.
lbcrypto::InterpolationMethod __interpolation_method_global = lbcrypto::HERMITE_INVALID;
uint32_t __eval_exp_degree_global                           = 58;
bool __complex_lut_global                                   = false;
bool __unsafe_global                                        = false;

TimingState g_timing;
PrecisionState g_precision;
RunMode g_runMode = RunMode::Benchmark;
DebugContext g_debug;

namespace {

const time_point_t& GetLutStageTime(const TimingState& timing) {
    auto it = timing.times.find("LUT1");
    if (it != timing.times.end())
        return it->second;
    return timing.times.at("LUT");
}

const time_point_t& GetExpStageTime(const TimingState& timing) {
    auto it = timing.times.find("Exp1");
    if (it != timing.times.end())
        return it->second;
    it = timing.times.find("Exp");
    if (it != timing.times.end())
        return it->second;
    return timing.times.at("Exp");
}

int64_t CenterOutputComponent(int64_t value, uint32_t modulus) {
    return (value > static_cast<int64_t>(modulus / 2)) ? value - static_cast<int64_t>(modulus) : value;
}

uint32_t ComponentError(int64_t actual, int64_t expected, uint32_t modulus) {
    return static_cast<uint32_t>(std::abs(actual - expected)) % modulus;
}

}  // namespace

// ---------------------------------------------------------------------------
// print_times: report benchmark timings after a run
// ---------------------------------------------------------------------------
void print_times() {
    auto elapsed_ms = [](time_point_t a, time_point_t b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    double overallMs = elapsed_ms(g_timing.times["Input"], g_timing.times["Output"]);
    double lutMs     = elapsed_ms(GetExpStageTime(g_timing), GetLutStageTime(g_timing));
    std::cout << "  Overall time: " << overallMs << " ms\n";
    std::cout << "  LUT time: " << lutMs << " ms\n";
    std::cout << "  Number of KeySwitch for EvalExp: " << g_timing.keySwitchCounterAtEvalExp << "\n";
    std::cout << "  Number of KeySwitch for LUT: " << g_timing.keySwitchCounterAtLUT << ", each with "
              << lutMs / g_timing.keySwitchCounterAtLUT << " ms\n";
}

// ---------------------------------------------------------------------------
// __debug: stage checkpoint callback invoked at each FBT step.
//   Always records timing.  In debug (non-BENCHMARK) mode also decrypts
//   the ciphertext and reports precision metrics for each stage.
// ---------------------------------------------------------------------------
double __debug(CiphertextT ct, std::string msg) {
    const bool isC2SStage = (msg.find("CoeffsToSlots") != std::string::npos);

    // Timing — always active.
    if (msg == "Input") {
        g_timing.times.clear();
        g_precision                        = {};
        g_timing.keySwitchCounterAtEvalExp = 0;
        g_timing.keySwitchCounterAtLUT     = 0;
    }
    auto timingMsg            = msg;
    g_timing.times[timingMsg] = std::chrono::steady_clock::now();
    if (msg == "Exp" || msg == "Exp1")
        g_timing.keySwitchCounterAtEvalExp = KeySwitchCounter;
    if (msg == "LUT" || msg == "LUT1")
        g_timing.keySwitchCounterAtLUT = KeySwitchCounter;

    if (g_runMode != RunMode::Precision)
        return 0;

    // -----------------------------------------------------------------------
    // Precision analysis: decrypt and compute stage metrics.
    // -----------------------------------------------------------------------
    auto& dbg = g_debug;

    // Decrypt to coefficient representation.
    auto b = DecryptCore(ct->GetElements(), dbg.sk);
    b.SetFormat(Format::COEFFICIENT);
    auto bigB       = b.CRTInterpolate();
    auto& bigCoeffs = bigB.GetValues();

    auto sf          = ct->GetScalingFactor();
    const auto& q    = b.GetParams()->GetModulus();
    const auto& half = q >> 1;
    const size_t N   = b.GetRingDimension();
    auto safe_stride = [&](size_t denom) {
        return std::max<size_t>(1, N / denom);
    };

    // Extract the relevant slot/coefficient values based on the stage.
    std::vector<double> values;
    for (size_t i = 0; i < N; ++i) {
        auto x   = bigCoeffs[i];
        bool neg = (x > half);
        if (neg)
            x = q - x;
        double value = x.ConvertToDouble() / sf;
        if (neg)
            value = -value;

        if (msg == "CoeffsToSlots" || msg == "Exp" || msg == "LUT" || msg == "PreEvalExp") {
            if (i % safe_stride(dbg.slots * 4) == 0)
                values.push_back(value);
        }
        else if (msg == "CoeffsToSlots0" || msg == "Exp0" || msg == "LUT0" || msg == "CoeffsToSlots1" ||
                 msg == "Exp1" || msg == "LUT1" || msg == "PreEvalExp0" || msg == "PreEvalExp1") {
            values.push_back(value);
        }
        else if ((msg == "ModRaise" || msg == "Input" || msg == "Output")) {
            if (i % safe_stride(dbg.slots * 2) == 0)
                values.push_back(value);
        }
    }

    // Stage-specific scaling corrections.
    for (auto& v : values) {
        if (isC2SStage)
            v *= (dbg.skDist == lbcrypto::SPARSE_TERNARY) ? 25.0 : 16.0;
    }

    // Store per-stage values.
    if (msg == "Input") {
        dbg.encrypted.input = values;
        dbg.precise.input   = get_precise_array(values, dbg.pIn);
    }
    else if (msg == "ModRaise") {
        dbg.encrypted.modraise = values;
        dbg.precise.modraise   = get_precise_array(values, dbg.pIn);
    }
    else if (msg == "Output") {
        dbg.encrypted.output        = values;
        dbg.precise.output          = get_precise_array(values, dbg.pOut);
        g_precision.outputPrecision = get_precision_array(values, dbg.pOut);
    }

    // --- CoeffsToSlots and LUT stages ---
    if (isC2SStage || msg.find("LUT") != std::string::npos) {
        auto complexValues = multiplyByUComplex(values);

        if (isC2SStage) {
            if (msg != "CoeffsToSlots1") {
                dbg.encrypted.c2s = complexValues;
                dbg.precise.c2s   = get_precise_array_complex(complexValues, dbg.pIn);
            }
            else {
                dbg.encrypted.c2s.insert(dbg.encrypted.c2s.end(), complexValues.begin(), complexValues.end());
                auto p = get_precise_array_complex(complexValues, dbg.pIn);
                dbg.precise.c2s.insert(dbg.precise.c2s.end(), p.begin(), p.end());
            }

            g_precision.inputPrecision = get_precision_array_complex(dbg.encrypted.c2s, dbg.pIn);
        }

        if (msg.find("LUT") != std::string::npos) {
            if (msg != "LUT1") {
                dbg.encrypted.lut = complexValues;
                dbg.precise.lut   = get_precise_array_complex(complexValues, dbg.pOut);
            }
            else {
                dbg.encrypted.lut.insert(dbg.encrypted.lut.end(), complexValues.begin(), complexValues.end());
                auto p = get_precise_array_complex(complexValues, dbg.pOut);
                dbg.precise.lut.insert(dbg.precise.lut.end(), p.begin(), p.end());
            }

            g_precision.lutPrecision = get_precision_array_complex(dbg.encrypted.lut, dbg.pOut);
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------
// ArbitraryLUT: run functional bootstrapping for a given LUT function and
//   set of scheme parameters, then report correctness and timing.
// ---------------------------------------------------------------------------
static std::vector<std::complex<double>> compute_coefficients(
    InterpolationMethod method, bool binaryLUT, std::function<int64_t(int64_t)> func, BigInteger PInput, size_t order,
    uint64_t scaleTHI, std::vector<int64_t>& coeffint,
    std::function<std::complex<double>(int64_t)> complexFunc = nullptr, bool complexLUT = false) {
    std::vector<std::complex<double>> coeffcomp;
    if (complexLUT) {
        if (!complexFunc)
            OPENFHE_THROW("A complex-valued LUT requires a complex-valued function.");
        if (method == lbcrypto::HERMITE_AKP25) {
            coeffcomp = GetHermiteTrigCoefficientsForComplexLUT(complexFunc, PInput.ConvertToInt(), order, scaleTHI);
        }
        else if (method == lbcrypto::HERMITE_SPARSE_THI) {
            coeffcomp =
                GetHermiteTrigCoefficientsSparseTHIForComplexLUT(complexFunc, PInput.ConvertToInt(), order, scaleTHI);
        }
        else if (method == lbcrypto::HERMITE_FULL_THI) {
            coeffcomp =
                GetHermiteTrigCoefficientsFullTHIForComplexLUT(complexFunc, PInput.ConvertToInt(), order, scaleTHI);
        }
        else if (method == lbcrypto::HERMITE_BKSS24) {
            coeffcomp = GetHermiteTrigCoefficientsBKSSForComplexLUT(complexFunc, PInput.ConvertToInt(), scaleTHI);
        }
        else {
            OPENFHE_THROW("Complex-valued LUTs not supported for the selected interpolation method.");
        }
        return coeffcomp;
    }

    if (method == lbcrypto::HERMITE_AKP25) {
        if (binaryLUT) {
            // Coefficients for [1, cos^2(pi x)], not [1, cos(2pi x)].
            coeffint = {func(1), func(0) - func(1)};
        }
        else {
            coeffcomp = GetHermiteTrigCoefficients(func, PInput.ConvertToInt(), order, scaleTHI);
        }
    }
    else if (method == lbcrypto::HERMITE_FULL_THI) {
        coeffcomp = GetHermiteTrigCoefficientsFullTHI(func, PInput.ConvertToInt(), order, scaleTHI);
    }
    else if (method == lbcrypto::HERMITE_SPARSE_THI) {
        coeffcomp = GetHermiteTrigCoefficientsSparseTHI(func, PInput.ConvertToInt(), order, scaleTHI);
    }
    else if (method == lbcrypto::HERMITE_BKSS24) {
        coeffcomp = GetHermiteTrigCoefficientsBKSS(func, PInput.ConvertToInt(), scaleTHI);
    }
    return coeffcomp;
}

double compute_threshold(bool complexLUT, InterpolationMethod method, std::function<int64_t(int64_t)> func,
                         std::function<std::complex<double>(int64_t)> complexFunc, BigInteger PInput, size_t order,
                         uint64_t scaleTHI) {
    double inputPrecThreshold = 0.0;
    if (complexLUT) {
        switch (method) {
            case lbcrypto::HERMITE_AKP25:
                inputPrecThreshold =
                    GetHermiteTrigAKPThresholdForComplexLUT(complexFunc, PInput.ConvertToInt(), order, scaleTHI);
                break;
            case lbcrypto::HERMITE_SPARSE_THI:
                inputPrecThreshold =
                    GetHermiteTrigSparseTHIThresholdForComplexLUT(complexFunc, PInput.ConvertToInt(), order, scaleTHI);
                break;
            case lbcrypto::HERMITE_FULL_THI:
                inputPrecThreshold =
                    GetHermiteTrigFullTHIThreshold(complexFunc, PInput.ConvertToInt(), order, scaleTHI);
                break;
            default:
                break;
        }
    }
    else {
        switch (method) {
            case lbcrypto::HERMITE_AKP25:
            case lbcrypto::HERMITE_BKSS24:
                inputPrecThreshold = GetHermiteTrigAKPThreshold(func, PInput.ConvertToInt(), order, scaleTHI);
                break;
            case lbcrypto::HERMITE_FULL_THI:
                inputPrecThreshold = GetHermiteTrigFullTHIThreshold(func, PInput.ConvertToInt(), order, scaleTHI);
                break;
            case lbcrypto::HERMITE_SPARSE_THI:
                inputPrecThreshold = GetHermiteTrigSparseTHIThreshold(func, PInput.ConvertToInt(), order, scaleTHI);
                break;
            default:
                break;
        }
    }
    return inputPrecThreshold;
}

std::vector<ArbitraryLUTRunResult> ArbitraryLUT(
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput, BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, bool allInputs, NoiseInjectionStage noiseInjectionStage,
    lbcrypto::ScalingTechnique scalingTechnique, lbcrypto::SecretKeyDist secretKeyDist,
    const std::vector<uint32_t>& levelBudget, std::function<int64_t(int64_t)> func,
    const std::vector<ArbitraryLUTRunRequest>& runRequests, InterpolationMethod method,
    std::function<std::complex<double>(int64_t)> complexFunc, bool complexLUT) {
    /* 1. Determine packing mode.
     * If numSlots <= ringDim/2 we use sparse packing; otherwise full packing. */
    bool flagSP       = (numSlots <= ringDim / 2);
    auto numSlotsCKKS = flagSP ? numSlots : numSlots / 2;
    if (complexLUT && !flagSP)
        OPENFHE_THROW("Complex-valued LUT benchmark currently supports sparse packing only.");

    /* 2. Build the input vector. */
    std::vector<int64_t> x;
    if (allInputs) {
        auto pInputInt = PInput.ConvertToInt<uint32_t>();
        if (numSlots < pInputInt)
            OPENFHE_THROW("all-inputs requires slots >= p_input");
        x.reserve(numSlots);
        for (uint32_t i = 0; i < pInputInt; ++i)
            x.push_back(static_cast<int64_t>(i));
        while (x.size() < numSlots)
            x.push_back(static_cast<int64_t>(x.size() % pInputInt));
    }
    else {
        x.reserve(numSlots);
        for (uint32_t i = 0; i < numSlots; ++i)
            x.push_back(static_cast<int64_t>(i % PInput.ConvertToInt()));
    }

    /* 3. Compute LUT coefficients.
     * Boolean LUTs with first-order Hermite interpolation can use real (integer)
     * coefficients; the general case uses complex coefficients. */
    bool binaryLUT = !complexLUT && (PInput.ConvertToInt() == 2) && (order == 1);
    std::vector<int64_t> coeffint;
    auto coeffcomp =
        compute_coefficients(method, binaryLUT, func, PInput, order, scaleTHI, coeffint, complexFunc, complexLUT);

    g_debug.coeffComp = coeffcomp;

    /* 4. Build CKKS crypto parameters. */
    uint32_t dcrtBits                       = Bigq.GetMSB() - 1;
    uint32_t firstMod                       = Bigq.GetMSB() - 1;
    uint32_t levelsAvailableAfterBootstrap  = 0;
    uint32_t levelsAvailableBeforeBootstrap = 0;
    uint32_t dnum                           = 3;
    std::vector<uint32_t> lvlb              = levelBudget;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(__unsafe_global ? lbcrypto::HEStd_NotSet : lbcrypto::HEStd_128_classic);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(scalingTechnique);
    parameters.SetFirstModSize(firstMod);
    parameters.SetNumLargeDigits(dnum);
    parameters.SetBatchSize(numSlotsCKKS);
    parameters.SetRingDim(ringDim);

    uint32_t depth = levelsAvailableAfterBootstrap;
    if (binaryLUT)
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffint, PInput, order, secretKeyDist, method);
    else
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffcomp, PInput, order, secretKeyDist, method);
    parameters.SetMultiplicativeDepth(depth);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    /* 5. Generate keys and populate debug context. */
    auto keyPair = cc->KeyGen();

    g_debug.slots      = numSlots;
    g_debug.pIn        = PInput.ConvertToInt();
    g_debug.pOut       = POutput.ConvertToInt();
    g_debug.cc         = cc;
    g_debug.sk         = keyPair.secretKey;
    g_debug.pInput     = PInput;
    g_debug.skDist     = secretKeyDist;
    g_debug.complexLUT = complexLUT;
    g_debug.order      = order;

    if (binaryLUT)
        cc->EvalFBTSetup(coeffint, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, 0, order);
    else
        cc->EvalFBTSetup(coeffcomp, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, 0, order);

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::vector<int64_t> exact;
    std::vector<std::complex<int64_t>> exactComplex;
    if (complexLUT) {
        exactComplex.resize(x.size());
        std::transform(x.begin(), x.end(), exactComplex.begin(), [&](int64_t elem) {
            auto value = complexFunc(elem);
            return std::complex<int64_t>(
                CenterOutputComponent(static_cast<int64_t>(std::llround(value.real())), POutput.ConvertToInt()),
                CenterOutputComponent(static_cast<int64_t>(std::llround(value.imag())), POutput.ConvertToInt()));
        });
    }
    else {
        exact = x;
        std::transform(x.begin(), x.end(), exact.begin(),
                       [&](int64_t elem) { return CenterOutputComponent(func(elem), POutput.ConvertToInt()); });
    }

    double inputPrecThreshold = 0.0;
    if (g_runMode == RunMode::Precision)
        inputPrecThreshold = compute_threshold(complexLUT, method, func, complexFunc, PInput, order, scaleTHI);
    auto elapsed_ms = [](time_point_t a, time_point_t b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    std::vector<ArbitraryLUTRunResult> results;
    results.reserve(runRequests.size());

    for (const auto& request : runRequests) {
        double noiseValue = std::ldexp(1.0, static_cast<int>(request.addNoiseBase)) / Bigq.ConvertToDouble();
        FHECKKSRNS::SetFBTPreEvalExpNoise(noiseValue);

        /* 6. Encrypt in RLWE with a larger initial modulus, then mod-switch down. */
        auto ctxtBFV = SchemeletRLWEMP::EncryptCoeff(x, QBFVInit, PInput, keyPair.secretKey, ep);
        SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVInit);

        /* 7. Convert RLWE -> CKKS (same secret key). */
        auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, ctxtBFV, keyPair.publicKey, Bigq, numSlotsCKKS,
                                                       depth - (levelsAvailableBeforeBootstrap > 0));

        auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersRNS>(cc->GetCryptoParameters());
        if (cryptoParams->GetScalingTechnique() == FLEXIBLEMANUAL) {
            ctxt->SetScalingFactor(Bigq.ConvertToDouble());
        }

        __debug(ctxt, "Input");

        /* 8. Evaluate the LUT. */
        Ciphertext<DCRTPoly> ctxtAfterFBT;
        if (binaryLUT)
            ctxtAfterFBT = cc->EvalFBT(ctxt, coeffint, PInput.GetMSB() - 1, ep->GetModulus(), scaleTHI, 0, order);
        else
            ctxtAfterFBT = cc->EvalFBT(ctxt, coeffcomp, PInput.GetMSB() - 1, ep->GetModulus(), scaleTHI, 0, order);

        __debug(ctxtAfterFBT, "Output");

        /* 9. Convert CKKS -> RLWE and decrypt. */
        std::vector<int64_t> computed;
        std::vector<std::complex<int64_t>> computedComplex;
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT, Q);
        if (complexLUT) {
            computed =
                SchemeletRLWEMP::DecryptCoeff(polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots * 2);
            for (size_t i = 0; i != computed.size() / 2; ++i) {
                computedComplex.emplace_back(computed[i], computed[i + computed.size() / 2]);
            }
        }
        else {
            computed = SchemeletRLWEMP::DecryptCoeff(polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        }

        if (g_runMode == RunMode::Precision) {
            g_precision.predictedOutputPrecision =
                PredictHermiteTrigNoise(inputPrecThreshold, order, g_precision.inputPrecision);
            g_precision.predictionGapBits = g_precision.outputPrecision - g_precision.predictedOutputPrecision;
        }

        /* 10. Compute max absolute error against exact output. */
        uint32_t maxError = 0;
        if (complexLUT) {
            for (size_t i = 0; i < exactComplex.size(); ++i) {
                maxError = std::max(maxError, ComponentError(computedComplex[i].real(), exactComplex[i].real(),
                                                             POutput.ConvertToInt()));
                maxError = std::max(maxError, ComponentError(computedComplex[i].imag(), exactComplex[i].imag(),
                                                             POutput.ConvertToInt()));
            }
        }
        else {
            for (size_t i = 0; i < exact.size(); ++i)
                maxError = std::max(maxError, ComponentError(computed[i], exact[i], POutput.ConvertToInt()));
        }

        double lutMs = std::numeric_limits<double>::quiet_NaN();
        if (g_runMode == RunMode::Benchmark)
            lutMs = elapsed_ms(GetExpStageTime(g_timing), GetLutStageTime(g_timing));

        results.push_back({request.runNumber, request.addNoiseBase, request.warmup, maxError,
                           elapsed_ms(g_timing.times["Input"], g_timing.times["Output"]), lutMs,
                           g_timing.keySwitchCounterAtEvalExp, g_timing.keySwitchCounterAtLUT, g_precision});
    }

    /* 11. Cleanup CryptoContext */
    cc->ClearStaticMapsAndVectors();

    return results;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
struct LutConfig {
    std::string name;
    std::function<int64_t(int64_t)> func;
    std::function<std::complex<double>(int64_t)> complexFunc;
};

static uint32_t ComputeConfiguredMulDepth(const BenchmarkOptions& options, const std::vector<LutConfig>& lutConfigs) {
    uint32_t maxDepth           = 0;
    SecretKeyDist secretKeyDist = options.secretKeyDist;
    std::vector<uint32_t> lvlb  = options.levelBudget;

    for (const auto& lutConfig : lutConfigs) {
        for (auto method : options.methods) {
            bool binaryLUT = !options.complexLUT && (options.pInput == 2) && (options.order == 1);
            std::vector<int64_t> coeffint;
            auto coeffcomp =
                compute_coefficients(method, binaryLUT, lutConfig.func, BigInteger(options.pInput), options.order,
                                     options.pOutput, coeffint, lutConfig.complexFunc, options.complexLUT);

            uint32_t depth = 0;
            if (binaryLUT)
                depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffint, BigInteger(options.pInput), options.order,
                                                 secretKeyDist, method);
            else
                depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffcomp, BigInteger(options.pInput), options.order,
                                                 secretKeyDist, method);
            maxDepth = std::max(maxDepth, depth);
        }
    }

    return maxDepth;
}

int main(int argc, char* argv[]) {
    constexpr int kWarmupRuns   = 1;
    constexpr int kMeasuredRuns = 5;

    BenchmarkOptions options;
    try {
        options = ParseBenchmarkOptions(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
    g_runMode                = options.runMode;
    __eval_exp_degree_global = options.evalExpDegree;
    __complex_lut_global     = options.complexLUT;
    __unsafe_global          = options.unsafe;

    std::vector<LutConfig> lutConfigs;
    if (options.complexLUT) {
        auto complexIdFunc = [pInput = options.pInput, pOutput = options.pOutput](int64_t x) {
            auto value = x % static_cast<int64_t>(pInput);
            value      = (value > static_cast<int64_t>(pOutput / 2)) ? value - static_cast<int64_t>(pOutput) : value;
            return std::complex<double>(static_cast<double>(value), static_cast<double>(value));
        };
        lutConfigs.push_back(
            {"ComplexID",
             [complexIdFunc](int64_t x) { return static_cast<int64_t>(std::llround(complexIdFunc(x).real())); },
             complexIdFunc});
#if 0
        auto complexSymbol = [pOutput = options.pOutput](int64_t x) {
            const double half = static_cast<double>(pOutput / 2);
            switch (static_cast<uint32_t>(x) & 3U) {
                case 0:
                    return std::complex<double>(0.0, 0.0);
                case 1:
                    return std::complex<double>(half, 0.0);
                case 2:
                    return std::complex<double>(0.0, half);
                default:
                    return std::complex<double>(half, half);
            }
        };
        lutConfigs.push_back(
            {"ComplexSymbol",
             [complexSymbol](int64_t x) { return static_cast<int64_t>(std::llround(complexSymbol(x).real())); },
             complexSymbol});
#endif
    }
    else {
        auto idFunc = [pInput = options.pInput, pOutput = options.pOutput](int64_t x) {
            auto value = x % static_cast<int64_t>(pInput);
            return (value > static_cast<int64_t>(pOutput / 2)) ? value - static_cast<int64_t>(pOutput) : value;
        };
        [[maybe_unused]] auto aesFunc = [pInput = options.pInput, pOutput = options.pOutput](int64_t x) {
            auto idx   = static_cast<size_t>(x % static_cast<int64_t>(pInput));
            auto value = static_cast<int64_t>(aesSBox[idx % aesSBox.size()] % pOutput);
            return (value > static_cast<int64_t>(pOutput / 2)) ? value - static_cast<int64_t>(pOutput) : value;
        };
        lutConfigs = {
            {"ID", idFunc,
             [idFunc](int64_t x) {
                 return std::complex<double>(static_cast<double>(idFunc(x)), 0.0);
             }},
#if 0
            {"AES_SBox", aesFunc, [aesFunc](int64_t x) {
                 return std::complex<double>(static_cast<double>(aesFunc(x)), 0.0);
             }}
#endif
        };
    }

    uint32_t mulDepth = ComputeConfiguredMulDepth(options, lutConfigs);

    std::vector<MethodResult> results;
    std::vector<RunResult> runResults;
    const int totalRuns   = (g_runMode == RunMode::Benchmark) ? (kWarmupRuns + kMeasuredRuns) : 1;
    size_t verifyFailures = 0;
    size_t verifyChecks   = 0;
    PrintBenchmarkHeader(options, options.secretKeyDist, options.slots, mulDepth);
    const std::string benchmarkLutName = options.complexLUT ? "ComplexID" : "ID";

    for (const auto& lutConfig : lutConfigs) {
        if (g_runMode == RunMode::Benchmark && lutConfig.name != benchmarkLutName)
            continue;
        for (auto method : options.methods) {
            __interpolation_method_global = method;
            std::vector<double> overallTimings;
            std::vector<double> lutTimings;
            uint32_t evalExpKeySwitchCount = 0;
            uint32_t lutKeySwitchCount     = 0;

            std::vector<ArbitraryLUTRunRequest> lutRunRequests;
            for (auto addNoiseBase : options.addNoiseBases) {
                for (int run = 0; run < totalRuns; ++run) {
                    bool warmup  = (g_runMode == RunMode::Benchmark) && (run < kWarmupRuns);
                    int runIndex = (g_runMode == RunMode::Benchmark) ? (run - kWarmupRuns + 1) : 1;
                    lutRunRequests.push_back({runIndex, addNoiseBase, warmup});
                }
            }

            auto lutRunResults = ArbitraryLUT(
                QBFVINIT, BigInteger(options.pInput), BigInteger(options.pOutput),
                (BigInteger(1) << options.scalingFactor), (BigInteger(1) << options.scalingFactor), options.pOutput,
                options.order, options.slots, options.ringDim, options.allInputs, options.noiseInjectionStage,
                options.scalingTechnique, options.secretKeyDist, options.levelBudget, lutConfig.func, lutRunRequests,
                method, lutConfig.complexFunc, options.complexLUT);

            for (const auto& runResult : lutRunResults) {
                if (g_runMode == RunMode::Benchmark && !runResult.warmup) {
                    evalExpKeySwitchCount = runResult.evalExpKeySwitchCount;
                    lutKeySwitchCount     = runResult.lutKeySwitchCount;
                    overallTimings.push_back(runResult.overallMs);
                    lutTimings.push_back(runResult.lutMs);
                    PrintBenchmarkRunLine(lutConfig.name, method, runResult.runNumber, runResult.addNoiseBase,
                                          runResult.maxError, runResult.overallMs, runResult.lutMs,
                                          evalExpKeySwitchCount, lutKeySwitchCount);
                }
                else if (g_runMode == RunMode::Precision) {
                    runResults.push_back({lutConfig.name, method, runResult.runNumber, runResult.addNoiseBase,
                                          runResult.maxError, runResult.lutMs, runResult.evalExpKeySwitchCount,
                                          runResult.lutKeySwitchCount, runResult.precision.inputPrecision,
                                          runResult.precision.lutPrecision, runResult.precision.lutMarginBits,
                                          runResult.precision.outputPrecision,
                                          runResult.precision.predictedOutputPrecision,
                                          runResult.precision.predictionGapBits, runResult.precision.outputMarginBits});
                    PrintPrecisionRunLine(lutConfig.name, method, runResult.addNoiseBase, runResult.maxError,
                                          runResult.precision, runResult.evalExpKeySwitchCount,
                                          runResult.lutKeySwitchCount);
                }
                else if (g_runMode == RunMode::Verify) {
                    ++verifyChecks;
                    if (runResult.maxError != 0)
                        ++verifyFailures;
                    PrintVerifyLine(lutConfig.name, method, runResult.maxError);
                }
            }

            if (g_runMode == RunMode::Benchmark) {
                double avgOverall =
                    std::accumulate(overallTimings.begin(), overallTimings.end(), 0.0) / overallTimings.size();
                double maxOverall = *std::max_element(overallTimings.begin(), overallTimings.end());
                double avgLut     = std::accumulate(lutTimings.begin(), lutTimings.end(), 0.0) / lutTimings.size();
                double maxLut     = *std::max_element(lutTimings.begin(), lutTimings.end());
                MethodResult result{lutConfig.name,        method,           avgOverall, maxOverall, avgLut, maxLut,
                                    evalExpKeySwitchCount, lutKeySwitchCount};
                PrintBenchmarkSummaryLine(result);
                results.push_back(result);
            }
        }
    }

    if (g_runMode == RunMode::Precision) {
        if (options.emitCsv) {
            std::cout << "\n";
            PrintPrecisionCsv(runResults);
        }
    }
    else if (g_runMode == RunMode::Verify) {
        PrintVerifySummary(verifyChecks - verifyFailures, verifyChecks);
        return verifyFailures == 0 ? 0 : 2;
    }

    return 0;
}
