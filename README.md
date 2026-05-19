Installation and Running instructions
=====================================

## Installation

```bash
# assuming the repo has been downloaded and we are in the root folder of the repo
mkdir build; cd build
CC=clang CXX=clang++ cmake -DWITH_INTEL_HEXL=ON -DINTEL_HEXL_HINT_DIR=$PWD/install -DCMAKE_INSTALL_PREFIX=$PWD/install -DMATHBACKEND=6 -DWITH_NTL=ON -DWITH_TCM=ON ..
make -j tcm
make -j fbt-benchmark
```

## Precision

The following commands are quick to run and verify:

```bash
./bin/examples/pke/fbt-benchmark --precision --methods=AKP25,BKSS24,FULLTHI,SPARSETHI --p-input=16 --p-output=16 --order=1 --slots=256 --ring-dim=256 --sf=59 --unsafe --noise-base=20,30,40
# Mode=precision ring_dim=256 key_dist=SPARSE_TERNARY slots=256 sf=59 eval_exp_degree=58 order=3 noise_stage=PRE_EVALEXP scaling=FIXEDMANUAL mul_depth=18 p_input=16 p_output=16 all_inputs=0 complex_lut=0 level_budget=1,1
# ID          AKP25                        noise=20 PASS max_error=0 in=-38.868 lut=-17.984 lut_ks=18
# ID          AKP25                        noise=30 PASS max_error=0 in=-29.163 lut=-17.984 lut_ks=18
# ID          AKP25                        noise=40 PASS max_error=0 in=-19.163 lut=-17.722 lut_ks=18
# ID          BKSS24                       noise=20 PASS max_error=0 in=-38.797 lut=-34.041 lut_ks=28
# ID          BKSS24                       noise=30 PASS max_error=0 in=-29.163 lut=-34.366 lut_ks=28
# ID          BKSS24                       noise=40 PASS max_error=0 in=-19.163 lut=-28.821 lut_ks=28
# ID          FULL THI                     noise=20 PASS max_error=0 in=-38.826 lut=-31.507 lut_ks=22
# ID          FULL THI                     noise=30 PASS max_error=0 in=-29.163 lut=-31.776 lut_ks=22
# ID          FULL THI                     noise=40 PASS max_error=0 in=-19.163 lut=-25.910 lut_ks=22
# ID          Sparse THI                   noise=20 PASS max_error=0 in=-38.868 lut=-33.709 lut_ks=22
# ID          Sparse THI                   noise=30 PASS max_error=0 in=-29.163 lut=-33.830 lut_ks=22
# ID          Sparse THI                   noise=40 PASS max_error=0 in=-19.163 lut=-27.871 lut_ks=22
./bin/examples/pke/fbt-benchmark --precision --methods=AKP25,BKSS24,FULLTHI,SPARSETHI --p-input=16 --p-output=16 --order=3 --slots=256 --ring-dim=256 --sf=59 --unsafe --noise-base=20,30,40
# ID          AKP25                        noise=20 PASS max_error=0 in=-38.866 lut=-18.358 lut_ks=24
# ID          AKP25                        noise=30 PASS max_error=0 in=-29.163 lut=-18.347 lut_ks=24
# ID          AKP25                        noise=40 PASS max_error=0 in=-19.163 lut=-18.610 lut_ks=24
# ID          BKSS24                       noise=20 PASS max_error=0 in=-38.808 lut=-34.000 lut_ks=28
# ID          BKSS24                       noise=30 PASS max_error=0 in=-29.163 lut=-34.416 lut_ks=28
# ID          BKSS24                       noise=40 PASS max_error=0 in=-19.163 lut=-28.601 lut_ks=28
# ID          FULL THI                     noise=20 PASS max_error=0 in=-38.886 lut=-39.712 lut_ks=32
# ID          FULL THI                     noise=30 PASS max_error=0 in=-29.163 lut=-39.712 lut_ks=32
# ID          FULL THI                     noise=40 PASS max_error=0 in=-19.163 lut=-39.712 lut_ks=32
# ID          Sparse THI                   noise=20 PASS max_error=0 in=-38.897 lut=-40.276 lut_ks=34
# ID          Sparse THI                   noise=30 PASS max_error=0 in=-29.163 lut=-40.273 lut_ks=34
# ID          Sparse THI                   noise=40 PASS max_error=0 in=-19.163 lut=-40.274 lut_ks=34
```

To obtain the real experimental data, use

```bash
./bin/examples/pke/fbt-benchmark --precision --methods=AKP25,BKSS24,FULLTHI,SPARSETHI --p-input=16 --p-output=16 --order=1 --slots=65536 --ring-dim=65536 --sf=59 --level-budget=3,2 --noise-base=20,30,40
# possibly use a much denser noise base
```

## Performance

The following commands are quick to run and verify:

```bash
./bin/examples/pke/fbt-benchmark --mode=benchmark --methods=AKP25,BKSS24,FULLTHI,SPARSETHI --p-input=16 --p-output=16 --order=1 --slots=16 --ring-dim=2048 --sf=59 --unsafe --level-budget=3,2
./bin/examples/pke/fbt-benchmark --mode=benchmark --methods=AKP25,BKSS24,FULLTHI,SPARSETHI --p-input=256 --p-output=256 --order=1 --slots=2048 --ring-dim=2048 --sf=59 --unsafe --level-budget=3,2
```

To obtain the real experimental data, use commands like

```bash
OMP_NUM_THREADS=1 taskset -c 0 ./bin/examples/pke/fbt-benchmark --mode=benchmark --methods=AKP25,BKSS24,FULLTHI,SPARSETHI --p-input=16 --p-output=16 --order=1 --slots=65536 --ring-dim=65536 --sf=59 --level-budget=3,2
```

OpenFHE - Open-Source Fully Homomorphic Encryption Library
=====================================

Fully Homomorphic Encryption (FHE) is a powerful cryptographic primitive that enables performing computations over encrypted data without having access to the secret key.
OpenFHE is an open-source FHE library that includes efficient implementations of all common FHE schemes:

* Brakerski/Fan-Vercauteren (BFV) scheme for integer arithmetic
* Brakerski-Gentry-Vaikuntanathan (BGV) scheme for integer arithmetic
* Cheon-Kim-Kim-Song (CKKS) scheme for real-number arithmetic (includes approximate bootstrapping)
* Ducas-Micciancio (DM/FHEW), Chillotti-Gama-Georgieva-Izabachene (CGGI/TFHE), and Lee-Micciancio-Kim-Choi-Deryabin-Eom-Yoo (LMKCDEY) schemes for evaluating Boolean circuits and arbitrary functions over larger plaintext spaces using lookup tables

OpenFHE also supports hybrid vectorized schemes, with the goal of enabling the FHEW/TFHE-like functional bootstrapping capability for schemes such as CKKS and BFV. In particular, OpenFHE supports

* Switching between CKKS and FHEW/TFHE to evaluate non-smooth functions, e.g., comparison, using (scalar) FHEW/TFHE functional bootstrapping
* Switching between RLWE (a scheme equivalent to the coefficient-encoded additive BFV scheme) and CKKS to evaluate arbitrary lookup tables over vectors of integers, e.g., modular reduction, comparison or S-box, using vectorized functional bootstrapping implemented in CKKS

OpenFHE also supports partial schemes, called schemelets, such as RLWE which is equivalent to the coefficient-encoded additive BFV scheme. In OpenFHE, the RLWE schemelet is the starting point for the vectorized functional bootstrapping capability, which allows the evaluation of arbitrary lookup tables over vectors of integers, e.g., modular reduction, comparison or Sbox, using CKKS in an intermediate step.

OpenFHE also includes the following multiparty extensions of FHE:

* Threshold FHE for BGV, BFV, and CKKS schemes
* Interactive bootstrapping for Threshold CKKS
* Proxy Re-Encryption for BGV, BFV, and CKKS schemes

OpenFHE supports any GNU C++ compiler version 9 or above and clang C++ compiler version 10 or above. To achieve the best runtime performance, we recommend following the
guidelines outlined in [building OpenFHE for best performance](https://github.com/openfheorg/openfhe-development/blob/main/docs/static_docs/Best_Performance.md).

## Links and Resources

* [OpenFHE documentation](https://openfhe-development.readthedocs.io/en/latest/)
* [Design paper for OpenFHE](https://eprint.iacr.org/2022/915)
* [OpenFHE website](https://openfhe.org)
* [Community forum for OpenFHE](https://openfhe.discourse.group/)
* [OpenFHE Release Notes](https://github.com/openfheorg/openfhe-development/blob/main/docs/static_docs/Release_Notes.md)
* [Quickstart](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/quickstart.html)
* [BSD 2-Clause License](LICENSE)
* [Contributing to OpenFHE](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/contributing/contributing.html)
* [OpenFHE Governance](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/misc/governance.html)
* [Openfhe-development Github Issues](https://github.com/openfheorg/openfhe-development/issues)
* To report security vulnerabilities, please email us at <contact@openfhe.org>

## Installation

Refer to our General Installation Information: [readthedocs](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html) for more information

Or refer to the following for your specific operating system:

* [Linux](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/linux.html)

* [MacOS](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/macos.html)

* [Windows](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/windows.html)

## Code Examples

To get familiar with the main API of OpenFHE, we recommend looking at the code of the following examples:

   1. FHE for arithmetic over integers (BFV):
       1. [Simple Code Example](src/pke/examples/simple-integers.cpp)
       2. [Simple Code Example with Serialization](src/pke/examples/simple-integers-serial.cpp)
   1. FHE for arithmetic over integers (BGV):
       1. [Simple Code Example](src/pke/examples/simple-integers-bgvrns.cpp)
       2. [Simple Code Example with Serialization](src/pke/examples/simple-integers-serial-bgvrns.cpp)
   1. FHE for arithmetic over real numbers (CKKS):
       1. [Simple Code Example](src/pke/examples/simple-real-numbers.cpp)
       2. [Advanced Code Example](src/pke/examples/advanced-real-numbers.cpp)
       3. [Advanced Code Example for High-Precision CKKS](src/pke/examples/advanced-real-numbers-128.cpp)
       4. [Arbitrary Smooth Function Evaluation](src/pke/examples/function-evaluation.cpp)
       5. [Simple CKKS Bootstrapping Example](src/pke/examples/simple-ckks-bootstrapping.cpp)
       6. [Advanced CKKS Bootstrapping Example](src/pke/examples/advanced-ckks-bootstrapping.cpp)
       7. [Double-Precision (Iterative) Bootstrapping Example](src/pke/examples/iterative-ckks-bootstrapping.cpp)
       8. [Basic CKKS Arithmetic in the CKKS Composite Scaling Mode](src/pke/examples/simple-real-numbers-composite-scaling.cpp)
   1. FHE for arithmetic over complex numbers (CKKS):
       1. [Leveled and Boostrapping Code Examples](src/pke/examples/simple-complex-numbers.cpp)
   1. FHE for Boolean circuits and larger plaintext spaces (FHEW/TFHE):
       1. [Simple Code Example with Symmetric Encryption](src/binfhe/examples/boolean.cpp)
       2. [Simple Code Example with PKE](src/binfhe/examples/pke/boolean-pke.cpp)
       3. [Evaluation of Multi-Input Gates](src/binfhe/examples/boolean-multi-input.cpp)
       4. [Code with JSON serialization](src/binfhe/examples/boolean-serial-json.cpp)
       5. [Code with Binary Serialization](src/binfhe/examples/boolean-serial-binary.cpp)
       6. [Large-Precision Comparison](src/binfhe/examples/eval-sign.cpp)
       7. [Small-Precison Arbitrary Function Evaluation](src/binfhe/examples/eval-function.cpp)
   1. Scheme Switching:
       1. [Examples with Scheme Switching between CKKS and FHEW/TFHE](src/pke/examples/scheme-switching.cpp)
   1. Functional Bootstrapping over integers (RLWE and CKKS):
       1. [Examples with Functional Bootstrapping using CKKS for lookup table evaluation, sign extraction and multivalue bootstrapping](src/pke/examples/functional-bootstrapping-ckks.cpp)
   1. Threshold FHE:
       1. [Code Example for BGV, BFV, and CKKS](src/pke/examples/threshold-fhe.cpp)
       1. [2-party Interactive Bootstrapping Examples](src/pke/examples/interactive-bootstrapping.cpp)
       1. [Simple n-party Interactive Bootstrapping Example](src/pke/examples/tckks-interactive-mp-bootstrapping.cpp)
       1. [n-party Interactive Bootstrapping after Chebyshev Approximation](src/pke/examples/tckks-interactive-mp-bootstrapping-Chebyshev.cpp)
       1. [Code Example for BFV with 5 parties](src/pke/examples/threshold-fhe-5p.cpp)

## Main API

* [PKE CryptoContext API (BGV/BFV/CKKS)](https://openfhe-development.readthedocs.io/en/latest/api/classlbcrypto_1_1CryptoContextImpl.html)
* [Description of CryptoContext Parameters for BGV, BFV, and CKKS](https://github.com/openfheorg/openfhe-development/blob/main/src/pke/examples#description-of-the-cryptocontext-parameters-and-their-restrictions)

* [BinFHE Context API (FHEW/TFHE)](https://openfhe-development.readthedocs.io/en/latest/api/classlbcrypto_1_1BinFHEContext.html)

## Code of Conduct

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making
participation in our project and our community a harassment-free experience for everyone, regardless of age, body size,
disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education,
socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

OpenFHE is a community-driven open source project developed by a diverse group of
[contributors](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/misc/contributors.html). The OpenFHE leadership has made a strong commitment to creating an open,
inclusive, and positive community. Please read our
[Code of Conduct](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/misc/code_of_conduct.html?highlight=code%20of%20) for guidance on how to interact with others in a way that
makes our community thrive.

## Call for Contributions

We welcome all contributions including but not limited to:

* [reporting issues](https://github.com/openfheorg/openfhe-development/issues)
* addressing [bugs](https://github.com/openfheorg/openfhe-development/issues) big or small. We label issues to help you filter them to your skill level.
* documentation changes
* talks and seminars using OpenFHE

## How to Cite OpenFHE

To cite OpenFHE in academic papers, please use the following BibTeX entry (updated version)

```
@misc{OpenFHE,
      author = {Ahmad Al Badawi and Andreea Alexandru and Jack Bates and Flavio Bergamaschi and David Bruce Cousins and Saroja Erabelli and Nicholas Genise and Shai Halevi and Hamish Hunt and Andrey Kim and Yongwoo Lee and Zeyu Liu and Daniele Micciancio and Carlo Pascoe and Yuriy Polyakov and Ian Quah and Saraswathy R.V. and Kurt Rohloff and Jonathan Saylor and Dmitriy Suponitsky and Matthew Triplett and Vinod Vaikuntanathan and Vincent Zucca},
      title = {{OpenFHE}: Open-Source Fully Homomorphic Encryption Library},
      howpublished = {Cryptology ePrint Archive, Paper 2022/915},
      year = {2022},
      note = {\url{https://eprint.iacr.org/2022/915}},
      url = {https://eprint.iacr.org/2022/915}
}
```

or, alternatively (original WAHC@CCS'22 version),

```
@inproceedings{10.1145/3560827.3563379,
      author = {Al Badawi, Ahmad and Bates, Jack and Bergamaschi, Flavio and Cousins, David Bruce and Erabelli, Saroja and Genise, Nicholas and Halevi, Shai and Hunt, Hamish and Kim, Andrey and Lee, Yongwoo and Liu, Zeyu and Micciancio, Daniele and Quah, Ian and Polyakov, Yuriy and R.V., Saraswathy and Rohloff, Kurt and Saylor, Jonathan and Suponitsky, Dmitriy and Triplett, Matthew and Vaikuntanathan, Vinod and Zucca, Vincent},
      title = {OpenFHE: Open-Source Fully Homomorphic Encryption Library},
      year = {2022},
      publisher = {Association for Computing Machinery},
      address = {New York, NY, USA},
      url = {https://doi.org/10.1145/3560827.3563379},
      doi = {10.1145/3560827.3563379},
      booktitle = {Proceedings of the 10th Workshop on Encrypted Computing \& Applied Homomorphic Cryptography},
      pages = {53-63},
      numpages = {11},
      location = {Los Angeles, CA, USA},
      series = {WAHC'22}
}
```

## Acknowledgments ##

Distribution Statement "A" (Approved for Public Release, Distribution Unlimited). This work is supported in part by DARPA through HR0011-21-9-0003 and HR0011-20-9-0102. The views, opinions, and/or findings expressed are those of the author(s) and should not be interpreted as representing the official views or policies of the Department of Defense or the U.S. Government.
