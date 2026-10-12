#!/usr/bin/env bash
# Reproduce everything: data, build, tests, samples, validation reports.
# Needs ROOT (root-config on PATH or ROOT_DIR) and a Python with numpy, scipy, matplotlib
# (set PY, e.g. PY="micromamba run -n nufitval python").
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
data/fetch_nufit.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release ${ROOTSYS:+-DCMAKE_PREFIX_PATH=$ROOTSYS}
cmake --build build -j
ctest --test-dir build --output-on-failure
mkdir -p results
for v in TBoff TByes; do
  for s in factorized coupled; do
    build/sample_direct --variant $v --surrogate $s --n 1000000 --out results/direct_${v}_${s}
  done
  build/sample_pt --variant $v --out results/pt_${v}_factorized
done
build/sample_pt --variant TBoff --surrogate coupled --out results/pt_TBoff_coupled
build/sample_pt --variant TBoff --jumps 0 --out results/pt_TBoff_factorized_nojumps
build/sample_direct --variant TBoff --angle-prior angle --dm21-prior log --n 1000000 \
  --out results/direct_TBoff_factorized_angleprior
# student's chain, for comparison (needs PyROOT)
python3 python/convert_legacy.py MCMCResults/vectorFile.root results/legacy_samples.npy || true
cd python
for cfg in "TBoff factorized" "TBoff coupled" "TByes factorized"; do
  set -- $cfg
  $PY -m nufitval.report --variant $1 --surrogate $2 --results ../results --data ../data/nufit
done
