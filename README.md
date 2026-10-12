# NuFitSampler

Samples the 3-flavour oscillation parameters from the NuFIT 6.0 global fit. The intended
use is as the "world knowledge" input to a model-independent consistency (goodness-of-fit)
test of the MINOS+ 6–12 GeV νe-appearance data.

The student's original code (`NuFuncs.*`, `libra.*`, `test.cpp`, `Nuprocess.py`,
`MCMCResults/`) is kept unchanged for reference. Its problems are listed under
[Review of the original implementation](#review-of-the-original-implementation).

## Quick start

```bash
data/fetch_nufit.sh                    # NuFIT 6.0 tables -> data/nufit (sha256-checked)
cmake -S . -B build -DCMAKE_PREFIX_PATH=$ROOTSYS && cmake --build build -j
ctest --test-dir build                 # interpolation, exact sampler, PT toy, NuFIT checks
build/sample_direct --variant TBoff --n 1000000 --out results/direct_TBoff_factorized
build/sample_pt     --variant TBoff --out results/pt_TBoff_factorized
cd python && python -m nufitval.report --variant TBoff --surrogate factorized
```

`run_all.sh` reproduces every result below. Each takes a few seconds: 1M exact samples
take about 2 s, and a 4-replica × 8-rung × 550k-iteration PT run takes about 3 s.

The Python suite needs numpy, scipy and matplotlib. On this machine they live in the
micromamba env `nufitval`:
`PY="micromamba run -n nufitval python" ./run_all.sh`.

### Feeding the p-value machinery

`build/sample_direct` writes **independent** draws, so no thinning is needed:

* `<out>_samples.npy` holds columns `s12 s13 s23 dcp dm3l dm21 logL logPrior`.
  * `dcp` is in degrees.
  * `dm3l` is in 1e-3 eV²: it is Δm²31 for NO (> 0) and Δm²32 for IO (< 0), the NuFIT convention.
  * `dm21` is in 1e-5 eV².
* `<out>.root` holds TTree `samples` with the same columns plus `dm31_eV2`, `dm32_eV2` and
  `ordering` (0 = NO, 1 = IO).

Options:

* `--variant TBoff|TByes`: NuFIT "IC19 w/o SK-atm" or "IC24 with SK-atm".
* `--surrogate factorized|coupled`: see below.
* `--angle-prior sin2|angle`: flat in sin²θ or in θ.
* `--dm21-prior linear|log`.

## What is sampled

NuFIT publishes only projections of its 6D χ²: one 3D table (sin²θ23, Δm²3ℓ, δ), fifteen 2D
tables and six 1D tables, separately for NO and IO. All values are relative to the global
minimum over both orderings. The surrogate likelihood is `exp(-Δχ²/2)` with

* **factorized** (the student's construction):
  `χ²_3D(s23, Δm²3ℓ, δ) + [χ²_2D(s12, log Δm²21) − m_O] + [χ²_1D(s13) − m_O]`
* **coupled**: the θ13 term is replaced by `χ²_{T13/T23}(s13, s23) − χ²_{T23}(s23)`, which keeps
  the θ13–θ23 correlation.

Here `m_O` is the minimum within ordering O. Priors are flat on the NuFIT table windows, and
both orderings get equal prior mass. Tables are interpolated N-linearly on their own grids
(`include/GridTable.h`).

**Samplers.**

* `DirectSampler` draws **exact i.i.d.** samples by cell-wise rejection. An N-linear χ² is
  extremal at cell corners, so the corner maximum of `exp(-χ²/2)` is a strict envelope.
  Acceptance is about 55%.
* `PTSampler` is a general-purpose parallel-tempering MCMC:
  * tempered likelihood (the prior is not tempered);
  * adaptive Haario random walk;
  * a ladder adapted during burn-in only (Vousden et al. 2016), then frozen;
  * non-reversible even/odd swaps (Syed et al. 2021);
  * involutive mode jumps: ordering mirror, octant flip, δ → 180° − δ.

  PT isn't needed for the current surrogate. It is built and validated against the exact
  sampler so that it is ready for likelihoods that don't factorize (for example a joint
  NuFIT + MINOS+ fit).

## Results (all in `results/report_*/report.pdf`)

| check | TBoff factorized | TBoff coupled | TByes factorized |
|---|---|---|---|
| Python mirror vs C++ log L | 4e-15 | 4e-15 | 5e-15 |
| P(IO): quadrature / exact sampler | 0.4694 / 0.4692 | 0.4701 / 0.4702 | 0.0293 / 0.0293 |
| P(IO): PT ± MCSE | 0.4688 ± 0.0027 | 0.4698 ± 0.0027 | 0.0298 ± 0.0005 |
| PT max R-hat / min ESS (bulk, tail) | 1.0001 / 33k | 1.0002 / 35k | 1.0000 / 101k |
| P(upper octant) | 0.713 | 0.714 | 0.429 |

The exact sampler was also checked with 2×10⁷ draws. All twelve per-ordering marginals pass
a χ² test against quadrature, and P(IO) agrees to 1.2σ at a precision of 10⁻⁴.

**Findings that matter for the physics analysis**

1. **The factorized surrogate drops θ13 correlations at the Δχ² ≈ 2–3 level.** Compare its
   profiles with the 2D NuFIT tables it was not built from. Inside 3σ, the largest residuals
   are T13/DCP 2.5, T13/T23 2.0 (3.0 for IO) and T13/DMA 1.4. Every solar-sector projection
   agrees to ≤ 0.06. The coupled surrogate makes T13/T23 exact but leaves T13/DCP at about 2.7.
   The 1D marginals barely change, but for a νe-appearance test the surrogate choice should be
   carried as a systematic.
2. **The NuFIT variant matters far more than the priors.** TBoff to TByes moves P(IO) from 0.47
   to 0.03 and P(upper octant) from 0.71 to 0.43, and reshapes δCP. Switching to flat-in-angle
   and log Δm²21 priors changes nothing visible. Quote the p-value for both variants.
3. **Interpolation.** Linear interpolation is exact at nodes and cannot overshoot. In a
   decimated-grid (leave-one-out) test at *twice* the real spacing, its rms error inside 3σ is
   0.15 (3D table), 0.16 (solar table) and 0.03 (θ13). The real error is about 4× smaller, so
   below the surrogate error above. PCHIP is about 4× better still. Switching is an option, but
   it is not the limiting approximation.
4. **Without mode jumps, PT cannot sample the orderings.** The support is disconnected in
   Δm²3ℓ, so swaps can only permute the K states, and each ordering's population is fixed at
   initialization. The no-jump run gives P(IO) = 0.30–0.68 per replica, against 0.47 exactly.
   The binary indicator's R-hat (1.05) only weakly flags this. The jump moves fix it (pull −0.3σ).

## Statistical approach (assessment)

The proposed procedure is sound. Draw θ from the NuFIT posterior, simulate MINOS+
pseudo-data at θ, and compare the data's test statistic with the pseudo-data distribution.
This is a prior-predictive (Box 1980), or Cousins–Highland hybrid, p-value: "if 3-flavour
oscillation is true, with parameters as the world data say, how often is a result at least
this discrepant?" Requirements:

* **Use one fixed statistic for data and pseudo-data.** Two good choices:
  * a constrained-fit deviance:
    `T = min_θ[−2 ln L_MINOS+(x|θ) + Δχ²_NuFIT(θ)] + 2 ln L_sat(x)`;
  * a prior-predictive deviance: `T = −2 ln[(1/N) Σ_j L(x|θ_j)] + 2 ln L_sat(x)`.

  Do not use the likelihood at the generating θ_i for pseudo-data. Subtracting the saturated
  Poisson likelihood keeps the normalization terms from dominating T.
* **Fluctuate the systematics in the pseudo-experiments**, and include them in T.
* **Choose before unblinding:** an omnibus GoF statistic, or a one-sided excess statistic. The
  excess statistic is still model-independent but more powerful for an excess.
* **Data overlap:** NuFIT 6.0 contains MINOS νμ disappearance and νe appearance
  (10.6e20 POT, 5 + 5 bins). A MINOS+-era sample is independent of these. Any overlap makes
  the test slightly conservative.
* **Pseudo-experiment count:** about 10⁵ independent θ draws gives a 10% error on p ≈ 10⁻³.
* **Robustness:** report the result for both NuFIT variants and both surrogates.

## Review of the original implementation

| where | problem | effect |
|---|---|---|
| `libra.cpp` `MH`, `MonteCarlo` | acceptance uses `exp(−Δχ²)` instead of `exp(−Δχ²/2)` | posterior √2 too narrow: σ(s13) 0.00035, true 0.00060 |
| `libra.cpp` `MonteCarlo` | "swaps" between equal-temperature chains, accepted on a χ² ratio, with only chain 0 recorded and `chis[]` not swapped | not a valid MCMC; biased toward low χ² |
| (structural) | local moves cannot change the mass ordering | ordering weights set by the swap heuristic: P(IO) 0.65, true 0.47; P(UO) 0.93, true 0.71 |
| `NuFuncs.cpp` NO Δm² grid | hard-coded grid rounded to 2 decimals (2.2, 2.2, 2.21, ...) | half resolution, misaligned cells |
| `NuFuncs.cpp` `threeDI` | δ ∈ (175°, 180°) interpolated between 175° and −180° using a width of −355° | wrong values and a discontinuity at 180° |
| `NuFuncs.cpp` `threeDI` | `lower_bound` result dereferenced at `end()` | undefined behaviour at upper edges |
| `libra.cpp` IO branch | snaps `dm231S` against the NO window instead of `dm232S` | ineffective edge handling |
| `NuFuncs.cpp` `twoDI` | `TGraph2D::Interpolate` (Delaunay) returns 0 (= best fit) outside the hull | latent bug; slow |
| `test.cpp` | "Δm²32" histogram filled with Δm²21 − Δm²31 = −Δm²32 | the IO plot actually shows NO samples |

The trilinear weights themselves were correct.

## Layout

```
include/, src/   GridTable (N-linear tables), NuFitModel (parser + surrogate), DirectSampler,
                 PTSampler (+ Target interface), NuFitTarget (jump moves), Output (npy/json/ROOT)
apps/            sample_direct, sample_pt
tests/           test_core (no data needed), test_nufit (skipped without data)
python/nufitval  tables, surrogate (exact profiles / quadrature marginals by variable
                 elimination), diagnostics (R-hat, ESS, ACF), report
python/convert_legacy.py   student's ROOT output -> npy (needs PyROOT)
```
