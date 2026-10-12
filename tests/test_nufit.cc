// Tests against the NuFIT 6.0 tables (skipped, exit 77, when data/nufit is absent).

#include "DirectSampler.h"
#include "NuFitModel.h"

#include <cmath>
#include <fstream>
#include <iostream>

using namespace nufit;

static int failures = 0;
#define CHECK_NEAR(a, b, tol)                                                                \
  do {                                                                                       \
    const double va = (a), vb = (b);                                                         \
    if (!(std::abs(va - vb) <= (tol))) {                                                     \
      std::cerr << __FILE__ << ":" << __LINE__ << ": " #a " = " << va << " vs " #b " = " << vb \
                << "\n";                                                                     \
      ++failures;                                                                            \
    }                                                                                        \
  } while (0)

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : "data/nufit";
  if (!std::ifstream(dir + "/v60.release-TBoff-NO.txt")) {
    std::cout << "NuFIT tables not found in " << dir << "; skipping\n";
    return 77;
  }
  for (const std::string variant : {"TBoff", "TByes"}) {
    ModelConfig cfg;
    cfg.dataDir = dir;
    cfg.variant = variant;
    const NuFitModel m(cfg);
    // minima found independently when inspecting the release files
    if (variant == "TBoff") {
      CHECK_NEAR(m.chi2Min(kNO), 0.591255, 1e-6);
      CHECK_NEAR(m.chi2Min(kIO), 0.0, 1e-9);
      const double bfIO[kNPar] = {0.31, 0.0222, 0.56, -75, -2.51, 7.4989420933245572}; // log10 = -4.125
      CHECK_NEAR(m.deltaChi2(bfIO), 0.0, 1e-6);
    } else {
      CHECK_NEAR(m.chi2Min(kNO), 0.0, 1e-9);
      CHECK_NEAR(m.chi2Min(kIO), 6.075067, 1e-6);
      const double bfNO[kNPar] = {0.31, 0.0221, 0.47, -150, 2.515, 7.4989420933245572};
      CHECK_NEAR(m.deltaChi2(bfNO), 0.0, 1e-6);
    }
    // delta_CP periodicity: continuous across +-180 and invariant under +-360
    for (double s23 : {0.45, 0.5523, 0.6}) {
      const double a[kNPar] = {0.3, 0.022, s23, 179.999999, 2.51, 7.4};
      const double b[kNPar] = {0.3, 0.022, s23, -179.999999, 2.51, 7.4};
      const double c[kNPar] = {0.3, 0.022, s23, 179.999999 + 720, 2.51, 7.4};
      CHECK_NEAR(m.deltaChi2(a), m.deltaChi2(b), 1e-4);
      CHECK_NEAR(m.deltaChi2(a), m.deltaChi2(c), 1e-9);
    }
    // outside support
    const double gap[kNPar] = {0.3, 0.022, 0.5, 0, 0.1, 7.4};
    if (!std::isnan(m.deltaChi2(gap)) || std::isfinite(m.logPrior(gap))) {
      std::cerr << "point in the dm3l gap not rejected\n";
      ++failures;
    }
    // the exact sampler runs and gives a sane IO fraction
    for (bool coupled : {false, true}) {
      cfg.coupledT13 = coupled;
      const NuFitModel mc(cfg);
      DirectSampler ds(mc);
      std::mt19937_64 rng(5);
      double p[kNPar];
      int nIO = 0;
      const int N = 20000;
      for (int i = 0; i < N; ++i) {
        ds.sample(rng, p);
        nIO += p[kDM3L] < 0;
        if (!std::isfinite(mc.deltaChi2(p))) {
          std::cerr << "direct sample outside support\n";
          ++failures;
          break;
        }
      }
      std::cout << variant << (coupled ? " coupled" : " factorized") << ": acceptance "
                << double(ds.nAccepted()) / ds.nProposed() << ", P(IO) " << double(nIO) / N << "\n";
      if (variant == "TByes" && nIO > 0.2 * N) {
        std::cerr << "TByes should strongly prefer NO\n";
        ++failures;
      }
    }
  }
  if (failures) {
    std::cerr << failures << " failure(s)\n";
    return 1;
  }
  std::cout << "all NuFIT tests passed\n";
  return 0;
}
