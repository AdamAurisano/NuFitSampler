#include "NuFitTarget.h"

#include <cmath>
#include <limits>

namespace nufit {

std::vector<std::string> NuFitTarget::names() const {
  std::vector<std::string> v;
  for (int p = 0; p < kNPar; ++p) v.push_back(parName(p));
  return v;
}

void NuFitTarget::evaluate(const double* x, double& logL, double& logPrior) const {
  logPrior = model_.logPrior(x);
  if (!std::isfinite(logPrior)) {
    logPrior = -std::numeric_limits<double>::infinity();
    logL = 0;
    return;
  }
  logL = -0.5 * model_.deltaChi2(x);
}

std::string NuFitTarget::jumpName(int k) const {
  static const char* n[] = {"mo", "octant", "dcp", "mo+octant"};
  return n[k];
}

void NuFitTarget::jump(int k, const double* x, double* xp) const {
  for (int d = 0; d < kNPar; ++d) xp[d] = x[d];
  if (k == 0 || k == 3) xp[kDM3L] = -x[kDM3L];
  if (k == 1 || k == 3) xp[kS23] = 1.0 - x[kS23];
  if (k == 2) xp[kDCP] = NuFitModel::wrapDeg(180.0 - x[kDCP]);
}

void NuFitTarget::initialScale(double* s) const {
  s[kS12] = 0.012;
  s[kS13] = 0.0006;
  s[kS23] = 0.02;
  s[kDCP] = 30;
  s[kDM3L] = 0.03;
  s[kDM21] = 0.2;
}

void NuFitTarget::initialPoint(std::mt19937_64& rng, double* x) const {
  // Over-dispersed relative to the posterior, random ordering and octant.
  std::uniform_real_distribution<double> u(0.0, 1.0);
  x[kS12] = 0.26 + 0.10 * u(rng);
  x[kS13] = 0.019 + 0.006 * u(rng);
  x[kS23] = 0.38 + 0.26 * u(rng);
  x[kDCP] = -180 + 360 * u(rng);
  x[kDM3L] = (u(rng) < 0.5 ? 1 : -1) * (2.3 + 0.4 * u(rng));
  x[kDM21] = 6.5 + 2.0 * u(rng);
}

} // namespace nufit
