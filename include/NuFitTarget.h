#pragma once
// PT target for the NuFIT surrogate posterior, with mode-jump moves.
//
// The support is disconnected in dm3l (|dm3l| >= 0.2e-3 eV^2 for both orderings), so no
// local move - at any temperature - can change the mass ordering. The involutions below
// (each its own inverse, unit Jacobian) connect the known degenerate modes:
//   0 "mo":        dm3l -> -dm3l                     (NO Dm31 <-> IO Dm32 of equal size)
//   1 "octant":    s23 -> 1 - s23
//   2 "dcp":       dcp -> 180 - dcp                  (the sin(dcp) - preserving reflection)
//   3 "mo+octant": both 0 and 1

#include "NuFitModel.h"
#include "PTSampler.h"

namespace nufit {

class NuFitTarget : public Target {
public:
  explicit NuFitTarget(const NuFitModel& model) : model_(model) {}
  int dim() const override { return kNPar; }
  std::vector<std::string> names() const override;
  void evaluate(const double* x, double& logL, double& logPrior) const override;
  double period(int d) const override { return d == kDCP ? 360.0 : 0.0; }
  double periodLo(int d) const override { return d == kDCP ? -180.0 : 0.0; }
  bool signFolded(int d) const override { return d == kDM3L; }
  int nJumps() const override { return 4; }
  std::string jumpName(int k) const override;
  void jump(int k, const double* x, double* xp) const override;
  void initialScale(double* s) const override;
  void initialPoint(std::mt19937_64& rng, double* x) const override;

private:
  const NuFitModel& model_;
};

} // namespace nufit
