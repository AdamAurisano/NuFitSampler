#pragma once
// NuFIT 6.0 Delta chi^2 tables and the surrogate 6D likelihood built from them.
//
// Parameter vector (index enum Par):
//   s12 = sin^2(theta12), s13 = sin^2(theta13), s23 = sin^2(theta23),
//   dcp  [deg, any real value; wrapped into [-180,180]],
//   dm3l [1e-3 eV^2]: Delta m^2_31 > 0 for NO, Delta m^2_32 < 0 for IO (NuFIT convention),
//   dm21 [1e-5 eV^2].
// The mass ordering is sign(dm3l).
//
// Surrogate Delta chi^2 (relative to the global minimum over both orderings, m_O = minimum
// of ordering O, which every NuFIT table of that ordering carries as an offset):
//   factorized: chi2_3D(s23,dm3l,dcp) + [chi2_2D(s12,log10 dm21) - m_O] + [chi2_1D(s13) - m_O]
//   coupled:    chi2_3D(s23,dm3l,dcp) + [chi2_2D(s12,log10 dm21) - m_O]
//               + [chi2_{T13/T23}(s13,s23) - chi2_{T23}(s23)]
// The coupled form keeps the theta13-theta23 correlation of the global fit; both are
// approximations to the (unpublished) 6D chi^2.

#include "GridTable.h"

#include <map>
#include <string>

namespace nufit {

enum Par { kS12 = 0, kS13, kS23, kDCP, kDM3L, kDM21, kNPar };
enum Ordering { kNO = 0, kIO = 1 };

const char* parName(int p);

// All sections of one NuFIT release file, keyed by tag ("T23/DMA/DCP", "T13", ...).
struct NuFitFile {
  std::string header; // e.g. "Normal Ordering"
  std::map<std::string, GridTable> tables;
  static NuFitFile read(const std::string& path);
  const GridTable& at(const std::string& tag) const;
};

enum class AnglePrior { kSin2, kAngle }; // flat in sin^2(theta) or in theta
enum class Dm21Prior { kLinear, kLog };   // flat in dm21 or in log(dm21)

struct ModelConfig {
  std::string dataDir = "data/nufit";
  std::string variant = "TBoff"; // TBoff = IC19 w/o SK-atm, TByes = IC24 with SK-atm
  bool coupledT13 = false;
  AnglePrior anglePrior = AnglePrior::kSin2;
  Dm21Prior dm21Prior = Dm21Prior::kLinear;
  std::string describe() const;
};

class NuFitModel {
public:
  explicit NuFitModel(const ModelConfig& cfg);

  const ModelConfig& config() const { return cfg_; }

  // Surrogate Delta chi^2; NaN outside the support (the table ranges).
  double deltaChi2(const double* p) const;
  // Log prior density in the sampled coordinates (unnormalized, identical
  // normalization for both orderings); -inf outside the support.
  double logPrior(const double* p) const;
  bool inSupport(const double* p) const;

  // Support: lower/upper bound of each parameter (for dm3l, of |dm3l|).
  double lo(int par) const { return lo_[par]; }
  double hi(int par) const { return hi_[par]; }

  static double wrapDeg(double dcp); // into [-180, 180)
  static Ordering ordering(const double* p) { return p[kDM3L] > 0 ? kNO : kIO; }

  // Table access for the direct sampler and validation.
  const GridTable& t3(Ordering o) const { return t3_[o]; }   // (s23, dm3l, dcp)
  const GridTable& t2(Ordering o) const { return t2_[o]; }   // (s12, log10 dm21/eV^2)
  const GridTable& t1(Ordering o) const { return t1_[o]; }   // (s13)
  const GridTable& t1323(Ordering o) const { return t1323_[o]; } // (s13, s23)
  const GridTable& t23(Ordering o) const { return t23_[o]; }     // (s23)
  double chi2Min(Ordering o) const { return min_[o]; }

  // Log prior density of a single coordinate, expressed in *table* coordinates
  // (s for mixing angles, log10(dm21/eV^2) for dm21, identity otherwise). Each is
  // monotone or convex on any interval, so its maximum over a cell is at an end point.
  double logPriorAngle(double s) const;
  double logPriorLog10Dm21(double l) const;

private:
  ModelConfig cfg_;
  GridTable t3_[2], t2_[2], t1_[2], t1323_[2], t23_[2];
  double min_[2];
  double lo_[kNPar], hi_[kNPar];
};

} // namespace nufit
