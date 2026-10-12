#pragma once
// Exact i.i.d. sampler for the NuFIT surrogate posterior.
//
// The surrogate density is a product of "blocks", each of the form
//     f_b(x) = exp(-(chi2_b(x) - off_b)/2) * prod_d pi_d(x_d),
// with chi2_b N-linear on a grid and each per-axis prior pi_d monotone or convex.
// On any grid cell c both factors are maximal at a cell corner, which gives a strict
// envelope M_c >= f_b on c. Sampling a cell with probability M_c V_c / W_b, a point
// uniformly inside it, and accepting with prob. prod_b f_b(x_b)/M_{c_b} is an exact
// rejection sampler for the product (the mass ordering is chosen with probability
// proportional to prod_b W_b, which completes the envelope).
//
// A block may also be *conditional*: chi2_b depends on axes that are sampled by an
// earlier block (s23 in the theta13-coupled surrogate). The envelope of such a block
// takes the maximum over all nodes of the conditioning axes, which remains a valid
// bound because N-linear functions are extremal at nodes.

#include "NuFitModel.h"

#include <cstdint>
#include <functional>
#include <random>
#include <vector>

namespace nufit {

class GridBlock {
public:
  // sampled[d]: true if axis d is sampled by this block, false if it is supplied (conditioned).
  // offset: constant subtracted from chi2; offsetTable (optional, 1D on the single
  // conditioning axis) is subtracted in addition.
  // logPrior[d]: per-axis log prior density for sampled axes (ignored otherwise).
  GridBlock(const GridTable* table, std::vector<bool> sampled, double offset,
            const GridTable* offsetTable, std::vector<std::function<double(double)>> logPrior);

  double totalMass() const { return cumMass_.empty() ? 0 : cumMass_.back(); }
  // Draw a proposal for the sampled axes into x (conditioned axes must already be set);
  // returns the acceptance probability f(x)/M_c.
  double propose(std::mt19937_64& rng, double* x) const;
  // Log of f_b at a full point x (all axes set); used by tests and validation.
  double logDensity(const double* x) const;

private:
  const GridTable* table_;
  std::vector<bool> sampled_;
  double offset_;
  const GridTable* offsetTable_;
  std::vector<std::function<double(double)>> logPrior_;
  int condAxis_ = -1;
  std::vector<int> sAxes_;                    // sampled axes
  std::vector<double> cumMass_;               // cumulative M_c V_c over sampled-axis cells
  std::vector<double> logEnv_;                // log M_c per cell
  std::vector<int> cellShape_;                // cells per sampled axis
};

class DirectSampler {
public:
  explicit DirectSampler(const NuFitModel& model);

  // Draw one exact sample of the parameter vector p[kNPar].
  void sample(std::mt19937_64& rng, double* p);

  uint64_t nProposed() const { return nProposed_; }
  uint64_t nAccepted() const { return nAccepted_; }
  // Envelope-mass fraction of IO (not the posterior P(IO); that is the accepted fraction).
  double envelopeFractionIO() const { return orderingCut_ < 1 ? 1 - orderingCut_ : 0; }

private:
  const NuFitModel& model_;
  // per ordering: blocks over (s23,dm3l,dcp), (s12,l21), (s13) or (s13|s23)
  std::vector<GridBlock> b3_, b2_, b1_;
  double orderingCut_ = 0; // P(propose NO)
  uint64_t nProposed_ = 0, nAccepted_ = 0;
};

} // namespace nufit
