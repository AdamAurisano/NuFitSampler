#pragma once
// Parallel-tempering MCMC.
//
// Rung k targets pi_k(x) ~ prior(x) * L(x)^beta_k with beta_0 = 1 > beta_1 > ... > beta_{K-1}
// (only the likelihood is tempered). Each iteration:
//   1. every rung makes one adaptive random-walk Metropolis move and, with probability
//      pJump, one deterministic involutive "mode jump" (unit Jacobian, so the MH ratio is
//      the density ratio), e.g. the mass-ordering mirror dm3l -> -dm3l;
//   2. adjacent rungs attempt state swaps, alternating even and odd pairs (the
//      non-reversible DEO scheme of Syed et al. 2021), accepted with
//      min(1, exp[(beta_i - beta_j)(logL_j - logL_i)]).
// During burn-in (only) the random-walk covariances (Haario et al. 2001 with Robbins-Monro
// scale targeting acceptance 0.234) and the temperature ladder (Vousden, Farr & Mandel 2016,
// equalizing adjacent swap rates with fixed T_max) are adapted. They are then frozen, so
// the post-burn-in chain is a time-homogeneous Markov chain with the correct stationary law.
// Several independent ensembles ("replicas") run in parallel threads for R-hat.

#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace nufit {

class Target {
public:
  virtual ~Target() = default;
  virtual int dim() const = 0;
  virtual std::vector<std::string> names() const = 0;
  // logPrior = -inf outside the support (logL is then not used).
  virtual void evaluate(const double* x, double& logL, double& logPrior) const = 0;
  // Periodic coordinates: period > 0 (the RW step is wrapped into [lo, lo+period)).
  virtual double period(int /*d*/) const { return 0; }
  virtual double periodLo(int /*d*/) const { return 0; }
  // Sign-folded coordinates: the RW acts on |x| and keeps the sign (moves that would
  // cross zero are rejected); the covariance is learned on |x|.
  virtual bool signFolded(int /*d*/) const { return false; }
  // Involutive jump moves with unit Jacobian.
  virtual int nJumps() const { return 0; }
  virtual std::string jumpName(int /*k*/) const { return ""; }
  virtual void jump(int /*k*/, const double* x, double* xp) const {
    for (int d = 0; d < dim(); ++d) xp[d] = x[d];
  }
  virtual void initialScale(double* s) const = 0;
  virtual void initialPoint(std::mt19937_64& rng, double* x) const = 0;
};

struct PTConfig {
  int nReplicas = 4;
  int nRungs = 8;
  double tMax = 50;
  uint64_t nBurn = 50000;
  uint64_t nSteps = 500000;
  int thin = 10;
  double pJump = 0.1;
  bool useJumps = true;
  bool adaptLadder = true;
  int adaptInterval = 1000; // burn-in iterations between covariance updates
  uint64_t seed = 12345;
};

struct PTReplicaResult {
  // thinned cold-chain samples: per row x[dim], logL, logPrior
  std::vector<double> cold;
  // thinned per-rung logL: per row K values
  std::vector<double> rungLogL;
  // thinned rung position of each label (state copy): per row K values
  std::vector<double> labelRung;
  // ladder snapshots during burn-in: per row iteration, beta[K]
  std::vector<double> ladderHistory;
  std::vector<double> betas;                         // final ladder
  std::vector<double> rwAccept;                      // per rung (post burn-in)
  std::vector<std::vector<double>> jumpAccept;       // [rung][jump]
  std::vector<double> swapAccept;                    // per adjacent pair (post burn-in)
  std::vector<double> rwScale;                       // final RW scale factor per rung
  uint64_t roundTrips = 0;
  double seconds = 0;
};

class PTSampler {
public:
  PTSampler(const Target& target, const PTConfig& cfg) : target_(target), cfg_(cfg) {}
  std::vector<PTReplicaResult> run() const;
  PTReplicaResult runReplica(int replica) const;

private:
  const Target& target_;
  PTConfig cfg_;
};

} // namespace nufit
