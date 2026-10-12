#include "PTSampler.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>

namespace nufit {

namespace {

constexpr double kTargetAccept = 0.234;

// Lower-triangular Cholesky factor of a (n x n, row-major); false if not positive definite.
bool cholesky(const std::vector<double>& a, int n, std::vector<double>& l) {
  l.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = a[i * n + j];
      for (int k = 0; k < j; ++k) s -= l[i * n + k] * l[j * n + k];
      if (i == j) {
        if (!(s > 0)) return false;
        l[i * n + i] = std::sqrt(s);
      } else {
        l[i * n + j] = s / l[j * n + j];
      }
    }
  }
  return true;
}

double wrapInto(double x, double lo, double period) {
  double w = std::fmod(x - lo, period);
  if (w < 0) w += period;
  return lo + w;
}

struct Rung {
  std::vector<double> x;
  double logL = 0, logPrior = 0;
};

} // namespace

std::vector<PTReplicaResult> PTSampler::run() const {
  std::vector<PTReplicaResult> out(cfg_.nReplicas);
  std::vector<std::thread> threads;
  for (int r = 0; r < cfg_.nReplicas; ++r)
    threads.emplace_back([this, r, &out]() { out[r] = runReplica(r); });
  for (auto& t : threads) t.join();
  return out;
}

PTReplicaResult PTSampler::runReplica(int replica) const {
  const auto t0 = std::chrono::steady_clock::now();
  const int n = target_.dim();
  const int K = cfg_.nRungs;
  const int nJ = cfg_.useJumps ? target_.nJumps() : 0;
  std::seed_seq seq{static_cast<uint64_t>(cfg_.seed), static_cast<uint64_t>(replica), uint64_t{0x5eed}};
  std::mt19937_64 rng(seq);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::normal_distribution<double> gauss(0.0, 1.0);

  // Temperature ladder: log T equally spaced in [0, log Tmax].
  std::vector<double> gaps(std::max(K - 1, 0), K > 1 ? std::log(cfg_.tMax) / (K - 1) : 0.0);
  std::vector<double> beta(K, 1.0);
  auto setBetas = [&]() {
    double logT = 0;
    beta[0] = 1;
    for (int k = 1; k < K; ++k) {
      logT += gaps[k - 1];
      beta[k] = std::exp(-logT);
    }
  };
  setBetas();

  // Per-rung random-walk proposal: step = exp(logScale) * chol * z.
  std::vector<double> scale0(n);
  target_.initialScale(scale0.data());
  std::vector<std::vector<double>> chol(K, std::vector<double>(n * n, 0.0));
  std::vector<double> logScale(K, 0.0);
  for (int k = 0; k < K; ++k)
    for (int d = 0; d < n; ++d) chol[k][d * n + d] = scale0[d] * std::sqrt(1.0 / beta[k]);

  // Initial states.
  std::vector<Rung> rung(K);
  for (auto& s : rung) {
    s.x.resize(n);
    for (int tries = 0;; ++tries) {
      target_.initialPoint(rng, s.x.data());
      target_.evaluate(s.x.data(), s.logL, s.logPrior);
      if (std::isfinite(s.logPrior) && std::isfinite(s.logL)) break;
      if (tries > 10000) throw std::runtime_error("PTSampler: cannot find an initial point in the support");
    }
  }
  std::vector<int> rungOfLabel(K), labelOfRung(K), dir(K, 0);
  for (int k = 0; k < K; ++k) rungOfLabel[k] = labelOfRung[k] = k;

  PTReplicaResult res;
  std::vector<uint64_t> rwTry(K, 0), rwAcc(K, 0), swTry(std::max(K - 1, 1), 0), swAcc(std::max(K - 1, 1), 0);
  std::vector<std::vector<uint64_t>> jTry(K, std::vector<uint64_t>(std::max(nJ, 1), 0)),
      jAcc(K, std::vector<uint64_t>(std::max(nJ, 1), 0));
  std::vector<double> swapRateEma(std::max(K - 1, 1), 0.5);

  // Burn-in history in covariance coordinates, per rung.
  std::vector<std::vector<double>> hist(K);
  for (auto& h : hist) h.reserve(static_cast<size_t>(cfg_.nBurn) * n);

  std::vector<double> xp(n), step(n), z(n), u(n);

  auto mh = [&](Rung& s, double b, const std::vector<double>& prop) {
    double lL, lP;
    target_.evaluate(prop.data(), lL, lP);
    if (!std::isfinite(lP) || !std::isfinite(lL)) return false;
    const double logA = b * (lL - s.logL) + (lP - s.logPrior);
    if (logA >= 0 || std::log(u01(rng)) < logA) {
      s.x = prop;
      s.logL = lL;
      s.logPrior = lP;
      return true;
    }
    return false;
  };

  auto updateCovariance = [&](int k, uint64_t t) {
    // samples from the second half of the burn-in so far
    const size_t rows = hist[k].size() / n;
    const size_t first = rows / 2;
    const size_t m = rows - first;
    if (m < static_cast<size_t>(10 * n)) return;
    // circular centring of periodic coordinates
    std::vector<double> centre(n, 0.0);
    for (int d = 0; d < n; ++d) {
      const double P = target_.period(d);
      if (P <= 0) continue;
      double sn = 0, cs = 0;
      for (size_t r = first; r < rows; ++r) {
        const double a = 2 * M_PI * hist[k][r * n + d] / P;
        sn += std::sin(a);
        cs += std::cos(a);
      }
      centre[d] = std::atan2(sn, cs) * P / (2 * M_PI);
    }
    std::vector<double> mean(n, 0.0), cov(n * n, 0.0);
    auto coord = [&](size_t r, int d) {
      double v = hist[k][r * n + d];
      const double P = target_.period(d);
      if (P > 0) v = wrapInto(v - centre[d], -0.5 * P, P);
      return v;
    };
    for (size_t r = first; r < rows; ++r)
      for (int d = 0; d < n; ++d) mean[d] += coord(r, d) / m;
    for (size_t r = first; r < rows; ++r)
      for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j)
          cov[i * n + j] += (coord(r, i) - mean[i]) * (coord(r, j) - mean[j]) / (m - 1);
    const double sd = 2.38 * 2.38 / n;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j <= i; ++j) {
        cov[i * n + j] *= sd;
        cov[j * n + i] = cov[i * n + j];
      }
    for (int i = 0; i < n; ++i) cov[i * n + i] += 1e-8 * scale0[i] * scale0[i];
    std::vector<double> l;
    if (cholesky(cov, n, l)) chol[k] = l;
    (void)t;
  };

  const uint64_t nTot = cfg_.nBurn + cfg_.nSteps;
  for (uint64_t t = 0; t < nTot; ++t) {
    const bool burn = t < cfg_.nBurn;

    // 1. local moves
    for (int k = 0; k < K; ++k) {
      Rung& s = rung[k];
      for (int d = 0; d < n; ++d) z[d] = gauss(rng);
      const double sc = std::exp(logScale[k]);
      bool valid = true;
      for (int i = 0; i < n; ++i) {
        double st = 0;
        for (int j = 0; j <= i; ++j) st += chol[k][i * n + j] * z[j];
        st *= sc;
        const double P = target_.period(i);
        if (target_.signFolded(i)) {
          const double a = std::abs(s.x[i]) + st;
          if (a <= 0) valid = false;
          xp[i] = std::copysign(a, s.x[i]);
        } else if (P > 0) {
          xp[i] = wrapInto(s.x[i] + st, target_.periodLo(i), P);
        } else {
          xp[i] = s.x[i] + st;
        }
      }
      const bool acc = valid && mh(s, beta[k], xp);
      if (burn) {
        const double g = std::pow(1.0 + t / 100.0, -0.6);
        logScale[k] += g * ((acc ? 1.0 : 0.0) - kTargetAccept);
      } else {
        ++rwTry[k];
        rwAcc[k] += acc;
      }
      if (nJ > 0 && u01(rng) < cfg_.pJump) {
        const int j = std::min(static_cast<int>(u01(rng) * nJ), nJ - 1);
        target_.jump(j, s.x.data(), xp.data());
        const bool ja = mh(s, beta[k], xp);
        if (!burn) {
          ++jTry[k][j];
          jAcc[k][j] += ja;
        }
      }
      if (burn) {
        for (int d = 0; d < n; ++d) hist[k].push_back(target_.signFolded(d) ? std::abs(s.x[d]) : s.x[d]);
      }
    }

    // 2. swaps between adjacent rungs, alternating even / odd pairs (DEO)
    for (int i = static_cast<int>(t % 2); i + 1 < K; i += 2) {
      const double logA = (beta[i] - beta[i + 1]) * (rung[i + 1].logL - rung[i].logL);
      const bool acc = logA >= 0 || std::log(u01(rng)) < logA;
      if (acc) {
        std::swap(rung[i], rung[i + 1]);
        std::swap(labelOfRung[i], labelOfRung[i + 1]);
        rungOfLabel[labelOfRung[i]] = i;
        rungOfLabel[labelOfRung[i + 1]] = i + 1;
      }
      if (burn) {
        swapRateEma[i] += 0.01 * ((acc ? 1.0 : 0.0) - swapRateEma[i]);
      } else {
        ++swTry[i];
        swAcc[i] += acc;
      }
    }
    if (K > 1 && !burn) {
      const int lb = labelOfRung[0], lt = labelOfRung[K - 1];
      if (dir[lb] == -1) ++res.roundTrips;
      dir[lb] = 1;
      if (dir[lt] == 1) dir[lt] = -1;
    }

    // 3. adaptation (burn-in only)
    if (burn) {
      if (cfg_.adaptLadder && K > 2 && t % 10 == 9) {
        const double kappa = 0.2 * 1e4 / (t + 1e4);
        double meanRate = 0;
        for (int i = 0; i < K - 1; ++i) meanRate += swapRateEma[i] / (K - 1);
        double sum = 0;
        for (int i = 0; i < K - 1; ++i) {
          gaps[i] *= std::exp(kappa * (swapRateEma[i] - meanRate));
          sum += gaps[i];
        }
        for (auto& g : gaps) g *= std::log(cfg_.tMax) / sum;
        setBetas();
      }
      if ((t + 1) % cfg_.adaptInterval == 0) {
        for (int k = 0; k < K; ++k) updateCovariance(k, t);
        res.ladderHistory.push_back(static_cast<double>(t + 1));
        res.ladderHistory.insert(res.ladderHistory.end(), beta.begin(), beta.end());
      }
      if (t + 1 == cfg_.nBurn)
        for (auto& h : hist) std::vector<double>().swap(h);
    } else if ((t - cfg_.nBurn) % cfg_.thin == 0) {
      res.cold.insert(res.cold.end(), rung[0].x.begin(), rung[0].x.end());
      res.cold.push_back(rung[0].logL);
      res.cold.push_back(rung[0].logPrior);
      for (int k = 0; k < K; ++k) res.rungLogL.push_back(rung[k].logL);
      for (int l = 0; l < K; ++l) res.labelRung.push_back(rungOfLabel[l]);
    }
  }

  res.betas = beta;
  for (int k = 0; k < K; ++k) {
    res.rwAccept.push_back(rwTry[k] ? double(rwAcc[k]) / rwTry[k] : 0);
    res.rwScale.push_back(std::exp(logScale[k]));
    std::vector<double> ja;
    for (int j = 0; j < nJ; ++j) ja.push_back(jTry[k][j] ? double(jAcc[k][j]) / jTry[k][j] : 0);
    res.jumpAccept.push_back(ja);
  }
  for (int i = 0; i + 1 < K; ++i) res.swapAccept.push_back(swTry[i] ? double(swAcc[i]) / swTry[i] : 0);
  res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return res;
}

} // namespace nufit
